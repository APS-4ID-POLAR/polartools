""" Vortex with DXP """

from ophyd.mca import SaturnDXP, EpicsMCARecord
from ophyd import (
    Staged,
    Device,
    DynamicDeviceComponent,
    Component,
    EpicsSignal,
    EpicsSignalRO,
    EpicsSignalWithRBV,
    SignalRO,
)
from ophyd.status import DeviceStatus
from collections import OrderedDict
from ..utils._logging_setup import logger
logger.info(__file__)

MAX_ROIS = 32


class MyDXP(SaturnDXP):
    live_time_output = None
    trigger_output = None


class MyMCA(EpicsMCARecord):
    check_acquiring = Component(
        EpicsSignal, '.ACQG', kind='omitted', string=False
    )


class SingleTrigger(Device):
    """
    This trigger mixin class takes one acquisition per trigger.
    Examples
    --------
    >>> class SimDetector(SingleTrigger):
    ...     pass
    >>> det = SimDetector('..pv..')
    # optionally, customize name of image
    >>> det = SimDetector('..pv..', image_name='fast_detector_image')
    """

    _status_type = DeviceStatus

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._acquisition_signal = self.erase_start
        self._status_signal = self.status

    def stage(self):
        self._status_signal.subscribe(self._acquire_changed)
        super().stage()

    def unstage(self):
        super().unstage()
        self._status_signal.clear_sub(self._acquire_changed)

    def trigger(self):
        "Trigger one acquisition."
        if self._staged != Staged.yes:
            raise RuntimeError("This detector is not ready to trigger."
                               "Call the stage() method before triggering.")

        self._status = self._status_type(self)
        self._acquisition_signal.put(1, wait=False)
        return self._status

    def _acquire_changed(self, value=None, old_value=None, **kwargs):
        "This is called when the 'acquire' signal changes."
        if self._status is None:
            return
        if (old_value == 1) and (value == 0):
            # Negative-going edge means an acquisition just finished.
            self._status.set_finished()
            self._status = None


class TotalCorrectedSignal(SignalRO):
    """ Signal that returns the deadtime corrected total counts """

    def __init__(self, prefix, roi_index=0, **kwargs):
        self.roi_index = roi_index
        super().__init__(**kwargs)

    def get(self, **kwargs):
        value = 0
        for ch_num in range(1, 4+1):
            roi = getattr(self.root, f'mca{ch_num}.rois.roi{self.roi_index}')
            dxp = getattr(self.root, f"dxp{ch_num}")
            _ocr = dxp.output_count_rate.get(**kwargs)
            correction = 1.0 if _ocr == 0 else dxp.input_count_rate.get(**kwargs)/_ocr
            value += roi.count.get(**kwargs) * correction
        return value


def _totals(attr_fix, id_range):
    defn = OrderedDict()
    for k in id_range:
        _kind = "normal" if k == 0 else "omitted"
        defn['{}{:d}'.format(attr_fix, k)] = (
            TotalCorrectedSignal, '', {'roi_index': k, 'kind': _kind}
        )
    return defn


class MyXMAP(SingleTrigger):

    # Buttons
    start = Component(EpicsSignal, "StartAll", kind="omitted")
    stop_ = Component(EpicsSignal, "StopAll", kind="omitted")
    erase_start = Component(EpicsSignal, "EraseStart", kind="omitted")
    erase = Component(EpicsSignal, "EraseAll", kind="omitted")

    # Status and configs
    status = Component(EpicsSignal, "Acquiring", kind="config")
    collection_mode = Component(EpicsSignal, "CollectMode", kind="config")
    preset_mode = Component(EpicsSignal, "PresetMode", kind="config")
    instant_deadtime = Component(EpicsSignalRO, "IDeadTime", kind="normal")
    average_deadtime = Component(EpicsSignalRO, "DeadTime", kind="normal")
    poll_time = Component(EpicsSignalWithRBV, "PollTime", kind="config")

    # Times
    real_preset = Component(EpicsSignal, "PresetReal", kind="config")
    live_preset = Component(EpicsSignal, "PresetLive", kind="config")
    real_elapsed = Component(EpicsSignal, "ElapsedReal", kind="normal")
    live_elapsed = Component(EpicsSignal, "ElapsedLive", kind="normal")

    events_preset = Component(EpicsSignal, "PresetEvents", kind="config")
    triggers_preset = Component(EpicsSignal, "PresetTriggers", kind="config")

    total = DynamicDeviceComponent(_totals('roi', range(MAX_ROIS)))

    # MCAs
    mca1 = Component(MyMCA, "mca1", kind="config")
    mca2 = Component(MyMCA, "mca2", kind="config")
    mca3 = Component(MyMCA, "mca3", kind="config")
    mca4 = Component(MyMCA, "mca4", kind="config")

    # DXPs
    dxp1 = Component(MyDXP, "dxp1:", kind="config")
    dxp2 = Component(MyDXP, "dxp2:", kind="config")
    dxp3 = Component(MyDXP, "dxp3:", kind="config")
    dxp4 = Component(MyDXP, "dxp4:", kind="config")

    _read_rois = [1]

    @property
    def preset_monitor(self):
        return self.real_preset

    def default_kinds(self):

        # TODO: This is setting A LOT of stuff as "configuration_attrs", should
        # be revised at some point.

        # self.mca1.configuration_attrs += [
        #     item for item in self.mca1.component_names
        # ]

        # self.dxp.configuration_attrs += [
        #     item for item in self.dxp.component_names
        # ]

        self.mca1.read_attrs = [
            "preset_real_time",
            "preset_live_time",
            "elapsed_real_time",
            "elapsed_live_time",
            "rois.roi0",
            "rois.roi1",
        ]

    def default_settings(self):
        self.stage_sigs['stop_'] = 1
        self.stage_sigs['erase'] = 1
        self.stage_sigs['preset_mode'] = "Real time"

    @property
    def read_rois(self):
        return self._read_rois

    @read_rois.setter
    def read_rois(self, rois):
        self._read_rois = list(rois)

    def select_roi(self, rois):

        for i in range(MAX_ROIS):
            k = (
                "hinted" if i in rois else
                "normal" if i in self.read_rois else
                "omitted"
            )

            getattr(self.total, f"roi{i}").kind = k

            if k == "hinted" and i not in self.read_rois:
                self.read_rois.append(i)

    def plot_roi0(self):
        self.select_roi([0])

    def plot_roi1(self):
        self.select_roi([1])

    def plot_roi2(self):
        self.select_roi([2])

    def plot_roi3(self):
        self.select_roi([3])

    def plot_roi4(self):
        self.select_roi([4])

    @property
    def label_option_map(self):
        return {f"ROI{i} Total": i for i in range(1, 8+1)}

    @property
    def plot_options(self):
        # Return all named scaler channels
        return list(self.label_option_map.keys())

    def select_plot(self, channels):
        chans = [self.label_option_map[i] for i in channels]
        self.select_roi(chans)


vortex = MyXMAP("dxpXMAPDP2:", name="vortex")
