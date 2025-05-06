""" Eiger 1M setup """

from ophyd import (
    ADComponent, Staged, SignalRO, DynamicDeviceComponent
)
from ophyd.mca import EpicsMCARecord
from ophyd.areadetector import DetectorBase
from ophyd.areadetector.trigger_mixins import ADTriggerStatus
from pathlib import Path
from collections import OrderedDict
from time import time as ttime
from .ad_mixins import TriggerBase
from .vortex_dante_parts import DanteCAM, DanteHDF1Plugin, DanteSCA
from ..utils.config import iconfig
from ..utils._logging_setup import logger
logger.info(__file__)

__all__ = ["vortex"]

# Bluesky and IOC have the same path root.
# IOC_FILES_ROOT = Path(iconfig["AREA_DETECTOR"]["VORTEX"]["IOC_FILES_ROOT"])
IOC_FILES_ROOT = Path("")

DEFAULT_FOLDER = Path(iconfig["AREA_DETECTOR"]["VORTEX"]["DEFAULT_FOLDER"])

HDF1_NAME_TEMPLATE = iconfig["AREA_DETECTOR"]["HDF5_FILE_TEMPLATE"]
HDF1_FILE_EXTENSION = iconfig["AREA_DETECTOR"]["HDF5_FILE_EXTENSION"]
HDF1_NAME_FORMAT = HDF1_NAME_TEMPLATE + "." + HDF1_FILE_EXTENSION

MAX_TIME = 60*60  # time used in align mode
MAX_ROIS = 32


class Trigger(TriggerBase):
    """
    This trigger mixin class takes one acquisition per trigger.
    """
    _status_type = ADTriggerStatus
    # _status_type = DeviceStatus

    def __init__(self, *args, image_name=None, **kwargs):
        super().__init__(
            *args,
            acquisition_signal_dev="cam.acquire_start",
            acquire_busy_signal_dev="cam.acquire_busy",
            **kwargs
        )

        if image_name is None:
            image_name = '_'.join([self.name, 'image'])
        self._image_name = image_name
        # self._acquisition_signal = self.cam.acquire_start
        self._acquisition_signal_stop = self.cam.acquire_stop
        # self._acquire_busy_signal = self.cam.acquire_busy

        # self._acquisition_signal_pv = "cam.acquire_start"
        # self._acquire_busy_signal_pv = "cam.acquire_busy"
        self._flysetup = False
        self._status = None

        self.stage_sigs = OrderedDict(
            [
                ("cam.acquire_stop", 1),  # If acquiring, stop
                ("cam.mca_mode", "MCA Mapping"),  # 'Multiple' mode
                ("cam.mca_mapping_points", 1),
            ]
        )

    def setup_manual_trigger(self):
        # Stage signals
        self.cam.stage_sigs["wait_for_plugins"] = "Yes"

    def setup_external_trigger(self):
        # Stage signals
        # self.cam.stage_sigs["trigger_mode"] = "TTL Veto Only"
        # self.cam.stage_sigs["num_images"] = MAX_IMAGES
        # self.cam.stage_sigs["wait_for_plugins"] = "No"
        raise NotImplementedError("This detector cannot be used in flyscans")

    def stage(self):

        if self._flysetup:
            self.setup_external_trigger()

        # Make sure that detector not running
        self._acquire_busy_signal.subscribe(self._acquire_changed)

        super().stage()

        if self._flysetup:
            self._acquisition_signal.set(1).wait(timeout=10)

    def unstage(self):
        super().unstage()
        self._acquisition_signal_stop.set(1).wait(timeout=10)
        self._flysetup = False
        self._acquire_busy_signal.clear_sub(self._acquire_changed)
        self._collect_image = False
        self.setup_manual_trigger()

    def trigger(self):
        if self._staged != Staged.yes:
            raise RuntimeError("This detector is not ready to trigger."
                               "Call the stage() method before triggering.")

        # Click the Acquire_button
        self._status = self._status_type(self)
        self._acquisition_signal.put(1, wait=False)
        if self.hdf1.enable.get() in (True, 1, "on", "Enable"):
            self.generate_datum(self._image_name, ttime(), {})

        return self._status

    def _acquire_changed(self, value=None, old_value=None, **kwargs):
        "This is called when the 'acquire_busy' signal changes."

        if self._status is None:
            return
        if (old_value != 0) and (value == 0):
            # Negative-going edge means an acquisition just finished.
            # sleep(self._delay)
            self._status.set_finished()
            self._status = None


class TotalCorrectedSignal(SignalRO):
    """ Signal that returns the deadtime corrected total counts """

    def __init__(self, prefix, roi_index=0, **kwargs):
        self.roi_index = roi_index
        super().__init__(**kwargs)

    def get(self, **kwargs):
        value = 0
        for ch_num in range(1, self.root._num_channels+1):
            roi = getattr(self.root.mcas, f'mca{ch_num}.rois.roi{self.roi_index}')
            sca = getattr(self.root.scas, f"sca{ch_num}")
            _ocr = sca.ocr.get(**kwargs)
            correction = 1.0 if _ocr == 0 else sca.icr.get(**kwargs)/_ocr
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


def _mcas(num_channels):
    defn = OrderedDict()
    for k in range(1, num_channels+1):
        defn[f'mca{k}'] = (EpicsMCARecord, f'mca{k}', {})
    return defn


def _scas(num_channels):
    defn = OrderedDict()
    for k in range(1, num_channels+1):
        defn[f'sca{k}'] = (DanteSCA, f'dante{k}:', {})
    return defn


class DanteDetector(Trigger, DetectorBase):

    _default_configuration_attrs = ('cam',)
    _default_read_attrs = (
        'hdf1',
        'mcas',
        'scas',
        'total'
    )

    _read_rois = [1]
    _num_channels = 4
    _mca_rois_read_attrs = ("count",)

    cam = ADComponent(DanteCAM, "dante:")

    mcas = DynamicDeviceComponent(_mcas(_num_channels))
    scas = DynamicDeviceComponent(_scas(_num_channels))

    total = DynamicDeviceComponent(_totals('roi', range(MAX_ROIS)))

    hdf1 = ADComponent(
        DanteHDF1Plugin,
        "HDF1:",
        ioc_path_root=IOC_FILES_ROOT,
    )

    # TODO: TEMPORARY!
    _local_folder = "/local/home/dpuser/sector4/"

    # Make this compatible with other detectors
    @property
    def preset_monitor(self):
        return self.cam.real_time_preset

    @property
    def num_channels(self):
        return self._num_channels

    def align_on(self):
        """Start detector in alignment mode"""
        self.save_images_off()
        self.cam.mca_mode.set("MCA").wait(timeout=10)
        self.preset_monitor.set(MAX_TIME).wait(timeout=10)
        self.cam.acquire_start.set(1).wait(timeout=10)

    def align_off(self):
        """Stop detector"""
        self.cam.acquire_stop.set(1).wait(timeout=10)

    def save_images_on(self):
        self.hdf1.enable.set("Enable").wait(timeout=10)

    def save_images_off(self):
        self.hdf1.enable.set("Disable").wait(timeout=10)

    def auto_save_on(self):
        self.hdf1.autosave.put("on")

    def auto_save_off(self):
        self.hdf1.autosave.put("off")

    def default_settings(self):

        self.hdf1.file_template.put(HDF1_NAME_FORMAT)
        self.hdf1.file_path.put(str(DEFAULT_FOLDER))
        self.hdf1.num_capture.put(0)

        if "enable" in self.hdf1.stage_sigs.keys():
            self.hdf1.stage_sigs.pop("enable")
        self.hdf1.stage_sigs["num_capture"] = 0
        self.hdf1.stage_sigs["capture"] = 1

        self.setup_manual_trigger()
        self.save_images_off()
        self.read_rois = [0]
        self.plot_roi0()

        for item in "mca1 mca2 mca3 mca4".split():
            mca = getattr(self.mcas, item)
            mca.preset_real_time.kind = "omitted"
            mca.elapsed_real_time.kind = "omitted"
            mca.spectrum.kind = "omitted"
            for roi in mca.rois.component_names:
                d = getattr(mca, f"rois.{roi}")
                for c in d.component_names:
                    k = (
                        "normal"
                        if c in self._mca_rois_read_attrs
                        else "config"
                    )
                    getattr(d, c).kind = k

    # TODO: Probably need to take another look at read_rois.setter and
    # select_roi
    # @property
    # def read_rois(self):
    #     return self._read_rois

    # @read_rois.setter
    # def read_rois(self, rois):
    #     self._read_rois = list(rois)

    # def select_roi(self, rois):

    #     for i in range(MAX_ROIS):
    #         kh = "hinted" if i in rois else "normal"
    #         getattr(self.total, f"roi{i}").kind = kh

    #         if kh == "hinted" and i not in self.read_rois:
    #             self.read_rois.append(i)

    #         kr = "normal" if i in self.read_rois else "omitted"
    #         getattr(self.total, f"roi{i}").kind = kr
    @property
    def read_rois(self):
        return self._read_rois

    @read_rois.setter
    def read_rois(self, rois):
        # Change total kinds
        for i in range(MAX_ROIS):
            if i in rois:
                ktot = getattr(self.total, f"roi{i}").kind.name
                if ktot == "omitted":
                    getattr(self.total, f"roi{i}").kind = "normal"
            else:
                getattr(self.total, f"roi{i}").kind = "omitted"

        # change ROISTAT kinds
        for pixel in range(1, self.num_channels+1):
            pix = getattr(self.mcas, f"mca{pixel}")
            for i in range(MAX_ROIS):
                k = "normal" if i in rois else "omitted"
                getattr(pix.rois, f"roi{i}").kind = k

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

    def setup_images(
            self, base_folder, file_name_base, file_number, flyscan=False
    ):

        self.hdf1.file_name.set(file_name_base).wait(timeout=10)
        self.hdf1.file_number.set(file_number).wait(timeout=10)
        self.auto_save_on()
        self._flysetup = flyscan

        base_folder = str(base_folder) + f"/{self.name}/"
        # self.hdf1.file_path.set(base_folder).wait(timeout=10)
        # TODO: need to temporarily change the saving folder.
        base_folder2 = self._local_folder
        self.hdf1.file_path.set(base_folder2).wait(timeout=10)

        _, full_path, relative_path = self.hdf1.make_write_read_paths(
            base_folder
        )

        return Path(full_path), Path(relative_path)

    @property
    def save_image_flag(self):
        _hdf1_auto = True if self.hdf1.autosave.get() == "on" else False
        _hdf1_on = True if self.hdf1.enable.get() == "Enable" else False
        return _hdf1_on or _hdf1_auto


vortex = DanteDetector("dp_dante8_xrd4:", name="vortex", labels=("detector",))
