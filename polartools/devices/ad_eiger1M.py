""" Eiger 1M setup """

from ophyd import ADComponent, Staged
from ophyd.status import wait as status_wait, SubscriptionStatus
from ophyd.areadetector import DetectorBase
from ophyd.areadetector.trigger_mixins import TriggerBase, ADTriggerStatus
from apstools.utils import run_in_thread
from pathlib import Path
from time import time as ttime, sleep
from .ad_mixins import (
    EigerDetectorCam,
    CodecPlugin,
    ImagePlugin,
    ROIPlugin,
    StatsPlugin,
    PolarHDF5Plugin,
    ProcessPlugin,
    TransformPlugin
)


class TriggerTime(TriggerBase):
    """
    This trigger mixin class takes one acquisition per trigger.
    """
    _status_type = ADTriggerStatus

    def __init__(self, *args, image_name=None, min_period=0.2, **kwargs):
        super().__init__(*args, **kwargs)
        if image_name is None:
            image_name = '_'.join([self.name, 'image'])
        self._image_name = image_name
        self._acquisition_signal_pv = "cam.special_trigger_button"
        self._min_period = min_period
        self._flysetup = False

    @property
    def acquisition_signal(self):
        return getattr(self, self._acquisition_signal_pv)

    @property
    def min_period(self):
        return self._min_period

    @min_period.setter
    def min_period(self, value):
        try:
            self._min_period = float(value)
        except ValueError:
            raise ValueError("min_period must be a number.")

    def setup_manual_trigger(self):
        # Stage signals
        self.cam.stage_sigs["trigger_mode"] = "Continuous"
        self.cam.stage_sigs["manual_trigger"] = "Enable"
        self.cam.stage_sigs["num_images"] = 1
        self.cam.stage_sigs["num_exposures"] = 1
        # TODO: I don't like this too much, would prefer that we set this for
        # each scan.
        self.cam.stage_sigs["num_triggers"] = int(1e5)

    def setup_external_trigger(self, trigger_type="gate"):
        if trigger_type not in "gate rising_edge".split():
            raise ValueError(
                "trigger_type must be either 'gate' or 'rising_edge', but"
                f"{trigger_type} was entered."
            )

        if trigger_type == "rising_edge":
            # Stage signals
            self.cam.stage_sigs["trigger_mode"] = "External Enable"
            self.cam.stage_sigs["manual_trigger"] = "Disable"
            self.cam.stage_sigs["num_images"] = 1
            self.cam.stage_sigs["num_exposures"] = 1
            # TODO: We may not need this.
            self.cam.stage_sigs["num_triggers"] = self.max_num_images

        elif trigger_type == "gate":

            # Stage signals
            self.cam.stage_sigs["num_triggers"] = 1
            # The num_triggers need to be the first in the Ordered dict! This is
            # because in EPICS, if trigger_mode = External Gate, then cannot
            # change the num_triggers.
            self.cam.stage_sigs.move_to_end("num_triggers", last=False)

            self.cam.stage_sigs["trigger_mode"] = "External Gate"
            self.cam.stage_sigs["manual_trigger"] = "Disable"
            self.cam.stage_sigs["num_images"] = self.max_num_images
            self.cam.stage_sigs["num_exposures"] = 1

    def stage(self):
        if self._flysetup:
            self.setup_external_trigger()

        # Make sure that detector is not armed.
        self.cam.acquire.set(0).wait(timeout=10)
        super().stage()
        self.cam.acquire.set(1).wait(timeout=10)

    def unstage(self):
        super().unstage()
        self.cam.acquire.set(0).wait(timeout=10)

        def check_value(*, old_value, value, **kwargs):
            "Return True when detector is done"
            return (value == "Ready" or value == "Acquisition aborted")

        # When stopping the detector, it may take some time processing the
        # images. This will block until it's done.
        status_wait(
            SubscriptionStatus(
                self.cam.status_message, check_value, timeout=10
            )
        )
        self._flysetup = False
        self.setup_manual_trigger()

    def trigger(self):
        "Trigger one acquisition."
        if self._staged != Staged.yes:
            raise RuntimeError("This detector is not ready to trigger."
                               "Call the stage() method before triggering.")

        @run_in_thread
        def add_delay(status_obj, min_period):
            count_time = self.cam.acquire_time.get()
            total_sleep = count_time if count_time > min_period else min_period
            sleep(total_sleep)
            status_obj.set_finished()

        self._status = self._status_type(self)
        self.acquisition_signal.put(1, wait=False)
        if self.hdf1.enable.get() in (True, 1, "on", "Enable"):
            self.generate_datum(self._image_name, ttime(), {})
        add_delay(self._status, self._min_period)
        return self._status


class Eiger1MDetector(TriggerTime, DetectorBase):

    _default_configuration_attrs = (
        'roi1', 'roi2', 'roi3', 'roi4', 'codec', 'image',
    )
    _default_read_attrs = (
        'cam', 'hdf1', 'stats1', 'stats2', 'stats3', 'stats4', 'stats5'
    )

    cam = ADComponent(EigerDetectorCam, "cam1:")
    codec = ADComponent(CodecPlugin, "Codec1:")
    proc = ADComponent(ProcessPlugin, "Proc1:")
    trans = ADComponent(TransformPlugin, "Trans1:")
    image = ADComponent(ImagePlugin, "image1:")
    hdf1 = ADComponent(PolarHDF5Plugin, "HDF1:")

    roi1 = ADComponent(ROIPlugin, "ROI1:")
    roi2 = ADComponent(ROIPlugin, "ROI2:")
    roi3 = ADComponent(ROIPlugin, "ROI3:")
    roi4 = ADComponent(ROIPlugin, "ROI4:")
    stats1 = ADComponent(StatsPlugin, "Stats1:", kind="normal")
    stats2 = ADComponent(StatsPlugin, "Stats2:", kind="normal")
    stats3 = ADComponent(StatsPlugin, "Stats3:", kind="normal")
    stats4 = ADComponent(StatsPlugin, "Stats4:", kind="normal")
    stats5 = ADComponent(StatsPlugin, "Stats5:", kind="normal")

    def __init__(
        self,
        *args,
        default_folder="",
        hdf1_name_template="%s/%s_%6.6d",
        hdf1_file_extension="h5",
        max_num_images=600000,
        **kwargs
    ):
        self.default_folder = default_folder
        self.hdf1_name_format = (
            hdf1_name_template + "." + hdf1_file_extension
        )
        self.max_num_images = max_num_images
        super().__init__(*args, **kwargs)

    # Make this compatible with other detectors
    @property
    def preset_monitor(self):
        return self.cam.acquire_time

    def align_on(self, time=0.1):
        """Start detector in alignment mode"""
        self.save_images_off()
        self.cam.manual_trigger.set("Disable").wait(timeout=10)
        self.cam.num_triggers.set(int(1e9)).wait(timeout=10)
        self.cam.trigger_mode.set("Continuous").wait(timeout=10)
        self.preset_monitor.set(time).wait(timeout=10)
        self.cam.acquire.set(1).wait(timeout=10)

    def align_off(self):
        """Stop detector"""
        self.cam.acquire.set(0).wait(timeout=10)

    def save_images_on(self):
        self.hdf1.enable.set("Enable").wait(timeout=10)

    def save_images_off(self):
        self.hdf1.enable.set("Disable").wait(timeout=10)

    def auto_save_on(self):
        self.hdf1.autosave.put("on")

    def auto_save_off(self):
        self.hdf1.autosave.put("off")

    def default_settings(self):

        self.cam.num_triggers.put(1)
        self.cam.manual_trigger.put("Disable")
        self.cam.trigger_mode.put("Internal Enable")
        self.cam.acquire.put(0)

        self.hdf1.file_template.put(self.hdf1_name_format)
        self.hdf1.file_path.put(str(self.default_folder))
        self.hdf1.num_capture.put(0)

        self.hdf1.stage_sigs.pop("enable")
        self.hdf1.stage_sigs["num_capture"] = 0
        self.hdf1.stage_sigs["capture"] = 1

        self.setup_manual_trigger()
        self.save_images_off()
        self.plot_stats1()

    def plot_all(self):
        self.plot_select([1, 2, 3, 4, 5])

    def plot_stats1(self):
        self.plot_select([1])

    def plot_stats2(self):
        self.plot_select([2])

    def plot_stats3(self):
        self.plot_select([3])

    def plot_stats4(self):
        self.plot_select([4])

    def plot_stats5(self):
        self.plot_select([5])

    def setup_images(
            self, base_path, name_template, file_number, flyscan=False
    ):

        self.hdf1.file_number.set(file_number).wait(timeout=10)
        self.hdf1.file_name.set(name_template).wait(timeout=10)
        # Make sure eiger will save image
        self.auto_save_on()
        # Changes the stage_sigs to the external trigger mode
        self._flysetup = flyscan

        base_path = str(base_path) + f"/{self.name}/"
        self.hdf1.file_path.set(base_path).wait(timeout=10)

        _, full_path, relative_path = self.hdf1.make_write_read_paths(base_path)

        return Path(full_path), Path(relative_path)

    def plot_select(self, stats):
        """
        Selects which stats will be plotted. All are being read.

        This assumes that 5 stats are setup in Bluesky.

        PARAMETERS
        ----------
        stats : iterable of ints
            List with the stats numbers to be plotted.
        """

        for i in range(1, 5 + 1):
            getattr(self, f"stats{i}").total.kind = (
                "hinted" if i in stats else "normal"
            )

    @property
    def save_image_flag(self):
        _hdf1_auto = True if self.hdf1.autosave.get() == "on" else False
        _hdf1_on = True if self.hdf1.enable.get() == "Enable" else False
        return _hdf1_on or _hdf1_auto

    @property
    def label_option_map(self):
        return {f"Stats{i} Total": i for i in range(1, 5 + 1)}

    @property
    def plot_options(self):
        # Return all named scaler channels
        return list(self.label_option_map.keys())

    def select_plot(self, channels):
        chans = [self.label_option_map[i] for i in channels]
        self.plot_select(chans)
