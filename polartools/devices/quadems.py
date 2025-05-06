"""
QuadEMs for POLAR
"""

from ophyd import Component, QuadEM, EpicsSignalRO, Device
from ophyd.quadem import QuadEMPort
from collections import OrderedDict
from .ad_mixins import ImagePlugin, StatsPlugin
from ..utils._logging_setup import logger
logger.info(__file__)


class StatsPluginQuadEM(StatsPlugin):
    # Remove subscriptions from StatsPlugin
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stop_auto_kind()
        self.kind = "config"


class QuadEMPOLAR(QuadEM):
    image = Component(ImagePlugin, "image1:")
    current1 = Component(StatsPluginQuadEM, "Current1:")
    current2 = Component(StatsPluginQuadEM, "Current2:")
    current3 = Component(StatsPluginQuadEM, "Current3:")
    current4 = Component(StatsPluginQuadEM, "Current4:")

    sum_all = Component(StatsPluginQuadEM, "SumAll:")

    # The way the QuadEM support computes things is a bit complicated, so
    # will expose main screen here eventhough it will be a duplicate of the
    # stats, and will leave the stats as config.

    sumall_mean = Component(EpicsSignalRO, "SumAll:MeanValue_RBV")
    sumall_fast = Component(EpicsSignalRO, "SumAllAve")
    sumall_sigma = Component(EpicsSignalRO, "SumAll:Sigma_RBV")

    sumx_mean = Component(EpicsSignalRO, "SumX:MeanValue_RBV")
    sumx_fast = Component(EpicsSignalRO, "SumXAve")
    sumx_sigma = Component(EpicsSignalRO, "SumX:Sigma_RBV")

    sumy_mean = Component(EpicsSignalRO, "SumY:MeanValue_RBV")
    sumy_fast = Component(EpicsSignalRO, "SumYAve")
    sumy_sigma = Component(EpicsSignalRO, "SumY:Sigma_RBV")

    diffx_mean = Component(EpicsSignalRO, "DiffX:MeanValue_RBV")
    diffx_fast = Component(EpicsSignalRO, "DiffXAve")
    diffx_sigma = Component(EpicsSignalRO, "DiffX:Sigma_RBV")

    diffy_mean = Component(EpicsSignalRO, "DiffY:MeanValue_RBV")
    diffy_fast = Component(EpicsSignalRO, "DiffYAve")
    diffy_sigma = Component(EpicsSignalRO, "DiffY:Sigma_RBV")

    posx_mean = Component(EpicsSignalRO, "PosX:MeanValue_RBV")
    posx_fast = Component(EpicsSignalRO, "PositionXAve")
    posx_sigma = Component(EpicsSignalRO, "PosX:Sigma_RBV")

    posy_mean = Component(EpicsSignalRO, "PosY:MeanValue_RBV")
    posy_fast = Component(EpicsSignalRO, "PositionYAve")
    posy_sigma = Component(EpicsSignalRO, "PosY:Sigma_RBV")

    @property
    def preset_monitor(self):
        return self.averaging_time


class TetrAMM(QuadEMPOLAR):
    conf = Component(QuadEMPort, port_name="TetrAMM")

    # TODO: If we ever want to trigger the TetrAMM, we need
    # to check if changes are needed to the trigger procedure.


class QuadEMRO_mixins:
    # Disables preset_monitor and trigger

    def trigger(self):
        self._status = self._status_type(self)
        self._status.set_finished()
        return self._status

    @property
    def preset_monitor(self, value):
        pass

    def stage(self):
        Device.stage(self)

    def unstage(self):
        Device.unstage(self)


class SydorEMRO(QuadEMRO_mixins, QuadEMPOLAR):

    conf = Component(QuadEMPort, port_name="T4U_BPM")

    # These are TetrAMM specific!
    num_channels = None
    read_format = None
    trigger_mode = None
    bias_interlock = None
    bias_state = None
    bias_voltage = None
    hvi_readback = None
    hvs_readback = None
    hvv_readback = None

    image = None

    def default_settings(self):
        # Remove all these from read_attrs
        for item in (
            'conf',
            'current_names',
            'current_names.ch1',
            'current_names.ch2',
            'current_names.ch3',
            'current_names.ch4',
            'current_offsets',
            'current_offsets.ch1',
            'current_offsets.ch2',
            'current_offsets.ch3',
            'current_offsets.ch4',
            'current_offset_calcs',
            'current_offset_calcs.ch1',
            'current_offset_calcs.ch2',
            'current_offset_calcs.ch3',
            'current_offset_calcs.ch4',
            'current_scales',
            'current_scales.ch1',
            'current_scales.ch2',
            'current_scales.ch3',
            'current_scales.ch4',
            'position_offset_x',
            'position_offset_y',
            'position_offset_calc_x',
            'position_offset_calc_y',
            'position_scale_x',
            'position_scale_y'
        ):
            getattr(self, item).kind = "config"

        self.stage_sigs = OrderedDict()


class TetrAMMRO(QuadEMRO_mixins, TetrAMM):

    def default_settings(self):
        # Remove all these from read_attrs
        for item in (
            'conf',
            'current_names',
            'current_names.ch1',
            'current_names.ch2',
            'current_names.ch3',
            'current_names.ch4',
            'current_offsets',
            'current_offsets.ch1',
            'current_offsets.ch2',
            'current_offsets.ch3',
            'current_offsets.ch4',
            'current_offset_calcs',
            'current_offset_calcs.ch1',
            'current_offset_calcs.ch2',
            'current_offset_calcs.ch3',
            'current_offset_calcs.ch4',
            'current_scales',
            'current_scales.ch1',
            'current_scales.ch2',
            'current_scales.ch3',
            'current_scales.ch4',
            'position_offset_x',
            'position_offset_y',
            'position_offset_calc_x',
            'position_offset_calc_y',
            'position_scale_x',
            'position_scale_y'
        ):
            getattr(self, item).kind = "config"

        self.stage_sigs = OrderedDict()
