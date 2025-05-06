"""
Dante CAM
"""

from ophyd import (
    EpicsSignalRO, EpicsSignal, DynamicDeviceComponent
)
from ophyd.areadetector import ADBase, ADComponent, EpicsSignalWithRBV, ad_group
from collections import OrderedDict
from time import sleep
from .ad_mixins import PolarHDF5Plugin


class DanteCAM(ADBase):

    _default_configuration_attrs = (
        'port_name',
        'manufacturer',
        'model',
        'firmware',
        'sdk_version',
        'driver_version',
        'adcore_version',
        'connected',
        'array_callbacks',
        'mca_mode',
        'mca_channels',
        'mca_mapping_points',
        'mca_gatting',
        'mca_list_buffer_size',
        'snl_connected'
    )

    _default_read_attrs = (
        'real_time_preset',
        'real_time_elapsed',
        'real_time_live',
        'instant_deadtime',
        'average_deadtime',
    )

    # Setup
    port_name = ADComponent(EpicsSignalRO, "PortName_RBV")
    manufacturer = ADComponent(EpicsSignalRO, "Manufacturer_RBV")
    model = ADComponent(EpicsSignalRO, "Model_RBV")
    firmware = ADComponent(EpicsSignalRO, "FirmwareVersion_RBV")
    sdk_version = ADComponent(EpicsSignalRO, "SDKVersion_RBV")
    driver_version = ADComponent(EpicsSignalRO, "DriverVersion_RBV")
    adcore_version = ADComponent(EpicsSignalRO, "ADCoreVersion_RBV")
    connected = ADComponent(EpicsSignal, "AsynIO.CNCT")

    array_size = DynamicDeviceComponent(
        ad_group(
            EpicsSignalRO,
            (
                ("array_size_x", "ArraySizeX_RBV"),
                ("array_size_y", "ArraySizeY_RBV"),
                ("array_size_z", "ArraySizeZ_RBV")
            )
        )
    )

    color_mode = ADComponent(EpicsSignalWithRBV, "ColorMode")
    data_type = ADComponent(EpicsSignalWithRBV, "DataType")

    # Acquire
    acquire_start = ADComponent(EpicsSignal, "EraseStart")
    acquire_stop = ADComponent(EpicsSignal, "StopAll")
    acquire_status = ADComponent(EpicsSignalRO, "MCAAcquiring")
    acquire_busy = ADComponent(EpicsSignalRO, "AcquireBusy")

    real_time_preset = ADComponent(EpicsSignal, "PresetReal")
    real_time_elapsed = ADComponent(EpicsSignalRO, "ElapsedReal")
    real_time_live = ADComponent(EpicsSignalRO, "ElapsedLive")

    instant_deadtime = ADComponent(EpicsSignalRO, "IDeadTime")
    average_deadtime = ADComponent(EpicsSignalRO, "DeadTime")

    current_pixel = ADComponent(EpicsSignalRO, "CurrentPixel")
    poll_time = ADComponent(EpicsSignalWithRBV, "PollTime")

    read_rate = ADComponent(EpicsSignal, "ReadAll.SCAN")
    read_rate_button = ADComponent(EpicsSignal, "ReadAllOnce.PROC")

    queued_arrays = ADComponent(EpicsSignalRO, "NumQueuedArrays")

    wait_for_plugins = ADComponent(EpicsSignal, "WaitForPlugins")

    array_counter = ADComponent(EpicsSignalWithRBV, "ArrayCounter")
    image_rate = ADComponent(EpicsSignalRO, "ArrayRate_RBV")

    array_callbacks = ADComponent(EpicsSignal, "ArrayCallbacks")

    # MCA setup
    mca_mode = ADComponent(EpicsSignalWithRBV, "CollectMode")
    mca_channels = ADComponent(EpicsSignalWithRBV, "NumMCAChannels")
    mca_mapping_points = ADComponent(EpicsSignalWithRBV, "MappingPoints")
    num_images = ADComponent(EpicsSignalWithRBV, "MappingPoints")
    mca_gatting = ADComponent(EpicsSignalWithRBV, "GatingMode")
    mca_list_buffer_size = ADComponent(EpicsSignalWithRBV, "ListBufferSize")

    # Multi Channel
    snl_connected = ADComponent(EpicsSignalRO, "SNL_Connected")


class DanteSCA(ADBase):

    _default_read_attrs = ("icr", "ocr", "f1_deadtime")

    _default_configuration_attrs = (
        'enable',
        'fast_peaking_time',
        'fast_threshold',
        'fast_flat_top_time',
        'peaking_time',
        'max_peaking_time',
        'energy_threshold',
        'baseline_threshold',
        'max_rise_time',
        'reset_recovery_time',
        'zero_peak_frequency',
        'baseline_samples',
        'gain',
        'input_mode',
        'input_polarity',
        'analog_offset',
        'base_offset',
        'reset_threshold',
        'time_constant',
        'max_energy'
    )

    # Statistics
    real_time = ADComponent(EpicsSignalRO, "ElapsedRealTime")
    live_time = ADComponent(EpicsSignalRO, "ElapsedLiveTime")
    icr = ADComponent(EpicsSignalRO, "InputCountRate")
    ocr = ADComponent(EpicsSignalRO, "OutputCountRate")
    triggers = ADComponent(EpicsSignalRO, "Triggers")
    events = ADComponent(EpicsSignalRO, "Events")
    fast_deadtime = ADComponent(EpicsSignalRO, "FastDeadTime")
    f1_deadtime = ADComponent(EpicsSignalRO, "F1DeadTime")
    zero_counts = ADComponent(EpicsSignalRO, "ZeroCounts")
    baseline_counts = ADComponent(EpicsSignalRO, "BaselineCount")
    pileup = ADComponent(EpicsSignalRO, "PileUp")
    f1_pileup = ADComponent(EpicsSignalRO, "F1PileUp")
    not_f1_pileup = ADComponent(EpicsSignalRO, "NotF1PileUp")
    reset_counts = ADComponent(EpicsSignalRO, "ResetCounts")

    # Parameters
    enable = ADComponent(EpicsSignal, "EnableBoard")
    fast_peaking_time = ADComponent(EpicsSignalWithRBV, "FastPeakingTime")
    fast_threshold = ADComponent(EpicsSignalWithRBV, "FastThreshold")
    fast_flat_top_time = ADComponent(EpicsSignalWithRBV, "FastFlatTopTime")
    peaking_time = ADComponent(EpicsSignalWithRBV, "PeakingTime")
    max_peaking_time = ADComponent(EpicsSignalWithRBV, "MaxPeakingTime")
    energy_threshold = ADComponent(EpicsSignalWithRBV, "EnergyThreshold")
    fast_peaking_time = ADComponent(EpicsSignalWithRBV, "FastPeakingTime")
    baseline_threshold = ADComponent(EpicsSignalWithRBV, "BaselineThreshold")
    max_rise_time = ADComponent(EpicsSignalWithRBV, "MaxRiseTime")
    reset_recovery_time = ADComponent(EpicsSignalWithRBV, "ResetRecoveryTime")
    zero_peak_frequency = ADComponent(EpicsSignalWithRBV, "ZeroPeakFreq")
    baseline_samples = ADComponent(EpicsSignalWithRBV, "BaselineSamples")
    gain = ADComponent(EpicsSignalWithRBV, "Gain")
    input_mode = ADComponent(EpicsSignalWithRBV, "InputMode")
    input_polarity = ADComponent(EpicsSignalWithRBV, "InputPolarity")
    analog_offset = ADComponent(EpicsSignalWithRBV, "AnalogOffset")
    base_offset = ADComponent(EpicsSignalWithRBV, "BaseOffset")
    reset_threshold = ADComponent(EpicsSignalWithRBV, "ResetThreshold")
    time_constant = ADComponent(EpicsSignalWithRBV, "TimeConstant")
    max_energy = ADComponent(EpicsSignalWithRBV, "MaxEnergy")


class DanteHDF1Plugin(PolarHDF5Plugin):
    # The array counter readback pv is different...
    array_counter = ADComponent(EpicsSignal, "ArrayCounter", kind="config")
    array_counter_readback = ADComponent(
        EpicsSignalRO, "ArrayCounter_RBV", kind="config"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._num_images_device = "cam.mca_mapping_points"

    def warmup(self):
        sigs = OrderedDict(
            [
                (self.enable, 1),
                (self.parent.cam.array_callbacks, 1),  # set by number
                (self.parent.cam.mca_mode, "MCA Mapping"),
                (self.parent.cam.mca_mapping_points, 1),
                (self.parent.preset_monitor, 1),
                (self.parent.cam.acquire_start, 1),  # set by number
            ]
        )

        original_vals = {sig: sig.get() for sig in sigs}

        for sig, val in sigs.items():
            sleep(0.1)  # abundance of caution
            sig.set(val).wait()

        sleep(2)  # wait for acquisition

        for sig, val in reversed(list(original_vals.items())):
            sleep(0.1)
            sig.set(val).wait()
