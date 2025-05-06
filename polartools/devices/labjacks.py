"""
Labjacks
"""

__all__ = [
    "labjack_4idb",
    "labjack_4ida"
]

from apstools.devices import LabJackT7
from apstools.devices.labjack import (
    KIND_CONFIG_OR_NORMAL, DigitalIO, Output
)
from ophyd import DynamicDeviceComponent, EpicsSignalRO, Component, EpicsSignal
from ..utils._logging_setup import logger
logger.info(__file__)


class AnalogOutput(Output):
    description = Component(EpicsSignal, ".DESC", kind="config")
    value = Component(EpicsSignal, "", kind="hinted")

    low_limit = Component(EpicsSignal, ".DRVL", kind="config")
    high_limit = Component(EpicsSignal, ".DRVH", kind="config")

    readback_value = None
    desired_value = None


def make_analog_outputs(num_aos: int):
    """Create a dictionary with analog output device definitions.

    For use with an ophyd DynamicDeviceComponent.

    Parameters
    ==========
    num_aos
      How many analog outputs to create.

    """
    defn = {}
    for n in range(num_aos):
        defn[f"ao{n}"] = (
            AnalogOutput, f"Ao{n}", dict(kind=KIND_CONFIG_OR_NORMAL)
        )
    return defn


def make_digital_ios(channels_list: list):
    """Create a dictionary with digital I/O device definitions.

    For use with an ophyd DynamicDeviceComponent.

    Parameters
    ==========
    num_dios
      How many digital I/Os to create.
    """
    defn = {}
    for n in channels_list:
        defn[f"dio{n}"] = (
            DigitalIO, "", dict(ch_num=n, kind=KIND_CONFIG_OR_NORMAL)
        )

    # Add the digital word outputs
    defn["dio"] = (EpicsSignalRO, "DIOIn", dict(kind=KIND_CONFIG_OR_NORMAL))
    defn["fio"] = (EpicsSignalRO, "FIOIn", dict(kind=KIND_CONFIG_OR_NORMAL))
    defn["eio"] = (EpicsSignalRO, "EIOIn", dict(kind=KIND_CONFIG_OR_NORMAL))
    defn["cio"] = (EpicsSignalRO, "CIOIn", dict(kind=KIND_CONFIG_OR_NORMAL))
    defn["mio"] = (EpicsSignalRO, "MIOIn", dict(kind=KIND_CONFIG_OR_NORMAL))
    return defn


class CustomLabJackT7(LabJackT7):
    # In the "default" BCDA setup, four IO channels (all CIO, #16-19) are
    # converted into analog outputs (thus now 6 DACs)

    analog_outputs = DynamicDeviceComponent(
        make_analog_outputs(6), kind=KIND_CONFIG_OR_NORMAL
    )

    digital_ios = DynamicDeviceComponent(
        make_digital_ios(list(range(0, 16)) + list(range(20, 23))),
        kind=KIND_CONFIG_OR_NORMAL
    )


labjack_4idb = CustomLabJackT7(
    "4idLabJackT7_1:", name="labjack_4idb", labels=("4idb",)
)

labjack_4ida = CustomLabJackT7(
    "4idaSoft:LJ:", name="labjack_4ida", labels=("4ida",)
)
labjack_4ida.analog_outputs.kind = "normal"
labjack_4ida.waveform_digitizer.kind = "omitted"
labjack_4ida.digital_ios.kind = "omitted"
labjack_4ida.analog_inputs.kind = "omitted"
for i in range(4):
    getattr(labjack_4ida.analog_outputs, f"ao{i}").kind = "normal"
