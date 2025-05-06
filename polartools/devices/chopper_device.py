"""
Chopper
"""

__all__ = ["chopper"]


from ophyd import (
    Component,
    FormattedComponent,
    Device,
    EpicsSignal,
    EpicsSignalRO,
    EpicsMotor
)
from ..utils._logging_setup import logger
logger.info(__file__)


class ChopperDevice(Device):

    # In/out motor
    translation = FormattedComponent(
        EpicsMotor, "4idbSoft:m14", labels=("motor",)
    )

    # Newport controller
    status = Component(EpicsSignalRO, "Connected", kind="config")

    frequency_readback = Component(EpicsSignalRO, "FreqSync")
    frequency_setpoint = Component(EpicsSignal, "FrequencySet")

    phase = Component(EpicsSignal, "PhaseDelaySet", kind="config")

    wheel = Component(EpicsSignal, "WheelSet", string=True, kind="config")
    sync = Component(EpicsSignal, "SyncSourceSet", string=True, kind="config")
    mode = Component(EpicsSignal, "ModeSet", string=True, kind="config")


chopper = ChopperDevice("4idChopper:", name="chopper", labels=("4idb",))
