"""
Undulator support
"""

from apstools.devices import STI_Undulator, TrackingSignal
from apstools.devices.aps_undulator import UndulatorPositioner
from ophyd import Component, Device, Signal, EpicsSignal, EpicsSignalRO
from ophyd.status import Status, StatusBase
from typing import Any, Callable
from numpy import abs
from ..utils._logging_setup import logger
logger.info(__file__)


class PolarUndulatorPositioner(UndulatorPositioner):

    def set(
        self,
        new_position: Any,
        *,
        timeout: float = None,
        moved_cb: Callable = None,
        wait: bool = False,
    ) -> StatusBase:
        # If position is within the deadband --> do nothing.
        if (
            abs(new_position - self.readback.get()) <
            self.parent.energy_deadband.get()*self.parent.harmonic_value.get()
        ):
            _status = Status()
            _status.set_finished()
        else:
            _status = super().set(
                new_position, timeout=timeout, moved_cb=moved_cb, wait=wait
            )
        return _status


class PolarUndulator(STI_Undulator):
    tracking = Component(TrackingSignal, value=False, kind='config')
    offset = Component(Signal, value=0, kind='config')
    energy_deadband = Component(Signal, value=0.003, kind='config')
    energy = Component(PolarUndulatorPositioner, "Energy")
    version_hpmu = None


class PhaseShifterDevice(Device):
    gap = Component(UndulatorPositioner, "Gap")

    start_button = Component(EpicsSignal, "StartC.VAL")
    stop_button = Component(EpicsSignal, "StopC.VAL")
    done = Component(EpicsSignalRO, "BusyM.VAL", kind="omitted")

    gap_deadband = Component(EpicsSignal, "DeadbandGapC")
    device_limit = Component(EpicsSignal, "DeviceLimitM.VAL")
    device = Component(EpicsSignalRO, "DeviceM", kind="config")
    location = Component(EpicsSignalRO, "LocationM", kind="config")
    message1 = Component(
        EpicsSignalRO, "Message1M.VAL", kind="config", string=True
    )
    message2 = Component(
        EpicsSignalRO, "Message2M.VAL", kind="config", string=True
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.gap.done_value = 0


class PolarUndulatorPair(Device):
    us = Component(PolarUndulator, "USID:", labels=("track_energy",))
    ds = Component(PolarUndulator, "DSID:", labels=("track_energy",))
    phase_shifter = Component(PhaseShifterDevice, "ILPS:")


undulators = PolarUndulatorPair(
    "S04ID:", name="undulators", labels=("energy", "source")
)
