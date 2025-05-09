
"""
Phase Shifter
"""

from ophyd import Device, Component, EpicsSignalRO, EpicsSignal
from ophyd.status import DeviceStatus
from apstools.devices.aps_undulator import UndulatorPositioner


class PSPositioner(UndulatorPositioner):

    def move(self, position, wait=True, timeout=None, moved_cb=None):
        if (
            abs(position - self.readback.get(use_monitor=False))
            < self.parent.gap_deadband.get(use_monitor=False) / 1000
            # deadband in microns, gap in mm
        ):
            status = DeviceStatus(self)
            status.set_finished()
        else:
            status = super().move(
                position, wait=wait, timeout=timeout, moved_cb=moved_cb
            )

        return status


class PhaseShifterDevice(Device):
    gap = Component(PSPositioner, "Gap")

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
