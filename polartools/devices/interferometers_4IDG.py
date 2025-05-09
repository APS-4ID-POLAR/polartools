""" Interferometer setup """

from ophyd import Device, Component, EpicsSignalRO


class InterferometerDevice(Device):
    mhor_up = Component(EpicsSignalRO, "pixelTrig-1_POS1")

    mhor_down = Component(EpicsSignalRO, "pixelTrig-1_POS2")

    shor = Component(EpicsSignalRO, "pixelTrig-2_POS1")
    mvert_up = Component(EpicsSignalRO, "pixelTrig-2_POS2")

    mvert_down = Component(EpicsSignalRO, "pixelTrig-3_POS1")
    svert = Component(EpicsSignalRO, "pixelTrig-3_POS2")

    def plot_first_pos1(self):
        self.mhor_up.kind = "hinted"
        self.mhor_down.kind = "normal"
        self.shor.kind = "normal"
        self.mvert_up.kind = "normal"
        self.mvert_down.kind = "normal"
        self.svert.kind = "normal"

    def plot_all(self):
        self.mhor_up.kind = "hinted"
        self.mhor_down.kind = "hinted"
        self.shor.kind = "hinted"
        self.mvert_up.kind = "hinted"
        self.mvert_down.kind = "hinted"
        self.svert.kind = "hinted"
