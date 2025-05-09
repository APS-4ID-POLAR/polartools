
"""
Electromagnet
"""

from ophyd import Component, Device, EpicsMotor
from .magnet_kepco_4idb import KepcoController


class Magnet2T(Device):
    # tablex = Component(EpicsMotor, "4idbSoft:m15")
    # tabley = Component(EpicsMotor, "4idbSoft:m16")
    samplex = Component(EpicsMotor, "4idb:m25", labels=("motor",))
    sampley = Component(EpicsMotor, "4idb:m17", labels=("motor",))
    kepco = Component(KepcoController, '4idbSoft:BOP:PS1:', labels=("magnet",))

    def default_settings(self):
        self.kepco.mode_change(value=self.kepco.mode.get())
