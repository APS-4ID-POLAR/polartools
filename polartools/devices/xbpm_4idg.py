"""
4idg XBPM
"""

from ophyd import Device, Component, EpicsMotor


class XBPM(Device):
    x = Component(EpicsMotor, "m48", labels=("motor",))
    y = Component(EpicsMotor, "m47", labels=("motor",))
    # detector
