"""
4-ID-B diamond window motors
"""

from ophyd import Component, Device, EpicsMotor


class WindowStages(Device):
    x = Component(EpicsMotor, "m1", labels=("motor",))
    y = Component(EpicsMotor, "m2", labels=("motor",))
