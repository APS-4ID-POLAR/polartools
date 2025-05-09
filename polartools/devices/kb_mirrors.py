"""
KB mirror
"""

from ophyd import Component, Device, EpicsMotor


class KBMirror(Device):
    # Overal motors
    x = Component(EpicsMotor, "m16", labels=("motor",))
    rot = Component(EpicsMotor, "m15", labels=("motor",))

    # KB mirror setup
    # TODO
