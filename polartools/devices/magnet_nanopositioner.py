'''
Magnet Nanopositioner motors
'''

from ophyd import Component, MotorBundle, EpicsMotor


class NanoPositioner(MotorBundle):
    x = Component(EpicsMotor, 'm1', labels=('motor',))
    y = Component(EpicsMotor, 'm2', labels=('motor',))
    z = Component(EpicsMotor, 'm3', labels=('motor',))
