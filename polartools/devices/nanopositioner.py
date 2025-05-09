'''
Nanopositioner motors
'''

from ophyd import Component, MotorBundle, EpicsMotor


class MyEpicsMotor(EpicsMotor):
    def unstage(self):
        try:
            self.stage_sigs.pop("velocity")
        except KeyError:
            pass
        return super().unstage()


class NanoPositioner(MotorBundle):
    nanoy = Component(MyEpicsMotor, 'm1')
    nanox = Component(MyEpicsMotor, 'm2')
    nanoz = Component(MyEpicsMotor, 'm3')
