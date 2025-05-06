'''
Nanopositioner motors
'''

__all__ = ['diff_nano']

from ophyd import Component, MotorBundle, EpicsMotor
from ..utils import logger
logger.info(__file__)


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


diff_nano = NanoPositioner(
    '4idIF:', name='diff_nano', labels=('motor', 'nanopositioner')
)
