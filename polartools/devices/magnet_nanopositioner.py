'''
Magnet Nanopositioner motors
'''

__all__ = ['magnet_nano']

from ophyd import Component, MotorBundle, EpicsMotor
from ..utils import logger
logger.info(__file__)


class NanoPositioner(MotorBundle):
    x = Component(EpicsMotor, 'm1', labels=('motor',))
    y = Component(EpicsMotor, 'm2', labels=('motor',))
    z = Component(EpicsMotor, 'm3', labels=('motor',))


magnet_nano = NanoPositioner(
    'cpscIOC:', name='magnet_nano', labels=('nanopositioner',)
)
