'''
Sample motors
'''

__all__ = "sx sy sz".split()

from ophyd import EpicsMotor
from ..utils import logger
logger.info(__file__)


sx = EpicsMotor('4tst:m1', labels=('motor'), name="sx")
sy = EpicsMotor('4tst:m2', labels=('motor'), name="sy")
sz = EpicsMotor('4tst:m3', labels=('motor'), name="sz")
