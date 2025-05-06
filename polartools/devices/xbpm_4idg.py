"""
4idg XBPM
"""

__all__ = ['gxbpm']

from ophyd import Device, Component, EpicsMotor
from ..utils._logging_setup import logger

logger.info(__file__)


class XBPM(Device):
    x = Component(EpicsMotor, "m48", labels=("motor",))
    y = Component(EpicsMotor, "m47", labels=("motor",))
    # detector


gxbpm = XBPM("4idgSoft:", name="gxbpm", labels=("4idg",))
