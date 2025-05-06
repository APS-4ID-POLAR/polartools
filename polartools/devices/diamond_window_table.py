"""
4-ID-B diamond window motors
"""

__all__ = [
    'diamond_window'
]

from ophyd import Component, Device, EpicsMotor
from ..utils._logging_setup import logger
logger.info(__file__)


class WindowStages(Device):
    x = Component(EpicsMotor, "m1", labels=("motor",))
    y = Component(EpicsMotor, "m2", labels=("motor",))


diamond_window = WindowStages(
    "4idbSoft:", name="diamond_window", labels=("4ida",)
)
