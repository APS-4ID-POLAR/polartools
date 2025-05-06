"""
KB mirror
"""

__all__ = ["kbb"]

from ophyd import Component, Device, EpicsMotor
from ..utils._logging_setup import logger
logger.info(__file__)


class KBMirror(Device):
    # Overal motors
    x = Component(EpicsMotor, "m16", labels=("motor",))
    rot = Component(EpicsMotor, "m15", labels=("motor",))

    # KB mirror setup
    # TODO


kbb = KBMirror("4idbSoft:", name="kbb", labels=("4idb", "optics"))
