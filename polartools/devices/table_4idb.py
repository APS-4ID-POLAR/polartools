"""
Table in middle of 4idb
"""

__all__ = ['midtable_4idb']

from ophyd import Device, Component, EpicsMotor
from ..utils._logging_setup import logger

logger.info(__file__)


class Table4idb(Device):
    x_us = Component(EpicsMotor, "m5", labels=("motor",))
    x_ds = Component(EpicsMotor, "m8", labels=("motor",))

    y_us = Component(EpicsMotor, "m4", labels=("motor",))
    y_ds_in = Component(EpicsMotor, "m7", labels=("motor",))
    y_ds_out = Component(EpicsMotor, "m6", labels=("motor",))


midtable_4idb = Table4idb("4idbSoft:", name="midtable_4idb", labels=("4idb",))
