"""
Table in middle of 4idb
"""

from ophyd import Device, Component, EpicsMotor


class Table4idb(Device):
    x_us = Component(EpicsMotor, "m5", labels=("motor",))
    x_ds = Component(EpicsMotor, "m8", labels=("motor",))

    y_us = Component(EpicsMotor, "m4", labels=("motor",))
    y_ds_in = Component(EpicsMotor, "m7", labels=("motor",))
    y_ds_out = Component(EpicsMotor, "m6", labels=("motor",))

    # TODO: add the combined motion pseudomotors.
