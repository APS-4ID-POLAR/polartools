'''
Ventus Laser
'''

__all__ = ["laser"]

from ophyd import EpicsSignal, Component, Device
from ..utils import logger
logger.info(__file__)


class VentusLaser(Device):

    # TODO: Should power and current be PVPositioners?
    power_setpoint = Component(EpicsSignal, "PowerOut", kind="normal")
    power_readback = Component(EpicsSignal, "PowerIn", kind="normal")

    current_setpoint = Component(EpicsSignal, "CurrentOut", kind="normal")
    current_readback = Component(EpicsSignal, "CurrentIn", kind="normal")

    laser_on = Component(EpicsSignal, "On", kind="omitted")
    laser_off = Component(EpicsSignal, "Off", kind="omitted")

    laser_temperature = Component(EpicsSignal, "LaserTemp", kind="config")
    psu_temperature = Component(EpicsSignal, "PSUTemp", kind="config")

    status = Component(EpicsSignal, "Status", kind="config")

    scan = Component(EpicsSignal, "Read.SCAN", kind="config")

    control_mode = Component(EpicsSignal, "ControlMode", kind="config")

    psu_time = Component(EpicsSignal, "PSUTime", kind="config")
    laser_enabled_time = Component(EpicsSignal, "LETime", kind="config")
    laser_operation_time = Component(EpicsSignal, "LOTime", kind="config")


laser = VentusLaser("4tst:LQE1:", name="laser", labels=("raman",))
