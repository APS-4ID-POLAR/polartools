 
"""
APS status
"""
__all__ = ["status_aps"]

from ophyd import (
    Component, EpicsSignalRO, Device
)
from ..utils._logging_setup import logger
logger.info(__file__)


class StatusAPS(Device):
    current = Component(EpicsSignalRO, "S-DCCT:CurrentM")
    machine_status = Component(EpicsSignalRO, "S:DesiredMode", string=True)
    operating_mode = Component(EpicsSignalRO, "S:ActualMode", string=True)
    shutter_status = Component(
        EpicsSignalRO, "RF-ACIS:FePermit:Sect1To35IdM", string=True
    )


status_aps = StatusAPS("", name="status_aps", labels=("source",))
