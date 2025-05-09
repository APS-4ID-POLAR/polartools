 
"""
APS status
"""

from ophyd import (
    Component, EpicsSignalRO, Device
)


class StatusAPS(Device):
    current = Component(EpicsSignalRO, "S-DCCT:CurrentM")
    machine_status = Component(EpicsSignalRO, "S:DesiredMode", string=True)
    operating_mode = Component(EpicsSignalRO, "S:ActualMode", string=True)
    shutter_status = Component(
        EpicsSignalRO, "RF-ACIS:FePermit:Sect1To35IdM", string=True
    )
