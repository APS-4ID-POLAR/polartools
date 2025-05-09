"""
Ring XBPM support
"""

__all__ = ["aps_xbpm"]

from ophyd import Component, Device, EpicsSignalRO


class CMsDevice(Device):
    current1 = Component(EpicsSignalRO, "Current1:MeanValue_NormM")
    current2 = Component(EpicsSignalRO, "Current2:MeanValue_NormM")
    current3 = Component(EpicsSignalRO, "Current3:MeanValue_NormM")
    current4 = Component(EpicsSignalRO, "Current4:MeanValue_NormM")


class MyXBPM(Device):
    vertical_position = Component(EpicsSignalRO, ":ID:SrcPt:VPositionM")
    vertical_angle = Component(EpicsSignalRO, ":ID:SrcPt:VAngleM")
    horizontal_position = Component(EpicsSignalRO, ":ID:SrcPt:HPositionM")
    horizontal_angle = Component(EpicsSignalRO, ":ID:SrcPt:HAngleM")

    x_axis = Component(EpicsSignalRO, "IDFE-XBPM:P1ds:XAxisNorm_")
    y_top = Component(EpicsSignalRO, "IDFE-XBPM:P1ds:YDataTopNorm_")
    y_bot = Component(EpicsSignalRO, "IDFE-XBPM:P1ds:YDataBottomNorm_")

    cm1 = Component(CMsDevice, "IDFE-XBPM:CM1ds:")
    cm2 = Component(CMsDevice, "IDFE-XBPM:CM2ds:")
    cm3 = Component(CMsDevice, "IDFE-XBPM:CM3ds:")


aps_xbpm = MyXBPM("S04", name="aps_xbpm", labels=("source",))
