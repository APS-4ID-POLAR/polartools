
"""
Polar status
"""


from ophyd import (
    Component, FormattedComponent, EpicsSignalRO, Device
)


class GHStatus(Device):

    user_enable = FormattedComponent(
        EpicsSignalRO, "{prefix}{_hutch}_USER_KEY", string=True
    )

    aps_enable = FormattedComponent(
        EpicsSignalRO, "{prefix}{_hutch}_APS_KEY", string=True
    )

    searched = FormattedComponent(
        EpicsSignalRO, "{prefix}{_hutch}_SEARCHED", string=True
    )

    beam_present = FormattedComponent(
        EpicsSignalRO, "{prefix}{_hutch}_BEAM_PRESENT", string=True
    )

    manual_shutter = FormattedComponent(
        EpicsSignalRO, "{prefix}S{_hutch}S_BLOCKING_BEAM.VAL", string=True
    )

    def __init__(self, prefix, hutch=None, **kwargs):
        self._hutch = hutch
        super().__init__(prefix, **kwargs)


class ABStatus(GHStatus):
    beam_ready = FormattedComponent(
        EpicsSignalRO, "{prefix}{self._beam_ready}", string=True
    )

    manual_shutter = None

    def __init__(self, prefix, hutch=None, **kwargs):
        if hutch == "A":
            self._beam_ready = "FES_PERMIT.VAL"
        elif hutch == "B":
            self._beam_ready = "SBS_PERMIT.VAL"
        super().__init__(prefix, hutch=hutch, **kwargs)


class Status4ID(Device):

    online = Component(
        EpicsSignalRO, "ACIS_GLOBAL_ONLINE.VAL", string=True
    )
    acis = Component(
        EpicsSignalRO, "ACIS_FES_PERMIT.VAL", string=True
    )

    a_hutch = Component(ABStatus, "", hutch="A", labels=("4ida",))
    b_hutch = Component(ABStatus, "", hutch="B", labels=("4idb",))
    g_hutch = Component(GHStatus, "", hutch="G", labels=("4idg",))
    h_hutch = Component(GHStatus, "", hutch="H", labels=("4idh",))
