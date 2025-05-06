"""
APS filter support
"""

__all__ = ["gfilter"]


from ophyd import (
    Component,
    DynamicDeviceComponent,
    Device,
    EpicsSignal,
    EpicsSignalRO,
    EpicsSignalWithRBV
)
from ..utils._logging_setup import logger
logger.info(__file__)

NUM_FILTERS = 12


class FilterSlot(Device):
    setpoint = Component(EpicsSignal, "", string=True)
    readback = Component(EpicsSignalRO, "_RBV", string=True)

    lock = Component(EpicsSignal, "_Lock", string=True, kind="config")
    enable = Component(EpicsSignal, "_Enable", string=True, kind="config")
    material = Component(EpicsSignal, "_material", string=True, kind="config")
    thickness = Component(EpicsSignal, "_thickness", string=True, kind="config")

    control_pv = Component(EpicsSignalRO, "_calc.OUT", kind="config")
    control_in_value = Component(EpicsSignalRO, "_calc.B", kind="config")

    readback_pv = Component(EpicsSignalRO, "_RBV_calc.INPA", kind="config")
    readback_in_value = Component(EpicsSignalRO, "_RBV_calc.B", kind="config")


def make_filter_slots(num: int):
    defn = {}
    for n in range(1, num+1):
        defn[f"f{n}"] = (
            FilterSlot, f"filter{n :02d}", dict(kind="config")
        )
    return defn


class APSFilter(Device):

    # Status and information

    energy_select = Component(
        EpicsSignal, "EnergySelect", string=True, kind="config"
    )
    mono_energy = Component(EpicsSignalRO, "EnergyBeamline", kind="config")
    local_energy = Component(EpicsSignal, "EnergyLocal", kind="config")

    energy = Component(EpicsSignalRO, "energy_RBV")

    status = Component(EpicsSignalRO, "filterBusy", string=True)

    attenuation_setpoint = Component(EpicsSignal, "attenuation")
    attenuation_readback = Component(
        EpicsSignalRO, "attenuation_actual", kind="config"
    )

    sorted_index = Component(EpicsSignalWithRBV, "sortedIndex")

    attenuation_2e_harmonic = Component(EpicsSignalRO, "attenuation_2E_actual")
    attenuation_3e_harmonic = Component(EpicsSignalRO, "attenuation_3E_actual")

    # Configuration
    inter_filter_delay = Component(
        EpicsSignal, "interFilterDelay", kind="config"
    )

    filters = DynamicDeviceComponent(make_filter_slots(NUM_FILTERS))


gfilter = APSFilter(
    "4idPyFilter:FL1:", name="gfilter", labels=("4idg", "filter")
)
