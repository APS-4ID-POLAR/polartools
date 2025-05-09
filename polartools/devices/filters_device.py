"""
APS filter support
"""


from ophyd import (
    Component, DynamicDeviceComponent, Device, EpicsSignal, EpicsSignalRO
)

NUM_FILTERS = 12


class FilterSlot(Device):
    status = Component(EpicsSignal, "Set", string=True)
    lock = Component(EpicsSignal, "Lock", string=True)
    material = Component(EpicsSignal, "Material", string=True)
    thickness = Component(EpicsSignal, "Thickness")
    enable = Component(EpicsSignal, "Enable", string=True)
    transmission = Component(EpicsSignalRO, "Transmission")


def make_filter_slots(num: int):
    defn = {}
    for n in range(1, num + 1):
        defn[f"f{n}"] = (
            FilterSlot, f"Fi{n}:", dict(kind="config")
        )
    return defn


class APSFilter(Device):

    # Status and information

    energy_select = Component(
        EpicsSignal, "EnergySelect", string=True, kind="config"
    )
    mono_energy = Component(EpicsSignalRO, "EnergyBeamline", kind="config")
    local_energy = Component(EpicsSignal, "EnergyLocal", kind="config")

    status = Component(EpicsSignalRO, "Status", string=True, kind="config")

    transmission_readback = Component(EpicsSignalRO, "Transmission")
    transmission_setpoint = Component(
        EpicsSignal, "TransmissionSetpoint", kind="config"
    )
    transmission_factor = Component(
        EpicsSignal, "TransmissionFactor", kind="config"
    )

    mask_readback = Component(EpicsSignalRO, "FilterMask", kind="config")
    mask_setpoint = Component(
        EpicsSignalRO, "FilterMaskSetpoint", kind="config"
    )

    message = Component(EpicsSignalRO, "Message", kind="config")

    slots = DynamicDeviceComponent(make_filter_slots(NUM_FILTERS))

    # Configuration
    wait_time = Component(EpicsSignal, "WaitTime", kind="config")
    debug_level = Component(EpicsSignal, "Debug", kind="config")
