
"""
HHL mirror
"""

from ophyd import (
    Component,
    FormattedComponent,
    Device,
    EpicsMotor,
    EpicsSignal,
    EpicsSignalRO
)


class ToroidalMirror(Device):
    """ Beamline toroidal mirror components. """

    # Motors
    y = Component(EpicsMotor, 'm1', labels=('motor',))
    x1 = Component(EpicsMotor, 'm2', labels=('motor',))
    x2 = Component(EpicsMotor, 'm3', labels=('motor',))
    us_bend = Component(EpicsMotor, 'm4', labels=('motor',))
    ds_bend = Component(EpicsMotor, 'm5', labels=('motor',))

    # Combined motions
    x = Component(EpicsMotor, 'pm1', labels=('motor',))
    pitch = Component(EpicsMotor, 'pm2', labels=('motor',))
    fine_pitch = FormattedComponent(
        EpicsMotor, '4idaSoft:m1', labels=('motor',)
    )
    curvature = Component(EpicsMotor, 'pm3', labels=('motor',))
    elipticity = Component(EpicsMotor, 'pm4', labels=('motor',))

    # Other parameters
    stripe = Component(EpicsSignal, 'stripe', string=True)
    radius_target = Component(EpicsSignalRO, 'EstimatedRoC')
    critical_energy = Component(EpicsSignalRO, 'Ecritical')
    beam_offset = Component(EpicsSignalRO, 'beam_offset')
    alpha = Component(EpicsSignalRO, 'alpha')
