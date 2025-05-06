"""
Monochromator with energy controller by bluesky
"""

__all__ = ["mono"]

from ophyd import (
    Component,
    FormattedComponent,
    EpicsMotor,
    EpicsSignal,
    PseudoPositioner,
    PseudoSingle,
)
from ophyd.pseudopos import pseudo_position_argument, real_position_argument
from scipy.constants import speed_of_light, Planck
from numpy import arcsin, pi, sin, cos
from .labjacks import AnalogOutput
from ..utils._logging_setup import logger

logger.info(__file__)


class MonoDevice(PseudoPositioner):

    energy = Component(PseudoSingle, limits=(2.6, 32))
    th = Component(EpicsMotor, 'm1', labels=('motor',))

    y2 = Component(EpicsMotor, 'm3', labels=('motor',))

    # Explicitly selects the real motors
    _real = ['th', 'y2']

    # Other motors
    crystal_select = Component(EpicsMotor, 'm2', labels=('motor',), kind="config")
    thf2 = Component(EpicsMotor, 'm4', labels=('motor',))
    chi2 = Component(EpicsMotor, 'm5', labels=('motor',))

    # PZTs from labjack
    pzt_thf2 = FormattedComponent(AnalogOutput, "4idaSoft:LJ:Ao5")
    pzt_chi2 = FormattedComponent(AnalogOutput, "4idaSoft:LJ:Ao3")

    # Parameters
    y_offset = Component(EpicsSignal, "Kohzu_yOffsetAO.VAL", kind="config")
    crystal_h = Component(EpicsSignal, "BraggHAO.VAL", kind="config")
    crystal_k = Component(EpicsSignal, "BraggKAO.VAL", kind="config")
    crystal_l = Component(EpicsSignal, "BraggLAO.VAL", kind="config")
    crystal_a = Component(EpicsSignal, "BraggAAO.VAL", kind="config")
    crystal_2d = Component(EpicsSignal, "Bragg2dSpacingAO", kind="config")
    crystal_type = Component(
        EpicsSignal, "BraggTypeMO", string=True, kind="config"
    )

    def convert_energy_to_theta(self, energy):
        # lambda in angstroms, theta in degrees, energy in keV
        lamb = speed_of_light*Planck*6.241509e15*1e10/energy
        theta = arcsin(lamb/self.crystal_2d.get())*180./pi
        return theta

    def convert_energy_to_y(self, energy):
        # lambda in angstroms, theta in degrees, energy in keV
        theta = self.convert_energy_to_theta(energy)
        return self.y_offset.get()/(2*cos(theta*pi/180))

    def convert_theta_to_energy(self, theta):
        # lambda in angstroms, theta in degrees, energy in keV
        lamb = self.crystal_2d.get()*sin(theta*pi/180)
        energy = speed_of_light*Planck*6.241509e15*1e10/lamb
        return energy

    @pseudo_position_argument
    def forward(self, pseudo_pos):
        '''Run a forward (pseudo -> real) calculation'''
        return self.RealPosition(
            th=self.convert_energy_to_theta(pseudo_pos.energy),
            y2=self.convert_energy_to_y(pseudo_pos.energy)
        )

    @real_position_argument
    def inverse(self, real_pos):
        '''Run an inverse (real -> pseudo) calculation'''
        # Changing y does not change the energy.
        return self.PseudoPosition(
            energy=self.convert_theta_to_energy(real_pos.th)
        )

    def set_energy(self, energy):
        # energy in keV, theta in degrees.
        theta = self.convert_energy_to_theta(energy)
        self.th.set_current_position(theta)


mono = MonoDevice(
    "4idVDCM:", name="mono", labels=("monochromator", "energy", "4ida")
)

mono._sub_devices.remove("crystal_select")
