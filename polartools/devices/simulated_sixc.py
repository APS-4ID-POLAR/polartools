"""
Simulated psic
"""

__all__ = ["psic"]

from ophyd import Component, PseudoSingle, Kind, Signal, EpicsMotor
from ..utils.run_engine import sd
import gi

gi.require_version("Hkl", "5.0")
# MUST come before `import hkl`
from hkl.geometries import E6C
from hkl.user import select_diffractometer
from ..utils import logger

logger.info(__file__)


class SixCircleDiffractometer(E6C):
    """
    E6C: huber diffractometer in 6-circle geometry with energy.

    4-ID-D setup.
    """

    # HKL and 4C motors
    h = Component(PseudoSingle, "", labels=("hkl", "psic"))
    k = Component(PseudoSingle, "", labels=("hkl", "psic"))
    l = Component(PseudoSingle, "", labels=("hkl", "psic"))

    mu = Component(EpicsMotor, "m9", labels=("motor", "psic"))
    omega = Component(EpicsMotor, "m10", labels=("motor", "psic"))
    chi = Component(EpicsMotor, "m11", labels=("motor", "psic"))
    phi = Component(EpicsMotor, "m12", labels=("motor", "psic"))
    delta = Component(EpicsMotor, "m13", labels=("motor", "psic"))
    gamma = Component(EpicsMotor, "m14", labels=("motor", "psic"))

    # Explicitly selects the real motors
    # _real = ['theta', 'chi', 'phi', 'tth']
    _real = "mu omega chi phi delta gamma".split()

    # Energy
    energy = Component(Signal, value=8)
    energy_update_calc_flag = Component(Signal, value=1)
    energy_offset = Component(Signal, value=0)

    # TODO: This is needed to prevent busy plotting.
    @property
    def hints_test(self):
        fields = []
        for _, component in self._get_components_of_kind(Kind.hinted):
            if (~Kind.normal & Kind.hinted) & component.kind:
                c_hints = component.hints
                fields.extend(c_hints.get("fields", []))
        return {"fields": fields}


psic = SixCircleDiffractometer("4idsoftmotors:", name="psic")


class SixcPSI(E6C):
    """
    Our 6-circle.  Eulerian.
    """

    # the reciprocal axes are called "pseudo" in hklpy
    psi = Component(PseudoSingle, "")

    # the motor axes are called "real" in hklpy
    mu = Component(EpicsMotor, "m9", labels=("motor", "sixcpsi"))
    omega = Component(EpicsMotor, "m10", labels=("motor", "sixcpsi"))
    chi = Component(EpicsMotor, "m11", labels=("motor", "sixcpsi"))
    phi = Component(EpicsMotor, "m12", labels=("motor", "sixcpsi"))
    delta = Component(EpicsMotor, "m13", labels=("motor", "sixcpsi"))
    gamma = Component(EpicsMotor, "m14", labels=("motor", "sixcpsi"))


sixcpsi = SixcPSI("", name="sixcpsi", engine="psi")

select_diffractometer(psic)
# sd.baseline.append(psic)
