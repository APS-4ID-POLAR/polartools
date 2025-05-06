"""
Simulated diffract
"""

__all__ = ['diffract']

from ophyd import Component, PseudoSingle, Kind, Signal, EpicsMotor
from ..utils.run_engine import sd
import gi
gi.require_version('Hkl', '5.0')
# MUST come before `import hkl`
from hkl.geometries import Petra3_p09_eh2
from hkl.user import select_diffractometer
from ..utils import logger
logger.info(__file__)


class FourCircleDiffractometer(Petra3_p09_eh2):
    """
    E4CV: huber diffractometer in 4-circle vertical geometry with energy.

    4-ID-D setup.
    """

    # HKL and 4C motors
    h = Component(PseudoSingle, '', labels=("hkl", "diffract"))
    k = Component(PseudoSingle, '', labels=("hkl", "diffract"))
    l = Component(PseudoSingle, '', labels=("hkl", "diffract"))

    mu = Component(EpicsMotor, "m9", labels=("motor", "diffract"))
    omega = Component(EpicsMotor, "m10", labels=("motor", "diffract"))
    chi = Component(EpicsMotor, "m11", labels=("motor", "diffract"))
    phi = Component(EpicsMotor, "m12", labels=("motor", "diffract"))
    delta = Component(EpicsMotor, "m13", labels=("motor", "diffract"))
    gamma = Component(EpicsMotor, "m14", labels=("motor", "diffract"))

    # Explicitly selects the real motors
    # _real = ['theta', 'chi', 'phi', 'tth']
    _real = "mu omega chi phi delta gamma".split()

    # Energy
    energy = Component(Signal, value=8)
    energy_update_calc_flag = Component(Signal, value=1)
    energy_offset = Component(Signal, value=0)

    # TODO: This is needed to prevent busy plotting.
    @property
    def hints(self):
        fields = []
        for _, component in self._get_components_of_kind(Kind.hinted):
            if (~Kind.normal & Kind.hinted) & component.kind:
                c_hints = component.hints
                fields.extend(c_hints.get('fields', []))
        return {'fields': fields}


diffract = FourCircleDiffractometer("4idsoftmotors:", name='diffract')
select_diffractometer(diffract)
# sd.baseline.append(diffract)
