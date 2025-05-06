"""
Simulated polar
"""

__all__ = ["huber_euler", "huber_hp", "huber_euler_psi", "huber_hp_psi"]

from ophyd import (
    Component,
    FormattedComponent,
    PseudoSingle,
    PseudoPositioner,
    Kind,
    Signal,
    EpicsMotor,
    EpicsSignal,
    EpicsSignalRO,
)
from ophyd.pseudopos import pseudo_position_argument, real_position_argument
from scipy.constants import speed_of_light, Planck
from numpy import arcsin, pi, sin, tan
from .jj_slits import SlitDevice
from .huber_filter import HuberFilter
from ..utils.analyzer_utils import check_structure_factor, calcdhkl

from ..utils import logger
import gi

gi.require_version("Hkl", "5.0")
# MUST come before `import hkl`
from hkl.geometries import ApsPolar
from hkl.user import select_diffractometer
import math

logger.info(__file__)

# Constants
WAVELENGTH_CONSTANT = 12.39
PTTH_MIN_DEGREES = 79
PTTH_MAX_DEGREES = 101
PTH_MIN_DEGREES = 39
PTH_MAX_DEGREES = 51
# TODO: Find a better way to get to this
ANALYZER_LIST_PATH = (
    "/home/beams17/POLAR/joerg/polar_instrument/src/instrument/devices/"
    "analyzerlist.dat"
)


class AnalyzerDevice(PseudoPositioner):

    energy = Component(PseudoSingle, limits=(2.6, 34))
    th = Component(EpicsMotor, "pmth", labels=("motor",))
    tth_trans = Component(EpicsMotor, "m25", labels=("motor",))

    _real = ["th", "tth_trans"]

    # Analyzer motors
    th_motor = Component(EpicsMotor, "m24", labels=("motor",))
    tth = Component(EpicsMotor, "pm2th", labels=("motor",))
    eta = Component(EpicsMotor, "m23", labels=("motor",))
    chi = Component(EpicsMotor, "m26", labels=("motor",))

    d_spacing = Component(Signal, value=1e4, kind="config")
    crystal = Component(Signal, value="None", kind="config")
    tth_detector_distance = Component(
        EpicsSignal, "TWTH:Drive.C", kind="config"
    )
    tracking = Component(Signal, value=False, kind="config")

    def move_single(self, pseudo, position, **kwargs):
        if self.d_spacing.get() == 1e4:
            raise RuntimeError(
                "The analyzer has not been setup, please run the .setup() "
                "function before moving the energy"
            )
        return super().move_single(pseudo, position, **kwargs)

    # These assume that the analyzer is part of the diffractometer.
    @property
    def beamline_wavelength(self):
        return self.parent.calc.wavelength

    @property
    def beamline_energy(self):
        return self.parent.energy.get()

    def convert_energy_to_theta(self, energy):
        # lambda in angstroms, theta in degrees, energy in keV
        lamb = speed_of_light * Planck * 6.241509e15 * 1e10 / energy
        theta = arcsin(lamb / 2 / self.d_spacing.get()) * 180.0 / pi
        return theta

    def convert_energy_to_tth_trans(self, energy):
        # lambda in angstroms, theta in degrees, energy in keV
        th = self.convert_energy_to_theta(energy)
        tth = 2 * th
        tth_trans = self.tth_detector_distance.get() * tan(
            (tth - 45 - th) * pi / 180.0
        )
        return tth_trans

    def convert_theta_to_energy(self, theta):
        # lambda in angstroms, theta in degrees, energy in keV
        lamb = 2 * self.d_spacing.get() * sin(theta * pi / 180)
        energy = speed_of_light * Planck * 6.241509e15 * 1e10 / lamb
        return energy

    @pseudo_position_argument
    def forward(self, pseudo_pos):
        """Run a forward (pseudo -> real) calculation"""
        return self.RealPosition(
            th=self.convert_energy_to_theta(pseudo_pos.energy),
            tth_trans=self.convert_energy_to_tth_trans(pseudo_pos.energy),
        )

    @real_position_argument
    def inverse(self, real_pos):
        """Run an inverse (real -> pseudo) calculation"""
        # Changing y does not change the energy.
        return self.PseudoPosition(
            energy=self.convert_theta_to_energy(real_pos.th)
        )

    def set_energy(self, energy):
        # energy in keV, theta in degrees.
        theta = self.convert_energy_to_theta(energy)
        self.th.set_current_position(theta)

    def calc(self):
        d_ana = self.d_spacing.get()
        if d_ana == 1e4:
            self.setup()
            d_ana = self.d_spacing.get()
        wavelength = self.beamline_wavelength
        cryst = self.crystal.get()
        th_angle = math.degrees(math.asin(wavelength / (2 * d_ana)))
        tth_angle = 2 * th_angle
        print(
            f"[ath, atth] = [{th_angle:6.2f}, {tth_angle:6.2f}] for {cryst} "
            "analyzer"
        )

    def setup(
        self, analyzer_energy=None, analyzer_list_path=ANALYZER_LIST_PATH
    ):
        if not analyzer_energy:
            energy = self.beamline_energy
            wavelength = self.beamline_wavelength
        else:
            energy = analyzer_energy
            wavelength = WAVELENGTH_CONSTANT / energy
        ptthmin = math.radians(PTTH_MIN_DEGREES)
        ptthmax = math.radians(PTTH_MAX_DEGREES)

        d_dict = {}
        print(
            f"{'#':>2}{'Crystal':>8}{'Refl.':>11}{'Range (keV)':>15}"
            f"{'d_hkl':>10}{'2th (deg)':>11}"
        )
        with open(analyzer_list_path, "r") as f:
            for num, line in enumerate(f):
                split = line.split()
                if "##REFERENCES" in line:
                    break
                (
                    analyzer,
                    hh,
                    kk,
                    ll,
                    a,
                    b,
                    c,
                    alpha,
                    beta,
                    gamma,
                    symmetry,
                    spacegroupnumber,
                    special,
                ) = split[:13]
                hh, kk, ll = map(int, [hh, kk, ll])
                a, b, c = map(float, [a, b, c])
                alpha, beta, gamma = map(
                    lambda x: math.radians(float(x)), [alpha, beta, gamma]
                )
                spacegroupnumber = int(spacegroupnumber)

                for i in range(
                    1, 100
                ):  # Arbitrary limit to prevent infinite loop
                    hhh, kkk, lll = hh * i, kk * i, ll * i
                    dhkl = calcdhkl(
                        hhh, kkk, lll, alpha, beta, gamma, symmetry, a, b, c
                    )
                    if check_structure_factor(
                        hhh, kkk, lll, spacegroupnumber, special
                    ):
                        ana_emax = WAVELENGTH_CONSTANT / (
                            2 * dhkl * math.sin(ptthmin / 2)
                        )
                        ana_emin = WAVELENGTH_CONSTANT / (
                            2 * dhkl * math.sin(ptthmax / 2)
                        )

                        if ana_emin <= energy <= ana_emax:
                            tt_angle = 2 * math.degrees(
                                math.asin(wavelength / (2 * dhkl))
                            )
                            print(
                                f"{num:>2} {analyzer:<9} {hhh:>2}{kkk:>3}"
                                f"{lll:>3}  [{ana_emin:5.2f},{ana_emax:5.2f}]"
                                f"{dhkl:10.3f}{tt_angle:11.2f}"
                            )
                            d_dict[num] = [
                                analyzer,
                                hhh,
                                kkk,
                                lll,
                                dhkl,
                                tt_angle,
                            ]
                            break
                        elif ana_emin > energy:
                            break
        ttdiff = max(abs(ptthmax - 90), abs(ptthmin - 90))
        d_best = {}
        for key, value in d_dict.items():
            if abs(value[5] - 90) < ttdiff:
                ttdiff = abs(value[5] - 90)
                d_best = [key, value]
        if analyzer_energy:
            print(
                f"Best analyzer to use at {energy}: {d_best[1][0]}_"
                f"{d_best[1][1]}{d_best[1][2]}{d_best[1][3]}"
            )
        else:
            anum = (
                input(f"Choice of polarization analyzer [{d_best[0]}]: ")
                or d_best[0]
            )
            anum = int(anum) if isinstance(anum, str) else anum
            if anum in d_dict:
                ana = d_dict[anum]
                cryst = f"{ana[0]}_{ana[1]}{ana[2]}{ana[3]}"
            else:
                ana = d_best[1:][0]
                cryst = f"{ana[0]}_{ana[1]}{ana[2]}{ana[3]}"
                print(f"Choice not possible, using {cryst}")

            d_ana = ana[4]
            self.d_spacing.put(d_ana)
            self.crystal.put(cryst)


class SixCircleDiffractometer(ApsPolar):
    """
    ApsPolar: Huber diffractometer in 6-circle horizontal geometry with energy.

    HKL engine.
    """

    # HKL and 6C motors
    h = Component(PseudoSingle, "", labels=("hkl",))
    k = Component(PseudoSingle, "", labels=("hkl",))
    l = Component(PseudoSingle, "", labels=("hkl",))

    # 03/16/2025 - Tau is the whole diffractometer "theta" angle, but
    # it is not currently setup. m73 is a simulated motor.
    tau = Component(EpicsMotor, "m73", labels=("motor",))
    mu = Component(EpicsMotor, "m4", labels=("motor",))
    gamma = Component(EpicsMotor, "m19", labels=("motor",))
    delta = Component(EpicsMotor, "m20", labels=("motor",))

    # Explicitly selects the real motors
    _real = "tau mu chi phi gamma delta".split()

    # Table vertical/horizontal
    tablex = Component(EpicsMotor, "m3", labels=("motor",))
    tabley = Component(EpicsMotor, "m1", labels=("motor",))

    # Area detector motors
    pad_rail = Component(EpicsMotor, "m21", labels=("motor",))
    point_rail = Component(EpicsMotor, "m22", labels=("motor",))

    # # Guard slit
    # guardslt  = ...

    # Filters
    filter = Component(HuberFilter, "atten:", labels=("filter",))

    # Detector JJ slit
    detslt = Component(
        SlitDevice,
        "",
        motorsDict={"top": "m31", "bot": "m32", "out": "m34", "inb": "m33"},
        slitnum=2,
        labels=("slit",),
    )

    # Analyzer
    ana = Component(AnalyzerDevice, "", labels=("track_energy",))

    # Energy
    energy = FormattedComponent(
        EpicsSignalRO, "4idVDCM:BraggERdbkAO", kind="config"
    )
    energy_update_calc_flag = Component(Signal, value=1, kind="config")
    energy_offset = Component(Signal, value=0, kind="config")

    # TODO: This is needed to prevent busy plotting.
    @property
    def hints(self):
        fields = []
        for _, component in self._get_components_of_kind(Kind.hinted):
            if (~Kind.normal & Kind.hinted) & component.kind:
                c_hints = component.hints
                fields.extend(c_hints.get("fields", []))
        return {"fields": fields}


class CradleDiffractometer(SixCircleDiffractometer):
    chi = Component(EpicsMotor, "m37", labels=("motor",))
    phi = Component(EpicsMotor, "m38", labels=("motor",))

    x = Component(EpicsMotor, "m40", labels=("motor",))
    y = Component(EpicsMotor, "m41", labels=("motor",))
    z = Component(EpicsMotor, "m42", labels=("motor",))


class HPDiffractometer(SixCircleDiffractometer):
    chi = Component(EpicsMotor, "m5", labels=("motor",))
    phi = Component(EpicsMotor, "m6", labels=("motor",))

    basex = Component(EpicsMotor, "m7", labels=("motor",))
    basey = Component(EpicsMotor, "m9", labels=("motor",))
    basez = Component(EpicsMotor, "m8", labels=("motor",))

    sample_tilt = Component(EpicsMotor, "m11", labels=("motor",))

    x = Component(EpicsMotor, "m12", labels=("motor",))
    y = Component(EpicsMotor, "m14", labels=("motor",))
    z = Component(EpicsMotor, "m13", labels=("motor",))

    nanox = FormattedComponent(
        EpicsMotor, "4idgSoftX:jena:m1", labels=("motors",)
    )
    nanoy = FormattedComponent(
        EpicsMotor, "4idgSoftX:jena:m2", labels=("motors",)
    )
    nanoz = FormattedComponent(
        EpicsMotor, "4idgSoftX:jena:m3", labels=("motors",)
    )


class PolarPSI(ApsPolar):
    """
    ApsPolar: Huber diffractometer in 6-circle horizontal geometry with energy.

    Psi engine.
    """

    # the reciprocal axes are called "pseudo" in hklpy
    psi = Component(PseudoSingle, "")

    # 03/16/2025 - Tau is the whole diffractometer "theta" angle, but
    # it is not currently setup. m73 is a simulated motor.
    tau = Component(EpicsMotor, "m73", labels=("motor",))
    mu = Component(EpicsMotor, "m4", labels=("motor",))
    gamma = Component(EpicsMotor, "m19", labels=("motor",))
    delta = Component(EpicsMotor, "m20", labels=("motor",))


class CradlePSI(PolarPSI):
    chi = Component(EpicsMotor, "m37", labels=("motor",))
    phi = Component(EpicsMotor, "m38", labels=("motor",))


class HPPSI(PolarPSI):
    chi = Component(EpicsMotor, "m5", labels=("motor",))
    phi = Component(EpicsMotor, "m6", labels=("motor",))


huber_euler = CradleDiffractometer(
    "4idgSoft:",
    name="huber_euler",
    labels=(
        "4idg",
        "diffractometer",
    ),
)


huber_hp = HPDiffractometer(
    "4idgSoft:",
    name="huber_hp",
    labels=(
        "4idg",
        "diffractometer",
    ),
)

huber_euler_psi = CradlePSI(
    "4idgSoft:",
    name="huber_euler_psi",
    engine="psi",
    labels=(
        "4idg",
        "diffractometer",
    ),
)

huber_hp_psi = CradlePSI(
    "4idgSoft:",
    name="huber_hp_psi",
    engine="psi",
    labels=(
        "4idg",
        "diffractometer",
    ),
)

select_diffractometer(huber_euler)
huber_euler._update_calc_energy()
huber_euler_psi._update_calc_energy()
