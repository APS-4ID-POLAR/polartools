"""
Transfocator
"""

__all__ = ['transfocator']

from ophyd import (
    Device,
    Component,
    DynamicDeviceComponent,
    FormattedComponent,
    PVPositioner,
    EpicsMotor,
    Signal,
    EpicsSignal,
    EpicsSignalRO,
    EpicsSignalWithRBV,
    DeviceStatus
)
from ophyd.status import AndStatus
from bluesky.plan_stubs import mv, sleep as bps_sleep, abs_set
from apstools.devices import TrackingSignal
from toolz import partition
from numpy import poly1d, loadtxt
from scipy.interpolate import interp1d
from time import sleep as tsleep
from .monochromator import mono
from ..utils._logging_setup import logger
from ..utils.transfocator_calculation_new import transfocator_calculation

logger.info(__file__)

MOTORS_IOC = "4idgSoft:"
EPICS_ENERGY_SLEEP = 0.15


def _make_lenses_motors(motors: list):
    defn = {}
    for n, mot in enumerate(motors):
        defn[f"l{n}"] = (
            EpicsMotor, f"{mot}", dict(kind="config", labels=("motor",))
        )
    return defn


class PyCRLSingleLens(PVPositioner):
    readback = Component(EpicsSignalRO, "_RBV")
    setpoint = Component(EpicsSignal, "", put_complete=True)

    done = Component(EpicsSignal, "_matchCalc.C")
    done_value = 1

    # Settings
    num_lenses = Component(EpicsSignal, "_NumLens", kind="config")
    radius = Component(EpicsSignal, "_LensRadius", kind="config")
    location = Component(EpicsSignal, "_Location", kind="config")
    material = Component(EpicsSignal, "_Material", kind="config")
    thickness_error = Component(EpicsSignal, "_ThickErr", kind="config")
    in_limit = Component(EpicsSignal, "_RBV_calc.CC", kind="config")

    def set(
        self,
        new_position,
        *,
        timeout: float = None,
        moved_cb=None,
        wait: bool = False,
    ):
        if self.readback.get() == new_position:
            _status = DeviceStatus(self)
            _status.set_finished()
            return _status
        else:
            return super().set(
                new_position, timeout=timeout, moved_cb=moved_cb, wait=wait
            )


class PyCRLSignal(EpicsSignal):
    value = Component(EpicsSignal, "")
    egu = Component(EpicsSignalRO, ".EGU")


class PyCRL(Device):

    # Energy
    energy_mono = Component(PyCRLSignal, "EnergyBeamline", kind="config")
    energy_local = Component(PyCRLSignal, "EnergyLocal", kind="config")
    energy_select = Component(PyCRLSignal, "EnergySelect", kind="config")

    # Slits
    slit_hor_size = Component(PyCRLSignal, "1:slitSize_H_RBV", kind="config")
    slit_hor_pv = Component(
        EpicsSignal, "1:slitSize_H.DOL", string=True, kind="config"
    )
    slit_vert_size = Component(PyCRLSignal, "1:slitSize_V_RBV", kind="config")
    slit_vert_pv = Component(
        EpicsSignal, "1:slitSize_V.DOL", string=True, kind="config"
    )

    # Focus info/control
    focal_size_setpoint = Component(EpicsSignal, "focalSize")
    focal_size_readback = Component(EpicsSignalRO, "fSize_actual")
    focal_power_index = Component(EpicsSignalWithRBV, "1:sortedIndex")
    # focal_power_index_readback = Component(EpicsSignal, "1:sortedIndex_RBV")
    focal_sizes = Component(EpicsSignal, "fSizes", kind="omitted")
    minimize_button = Component(EpicsSignal, "minimizeFsize.PROC", kind="omitted")
    system_done = Component(EpicsSignalRO, "sysBusy", kind="omitted")

    # Parameters readbacks
    dq = Component(PyCRLSignal, "dq", kind="config")
    q = Component(PyCRLSignal, "q", kind="config")
    z_offset = Component(PyCRLSignal, "1:oePositionOffset_RBV", kind="config")
    z_offset_pv = Component(
        EpicsSignal, "1:oePositionOffset.DOL", kind="config"
    )
    z_from_source = Component(PyCRLSignal, "1:oePosition_RBV", kind="config")
    sample_offset = Component(
        PyCRLSignal, "samplePositionOffset_RBV", kind="config"
    )
    sample_offset_pv = Component(
        EpicsSignal, "samplePositionOffset.DOL", kind="config"
    )
    sample = Component(PyCRLSignal, "samplePosition_RBV", kind="config")

    # Lenses indices
    binary = Component(EpicsSignalRO, "1:lenses", kind="config")
    ind_control = Component(EpicsSignalRO, "1:lensConfig_BW", kind="config")
    readbacks = Component(EpicsSignalRO, "1:lensConfig_RBV", kind="config")

    # Other options
    preview_index = Component(EpicsSignal, "previewIndex", kind="config")
    focal_size_preview = Component(
        EpicsSignalRO, "fSize_preview", kind="config"
    )
    inter_lens_delay = Component(EpicsSignal, "1:interLensDelay", kind="config")
    verbose_console = Component(EpicsSignal, "verbosity", kind="config")
    thickness_error_flag = Component(
        EpicsSignal, "thickerr_flag", kind="config"
    )
    beam_mode = Component(EpicsSignalWithRBV, "beamMode", kind="config")

    # Lenses
    lens1 = Component(PyCRLSingleLens, "1:stack01")
    lens2 = Component(PyCRLSingleLens, "1:stack02")
    lens3 = Component(PyCRLSingleLens, "1:stack03")
    lens4 = Component(PyCRLSingleLens, "1:stack04")
    lens5 = Component(PyCRLSingleLens, "1:stack05")
    lens6 = Component(PyCRLSingleLens, "1:stack06")
    lens7 = Component(PyCRLSingleLens, "1:stack07")
    lens8 = Component(PyCRLSingleLens, "1:stack08")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._status = None
        self.system_done.subscribe(self._update_status_subscription, run=False)
    
    def _update_status_subscription(self, value, old_value, **kwarg):
        if (
            (self._status is not None) and
            (value in ["Done", 0]) and
            (old_value in ["Changing", 1])
        ):
            self._status.set_finished()
            self._status = None

    # def minimize_beam(self, value=1, **kwargs):
    #     _button_status = self.minimize_button.set(value)
    #     self._status = DeviceStatus(self)
    #     return AndStatus(_button_status, self._status)

    def set(self, value, **kwargs):
        _st = DeviceStatus(self)

        if self.system_done.get() in ["Done", 0]:
            _st.set_finished()
        else:
            self._status = _st
        
        return _st


class EnergySignal(Signal):

    _epics_sleep = EPICS_ENERGY_SLEEP

    def put(self, *args, **kwargs):
        raise NotImplementedError("put operation not setup in this signal.")

    def set(self, value, **kwargs):

        self._readback = value

        if self.parent.energy_select.get() != 1:
            self.parent.energy_select.set(1).wait(1)

        self.parent.energy_local.set(value).wait(1)
        tsleep(self._epics_sleep)
        # this is needed because the scan of the transfocator is 0.1 s

        zpos = self.parent.z.user_readback.get() - self.parent.dq.get()*1000.
        # dq in meters

        return self.parent.z.set(zpos, **kwargs)


class ZMotor(EpicsMotor):
    def set(self, new_position, **kwargs):

        zstatus = super().set(new_position, **kwargs)

        if self.parent.trackxy.get():
            if self.parent._x_interpolation is None:
                raise ValueError(
                    "The reference data for X tracking has not been entered. "
                    "Cannot track the X motion."
                )

            if self.parent._y_interpolation is None:
                raise ValueError(
                    "The reference data for Y tracking has not been entered. "
                    "Cannot track the Y motion."
                )

            xpos = (
                self.parent._x_interpolation(new_position) + self.parent.deltax.get()
            )
            ypos = (
                self.parent._y_interpolation(new_position) + self.parent.deltay.get()
            )

            xystatus = AndStatus(
                self.parent.x.set(xpos),
                self.parent.y.set(ypos)
            )

            return AndStatus(zstatus, xystatus)
        else:
            return zstatus

    def stop(self, *, success=False):
        super().stop(success=success)
        if self.parent.trackxy.get():
            self.parent.x.stop(success=success)
            self.parent.y.stop(success=success)


class TransfocatorClass(PyCRL):

    energy = Component(EnergySignal)
    tracking = Component(TrackingSignal, value=False, kind="config")

    # Motors -- setup in 4idgSoft
    x = FormattedComponent(EpicsMotor, "{_motors_IOC}m58", labels=("motor",))
    y = FormattedComponent(EpicsMotor, "{_motors_IOC}m57", labels=("motor",))
    z = FormattedComponent(ZMotor, "{_motors_IOC}m61", labels=("motor",))
    pitch = FormattedComponent(
        EpicsMotor, "{_motors_IOC}m60", labels=("motor",)
    )
    yaw = FormattedComponent(EpicsMotor, "{_motors_IOC}m59", labels=("motor",))

    lens_motors = DynamicDeviceComponent(
        _make_lenses_motors(
            [
                f"{MOTORS_IOC}m69",
                f"{MOTORS_IOC}m68",
                f"{MOTORS_IOC}m67",
                f"{MOTORS_IOC}m66",
                f"{MOTORS_IOC}m65",
                f"{MOTORS_IOC}m64",
                f"{MOTORS_IOC}m63",
                f"{MOTORS_IOC}m62"
            ]
        ),
        component_class=FormattedComponent
    )

    reference_data_x = Component(Signal, kind="config")
    reference_data_y = Component(Signal, kind="config")
    deltax = Component(Signal, value=0, kind="config")
    deltay = Component(Signal, value=0, kind="config")
    trackxy = Component(TrackingSignal, value=False, kind="config")

    def __init__(
            self,
            *args,
            lens_pos=30,
            default_distance=2591,
            # reference_x=0,
            # reference_y=0,
            # x_polynomial=[0],
            # y_polynomial=[0],
            **kwargs
    ):
        self._motors_IOC = MOTORS_IOC
        PyCRL.__init__(self, *args, **kwargs)
        self._lens_pos = lens_pos
        self._default_distance = default_distance  # mm
        # self.reference_x.put(reference_x)
        # self.reference_y.put(reference_y)
        # self.polynomial_x.put(x_polynomial)
        # self.polynomial_y.put(y_polynomial)
        self._x_interpolation = None
        self._y_interpolation = None
        self.reference_data_x.subscribe(self._update_interpolation_x, run=False)
        self.reference_data_y.subscribe(self._update_interpolation_y, run=False)

    def load_reference_data(self, fname, axis):
        if axis not in "x y".split():
            raise ValueError(f"axis must be x or y. {axis} is not valid.")
        # x, y = loadtxt(fname, unpack=True)
        getattr(self, f"reference_data_{axis}").put(loadtxt(fname))

    def _update_interpolation_x(self, value, **kwargs):
        z = value[:, 0]
        x = value[:, 1]
        self._x_interpolation = interp1d(z, x)

    def _update_interpolation_y(self, value, **kwargs):
        z = value[:, 0]
        y = value[:, 1]
        self._y_interpolation = interp1d(z, y)

    def lens_status(self, i):
        return getattr(self, f"lens{i}").readback.get(as_string=True)

    @property
    def lenses_in(self):
        selected = []
        for i in range(1, 9):
            _status = self.lens_status(i)
            if _status == "In":
                selected.append(i)
            elif _status == "Both out":
                pass
                # logger.info(f"WARNING: the status of lens #{i} is unknown.")
        return selected

    def _setup_lenses_move(self, lenses_in: list = []):
        """
        Adjust lenses

        PARAMETERS
        ----------
        lenses_in : list or iterable
            Index of the lenses that will be inserted. The ones not in this list
            will be removed.
        type : "plan" or "noplan"
            Determines how the lenses will be used, using a bluesky plan
            ("plan" option), or "noplan".
        """

        for i in lenses_in:
            if (i > 8) or (i < 1):
                raise ValueError("Lens index must be from 1 to 8.")

        # Positive/negative step moves lens in/out respectively.
        # We want to move it to the hard limit.

        args = []
        for lens in range(1, 9):
            step = 1 if lens in lenses_in else 0
            args += [
                getattr(self, f"lens{lens}"), step
            ]

        return args

    def set_lenses(self, selected_lenses: list):
        args = self._setup_lenses_move(selected_lenses)
        for dev, pos in partition(2, args):
            dev.setpoint.put(pos)

    def set_lenses_plan(self, selected_lenses: list):
        args = self._setup_lenses_move(selected_lenses)
        return (yield from mv(*args))

    def _check_z_lims(self, position):
        if (
            (position > self.z.low_limit_travel.get()) &
            (position < self.z.high_limit_travel.get())
        ):
            return True
        else:
            return False

    def _setup_optimize_distance(self):

        if self.energy_select.get() in (1, "Local"):
            logger.info("WARNING: transfocator in 'Local' energy mode")

        distance = (
            self.z.user_readback.get() - self.dq.get()*1000
        )

        if not self._check_z_lims(distance):
            raise ValueError(
                f"The distance {distance} is outsize the Z travel range. No"
                "motion will occur."
            )

        return distance

    def optimize_lenses(self,):

        self.focal_power_index.set(
            self.focal_sizes.get().argmin()
        ).wait()

        self.z.move(
            self._setup_optimize_distance()
        ).wait()

        self.set(1).wait()

    def optimize_lenses_plan(self):

        def _moves():
            yield from mv(
                self.focal_power_index,
                self.focal_sizes.get().argmin()
            )
            yield from mv(
                self.z, self._setup_optimize_distance(),
                self, 1
            )

        return (yield from _moves())

    def optimize_distance(self):
        self.z.move(
            self._setup_optimize_distance()
        ).wait()

    def optimize_distance_plan(self):
        return (
            yield from mv(
                self.z,
                self._setup_optimize_distance()
            )
        )

    def calc(
        self,
        optimize_position=None,
        reference_distance=None,
        energy=None,
        experiment="diffractometer",
        distance_only=False,
        selected_lenses=None,
        verbose=True
    ):
        if energy is None:
            energy = mono.energy.get()

        if selected_lenses is None:
            selected_lenses = self.lenses_in

        if reference_distance is None:
            reference_distance = self._default_distance

        if optimize_position is None:
            optimize_position = 0

        return transfocator_calculation(
            energy,
            optimize_position=optimize_position,
            reference_distance=reference_distance,
            experiment=experiment,
            distance_only=distance_only,
            selected_lenses=selected_lenses,
            verbose=verbose
        )

    def move_z_correct_xy_plan(self, zpos):
        xpos = (
            self.reference_x.get() +
            poly1d(self.polynomial_x.get())(zpos)
        )
        ypos = (
            self.reference_y.get() +
            poly1d(self.polynomial_y.get())(zpos)
        )

        yield from mv(
            self.x, xpos,
            self.y, ypos,
            self.z, zpos
        )


transfocator = TransfocatorClass(
    "4idPyCRL:CRL4ID:",
    name="transfocator",
    labels=("4idg", "optics", "track_energy")
)

transfocator.stage_sigs["energy_select"] = 1
