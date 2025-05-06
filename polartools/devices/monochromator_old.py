"""
Monochromator motors
"""

__all__ = ['mono']

from apstools.devices import KohzuSeqCtl_Monochromator
from ophyd import (
    Component,
    Device,
    FormattedComponent,
    EpicsMotor,
    EpicsSignal,
    EpicsSignalRO,
    PVPositioner
)
from ophyd.status import Status
from ..utils._logging_setup import logger
logger.info(__file__)

KOHZU_SETTLE_TIME = 0.1


class MonoFeedback(Device):
    """ Mono feedback reading """

    readback = Component(EpicsSignalRO, 'mono_pid2.CVAL', kind='config')
    setpoint = Component(
        EpicsSignal,
        'mono_pid2.VAL',
        kind='config',
        put_complete=True
    )
    onoff = Component(
        EpicsSignal,
        'mono_pid2.FBON',
        kind='config',
        put_complete=True
    )


class KohzuPositioner(PVPositioner):

    readback = FormattedComponent(
        EpicsSignalRO,
        "{prefix}{_readback_pv}",
        kind="hinted",
        auto_monitor=True
    )
    setpoint = FormattedComponent(
        EpicsSignal,
        "{prefix}{_setpoint_pv}",
        kind="normal",
        put_complete=True
    )

    stop_theta = FormattedComponent(
        EpicsSignal, "{_theta_pv}.STOP", kind="omitted"
    )
    stop_y = FormattedComponent(EpicsSignal, "{_y_pv}.STOP", kind="omitted")

    actuate = Component(
        EpicsSignal, "KohzuPutBO", put_complete=True, kind="omitted"
    )
    actuate_value = 1

    done = Component(EpicsSignalRO, "KohzuMoving", kind="omitted")
    done_value = 0

    def __init__(
        self,
        prefix,
        *,
        limits=None,
        readback_pv="",
        setpoint_pv="",
        name=None,
        read_attrs=None,
        configuration_attrs=None,
        parent=None,
        egu="",
        **kwargs
    ):

        self._setpoint_pv = setpoint_pv
        self._readback_pv = readback_pv

        def get_motor_pv(label):
            _pv_signal = EpicsSignalRO(f"{prefix}Kohzu{label}PvSI", name="tmp")
            _pv_signal.wait_for_connection()
            return _pv_signal.get(as_string=True)

        self._theta_pv = get_motor_pv("Theta")
        self._y_pv = get_motor_pv("Y")

        super().__init__(
            prefix,
            limits=limits, name=name, read_attrs=read_attrs,
            configuration_attrs=configuration_attrs, parent=parent, egu=egu,
            **kwargs
        )

        # Make the default alias for the readback the name of the
        # positioner itself as in EpicsMotor.
        self.readback.name = self.name

    def _setup_move(self, position):
        '''Move and do not wait until motion is complete (asynchronous)'''
        self.log.debug('%s.setpoint = %s', self.name, position)

        # If wait is needed...
        # self.setpoint.set(position).wait(timeout=10)

        # If not...
        self.setpoint.put(position, wait=False)

        if self.actuate is not None:
            self.log.debug('%s.actuate = %s', self.name, self.actuate_value)
            self.actuate.put(self.actuate_value, wait=False)

    def stop(self, *, success=False):
        for motor in ["theta", "y"]:
            getattr(self, f"stop_{motor}").put(1, wait=False)
        super().stop(success=success)

    def move(self, position, wait=True, timeout=None, moved_cb=None):
        if abs(position - self.setpoint.get()) <= self.setpoint.tolerance:
            status = Status()
            status.set_finished()
        else:
            status = super().move(position, wait=wait, timeout=timeout,
                                  moved_cb=moved_cb)
        return status


class Monochromator(KohzuSeqCtl_Monochromator):
    """ Tweaks from apstools mono """

    wavelength = Component(
        KohzuPositioner,
        "",
        readback_pv="BraggLambdaRdbkAO",
        setpoint_pv="BraggLambdaAO",
        settle_time=KOHZU_SETTLE_TIME
    )

    energy = Component(
        KohzuPositioner,
        "",
        readback_pv="BraggERdbkAO",
        setpoint_pv="BraggEAO",
        settle_time=KOHZU_SETTLE_TIME
    )

    theta = Component(EpicsMotor, 'm1', labels=('motor',))

    theta_kohzu_screen = Component(
        KohzuPositioner,
        "",
        readback_pv="BraggThetaRdbkAO",
        setpoint_pv="BraggThetaAO",
        settle_time=KOHZU_SETTLE_TIME
    )

    y1 = None
    crystal_select = Component(EpicsMotor, 'm2', labels=('motor',))

    x2 = None
    y2 = Component(EpicsMotor, 'm3', labels=('motor',))
    z2 = Component(EpicsSignalRO, 'Zdummy', labels=('motor',))

    thf2 = Component(EpicsMotor, 'm4', labels=('motor',))
    chi2 = Component(EpicsMotor, 'm5', labels=('motor',))

    # TODO: these are offline for some reason
    # temp1 = FormattedComponent(
    #     LakeShore336Device, "4idaSoft:LS336:TC1:", labels=("lakeshore",)
    # )
    # temp2 = FormattedComponent(
    #     LakeShore336Device, "4idaSoft:LS336:TC2:", labels=("lakeshore",)
    # )

    # feedback = FormattedComponent(MonoFeedback, '4id:')

    def calibrate_energy(self, value):
        """Calibrate the mono energy.
        Parameters
        ----------
        value: float
            New energy for the current monochromator position.
        """
        self.use_set.put('Set', use_complete=True)
        self.energy.setpoint.put(value)
        self.use_set.put('Use', use_complete=True)


mono = Monochromator(
    '4idVDCM:', name='mono', labels=("monochromator", "energy")
)
