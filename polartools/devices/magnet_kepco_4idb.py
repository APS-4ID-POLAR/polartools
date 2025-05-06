"""
Kepko power supply
"""

# __all__ = ['kepco']

from ophyd import Component, FormattedComponent, Device, Kind
from ophyd import EpicsSignal, EpicsSignalRO
from apstools.devices import PVPositionerSoftDoneWithStop

from ..utils._logging_setup import logger
logger.info(__file__)


class LocalPositioner(PVPositionerSoftDoneWithStop):
    """ Voltage/Current positioner """

    readback = FormattedComponent(
        EpicsSignalRO, '{prefix}d{_type}', kind='hinted',
    )

    setpoint = FormattedComponent(
        EpicsSignal, "{prefix}{_type}", write_pv="{prefix}set{_type}",
    )

    def __init__(self, *args, progtype, **kwargs):
        self._type = progtype
        super().__init__(*args, readback_pv="1", **kwargs)


class KepcoController(Device):

    voltage = Component(LocalPositioner, '', progtype='V', tolerance=0.02)
    current = Component(LocalPositioner, '', progtype='C', tolerance=0.03)

    mode = Component(
        EpicsSignal, 'setMode', kind='config', string=True, auto_monitor=True
    )

    remote = Component(
        EpicsSignal, 'setRemote', kind='config', string=True, auto_monitor=True
    )

    enable = Component(EpicsSignal, 'Enable.VAL', kind='omitted', string=True)

    @mode.sub_value
    def mode_change(self, value=None, **kwargs):
        if value == 'Current':
            self.current.readback.kind = Kind.hinted
            self.voltage.readback.kind = Kind.normal

        if value == 'Voltage':
            self.current.readback.kind = Kind.normal
            self.voltage.readback.kind = Kind.hinted


# kepco = KepcoController('4idd:BOP:PS1:', name='kepco', labels=("magnet",))
# kepco.mode_change(value=kepco.mode.get())
