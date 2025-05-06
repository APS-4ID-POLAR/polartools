"""
WB Slit
"""

__all__ = [
    'wbslt'
]

from ophyd import Device, FormattedComponent, EpicsMotor
from ..utils import logger
logger.info(__file__)


class SlitDevice(Device):

    # Setting motors
    horizontal = FormattedComponent(
        EpicsMotor, '{prefix}{_motorsDict[hor]}', labels=('motor',)
    )

    diagonal = FormattedComponent(
        EpicsMotor, '{prefix}{_motorsDict[diag]}', labels=('motor',)
    )

    pitch = FormattedComponent(
        EpicsMotor, '{prefix}{_motorsDict[pitch]}', labels=('motor',)
    )

    yaw = FormattedComponent(
        EpicsMotor, '{prefix}{_motorsDict[yaw]}', labels=('motor',)
    )

    # Setting pseudo positioners
    vcen = FormattedComponent(
        EpicsMotor, '{prefix}{_slit_prefix}vCenter', labels=('motor',)
    )

    vsize = FormattedComponent(
        EpicsMotor, '{prefix}{_slit_prefix}vSize', labels=('motor',)
    )

    hcen = FormattedComponent(
        EpicsMotor, '{prefix}{_slit_prefix}hCenter', labels=('motor',)
    )

    hsize = FormattedComponent(
        EpicsMotor, '{prefix}{_slit_prefix}hSize', labels=('motor',)
    )

    def __init__(self, PV, name, motorsDict, slitnum, **kwargs):

        self._motorsDict = motorsDict
        self._slit_prefix = f'Slit{slitnum}:'

        super().__init__(prefix=PV, name=name, **kwargs)


# White beam slit
wbslt = SlitDevice(
    '4idVDCM:',
    'wbslt',
    {'hor': 'm9', 'diag': 'm10', 'pitch': 'm11', 'yaw': 'm12'},
    1,
    labels=("4ida", 'slit',)
)
