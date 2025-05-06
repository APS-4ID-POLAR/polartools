"""
Polarization analyzer
"""

__all__ = ["pol"]

from ophyd import Device, Component, EpicsMotor
from .preamps import LocalPreAmp
from ..utils._logging_setup import logger
logger.info(__file__)


class PolAnalyzer(Device):
    y = Component(EpicsMotor, "m17", labels=("motor",))
    th = Component(EpicsMotor, "m9", labels=("motor",))

    vertical_preamp = Component(
        LocalPreAmp, 'A1', labels=('preamp', 'detector',), kind="config"
    )
    horizontal_preamp = Component(
        LocalPreAmp, 'A2', labels=('preamp', 'detector',), kind="config"
    )


pol = PolAnalyzer("4idbSoft:", name="pol", labels=("4idb",))

for pa in [pol.vertical_preamp, pol.horizontal_preamp]:
    pa.offset_fine._string = False
    for item in (
        "offset_fine set_all offset_value offset_unit offset_fine"
    ).split():
        getattr(pa, item).put_complete = True
