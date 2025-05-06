"""
USB-CTR8
"""

__all__ = ["ctr8"]

from apstools.devices import MeasCompCtr
from ..utils import logger
logger.info(__file__)


class CustomMeasCompCtr(MeasCompCtr):

    # TODO: We don't have everything setup yet

    long_in = None

    binary_in_1 = None
    binary_in_2 = None
    binary_in_3 = None
    binary_in_4 = None
    binary_in_5 = None
    binary_in_6 = None
    binary_in_7 = None
    binary_in_8 = None

    long_out = None

    binary_out_1 = None
    binary_out_2 = None
    binary_out_3 = None
    binary_out_4 = None
    binary_out_5 = None
    binary_out_6 = None
    binary_out_7 = None
    binary_out_8 = None

    binary_direction_1 = None
    binary_direction_2 = None
    binary_direction_3 = None
    binary_direction_4 = None
    binary_direction_5 = None
    binary_direction_6 = None
    binary_direction_7 = None
    binary_direction_8 = None


ctr8 = CustomMeasCompCtr(
    "4idCTR8_1:", name="ctr8", labels=("detector", "scaler")
)
