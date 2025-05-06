
"""
Flags
"""

__all__ = [
    'flag_4ida_up',
    'flag_4ida_down',
    'flag_4idb'
]

from ophyd import EpicsMotor
from ..utils._logging_setup import logger
logger.info(__file__)

flag_4ida_up = EpicsMotor(
    "4idVDCM:m6", name="flag_4ida_up", labels=("4ida", "motor", "flag")
)

flag_4ida_down = EpicsMotor(
    "4idVDCM:m7", name="flag_4ida_down", labels=("4ida", "motor", "flag")
)

flag_4idb = EpicsMotor(
    "4idbSoft:m3", name="flag_4idb", labels=("4idb", "motor", "flag")
)