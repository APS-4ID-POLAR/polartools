"""
TetrAMMs
"""

__all__ = ["btetramm"]

from .quadems import TetrAMMRO
from ..utils._logging_setup import logger
logger.info(__file__)


btetramm = TetrAMMRO(
    "4idbSoft:TetrAMM:", name="btetramm", labels=("detector",)
)
