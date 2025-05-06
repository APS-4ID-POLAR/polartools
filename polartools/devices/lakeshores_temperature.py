"""
Lakeshore Temperature Controllers
"""

from apstools.devices import LakeShore336Device
from ..utils._logging_setup import logger
logger.info(__file__)

temperature_4idg = LakeShore336Device(
    "4idgSoft:LS336:cryo:", name="temperature_4idg"
)
