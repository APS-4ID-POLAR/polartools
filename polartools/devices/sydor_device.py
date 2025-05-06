
__all__ = ["gsydor"]

from .quadems import SydorEMRO
from ..utils._logging_setup import logger
logger.info(__file__)

gsydor = SydorEMRO(
    "4idgSydor:T4U_BPM:", name="gsydor", labels=("detector",)
)
