"""
X-ray Eye
"""

__all__ = ["xrayeye"]

from instrument.devices.ad_vimba import VimbaDetector


xrayeye = VimbaDetector(
    "4idXrayEye:",
    name="xrayeye",
    labels=("4idg", "camera", "detector")
)
