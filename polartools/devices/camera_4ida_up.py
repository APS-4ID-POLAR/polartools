"""
Camera Upstream
"""

__all__ = ["flag_camera_4ida_up"]

from .ad_vimba import VimbaDetector


flag_camera_4ida_up = VimbaDetector(
    "4idaPostMirrBeam:",
    name="flag_camera_4ida_up",
    labels=("4ida", "camera", "detector", "flag")
)
