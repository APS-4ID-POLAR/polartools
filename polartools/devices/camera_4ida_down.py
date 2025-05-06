"""
Camera Upstream
"""

__all__ = ["flag_camera_4ida_down"]

from .ad_vimba import VimbaDetector


flag_camera_4ida_down = VimbaDetector(
    "4idaPostMonoBeam:",
    name="flag_camera_4ida_down",
    labels=("4ida", "camera", "detector", "flag")
)
