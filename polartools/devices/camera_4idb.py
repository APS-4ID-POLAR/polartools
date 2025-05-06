"""
Camera Upstream
"""

__all__ = ["flag_camera_4idb"]

from .ad_vimba import VimbaDetector


flag_camera_4idb = VimbaDetector(
    "4idbPostToroBeam:",
    name="flag_camera_4idb",
    labels=("4idb", "camera", "detector", "flag")
)
