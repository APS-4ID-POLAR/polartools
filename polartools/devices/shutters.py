
"""
Shutters
"""
__all__ = ["ashutter", "bshutter"]

from apstools.devices import ApsPssShutterWithStatus
from time import sleep
from ..utils._logging_setup import logger
logger.info(__file__)


class PolarFEShutter(ApsPssShutterWithStatus):

    sleep_time = 5

    def _auto_shutter_subs(
        self, value, **kwargs
    ):
        if value == 0:
            while True:
                if self.pss_state.get() == 0:
                    self.open_signal.set(1)
                    sleep(self.sleep_time)
                else:
                    break

    def start_auto_shutter(self):
        self.pss_state.subscribe(
            self._auto_shutter_subs
        )

    def stop_auto_shutter(self):
        self.pss_state.unsubscribe_all()


ashutter = PolarFEShutter(
    "",
    "PA:04ID:A_BEAM_PRESENT",
    open_pv="PC:04ID:FES_OPEN_REQUEST",
    close_pv="PC:04ID:FES_CLOSE_REQUEST",
    name="ashutter"
)

bshutter = ApsPssShutterWithStatus(
    "",
    "PA:04ID:B_BEAM_PRESENT",
    open_pv="PC:04ID:SBS_OPEN_REQUEST",
    close_pv="PC:04ID:SBS_CLOSE_REQUEST",
    name="bshutter"
)
