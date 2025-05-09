
"""
Shutters
"""

from apstools.devices import ApsPssShutterWithStatus
from time import sleep


class PolarShutter(ApsPssShutterWithStatus):

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
