"""
Caen power supply
"""

# WARNING: THIS IS A TEMPORARY SETUP WHILE WE DON'T HAVE THE EPICS SUPPORT #

__all__ = ["caenps"]

try:
    from caen_libs import caenhvwrapper as hv
except RuntimeError as excerror:
    print(
        "WARNING: could not find the CAEN library, please add the path to the "
        "LD_LIBRARY_PATH to the environment variables."
    )
    raise RuntimeError(excerror)

from ophyd import Device, Component, Signal, SignalRO
from ophyd.status import Status
from apstools.devices import PVPositionerSoftDoneWithStop
import threading
from time import sleep

DEVICE = hv.Device.open(
    hv.SystemType["SMARTHV"], hv.LinkType["TCPIP"], "10.54.115.56"
)


class StoppableThread(threading.Thread):
    def __init__(self, device, sleep_time=0.2):
        super().__init__()
        self._stop_event = threading.Event()
        self.device = device
        self.sleep_time = sleep_time

    def run(self):
        while not self._stop_event.is_set():
            self.device.cb_readback()
            sleep(self.sleep_time)

    def stop(self):
        self._stop_event.set()


class CaenSignal(Signal):
    def __init__(self, *args, channel=0, param_name="VSet", **kwargs):
        super().__init__(*args, **kwargs)
        self._channel = channel
        self._param = param_name

    def get(self, **kwargs):
        self._readback = DEVICE.get_ch_param(0, [self._channel], self._param)[0]
        return self._readback

    def put(self, value, **kwargs):
        if not isinstance(value, (int, float)):
            raise ValueError(
                f"file_path needs to be a number, but {type(value)} was "
                "entered."
            )
        super().put(value, **kwargs)
        DEVICE.set_ch_param(0, [self._channel], self._param, value)

    def set(self, value, **kwargs):
        self.put(value, **kwargs)
        # Do not check completion.
        st = Status()
        st.set_finished()
        return st


class CaenSignalRO(SignalRO):
    def __init__(self, *args, channel=0, param_name="VMon", **kwargs):
        super().__init__(*args, **kwargs)
        self._channel = channel
        self._param = param_name

    def get(self, **kwargs):
        self._readback = DEVICE.get_ch_param(0, [self._channel], self._param)[0]
        return self._readback


class CaenPositioner(PVPositionerSoftDoneWithStop):
    readback = Component(CaenSignalRO)
    setpoint = Component(CaenSignal)

    def __init__(self, *args, channel=0, **kwargs):
        super().__init__(
            *args, readback_pv="1", tolerance=0.5, use_target=True, **kwargs
        )
        self.readback._channel = channel
        self.setpoint._channel = channel
        self.timeout = 120
        self.settle_time = 5
        self.setpoint.subscribe(self.cb_update_target)
        self.thread = None

    def cb_update_target(self, value, *args, **kwargs):
        self.target.put(value)

    def set(self, new_position, *, timeout=None, moved_cb=None, wait=False):
        if self.thread is not None and self.thread.is_alive():
            self.thread.stop()
        self.thread = StoppableThread(self)
        self.thread.start()
        return super().set(
            new_position, timeout=timeout, moved_cb=moved_cb, wait=wait
        )

    def _done_moving(self, **kwargs):
        super()._done_moving(**kwargs)
        if self.thread is not None and self.thread.is_alive():
            self.thread.stop()
            self.thread = None


class CaenDevice(Device):
    ch1 = Component(CaenPositioner, "", channel=0)
    ch2 = Component(CaenPositioner, "", channel=1)
    ch3 = Component(CaenPositioner, "", channel=2)
    ch4 = Component(CaenPositioner, "", channel=3)
    ch5 = Component(CaenPositioner, "", channel=4)
    ch6 = Component(CaenPositioner, "", channel=5)
    ch7 = Component(CaenPositioner, "", channel=6)
    ch8 = Component(CaenPositioner, "", channel=7)

    device = DEVICE


caenps = CaenDevice("", name="caenps", labels=("detector",))
