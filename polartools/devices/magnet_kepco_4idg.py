"""
Diffractometer magnet
"""

from apstools.devices import PVPositionerSoftDone
from ophyd import Component, EpicsSignal, EpicsSignalRO, Device
from ..utils._logging_setup import logger
logger.info(__file__)


class KepcoDevice(Device):

    # Info and status
    manufacturer = Component(EpicsSignalRO, "manufacturer", kind="omitted")
    model = Component(EpicsSignalRO, "model", kind="config")
    model2 = Component(EpicsSignalRO, "model2", kind="config")
    calibration = Component(EpicsSignalRO, "calibration", kind="config")
    serial = Component(EpicsSignalRO, "serial", kind="omitted")
    firmware = Component(EpicsSignalRO, "firmware", kind="config")

    status = Component(EpicsSignalRO, "meas", kind="config")
    status_output = Component(EpicsSignalRO, "D0", string=True, kind="omitted")
    status_list = Component(EpicsSignalRO, "D1", string=True, kind="omitted")
    status_errors = Component(EpicsSignalRO, "D2", string=True, kind="omitted")
    status_mode = Component(EpicsSignalRO, "D3", string=True, kind="omitted")
    status_protect = Component(EpicsSignalRO, "D4", string=True, kind="omitted")
    status_message = Component(EpicsSignalRO, "D5", string=True, kind="omitted")
    last_error_message = Component(
        EpicsSignalRO, "systerr", string=True, kind="omitted"
    )

    # Field control
    field = Component(
        PVPositionerSoftDone,
        "",
        readback_pv="interpAtoT_y1.VAL",
        setpoint_pv="interpTtoA_x.VAL",
        tolerance=0.01
    )

    voltage = Component(EpicsSignalRO, "measvolt", kind="config")

    current = Component(
        PVPositionerSoftDone,
        "",
        readback_pv="meascurr",
        setpoint_pv="setcurr",
        tolerance=0.1
    )

    mode = Component(
        EpicsSignal, "funcmode", write_pv="setfuncmode", kind="config"
    )
    _auto_mode_subs = []

    enable = Component(EpicsSignal, "outp", write_pv="setoutp", kind="config")

    # Configs
    current_protection_positive = Component(
        EpicsSignalRO, "currprotpos", kind="omitted"
    )
    current_protection_negative = Component(
        EpicsSignalRO, "currprotneg", kind="omitted"
    )

    voltage_protection_positive = Component(
        EpicsSignal, "voltprotpos", write_pv="setvoltprotpos", kind="config"
    )
    voltage_protection_negative = Component(
        EpicsSignal, "voltprotneg", write_pv="setvoltprotneg", kind="config"
    )

    current_limit_positive = Component(
        EpicsSignal, "currlimpos", write_pv="setcurrlimpos", kind="config"
    )
    current_limit_negative = Component(
        EpicsSignal, "currlimneg", write_pv="setcurrlimneg", kind="config"
    )

    voltage_limit_positive = Component(
        EpicsSignalRO, "voltlimpos", kind="omitted"
    )
    voltage_limit_negative = Component(
        EpicsSignalRO, "voltlimneg", kind="omitted"
    )

    # Temperatures
    t1 = Component(EpicsSignalRO, "T1", kind="config")
    t2 = Component(EpicsSignalRO, "T2", kind="config")
    t3 = Component(EpicsSignalRO, "T3", kind="config")
    status_temperature = Component(EpicsSignalRO, "hts_status", kind="config")

    def default_settings(self):
        self.start_auto_mode()

    def _auto_mode(self, value, **kwargs):
        if value in ("VOLT", 0):
            logger.warning(
                "Cannot use the voltage mode. If you want to change this "
                "behavior run the .stop_auto_mode() command."
            )
            self.mode.set("CURRENT").wait()

    def start_auto_mode(self):
        self._auto_mode_subs.append(
            self.mode.subscribe(self._auto_mode, run=True)
        )

    def stop_auto_mode(self):
        for _sub in self._auto_mode_subs:
            self.unsubscribe(_sub)


magnet2t = KepcoDevice("4idkepco:", name="magnet2t", labels=("magnet", "4idg"))
