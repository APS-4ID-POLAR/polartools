'''
SoftGlueZynq
'''

from ophyd import (
    Component, Device, EpicsSignal, EpicsSignalRO, DynamicDeviceComponent
)
from collections import OrderedDict
from bluesky.plan_stubs import mv
from logging import getLogger

logger = getLogger(__name__)


def _io_fields(num=16):
    defn = OrderedDict()
    for i in range(1, num + 1):
        defn[f"fi{i}"] = (SoftGlueSignal, f"SG:FI{i}", {"kind": "config"})
        defn[f"fo{i}"] = (SoftGlueSignal, f"SG:FO{i}", {"kind": "config"})
    return defn


def _buffer_fields(num=4):
    defn = OrderedDict()
    for i in range(1, num + 1):
        defn[f"in{i}"] = (
            SoftGlueSignal, f"SG:BUFFER-{i}_IN", {"kind": "config"}
        )
        defn[f"out{i}"] = (
            SoftGlueSignal, f"SG:BUFFER-{i}_OUT", {"kind": "config"}
        )
    return defn


def _dma_fields(num=8, first_letter="I"):
    defn = OrderedDict()
    defn["enable"] = (EpicsSignal, "1acquireDmaEnable", {"kind": "config"})
    defn["scan"] = (EpicsSignal, "1acquireDma.SCAN", {"kind": "config"})
    defn["read_button"] = (EpicsSignal, "1acquireDma.PROC", {"kind": "omitted"})
    defn["clear_button"] = (EpicsSignal, "1acquireDma.D", {"kind": "omitted"})
    defn["clear_buffer"] = (EpicsSignal, "1acquireDma.F", {"kind": "omitted"})
    defn["words_in_buffer"] = (
        EpicsSignalRO, "1acquireDma.VALJ", {"kind": "config"}
    )
    defn["events"] = (EpicsSignalRO, "1acquireDma.VALI", {"kind": "config"})
    for i in range(1, num + 1):
        defn[f"channel_{i}_name"] = (
            EpicsSignal, f"1s{i}name", {"kind": "config"}
        )
        defn[f"channel_{i}_scale"] = (
            EpicsSignal,
            f"1acquireDma.{chr(ord(first_letter)+i-1)}",
            {"kind": "config"}
        )
    return defn


class SoftGlueSignal(Device):
    signal = Component(EpicsSignal, "_Signal", kind="config")
    bi = Component(EpicsSignal, "_BI", kind="config")


class SoftGlueZynqDevideByN(Device):
    enable = Component(SoftGlueSignal, "ENABLE", kind="config")
    clock = Component(SoftGlueSignal, "CLOCK", kind="config")
    reset = Component(SoftGlueSignal, "RESET", kind="config")
    out = Component(SoftGlueSignal, "OUT", kind="config")
    n = Component(EpicsSignal, "N", kind="config")


class SoftGlueZynqUpCounter(Device):
    enable = Component(SoftGlueSignal, "ENABLE", kind="config")
    clock = Component(SoftGlueSignal, "CLOCK", kind="config")
    reset = Component(SoftGlueSignal, "CLEAR", kind="config")
    counts = Component(EpicsSignalRO, "COUNTS", kind="config")


class SoftGlueZynqGateDly(Device):
    input = Component(SoftGlueSignal, "IN", kind="config")
    clock = Component(SoftGlueSignal, "CLK", kind="config")
    delay = Component(EpicsSignal, "DLY", kind="config")
    width = Component(EpicsSignal, "WIDTH", kind="config")
    out = Component(SoftGlueSignal, "OUT", kind="config")


class SoftGlueScalToStream(Device):
    reset = Component(SoftGlueSignal, "RESET", kind="config")
    chadv = Component(SoftGlueSignal, "CHADV", kind="config")
    imtrig = Component(SoftGlueSignal, "IMTRIG", kind="config")
    flush = Component(SoftGlueSignal, "FLUSH", kind="config")
    full = Component(SoftGlueSignal, "FULL", kind="config")
    advdone = Component(SoftGlueSignal, "ADVDONE", kind="config")
    imdone = Component(SoftGlueSignal, "IMDONE", kind="config")
    fifo = Component(EpicsSignalRO, "FIFO", kind="config")
    dmawords = Component(EpicsSignal, "DMAWORDS", kind="config")


class SoftGlueClocks(Device):
    clock_10MHz = Component(SoftGlueSignal, "10MHZ_CLOCK", kind="config")
    clock_20MHz = Component(SoftGlueSignal, "20MHZ_CLOCK", kind="config")
    clock_50MHz = Component(SoftGlueSignal, "50MHZ_CLOCK", kind="config")
    clock_variable = Component(SoftGlueSignal, "VAR_CLOCK", kind="config")


class SoftGlueZynqDevice(Device):
    dma = DynamicDeviceComponent(_dma_fields())

    # Buffer 1 --> general enable
    # Buffer 2 --> detector enable
    # Buffer 3 --> reset counters
    # Buffer 4 --> reset DMA/FIFO
    buffers = DynamicDeviceComponent(_buffer_fields())

    # I/O fields
    io = DynamicDeviceComponent(_io_fields())

    # Using channel #4 to count when the gate is off.
    up_counter_count = Component(
        SoftGlueZynqUpCounter, "SG:UpCntr-1_", kind="config"
    )
    up_counter_trigger = Component(
        SoftGlueZynqUpCounter, "SG:UpCntr-2_", kind="config"
    )
    up_counter_gate_on = Component(
        SoftGlueZynqUpCounter, "SG:UpCntr-3_", kind="config"
    )
    up_counter_gate_off = Component(
        SoftGlueZynqUpCounter, "SG:UpCntr-4_", kind="config"
    )

    # Setup the frequency of the count and trigger based on 10 MHz clock.
    div_by_n_count = Component(
        SoftGlueZynqDevideByN, "SG:DivByN-1_", kind="config"
    )
    div_by_n_trigger = Component(
        SoftGlueZynqDevideByN, "SG:DivByN-2_", kind="config"
    )
    div_by_n_interrupt = Component(
        SoftGlueZynqDevideByN, "SG:DivByN-3_", kind="config"
    )

    # Create a gate pulse
    gate_trigger = Component(
        SoftGlueZynqGateDly, "SG:GateDly-1_", kind="config"
    )

    # Send data to DMA
    scaltostream = Component(SoftGlueScalToStream, "SG:scalToStream-1_")

    # Clocks
    clocks = Component(SoftGlueClocks, "SG:", kind="config")

    def __init__(
        self, *args, reset_sleep_time=0.2, reference_clock=1e7, **kwargs
    ):
        super().__init__(*args, **kwargs)
        self._reset_sleep_time = reset_sleep_time
        self._reference_clock = reference_clock

    def start_softglue(self):
        yield from mv(self.buffers.in1.signal, "1")

    def stop_softglue(self):
        yield from mv(self.buffers.in1.signal, "0")

    def start_detectors(self):
        yield from mv(self.buffers.in2.signal, "1")

    def stop_detectors(self):
        yield from mv(self.buffers.in2.signal, "0")

    def reset_plan(self):
        yield from mv(
            self.buffers.in3.signal, "1!", self.buffers.in4.signal, "1!"
        )

    def clear_enable_dma(self):
        yield from mv(self.dma.clear_button, 1, self.dma.clear_buffer, 1)
        yield from mv(self.dma.enable, 1)

    def clear_disable_dma(self):
        yield from mv(self.dma.clear_button, 1, self.dma.clear_buffer, 1)
        yield from mv(self.dma.enable, 0)

    def setup_trigger_plan(
        self, period_time, pulse_width_time, pulse_delay_time=0
    ):
        yield from mv(
            self.div_by_n_trigger.n, self._reference_clock * period_time,
            self.gate_trigger.delay, self._reference_clock * pulse_delay_time,
            self.gate_trigger.width, self._reference_clock * pulse_width_time
        )

    def setup_count_plan(self, time):
        yield from mv(self.div_by_n_count.n, self._reference_clock * time)

    def default_settings(self, timeout=10):

        logger.info("Setting up clocks.")

        self.clocks.clock_10MHz.signal.set("ck10").wait(timeout)

        self.div_by_n_count.enable.signal.set("enable").wait(timeout)
        self.div_by_n_count.clock.signal.set("ck10").wait(timeout)
        self.div_by_n_count.reset.signal.set("1").wait(timeout)
        self.div_by_n_count.reset.signal.set("0").wait(timeout)
        self.div_by_n_count.out.signal.set("ckUser").wait(timeout)
        self.div_by_n_count.n.set(10000).wait(timeout)

        self.div_by_n_trigger.enable.signal.set("enableDet").wait(timeout)
        self.div_by_n_trigger.clock.signal.set("ck10").wait(timeout)
        self.div_by_n_trigger.reset.signal.set("1").wait(timeout)
        self.div_by_n_trigger.reset.signal.set("0").wait(timeout)
        self.div_by_n_trigger.out.signal.set("ckDet").wait(timeout)
        self.div_by_n_trigger.n.set(1000000).wait(timeout)

        self.div_by_n_interrupt.enable.signal.set("enable").wait(timeout)
        self.div_by_n_interrupt.clock.signal.set("ck10").wait(timeout)
        self.div_by_n_interrupt.reset.signal.set("1").wait(timeout)
        self.div_by_n_interrupt.reset.signal.set("0").wait(timeout)
        self.div_by_n_interrupt.out.signal.set("ckInt").wait(timeout)
        self.div_by_n_interrupt.n.set(100000).wait(timeout)

        logger.info("Setting up buffers.")

        for i in range(1, 5):
            getattr(self.buffers, f"in{i}").signal.set("1!").wait(timeout)

        self.buffers.out1.signal.set("enable").wait(timeout)
        self.buffers.out2.signal.set("enableDet").wait(timeout)
        self.buffers.out3.signal.set("resetCnters").wait(timeout)
        self.buffers.out4.signal.set("reset").wait(timeout)

        logger.info("Setting up counters.")

        self.up_counter_count.enable.signal.set("enable").wait(timeout)
        self.up_counter_count.clock.signal.set("ckUser").wait(timeout)
        self.up_counter_count.reset.signal.set("resetCnters").wait(timeout)

        self.up_counter_trigger.enable.signal.set("enable").wait(timeout)
        self.up_counter_trigger.clock.signal.set("ckDet").wait(timeout)
        self.up_counter_trigger.reset.signal.set("resetCnters").wait(timeout)

        self.up_counter_gate_on.enable.signal.set("enable").wait(timeout)
        self.up_counter_gate_on.clock.signal.set("gateTrigger").wait(timeout)
        self.up_counter_gate_on.reset.signal.set("resetCnters").wait(timeout)

        self.up_counter_gate_on.enable.signal.set("enable").wait(timeout)
        self.up_counter_gate_on.clock.signal.set("gateTrigger*").wait(timeout)
        self.up_counter_gate_on.reset.signal.set("resetCnters").wait(timeout)

        logger.info("Setting up outputs.")

        self.io.fo1.signal.set("gateTrigger").wait(timeout)
        self.io.fo15.signal.set("ckInt").wait(timeout)
        self.io.fo16.signal.set("reset").wait(timeout)

        logger.info("Setting up DMA transfer.")

        self.scaltostream.reset.signal.set("reset*").wait(timeout)
        self.scaltostream.chadv.signal.set("ckUser").wait(timeout)
        self.scaltostream.dmawords.set(16384).wait(timeout)
        self.dma.scan.set(6).wait(timeout)

        logger.info("Setting up trigger transfer.")

        self.gate_trigger.input.signal.set("ckDet").wait(timeout)
        self.gate_trigger.clock.signal.set("ck10").wait(timeout)
        self.gate_trigger.out.signal.set("gateTrigger").wait(timeout)
        self.gate_trigger.delay.set(0).wait(timeout)
        self.gate_trigger.width.set(500000).wait(timeout)
