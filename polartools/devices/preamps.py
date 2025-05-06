"""
SRS 570 pre-amplifiers
"""

# __all__ = ["preamp1", "preamp2"]

from apstools.devices import SRS570_PreAmplifier
from pint import Quantity
from pandas import DataFrame
from bluesky.plan_stubs import mv, rd, trigger, checkpoint, sleep
from numpy import array, where, round, linspace, polyfit, poly1d
from ..utils._logging_setup import logger
logger.info(__file__)


class LocalPreAmp(SRS570_PreAmplifier):

    def __init__(self, *args, scaler_channel=None, shutter=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._scaler_channel = scaler_channel
        self._shutter = shutter
        self._shutter_vals = dict(on=1, off=0)

    @property
    def shutter(self):
        return self._shutter

    @shutter.setter
    def shutter(self, device):
        if not hasattr(device, "set"):
            raise ValueError("device must have a 'set' attribute.")
        self._shutter = device

    @property
    def _offset_current_table(self):
        convert = dict(vals=[], units=[], mags=[])
        for val in self.offset_value.enum_strs:
            for unit in self.offset_unit.enum_strs:
                convert["units"].append(unit)
                convert["vals"].append(val)
                convert["mags"].append(round(
                    Quantity(float(val), unit).to("A").magnitude, decimals=12
                ))
        return DataFrame(convert).set_index("mags").sort_index()

    @property
    def _sensitivity_table(self):
        convert = dict(vals=[], units=[], mags=[])
        for val in self.sensitivity_value.enum_strs:
            for unit in self.sensitivity_unit.enum_strs:
                convert["units"].append(unit)
                convert["vals"].append(val)
                convert["mags"].append(round(
                    Quantity(float(val), unit).to("A/V").magnitude, decimals=12
                ))
        return DataFrame(convert).set_index("mags").sort_index()

    def opt_sens_plan(self, scaler_channel=None, time=0.1, delay=1):

        table = self._sensitivity_table
        if scaler_channel is None:
            scaler_channel = self._scaler_channel
        self._scaler_channel.root.monitor = "Time"
        yield from mv(scaler_channel.root.preset_monitor, time)

        # Count the scaler
        yield from mv(self.set_all, 1)
        yield from sleep(delay)
        yield from trigger(scaler_channel.root, wait=True)
        value = yield from rd(scaler_channel.s)
        gain = round(self.computed_gain, 12)

        # print(value, gain)

        best = [None, None, 1e10, False]
        if (value/time > 1000) and (value/time < 950000):
            optimal_gain = value/time/6e5*gain
            opt_list = array(table.index)/optimal_gain
            i = where(opt_list == opt_list[opt_list > 1].min())[0][0]
            best = [table.iloc[i]["vals"], table.iloc[i]["units"], value, True]
        else:
            direction = 1 if value/time > 4e5 else -1
            start = table.index.get_loc(gain)
            end = 0 if direction == -1 else table.shape[0]

            for i in range(start, end, direction):
                yield from mv(
                    self.sensitivity_value, table.iloc[i]["vals"],
                    self.sensitivity_unit, table.iloc[i]["units"]
                )
                yield from mv(self.set_all, 1)
                yield from sleep(delay)
                yield from trigger(scaler_channel.root, wait=True)
                value = yield from rd(scaler_channel.s)

                # print(value, best)
                if (
                    (abs(value/time - 5e5) < abs(best[2]/time - 5e5))
                    & (value/time < 6e5)
                ):
                    best = [
                        table.iloc[i]["vals"],
                        table.iloc[i]["units"],
                        value,
                        True
                    ]
                else:
                    if best[-1] is True:
                        break

        if best[0] is not None:
            yield from mv(
                    self.sensitivity_value, best[0],
                    self.sensitivity_unit, best[1]
                )
            yield from mv(self.set_all, 1)

    def opt_offset_plan(self, scaler_channel=None, time=0.1, delay=1):

        gain_pv_conversion = self._offset_current_table
        if scaler_channel is None:
            scaler_channel = self._scaler_channel
        scaler_channel.root.monitor = "Time"
        yield from mv(scaler_channel.root.preset_monitor, time)

        sign_str = yield from rd(self.offset_sign)
        current_sign = +1 if sign_str == "+" else -1

        def _offset_scan(
            rng,
            best=dict(vals=None, units=None, count=1e10, done=False, fine=500)
        ):
            fine = yield from rd(self.offset_fine)
            for i in rng:
                yield from mv(
                    self.offset_value, gain_pv_conversion.iloc[i]["vals"],
                    self.offset_unit, gain_pv_conversion.iloc[i]["units"]
                )
                yield from mv(self.set_all, 1)
                yield from sleep(delay)
                yield from trigger(scaler_channel.root, wait=True)
                value = yield from rd(scaler_channel.s)
                value /= time

                # print(value, best)
                # If value is better than previous one, then update.
                if (
                    (abs(value - 200) < abs(best["count"] - 200)) &
                    (value*time > 2)
                ):
                    best["vals"] = gain_pv_conversion.iloc[i]["vals"]
                    best["units"] = gain_pv_conversion.iloc[i]["units"]
                    best["count"] = value
                    best["fine"] = fine
                # If this is a good value.
                elif (best["count"] > 50) & (best["count"] < 400):
                    best["done"] = True
                    break
            return best

        # Keep the offset sign, start with the same "number" as the sensitivity
        yield from mv(self.offset_fine, current_sign*500)
        yield from sleep(delay)
        start = gain_pv_conversion.index.get_loc(
            round(self.computed_gain, 12)
        )
        best = yield from _offset_scan(range(start, -1, -1))
        # print(best)
        if not best["done"]:
            # Change sign
            yield from mv(self.offset_fine, -1*current_sign*500)
            yield from sleep(delay)
            best = yield from _offset_scan(range(0, start, 1), best=best)

        yield from mv(
                self.offset_fine, best["fine"],
                self.offset_value, best["vals"],
                self.offset_unit, best["units"]
            )
        yield from mv(self.set_all, 1)

    def opt_fine_plan(
            self,
            scaler_channel=None,
            start=None,
            end=None,
            steps=11,
            time=0.1,
            delay=1
    ):

        if scaler_channel is None:
            scaler_channel = self._scaler_channel
        scaler_channel.root.monitor = "Time"
        yield from mv(scaler_channel.root.preset_monitor, time)

        sign = yield from rd(self.offset_sign)
        factor = +1 if sign == "+" else -1

        if start is None:
            start = factor*1

        if end is None:
            end = factor*1000

        yield from mv(self.offset_fine, start)
        yield from mv(self.set_all, 1)
        yield from sleep(delay)

        fines = linspace(start, end, steps)
        counts = []
        for fine in fines:
            yield from mv(self.offset_fine, fine)
            yield from mv(self.set_all, 1)
            yield from sleep(delay)
            yield from trigger(scaler_channel.root, wait=True)
            value = yield from rd(scaler_channel.s)
            counts.append(value/time)

        pos = int(poly1d(polyfit(counts, fines, 1))(200))
        pos = pos if pos < 1000 else 1000
        pos = pos if pos > -1000 else -1000

        yield from mv(self.offset_fine, pos)
        yield from mv(self.set_all, 1)

    def optimize_plan(self, scaler_channel=None, time=0.1, delay=1):

        if self._shutter is None:
            raise ValueError("Shutter is not configured.")

        if scaler_channel is None:
            scaler_channel = self._scaler_channel
        scaler_channel.root.monitor = "Time"
        yield from mv(scaler_channel.root.preset_monitor, time)

        yield from mv(self.set_all, 1)
        yield from sleep(delay)
        # normally this would be the beam shutter, and not a argument.
        yield from mv(self.shutter, self._shutter_vals["off"])

        # If initial dark current is too high, then get it roughly right
        yield from trigger(scaler_channel.root, wait=True)
        zero = yield from rd(scaler_channel.s)

        if zero/time > 5e4:
            yield from self.opt_offset_plan(scaler_channel, time, delay)

        yield from checkpoint()

        yield from mv(self.shutter, self._shutter_vals["on"])
        yield from self.opt_sens_plan(scaler_channel, time, delay)

        yield from checkpoint()

        yield from mv(self.shutter, self._shutter_vals["off"])
        for func in [self.opt_offset_plan, self.opt_fine_plan]:
            yield from func(scaler_channel, time=time, delay=delay)
            yield from checkpoint()

        yield from mv(self.shutter, self._shutter_vals["on"])
        yield from sleep(delay)
        yield from trigger(scaler_channel.root, wait=True)
        value = yield from rd(scaler_channel.s)

        if (value/time > 6e5) or (value/time < 3e5):

            yield from mv(self.shutter, self._shutter_vals["on"])
            yield from self.opt_sens_plan(scaler_channel, time, delay)
            yield from checkpoint()

            yield from mv(self.shutter, self._shutter_vals["off"])
            for func in [self.opt_offset_plan, self.opt_fine_plan]:
                yield from func(scaler_channel, time=time, delay=delay)
                yield from checkpoint()

        yield from mv(self.shutter, self._shutter_vals["on"])
        yield from mv(self.set_all, 1)


# preamp1 = LocalPreAmp(
#     '4tst:A1', name="preamp1", labels=('preamp', 'detector',)
# )
# preamp2 = LocalPreAmp(
#     '4tst:A2', name="preamp2", labels=('preamp', 'detector',)
# )

# preamp1.offset_fine._string = False
# preamp2.offset_fine._string = False

# for pa in [preamp1, preamp2]:
#     for item in (
#         "offset_fine set_all offset_value offset_unit offset_fine"
#     ).split():
#         getattr(pa, item).put_complete = True
