
"""
Scalers
"""

__all__ = ['scaler_4tst']

from ophyd.scaler import ScalerCH
from ophyd.signal import Signal
from ophyd import Kind, Component
import time

from ..utils import logger
logger.info(__file__)


class PresetMonitorSignal(Signal):
    """ Signal that control the selected monitor channel """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._readback = 0

    def get(self, **kwargs):
        self._readback = self.parent._monitor.preset.get()
        if self.parent._monitor.s.name == 'Time':
            self._readback /= 1e7  # convert to seconds
        return self._readback

    def put(self, value, *, timestamp=None, force=False, metadata=None):
        """Put updates the internal readback value.

        The value is optionally checked first, depending on the value of force.
        In addition, VALUE subscriptions are run.
        Extra kwargs are ignored (for API compatibility with EpicsSignal kwargs
        pass through).
        Parameters
        ----------
        value : any
            Value to set
        timestamp : float, optional
            The timestamp associated with the value, defaults to time.time()
        metadata : dict, optional
            Further associated metadata with the value (such as alarm status,
            severity, etc.)
        force : bool, optional
            Check the value prior to setting it, defaults to False
        """
        self.log.debug(
            'put(value=%s, timestamp=%s, force=%s, metadata=%s)',
            value, timestamp, force, metadata
        )

        if float(value) <= 0:
            raise ValueError('preset_value has to be > 0.')

        # if self.parent._monitor.s.name == 'Time':
        if "chan01" in self.parent._monitor.name:
            value_put = 1e7*value  # convert to seconds
        else:
            value_put = value

        old_value = self._readback
        self.parent._monitor.preset.put(value_put)
        self._readback = value

        if metadata is None:
            metadata = {}

        if timestamp is None:
            timestamp = metadata.get('timestamp', time.time())

        metadata = metadata.copy()
        metadata['timestamp'] = timestamp
        self._metadata.update(**metadata)

        md_for_callback = {key: metadata[key]
                           for key in self._metadata_keys
                           if key in metadata}

        if 'timestamp' not in self._metadata_keys:
            md_for_callback['timestamp'] = timestamp

        self._run_subs(sub_type=self.SUB_VALUE, old_value=old_value,
                       value=value, **md_for_callback)


class LocalScalerCH(ScalerCH):

    preset_time = None
    preset_monitor = Component(PresetMonitorSignal, kind=Kind.config)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._monitor = self.channels.chan01  # Time is the default monitor.

    @property
    def channels_name_map(self):
        name_map = {}
        for channel in self.channels.component_names:
            # as defined in self.match_names()
            name = getattr(self.channels, channel).s.name
            if len(name) > 0:
                name_map[name] = channel
        return name_map

    def select_plot_channels(self, chan_names=None):

        self.match_names()
        name_map = self.channels_name_map

        if not chan_names:
            chan_names = name_map.keys()

        for ch in name_map.keys():
            try:
                channel = getattr(self.channels, name_map[ch])
            except KeyError:
                raise RuntimeError("The channel {} is not configured "
                                   "on the scaler.  The named channels are "
                                   "{}".format(ch, tuple(name_map)))
            if ch in chan_names:
                channel.s.kind = Kind.hinted
            else:
                if channel.kind.value != 0:
                    channel.s.kind = Kind.normal

    def select_read_channels(self, chan_names=None):
        """Select channels based on the EPICS name PV.

        Parameters
        ----------
        chan_names : Iterable[str] or None

            The names (as reported by the channel.chname signal)
            of the channels to select.
            If *None*, select all channels named in the EPICS scaler.
        """
        self.match_names()
        name_map = self.channels_name_map

        if chan_names is None:
            chan_names = name_map.keys()

        read_attrs = ['chan01']  # always include time
        for ch in chan_names:
            try:
                read_attrs.append(name_map[ch])
            except KeyError:
                raise RuntimeError("The channel {} is not configured "
                                   "on the scaler.  The named channels are "
                                   "{}".format(ch, tuple(name_map)))

        self.channels.kind = Kind.normal
        self.channels.read_attrs = list(read_attrs)
        self.channels.configuration_attrs = list(read_attrs)
        if len(self.hints['fields']) == 0:
            self.select_plot_channels(chan_names)

    @property
    def monitor(self):
        return self._monitor.s.name

    @monitor.setter
    def monitor(self, value):
        """
        Selects the monitor channel.

        Parameters
        ----------
        value : str
            Can be either the name of the component channel (like 'chan01'),
            or the name of that channel (like 'Ion Ch 1').
        """

        # Check that value is a valid name.
        name_map = self.channels_name_map
        if value not in (set(name_map.keys()) | set(name_map.values())):
            raise ValueError(f"Monitor must be either a channel name or the "
                             "channel component. Valid entries are one of "
                             f"these: {name_map.keys()}, or these: "
                             f"{name_map.keys()}.")

        # Changes value to the channel number if needed. From here on,
        # value will always be something like 'chan01'.
        if value in name_map.keys():
            value = name_map[value]

        # Checks/modifies the channel Kind.
        channel = getattr(self.channels, value)
        if channel.kind == Kind.omitted:
            channel.kind = Kind.normal

        # Adjust gates
        for channel_name in self.channels.component_names:
            chan = getattr(self.channels, channel_name)
            target = 'Y' if chan == channel else 'N'
            chan.gate.put(target, use_complete=True)

        self._monitor = channel


scaler_4tst = LocalScalerCH(
    '4tst:scaler1', name='scaler_4tst', labels=('detector', 'scaler')
)
scaler_4tst.monitor = 'chan01'
scaler_4tst.select_read_channels()
scaler_4tst.select_plot_channels()
