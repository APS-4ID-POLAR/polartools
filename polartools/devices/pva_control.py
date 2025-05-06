
"""
Device to control the PositionerStream
"""

__all__ = ["positioner_stream"]

from pvapy import Channel
from ophyd import Device, Signal, Component
from ophyd.status import Status
from pathlib import Path
from ..utils.config import iconfig
from ..utils import logger
logger.info(__file__)

HDF1_NAME_TEMPLATE = iconfig["AREA_DETECTOR"]["HDF5_FILE_TEMPLATE"]
HDF1_FILE_EXTENSION = iconfig["AREA_DETECTOR"]["HDF5_FILE_EXTENSION"]
HDF1_NAME = Path(HDF1_NAME_TEMPLATE + "." + HDF1_FILE_EXTENSION)


class PVASignal(Signal):
    def __init__(self, *args, pva_channel="", pva_label="", **kwargs):
        super().__init__(*args, **kwargs)
        self._pva = Channel(pva_channel)
        self._pva_label = pva_label

    def get(self, **kwargs):
        return self._pva.get().toDict()[self._pva_label]

    def put(self, value, **kwargs):
        if not isinstance(value, str):
            raise ValueError(
                f"file_path needs to be a string, but {type(value)} was "
                "entered."
            )
        self._pva.putString(value, self._pva_label)

    def set(self, value, **kwargs):
        self.set(value, **kwargs)
        # Do not check completion.
        st = Status()
        st.set_finished()
        return st


class PositionerStream(Device):
    file_pva = Channel("4idSoftGluePVA:outputFile")
    status_pva = Channel("4idSoftGluePVA:status")
    start_pva = Channel("4idSoftGluePVA:start")
    stop_pva = Channel("4idSoftGluePVA:stop")

    # These will be signals that Bluesky can read and save in the catalog.
    file_path = Component(
        PVASignal,
        pva_channel="4idSoftGluePVA:outputFile",
        pva_label="filePath",
        kind="normal"
    )

    file_name = Component(
        PVASignal,
        pva_channel="4idSoftGluePVA:outputFile",
        pva_label="fileName",
        kind="normal"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    _status_obj = None

    @property
    def status(self):
        return self.status_pva.get().toDict()["value"]

    def start_signal(self):
        self.start_pva.putInt(1)

    def stop_signal(self):
        self.stop_pva.putInt(1)

    def start_stream(self):
        def _status_sub(inp):
            if inp["value"] == "Acquiring":
                self._status_obj.set_finished()
                self.status_pva.stopMonitor()

        self._status_obj = Status()

        if self.status != "Acquiring":
            self.start_pva.stopMonitor()
            self.start_signal()
            self.status_pva.monitor(
                _status_sub, "field(value, alarm, timeStamp)"
            )
        else:
            self._status_obj.set_finished()

        return self._status_obj

    def stop_stream(self):
        def _status_sub(inp):
            if inp["value"] == "Idle":
                self._status_obj.set_finished()
                self.status_pva.stopMonitor()

        self._status_obj = Status()

        if self.status != 'Idle':
            self.start_pva.stopMonitor()
            self.stop_signal()
            self.status_pva.monitor(
                _status_sub, "field(value, alarm, timeStamp)"
            )
        else:
            self._status_obj.set_finished()

        return self._status_obj

    def set(self, value, **kwargs):

        if value not in [1, 0]:
            raise ValueError("Value must be 1 or 0.")

        return self.start_stream() if value == 1 else self.stop_stream()

    def stop(self, **kwargs):
        super().stop(**kwargs)
        self.stop_signal()

    def setup_file_path_name(self, path, name_base, file_number):

        if path is None:
            path = Path(self.file_path.get())

        # Add the sample name from metadata to the folder.
        # path /= experiment.sample
        # Add the name of the device
        path /= self.name

        full_path = str(HDF1_NAME) % (
            str(path), name_base, file_number
        )

        relative_path = str(HDF1_NAME) % (
            self.name, name_base, file_number
        )

        return path, full_path, relative_path

    def setup_images(
            self, path, name_base, file_number, flyscan=False
    ):

        folder, full_path, relative_path = self.setup_file_path_name(
            path, name_base, file_number
        )

        # Setup positioner stream
        if not folder.is_dir():
            folder.mkdir()

        _ps_fname = Path(full_path).relative_to(folder)

        # Setup path and file name in positioner_stream
        self.file_path.put(str(folder))
        self.file_name.put(str(_ps_fname))

        return Path(full_path), Path(relative_path)


positioner_stream = PositionerStream(
    "", name="positioner_stream", labels=("detector",)
)
