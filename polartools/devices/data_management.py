from apstools.devices import DM_WorkflowConnector
from apstools.utils import dm_api_proc
from ophyd import Signal
from ..utils import logger
logger.info(__file__)

__all__ = """
    dm_experiment
    dm_workflow
""".split()


dm_workflow = DM_WorkflowConnector(name="dm_workflow", labels=("dm",))
dm_workflow.owner.put(dm_api_proc().username)

# TODO: make this an EpicsSignal instead
dm_experiment = Signal(name="dm_experiment", value="", labels=("dm",))
