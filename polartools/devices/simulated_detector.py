
__all__ = ['simdet']

from ophyd.signal import SignalRO
from numpy.random import default_rng
from time import sleep
from ..utils import logger
logger.info(__file__)


class RandomDet(SignalRO):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._gen = default_rng()
        self._count_time = 1

    def get(self, **kwargs):
        sleep(self._count_time)
        return self._gen.random()

    @property
    def count_time(self):
        return self._count_time

    @count_time.setter
    def count_time(self, value):
        if isinstance(value, (int, float)):
            self._count_time = value
        else:
            raise ValueError("Count time must be a number.")


simdet = RandomDet(name="simdet")
