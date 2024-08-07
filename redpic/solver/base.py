import logging
from abc import ABC, abstractmethod
from collections.abc import Iterable

from redpic import constants as const
from redpic.accelerator import Accelerator
from redpic.beam.base import BaseBeam

module_logger = logging.getLogger(__name__)


class BaseSimulation(ABC):
    def __init__(
        self,
        beam: BaseBeam | Iterable[BaseBeam],
        accelerator: Accelerator,
        *,
        t_start: float = 0.0e0,
        t_stop: float = 0.0e0,
        dt: float = 0.0e0,
    ):
        self.beams = beam if isinstance(beam, Iterable) else [beam]
        self.acc = accelerator
        self.result = {}
        self.t_start = t_start
        self.t_stop = t_stop if t_stop else (self.acc.z_stop - self.acc.z_start) / const.c
        self.dt = dt if dt else self.acc.dz / const.c

    @abstractmethod
    def _track(self, *, n_files: int) -> None:
        pass

    def track(self, *, n_files: int = 20) -> None:
        if n_files <= 0:
            raise ValueError("The number of files (n_files) must be a positive number!")

        try:
            self._track(n_files=n_files)
        except KeyboardInterrupt:
            module_logger.info("Simulation interrupted by user")
