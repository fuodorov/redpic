import logging
from abc import ABC, abstractmethod

from redpic.accelerator import Accelerator
from redpic.beam.base import BaseBeam

module_logger = logging.getLogger(__name__)


class BaseSimulation(ABC):
    def __init__(self, beam: BaseBeam, accelerator: Accelerator):
        self.beam = beam
        self.acc = accelerator
        self.result = {}

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
