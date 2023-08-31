from abc import ABC, abstractmethod

from redpic.accelerator import Accelerator
from redpic.beam.base import BaseBeam


class BaseSimulation(ABC):
    def __init__(self, beam: BaseBeam, accelerator: Accelerator):
        self.beam = beam
        self.acc = accelerator
        self.result = {}

    @abstractmethod
    def track(self, *, n_files: int = 20) -> None:
        pass
