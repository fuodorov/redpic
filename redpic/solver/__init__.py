from redpic.solver.base import BaseSimulation
from redpic.solver.red import REDSimulation
from redpic.solver.utils import (
    get_field_accelerator,
    get_field_beam,
    sum_field_particles,
)

Simulation = REDSimulation

__all__ = [
    "BaseSimulation",
    "Simulation",
    "REDSimulation",
    "get_field_accelerator",
    "get_field_beam",
    "sum_field_particles",
]
