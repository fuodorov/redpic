"""
Relativictic Difference Scheme Particle-in-Cell code (REDPIC) init file.

"""

from .accelerator import *
from .beam import *
from .constants import *
from .solver import *

__version__ = "1.0.0.1"
__doc__ = """Relativictic Difference Scheme Particle-in-Cell code (REDPIC)"""

__all__ = constants.__all__ + beam.__all__ + solver.__all__ + accelerator.__all__
