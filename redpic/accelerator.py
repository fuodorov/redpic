"""
Relativictic Difference Scheme Particle-in-Cell code (REDPIC) accelerator file.

"""

import numpy as np
from scipy import interpolate, integrate
from .constants import *
from .beam import *
from kenv.accelerator import *

__all__ = ['Element',
           'Accelerator',
           'read_fields',
           'read_offsets']
