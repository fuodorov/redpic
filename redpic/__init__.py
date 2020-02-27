'''
Relativictic Difference Scheme Particle-in-Cell code (REDPIC) init file.

'''

from .constants import *
from .beam import *
from .accelerator import *
from .solver import *

__version__ = '0.1.0'
__doc__ = '''Relativictic Difference Scheme Particle-in-Cell code (REDPIC)'''

__all__ = constants.__all__ + beam.__all__ + solver.__all__ + accelerator.__all__
