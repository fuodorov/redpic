from redpic.accelerator import *
from redpic.beam import *
from redpic.core import config
from redpic.solver import *

__version__ = config.PROJECT_VERSION

__doc__ = config.PROJECT_DOC

__all__ = accelerator.__all__ + beam.__all__ + solver.__all__
