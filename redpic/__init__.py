from redpic.accelerator import __all__ as _accelerator_all
from redpic.beam import __all__ as _beam_all
from redpic.constants import __all__ as _constants_all
from redpic.core import config
from redpic.solver import __all__ as _solver_all

__version__ = config.PROJECT_VERSION

__doc__ = config.PROJECT_DOC

__all__ = _constants_all + _accelerator_all + _beam_all + _solver_all
