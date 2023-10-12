from logging import config as logging_config

from redpic.core.logger import LOGGING

logging_config.dictConfig(LOGGING)

NAME = "redpic"

VERSION = "0.8.1"

AUTHOR = "Vyacheslav Fedorov"

AUTHOR_EMAIL = "slava@fuodorov.ru"

DOC = "Relativistic Difference Scheme Particle-in-Cell code (REDPIC)"

URL = "https://github.com/fuodorov/redpic"

DISABLE_JIT = False

DISABLE_NOPYTHON = False

DISABLE_PARALLEL = False

DISABLE_FASTMATH = False

DISABLE_CACHE = False

DISABLE_NOGIL = False

ENABLE_CUDA = False

CUDA_THREADS_PER_BLOCK = 1024
