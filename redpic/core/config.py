from logging import config as logging_config

from .logger import LOGGING

logging_config.dictConfig(LOGGING)

PROJECT_NAME = "redpic"

PROJECT_VERSION = "0.7.29"

PROJECT_AUTHOR = "Vyacheslav Fedorov"

PROJECT_AUTHOR_EMAIL = "slava@fuodorov.ru"

PROJECT_DOC = "Relativistic Difference Scheme Particle-in-Cell code (REDPIC)"

PROJECT_URL = "https://github.com/fuodorov/redpic"
