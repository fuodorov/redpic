from numba import jit as numba_jit
from numba.cuda import jit as numba_cuda_jit

from redpic.core import config as cfg


def jit(func):
    def wrapper(*args, **kwargs):
        if cfg.DISABLE_JIT:
            return func(*args, **kwargs)

        if cfg.ENABLE_CUDA:
            return numba_cuda_jit(
                fastmath=not cfg.DISABLE_FASTMATH,
                cache=not cfg.DISABLE_CACHE,
            )(
                func
            )(*args, **kwargs)

        return numba_jit(
            nopython=not cfg.DISABLE_NOPYTHON,
            parallel=not cfg.DISABLE_PARALLEL,
            fastmath=not cfg.DISABLE_FASTMATH,
            cache=not cfg.DISABLE_CACHE,
            nogil=not cfg.DISABLE_NOGIL,
        )(func)(*args, **kwargs)

    return wrapper
