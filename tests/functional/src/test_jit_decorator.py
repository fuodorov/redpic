import timeit

import numpy as np
import pytest

from redpic.core import config as cfg
from redpic.utils import jit

test_arr = np.random.randn(1_000_000)


def expensive_function(arr):
    ans = 0.0
    for a in arr:
        ans += a
    return ans


def test_performance_jit():
    jit_function = jit(expensive_function)
    t1 = timeit.timeit(lambda: expensive_function(test_arr), number=10)
    t2 = timeit.timeit(lambda: jit_function(test_arr), number=10)
    assert t2 != pytest.approx(t1, t1 / 10)


def test_disable_jit():
    cfg.DISABLE_JIT = True
    jit_function = jit(expensive_function)
    t1 = timeit.timeit(lambda: expensive_function(test_arr), number=10)
    t2 = timeit.timeit(lambda: jit_function(test_arr), number=10)
    cfg.DISABLE_JIT = False
    assert t2 == pytest.approx(t1, t1 / 10)
