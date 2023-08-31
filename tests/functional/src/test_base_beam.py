import random

import pytest

from redpic import constants as const
from redpic.beam import BaseBeam


class BaseBeamTest(BaseBeam):
    def _generate(self, n: int = 0, *, file_name: str = "", **kwargs) -> None:
        pass


@pytest.mark.parametrize(
    "kwargs",
    [
        {
            "type": const.electron,
            "current": random.random(),
            "energy": random.random(),
            "radius": random.random(),
            "radius_z": random.random(),
            "rp": random.random(),
            "normalized_emittance": random.random(),
            "x": random.random(),
            "y": random.random(),
            "z": random.random(),
            "xp": random.random(),
            "yp": random.random(),
        }
        for _ in range(10)
    ],
)
def test_base_beam_init(kwargs):
    beam = BaseBeamTest(**kwargs)
    assert beam.type == pytest.approx(kwargs["type"])
    assert beam.current == pytest.approx(kwargs["current"])
    assert beam.energy == pytest.approx(kwargs["energy"])
    assert beam.radius == pytest.approx(kwargs["radius"])
    assert beam.radius_z == pytest.approx(kwargs["radius_z"])
    assert beam.rp == pytest.approx(kwargs["rp"])
    assert beam.normalized_emittance == pytest.approx(kwargs["normalized_emittance"])
    assert beam.x == pytest.approx(kwargs["x"])
    assert beam.y == pytest.approx(kwargs["y"])
    assert beam.z == pytest.approx(kwargs["z"])
    assert beam.xp == pytest.approx(kwargs["xp"])
    assert beam.yp == pytest.approx(kwargs["yp"])
    assert beam.n == 0.0
    assert beam.df.empty


def test_base_beam_generate():
    beam = BaseBeamTest()
    with pytest.raises(ValueError):
        beam.generate(n=-1)
