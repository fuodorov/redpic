import random

import pytest

from redpic import constants as const
from redpic.beam import UniformBeam


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
            "n": random.randint(1, 100),
        }
        for _ in range(10)
    ],
)
def test_uniform_beam_generate(kwargs):
    beam_kwargs = kwargs.copy()
    del beam_kwargs["n"]
    beam = UniformBeam(**beam_kwargs)
    beam.generate(n=kwargs["n"])
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
    assert beam.n == kwargs["n"]
    assert beam.df.shape[0] == kwargs["n"]
    assert beam.df.shape[1] == 6
    assert beam.df["x"].min() >= -beam.radius_x - beam.x
    assert beam.df["x"].max() <= beam.radius_x + beam.x
    assert beam.df["y"].min() >= -beam.radius_y - beam.y
    assert beam.df["y"].max() <= beam.radius_y + beam.y
    assert beam.df["z"].min() >= -beam.radius_z + beam.z
    assert beam.df["z"].max() <= beam.radius_z + beam.z
