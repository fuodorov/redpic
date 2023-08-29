import pytest

from redpic.beam import AstraBeam
from tests.functional.utils.astra import read_beam_astra_particles


@pytest.mark.parametrize(
    "file_name",
    [
        "tests/functional/test_data/data/astra/beam.ini",
    ],
)
def test_astra_beam_generate(file_name):
    beam = AstraBeam()
    beam.generate(file_name=file_name)
    df = read_beam_astra_particles(file_name)
    assert beam.df is not None
    assert beam.df["x"].values.tolist() == df["x"].values.tolist()
    assert beam.df["y"].values.tolist() == df["y"].values.tolist()
    assert beam.df["z"].values.tolist() == df["z"].values.tolist()
    assert beam.df["px"].values.tolist() == df["px"].values.tolist()
    assert beam.df["py"].values.tolist() == df["py"].values.tolist()
    assert beam.df["pz"].values.tolist() == df["pz"].values.tolist()
