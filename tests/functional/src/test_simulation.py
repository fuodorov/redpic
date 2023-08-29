import pytest
from scipy import stats

from redpic import constants as const
from redpic.accelerator import Accelerator
from redpic.beam import AstraBeam
from redpic.solver import Simulation
from tests.functional.utils.astra import read_track_astra_particles


@pytest.fixture
def test_accelerator():
    acc = Accelerator(0.7, 1.7, 0.01)
    acc.add_accel("Acc. 1", 4.096, -1.1, "tests/functional/test_data/data/fields/Ez.dat")
    acc.add_accel("Acc. 2", 5.944, -1.1, "tests/functional/test_data/data/fields/Ez.dat")
    acc.add_solenoid("Sol. 1", 0.450, -0.0580, "tests/functional/test_data/data/fields/Bz.dat")
    acc.add_solenoid("Sol. 2", 0.957, 0.0390, "tests/functional/test_data/data/fields/Bz.dat")
    acc.add_solenoid("Sol. 3", 2.107, 0.0250, "tests/functional/test_data/data/fields/Bz.dat")
    acc.add_solenoid("Sol. 4", 2.907, 0.0440, "tests/functional/test_data/data/fields/Bz.dat")
    acc.add_solenoid("Sol. 5", 3.670, 0.0400, "tests/functional/test_data/data/fields/Bz.dat")
    acc.add_solenoid("Sol. 6", 4.570, 0.0595, "tests/functional/test_data/data/fields/Bz.dat")
    acc.add_solenoid("Sol. 7", 5.470, 0.0590, "tests/functional/test_data/data/fields/Bz.dat")
    acc.compile()
    return acc


@pytest.fixture
def test_beam():
    beam = AstraBeam(
        type=const.electron,
        energy=1.32,  # MeV
        current=0.5e3,  # A
        radius_x=48e-3,  # initial r (m)
        radius_y=48e-3,  # initial r (m)
        radius_z=3.5,
        radius_xp=2 * 35.0e-3,  # initial r' (rad)
        radius_yp=2 * 35.0e-3,  # initial r' (rad)
        x=0.0e-3,  # horizontal centroid position (m)
        xp=0.0e-3,  # horizontal centroid angle (rad)
        y=0,  # vertical centroid position (m)
        normalized_emittance=200e-6,  # m*rad
    )
    beam.generate(file_name="tests/functional/test_data/data/astra/beam.ini")
    return beam


def test_simulation_vs_astra(test_beam, test_accelerator):
    redpic_sim = Simulation(test_beam, test_accelerator)
    redpic_sim.track()
    redpic_df = redpic_sim.result[list(redpic_sim.result.keys())[-1]]
    astra_df = read_track_astra_particles("tests/functional/test_data/data/astra/track.data")

    redpic_df = redpic_df[test_accelerator.z_start < redpic_df["z"]]
    redpic_df = redpic_df[redpic_df["z"] < test_accelerator.z_stop]
    astra_df = astra_df[test_accelerator.z_start < astra_df["z"]]
    astra_df = astra_df[astra_df["z"] < test_accelerator.z_stop]

    assert stats.ks_2samp(redpic_df["x"], astra_df["x"]).pvalue > 0.01
    assert stats.ks_2samp(redpic_df["y"], astra_df["y"]).pvalue > 0.01
