import numpy as np
import pytest
from scipy import constants as scipy_const

from redpic import constants as redpic_const


@pytest.mark.parametrize(
    "name, value",
    [
        ("c", scipy_const.c),
        ("speed_of_light", scipy_const.c),
        ("epsilon_0", scipy_const.epsilon_0),
        ("ep_0", scipy_const.epsilon_0),
        ("mu_0", scipy_const.mu_0),
        ("ke", 1 / (4 * np.pi * scipy_const.epsilon_0)),
        ("km", scipy_const.mu_0 / (4 * np.pi)),
        ("h", scipy_const.h),
        ("hbar", scipy_const.hbar),
        ("elementary_charge", scipy_const.e),
        ("e", scipy_const.e),
        ("electron_mass", scipy_const.m_e),
        ("m_e", scipy_const.m_e),
        ("electron_mass_energy", scipy_const.physical_constants["electron mass energy equivalent in MeV"][0]),
        ("mc", scipy_const.physical_constants["electron mass energy equivalent in MeV"][0]),
        ("electron_radius", scipy_const.physical_constants["classical electron radius"][0]),
        ("r_0", scipy_const.physical_constants["classical electron radius"][0]),
        ("proton_mass", scipy_const.m_p),
        ("m_p", scipy_const.m_p),
        ("neutron_mass", scipy_const.m_n),
        ("m_n", scipy_const.m_n),
        ("atomic_constant_mass", scipy_const.physical_constants["atomic mass constant"][0]),
        ("m_u", scipy_const.physical_constants["atomic mass constant"][0]),
        ("u", scipy_const.physical_constants["atomic mass constant"][0]),
    ],
)
def test_constants(name, value):
    assert getattr(redpic_const, name) == value


@pytest.mark.parametrize(
    "name, symbol, mass, charge",
    [
        ("electron", "e", scipy_const.m_e, -scipy_const.e),
        ("positron", "e+", scipy_const.m_e, scipy_const.e),
        ("proton", "p", scipy_const.m_p, scipy_const.e),
        ("antiproton", "p-", scipy_const.m_p, -scipy_const.e),
        ("neutron", "n", scipy_const.m_n, 0),
        ("antineutron", "n", scipy_const.m_n, 0),
    ],
)
def test_particles(name, symbol, mass, charge):
    particle = getattr(redpic_const, name)
    assert particle.name == name
    assert particle.symbol == symbol
    assert particle.mass == mass
    assert particle.charge == charge


def test_element():
    element = redpic_const.Element(name="test", symbol="T", mass=1, charge=2)
    assert element.name == "test"
    assert element.symbol == "T"
    assert element.mass == 1
    assert element.charge == 2


def test_particle():
    particle = redpic_const.Particle(name="test", symbol="T", mass=1, charge=2)
    assert particle.name == "test"
    assert particle.symbol == "T"
    assert particle.mass == 1
    assert particle.charge == 2
