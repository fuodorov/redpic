from collections import namedtuple

import numpy as np
import periodictable
from scipy import constants as const

__all__ = [
    "speed_of_light",
    "c",
    "epsilon_0",
    "ep_0",
    "mu_0",
    "ke",
    "km",
    "h",
    "hbar",
    "elementary_charge",
    "e",
    "electron_mass",
    "m_e",
    "electron_mass_energy",
    "mc",
    "electron_radius",
    "r_0",
    "proton_mass",
    "m_p",
    "neutron_mass",
    "m_n",
    "atomic_constant_mass",
    "m_u",
    "u",
    "electron",
    "positron",
    "proton",
    "antiproton",
    "neutron",
    "antineutron",
    "Element",
    "Particle",
]

# Constants
c = speed_of_light = const.c

ep_0 = epsilon_0 = const.epsilon_0
mu_0 = const.mu_0
h = const.h
hbar = const.hbar

ke = 1 / (4 * np.pi * ep_0)
km = mu_0 / (4 * np.pi)

e = elementary_charge = const.e

m_e = electron_mass = const.m_e
mc = electron_mass_energy = const.physical_constants["electron mass energy equivalent in MeV"][0]
r_0 = electron_radius = const.physical_constants["classical electron radius"][0]
m_p = proton_mass = const.m_p
m_n = neutron_mass = const.m_n
u = m_u = atomic_constant_mass = const.physical_constants["atomic mass constant"][0]

Element = namedtuple("Element", ["name", "symbol", "mass", "charge"])

electron = Element(name="electron", symbol="e", mass=m_e, charge=-e)
positron = Element(name="positron", symbol="e+", mass=m_e, charge=e)
proton = Element(name="proton", symbol="p", mass=m_p, charge=e)
antiproton = Element(name="antiproton", symbol="p-", mass=m_p, charge=-e)
neutron = Element(name="neutron", symbol="n", mass=m_n, charge=0)
antineutron = Element(name="antineutron", symbol="n", mass=m_n, charge=0)
Particle = Element

# Get mass of each element from periodictable
elements = periodictable.core.default_table()
__all__ += periodictable.core.define_elements(elements, globals())
