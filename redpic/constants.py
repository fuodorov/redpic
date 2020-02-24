"""
Relativictic Difference Scheme
Particle-in-Cell code (REDPIC) constants file.

"""

import periodictable
import numpy as np
from collections import namedtuple
from scipy import constants

__all__ = [ 'speed_of_light',
            'c',
            'epsilon_0',
            'ep_0',
            'mu_0',
            'h',
            'hbar',
            'elementary_charge',
            'e',
            'electron_mass',
            'm_e',
            'electron_mass_energy',
            'mc',
            'electron_radius',
            'r_0',
            'proton_mass',
            'm_p',
            'neutron_mass',
            'm_n',
            'atomic_constant_mass',
            'm_u',
            'u',
            'electron',
            'positron',
            'proton',
            'antiproton',
            'neutron',
            'antineutron',
            'Element',
            'very_few',
            'few',
            'average',
            'many',
            'very_many',
            'E_0', 'B_0' ]

# Constants
c = speed_of_light = constants.c

ep_0 = epsilon_0 = constants.epsilon_0
mu_0 = constants.mu_0
h = constants.h
hbar = constants.hbar

e = elementary_charge = constants.e

m_e = electron_mass = constants.m_e
mc = electron_mass_energy = constants.physical_constants[ 'electron mass energy equivalent in MeV' ][0]
r_0 = electron_radius = constants.physical_constants[ 'classical electron radius' ][0]
m_p = proton_mass = constants.m_p
m_n = neutron_mass = constants.m_n
u = m_u = atomic_constant_mass = constants.physical_constants[ 'atomic mass constant' ][0]

Element = namedtuple('Element', [ 'name', 'symbol', 'mass', 'charge' ])

electron = Element(name='electron', symbol='e', mass=m_e, charge=-e)
positron = Element(name='positron', symbol='e+', mass=m_e, charge=e)
proton = Element(name='proton', symbol='p', mass=m_p, charge=e)
antiproton = Element(name='antiproton', symbol='p-', mass=m_p, charge=-e)
neutron = Element(name='neutron', symbol='n', mass=m_n, charge=0)
antineutron = Element(name='antineutron', symbol='n', mass=m_n, charge=0)

##########!!!########
E_0_cgs = 4.803e-10/(r_0*r_0*1e4) # statvolt/cm
E_0 = E_0_cgs*299.792458*100 # V/m
B_0 = E_0_cgs/10_000 # Tesla

# Get mass of each element from periodictable
elements = periodictable.core.default_table()
__all__  += periodictable.core.define_elements(elements, globals())

# Size
very_few = 1_000
few = 10_000
average = 100_000
many = 1_000_000
very_many = 10_000_000
