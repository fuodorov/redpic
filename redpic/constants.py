"""
Relativictic Difference Scheme
Particle-in-Cell code (REDPIC) constants file.

"""

import periodictable

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
            'anti_proton',
            'neutron',
            'anti_neutron' ]

c = speed_of_light = constants.c

ep_0 = epsilon_0 = constants.epsilon_0
mu_0 = constants.mu_0
h = constants.h
hbar = constants.hbar

e = elementary_charge = constants.e

m_e = electron_mass = constants.m_e
m_p = proton_mass = constants.m_p
m_n = neutron_mass = constants.m_n
u = m_u = atomic_constant_mass = constants.physical_constants[ 'atomic mass constant' ][0]

particle = namedtuple('particle', [ 'symbol', 'mass', 'charge' ])

electron = particle(symbol='e', mass=electron_mass, charge=-elementary_charge)
positron = particle(symbol='e+', mass=electron_mass, charge=elementary_charge)
proton = particle(symbol='p', mass=proton_mass, charge=elementary_charge)
anti_proton = particle(symbol='p-', mass=proton_mass, charge=-elementary_charge)
neutron = particle(symbol='n', mass=neutron_mass, charge=0)
anti_neutron = particle(symbol='n', mass=neutron_mass, charge=0)

# Get mass of each element from periodictable
elements = periodictable.core.default_table()
__all__  += periodictable.core.define_elements(elements, globals())
