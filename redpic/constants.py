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

element = namedtuple('element', [ 'name', 'symbol', 'mass', 'charge' ])

electron = element(name='electron', symbol='e', mass=m_e, charge=-e)
positron = element(name='positron', symbol='e+', mass=m_e, charge=e)
proton = element(name='proton', symbol='p', mass=m_p, charge=e)
anti_proton = element(name='antiproton', symbol='p-', mass=m_p, charge=-e)
neutron = element(name='neutron', symbol='n', mass=m_n, charge=0)
anti_neutron = element(name='antineutron', symbol='n', mass=m_n, charge=0)

# Get mass of each element from periodictable
elements = periodictable.core.default_table()
__all__  += periodictable.core.define_elements(elements, globals())
