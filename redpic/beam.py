'''
Relativictic Difference Scheme
Particle-in-Cell code (REDPIC) beam file.

'''

import numpy as np
import pandas as pd
from collections import namedtuple
from .constants import *

__all__ = [ 'Distribution',
            'Beam' ]

Distribution = namedtuple('Distribution', [ 'name', 'x', 'y', 'z', 'px', 'py', 'pz' ])

class Beam:
    '''Creates a beam of selected type particles

    '''

    def __init__(self, type: Element, *, n: int, macro_charge: float=0.0) -> None:

        assert n > 0, 'The number of particles (n) must be a positive number!'

        self.type = type
        self.n = n
        self.x = np.zeros(int(n))
        self.y = np.zeros(int(n))
        self.z = np.zeros(int(n))
        self.px = np.zeros(int(n))
        self.py = np.zeros(int(n))
        self.pz = np.zeros(int(n))
        self.data = np.zeros((6, int(n)))
        if macro_charge == 0.0:
            self.macro_charge = type.charge
        else:
            self.macro_charge = macro_charge
        self.df = pd.DataFrame(np.transpose(self.data), columns=['x','y','z','px','py','pz'])

    def generate(self, distribution: Distribution, *,
                 x_off: float=0.0, y_off: float=0.0, z_off: float=0.0, sig_pz: float=0.1) -> None:
        '''Beam generator

        This function generates a beam with a given distribution and initial beam displacement.
        '''
        assert type(distribution.name) == str, 'Distribution name must be a string!'
        condition = ((distribution.name == 'KV' or distribution.name == 'Uniform') or
                    (distribution.name == 'Gauss' or distribution.name == 'GA'))
        assert condition, 'Сheck distribution name!'

        if distribution.name == 'KV' or distribution.name == 'Uniform':
            phi = np.random.uniform(0, 2*np.pi, int(self.n))
            theta = np.random.uniform(0, np.pi, int(self.n))
            self.x = distribution.x * np.sin(theta) * np.cos(phi)
            self.y = distribution.y * np.sin(theta) * np.sin(phi)
            self.z = distribution.z * np.cos(theta)

        if distribution.name == 'Gauss' or distribution.name == 'GA':
            self.x = np.random.normal(x_off, distribution.x, int(self.n))
            self.y = np.random.normal(y_off, distribution.y, int(self.n))
            self.z = np.random.normal(z_off, distribution.z, int(self.n))

        self.px = np.random.normal(0, distribution.px, int(self.n))
        self.py = np.random.normal(0, distribution.py, int(self.n))
        self.pz = np.random.normal(distribution.pz, distribution.pz*sig_pz, int(self.n))

        self.data = np.row_stack((self.x, self.y, self.z, self.px, self.py, self.pz))
        self.df = pd.DataFrame(np.transpose(self.data), columns=['x','y','z','px','py','pz'])

    def __str__(self):
        return str(self.df)

    def _ipython_display_(self):
        print(str(self.df))
