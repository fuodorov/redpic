'''
Relativictic Difference Scheme Particle-in-Cell code (REDPIC) beam file.

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

    def __init__(self, type: Element, *, charge: float=0.0) -> None:
        self.type = type        # particles type
        self.charge = charge
        self.n = 0.0            # quantity
        self.da = np.array      # data array
        self.df = pd.DataFrame  # data frame

    def generate(self, distribution: Distribution, *, n: float,
                 x_off: float=0.0, y_off: float=0.0, z_off: float=0.0, sig_pz: float=0.01) -> None:
        '''Beam generator

        This function generates a beam with a given distribution and initial beam displacement.
        '''
        assert n > 0, 'The number of particles (n) must be a positive number!'
        assert type(distribution.name) == str, 'Distribution name must be a string!'
        condition = ((distribution.name == 'KV' or distribution.name == 'Uniform') or
                    (distribution.name == 'Gauss' or distribution.name == 'GA'))
        assert condition, 'Сheck distribution name!'

        self.n = n

        if distribution.name == 'Uniform' or distribution.name == 'KV':
            phi = np.random.uniform(0, 2*np.pi, int(self.n))
            x = distribution.x * np.sqrt(np.random.uniform(0, 1, int(self.n))) * np.cos(phi) + x_off
            y = distribution.y * np.sqrt(np.random.uniform(0, 1, int(self.n))) * np.sin(phi) + y_off
            px = distribution.px * np.sqrt(np.random.uniform(0, 1, int(self.n))) * np.cos(phi)
            py = distribution.py * np.sqrt(np.random.uniform(0, 1, int(self.n))) * np.sin(phi)
        if distribution.name == 'Gauss' or distribution.name == 'GA':
            x = np.random.normal(x_off, distribution.x, int(self.n))
            y = np.random.normal(y_off, distribution.y, int(self.n))
            px = np.random.normal(0, distribution.px, int(self.n))
            py = np.random.normal(0, distribution.py, int(self.n))

        z = np.random.uniform(-distribution.z+z_off, distribution.z+z_off, int(self.n)) / 2
        pz = np.random.normal(distribution.pz, distribution.pz * sig_pz, int(self.n))

        self.da = np.row_stack((x, y, z, px, py, pz))
        self.df = pd.DataFrame(np.transpose(self.da), columns=['x','y','z','px','py','pz'])
        self.df.to_csv(self.type.symbol + 'Beam.csv')

    def upload_particles(self, x: np.array, y: np.array, z: np.array,
                         px: np.array, py: np.array, pz: np.array) -> None:
        ''' Particle loading

        Y = np.array[[ ... ] # x[m]
                     [ ... ] # y[m]
                     [ ... ] # z[m]
                     [ ... ] # px[MeV/c]
                     [ ... ] # py[MeV/c]
                     [ ... ] # pz[MeV/c]
                     ]
        '''
        self.da = np.row_stack((x, y, z, px, py, pz))
        self.df = pd.DataFrame(np.transpose(self.da), columns=['x','y','z','px','py','pz'])
        self.df.to_csv(self.type.symbol + 'Beam.csv')

    def __str__(self):
        return str(self.df)

    def _ipython_display_(self):
        print(str(self.df))
