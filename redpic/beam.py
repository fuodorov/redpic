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

    def __init__(self, type: Element) -> None:
        self.type = type        # particles type
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
            theta = np.random.uniform(0, np.pi, int(self.n))
            x = distribution.x * np.sin(theta) * np.cos(phi)
            y = distribution.y * np.sin(theta) * np.sin(phi)
            z = distribution.z * np.cos(theta)

        if distribution.name == 'Gauss' or distribution.name == 'GA':
            x = np.random.normal(x_off, distribution.x, int(self.n))
            y = np.random.normal(y_off, distribution.y, int(self.n))
            z = np.random.normal(z_off, distribution.z, int(self.n))

        px = np.random.normal(0, distribution.px, int(self.n))
        py = np.random.normal(0, distribution.py, int(self.n))
        pz = np.random.normal(distribution.pz, distribution.pz * sig_pz, int(self.n))

        self.da = np.row_stack((x, y, z, px, py, pz))
        self.df = pd.DataFrame(np.transpose(self.da), columns=['x','y','z','px','py','pz'])
        self.df.to_csv(self.type.symbol + 'Beam.csv')

    def upload(self, Y: np.array) -> None:
        ''' Particle loading

        Y = np.array[[ ... ] # x[m]
                     [ ... ] # y[m]
                     [ ... ] # z[m]
                     [ ... ] # px[MeV/c]
                     [ ... ] # py[MeV/c]
                     [ ... ] # pz[MeV/c]
                     ]
        '''
        self.da = Y
        self.df = pd.DataFrame(np.transpose(self.da), columns=['x','y','z','px','py','pz'])

    def __str__(self):
        return str(self.df)

    def _ipython_display_(self):
        print(str(self.df))
