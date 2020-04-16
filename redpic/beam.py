'''
Relativictic Difference Scheme Particle-in-Cell code (REDPIC) beam file.

'''

import numpy as np
import pandas as pd
from collections import namedtuple
from .constants import *

__all__ = [ 'Distribution',
            'Beam' ]

def read_distribution_file(fname):
    ''' Read distribution from file .csv

    '''
    cols = ['x', 'y', 'z', 'px', 'py', 'pz']
    #        m    m    m    MeV/c  MeV/c  MeV/c

    df = pd.read_csv(fname, dtype='float32')
    return df['x'], df['y'], df['z'], df['px'], df['py'], df['pz']


def read_distribution_astra(fname):
    ''' Read distribution from Astra file

    '''
    cols = ['x', 'y', 'z', 'px', 'py', 'pz', 'clock', 'charge', 'id', 'flag']
    #        m    m    m    eV/c  eV/c  eV/c  ns       nC
    df = pd.read_csv(fname, header=None, delim_whitespace=True, names=cols, dtype='float32')

    df = df[df.flag != -15] # ignore the lost particles
    df['px'] = df['px']/1e6 # MeV/c
    df['py'] = df['py']/1e6 # MeV/c

    # remove the reference particle
    df0 = df.head(1)
    df  = df.drop(df0.index)

    z0  = df0.z.values[0]
    pz0 = df0.pz.values[0]

    # Recalculate z and pz:
    df['z'] = z0 + df['clock']*1e-9*c # m
    df['pz'] = (pz0 + df['pz'])/1e6 # MeV/c
    return df['x'], df['y'], df['z'], df['px'], df['py'], df['pz']

Distribution = namedtuple('Distribution', [ 'name', 'x', 'y', 'z', 'px', 'py', 'pz' ])

class Beam:
    '''Creates a beam of selected type particles

    '''

    def __init__(self, type: Element, *, charge: float=0.0) -> None:
        self.type = type        # particles type
        self.charge = charge    # beam charge
        self.n = 0.0            # quantity
        self.da = np.array      # data array
        self.df = pd.DataFrame  # data frame

    def generate(self, distribution: Distribution, *, n: float,
                 x_off: float=0.0, y_off: float=0.0, z_off: float=0.0, sig_pz: float=0.001, path: str='') -> None:
        '''Beam generator

        This function generates a beam with a given distribution and initial beam displacement.
        '''
        condition = ((distribution.name == 'KV' or distribution.name == 'Uniform') or
                    (distribution.name == 'Gauss' or distribution.name == 'GA'))

        assert n > 0, 'The number of particles (n) must be a positive number!'
        assert type(distribution.name) == str, 'Distribution name must be a string!'
        assert condition, 'Сheck distribution name!'

        self.n = n

        if distribution.name == 'Uniform' or distribution.name == 'KV':
            s = np.random.normal(0, distribution.x, int(self.n))
            t = np.random.normal(0, distribution.y, int(self.n))
            u = np.random.normal(0, distribution.x, int(self.n))
            v = np.random.normal(0, distribution.y, int(self.n))
            norm = (s*s + t*t + u*u + v*v)**0.5
            (x,y) = (u+x_off, v+y_off) / norm

            s = np.random.normal(0, distribution.px, int(self.n))
            t = np.random.normal(0, distribution.py, int(self.n))
            u = np.random.normal(0, distribution.px, int(self.n))
            v = np.random.normal(0, distribution.py, int(self.n))
            norm = (s*s + t*t + u*u + v*v)**0.5
            (px,py) = (u,v) / norm
        if distribution.name == 'Gauss' or distribution.name == 'GA':
            x = np.random.normal(x_off, distribution.x, int(self.n))
            y = np.random.normal(y_off, distribution.y, int(self.n))
            px = np.random.normal(0, distribution.px, int(self.n))
            py = np.random.normal(0, distribution.py, int(self.n))

        z = np.random.uniform(-distribution.z+z_off, distribution.z+z_off, int(self.n)) / 2
        pz = np.random.normal(distribution.pz, distribution.pz * sig_pz, int(self.n))

        self.da = np.row_stack((x, y, z, px, py, pz))
        self.df = pd.DataFrame(np.transpose(self.da), columns=['x','y','z','px','py','pz'])
        self.df.to_csv(path + self.type.symbol + 'Beam.csv')

    def upload(self, file_name: str, *, path: str=''):
        ''' Particle loading

        '''
        file_extension = file_name.split('.')[-1]

        if file_extension == 'ini':
            x, y, z, px, py, pz = read_distribution_astra(file_name)
            self.n = int(len(x))
            self.da = np.row_stack((x, y, z, px, py, pz))
            self.df = pd.DataFrame(np.transpose(self.da), columns=['x','y','z','px','py','pz'])
            self.df.to_csv(path + self.type.symbol + 'Beam.csv')
        if file_extension == 'csv':
            x, y, z, px, py, pz = read_distribution_file(file_name)
            self.n = int(len(x))
            self.da = np.row_stack((x, y, z, px, py, pz))
            self.df = pd.DataFrame(np.transpose(self.da), columns=['x','y','z','px','py','pz'])
            self.df.to_csv(path + self.type.symbol + 'Beam.csv')

    def __str__(self):
        return str(self.df)

    def _ipython_display_(self):
        print(str(self.df))
