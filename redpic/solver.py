'''
Relativictic Difference Scheme
Particle-in-Cell code (REDPIC) solver file.

'''

import numpy as np
import pandas as pd
from .constants import *
from .beam import *

__all__ = [ 'Simulation',
            'Sim' ]

def F(Y,t=0):
    x  = Y[0]
    y  = Y[1]
    z  = Y[2]
    px = Y[3]
    py = Y[4]
    pz = Y[5]

    p = np.sqrt(px*px + py*py + pz*pz)

    gamma = np.sqrt(1 + p*p)

    vx = px / gamma
    vy = py / gamma
    vz = pz / gamma

    Ex = 0 # 100e6/E0
    Ey = 0
    Ez = 0

    Bx = 0
    By = 0
    Bz = 1.0/B_0 # 1.0 Tesla

    Fx = Ex + vy*Bz - vz*By
    Fy = Ey + vz*Bx - vx*Bz
    Fz = Ez + vx*By - vy*Bx

    return np.array([vx, vy, vz, Fx, Fy, Fz])

class Simulation:
    def __init__(self, beam):
        self.beam = beam

    beam_track = {}

    def track(self):
        Y = self.beam.particles

        cdt = 0.1 # m
        dt = cdt/r_0

        ct_max = 50 # m
        t_max = ct_max/r_0

        for t in np.arange(0,t_max,dt):
            F0 = F(Y,t)

            K1 = F0*dt
            K2 = F(Y + K1/2, t + dt/2)*dt
            K3 = F(Y + K2/2, t + dt/2)*dt
            K4 = F(Y + K3  , t + dt)*dt

            Y  = Y + (K1 + 2*K2 + 2*K3 + K4)/6

            # output:
            df = pd.DataFrame(np.transpose(Y), columns=['x','y','z','px','py','pz'])

            self.beam_track[t*r_0] = df


Sim = Simulation
