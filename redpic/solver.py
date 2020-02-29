'''
Relativictic Difference Scheme Particle-in-Cell code (REDPIC) solver file.

'''

import numba
import numpy as np
import pandas as pd
from scipy import misc
from .constants import *
from .beam import *
from .accelerator import *

__all__ = [ 'Simulation',
            'Sim'
             ]

def get_field_accelerator(acc: Accelerator, type: str,
                          x: np.array, y: np.array, z: np.array) -> (np.array, np.array, np.array):
    ''' Get an electric or magnetic field at a specific location in the accelerator

    '''
    assert type == 'E' or type == 'B', 'Check type field! (E or B)'

    dz = acc.dz
    offset_x = acc.Dx(z)
    offset_xp = acc.Dxp(z)
    offset_y = acc.Dy(z)
    offset_yp = acc.Dyp(z)
    x_corr = x - offset_x
    y_corr = y - offset_y
    r_corr = np.sqrt(x_corr*x_corr + y_corr*y_corr)

    if type == 'E':
        Ez = acc.Ez(z)
        dEzdz = acc.dEzdz(z)
        d2Ezdz2 = misc.derivative(acc.dEzdz, z, dx=dz, n=1)
        Ez = (Ez - d2Ezdz2*r_corr*r_corr/4 - d2Ezdz2*r_corr*r_corr/4)           # row remainder
        Ex = (- dEzdz*x_corr/2 - dEzdz*x_corr/2 + Ez*offset_xp)                 # row remainder
        Ey = (- dEzdz*y_corr/2 - dEzdz*y_corr/2 + Ez*offset_yp)                 # row remainder
        return Ex, Ey, Ez
    if type == 'B':
        Bx = acc.Bx(z)
        By = acc.By(z)
        Bz = acc.Bz(z)
        dBzdz = acc.dBzdz(z)
        d2Bzdz2 = misc.derivative(acc.dBzdz, z, dx=dz, n=1)
        Gz = acc.Gz(z)
        Bz = (Bz - d2Bzdz2*r_corr*r_corr/4 - d2Bzdz2*r_corr*r_corr/4)              # row remainder
        Bx = (Bx + Gz*y_corr - dBzdz*x_corr/2 - dBzdz*x_corr/2 + Bz*offset_xp)     # row remainder
        By = (By + Gz*x_corr - dBzdz*y_corr/2 - dBzdz*y_corr/2 + Bz*offset_yp)     # row remainder
        return Bx, By, Bz

@numba.jit(nopython=True, parallel=True)
def sum_field_particle(n: int, Fx: np.array, Fy: np.array, Fz: np.array,
                       x: np.array, y: np.array, z: np.array) -> (np.array, np.array, np.array):
    ''' Sum particles fields

    '''
    rrr = np.zeros(n)
    for i in np.arange(int(n)):
        rrr = np.sqrt(((x - x[i])*(x - x[i]) + (y - y[i])*(y - y[i]) + (z - z[i])*(z - z[i]))*
                      ((x - x[i])*(x - x[i]) + (y - y[i])*(y - y[i]) + (z - z[i])*(z - z[i]))*
                      ((x - x[i])*(x - x[i]) + (y - y[i])*(y - y[i]) + (z - z[i])*(z - z[i])))
        rrr[i] = float(np.inf)
        Fx = Fx + (x - x[i])/rrr 
        Fy = Fy + (y - y[i])/rrr
        Fz = Fz + (z - z[i])/rrr
    return Fx, Fy, Fz

def get_field_beam(beam: Beam, type: str,
                   x: np.array, y: np.array, z: np.array) -> (np.array, np.array, np.array):
    ''' Get an electric [MV/m] or magnetic [T] field space charge at a specific location in the accelerator

    '''
    assert type == 'E' or type == 'B', 'Check type field! (E or B)'
    Q = beam.macro_charge
    n = int(beam.n)
    ke = 1/(4*np.pi*ep_0*1e6)
    kb = mu_0/(4*np.pi)
    Ex, Ey, Ez = np.zeros(n), np.zeros(n), np.zeros(n)
    Bx, By, Bz = np.zeros(n), np.zeros(n), np.zeros(n)
    if type == 'E':
        Ex, Ey, Ez = sum_field_particle(n, Ex, Ey, Ez, x, y, z)
        return ke*Q*Ex, ke*Q*Ey, ke*Q*Ez
    if type == 'B':
        Bx, By, Bz = sum_field_particle(n, Bx, By, Bz, x, y, z)
        return kb*Q*Bx, -kb*Q*By, 0*Bz

class Simulation:
    ''' Simulation of the beam tracking in the accelerator

    '''
    def __init__(self, beam, accelerator):
        self.beam = beam
        self.acc = accelerator

    def track(self):
        ''' Tracking!

        '''
        # Constants
        m = self.beam.type.mass
        q = self.beam.type.charge
        Q = self.beam.macro_charge
        dz = self.acc.dz
        dt = dz/c
        t_max = (self.acc.z_stop-self.acc.z_start)/c
        Y = self.beam.da
        Y[2] = Y[2] + self.acc.z_start
        P0 = m*c*c / (e*1e6)
        E0 = m*c / (q*dt*1e6)
        B0 = m / (q*dt)

        # RED
        for t in np.arange(0, t_max, 2*dt):

            # into the internal system
            Y[0], Y[1], Y[2] = Y[0]/dz, Y[1]/dz, Y[2]/dz
            Y[3], Y[4], Y[5] = Y[3]/P0, Y[4]/P0, Y[5]/P0

            # get electric field from accelerator
            Ex, Ey, Ez = get_field_accelerator(self.acc, 'E',
                                               Y[0]*dz, Y[1]*dz, Y[2]*dz)
            Ex, Ey, Ez = Ex/E0, Ey/E0, Ez/E0

            # get electric field from beam
            if self.beam.macro_charge:
                ex, ey, ez = get_field_beam(self.beam, 'E',
                                            Y[0]*dz, Y[1]*dz, Y[2]*dz)
                ex, ey, ez = ex/E0, ey/E0, ez/E0
                Ex, Ey, Ez = Ex + ex, Ey + ey, Ez + ez

            # first step in RED
            Y[3], Y[4], Y[5] = Y[3] + 2*Ex, Y[4] + 2*Ey, Y[5] + 2*Ez
            gamma = np.sqrt(1 + Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])
            vx, vy, vz = Y[3]/gamma, Y[4]/gamma, Y[5]/gamma
            Y[0], Y[1], Y[2] = Y[0] + vx, Y[1] + vy, Y[2] + vz

            # get magnetic field from accelerator
            Bx, By, Bz = get_field_accelerator(self.acc, 'B',
                                               Y[0]*dz, Y[1]*dz, Y[2]*dz)
            Bx, By, Bz = Bx/B0/gamma, By/B0/gamma, Bz/B0/gamma

            # get magnetic field from beam
            if self.beam.macro_charge:
                bx, by, bz = get_field_beam(self.beam, 'B',
                                            Y[0]*dz, Y[1]*dz, Y[2]*dz)
                bx, by, bz = bx*vz/B0, by*vz/B0, bz*vz/B0
                Bx, By, Bz = Bx + bx, By + by, Bz + bz

            # second step in RED
            b2 = 1 + Bx*Bx + By*By + Bz*Bz
            b1 = 2 - b2
            b3 = 2 * (vx*Bx + vy*By + vz*Bz)
            fx = 2 * (vy*Bz - vz*By)
            fy = 2 * (vz*Bx - vx*Bz)
            fz = 2 * (vx*By - vy*Bx)
            vx = (vx*b1 + fx + Bx*b3)/b2
            vy = (vy*b1 + fy + By*b3)/b2
            vz = (vz*b1 + fz + Bz*b3)/b2
            Y[0], Y[1], Y[2] = Y[0] + vx, Y[1] + vy, Y[2] + vz
            Y[3], Y[4], Y[5] = vx*gamma, vy*gamma, vz*gamma

            # into the SI
            Y[0], Y[1], Y[2] = Y[0]*dz, Y[1]*dz, Y[2]*dz
            Y[3], Y[4], Y[5] = Y[3]*P0, Y[4]*P0, Y[5]*P0
            Bx, By, Bz = Bx*gamma*B0, By*gamma*B0, Bz*gamma*B0
            Ex, Ey, Ez = Ex*E0, Ey*E0, Ez*E0

            # output in files (symbolBeam.XXXX.csv)
            progress = t/t_max*100
            meters = self.acc.z_start+t*c
            X = np.row_stack((Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Bx, By, Bz, Ex, Ey, Ez))
            Xt = np.transpose(X)
            df = pd.DataFrame(Xt, columns=['x','y','z','px','py','pz',
                             'Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez' ])
            if progress%2 <= dt/t_max*100:
                df.to_csv(self.beam.type.symbol + 'Beam.'+ '%04.0f' % (meters*100) +'.csv')
            print( '\rz = %.2f m (%.1f %%) ' % (meters, progress), end='')

Sim = Simulation
