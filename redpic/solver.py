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
                          x: np.array, y: np.array, z: np.array) -> (np.array,
                                                                     np.array,
                                                                     np.array):
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
        Ez = Ez - d2Ezdz2*r_corr*r_corr/4 - d2Ezdz2*r_corr*r_corr/4           # row remainder
        Ex = - dEzdz*x_corr/2 - dEzdz*x_corr/2 + Ez*offset_xp                # row remainder
        Ey = - dEzdz*y_corr/2 - dEzdz*y_corr/2 + Ez*offset_yp                # row remainder
        return Ex, Ey, Ez
    if type == 'B':
        Bx = acc.Bx(z)
        By = acc.By(z)
        Bz = acc.Bz(z)
        dBzdz = acc.dBzdz(z)
        d2Bzdz2 = misc.derivative(acc.dBzdz, z, dx=dz, n=1)
        Gz = acc.Gz(z)
        Bz = Bz - d2Bzdz2*r_corr*r_corr/4 - d2Bzdz2*r_corr*r_corr/4                # row remainder
        Bx = Bx + Gz*y_corr - dBzdz*x_corr/2 - dBzdz*x_corr/2 + Bz*offset_xp     # row remainder
        By = By + Gz*x_corr - dBzdz*y_corr/2 - dBzdz*y_corr/2 + Bz*offset_yp     # row remainder
        return Bx, By, Bz

@numba.jit(nopython=True, parallel=True, fastmath=True, nogil=True, cache=True)
def sum_field_particle(x: np.array, y: np.array, z: np.array,
                       z_start: float=0.0) -> (np.array, np.array, np.array):
    ''' Sum particles fields

    '''
    n = int(len(x))
    Fx, Fy, Fz = np.zeros(n), np.zeros(n), np.zeros(n)
    r3 = np.zeros(n)
    for i in np.arange(int(n)):
        if z[i] >= z_start:
            r3 = np.sqrt(((x - x[i])*(x - x[i]) + (y - y[i])*(y - y[i]) + (z - z[i])*(z - z[i]))*
                          ((x - x[i])*(x - x[i]) + (y - y[i])*(y - y[i]) + (z - z[i])*(z - z[i]))*
                          ((x - x[i])*(x - x[i]) + (y - y[i])*(y - y[i]) + (z - z[i])*(z - z[i])))
            r3[i] = float(np.inf)
            Fx = Fx + (x - x[i])/r3
            Fy = Fy + (y - y[i])/r3
            Fz = Fz + (z - z[i])/r3
    return Fx, Fy, Fz

def get_field_beam(beam: Beam, acc: Accelerator, type: str,
                   x: np.array, y: np.array, z: np.array) -> (np.array,
                                                              np.array,
                                                              np.array):
    ''' Get an electric [MV/m] or magnetic [T] field space charge at a specific location in the accelerator

    '''
    assert type == 'E' or type == 'B', 'Check type field! (E or B)'
    Q = beam.charge
    q = beam.charge / beam.n
    ke = 1 / (4*np.pi*ep_0*1e6) / 33.3564
    km = mu_0 / (4*np.pi) / 33.3564

    if type == 'E':
        Ex, Ey, Ez = sum_field_particle(x, y, z, acc.z_start)
        return ke*q*Ex, ke*q*Ey, ke*q*Ez
    if type == 'B':
        Bx, By, Bz = sum_field_particle(x, y, z, acc.z_start)
        return km*q*Bx, -km*q*By, 0*Bz

@numba.jit(nopython=True, parallel=True, fastmath=True, nogil=True, cache=True)
def first_step_red(x: np.array, y: np.array, z: np.array,
                   px: np.array, py: np.array, pz: np.array,
                   Fx: np.array, Fy: np.array, Fz: np.array,
                   z_start: float=0.0) -> (np.array, np.array, np.array,
                                           np.array, np.array, np.array):
    ''' First step for Relativictic Difference Scheme

    '''
    n = int(len(x))
    for i in np.arange(int(n)):
        if z[i] >= z_start:
            px[i], py[i], pz[i] = px[i] + 2*Fx[i], py[i] + 2*Fy[i], pz[i] + 2*Fz[i]
            gamma = np.sqrt(1 + px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i])
            vx, vy, vz = px[i]/gamma, py[i]/gamma, pz[i]/gamma
            x[i], y[i], z[i] = x[i] + vx, y[i] + vy, z[i] + vz
        else:
            gamma = np.sqrt(1 + px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i])
            vz = pz[i]/gamma
            x[i], y[i], z[i] = x[i], y[i], z[i] + vz
    return x, y, z, px, py, pz


@numba.jit(nopython=True, parallel=True, fastmath=True, nogil=True, cache=True)
def second_step_red(x: np.array, y: np.array, z: np.array,
                    px: np.array, py: np.array, pz: np.array,
                    Fx: np.array, Fy: np.array, Fz: np.array,
                    z_start: float=0.0) -> (np.array, np.array, np.array,
                                            np.array, np.array, np.array):
    ''' Second step for Relativictic Difference Scheme

    '''
    n = int(len(x))
    for i in np.arange(int(n)):
        if z[i] >= z_start:
            gamma = np.sqrt(1 + px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i])
            vx, vy, vz = px[i]/gamma, py[i]/gamma, pz[i]/gamma
            b2 = 1 + Fx[i]*Fx[i] + Fy[i]*Fy[i] + Fz[i]*Fz[i]
            b1 = 2 - b2
            b3 = 2 * (vx*Fx[i] + vy*Fy[i] + vz*Fz[i])
            fx = 2 * (vy*Fz[i] - vz*Fy[i])
            fy = 2 * (vz*Fx[i] - vx*Fz[i])
            fz = 2 * (vx*Fy[i] - vy*Fx[i])
            vx = (vx*b1 + fx + Fx[i]*b3)/b2
            vy = (vy*b1 + fy + Fy[i]*b3)/b2
            vz = (vz*b1 + fz + Fz[i]*b3)/b2
            x[i], y[i], z[i] = x[i] + vx, y[i] + vy, z[i] + vz
            px[i], py[i], pz[i] = vx*gamma, vy*gamma, vz*gamma
        else:
            gamma = np.sqrt(1 + px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i])
            vz = pz[i]/gamma
            x[i], y[i], z[i] = x[i], y[i], z[i] + vz
    return x, y, z, px, py, pz

class Simulation:
    ''' Simulation of the beam tracking in the accelerator

    '''
    def __init__(self, beam, accelerator):
        self.beam = beam
        self.acc = accelerator

    def track(self, *, n_files: int=30) -> None:
        ''' Tracking!

        '''
        assert n_files > 0, 'The number of files (n_files) must be a positive number!'

        # Constants
        Y = self.beam.da
        Y[2] = Y[2] + self.acc.z_start - max(Y[2]) # set initial beam position
        m = self.beam.type.mass
        q = self.beam.type.charge
        Q = self.beam.charge
        dz = self.acc.dz
        dt = dz/c
        t_max = (self.acc.z_stop-self.acc.z_start)/c + (max(Y[2])-min(Y[2]))/c
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
            if self.beam.charge:
                ex, ey, ez = get_field_beam(self.beam, self.acc, 'E',
                                            Y[0]*dz, Y[1]*dz, Y[2]*dz)
                ex, ey, ez = ex/E0, ey/E0, ez/E0
                Ex, Ey, Ez = Ex + ex, Ey + ey, Ez + ez

            # first step RED
            Y[0], Y[1], Y[2], Y[3], Y[4], Y[5] = first_step_red(Y[0], Y[1], Y[2],
                                                                Y[3], Y[4], Y[5],
                                                                Ex, Ey, Ez,
                                                                self.acc.z_start/dz)
            gamma = np.sqrt(1 + Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])
            vx, vy, vz = Y[3]/gamma, Y[4]/gamma, Y[5]/gamma

            # get magnetic field from accelerator
            Bx, By, Bz = get_field_accelerator(self.acc, 'B',
                                               Y[0]*dz, Y[1]*dz, Y[2]*dz)
            Bx, By, Bz = Bx/B0/gamma, By/B0/gamma, Bz/B0/gamma

            # get magnetic field from beam
            if self.beam.charge:
                bx, by, bz = get_field_beam(self.beam, self.acc, 'B',
                                            Y[0]*dz, Y[1]*dz, Y[2]*dz)
                bx, by, bz = bx*vz/B0, by*vz/B0, bz*vz/B0
                Bx, By, Bz = Bx + bx, By + by, Bz + bz

            # second step RED
            Y[0], Y[1], Y[2], Y[3], Y[4], Y[5] = second_step_red(Y[0], Y[1], Y[2],
                                                                 Y[3], Y[4], Y[5],
                                                                 Bx, By, Bz,
                                                                 self.acc.z_start/dz)
            gamma = np.sqrt(1 + Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])
            vx, vy, vz = Y[3]/gamma, Y[4]/gamma, Y[5]/gamma

            # into the SI
            Y[0], Y[1], Y[2] = Y[0]*dz, Y[1]*dz, Y[2]*dz
            Y[3], Y[4], Y[5] = Y[3]*P0, Y[4]*P0, Y[5]*P0
            Bx, By, Bz = Bx*gamma*B0, By*gamma*B0, Bz*gamma*B0
            Ex, Ey, Ez = Ex*E0, Ey*E0, Ez*E0

            # output in files (symbolBeam.XXXX.csv)
            progress = t / t_max * 100
            meters = self.acc.z_start+t*c
            X = np.row_stack((Y[0], Y[1], Y[2], Y[3], Y[4], Y[5],
                              Bx, By, Bz, Ex, Ey, Ez))
            Xt = np.transpose(X)
            df = pd.DataFrame(Xt, columns=[ 'x','y','z','px','py','pz',
                             'Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez' ])
            if progress % (100 // n_files) <= dt / t_max * 100:
                df.to_csv(self.beam.type.symbol + 'Beam.'+ '%04.0f' % (meters * 100) +'.csv')
            print( '\rz = %.2f m (%.1f %%) ' % (meters, progress), end='')

Sim = Simulation
