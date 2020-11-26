"""
Relativictic Difference Scheme Particle-in-Cell code (REDPIC) solver file.

"""

import numpy as np
import pandas as pd
from numba import jit, prange
from multiprocessing import Process
from scipy import misc
from .constants import *
from .beam import *
from .accelerator import *

__all__ = [ 'Simulation',
            'Sim']

def get_field_accelerator(acc: Accelerator, type: str,
                          x: np.array, y: np.array, z: np.array) -> (np.array,
                                                                     np.array,
                                                                     np.array):
    """ Get an electric or magnetic field at a specific location in the accelerator

    """
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
        Ez = Ez - d2Ezdz2*r_corr*r_corr/4         # row remainder
        Ex = - dEzdz*x_corr/2 + Ez*offset_xp                # row remainder
        Ey = - dEzdz*y_corr/2 + Ez*offset_yp                # row remainder
        return Ex, Ey, Ez
    if type == 'B':
        Bx = acc.Bx(z)
        By = acc.By(z)
        Bz = acc.Bz(z)
        dBzdz = acc.dBzdz(z)
        d2Bzdz2 = misc.derivative(acc.dBzdz, z, dx=dz, n=1)
        Gz = acc.Gz(z)
        Bz = Bz - d2Bzdz2*r_corr*r_corr/4              # row remainder
        Bx = Bx + Gz*y_corr - dBzdz*x_corr/2 + Bz*offset_xp     # row remainder
        By = By + Gz*x_corr - dBzdz*y_corr/2 + Bz*offset_yp     # row remainder
        return Bx, By, Bz

@jit(nopython=True, parallel=True, fastmath=True, cache=True, nogil=True)
def sum_field_particles(x: np.array, y: np.array, z: np.array,
                        z_start: float, z_stop: float) -> (np.array, np.array, np.array):
    """ Sum particles fields

    """
    n = int(len(x))
    Fx = np.zeros(n)
    Fy = np.zeros(n)
    Fz = np.zeros(n)
    for i in prange(n):
        if z[i] >= z_start and z[i] <= z_stop:
            for j in prange(n):
                if z[j] >= z_start and z[j] <= z_stop and i != j:
                    Fx[i] += (x[i]-x[j]) / ((x[j] - x[i])**2 + (y[j] - y[i])**2 + (z[j] - z[i])**2)**(3/2)
                    Fy[i] += (y[i]-y[j]) / ((x[j] - x[i])**2 + (y[j] - y[i])**2 + (z[j] - z[i])**2)**(3/2)
                    Fz[i] += (z[i]-z[j]) / ((x[j] - x[i])**2 + (y[j] - y[i])**2 + (z[j] - z[i])**2)**(3/2)
    return Fx, Fy, Fz

def get_field_beam(beam: Beam, acc: Accelerator, type: str,
                   x: np.array, y: np.array, z: np.array) -> (np.array, np.array, np.array):
    """ Get an electric [MV/m] or magnetic [T] field space charge at a specific location in the accelerator

    """
    assert type == 'E' or type == 'B', 'Check type field! (E or B)'
    Q = beam.total_charge
    q = Q / beam.n
    if type == 'E':
        Ex, Ey, Ez = sum_field_particles(x, y, z, acc.z_start, acc.z_stop)
        return ke/(4*np.pi)*q*Ex/1e6, ke/(4*np.pi)*q*Ey/1e6, ke/(4*np.pi)*q*Ez/1e6
    if type == 'B':
        Bx, By, Bz = sum_field_particles(x, y, z, acc.z_start, acc.z_stop)
        return km/(4*np.pi)*q*Bx, -km/(4*np.pi)*q*By, 0*km/(4*np.pi)*q*Bz

@jit(nopython=True, parallel=True, fastmath=True, cache=True, nogil=True)
def first_step_red(x: np.array, y: np.array, z: np.array,
                   px: np.array, py: np.array, pz: np.array,
                   Fx: np.array, Fy: np.array, Fz: np.array,
                   z_start: float, z_stop: float) -> (np.array, np.array, np.array,
                                                      np.array, np.array, np.array):
    """ First step for Relativictic Difference Scheme

    """
    n = int(len(x))
    for i in prange(n):
        if z[i] >= z_start and z[i] <= z_stop:
            px[i] += 2*Fx[i]
            py[i] += 2*Fy[i]
            pz[i] += 2*Fz[i]
            gamma = (1 + px[i]**2 + py[i]**2 + pz[i]**2)**(1/2)
            vx = px[i]/gamma
            vy = py[i]/gamma
            vz = pz[i]/gamma
            x[i] += vx
            y[i] += vy
            z[i] += vz
        else:
            gamma = (1 + px[i]**2 + py[i]**2 + pz[i]**2)**(1/2)
            vz = pz[i]/gamma
            z[i] += vz
    return x, y, z, px, py, pz

@jit(nopython=True, parallel=True, fastmath=True, cache=True, nogil=True)
def second_step_red(x: np.array, y: np.array, z: np.array,
                    px: np.array, py: np.array, pz: np.array,
                    Fx: np.array, Fy: np.array, Fz: np.array,
                    z_start: float, z_stop: float) -> (np.array, np.array, np.array,
                                                       np.array, np.array, np.array):
    """ Second step for Relativictic Difference Scheme

    """
    n = int(len(x))
    for i in prange(n):
        if z[i] >= z_start and z[i] <= z_stop:
            gamma = (1 + px[i]**2 + py[i]**2 + pz[i]**2)**(1/2)
            vx = px[i]/gamma
            vy = py[i]/gamma
            vz = pz[i]/gamma
            b2 = 1 + Fx[i]**2 + Fy[i]**2 + Fz[i]**2
            b1 = 2 - b2
            b3 = 2 * (vx*Fx[i] + vy*Fy[i] + vz*Fz[i])
            fx = 2 * (vy*Fz[i] - vz*Fy[i])
            fy = 2 * (vz*Fx[i] - vx*Fz[i])
            fz = 2 * (vx*Fy[i] - vy*Fx[i])
            vx = (vx*b1 + fx + Fx[i]*b3)/b2
            vy = (vy*b1 + fy + Fy[i]*b3)/b2
            vz = (vz*b1 + fz + Fz[i]*b3)/b2
            x[i] += vx
            y[i] += vy
            z[i] += vz
            px[i] = vx*gamma
            py[i] = vy*gamma
            pz[i] = vz*gamma
        else:
            gamma = (1 + px[i]**2 + py[i]**2 + pz[i]**2)**(1/2)
            vz = pz[i]/gamma
            z[i] += vz
    return x, y, z, px, py, pz

class Simulation:
    """ Simulation of the beam tracking in the accelerator

    """
    def __init__(self, beam, accelerator):
        self.beam = beam
        self.acc = accelerator
        self.result = {}

    def track(self, *, n_files: int=20) -> None:
        """ Tracking!

        """
        assert n_files > 0, 'The number of files (n_files) must be a positive number!'

        # Init parameterss
        Y = self.beam.da
        Y[2] = Y[2] + self.acc.z_start - max(Y[2]) # set initial beam position

        z_start = self.acc.z_start
        z_stop = self.acc.z_stop
        dz = self.acc.dz
        dt = dz / c
        len_beam = max(Y[2]) - min(Y[2])
        t_max = (z_stop - z_start) / c

        m = self.beam.type.mass
        q = self.beam.type.charge

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
            if self.beam.total_charge:
                ex, ey, ez = get_field_beam(self.beam, self.acc, 'E',
                                            Y[0]*dz, Y[1]*dz, Y[2]*dz)
                ex, ey, ez = ex/E0, ey/E0, ez/E0
                Ex, Ey, Ez = Ex + ex, Ey + ey, Ez + ez

            # first step RED
            Y[0], Y[1], Y[2], Y[3], Y[4], Y[5] = first_step_red(Y[0], Y[1], Y[2],
                                                                Y[3], Y[4], Y[5],
                                                                Ex, Ey, Ez,
                                                                z_start/dz, z_stop/dz)
            gamma = np.sqrt(1 + Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])
            vx, vy, vz = Y[3]/gamma, Y[4]/gamma, Y[5]/gamma

            # get magnetic field from accelerator
            Bx, By, Bz = get_field_accelerator(self.acc, 'B',
                                               Y[0]*dz, Y[1]*dz, Y[2]*dz)
            Bx, By, Bz = Bx/B0/gamma, By/B0/gamma, Bz/B0/gamma

            # get magnetic field from beam
            if self.beam.total_charge:
                bx, by, bz = get_field_beam(self.beam, self.acc, 'B',
                                            Y[0]*dz, Y[1]*dz, Y[2]*dz)
                bx, by, bz = bx*vz/B0, by*vz/B0, bz*vz/B0
                Bx, By, Bz = Bx + bx, By + by, Bz + bz

            # second step RED
            Y[0], Y[1], Y[2], Y[3], Y[4], Y[5] = second_step_red(Y[0], Y[1], Y[2],
                                                                 Y[3], Y[4], Y[5],
                                                                 Bx, By, Bz,
                                                                 z_start/dz, z_stop/dz)
            gamma = np.sqrt(1 + Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])
            vx, vy, vz = Y[3]/gamma, Y[4]/gamma, Y[5]/gamma

            # into the SI
            Y[0], Y[1], Y[2] = Y[0]*dz, Y[1]*dz, Y[2]*dz
            Y[3], Y[4], Y[5] = Y[3]*P0, Y[4]*P0, Y[5]*P0
            Bx, By, Bz = Bx*gamma*B0, By*gamma*B0, Bz*gamma*B0
            Ex, Ey, Ez = Ex*E0, Ey*E0, Ez*E0

            progress = t / t_max * 100
            meters = z_start + t * c
            X = np.row_stack((Y[0], Y[1], Y[2], Y[3], Y[4], Y[5],
                              Bx, By, Bz, Ex, Ey, Ez))
            Xt = np.transpose(X)
            df = pd.DataFrame(Xt, columns=[ 'x','y','z','px','py','pz',
                             'Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez' ])
            if progress % (100 // n_files) < 2*dt / t_max * 100:
                self.result.update({round(meters, 3): df})
            print( '\rz = %.2f m (%.1f %%) ' % (meters, progress), end='')

Sim = Simulation
