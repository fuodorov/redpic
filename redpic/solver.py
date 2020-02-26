'''
Relativictic Difference Scheme
Particle-in-Cell code (REDPIC) solver file.

'''

import numpy as np
import pandas as pd
from scipy import misc
from .constants import *
from .beam import *

__all__ = [ 'Simulation',
            'Sim' ]

def red(Y,t=0):
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
    ''' Simulation of the envelope beam in the accelerator

    '''
    def __init__(self, beam, accelerator):
        self.beam = beam
        self.acc = accelerator

    def track(self):
        m = self.beam.type.mass
        q = self.beam.type.charge
        dz = self.acc.dz
        dt = dz/c
        t_max = (self.acc.z_stop-self.acc.z_start)/c
        Y = self.beam.da

        Y[0], Y[1], Y[2] = Y[0]/dz, Y[1]/dz, Y[2]/dz
        Y[3], Y[4], Y[5] = Y[3]*1e6*e/m/c**2, Y[4]*1e6*e/m/c**2, Y[5]*1e6*e/m/c**2
        # RED schema
        for t in np.arange(0, t_max, 2*dt):

            offset_x = self.acc.Dx(Y[2]*dz)
            offset_xp = self.acc.Dxp(Y[2]*dz)
            offset_y = self.acc.Dy(Y[2]*dz)
            offset_yp = self.acc.Dyp(Y[2]*dz)
            x_corr = Y[0]*dz - offset_x
            y_corr = Y[1]*dz - offset_y
            r_corr = np.sqrt((x_corr)**2 + (y_corr)**2)

            Ez = self.acc.Ez(Y[2]*dz)*1e6
            dEzdz = self.acc.dEzdz(Y[2]*dz)*1e6
            d2Ezdz2 = misc.derivative(self.acc.dEzdz, Y[2]*dz, dx=dz, n=1)*1e6
            Ez = Ez - d2Ezdz2*r_corr**2/4 - d2Ezdz2*r_corr**2/4            # row remainder
            Ex = - dEzdz*x_corr/2 - dEzdz*x_corr/2 + Ez*offset_xp          # row remainder
            Ey = - dEzdz*y_corr/2 - dEzdz*y_corr/2 + Ez*offset_yp          # row remainder
            Ez = Ez*(e*dt/m/c/30_000)
            Ex = Ex*(e*dt/m/c/30_000)
            Ey = Ey*(e*dt/m/c/30_000)

            Y[3], Y[4], Y[5] = Y[3] + 2*Ex, Y[4] + 2*Ey, Y[5] + 2*Ez
            gamma = np.sqrt(1 + Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])
            vx, vy, vz = Y[3]/gamma, Y[4]/gamma, Y[5]/gamma
            Y[0], Y[1], Y[2] = Y[0] + vx, Y[1] + vy, Y[2] + vz

            Bx = self.acc.Bx(Y[2]*dz)
            By = self.acc.By(Y[2]*dz)
            Bz = self.acc.Bz(Y[2]*dz)
            dBzdz = self.acc.dBzdz(Y[2]*dz)
            d2Bzdz2 = misc.derivative(self.acc.dBzdz, Y[2]*dz, dx=dz, n=1)
            Gz = self.acc.Gz(Y[2]*dz)
            Bz = Bz - d2Bzdz2*r_corr**2/4 - d2Bzdz2*r_corr**2/4              # row remainder
            Bx = Bx + Gz*y_corr - dBzdz*x_corr/2 - dBzdz*x_corr/2 + Bz*offset_xp     # row remainder
            By = By + Gz*x_corr - dBzdz*y_corr/2 - dBzdz*y_corr/2 + Bz*offset_yp     # row remainder
            Bz = Bz*(e*dt/m/c*1.2566e-2)
            Bx = Bx*(e*dt/m/c*1.2566e-2)
            By = By*(e*dt/m/c*1.2566e-2)
            Bx, By, Bz = Bx/gamma, By/gamma, Bz/gamma

            b2 = 1 + Bx*Bx + By*By + Bz*Bz
            b1 = 2 - b2
            b3 = 2*(vx*Bx + vy*By + vz*Bz)
            fx = 2*(vy*Bz - vz*By)
            fy = 2*(vz*Bx - vx*Bz)
            fz = 2*(vx*By - vy*Bx)

            vx = (vx*b1 + fx + Bx*b3)/b2
            vy = (vy*b1 + fy + By*b3)/b2
            vz = (vz*b1 + fz + Bz*b3)/b2

            Y[0], Y[1], Y[2] = Y[0] + vx, Y[1] + vy, Y[2] + vz
            Y[3], Y[4], Y[5] = vx*gamma, vy*gamma, vz*gamma
            # output:
            df = pd.DataFrame(np.transpose(Y),
                              columns=['x','y','z','px','py','pz'])
            df['x'] = df['x']*dz
            df['y'] = df['y']*dz
            df['z'] = df['z']*dz
            df['px'] = df['px']*m*c**2/e/1e6
            df['py'] = df['py']*m*c**2/e/1e6
            df['pz'] = df['pz']*m*c**2/e/1e6
            #print(df)
            print( "\rz = %g m " % (self.acc.z_start+t*c), end="")

Sim = Simulation
