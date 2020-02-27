'''
Relativictic Difference Scheme Particle-in-Cell code (REDPIC) solver file.

'''

import numpy as np
import pandas as pd
from scipy import misc
from .constants import *
from .beam import *

__all__ = [ 'Simulation',
            'Sim'
             ]

class Simulation:
    ''' Simulation of the envelope beam in the accelerator

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
        dz = self.acc.dz
        dt = dz/c
        t_max = (self.acc.z_stop-self.acc.z_start)/c
        Y = self.beam.da
        P0 = m*c*c / (e*1e6)
        E0 = m*c / (e*dt*1e6)
        B0 = m / (e*dt)

        # RED
        for t in np.arange(0, t_max, 2*dt):
            Y[0], Y[1], Y[2] = Y[0]/dz, Y[1]/dz, Y[2]/dz
            Y[3], Y[4], Y[5] = Y[3]/P0, Y[4]/P0, Y[5]/P0


            offset_x = self.acc.Dx(Y[2]*dz)
            offset_xp = self.acc.Dxp(Y[2]*dz)
            offset_y = self.acc.Dy(Y[2]*dz)
            offset_yp = self.acc.Dyp(Y[2]*dz)
            x_corr = Y[0]*dz - offset_x
            y_corr = Y[1]*dz - offset_y
            r_corr = np.sqrt(x_corr*x_corr + y_corr*y_corr)
            Ez = self.acc.Ez(Y[2]*dz)
            dEzdz = self.acc.dEzdz(Y[2]*dz)
            d2Ezdz2 = misc.derivative(self.acc.dEzdz, Y[2]*dz, dx=dz, n=1)
            Ez = (Ez - d2Ezdz2*r_corr*r_corr/4 - d2Ezdz2*r_corr*r_corr/4)/E0            # row remainder
            Ex = (- dEzdz*x_corr/2 - dEzdz*x_corr/2 + Ez*offset_xp)/E0          # row remainder
            Ey = (- dEzdz*y_corr/2 - dEzdz*y_corr/2 + Ez*offset_yp)/E0          # row remainder


            Y[3], Y[4], Y[5] = Y[3] - 2*Ex, Y[4] - 2*Ey, Y[5] - 2*Ez
            gamma = np.sqrt(1 + Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])
            vx, vy, vz = Y[3]/gamma, Y[4]/gamma, Y[5]/gamma
            Y[0], Y[1], Y[2] = Y[0] + vx, Y[1] + vy, Y[2] + vz


            offset_x = self.acc.Dx(Y[2]*dz)
            offset_xp = self.acc.Dxp(Y[2]*dz)
            offset_y = self.acc.Dy(Y[2]*dz)
            offset_yp = self.acc.Dyp(Y[2]*dz)
            x_corr = Y[0]*dz - offset_x
            y_corr = Y[1]*dz - offset_y
            r_corr = np.sqrt(x_corr*x_corr + y_corr*y_corr)
            Bx = self.acc.Bx(Y[2]*dz)
            By = self.acc.By(Y[2]*dz)
            Bz = self.acc.Bz(Y[2]*dz)
            dBzdz = self.acc.dBzdz(Y[2]*dz)
            d2Bzdz2 = misc.derivative(self.acc.dBzdz, Y[2]*dz, dx=dz, n=1)
            Gz = self.acc.Gz(Y[2]*dz)
            Bz = (Bz - d2Bzdz2*r_corr*r_corr/4 - d2Bzdz2*r_corr*r_corr/4)/B0/gamma              # row remainder
            Bx = (Bx + Gz*y_corr - dBzdz*x_corr/2 - dBzdz*x_corr/2 + Bz*offset_xp)/B0/gamma     # row remainder
            By = (By + Gz*x_corr - dBzdz*y_corr/2 - dBzdz*y_corr/2 + Bz*offset_yp)/B0/gamma     # row remainder


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

            Y[0], Y[1], Y[2] = Y[0]*dz, Y[1]*dz, Y[2]*dz
            Y[3], Y[4], Y[5] = Y[3]*P0, Y[4]*P0, Y[5]*P0
            Bx, By, Bz = Bx*gamma*B0, By*gamma*B0, Bz*gamma*B0
            Ex, Ey, Ez = Ex*E0, Ey*E0, Ez*E0

            X = np.row_stack((Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Bx, By, Bz, Ex, Ey, Ez))
            Xt = np.transpose(X)
            df = pd.DataFrame(Xt, columns=['x','y','z','px','py','pz',
                             'Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez' ])
            progress = t/t_max*100
            meters = self.acc.z_start+t*c
            print( '\rz = %.2f m (%.1f %%) ' % (meters, progress), end='')
            if progress%3 == 0 or progress%4 == 0:
                df.to_csv(self.beam.type.symbol + 'Beam.'+ '%04.0f' % (meters*100) +'.csv')

Sim = Simulation
