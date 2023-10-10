import numpy as np
import pandas as pd
from numba import cuda, prange

from redpic import constants as const
from redpic.core import config as cfg
from redpic.solver.base import BaseSimulation
from redpic.utils.field import get_field_accelerator, get_field_beam
from redpic.utils.jit import jit


class REDSimulation(BaseSimulation):
    @staticmethod
    @jit
    def _first_step_red(
        x: np.array,
        y: np.array,
        z: np.array,
        px: np.array,
        py: np.array,
        pz: np.array,
        Fx: np.array,
        Fy: np.array,
        Fz: np.array,
        z_start: float,
        z_stop: float,
    ) -> None:
        for i in prange(int(len(x))):  # pylint: disable=E1133
            if z_start <= z[i] <= z_stop:
                px[i] += 2 * Fx[i]
                py[i] += 2 * Fy[i]
                pz[i] += 2 * Fz[i]
                gamma = (1 + px[i] ** 2 + py[i] ** 2 + pz[i] ** 2) ** (1 / 2)
                vx = px[i] / gamma
                vy = py[i] / gamma
                vz = pz[i] / gamma
                x[i] += vx
                y[i] += vy
                z[i] += vz
            else:
                gamma = (1 + px[i] ** 2 + py[i] ** 2 + pz[i] ** 2) ** (1 / 2)
                vz = pz[i] / gamma
                z[i] += vz

    @staticmethod
    @cuda.jit
    def _first_step_red_cuda(
        x: np.array,
        y: np.array,
        z: np.array,
        px: np.array,
        py: np.array,
        pz: np.array,
        Fx: np.array,
        Fy: np.array,
        Fz: np.array,
        z_start: float,
        z_stop: float,
    ) -> None:
        i = cuda.grid(1)  # pylint: disable=no-value-for-parameter
        if z_start <= z[i] <= z_stop and i < x.size:
            px[i] += 2 * Fx[i]
            py[i] += 2 * Fy[i]
            pz[i] += 2 * Fz[i]
            gamma = (1 + px[i] ** 2 + py[i] ** 2 + pz[i] ** 2) ** (1 / 2)
            vx = px[i] / gamma
            vy = py[i] / gamma
            vz = pz[i] / gamma
            x[i] += vx
            y[i] += vy
            z[i] += vz
        else:
            gamma = (1 + px[i] ** 2 + py[i] ** 2 + pz[i] ** 2) ** (1 / 2)
            vz = pz[i] / gamma
            z[i] += vz

    @staticmethod
    @jit
    def _second_step_red(
        x: np.array,
        y: np.array,
        z: np.array,
        px: np.array,
        py: np.array,
        pz: np.array,
        Fx: np.array,
        Fy: np.array,
        Fz: np.array,
        z_start: float,
        z_stop: float,
    ) -> None:
        for i in prange(int(len(x))):  # pylint: disable=E1133
            if z_start <= z[i] <= z_stop:
                gamma = (1 + px[i] ** 2 + py[i] ** 2 + pz[i] ** 2) ** (1 / 2)
                vx = px[i] / gamma
                vy = py[i] / gamma
                vz = pz[i] / gamma
                b2 = 1 + Fx[i] ** 2 + Fy[i] ** 2 + Fz[i] ** 2
                b1 = 2 - b2
                b3 = 2 * (vx * Fx[i] + vy * Fy[i] + vz * Fz[i])
                fx = 2 * (vy * Fz[i] - vz * Fy[i])
                fy = 2 * (vz * Fx[i] - vx * Fz[i])
                fz = 2 * (vx * Fy[i] - vy * Fx[i])
                vx = (vx * b1 + fx + Fx[i] * b3) / b2
                vy = (vy * b1 + fy + Fy[i] * b3) / b2
                vz = (vz * b1 + fz + Fz[i] * b3) / b2
                x[i] += vx
                y[i] += vy
                z[i] += vz
                px[i] = vx * gamma
                py[i] = vy * gamma
                pz[i] = vz * gamma
            else:
                gamma = (1 + px[i] ** 2 + py[i] ** 2 + pz[i] ** 2) ** (1 / 2)
                vz = pz[i] / gamma
                z[i] += vz

    @staticmethod
    @cuda.jit
    def _second_step_red_cuda(
        x: np.array,
        y: np.array,
        z: np.array,
        px: np.array,
        py: np.array,
        pz: np.array,
        Fx: np.array,
        Fy: np.array,
        Fz: np.array,
        z_start: float,
        z_stop: float,
    ) -> None:
        i = cuda.grid(1)  # pylint: disable=no-value-for-parameter
        if z_start <= z[i] <= z_stop and i < x.size:
            gamma = (1 + px[i] ** 2 + py[i] ** 2 + pz[i] ** 2) ** (1 / 2)
            vx = px[i] / gamma
            vy = py[i] / gamma
            vz = pz[i] / gamma
            b2 = 1 + Fx[i] ** 2 + Fy[i] ** 2 + Fz[i] ** 2
            b1 = 2 - b2
            b3 = 2 * (vx * Fx[i] + vy * Fy[i] + vz * Fz[i])
            fx = 2 * (vy * Fz[i] - vz * Fy[i])
            fy = 2 * (vz * Fx[i] - vx * Fz[i])
            fz = 2 * (vx * Fy[i] - vy * Fx[i])
            vx = (vx * b1 + fx + Fx[i] * b3) / b2
            vy = (vy * b1 + fy + Fy[i] * b3) / b2
            vz = (vz * b1 + fz + Fz[i] * b3) / b2
            x[i] += vx
            y[i] += vy
            z[i] += vz
            px[i] = vx * gamma
            py[i] = vy * gamma
            pz[i] = vz * gamma
        else:
            gamma = (1 + px[i] ** 2 + py[i] ** 2 + pz[i] ** 2) ** (1 / 2)
            vz = pz[i] / gamma
            z[i] += vz

    def _track(self, *, n_files: int = 20) -> None:
        # Init parameterss
        Y = self.beam.da
        Y[2] = Y[2] + self.acc.z_start - max(Y[2])  # set initial beam position

        z_start = self.acc.z_start
        z_stop = self.acc.z_stop
        dz = self.acc.dz
        dt = dz / const.c
        t_max = (z_stop - z_start) / const.c

        m = self.beam.type.mass
        q = self.beam.type.charge

        P0 = m * const.c * const.c / (const.e * 1e6)
        E0 = m * const.c / (q * dt * 1e6)
        B0 = m / (q * dt)

        # RED
        for t in np.arange(0, t_max, 2 * dt):
            # into the internal system
            Y[0], Y[1], Y[2] = Y[0] / dz, Y[1] / dz, Y[2] / dz
            Y[3], Y[4], Y[5] = Y[3] / P0, Y[4] / P0, Y[5] / P0

            # get electric field from accelerator
            Ex, Ey, Ez = get_field_accelerator(self.acc, "E", Y[0] * dz, Y[1] * dz, Y[2] * dz)
            Ex, Ey, Ez = Ex / E0, Ey / E0, Ez / E0

            # get electric field from beam
            if self.beam.total_charge:
                ex, ey, ez = get_field_beam(self.beam, self.acc, "E", Y[0] * dz, Y[1] * dz, Y[2] * dz)
                ex, ey, ez = ex / E0, ey / E0, ez / E0
                Ex, Ey, Ez = Ex + ex, Ey + ey, Ez + ez

            # first step RED
            if cfg.ENABLE_CUDA:
                self._first_step_red_cuda[self.beam.n // cfg.CUDA_THREADS_PER_BLOCK + 1, cfg.CUDA_THREADS_PER_BLOCK](
                    Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Ex, Ey, Ez, z_start / dz, z_stop / dz
                )
            else:
                self._first_step_red(Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Ex, Ey, Ez, z_start / dz, z_stop / dz)
            gamma = np.sqrt(1 + Y[3] * Y[3] + Y[4] * Y[4] + Y[5] * Y[5])
            vz = Y[5] / gamma

            # get magnetic field from accelerator
            Bx, By, Bz = get_field_accelerator(self.acc, "B", Y[0] * dz, Y[1] * dz, Y[2] * dz)
            Bx, By, Bz = Bx / B0 / gamma, By / B0 / gamma, Bz / B0 / gamma

            # get magnetic field from beam
            if self.beam.total_charge:
                bx, by, bz = get_field_beam(self.beam, self.acc, "B", Y[0] * dz, Y[1] * dz, Y[2] * dz)
                bx, by, bz = bx * vz / B0, by * vz / B0, bz * vz / B0
                Bx, By, Bz = Bx + bx, By + by, Bz + bz

            # second step RED
            if cfg.ENABLE_CUDA:
                self._second_step_red_cuda[self.beam.n // cfg.CUDA_THREADS_PER_BLOCK + 1, cfg.CUDA_THREADS_PER_BLOCK](
                    Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Bx, By, Bz, z_start / dz, z_stop / dz
                )
            else:
                self._second_step_red(Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Bx, By, Bz, z_start / dz, z_stop / dz)
            gamma = np.sqrt(1 + Y[3] * Y[3] + Y[4] * Y[4] + Y[5] * Y[5])
            vz = Y[5] / gamma

            # into the SI
            Y[0], Y[1], Y[2] = Y[0] * dz, Y[1] * dz, Y[2] * dz
            Y[3], Y[4], Y[5] = Y[3] * P0, Y[4] * P0, Y[5] * P0
            Bx, By, Bz = Bx * gamma * B0, By * gamma * B0, Bz * gamma * B0
            Ex, Ey, Ez = Ex * E0, Ey * E0, Ez * E0

            progress = t / t_max * 100
            meters = z_start + t * const.c
            X = np.row_stack((Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Bx, By, Bz, Ex, Ey, Ez))
            Xt = np.transpose(X)
            df = pd.DataFrame(Xt, columns=["x", "y", "z", "px", "py", "pz", "Bx", "By", "Bz", "Ex", "Ey", "Ez"])
            if progress % (100 // n_files) < 2 * dt / t_max * 100:
                self.result.update({round(meters, 3): df})
            print("\rz = %.2f m (%.1f %%) " % (meters, progress), end="")
