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

    def _track(self, *, n_files: int = cfg.DEFAULT_TRACK_SAVE_N_FILES) -> None:
        # Init parameters
        z_start, z_stop, dz = self.acc.z_start, self.acc.z_stop, self.acc.dz
        t_start, t_stop, dt = self.t_start, self.t_stop, self.dt

        Y = np.transpose(np.concatenate([beam.df.to_numpy() for beam in self.beams]))
        Y[2] = Y[2] + z_start - max(Y[2])  # set initial beam position

        m = np.transpose(np.concatenate([np.full(beam.n, beam.type.mass) for beam in self.beams]))
        phys_q = np.transpose(np.concatenate([np.full(beam.n, beam.type.charge) for beam in self.beams]))
        macro_q = np.transpose(np.concatenate([np.full(beam.n, beam.total_charge / beam.n) for beam in self.beams]))

        P0 = m * const.c * const.c / (const.e * const.mega)
        E0 = m * const.c / (phys_q * dt * const.mega)
        B0 = m / (phys_q * dt)

        # RED
        for t in np.arange(t_start, t_stop, 2 * dt):
            # into the internal system
            Y[0], Y[1], Y[2] = Y[0] / dz, Y[1] / dz, Y[2] / dz
            Y[3], Y[4], Y[5] = Y[3] / P0, Y[4] / P0, Y[5] / P0

            # get electric field from accelerator
            Ex, Ey, Ez = get_field_accelerator(self.acc, "E", Y[0] * dz, Y[1] * dz, Y[2] * dz)
            Ex, Ey, Ez = Ex / E0, Ey / E0, Ez / E0

            # get electric field from beam
            if bool(beam.total_charge for beam in self.beams):
                ex, ey, ez = get_field_beam(macro_q, self.acc, "E", Y[0] * dz, Y[1] * dz, Y[2] * dz)
                ex, ey, ez = ex / E0, ey / E0, ez / E0
                Ex, Ey, Ez = Ex + ex, Ey + ey, Ez + ez

            # first step RED
            if cfg.ENABLE_CUDA:
                self._first_step_red_cuda[
                    min(beam.n for beam in self.beams) // cfg.CUDA_THREADS_PER_BLOCK + 1, cfg.CUDA_THREADS_PER_BLOCK
                ](Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Ex, Ey, Ez, z_start / dz, z_stop / dz)
            else:
                self._first_step_red(Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Ex, Ey, Ez, z_start / dz, z_stop / dz)
            gamma = np.sqrt(1 + Y[3] * Y[3] + Y[4] * Y[4] + Y[5] * Y[5])
            vz = Y[5] / gamma

            # get magnetic field from accelerator
            Bx, By, Bz = get_field_accelerator(self.acc, "B", Y[0] * dz, Y[1] * dz, Y[2] * dz)
            Bx, By, Bz = Bx / B0 / gamma, By / B0 / gamma, Bz / B0 / gamma

            # get magnetic field from beam
            if bool(beam.total_charge for beam in self.beams):
                bx, by, bz = get_field_beam(macro_q, self.acc, "B", Y[0] * dz, Y[1] * dz, Y[2] * dz)
                bx, by, bz = bx * vz / B0, by * vz / B0, bz * vz / B0
                Bx, By, Bz = Bx + bx, By + by, Bz + bz

            # second step RED
            if cfg.ENABLE_CUDA:
                self._second_step_red_cuda[
                    min(beam.n for beam in self.beams) // cfg.CUDA_THREADS_PER_BLOCK + 1, cfg.CUDA_THREADS_PER_BLOCK
                ](Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Bx, By, Bz, z_start / dz, z_stop / dz)
            else:
                self._second_step_red(Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Bx, By, Bz, z_start / dz, z_stop / dz)
            gamma = np.sqrt(1 + Y[3] * Y[3] + Y[4] * Y[4] + Y[5] * Y[5])
            vz = Y[5] / gamma

            # into the SI
            Y[0], Y[1], Y[2] = Y[0] * dz, Y[1] * dz, Y[2] * dz
            Y[3], Y[4], Y[5] = Y[3] * P0, Y[4] * P0, Y[5] * P0
            Bx, By, Bz = Bx * gamma * B0, By * gamma * B0, Bz * gamma * B0
            Ex, Ey, Ez = Ex * E0, Ey * E0, Ez * E0

            progress, meters = t / t_stop * 100, z_start + t * const.c
            if progress % (100 // n_files) < 2 * dt / t_stop * 100:
                self.result.update(
                    {
                        round(meters, 3): pd.DataFrame(
                            np.transpose(
                                np.vstack((Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Bx, By, Bz, Ex, Ey, Ez, macro_q))
                            ),
                            columns=["x", "y", "z", "px", "py", "pz", "Bx", "By", "Bz", "Ex", "Ey", "Ez", "charge"],
                        )
                    }
                )
            print("\rz = %.2f m (%.1f %%) " % (meters, progress), end="")
