import numpy as np
from numba import cuda, prange
from scipy import misc

from redpic import constants as const
from redpic.accelerator import Accelerator
from redpic.beam.base import BaseBeam
from redpic.core import config as cfg
from redpic.utils.jit import jit


def get_field_accelerator(
    acc: Accelerator, type: str, x: np.array, y: np.array, z: np.array, t: np.float64
) -> (np.array, np.array, np.array):
    dz = acc.dz
    offset_x = acc.Dx(z)
    offset_xp = acc.Dxp(z)
    offset_y = acc.Dy(z)
    offset_yp = acc.Dyp(z)
    x_corr = x - offset_x
    y_corr = y - offset_y
    r_corr = np.sqrt(x_corr * x_corr + y_corr * y_corr)
    if type == "E":
        Ez = acc.Ez(z)
        dEzdz = acc.dEzdz(z)
        d2Ezdz2 = misc.derivative(acc.dEzdz, z, dx=dz, n=1)
        E_freq, E_phase = acc.E_freq(z) * 1e9, acc.E_phase(z)
        Ez = (Ez - (d2Ezdz2 + Ez * (2 * np.pi * E_freq) ** 2 / const.c**2) * r_corr * r_corr / 4) * np.sin(
            2 * np.pi * E_freq * t + E_phase
        )
        Ex = (-dEzdz * x_corr / 2 + Ez * offset_xp) * np.sin(2 * np.pi * E_freq * t + E_phase)
        Ey = (-dEzdz * y_corr / 2 + Ez * offset_yp) * np.sin(2 * np.pi * E_freq * t + E_phase)
        return Ex, Ey, Ez
    if type == "B":
        Bx = acc.Bx(z)
        By = acc.By(z)
        Bz = acc.Bz(z)
        dBzdz = acc.dBzdz(z)
        d2Bzdz2 = misc.derivative(acc.dBzdz, z, dx=dz, n=1)
        B_freq, B_phase = acc.B_freq(z) * 1e9, acc.B_phase(z)
        Gz = acc.Gz(z)
        Bz = (Bz - (d2Bzdz2 + Bz * (2 * np.pi * B_freq) ** 2 / const.c**2) * r_corr * r_corr / 4) * np.cos(
            2 * np.pi * B_freq * t + B_phase
        )
        Bx = (Bx - dBzdz * x_corr / 2 + Bz * offset_xp) * np.cos(2 * np.pi * B_freq * t + B_phase) + Gz * y_corr
        By = (By - dBzdz * y_corr / 2 + Bz * offset_yp) * np.cos(2 * np.pi * B_freq * t + B_phase) + Gz * x_corr
        return Bx, By, Bz


@jit
def sum_field_particles(
    x: np.array, y: np.array, z: np.array, z_start: float, z_stop: float, Fx: np.array, Fy: np.array, Fz: np.array
) -> None:
    for i in prange(int(len(x))):  # pylint: disable=E1133
        if z_start <= z[i] <= z_stop:
            for j in range(int(len(x))):
                if z_start <= z[j] <= z_stop and i != j:
                    Fx[i] += (x[i] - x[j]) / ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2) ** (3 / 2)
                    Fy[i] += (y[i] - y[j]) / ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2) ** (3 / 2)
                    Fz[i] += (z[i] - z[j]) / ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2) ** (3 / 2)


@cuda.jit
def sum_field_particles_cuda(
    x: np.array, y: np.array, z: np.array, z_start: float, z_stop: float, Fx: np.array, Fy: np.array, Fz: np.array
) -> None:
    i = cuda.grid(1)  # pylint: disable=no-value-for-parameter

    if i >= x.size:
        return

    for j in range(x.size):
        if z_start <= z[i] <= z_stop and z_start <= z[j] <= z_stop and i != j:
            Fx[i] += (x[i] - x[j]) / ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2) ** (3 / 2)
            Fy[i] += (y[i] - y[j]) / ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2) ** (3 / 2)
            Fz[i] += (z[i] - z[j]) / ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2) ** (3 / 2)


def get_field_beam(
    beam: BaseBeam, acc: Accelerator, type: str, x: np.array, y: np.array, z: np.array
) -> (np.array, np.array, np.array):
    q = beam.total_charge / beam.n
    Fx, Fy, Fz = np.zeros(beam.n), np.zeros(beam.n), np.zeros(beam.n)

    if cfg.ENABLE_CUDA:
        sum_field_particles_cuda[beam.n // cfg.CUDA_THREADS_PER_BLOCK + 1, cfg.CUDA_THREADS_PER_BLOCK](
            x, y, z, acc.z_start, acc.z_stop, Fx, Fy, Fz
        )
    else:
        sum_field_particles(x, y, z, acc.z_start, acc.z_stop, Fx, Fy, Fz)

    if type == "E":
        k = const.ke / (4 * np.pi) * q
        return k * Fx * 1e-6, k * Fy * 1e-6, k * Fz * 1e-6
    if type == "B":
        k = const.km / (4 * np.pi) * q
        return k * Fx, -k * Fy, 0 * k * Fz
