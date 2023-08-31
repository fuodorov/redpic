import numpy as np
from numba import prange
from scipy import misc

from redpic import constants as const
from redpic.accelerator import Accelerator
from redpic.beam.base import BaseBeam
from redpic.utils.jit import jit


def get_field_accelerator(
    acc: Accelerator, type: str, x: np.array, y: np.array, z: np.array
) -> (np.array, np.array, np.array):
    assert type in ("E", "B"), "Check type field! (E or B)"

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
        Ez = Ez - d2Ezdz2 * r_corr * r_corr / 4  # row remainder
        Ex = -dEzdz * x_corr / 2 + Ez * offset_xp  # row remainder
        Ey = -dEzdz * y_corr / 2 + Ez * offset_yp  # row remainder
        return Ex, Ey, Ez
    if type == "B":
        Bx = acc.Bx(z)
        By = acc.By(z)
        Bz = acc.Bz(z)
        dBzdz = acc.dBzdz(z)
        d2Bzdz2 = misc.derivative(acc.dBzdz, z, dx=dz, n=1)
        Gz = acc.Gz(z)
        Bz = Bz - d2Bzdz2 * r_corr * r_corr / 4  # row remainder
        Bx = Bx + Gz * y_corr - dBzdz * x_corr / 2 + Bz * offset_xp  # row remainder
        By = By + Gz * x_corr - dBzdz * y_corr / 2 + Bz * offset_yp  # row remainder
        return Bx, By, Bz


@jit
def sum_field_particles(
    x: np.array, y: np.array, z: np.array, z_start: float, z_stop: float
) -> (np.array, np.array, np.array):
    n = int(len(x))
    Fx = np.zeros(n)
    Fy = np.zeros(n)
    Fz = np.zeros(n)
    for i in prange(n):  # pylint: disable=E1133
        if z_start <= z[i] <= z_stop:
            for j in prange(n):  # pylint: disable=E1133
                if z_start <= z[j] <= z_stop and i != j:
                    Fx[i] += (x[i] - x[j]) / ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2) ** (3 / 2)
                    Fy[i] += (y[i] - y[j]) / ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2) ** (3 / 2)
                    Fz[i] += (z[i] - z[j]) / ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2) ** (3 / 2)
    return Fx, Fy, Fz


def get_field_beam(
    beam: BaseBeam, acc: Accelerator, type: str, x: np.array, y: np.array, z: np.array
) -> (np.array, np.array, np.array):
    assert type in ("E", "B"), "Check type field! (E or B)"
    Q = beam.total_charge
    q = Q / beam.n
    if type == "E":
        Ex, Ey, Ez = sum_field_particles(x, y, z, acc.z_start, acc.z_stop)
        return (
            const.ke / (4 * np.pi) * q * Ex / 1e6,
            const.ke / (4 * np.pi) * q * Ey / 1e6,
            const.ke / (4 * np.pi) * q * Ez / 1e6,
        )
    if type == "B":
        Bx, By, Bz = sum_field_particles(x, y, z, acc.z_start, acc.z_stop)
        return const.km / (4 * np.pi) * q * Bx, -const.km / (4 * np.pi) * q * By, 0 * const.km / (4 * np.pi) * q * Bz
