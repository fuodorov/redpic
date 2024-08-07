import logging
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd

from redpic import constants as const

module_logger = logging.getLogger(__name__)


class BeamDistributionType:
    NONE = "none"
    UNIFORM = "uniform"
    RADIAL_UNIFORM = "radial-uniform"
    GAUSSIAN = "gaussian"
    ASTRA = "astra"


class BaseBeam(ABC):
    distribution: BeamDistributionType = BeamDistributionType.NONE

    def __init__(
        self,
        *,
        type: const.Element = const.electron,
        current: float = 0.0e0,
        energy: float = 0.0e0,
        radius: float = 0.0e0,
        radius_x: float = 0.0e0,
        radius_y: float = 0.0e0,
        radius_z: float = 0.0e0,
        rp: float = 0.0e0,
        radius_xp: float = 0.0e0,
        radius_yp: float = 0.0e0,
        normalized_emittance: float = 0.0e0,
        normalized_emittance_x: float = 0.0e0,
        normalized_emittance_y: float = 0.0e0,
        x: float = 0.0e0,
        y: float = 0.0e0,
        z: float = 0.0e0,
        xp: float = 0.0e0,
        yp: float = 0.0e0,
    ):
        """
        Initialization of an electron beam.

        Parameters
        ----------
        current: float
            Beam current, [A]
        energy: float
            Beam energy, [MeV]
        radius: float
            Beam radius, [m]
        radius_x: float, optional
            Elliptical beam radius along the x-axis, [m]
        radius_y: float, optional
            Elliptical beam radius along the y-axis, [m]
        rp: float
            Beam radius prime, [rad]
        radius_xp: float, optional
            Elliptical beam radius prime along the x-axis, [rad]
        radius_yp: float, optional
            Elliptical beam radius prime along the y-axis, [rad]
        normalized_emittance: float
            Beam normalized emittance, [m*rad]
        normalized_emittance_x: float, optional
            Elliptical beam normalized emittance along the x-axis, [m*rad]
        normalized_emittance_y: float, optional
            Elliptical beam normalized emittance along the y-axis, [m*rad]

        x: float, optional
            Offset of the centroid along the x-axis, [m]
        xp: float, optional
            Centroid rotation in the z-x plane, [rad]
        y: float, optional
            Offset of the centroid along the y-axis, [m]
        yp: float, optional
            Centroid rotation in the z-y plane, [rad]
        """
        self.current = current
        self.energy = energy
        self.radius = radius
        self.rp = rp
        self.radius_x = radius_x
        self.radius_y = radius_y
        self.radius_xp = radius_xp
        self.radius_yp = radius_yp
        self.normalized_emittance = normalized_emittance
        self.normalized_emittance_x = normalized_emittance_x
        self.normalized_emittance_y = normalized_emittance_y
        if radius != 0.0e0:
            self.radius_x = radius
            self.radius_y = radius
        if rp != 0.0e0:
            self.radius_xp = rp
            self.radius_yp = rp
        if normalized_emittance != 0.0e0:
            self.normalized_emittance_x = normalized_emittance
            self.normalized_emittance_y = normalized_emittance

        self.x = x
        self.y = y
        self.xp = xp
        self.yp = yp

        self.gamma = gamma = self.energy / const.electron_mass_energy + 1
        self.beta = beta = np.sqrt(1 - 1 / (gamma * gamma))

        self.p = self.momentum = gamma * beta * const.electron_mass_energy
        self.px = self.p * self.radius_xp
        self.py = self.p * self.radius_yp
        self.pz = self.p
        self.description = ""

        self.type = type  # particles type
        self.n = 0.0  # quantity
        self.df = pd.DataFrame  # data frame
        self.radius_z = radius_z
        self.z = z
        self.total_charge = (
            (2 * type.charge / abs(type.charge) * self.current * self.radius_z / const.c / self.beta)
            if self.beta != 0.0e0
            else 0.0e0
        )

    @abstractmethod
    def _generate(self, n: int = 0, *, file_name: str = "", **kwargs) -> None:
        pass

    def generate(self, n: int = 0, *, file_name: str = "", **kwargs) -> None:
        if n < 0:
            raise ValueError("The number of particles (n) must be a positive number!")
        module_logger.info("Generate a %s beam with %d particles", self.distribution, n)
        self._generate(n, file_name=file_name, **kwargs)

    def __str__(self):
        return f"""Beam parameters:
            Type\t{self.type.name}
            Distribution\t{self.distribution}
            Particles\t{self.n}
            Current\t{self.current} A
            Energy\t{self.energy} MeV
            Total momentum\t{self.momentum} MeV/c
            Rel. factor\t{self.gamma}
            Radius x\t{self.radius_x / const.milli} mm
            Radius y\t{self.radius_y / const.milli} mm
            Radius z\t{self.radius_z} m
            Radius x prime\t{self.radius_xp / const.milli} mrad
            Radius y prime\t{self.radius_yp / const.milli} mrad
            Horizontal centroid position\t{self.x / const.milli} mm
            Vertical centroid position\t{self.y / const.milli} mm
            Horizontal centroid angle\t{self.xp / const.milli} mrad
            Vertical centroid angle\t{self.yp / const.milli} mrad
            Normalized emittance x\t{self.normalized_emittance_x / const.milli / const.milli} mm*mrad
            Normalized emittance y\t{self.normalized_emittance_y / const.milli / const.milli} mm*mrad
        """
