import logging
from abc import ABC, abstractmethod

import kenv as kv
import numpy as np
import pandas as pd

from redpic import constants as const

module_logger = logging.getLogger(__name__)


class BeamDistributionType:
    NONE = "none"
    UNIFORM = "uniform"
    GAUSSIAN = "gaussian"
    ASTRA = "astra"


class BaseBeam(ABC, kv.Beam):
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
        super().__init__(
            current=current,
            energy=energy,
            radius=radius,
            radius_x=radius_x,
            radius_y=radius_y,
            rp=rp,
            radius_xp=radius_xp,
            radius_yp=radius_yp,
            normalized_emittance=normalized_emittance,
            normalized_emittance_x=normalized_emittance_x,
            normalized_emittance_y=normalized_emittance_y,
            x=x,
            y=y,
            xp=xp,
            yp=yp,
        )
        self.type = type  # particles type
        self.n = 0.0  # quantity
        self.df = pd.DataFrame  # data frame
        self.da = np.array
        self.radius_z = radius_z
        self.z = z
        try:
            self.total_charge = (
                2 * type.charge / abs(type.charge) * self.current * self.radius_z / const.c / self.beta
            )  # beam charge
        except ZeroDivisionError:
            self.total_charge = 0.0e0

    @abstractmethod
    def generate(self, n: int = 0, *, file_name: str = "", **kwargs) -> None:
        if n < 0:
            raise ValueError("The number of particles (n) must be a positive number!")

    def __str__(self):
        return f"""Beam parameters:
            Type\t{self.type.name}
            Distribution\t{self.distribution}
            Particles\t{self.n}
            Current\t{self.current} A
            Energy\t{self.energy} MeV
            Total momentum\t{self.momentum} MeV/c
            Rel. factor\t{self.gamma}
            Radius x\t{self.radius_x * 1e3} mm
            Radius y\t{self.radius_y * 1e3} mm
            Radius z\t{self.radius_z} m
            Radius x prime\t{self.radius_xp * 1e3} mrad
            Radius y prime\t{self.radius_yp * 1e3} mrad
            Horizontal centroid position\t{self.x * 1e3} mm
            Vertical centroid position\t{self.y * 1e3} mm
            Horizontal centroid angle\t{self.xp * 1e3} mrad
            Vertical centroid angle\t{self.yp * 1e3} mrad
            Normalized emittance x\t{self.normalized_emittance_x * 1e6} mm*mrad
            Normalized emittance y\t{self.normalized_emittance_y * 1e6} mm*mrad
        """
