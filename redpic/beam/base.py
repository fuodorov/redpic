import logging
from abc import ABC, abstractmethod

import kenv as kv
import numpy as np
import pandas as pd

from redpic import constants as const

module_logger = logging.getLogger(__name__)


class BaseBeam(ABC, kv.Beam):
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
        self.distribution = ""
        self.n = 0.0  # quantity
        self.df = pd.DataFrame  # data frame
        self.da = np.array
        self.radius_z = radius_z
        self.z = z
        self.total_charge = (
            2 * type.charge / abs(type.charge) * self.current * self.radius_z / const.c / self.beta
        )  # beam charge

    @abstractmethod
    def generate(self, n: int, **kwargs) -> None:
        pass

    def _read_distribution_file(self, fname):
        module_logger.info("Read distribution from file")
        df = pd.read_csv(fname, dtype="float32")
        return df["x"], df["y"], df["z"], df["px"], df["py"], df["pz"]

    def _read_distribution_astra(self, fname):
        module_logger.info("Read distribution from Astra file")

        cols = ["x", "y", "z", "px", "py", "pz", "clock", "charge", "id", "flag"]
        #        m    m    m    eV/c  eV/c  eV/c  ns       nC
        df = pd.read_csv(fname, header=None, delim_whitespace=True, names=cols, dtype="float32")

        df = df[df.flag != -15]  # ignore the lost particles
        df["px"] = df["px"] / 1e6  # MeV/c
        df["py"] = df["py"] / 1e6  # MeV/c

        # remove the reference particle
        df0 = df.head(1)
        df = df.drop(df0.index)

        z0 = df0.z.values[0]
        pz0 = df0.pz.values[0]

        # Recalculate z and pz:
        df["z"] = z0 + df["clock"] * 1e-9 * const.c  # m
        df["pz"] = (pz0 + df["pz"]) / 1e6  # MeV/c
        return df["x"], df["y"], df["z"], df["px"], df["py"], df["pz"]

    def upload(self, file_name: str):
        module_logger.info("Particle loading")
        file_extension = file_name.split(".")[-1]

        if file_extension == "ini":
            x, y, z, px, py, pz = self._read_distribution_astra(file_name)
            self.n = int(len(x))
            self.da = np.row_stack((x, y, z, px, py, pz))
            self.df = pd.DataFrame(np.transpose(self.da), columns=["x", "y", "z", "px", "py", "pz"])
        if file_extension == "csv":
            x, y, z, px, py, pz = self._read_distribution_file(file_name)
            self.n = int(len(x))
            self.da = np.row_stack((x, y, z, px, py, pz))
            self.df = pd.DataFrame(np.transpose(self.da), columns=["x", "y", "z", "px", "py", "pz"])

    def __str__(self):
        self.description = (
            "Beam parameters:"
            + "\n"
            + "\tType\t"
            + (self.type.name)
            + "\n"
            + "\tDistribution\t"
            + (self.distribution)
            + "\n"
            + "\tParticles\t%0.0f" % (self.n)
            + "\n"
            + "\tCurrent\t%0.0f A" % (self.current)
            + "\n"
            + "\tEnergy\t%0.3f MeV" % (self.energy)
            + "\n"
            + "\tTotal momentum\t%0.3f MeV/c" % (self.momentum)
            + "\n"
            + "\tRel. factor\t%0.3f" % (self.gamma)
            + "\n"
            + "\tRadius x\t%0.1f mm" % (self.radius_x * 1e3)
            + "\n"
            + "\tRadius y\t%0.1f mm" % (self.radius_y * 1e3)
            + "\n"
            + "\tRadius z\t%0.1f m" % (self.radius_z)
            + "\n"
            + "\tRadius x prime\t%0.1f mrad" % (self.radius_xp * 1e3)
            + "\n"
            + "\tRadius y prime\t%0.1f mrad" % (self.radius_yp * 1e3)
            + "\n"
            + "\tHorizontal centroid position\t%0.1f mm" % (self.x * 1e3)
            + "\n"
            + "\tVertical centroid position\t%0.1f mm" % (self.y * 1e3)
            + "\n"
            + "\tHorizontal centroid angle\t%0.1f mrad" % (self.xp * 1e3)
            + "\n"
            + "\tVertical centroid angle\t%0.1f mrad" % (self.yp * 1e3)
            + "\n"
            + "\tNormalized emittance x\t%0.1f mm*mrad" % (self.normalized_emittance_x * 1e6)
            + "\n"
            + "\tNormalized emittance y\t%0.1f mm*mrad" % (self.normalized_emittance_y * 1e6)
            + "\n"
        )
        return self.description
