import logging

import numpy as np
import pandas as pd

from redpic import constants as const
from redpic.beam.base import BaseBeam, BeamDistributionType

module_logger = logging.getLogger(__name__)


class AstraBeam(BaseBeam):
    distribution = BeamDistributionType.ASTRA

    def _read_astra_particles(self, file_name: str) -> tuple:
        module_logger.info("Read particles from Astra file")
        cols = [
            "x",
            "y",
            "z",
            "px",
            "py",
            "pz",
            "clock",
            "charge",
            "id",
            "flag",
        ]  # m    m    m    eV/c  eV/c  eV/c  ns       nC
        df = pd.read_csv(file_name, header=None, sep=r"\s+", names=cols, dtype="float32")
        df = df[df.flag != -15]  # ignore the lost particles
        df["px"] = df["px"] / const.mega  # MeV/c
        df["py"] = df["py"] / const.mega  # MeV/c
        df0 = df.head(1)  # remove the reference particle
        df = df.drop(df0.index)
        z0 = df0.z.values[0]
        pz0 = df0.pz.values[0]
        df["z"] = z0 + df["clock"] * 1e-9 * const.c  # m
        df["pz"] = (pz0 + df["pz"]) / const.mega  # MeV/c
        return df["x"], df["y"], df["z"], df["px"], df["py"], df["pz"]

    def _generate(self, n: int = 0, *, file_name: str = "", **kwargs) -> None:
        module_logger.info("Particle loading from Astra %s file", file_name)
        x, y, z, px, py, pz = self._read_astra_particles(file_name)
        self.n = int(len(x))
        self.df = pd.DataFrame(
            np.transpose(np.vstack((x, y, z, px, py, pz))), columns=["x", "y", "z", "px", "py", "pz"]
        )
