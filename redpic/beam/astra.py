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
        df = pd.read_csv(file_name, header=None, delim_whitespace=True, names=cols, dtype="float32")
        df = df[df.flag != -15]  # ignore the lost particles
        df["px"] = df["px"] / 1e6  # MeV/c
        df["py"] = df["py"] / 1e6  # MeV/c
        df0 = df.head(1)  # remove the reference particle
        df = df.drop(df0.index)
        z0 = df0.z.values[0]
        pz0 = df0.pz.values[0]
        df["z"] = z0 + df["clock"] * 1e-9 * const.c  # m
        df["pz"] = (pz0 + df["pz"]) / 1e6  # MeV/c
        return df["x"], df["y"], df["z"], df["px"], df["py"], df["pz"]

    def generate(self, n: int = 0, *, file_name: str = "", **kwargs) -> None:
        super().generate(n, file_name=file_name, **kwargs)
        if file_name.split(".")[-1] != "ini":
            raise ValueError("The file format is not supported! Astra file should be .ini")

        module_logger.info("Particle loading from Astra %s file", file_name)
        x, y, z, px, py, pz = self._read_astra_particles(file_name)
        self.n = int(len(x))
        self.da = np.row_stack((x, y, z, px, py, pz))
        self.df = pd.DataFrame(np.transpose(self.da), columns=["x", "y", "z", "px", "py", "pz"])
