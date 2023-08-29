import logging

import numpy as np
import pandas as pd

from redpic.beam.base import BaseBeam, BeamDistributionType

module_logger = logging.getLogger(__name__)


class UniformBeam(BaseBeam):
    distribution = BeamDistributionType.UNIFORM

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def generate(self, n: int = 0, *, file_name: str = "", **kwargs) -> None:
        super().generate(n, file_name=file_name, **kwargs)
        module_logger.info("Generate a %s beam with %d particles", self.distribution, n)
        dpz = kwargs.get("dpz", 0.01)
        self.n = n
        s = np.random.normal(0, 1, int(self.n))
        t = np.random.normal(0, 1, int(self.n))
        u = np.random.normal(0, 1, int(self.n))  # pylint: disable=W0621
        v = np.random.normal(0, 1, int(self.n))
        norm = (s * s + t * t + u * u + v * v) ** 0.5
        (s, t, u, v) = (s, t, u, v) / norm
        (x, y) = (self.radius_x * s, self.radius_y * t)
        u = 2 * (self.px * self.normalized_emittance_x * u + self.radius_xp * x) / self.radius_x
        v = 2 * (self.py * self.normalized_emittance_y * v + self.radius_yp * y) / self.radius_y
        (x, y) = (x + self.x, y + self.y)
        (px, py) = (u + self.p * self.xp, v + self.p * self.yp)

        z = np.random.uniform(-self.radius_z + self.z, self.radius_z + self.z, int(self.n))
        pz = np.random.uniform((1 - dpz) * self.pz, (1 + dpz) * self.pz, int(self.n))

        self.da = np.row_stack((x, y, z, px, py, pz))
        self.df = pd.DataFrame(np.transpose(self.da), columns=["x", "y", "z", "px", "py", "pz"])
