import logging

import numpy as np
import pandas as pd

from redpic.beam.base import BaseBeam, BeamDistributionType
from redpic.beam.uniform import get_4_sphere_normal

module_logger = logging.getLogger(__name__)


class RadialUniformBeam(BaseBeam):
    distribution = BeamDistributionType.RADIAL_UNIFORM

    def _generate(self, n: int = 0, *, file_name: str = "", **kwargs) -> None:
        dpz = kwargs.get("dpz", 0.01)
        self.n = n
        (s, t, u, v) = get_4_sphere_normal(self.n)
        (x, y) = (self.radius_x * s, self.radius_y * t)
        u = 2 * (self.px * self.normalized_emittance_x * u + self.radius_xp * x) / self.radius_x
        v = 2 * (self.py * self.normalized_emittance_y * v + self.radius_yp * y) / self.radius_y
        (x, y) = (x + self.x, y + self.y)
        (px, py) = (u + self.p * self.xp, v + self.p * self.yp)

        z = np.random.uniform(-self.radius_z + self.z, self.radius_z + self.z, int(self.n))
        pz = np.random.uniform((1 - dpz) * self.pz, (1 + dpz) * self.pz, int(self.n))

        self.da = np.row_stack((x, y, z, px, py, pz))
        self.df = pd.DataFrame(np.transpose(self.da), columns=["x", "y", "z", "px", "py", "pz"])
