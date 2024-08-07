import logging

import numpy as np
import pandas as pd

from redpic.beam.base import BaseBeam, BeamDistributionType

module_logger = logging.getLogger(__name__)


@staticmethod
def get_4_sphere_normal(n):
    s = np.random.normal(0, 1, int(n))
    t = np.random.normal(0, 1, int(n))
    u = np.random.normal(0, 1, int(n))  # pylint: disable=W0621
    v = np.random.normal(0, 1, int(n))
    norm = (s * s + t * t + u * u + v * v) ** 0.5
    return (s, t, u, v) / norm


class UniformBeam(BaseBeam):
    distribution = BeamDistributionType.UNIFORM

    def _generate(self, n: int = 0, *, file_name: str = "", **kwargs) -> None:
        module_logger.info("Generate a %s beam with %d particles", self.distribution, n)
        dpz = kwargs.get("dpz", 0.01)
        self.n = n
        (s, t, u, v) = get_4_sphere_normal(self.n)
        (x, y) = (self.radius_x * s + self.x, self.radius_y * t + self.y)
        (px, py) = (self.px * u + self.p * self.xp, v * self.py + self.p * self.yp)

        z = np.random.uniform(-self.radius_z + self.z, self.radius_z + self.z, int(self.n))
        pz = np.random.uniform((1 - dpz) * self.pz, (1 + dpz) * self.pz, int(self.n))

        self.df = pd.DataFrame(
            np.transpose(np.vstack((x, y, z, px, py, pz))), columns=["x", "y", "z", "px", "py", "pz"]
        )
