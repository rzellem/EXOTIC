"""Scalar proposal geometry for UltraNest's one-parameter fits."""

import numpy as np
from ultranest.mlfriends import RobustEllipsoidRegion


class OneDimensionalRegion(RobustEllipsoidRegion):
    """Bootstrap an enclosing interval without scalar covariance inversion.

    An interval is the one-dimensional ellipsoid. Bootstrap rounds compare
    held-out distances with the selected interval's radius. A round with no
    held-out points, or zero held-out distance, needs no additional expansion.
    The resulting region always encloses every live point and respects the
    sampler's remaining-volume estimate; it changes only proposal geometry.
    """

    @staticmethod
    def _interval(values, minvol):
        center = float(np.mean(values))
        radius = max(float(np.max(np.abs(values - center))), 0.5 * minvol)
        if radius == 0.0:
            # Coincident live points give no interval width. Keep the unit
            # prior cube accessible when no remaining volume was supplied.
            radius = max(center, 1.0 - center)
        return center, radius

    def compute_enlargement(self, nbootstraps=50, minvol=0.0, rng=np.random):
        values = self.u[:, 0]
        enlargement = 1.0
        selected = np.empty(len(values), dtype=bool)
        for _ in range(nbootstraps):
            selected[:] = False
            selected[rng.randint(len(values), size=len(values))] = True
            if selected.all():
                continue
            center, radius = self._interval(values[selected], minvol)
            distance = float(np.max(np.abs(values[~selected] - center)))
            enlargement = max(enlargement, (distance / radius) ** 2)
        # RobustEllipsoidRegion uses only the enclosing ellipsoid, so its
        # MLFriends neighbourhood radius is deliberately unbounded upstream.
        return 1e300, enlargement

    def create_ellipsoid(self, minvol=0.0):
        center, radius = self._interval(self.u[:, 0], minvol)
        self.ellipsoid_center = np.array([center])
        self.ellipsoid_cov = np.array([[radius ** 2]])
        self.ellipsoid_invcov = np.array([[1.0 / radius ** 2]])
        self.ellipsoid_axlens = np.array([radius])
        self.ellipsoid_axes = np.array([[radius]])
        self.ellipsoid_axes_T = self.ellipsoid_axes.T
        self.ellipsoid_inv_axlens = np.array([1.0 / radius])
        self.ellipsoid_inv_axes = np.array([[1.0 / radius]])
