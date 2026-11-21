"""limbdarkfast: transit light curves for the hybrid limb-darkening laws (NumPy).

The NumPy version, and the reference for the others. It gives the flux
of a star with hybrid2, hybrid4 or hybrid5 limb darkening, occulted by a
spherical or an oblate planet, in closed form and in float64. There are no
orbits, exposures or derivatives: give it the planet's position on the sky
and it returns the normalised flux.

If you use this code, please cite

    Kipping, D. 2027, "Limb Darkened Transits of Oblate Planets as
    Elementary Functions", AJ, submitted.

Copyright (c) 2026 David Kipping. Released under the MIT License: see
the LICENSE file, or https://opensource.org/licenses/MIT.

    F = hybrid_flux(z, r, w, law)                        spherical planet
    F = hybrid_flux_oblate(X, Y, r, f, theta, w, law)    oblate planet

Lengths are in stellar radii. ``z`` is the sky separation of the planet's
centre from the star's; ``(X, Y)`` is its sky position, X along the
direction of motion. ``r`` is the planet's radius, for an oblate planet its
area-equivalent radius sqrt(A B); ``f = 1 - B/A`` is the projected
flattening and ``theta`` (radians) the angle of the long axis from X.
``w`` are the law's shape weights and ``law`` is "hybrid2", "hybrid4" or
"hybrid5". Every argument but ``w`` and ``law`` may be an array; they
broadcast together.

The laws are written in a shape basis,

    I(mu) / I(1) = 1 - sum_j w_j S_j(mu),

each shape falling from 1 at the limb to 0 at the centre. They combine
mu^0, mu^2, mu^4 with double poles (mu^2 + eps)^-2, chosen so that the
occulted flux of every term is elementary for any elliptical planet.

Method. Green's theorem with the radial potential h(s) = (1/2s) int_0^s I,
s = x^2 + y^2, turns the occulted flux into a line integral around the
overlap of planet and star: h(1) times the limb arc inside the planet,
plus the integral of h(s) along the planet's outline inside the star.

* Spherical planet: everything reduces to the lens geometry of two
  circles. mu^0, mu^2, mu^4 come from the even recursion of Agol, Luger &
  Foreman-Mackey (2020); a pole is an arctangent or logarithm of rational
  functions of z and r.
* Oblate planet: along the ellipse the even terms are finite Fourier sums,
  and each pole is a partial-fraction sum over the roots of a quartic,
  integrated as logarithms along the arc. Near-double roots, which every
  transit meets twice per pole, are integrated as one quadratic factor.

The oblate planet must satisfy r <= (1 - f)^(3/2), so that its outline
crosses the stellar limb at most twice.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = ["hybrid_flux", "hybrid_flux_oblate", "intensity", "get_law",
           "LAWS", "HYBRID2", "HYBRID4", "HYBRID5", "ladder",
           "hybrid4_basis"]

_PI = math.pi
_TWO_PI = 2.0 * math.pi
_EPS = 2.220446049250313e-16
_TINY = 1e-280

# spherical pole column: half-width in y of the series across Q = 0
_Y_SW = 1e-3
# oblate planet
_F_SW = 1e-14          # below this flattening, the spherical forms (error ~ 200 f)
_Z_TOL = 1e-4          # |z| ~ 1: a candidate limb crossing
_G_TOL = 1e-12         # accepted crossing residual, relative
_PAIR_TOL = 1e-6       # crossings closer than this (rad) are a tangency
_PAIR_RATIO = 0.1      # near-double pole roots: |Delta| <= 0.1 d^2
_N_SERIES = 18         # terms of the pair series J
_N_ABERTH = 8
_N_BAIRSTOW = 2
_QR_SWEEPS = (6, 5)
_QR_BALANCE = 4
_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


# ===========================================================================
# the laws
# ===========================================================================

def ladder(K: int, alpha: float = 1.0) -> tuple:
    """Analytic pole ladder eps_k = exp[(2 pi / sqrt alpha)
    (sqrt(k + 1/2) - sqrt K)], k = 0..K-1."""
    return tuple(math.exp(2.0 * _PI / math.sqrt(alpha)
                          * (math.sqrt(k + 0.5) - math.sqrt(K)))
                 for k in range(K))


def hybrid4_basis(eps) -> np.ndarray:
    """hybrid4's curvature basis: rows are its shapes S1..S4 as
    combinations of the terms (1 - mu^4, (1 - mu^2)^2, Pi_eps0, Pi_eps1).
    With x = mu^2 and g = dI/dx, S3 has I(0) = 0, g(1) = g'(1) = g''(1) = 0
    and S4 has I(0) = 0, g(1) = g'(1) = 0, g''(0) = 0."""
    e0, e1 = (float(e) for e in eps)
    a0, a1 = (1.0 / (e ** -2 - (1.0 + e) ** -2) for e in (e0, e1))

    def g(x):
        return [2 * x, 2 * (1 - x), 2 * a0 * (x + e0) ** -3,
                2 * a1 * (x + e1) ** -3]

    def g1(x):
        return [2.0, -2.0, -6 * a0 * (x + e0) ** -4, -6 * a1 * (x + e1) ** -4]

    def g2(x):
        return [0.0, 0.0, 24 * a0 * (x + e0) ** -5, 24 * a1 * (x + e1) ** -5]

    rhs = [1.0, 0.0, 0.0, 0.0]
    S3 = np.linalg.solve([[1.0] * 4, g(1.0), g1(1.0), g2(1.0)], rhs)
    S4 = np.linalg.solve([[1.0] * 4, g(1.0), g1(1.0), g2(0.0)], rhs)
    return np.array([[0.5, 0.5, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], S3, S4])


class Law:
    """A hybrid law: its poles ``eps`` and the matrix ``G`` (n_w, 3 + K)
    giving each shape S_j as a combination of the generator profiles
    [mu^0, mu^2, mu^4, (mu^2 + eps_0)^-2, ...]. ``norms[j]`` is the
    unocculted flux of S_j."""

    def __init__(self, name, eps, shapes, mix=None):
        self.name = name
        self.eps = tuple(float(e) for e in eps)
        K = len(self.eps)
        n_even = len(shapes)
        self.n_w = n_even + K
        G = np.zeros((self.n_w, 3 + K))
        for j, a in enumerate(shapes):
            G[j, :3] = a
        for k, e in enumerate(self.eps):
            p = 1.0 + e
            N = e ** -2 - p ** -2          # Pi_eps = (P_eps - p^-2) / N
            G[n_even + k, 0] = -(p ** -2) / N
            G[n_even + k, 3 + k] = 1.0 / N
        if mix is not None:
            G = np.asarray(mix, dtype=np.float64) @ G
        self.G = G
        gen_norms = np.array([_PI, _PI / 2.0, _PI / 3.0]
                             + [_PI / (e * (1.0 + e)) for e in self.eps])
        self.norms = G @ gen_norms

    def __repr__(self):
        return f"Law({self.name!r}, eps={self.eps})"


HYBRID2_EPS = (0.36,)
HYBRID4_EPS = ladder(2, 25.0 / 16.0)          # 0.028601, 0.385824
HYBRID5_EPS = ladder(3)                       # 0.001597, 0.041275, 0.387435

HYBRID2 = Law("hybrid2", HYBRID2_EPS, [(1.0, -1.0, 0.0)])
HYBRID4 = Law("hybrid4", HYBRID4_EPS, [(1.0, 0.0, -1.0), (1.0, -2.0, 1.0)],
              mix=hybrid4_basis(HYBRID4_EPS))
HYBRID5 = Law("hybrid5", HYBRID5_EPS, [(1.0, 0.0, -1.0), (1.0, -2.0, 1.0)])
LAWS = {law.name: law for law in (HYBRID2, HYBRID4, HYBRID5)}


def get_law(law) -> Law:
    if isinstance(law, Law):
        return law
    try:
        return LAWS[law]
    except KeyError:
        raise ValueError(f"unknown law {law!r}; expected one of "
                         f"{sorted(LAWS)}") from None


def _weights(w, law):
    w = np.asarray(w, dtype=np.float64).reshape(-1)
    if w.size != law.n_w:
        raise ValueError(f"{law.name} takes {law.n_w} weights; got {w.size}")
    return w


def intensity(mu, w, law):
    """I(mu) / I(1) = 1 - sum_j w_j S_j(mu)."""
    law = get_law(law)
    w = _weights(w, law)
    m2 = np.asarray(mu, dtype=np.float64) ** 2
    gen = [np.ones_like(m2), m2, m2 * m2] + [(m2 + e) ** -2 for e in law.eps]
    coef = w @ law.G
    return 1.0 - sum(c * g for c, g in zip(coef, gen))


def _flux_from_occ(occ, w, law):
    """F = 1 - (O_0 - sum_j w_j O_Sj) / (pi - sum_j w_j N_j), from the
    occulted flux O of each generator (trailing axis)."""
    coef = w @ law.G                          # weight of each generator
    blocked = occ[..., 0] - occ @ coef
    return 1.0 - blocked / (_PI - w @ law.norms)


# ===========================================================================
# spherical planet
# ===========================================================================

def _kite_sqarea(z, r):
    """16 x (area of the triangle with sides 1, z, r)^2, sorted (Kahan)
    form; negative when no such triangle exists."""
    a = np.maximum(z, r)
    c = np.minimum(z, r)
    b = np.minimum(a, 1.0)
    a = np.maximum(a, 1.0)
    b, c = np.maximum(b, c), np.minimum(b, c)
    return (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c))


def _sphere_occ(z, r, eps):
    """Occulted flux of [mu^0, mu^2, mu^4, P_eps...] for a circular planet
    of radius r at separation z (broadcast), stacked on a trailing axis."""
    z = np.abs(np.asarray(z, dtype=np.float64))
    r = np.asarray(r, dtype=np.float64)
    z, r = np.broadcast_arrays(z, r)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        r2, z2 = r * r, z * z
        m_none = z >= 1.0 + r
        m_comp = (z <= 1.0 - r) & ~m_none
        m_part = ~(m_none | m_comp)
        m_tot = z <= r - 1.0                  # total occultation (r > 1)
        m_lens = m_part & ~m_tot

        A = np.maximum((r + 1.0 - z) * (1.0 - r + z), _TINY)   # 1 - (z-r)^2
        # 1 - (z+r)^2, the larger of z and r taken from 1 first (exactly,
        # when it is near 1), so that the smaller is not rounded against 1
        Bp = ((1.0 - np.maximum(z, r)) - np.minimum(z, r)) * (1.0 + z + r)
        sqarea = _kite_sqarea(z, r)
        # sqarea > 0 inside the lens, to a few ulps (the floor only keeps
        # kite^2 representable); r^2 - 1 is factored, as z^2 is lost against
        # it when r = 1 and z is small
        kite = np.where(m_lens, np.sqrt(np.maximum(sqarea, _TINY)), 0.0)
        rm1 = (r - 1.0) * (r + 1.0)              # r^2 - 1
        kap0 = np.arctan2(kite, rm1 + z2)        # half-angle at the planet
        kap1 = np.arctan2(kite, z2 - rm1)        # half-angle at the star

        # ---- even terms: s0, s2 and the recursion M0, M2 -> M4 ----------
        alpha = -rm1 - z2                        # 1 - r^2 - z^2
        s0_part = -(kap1 + r2 * kap0 - 0.5 * kite)
        s0 = np.where(m_part, s0_part, np.where(m_comp, -_PI * r2, 0.0))
        eta2 = r2 * (r2 + 2.0 * z2)
        s2_part = 2.0 * s0_part + 2.0 * (
            kap1 + eta2 * kap0 - 0.25 * kite * (1.0 + 5.0 * r2 + z2))
        s2 = np.where(m_part, s2_part,
                      np.where(m_comp, _TWO_PI * r2 * (r2 + 2.0 * z2 - 1.0),
                               0.0))
        M0 = np.where(m_comp, _PI, kap0)
        M2 = np.where(m_comp, _PI * alpha, kap0 * alpha + kite)
        M4 = (6.0 * alpha * M2 + 2.0 * sqarea * M0) / 4.0
        s4 = -(2.0 * r2 * M4 - (2.0 / 3.0) * (alpha * M4 + sqarea * M2))
        s4 = np.where(m_none, 0.0, s4)
        cols = [-s0, -(0.5 * s0 + 0.25 * s2), -(s0 / 3.0 + s2 / 6.0 + s4 / 6.0)]

        # ---- poles -------------------------------------------------------
        # the smaller of A and -Bp is taken from kite^2 = A (-Bp), exact
        inner = m_lens & (A >= -Bp)
        nBp = np.where(inner, kite * kite / A,
                       np.where(-Bp > 0.0, -Bp, kite * kite / A))
        nBp = np.where(m_part, nBp, 1.0)
        for e in eps:
            p = 1.0 + e
            a_ = p - z2 - r2
            K = p + r2 - z2
            Q = (e + A) * (e + Bp)
            # complete transit (Q > 0 there)
            sqQ = np.sqrt(np.where(m_comp, Q, 1.0))
            occ_c = _TWO_PI * r2 / (sqQ * (K + sqQ))
            # partial: J = int_{-kap0}^{kap0} dpsi / (a_ + 2 z r cos psi)
            apb = e + A
            U = apb * nBp
            ratio = (e + Bp) / nBp
            y = ratio * kite * kite / U
            m_atan = m_part & (y > _Y_SW)
            m_log = m_part & (y < -_Y_SW)
            sqU = np.sqrt(U)
            J_atan = (4.0 / np.sqrt(np.abs(Q))
                      * np.arctan2(kite * np.sqrt(np.abs(ratio)), sqU))
            J_log = (4.0 / np.sqrt(np.abs(Q))
                     * np.log((sqU + kite * np.sqrt(np.abs(ratio)))
                              / (2.0 * np.sqrt(np.maximum(z * r * e, _TINY)))))
            T = kite / nBp                    # tan(kap0 / 2)
            y2 = y * y
            G = (1.0 - y / 3.0 + y2 / 5.0 - y * y2 / 7.0 + y2 * y2 / 9.0
                 - y * y2 * y2 / 11.0)
            J_ser = 4.0 * T * G / apb
            J = np.where(m_atan, J_atan, np.where(m_log, J_log, J_ser))
            occ_p = kap1 / (p * e) + (K * J - 2.0 * kap0) / (4.0 * p)
            cols.append(np.where(m_part, occ_p, np.where(m_comp, occ_c, 0.0)))
    return np.stack(cols, axis=-1)


def hybrid_flux(z, r, w, law="hybrid4"):
    """Normalised flux of a star with a hybrid law, occulted by a spherical
    planet of radius r (stellar radii) at sky separation z.

    z, r: scalars or arrays (broadcast); w: the law's shape weights
    (2, 4 or 5 of them); law: "hybrid2", "hybrid4" or "hybrid5".
    Exact for r < 1; for a companion larger than the star the even terms
    lose about r^5 x 1e-16 (tested to r = 2; 2e-12 at r = 10, 2e-7 at
    r = 100).
    """
    law = get_law(law)
    w = _weights(w, law)
    r = np.asarray(r, dtype=np.float64)
    if np.any(r < 0.0):
        raise ValueError("r must be >= 0")
    return _flux_from_occ(_sphere_occ(z, r, law.eps), w, law)


# ===========================================================================
# oblate planet: the ellipse x = x0 + A cos phi, y = y0 + B sin phi
# ===========================================================================
#
# Along it, with z = e^{i phi},
#   s(phi) = x^2 + y^2 = S0 + 2 Re(S1 z) + 2 S2 cos 2 phi,
#   (x dy - y dx)/dphi = K0 + 2 Re(K1 z),
# and the points where s = p are the roots of the self-inversive quartic
#   D(z) = S2 z^4 + S1 z^3 + (S0 - p) z^2 + conj(S1) z + S2.

def _coeffs(x0, y0, a, b):
    S0 = x0 * x0 + y0 * y0 + 0.5 * (a * a + b * b)
    S1 = a * x0 - 1j * (b * y0)
    S2 = 0.25 * (a * a - b * b)
    K0 = a * b
    K1 = 0.5 * (b * x0 - 1j * (a * y0))
    return S0, S1, S2, K0, K1


def _D(z, S0, S1, S2, p):
    return (((S2 * z + S1) * z + (S0 - p)) * z + np.conj(S1)) * z + S2


def _dD(z, S0, S1, S2, p):
    return ((4.0 * S2 * z + 3.0 * S1) * z + 2.0 * (S0 - p)) * z + np.conj(S1)


def _N(z, K0, K1):
    return np.conj(K1) + K0 * z + K1 * z * z


# ---- the quartic's roots: a 4x4 eigensolver, then Aberth ------------------
# The companion matrix is upper Hessenberg already. Balanced, then
# Wilkinson-shifted complex QR with a fixed number of sweeps per deflation
# stage; the last 2x2 block is solved directly. Every operation is
# elementwise over the leading (points) axis.

def _givens(x, y):
    """(c, s), c real: [[c, s], [-conj(s), c]] @ [x, y] = [rho, 0]."""
    xa = np.abs(x)
    rho = np.sqrt(xa * xa + np.abs(y) ** 2)
    ph = np.where(xa > 0, x / np.where(xa > 0, xa, 1.0), 1.0)
    ok = rho > 0
    rs = np.where(ok, rho, 1.0)
    c = np.where(ok, xa / rs, 1.0)
    s = np.where(ok, ph * np.conj(y) / rs, 0.0)
    return c, s


def _wilkinson(a, b, c, d):
    """The eigenvalue of [[a, b], [c, d]] nearer d."""
    h = 0.5 * (a - d)
    disc = np.sqrt(h * h + b * c)
    den = np.where(np.abs(h + disc) >= np.abs(h - disc), h + disc, h - disc)
    ok = den != 0
    return d - np.where(ok, b * c / np.where(ok, den, 1.0), 0.0)


def _qr_sweep(H, m):
    """One shifted QR step on the leading m x m block of H (n, 4, 4)."""
    sig = _wilkinson(H[:, m - 2, m - 2], H[:, m - 2, m - 1],
                     H[:, m - 1, m - 2], H[:, m - 1, m - 1])
    for i in range(m):
        H[:, i, i] -= sig
    rots = []
    for k in range(m - 1):
        c, s = _givens(H[:, k, k], H[:, k + 1, k])
        rots.append((c, s))
        for j in range(k, m):
            x, y = H[:, k, j].copy(), H[:, k + 1, j].copy()
            H[:, k, j] = c * x + s * y
            H[:, k + 1, j] = c * y - np.conj(s) * x
    for k, (c, s) in enumerate(rots):
        for i in range(min(k + 2, m - 1) + 1):
            x, y = H[:, i, k].copy(), H[:, i, k + 1].copy()
            H[:, i, k] = c * x + np.conj(s) * y
            H[:, i, k + 1] = c * y - s * x
    for i in range(m):
        H[:, i, i] += sig


def _qr_seeds(S0, S1, S2, p):
    n = S0.size
    H = np.zeros((n, 4, 4), dtype=np.complex128)
    inv = 1.0 / S2
    H[:, 0, 0] = -S1 * inv
    H[:, 0, 1] = -(S0 - p) * inv
    H[:, 0, 2] = -np.conj(S1) * inv
    H[:, 0, 3] = -1.0
    H[:, 1, 0] = H[:, 2, 1] = H[:, 3, 2] = 1.0
    for _ in range(_QR_BALANCE):              # diagonal similarity
        for i in range(4):
            col = np.abs(H[:, :, i]).sum(1) - np.abs(H[:, i, i])
            row = np.abs(H[:, i, :]).sum(1) - np.abs(H[:, i, i])
            ok = (col > 0) & (row > 0)
            f = np.where(ok, np.sqrt(np.where(ok, row, 1.0)
                                     / np.where(ok, col, 1.0)), 1.0)
            H[:, :, i] *= f[:, None]
            H[:, i, :] /= f[:, None]
    for _ in range(_QR_SWEEPS[0]):
        _qr_sweep(H, 4)
    for _ in range(_QR_SWEEPS[1]):
        _qr_sweep(H, 3)
    a, b, c, d = H[:, 0, 0], H[:, 0, 1], H[:, 1, 0], H[:, 1, 1]
    h = 0.5 * (a - d)
    disc = np.sqrt(h * h + b * c)
    mid = 0.5 * (a + d)
    return np.stack([mid + disc, mid - disc, H[:, 2, 2], H[:, 3, 3]], axis=1)


def _roots(S0, S1, S2, p):
    """The four roots (n, 4) of D at level p: eigenvalue seeds, then
    Aberth iterations (at most _N_ABERTH, stopping once every correction
    is at rounding) and two Newton polishes; a correction that is not
    finite is skipped. A near-double pair converges only to about
    sqrt(eps); the pair formula then uses the pair's sum and product."""
    z = _qr_seeds(S0, S1, S2, p)
    c = [v[:, None] for v in (S0, S1, S2)]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for _ in range(_N_ABERTH):
            rat = _D(z, *c, p) / _dD(z, *c, p)
            diff = z[:, :, None] - z[:, None, :]
            inv = np.where(diff != 0, 1.0 / np.where(diff != 0, diff, 1.0), 0.0)
            s = inv.sum(axis=2)
            step = rat / (1.0 - rat * s)
            ok = np.isfinite(step)
            z = np.where(ok, z - step, z)
            if not np.any(ok & (np.abs(step) > 4.0 * _EPS * np.abs(z))):
                break
        for _ in range(2):
            new = z - _D(z, *c, p) / _dD(z, *c, p)
            z = np.where(np.isfinite(new), new, z)
    return z


# ---- limb crossings and the arcs -------------------------------------------

def _mod2pi(t):
    return t - _TWO_PI * np.floor(t / _TWO_PI)


def _g_dg(phi, x0, y0, a, b):
    s, c = np.sin(phi), np.cos(phi)
    x, y = x0 + a * c, y0 + b * s
    return x * x + y * y - 1.0, 2.0 * (-x * a * s + y * b * c)


def _crossings(x0, y0, a, b):
    """(count, phi_a <= phi_b): the limb crossings on the ellipse
    parameter. Candidates are roots of D at level 1 with |z| ~ 1, polished
    by Newton on g = s - 1 and accepted on the residual; crossings closer
    than _PAIR_TOL are a tangency and dropped together."""
    S0, S1, S2, _, _ = _coeffs(x0, y0, a, b)
    z = _roots(S0, S1, S2, 1.0)
    X0, Y0, A, B = (v[:, None] for v in (x0, y0, a, b))
    sc = 1.0 + (np.abs(X0) + A) ** 2 + (np.abs(Y0) + B) ** 2
    cand = np.abs(np.abs(z) - 1.0) < _Z_TOL
    phi = np.where(cand, np.angle(z), 0.0)
    g, dg = _g_dg(phi, X0, Y0, A, B)
    for _ in range(3):
        with np.errstate(divide="ignore", invalid="ignore"):
            trial = phi - np.where(dg != 0, g / np.where(dg != 0, dg, 1.0), 0.0)
        gt, dgt = _g_dg(trial, X0, Y0, A, B)
        better = np.abs(gt) < np.abs(g)
        phi = np.where(better, trial, phi)
        g = np.where(better, gt, g)
        dg = np.where(better, dgt, dg)
    valid = cand & (np.abs(g) <= _G_TOL * sc)
    phi = _mod2pi(phi)
    drop = np.zeros_like(valid)
    for i in range(4):
        for j in range(i + 1, 4):
            d = np.abs(_mod2pi(phi[:, i] - phi[:, j] + _PI) - _PI)
            close = valid[:, i] & valid[:, j] & (d < _PAIR_TOL)
            drop[:, i] |= close
            drop[:, j] |= close
    valid &= ~drop
    key = np.sort(np.where(valid, phi, np.inf), axis=1)
    return valid.sum(axis=1), key[:, 0], key[:, 1]


def _arcs(x0, y0, a, b, phi_a, phi_b):
    """(lo, hi) of the ellipse arc inside the star, and the length of the
    limb arc inside the ellipse, for two crossings. The ellipse arc is
    chosen at the midpoint of the longer of its two arcs, far from both
    crossings; the limb arc inside the ellipse is always the shorter one
    (the planet lies within a circle of radius A < 1)."""
    long_first = (phi_b - phi_a) >= _PI
    mid = np.where(long_first, 0.5 * (phi_a + phi_b),
                   0.5 * (phi_a + phi_b) + _PI)
    xm, ym = x0 + a * np.cos(mid), y0 + b * np.sin(mid)
    mid_in = xm * xm + ym * ym < 1.0
    first = np.where(long_first, mid_in, ~mid_in)  # [phi_a, phi_b] inside
    lo = np.where(first, phi_a, phi_b)
    hi = np.where(first, phi_b, phi_a + _TWO_PI)
    ta = np.arctan2(y0 + b * np.sin(phi_a), x0 + a * np.cos(phi_a))
    tb = np.arctan2(y0 + b * np.sin(phi_b), x0 + a * np.cos(phi_b))
    dt = np.abs(_mod2pi(ta) - _mod2pi(tb))
    return lo, hi, np.minimum(dt, _TWO_PI - dt)


# ---- even terms: Fourier coefficients of h_n(s) K ---------------------------

def _q_coeffs(S0, S1, S2, K0, K1):
    """q[:, n, m], n = 0..2, m = 0..5: h_n(s) K = Re sum_m q_m e^{i m phi}
    (positive frequencies), with h_n(s) = sum_p alpha[n, p] s^p,
    alpha[n, p] = (-1)^p C(n, p) / (2 (p + 1))."""
    n_max, H = 2, 5
    W = 2 * H + 1
    T = np.zeros((S0.size, W), dtype=np.complex128)
    T[:, H - 1], T[:, H], T[:, H + 1] = np.conj(K1), K0, K1
    Sk = [S2 + 0j, np.conj(S1), S0 + 0j, S1, S2 + 0j]       # j = -2..2
    q = np.zeros((S0.size, n_max + 1, H + 1), dtype=np.complex128)
    for p in range(n_max + 1):
        if p > 0:                                          # T <- S * T
            Tn = np.zeros_like(T)
            for jj, j in enumerate(range(-2, 3)):
                if j >= 0:
                    Tn[:, j:] += Sk[jj][:, None] * T[:, :W - j]
                else:
                    Tn[:, :j] += Sk[jj][:, None] * T[:, -j:]
            T = Tn
        for n in range(p, n_max + 1):
            alpha = (-1) ** p * math.comb(n, p) / (2.0 * (p + 1))
            q[:, n, :] += alpha * T[:, H:]
    return q


def _even_arc(q, lo, hi):
    """Psi_n(hi) - Psi_n(lo), n = 0..2, Psi_n the exact antiderivative."""
    mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
    m = np.arange(1, q.shape[-1])
    osc = (4.0 / m) * np.real(q[:, :, 1:] * np.exp(1j * m * mid[:, None, None])) \
        * np.sin(m * half[:, None, None])
    return q[:, :, 0].real * (hi - lo)[:, None] + osc.sum(axis=-1)


# ---- poles: partial fractions over the roots of D ----------------------------

def _bairstow(S0, S1, S2, p, sig, pi_):
    """Newton steps on (sigma, pi) for the quadratic factor
    z^2 - sigma z + pi of D (the remainder of D mod it driven to 0)."""
    dc = (S2 + 0j, np.conj(S1), (S0 - p) + 0j, S1, S2 + 0j)   # low -> high
    for _ in range(_N_BAIRSTOW):
        zero = np.zeros_like(sig)
        al, be = zero, zero + 1.0          # z^k mod q = al z + be
        das = dbs = dap = dbp = zero
        r1, r0 = dc[0] * al, dc[0] * be
        J11 = J12 = J21 = J22 = zero
        for k in range(1, 5):
            al, be, das, dbs, dap, dbp = (al * sig + be, -al * pi_,
                                          das * sig + al + dbs, -das * pi_,
                                          dap * sig + dbp, -dap * pi_ - al)
            r1 = r1 + dc[k] * al
            r0 = r0 + dc[k] * be
            J11, J12 = J11 + dc[k] * das, J12 + dc[k] * dap
            J21, J22 = J21 + dc[k] * dbs, J22 + dc[k] * dbp
        det = J11 * J22 - J12 * J21
        ok = det != 0
        det = np.where(ok, det, 1.0)
        sig = np.where(ok, sig - (J22 * r1 - J12 * r0) / det, sig)
        pi_ = np.where(ok, pi_ - (-J21 * r1 + J11 * r0) / det, pi_)
    return sig, pi_


def _pair_terms(S0, S1, S2, K0, K1, p, sig, pi_):
    """For the quadratic factor q = z^2 - sigma z + pi of D, with
    N = aN z + bN and Q = D / q = aQ z + bQ (mod q): the coefficients of
    -N / (i D) = (A z + B) / q + ..., as (A, C) with C = B + A sigma / 2,
    and their common denominator den. Symmetric in the two roots, so no
    division by their separation."""
    Q2 = S2 + 0j
    Q1 = S1 + sig * Q2
    Q0 = sig * Q1 + (S0 - p) - pi_ * Q2
    aN, bN = K1 * sig + K0, np.conj(K1) - K1 * pi_
    aQ, bQ = Q2 * sig + Q1, Q0 - Q2 * pi_
    den = aQ * aQ * pi_ + aQ * bQ * sig + bQ * bQ
    den = np.where(den != 0, den, 1.0)
    num_A = aN * bQ - bN * aQ          # = den * sum over the pair of N / D'
    num_C = 0.5 * (2.0 * aN * aQ * pi_ + (aN * bQ + bN * aQ) * sig
                   + 2.0 * bN * bQ)
    return num_A / den, num_C / den


def _pole_inside(S0, S1, S2, K0, K1, p, z):
    """Occulted flux of (mu^2 + eps)^-2 for an ellipse wholly on the disc:
    the whole contour, -(pi/p) Re sum_{|z_r| < 1} N(z_r) / D'(z_r). No root
    lies on the unit circle and they pair as (z, 1/conj z), so exactly two
    lie inside; their sum is taken as one quadratic factor."""
    inner = np.abs(z) < 1.0
    sig = np.where(inner, z, 0.0).sum(axis=1)
    pi_ = np.where(inner, z, 1.0).prod(axis=1)
    sig, pi_ = _bairstow(S0, S1, S2, p, sig, pi_)
    pair_sum, _ = _pair_terms(S0, S1, S2, K0, K1, p, sig, pi_)
    return -_PI * pair_sum.real / p


def _dlog(z, lo, hi, inner, w1, w2):
    """Change of log(e^{i phi} - z) for phi from lo to hi, continuous along
    the arc (which never meets z)."""
    dm = np.log(np.abs(w2 - z)) - np.log(np.abs(w1 - z))
    zi = np.where(inner, z, 0.5)
    zo = np.where(inner, 2.0, z)
    d_in = (hi - lo) + np.angle(1.0 - zi / w2) - np.angle(1.0 - zi / w1)
    d_out = np.angle(1.0 - w2 / zo) - np.angle(1.0 - w1 / zo)
    return dm + 1j * np.where(inner, d_in, d_out)


def _J(u, Delta):
    """int du / (u^2 - Delta) = -(1/u) sum_k (Delta/u^2)^k / (2k + 1)."""
    x = Delta / (u * u)
    acc = np.zeros_like(u)
    for k in range(_N_SERIES - 1, -1, -1):
        acc = acc * x + 1.0 / (2 * k + 1)
    return -acc / u


def _pole_arc(S0, S1, S2, K0, K1, p, z, lo, hi):
    """(1/2p) int_lo^hi K / (p - s) dphi along one ellipse arc inside the
    star, as partial fractions over the roots z (n, 4) of D:
    (1/2p) Re sum_r R_r [log(e^{i phi} - z_r)], R_r = -N(z_r) / (i D'(z_r))."""
    n = S0.size
    rows = np.arange(n)
    w1, w2 = np.exp(1j * lo), np.exp(1j * hi)
    # near-double pairs: up to two disjoint, scored by |Delta| / d^2, d the
    # distance from the pair's centre to the arc
    used = np.zeros((n, 4), dtype=bool)
    pairs = []
    PJ = np.array([j for j, _ in _PAIRS])
    PL = np.array([l for _, l in _PAIRS])
    for _ in range(2):
        score = np.full((n, len(_PAIRS)), np.inf)
        for m, (j, l) in enumerate(_PAIRS):
            ok = ~used[:, j] & ~used[:, l]
            Dl = 0.25 * (z[:, j] - z[:, l]) ** 2
            c = 0.5 * (z[:, j] + z[:, l])
            ang = _mod2pi(np.angle(c))
            on = (((ang >= lo) & (ang <= hi))
                  | ((ang + _TWO_PI >= lo) & (ang + _TWO_PI <= hi)))
            d_end = np.minimum(np.abs(w1 - c), np.abs(w2 - c))
            d = np.where(on, np.abs(np.abs(c) - 1.0), d_end)
            score[:, m] = np.where(ok, np.abs(Dl) / np.maximum(d * d, 1e-300),
                                   np.inf)
        best = np.argmin(score, axis=1)
        good = score[rows, best] <= _PAIR_RATIO
        jj, ll = PJ[best], PL[best]
        used[rows, jj] |= good
        used[rows, ll] |= good
        pairs.append((jj, ll, good))

    ell = np.zeros(n, dtype=np.complex128)
    c = (S0, S1, S2)
    # isolated roots. A root beyond |z| = 2 is taken through its partner
    # zeta = 1/conj(z) inside the disc: D and N are self-reciprocal, so
    # R(z) = conj(R(zeta)) and the log's change is that of
    # log(1 - e^{i phi} conj(zeta)); a nearly circular outline's far root
    # is ~1/S2, whose log would otherwise lose digits.
    for k in range(4):
        zk = z[:, k]
        a2 = np.abs(zk) ** 2
        far = a2 > 4.0
        zeta = np.where(far, zk / np.where(far, a2, 1.0), zk)
        new = zeta - _D(zeta, *c, p) / _dD(zeta, *c, p)
        zeta = np.where(np.isfinite(new) & (np.abs(_D(new, *c, p))
                                            < np.abs(_D(zeta, *c, p))),
                        new, zeta)
        dD = _dD(zeta, *c, p)
        live = ~used[:, k] & (dD != 0)
        R = -_N(zeta, K0, K1) / (1j * np.where(live, dD, 1.0))
        b2, b1 = 1.0 - w2 * np.conj(zeta), 1.0 - w1 * np.conj(zeta)
        dlw = np.log(np.abs(b2)) - np.log(np.abs(b1)) \
            + 1j * (np.angle(b2) - np.angle(b1))
        dln = _dlog(zeta, lo, hi, np.abs(zeta) < 1.0, w1, w2)
        ell += np.where(live & ~far, R * dln, 0.0) \
            + np.where(live & far, np.conj(R) * dlw, 0.0)

    # near-double pairs, each as one quadratic factor q = z^2 - sigma z + pi:
    # int (A z + B)/q dz = (A/2) log q + C J(z - sigma/2), with the winding
    # of log q taken from the two roots' own continued logs
    for jj, ll, good in pairs:
        g = np.nonzero(good)[0]
        if g.size == 0:
            continue
        z1, z2 = z[g, jj[g]], z[g, ll[g]]
        s0, s1, s2, k0, k1 = S0[g], S1[g], S2[g], K0[g], K1[g]
        sig, pi_ = _bairstow(s0, s1, s2, p, z1 + z2, z1 * z2)
        pair_sum, Ccoef = _pair_terms(s0, s1, s2, k0, k1, p, sig, pi_)
        Acoef = 1j * pair_sum
        Ccoef = 1j * Ccoef
        lo_, hi_, a1, a2_ = lo[g], hi[g], w1[g], w2[g]
        qh = a2_ * a2_ - sig * a2_ + pi_
        ql = a1 * a1 - sig * a1 + pi_
        darg = np.angle(qh / ql)
        rough = (_dlog(z1, lo_, hi_, np.abs(z1) < 1.0, a1, a2_).imag
                 + _dlog(z2, lo_, hi_, np.abs(z2) < 1.0, a1, a2_).imag)
        wind = np.round((rough - darg) / _TWO_PI)
        LQ = np.log(np.abs(qh)) - np.log(np.abs(ql)) \
            + 1j * (darg + _TWO_PI * wind)
        Delta = 0.25 * sig * sig - pi_
        cen = 0.5 * sig
        Jd = _J(a2_ - cen, Delta) - _J(a1 - cen, Delta)
        ell[g] += 0.5 * Acoef * LQ + Ccoef * Jd
    return ell.real / (2.0 * p)


def _ellipse_occ(x0, y0, a, b, eps):
    """Occulted flux of [mu^0, mu^2, mu^4, P_eps...] for the ellipse with
    centre (x0, y0) in its principal frame and semi-axes a, b (1-D arrays,
    a > b, at most two limb crossings), shape (n, 3 + K)."""
    n = x0.size
    occ = np.zeros((n, 3 + len(eps)))
    d = np.hypot(x0, y0)
    fast_in = d + a < 1.0
    general = ~fast_in & (d - a < 1.0)
    cnt = np.zeros(n, dtype=int)
    phi_a, phi_b = np.zeros(n), np.zeros(n)
    gi = np.nonzero(general)[0]
    if gi.size:
        cnt[gi], phi_a[gi], phi_b[gi] = _crossings(x0[gi], y0[gi], a[gi], b[gi])
    partial = general & (cnt == 2)
    inside = fast_in | (general & (cnt != 2) & (d < 1.0))

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        ii = np.nonzero(inside)[0]
        if ii.size:
            S0, S1, S2, K0, K1 = _coeffs(x0[ii], y0[ii], a[ii], b[ii])
            q = _q_coeffs(S0, S1, S2, K0, K1)
            occ[ii, :3] = _TWO_PI * q[:, :, 0].real
            for k, e in enumerate(eps):
                z = _roots(S0, S1, S2, 1.0 + e)
                occ[ii, 3 + k] = _pole_inside(S0, S1, S2, K0, K1, 1.0 + e, z)

        pi = np.nonzero(partial)[0]
        if pi.size:
            px0, py0, pa, pb = x0[pi], y0[pi], a[pi], b[pi]
            lo, hi, dth = _arcs(px0, py0, pa, pb, phi_a[pi], phi_b[pi])
            S0, S1, S2, K0, K1 = _coeffs(px0, py0, pa, pb)
            q = _q_coeffs(S0, S1, S2, K0, K1)
            occ[pi, :3] = _even_arc(q, lo, hi) \
                + dth[:, None] / (2.0 * np.arange(1, 4))
            for k, e in enumerate(eps):
                p = 1.0 + e
                z = _roots(S0, S1, S2, p)
                occ[pi, 3 + k] = _pole_arc(S0, S1, S2, K0, K1, p, z, lo, hi) \
                    + dth / (2.0 * p * e)
    return occ


def _oblate_occ(x0, y0, r, f, eps):
    """Occulted generator fluxes for an oblate planet, broadcast over
    (x0, y0, r, f); trailing axis [mu^0, mu^2, mu^4, P_eps...]."""
    x0, y0, r, f = np.broadcast_arrays(*(np.asarray(v, dtype=np.float64)
                                         for v in (x0, y0, r, f)))
    shape = x0.shape
    x0, y0, r, f = (v.ravel() for v in (x0, y0, r, f))
    occ = np.empty((x0.size, 3 + len(eps)))
    small = f < _F_SW
    if np.any(small):
        occ[small] = _sphere_occ(np.hypot(x0[small], y0[small]), r[small], eps)
    big = ~small
    if np.any(big):
        a = r[big] / np.sqrt(1.0 - f[big])
        occ[big] = _ellipse_occ(x0[big], y0[big], a, a * (1.0 - f[big]), eps)
    return occ.reshape(shape + (3 + len(eps),))


def hybrid_flux_oblate(X, Y, r, f, theta, w, law="hybrid4"):
    """Normalised flux of a star with a hybrid law, occulted by an oblate
    planet whose projected outline is an ellipse.

    X, Y: the planet's sky position (stellar radii), X along its direction
    of motion. r: area-equivalent radius sqrt(A B), with semi-axes
    A = r / sqrt(1 - f) and B = A (1 - f). f: projected flattening, in
    [0, 1). theta: angle of the long axis from X, in radians. All of these
    broadcast. w: the law's shape weights; law: "hybrid2", "hybrid4" or
    "hybrid5". Requires r <= (1 - f)^(3/2).
    """
    law = get_law(law)
    w = _weights(w, law)
    r = np.asarray(r, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    if np.any(f < 0.0) or np.any(f >= 1.0):
        raise ValueError("f must be in [0, 1)")
    if np.any(r < 0.0) or np.any(r > (1.0 - f) ** 1.5 * (1.0 + 1e-12)):
        raise ValueError("need 0 <= r <= (1 - f)^(3/2): the planet must cross "
                         "the stellar limb at most twice")
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    c, s = np.cos(theta), np.sin(theta)
    x0, y0 = X * c + Y * s, -X * s + Y * c          # principal frame
    return _flux_from_occ(_oblate_occ(x0, y0, r, f, law.eps), w, law)

