"""Native bounded Gaussian fitting and batched exact/center aperture kernels.

Only small star cutouts are converted to native contiguous doubles. FITS images
are neither copied in full nor reopened. The solver uses the installed SciPy
LAPACK SVD through one callback per iteration; the TRF loop and derivatives run
in C++. Callback failures propagate to the caller.
"""
import ctypes
from functools import lru_cache

import numpy as np
from scipy.linalg import svd
from scipy.optimize import OptimizeResult

from .limbdarkfast_cpp import _library

_pointer = ctypes.POINTER(ctypes.c_double)
_svd_type = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_size_t, _pointer, _pointer,
                            _pointer, _pointer, _pointer)


@lru_cache(maxsize=1)
def prepare():
    lib = _library()
    array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
    lib.exotic_psf_fit.argtypes = [ctypes.c_size_t, array, array, array, array,
                                 array, array, ctypes.c_int, _svd_type, array,
                                 array, ctypes.c_char_p, ctypes.c_size_t]
    lib.exotic_psf_fit.restype = ctypes.c_int
    lib.exotic_aperture_grid.argtypes = [array, ctypes.c_int, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
        ctypes.c_size_t, array, ctypes.c_size_t, array, array, ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_int, array, array,
        ctypes.c_char_p, ctypes.c_size_t]
    lib.exotic_aperture_grid.restype = ctypes.c_int
    return lib


def _flat(a):
    return np.ascontiguousarray(a, dtype=np.float64).reshape(-1)


def fit_gaussian(x, y, data, initial, bounds, method='2-point'):
    x, y, data = _flat(x), _flat(y), _flat(data)
    initial, lo, hi = _flat(initial), _flat(bounds[0]), _flat(bounds[1])
    if initial.size != 7 or lo.size != 7 or hi.size != 7:
        raise ValueError('Gaussian PSF fitting requires seven parameters and bounds.')
    if x.size != y.size or x.size != data.size or not x.size:
        raise ValueError('PSF cutout and coordinate sizes differ or are empty.')
    if method not in ('2-point', '3-point'):
        raise ValueError('Native PSF Jacobian must be 2-point or 3-point.')
    failure = []

    @_svd_type
    def decompose(n, matrix, residuals, s, vectors, uf):
        try:
            a = np.ctypeslib.as_array(matrix, shape=(n * 7,)).reshape(n, 7)
            f = np.ctypeslib.as_array(residuals, shape=(n,))
            u, singular, vt = svd(a, full_matrices=False)
            np.ctypeslib.as_array(s, shape=(7,))[:] = singular
            np.ctypeslib.as_array(vectors, shape=(49,))[:] = vt.T.reshape(-1)
            np.ctypeslib.as_array(uf, shape=(7,))[:] = u.T.dot(f)
            return 0
        except BaseException as exc:
            failure.append(exc)
            return 1

    output, info = np.empty(7), np.empty(5)
    error = ctypes.create_string_buffer(512)
    status = prepare().exotic_psf_fit(data.size, x, y, data, initial, lo, hi,
        int(method[0]), decompose, output, info, error, len(error))
    if failure:
        raise failure[0]
    if status:
        raise ValueError(error.value.decode('utf-8', errors='replace'))
    return OptimizeResult(x=output, cost=info[0], nfev=int(info[1]),
                          njev=int(info[2]), status=int(info[3]),
                          optimality=info[4], success=info[3] > 0)


def aperture_grid(data, xc, yc, radii, inner, outer, *, center=False,
                  ptol=99., sigma=3., max_iters=3):
    """Return raw sums/areas/intersection flags and sky median/scatter/area.

    inner/outer have shape (number of apertures, number of annuli). A zero
    outer radius requests no sky subtraction; NaNs retain an invalid cell.
    """
    radii = np.asarray(radii, dtype=float).reshape(-1)
    inner, outer = np.asarray(inner, dtype=float), np.asarray(outer, dtype=float)
    if inner.shape != outer.shape or inner.ndim != 2 or inner.shape[0] != radii.size:
        raise ValueError('Sky geometry must have aperture-by-annulus shape.')
    candidates = np.concatenate((radii, outer.reshape(-1)))
    valid = candidates[np.isfinite(candidates) & (candidates > 0)]
    extent = float(np.max(valid)) if valid.size else 0.
    height, width = data.shape
    # Leave sky weights/pixel positions in original coordinates. Empty crops
    # retain the reference out-of-frame behavior without reading an image.
    x0 = min(width, max(0, int(np.floor(xc - extent + .5))))
    y0 = min(height, max(0, int(np.floor(yc - extent + .5))))
    x1 = min(width, max(0, int(np.ceil(xc + extent + .5))))
    y1 = min(height, max(0, int(np.ceil(yc + extent + .5))))
    image = np.ascontiguousarray(data[y0:y1, x0:x1], dtype=np.float64)
    ap, bg = np.empty((radii.size, 3)), np.empty((*inner.shape, 3))
    error = ctypes.create_string_buffer(512)
    status = prepare().exotic_aperture_grid(image.reshape(-1), image.shape[1], image.shape[0],
        x0, y0, float(xc), float(yc), radii.size, _flat(radii), inner.shape[1],
        _flat(inner), _flat(outer), int(center), float(ptol), float(sigma),
        int(max_iters), ap.reshape(-1), bg.reshape(-1), error, len(error))
    if status:
        raise ValueError(error.value.decode('utf-8', errors='replace'))
    return ap, bg
