"""NumPy/ctypes interface to EXOTIC's native Claret and hybrid4 models.

The bundled library is built automatically when missing or older than its sources.
No Python extension installation or additional Python dependency is needed.
An explicitly requested C++ evaluation raises if its library cannot be loaded;
it never substitutes the NumPy implementation.
"""
import ctypes
from functools import lru_cache

import numpy as np
from .native_build import ABI_VERSION, NATIVE, ensure_hybrid4_native, library_filename


def library_path():
    return NATIVE / library_filename()


@lru_cache(maxsize=1)
def _library():
    path = library_path()
    ensure_hybrid4_native(output_dir=path.parent)
    try:
        library = ctypes.CDLL(str(path))
    except OSError as exc:
        raise RuntimeError(
            f'Cannot load C++ hybrid4 library {path}. Build it from the repository '
            'console with: python scripts/build_hybrid4_cpp.py'
        ) from exc
    library.exotic_hybrid4_abi_version.argtypes = []
    library.exotic_hybrid4_abi_version.restype = ctypes.c_int
    if library.exotic_hybrid4_abi_version() != ABI_VERSION:
        raise RuntimeError('C++ hybrid4 ABI differs from this checkout; rebuild with scripts/build_hybrid4_cpp.py.')
    array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
    library.exotic_hybrid4_flux.argtypes = [
        ctypes.c_size_t, array, array, array, array, ctypes.c_char_p, ctypes.c_size_t,
    ]
    library.exotic_hybrid4_flux.restype = ctypes.c_int
    library.exotic_hybrid4_transit.argtypes = [
        ctypes.c_size_t, array, array, ctypes.c_size_t, array, array,
        ctypes.c_size_t, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        array, ctypes.c_char_p, ctypes.c_size_t,
    ]
    library.exotic_hybrid4_transit.restype = ctypes.c_int
    library.exotic_claret_transit.argtypes = library.exotic_hybrid4_transit.argtypes
    library.exotic_claret_transit.restype = ctypes.c_int
    byte_array = np.ctypeslib.ndpointer(dtype=np.uint8, ndim=1, flags='C_CONTIGUOUS')
    library.exotic_orbital_replay.argtypes = [
        ctypes.c_size_t, ctypes.c_size_t, array, array, array, ctypes.c_size_t,
        array, array, array, array, byte_array, ctypes.c_int, ctypes.c_int,
        ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_size_t, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, array, array, array, array,
        ctypes.c_char_p, ctypes.c_size_t,
    ]
    library.exotic_orbital_replay.restype = ctypes.c_int
    integer_array = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
    library.exotic_loglikelihood.argtypes = [
        ctypes.c_size_t, ctypes.c_size_t, array, ctypes.c_size_t, array, array,
        ctypes.c_size_t, array, integer_array, array, array, array, array, array,
        byte_array, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
        ctypes.c_size_t, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        array, ctypes.c_char_p, ctypes.c_size_t,
    ]
    library.exotic_loglikelihood.restype = ctypes.c_int
    library.exotic_low_pixel_mask.argtypes = [
        ctypes.c_size_t, array, ctypes.c_int64, ctypes.c_double, byte_array,
        ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_size_t,
    ]
    library.exotic_low_pixel_mask.restype = ctypes.c_int
    library.exotic_repair_pixels.argtypes = [
        ctypes.c_size_t, ctypes.c_size_t, array, byte_array, ctypes.c_int, byte_array,
        ctypes.POINTER(ctypes.c_size_t), ctypes.c_char_p, ctypes.c_size_t,
    ]
    library.exotic_repair_pixels.restype = ctypes.c_int
    return library


def hybrid4_flux(separation, radius, weights):
    """Spherical hybrid4 flux; separation and radius broadcast like NumPy."""
    z, r = np.broadcast_arrays(np.asarray(separation, dtype=np.float64),
                               np.asarray(radius, dtype=np.float64))
    shape = z.shape
    z = np.ascontiguousarray(z).reshape(-1)
    r = np.ascontiguousarray(r).reshape(-1)
    w = np.ascontiguousarray(weights, dtype=np.float64).reshape(-1)
    if w.size != 4:
        raise ValueError('hybrid4 requires four shape weights.')
    output = np.empty(z.size, dtype=np.float64)
    error = ctypes.create_string_buffer(512)
    status = _library().exotic_hybrid4_flux(z.size, z, r, w, output, error, len(error))
    if status:
        raise ValueError(error.value.decode('utf-8', errors='replace'))
    return output.reshape(shape)


def _native_transit(times, values, law, *, samples=7, tolerance=1.e-5,
                    minimum_exposure_days=1.0/86400.0, padding=2.0):
    """One native call for orbit, occulted flux and midpoint exposure integration.

    API normalization follows EXOTIC: a scalar exposure broadcasts; a mismatched
    or nonnumeric exposure means instantaneous flux. Native failures propagate.
    """
    time_array = np.asarray(times, dtype=np.float64)
    shape = time_array.shape
    t = np.ascontiguousarray(time_array).reshape(-1)
    exposures = np.empty(0, dtype=np.float64)
    raw = values.get('_exposure_time_days')
    if raw is not None:
        try:
            candidate = np.asarray(raw, dtype=np.float64)
            if candidate.shape == () or candidate.shape == shape:
                exposures = np.ascontiguousarray(candidate).reshape(-1)
        except (TypeError, ValueError):
            pass
    parameters = np.array([values[key] for key in ('rprs','per','ars','ecc','inc','omega','tmid')], dtype=np.float64)
    prefix = 'h' if law == 'hybrid4' else 'u'
    weights = np.array([values[f'{prefix}{i}'] for i in range(4)], dtype=np.float64)
    output = np.empty(t.size, dtype=np.float64)
    error = ctypes.create_string_buffer(512)
    function = getattr(_library(), f'exotic_{law}_transit')
    status = function(
        t.size, t, exposures, exposures.size, parameters, weights,
        samples, tolerance, minimum_exposure_days, padding, output, error, len(error),
    )
    if status:
        raise ValueError(error.value.decode('utf-8', errors='replace'))
    return output.reshape(shape)


def hybrid4_transit(times, values, **settings):
    return _native_transit(times, values, 'hybrid4', **settings)


def claret_transit(times, values, **settings):
    """Existing precision=3 Claret model, orbit and exposure integration."""
    return _native_transit(times, values, 'claret', **settings)
