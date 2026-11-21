"""Native detector histogram and repair, with existing SciPy EDT tie handling."""
import ctypes

import numpy as np
from scipy.ndimage import distance_transform_edt

from .limbdarkfast_cpp import _library


def low_pixel_mask(image, gap_width=17, population_threshold=25):
    values = np.asarray(image, dtype=float)
    shape = values.shape
    values = np.ascontiguousarray(values).reshape(-1)
    mask = np.empty(shape, dtype=np.uint8)
    threshold = ctypes.c_double()
    error = ctypes.create_string_buffer(512)
    status = _library().exotic_low_pixel_mask(
        values.size, values.reshape(-1), gap_width, population_threshold, mask.reshape(-1),
        ctypes.byref(threshold), error, len(error),
    )
    if status:
        raise ValueError(error.value.decode('utf-8', errors='replace'))
    return mask.view(np.bool_), threshold.value if np.isfinite(threshold.value) else None


def repair_frame(image, reference, *, copy=True):
    # The public default owns a new result. Frame loaders may reuse their own
    # writable float64 buffer; dtype conversion remains the original float64.
    values = np.asarray(image)
    if copy or values.dtype != np.float64 or not values.flags.c_contiguous or not values.flags.writeable:
        values = np.array(values, dtype=float, order='C', copy=True)
    source_mask = np.asarray(reference['mask'], dtype=bool)
    if source_mask.shape != values.shape:
        raise ValueError('Bad-pixel mask shape differs from current frame')
    source_mask = np.ascontiguousarray(source_mask).view(np.uint8)
    mask = np.empty(values.shape, dtype=np.uint8)
    unresolved = ctypes.c_size_t()
    error = ctypes.create_string_buffer(512)
    status = _library().exotic_repair_pixels(
        values.shape[0], values.shape[1], values.reshape(-1), source_mask.reshape(-1),
        bool(reference.get('detect_low')), mask.reshape(-1), ctypes.byref(unresolved), error, len(error),
    )
    if status:
        raise ValueError(error.value.decode('utf-8', errors='replace'))
    if unresolved.value:
        # This is the same compiled nearest-original-good-pixel operation used
        # by the reference. Mask values 1/2 both remain bad for EDT, so previous
        # local repairs never become sources, including for tied distances.
        nearest = distance_transform_edt(mask, return_distances=False, return_indices=True)
        y, x = np.nonzero(mask == 2)
        values[y, x] = values[tuple(nearest[:, y, x])]
    return values
