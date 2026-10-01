"""Flux-preserving Bayer superpixels and their FITS geometry."""

import re
import numpy as np

from .bad_pixels import repair_frame


def bin2x2(image):
    """Sum R, G1, G2, B without interpolation; discard incomplete edge tiles."""
    image = np.asarray(image, dtype=np.float64)
    if image.ndim != 2 or min(image.shape) < 2:
        raise ValueError('bin2x2 requires a 2-D Bayer frame at least 2x2 pixels')
    height, width = (size // 2 * 2 for size in image.shape)
    return image[:height, :width].reshape(height // 2, 2, width // 2, 2).sum(axis=(1, 3))


def repair_bayer_pixels(image, reference):
    """Repair each Bayer phase using neighbours of the same colour."""
    if reference is None:
        return image
    result = np.array(image, dtype=float, copy=True)
    for y in range(2):
        for x in range(2):
            phase = dict(reference, mask=reference['mask'][y::2, x::2])
            result[y::2, x::2] = repair_frame(result[y::2, x::2], phase)
    return result


def binned_header(header):
    """Map zero-based pixel centres with p_out=(p_in-0.5)/2, including SIP."""
    header = header.copy()
    native_bins = None
    for xkey, ykey in (('XBINNING', 'YBINNING'), ('XBINING', 'YBINING'), ('CCDXBIN', 'CCDYBIN')):
        if xkey in header and ykey in header:
            native_bins = (int(header[xkey]), int(header[ykey]))
            break
    if native_bins is None:
        for key in ('CCDSUM', 'BINNING'):
            parts = re.split(r'[xX,\s]+', str(header.get(key, '')).strip())
            if len(parts) == 2:
                native_bins = tuple(int(p) for p in parts)
                break
    native_bins = native_bins or (1, 1)
    for key in list(header):
        if re.fullmatch(r'CRPIX[12][A-Z]?', key):
            header[key] = (float(header[key]) + 0.5) / 2
        elif re.fullmatch(r'(CD[12]_[12]|CDELT[12])[A-Z]?', key):
            header[key] = float(header[key]) * 2
        elif re.fullmatch(r'(A|B|AP|BP)_\d+_\d+', key):
            i, j = map(int, key.split('_')[-2:])
            header[key] = float(header[key]) * 2 ** (i + j - 1)
    for key in ('PIXSCALE', 'PIXSCAL1', 'PIXSCAL2', 'SECPIX', 'SECPIX1', 'SECPIX2', 'XPIXSZ', 'YPIXSZ'):
        if key in header:
            header[key] = float(header[key]) * 2
    for key in ('XBINNING', 'YBINNING', 'XBINING', 'YBINING', 'CCDXBIN', 'CCDYBIN'):
        if key in header:
            header[key] = int(header[key]) * 2
    for key in ('BINNING', 'CCDSUM'):
        if key in header:
            parts = re.split(r'[xX,\s]+', str(header[key]).strip())
            if len(parts) == 2:
                header[key] = ('x' if key == 'BINNING' else ' ').join(str(int(p) * 2) for p in parts)
    header['XBINNING'] = native_bins[0] * 2
    header['YBINNING'] = native_bins[1] * 2
    for key in ('SATURATE', 'RDNOISE', 'READNOIS'):
        if key in header:
            header[key] = float(header[key]) * (4 if key == 'SATURATE' else 2)
    for key in ('BAYERPAT', 'XBAYROFF', 'YBAYROFF', 'CHECKSUM', 'DATASUM'):
        header.remove(key, ignore_missing=True, remove_all=True)
    header['EXOBIN'] = ('2x2', 'Calibrated Bayer superpixel sum: R+G1+G2+B')
    header.add_history('EXOTIC: bias/dark/flat and same-colour bad-pixel repair before bin2x2.')
    return header
