"""Detector masks and conservative low-tail rejection before photometry."""

import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import distance_transform_edt, median_filter


def enabled(value):
    return value is True or str(value).strip().lower() in {'y', 'yes', 'true', '1'}


def low_pixel_mask(image, gap_width=17, population_threshold=25):
    """OBS principle: find the first quiet histogram run below the mode.

    Use the whole detector, integer ADU bins, and fewer than 25 pixels per
    bin for 17 consecutive bins. No gap means no rejection. This is a lower
    histogram-tail test, not symmetric clipping of stars and background.
    """
    values = np.asarray(image, dtype=float)
    finite = np.isfinite(values)
    mask = np.zeros(values.shape, dtype=bool)
    if not finite.any():
        return mask, None
    bins, counts = np.unique(np.trunc(values[finite]), return_counts=True)
    mode = bins[counts.argmax()]
    population = dict(zip(bins, counts))
    for offset in range(1, 10000):
        threshold = mode - offset
        if all(population.get(threshold - step, 0) < population_threshold
               for step in range(gap_width)):
            return finite & (values < threshold), float(threshold)
    return mask, None


def robust_sigma(values):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if not values.size:
        return 0.0
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def load_mask(path, shape):
    """FITS mask convention: zero good, nonzero/nonfinite bad."""
    mask = np.asarray(fits.getdata(Path(path).expanduser()))
    if mask.shape != tuple(shape):
        raise ValueError(f'Bad-pixel map shape {mask.shape} does not match image shape {tuple(shape)}')
    return ~np.isfinite(mask) | (mask != 0)


def dark_masks(files, shape, hot_sigma=10.0, hot_floor=50.0, variability_sigma=8.0):
    """Flag extreme dark defects; compare variability within exposure groups.

    Hot/cold defects in any component can contaminate a combined master.
    They are candidates, not a claim of persistence. Variability needs five
    independent frames with matching positive exposures. Streaming moments
    avoid retaining the dark stack in RAM. Frame offsets are removed first.
    """
    hot = np.zeros(shape, bool)
    cold = np.zeros(shape, bool)
    invalid = np.zeros(shape, bool)
    variable = np.zeros(shape, bool)
    groups = {}
    files = list(dict.fromkeys(str(Path(file).resolve()) for file in (files or [])))
    for file in files:
        header = fits.getheader(file)
        exposure = header.get('EXPTIME', header.get('EXPOSURE'))
        try:
            exposure = float(exposure)
        except (TypeError, ValueError):
            exposure = None
        if exposure is not None and (not np.isfinite(exposure) or exposure <= 0):
            exposure = None
        # Never count a combined master as independent dark exposures.
        try:
            combined_count = int(header.get('NCOMBINE', 1))
        except (TypeError, ValueError):
            combined_count = 1
        is_master = 'master' in Path(file).stem.lower() or combined_count > 1
        key = exposure if exposure is not None and not is_master else None
        groups.setdefault(key, []).append(str(file))
    group_reports = []
    for exposure, paths in groups.items():
        assess_variability = exposure is not None and len(paths) >= 5
        count = np.zeros(shape, np.int32) if assess_variability else None
        mean = np.zeros(shape, float) if assess_variability else None
        m2 = np.zeros(shape, float) if assess_variability else None
        for file in paths:
            frame = np.asarray(fits.getdata(file), dtype=float)
            if frame.shape != tuple(shape):
                raise ValueError(f'Dark {file} shape {frame.shape} does not match detector {tuple(shape)}')
            finite = np.isfinite(frame)
            invalid |= ~finite
            if not finite.any():
                raise ValueError(f'Dark {file} has no finite pixels')
            center = float(np.median(frame[finite]))
            working = np.where(finite, frame, center)
            excess = working - median_filter(working, size=5, mode='mirror')
            limit = max(float(hot_floor), float(hot_sigma) * robust_sigma(excess[finite]))
            hot |= finite & (excess > limit)
            low, _ = low_pixel_mask(frame)
            cold |= low
            if assess_variability:
                centered = working - center
                count += finite
                delta = centered - mean
                mean += np.divide(delta, count, out=np.zeros_like(mean), where=finite & (count > 0))
                m2 += np.where(finite, delta * (centered - mean), 0.0)
        report = {'exposure_seconds': exposure, 'frame_count': len(paths),
                  'variability_evaluated': assess_variability}
        if assess_variability:
            usable = count >= 5
            scatter = np.sqrt(np.maximum(np.divide(m2, count - 1,
                out=np.zeros_like(m2), where=usable), 0.0))
            typical = float(np.median(scatter[usable])) if usable.any() else 0.0
            limit = max(3.0 * typical, typical + variability_sigma * robust_sigma(scatter[usable]))
            flagged = usable & (scatter > limit)
            variable |= flagged
            report.update(typical_scatter_adu=typical, variability_threshold_adu=limit,
                          variable_pixel_count=int(flagged.sum()))
        else:
            report['variability_skip_reason'] = 'Need at least 5 individual darks with matching known exposure'
        group_reports.append(report)
    return {'dark_hot': hot, 'dark_low': cold, 'dark_nonfinite': invalid,
            'dark_variable': variable}, group_reports


def prepare_reference(info, shape, existing=None, dark_files=None):
    """Combine configured masks and write auditable detector-coordinate maps."""
    external = info.get('bad_pixel_map')
    use_darks = enabled(info.get('detect_bad_pixels_from_darks', True))
    use_low = enabled(info.get('detect_low_pixels_before_photometry', True))
    dark_source = info.get('bad_pixel_dark_source')
    if use_darks and dark_source:
        if isinstance(dark_source, (list, tuple)):
            dark_files = list(dark_source)
        else:
            source_path = Path(dark_source).expanduser()
            if source_path.is_dir():
                extensions = ('.fits', '.fit', '.fts', '.fz', '.fits.gz', '.fit.gz')
                dark_files = sorted(p for p in source_path.iterdir()
                                    if p.is_file() and p.name.lower().endswith(extensions))
                if not dark_files:
                    raise ValueError(f'No dark FITS frames found for detector mask in {source_path}')
            elif source_path.is_file():
                dark_files = [source_path]
            else:
                raise ValueError(f'Dark source for detector mask does not exist: {source_path}')
    if not external and not use_low and not (use_darks and dark_files):
        return existing
    layers = {}
    if existing is not None:
        old = np.asarray(existing['mask'], dtype=bool)
        if old.shape != tuple(shape):
            raise ValueError('Science-derived bad-pixel mask does not match detector shape')
        layers['science_persistent'] = old
    if external:
        layers['supplied'] = load_mask(external, shape)
    groups = []
    if use_darks and dark_files:
        layers_from_darks, groups = dark_masks(dark_files, shape)
        layers.update(layers_from_darks)
    mask = np.zeros(shape, bool)
    for layer in layers.values():
        mask |= layer
    summary = {'external_map': str(external) if external else None,
               'low_tail_per_frame': use_low, 'low_tail_gap_bins': 17,
               'low_tail_population_threshold': 25,
               'dark_hot_sigma': 10.0, 'dark_hot_floor_adu': 50.0,
               'dark_variability_sigma': 8.0, 'dark_variability_min_frames': 5,
               'dark_groups': groups, 'static_pixel_count': int(mask.sum()),
               'layer_counts': {name: int(layer.sum()) for name, layer in layers.items()},
               'repair': 'NaN mask, median of valid 8 neighbours; nearest valid pixel for enclosed clusters'}
    destination = Path(info['save']) / 'working_artifacts'
    destination.mkdir(parents=True, exist_ok=True)
    for name, layer in layers.items():
        fits.writeto(destination / f'BadPixelMask_{name}.fits', layer.astype(np.uint8), overwrite=True)
    path = destination / 'BadPixelMask.fits'
    fits.writeto(path, mask.astype(np.uint8), overwrite=True)
    (destination / 'BadPixelSummary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    result = dict(existing or {})
    y, x = np.nonzero(mask)
    result.update(mask=mask, coord_y=y, coord_x=x, mask_path=path,
                  detect_low=use_low, summary=summary)
    return result


def repair_frame(image, reference):
    """Mask simultaneously so no defective neighbour is used for repair."""
    if reference is None:
        return image
    repaired = np.array(image, dtype=float, copy=True)
    mask = np.asarray(reference['mask'], dtype=bool).copy()
    if mask.shape != repaired.shape:
        raise ValueError('Bad-pixel mask shape differs from current frame')
    mask |= ~np.isfinite(repaired)
    if reference.get('detect_low'):
        low, _ = low_pixel_mask(repaired)
        mask |= low
    if not mask.any():
        return repaired
    if mask.all():
        raise ValueError('Cannot infill a frame with no valid pixels')
    repaired[mask] = np.nan
    y, x = np.nonzero(mask)
    padded = np.pad(repaired, 1, mode='constant', constant_values=np.nan)
    neighbours = np.stack([padded[y + dy + 1, x + dx + 1]
                           for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dy or dx])
    has_neighbours = np.isfinite(neighbours).any(axis=0)
    repaired[y[has_neighbours], x[has_neighbours]] = np.nanmedian(neighbours[:, has_neighbours], axis=0)
    if not has_neighbours.all():
        # Use original good pixels, not interpolated values, for wider defects.
        nearest = distance_transform_edt(mask, return_distances=False, return_indices=True)
        original = np.asarray(image)
        missing_y, missing_x = y[~has_neighbours], x[~has_neighbours]
        repaired[missing_y, missing_x] = original[tuple(nearest[:, missing_y, missing_x])]
    return repaired
