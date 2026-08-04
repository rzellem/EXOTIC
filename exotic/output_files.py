from json import dump, dumps
import shutil
from numpy import mean, std
from pathlib import Path
import numpy as np

try:
    from utils import (
        aavso_output_directory,
        MAGNITUDE_DECIMAL_PLACES,
        filename_date_token,
        format_magnitude_error,
        format_magnitude,
        magnitude_text,
        normalized_magnitude_error,
        round_to_2,
        rounded_magnitude_error,
        rounded_magnitude_value,
        safe_output_filename,
    )
except ImportError:
    from .utils import (
        aavso_output_directory,
        MAGNITUDE_DECIMAL_PLACES,
        filename_date_token,
        format_magnitude_error,
        format_magnitude,
        magnitude_text,
        normalized_magnitude_error,
        round_to_2,
        rounded_magnitude_error,
        rounded_magnitude_value,
        safe_output_filename,
    )
try:
    from version import __version__
except ImportError:
    from .version import __version__
try:
    from plate_status import PlateStatus
except ImportError:
    from .plate_status import PlateStatus
try:
    from transit_depth import (
        AREA_DEPTH_LABEL,
        OBSERVABLE_DEPTH_DELTA_LABEL,
        OBSERVABLE_DEPTH_LABEL,
        PRIOR_OBSERVABLE_DEPTH_LABEL,
        fit_transit_depth_summary,
        planet_dict_transit_errors,
        planet_dict_transit_parameters,
    )
except ImportError:
    from .transit_depth import (
        AREA_DEPTH_LABEL,
        OBSERVABLE_DEPTH_DELTA_LABEL,
        OBSERVABLE_DEPTH_LABEL,
        PRIOR_OBSERVABLE_DEPTH_LABEL,
        fit_transit_depth_summary,
        planet_dict_transit_errors,
        planet_dict_transit_parameters,
    )


AAVSO_FINDER_STRETCH_NAMES = (
    'LinearStretch',
    'SquaredStretch',
    'SqrtStretch',
    'LogStretch',
)


def copy_aavso_supporting_artifacts(save, target_name, observation_date):
    """Copy final lightcurve, finder, triangle, and QC products into ``AAVSO_Files``."""

    output_dir = Path(save)
    working_artifacts_dir = output_dir / 'working_artifacts'
    diagnostics_dir = output_dir / 'Diagnostics'
    date_token = filename_date_token(observation_date)
    source_paths = [
        output_dir / safe_output_filename(
            'FinalLightCurve', target_name, date_token, extension=extension
        )
        for extension in ('png', 'pdf')
    ]
    source_paths.append(
        working_artifacts_dir / safe_output_filename(
            'FinalLightCurve', target_name, date_token, extension='csv'
        )
    )
    for stretch_name in AAVSO_FINDER_STRETCH_NAMES:
        source_paths.extend(
            working_artifacts_dir / safe_output_filename(
                'FOV', target_name, stretch_name, date_token, extension=extension
            )
            for extension in ('png', 'pdf')
        )
    for prefix, extensions in (
        ('FinalTriangle', ('png',)),
        ('Triangle', ('png',)),
        ('ZoomedTrianglePlot', ('png',)),
        ('KTMF_QC', ('png', 'pdf')),
        ('PriorPosteriorComparison', ('png', 'pdf')),
    ):
        source_paths.extend(
            diagnostics_dir / safe_output_filename(
                prefix, target_name, date_token, extension=extension
            )
            for extension in extensions
        )

    aavso_dir = aavso_output_directory(output_dir)
    copied_paths = []
    for source_path in source_paths:
        if not source_path.is_file():
            continue
        destination_path = aavso_dir / source_path.name
        shutil.copy2(source_path, destination_path)
        copied_paths.append(destination_path)
    return copied_paths


def aavso_airmass_results(fit):
    if getattr(fit, 'airmass_fit_skipped', False):
        return (
            ('Am1', '0', '0'),
            ('Am2', '0', '0'),
        )

    if 'a0' in fit.parameters:
        first_result = (
            'A0',
            str(round_to_2(fit.parameters['a0'], fit.errors['a0'])),
            str(round_to_2(fit.errors['a0'])),
        )
    else:
        first_result = (
            'Am1',
            str(round_to_2(fit.parameters['a1'], fit.errors['a1'])),
            str(round_to_2(fit.errors['a1'])),
        )

    return (
        first_result,
        (
            'Am2',
            str(round_to_2(fit.parameters.get('a2', 0), fit.errors.get('a2', 0))),
            str(round_to_2(fit.errors.get('a2', 0))),
        ),
    )


def aavso_detrend_model(fit):
    if getattr(fit, 'airmass_fit_skipped', False):
        return np.ones(len(fit.time), dtype=float)
    return np.asarray(fit.airmass_model, dtype=float)


def finite_float(value, default=np.nan):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return default
    return value if np.isfinite(value) else default


def apparent_magnitude_calibration_from_vsp_params(vsp_params):
    rows = []
    zero_points = []
    zero_point_errors = []
    for vsp_p in vsp_params or []:
        mag = finite_float(vsp_p.get('mag'))
        mag_err = normalized_magnitude_error(vsp_p.get('mag_err'))
        if not np.isfinite(mag) or mag_err is None:
            continue
        rows.append((mag, mag_err, vsp_p.get('mag_band') or 'V'))

        differential_mag, differential_err = differential_magnitude_from_vsp_param(
            vsp_p
        )
        if np.isfinite(differential_mag):
            comparison_mag = finite_float(vsp_p.get('cmag'))
            comparison_mag_err = normalized_magnitude_error(vsp_p.get('cmag_err'))
            zero_points.append(
                comparison_mag if np.isfinite(comparison_mag) else mag - differential_mag
            )
            if comparison_mag_err is not None:
                zero_point_errors.append(comparison_mag_err)
            elif np.isfinite(differential_err):
                zero_point_errors.append(
                    np.sqrt(max(mag_err ** 2 - differential_err ** 2, 0.0))
                )
            else:
                zero_point_errors.append(mag_err)
            continue

        comparison_mag = finite_float(vsp_p.get('cmag'))
        comparison_mag_err = normalized_magnitude_error(vsp_p.get('cmag_err'))
        if np.isfinite(comparison_mag):
            zero_points.append(comparison_mag)
            zero_point_errors.append(
                comparison_mag_err if comparison_mag_err is not None else mag_err
            )

    if not rows:
        return None

    magnitudes = np.array([row[0] for row in rows], dtype=float)
    magnitude_errors = np.array([row[1] for row in rows], dtype=float)
    return {
        'baseline_magnitude': float(np.nanmedian(magnitudes)),
        'baseline_error': float(np.nanmedian(magnitude_errors)),
        'band': rows[0][2],
        'zero_point_magnitude': (
            float(np.nanmedian(zero_points)) if zero_points else np.nan
        ),
        'zero_point_error': (
            float(np.nanmedian(zero_point_errors)) if zero_point_errors else np.nan
        ),
    }


def differential_magnitude_from_vsp_param(vsp_param):
    """Return target-minus-reference magnitude and its flux-only uncertainty."""
    differential_mag = finite_float(
        vsp_param.get(
            'differential_mag',
            vsp_param.get('differential_magnitude'),
        )
    )
    differential_err = finite_float(
        vsp_param.get(
            'differential_mag_err',
            vsp_param.get('differential_magnitude_error'),
        )
    )
    if np.isfinite(differential_err) and differential_err < 0:
        differential_err = np.nan

    if np.isfinite(differential_mag):
        return differential_mag, differential_err

    apparent_mag = finite_float(vsp_param.get('mag'))
    comparison_mag = finite_float(vsp_param.get('cmag'))
    if not (np.isfinite(apparent_mag) and np.isfinite(comparison_mag)):
        return np.nan, np.nan

    apparent_err = normalized_magnitude_error(vsp_param.get('mag_err'))
    comparison_err = normalized_magnitude_error(vsp_param.get('cmag_err'))
    if apparent_err is not None and comparison_err is not None:
        differential_err = np.sqrt(max(apparent_err ** 2 - comparison_err ** 2, 0.0))
    elif apparent_err is not None:
        differential_err = apparent_err
    else:
        differential_err = np.nan
    return apparent_mag - comparison_mag, differential_err


def differential_magnitude_series_from_fit(fit, out_of_transit_only=False,
                                           apply_airmass_correction=None):
    """Return target-minus-reference instrumental magnitudes.

    Unlike apparent magnitudes, this series needs no catalogue magnitude.  Raw
    target/reference fluxes are preferred so the instrumental zero point is
    preserved; pre-reduced light curves fall back to their relative flux.
    Stellar-variability fits intentionally remain uncorrected for airmass so a
    real time-dependent stellar signal is not fitted away.
    """
    fit_data = np.asarray(
        getattr(fit, 'data', getattr(fit, 'detrended', [])),
        dtype=float,
    ).reshape(-1)
    fit_times = np.asarray(
        getattr(fit, 'time', getattr(fit, 'jd_times', [])),
        dtype=float,
    ).reshape(-1)
    if fit_data.size == 0 or fit_times.shape != fit_data.shape:
        return None

    target_flux = np.asarray(
        getattr(
            fit,
            'differential_magnitude_target_flux',
            getattr(fit, 'stellar_variability_target_flux', []),
        ),
        dtype=float,
    ).reshape(-1)
    reference_flux = np.asarray(
        getattr(
            fit,
            'differential_magnitude_reference_flux',
            getattr(fit, 'stellar_variability_comp_flux', []),
        ),
        dtype=float,
    ).reshape(-1)
    target_error = np.asarray(
        getattr(
            fit,
            'differential_magnitude_target_flux_error',
            getattr(fit, 'stellar_variability_target_flux_error', []),
        ),
        dtype=float,
    ).reshape(-1)
    reference_error = np.asarray(
        getattr(
            fit,
            'differential_magnitude_reference_flux_error',
            getattr(fit, 'stellar_variability_comp_flux_error', []),
        ),
        dtype=float,
    ).reshape(-1)

    has_raw_photometry = (
        target_flux.shape == fit_data.shape
        and reference_flux.shape == fit_data.shape
    )
    if not has_raw_photometry:
        target_flux = np.asarray(getattr(fit, 'detrended', fit_data), dtype=float).reshape(-1)
        if target_flux.shape != fit_data.shape:
            target_flux = fit_data.copy()
        reference_flux = np.ones(fit_data.shape, dtype=float)
        target_error = np.asarray(
            getattr(fit, 'detrendederr', getattr(fit, 'dataerr', [])),
            dtype=float,
        ).reshape(-1)
        reference_error = np.zeros(fit_data.shape, dtype=float)

    if target_error.shape != fit_data.shape:
        target_error = np.full(fit_data.shape, np.nan, dtype=float)
    if reference_error.shape != fit_data.shape:
        reference_error = np.full(fit_data.shape, np.nan, dtype=float)

    if apply_airmass_correction is None:
        apply_airmass_correction = not bool(
            getattr(fit, 'stellar_variability_only', False)
        )
    # ``fit.detrended`` is already corrected.  Only divide an explicitly
    # retained raw target/reference ratio by the fitted airmass model.
    apply_airmass_correction = bool(apply_airmass_correction and has_raw_photometry)
    relative_airmass_model = np.ones(fit_data.shape, dtype=float)
    if apply_airmass_correction:
        airmass_model = np.asarray(
            getattr(fit, 'airmass_model', np.ones(fit_data.shape)),
            dtype=float,
        ).reshape(-1)
        if airmass_model.shape != fit_data.shape:
            airmass_model = np.ones(fit_data.shape, dtype=float)
        valid_airmass_model = np.isfinite(airmass_model) & (airmass_model > 0)
        airmass_reference = (
            float(np.nanmedian(airmass_model[valid_airmass_model]))
            if np.any(valid_airmass_model)
            else 1.0
        )
        if not np.isfinite(airmass_reference) or airmass_reference <= 0:
            airmass_reference = 1.0
        relative_airmass_model = np.divide(
            airmass_model,
            airmass_reference,
            out=np.full(fit_data.shape, np.nan, dtype=float),
            where=valid_airmass_model,
        )

    with np.errstate(divide='ignore', invalid='ignore'):
        raw_ratio = np.divide(target_flux, reference_flux)
        corrected_ratio = np.divide(raw_ratio, relative_airmass_model)
        differential_magnitude = -2.5 * np.log10(corrected_ratio)
        magnitude_factor = 2.5 / np.log(10.0)
        explicit_error = magnitude_factor * np.sqrt(
            (target_error / target_flux) ** 2
            + (reference_error / reference_flux) ** 2
        )

    fit_error = np.asarray(getattr(fit, 'dataerr', []), dtype=float).reshape(-1)
    if fit_error.shape == fit_data.shape:
        with np.errstate(divide='ignore', invalid='ignore'):
            fallback_error = magnitude_factor * np.abs(fit_error / fit_data)
    else:
        fallback_error = np.full(fit_data.shape, np.nan, dtype=float)
    differential_error = np.where(
        np.isfinite(explicit_error) & (explicit_error >= 0),
        explicit_error,
        fallback_error,
    )

    airmass = np.asarray(
        getattr(fit, 'airmass', np.full(fit_data.shape, np.nan)),
        dtype=float,
    ).reshape(-1)
    if airmass.shape != fit_data.shape:
        airmass = np.full(fit_data.shape, np.nan, dtype=float)

    keep = (
        np.isfinite(fit_times)
        & np.isfinite(corrected_ratio)
        & (corrected_ratio > 0)
        & np.isfinite(differential_magnitude)
    )
    if out_of_transit_only:
        transit_model = np.asarray(
            getattr(fit, 'transit', np.ones(fit_data.shape)),
            dtype=float,
        ).reshape(-1)
        if transit_model.shape == fit_data.shape and np.any(transit_model == 1):
            keep &= transit_model == 1

    if not np.any(keep):
        return None
    return {
        'time': fit_times[keep],
        'airmass': airmass[keep],
        'magnitude': differential_magnitude[keep],
        'magnitude_error': differential_error[keep],
        'source_mask': keep,
        'airmass_corrected': apply_airmass_correction,
        'has_raw_photometry': has_raw_photometry,
    }


def magnitude_series_from_fit(fit, out_of_transit_only=False,
                              apply_airmass_correction=None):
    """Return aligned differential and, when calibrated, apparent magnitudes."""
    fit_data = np.asarray(
        getattr(fit, 'data', getattr(fit, 'detrended', [])),
        dtype=float,
    ).reshape(-1)
    result = {
        'differential_magnitude': np.full(fit_data.shape, np.nan, dtype=float),
        'differential_magnitude_error': np.full(fit_data.shape, np.nan, dtype=float),
        'apparent_magnitude': np.full(fit_data.shape, np.nan, dtype=float),
        'apparent_magnitude_error': np.full(fit_data.shape, np.nan, dtype=float),
        'band': None,
        'airmass_corrected': False,
        'has_raw_photometry': False,
        'apparent_calibrated': False,
    }
    series = differential_magnitude_series_from_fit(
        fit,
        out_of_transit_only=out_of_transit_only,
        apply_airmass_correction=apply_airmass_correction,
    )
    if series is None:
        return result

    source_mask = np.asarray(series['source_mask'], dtype=bool)
    if source_mask.shape != fit_data.shape:
        return result
    result['differential_magnitude'][source_mask] = series['magnitude']
    result['differential_magnitude_error'][source_mask] = series['magnitude_error']
    result['airmass_corrected'] = bool(series['airmass_corrected'])
    result['has_raw_photometry'] = bool(series['has_raw_photometry'])

    calibration = apparent_magnitude_calibration_from_vsp_params(
        getattr(fit, 'stellar_variability_params', None)
    )
    if calibration is None:
        return result

    # A calibrated comparison ensemble already carries its independently
    # derived apparent-magnitude series.  Do not reconstruct that series by
    # adding a constant to the raw instrumental differential magnitudes: the
    # calibrated ensemble and the raw median-scaled ensemble intentionally use
    # different reference constructions and can have different time trends.
    calibrated_magnitude = np.asarray(
        getattr(fit, 'stellar_variability_ensemble_magnitudes', []),
        dtype=float,
    ).reshape(-1)
    calibrated_error = np.asarray(
        getattr(fit, 'stellar_variability_ensemble_magnitude_errors', []),
        dtype=float,
    ).reshape(-1)
    if calibrated_magnitude.shape == fit_data.shape:
        calibrated_mask = source_mask & np.isfinite(calibrated_magnitude)
        result['apparent_magnitude'][calibrated_mask] = calibrated_magnitude[calibrated_mask]
        if calibrated_error.shape == fit_data.shape:
            calibrated_error_mask = calibrated_mask & np.isfinite(calibrated_error)
            result['apparent_magnitude_error'][calibrated_error_mask] = (
                calibrated_error[calibrated_error_mask]
            )
        result['band'] = calibration['band']
        result['apparent_calibrated'] = bool(np.any(calibrated_mask))
        return result

    magnitude_offset = np.nan
    calibration_error = calibration['baseline_error']
    if np.isfinite(calibration['zero_point_error']):
        calibration_error = calibration['zero_point_error']

    differential_magnitude = result['differential_magnitude']
    differential_error = result['differential_magnitude_error']
    fit_times = np.asarray(
        getattr(fit, 'time', getattr(fit, 'jd_times', [])),
        dtype=float,
    ).reshape(-1)
    matched_offsets = []
    if result['has_raw_photometry'] and fit_times.shape == differential_magnitude.shape:
        for vsp_param in getattr(fit, 'stellar_variability_params', None) or []:
            apparent_mag = finite_float(vsp_param.get('mag'))
            apparent_time = finite_float(vsp_param.get('time'))
            if not (np.isfinite(apparent_mag) and np.isfinite(apparent_time)):
                continue
            matches = np.flatnonzero(np.isclose(
                fit_times,
                apparent_time,
                rtol=0.0,
                atol=1.0e-7,
            ))
            if matches.size == 0:
                continue
            matched_differential = differential_magnitude[matches[0]]
            if np.isfinite(matched_differential):
                matched_offsets.append(apparent_mag - matched_differential)
    if not result['has_raw_photometry']:
        magnitude_offset = calibration['baseline_magnitude']
    elif matched_offsets:
        magnitude_offset = float(np.nanmedian(matched_offsets))
    else:
        reference_mask = np.isfinite(differential_magnitude)
        transit_model = np.asarray(
            getattr(fit, 'transit', np.ones(fit_data.shape)),
            dtype=float,
        ).reshape(-1)
        if transit_model.shape == fit_data.shape and np.any(transit_model == 1):
            reference_mask &= transit_model == 1
        if np.any(reference_mask):
            magnitude_offset = (
                calibration['baseline_magnitude']
                - float(np.nanmedian(differential_magnitude[reference_mask]))
            )
    if not np.isfinite(magnitude_offset):
        return result

    finite_differential = np.isfinite(differential_magnitude)
    result['apparent_magnitude'][finite_differential] = (
        differential_magnitude[finite_differential] + magnitude_offset
    )
    if np.isfinite(calibration_error):
        finite_error = finite_differential & np.isfinite(differential_error)
        result['apparent_magnitude_error'][finite_error] = np.hypot(
            differential_error[finite_error],
            calibration_error,
        )
        missing_error = finite_differential & ~np.isfinite(differential_error)
        result['apparent_magnitude_error'][missing_error] = calibration_error
    else:
        result['apparent_magnitude_error'][finite_differential] = (
            differential_error[finite_differential]
        )
    result['band'] = calibration['band']
    result['apparent_calibrated'] = True
    return result


def write_differential_magnitude_csv(fit, save, target_name, observation_date=None,
                                     observed_filter=None, out_of_transit_only=False,
                                     apply_airmass_correction=None,
                                     filename_prefix='DifferentialMagnitude'):
    series = differential_magnitude_series_from_fit(
        fit,
        out_of_transit_only=out_of_transit_only,
        apply_airmass_correction=apply_airmass_correction,
    )
    if series is None:
        return None

    output_dir = Path(save)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / safe_output_filename(
        filename_prefix,
        target_name,
        filename_date_token(observation_date) if observation_date else 'undated',
        extension='csv',
    )
    comparison = (
        getattr(fit, 'stellar_variability_reference_label', None)
        or getattr(fit, 'differential_magnitude_reference_label', None)
        or 'selected comparison reference'
    )
    with output_path.open('w', encoding='utf-8') as handle:
        handle.write(
            '# AIRMASS_CORRECTION='
            f"{'YES' if series['airmass_corrected'] else 'NO'}\n"
        )
        handle.write(
            '# BJD_TDB,Airmass,Differential Magnitude,'
            'Differential Magnitude Uncertainty,Filter,Comparison\n'
        )
        for time_value, airmass, magnitude, magnitude_error in zip(
            series['time'],
            series['airmass'],
            series['magnitude'],
            series['magnitude_error'],
        ):
            airmass_text = f"{airmass}" if np.isfinite(airmass) else 'na'
            error_text = (
                f"{magnitude_error:.{MAGNITUDE_DECIMAL_PLACES}f}"
                if np.isfinite(magnitude_error)
                else 'na'
            )
            handle.write(
                f"{time_value}, {airmass_text}, "
                f"{magnitude:.{MAGNITUDE_DECIMAL_PLACES}f}, {error_text}, "
                f"{observed_filter or 'na'}, {comparison}\n"
            )
    return output_path


def aavso_json_safe(value):
    if isinstance(value, dict):
        return {str(key): aavso_json_safe(subvalue) for key, subvalue in value.items()}
    if isinstance(value, (list, tuple)):
        return [aavso_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return aavso_json_safe(value.item())
        return [aavso_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return aavso_json_safe(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, float):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, int):
        return int(value)
    return value


def format_optional_float(value, digits=7):
    value = finite_float(value)
    if np.isfinite(value):
        return f"{value:.{digits}f}"
    return "na"


def stellar_variability_reference_summary(vsp_param):
    if not vsp_param:
        return "na"

    cname = vsp_param.get('cname', 'na')
    band = vsp_param.get('catalog_mag_band') or vsp_param.get('mag_band') or 'V'
    if vsp_param.get('ensemble_reference'):
        member_count = int(vsp_param.get('ensemble_member_count', 0) or 0)
        labels = vsp_param.get('ensemble_member_labels') or []
        label_text = ", ".join(str(label) for label in labels)
        return (
            f"Calibrated comparison-star ensemble ({member_count} stars)"
            + (f": {label_text}" if label_text else "")
        )
    if vsp_param.get('is_aavso_vsp', True):
        return f"AAVSO Label: {cname}, Position: {vsp_param.get('pos')}"

    source = vsp_param.get('catalog_source') or 'NextAstro photometry catalog'
    comp_ra = format_optional_float(vsp_param.get('comp_ra'))
    comp_dec = format_optional_float(vsp_param.get('comp_dec'))
    details = [f"{source}: RA={comp_ra}", f"Dec={comp_dec}"]
    mag_text = magnitude_text(band, vsp_param.get('cmag'), vsp_param.get('cmag_err'))
    if mag_text is not None:
        details.append(mag_text)
    return ", ".join(details)


def stellar_variability_measurement_summary(vsp_params, transit_fit_comp_star=None):
    point_count = len(vsp_params or [])
    if point_count == 0:
        return None

    if transit_fit_comp_star is None:
        return None

    if vsp_params[0].get('ensemble_reference'):
        return (
            f"Combined {point_count} out-of-transit target measurements against the calibrated "
            "comparison-star ensemble; AID rows list the BJD_TDB timestamps used."
        )

    return (
        f"Remeasured {point_count} out-of-transit target/reference point(s) against the transit-fit "
        f"{'derived ' if vsp_params[0].get('derived_catalog_reference') else ''}catalog reference "
        "for AID magnitudes; AID rows list the BJD_TDB timestamps used."
    )


def aid_comparison_metadata(vsp_param):
    if not vsp_param:
        return {}
    anchor_labels = vsp_param.get('derived_reference_anchor_labels')
    anchor_label_sample = None
    if isinstance(anchor_labels, np.ndarray):
        anchor_labels = anchor_labels.tolist()
    if isinstance(anchor_labels, (list, tuple)):
        anchor_labels = list(anchor_labels)
        if len(anchor_labels) > 10:
            anchor_label_sample = anchor_labels[:10]
            anchor_labels = None
    metadata = {
        'source': vsp_param.get('catalog_source', 'AAVSO VSP'),
        'is_aavso_vsp': bool(vsp_param.get('is_aavso_vsp', True)),
        'comparison_name': vsp_param.get('cname'),
        'comparison_position_pixels': vsp_param.get('pos'),
        'comparison_ra_deg': vsp_param.get('comp_ra'),
        'comparison_dec_deg': vsp_param.get('comp_dec'),
        'catalog_ra_deg': vsp_param.get('catalog_ra'),
        'catalog_dec_deg': vsp_param.get('catalog_dec'),
        'catalog_source_id': vsp_param.get('source_id'),
        'catalog_id': vsp_param.get('catalog_id'),
        'catalog_match_separation_arcsec': vsp_param.get('separation_arcsec'),
        'derived_catalog_reference': bool(vsp_param.get('derived_catalog_reference', False)),
        'derived_reference_anchor_count': vsp_param.get('derived_reference_anchor_count'),
        'magnitude_band': vsp_param.get('catalog_mag_band') or vsp_param.get('mag_band'),
        'reported_measurement_band': vsp_param.get('mag_band'),
        'apparent_magnitude': rounded_magnitude_value(vsp_param.get('cmag')),
        'apparent_magnitude_error': rounded_magnitude_error(vsp_param.get('cmag_err')),
        'ensemble_reference': bool(vsp_param.get('ensemble_reference', False)),
        'ensemble_member_count': vsp_param.get('ensemble_member_count'),
        'ensemble_member_labels': vsp_param.get('ensemble_member_labels'),
        'ensemble_member_positions': vsp_param.get('ensemble_member_positions'),
        'ensemble_member_catalog_magnitudes': vsp_param.get('ensemble_member_catalog_magnitudes'),
        'ensemble_member_catalog_errors': vsp_param.get('ensemble_member_catalog_errors'),
        'ensemble_member_catalog_sources': vsp_param.get('ensemble_member_catalog_sources'),
        'ensemble_member_ra_degs': vsp_param.get('ensemble_member_ra_degs'),
        'ensemble_member_dec_degs': vsp_param.get('ensemble_member_dec_degs'),
        'ensemble_members': vsp_param.get('ensemble_members'),
        'ensemble_member_catalog_colors': vsp_param.get('ensemble_member_catalog_colors'),
        'ensemble_member_catalog_color_labels': vsp_param.get('ensemble_member_catalog_color_labels'),
        'ensemble_member_color_deltas': vsp_param.get('ensemble_member_color_deltas'),
        'ensemble_member_magnitude_deltas': vsp_param.get('ensemble_member_magnitude_deltas'),
        'ensemble_member_similarity_scores': vsp_param.get('ensemble_member_similarity_scores'),
    }
    if anchor_labels is not None:
        metadata['derived_reference_anchor_labels'] = anchor_labels
    if anchor_label_sample is not None:
        metadata['derived_reference_anchor_label_sample'] = anchor_label_sample
    return aavso_json_safe(metadata)


def aid_comparison_coordinate_headers(vsp_params, indexed=False):
    """Return standards-safe comparison coordinates with RA and Dec on separate lines."""
    if isinstance(vsp_params, dict):
        vsp_params = [vsp_params]

    coordinates = []
    seen = set()
    for vsp_param in vsp_params or []:
        if not isinstance(vsp_param, dict) or vsp_param.get('ensemble_reference'):
            continue
        comp_ra = finite_float(vsp_param.get('comp_ra'))
        comp_dec = finite_float(vsp_param.get('comp_dec'))
        if not np.isfinite(comp_ra) or not np.isfinite(comp_dec):
            continue
        comparison_name = format_aavso_header_value(vsp_param.get('cname'))
        identity = (comparison_name, round(float(comp_ra), 10), round(float(comp_dec), 10))
        if identity in seen:
            continue
        seen.add(identity)
        coordinates.append((comparison_name, float(comp_ra), float(comp_dec)))

    if not coordinates:
        return ""
    if not indexed and len(coordinates) == 1:
        _, comp_ra, comp_dec = coordinates[0]
        return f"#COMPARISON_RA={comp_ra:.7f}\n#COMPARISON_DEC={comp_dec:.7f}\n"

    headers = []
    for index, (comparison_name, comp_ra, comp_dec) in enumerate(coordinates, start=1):
        if comparison_name:
            headers.append(f"#COMPARISON_{index}_NAME={comparison_name}")
        headers.append(f"#COMPARISON_{index}_RA={comp_ra:.7f}")
        headers.append(f"#COMPARISON_{index}_DEC={comp_dec:.7f}")
    return "\n".join(headers) + "\n"


def aid_ensemble_comparison_metadata(vsp_param):
    if not vsp_param or not vsp_param.get('ensemble_reference'):
        return {}

    members = vsp_param.get('ensemble_members')
    if isinstance(members, np.ndarray):
        members = members.tolist()
    if not isinstance(members, (list, tuple)) or not members:
        labels = list(vsp_param.get('ensemble_member_labels') or [])
        positions = list(vsp_param.get('ensemble_member_positions') or [])
        ra_degs = list(vsp_param.get('ensemble_member_ra_degs') or [])
        dec_degs = list(vsp_param.get('ensemble_member_dec_degs') or [])
        magnitudes = list(vsp_param.get('ensemble_member_catalog_magnitudes') or [])
        magnitude_errors = list(vsp_param.get('ensemble_member_catalog_errors') or [])
        catalog_sources = list(vsp_param.get('ensemble_member_catalog_sources') or [])
        member_count = max(
            int(vsp_param.get('ensemble_member_count', 0) or 0),
            len(labels),
            len(ra_degs),
            len(dec_degs),
        )

        def value_at(values, index):
            return values[index] if index < len(values) else None

        members = [
            {
                'label': value_at(labels, index),
                'ra_deg': value_at(ra_degs, index),
                'dec_deg': value_at(dec_degs, index),
                'pixel_position': value_at(positions, index),
                'catalog_magnitude': value_at(magnitudes, index),
                'catalog_magnitude_error': value_at(magnitude_errors, index),
                'catalog_source': value_at(catalog_sources, index),
            }
            for index in range(member_count)
        ]
    else:
        members = list(members)

    return prune_aavso_metadata({
        'member_count': int(vsp_param.get('ensemble_member_count', len(members)) or len(members)),
        'members': members,
    })


def prune_aavso_metadata(value):
    if isinstance(value, dict):
        pruned = {}
        for key, subvalue in value.items():
            cleaned = prune_aavso_metadata(subvalue)
            if cleaned is None or cleaned == "" or cleaned == [] or cleaned == {}:
                continue
            pruned[key] = cleaned
        return pruned
    if isinstance(value, np.ndarray):
        return prune_aavso_metadata(value.tolist())
    if isinstance(value, (list, tuple)):
        return [
            cleaned for cleaned in (prune_aavso_metadata(item) for item in value)
            if cleaned is not None and cleaned != "" and cleaned != [] and cleaned != {}
        ]
    return aavso_json_safe(value)


def format_aavso_json_header(name, payload, preserve_nulls=False):
    payload = aavso_json_safe(payload) if preserve_nulls else prune_aavso_metadata(payload)
    if not payload:
        return ""
    return f"#{name}={dumps(payload, sort_keys=True)}\n"


def aavso_result_entry(value, uncertainty=None, units=None):
    value = finite_float(value)
    uncertainty = finite_float(uncertainty)
    if not np.isfinite(value):
        return None

    entry = {
        'value': str(round_to_2(value, uncertainty)) if np.isfinite(uncertainty) else str(round_to_2(value)),
    }
    if np.isfinite(uncertainty):
        entry['uncertainty'] = str(round_to_2(uncertainty))
    if units:
        entry['units'] = units
    return entry


def numeric_series_summary(values):
    if values is None:
        return {}

    try:
        series = np.asarray(values, dtype=float).reshape(-1)
    except (TypeError, ValueError):
        return {}

    finite = series[np.isfinite(series)]
    if finite.size == 0:
        return {}

    return {
        'count': int(finite.size),
        'median': float(np.nanmedian(finite)),
        'std': float(np.nanstd(finite)),
        'min': float(np.nanmin(finite)),
        'max': float(np.nanmax(finite)),
    }


def path_name(value):
    if value is None:
        return None
    return Path(str(value)).name


def file_list_summary(files, limit=10):
    files = list(files or [])
    return {
        'count': len(files),
        'files': [path_name(file_name) for file_name in files[:limit]],
        'omitted_file_count': max(0, len(files) - limit),
    }


def residual_scatter_fraction(fit):
    transit_qc = getattr(fit, 'transit_qc', None)
    if isinstance(transit_qc, dict):
        qc_residual_scatter = finite_float(transit_qc.get('residual_scatter'))
        if np.isfinite(qc_residual_scatter):
            return qc_residual_scatter

    residuals = np.asarray(getattr(fit, 'residuals', np.array([])), dtype=float)
    data = np.asarray(getattr(fit, 'data', np.array([])), dtype=float)
    if residuals.size == 0 or data.size == 0:
        return np.nan

    median_flux = np.nanmedian(data)
    if not np.isfinite(median_flux) or median_flux == 0:
        return np.nan

    if residuals.shape == data.shape:
        return float(np.nanstd(residuals) / median_flux)
    if residuals.size == 1:
        return float(abs(residuals.reshape(-1)[0]) / median_flux)
    return np.nan


def fit_data_model_uncertainty(fit):
    data = np.asarray(getattr(fit, 'data', np.array([])), dtype=float)
    if data.ndim != 1 or data.size == 0:
        return None, None, None

    model = getattr(fit, 'model', None)
    if model is None:
        residuals = np.asarray(getattr(fit, 'residuals', np.array([])), dtype=float)
        if residuals.shape == data.shape:
            model = data - residuals
    if model is None:
        transit_model = getattr(fit, 'transit', None)
        systematics_model = getattr(fit, 'airmass_model', None)
        if transit_model is not None and systematics_model is not None:
            model = np.asarray(transit_model, dtype=float) * np.asarray(systematics_model, dtype=float)
    if model is None:
        return data, None, None

    model = np.asarray(model, dtype=float)
    if model.shape != data.shape:
        return data, None, None

    uncertainty = getattr(fit, 'dataerr', None)
    if uncertainty is not None:
        uncertainty = np.asarray(uncertainty, dtype=float)
        if uncertainty.shape != data.shape:
            uncertainty = None

    return data, model, uncertainty


def infer_fit_quality_parameter_count(fit):
    transit_qc = getattr(fit, 'transit_qc', None)
    if isinstance(transit_qc, dict):
        parameter_count = finite_float(transit_qc.get('transit_parameter_count'))
        if np.isfinite(parameter_count) and parameter_count > 0:
            return int(parameter_count)

    bounds = getattr(fit, 'bounds', None)
    if isinstance(bounds, dict) and bounds:
        return len(bounds)

    parameters = getattr(fit, 'parameters', None)
    if isinstance(parameters, dict) and parameters:
        return len(parameters)

    return 0


def build_fit_quality_metadata(fit):
    data, model, uncertainty = fit_data_model_uncertainty(fit)
    if data is None or model is None:
        return {}

    residuals = data - model
    finite_mask = np.isfinite(data) & np.isfinite(model) & np.isfinite(residuals)
    if not np.any(finite_mask):
        return {}

    finite_residuals = residuals[finite_mask]
    median_flux = np.nanmedian(data[finite_mask])
    rms_residual = float(np.sqrt(np.nanmean(finite_residuals ** 2)))
    mad_residual = float(np.nanmedian(np.abs(finite_residuals)))
    residual_scatter = (
        float(np.nanstd(finite_residuals) / median_flux)
        if np.isfinite(median_flux) and median_flux != 0
        else np.nan
    )
    point_count = int(np.count_nonzero(finite_mask))
    parameter_count = infer_fit_quality_parameter_count(fit)
    degrees_of_freedom = point_count - parameter_count

    payload = {
        'point_count': point_count,
        'parameter_count': parameter_count,
        'degrees_of_freedom': degrees_of_freedom,
        'rms_residual': rms_residual,
        'rms_residual_percent': (
            100.0 * rms_residual / median_flux
            if np.isfinite(median_flux) and median_flux != 0
            else np.nan
        ),
        'median_absolute_residual': mad_residual,
        'median_flux': float(median_flux) if np.isfinite(median_flux) else np.nan,
        'residual_scatter': residual_scatter,
        'residual_scatter_percent': 100.0 * residual_scatter if np.isfinite(residual_scatter) else np.nan,
        'uses_uncertainties': False,
    }

    if uncertainty is None:
        return payload

    uncertainty_mask = finite_mask & np.isfinite(uncertainty) & (uncertainty > 0)
    if not np.any(uncertainty_mask):
        return payload

    weighted_residuals = residuals[uncertainty_mask]
    weighted_uncertainties = uncertainty[uncertainty_mask]
    normalized_residuals = weighted_residuals / weighted_uncertainties
    chi_square = float(np.sum(normalized_residuals ** 2))
    weighted_point_count = int(np.count_nonzero(uncertainty_mask))
    weighted_degrees_of_freedom = weighted_point_count - parameter_count
    median_uncertainty = float(np.nanmedian(weighted_uncertainties))

    payload.update({
        'uses_uncertainties': True,
        'weighted_point_count': weighted_point_count,
        'degrees_of_freedom': weighted_degrees_of_freedom,
        'chi_square': chi_square,
        'reduced_chi_square': (
            chi_square / weighted_degrees_of_freedom
            if weighted_degrees_of_freedom > 0
            else np.nan
        ),
        'rms_normalized_residual': float(np.sqrt(np.nanmean(normalized_residuals ** 2))),
        'median_absolute_normalized_residual': float(np.nanmedian(np.abs(normalized_residuals))),
        'max_absolute_normalized_residual': float(np.nanmax(np.abs(normalized_residuals))),
        'median_uncertainty': median_uncertainty,
        'rms_residual_to_median_uncertainty': (
            rms_residual / median_uncertainty
            if np.isfinite(median_uncertainty) and median_uncertainty > 0
            else np.nan
        ),
    })
    return payload


def red_noise_beta_factor(residual_fraction, coordinates=None, min_bin_size=2, max_bin_size=None):
    residual_fraction = np.asarray(residual_fraction, dtype=float)
    finite_mask = np.isfinite(residual_fraction)

    coordinates_array = None
    if coordinates is not None:
        coordinates_array = np.asarray(coordinates, dtype=float)
        if coordinates_array.shape == residual_fraction.shape:
            finite_mask &= np.isfinite(coordinates_array)
        else:
            coordinates_array = None

    residual_fraction = residual_fraction[finite_mask]
    if coordinates_array is not None:
        coordinates_array = coordinates_array[finite_mask]
        sort_index = np.argsort(coordinates_array)
        residual_fraction = residual_fraction[sort_index]

    point_count = int(residual_fraction.size)
    payload = {
        'factor': 1.0,
        'point_count': point_count,
        'bin_sizes': [],
        'beta_by_bin': {},
        'max_bin_size': np.nan,
    }
    min_bin_size = int(max(2, finite_float(min_bin_size, 2)))
    if point_count < min_bin_size * 2:
        return payload

    residual_fraction = residual_fraction - np.nanmedian(residual_fraction)
    unbinned_rms = finite_float(np.nanstd(residual_fraction, ddof=1))
    if not np.isfinite(unbinned_rms) or unbinned_rms <= 0:
        return payload

    if max_bin_size is None:
        max_bin_size = min(10, max(min_bin_size, point_count // 4))
    max_bin_size = int(max(min_bin_size, finite_float(max_bin_size, min_bin_size)))
    max_bin_size = min(max_bin_size, point_count // 2)
    if max_bin_size < min_bin_size:
        return payload

    beta_values = []
    for bin_size in range(min_bin_size, max_bin_size + 1):
        bin_count = point_count // bin_size
        if bin_count < 2:
            continue
        trimmed = residual_fraction[:bin_count * bin_size]
        binned_means = np.nanmean(trimmed.reshape(bin_count, bin_size), axis=1)
        binned_rms = finite_float(np.nanstd(binned_means, ddof=1))
        expected_rms = (
            unbinned_rms
            / np.sqrt(bin_size)
            * np.sqrt(bin_count / (bin_count - 1.0))
        )
        if not np.isfinite(binned_rms) or not np.isfinite(expected_rms) or expected_rms <= 0:
            continue
        beta = max(1.0, float(binned_rms / expected_rms))
        payload['bin_sizes'].append(bin_size)
        payload['beta_by_bin'][str(bin_size)] = beta
        beta_values.append(beta)

    if beta_values:
        payload['factor'] = float(max(beta_values))
        payload['max_bin_size'] = max(payload['bin_sizes'])
    return payload


def _fit_uncertainty_time_coordinates(fit, expected_shape):
    for name in ('time', 'phase'):
        values = getattr(fit, name, None)
        if values is None:
            continue
        try:
            values = np.asarray(values, dtype=float)
        except (TypeError, ValueError):
            continue
        if values.shape == expected_shape:
            return values
    return None


def fit_empirical_transit_uncertainty(fit, fit_quality=None, transit_depth_threshold_fraction=0.05):
    data, model, _ = fit_data_model_uncertainty(fit)
    if data is None or model is None:
        return {}

    data = np.asarray(data, dtype=float)
    model = np.asarray(model, dtype=float)
    if data.shape != model.shape or data.ndim != 1:
        return {}

    residuals = data - model
    finite_mask = np.isfinite(data) & np.isfinite(model) & np.isfinite(residuals)
    if not np.any(finite_mask):
        return {}

    fit_quality = fit_quality or {}
    median_flux = np.nanmedian(data[finite_mask])
    residual_scatter = finite_float(fit_quality.get('residual_scatter'))
    if not np.isfinite(residual_scatter):
        if np.isfinite(median_flux) and median_flux != 0:
            residual_scatter = float(np.nanstd(residuals[finite_mask]) / median_flux)
    if not np.isfinite(residual_scatter) or residual_scatter < 0:
        return {}

    transit = np.asarray(getattr(fit, 'transit', np.array([])), dtype=float)
    if transit.shape != data.shape:
        return {}

    transit_mask = finite_mask & np.isfinite(transit)
    if np.count_nonzero(transit_mask) < 2:
        return {}

    transit_values = transit[transit_mask]
    baseline = finite_float(np.nanpercentile(transit_values, 95))
    if not np.isfinite(baseline):
        return {}

    transit_depth_profile = baseline - transit
    max_depth = finite_float(np.nanmax(transit_depth_profile[transit_mask]))
    if not np.isfinite(max_depth) or max_depth <= 0:
        return {}

    threshold_fraction = finite_float(transit_depth_threshold_fraction, 0.05)
    if not np.isfinite(threshold_fraction) or threshold_fraction <= 0:
        threshold_fraction = 0.05
    depth_threshold = max_depth * threshold_fraction
    in_transit_mask = transit_mask & (transit_depth_profile >= depth_threshold)
    out_of_transit_mask = transit_mask & (transit_depth_profile < depth_threshold)

    in_transit_count = int(np.count_nonzero(in_transit_mask))
    out_of_transit_count = int(np.count_nonzero(out_of_transit_mask))
    if in_transit_count <= 0:
        return {}

    sample_term = 1.0 / in_transit_count
    if out_of_transit_count > 0:
        sample_term += 1.0 / out_of_transit_count

    depth_standard_error_fraction = float(residual_scatter * np.sqrt(sample_term))
    if not np.isfinite(depth_standard_error_fraction) or depth_standard_error_fraction < 0:
        return {}
    if out_of_transit_count > 0:
        baseline_standard_error_fraction = float(residual_scatter / np.sqrt(out_of_transit_count))
    else:
        baseline_standard_error_fraction = float(depth_standard_error_fraction)

    residual_fraction = residuals
    if np.isfinite(median_flux) and median_flux != 0:
        residual_fraction = residuals / median_flux
    coordinates = _fit_uncertainty_time_coordinates(fit, data.shape)
    max_beta_bin_size = min(10, max(2, in_transit_count // 4))
    beta_payload = red_noise_beta_factor(
        residual_fraction[finite_mask],
        coordinates=coordinates[finite_mask] if coordinates is not None else None,
        min_bin_size=2,
        max_bin_size=max_beta_bin_size,
    )
    red_noise_beta = finite_float(beta_payload.get('factor'), 1.0)
    if not np.isfinite(red_noise_beta) or red_noise_beta < 1.0:
        red_noise_beta = 1.0
    depth_uncertainty_fraction = float(depth_standard_error_fraction * red_noise_beta)
    baseline_uncertainty_fraction = float(baseline_standard_error_fraction * red_noise_beta)
    depth_flux_scatter_fraction = float(residual_scatter)

    parameters = getattr(fit, 'parameters', {}) or {}
    errors = getattr(fit, 'errors', {}) or {}
    rprs = finite_float(parameters.get('rprs'))
    model_rprs_uncertainty = finite_float(errors.get('rprs'))
    rprs_prior_fallback = bool(getattr(fit, 'rprs_prior_fallback_applied', False))
    if rprs_prior_fallback:
        model_rprs_uncertainty = np.nan
    data_rprs_uncertainty = np.nan
    data_rprs_standard_error = np.nan
    data_rprs_flux_scatter_uncertainty = np.nan
    combined_rprs_uncertainty = np.nan
    combined_rprs_standard_error = np.nan
    conservative_rprs_uncertainty = np.nan
    if np.isfinite(rprs) and rprs > 0:
        data_rprs_uncertainty = float(depth_uncertainty_fraction / (2.0 * rprs))
        data_rprs_standard_error = float(depth_standard_error_fraction / (2.0 * rprs))
        data_rprs_flux_scatter_uncertainty = float(depth_flux_scatter_fraction / (2.0 * rprs))
        if rprs_prior_fallback:
            combined_rprs_uncertainty = data_rprs_uncertainty
            combined_rprs_standard_error = data_rprs_standard_error
            conservative_rprs_uncertainty = data_rprs_uncertainty
        elif np.isfinite(model_rprs_uncertainty) and model_rprs_uncertainty >= 0:
            combined_rprs_uncertainty = float(
                np.sqrt(model_rprs_uncertainty ** 2 + data_rprs_uncertainty ** 2)
            )
            combined_rprs_standard_error = float(
                np.sqrt(model_rprs_uncertainty ** 2 + data_rprs_standard_error ** 2)
            )
            conservative_rprs_uncertainty = float(
                max(model_rprs_uncertainty, data_rprs_uncertainty)
            )

    return {
        'available': True,
        'residual_scatter': residual_scatter,
        'residual_scatter_percent': residual_scatter * 100.0,
        'in_transit_point_count': in_transit_count,
        'out_of_transit_point_count': out_of_transit_count,
        'transit_depth_threshold_fraction': threshold_fraction,
        'model_depth_fraction': max_depth,
        'model_depth_percent': max_depth * 100.0,
        'depth_uncertainty_fraction': depth_uncertainty_fraction,
        'depth_uncertainty_percent': depth_uncertainty_fraction * 100.0,
        'depth_red_noise_uncertainty_fraction': depth_uncertainty_fraction,
        'depth_red_noise_uncertainty_percent': depth_uncertainty_fraction * 100.0,
        'baseline_standard_error_fraction': baseline_standard_error_fraction,
        'baseline_standard_error_percent': baseline_standard_error_fraction * 100.0,
        'baseline_red_noise_uncertainty_fraction': baseline_uncertainty_fraction,
        'baseline_red_noise_uncertainty_percent': baseline_uncertainty_fraction * 100.0,
        'depth_flux_scatter_fraction': depth_flux_scatter_fraction,
        'depth_flux_scatter_percent': depth_flux_scatter_fraction * 100.0,
        'depth_standard_error_fraction': depth_standard_error_fraction,
        'depth_standard_error_percent': depth_standard_error_fraction * 100.0,
        'red_noise_beta_factor': red_noise_beta,
        'red_noise_beta_bin_sizes': beta_payload.get('bin_sizes', []),
        'red_noise_beta_by_bin': beta_payload.get('beta_by_bin', {}),
        'red_noise_beta_max_bin_size': beta_payload.get('max_bin_size', np.nan),
        'rprs': rprs,
        'model_rprs_uncertainty': model_rprs_uncertainty,
        'data_rprs_uncertainty': data_rprs_uncertainty,
        'data_rprs_red_noise_uncertainty': data_rprs_uncertainty,
        'data_rprs_standard_error': data_rprs_standard_error,
        'data_rprs_flux_scatter_uncertainty': data_rprs_flux_scatter_uncertainty,
        'combined_rprs_uncertainty': combined_rprs_uncertainty,
        'combined_rprs_red_noise_uncertainty': combined_rprs_uncertainty,
        'combined_rprs_standard_error': combined_rprs_standard_error,
        'conservative_rprs_uncertainty': conservative_rprs_uncertainty,
        'rprs_uncertainty_basis': (
            'prior_assumed_data_only'
            if rprs_prior_fallback
            else 'model_plus_red_noise'
        ),
        'rprs_prior_fallback_applied': rprs_prior_fallback,
        'rprs_prior_fallback_prior_value': finite_float(
            getattr(fit, 'rprs_prior_fallback_prior_value', np.nan)
        ),
        'rprs_prior_fallback_original_fit_value': finite_float(
            getattr(fit, 'rprs_prior_fallback_original_fit_value', np.nan)
        ),
        'rprs_prior_fallback_data_uncertainty': finite_float(
            getattr(fit, 'rprs_prior_fallback_data_uncertainty', np.nan)
        ),
        'rprs_prior_fallback_note': getattr(fit, 'rprs_prior_fallback_note', None),
    }


def photometry_method_from_info(photometry_info):
    if not isinstance(photometry_info, dict):
        return None

    min_aperture = photometry_info.get('min_aperture')
    min_aperture = finite_float(min_aperture)
    if not np.isfinite(min_aperture):
        return None
    if min_aperture == 0:
        return "PSF photometry"
    if min_aperture < 0:
        return "Aperture photometry without comparison star"
    return "Aperture photometry"


def build_aavso_qc_metadata(fit):
    transit_qc = getattr(fit, 'transit_qc', None)
    if not isinstance(transit_qc, dict):
        return {}

    fields = (
        'computed', 'status', 'summary', 'preferred_model', 'point_count',
        'transit_chi2', 'flat_chi2', 'delta_chi2', 'transit_bic', 'flat_bic',
        'delta_bic', 'transit_parameter_count', 'flat_parameter_count',
        'flat_baseline', 'flat_a2', 'flat_model_note', 'residual_scatter',
        'transit_depth_for_residual_scatter', 'residual_scatter_to_depth_ratio',
        'residual_flatness_score', 'residual_flatness_trend_strength',
        'residual_flatness_curve_strength', 'residual_flatness_scatter_ratio',
        'residual_flatness_zero_offset_strength',
        'residual_flatness_sign_imbalance', 'residual_flatness_zero_bias_score',
        'residual_flatness_trend_score', 'residual_flatness_curve_score',
        'residual_flatness_scatter_stability_score',
        'residual_flatness_dominant_metric', 'residual_flatness_detail',
        'tmid_gaussianity_score', 'tmid_gaussianity_score_uncertainty',
        'tmid_gaussianity_gaussian_distance', 'tmid_gaussianity_flat_distance',
        'tmid_gaussianity_effective_sample_count', 'tmid_gaussianity_sample_count',
        'tmid_gaussianity_robust_sigma', 'tmid_gaussianity_detail',
        'rprs_sigma', 'duration_ratio', 'eebls_depth_snr',
        'sampling_score', 'sampling_detail', 'sampling_ingress_count',
        'sampling_egress_count', 'sampling_in_transit_count',
        'sampling_pre_baseline_count', 'sampling_post_baseline_count',
        'sampling_total_duration', 'sampling_ingress_duration',
        'use_deviation_from_expected_transit_in_qc', 'deviation_sigma_threshold',
        'expected_tmid', 'expected_tmid_unc', 'fitted_tmid',
        'expected_rprs', 'expected_rprs_unc', 'fitted_rprs', 'fitted_rprs_unc',
        'rprs_deviation_fit_unc', 'rprs_deviation_model_fit_unc',
        'rprs_deviation_data_fit_unc', 'rprs_deviation_combined_fit_unc',
        'rprs_deviation_expected_unc',
        'rprs_deviation_systematic_floor', 'rprs_deviation_unc',
        'rprs_deviation_sigma', 'rprs_deviation_score',
        'rprs_prior_assumed', 'rprs_prior_assumed_note',
        'deviation_from_expected_value', 'ktmf_metric', 'ktmf_contributions',
        'notes',
    )
    return {field: transit_qc.get(field) for field in fields if field in transit_qc}


def compact_ktmf_contributions(contributions):
    compact = []
    for contribution in contributions or []:
        if not isinstance(contribution, dict):
            continue
        compact.append({
            'label': contribution.get('label'),
            'available': contribution.get('available'),
            'points': contribution.get('points'),
            'max_points': contribution.get('max_points'),
            'score': contribution.get('score'),
            'score_uncertainty': contribution.get('score_uncertainty'),
            'detail': contribution.get('detail'),
        })
    return compact


def compact_comparison_attempt_decision(attempt):
    if not isinstance(attempt, dict):
        return {}

    comp_index = attempt.get('comp_index')
    try:
        comp_number = int(comp_index) + 1
    except (TypeError, ValueError):
        comp_number = None

    return {
        'rank': attempt.get('rank'),
        'comparison_star': comp_number,
        'label': attempt.get('label'),
        'selected': attempt.get('selected'),
        'selection_reason': attempt.get('selection_reason'),
        'selection_pass_ktmf_metric': attempt.get('selection_pass_ktmf_metric'),
        'selection_pass_transit_delta_bic': attempt.get('selection_pass_transit_delta_bic'),
        'selection_pass_eebls_snr': attempt.get('selection_pass_eebls_snr'),
        'selection_pass_residual_scatter': attempt.get('selection_pass_residual_scatter'),
        'target_model_scatter_basis': attempt.get('target_model_scatter_basis'),
        'projected_full_residual_scatter': attempt.get('projected_full_residual_scatter'),
        'selection_scatter': attempt.get('selection_scatter'),
        'selection_scatter_basis': attempt.get('selection_scatter_basis'),
        'target_comp_scatter': attempt.get('target_comp_scatter'),
        'selection_pass_target_comp_scatter': attempt.get('selection_pass_target_comp_scatter'),
        'selection_pass_transit_qc_status': attempt.get('selection_pass_transit_qc_status'),
        'selection_pass_transit_qc_summary': attempt.get('selection_pass_transit_qc_summary'),
        'scatter_gate_passed': attempt.get('scatter_gate_passed'),
        'scatter_gate_lowest_residual_scatter': attempt.get('scatter_gate_lowest_residual_scatter'),
        'scatter_gate_threshold': attempt.get('scatter_gate_threshold'),
        'scatter_adjusted_ktmf_metric': attempt.get('scatter_adjusted_ktmf_metric'),
        'combined_quality_ktmf_metric': attempt.get('combined_quality_ktmf_metric'),
        'combined_quality_best_residual_scatter': attempt.get('combined_quality_best_residual_scatter'),
        'combined_quality_best_target_comp_scatter': attempt.get('combined_quality_best_target_comp_scatter'),
        'combined_quality_best_comp_stability': attempt.get('combined_quality_best_comp_stability'),
        'final_refit_metric_note': attempt.get('final_refit_metric_note'),
        'ktmf_metric': attempt.get('ktmf_metric'),
        'ktmf_contributions': compact_ktmf_contributions(attempt.get('ktmf_contributions')),
        'transit_delta_bic': attempt.get('transit_delta_bic'),
        'eebls_snr': attempt.get('eebls_snr'),
        'residual_scatter': attempt.get('residual_scatter'),
        'fit_point_count': attempt.get('fit_point_count'),
        'transit_qc_status': attempt.get('transit_qc_status'),
        'transit_qc_summary': attempt.get('transit_qc_summary'),
        'rejected_by_transit_qc': attempt.get('rejected_by_transit_qc'),
        'failure_reason': attempt.get('failure_reason'),
    }


def compact_comparison_attempt_decisions(attempts, limit=10):
    attempts = list(attempts or [])
    return {
        'candidate_count': len(attempts),
        'candidates': [
            compact_comparison_attempt_decision(attempt)
            for attempt in attempts[:limit]
        ],
        'omitted_candidate_count': max(0, len(attempts) - limit),
    }


def build_ktmf_decision_metadata(fit, photometry_info=None):
    transit_qc = getattr(fit, 'transit_qc', None)
    payload = {}
    if isinstance(transit_qc, dict):
        payload['target_fit'] = {
            'status': transit_qc.get('status'),
            'summary': transit_qc.get('summary'),
            'ktmf_metric': transit_qc.get('ktmf_metric'),
            'ktmf_contributions': compact_ktmf_contributions(transit_qc.get('ktmf_contributions')),
            'delta_bic': transit_qc.get('delta_bic'),
            'delta_chi2': transit_qc.get('delta_chi2'),
            'eebls_depth_snr': transit_qc.get('eebls_depth_snr'),
            'residual_scatter': transit_qc.get('residual_scatter'),
            'transit_depth_for_residual_scatter': transit_qc.get('transit_depth_for_residual_scatter'),
            'residual_scatter_to_depth_ratio': transit_qc.get('residual_scatter_to_depth_ratio'),
            'residual_flatness_score': transit_qc.get('residual_flatness_score'),
            'residual_flatness_detail': transit_qc.get('residual_flatness_detail'),
            'residual_flatness_trend_strength': transit_qc.get('residual_flatness_trend_strength'),
            'residual_flatness_curve_strength': transit_qc.get('residual_flatness_curve_strength'),
            'residual_flatness_scatter_ratio': transit_qc.get('residual_flatness_scatter_ratio'),
            'residual_flatness_zero_offset_strength': transit_qc.get('residual_flatness_zero_offset_strength'),
            'residual_flatness_sign_imbalance': transit_qc.get('residual_flatness_sign_imbalance'),
            'residual_flatness_zero_bias_score': transit_qc.get('residual_flatness_zero_bias_score'),
            'sampling_score': transit_qc.get('sampling_score'),
            'sampling_detail': transit_qc.get('sampling_detail'),
            'sampling_ingress_count': transit_qc.get('sampling_ingress_count'),
            'sampling_egress_count': transit_qc.get('sampling_egress_count'),
            'sampling_in_transit_count': transit_qc.get('sampling_in_transit_count'),
            'sampling_pre_baseline_count': transit_qc.get('sampling_pre_baseline_count'),
            'sampling_post_baseline_count': transit_qc.get('sampling_post_baseline_count'),
            'deviation_from_expected_value': transit_qc.get('deviation_from_expected_value'),
            'rprs_deviation_fit_unc': transit_qc.get('rprs_deviation_fit_unc'),
            'rprs_deviation_model_fit_unc': transit_qc.get('rprs_deviation_model_fit_unc'),
            'rprs_deviation_data_fit_unc': transit_qc.get('rprs_deviation_data_fit_unc'),
            'rprs_deviation_expected_unc': transit_qc.get('rprs_deviation_expected_unc'),
            'rprs_deviation_unc': transit_qc.get('rprs_deviation_unc'),
        }

    if isinstance(photometry_info, dict):
        selected_attempt = photometry_info.get('selected_comparison_attempt')
        selected_payload = compact_comparison_attempt_decision(selected_attempt)
        if not selected_payload:
            selected_payload = {
                'comparison_star': photometry_info.get('comp_star_num'),
                'selected': photometry_info.get('comp_star_num') is not None,
                'selection_reason': photometry_info.get('selected_comparison_selection_reason'),
                'ktmf_metric': photometry_info.get('comparison_ktmf_metric'),
                'ktmf_contributions': compact_ktmf_contributions(
                    photometry_info.get('selected_comparison_ktmf_contributions')
                ),
                'transit_delta_bic': photometry_info.get('comparison_transit_delta_bic'),
                'eebls_snr': photometry_info.get('comparison_eebls_snr'),
                'fit_point_count': photometry_info.get('selected_comparison_fit_point_count'),
                'transit_qc_status': photometry_info.get('selected_comparison_transit_qc_status'),
                'transit_qc_summary': photometry_info.get('selected_comparison_transit_qc_summary'),
            }

        payload['comparison_selection'] = {
            'basis': photometry_info.get('selection_basis'),
            'metric': photometry_info.get('selection_metric'),
            'field_score': photometry_info.get('calibration_field_score'),
            'selected': selected_payload,
        }

        attempt_summary = compact_comparison_attempt_decisions(
            photometry_info.get('comparison_fit_attempt_summaries')
        )
        if attempt_summary['candidate_count']:
            payload['comparison_selection'].update(attempt_summary)

    return payload


def format_ktmf_metric(value):
    value = finite_float(value)
    return f"{value:.2f} / 5.00" if np.isfinite(value) else "n/a"


def format_ktmf_status(value):
    value = finite_float(value)
    if not np.isfinite(value):
        return "UNKNOWN"
    if value >= 4.0:
        return "PASS"
    if value >= 3.0:
        return "MARGINAL"
    return "FAIL"


def format_transit_qc_headline_final_params(transit_qc):
    params = {}
    if not isinstance(transit_qc, dict) or not transit_qc:
        return params

    qc_status = transit_qc.get('status')
    if qc_status:
        params["Transit detection QC"] = str(qc_status).upper()

    qc_ktmf = finite_float(transit_qc.get('ktmf_metric'))
    if np.isfinite(qc_ktmf):
        params["KTMF"] = f"{qc_ktmf:.2f} / 5.00"

    return params


def format_optional_metric(label, value, precision=2):
    value = finite_float(value)
    if not np.isfinite(value):
        return None
    return f"{label}={value:.{precision}f}"


def format_transit_delta_bic(value):
    value = finite_float(value)
    return f"{value:.2f}" if np.isfinite(value) else "n/a"


def format_percent_metric(label, value, precision=4):
    value = finite_float(value)
    if not np.isfinite(value):
        return None
    return f"{label}={value * 100.0:.{precision}f}%"


def metric_values_differ(first_value, second_value, tolerance=5.0e-3):
    first_value = finite_float(first_value)
    second_value = finite_float(second_value)
    return (
        np.isfinite(first_value)
        and np.isfinite(second_value)
        and abs(first_value - second_value) > tolerance
    )


def format_ktmf_candidate_decision(attempt):
    attempt = compact_comparison_attempt_decision(attempt)
    label = attempt.get('label') or (
        f"Comp {attempt['comparison_star']}" if attempt.get('comparison_star') is not None else "Comparison candidate"
    )
    selected_text = " [selected]" if attempt.get('selected') else ""
    ktmf_text = f"KTMF={format_ktmf_metric(attempt.get('ktmf_metric'))}"
    if metric_values_differ(
        attempt.get('selection_pass_ktmf_metric'),
        attempt.get('ktmf_metric'),
    ):
        ktmf_text += (
            f" (selection-pass {format_ktmf_metric(attempt.get('selection_pass_ktmf_metric'))})"
        )
    parts = [f"{label}{selected_text}: {ktmf_text}"]
    for metric_text in (
        format_optional_metric("Delta BIC", attempt.get('transit_delta_bic')),
        format_optional_metric("EEBLS SNR", attempt.get('eebls_snr')),
        format_percent_metric("Target/comp scatter", attempt.get('target_comp_scatter')),
        format_percent_metric("Target model scatter", attempt.get('residual_scatter')),
        format_percent_metric("Selection scatter", attempt.get('selection_scatter')),
        format_optional_metric("KTMF/projected-scatter score", attempt.get('combined_quality_ktmf_metric')),
    ):
        if metric_text:
            parts.append(metric_text)
    target_model_basis = attempt.get('target_model_scatter_basis')
    selection_scatter_basis = attempt.get('selection_scatter_basis')
    if target_model_basis:
        parts.append(f"target model scatter basis={target_model_basis}")
    if selection_scatter_basis:
        parts.append(f"selection scatter basis={selection_scatter_basis}")
    selection_pass_target_comp_scatter = finite_float(attempt.get('selection_pass_target_comp_scatter'))
    target_comp_scatter = finite_float(attempt.get('target_comp_scatter'))
    if (
        np.isfinite(selection_pass_target_comp_scatter)
        and (
            not np.isfinite(target_comp_scatter)
            or abs(selection_pass_target_comp_scatter - target_comp_scatter) > 1.0e-5
        )
    ):
        parts.append(
            "selection-pass target/comp scatter="
            f"{selection_pass_target_comp_scatter * 100.0:.4f}%"
        )
    if metric_values_differ(
        attempt.get('selection_pass_residual_scatter'),
        attempt.get('residual_scatter'),
        tolerance=1.0e-5,
    ):
        parts.append(
            "selection-pass target model scatter="
            f"{finite_float(attempt.get('selection_pass_residual_scatter')) * 100.0:.4f}%"
        )
    if metric_values_differ(
        attempt.get('selection_pass_transit_delta_bic'),
        attempt.get('transit_delta_bic'),
        tolerance=1.0e-2,
    ):
        parts.append(
            "selection-pass Delta BIC="
            f"{format_transit_delta_bic(attempt.get('selection_pass_transit_delta_bic'))}"
        )
    if metric_values_differ(
        attempt.get('selection_pass_eebls_snr'),
        attempt.get('eebls_snr'),
        tolerance=1.0e-2,
    ):
        parts.append(
            "selection-pass EEBLS SNR="
            f"{finite_float(attempt.get('selection_pass_eebls_snr')):.2f}"
        )
    qc_status = attempt.get('transit_qc_status')
    if qc_status:
        parts.append(f"QC={str(qc_status).upper()}")
    reason = attempt.get('selection_reason') or attempt.get('failure_reason')
    if reason:
        parts.append(f"reason={reason}")
    return ", ".join(parts)


def format_ktmf_decision_final_params(fit, photometry_info=None):
    params = {}

    transit_qc = getattr(fit, 'transit_qc', None)
    if isinstance(transit_qc, dict):
        ktmf_metric = finite_float(transit_qc.get('ktmf_metric'))
        if np.isfinite(ktmf_metric):
            target_status = format_ktmf_status(ktmf_metric)
            params["KTMF target-fit decision"] = (
                f"{target_status}: KTMF={format_ktmf_metric(ktmf_metric)}"
            )
        for contribution_index, contribution in enumerate(
            compact_ktmf_contributions(transit_qc.get('ktmf_contributions')),
            start=1,
        ):
            label = contribution.get('label', f'Component {contribution_index}')
            available = bool(contribution.get('available'))
            points = finite_float(contribution.get('points'), 0.0)
            max_points = finite_float(contribution.get('max_points'), 0.0)
            score = finite_float(contribution.get('score'))
            detail = contribution.get('detail') or 'n/a'
            if available and np.isfinite(score):
                params[f"KTMF target contribution {contribution_index}"] = (
                    f"{label}: +{points:.2f}/{max_points:.2f} (score={score:.2f}; {detail})"
                )
            else:
                params[f"KTMF target contribution {contribution_index}"] = (
                    f"{label}: +0.00/0.00 (unavailable; {detail})"
                )

    if not isinstance(photometry_info, dict):
        return params

    basis = photometry_info.get('selection_basis')
    metric = photometry_info.get('selection_metric')
    if basis or metric:
        params["KTMF comparison selection mode"] = (
            f"basis={basis or 'n/a'}, metric={metric or 'n/a'}"
        )

    selected_attempt = photometry_info.get('selected_comparison_attempt')
    if selected_attempt:
        params["KTMF selected comparison decision"] = format_ktmf_candidate_decision(selected_attempt)
    elif photometry_info.get('comp_star_num') is not None:
        selected_payload = {
            'label': f"Comp {photometry_info.get('comp_star_num')}",
            'selected': True,
            'selection_reason': photometry_info.get('selected_comparison_selection_reason'),
            'ktmf_metric': photometry_info.get('comparison_ktmf_metric'),
            'ktmf_contributions': photometry_info.get('selected_comparison_ktmf_contributions'),
            'transit_delta_bic': photometry_info.get('comparison_transit_delta_bic'),
            'eebls_snr': photometry_info.get('comparison_eebls_snr'),
            'transit_qc_status': photometry_info.get('selected_comparison_transit_qc_status'),
        }
        params["KTMF selected comparison decision"] = format_ktmf_candidate_decision(selected_payload)

    for contribution_index, contribution in enumerate(
        compact_ktmf_contributions(photometry_info.get('selected_comparison_ktmf_contributions')),
        start=1,
    ):
        label = contribution.get('label', f'Component {contribution_index}')
        available = bool(contribution.get('available'))
        points = finite_float(contribution.get('points'), 0.0)
        max_points = finite_float(contribution.get('max_points'), 0.0)
        score = finite_float(contribution.get('score'))
        detail = contribution.get('detail') or 'n/a'
        if available and np.isfinite(score):
            params[f"KTMF selected comparison contribution {contribution_index}"] = (
                f"{label}: +{points:.2f}/{max_points:.2f} (score={score:.2f}; {detail})"
            )
        else:
            params[f"KTMF selected comparison contribution {contribution_index}"] = (
                f"{label}: +0.00/0.00 (unavailable; {detail})"
            )

    for attempt_index, attempt in enumerate(
        (photometry_info.get('comparison_fit_attempt_summaries') or [])[:10],
        start=1,
    ):
        params[f"KTMF comparison candidate {attempt_index}"] = format_ktmf_candidate_decision(attempt)

    return params


def format_fit_quality_final_params(fit_quality):
    fit_quality = fit_quality or {}
    params = {}

    reduced_chi_square = finite_float(fit_quality.get('reduced_chi_square'))
    if np.isfinite(reduced_chi_square):
        params["Fit quality reduced chi-square"] = f"{reduced_chi_square:.3f}"

    chi_square = finite_float(fit_quality.get('chi_square'))
    if np.isfinite(chi_square):
        params["Fit quality chi-square"] = f"{chi_square:.2f}"

    degrees_of_freedom = fit_quality.get('degrees_of_freedom')
    try:
        degrees_of_freedom = int(degrees_of_freedom)
    except (TypeError, ValueError):
        degrees_of_freedom = None
    if degrees_of_freedom is not None:
        params["Fit quality degrees of freedom"] = str(degrees_of_freedom)

    rms_residual_percent = finite_float(fit_quality.get('rms_residual_percent'))
    if np.isfinite(rms_residual_percent):
        params["Fit quality RMS residual"] = f"{rms_residual_percent:.4f} %"

    median_abs_normalized_residual = finite_float(
        fit_quality.get('median_absolute_normalized_residual')
    )
    if np.isfinite(median_abs_normalized_residual):
        params["Fit quality median absolute normalized residual"] = (
            f"{median_abs_normalized_residual:.2f} sigma"
        )

    rms_uncertainty_ratio = finite_float(fit_quality.get('rms_residual_to_median_uncertainty'))
    if np.isfinite(rms_uncertainty_ratio):
        params["Fit quality RMS residual / median uncertainty"] = f"{rms_uncertainty_ratio:.2f}"

    point_count = fit_quality.get('weighted_point_count', fit_quality.get('point_count'))
    try:
        point_count = int(point_count)
    except (TypeError, ValueError):
        point_count = None
    if point_count is not None:
        params["Fit quality point count"] = str(point_count)

    if fit_quality and not fit_quality.get('uses_uncertainties'):
        params["Fit quality note"] = "Per-point uncertainties unavailable; chi-square metrics not reported."

    return params


def format_empirical_transit_uncertainty_final_params(empirical_uncertainty):
    empirical_uncertainty = empirical_uncertainty or {}
    if not empirical_uncertainty.get('available'):
        return {}

    params = {}
    rprs = finite_float(empirical_uncertainty.get('rprs'))
    model_rprs_uncertainty = finite_float(empirical_uncertainty.get('model_rprs_uncertainty'))
    data_rprs_uncertainty = finite_float(empirical_uncertainty.get('data_rprs_uncertainty'))
    combined_rprs_uncertainty = finite_float(empirical_uncertainty.get('combined_rprs_uncertainty'))
    conservative_rprs_uncertainty = finite_float(
        empirical_uncertainty.get('conservative_rprs_uncertainty')
    )
    data_rprs_standard_error = finite_float(empirical_uncertainty.get('data_rprs_standard_error'))
    data_rprs_flux_scatter_uncertainty = finite_float(
        empirical_uncertainty.get('data_rprs_flux_scatter_uncertainty')
    )
    combined_rprs_standard_error = finite_float(
        empirical_uncertainty.get('combined_rprs_standard_error')
    )
    depth_uncertainty_percent = finite_float(empirical_uncertainty.get('depth_uncertainty_percent'))
    depth_flux_scatter_percent = finite_float(
        empirical_uncertainty.get('depth_flux_scatter_percent')
    )
    depth_standard_error_percent = finite_float(
        empirical_uncertainty.get('depth_standard_error_percent')
    )
    baseline_red_noise_percent = finite_float(
        empirical_uncertainty.get('baseline_red_noise_uncertainty_percent')
    )
    baseline_standard_error_percent = finite_float(
        empirical_uncertainty.get('baseline_standard_error_percent')
    )
    residual_scatter_percent = finite_float(empirical_uncertainty.get('residual_scatter_percent'))
    red_noise_beta = finite_float(empirical_uncertainty.get('red_noise_beta_factor'))
    rprs_prior_fallback = bool(empirical_uncertainty.get('rprs_prior_fallback_applied'))
    rprs_uncertainty_basis = empirical_uncertainty.get('rprs_uncertainty_basis')

    if (
        not rprs_prior_fallback
        and np.isfinite(rprs)
        and np.isfinite(model_rprs_uncertainty)
        and model_rprs_uncertainty >= 0
    ):
        params["Ratio of Planet to Stellar Radius (Rp/R*) model-fit uncertainty"] = (
            f"{round_to_2(rprs, model_rprs_uncertainty)} +/- {round_to_2(model_rprs_uncertainty)}"
        )
    if np.isfinite(rprs) and np.isfinite(data_rprs_uncertainty) and data_rprs_uncertainty >= 0:
        if rprs_prior_fallback:
            params["Ratio of Planet to Stellar Radius (Rp/R*) prior-assumed data-only uncertainty"] = (
                f"{round_to_2(rprs, data_rprs_uncertainty)} +/- {round_to_2(data_rprs_uncertainty)}"
            )
        else:
            params["Ratio of Planet to Stellar Radius (Rp/R*) data-fit red-noise uncertainty"] = (
                f"{round_to_2(rprs, data_rprs_uncertainty)} +/- {round_to_2(data_rprs_uncertainty)}"
            )
    if np.isfinite(rprs) and np.isfinite(combined_rprs_uncertainty) and combined_rprs_uncertainty >= 0:
        if rprs_prior_fallback:
            params["Ratio of Planet to Stellar Radius (Rp/R*) data-only uncertainty used for primary value"] = (
                f"{round_to_2(rprs, combined_rprs_uncertainty)} +/- "
                f"{round_to_2(combined_rprs_uncertainty)}"
            )
        else:
            params["Ratio of Planet to Stellar Radius (Rp/R*) model+red-noise uncertainty"] = (
                f"{round_to_2(rprs, combined_rprs_uncertainty)} +/- "
                f"{round_to_2(combined_rprs_uncertainty)}"
            )
    if np.isfinite(conservative_rprs_uncertainty):
        params["Conservative Rp/R* uncertainty to quote"] = (
            f"+/- {round_to_2(conservative_rprs_uncertainty)}"
        )
    if np.isfinite(rprs) and np.isfinite(data_rprs_standard_error) and data_rprs_standard_error >= 0:
        params["Ratio of Planet to Stellar Radius (Rp/R*) data-fit standard-error estimate"] = (
            f"{round_to_2(rprs, data_rprs_standard_error)} +/- {round_to_2(data_rprs_standard_error)}"
        )
    if (
        not rprs_prior_fallback
        and np.isfinite(rprs)
        and np.isfinite(combined_rprs_standard_error)
        and combined_rprs_standard_error >= 0
    ):
        params["Ratio of Planet to Stellar Radius (Rp/R*) model+standard-error estimate"] = (
            f"{round_to_2(rprs, combined_rprs_standard_error)} +/- "
            f"{round_to_2(combined_rprs_standard_error)}"
        )
    if (
        np.isfinite(rprs)
        and np.isfinite(data_rprs_flux_scatter_uncertainty)
        and data_rprs_flux_scatter_uncertainty >= 0
    ):
        params["Ratio of Planet to Stellar Radius (Rp/R*) flux-scatter equivalent"] = (
            f"{round_to_2(rprs, data_rprs_flux_scatter_uncertainty)} +/- "
            f"{round_to_2(data_rprs_flux_scatter_uncertainty)}"
        )
    if np.isfinite(depth_uncertainty_percent):
        params["Transit depth red-noise uncertainty"] = (
            f"+/- {depth_uncertainty_percent:.4f} %"
        )
    if np.isfinite(depth_flux_scatter_percent):
        params["Transit depth flux-scatter equivalent"] = (
            f"+/- {depth_flux_scatter_percent:.4f} %"
        )
    if np.isfinite(depth_standard_error_percent):
        params["Transit depth data-fit standard-error estimate"] = (
            f"+/- {depth_standard_error_percent:.4f} %"
        )
    if np.isfinite(baseline_red_noise_percent):
        params["Flux baseline red-noise uncertainty"] = (
            f"+/- {baseline_red_noise_percent:.4f} %"
        )
    if np.isfinite(baseline_standard_error_percent):
        params["Flux baseline standard-error estimate"] = (
            f"+/- {baseline_standard_error_percent:.4f} %"
        )
    if np.isfinite(red_noise_beta):
        params["Red-noise beta factor"] = f"{red_noise_beta:.3f}"
    if rprs_uncertainty_basis:
        params["Rp/R* uncertainty basis"] = str(rprs_uncertainty_basis)
    fallback_note = empirical_uncertainty.get('rprs_prior_fallback_note')
    if fallback_note:
        params["Rp/R* prior fallback note"] = str(fallback_note)

    beta_bins = empirical_uncertainty.get('red_noise_beta_bin_sizes')
    if beta_bins:
        params["Red-noise beta bin sizes"] = ", ".join(str(int(item)) for item in beta_bins)

    in_count = empirical_uncertainty.get('in_transit_point_count')
    out_count = empirical_uncertainty.get('out_of_transit_point_count')
    try:
        in_count = int(in_count)
        out_count = int(out_count)
    except (TypeError, ValueError):
        in_count = None
        out_count = None
    if in_count is not None and out_count is not None:
        params["Data-fit uncertainty point counts"] = (
            f"{in_count} in transit, {out_count} out of transit"
        )

    if np.isfinite(residual_scatter_percent):
        params["Flux residual scatter around model"] = (
            f"{residual_scatter_percent:.4f} %"
        )
    if rprs_prior_fallback:
        params["Uncertainty interpretation note"] = (
            "The primary Rp/R* and radius-ratio area-depth uncertainties use the input "
            "prior Rp/R* value with a data-only red-noise uncertainty because the sampled "
            "Rp/R* posterior was pinned against a search bound and automatic Rp/R* "
            "posterior expansion was disabled. Tmid, a/Rs, inclination, and impact "
            "parameter remain fitted parameters; their primary uncertainties apply the "
            "same residual time-binning beta factor to the posterior uncertainties. "
            "The flux baseline red-noise uncertainty is the out-of-transit baseline "
            "component used for the final-plot baseline band."
        )
    else:
        params["Uncertainty interpretation note"] = (
            "The primary Rp/R* and radius-ratio area-depth uncertainties use model+red-noise "
            "when available. The red-noise uncertainty inflates the data-fit standard error "
            "by a residual time-binning beta factor before combining it with the model "
            "posterior. The flux baseline red-noise uncertainty is the out-of-transit "
            "baseline component used for the final-plot baseline band; the transit-depth "
            "red-noise uncertainty already includes both the in-transit and baseline terms. "
            "The same residual time-binning beta factor is applied to the posterior "
            "uncertainties for Tmid, a/Rs, inclination, and impact parameter when reporting "
            "their primary model+red-noise uncertainties. "
            "The flux-scatter equivalent is also shown as a diagnostic of the full residual "
            "scatter around the model."
        )
    return params


def build_aavso_photometry_metadata(photometry_info):
    if not isinstance(photometry_info, dict):
        return {}

    selected_source_indices = photometry_info.get('selected_source_indices')
    selected_source_count = None
    if selected_source_indices is not None:
        try:
            selected_source_count = int(np.asarray(selected_source_indices).size)
        except (TypeError, ValueError):
            selected_source_count = None

    selected_times = numeric_series_summary(photometry_info.get('selected_fit_good_times'))
    return {
        'method': photometry_method_from_info(photometry_info),
        'selected_comparison_star': photometry_info.get('comp_star_num'),
        'selected_comparison_coordinates': photometry_info.get('comp_star_coords'),
        'noise_budget_summary': photometry_info.get('noise_budget_summary'),
        'noise_budget_terms': photometry_info.get('noise_budget_terms'),
        'comparison_selection_basis': photometry_info.get('selection_basis'),
        'comparison_selection_metric': photometry_info.get('selection_metric'),
        'comparison_field_score': photometry_info.get('calibration_field_score'),
        'comparison_field_score_percent': (
            100.0 * finite_float(photometry_info.get('calibration_field_score'))
            if np.isfinite(finite_float(photometry_info.get('calibration_field_score')))
            else np.nan
        ),
        'selected_comparison_ktmf': photometry_info.get('comparison_ktmf_metric'),
        'selected_comparison_eebls_snr': photometry_info.get('comparison_eebls_snr'),
        'selected_comparison_transit_delta_bic': photometry_info.get('comparison_transit_delta_bic'),
        'reused_selected_full_reduction_fit': photometry_info.get('reuse_selected_full_reduction_fit'),
        'selected_source_point_count': selected_source_count,
        'selected_fit_time_range': selected_times,
    }


def build_aavso_aperture_metadata(photometry_info):
    if not isinstance(photometry_info, dict):
        return {}

    adaptive_summary = photometry_info.get('adaptive_summary')
    payload = {
        'method': photometry_method_from_info(photometry_info),
        'aperture_index': photometry_info.get('aperture_index'),
        'annulus_index': photometry_info.get('annulus_index'),
        'configured_aperture_px': photometry_info.get('min_aperture'),
        'configured_annulus_px': photometry_info.get('min_annulus'),
        'adaptive': adaptive_summary is not None,
    }
    if not isinstance(adaptive_summary, dict):
        return payload

    payload.update({
        'aperture_sigma': adaptive_summary.get('aperture_sigma'),
        'annulus_sigma': adaptive_summary.get('annulus_sigma'),
        'aperture_px': {
            'median': adaptive_summary.get('aperture_median'),
            'std': adaptive_summary.get('aperture_std'),
            'min': adaptive_summary.get('aperture_min'),
            'max': adaptive_summary.get('aperture_max'),
        },
        'annulus_px': {
            'median': adaptive_summary.get('annulus_median'),
            'std': adaptive_summary.get('annulus_std'),
            'min': adaptive_summary.get('annulus_min'),
            'max': adaptive_summary.get('annulus_max'),
        },
        'fwhm_px': numeric_series_summary(adaptive_summary.get('fwhm_series')),
        'frame_sigma_px': numeric_series_summary(adaptive_summary.get('frame_sigma')),
        'sky_inner_px': numeric_series_summary(adaptive_summary.get('sky_inner_series')),
        'sky_outer_px': numeric_series_summary(adaptive_summary.get('sky_outer_series')),
        'sky_pixels': numeric_series_summary(adaptive_summary.get('sky_pixel_series')),
    })
    return payload


def build_aavso_frame_filtering_metadata(fit, frame_filtering_info):
    payload = dict(frame_filtering_info or {})

    for source_key, target_key in (
        ('dropped_missing_wcs_files', 'missing_wcs_rejections'),
        ('dropped_target_wcs_files', 'target_wcs_rejections'),
        ('dropped_pointing_files', 'pointing_rejections'),
    ):
        if source_key in payload:
            payload[target_key] = file_list_summary(payload.pop(source_key))

    diagnostics = getattr(fit, 'frame_filter_diagnostics', None)
    if diagnostics:
        payload['lightcurve_filter_diagnostics'] = diagnostics
        payload['lightcurve_dropped_point_count'] = sum(
            int((diagnostic or {}).get('dropped_point_count', 0))
            for diagnostic in diagnostics
        )
    final_residual_rejection = getattr(fit, 'final_residual_rejection', None)
    if isinstance(final_residual_rejection, dict) and final_residual_rejection:
        payload['final_residual_rejection'] = aavso_json_safe(final_residual_rejection)
    return payload


def build_aavso_astrometry_metadata(astrometry_info, comp_star):
    payload = dict(astrometry_info or {})
    if payload.get('wcs_file'):
        payload['wcs_file'] = path_name(payload['wcs_file'])
    if comp_star:
        payload['comparison_star_aavso_header'] = comp_star
    return payload


def build_aavso_bad_pixel_metadata(bad_pixel_info):
    if not isinstance(bad_pixel_info, dict):
        return {}

    payload = dict(bad_pixel_info)
    for key in ('counts_path', 'mask_path'):
        if payload.get(key):
            payload[key] = path_name(payload[key])
    return payload


def format_parameter_with_error(value, error):
    value = finite_float(value)
    error = finite_float(error)
    if not np.isfinite(value):
        return None
    if np.isfinite(error) and error >= 0:
        return f"{round_to_2(value, error)} +/- {round_to_2(error)}"
    return f"{round_to_2(value)} +/- n/a"


def format_percent_parameter_with_error(value, error):
    text = format_parameter_with_error(value, error)
    return f"{text} [%]" if text is not None else None


def formatted_transit_depth_parameters(fit, planet_dict=None, limb_darkening=None, empirical_uncertainty=None):
    prior_parameters = planet_dict_transit_parameters(
        planet_dict,
        limb_darkening=limb_darkening,
        fallback=getattr(fit, 'prior', None),
    )
    prior_errors = planet_dict_transit_errors(planet_dict, limb_darkening=limb_darkening)
    summary = fit_transit_depth_summary(
        fit,
        prior_parameters=prior_parameters,
        prior_errors=prior_errors,
    )

    entries = {}
    for label, value_key, error_key in (
        (AREA_DEPTH_LABEL, 'area_depth', 'area_depth_error'),
        (OBSERVABLE_DEPTH_LABEL, 'observable_depth', 'observable_depth_error'),
        (PRIOR_OBSERVABLE_DEPTH_LABEL, 'prior_observable_depth', 'prior_observable_depth_error'),
        (OBSERVABLE_DEPTH_DELTA_LABEL, 'observable_depth_prior_delta', 'observable_depth_prior_delta_error'),
    ):
        text = format_percent_parameter_with_error(summary.get(value_key), summary.get(error_key))
        if text is not None:
            entries[label] = text

    empirical_uncertainty = empirical_uncertainty or {}
    if empirical_uncertainty.get('available'):
        area_depth = finite_float(summary.get('area_depth'))
        rprs = finite_float(empirical_uncertainty.get('rprs'))
        rprs_prior_fallback = bool(empirical_uncertainty.get('rprs_prior_fallback_applied'))

        def area_error_percent_from_rprs_error(error):
            error = finite_float(error)
            if np.isfinite(area_depth) and np.isfinite(rprs) and np.isfinite(error) and error >= 0:
                return float(200.0 * abs(rprs) * error)
            return np.nan

        model_area_error = area_error_percent_from_rprs_error(
            empirical_uncertainty.get('model_rprs_uncertainty')
        )
        data_area_error = area_error_percent_from_rprs_error(
            empirical_uncertainty.get('data_rprs_uncertainty')
        )
        combined_area_error = area_error_percent_from_rprs_error(
            empirical_uncertainty.get('combined_rprs_uncertainty')
        )
        flux_scatter_area_error = area_error_percent_from_rprs_error(
            empirical_uncertainty.get('data_rprs_flux_scatter_uncertainty')
        )
        if np.isfinite(area_depth) and np.isfinite(combined_area_error):
            combined_text = format_percent_parameter_with_error(area_depth, combined_area_error)
            if combined_text is not None:
                entries[AREA_DEPTH_LABEL] = combined_text
                if rprs_prior_fallback:
                    entries[f"{AREA_DEPTH_LABEL} prior-assumed data-only uncertainty"] = combined_text
                else:
                    entries[f"{AREA_DEPTH_LABEL} model+red-noise uncertainty"] = combined_text
        if not rprs_prior_fallback and np.isfinite(area_depth) and np.isfinite(model_area_error):
            text = format_percent_parameter_with_error(area_depth, model_area_error)
            if text is not None:
                entries[f"{AREA_DEPTH_LABEL} model-fit uncertainty"] = text
        if np.isfinite(area_depth) and np.isfinite(data_area_error):
            text = format_percent_parameter_with_error(area_depth, data_area_error)
            if text is not None:
                if rprs_prior_fallback:
                    entries[f"{AREA_DEPTH_LABEL} prior-assumed data-fit uncertainty"] = text
                else:
                    entries[f"{AREA_DEPTH_LABEL} data-fit red-noise uncertainty"] = text
        if np.isfinite(area_depth) and np.isfinite(flux_scatter_area_error):
            text = format_percent_parameter_with_error(area_depth, flux_scatter_area_error)
            if text is not None:
                entries[f"{AREA_DEPTH_LABEL} flux-scatter equivalent"] = text
    return entries


def empirical_red_noise_error_scale(empirical_uncertainty):
    empirical_uncertainty = empirical_uncertainty or {}
    if not empirical_uncertainty.get('available'):
        return 1.0

    beta = finite_float(empirical_uncertainty.get('red_noise_beta_factor'))
    if np.isfinite(beta) and beta > 1.0:
        return float(beta)
    return 1.0


def fit_parameter_model_data_uncertainty(fit, parameter_name, empirical_uncertainty=None):
    errors = getattr(fit, 'errors', {}) or {}
    model_error = finite_float(errors.get(parameter_name))
    if not np.isfinite(model_error) or model_error < 0:
        return np.nan

    return float(model_error * empirical_red_noise_error_scale(empirical_uncertainty))


def fit_rprs_report_error(fit, empirical_uncertainty=None):
    if empirical_uncertainty is None:
        empirical_uncertainty = fit_empirical_transit_uncertainty(fit)
    if not isinstance(empirical_uncertainty, dict):
        empirical_uncertainty = {}

    report_error = finite_float(empirical_uncertainty.get('combined_rprs_uncertainty'))
    if np.isfinite(report_error) and report_error >= 0:
        return report_error

    errors = getattr(fit, 'errors', {}) or {}
    report_error = finite_float(errors.get('rprs'))
    if np.isfinite(report_error) and report_error >= 0:
        return report_error

    report_error = finite_float(getattr(fit, 'rprs_prior_fallback_data_uncertainty', np.nan))
    if np.isfinite(report_error) and report_error >= 0:
        return report_error

    return np.nan


def fit_impact_parameter_value_error(fit, errors_override=None):
    parameters = getattr(fit, 'parameters', {}) or {}
    errors = getattr(fit, 'errors', {}) or {}
    errors_override = errors_override or {}
    sample_parameters = getattr(fit, 'sample_parameters', {}) or {}
    sample_errors = getattr(fit, 'sample_errors', {}) or {}

    if 'b' in sample_parameters:
        impact_parameter = finite_float(sample_parameters.get('b'))
        impact_error = finite_float(errors_override.get('b'), finite_float(sample_errors.get('b')))
        if np.isfinite(impact_parameter):
            return impact_parameter, impact_error

    if 'b' in parameters:
        impact_parameter = finite_float(parameters.get('b'))
        impact_error = finite_float(errors_override.get('b'), finite_float(errors.get('b')))
        if np.isfinite(impact_parameter):
            return impact_parameter, impact_error

    ars = finite_float(parameters.get('ars'))
    inc = finite_float(parameters.get('inc'))
    if not np.isfinite(ars) or not np.isfinite(inc):
        return np.nan, np.nan

    ecc = finite_float(parameters.get('ecc'), 0.0)
    omega = np.deg2rad(finite_float(parameters.get('omega'), 0.0))
    denominator = 1.0 + ecc * np.sin(omega)
    if not np.isfinite(denominator) or np.isclose(denominator, 0.0):
        return np.nan, np.nan

    scale_factor = (1.0 - ecc ** 2) / denominator
    inc_rad = np.deg2rad(inc)
    impact_parameter = scale_factor * ars * np.cos(inc_rad)

    ars_error = finite_float(errors_override.get('ars'), finite_float(errors.get('ars')))
    inc_error = finite_float(errors_override.get('inc'), finite_float(errors.get('inc')))
    if np.isfinite(ars_error) and np.isfinite(inc_error):
        impact_error = np.hypot(
            scale_factor * np.cos(inc_rad) * ars_error,
            scale_factor * ars * np.sin(inc_rad) * np.deg2rad(inc_error),
        )
    else:
        impact_error = np.nan

    return float(impact_parameter), float(impact_error) if np.isfinite(impact_error) else np.nan


class OutputFiles:
    def __init__(self, fit, p_dict, i_dict, durs):
        self.fit = fit
        self.p_dict = p_dict
        self.i_dict = i_dict
        self.durs = durs
        self.dir = Path(self.i_dict['save'])

    def final_lightcurve(self, phase):
        params_file = self.dir / "working_artifacts" / safe_output_filename(
            "FinalLightCurve",
            self.p_dict['pName'],
            filename_date_token(self.i_dict['date']),
            extension="csv",
        )

        if getattr(self.fit, 'stellar_variability_only', False):
            vsp_params = getattr(self.fit, 'stellar_variability_params', None) or []
            with params_file.open('w') as f:
                target_name = self.p_dict.get('sName', self.p_dict['pName'])
                f.write(f"# FINAL STELLAR VARIABILITY TIMESERIES OF {target_name}\n")
                reference_label = (
                    getattr(self.fit, 'stellar_variability_reference_label', None)
                    or (vsp_params[0].get('cname') if vsp_params else None)
                    or 'selected comparison reference'
                )
                f.write(f"# DIFFERENTIAL_MAGNITUDE_REFERENCE={reference_label}\n")
                f.write("# DIFFERENTIAL_MAGNITUDE_AIRMASS_CORRECTED=NO\n")
                f.write(
                    "# BJD_TDB,Apparent Magnitude,Apparent Magnitude Uncertainty,"
                    "Differential Magnitude,Differential Magnitude Uncertainty,Band,Airmass\n"
                )
                for vsp_p in vsp_params:
                    time_value = finite_float(vsp_p.get('time'))
                    mag_value = format_magnitude(vsp_p.get('mag'), default=None)
                    mag_error = format_magnitude_error(vsp_p.get('mag_err'), default=None)
                    if not np.isfinite(time_value) or mag_value is None or mag_error is None:
                        continue
                    differential_mag, differential_error = differential_magnitude_from_vsp_param(
                        vsp_p
                    )
                    differential_mag_text = format_magnitude(
                        differential_mag,
                        default="na",
                        digits=MAGNITUDE_DECIMAL_PLACES,
                    )
                    differential_error_text = format_magnitude_error(
                        differential_error,
                        default="na",
                        digits=MAGNITUDE_DECIMAL_PLACES,
                    )
                    band = vsp_p.get('mag_band') or self.i_dict.get('filter') or 'V'
                    airmass = finite_float(vsp_p.get('airmass'))
                    airmass_text = f"{airmass}" if np.isfinite(airmass) else "na"
                    f.write(
                        f"{time_value}, {mag_value}, {mag_error}, {differential_mag_text}, "
                        f"{differential_error_text}, {band}, {airmass_text}\n"
                    )
            return

        magnitude_series = magnitude_series_from_fit(self.fit)
        band = magnitude_series['band'] or self.i_dict.get('filter') or 'na'

        with params_file.open('w') as f:
            f.write(f"# FINAL TIMESERIES OF {self.p_dict['pName']}\n")
            reference_label = (
                getattr(self.fit, 'differential_magnitude_reference_label', None)
                or getattr(self.fit, 'stellar_variability_reference_label', None)
                or 'selected comparison reference'
            )
            f.write(f"# DIFFERENTIAL_MAGNITUDE_REFERENCE={reference_label}\n")
            f.write(
                "# DIFFERENTIAL_MAGNITUDE_AIRMASS_CORRECTED="
                f"{'YES' if magnitude_series['airmass_corrected'] else 'NO'}\n"
            )
            f.write(
                "# BJD_TDB,Orbital Phase,Flux,Uncertainty,Model,Airmass,"
                "Differential Magnitude,Differential Magnitude Uncertainty,"
                "Apparent Magnitude,Apparent Magnitude Uncertainty,Band\n"
            )

            for row_index, (bjd, phase, flux, fluxerr, model, am) in enumerate(zip(
                    self.fit.time,
                    phase,
                    self.fit.detrended,
                    self.fit.dataerr / self.fit.airmass_model,
                    self.fit.transit,
                    self.fit.airmass_model)):
                row = f"{bjd}, {phase}, {flux}, {fluxerr}, {model}, {am}"
                differential_mag_text = format_magnitude(
                    magnitude_series['differential_magnitude'][row_index],
                    default="na",
                    digits=MAGNITUDE_DECIMAL_PLACES,
                )
                differential_error_text = format_magnitude_error(
                    magnitude_series['differential_magnitude_error'][row_index],
                    default="na",
                    digits=MAGNITUDE_DECIMAL_PLACES,
                )
                apparent_mag_text = format_magnitude(
                    magnitude_series['apparent_magnitude'][row_index],
                    default="na",
                )
                apparent_error_text = format_magnitude_error(
                    magnitude_series['apparent_magnitude_error'][row_index],
                    default="na",
                )
                row = (
                    f"{row}, {differential_mag_text}, {differential_error_text}, "
                    f"{apparent_mag_text}, {apparent_error_text}, {band}"
                )
                f.write(f"{row}\n")

    def differential_magnitude(self):
        target_name = self.p_dict.get('sName') or self.p_dict.get('pName') or 'target'
        return write_differential_magnitude_csv(
            self.fit,
            self.dir,
            target_name,
            observation_date=self.i_dict.get('date'),
            observed_filter=(
                self.i_dict.get('observed_filter')
                or self.i_dict.get('filter')
            ),
            out_of_transit_only=False,
        )

    def stellar_variability_differential_magnitude(self):
        """Write the raw, out-of-transit stellar-variability counterpart."""
        if getattr(self.fit, 'stellar_variability_only', False):
            return None
        fit_shape = np.asarray(getattr(self.fit, 'data', []), dtype=float).shape
        target_shape = np.asarray(
            getattr(self.fit, 'stellar_variability_target_flux', []),
            dtype=float,
        ).shape
        reference_shape = np.asarray(
            getattr(self.fit, 'stellar_variability_comp_flux', []),
            dtype=float,
        ).shape
        if target_shape != fit_shape or reference_shape != fit_shape:
            return None
        target_name = self.p_dict.get('sName') or self.p_dict.get('pName') or 'target'
        return write_differential_magnitude_csv(
            self.fit,
            self.dir,
            target_name,
            observation_date=self.i_dict.get('date'),
            observed_filter=(
                self.i_dict.get('observed_filter')
                or self.i_dict.get('filter')
            ),
            out_of_transit_only=True,
            apply_airmass_correction=False,
            filename_prefix='StellarVariabilityDifferentialMagnitude',
        )

    def final_planetary_params(self, phot_opt, vsp_params, comp_star=None, comp_coords=None, min_aper=None,
                               min_annul=None, adaptive_summary=None, photometry_info=None,
                               publish_to_root=False):
        params_file = self.dir / "working_artifacts" / safe_output_filename(
            "FinalParams",
            self.p_dict['pName'],
            filename_date_token(self.i_dict['date']),
            extension="json",
        )

        if getattr(self.fit, 'stellar_variability_only', False):
            exclusion = getattr(self.fit, 'stellar_variability_transit_exclusion', {}) or {}
            scatter = getattr(self.fit, 'stellar_variability_scatter', np.nan)
            params_num = {
                "Analysis Mode": "Stellar variability only",
                "Transit model fitting": "Skipped",
                "Out-of-transit lightcurve point count": str(len(getattr(self.fit, 'time', []))),
                "Predicted in-transit points excluded": str(exclusion.get('rejected_point_count', 0)),
            }
            duration = exclusion.get('duration_days', np.nan)
            if np.isfinite(duration):
                params_num["Excluded transit-window duration (day)"] = f"{duration:.8f}"
            if np.isfinite(scatter):
                params_num["Residual scatter around flat stellar-variability model"] = f"{scatter * 100.0:.4f} %"
            note = exclusion.get('note')
            if note:
                params_num["Transit-window exclusion note"] = str(note)
            if getattr(self.fit, 'airmass_fit_skipped', False):
                params_num["Airmass correction"] = getattr(
                    self.fit,
                    'airmass_correction_note',
                    "Skipped; no airmass correction applied.",
                )
            if isinstance(photometry_info, dict) and photometry_info.get('noise_budget_summary'):
                params_num["Photometry noise budget"] = str(photometry_info.get('noise_budget_summary'))
                if photometry_info.get('noise_budget_terms'):
                    params_num["Photometry noise budget terms"] = ", ".join(
                        str(term) for term in photometry_info.get('noise_budget_terms')
                    )

            if vsp_params:
                params_num["Variable Reference Star"] = stellar_variability_reference_summary(vsp_params[0])
                if vsp_params[0].get('ensemble_reference'):
                    params_num["Variable Reference Measurement"] = (
                        f"Combined {len(vsp_params)} out-of-transit target measurements against the "
                        "calibrated comparison-star ensemble; AID rows list the BJD_TDB timestamps used."
                    )
                else:
                    params_num["Variable Reference Measurement"] = (
                        f"Remeasured {len(vsp_params)} out-of-transit target/reference point(s) "
                        "against the stellar-variability reference catalog star; AID rows list the "
                        "BJD_TDB timestamps used."
                    )

            if phot_opt:
                if comp_star == 'ensemble':
                    reference_text = "ensemble"
                else:
                    reference_text = (
                        f"#{comp_star} - {comp_coords}"
                        if comp_star is not None and min_aper is not None and min_aper >= 0
                        else str(comp_star)
                    )
                params_num["Stellar Variability Reference Star"] = reference_text
                if min_aper == 0:
                    params_num["Optimal Method"] = "PSF photometry"
                else:
                    if adaptive_summary:
                        params_num["Adaptive Aperture Scale"] = f"{adaptive_summary['aperture_sigma']:.2f} sigma"
                        params_num["Adaptive Annulus Scale"] = f"{adaptive_summary['annulus_sigma']:.2f} sigma"
                        params_num["Optimal Aperture"] = (
                            f"{adaptive_summary['aperture_median']:.2f} +/- "
                            f"{adaptive_summary['aperture_std']:.2f} px"
                        )
                        params_num["Aperture Range"] = (
                            f"{adaptive_summary['aperture_min']:.2f} to "
                            f"{adaptive_summary['aperture_max']:.2f} px"
                        )
                        params_num["Optimal Annulus"] = (
                            f"{adaptive_summary['annulus_median']:.2f} +/- "
                            f"{adaptive_summary['annulus_std']:.2f} px"
                        )
                        params_num["Annulus Range"] = (
                            f"{adaptive_summary['annulus_min']:.2f} to "
                            f"{adaptive_summary['annulus_max']:.2f} px"
                        )
                    else:
                        params_num["Optimal Aperture"] = f"{abs(min_aper)}"
                        params_num["Optimal Annulus"] = f"{min_annul}"

            final_params = {'FINAL STELLAR VARIABILITY PARAMETERS': params_num}
            with params_file.open('w') as f:
                dump(final_params, f, indent=4)
            if publish_to_root:
                root_params_file = self.dir / params_file.name
                if root_params_file != params_file:
                    root_params_file.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(params_file, root_params_file)
            return

        transit_qc = getattr(self.fit, 'transit_qc', None)
        fit_quality = build_fit_quality_metadata(self.fit)
        empirical_uncertainty = fit_empirical_transit_uncertainty(self.fit, fit_quality=fit_quality)
        qc_residual_scatter = np.nan
        if isinstance(transit_qc, dict):
            qc_residual_scatter = transit_qc.get('residual_scatter', np.nan)
        if not np.isfinite(qc_residual_scatter):
            residuals = np.asarray(getattr(self.fit, 'residuals', np.array([])), dtype=float)
            data = np.asarray(getattr(self.fit, 'data', np.array([])), dtype=float)
            if residuals.size and data.size:
                if residuals.shape == data.shape:
                    median_flux = np.nanmedian(data)
                    if np.isfinite(median_flux) and median_flux != 0:
                        qc_residual_scatter = float(np.std(residuals) / median_flux)
                elif residuals.size == 1:
                    median_flux = np.nanmedian(data)
                    if np.isfinite(median_flux) and median_flux != 0:
                        qc_residual_scatter = float(abs(residuals.reshape(-1)[0]) / median_flux)

        headline_params = format_transit_qc_headline_final_params(transit_qc)
        rprs_report_error = fit_rprs_report_error(
            self.fit,
            empirical_uncertainty=empirical_uncertainty,
        )
        if not np.isfinite(rprs_report_error) or rprs_report_error < 0:
            rprs_report_error = finite_float(self.p_dict.get('rprsUnc'))
        tmid_report_error = fit_parameter_model_data_uncertainty(
            self.fit,
            'tmid',
            empirical_uncertainty=empirical_uncertainty,
        )
        if not np.isfinite(tmid_report_error) or tmid_report_error < 0:
            tmid_report_error = self.fit.errors['tmid']
        inc_report_error = fit_parameter_model_data_uncertainty(
            self.fit,
            'inc',
            empirical_uncertainty=empirical_uncertainty,
        )
        if not np.isfinite(inc_report_error) or inc_report_error < 0:
            inc_report_error = self.fit.errors['inc']
        ars_report_error = fit_parameter_model_data_uncertainty(
            self.fit,
            'ars',
            empirical_uncertainty=empirical_uncertainty,
        )
        if not np.isfinite(ars_report_error) or ars_report_error < 0:
            ars_report_error = self.fit.errors.get('ars', np.nan)
        core_params = {
            "Mid-Transit Time (Tmid)": f"{round_to_2(self.fit.parameters['tmid'], tmid_report_error)} +/- "
                                       f"{round_to_2(tmid_report_error)} BJD_TDB",
            "Ratio of Planet to Stellar Radius (Rp/R*)": f"{round_to_2(self.fit.parameters['rprs'], rprs_report_error)} +/- "
                                                         f"{round_to_2(rprs_report_error)}",
            "Orbital Inclination (inc)": f"{round_to_2(self.fit.parameters['inc'], inc_report_error)} +/- "
                                                    f"{round_to_2(inc_report_error)} ",
        }
        depth_params = formatted_transit_depth_parameters(
            self.fit,
            self.p_dict,
            empirical_uncertainty=empirical_uncertainty,
        )
        params_num = {
            **headline_params,
            "Mid-Transit Time (Tmid)": core_params["Mid-Transit Time (Tmid)"],
            "Ratio of Planet to Stellar Radius (Rp/R*)": core_params["Ratio of Planet to Stellar Radius (Rp/R*)"],
            **depth_params,
            "Orbital Inclination (inc)": core_params["Orbital Inclination (inc)"],
        }
        ars_text = format_parameter_with_error(
            self.fit.parameters.get('ars'),
            ars_report_error,
        )
        if ars_text is not None:
            params_num["Ratio of Distance to Stellar Radius (a/Rs)"] = ars_text
        model_impact_parameter, model_impact_error = fit_impact_parameter_value_error(self.fit)
        impact_error_overrides = {
            'ars': ars_report_error,
            'inc': inc_report_error,
        }
        if (
            np.isfinite(model_impact_error)
            and ('b' in (getattr(self.fit, 'parameters', {}) or {})
                 or 'b' in (getattr(self.fit, 'sample_parameters', {}) or {}))
        ):
            impact_error_overrides['b'] = model_impact_error * empirical_red_noise_error_scale(
                empirical_uncertainty
            )
        impact_parameter, impact_error = fit_impact_parameter_value_error(
            self.fit,
            errors_override=impact_error_overrides,
        )
        impact_text = format_parameter_with_error(impact_parameter, impact_error)
        if impact_text is not None:
            params_num["Impact Parameter (b)"] = impact_text
        if empirical_uncertainty.get('available'):
            tmid_model_text = (
                f"{round_to_2(self.fit.parameters['tmid'], self.fit.errors['tmid'])} +/- "
                f"{round_to_2(self.fit.errors['tmid'])} BJD_TDB"
            )
            tmid_combined_text = core_params["Mid-Transit Time (Tmid)"]
            params_num["Mid-Transit Time (Tmid) model-fit uncertainty"] = tmid_model_text
            params_num["Mid-Transit Time (Tmid) model+red-noise uncertainty"] = tmid_combined_text

            inc_model_text = (
                f"{round_to_2(self.fit.parameters['inc'], self.fit.errors['inc'])} +/- "
                f"{round_to_2(self.fit.errors['inc'])} "
            )
            params_num["Orbital Inclination (inc) model-fit uncertainty"] = inc_model_text
            params_num["Orbital Inclination (inc) model+red-noise uncertainty"] = core_params[
                "Orbital Inclination (inc)"
            ]

            ars_model_text = format_parameter_with_error(
                self.fit.parameters.get('ars'),
                self.fit.errors.get('ars'),
            )
            if ars_model_text is not None:
                params_num["Ratio of Distance to Stellar Radius (a/Rs) model-fit uncertainty"] = ars_model_text
            if ars_text is not None:
                params_num["Ratio of Distance to Stellar Radius (a/Rs) model+red-noise uncertainty"] = ars_text

            impact_model_text = format_parameter_with_error(model_impact_parameter, model_impact_error)
            if impact_model_text is not None:
                params_num["Impact Parameter (b) model-fit uncertainty"] = impact_model_text
            if impact_text is not None:
                params_num["Impact Parameter (b) model+red-noise uncertainty"] = impact_text
        if getattr(self.fit, 'ns_type', None) is not None:
            params_num["Fit parameter point estimate"] = (
                "Best-fit likelihood point; uncertainties are posterior spread."
            )
        prefit_refinement_note = getattr(self.fit, 'prefit_refinement_note', None)
        if prefit_refinement_note:
            params_num["Prefit refinement note"] = str(prefit_refinement_note)
        geometry_prior_note = getattr(self.fit, 'partial_transit_geometry_prior_assumption_note', None)
        if geometry_prior_note:
            params_num["Prior-assumed partial-transit geometry note"] = str(geometry_prior_note)
        ars_prior_fallback_note = getattr(self.fit, 'ars_prior_fallback_note', None)
        if ars_prior_fallback_note:
            params_num["a/Rs prior fallback note"] = str(ars_prior_fallback_note)
        oot_baseline_parameter_note = getattr(self.fit, 'oot_baseline_parameter_fit_note', None)
        if oot_baseline_parameter_note:
            params_num["Out-of-transit baseline parameter-fit note"] = str(oot_baseline_parameter_note)
        oot_baseline_note = getattr(self.fit, 'oot_baseline_detrending_note', None)
        if oot_baseline_note:
            params_num["Out-of-transit baseline detrending note"] = str(oot_baseline_note)
        sparse_posterior_note = getattr(self.fit, 'sparse_posterior_live_point_extension_note', None)
        if sparse_posterior_note:
            params_num["Sparse posterior live-point extension note"] = str(sparse_posterior_note)
        final_residual_note = getattr(self.fit, 'final_residual_rejection_note', None)
        if final_residual_note:
            params_num["Final residual rejection note"] = str(final_residual_note)
        if np.isfinite(qc_residual_scatter):
            params_num["Residual scatter around full model fit"] = f"{qc_residual_scatter * 100.0:.4f} %"
        params_num.update(format_fit_quality_final_params(fit_quality))
        params_num.update(format_empirical_transit_uncertainty_final_params(empirical_uncertainty))
        params_num.update(format_ktmf_decision_final_params(self.fit, photometry_info))
        if isinstance(photometry_info, dict) and photometry_info.get('noise_budget_summary'):
            params_num["Photometry noise budget"] = str(photometry_info.get('noise_budget_summary'))
            if photometry_info.get('noise_budget_terms'):
                params_num["Photometry noise budget terms"] = ", ".join(
                    str(term) for term in photometry_info.get('noise_budget_terms')
                )
        if getattr(self.fit, 'airmass_fit_skipped', False):
            params_num["Airmass correction"] = getattr(
                self.fit,
                'airmass_correction_note',
                "Skipped; no airmass correction applied.",
            )
        else:
            if 'a0' in self.fit.parameters:
                a0_error = self.fit.errors.get('a0') if isinstance(self.fit.errors, dict) else None
                if a0_error is not None and np.isfinite(a0_error):
                    params_num["Baseline flux (a0)"] = (
                        f"{round_to_2(self.fit.parameters['a0'], a0_error)} +/- "
                        f"{round_to_2(a0_error)}"
                    )
                else:
                    params_num["Baseline flux (a0)"] = (
                        f"{round_to_2(self.fit.parameters['a0'])} (fixed; uncertainty unavailable)"
                    )
            else:
                a1_error = self.fit.errors.get('a1') if isinstance(self.fit.errors, dict) else None
                if a1_error is not None and np.isfinite(a1_error):
                    params_num["Flux normalization (a1)"] = (
                        f"{round_to_2(self.fit.parameters['a1'], a1_error)} +/- "
                        f"{round_to_2(a1_error)}"
                    )
                else:
                    params_num["Flux normalization (a1)"] = (
                        f"{round_to_2(self.fit.parameters['a1'])} (fixed; uncertainty unavailable)"
                    )
            if 'a2' in self.fit.parameters:
                a2_error = self.fit.errors.get('a2') if isinstance(self.fit.errors, dict) else None
                if a2_error is not None and np.isfinite(a2_error):
                    params_num["Airmass coefficient 2 (a2)"] = (
                        f"{round_to_2(self.fit.parameters['a2'], a2_error)} +/- "
                        f"{round_to_2(a2_error)}"
                    )
                else:
                    params_num["Airmass coefficient 2 (a2)"] = (
                        f"{round_to_2(self.fit.parameters['a2'])} (fixed; uncertainty unavailable)"
                    )

        if isinstance(transit_qc, dict) and transit_qc:
            qc_status = transit_qc.get('status')
            qc_summary = transit_qc.get('summary')
            qc_notes = transit_qc.get('notes') or []
            qc_delta_bic = transit_qc.get('delta_bic', np.nan)
            qc_delta_chi2 = transit_qc.get('delta_chi2', np.nan)
            qc_rprs_sigma = transit_qc.get('rprs_sigma', np.nan)
            qc_duration_ratio = transit_qc.get('duration_ratio', np.nan)
            qc_eebls_depth_snr = transit_qc.get('eebls_depth_snr', np.nan)
            qc_deviation_metric = transit_qc.get('deviation_from_expected_value', np.nan)
            qc_rprs_deviation_sigma = transit_qc.get('rprs_deviation_sigma', np.nan)
            qc_rprs_deviation_fit_unc = transit_qc.get('rprs_deviation_fit_unc', np.nan)
            qc_rprs_deviation_model_unc = transit_qc.get('rprs_deviation_model_fit_unc', np.nan)
            qc_rprs_deviation_data_unc = transit_qc.get('rprs_deviation_data_fit_unc', np.nan)
            qc_rprs_deviation_expected_unc = transit_qc.get('rprs_deviation_expected_unc', np.nan)
            qc_rprs_deviation_comparison_unc = transit_qc.get('rprs_deviation_unc', np.nan)
            qc_sigma_threshold = transit_qc.get('deviation_sigma_threshold', np.nan)
            qc_ktmf = transit_qc.get('ktmf_metric', np.nan)
            qc_ktmf_contributions = transit_qc.get('ktmf_contributions') or []

            if qc_status:
                params_num["Transit detection QC"] = str(qc_status).upper()
            if qc_summary:
                params_num["Transit vs flat model"] = qc_summary
            if np.isfinite(qc_delta_bic):
                params_num["Transit vs flat Delta BIC"] = f"{qc_delta_bic:.2f}"
            if np.isfinite(qc_delta_chi2):
                params_num["Transit vs flat Delta chi2"] = f"{qc_delta_chi2:.2f}"
            if np.isfinite(qc_rprs_sigma):
                params_num["Transit depth significance"] = f"{qc_rprs_sigma:.2f} sigma"
            if np.isfinite(qc_duration_ratio):
                params_num["Transit duration consistency"] = f"{qc_duration_ratio:.2f}x modeled duration"
            if np.isfinite(qc_eebls_depth_snr):
                params_num["EEBLS depth SNR"] = f"{qc_eebls_depth_snr:.2f}"
            if np.isfinite(qc_deviation_metric):
                params_num["Deviation From Expected Value"] = f"{qc_deviation_metric:.2f} / 1.00"
            if np.isfinite(qc_sigma_threshold):
                params_num["Expected-value QC threshold"] = f"{qc_sigma_threshold:.2f} sigma"
            if np.isfinite(qc_rprs_deviation_sigma):
                params_num["Expected-value Rp/R* deviation"] = f"{qc_rprs_deviation_sigma:.2f} sigma"
            if np.isfinite(qc_rprs_deviation_fit_unc):
                params_num["Expected-value Rp/R* fit uncertainty used"] = (
                    f"+/- {round_to_2(qc_rprs_deviation_fit_unc)}"
                )
            if np.isfinite(qc_rprs_deviation_model_unc):
                params_num["Expected-value Rp/R* model-fit uncertainty"] = (
                    f"+/- {round_to_2(qc_rprs_deviation_model_unc)}"
                )
            if np.isfinite(qc_rprs_deviation_data_unc):
                params_num["Expected-value Rp/R* data-fit red-noise uncertainty"] = (
                    f"+/- {round_to_2(qc_rprs_deviation_data_unc)}"
                )
            if np.isfinite(qc_rprs_deviation_expected_unc):
                params_num["Expected-value Rp/R* prior uncertainty"] = (
                    f"+/- {round_to_2(qc_rprs_deviation_expected_unc)}"
                )
            if np.isfinite(qc_rprs_deviation_comparison_unc):
                params_num["Expected-value Rp/R* total comparison uncertainty"] = (
                    f"+/- {round_to_2(qc_rprs_deviation_comparison_unc)}"
                )
            if np.isfinite(qc_ktmf):
                params_num["KTMF"] = f"{qc_ktmf:.2f} / 5.00"
            for contribution_index, contribution in enumerate(qc_ktmf_contributions, start=1):
                label = contribution.get('label', f'Component {contribution_index}')
                detail = contribution.get('detail') or 'n/a'
                available = bool(contribution.get('available'))
                points = float(contribution.get('points', 0.0) or 0.0)
                max_points = float(contribution.get('max_points', 0.0) or 0.0)
                score = contribution.get('score', np.nan)
                if available and np.isfinite(score):
                    params_num[f"KTMF contribution {contribution_index}"] = (
                        f"{label}: +{points:.2f}/{max_points:.2f} (score={score:.2f}; {detail})"
                    )
                else:
                    params_num[f"KTMF contribution {contribution_index}"] = (
                        f"{label}: +0.00/0.00 (unavailable; {detail})"
                    )
            if qc_notes:
                params_num["Transit QC notes"] = " ".join(str(note) for note in qc_notes)

        if vsp_params and not (phot_opt and comp_star is None):
            params_num["Variable Reference Star"] = stellar_variability_reference_summary(vsp_params[0])
            measurement_summary = stellar_variability_measurement_summary(vsp_params, comp_star)
            if measurement_summary:
                params_num["Variable Reference Measurement"] = measurement_summary

        if phot_opt:
            if comp_star == 'ensemble':
                transit_fit_comp_text = "ensemble"
            else:
                transit_fit_comp_text = (
                    f"#{comp_star} - {comp_coords}"
                    if comp_star is not None and min_aper >= 0
                    else str(comp_star)
                )
            phot_ext = {
                "Transit Fit Comparison Star": transit_fit_comp_text
            }
            if min_aper == 0:
                phot_ext["Optimal Method"] = "PSF photometry"
            else:
                if adaptive_summary:
                    phot_ext["Adaptive Aperture Scale"] = f"{adaptive_summary['aperture_sigma']:.2f} sigma"
                    phot_ext["Adaptive Annulus Scale"] = f"{adaptive_summary['annulus_sigma']:.2f} sigma"
                    phot_ext["Optimal Aperture"] = (
                        f"{adaptive_summary['aperture_median']:.2f} +/- {adaptive_summary['aperture_std']:.2f} px"
                    )
                    phot_ext["Aperture Range"] = (
                        f"{adaptive_summary['aperture_min']:.2f} to {adaptive_summary['aperture_max']:.2f} px"
                    )
                    phot_ext["Optimal Annulus"] = (
                        f"{adaptive_summary['annulus_median']:.2f} +/- {adaptive_summary['annulus_std']:.2f} px"
                    )
                    phot_ext["Annulus Range"] = (
                        f"{adaptive_summary['annulus_min']:.2f} to {adaptive_summary['annulus_max']:.2f} px"
                    )
                else:
                    phot_ext["Optimal Aperture"] = f"{abs(min_aper)}"
                    phot_ext["Optimal Annulus"] = f"{min_annul}"
            params_num.update(phot_ext)

        params_num["Transit Duration (day)"] = (f"{round_to_2(mean(self.durs), std(self.durs))} +/- "
                                                f"{round_to_2(std(self.durs))}")
        final_params = {'FINAL PLANETARY PARAMETERS': params_num}

        with params_file.open('w') as f:
            dump(final_params, f, indent=4)
        if publish_to_root:
            root_params_file = self.dir / params_file.name
            if root_params_file != params_file:
                root_params_file.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(params_file, root_params_file)

    def aavso(self, comp_star, airmasses, ld0, ld1, ld2, ld3, epw_md5,
              photometry_info=None, astrometry_info=None, frame_filtering_info=None,
              bad_pixel_info=None):
        priors_dict, filter_dict, results_dict = aavso_dicts(self.p_dict, self.fit, self.i_dict, self.durs,
                                                             ld0, ld1, ld2, ld3)
        aavso_airmass_terms = aavso_airmass_results(self.fit)
        detrend_model = aavso_detrend_model(self.fit)
        qc_metadata = build_aavso_qc_metadata(self.fit)
        fit_quality_metadata = build_fit_quality_metadata(self.fit)
        rprs_report_error = fit_rprs_report_error(self.fit)
        if not np.isfinite(rprs_report_error) or rprs_report_error < 0:
            rprs_report_error = finite_float(self.p_dict.get('rprsUnc'))
        ktmf_decision_metadata = build_ktmf_decision_metadata(self.fit, photometry_info)
        photometry_metadata = build_aavso_photometry_metadata(photometry_info)
        aperture_metadata = build_aavso_aperture_metadata(photometry_info)
        frame_filtering_metadata = build_aavso_frame_filtering_metadata(self.fit, frame_filtering_info)
        astrometry_metadata = build_aavso_astrometry_metadata(astrometry_info, comp_star)
        bad_pixel_metadata = build_aavso_bad_pixel_metadata(bad_pixel_info)
        magnitude_series = magnitude_series_from_fit(self.fit)
        obs_name = format_aavso_header_value(self.i_dict.get('obs_name'))
        obs_name_header = f"#OBSNAME={obs_name}\n" if obs_name else ""
        gaia_dist = format_aavso_header_value(self.p_dict.get('dist'))
        gaia_pmra = format_aavso_header_value(self.p_dict.get('pm_ra'))
        gaia_pmdec = format_aavso_header_value(self.p_dict.get('pm_dec'))
        gaia_dist_header = f"#GAIADIST={gaia_dist}\n" if gaia_dist else ""
        gaia_pmra_header = f"#GAIAPMRA={gaia_pmra}\n" if gaia_pmra else ""
        gaia_pmdec_header = f"#GAIAPMDEC={gaia_pmdec}\n" if gaia_pmdec else ""

        params_file = aavso_output_directory(self.dir) / safe_output_filename(
            "AAVSO",
            self.p_dict['pName'],
            filename_date_token(self.i_dict['date']),
            extension="txt",
        )

        with params_file.open('w', encoding="utf-8") as f:
            f.write("#TYPE=EXOPLANET\n"  # fixed
                    f"#OBSCODE={self.i_dict['aavso_num']}\n"  # UI
                    f"#SECONDARY_OBSCODES={self.i_dict['second_obs']}\n"  # UI
                    f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                    "#DELIM=,\n"  # fixed
                    "#DATE_TYPE=BJD_TDB\n"  # fixed
                    f"#OBSDATE={format_aavso_header_value(self.i_dict.get('date'))}\n"
                    f"{obs_name_header}"
                    f"#OBSTYPE={self.i_dict['camera']}\n"
                    f"#STAR_NAME={self.p_dict['sName']}\n"  # code yields
                    f"#EXOPLANET_NAME={self.p_dict['pName']}\n"  # code yields
                    f"#BINNING={self.i_dict['pixel_bin']}\n"  # user input
                    f"#EXPOSURE_TIME={self.i_dict.get('exposure', -1)}\n"  # UI
                    f"#OBSLAT={format_aavso_header_value(self.i_dict.get('lat'))}\n"
                    f"#OBSLON={format_aavso_header_value(self.i_dict.get('long'))}\n"
                    f"#OBSELEV={format_aavso_header_value(self.i_dict.get('elev'))}\n"
                    f"{gaia_dist_header}"
                    f"{gaia_pmra_header}"
                    f"{gaia_pmdec_header}"
                    f"#COMP_STAR-XC={dumps(comp_star)}\n"
                    f"#NOTES={self.i_dict['notes']}\n"
                    "#DETREND_PARAMETERS=AIRMASS, AIRMASS CORRECTION FUNCTION\n"  # fixed
                    "#MEASUREMENT_TYPE=Rnflux\n"  # fixed
                    f"#FILTER={self.i_dict['filter']}\n"
                    f"#FILTER-XC={dumps(filter_dict)}\n"
                    f"#PRIORS=Period={round_to_2(self.p_dict['pPer'], self.p_dict['pPerUnc'])} +/- {round_to_2(self.p_dict['pPerUnc'])}"
                    f",Rp/R*={round_to_2(self.p_dict['rprs'], self.p_dict['rprsUnc'])} +/- {round_to_2(self.p_dict['rprsUnc'])}"
                    f",a/R*={round_to_2(self.p_dict['aRs'], self.p_dict['aRsUnc'])} +/- {round_to_2(self.p_dict['aRsUnc'])}"
                    f",inc={round_to_2(self.p_dict['inc'], self.p_dict['incUnc'])} +/- {round_to_2(self.p_dict['incUnc'])}"
                    f",ecc={round_to_2(self.p_dict['ecc'])}"
                    f",u0={round_to_2(ld0[0], ld0[1])} +/- {round_to_2(ld0[1])}"
                    f",u1={round_to_2(ld1[0], ld1[1])} +/- {round_to_2(ld1[1])}"
                    f",u2={round_to_2(ld2[0], ld2[1])} +/- {round_to_2(ld2[1])}"
                    f",u3={round_to_2(ld3[0], ld3[1])} +/- {round_to_2(ld3[1])}\n"
                    f"#PRIORS-XC={dumps(priors_dict)}\n"  # code yields
                    f"#RESULTS=Tc={round_to_2(self.fit.parameters['tmid'], self.fit.errors['tmid'])} +/- {round_to_2(self.fit.errors['tmid'])}"
                    f",Rp/R*={round_to_2(self.fit.parameters['rprs'], rprs_report_error)} +/- {round_to_2(rprs_report_error)}"
                    f",inc={round_to_2(self.fit.parameters['inc'], self.fit.errors['inc'])} +/- {round_to_2(self.fit.errors['inc'])}"
                    f",{aavso_airmass_terms[0][0]}={aavso_airmass_terms[0][1]} +/- {aavso_airmass_terms[0][2]}"
                    f",{aavso_airmass_terms[1][0]}={aavso_airmass_terms[1][1]} +/- {aavso_airmass_terms[1][2]}\n"
                    f"#RESULTS-XC={dumps(results_dict)}\n")  # code yields
            f.write(format_aavso_json_header("QC-XC", qc_metadata))
            f.write(format_aavso_json_header("FIT_QUALITY-XC", fit_quality_metadata))
            f.write(format_aavso_json_header("KTMF_DECISION-XC", ktmf_decision_metadata))
            f.write(format_aavso_json_header("PHOTOMETRY-XC", photometry_metadata))
            f.write(format_aavso_json_header("APERTURE-XC", aperture_metadata))
            f.write(format_aavso_json_header("FRAME_FILTERING-XC", frame_filtering_metadata))
            f.write(format_aavso_json_header("ASTROMETRY-XC", astrometry_metadata))
            f.write(format_aavso_json_header("BAD_PIXEL-XC", bad_pixel_metadata))
            f.write(format_aavso_json_header("MAGNITUDE_FIELDS-XC", {
                'apparent_magnitude': 'catalogue-calibrated target magnitude',
                'apparent_magnitude_error': 'flux and catalogue calibration uncertainty',
                'differential_magnitude': (
                    'target minus selected comparison reference; '
                    '-2.5 log10(target_flux/reference_flux)'
                ),
                'differential_magnitude_error': 'flux-only uncertainty',
                'per_point_header': 'MAGNITUDE-XC',
                'band': magnitude_series['band'] or self.i_dict.get('filter'),
                'airmass_corrected': magnitude_series['airmass_corrected'],
                'apparent_calibrated': magnitude_series['apparent_calibrated'],
                'comparison_reference': (
                    getattr(self.fit, 'differential_magnitude_reference_label', None)
                    or getattr(self.fit, 'stellar_variability_reference_label', None)
                    or 'selected comparison reference'
                ),
            }))
            for magnitude_index in range(0, len(self.fit.time)):
                f.write(format_aavso_json_header("MAGNITUDE-XC", {
                    'date_bjd_tdb': finite_float(self.fit.time[magnitude_index]),
                    'differential_magnitude': rounded_magnitude_value(
                        magnitude_series['differential_magnitude'][magnitude_index],
                        digits=MAGNITUDE_DECIMAL_PLACES,
                    ),
                    'differential_magnitude_error': rounded_magnitude_error(
                        magnitude_series['differential_magnitude_error'][magnitude_index],
                        digits=MAGNITUDE_DECIMAL_PLACES,
                    ),
                    'apparent_magnitude': rounded_magnitude_value(
                        magnitude_series['apparent_magnitude'][magnitude_index],
                        digits=MAGNITUDE_DECIMAL_PLACES,
                    ),
                    'apparent_magnitude_error': rounded_magnitude_error(
                        magnitude_series['apparent_magnitude_error'][magnitude_index],
                        digits=MAGNITUDE_DECIMAL_PLACES,
                    ),
                    'band': magnitude_series['band'] or self.i_dict.get('filter'),
                }, preserve_nulls=True))

            if epw_md5:
                f.write(f"#EPW_MD5-XC={dumps({'epw_checkout_md5': epw_md5})}\n")

            f.write(
                "# EXOTIC is developed by Exoplanet Watch (exoplanets.nasa.gov/exoplanet-watch/), a citizen science "
                "project managed by NASA's Jet Propulsion Laboratory on behalf of NASA's Universe of Learning. "
                "This work is supported by NASA under award number NNX16AC65A to the "
                "Space Telescope Science Institute.\n"
                "# Use of this data is governed by the AAVSO Data Usage Guidelines: "
                "aavso.org/data-usage-guidelines\n")

            f.write("#DATE,DIFF,ERR,DETREND_1,DETREND_2\n")
            for aavsoC in range(0, len(self.fit.time)):
                # f.write(f"{round(self.fit.time[aavsoC], 8)},{round(self.fit.data[aavsoC] / self.fit.parameters['a1'], 7)},"
                #         f"{round(self.fit.dataerr[aavsoC] / self.fit.parameters['a1'], 7)},{round(airmasses[aavsoC], 7)},"
                #         f"{round(self.fit.airmass_model[aavsoC] / self.fit.parameters['a1'], 7)}\n")
                f.write(f"{round(self.fit.time[aavsoC], 8)},{round(self.fit.data[aavsoC], 7)},"
                        f"{round(self.fit.dataerr[aavsoC], 7)},{round(airmasses[aavsoC], 7)},"
                        f"{round(detrend_model[aavsoC], 7)}\n")
        copy_aavso_supporting_artifacts(
            self.dir,
            self.p_dict['pName'],
            self.i_dict['date'],
        )

    def plate_status(self, plate_status: PlateStatus):
        plate_status_file = self.dir / "working_artifacts" / safe_output_filename(
            "PlateStatus",
            self.p_dict['pName'],
            filename_date_token(self.i_dict['date']),
            extension="csv",
        )
        plate_status.writePlateStatus(plate_status_file)

class AIDOutputFiles:
    def __init__(self, fit, p_dict, i_dict, auid, chart_id, vsp_params):
        self.fit = fit
        self.auid = auid
        self.chart_id = chart_id
        self.p_dict = p_dict
        self.i_dict = i_dict
        self.dir = Path(self.i_dict['save'])
        self.vsp_params = vsp_params

    def _aavso_path(self):
        return aavso_output_directory(self.dir) / safe_output_filename(
            "AID_AAVSO",
            self.p_dict['sName'],
            filename_date_token(self.i_dict['date']),
            extension="txt",
        )

    def _write_aavso(self, params_file, use_row_names=False, include_comparison_metadata=True):
        first_vsp_param = self.vsp_params[0] if self.vsp_params else {}
        comparison_metadata = aid_comparison_metadata(first_vsp_param)
        ensemble_comparison_metadata = aid_ensemble_comparison_metadata(first_vsp_param)
        fallback_differential_series = differential_magnitude_series_from_fit(
            self.fit,
            apply_airmass_correction=False,
        )
        fallback_times = np.asarray(
            [] if fallback_differential_series is None else fallback_differential_series['time'],
            dtype=float,
        )
        comparison_coordinate_headers = aid_comparison_coordinate_headers(
            self.vsp_params,
            indexed=use_row_names,
        )
        default_variable_name = self.auid or self.p_dict.get('sName') or self.p_dict.get('pName')

        with params_file.open('w', encoding="utf-8") as f:
            f.write("#TYPE=EXTENDED\n"  # fixed
                    f"#OBSCODE={self.i_dict['aavso_num']}\n"  # UI
                    f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                    "#DELIM=,\n"  # fixed
                    "#DATE=BJD_TDB\n"  # fixed
                    f"#OBSDATE={format_aavso_header_value(self.i_dict.get('date'))}\n"
                    f"#OBSTYPE={self.i_dict['camera']}\n"
                    f"#OBSLAT={format_aavso_header_value(self.i_dict.get('lat'))}\n"
                    f"#OBSLON={format_aavso_header_value(self.i_dict.get('long'))}\n"
                    f"#OBSELEV={format_aavso_header_value(self.i_dict.get('elev'))}\n")
            f.write(
                "# EXOTIC is developed by Exoplanet Watch (exoplanets.nasa.gov/exoplanet-watch/), a citizen science "
                "project managed by NASA's Jet Propulsion Laboratory on behalf of NASA's Universe of Learning. "
                "This work is supported by NASA under award number NNX16AC65A to the "
                "Space Telescope Science Institute.\n"
                "# Use of this data is governed by the AAVSO Data Usage Guidelines: "
                "aavso.org/data-usage-guidelines\n")
            if include_comparison_metadata and comparison_metadata:
                f.write(f"#COMPARISON-CATALOG-XC={dumps(comparison_metadata, sort_keys=True)}\n")
            if comparison_coordinate_headers:
                f.write(comparison_coordinate_headers)
            if include_comparison_metadata and ensemble_comparison_metadata:
                f.write(format_aavso_json_header(
                    "ENSEMBLE-COMPARISONS-XC",
                    ensemble_comparison_metadata,
                ))
            f.write(format_aavso_json_header("MAGNITUDE_FIELDS-XC", {
                'apparent_magnitude': 'MAG',
                'apparent_magnitude_error': 'MERR',
                'differential_magnitude': 'DIFFMAG',
                'differential_magnitude_error': 'DIFFERR',
                'differential_magnitude_definition': (
                    'target minus selected comparison reference; '
                    '-2.5 log10(target_flux/reference_flux)'
                ),
                'airmass_corrected': False,
            }))

            f.write(
                "#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,"
                "AMASS,GROUP,CHART,NOTES,DIFFMAG,DIFFERR\n"
            )
            for vsp_p in self.vsp_params:
                variable_name = default_variable_name
                if use_row_names:
                    variable_name = vsp_p.get('_aid_name') or variable_name
                mag = format_magnitude(
                    vsp_p.get('mag'),
                    default=None,
                    digits=MAGNITUDE_DECIMAL_PLACES,
                )
                if mag is None:
                    continue
                mag_err = format_magnitude_error(
                    vsp_p.get('mag_err'),
                    digits=MAGNITUDE_DECIMAL_PLACES,
                )
                cmag = format_magnitude(vsp_p.get('cmag'))
                chart_id = self.chart_id or vsp_p.get('chart_id') or 'na'
                differential_mag, differential_error = differential_magnitude_from_vsp_param(
                    vsp_p
                )
                if not np.isfinite(differential_mag) and fallback_times.size:
                    row_time = finite_float(vsp_p.get('time'))
                    time_matches = np.flatnonzero(np.isclose(
                        fallback_times,
                        row_time,
                        rtol=0.0,
                        atol=5.0e-5,
                    ))
                    if time_matches.size:
                        fallback_index = time_matches[0]
                        differential_mag = fallback_differential_series['magnitude'][fallback_index]
                        differential_error = (
                            fallback_differential_series['magnitude_error'][fallback_index]
                        )
                differential_mag_text = format_magnitude(
                    differential_mag,
                    default=None,
                    digits=MAGNITUDE_DECIMAL_PLACES,
                )
                differential_error_text = format_magnitude_error(
                    differential_error,
                    default=None,
                    digits=MAGNITUDE_DECIMAL_PLACES,
                )
                differential_mag_text = differential_mag_text or 'na'
                differential_error_text = differential_error_text or 'na'
                f.write(f"{variable_name},{round(vsp_p['time'], 5)},{mag},{mag_err},"
                        f"{self.i_dict['filter']},NO,STD,{vsp_p['cname']},{cmag},na,na,"
                        f"{round(vsp_p['airmass'], 7)},na,{chart_id},na,"
                        f"{differential_mag_text},{differential_error_text}\n")
        return params_file

    def aavso(self):
        params_file = self._write_aavso(self._aavso_path())
        copy_aavso_supporting_artifacts(
            self.dir,
            self.p_dict.get('pName') or self.p_dict.get('sName'),
            self.i_dict['date'],
        )
        return params_file

    def combined_aavso(self):
        """Write one AID file containing rows for multiple named variables."""
        params_file = self._write_aavso(
            self._aavso_path(),
            use_row_names=True,
            include_comparison_metadata=False,
        )
        copy_aavso_supporting_artifacts(
            self.dir,
            self.p_dict.get('pName') or self.p_dict.get('sName'),
            self.i_dict['date'],
        )
        return params_file


def aavso_dicts(planet_dict, fit, info_dict, durs, ld0, ld1, ld2, ld3):
    aavso_airmass_terms = aavso_airmass_results(fit)
    rprs_report_error = fit_rprs_report_error(fit)
    if not np.isfinite(rprs_report_error) or rprs_report_error < 0:
        rprs_report_error = finite_float(planet_dict.get('rprsUnc'))
    priors = {
        'Period': {
            'value': str(round_to_2(planet_dict['pPer'], planet_dict['pPerUnc'])),
            'uncertainty': str(round_to_2(planet_dict['pPerUnc'])) if planet_dict['pPerUnc'] else planet_dict['pPerUnc'],
            'units': "days"
        },
        'Rp/R*': {
            'value': str(round_to_2(planet_dict['rprs'], planet_dict['rprsUnc'])),
            'uncertainty': str(round_to_2(planet_dict['rprsUnc'])) if planet_dict['rprsUnc'] else planet_dict['rprsUnc'],
        },
        'a/R*': {
            'value': str(round_to_2(planet_dict['aRs'], planet_dict['aRsUnc'])),
            'uncertainty': str(round_to_2(planet_dict['aRsUnc'])) if planet_dict['aRsUnc'] else planet_dict['aRsUnc'],
        },
        'inc': {
            'value': str(round_to_2(planet_dict['inc'], planet_dict['incUnc'])),
            'uncertainty': str(round_to_2(planet_dict['incUnc'])) if planet_dict['incUnc'] else planet_dict['incUnc'],
            'units': "degrees"
        },
        'ecc': {
            'value': str(round_to_2(planet_dict['ecc'])),
            'uncertainty': None,
        },
        'u0': {
            'value': str(round_to_2(ld0[0], ld0[1])),
            'uncertainty': str(round_to_2(ld0[1]))
        },
        'u1': {
            'value': str(round_to_2(ld1[0], ld1[1])),
            'uncertainty': str(round_to_2(ld1[1]))
        },
        'u2': {
            'value': str(round_to_2(ld2[0], ld2[1])),
            'uncertainty': str(round_to_2(ld2[1]))
        },
        'u3': {
            'value': str(round_to_2(ld3[0], ld3[1])),
            'uncertainty': str(round_to_2(ld3[1]))
        }
    }

    filter_type = {
        'name': info_dict['filter'],
        'desc': info_dict['filter_desc'],
        'filter_width': {
            'left_side_wavelength': {
                'value': str(info_dict['wl_min']) if info_dict['wl_min'] else info_dict['wl_min'],
                'units': "nm"
            },
            'right_side_wavelength': {
                'value': str(info_dict['wl_max']) if info_dict['wl_max'] else info_dict['wl_max'],
                'units': "nm"
            }
        },
    }

    results = {
        'Tc': {
            'value': str(round_to_2(fit.parameters['tmid'], fit.errors['tmid'])),
            'uncertainty': str(round_to_2(fit.errors['tmid'])),
            'units': "BJD_TDB"
        },
        'Rp/R*': {
            'value': str(round_to_2(fit.parameters['rprs'], rprs_report_error)),
            'uncertainty': str(round_to_2(rprs_report_error))
        },
        'inc': {
            'value': str(round_to_2(fit.parameters['inc'], fit.errors['inc'])),
            'uncertainty': str(round_to_2(fit.errors['inc'])),
        },
        'Am2': {
            'value': aavso_airmass_terms[1][1],
            'uncertainty': aavso_airmass_terms[1][2]
        },
        'Duration': {
            'value': str(round_to_2(mean(durs))),
            'uncertainty': str(round_to_2(std(durs))),
            'units': "days"
        }
    }

    results[aavso_airmass_terms[0][0]] = {
        'value': aavso_airmass_terms[0][1],
        'uncertainty': aavso_airmass_terms[0][2]
    }
    optional_results = {
        'a/R*': aavso_result_entry(
            fit.parameters.get('ars'),
            fit.errors.get('ars'),
        ),
    }
    impact_parameter, impact_error = fit_impact_parameter_value_error(fit)
    optional_results['Impact Parameter (b)'] = aavso_result_entry(impact_parameter, impact_error)

    limb_darkening = (ld0, ld1, ld2, ld3)
    prior_parameters = planet_dict_transit_parameters(
        planet_dict,
        limb_darkening=limb_darkening,
        fallback=getattr(fit, 'prior', None),
    )
    prior_errors = planet_dict_transit_errors(planet_dict, limb_darkening=limb_darkening)
    depth_summary = fit_transit_depth_summary(
        fit,
        prior_parameters=prior_parameters,
        prior_errors=prior_errors,
    )
    for label, value_key, error_key in (
        (AREA_DEPTH_LABEL, 'area_depth', 'area_depth_error'),
        (OBSERVABLE_DEPTH_LABEL, 'observable_depth', 'observable_depth_error'),
        (PRIOR_OBSERVABLE_DEPTH_LABEL, 'prior_observable_depth', 'prior_observable_depth_error'),
        (OBSERVABLE_DEPTH_DELTA_LABEL, 'observable_depth_prior_delta', 'observable_depth_prior_delta_error'),
    ):
        optional_results[label] = aavso_result_entry(
            depth_summary.get(value_key),
            depth_summary.get(error_key),
            units="percent",
        )

    scatter = residual_scatter_fraction(fit)
    optional_results['Residual scatter around full model fit'] = aavso_result_entry(
        100.0 * scatter if np.isfinite(scatter) else np.nan,
        units="percent",
    )
    if 'a0' in fit.parameters:
        optional_results['a0'] = aavso_result_entry(fit.parameters.get('a0'), fit.errors.get('a0'))
    elif 'a1' in fit.parameters:
        optional_results['a1'] = aavso_result_entry(fit.parameters.get('a1'), fit.errors.get('a1'))

    results.update({
        key: value for key, value in optional_results.items()
        if value is not None
    })

    return priors, filter_type, results


def format_aavso_header_value(value):
    if value is None:
        return ""
    if isinstance(value, str):
        stripped = value.strip()
        return "" if stripped.lower() in ('', 'n/a', 'na', 'null', 'none') else stripped
    return str(value)


def save_comp_star_calibration_summary(save_dir, target_name, date, method_label, field_score,
                                       comp_summaries, best_comp_index):
    temp_dir = Path(save_dir) / "working_artifacts"
    temp_dir.mkdir(parents=True, exist_ok=True)
    summary_file = temp_dir / safe_output_filename(
        "CompStarCalibrationSummary",
        target_name,
        filename_date_token(date),
        extension="csv",
    )

    with summary_file.open('w') as handle:
        handle.write(f"# Comparison-star calibration summary for {target_name}\n")
        handle.write(f"# Method,{method_label}\n")
        if field_score is not None and field_score == field_score:
            handle.write(f"# Field suitability score,{field_score}\n")
        else:
            handle.write("# Field suitability score,\n")
        handle.write(f"# Selected comparison star,{'' if best_comp_index is None else best_comp_index + 1}\n")
        handle.write("comp_star,x_pixel,y_pixel,selected,suitability_score,ensemble_score,pairwise_median_score,"
                     "pairwise_max_score,self_score,valid_pair_count,coverage_count,coverage_peer_median,"
                     "coverage_min_required,coverage_rejected,suitability_outlier_rejected,"
                     "psf_quality_rejected_count,overexposure_rejected_count,ensemble_frame_rejected_count,"
                     "ensemble_frame_required_valid_pairs\n")

        for summary in comp_summaries:
            position = summary.get('position') or [None, None]
            values = [
                summary.get('label', ''),
                position[0],
                position[1],
                str(bool(summary.get('selected'))).lower(),
                summary.get('aggregate_score'),
                summary.get('ensemble_score'),
                summary.get('pairwise_median_score'),
                summary.get('pairwise_max_score'),
                summary.get('self_score'),
                summary.get('valid_pair_count'),
                summary.get('coverage_count'),
                summary.get('coverage_reference_count'),
                summary.get('coverage_min_required_count'),
                summary.get('coverage_rejected'),
                summary.get('suitability_outlier_rejected'),
                summary.get('psf_quality_rejected_count', 0),
                summary.get('overexposure_rejected_count', 0),
                summary.get('ensemble_frame_rejected_count', 0),
                summary.get('ensemble_frame_required_valid_pairs', 0),
            ]
            handle.write(",".join("" if value is None else str(value) for value in values) + "\n")

    return summary_file
