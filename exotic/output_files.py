from json import dump, dumps
import shutil
from numpy import mean, std
from pathlib import Path
import numpy as np

try:
    from utils import filename_date_token, round_to_2, safe_output_filename
except ImportError:
    from .utils import filename_date_token, round_to_2, safe_output_filename
try:
    from version import __version__
except ImportError:
    from .version import __version__
try:
    from plate_status import PlateStatus
except ImportError:
    from .plate_status import PlateStatus


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


def format_aavso_json_header(name, payload):
    payload = prune_aavso_metadata(payload)
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
        'rprs_sigma', 'duration_ratio', 'eebls_depth_snr',
        'use_deviation_from_expected_transit_in_qc', 'deviation_sigma_threshold',
        'expected_tmid', 'expected_tmid_unc', 'expected_tmid_unc_minutes',
        'fitted_tmid', 'expected_rprs', 'expected_rprs_unc',
        'tmid_deviation_days', 'tmid_deviation_minutes',
        'tmid_deviation_threshold_minutes', 'tmid_deviation_sigma',
        'rprs_deviation_sigma', 'tmid_deviation_score', 'rprs_deviation_score',
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
            'deviation_from_expected_value': transit_qc.get('deviation_from_expected_value'),
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


def format_optional_metric(label, value, precision=2):
    value = finite_float(value)
    if not np.isfinite(value):
        return None
    return f"{label}={value:.{precision}f}"


def format_ktmf_candidate_decision(attempt):
    attempt = compact_comparison_attempt_decision(attempt)
    label = attempt.get('label') or (
        f"Comp {attempt['comparison_star']}" if attempt.get('comparison_star') is not None else "Comparison candidate"
    )
    selected_text = " [selected]" if attempt.get('selected') else ""
    parts = [
        f"{label}{selected_text}: KTMF={format_ktmf_metric(attempt.get('ktmf_metric'))}",
    ]
    for metric_text in (
        format_optional_metric("Delta BIC", attempt.get('transit_delta_bic')),
        format_optional_metric("EEBLS SNR", attempt.get('eebls_snr')),
    ):
        if metric_text:
            parts.append(metric_text)
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
            target_status = str(transit_qc.get('status', 'unknown')).upper()
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


def fit_impact_parameter_value_error(fit):
    parameters = getattr(fit, 'parameters', {}) or {}
    errors = getattr(fit, 'errors', {}) or {}
    sample_parameters = getattr(fit, 'sample_parameters', {}) or {}
    sample_errors = getattr(fit, 'sample_errors', {}) or {}

    if 'b' in sample_parameters:
        impact_parameter = finite_float(sample_parameters.get('b'))
        impact_error = finite_float(sample_errors.get('b'))
        if np.isfinite(impact_parameter):
            return impact_parameter, impact_error

    if 'b' in parameters:
        impact_parameter = finite_float(parameters.get('b'))
        impact_error = finite_float(errors.get('b'))
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

    ars_error = finite_float(errors.get('ars'))
    inc_error = finite_float(errors.get('inc'))
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
        params_file = self.dir / "temp" / safe_output_filename(
            "FinalLightCurve",
            self.p_dict['pName'],
            filename_date_token(self.i_dict['date']),
            extension="csv",
        )

        with params_file.open('w') as f:
            f.write(f"# FINAL TIMESERIES OF {self.p_dict['pName']}\n")
            f.write("# BJD_TDB,Orbital Phase,Flux,Uncertainty,Model,Airmass\n")

            for bjd, phase, flux, fluxerr, model, am in zip(self.fit.time, phase, self.fit.detrended,
                                                            self.fit.dataerr / self.fit.airmass_model,
                                                            self.fit.transit, self.fit.airmass_model):
                f.write(f"{bjd}, {phase}, {flux}, {fluxerr}, {model}, {am}\n")

    def final_planetary_params(self, phot_opt, vsp_params, comp_star=None, comp_coords=None, min_aper=None,
                               min_annul=None, adaptive_summary=None, photometry_info=None,
                               publish_to_root=False):
        params_file = self.dir / "temp" / safe_output_filename(
            "FinalParams",
            self.p_dict['pName'],
            filename_date_token(self.i_dict['date']),
            extension="json",
        )

        transit_qc = getattr(self.fit, 'transit_qc', None)
        fit_quality = build_fit_quality_metadata(self.fit)
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

        params_num = {
            "Mid-Transit Time (Tmid)": f"{round_to_2(self.fit.parameters['tmid'], self.fit.errors['tmid'])} +/- "
                                       f"{round_to_2(self.fit.errors['tmid'])} BJD_TDB",
            "Ratio of Planet to Stellar Radius (Rp/R*)": f"{round_to_2(self.fit.parameters['rprs'], self.fit.errors['rprs'])} +/- "
                                                         f"{round_to_2(self.fit.errors['rprs'])}",
            "Transit depth (Rp/Rs)^2": f"{round_to_2(100. * (self.fit.parameters['rprs'] ** 2.))} +/- "
                                       f"{round_to_2(100. * 2. * self.fit.parameters['rprs'] * self.fit.errors['rprs'])} [%]",
            "Orbital Inclination (inc)": f"{round_to_2(self.fit.parameters['inc'], self.fit.errors['inc'])} +/- "
                                                   f"{round_to_2(self.fit.errors['inc'])} ",
        }
        ars_text = format_parameter_with_error(
            self.fit.parameters.get('ars'),
            self.fit.errors.get('ars'),
        )
        if ars_text is not None:
            params_num["Ratio of Distance to Stellar Radius (a/Rs)"] = ars_text
        impact_parameter, impact_error = fit_impact_parameter_value_error(self.fit)
        impact_text = format_parameter_with_error(impact_parameter, impact_error)
        if impact_text is not None:
            params_num["Impact Parameter (b)"] = impact_text
        if getattr(self.fit, 'ns_type', None) is not None:
            params_num["Fit parameter point estimate"] = (
                "Best-fit likelihood point; uncertainties are posterior spread."
            )
        prefit_refinement_note = getattr(self.fit, 'prefit_refinement_note', None)
        if prefit_refinement_note:
            params_num["Prefit refinement note"] = str(prefit_refinement_note)
        oot_baseline_parameter_note = getattr(self.fit, 'oot_baseline_parameter_fit_note', None)
        if oot_baseline_parameter_note:
            params_num["Out-of-transit baseline parameter-fit note"] = str(oot_baseline_parameter_note)
        oot_baseline_note = getattr(self.fit, 'oot_baseline_detrending_note', None)
        if oot_baseline_note:
            params_num["Out-of-transit baseline detrending note"] = str(oot_baseline_note)
        sparse_posterior_note = getattr(self.fit, 'sparse_posterior_live_point_extension_note', None)
        if sparse_posterior_note:
            params_num["Sparse posterior live-point extension note"] = str(sparse_posterior_note)
        if np.isfinite(qc_residual_scatter):
            params_num["Residual scatter around full model fit"] = f"{qc_residual_scatter * 100.0:.4f} %"
        params_num.update(format_fit_quality_final_params(fit_quality))
        params_num.update(format_ktmf_decision_final_params(self.fit, photometry_info))
        if getattr(self.fit, 'airmass_fit_skipped', False):
            params_num["Airmass correction"] = getattr(
                self.fit,
                'airmass_correction_note',
                "Skipped; no airmass correction applied.",
            )
        else:
            if 'a0' in self.fit.parameters:
                params_num["Baseline flux (a0)"] = (
                    f"{round_to_2(self.fit.parameters['a0'], self.fit.errors['a0'])} +/- "
                    f"{round_to_2(self.fit.errors['a0'])}"
                )
            else:
                params_num["Flux normalization (a1)"] = (
                    f"{round_to_2(self.fit.parameters['a1'], self.fit.errors['a1'])} +/- "
                    f"{round_to_2(self.fit.errors['a1'])}"
                )
            params_num["Airmass coefficient 2 (a2)"] = (
                f"{round_to_2(self.fit.parameters['a2'], self.fit.errors['a2'])} +/- "
                f"{round_to_2(self.fit.errors['a2'])}"
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
            qc_tmid_deviation_sigma = transit_qc.get('tmid_deviation_sigma', np.nan)
            qc_tmid_deviation_minutes = transit_qc.get('tmid_deviation_minutes', np.nan)
            qc_tmid_threshold_minutes = transit_qc.get('tmid_deviation_threshold_minutes', np.nan)
            qc_expected_tmid_unc_minutes = transit_qc.get('expected_tmid_unc_minutes', np.nan)
            qc_rprs_deviation_sigma = transit_qc.get('rprs_deviation_sigma', np.nan)
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
            if np.isfinite(qc_tmid_deviation_sigma):
                params_num["Expected-value Tmid deviation"] = f"{qc_tmid_deviation_sigma:.2f} sigma"
            if np.isfinite(qc_tmid_deviation_minutes):
                params_num["Expected-value Tmid offset"] = f"{qc_tmid_deviation_minutes:.2f} minutes"
            if np.isfinite(qc_expected_tmid_unc_minutes):
                params_num["Expected-value Tmid uncertainty"] = f"{qc_expected_tmid_unc_minutes:.2f} minutes"
            if np.isfinite(qc_tmid_threshold_minutes):
                params_num["Expected-value Tmid QC window"] = f"{qc_tmid_threshold_minutes:.2f} minutes"
            if np.isfinite(qc_rprs_deviation_sigma):
                params_num["Expected-value Rp/R* deviation"] = f"{qc_rprs_deviation_sigma:.2f} sigma"
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

        if vsp_params:
            params_num["Variable Reference Star"] = f"AAVSO Label: {vsp_params[0]['cname']}, " + \
                                                    f"Position: {vsp_params[0]['pos']}"

        if phot_opt:
            phot_ext = {"Best Comparison Star": f"#{comp_star} - {comp_coords}" if min_aper >= 0 else str(comp_star)}
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
        ktmf_decision_metadata = build_ktmf_decision_metadata(self.fit, photometry_info)
        photometry_metadata = build_aavso_photometry_metadata(photometry_info)
        aperture_metadata = build_aavso_aperture_metadata(photometry_info)
        frame_filtering_metadata = build_aavso_frame_filtering_metadata(self.fit, frame_filtering_info)
        astrometry_metadata = build_aavso_astrometry_metadata(astrometry_info, comp_star)
        bad_pixel_metadata = build_aavso_bad_pixel_metadata(bad_pixel_info)
        obs_name = format_aavso_header_value(self.i_dict.get('obs_name'))
        obs_name_header = f"#OBSNAME={obs_name}\n" if obs_name else ""
        gaia_dist = format_aavso_header_value(self.p_dict.get('dist'))
        gaia_pmra = format_aavso_header_value(self.p_dict.get('pm_ra'))
        gaia_pmdec = format_aavso_header_value(self.p_dict.get('pm_dec'))
        gaia_dist_header = f"#GAIADIST={gaia_dist}\n" if gaia_dist else ""
        gaia_pmra_header = f"#GAIAPMRA={gaia_pmra}\n" if gaia_pmra else ""
        gaia_pmdec_header = f"#GAIAPMDEC={gaia_pmdec}\n" if gaia_pmdec else ""

        params_file = self.dir / safe_output_filename(
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
                    f",Rp/R*={round_to_2(self.fit.parameters['rprs'], self.fit.errors['rprs'])} +/- {round_to_2(self.fit.errors['rprs'])}"
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
    def plate_status(self, plate_status: PlateStatus):
        plate_status_file = self.dir / "temp" / safe_output_filename(
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

    def aavso(self):
        params_file = self.dir / safe_output_filename(
            "AID_AAVSO",
            self.p_dict['sName'],
            filename_date_token(self.i_dict['date']),
            extension="txt",
        )
        with params_file.open('w', encoding="utf-8") as f:
            f.write("#TYPE=EXTENDED\n"  # fixed
                    f"#OBSCODE={self.i_dict['aavso_num']}\n"  # UI
                    f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                    "#DELIM=,\n"  # fixed
                    "#DATE=JD\n"  # fixed
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

            f.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")
            for vsp_p in self.vsp_params:
                f.write(f"{self.auid},{round(vsp_p['time'], 5)},{round(vsp_p['mag'], 5)},{round(vsp_p['mag_err'], 5)},"
                        f"{self.i_dict['filter']},NO,STD,{vsp_p['cname']},{round(vsp_p['cmag'], 5)},na,na," 
                        f"{round(vsp_p['airmass'], 7)},na,{self.chart_id},na\n")


def aavso_dicts(planet_dict, fit, info_dict, durs, ld0, ld1, ld2, ld3):
    aavso_airmass_terms = aavso_airmass_results(fit)
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
            'value': str(round_to_2(fit.parameters['rprs'], fit.errors['rprs'])),
            'uncertainty': str(round_to_2(fit.errors['rprs']))
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

    rprs = finite_float(fit.parameters.get('rprs'))
    rprs_error = finite_float(fit.errors.get('rprs'))
    if np.isfinite(rprs):
        optional_results['Transit depth (Rp/R*)^2'] = aavso_result_entry(
            100.0 * (rprs ** 2.0),
            100.0 * 2.0 * rprs * rprs_error if np.isfinite(rprs_error) else np.nan,
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
    temp_dir = Path(save_dir) / "temp"
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
                     "coverage_min_required,coverage_rejected,suitability_outlier_rejected\n")

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
            ]
            handle.write(",".join("" if value is None else str(value) for value in values) + "\n")

    return summary_file
