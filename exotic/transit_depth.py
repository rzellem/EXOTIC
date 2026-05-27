import math

import numpy as np


AREA_DEPTH_LABEL = "Radius-ratio area depth (Rp/R*)^2"
OBSERVABLE_DEPTH_LABEL = "Observable model transit depth"
PRIOR_OBSERVABLE_DEPTH_LABEL = "Prior observable model transit depth"
OBSERVABLE_DEPTH_DELTA_LABEL = "Observable model depth change from prior"

_DEPTH_ERROR_KEYS = ("rprs", "ars", "inc", "ecc", "omega", "u0", "u1", "u2", "u3")
_REQUIRED_TRANSIT_KEYS = ("rprs", "per", "ars", "inc", "ecc", "omega", "tmid", "u0", "u1", "u2", "u3")


def finite_float(value, default=np.nan):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return default
    return value if np.isfinite(value) else default


def radius_ratio_area_depth_percent(rprs, rprs_error=None):
    rprs = finite_float(rprs)
    if not np.isfinite(rprs) or rprs < 0:
        return np.nan, np.nan

    depth = 100.0 * rprs ** 2
    rprs_error = finite_float(rprs_error)
    if np.isfinite(rprs_error) and rprs_error >= 0:
        return float(depth), float(200.0 * abs(rprs) * rprs_error)
    return float(depth), np.nan


def planet_dict_transit_parameters(planet_dict, limb_darkening=None, fallback=None):
    values = dict(fallback or {})
    planet_dict = planet_dict or {}

    aliases = {
        "rprs": ("rprs", "pl_ratror"),
        "per": ("pPer", "pl_orbper", "per", "period"),
        "ars": ("aRs", "pl_ratdor", "ars"),
        "inc": ("inc", "pl_orbincl"),
        "ecc": ("ecc", "pl_orbeccen"),
        "omega": ("omega", "pl_orblper"),
        "tmid": ("midT", "pl_tranmid", "tmid"),
    }
    for target, names in aliases.items():
        for name in names:
            if name not in planet_dict:
                continue
            value = finite_float(planet_dict.get(name))
            if np.isfinite(value):
                values[target] = value
                break

    if limb_darkening is not None:
        for index, item in enumerate(limb_darkening):
            if index > 3:
                break
            if isinstance(item, (list, tuple, np.ndarray)):
                value = item[0] if len(item) else np.nan
            else:
                value = item
            value = finite_float(value)
            if np.isfinite(value):
                values[f"u{index}"] = value

    for key in ("u0", "u1", "u2", "u3"):
        values.setdefault(key, 0.0)
    values.setdefault("ecc", 0.0)
    values.setdefault("omega", 0.0)
    values.setdefault("tmid", 0.0)
    return values


def planet_dict_transit_errors(planet_dict, limb_darkening=None, fallback=None):
    errors = dict(fallback or {})
    planet_dict = planet_dict or {}

    aliases = {
        "rprs": ("rprsUnc", "pl_ratrorerr1"),
        "per": ("pPerUnc", "pl_orbpererr1"),
        "ars": ("aRsUnc", "pl_ratdorerr1"),
        "inc": ("incUnc", "pl_orbinclerr1"),
        "tmid": ("midTUnc", "pl_tranmiderr1"),
    }
    for target, names in aliases.items():
        for name in names:
            if name not in planet_dict:
                continue
            value = abs(finite_float(planet_dict.get(name)))
            if np.isfinite(value):
                errors[target] = value
                break

    if limb_darkening is not None:
        for index, item in enumerate(limb_darkening):
            if index > 3:
                break
            if not isinstance(item, (list, tuple, np.ndarray)) or len(item) < 2:
                continue
            value = abs(finite_float(item[1]))
            if np.isfinite(value):
                errors[f"u{index}"] = value
    return errors


def complete_transit_parameters(parameters):
    values = planet_dict_transit_parameters(parameters)
    if "per" not in values and "period" in values:
        values["per"] = values["period"]
    return values


def transit_duration_days(parameters):
    values = complete_transit_parameters(parameters)
    period = finite_float(values.get("per"))
    rprs = finite_float(values.get("rprs"))
    ars = finite_float(values.get("ars"))
    inc = finite_float(values.get("inc"))
    ecc = finite_float(values.get("ecc"), 0.0)
    omega = math.radians(finite_float(values.get("omega"), 0.0))

    if (
        not np.isfinite(period) or period <= 0
        or not np.isfinite(rprs) or rprs < 0
        or not np.isfinite(ars) or ars <= 0
        or not np.isfinite(inc)
        or not np.isfinite(ecc) or ecc < 0 or ecc >= 1
    ):
        return np.nan

    sin_inc = math.sin(math.radians(inc))
    if not np.isfinite(sin_inc) or sin_inc <= 0:
        return np.nan

    denominator = 1.0 + ecc * math.sin(omega)
    if not np.isfinite(denominator) or math.isclose(denominator, 0.0):
        return np.nan

    impact_scale = ars * (1.0 - ecc ** 2) / denominator
    impact_parameter = impact_scale * math.cos(math.radians(inc))
    chord_sq = (1.0 + rprs) ** 2 - impact_parameter ** 2
    if not np.isfinite(chord_sq) or chord_sq <= 0 or impact_scale <= 0:
        return np.nan

    argument = math.sqrt(chord_sq) / (impact_scale * sin_inc)
    argument = float(np.clip(argument, -1.0, 1.0))
    duration = (period / math.pi) * math.asin(argument)
    return float(duration) if np.isfinite(duration) and duration > 0 else np.nan


def impact_parameter(parameters):
    values = complete_transit_parameters(parameters)
    ars = finite_float(values.get("ars"))
    inc = finite_float(values.get("inc"))
    ecc = finite_float(values.get("ecc"), 0.0)
    omega = math.radians(finite_float(values.get("omega"), 0.0))
    if not np.isfinite(ars) or not np.isfinite(inc) or not np.isfinite(ecc):
        return np.nan
    denominator = 1.0 + ecc * math.sin(omega)
    if not np.isfinite(denominator) or math.isclose(denominator, 0.0):
        return np.nan
    return float(ars * (1.0 - ecc ** 2) * math.cos(math.radians(inc)) / denominator)


def geometric_observable_depth_fraction(parameters):
    values = complete_transit_parameters(parameters)
    rprs = finite_float(values.get("rprs"))
    b = abs(impact_parameter(values))
    if not np.isfinite(rprs) or rprs < 0 or not np.isfinite(b):
        return np.nan
    if b >= 1.0 + rprs:
        return 0.0
    if b <= abs(1.0 - rprs):
        return float(min(rprs ** 2, 1.0))
    if b <= 0:
        return float(min(rprs ** 2, 1.0))

    star_radius = 1.0
    planet_radius = rprs
    cos_star = np.clip(
        (b ** 2 + star_radius ** 2 - planet_radius ** 2) / (2.0 * b * star_radius),
        -1.0,
        1.0,
    )
    cos_planet = np.clip(
        (b ** 2 + planet_radius ** 2 - star_radius ** 2) / (2.0 * b * planet_radius),
        -1.0,
        1.0,
    )
    overlap = (
        star_radius ** 2 * math.acos(cos_star)
        + planet_radius ** 2 * math.acos(cos_planet)
        - 0.5 * math.sqrt(
            max(
                0.0,
                (-b + star_radius + planet_radius)
                * (b + star_radius - planet_radius)
                * (b - star_radius + planet_radius)
                * (b + star_radius + planet_radius),
            )
        )
    )
    return float(np.clip(overlap / math.pi, 0.0, 1.0))


def transit_depth_evaluation_times(parameters, sample_count=2000):
    values = complete_transit_parameters(parameters)
    tmid = finite_float(values.get("tmid"), 0.0)
    period = finite_float(values.get("per"))
    duration = transit_duration_days(values)
    if np.isfinite(duration) and duration > 0:
        half_window = duration
    elif np.isfinite(period) and period > 0:
        half_window = min(0.2, 0.1 * period)
    else:
        half_window = 0.2
    half_window = max(float(half_window), 1.0e-4)
    return np.linspace(tmid - half_window, tmid + half_window, int(sample_count))


def _load_transit_model():
    try:
        from .api.elca import transit
    except Exception:
        try:
            from api.elca import transit
        except Exception:
            return None
    return transit


def _depth_fraction_from_model_flux(model_flux):
    try:
        flux = np.asarray(model_flux, dtype=float).reshape(-1)
    except (TypeError, ValueError):
        return np.nan
    finite = flux[np.isfinite(flux)]
    if finite.size == 0:
        return np.nan
    return float(max(0.0, 1.0 - np.nanmin(finite)))


def observable_depth_fraction(parameters, model_flux=None, times=None):
    values = complete_transit_parameters(parameters)
    if not all(np.isfinite(finite_float(values.get(key))) for key in _REQUIRED_TRANSIT_KEYS):
        fallback_depth = _depth_fraction_from_model_flux(model_flux)
        return fallback_depth if np.isfinite(fallback_depth) else geometric_observable_depth_fraction(values)

    transit_model = _load_transit_model()
    if transit_model is not None:
        try:
            if times is None:
                times = transit_depth_evaluation_times(values)
            flux = transit_model(np.asarray(times, dtype=float), values)
            depth = _depth_fraction_from_model_flux(flux)
            if np.isfinite(depth):
                return depth
        except Exception:
            pass

    fallback_depth = _depth_fraction_from_model_flux(model_flux)
    if np.isfinite(fallback_depth):
        return fallback_depth
    return geometric_observable_depth_fraction(values)


def _perturbed_value(key, value):
    value = finite_float(value)
    if not np.isfinite(value):
        return np.nan
    if key in ("rprs", "ars", "per"):
        return max(value, np.finfo(float).eps)
    if key == "ecc":
        return float(np.clip(value, 0.0, 0.999999))
    if key == "inc":
        return float(np.clip(value, 0.0, 180.0))
    return value


def observable_depth_uncertainty_fraction(parameters, errors):
    values = complete_transit_parameters(parameters)
    errors = errors or {}
    contributions = []
    for key in _DEPTH_ERROR_KEYS:
        center = finite_float(values.get(key))
        error = abs(finite_float(errors.get(key)))
        if not np.isfinite(center) or not np.isfinite(error) or error <= 0:
            continue

        lower_values = dict(values)
        upper_values = dict(values)
        lower_values[key] = _perturbed_value(key, center - error)
        upper_values[key] = _perturbed_value(key, center + error)
        lower_depth = observable_depth_fraction(lower_values)
        upper_depth = observable_depth_fraction(upper_values)
        if np.isfinite(lower_depth) and np.isfinite(upper_depth):
            contributions.append(0.5 * abs(upper_depth - lower_depth))

    if not contributions:
        return np.nan
    return float(np.sqrt(np.sum(np.square(contributions))))


def observable_depth_percent(parameters, errors=None, model_flux=None, times=None):
    depth_fraction = observable_depth_fraction(parameters, model_flux=model_flux, times=times)
    if not np.isfinite(depth_fraction):
        return np.nan, np.nan
    error_fraction = observable_depth_uncertainty_fraction(parameters, errors or {})
    return (
        float(100.0 * depth_fraction),
        float(100.0 * error_fraction) if np.isfinite(error_fraction) else np.nan,
    )


def fit_transit_depth_summary(fit, prior_parameters=None, prior_errors=None):
    parameters = dict(getattr(fit, "parameters", {}) or {})
    errors = dict(getattr(fit, "errors", {}) or {})
    model_flux = getattr(fit, "transit_upsample", None)
    times = getattr(fit, "time_upsample", None)

    area_depth, area_error = radius_ratio_area_depth_percent(
        parameters.get("rprs"),
        errors.get("rprs"),
    )
    observable_depth, observable_error = observable_depth_percent(
        parameters,
        errors,
        model_flux=model_flux,
        times=times,
    )

    if prior_parameters is None:
        prior_parameters = getattr(fit, "prior", None)
    prior_depth = np.nan
    prior_error = np.nan
    if prior_parameters:
        prior_depth, prior_error = observable_depth_percent(prior_parameters, prior_errors or {})

    delta = np.nan
    delta_error = np.nan
    if np.isfinite(observable_depth) and np.isfinite(prior_depth):
        delta = float(observable_depth - prior_depth)
        if np.isfinite(observable_error) and np.isfinite(prior_error):
            delta_error = float(np.hypot(observable_error, prior_error))

    return {
        "area_depth": area_depth,
        "area_depth_error": area_error,
        "observable_depth": observable_depth,
        "observable_depth_error": observable_error,
        "prior_observable_depth": prior_depth,
        "prior_observable_depth_error": prior_error,
        "observable_depth_prior_delta": delta,
        "observable_depth_prior_delta_error": delta_error,
    }
