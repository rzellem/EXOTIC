"""Time-reference conversion and precision-led FITS exposure timing.

Precision is estimated from the stored digits, not a claim about shutter or
clock calibration. Ambiguous legacy standards remain explicit diagnostics.
"""
from datetime import timezone
from decimal import Decimal, InvalidOperation
import re

import dateutil.parser
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u


def time_standard(value):
    if value is None or isinstance(value, (float, np.floating, int)) and (not np.isfinite(value) or value == 0):
        return ''
    text = str(value or '').upper().strip().replace('(', '_').replace(')', '')
    text = re.sub(r'[\s/-]+', '_', text).strip('_')
    if text in ('UNKNOWN', 'NONE', 'NAN', 'NULL'):
        return ''
    return text.replace('BJDTDB', 'BJD_TDB').replace('BJDUTC', 'BJD_UTC')


def reference_to_bjd_tdb(value, standard, ra=None, dec=None):
    """Convert a declared reference once; BJD clocks need no light-time again.

    BJD_UTC clock conversion uses Astropy's clock correction, accurate at the
    tens-of-ms level for conventional published BJD_UTC (Eastman et al. 2010).
    HJD/JD references use the geocentre: a published event has no site attached.
    """
    standard = time_standard(standard)
    frame, _, scale = standard.partition('_')
    if frame not in ('BJD', 'HJD', 'JD', 'MJD') or scale.lower() not in Time.SCALES:
        return None
    value = float(value) + (2400000.5 if frame == 'MJD' else 0)
    instant = Time(value, format='jd', scale=scale.lower())
    if frame == 'BJD':
        return float(instant.tdb.jd)
    if ra is None or dec is None:
        return None
    units = (u.hourangle if isinstance(ra, str) and ':' in ra else u.deg, u.deg)
    target = SkyCoord(ra, dec, unit=units)
    location = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
    site = Time(value, format='jd', scale=scale.lower(), location=location)
    if frame == 'HJD':
        # Fixed-point inversion. Each update contracts by the Earth's v/c
        # (~1e-4); five iterations resolve a <510-s correction below float JD
        # precision, without requiring an unknown original observatory.
        for _ in range(5):
            light = site.light_travel_time(target, kind='heliocentric').to_value(u.day)
            site = Time(value-light, format='jd', scale=scale.lower(), location=location)
    return float((site.tdb + site.light_travel_time(target)).jd)


def normalise_ephemeris(parameters, archive_parameters=None, lookup=None):
    """Retain the source epoch, standard and conversion separately from midT.

    Old EXOTIC's BJD-UTC key labelled archive TDB values incorrectly. It is a
    hint, not sufficient evidence to shift a legacy value by a minute.
    """
    result = dict(parameters)
    if (result.get('ephemeris_timing') or {}).get('normalised'):
        return result
    original = float(result['midT'])
    standard = time_standard(result.get('midTStandard'))
    source = result.get('midTSource')
    origin = result.get('midTStandardOrigin', 'explicit_metadata')
    lookup_error = None
    if not standard:
        matched = archive_parameters if isinstance(archive_parameters, dict) else {}
        if matched.get('midT') != original or not matched.get('midTStandard'):
            matched = {}
            if callable(lookup):
                try:
                    matched = lookup(result.get('pName'), original, result.get('pPer')) or {}
                except Exception as exc:
                    lookup_error = f'{type(exc).__name__}: {exc}'
        if matched.get('midT') == original:
            standard = time_standard(matched.get('midTStandard'))
            source = matched.get('midTSource')
            if standard:
                origin = 'matching_archive_epoch'
    converted = reference_to_bjd_tdb(original, standard, result.get('ra'), result.get('dec'))
    verified = converted is not None
    result['midT'] = converted if verified else original
    result['midTStandard'] = 'BJD_TDB' if verified else (standard or 'UNKNOWN')
    result['midTSource'] = source
    result['ephemeris_timing'] = {
        'normalised': True, 'source_epoch': original,
        'source_standard': standard or 'UNKNOWN', 'source_reference': source,
        'standard_origin': origin if standard else 'unverified_legacy_input',
        'legacy_label_hint': result.get('midTLegacyLabel'),
        'model_epoch': result['midT'], 'model_standard': 'BJD_TDB',
        'status': 'verified' if verified else 'unverified_assumed_model_standard',
        'conversion_applied': bool(verified and standard != 'BJD_TDB'),
        'conversion_seconds': (result['midT']-original)*86400 if verified else None,
        'lookup_error': lookup_error,
        'note': ('No light-time correction is applied again to an existing BJD. '
                 'Quoted epoch/period errors are retained; clock conversion is deterministic. '
                 'A period originally fitted in UTC across leap seconds may require its original timings.'
                 if verified else
                 'Time standard could not be verified. Numeric epoch retained under the existing '
                 'BJD_TDB model convention; catalogue timing offsets remain conditional on that assumption.'),
    }
    return result


def input_ephemeris_metadata(values):
    standard = values.get('Published Mid-Transit Time Standard', values.get('midTStandard'))
    legacy = None
    # New explicit keys are trustworthy; the old BJD-UTC key is ambiguous.
    for key in values:
        match = re.fullmatch(r'Published Mid-Transit Time \(([^)]+)\)', key)
        if not match or values[key] is None:
            continue
        declared = time_standard(match[1])
        if key == 'Published Mid-Transit Time (BJD-UTC)':
            legacy = key
        elif not standard:
            standard = declared
    return {'midTStandard': standard, 'midTSource': values.get('Published Mid-Transit Time Reference'),
            'midTLegacyLabel': legacy}


def _stored_number(header, key):
    value = header[key]
    if isinstance(value, tuple):
        value = value[0]
    text = str(value).strip()
    # FITS card spelling preserves trailing zeros (float repr does not).
    try:
        card_value = header.cards[key].image.split('=', 1)[1].split('/', 1)[0].strip()
        if not card_value.startswith("'"):
            text = card_value
    except (AttributeError, KeyError, IndexError):
        pass
    number = Decimal(text.replace('D', 'E'))
    if not number.is_finite():
        raise ValueError('non-finite timestamp')
    return float(number), float(Decimal(10) ** number.as_tuple().exponent)


def _date_candidate(header, key, scale):
    text = str(header[key]).strip()
    if key == 'DATE-OBS' and 'T' not in text and 'TIME-OBS' in header:
        text += 'T' + str(header['TIME-OBS']).strip()
    if key == 'UT-OBS' and 'T' not in text and 'DATE-OBS' in header:
        text = str(header['DATE-OBS']).split('T')[0] + 'T' + text
    stamp = dateutil.parser.parse(text)
    zone = re.search(r'(Z|[+-]\d\d:?\d\d)$', text) if 'T' in text else None
    if zone:
        offset = stamp.utcoffset().total_seconds() if stamp.utcoffset() is not None else 0.
        instant = Time(text[:zone.start()], scale=scale)-offset*u.s
    elif stamp.tzinfo is not None:
        instant = Time(stamp.astimezone(timezone.utc), scale=scale)
    # Time parses the text directly to retain sub-microsecond fractions.
    elif stamp.tzinfo is None and 'T' in text:
        try:
            instant = Time(text, scale=scale)
        except ValueError:
            instant = Time(stamp, scale=scale)
    else:
        instant = Time(stamp, scale=scale)
    match = re.search(r'\d\d:\d\d:\d\d(?:\.(\d+))?', text)
    precision = 10.0**(-len(match[1])) if match and match[1] else (1.0 if match else 60.0 if ':' in text else 86400.0)
    return instant, precision


def select_exposure_timestamp(header, exposure_seconds):
    """Rank all usable midpoint representations by stored precision.

    Declared reference/scale is preferred over an ambiguous generic BJD/TDB.
    Ties favour an explicit midpoint, then start/end, then start+EXPTIME.
    No disagreement cutoff removes otherwise usable timestamps.
    """
    declared_scale = str(header.get('TIMESYS', 'UTC')).strip().lower()
    raw_scale = declared_scale if declared_scale in Time.SCALES else 'utc'
    candidates, starts, ends, issues = [], [], [], []
    try:
        pixel = float(header.get('TIMEPIXR', 0.0))
        if not np.isfinite(pixel):
            raise ValueError('non-finite TIMEPIXR')
    except (TypeError, ValueError) as exc:
        issues.append({'source': 'TIMEPIXR', 'status': 'unusable', 'reason': str(exc),
                       'assumption': 'Existing start-of-exposure convention retained.'})
        pixel = 0.0
    try:
        _, exposure_precision = _stored_number(header, 'EXPTIME')
    except (KeyError, ValueError, InvalidOperation):
        exposure_precision = 0.0

    def add(key, frame, scale, role, date=False, ambiguous=False):
        if key not in header:
            return
        try:
            if frame == 'BJD' and 'MID' not in key:
                try:
                    comment = str(header.comments[key]).lower()
                except (AttributeError, KeyError):
                    comment = ''
                if 'start' in comment or 'begin' in comment:
                    role = 'start'
                elif 'end of' in comment:
                    role = 'end'
            if date:
                instant, precision = _date_candidate(header, key, scale)
            else:
                value, step = _stored_number(header, key)
                if frame == 'MJD':
                    value += 2400000.5
                instant = Time(value, format='jd', scale=scale)
                precision = max(step*86400, abs(np.spacing(value))*86400)
            barycentric_header = str(header.get('TIMEREF', '')).upper() in ('SOLARSYSTEM', 'BARYCENTER', 'BARYCENTRIC')
            row = {'source': key, 'frame': 'BJD' if frame == 'BJD' or barycentric_header else 'JD',
                   'scale': scale.upper(), 'precision_seconds': precision,
                   'role': role, 'ambiguous_standard': ambiguous,
                   'value': str(header[key]), '_instant': instant}
            if role == 'start':
                starts.append(row)
            elif role == 'end':
                ends.append(row)
            else:
                candidates.append(row)
        except (ValueError, TypeError, InvalidOperation, OverflowError) as exc:
            issues.append({'source': key, 'status': 'unusable', 'reason': str(exc)})

    for scale in ('TDB', 'UTC', 'TT', 'TAI', 'TCB'):
        for key in (f'BJD_{scale}', f'BJD-{scale}'):
            add(key, 'BJD', scale.lower(), 'mid')
    for key in ('BJD_TBD', 'TDB-MID'):
        add(key, 'BJD', 'tdb', 'mid')
    for key in ('BJD-MID', 'BJD', 'TDB'):
        scale = raw_scale if 'TIMESYS' in header else 'tdb'
        add(key, 'BJD', scale, 'mid', ambiguous='TIMESYS' not in header)
    for key in ('JD-MID', 'MJD-MID'):
        add(key, 'MJD' if key.startswith('MJD') else 'JD', raw_scale, 'mid')
    for key in ('DATE-AVG', 'DATE-MID'):
        add(key, 'JD', raw_scale, 'mid', date=True)
    for key in ('JD-START', 'JD', 'JULIAN', 'MJD-OBS', 'MJD'):
        add(key, 'MJD' if key.startswith('MJD') else 'JD', raw_scale, 'start')
    for key in ('DATE-UTC', 'DATE-BEG', 'DATE-OBS', 'UT-OBS'):
        add(key, 'JD', 'utc' if key in ('DATE-UTC', 'UT-OBS') else raw_scale, 'start', date=True)
    for key in ('DATE-END', 'END-OBS'):
        add(key, 'JD', raw_scale, 'end', date=True)
    for start in starts:
        for end in ends:
            if end['_instant'] < start['_instant'] or pixel != 0 or end['frame'] != start['frame']:
                continue
            row = dict(start)
            row.update(source=f"{start['source']}/{end['source']}", role='start_end_mid',
                       precision_seconds=(start['precision_seconds']+end['precision_seconds'])/2,
                       _instant=start['_instant']+(end['_instant']-start['_instant'])/2)
            candidates.append(row)
        row = dict(start)
        row.update(source=f"{start['source']}+EXPTIME*(0.5-TIMEPIXR)", role='start_exposure_mid',
                   precision_seconds=start['precision_seconds']+abs(.5-pixel)*exposure_precision,
                   _instant=start['_instant']+(.5-pixel)*float(exposure_seconds)*u.s)
        candidates.append(row)
    for end in ends:
        row = dict(end)
        row.update(source=f"{end['source']}-EXPTIME/2", role='end_exposure_mid',
                   precision_seconds=end['precision_seconds']+.5*exposure_precision,
                   _instant=end['_instant']-.5*float(exposure_seconds)*u.s)
        candidates.append(row)
    if not candidates:
        return {'status': 'unavailable', 'selected': None, 'candidates': [], 'issues': issues}
    order = {'mid': 0, 'start_end_mid': 1, 'start_exposure_mid': 2, 'end_exposure_mid': 2}
    candidates.sort(key=lambda row: (row['ambiguous_standard'], row['precision_seconds'], order[row['role']]))
    chosen = candidates[0]
    for row in candidates:
        row['jd_in_source_scale'] = float(row['_instant'].jd)
        row['jd_utc'] = float(row['_instant'].utc.jd) if row['frame'] == 'JD' else None
    return {'status': 'selected', 'selected': chosen, 'candidates': candidates, 'issues': issues,
            'precision_basis': 'Stored decimal resolution; not a calibrated clock/shutter accuracy or a Tmid error bar.',
            'raw_timesys': declared_scale.upper(),
            'raw_timesys_assumed': 'TIMESYS' not in header or declared_scale not in Time.SCALES,
            'declared_clock_accuracy': {key:str(header[key]) for key in ('TIMSYER','TIMRDER','TIMEUNIT') if key in header}}


def public_timestamp_selection(selection):
    return {key: ([{k:v for k,v in row.items() if k != '_instant'} for row in value]
                  if key == 'candidates' else
                  {k:v for k,v in value.items() if k != '_instant'} if key == 'selected' and value else value)
            for key,value in selection.items()}


def calculate_timestamp_tmid_difference(fit, selections):
    """Paired deterministic fits to the same retained observations.

    Only timestamps change. Reuse the actual exposure-integrated transit,
    likelihood baseline and duration prior. These are joint point estimates,
    not new posterior medians or histogram/Gaussian centres.
    """
    from scipy.optimize import minimize
    context = getattr(fit, 'orbital_likelihood_context', None)
    unavailable = {'status': 'unavailable', 'production_fit_changed': False}
    if not context or not hasattr(fit, '_transit_model'):
        return dict(unavailable, reason='The production likelihood context is unavailable.')
    selected_times = np.asarray(context['time'])
    # Alternatives are converted through the same barycentric calculation in
    # image_timestamp_solution. Do not approximate their correction by merely
    # adding a UTC difference (the barycentric term varies with arrival time).
    converted_pairs = {record['bjd_tdb']:record['mjd_date_comparison'] for record in selections
                       if record.get('mjd_date_comparison')}
    if not all(float(t) in converted_pairs for t in selected_times):
        return dict(unavailable, reason='Converted MJD-OBS and DATE-OBS timestamps are needed for each retained observation.',
                    matched_point_count=sum(float(t) in converted_pairs for t in selected_times),
                    retained_point_count=len(selected_times))
    times = {kind:np.array([converted_pairs[float(t)][kind+'_bjd_tdb'] for t in selected_times])
             for kind in ('date_obs','mjd_obs')}
    keys = [key for key in fit.sampled_keys if key not in ('ecc', 'omega')]
    if 'tmid' not in keys:
        return dict(unavailable, reason='Tmid was held fixed in the production fit.')
    base = dict(context['base_physical'])
    base.update(fit.parameters)
    bounds = getattr(fit, 'sample_bounds', fit.bounds)
    lower = np.array([0. if key == 'b' else bounds[key][0] for key in keys])
    width = np.array([1. if key == 'b' else bounds[key][1]-bounds[key][0] for key in keys])

    def physical(point):
        params = dict(base)
        params.update({key:float(low+span*x) for key,low,span,x in zip(keys,lower,width,point) if key != 'b'})
        if 'b' in keys:
            f = (1-params['ecc']**2)/(1+params['ecc']*np.sin(np.radians(params['omega'])))
            upper = min(1+params['rprs'], params['ars']*f)
            # Transit shape is approximately linear in b**2 near a central
            # transit. This keeps the diagnostic optimizer conditioned at b=0.
            params['b'] = float(np.sqrt(point[keys.index('b')])*upper)
            params['inc'] = float(np.degrees(np.arccos(params['b']/(params['ars']*f))))
        return params

    x0 = np.array([(base[key]-low)/span if key != 'b' else 0. for key,low,span in zip(keys,lower,width)])
    if 'b' in keys:
        f = (1-base['ecc']**2)/(1+base['ecc']*np.sin(np.radians(base['omega'])))
        impact = base.get('b',base['ars']*f*np.cos(np.radians(base['inc'])))
        x0[keys.index('b')] = (impact/min(1+base['rprs'],base['ars']*f))**2
    initial_outside = bool(np.any((x0 < 0) | (x0 > 1)))
    x0 = np.clip(x0, 0, 1)

    def residual(point, timestamps):
        params = physical(point)
        model = np.asarray(fit._transit_model(timestamps, params),dtype=float)
        model *= np.exp(params.get('a2',0.)*context['centered_airmass'])
        baseline_key = context['free_flux_baseline_key']
        if baseline_key:
            model *= params[baseline_key]
        elif context['uses_fixed_flux_baseline']:
            model *= context['fixed_flux_baseline_value']
        elif context['has_free_flux_baseline']:
            if __package__:
                from .api.elca import get_flux_baseline
            else:  # The command-line launcher also supports flat imports.
                from api.elca import get_flux_baseline
            model *= get_flux_baseline(params)
        else:
            mask = context['baseline_static_mask'] & np.isfinite(model) & (model != 0)
            denom = np.sum(context['baseline_weights'][mask]*model[mask]**2)
            if not np.any(mask):
                baseline = 1.
            elif not np.isfinite(denom) or denom <= 0:
                ratios = context['data'][mask]/model[mask]
                ratios = ratios[np.isfinite(ratios)]
                baseline = float(np.median(ratios)) if ratios.size else 1.
            else:
                baseline = np.sum(context['baseline_weights'][mask]*context['data'][mask]*model[mask])/denom
            model *= baseline if np.isfinite(baseline) else 1.
        res = (context['data']-model)*context['inverse_dataerr']
        if context['duration_prior_valid']:
            if __package__:
                from .api.elca import transit_duration
            else:
                from api.elca import transit_duration
            duration = transit_duration(params)
            res = np.append(res, np.log(duration/context['expected_duration'])/context['sigma_log_duration'])
        return res

    fitted = {}
    for kind,timestamps in times.items():
        # Adaptive exposure integration can introduce tiny model discontinuities
        # at its existing tolerance. Use a derivative-free optimizer for both
        # series rather than finite differences of that production model.
        solution = minimize(lambda point: float(np.sum(residual(point, timestamps)**2)),
                            x0=x0.copy(), method='Powell', bounds=[(0.,1.)]*len(keys),
                            options={'ftol':1e-10, 'xtol':1e-8})
        params = physical(solution.x)
        fitted[kind] = {'tmid_bjd_tdb':params['tmid'], 'parameters':params,
                        'objective_chisquare_including_duration_prior':float(solution.fun),
                        'optimizer_success':bool(solution.success), 'optimizer_message':str(solution.message),
                        'optimizer_evaluations':int(solution.nfev)}
    return {'status':'calculated', 'method':'paired_deterministic_transit_refits_same_retained_photometry',
            'optimizer':'bounded_Powell_in_unit_cube_with_b_squared_geometry',
            'point_count':len(selected_times), 'date_obs':fitted['date_obs'], 'mjd_obs':fitted['mjd_obs'],
            'date_obs_minus_mjd_obs_seconds':(fitted['date_obs']['tmid_bjd_tdb']-fitted['mjd_obs']['tmid_bjd_tdb'])*86400,
            'free_parameters':keys, 'initial_point_projected_into_existing_bounds':initial_outside,
            'production_fit_changed':False,
            'note':'Calculated by fitting both timestamp series to identical retained fluxes, errors, airmass, '
                   'baseline masks, parameter bounds, exposure integration and duration prior. Eccentricity '
                   'and periastron stay fixed. This compares deterministic joint point estimates, not posterior '
                   'histogram peaks, Gaussian centres or medians; no theoretical timestamp average is substituted.'}


def timing_reporting_metadata(parameters, info, fit=None):
    comparison = getattr(fit, 'timestamp_tmid_comparison', None) if fit is not None else None
    if fit is not None and comparison is None and info.get('timestamp_selection_summary'):
        import json
        from pathlib import Path
        try:
            path = info['timestamp_selection_summary']['report_path']
            comparison = calculate_timestamp_tmid_difference(fit, json.loads(Path(path).read_text()))
        except Exception as exc:
            comparison = {'status':'failed', 'reason':f'{type(exc).__name__}: {exc}', 'production_fit_changed':False}
        fit.timestamp_tmid_comparison = comparison
    return {'model_time_standard': 'BJD_TDB',
            'ephemeris': parameters.get('ephemeris_timing'),
            'image_timestamps': info.get('timestamp_selection_summary'),
            'calculated_mjd_obs_vs_date_obs_tmid': comparison,
            'prereduced_input_time_standard': info.get('file_time')}
