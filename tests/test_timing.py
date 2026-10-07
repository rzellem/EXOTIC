import json
import numpy as np
import pytest
from astropy.io import fits
from astropy.time import Time
from astropy.utils import iers
from exotic.timing import (select_exposure_timestamp, public_timestamp_selection,
                           normalise_ephemeris, input_ephemeris_metadata,
                           reference_to_bjd_tdb)

iers.conf.auto_download = False


def test_precise_start_end_outweigh_rounded_mjd():
    hdr = fits.Header({'DATE-OBS': '2021-06-11T20:25:45.449-0700',
                       'DATE-END': '2021-06-11T20:26:47.403-0700',
                       'MJD-OBS': 59377.143, 'EXPTIME': 60.0})
    result = select_exposure_timestamp(hdr, 60)
    chosen = result['selected']
    assert chosen['source'] == 'DATE-OBS/DATE-END'
    assert chosen['precision_seconds'] == pytest.approx(.001)
    assert chosen['_instant'].utc.isot == '2021-06-12T03:26:16.426'
    coarse = next(row for row in result['candidates'] if row['source'] == 'MJD-OBS/DATE-END')
    assert (coarse['_instant']-chosen['_instant']).sec == pytest.approx(4.8755, abs=.00005)
    assert json.loads(json.dumps(public_timestamp_selection(result)))['selected']['source'] == chosen['source']


@pytest.mark.parametrize('key,value', [('JD-MID', 2459377.6), ('MJD-MID', 59377.1),
                                      ('DATE-AVG', '2021-06-12T03:00:00')])
def test_coarse_midpoint_does_not_override_precise_start(key, value):
    hdr = fits.Header({key: value, 'DATE-OBS': '2021-06-12T03:00:00.123456', 'EXPTIME': 30.0})
    chosen = select_exposure_timestamp(hdr, 30)['selected']
    assert chosen['source'].startswith('DATE-OBS+')


def test_precise_numeric_midpoint_can_win_over_date():
    hdr = fits.Header({'JD-MID': 2459377.600000123, 'DATE-OBS': '2021-06-12T03:00:00', 'EXPTIME': 30.0})
    assert select_exposure_timestamp(hdr, 30)['selected']['source'] == 'JD-MID'


def test_midpoint_is_not_offset_twice():
    hdr = fits.Header({'DATE-AVG': '2021-06-12T03:00:00.123456', 'EXPTIME': 30.0})
    result = select_exposure_timestamp(hdr, 30)['selected']
    assert result['_instant'].isot == '2021-06-12T03:00:00.123'
    assert result['source'] == 'DATE-AVG'


def test_end_only_timestamp_remains_usable():
    hdr=fits.Header({'DATE-END':'2021-06-12T03:00:30.123','EXPTIME':60.})
    row=select_exposure_timestamp(hdr,60)['selected']
    assert row['source']=='DATE-END-EXPTIME/2'
    assert row['_instant'].utc.isot=='2021-06-12T03:00:00.123'


def test_timepixr_midpoint_uses_no_half_exposure_offset():
    hdr = fits.Header({'DATE-OBS': '2021-06-12T03:00:00.123456', 'TIMEPIXR': .5, 'EXPTIME': 30.0})
    result = select_exposure_timestamp(hdr, 30)['selected']
    assert abs((result['_instant']-Time('2021-06-12T03:00:00.123456')).sec) < 1e-7


def test_timezone_date_preserves_nanosecond_digits():
    hdr=fits.Header({'DATE-AVG':'2021-06-11T20:00:00.123456789-0700'})
    row=select_exposure_timestamp(hdr,30)['selected']
    expected=Time('2021-06-12T03:00:00.123456789',scale='utc')
    assert abs((row['_instant']-expected).sec) < 1e-9
    assert row['precision_seconds']==pytest.approx(1e-9)


def test_bjd_tdb_is_preserved_and_generic_bjd_is_diagnosed():
    hdr = fits.Header({'BJD_TDB': 2459377.72331181})
    chosen = select_exposure_timestamp(hdr, 60)['selected']
    assert chosen['frame'] == 'BJD'
    assert chosen['scale'] == 'TDB'
    assert chosen['_instant'].tdb.jd == hdr['BJD_TDB']
    generic = select_exposure_timestamp(fits.Header({'BJD': 2459377.72331181}), 60)
    assert generic['selected']['ambiguous_standard'] is True


def test_barycentric_start_comment_gets_exposure_midpoint_offset():
    hdr=fits.Header()
    hdr['BJD_TDB']=(2459377.7,'Barycentric time at exposure start')
    hdr['EXPTIME']=60.
    row=select_exposure_timestamp(hdr,60)['selected']
    assert row['frame']=='BJD'
    assert (row['_instant'].tdb.jd-hdr['BJD_TDB'])*86400==pytest.approx(30,abs=.0001)


def test_verified_raw_timestamp_preferred_to_ambiguous_bjd():
    hdr = fits.Header({'BJD': 2459377.723311811, 'DATE-AVG': '2021-06-12T03:00:00.123'})
    assert select_exposure_timestamp(hdr, 60)['selected']['source'] == 'DATE-AVG'


def test_explicit_bjd_utc_converts_clock_without_repeating_light_time():
    hdr = fits.Header({'BJD_UTC': 2459377.72331181})
    row = select_exposure_timestamp(hdr, 60)['selected']
    assert row['frame'] == 'BJD'
    assert (row['_instant'].tdb.jd-row['_instant'].utc.jd)*86400 == pytest.approx(69.18466, abs=.0001)


def test_timeref_barycentric_date_has_no_site_correction():
    hdr = fits.Header({'DATE-AVG': '2021-06-12T03:00:00.123', 'TIMEREF': 'SOLARSYSTEM', 'TIMESYS': 'TDB'})
    row = select_exposure_timestamp(hdr, 60)['selected']
    assert row['frame'] == 'BJD'
    assert row['scale'] == 'TDB'


def test_raw_tt_is_normalised_to_utc_before_barycentric_conversion():
    hdr = fits.Header({'JD-MID': 2459377.72331181, 'TIMESYS': 'TT'})
    row = select_exposure_timestamp(hdr, 60)['selected']
    assert (row['jd_in_source_scale']-row['jd_utc'])*86400 == pytest.approx(69.184, abs=.0001)


def test_invalid_alternative_is_recorded_and_does_not_block_valid_time():
    hdr = fits.Header({'BJD_TDB': 'unknown', 'DATE-AVG': '2021-06-12T03:00:00.123'})
    result = select_exposure_timestamp(hdr, 60)
    assert result['selected']['source'] == 'DATE-AVG'
    assert result['issues'][0]['source'] == 'BJD_TDB'


def test_no_disagreement_rejection():
    hdr = fits.Header({'JD-MID': 2440000.5, 'DATE-AVG': '2021-06-12T03:00:00.123'})
    result = select_exposure_timestamp(hdr, 60)
    assert result['status'] == 'selected'
    assert len(result['candidates']) == 2


def test_old_utc_label_does_not_imply_conversion():
    fields = input_ephemeris_metadata({'Published Mid-Transit Time (BJD-UTC)': 2456355.79653})
    assert fields['midTStandard'] is None
    params = {'midT': 2456355.79653, 'pPer': 3.11860349, 'pName': 'WASP-16 b', **fields}
    called = []
    def lookup(*args):
        called.append(args)
        return {'midT': 2456355.79653, 'midTStandard': 'BJD-TDB', 'midTSource': 'Kokori 2023'}
    out = normalise_ephemeris(params, lookup=lookup)
    assert called == [('WASP-16 b', 2456355.79653, 3.11860349)]
    assert out['midT'] == params['midT']
    assert out['ephemeris_timing']['status'] == 'verified'
    assert out['ephemeris_timing']['conversion_applied'] is False


def test_explicit_utc_conversion_is_idempotent_and_errors_retained():
    params = {'midT': 2456355.79653, 'midTStandard': 'BJD_UTC', 'midTUnc': .00013, 'pPer': 3.11860349}
    out = normalise_ephemeris(params)
    assert out['ephemeris_timing']['conversion_seconds'] == pytest.approx(67.18543, abs=.0001)
    assert normalise_ephemeris(out) == out
    assert out['midTUnc'] == params['midTUnc']
    assert out['pPer'] == params['pPer']


def test_unknown_scale_does_not_block_or_invent_random_error():
    params = {'midT': 2456355.79653}
    def lookup(*_):
        raise RuntimeError('archive unavailable')
    out = normalise_ephemeris(params, lookup=lookup)
    assert out['midT'] == params['midT']
    assert out['ephemeris_timing']['status'] == 'unverified_assumed_model_standard'
    assert out['ephemeris_timing']['conversion_seconds'] is None
    assert 'archive unavailable' in out['ephemeris_timing']['lookup_error']


def test_unrelated_archive_epoch_does_not_supply_standard():
    out = normalise_ephemeris({'midT': 2456355.79653}, archive_parameters={'midT': 2450000., 'midTStandard': 'BJD_UTC'})
    assert out['midT'] == 2456355.79653
    assert out['ephemeris_timing']['status'].startswith('unverified')


@pytest.mark.parametrize('label,standard', [('BJD-TDB','BJD_TDB'),('BJD_UTC','BJD_UTC'),('HJD-UTC','HJD_UTC'),('JD-UTC','JD_UTC')])
def test_explicit_input_time_standards(label, standard):
    assert input_ephemeris_metadata({f'Published Mid-Transit Time ({label})': 2456355.8})['midTStandard'] == standard


def test_heliocentric_and_site_references_are_distinct():
    utc = 2459377.72331181
    jd = reference_to_bjd_tdb(utc, 'JD_UTC', 214.68301, -20.27544)
    hjd = reference_to_bjd_tdb(utc, 'HJD_UTC', 214.68301, -20.27544)
    assert 400 < (jd-utc)*86400 < 500
    assert 60 < (hjd-utc)*86400 < 80


def test_fits_trailing_zero_precision_is_retained():
    hdr = fits.Header.fromstring('JD-MID  =      2459377.600000000 / middle time'.ljust(80), sep='')
    selected = select_exposure_timestamp(hdr, 60)['selected']
    assert selected['precision_seconds'] == pytest.approx(.0000864)


@pytest.mark.parametrize('flat_import', [False, True])
@pytest.mark.parametrize('central_geometry', [False, True])
def test_calculated_tmid_difference_fits_both_series_without_changing_production(monkeypatch, flat_import, central_geometry):
    from types import SimpleNamespace
    import importlib.util
    import sys
    from pathlib import Path
    if flat_import:
        spec = importlib.util.spec_from_file_location('timing', Path(__file__).parents[1]/'exotic'/'timing.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    else:
        import exotic.timing as module
    calculate_timestamp_tmid_difference = module.calculate_timestamp_tmid_difference
    # Exercise the duration-prior import in both supported launcher modes.
    monkeypatch.setitem(sys.modules, 'api.elca' if flat_import else 'exotic.api.elca',
                        SimpleNamespace(transit_duration=lambda params: 1.))
    date = 2459377.7 + np.linspace(-.1,.1,81)
    tmid = 2459377.7
    def model(t,p):
        width = .025*(1-.5*p.get('b',0.)**2)
        return 1-.02*np.exp(-.5*((t-p['tmid'])/width)**2)
    params = {'tmid':tmid,'a2':0.}
    bounds = {'tmid':[tmid-.05,tmid+.05]}
    if central_geometry:
        params.update(b=.05,ars=10.,rprs=.1,ecc=0.,omega=90.,inc=89.)
        bounds['b'] = [0.,1.1]
    data = model(date,dict(params,b=0.))
    context={'time':date, 'base_physical':params.copy(), 'data':data,
             'inverse_dataerr':np.full(len(date),1000.),'centered_airmass':np.zeros(len(date)),
             'free_flux_baseline_key':None,'uses_fixed_flux_baseline':True,
             'fixed_flux_baseline_value':1.,'has_free_flux_baseline':False,'duration_prior_valid':True,
             'expected_duration':1.,'sigma_log_duration':.1}
    fit=SimpleNamespace(parameters=params.copy(),orbital_likelihood_context=context,
                        sampled_keys=list(bounds),sample_bounds=bounds,
                        bounds=bounds,_transit_model=model)
    selections=[{'bjd_tdb':float(t),'selected':{'frame':'JD','jd_utc':float(t)},
                 'candidates':[{'source':'DATE-OBS/DATE-END','jd_utc':float(t)},
                               {'source':'MJD-OBS/DATE-END','jd_utc':float(t+10/86400)}],
                 'mjd_date_comparison':{'date_obs_bjd_tdb':float(t),'mjd_obs_bjd_tdb':float(t+10/86400)}} for t in date]
    result=calculate_timestamp_tmid_difference(fit,selections)
    assert result['status']=='calculated'
    assert result['date_obs_minus_mjd_obs_seconds']==pytest.approx(-10,abs=.001)
    assert result['date_obs']['optimizer_success'] is True
    assert result['mjd_obs']['optimizer_success'] is True
    if central_geometry:
        assert result['date_obs']['parameters']['b'] < .001
        assert result['mjd_obs']['parameters']['b'] < .001
    assert fit.parameters == params
    np.testing.assert_array_equal(fit.orbital_likelihood_context['time'],date)
    assert result['production_fit_changed'] is False
