"""Regression coverage for epoch values travelling with their time standards."""
import json
import ast
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.time import Time

import exotic.exotic as reduction
from exotic.exotic_gui import ephemeris_export_fields
from exotic.inputs import Inputs, data_file_time, parse_aavso_time_format_from_metadata
from exotic.api.nea import NASAExoplanetArchive
from exotic.timing import (
    EPHEMERIS_EPOCH_KEYS, input_ephemeris_metadata, normalise_ephemeris, time_standard,
)
from exotic.transit_depth import planet_dict_transit_parameters


def parse_parameters(tmp_path, parameters, archive=None):
    path = tmp_path/'inits.json'
    path.write_text(json.dumps({'user_info': {}, 'optional_info': {}, 'planetary_parameters': parameters}))
    return Inputs(init_opt='y').comp_params(path, dict(archive or {}))


@pytest.mark.parametrize('key', EPHEMERIS_EPOCH_KEYS)
def test_epoch_alias_and_metadata_are_selected_together(tmp_path, key):
    parameters = {key:2457218.1101306}
    result = parse_parameters(tmp_path, parameters)
    assert result['midT'] == parameters[key]
    assert result['midTInputKey'] == key
    assert result['midTStandard'] == input_ephemeris_metadata(parameters)['midTStandard']
    assert planet_dict_transit_parameters(parameters)['tmid'] == parameters[key]


@pytest.mark.parametrize('reverse', [False, True])
def test_neutral_epoch_cannot_acquire_another_alias_standard(tmp_path, reverse):
    fields = [('Published Mid-Transit Time',2457218.1101306),
              ('Published Mid-Transit Time (BJD_UTC)',2450000.5)]
    result = parse_parameters(tmp_path,dict(reversed(fields) if reverse else fields))
    assert result['midT'] == 2457218.1101306
    assert result['midTStandard'] is None
    assert result['midTInputKey'] == 'Published Mid-Transit Time'


def test_epoch_metadata_is_not_erased_without_an_input_epoch(tmp_path):
    archive = {'midT':2457218.1101306, 'midTStandard':'BJD-TDB', 'midTSource':'Published epoch'}
    result = parse_parameters(tmp_path,{},archive)
    assert all(result[key] == value for key,value in archive.items())


def test_new_input_epoch_discards_old_normalisation_cache(tmp_path):
    archive = normalise_ephemeris({'midT':2456355.79653, 'midTStandard':'BJD_TDB'})
    result = parse_parameters(tmp_path,{'Published Mid-Transit Time (BJD_UTC)':2457218.1101306},archive)
    assert 'ephemeris_timing' not in result
    normalized = normalise_ephemeris(result)
    assert normalized['ephemeris_timing']['source_epoch'] == 2457218.1101306
    assert normalized['ephemeris_timing']['conversion_seconds'] > 60


@pytest.mark.parametrize('standard', ['BJD-TDB', 'BJD_UTC', None, float('nan')])
def test_gui_export_preserves_epoch_and_standard(tmp_path, standard):
    epoch = 2457218.1101306
    reference = '<a href="https://example.org/paper">Published epoch</a>'
    fields = ephemeris_export_fields({'midT':epoch, 'midTStandard':standard, 'midTSource':reference})
    assert 'Published Mid-Transit Time (BJD-UTC)' not in fields
    result = parse_parameters(tmp_path,fields)
    assert result['midT'] == epoch
    assert result['midTStandard'] == (time_standard(standard) or 'UNKNOWN')
    assert result['midTSource'] == reference


def archive_object(standard='BJD-TDB', reference='Published epoch'):
    obj = NASAExoplanetArchive('Qatar-2 b')
    obj._get_params({'pl_name':'Qatar-2 b','hostname':'Qatar-2','ra':207.65,'dec':-6.804,
                    'pl_ratdor':6.45,'pl_ratror':.18258,'pl_ratrorerr1':.0001,'pl_ratrorerr2':-.0001,
                    'pl_tranmid':2457218.1101306,
                    'pl_tsystemref':standard,'pl_refname':reference})
    return obj


@pytest.mark.parametrize('missing', [None, '', float('nan'), 'NaN'])
def test_archive_missing_metadata_exports_unknown_and_null(missing):
    obj = archive_object(missing,missing)
    result = json.loads(obj.planet_info(fancy=True))
    assert result['Published Mid-Transit Time Standard'] == 'UNKNOWN'
    assert result['Published Mid-Transit Time Reference'] is None


def test_archive_standard_lookup_recognizes_equivalent_spellings(monkeypatch):
    obj = archive_object()
    csv = ('pl_tranmid,pl_orbper,pl_tsystemref,pl_refname\n'
           '2457218.1101306,1.33711644,BJD-TDB,first\n'
           '2457218.1101306,1.33711644,BJD_TDB,second\n'
           '2457218.1101306,1.33711644,UNKNOWN,third\n')
    monkeypatch.setattr(obj,'_tap_query',lambda *_args, **_kwargs:csv)
    result = obj.lookup_ephemeris_time_standard('Qatar-2 b',2457218.1101306,1.33711644)
    assert result['midTStandard'] == 'BJD-TDB'
    assert result['midTSource'] == 'first'


def test_archive_csv_keeps_full_epoch_precision(monkeypatch):
    obj = archive_object()
    response = SimpleNamespace(status_code=200,headers={'content-type':'text/csv'},
                               text='pl_tranmid\n2459377.6483981977\n')
    monkeypatch.setattr('exotic.api.nea.requests.get',lambda *_args,**_kwargs:response)
    result = obj._tap_query('https://example.org/?query=',{'select':'pl_tranmid','from':'ps'})
    assert result['pl_tranmid'][0] == float('2459377.6483981977')


def test_adopting_archive_epoch_replaces_its_metadata(monkeypatch):
    obj = archive_object()
    source = obj.pl_dict.copy()
    for key in source:
        if isinstance(source[key],float) and np.isnan(source[key]):
            source[key] = .01
    supplied = source.copy()
    supplied.update(midT=2450000.5, midTStandard='BJD_UTC', midTSource='old',
                    ephemeris_timing={'normalised':True})
    calls=[]
    def answer(prompt,**kwargs):
        calls.append(prompt)
        return 1
    monkeypatch.setattr(reduction,'user_input',answer)
    reduction.get_planetary_parameters(False,supplied,source)
    assert len(calls) == 1
    assert supplied['midT'] == source['midT']
    assert supplied['midTStandard'] == source['midTStandard']
    assert supplied['midTSource'] == source['midTSource']
    assert 'ephemeris_timing' not in supplied


def test_numeric_repairs_preserve_null_and_structured_timing_metadata():
    timing = {'normalised':True, 'source_standard':'BJD_TDB', 'source_epoch':2457218.1101306}
    parameters = {'rprs':.18258, 'rprsUnc':0., 'omega':None,
                  'midTStandard':'BJD_TDB', 'midTSource':None, 'midTLegacyLabel':None,
                  'midTInputKey':None, 'ephemeris_timing':timing.copy()}
    reduction.repair_missing_planetary_numeric_values(parameters)
    assert parameters['rprsUnc'] == 1.
    assert parameters['omega'] == 0.
    assert parameters['midTSource'] is None
    assert parameters['midTLegacyLabel'] is None
    assert parameters['midTInputKey'] is None
    assert parameters['midTStandard'] == 'BJD_TDB'
    assert parameters['ephemeris_timing'] == timing


@pytest.mark.parametrize('name', ['BJD_TDB','bjd-tdb',' BJD (TDB) ','BJD_UTC','bjd-utc','JD-UTC','MJD-UTC'])
def test_prereduced_time_spellings_do_not_reprompt(monkeypatch,name):
    monkeypatch.setattr('exotic.inputs.user_input',lambda *_args,**_kwargs:pytest.fail('recognized time standard'))
    expected = time_standard(name)
    assert data_file_time(name) == expected
    assert parse_aavso_time_format_from_metadata({'DATE_TYPE':name}) == expected


@pytest.mark.parametrize('standard',['BJD-TDB','BJD_UTC'])
def test_existing_bjd_only_changes_clocks_without_site_or_target(monkeypatch,standard):
    monkeypatch.setattr(reduction,'convert_jd_to_bjd',lambda *_args:pytest.fail('no second barycentric correction'))
    values = np.array([2457218.1101306,2459377.72331181])
    result = reduction.convert_prereduced_to_bjd_tdb(values,{}, {'file_time':standard})
    expected = values if standard == 'BJD-TDB' else Time(values,format='jd',scale='utc').tdb.jd
    np.testing.assert_array_equal(result,expected)
    np.testing.assert_array_equal(values,np.array([2457218.1101306,2459377.72331181]))


@pytest.mark.parametrize('standard', ['JD_UTC','MJD_UTC'])
def test_jd_mjd_still_apply_one_barycentric_conversion(monkeypatch,standard):
    values = np.array([2457218.1101306])
    source = values-2400000.5 if standard == 'MJD_UTC' else values.copy()
    calls=[]
    def convert(times,p,info):
        calls.append(np.asarray(times))
        return np.asarray(times)+.005
    monkeypatch.setattr(reduction,'convert_jd_to_bjd',convert)
    result = reduction.convert_prereduced_to_bjd_tdb(source,{}, {'file_time':standard})
    assert len(calls) == 1
    np.testing.assert_array_equal(calls[0],values)
    np.testing.assert_array_equal(result,values+.005)


@pytest.mark.parametrize('script', ['for_exotic_py_candidate_inits_maker.py','for_toi_py_candidate_inits_maker.py'])
@pytest.mark.parametrize('standard', [None, float('nan'), 'BJD-TDB'])
def test_candidate_writers_preserve_declared_clock_without_inventing_utc(tmp_path,script,standard):
    path = Path(__file__).parents[1]/'examples'/'tess'/'candidates'/script
    tree = ast.parse(path.read_text())
    # These examples execute network and interactive work at module scope.
    # Load their actual writer functions without running those unrelated steps.
    nodes = [node for node in tree.body if isinstance(node,ast.FunctionDef) and
             node.name in ('extract_host_star_name','create_inits_file')]
    namespace = {'np':np,'os':os,'json':json}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),namespace)
    params = {'TIC ID':'123.01','RA (hms)':'13:50:37.318272','DEC (dms)':'-06:48:14.65704',
              'Per (days)':1.33711644,'Per_err (days)':3.2e-8,'Epoch (BJD)':2457218.1101306,
              'Epoch_err (BJD)':6.3e-6,'Epoch Time Standard':standard,'Epoch Reference':None,
              'Teff (K)':4645.,'Teff_err (K)':50.,'logg':4.601,'logg_err':.018}
    destination = tmp_path/'inits.json' if script.startswith('for_exotic_') else tmp_path
    result = namespace['create_inits_file'](params,str(destination))
    if script.startswith('for_exotic_'):
        fields = result['planetary_parameters']
        assert fields['Published Mid-Transit Time'] == params['Epoch (BJD)']
        assert fields['Published Mid-Transit Time Standard'] == (time_standard(standard) or 'UNKNOWN')
        assert fields['Published Mid-Transit Time Reference'] is None
        assert 'Published Mid-Transit Time (BJD-UTC)' not in fields
    else:
        # This writer saves a NASA-style parameter dictionary.
        fields = json.loads(next(tmp_path.glob('*.json')).read_text())
        assert fields['pl_tranmid'] == params['Epoch (BJD)']
        assert fields['pl_tsystemref'] == (time_standard(standard) or 'UNKNOWN')
        assert fields['pl_refname'] is None
