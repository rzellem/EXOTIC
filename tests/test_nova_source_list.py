import io
import json
from pathlib import Path

import numpy as np
from astropy.io.fits import getdata, getheader, writeto
from astropy.io import fits

from exotic.api.plate_solution import NextAstroPlateSolution, NovaSourceListPlateSolution, PlateSolution


class DummyResponse:
    def __init__(self, payload, status_code=200, content=b''):
        self._payload = payload
        self.status_code = status_code
        self.content = content

    def json(self):
        return self._payload


def _create_test_fits(tmp_path: Path) -> Path:
    image = np.zeros((100, 120), dtype=float)
    image[30, 25] = 10000.0
    image[70, 80] = 8000.0
    image[50, 60] = 7000.0
    fits_path = tmp_path / "image.fits"
    writeto(fits_path, image, overwrite=True)
    return fits_path


def _wcs_bytes():
    header = fits.Header()
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['CRVAL1'] = 150.0
    header['CRVAL2'] = -2.0
    header['CRPIX1'] = 60.0
    header['CRPIX2'] = 50.0
    header['CDELT1'] = -0.0004
    header['CDELT2'] = 0.0004
    buffer = io.BytesIO()
    fits.PrimaryHDU(header=header).writeto(buffer)
    return buffer.getvalue()


def test_source_list_is_shared_with_nextastro(tmp_path):
    fits_path = _create_test_fits(tmp_path)
    nova = NovaSourceListPlateSolution(file=fits_path, directory=tmp_path)._generate_source_list()
    nextastro = NextAstroPlateSolution(file=fits_path, directory=tmp_path)._generate_source_list()
    assert nova == nextastro
    assert nova["pixel_indexing"] == "0-based"


def test_xylist_is_one_based_with_image_dims_and_private_flags(tmp_path, monkeypatch):
    fits_path = _create_test_fits(tmp_path)
    solver = NovaSourceListPlateSolution(file=fits_path, directory=tmp_path, ra=150.123, dec=-2.456, pixel_scale=1.5)
    source_list = solver._generate_source_list()

    captured = {}

    def fake_post(url, files, data, timeout):
        captured['url'] = url
        captured['request_json'] = json.loads(data['request-json'])
        captured['upload_name'] = Path(files['file'].name).name
        captured['table'] = getdata(files['file'].name, ext=1)
        return DummyResponse({'status': 'success', 'subid': 42})

    monkeypatch.setattr('exotic.api.plate_solution.requests.post', fake_post)

    assert solver._write_source_list() is not None
    assert solver._upload(session='session-id') == 42

    assert captured['url'].endswith('/upload')
    assert captured['upload_name'] == 'nova_sources.xyls'
    table = captured['table']
    assert list(table.columns.names) == ['X', 'Y', 'FLUX']
    np.testing.assert_allclose(table['X'], np.asarray(source_list['x']) + 1.0)
    np.testing.assert_allclose(table['Y'], np.asarray(source_list['y']) + 1.0)
    np.testing.assert_allclose(table['FLUX'], source_list['flux'])

    request_json = captured['request_json']
    assert request_json['session'] == 'session-id'
    assert request_json['image_width'] == 120
    assert request_json['image_height'] == 100
    assert request_json['publicly_visible'] == 'n'
    assert request_json['allow_commercial_use'] == 'n'
    assert request_json['allow_modifications'] == 'n'
    assert request_json['center_ra'] == 150.123
    assert request_json['center_dec'] == -2.456
    assert request_json['scale_est'] == 1.5
    assert request_json['scale_units'] == 'arcsecperpix'


def test_image_upload_puts_visibility_flags_inside_request_json(tmp_path, monkeypatch):
    fits_path = _create_test_fits(tmp_path)
    captured = {}

    def fake_post(url, files, data, timeout):
        captured['form_keys'] = set(data)
        captured['request_json'] = json.loads(data['request-json'])
        captured['upload_name'] = Path(files['file'].name).name
        return DummyResponse({'status': 'success', 'subid': 7})

    monkeypatch.setattr('exotic.api.plate_solution.requests.post', fake_post)

    assert PlateSolution(file=fits_path, directory=tmp_path)._upload(session='s') == 7
    assert captured['upload_name'] == 'image.fits'
    assert captured['form_keys'] == {'request-json'}
    assert captured['request_json']['publicly_visible'] == 'n'
    assert 'image_width' not in captured['request_json']


def test_plate_solution_runs_the_inherited_nova_flow(tmp_path, monkeypatch):
    fits_path = _create_test_fits(tmp_path)
    calls = []

    def fake_post(url, data, timeout, files=None):
        calls.append(url)
        if url.endswith('/login'):
            assert json.loads(data['request-json']) == {'apikey': 'vfsyxlmdxfryhprq'}
            return DummyResponse({'status': 'success', 'session': 'sess'})
        if url.endswith('/upload'):
            assert Path(files['file'].name).name == 'nova_sources.xyls'
            assert json.loads(data['request-json'])['session'] == 'sess'
            return DummyResponse({'status': 'success', 'subid': 5})
        raise AssertionError(url)

    def fake_get(url, timeout):
        calls.append(url)
        if url.endswith('/submissions/5'):
            return DummyResponse({'job_calibrations': [[9, 1]], 'jobs': [9]})
        if url.endswith('/jobs/9'):
            return DummyResponse({'status': 'success'})
        if url.endswith('/wcs_file/9/'):
            return DummyResponse(None, content=_wcs_bytes())
        raise AssertionError(url)

    monkeypatch.setattr('exotic.api.plate_solution.requests.post', fake_post)
    monkeypatch.setattr('exotic.api.plate_solution.requests.get', fake_get)

    solver = NovaSourceListPlateSolution(file=fits_path, directory=tmp_path, suppress_fail_warning=True)
    wcs_file = solver.plate_solution()

    assert wcs_file == tmp_path / 'working_artifacts' / 'wcs.fits'
    assert (tmp_path / 'working_artifacts' / 'nova_sources.xyls').exists()
    assert getheader(wcs_file)['CTYPE1'] == 'RA---TAN'
    assert getdata(wcs_file).shape == (100, 120)
    assert [c.rsplit('/api/', 1)[-1] if '/api/' in c else c for c in calls] == [
        'login', 'upload', 'submissions/5', 'jobs/9', 'http://nova.astrometry.net/wcs_file/9/']
    assert solver.last_error_type is None


def test_blank_frame_fails_before_any_request(tmp_path, monkeypatch):
    fits_path = tmp_path / "blank.fits"
    writeto(fits_path, np.zeros((50, 50), dtype=float), overwrite=True)

    def no_network(*args, **kwargs):
        raise AssertionError('no request expected')

    monkeypatch.setattr('exotic.api.plate_solution.requests.post', no_network)
    monkeypatch.setattr('exotic.api.plate_solution.requests.get', no_network)

    solver = NovaSourceListPlateSolution(file=fits_path, directory=tmp_path, suppress_fail_warning=True)
    assert solver.plate_solution() is False
    assert solver.last_error_type == 'Source extraction'
