import gzip
import json
from pathlib import Path

import numpy as np
from astropy.io.fits import getheader, writeto
import pytest

from exotic.api.plate_solution import NextAstroPlateSolution, PlateSolution


class DummyResponse:
    def __init__(self, payload=None, status_code=200, text=None, json_error=None):
        self._payload = payload
        self.status_code = status_code
        self.text = text if text is not None else ("" if payload is None else str(payload))
        self._json_error = json_error

    def json(self):
        if self._json_error is not None:
            raise self._json_error
        return self._payload


def _decode_request_body(body, headers):
    encoding = headers["Content-Encoding"]
    if encoding == "gzip":
        return json.loads(gzip.decompress(body).decode("utf-8"))
    if encoding == "zstd":
        zstandard = pytest.importorskip("zstandard")
        return json.loads(zstandard.ZstdDecompressor().decompress(body).decode("utf-8"))
    raise AssertionError(f"Unexpected content encoding: {encoding}")


def _create_test_fits(tmp_path: Path) -> Path:
    image = np.zeros((100, 120), dtype=float)
    image[30, 25] = 10000.0
    image[70, 80] = 8000.0
    image[50, 60] = 7000.0
    fits_path = tmp_path / "image.fits"
    writeto(fits_path, image, overwrite=True)
    return fits_path


def test_generate_source_list(tmp_path):
    fits_path = _create_test_fits(tmp_path)
    solver = NextAstroPlateSolution(file=fits_path, directory=tmp_path)

    source_list = solver._generate_source_list()

    assert source_list is not None
    assert source_list["pixel_indexing"] == "0-based"
    assert source_list["origin"] == "exotic"
    assert len(source_list["x"]) > 0
    assert len(source_list["x"]) == len(source_list["y"]) == len(source_list["flux"])


def test_plate_solution_writes_wcs_file(tmp_path, monkeypatch):
    fits_path = _create_test_fits(tmp_path)
    (tmp_path / "working_artifacts").mkdir()

    def fake_post(url, data, headers, timeout):
        payload = _decode_request_body(data, headers)
        assert url.endswith('/solve')
        assert headers["Content-Type"] == "application/json"
        assert headers["Content-Encoding"] in {"gzip", "zstd"}
        assert headers["X-NextAstro-Software"].startswith("EXOTIC/")
        assert payload["sources"]["origin"] == "exotic"
        assert payload['image'] == {'width': 120, 'height': 100}
        assert payload['hints']['ra_deg'] == 210.8023
        assert payload['hints']['dec_deg'] == 54.3489
        assert payload['hints']['scale_arcsec_per_pix'] == 1.23
        assert payload['hints']['scale_tolerance_frac'] == 0.25
        return DummyResponse({'status': 'queued', 'request_id': 'abc123'})

    def fake_get(url, timeout):
        assert url.endswith('/status/abc123')
        return DummyResponse({
            'status': 'solved',
            'solution': {
                'wcs_header': {
                    'SIMPLE': True,
                    'BITPIX': -64,
                    'NAXIS': 2,
                    'NAXIS1': 120,
                    'NAXIS2': 100,
                    'CTYPE1': 'RA---TAN',
                    'CTYPE2': 'DEC--TAN',
                    'CRVAL1': 210.8,
                    'CRVAL2': 54.3,
                    'CRPIX1': 60.0,
                    'CRPIX2': 50.0,
                    'CD1_1': -0.00028,
                    'CD1_2': 0.0,
                    'CD2_1': 0.0,
                    'CD2_2': 0.00028,
                }
            }
        })

    monkeypatch.setattr('exotic.api.plate_solution.requests.post', fake_post)
    monkeypatch.setattr('exotic.api.plate_solution.requests.get', fake_get)
    monkeypatch.setattr('exotic.api.plate_solution.time.sleep', lambda _: None)

    solver = NextAstroPlateSolution(file=fits_path, directory=tmp_path, ra=210.8023, dec=54.3489, pixel_scale=1.23)
    wcs_file = solver.plate_solution()

    assert wcs_file == tmp_path / 'working_artifacts' / 'wcs.fits'
    header = getheader(wcs_file)
    assert header['CTYPE1'] == 'RA---TAN'
    assert header['CTYPE2'] == 'DEC--TAN'


def test_plate_solution_logs_json_via_message_logger_when_fail_warnings_suppressed(tmp_path, monkeypatch):
    fits_path = _create_test_fits(tmp_path)
    (tmp_path / "working_artifacts").mkdir()
    logged = []

    monkeypatch.setattr(
        'exotic.api.plate_solution.requests.post',
        lambda url, data, headers, timeout: DummyResponse({'status': 'queued', 'request_id': 'abc123'})
    )
    monkeypatch.setattr(
        'exotic.api.plate_solution.requests.get',
        lambda url, timeout: DummyResponse({
            'status': 'solved',
            'solution': {
                'wcs_header': {
                    'SIMPLE': True,
                    'BITPIX': -64,
                    'NAXIS': 2,
                    'NAXIS1': 120,
                    'NAXIS2': 100,
                    'CTYPE1': 'RA---TAN',
                    'CTYPE2': 'DEC--TAN',
                }
            }
        })
    )
    monkeypatch.setattr('exotic.api.plate_solution.time.sleep', lambda _: None)

    solver = NextAstroPlateSolution(
        file=fits_path,
        directory=tmp_path,
        suppress_fail_warning=True,
        message_logger=logged.append
    )

    wcs_file = solver.plate_solution()

    assert wcs_file == tmp_path / 'working_artifacts' / 'wcs.fits'
    assert any('NextAstro astrometry request JSON:' in message for message in logged)
    assert any('NextAstro astrometry request compression:' in message for message in logged)
    assert any('NextAstro astrometry submission response JSON:' in message for message in logged)
    assert any('NextAstro astrometry status response JSON (solved):' in message for message in logged)


def test_extract_astrometry_hints_with_scale_only(tmp_path):
    fits_path = _create_test_fits(tmp_path)
    solver = NextAstroPlateSolution(file=fits_path, directory=tmp_path, pixel_scale=2.0)

    hints = solver._extract_astrometry_hints()

    assert hints == {'scale_arcsec_per_pix': 2.0, 'scale_tolerance_frac': 0.25}


def test_poll_for_solution_accepts_case_insensitive_running_status(tmp_path, monkeypatch):
    fits_path = _create_test_fits(tmp_path)
    solver = NextAstroPlateSolution(file=fits_path, directory=tmp_path)

    responses = iter([
        DummyResponse({'status': 'RUNNING'}),
        DummyResponse({
            'status': 'solved',
            'solution': {
                'wcs_header': {
                    'SIMPLE': True,
                    'BITPIX': -64,
                    'NAXIS': 2,
                    'NAXIS1': 120,
                    'NAXIS2': 100,
                    'CTYPE1': 'RA---TAN',
                    'CTYPE2': 'DEC--TAN',
                }
            }
        })
    ])

    monkeypatch.setattr('exotic.api.plate_solution.requests.get', lambda url, timeout: next(responses))
    monkeypatch.setattr('exotic.api.plate_solution.time.sleep', lambda _: None)

    header = solver._poll_for_solution('abc123')

    assert header is not False
    assert header['CTYPE1'] == 'RA---TAN'


def test_poll_for_solution_logs_unexpected_status(tmp_path, monkeypatch, capsys):
    fits_path = _create_test_fits(tmp_path)
    solver = NextAstroPlateSolution(file=fits_path, directory=tmp_path)

    monkeypatch.setattr('exotic.api.plate_solution.requests.get',
                        lambda url, timeout: DummyResponse({'status': 'processing'}))

    header = solver._poll_for_solution('abc123')

    assert header is False
    output = capsys.readouterr().out
    assert "NextAstro astrometry status response JSON (unexpected)" in output
    assert '"status": "processing"' in output


def test_submit_solve_request_handles_non_json_response(tmp_path, monkeypatch, capsys):
    fits_path = _create_test_fits(tmp_path)
    solver = NextAstroPlateSolution(file=fits_path, directory=tmp_path)

    monkeypatch.setattr(
        'exotic.api.plate_solution.requests.post',
        lambda url, data, headers, timeout: DummyResponse(
            status_code=502,
            text='<html>bad gateway</html>',
            json_error=ValueError('not json')
        )
    )

    request_id = solver._submit_solve_request({
        'x': [25.0],
        'y': [30.0],
        'flux': [10000.0],
        'pixel_indexing': '0-based'
    })

    assert request_id is False
    output = capsys.readouterr().out
    assert "NextAstro astrometry request JSON:" in output
    assert "Solve response returned non-JSON response" not in output
    assert "bad gateway" not in output.lower()
    assert solver.last_http_status == 502


def test_poll_for_solution_handles_non_json_response(tmp_path, monkeypatch, capsys):
    fits_path = _create_test_fits(tmp_path)
    solver = NextAstroPlateSolution(file=fits_path, directory=tmp_path)

    monkeypatch.setattr(
        'exotic.api.plate_solution.requests.get',
        lambda url, timeout: DummyResponse(
            status_code=200,
            text='',
            json_error=ValueError('not json')
        )
    )

    header = solver._poll_for_solution('abc123')

    assert header is False
    output = capsys.readouterr().out
    assert "Status response returned non-JSON response" in output
    assert "<empty response body>" in output


def test_nova_upload_includes_astrometry_hints(tmp_path, monkeypatch):
    fits_path = _create_test_fits(tmp_path)

    captured_payload = {}

    def fake_post(url, files, data, timeout):
        captured_payload['url'] = url
        captured_payload['request_json'] = data['request-json']
        return DummyResponse({'status': 'success', 'subid': 42})

    monkeypatch.setattr('exotic.api.plate_solution.requests.post', fake_post)

    solver = PlateSolution(file=fits_path, directory=tmp_path, ra=150.123, dec=-2.456,
                           pixel_scale=1.5, radius=1.2, scale_err=30)
    sub_id = solver._upload(session='session-id')

    assert sub_id == 42
    assert captured_payload['url'].endswith('/upload')
    assert '"session": "session-id"' in captured_payload['request_json']
    assert '"center_ra": 150.123' in captured_payload['request_json']
    assert '"center_dec": -2.456' in captured_payload['request_json']
    assert '"radius": 1.2' in captured_payload['request_json']
    assert '"scale_units": "arcsecperpix"' in captured_payload['request_json']
    assert '"scale_type": "ev"' in captured_payload['request_json']
    assert '"scale_est": 1.5' in captured_payload['request_json']
    assert '"scale_err": 30.0' in captured_payload['request_json']
