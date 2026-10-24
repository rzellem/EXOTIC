"""Initialization-file pre-flight (``exotic -pf inits.json``).

Each check targets a mistake that -ov otherwise lets through until the
reduction has already been run: a template copy's ephemeris, a transit
outside the window, a seed pixel off the target, comparison stars off the
frame or saturated.
"""
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import pytest

from exotic.api import preflight
from exotic.api.preflight import (
    FAIL, LOOK, PASS, SKIP,
    PreflightReport,
    check_comparison_stars,
    check_target_pixel,
    check_transit_window,
    compare_archive_parameters,
    predicted_mid_transit,
)

WASP11 = {'pName': 'WASP-11 b', 'pPer': 3.72247975, 'midT': 2456933.615316, 'rprs': 0.1413, 'aRs': 12.3,
          'inc': 89.1, 'ecc': 0.0, 'omega': 0.0, 'ra': 47.368958, 'dec': 30.673389}
HATP32_SAMPLE_PERIOD = 2.1500082  # the shipped sample's period, the classic template-copy leftover


def statuses(report):
    return {label: status for status, label, _ in report.checks}


# --- archive agreement -----------------------------------------------------

def test_archive_agreement_passes_when_values_match():
    report = compare_archive_parameters(dict(WASP11), dict(WASP11), PreflightReport())
    assert report.ok
    assert set(statuses(report).values()) == {PASS}


def test_template_period_fails_loudly():
    mine = dict(WASP11, pPer=HATP32_SAMPLE_PERIOD)
    report = compare_archive_parameters(mine, dict(WASP11), PreflightReport())
    assert report.failures == ['Orbital Period (days)']


def test_geometry_tolerance_is_looser_than_period():
    mine = dict(WASP11, rprs=WASP11['rprs'] * 1.10, aRs=WASP11['aRs'] * 0.90, inc=WASP11['inc'] * 1.01)
    assert compare_archive_parameters(mine, dict(WASP11), PreflightReport()).ok
    mine = dict(WASP11, pPer=WASP11['pPer'] * 1.001)
    assert not compare_archive_parameters(mine, dict(WASP11), PreflightReport()).ok


def test_missing_archive_value_is_skipped_not_failed():
    archive = dict(WASP11, aRs=np.nan)
    report = compare_archive_parameters(dict(WASP11), archive, PreflightReport())
    assert report.ok
    assert statuses(report)['Ratio of Distance to Stellar Radius (a/Rs)'] == SKIP


def test_missing_inits_value_fails():
    mine = dict(WASP11, inc=None)
    report = compare_archive_parameters(mine, dict(WASP11), PreflightReport())
    assert report.failures == ['Orbital Inclination (deg)']


# --- timing ------------------------------------------------------------------

# MicroObservatory WASP-11 b night of 2026-09-05: 07:42-12:21 UT, Tmid 09:59:56 UT.
NIGHT_START, NIGHT_END = 2461288.8211, 2461289.0146
DURATION_DAYS = 2.6 / 24.0


def test_predicted_epoch_is_the_one_nearest_the_window():
    tmid, cycle = predicted_mid_transit(NIGHT_START, NIGHT_END, WASP11['midT'], WASP11['pPer'])
    assert cycle == 1170
    assert abs(tmid - 2461288.9166) < 1e-3


def test_full_transit_inside_window_passes_cleanly():
    report = check_transit_window(NIGHT_START, NIGHT_END, WASP11['midT'], WASP11['pPer'], DURATION_DAYS,
                                  PreflightReport())
    assert report.ok
    assert statuses(report)['full transit inside the window'] == PASS


def test_partial_transit_is_flagged_not_failed():
    # window ends 20 minutes after mid-transit: egress is cut
    report = check_transit_window(NIGHT_START, 2461288.9166 + 20 / 1440, WASP11['midT'], WASP11['pPer'],
                                  DURATION_DAYS, PreflightReport())
    assert report.ok
    assert statuses(report)['full transit inside the window'] == LOOK


def test_template_period_puts_transit_outside_window():
    report = check_transit_window(NIGHT_START, NIGHT_END, WASP11['midT'], HATP32_SAMPLE_PERIOD, DURATION_DAYS,
                                  PreflightReport())
    assert report.failures == ['predicted transit overlaps the window']


def test_no_duration_falls_back_to_mid_transit_in_window():
    report = check_transit_window(NIGHT_START, NIGHT_END, WASP11['midT'], WASP11['pPer'], np.nan,
                                  PreflightReport())
    assert report.ok
    assert 'predicted mid-transit inside the window' in statuses(report)


def test_unreadable_window_fails():
    report = check_transit_window(np.nan, NIGHT_END, WASP11['midT'], WASP11['pPer'], DURATION_DAYS,
                                  PreflightReport())
    assert report.failures == ['observing window from the FITS headers']


# --- pointing --------------------------------------------------------------------

def tan_wcs(crpix=(208.0, 342.0), scale_arcsec=5.0, ra=WASP11['ra'], dec=WASP11['dec']):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.crpix = [crpix[0] + 1, crpix[1] + 1]  # FITS 1-based
    wcs.wcs.cdelt = [-scale_arcsec / 3600.0, scale_arcsec / 3600.0]
    return wcs


def test_seed_pixel_on_target_passes():
    report = check_target_pixel(tan_wcs(), WASP11['ra'], WASP11['dec'], [208, 342], PreflightReport())
    assert report.ok


def test_seed_pixel_off_target_fails_with_distance():
    report = check_target_pixel(tan_wcs(), WASP11['ra'], WASP11['dec'], [215, 300], PreflightReport())
    assert len(report.failures) == 1
    _, _, detail = report.checks[0]
    assert '42.6 px' in detail


def test_seed_pixel_within_tolerance_passes():
    report = check_target_pixel(tan_wcs(), WASP11['ra'], WASP11['dec'], [212, 345], PreflightReport())
    assert report.ok


def test_missing_coordinates_skip_the_pointing_check():
    report = check_target_pixel(tan_wcs(), None, WASP11['dec'], [208, 342], PreflightReport())
    assert report.ok
    assert statuses(report)['target pixel on the target'] == SKIP


# --- comparison stars ------------------------------------------------------------

def frame(shape=(500, 650), background=420.0, noise=8.0, stars=(), seed=1):
    rng = np.random.default_rng(seed)
    image = background + noise * rng.standard_normal(shape)
    for x, y, peak in stars:
        image[y, x] = background + peak
    return image


TARGET = (208, 342, 900.0)
SATURATION = 4096.0


def test_good_comparisons_pass():
    image = frame(stars=[TARGET, (388, 471, 1500.0), (589, 211, 400.0)])
    report = check_comparison_stars(image, TARGET[:2], [[388, 471], [589, 211]], SATURATION, PreflightReport())
    assert report.ok
    assert LOOK not in statuses(report).values()


def test_comparison_off_frame_fails():
    image = frame(stars=[TARGET])
    report = check_comparison_stars(image, TARGET[:2], [[700, 100]], SATURATION, PreflightReport())
    assert report.failures == ['comp (700, 100) inside the frame']


def test_saturated_comparison_fails_at_overexposure_fraction():
    image = frame(stars=[TARGET, (388, 471, SATURATION * 0.92 - 420.0)])
    report = check_comparison_stars(image, TARGET[:2], [[388, 471]], SATURATION, PreflightReport(),
                                    overexposure_fraction=0.9)
    assert report.failures == ['comp (388, 471) unsaturated']
    report = check_comparison_stars(image, TARGET[:2], [[388, 471]], SATURATION, PreflightReport(),
                                    overexposure_fraction=0.95)
    assert report.ok


def test_blank_sky_comparison_is_flagged_for_a_look():
    image = frame(stars=[TARGET])
    report = check_comparison_stars(image, TARGET[:2], [[100, 100]], SATURATION, PreflightReport())
    assert report.ok
    assert statuses(report)['comp (100, 100) brightness'] == LOOK


def test_comparison_far_brighter_than_target_is_flagged():
    image = frame(stars=[TARGET, (388, 471, 900.0 * 49)])
    report = check_comparison_stars(image, TARGET[:2], [[388, 471]], 65535.0, PreflightReport())
    assert statuses(report)['comp (388, 471) brightness'] == LOOK


def test_target_off_frame_fails():
    report = check_comparison_stars(frame(), (900, 900), [[388, 471]], SATURATION, PreflightReport())
    assert 'target pixel inside the frame' in report.failures


def test_no_saturation_value_skips_saturation_only():
    image = frame(stars=[TARGET, (388, 471, 60000.0)])
    report = check_comparison_stars(image, TARGET[:2], [[388, 471]], None, PreflightReport())
    assert 'comp (388, 471) unsaturated' not in statuses(report)


# --- report and CLI --------------------------------------------------------------

def test_report_summary_counts_failures():
    report = PreflightReport()
    report.passed(True, 'a')
    report.passed(False, 'b')
    report.look(False, 'c')
    assert report.failures == ['b']
    assert report.summary().startswith('1 pre-flight check(s) failed')
    assert [line.strip()[:6] for line in report.lines()] == ['[ok  ]', '[FAIL]', '[look]']


def test_preflight_flag_parses(monkeypatch):
    from exotic import exotic as runtime
    monkeypatch.setattr('sys.argv', ['exotic', '-pf', 'inits.json'])
    assert runtime.parse_args().preflight == 'inits.json'
    monkeypatch.setattr('sys.argv', ['exotic', '-red', 'inits.json'])
    assert runtime.parse_args().preflight is None


# --- end to end on synthetic frames, no network -----------------------------------

def write_night(tmp_path, count=5, target_px=(208, 342), comps=((388, 471), (589, 211)), period=WASP11['pPer']):
    frames = tmp_path / 'frames'
    frames.mkdir()
    tmid, _ = predicted_mid_transit(NIGHT_START, NIGHT_END, WASP11['midT'], WASP11['pPer'])
    times = np.linspace(tmid - 0.08, tmid + 0.08, count)
    stars = [TARGET] + [(x, y, 1200.0) for x, y in comps if x < 650 and y < 500]
    for index, jd in enumerate(times):
        header = fits.Header()
        header['MJD-OBS'] = jd - 2400000.5
        header['EXPTIME'] = 60.0
        header['TELESCOP'] = 'Cecilia'
        header['IM_SCALE'] = 5.0
        fits.PrimaryHDU(data=frame(stars=stars, seed=index).astype(np.float32), header=header).writeto(
            frames / f'frame{index:02d}.fits')
    out = tmp_path / 'out'
    out.mkdir()
    inits = {
        'user_info': {
            'Directory with FITS files': str(frames), 'Directory to Save Plots': str(out),
            'Directory of Flats': None, 'Directory of Darks': None, 'Directory of Biases': None,
            'AAVSO Observer Code (blank if none)': '', 'Secondary Observer Codes (blank if none)': '',
            'Observation date': '2026-09-05', 'Obs. Latitude': '+31.68', 'Obs. Longitude': '-110.88',
            'Obs. Elevation (meters)': 1268, 'Camera Type (CCD or DSLR)': 'CCD', 'Pixel Binning': '2x2',
            'Filter Name (aavso.org/filters)': 'CV', 'Observing Notes': '', 'Plate Solution? (y/n)': 'n',
            'Add Comparison Stars from AAVSO? (y/n)': 'n',
            'Target Star X & Y Pixel': list(target_px),
            'Comparison Star(s) X & Y Pixel': [list(c) for c in comps],
        },
        'planetary_parameters': {
            'Target Star RA': '03:09:28.55', 'Target Star Dec': '+30:40:24.2',
            'Planet Name': 'WASP-11 b', 'Host Star Name': 'WASP-11',
            'Orbital Period (days)': period, 'Orbital Period Uncertainty': 1.3e-07,
            'Published Mid-Transit Time (BJD-UTC)': WASP11['midT'], 'Mid-Transit Time Uncertainty': 5e-05,
            'Ratio of Planet to Stellar Radius (Rp/Rs)': WASP11['rprs'],
            'Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty': 0.0004,
            'Ratio of Distance to Stellar Radius (a/Rs)': WASP11['aRs'],
            'Ratio of Distance to Stellar Radius (a/Rs) Uncertainty': 0.11,
            'Orbital Inclination (deg)': WASP11['inc'], 'Orbital Inclination (deg) Uncertainty': 0.5,
            'Orbital Eccentricity (0 if null)': 0.0, 'Argument of Periastron (deg)': 0.0,
            'Star Effective Temperature (K)': 4980.0, 'Star Effective Temperature (+) Uncertainty': 60.0,
            'Star Effective Temperature (-) Uncertainty': -60.0, 'Star Metallicity ([FE/H])': 0.0,
            'Star Metallicity (+) Uncertainty': 0.2, 'Star Metallicity (-) Uncertainty': -0.2,
            'Star Surface Gravity (log(g))': 4.45, 'Star Surface Gravity (+) Uncertainty': 0.2,
            'Star Surface Gravity (-) Uncertainty': -0.2, 'Star Distance (pc)': 124.73,
            'Star Proper Motion RA (mas/yr)': 3.85579, 'Star Proper Motion DEC (mas/yr)': -44.8259,
        },
        'optional_info': {'Exposure Time (s)': 60.0},
    }
    path = tmp_path / 'inits.json'
    path.write_text(json.dumps(inits))
    return path


@pytest.fixture
def offline_runtime(monkeypatch):
    from exotic import exotic as runtime

    class FakeArchive:
        def __init__(self, planet=None, **kwargs):
            self.planet = planet

        def planet_info(self):
            return self.planet, False, dict(WASP11, pPerUnc=1.3e-7, midTUnc=5e-5)

    monkeypatch.setattr(runtime, 'NASAExoplanetArchive', FakeArchive)
    monkeypatch.setattr(runtime, '_preflight_wcs_for_frame', lambda *args, **kwargs: (tan_wcs(), 'test WCS'))
    return runtime


def test_end_to_end_clean_inits_exits_zero(tmp_path, offline_runtime, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert offline_runtime.run_inits_preflight(str(write_night(tmp_path))) == 0


def test_end_to_end_template_period_exits_one(tmp_path, offline_runtime, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    code = offline_runtime.run_inits_preflight(str(write_night(tmp_path, period=HATP32_SAMPLE_PERIOD)))
    assert code == 1
    out = capsys.readouterr().out
    assert '[FAIL] Orbital Period (days)' in out
    assert '[FAIL] predicted transit overlaps the window' in out


def test_end_to_end_seed_off_target_and_comp_off_frame(tmp_path, offline_runtime, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    code = offline_runtime.run_inits_preflight(
        str(write_night(tmp_path, target_px=(215, 300), comps=((388, 471), (700, 100)))))
    assert code == 1
    out = capsys.readouterr().out
    assert '[FAIL] target pixel on the target' in out
    assert '[FAIL] comp (700, 100) inside the frame' in out
