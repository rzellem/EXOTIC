"""Comparison stars that drift off the frame in more than 10% of the images
are dropped as references; frames where a kept comp is off the frame are
excluded from its photometry.

Motivating case (Exoplanet Watch, 2026-08 MicroObservatory nights): pointing
drift of ~70-150 px per night walks an edge comparison star out of the frame
for the tail of a series. For those frames the centroid may still fit, but
the science aperture crosses the edge and photutils returns a truncated sum
that is finite and positive, so the flux-based coverage rule never sees it.
"""
import numpy as np

from exotic.exotic import (
    APERTURE_SIGMA_MAX,
    COMPARISON_STAR_MAX_OFF_FRAME_FRACTION,
    apply_comparison_star_off_frame_summary,
    comparison_star_off_frame_mask,
    comparison_star_off_frame_summary,
    comparison_star_out_of_frame_flags,
    psf_quality_mask_for_key,
)
from exotic.plate_status import PlateStatus

FRAME_SHAPE = (500, 650)  # (height, width) of a MicroObservatory frame
SIGMA = 2.0


def _rows(x_values, y_values, sigma=SIGMA):
    rows = np.zeros((len(x_values), 7), dtype=float)
    rows[:, 0] = x_values
    rows[:, 1] = y_values
    rows[:, 2] = 1000.0
    rows[:, 3] = sigma
    rows[:, 4] = sigma
    return rows


def _drifting_rows(frame_count, x_start, drift_per_frame, y=250.0):
    x_values = x_start + drift_per_frame * np.arange(frame_count)
    return _rows(x_values, np.full(frame_count, y))


def test_limit_is_ten_percent():
    assert COMPARISON_STAR_MAX_OFF_FRAME_FRACTION == 0.10


def test_mask_flags_frames_where_widest_aperture_crosses_an_edge():
    margin = APERTURE_SIGMA_MAX * SIGMA
    rows = _rows(
        [325.0, margin + 0.5, margin - 0.5, 650.0 - margin - 0.5, 650.0 - margin + 0.5, 325.0],
        [250.0, 250.0, 250.0, 250.0, 250.0, margin - 0.5],
    )
    mask = comparison_star_off_frame_mask(rows, FRAME_SHAPE)
    assert mask.tolist() == [False, False, True, False, True, True]


def test_mask_counts_plate_status_flags_but_not_other_fit_failures():
    rows = _rows([325.0, np.nan, np.nan], [250.0, np.nan, np.nan])
    flagged = np.array([False, True, False])
    mask = comparison_star_off_frame_mask(rows, FRAME_SHAPE, flagged_frames=flagged)
    # Frame 1: seed fell outside the image (flagged). Frame 2: fit failed for
    # some other reason and is not the drift rule's business.
    assert mask.tolist() == [False, True, False]


def test_comp_off_frame_in_more_than_ten_percent_is_rejected():
    frame_count = 100
    psf_data = {
        'target': _drifting_rows(frame_count, 325.0, 0.0),
        # Drifts 1 px/frame toward the right edge; with sigma 2 the widest
        # aperture is APERTURE_SIGMA_MAX*2 px, and the last 15 frames cross it.
        'comp1': _drifting_rows(frame_count, 650.0 - APERTURE_SIGMA_MAX * SIGMA - 84.5, 1.0),
        'comp2': _drifting_rows(frame_count, 200.0, 1.0),
    }
    summary = comparison_star_off_frame_summary(psf_data, ['comp1', 'comp2'], FRAME_SHAPE)
    assert summary['comp1']['off_frame_count'] == 15
    assert summary['comp1']['rejected']
    assert summary['comp2']['off_frame_count'] == 0
    assert not summary['comp2']['rejected']


def test_comp_off_frame_in_at_most_ten_percent_is_kept_with_those_frames_excluded():
    frame_count = 100
    psf_data = {
        'target': _drifting_rows(frame_count, 325.0, 0.0),
        'comp1': _drifting_rows(frame_count, 650.0 - APERTURE_SIGMA_MAX * SIGMA - 89.5, 1.0),
    }
    summary = comparison_star_off_frame_summary(psf_data, ['comp1'], FRAME_SHAPE)
    assert summary['comp1']['off_frame_count'] == 10
    assert not summary['comp1']['rejected']

    apply_comparison_star_off_frame_summary(psf_data, summary)
    quality = psf_quality_mask_for_key(psf_data, 'comp1', frame_count)
    assert quality[:90].all()
    assert not quality[90:].any()


def test_rejected_comp_has_every_frame_excluded_downstream():
    frame_count = 40
    psf_data = {
        'target': _drifting_rows(frame_count, 325.0, 0.0),
        'comp1': _drifting_rows(frame_count, 20.0, -1.0),
        'comp2': _drifting_rows(frame_count, 400.0, -1.0),
    }
    summary = comparison_star_off_frame_summary(psf_data, ['comp1', 'comp2'], FRAME_SHAPE)
    assert summary['comp1']['rejected']
    apply_comparison_star_off_frame_summary(psf_data, summary)
    assert not psf_quality_mask_for_key(psf_data, 'comp1', frame_count).any()
    # The other comp, the target and keys without a record are untouched.
    assert psf_quality_mask_for_key(psf_data, 'comp2', frame_count).all()
    assert psf_quality_mask_for_key(psf_data, 'target', frame_count).all()


def test_last_comparison_star_is_never_rejected():
    frame_count = 50
    psf_data = {
        'target': _drifting_rows(frame_count, 325.0, 0.0),
        'comp1': _drifting_rows(frame_count, 20.0, -1.0),   # off for most of the series
        'comp2': _drifting_rows(frame_count, 30.0, -1.0),   # off for a bit less of it
    }
    summary = comparison_star_off_frame_summary(psf_data, ['comp1', 'comp2'], FRAME_SHAPE)
    assert summary['comp1']['off_frame_fraction'] > COMPARISON_STAR_MAX_OFF_FRAME_FRACTION
    assert summary['comp2']['off_frame_fraction'] > COMPARISON_STAR_MAX_OFF_FRAME_FRACTION
    assert summary['comp1']['rejected']
    assert not summary['comp2']['rejected']
    assert summary['comp2']['kept_as_last_comparison']


def test_sole_comparison_star_is_kept_even_when_over_the_limit():
    frame_count = 40
    psf_data = {'comp1': _drifting_rows(frame_count, 20.0, -1.0)}
    summary = comparison_star_off_frame_summary(psf_data, ['comp1'], FRAME_SHAPE)
    assert summary['comp1']['off_frame_fraction'] > COMPARISON_STAR_MAX_OFF_FRAME_FRACTION
    assert not summary['comp1']['rejected']
    assert summary['comp1']['kept_as_last_comparison']
    apply_comparison_star_off_frame_summary(psf_data, summary)
    quality = psf_quality_mask_for_key(psf_data, 'comp1', frame_count)
    assert quality[:11].all() and not quality[11:].any()


def test_skip_rejection_keeps_masks_but_rejects_nothing():
    frame_count = 40
    psf_data = {'comp1': _drifting_rows(frame_count, 20.0, -1.0)}
    summary = comparison_star_off_frame_summary(psf_data, ['comp1'], FRAME_SHAPE, skip_rejection=True)
    assert summary['comp1']['off_frame_count'] > 4
    assert not summary['comp1']['rejected']


def test_plate_status_out_of_frame_flags_follow_file_order():
    status = PlateStatus(lambda *args, **kwargs: None)
    files = ['a.fits', 'b.fits', 'c.fits']
    status.initializeFilenames(files)
    status.initializeComparisonStarCount(2)
    status.setCurrentFilename('b.fits')
    status.outOfFrameWarning(2)
    flags = comparison_star_out_of_frame_flags(status, files, ['comp1', 'comp2', 'target'])
    assert set(flags) == {'comp1', 'comp2'}
    assert flags['comp1'].tolist() == [False, False, False]
    assert flags['comp2'].tolist() == [False, True, False]


def test_log_lines_name_the_reason(monkeypatch):
    import exotic.exotic as exotic_module

    lines = []
    monkeypatch.setattr(exotic_module, 'log_info', lambda message, **kwargs: lines.append(message))
    frame_count = 100
    psf_data = {
        'comp1': _drifting_rows(frame_count, 650.0 - APERTURE_SIGMA_MAX * SIGMA - 84.5, 1.0),
        'comp2': _drifting_rows(frame_count, 650.0 - APERTURE_SIGMA_MAX * SIGMA - 96.5, 1.0),
        'comp3': _drifting_rows(frame_count, 325.0, 0.0),
    }
    summary = comparison_star_off_frame_summary(psf_data, ['comp1', 'comp2', 'comp3'], FRAME_SHAPE)
    exotic_module.log_comparison_star_off_frame_summary(summary)
    assert lines == [
        "Comparison star #1 rejected as a comparison star: off the frame in 15 of 100 frame(s) "
        "(15.0%; limit 10.0%). Its photometry will not be used as a reference.",
        "Comparison star #2: off the frame in 3 of 100 frame(s) (3.0%; limit 10.0%); "
        "those frames are excluded from its photometry.",
    ]

    lines.clear()
    exotic_module.log_comparison_star_off_frame_summary(
        comparison_star_off_frame_summary(psf_data, ['comp3'], FRAME_SHAPE)
    )
    assert lines == [
        "Comparison-star drift check: no comparison star left the frame in any of 100 frame(s) (limit 10.0%).",
    ]
