# ########################################################################### #
#    Copyright (c) 2019-2020, California Institute of Technology.
#    All rights reserved.  Based on Government Sponsored Research under
#    contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.
#
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions
#    are met:
#      1. Redistributions of source code must retain the above copyright
#         notice, this list of conditions and the following disclaimer.
#      2. Redistributions in binary form must reproduce the above copyright
#         notice, this list of conditions and the following disclaimer in
#         the documentation and/or other materials provided with the
#         distribution.
#      3. Neither the name of the California Institute of
#         Technology (Caltech), its operating division the Jet Propulsion
#         Laboratory (JPL), the National Aeronautics and Space
#         Administration (NASA), nor the names of its contributors may be
#         used to endorse or promote products derived from this software
#         without specific prior written permission.
#
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CALIFORNIA
#    INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ########################################################################### #
#    EXOplanet Transit Interpretation Code (EXOTIC)
#    # NOTE: See companion file version.py for version info.
# ########################################################################### #
# -- IMPORTS START ------------------------------------------------------------
# ########## IMPORTS -- PRELOAD ANIMATION START ##########

try:  # animation
    from animate import *
except ImportError:
    from .animate import *


if __name__ == "__main__":
    print("Importing Python Packages - please wait.")
    animate_toggle(True)

# ########## IMPORTS -- PRELOAD ANIMATION END   ##########

# preload to limit import warnings
import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning
warnings.simplefilter('ignore', category=AstropyDeprecationWarning)

# standard imports
import argparse
import copy
import faulthandler
from functools import lru_cache
import inspect
import json
import hashlib
import multiprocessing
import os
import shutil
import sys
import threading
import traceback
from concurrent.futures import ProcessPoolExecutor as _ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from time import sleep, perf_counter
# Image alignment import
import astroalign as aa
aa.PIXEL_TOL = 1
# aa.NUM_NEAREST_NEIGHBORS=10
# astropy imports
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import BoxLeastSquares
from astropy.visualization import astropy_mpl_style
from astropy.wcs import WCS, FITSFixedWarning
# UTC to BJD converter import
from barycorrpy.utc_tdb import JDUTC_to_BJDTDB
# julian conversion imports
import dateutil.parser as dup
import imreg_dft as ird
from pathlib import Path
import logging
from logging.handlers import TimedRotatingFileHandler
from matplotlib.animation import FuncAnimation
# Pyplot imports
import bottleneck as bn
import matplotlib.pyplot as plt
import numpy as np
# photometry
from photutils.aperture import CircularAperture, CircularAnnulus
import re
import requests
# scipy imports
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
from scipy.ndimage import binary_erosion, gaussian_filter, maximum_filter, median_filter
from skimage.registration import phase_cross_correlation
from skimage.transform import SimilarityTransform
# error handling for scraper
from tenacity import RetryError, retry, retry_if_exception, stop_after_attempt, stop_after_delay, wait_fixed
# color, color_demosaicing
from colour_demosaicing import demosaicing_CFA_Bayer_bilinear
# ########## EXOTIC imports ##########
try:  # light curve numerics
    from .api.elca import lc_fitter, transit, get_phase
except ImportError:  # package import
    from api.elca import lc_fitter, transit, get_phase
try:  # output files
    from inputs import Inputs, comparison_star_coords
except ImportError:  # package import
    from .inputs import Inputs, comparison_star_coords
try:  # ld
    from .api.ld import LimbDarkening, ld_re_punct_p
except ImportError:  # package import
    from api.ld import LimbDarkening, ld_re_punct_p
try:  # plate solution
    from .api.plate_solution import NextAstroPlateSolution, PlateSolution
except ImportError:  # package import
    from api.plate_solution import NextAstroPlateSolution, PlateSolution
try:
    from .api.ultranest_utils import (
        get_mpi_status,
        suppress_inherited_tk_cleanup_in_worker,
        suppress_tk_cleanup_during_process_pool,
    )
except ImportError:
    from api.ultranest_utils import (
        get_mpi_status,
        suppress_inherited_tk_cleanup_in_worker,
        suppress_tk_cleanup_during_process_pool,
    )
try:
    from .api.http_compression import build_compressed_json_request
except ImportError:
    from api.http_compression import build_compressed_json_request
try:  # nea
    from .api.nea import NASAExoplanetArchive
except ImportError:  # package import
    from api.nea import NASAExoplanetArchive
try:  # output files
    from output_files import (
        OutputFiles,
        AIDOutputFiles,
        fit_impact_parameter_value_error,
        format_parameter_with_error,
        save_comp_star_calibration_summary,
    )
except ImportError:  # package import
    from .output_files import (
        OutputFiles,
        AIDOutputFiles,
        fit_impact_parameter_value_error,
        format_parameter_with_error,
        save_comp_star_calibration_summary,
    )
try:
    from plate_status import PlateStatus
except ImportError:
    from .plate_status import PlateStatus
try:  # plots
    from plots import plot_fov, plot_centroids, plot_obs_stats, plot_final_lightcurve, plot_flux, \
        plot_stellar_variability, plot_variable_residuals, plot_comp_star_pairwise_matrix, \
        plot_comp_star_calibration_series, plot_individual_comp_star_calibration_series, \
        plot_comp_star_candidate_lightcurve_fits, plot_comp_star_suitability, \
        plot_adaptive_aperture_diagnostics
except ImportError:  # package import
    from .plots import plot_fov, plot_centroids, plot_obs_stats, plot_final_lightcurve, plot_flux, \
        plot_stellar_variability, plot_variable_residuals, plot_comp_star_pairwise_matrix, \
        plot_comp_star_calibration_series, plot_individual_comp_star_calibration_series, \
        plot_comp_star_candidate_lightcurve_fits, plot_comp_star_suitability, \
        plot_adaptive_aperture_diagnostics
try:  # tools
    from utils import filename_date_token, round_to_2, safe_output_filename, user_input
except ImportError: # package import
    from .utils import filename_date_token, round_to_2, safe_output_filename, user_input
try:  # simple version
    from .version import __version__
except ImportError:  # package import
    from version import __version__

animate_toggle()  # CLOSE PRELOAD ANIMATION
# -- IMPORTS END --------------------------------------------------------------

# SETTINGS
plt.style.use(astropy_mpl_style)
# To increase memory allocation for EXOTIC; allows for more fits files
# resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
# ################### END PROPERTIES/SETTINGS ############################### #

# logging -- https://docs.python.org/3/library/logging.html
log = logging.getLogger(__name__)
_RUNTIME_LOGGING_CONFIGURED = False
_EXCEPTION_HOOKS_INSTALLED = False
_UNHANDLED_EXCEPTION_LOGGED = False
_BJD_FALLBACK_WARNING_LOGGED = False
_RUNTIME_FILE_HANDLER_NAME = "exotic-runtime-file"
_RUNTIME_CONSOLE_HANDLER_NAME = "exotic-runtime-console"
_RUNTIME_TRACEBACK_WATCHDOG_SECONDS_ENV = "EXOTIC_RUNTIME_TRACEBACK_WATCHDOG_SECONDS"
_RUNTIME_TRACEBACK_WATCHDOG_DEFAULT_SECONDS = 1800.0
_RUNTIME_TRACEBACK_WATCHDOG_ACTIVE = False
_mid_transit_warning_reported = False
RELATIVE_FLUX_MAX = 2.0  # Legacy threshold retained for compatibility; no longer used as a hard rejection cap.
AIRMASS_FLAT_RANGE_THRESHOLD = 0.05
LIGHTCURVE_MIN_VALID_POINTS = 5
COMPARISON_STAR_MIN_COVERAGE_FRACTION = 0.8
COMPARISON_STAR_MIN_VALID_FRAMES = 5
COMPARISON_STAR_COVERAGE_SIGMA = 3.0
COMPARISON_STAR_COVERAGE_MAX_ITERS = 10
COMPARISON_STAR_SUITABILITY_OUTLIER_SIGMA = 4.25
COMPARISON_STAR_SUITABILITY_MIN_CANDIDATES = 5
COMPARISON_STAR_SUITABILITY_MAX_ITERS = 10
COMPARISON_IMAGE_OUTLIER_SIGMA = COMPARISON_STAR_SUITABILITY_OUTLIER_SIGMA
COMPARISON_IMAGE_OUTLIER_MIN_ACTIVE_STARS = 3
COMPARISON_IMAGE_OUTLIER_MIN_VALID_PAIRS = 2
COMPARISON_IMAGE_OUTLIER_MIN_SCATTER = 1e-4
OUT_OF_TRANSIT_BASELINE_DEPTH_FRACTION = 0.05
FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT = 1.0
ULTRANEST_MIN_NUM_LIVE_POINTS_DEFAULT = 200
ULTRANEST_MIN_NUM_LIVE_POINTS_ENV = "EXOTIC_ULTRANEST_MIN_NUM_LIVE_POINTS"
FAST_ULTRANEST_BEFORE_FINAL_RUN_DEFAULT = True
FAST_ULTRANEST_MAX_BINNED_POINTS = 20
FAST_ULTRANEST_MIN_POINTS_TO_BIN = 60
COMPARISON_PREFLIGHT_FIELD_SCORE_RELATIVE_BAND = 0.25
COMPARISON_PREFLIGHT_FIELD_SCORE_ABSOLUTE_BAND = 2.5e-4
PARTIAL_COVERAGE_RPRS_POSTERIOR_MAX_RETRIES = 1
PARTIAL_COVERAGE_ARS_POSTERIOR_MAX_RETRIES = 0
PARTIAL_COVERAGE_IMPACT_PARAMETER_POSTERIOR_MAX_RETRIES = 0
PROMISING_PARTIAL_COMPARISON_KTMF_MIN = 3.0
SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_DEFAULT = True
SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_ENV = "EXOTIC_SPARSE_POSTERIOR_LIVE_POINT_RETRY"
SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT = 5
SPARSE_POSTERIOR_RETRY_PARAMETER_KEYS = ('rprs', 'tmid', 'ars')
SPARSE_POSTERIOR_MIN_EFFECTIVE_SAMPLES_FLOOR = 1000
SPARSE_POSTERIOR_MIN_EFFECTIVE_SAMPLES_PER_LIVE_POINT = 5.0
SPARSE_POSTERIOR_MIN_OCCUPIED_BINS = 8
SPARSE_POSTERIOR_MIN_OCCUPIED_BIN_FRACTION = 0.65
SPARSE_POSTERIOR_MIN_EFFECTIVE_SAMPLES_PER_OCCUPIED_BIN = 25.0
RPRS_POSTERIOR_MAX_RETRIES_DEFAULT = 5
RPRS_SEARCH_BOUND_MIN = 0.0
RPRS_SEARCH_BOUND_MAX_DEFAULT = 0.5
RPRS_SEARCH_BOUND_ABSOLUTE_MAX = 1.0
RPRS_SEARCH_BOUND_MAX = RPRS_SEARCH_BOUND_MAX_DEFAULT
RPRS_RETRY_MIN_HALF_WIDTH = 0.05
INITIAL_RPRS_BOUND_LOWER_SCALE = 0.0
INITIAL_RPRS_BOUND_UPPER_SCALE = 3.0
ARS_SEARCH_BOUND_MIN = 1e-6
ARS_SEARCH_BOUND_FALLBACK_MAX = 100.0
ARS_POSTERIOR_MAX_RETRIES_DEFAULT = 5
ARS_RETRY_MIN_HALF_WIDTH = 0.0
IMPACT_PARAMETER_POSTERIOR_MAX_RETRIES_DEFAULT = 5
INCLINATION_SEARCH_BOUND_MIN = 0.0
INCLINATION_SEARCH_BOUND_MAX = 90.0
INITIAL_ARS_BOUND_SIGMA_MULTIPLIER = 5.0
INITIAL_ARS_BOUND_FALLBACK_RELATIVE_HALF_WIDTH = 0.25
DURATION_PRIOR_MONTE_CARLO_SAMPLES = 256
DURATION_PRIOR_MIN_VALID_MONTE_CARLO_SAMPLES = 64
DURATION_PRIOR_MIN_RELATIVE_SIGMA = 0.05
DURATION_PRIOR_FALLBACK_RELATIVE_SIGMA = 0.15
FINAL_FIT_TMID_HALF_DURATION_MULTIPLIER = 0.5
EEBLS_DURATION_GRID_SIZE = 15
EEBLS_DURATION_MIN_FRACTION = 0.5
EEBLS_DURATION_MAX_FRACTION = 1.75
EEBLS_TMID_HALF_WIDTH_DURATION_MULTIPLIER = 1.5
EEBLS_MIN_VALID_POINTS = 10
EPHEMERIS_BRACKETED_TMID_HALF_WIDTH_DURATION_MULTIPLIER = 2.0
COMPARISON_STAR_DUPLICATE_DISTANCE_PIXELS = 15.0
ROBUST_FLUX_MIN_FRACTION_OF_MEDIAN = 0.02
ROBUST_FLUX_MIN_POINTS = 20
WCS_REFERENCE_GEOMETRY_TOLERANCE_PIXELS = 5.0
WCS_MIN_GEOMETRY_MATCH_FRACTION = 0.5
TIME_REJECTION_RANGE_DISPLAY_LIMIT = 6
TIME_REJECTION_GROUP_GAP_CADENCE_MULTIPLIER = 2.5
NEXTASTRO_VARIABILITY_MAX_RETRY_ATTEMPTS = 5
NEXTASTRO_VARIABILITY_RETRY_WAIT_SECONDS = 10
NEXTASTRO_VARIABILITY_RETRYABLE_HTTP_STATUS_CODES = {408, 425, 429, 500, 502, 503, 504}
NEXTASTRO_PHOTOMETRY_API_URL = 'https://photometry.nextastro.org'
NEXTASTRO_PHOTOMETRY_COLUMNS = (
    'id', 'source_id', 'ra', 'dec',
    'Bmag', 'err_Bmag', 'Vmag', 'err_Vmag',
    'umag', 'err_umag', 'g', 'dg', 'r', 'dr', 'i', 'di', 'z', 'dz',
)
NEXTASTRO_PHOTOMETRY_FIELD_PADDING_ARCSEC = 30.0
NEXTASTRO_PHOTOMETRY_MATCH_RADIUS_ARCSEC = 30.0
BAD_PIXEL_DETECTION_FRACTION = 0.30
BAD_PIXEL_PRECHECK_MIN_FRAMES = 5
BAD_PIXEL_PROGRESS_LOG_INTERVAL = 25
BAD_PIXEL_OUTLIER_SIGMA = 8.0
BAD_PIXEL_GLOBAL_SIGMA = 3.0
BAD_PIXEL_ISOLATION_SIGMA = 5.0
BAD_PIXEL_ISOLATION_RATIO = 2.0
MAX_MULTIPROCESS_BAD_PIXEL_WORKERS = 8
BAD_PIXEL_COUNTS_FILENAME = "BadPixelDetectionCounts.fits"
BAD_PIXEL_MASK_FILENAME = "BadPixelMask.fits"
BAD_PIXEL_NEIGHBOR_FOOTPRINT = np.array(
    [[1, 1, 1],
     [1, 0, 1],
     [1, 1, 1]],
    dtype=bool,
)
TRANSIT_QC_DELTA_BIC_FAIL_THRESHOLD = 6.0
TRANSIT_QC_DELTA_BIC_PASS_THRESHOLD = 10.0
TRANSIT_QC_MIN_RPRS_SIGMA = 3.0
TRANSIT_QC_MARGINAL_RPRS_SIGMA = 5.0
TRANSIT_QC_MIN_EEBLS_SNR = 4.0
TRANSIT_QC_DURATION_RATIO_MIN = 0.5
TRANSIT_QC_DURATION_RATIO_MAX = 2.0
TRANSIT_QC_DEFAULT_A2_BOUNDS = (-3.0, 3.0)
TRANSIT_QC_USE_DEVIATION_FROM_EXPECTED_DEFAULT = True
TRANSIT_QC_DEVIATION_SIGMA_DEFAULT = 5.0
TRANSIT_QC_KTMF_COMPONENT_MAX_POINTS = {
    'model_evidence': 0.8,
    'deviation_from_expected_value': 1.5,
    'residual_scatter': 0.7,
    'rprs_significance': 0.5,
    'duration_consistency': 0.75,
    'eebls_depth_snr': 0.75,
}


def airmass_span(airmass):
    values = np.asarray(airmass, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.nan
    return float(np.nanmax(finite) - np.nanmin(finite))


def should_skip_airmass_fit(airmass, max_span=AIRMASS_FLAT_RANGE_THRESHOLD):
    span = airmass_span(airmass)
    return np.isfinite(span) and span <= max_span


def annotate_airmass_fit(fit, airmass, skipped, max_span=AIRMASS_FLAT_RANGE_THRESHOLD, note=None):
    if fit is None:
        return

    span = airmass_span(airmass)
    fit.airmass_span = span
    fit.airmass_fit_threshold = max_span
    fit.airmass_fit_skipped = bool(skipped)
    fit.airmass_correction_note = None
    if fit.airmass_fit_skipped:
        if note:
            fit.airmass_correction_note = note
        elif np.isfinite(span):
            fit.airmass_correction_note = (
                f"Skipped (airmass span {span:.4f} <= {max_span:.2f}); no airmass correction applied."
            )
        else:
            fit.airmass_correction_note = "Skipped; no airmass correction applied."


def annotate_out_of_transit_baseline_detrending(
    fit,
    applied,
    note=None,
    slope=None,
    intercept=None,
    pre_points=0,
    post_points=0,
):
    if fit is None:
        return

    fit.oot_baseline_detrending_applied = bool(applied)
    fit.oot_baseline_detrending_note = note
    fit.oot_baseline_slope = slope
    fit.oot_baseline_intercept = intercept
    fit.oot_baseline_pre_points = int(pre_points) if pre_points is not None else 0
    fit.oot_baseline_post_points = int(post_points) if post_points is not None else 0


def annotate_out_of_transit_baseline_parameter_fit(
    fit,
    applied,
    note=None,
    pre_points=0,
    post_points=0,
    a0=None,
    a0_error=None,
    a2=None,
    a2_error=None,
):
    if fit is None:
        return

    fit.oot_baseline_parameter_fit_applied = bool(applied)
    fit.oot_baseline_parameter_fit_note = note
    fit.oot_baseline_parameter_fit_pre_points = int(pre_points) if pre_points is not None else 0
    fit.oot_baseline_parameter_fit_post_points = int(post_points) if post_points is not None else 0
    fit.oot_baseline_parameter_fit_a0 = a0
    fit.oot_baseline_parameter_fit_a0_error = a0_error
    fit.oot_baseline_parameter_fit_a2 = a2
    fit.oot_baseline_parameter_fit_a2_error = a2_error


def annotate_final_fit_prefit_refinement(
    fit,
    applied,
    note=None,
    baseline_duration_multiplier=FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT,
    duration=None,
    original_point_count=None,
    refined_point_count=None,
    trimmed_pre_points=0,
    trimmed_post_points=0,
    original_tmid_bounds=None,
    refined_tmid_bounds=None,
):
    if fit is None:
        return

    fit.prefit_refinement_applied = bool(applied)
    fit.prefit_refinement_note = note
    fit.prefit_refinement_baseline_duration_multiplier = float(baseline_duration_multiplier)
    fit.prefit_refinement_duration = duration
    fit.prefit_refinement_original_point_count = (
        int(original_point_count) if original_point_count is not None else None
    )
    fit.prefit_refinement_point_count = (
        int(refined_point_count) if refined_point_count is not None else None
    )
    fit.prefit_refinement_trimmed_pre_points = int(trimmed_pre_points or 0)
    fit.prefit_refinement_trimmed_post_points = int(trimmed_post_points or 0)
    fit.prefit_refinement_original_tmid_bounds = original_tmid_bounds
    fit.prefit_refinement_tmid_bounds = refined_tmid_bounds


def annotate_nested_tmid_refinement(
    fit,
    applied,
    note=None,
    original_tmid_bounds=None,
    refined_tmid_bounds=None,
):
    if fit is None:
        return

    fit.nested_tmid_refinement_applied = bool(applied)
    fit.nested_tmid_refinement_note = note
    fit.nested_tmid_refinement_original_tmid_bounds = original_tmid_bounds
    fit.nested_tmid_refinement_tmid_bounds = refined_tmid_bounds


def annotate_duration_prior(fit, duration_prior):
    if fit is None:
        return

    summary = duration_prior if isinstance(duration_prior, dict) else {}
    fit.duration_prior_applied = bool(summary.get('applied', False))
    fit.duration_prior_note = summary.get('note')
    fit.duration_prior_expected_duration = coerce_finite_transit_qc_scalar(
        summary.get('expected_duration', np.nan)
    )
    fit.duration_prior_sigma_log = coerce_finite_transit_qc_scalar(
        summary.get('sigma_log_duration', np.nan)
    )
    fit.duration_prior_relative_sigma = coerce_finite_transit_qc_scalar(
        summary.get('relative_sigma', np.nan)
    )
    fit.duration_prior_source = summary.get('source')


def annotate_pre_ultranest_transit_coverage(fit, assessment):
    if fit is None:
        return

    assessment = assessment if isinstance(assessment, dict) else {}
    fit.pre_ultranest_transit_coverage = dict(assessment)
    fit.pre_ultranest_transit_coverage_valid = bool(assessment.get('valid', False))
    fit.pre_ultranest_transit_coverage_status = assessment.get('success_label')
    fit.pre_ultranest_transit_coverage_chance = assessment.get('success_chance')
    fit.pre_ultranest_transit_coverage_expected_successful = bool(
        assessment.get('expected_successful', False)
    )
    fit.pre_ultranest_transit_coverage_note = assessment.get('note')


def annotate_lightcurve_filter_diagnostics(fit, diagnostics):
    if fit is None:
        return

    fit.frame_filter_diagnostics = [dict(diagnostic) for diagnostic in (diagnostics or [])]


def prepend_lightcurve_filter_diagnostic(fit, diagnostic):
    if fit is None or not diagnostic:
        return

    existing = getattr(fit, 'frame_filter_diagnostics', [])
    diagnostics = [dict(diagnostic)]
    if isinstance(existing, list):
        diagnostics.extend(dict(item) for item in existing if isinstance(item, dict))
    fit.frame_filter_diagnostics = diagnostics


def annotate_selected_photometry_debug(
    fit,
    times,
    target_flux,
    comp_flux,
    raw_ratio,
    initial_sigma_keep_mask,
    phase_clip_keep_mask_on_sigma_filtered=None,
):
    if fit is None:
        return

    sigma_keep_mask = np.asarray(initial_sigma_keep_mask, dtype=bool)
    sigma_kept_count = int(np.count_nonzero(sigma_keep_mask))
    if phase_clip_keep_mask_on_sigma_filtered is None:
        phase_keep_mask = np.ones(sigma_kept_count, dtype=bool)
    else:
        phase_keep_mask = np.asarray(phase_clip_keep_mask_on_sigma_filtered, dtype=bool)
        if phase_keep_mask.shape[0] != sigma_kept_count:
            phase_keep_mask = np.ones(sigma_kept_count, dtype=bool)

    fit.selected_photometry_debug = {
        'times': np.asarray(times, dtype=float).copy(),
        'target_flux': np.asarray(target_flux, dtype=float).copy(),
        'comp_flux': np.asarray(comp_flux, dtype=float).copy(),
        'raw_ratio': np.asarray(raw_ratio, dtype=float).copy(),
        'initial_sigma_keep_mask': sigma_keep_mask.copy(),
        'phase_clip_keep_mask_on_sigma_filtered': phase_keep_mask.copy(),
    }


def transit_qc_airmass_reference(airmass):
    values = np.asarray(airmass, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0
    return float(np.nanmean(finite))


def transit_qc_airmass_trend(a2, airmass, reference=None):
    values = np.asarray(airmass, dtype=float)
    if reference is None:
        reference = transit_qc_airmass_reference(values)
    return np.exp(float(a2) * (values - float(reference)))


def solve_transit_qc_flux_baseline(systematics, data, dataerr=None):
    systematics = np.asarray(systematics, dtype=float)
    data = np.asarray(data, dtype=float)
    if systematics.shape != data.shape:
        return np.nan

    weights = np.ones(systematics.shape, dtype=float)
    if dataerr is not None:
        dataerr = np.asarray(dataerr, dtype=float)
        if dataerr.shape != data.shape:
            dataerr = None
        else:
            weights = np.zeros(systematics.shape, dtype=float)
            valid_err = np.isfinite(dataerr) & (dataerr > 0)
            weights[valid_err] = 1.0 / (dataerr[valid_err] ** 2)

    mask = np.isfinite(systematics) & np.isfinite(data) & (systematics != 0)
    if dataerr is not None:
        mask &= np.isfinite(weights) & (weights > 0)

    if not np.any(mask):
        return np.nan

    masked_systematics = systematics[mask]
    masked_data = data[mask]
    masked_weights = weights[mask]
    denom = np.sum(masked_weights * masked_systematics ** 2)
    if np.isfinite(denom) and denom > 0:
        baseline = np.sum(masked_weights * masked_data * masked_systematics) / denom
        if np.isfinite(baseline):
            return float(baseline)

    ratio = masked_data / masked_systematics
    ratio = ratio[np.isfinite(ratio)]
    if ratio.size == 0:
        return np.nan
    return float(np.nanmedian(ratio))


def compute_transit_qc_model_chi2(data, model, dataerr=None):
    data = np.asarray(data, dtype=float)
    model = np.asarray(model, dtype=float)
    if data.shape != model.shape:
        return np.nan, 0

    mask = np.isfinite(data) & np.isfinite(model)
    if dataerr is not None:
        dataerr = np.asarray(dataerr, dtype=float)
        if dataerr.shape != data.shape:
            dataerr = None
        else:
            mask &= np.isfinite(dataerr) & (dataerr > 0)

    point_count = int(np.count_nonzero(mask))
    if point_count == 0:
        return np.nan, 0

    residuals = data[mask] - model[mask]
    if dataerr is not None:
        chi2 = np.sum((residuals / dataerr[mask]) ** 2)
    else:
        chi2 = np.sum(residuals ** 2)
    return float(chi2), point_count


def compute_transit_qc_bic(chi2, point_count, parameter_count):
    try:
        chi2 = float(chi2)
        point_count = int(point_count)
        parameter_count = int(parameter_count)
    except (TypeError, ValueError):
        return np.nan

    if not np.isfinite(chi2) or point_count <= 0 or parameter_count <= 0:
        return np.nan
    return float(chi2 + parameter_count * np.log(point_count))


def clip_unit_interval(value):
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return np.nan

    if not np.isfinite(numeric_value):
        return np.nan
    return float(np.clip(numeric_value, 0.0, 1.0))


def transit_qc_residual_scatter(data, model):
    data = np.asarray(data, dtype=float)
    model = np.asarray(model, dtype=float)
    if data.shape != model.shape or data.size == 0:
        return np.nan

    mask = np.isfinite(data) & np.isfinite(model)
    if not np.any(mask):
        return np.nan

    median_flux = np.nanmedian(data[mask])
    if not np.isfinite(median_flux) or median_flux == 0:
        return np.nan

    residuals = data[mask] - model[mask]
    return float(np.std(residuals) / median_flux)


def transit_qc_deviation_score_from_sigma(sigma_offset, sigma_threshold):
    try:
        sigma_offset = abs(float(sigma_offset))
        sigma_threshold = float(sigma_threshold)
    except (TypeError, ValueError):
        return np.nan

    if not np.isfinite(sigma_offset) or not np.isfinite(sigma_threshold) or sigma_threshold <= 0:
        return np.nan

    return float(max(0.0, 1.0 - sigma_offset / sigma_threshold))


def transit_qc_duration_score(duration_ratio):
    try:
        duration_ratio = float(duration_ratio)
    except (TypeError, ValueError):
        return np.nan

    if not np.isfinite(duration_ratio) or duration_ratio <= 0:
        return np.nan

    max_log_deviation = np.log(TRANSIT_QC_DURATION_RATIO_MAX)
    if not np.isfinite(max_log_deviation) or max_log_deviation <= 0:
        return np.nan

    score = 1.0 - abs(np.log(duration_ratio)) / max_log_deviation
    return float(np.clip(score, 0.0, 1.0))


def transit_qc_saturating_score(value, scale):
    try:
        value = float(value)
        scale = float(scale)
    except (TypeError, ValueError):
        return np.nan

    if not np.isfinite(value) or not np.isfinite(scale) or scale <= 0:
        return np.nan

    return float(np.clip(1.0 - np.exp(-max(value, 0.0) / scale), 0.0, 1.0))


def transit_qc_residual_scatter_score(residual_scatter, reference_scatter=0.005):
    try:
        residual_scatter = float(residual_scatter)
        reference_scatter = float(reference_scatter)
    except (TypeError, ValueError):
        return np.nan

    if not np.isfinite(residual_scatter) or residual_scatter < 0 or not np.isfinite(reference_scatter) or reference_scatter <= 0:
        return np.nan

    return float(np.clip(1.0 / (1.0 + residual_scatter / reference_scatter), 0.0, 1.0))


def transit_qc_mean_available_score(*scores):
    finite_scores = [float(score) for score in scores if np.isfinite(score)]
    if not finite_scores:
        return np.nan
    return float(np.clip(np.mean(finite_scores), 0.0, 1.0))


def coerce_finite_transit_qc_scalar(value):
    if value is None:
        return np.nan

    if isinstance(value, (list, tuple, np.ndarray)):
        array_value = np.asarray(value)
        if array_value.size != 1:
            return np.nan
        value = array_value.reshape(-1)[0]

    try:
        numeric_value = float(value.strip()) if isinstance(value, str) else float(value)
    except (AttributeError, TypeError, ValueError):
        return np.nan

    if not np.isfinite(numeric_value):
        return np.nan
    return float(numeric_value)


def fit_transit_qc_expected_context(fit):
    if fit is None:
        return {}

    return {
        'expected_tmid': coerce_finite_transit_qc_scalar(
            getattr(fit, 'transit_qc_expected_tmid', np.nan)
        ),
        'expected_tmid_unc': coerce_finite_transit_qc_scalar(
            getattr(fit, 'transit_qc_expected_tmid_unc', np.nan)
        ),
        'expected_rprs': coerce_finite_transit_qc_scalar(
            getattr(fit, 'transit_qc_expected_rprs', np.nan)
        ),
        'expected_rprs_unc': coerce_finite_transit_qc_scalar(
            getattr(fit, 'transit_qc_expected_rprs_unc', np.nan)
        ),
        'use_deviation_from_expected_transit_in_qc': should_use_deviation_from_expected_transit_in_qc(
            getattr(
                fit,
                'transit_qc_use_deviation_from_expected_transit_in_qc',
                TRANSIT_QC_USE_DEVIATION_FROM_EXPECTED_DEFAULT,
            )
        ),
        'deviation_sigma_threshold': parse_deviation_from_expected_transit_in_qc_sigma(
            getattr(
                fit,
                'transit_qc_deviation_sigma_threshold',
                TRANSIT_QC_DEVIATION_SIGMA_DEFAULT,
            )
        ),
    }


def annotate_transit_qc_expected_values(fit, planet_dict):
    if fit is None or not isinstance(planet_dict, dict):
        return

    expected_tmid = coerce_finite_transit_qc_scalar(
        getattr(fit, 'initial_tmid_search_tmid', np.nan)
    )
    expected_tmid_unc = coerce_finite_transit_qc_scalar(
        getattr(fit, 'initial_tmid_search_uncertainty', np.nan)
    )
    if not np.isfinite(expected_tmid):
        expected_tmid = coerce_finite_transit_qc_scalar(planet_dict.get('midT', np.nan))
    if not np.isfinite(expected_tmid_unc):
        expected_tmid_unc = coerce_finite_transit_qc_scalar(planet_dict.get('midTUnc', np.nan))

    fit.transit_qc_expected_tmid = expected_tmid
    fit.transit_qc_expected_tmid_unc = expected_tmid_unc
    fit.transit_qc_expected_rprs = coerce_finite_transit_qc_scalar(
        planet_dict.get('rprs', np.nan)
    )
    fit.transit_qc_expected_rprs_unc = coerce_finite_transit_qc_scalar(
        planet_dict.get('rprsUnc', np.nan)
    )
    fit.transit_qc_use_deviation_from_expected_transit_in_qc = should_use_deviation_from_expected_transit_in_qc(
        planet_dict.get(
            'use_deviation_from_expected_transit_in_qc',
            TRANSIT_QC_USE_DEVIATION_FROM_EXPECTED_DEFAULT,
        )
    )
    fit.transit_qc_deviation_sigma_threshold = parse_deviation_from_expected_transit_in_qc_sigma(
        planet_dict.get(
            'deviation_from_expected_transit_in_qc_sigma',
            TRANSIT_QC_DEVIATION_SIGMA_DEFAULT,
        )
    )


def annotate_transit_qc_fit_context(
    fit,
    planet_dict=None,
    tmid_search_summary=None,
    eebls_search_summary=None,
):
    if fit is None:
        return

    if tmid_search_summary is not None:
        annotate_lightcurve_tmid_search(fit, tmid_search_summary)
    if eebls_search_summary is not None:
        annotate_lightcurve_eebls_diagnostic(fit, eebls_search_summary)
    annotate_transit_qc_expected_values(fit, planet_dict)


def copy_transit_qc_expected_values(source_fit, target_fit):
    if source_fit is None or target_fit is None:
        return

    for attr in (
        'transit_qc_expected_tmid',
        'transit_qc_expected_tmid_unc',
        'transit_qc_expected_rprs',
        'transit_qc_expected_rprs_unc',
        'transit_qc_use_deviation_from_expected_transit_in_qc',
        'transit_qc_deviation_sigma_threshold',
    ):
        if hasattr(source_fit, attr):
            setattr(target_fit, attr, getattr(source_fit, attr))


def evaluate_transit_qc_expected_value_deviation(fit, sigma_threshold, enabled=True):
    summary = {
        'enabled': bool(enabled),
        'sigma_threshold': sigma_threshold,
        'expected_tmid': np.nan,
        'expected_tmid_unc': np.nan,
        'expected_tmid_unc_minutes': np.nan,
        'fitted_tmid': np.nan,
        'tmid_deviation_days': np.nan,
        'tmid_deviation_minutes': np.nan,
        'tmid_deviation_threshold_minutes': np.nan,
        'tmid_deviation_sigma': np.nan,
        'rprs_deviation_sigma': np.nan,
        'tmid_deviation_score': np.nan,
        'rprs_deviation_score': np.nan,
        'deviation_from_expected_value': np.nan,
        'available': False,
        'failed': False,
        'notes': [],
        'failure_reasons': [],
    }
    if fit is None:
        return summary

    parameters = getattr(fit, 'parameters', {}) or {}
    expected = fit_transit_qc_expected_context(fit)
    sigma_threshold = expected.get('deviation_sigma_threshold', sigma_threshold)
    summary['sigma_threshold'] = sigma_threshold

    expected_tmid = expected.get('expected_tmid', np.nan)
    expected_tmid_unc = expected.get('expected_tmid_unc', np.nan)
    fitted_tmid = parameters.get('tmid', np.nan)
    summary['expected_tmid'] = expected_tmid
    summary['expected_tmid_unc'] = expected_tmid_unc
    summary['fitted_tmid'] = fitted_tmid
    if not enabled:
        return summary

    if (
        np.isfinite(expected_tmid)
        and np.isfinite(expected_tmid_unc)
        and expected_tmid_unc > 0
        and np.isfinite(fitted_tmid)
    ):
        tmid_deviation_days = float(abs(fitted_tmid - expected_tmid))
        tmid_sigma = float(tmid_deviation_days / expected_tmid_unc)
        summary['tmid_deviation_days'] = tmid_deviation_days
        summary['tmid_deviation_minutes'] = tmid_deviation_days * 24.0 * 60.0
        summary['expected_tmid_unc_minutes'] = float(expected_tmid_unc) * 24.0 * 60.0
        if np.isfinite(sigma_threshold) and sigma_threshold > 0:
            summary['tmid_deviation_threshold_minutes'] = (
                float(sigma_threshold) * float(expected_tmid_unc) * 24.0 * 60.0
            )
        summary['tmid_deviation_sigma'] = tmid_sigma
        summary['tmid_deviation_score'] = transit_qc_deviation_score_from_sigma(tmid_sigma, sigma_threshold)

    expected_rprs = expected.get('expected_rprs', np.nan)
    expected_rprs_unc = expected.get('expected_rprs_unc', np.nan)
    fitted_rprs = parameters.get('rprs', np.nan)
    if (
        np.isfinite(expected_rprs)
        and np.isfinite(expected_rprs_unc)
        and expected_rprs_unc > 0
        and np.isfinite(fitted_rprs)
    ):
        rprs_sigma = float(abs(fitted_rprs - expected_rprs) / expected_rprs_unc)
        summary['rprs_deviation_sigma'] = rprs_sigma
        summary['rprs_deviation_score'] = transit_qc_deviation_score_from_sigma(rprs_sigma, sigma_threshold)

    if np.isfinite(summary['rprs_deviation_score']):
        summary['available'] = True
        summary['deviation_from_expected_value'] = float(summary['rprs_deviation_score'])

    if np.isfinite(summary['tmid_deviation_sigma']):
        summary['notes'].append(
            "Expected-value Tmid deviation: "
            f"{summary['tmid_deviation_minutes']:.2f} minutes "
            f"({summary['tmid_deviation_sigma']:.2f} sigma; "
            f"fit={summary['fitted_tmid']:.6f}, "
            f"ephemeris={summary['expected_tmid']:.6f} +/- "
            f"{summary['expected_tmid_unc_minutes']:.2f} minutes)."
        )
    if np.isfinite(summary['rprs_deviation_sigma']):
        summary['notes'].append(
            f"Expected-value Rp/R* deviation: {summary['rprs_deviation_sigma']:.2f} sigma."
        )

    rprs_sigma = summary['rprs_deviation_sigma']
    if np.isfinite(rprs_sigma) and np.isfinite(sigma_threshold) and sigma_threshold > 0 and rprs_sigma > sigma_threshold:
        summary['failed'] = True
        reason = f"Rp/R* differs from the expected value by more than {sigma_threshold:.2f} sigma"
        summary['notes'].append(reason + ".")
        summary['failure_reasons'].append(reason)

    return summary


def compute_transit_qc_ktmf(summary):
    if not isinstance(summary, dict):
        return np.nan, []

    delta_bic_score = transit_qc_saturating_score(
        summary.get('delta_bic', np.nan),
        TRANSIT_QC_DELTA_BIC_PASS_THRESHOLD,
    )
    delta_chi2_score = transit_qc_saturating_score(summary.get('delta_chi2', np.nan), 25.0)
    model_evidence_score = transit_qc_mean_available_score(delta_bic_score, delta_chi2_score)
    model_evidence_detail_parts = [
        f"Delta BIC={format_transit_delta_bic(summary.get('delta_bic', np.nan))}",
        f"Delta chi2={summary.get('delta_chi2', np.nan):.2f}"
        if np.isfinite(summary.get('delta_chi2', np.nan))
        else "Delta chi2=n/a",
    ]

    raw_components = [
        {
            'key': 'model_evidence',
            'label': 'Model Evidence',
            'score': model_evidence_score,
            'detail': ", ".join(model_evidence_detail_parts),
        },
        {
            'key': 'deviation_from_expected_value',
            'label': 'Deviation From Expected Value',
            'score': summary.get('deviation_from_expected_value', np.nan),
            'detail': (
                f"score={summary.get('deviation_from_expected_value', np.nan):.2f}, "
                f"Rp/R* sigma={summary.get('rprs_deviation_sigma', np.nan):.2f}"
                if np.isfinite(summary.get('deviation_from_expected_value', np.nan))
                else "expected-value deviation disabled or unavailable"
            ),
        },
        {
            'key': 'residual_scatter',
            'label': 'Residual Scatter Around Full Model Fit',
            'score': transit_qc_residual_scatter_score(summary.get('residual_scatter', np.nan)),
            'detail': (
                f"{summary.get('residual_scatter', np.nan) * 100.0:.4f}%"
                if np.isfinite(summary.get('residual_scatter', np.nan))
                else "n/a"
            ),
        },
        {
            'key': 'rprs_significance',
            'label': 'Rp/R* Significance',
            'score': transit_qc_saturating_score(summary.get('rprs_sigma', np.nan), TRANSIT_QC_MIN_RPRS_SIGMA),
            'detail': (
                f"{summary.get('rprs_sigma', np.nan):.2f} sigma"
                if np.isfinite(summary.get('rprs_sigma', np.nan))
                else "n/a"
            ),
        },
        {
            'key': 'duration_consistency',
            'label': 'Duration Consistency',
            'score': transit_qc_duration_score(summary.get('duration_ratio', np.nan)),
            'detail': (
                f"{summary.get('duration_ratio', np.nan):.2f}x expected duration"
                if np.isfinite(summary.get('duration_ratio', np.nan))
                else "n/a"
            ),
        },
        {
            'key': 'eebls_depth_snr',
            'label': 'EEBLS Depth SNR',
            'score': transit_qc_saturating_score(summary.get('eebls_depth_snr', np.nan), TRANSIT_QC_MIN_EEBLS_SNR),
            'detail': (
                f"{summary.get('eebls_depth_snr', np.nan):.2f}"
                if np.isfinite(summary.get('eebls_depth_snr', np.nan))
                else "n/a"
            ),
        },
    ]

    available_components = [
        component
        for component in raw_components
        if np.isfinite(component.get('score', np.nan))
        and component['key'] in TRANSIT_QC_KTMF_COMPONENT_MAX_POINTS
    ]
    if not available_components:
        return np.nan, []

    available_max_points = sum(TRANSIT_QC_KTMF_COMPONENT_MAX_POINTS[component['key']] for component in available_components)
    if not np.isfinite(available_max_points) or available_max_points <= 0:
        return np.nan, []

    scale_factor = 5.0 / available_max_points
    ktmf_contributions = []
    total_points = 0.0
    for component in raw_components:
        nominal_max_points = TRANSIT_QC_KTMF_COMPONENT_MAX_POINTS.get(component['key'], 0.0)
        score = component.get('score', np.nan)
        if np.isfinite(score) and nominal_max_points > 0:
            max_points = nominal_max_points * scale_factor
            points = float(np.clip(score, 0.0, 1.0) * max_points)
            total_points += points
            ktmf_contributions.append({
                'label': component['label'],
                'score': float(np.clip(score, 0.0, 1.0)),
                'max_points': float(max_points),
                'points': points,
                'detail': component.get('detail'),
                'available': True,
            })
        else:
            ktmf_contributions.append({
                'label': component['label'],
                'score': np.nan,
                'max_points': 0.0,
                'points': 0.0,
                'detail': component.get('detail'),
                'available': False,
            })

    return float(np.clip(total_points, 0.0, 5.0)), ktmf_contributions


def infer_transit_qc_parameter_count(fit, allow_airmass_term):
    bounds = getattr(fit, 'bounds', None)
    if isinstance(bounds, dict) and bounds:
        parameter_count = len(bounds)
        if 'a0' not in bounds and 'a1' not in bounds:
            parameter_count += 1
        return max(int(parameter_count), 1)

    parameters = getattr(fit, 'parameters', {}) or {}
    errors = getattr(fit, 'errors', {}) or {}
    available_keys = set(parameters.keys()) | set(errors.keys())

    parameter_count = 1  # profiled flux baseline
    if 'rprs' in available_keys:
        parameter_count += 1
    if 'tmid' in available_keys:
        parameter_count += 1
    if 'inc' in available_keys or 'b' in available_keys:
        parameter_count += 1
    if allow_airmass_term and 'a2' in available_keys:
        parameter_count += 1
    return max(parameter_count, 1)


def fit_profiled_flat_null_model(data, dataerr, airmass, initial_a2=0.0, a2_bounds=None, allow_airmass_term=True):
    data = np.asarray(data, dtype=float)
    dataerr = None if dataerr is None else np.asarray(dataerr, dtype=float)
    airmass = None if airmass is None else np.asarray(airmass, dtype=float)

    result = {
        'available': False,
        'used_airmass_term': False,
        'baseline': np.nan,
        'a2': 0.0,
        'model': np.full(data.shape, np.nan, dtype=float),
        'chi2': np.nan,
        'bic': np.nan,
        'point_count': 0,
        'param_count': 1,
        'note': 'Flat/null model was not evaluated.',
    }

    if data.ndim != 1 or data.size == 0:
        result['note'] = 'Flat/null model comparison unavailable: no 1D lightcurve data were provided.'
        return result

    if dataerr is not None and dataerr.shape != data.shape:
        dataerr = None
    if airmass is not None and airmass.shape != data.shape:
        airmass = None

    use_airmass_term = bool(
        allow_airmass_term
        and airmass is not None
        and airmass.ndim == 1
        and airmass.shape == data.shape
        and not should_skip_airmass_fit(airmass)
    )
    result['used_airmass_term'] = use_airmass_term
    result['param_count'] = 1 + int(use_airmass_term)

    if not use_airmass_term:
        baseline = solve_transit_qc_flux_baseline(np.ones_like(data, dtype=float), data, dataerr)
        if not np.isfinite(baseline):
            result['note'] = 'Flat/null model comparison unavailable: could not solve the baseline flux level.'
            return result
        model = np.full(data.shape, baseline, dtype=float)
        chi2, point_count = compute_transit_qc_model_chi2(data, model, dataerr)
        result.update({
            'available': np.isfinite(chi2),
            'baseline': float(baseline),
            'model': model,
            'chi2': chi2,
            'point_count': point_count,
            'bic': compute_transit_qc_bic(chi2, point_count, result['param_count']),
            'note': 'Compared against a flat baseline-only null model.',
        })
        return result

    lower, upper = TRANSIT_QC_DEFAULT_A2_BOUNDS
    if a2_bounds is not None:
        try:
            lower, upper = np.asarray(a2_bounds, dtype=float).reshape(-1)[:2]
        except (TypeError, ValueError, IndexError):
            lower, upper = TRANSIT_QC_DEFAULT_A2_BOUNDS
    if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
        lower, upper = TRANSIT_QC_DEFAULT_A2_BOUNDS

    if not np.isfinite(initial_a2):
        initial_a2 = 0.0
    initial_a2 = float(np.clip(initial_a2, lower + np.finfo(float).eps, upper - np.finfo(float).eps))
    reference = transit_qc_airmass_reference(airmass)

    def build_model(a2_value):
        systematics = transit_qc_airmass_trend(a2_value, airmass, reference=reference)
        baseline = solve_transit_qc_flux_baseline(systematics, data, dataerr)
        if not np.isfinite(baseline):
            return np.full(data.shape, np.nan, dtype=float), np.nan
        return baseline * systematics, baseline

    def residual_vector(params):
        model, baseline = build_model(params[0])
        if not np.isfinite(baseline):
            return np.full(max(1, data.size), 1e6, dtype=float)

        mask = np.isfinite(data) & np.isfinite(model)
        if dataerr is not None:
            mask &= np.isfinite(dataerr) & (dataerr > 0)
        if not np.any(mask):
            return np.full(max(1, data.size), 1e6, dtype=float)

        if dataerr is not None:
            return (data[mask] - model[mask]) / dataerr[mask]
        return data[mask] - model[mask]

    best_a2 = float(initial_a2)
    try:
        fit_result = least_squares(
            residual_vector,
            x0=np.array([initial_a2], dtype=float),
            bounds=([lower], [upper]),
        )
        if fit_result.x.size:
            best_a2 = float(fit_result.x[0])
    except Exception:
        pass

    model, baseline = build_model(best_a2)
    if not np.isfinite(baseline):
        result['note'] = 'Flat/null model comparison unavailable: the null-model fit did not converge.'
        return result

    chi2, point_count = compute_transit_qc_model_chi2(data, model, dataerr)
    result.update({
        'available': np.isfinite(chi2),
        'baseline': float(baseline),
        'a2': float(best_a2),
        'model': model,
        'chi2': chi2,
        'point_count': point_count,
        'bic': compute_transit_qc_bic(chi2, point_count, result['param_count']),
        'note': 'Compared against a flat null model with the same profiled baseline and airmass trend.',
    })
    return result


def evaluate_transit_detection_qc(fit):
    expected_context = fit_transit_qc_expected_context(fit)
    use_deviation_from_expected_transit_in_qc = expected_context.get(
        'use_deviation_from_expected_transit_in_qc',
        TRANSIT_QC_USE_DEVIATION_FROM_EXPECTED_DEFAULT,
    )
    deviation_sigma_threshold = expected_context.get(
        'deviation_sigma_threshold',
        TRANSIT_QC_DEVIATION_SIGMA_DEFAULT,
    )
    summary = {
        'computed': False,
        'status': 'unknown',
        'preferred_model': 'unknown',
        'summary': 'Transit QC unavailable: no fit result was provided.',
        'notes': [],
        'transit_chi2': np.nan,
        'flat_chi2': np.nan,
        'delta_chi2': np.nan,
        'transit_bic': np.nan,
        'flat_bic': np.nan,
        'delta_bic': np.nan,
        'transit_parameter_count': 0,
        'flat_parameter_count': 0,
        'flat_baseline': np.nan,
        'flat_a2': np.nan,
        'flat_model_note': None,
        'rprs_sigma': np.nan,
        'duration_ratio': np.nan,
        'eebls_depth_snr': np.nan,
        'residual_scatter': np.nan,
        'use_deviation_from_expected_transit_in_qc': bool(use_deviation_from_expected_transit_in_qc),
        'deviation_sigma_threshold': deviation_sigma_threshold,
        'expected_tmid': expected_context.get('expected_tmid', np.nan),
        'expected_tmid_unc': expected_context.get('expected_tmid_unc', np.nan),
        'expected_tmid_unc_minutes': np.nan,
        'fitted_tmid': np.nan,
        'expected_rprs': expected_context.get('expected_rprs', np.nan),
        'expected_rprs_unc': expected_context.get('expected_rprs_unc', np.nan),
        'tmid_deviation_days': np.nan,
        'tmid_deviation_minutes': np.nan,
        'tmid_deviation_threshold_minutes': np.nan,
        'tmid_deviation_sigma': np.nan,
        'rprs_deviation_sigma': np.nan,
        'tmid_deviation_score': np.nan,
        'rprs_deviation_score': np.nan,
        'deviation_from_expected_value': np.nan,
        'ktmf_metric': np.nan,
        'ktmf_contributions': [],
        'point_count': 0,
    }
    if fit is None:
        return summary

    data = np.asarray(getattr(fit, 'data', np.array([])), dtype=float)
    if data.ndim != 1 or data.size == 0:
        summary['summary'] = (
            "Transit QC unavailable: fit results do not expose the 1D lightcurve data needed for "
            "a transit-vs-flat comparison."
        )
        return summary

    dataerr_obj = getattr(fit, 'dataerr', None)
    dataerr = None if dataerr_obj is None else np.asarray(dataerr_obj, dtype=float)
    if dataerr is not None and dataerr.shape != data.shape:
        dataerr = None

    transit_model_obj = getattr(fit, 'model', None)
    if transit_model_obj is None and hasattr(fit, 'residuals'):
        residuals = np.asarray(getattr(fit, 'residuals'), dtype=float)
        if residuals.shape == data.shape:
            transit_model_obj = data - residuals
    if transit_model_obj is None:
        summary['summary'] = (
            "Transit QC unavailable: fit results do not expose the modeled transit lightcurve needed "
            "for a transit-vs-flat comparison."
        )
        return summary

    transit_model = np.asarray(transit_model_obj, dtype=float)
    if transit_model.shape != data.shape:
        summary['summary'] = (
            "Transit QC unavailable: the fitted transit model shape does not match the lightcurve data."
        )
        return summary

    airmass_obj = getattr(fit, 'airmass', None)
    airmass = None if airmass_obj is None else np.asarray(airmass_obj, dtype=float)
    if airmass is not None and airmass.shape != data.shape:
        airmass = None

    allow_airmass_term = bool(
        airmass is not None
        and airmass.ndim == 1
        and not getattr(fit, 'airmass_fit_skipped', False)
    )
    bounds = getattr(fit, 'bounds', None)
    a2_bounds = bounds.get('a2') if isinstance(bounds, dict) else None
    parameters = getattr(fit, 'parameters', {}) or {}
    errors = getattr(fit, 'errors', {}) or {}
    initial_a2 = parameters.get('a2', 0.0)

    flat_model = fit_profiled_flat_null_model(
        data,
        dataerr,
        airmass,
        initial_a2=initial_a2,
        a2_bounds=a2_bounds,
        allow_airmass_term=allow_airmass_term,
    )
    transit_chi2, point_count = compute_transit_qc_model_chi2(data, transit_model, dataerr)
    transit_parameter_count = infer_transit_qc_parameter_count(fit, flat_model.get('used_airmass_term', False))
    transit_bic = compute_transit_qc_bic(transit_chi2, point_count, transit_parameter_count)
    flat_bic = flat_model.get('bic', np.nan)
    delta_chi2 = flat_model.get('chi2', np.nan) - transit_chi2
    delta_bic = flat_bic - transit_bic

    summary.update({
        'computed': bool(bool(flat_model.get('available')) and np.isfinite(transit_chi2)),
        'transit_chi2': transit_chi2,
        'flat_chi2': flat_model.get('chi2', np.nan),
        'delta_chi2': delta_chi2,
        'transit_bic': transit_bic,
        'flat_bic': flat_bic,
        'delta_bic': delta_bic,
        'transit_parameter_count': int(transit_parameter_count),
        'flat_parameter_count': int(flat_model.get('param_count', 0)),
        'flat_baseline': flat_model.get('baseline', np.nan),
        'flat_a2': flat_model.get('a2', np.nan),
        'flat_model_note': flat_model.get('note'),
        'residual_scatter': transit_qc_residual_scatter(data, transit_model),
        'point_count': int(point_count),
    })

    if not summary['computed']:
        note = flat_model.get('note') or 'flat/null model comparison failed.'
        summary['summary'] = f"Transit QC unavailable: {note}"
        return summary

    if np.isfinite(delta_chi2):
        if delta_chi2 > 1e-12:
            summary['preferred_model'] = 'transit'
        elif delta_chi2 < -1e-12:
            summary['preferred_model'] = 'flat'
        else:
            summary['preferred_model'] = 'ambiguous'

    rprs = parameters.get('rprs', np.nan)
    rprs_err = errors.get('rprs', np.nan)
    if np.isfinite(rprs) and np.isfinite(rprs_err) and rprs_err > 0:
        summary['rprs_sigma'] = float(abs(rprs) / rprs_err)

    duration_expected = getattr(fit, 'duration_expected', np.nan)
    duration_measured = getattr(fit, 'duration_measured', np.nan)
    if (
        np.isfinite(duration_expected)
        and duration_expected > 0
        and np.isfinite(duration_measured)
        and duration_measured >= 0
    ):
        summary['duration_ratio'] = float(duration_measured / duration_expected)

    ensure_lightcurve_fit_eebls_diagnostic(fit)
    summary['eebls_depth_snr'] = extract_lightcurve_fit_eebls_snr(fit)
    deviation_summary = evaluate_transit_qc_expected_value_deviation(
        fit,
        deviation_sigma_threshold,
        enabled=use_deviation_from_expected_transit_in_qc,
    )
    summary.update({
        'expected_tmid': deviation_summary.get('expected_tmid', np.nan),
        'expected_tmid_unc': deviation_summary.get('expected_tmid_unc', np.nan),
        'expected_tmid_unc_minutes': deviation_summary.get('expected_tmid_unc_minutes', np.nan),
        'fitted_tmid': deviation_summary.get('fitted_tmid', np.nan),
        'tmid_deviation_days': deviation_summary.get('tmid_deviation_days', np.nan),
        'tmid_deviation_minutes': deviation_summary.get('tmid_deviation_minutes', np.nan),
        'tmid_deviation_threshold_minutes': deviation_summary.get('tmid_deviation_threshold_minutes', np.nan),
        'tmid_deviation_sigma': deviation_summary.get('tmid_deviation_sigma', np.nan),
        'rprs_deviation_sigma': deviation_summary.get('rprs_deviation_sigma', np.nan),
        'tmid_deviation_score': deviation_summary.get('tmid_deviation_score', np.nan),
        'rprs_deviation_score': deviation_summary.get('rprs_deviation_score', np.nan),
        'deviation_from_expected_value': deviation_summary.get('deviation_from_expected_value', np.nan),
    })

    notes = []
    failure_reasons = []
    status = 'pass'
    comparison_text = (
        f"Delta BIC={delta_bic:.2f}, Delta chi2={delta_chi2:.2f}"
        if np.isfinite(delta_bic) and np.isfinite(delta_chi2)
        else "model comparison unavailable"
    )

    if not np.isfinite(delta_bic) or not np.isfinite(delta_chi2):
        status = 'unknown'
        notes.append("Transit-vs-flat model comparison was not finite.")
    elif delta_chi2 <= 0:
        status = 'fail'
        notes.append("The flat/null model fits the lightcurve at least as well as the transit model.")
        failure_reasons.append("the flat/null model fits the lightcurve at least as well as the transit model")
    elif delta_bic < TRANSIT_QC_DELTA_BIC_FAIL_THRESHOLD:
        status = 'fail'
        notes.append(
            "The transit model does not beat the flat/null model strongly enough to claim a detection."
        )
        failure_reasons.append(
            "the transit model does not beat the flat/null model strongly enough to claim a detection"
        )
    elif delta_bic < TRANSIT_QC_DELTA_BIC_PASS_THRESHOLD:
        status = 'marginal'
        notes.append(
            "The transit model is preferred over the flat/null model, but the evidence is only moderate."
        )
    else:
        notes.append("The transit model is strongly preferred over the flat/null model.")

    if np.isfinite(summary['rprs_sigma']):
        if summary['rprs_sigma'] < TRANSIT_QC_MIN_RPRS_SIGMA:
            status = 'fail'
            notes.append(
                f"The fitted transit depth is only {summary['rprs_sigma']:.2f}-sigma."
            )
            failure_reasons.append(
                f"the fitted transit depth is only {summary['rprs_sigma']:.2f}-sigma"
            )
        elif summary['rprs_sigma'] < TRANSIT_QC_MARGINAL_RPRS_SIGMA and status == 'pass':
            status = 'marginal'
            notes.append(
                f"The fitted transit depth is only {summary['rprs_sigma']:.2f}-sigma."
            )

    if np.isfinite(summary['duration_ratio']):
        if (
            summary['duration_ratio'] < TRANSIT_QC_DURATION_RATIO_MIN
            or summary['duration_ratio'] > TRANSIT_QC_DURATION_RATIO_MAX
        ):
            if status == 'pass':
                status = 'marginal'
            notes.append(
                f"The measured transit duration is {summary['duration_ratio']:.2f}x the modeled duration."
            )

    if np.isfinite(summary['eebls_depth_snr']) and summary['eebls_depth_snr'] < TRANSIT_QC_MIN_EEBLS_SNR:
        if status == 'pass':
            status = 'marginal'
        notes.append(
            f"EEBLS only found a weak box-like event (depth SNR={summary['eebls_depth_snr']:.2f})."
        )

    if use_deviation_from_expected_transit_in_qc:
        notes.extend(deviation_summary.get('notes', []))
        if deviation_summary.get('failed'):
            status = 'fail'
            detailed_reasons = [
                reason for reason in deviation_summary.get('failure_reasons', [])
                if isinstance(reason, str) and reason.strip()
            ]
            if detailed_reasons:
                failure_reasons.extend(detailed_reasons)
            else:
                notes.append(
                    "The fit deviates too far from the expected published Rp/R* value."
                )
                failure_reasons.append(
                    "the fit deviates too far from the expected published Rp/R* value"
                )

    ktmf_metric, ktmf_contributions = compute_transit_qc_ktmf(summary)
    summary['ktmf_metric'] = ktmf_metric
    summary['ktmf_contributions'] = ktmf_contributions

    if status == 'pass':
        summary_text = f"Transit model strongly preferred over flat/null model ({comparison_text})."
    elif status == 'marginal':
        summary_text = f"Transit model preferred over flat/null model, but the detection is marginal ({comparison_text})."
    elif status == 'fail':
        if failure_reasons:
            flat_model_only_failure = all(
                "flat/null model" in reason or "does not beat the flat/null model" in reason
                for reason in failure_reasons
            )
            if flat_model_only_failure:
                summary_text = (
                    f"Transit detection not supported strongly enough against a flat/null model ({comparison_text})."
                )
            else:
                summary_text = (
                    "Transit model is preferred over the flat/null model, but QC rejected the fit because "
                    + "; ".join(failure_reasons)
                    + f" ({comparison_text})."
                )
        else:
            summary_text = f"Transit detection QC rejected this fit ({comparison_text})."
    else:
        summary_text = f"Transit QC unavailable ({comparison_text})."

    summary.update({
        'status': status,
        'summary': summary_text,
        'notes': notes,
    })
    return summary


def annotate_transit_detection_qc(fit, summary=None):
    if fit is None:
        return

    summary = evaluate_transit_detection_qc(fit) if summary is None else dict(summary)
    fit.transit_qc = summary
    fit.transit_qc_computed = bool(summary.get('computed'))
    fit.transit_qc_status = summary.get('status')
    fit.transit_qc_summary = summary.get('summary')
    fit.transit_qc_preferred_model = summary.get('preferred_model')
    fit.transit_qc_delta_bic = summary.get('delta_bic')
    fit.transit_qc_delta_chi2 = summary.get('delta_chi2')
    fit.transit_qc_rprs_sigma = summary.get('rprs_sigma')
    fit.transit_qc_duration_ratio = summary.get('duration_ratio')
    fit.transit_qc_eebls_depth_snr = summary.get('eebls_depth_snr')
    fit.transit_qc_residual_scatter = summary.get('residual_scatter')
    fit.transit_qc_deviation_from_expected_value = summary.get('deviation_from_expected_value')
    fit.transit_qc_deviation_sigma_threshold = summary.get('deviation_sigma_threshold')
    fit.transit_qc_expected_tmid_value = summary.get('expected_tmid')
    fit.transit_qc_expected_tmid_unc = summary.get('expected_tmid_unc')
    fit.transit_qc_expected_tmid_unc_minutes = summary.get('expected_tmid_unc_minutes')
    fit.transit_qc_fitted_tmid = summary.get('fitted_tmid')
    fit.transit_qc_tmid_deviation_days = summary.get('tmid_deviation_days')
    fit.transit_qc_tmid_deviation_minutes = summary.get('tmid_deviation_minutes')
    fit.transit_qc_tmid_deviation_threshold_minutes = summary.get('tmid_deviation_threshold_minutes')
    fit.transit_qc_tmid_deviation_sigma = summary.get('tmid_deviation_sigma')
    fit.transit_qc_rprs_deviation_sigma = summary.get('rprs_deviation_sigma')
    fit.transit_qc_expected_rprs_deviation_sigma = summary.get('rprs_deviation_sigma')
    fit.transit_qc_ktmf_metric = summary.get('ktmf_metric')
    fit.transit_qc_ktmf_contributions = summary.get('ktmf_contributions')


def lightcurve_fit_transit_qc_failure_reason(fit):
    if fit is None:
        return None

    transit_qc = getattr(fit, 'transit_qc', None)
    if not isinstance(transit_qc, dict):
        return None

    status = str(transit_qc.get('status', '')).strip().lower()
    if status != 'fail':
        return None

    summary = transit_qc.get('summary')
    if isinstance(summary, str) and summary.strip():
        return summary.strip()

    return "Transit detection QC flagged this fit as a poor transit candidate."


def lightcurve_fit_transit_qc_passed(fit):
    if fit is None:
        return False

    transit_qc = getattr(fit, 'transit_qc', None)
    if isinstance(transit_qc, dict):
        return str(transit_qc.get('status', '')).strip().lower() == 'pass'

    return str(getattr(fit, 'transit_qc_status', '')).strip().lower() == 'pass'


def make_json_safe(value):
    if isinstance(value, dict):
        return {str(key): make_json_safe(subvalue) for key, subvalue in value.items()}
    if isinstance(value, (list, tuple)):
        return [make_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [make_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def archive_exception_payload(action, exc):
    traceback_text = ''.join(
        traceback.format_exception(type(exc), exc, exc.__traceback__)
    ).strip()
    message = f"{action}: {type(exc).__name__}: {exc}"
    try:
        log_info(f"Warning: {message}", warn=True)
    except Exception:
        print(f"Warning: {message}", flush=True)
    log.debug("%s\n%s", message, traceback_text)
    return {
        'action': action,
        'error_type': type(exc).__name__,
        'message': str(exc),
        'traceback': traceback_text,
    }


def failed_comparison_archive_dir(save_dir, comp_index):
    base_dir = Path(save_dir)
    candidate_dir = base_dir / f"comp_{comp_index + 1}_failed"
    if not candidate_dir.exists():
        return candidate_dir

    suffix = 2
    while True:
        fallback_dir = base_dir / f"comp_{comp_index + 1}_failed_{suffix}"
        if not fallback_dir.exists():
            return fallback_dir
        suffix += 1


def comparison_candidate_output_dir(save_dir, comp_index):
    return Path(save_dir) / f"comp{comp_index + 1}"


def triangle_plot_output_path(save_dir, planet_name, observation_date):
    return (
        Path(save_dir)
        / "temp"
        / safe_output_filename("Triangle", planet_name, filename_date_token(observation_date), extension="png")
    )


def final_triangle_plot_output_path(save_dir, planet_name, observation_date):
    return (
        Path(save_dir)
        / safe_output_filename("FinalTriangle", planet_name, filename_date_token(observation_date), extension="png")
    )


def zoomed_final_triangle_plot_output_path(save_dir, planet_name, observation_date):
    return (
        Path(save_dir)
        / safe_output_filename("ZoomedTrianglePlot", planet_name, filename_date_token(observation_date), extension="png")
    )


def comparison_candidate_triangle_plot_output_path(save_dir, planet_name, observation_date, comp_index):
    return (
        Path(save_dir)
        / "temp"
        / safe_output_filename(
            f"Comp{int(comp_index) + 1}_Triangle",
            planet_name,
            filename_date_token(observation_date),
            extension="png",
        )
    )


def comparison_candidate_label_from_output_dir(output_dir):
    if output_dir is None:
        return None

    for part in reversed(Path(output_dir).parts):
        match = re.fullmatch(r"comp(\d+)", str(part), re.IGNORECASE)
        if match:
            return f"comparison candidate #{int(match.group(1))}"
    return None


def _plot_triangle_for_output(fit, plot_title=None, required_keywords=(), **plot_kwargs):
    plotter = getattr(fit, 'plot_triangle', None)
    if not callable(plotter):
        return None

    supported_kwargs = {}
    for keyword in required_keywords:
        if not callable_accepts_keyword(plotter, keyword):
            return None

    if plot_title and callable_accepts_keyword(plotter, 'plot_title'):
        supported_kwargs['plot_title'] = plot_title
    for keyword, value in plot_kwargs.items():
        if callable_accepts_keyword(plotter, keyword):
            supported_kwargs[keyword] = value

    return plotter(**supported_kwargs)


def _close_plot_figure(fig):
    try:
        plt.close(fig)
    except TypeError:
        pass


def save_final_triangle_plot(fit, save_dir, planet_name, observation_date, source_dir=None):
    output_path = final_triangle_plot_output_path(save_dir, planet_name, observation_date)
    zoomed_output_path = zoomed_final_triangle_plot_output_path(save_dir, planet_name, observation_date)
    compatibility_path = triangle_plot_output_path(save_dir, planet_name, observation_date)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    zoomed_output_path.parent.mkdir(parents=True, exist_ok=True)
    compatibility_path.parent.mkdir(parents=True, exist_ok=True)

    source_label = comparison_candidate_label_from_output_dir(source_dir)
    plot_title = "Final selected fit"
    if source_label:
        plot_title = f"{plot_title} ({source_label})"
    fig = _plot_triangle_for_output(fit, plot_title=plot_title)
    if fig is None:
        return None
    fig.savefig(output_path)
    if compatibility_path != output_path:
        try:
            shutil.copy2(output_path, compatibility_path)
        except Exception:
            fig.savefig(compatibility_path)
    _close_plot_figure(fig)

    zoomed_fig = None
    try:
        zoomed_title = f"{plot_title} (5-sigma zoom)"
        zoomed_fig = _plot_triangle_for_output(
            fit,
            plot_title=zoomed_title,
            required_keywords=('zoom_sigma',),
            zoom_sigma=5.0,
        )
        if zoomed_fig is not None:
            zoomed_fig.savefig(zoomed_output_path)
    except Exception as exc:
        try:
            log_info(f"Warning: Could not save zoomed final triangle plot: {exc}", warn=True)
        except Exception:
            pass
    finally:
        if zoomed_fig is not None:
            _close_plot_figure(zoomed_fig)
    return output_path


def estimate_transit_duration_samples_from_fit(fit, sample_count=1000, grid_size=1000):
    if fit is None or not hasattr(fit, 'parameters') or not hasattr(fit, 'errors'):
        return None, np.array([], dtype=float)

    fit_times = np.asarray(getattr(fit, 'time', []), dtype=float)
    fit_times = fit_times[np.isfinite(fit_times)]
    if fit_times.size < 2:
        return None, np.array([], dtype=float)

    parameters = getattr(fit, 'parameters', {}) or {}
    errors = getattr(fit, 'errors', {}) or {}
    transit_times = np.linspace(np.nanmin(fit_times), np.nanmax(fit_times), int(grid_size))
    if transit_times.size < 2:
        return None, np.array([], dtype=float)

    baseline_parameters = dict(parameters)
    baseline_model = transit(transit_times, baseline_parameters)
    dt = float(np.nanmean(np.diff(transit_times)))
    if not np.isfinite(dt) or dt <= 0:
        return baseline_model, np.array([], dtype=float)

    duration_samples = []
    sample_count = max(1, int(sample_count))
    for _ in range(sample_count):
        sampled_parameters = dict(parameters)
        for key, error_value in errors.items():
            parameter_value = parameters.get(key)
            if parameter_value is None:
                continue
            try:
                numeric_error = float(error_value)
                numeric_value = float(parameter_value)
            except (TypeError, ValueError):
                continue
            if not np.isfinite(numeric_error) or numeric_error <= 0 or not np.isfinite(numeric_value):
                continue
            sampled_parameters[key] = np.random.normal(numeric_value, numeric_error)

        sampled_model = transit(transit_times, sampled_parameters)
        in_transit_mask = np.asarray(sampled_model, dtype=float) < 1
        duration_samples.append(float(np.count_nonzero(in_transit_mask)) * dt)

    return baseline_model, np.asarray(duration_samples, dtype=float)


def build_comparison_candidate_adaptive_summary(comparison_calibration, psf_data,
                                                use_adaptive_apertures=False,
                                                adaptive_aperture_values=None,
                                                adaptive_annulus_values=None,
                                                fallback_sigma=np.nan):
    if (
        not use_adaptive_apertures
        or comparison_calibration is None
        or comparison_calibration.get('method') == 'psf'
        or adaptive_aperture_values is None
        or adaptive_annulus_values is None
    ):
        return None

    aperture_index = comparison_calibration.get('a')
    annulus_index = comparison_calibration.get('an')
    if aperture_index is None or annulus_index is None:
        return None

    aperture_grid = np.asarray(adaptive_aperture_values, dtype=float)
    annulus_grid = np.asarray(adaptive_annulus_values, dtype=float)
    if aperture_index >= aperture_grid.size or annulus_index >= annulus_grid.size:
        return None

    return summarize_adaptive_aperture_usage(
        psf_data['target'],
        aperture_grid[aperture_index],
        annulus_grid[annulus_index],
        fallback_sigma=fallback_sigma,
    )


def build_comparison_candidate_transit_prior(p_dict, ld):
    return {
        'rprs': p_dict['rprs'],
        'ars': p_dict['aRs'],
        'per': p_dict['pPer'],
        'inc': p_dict['inc'],
        'u0': ld[0], 'u1': ld[1], 'u2': ld[2], 'u3': ld[3],
        'ecc': p_dict['ecc'],
        'omega': p_dict['omega'],
        'tmid': p_dict['midT'],
        'a2': 0,
    }


def prepare_comparison_candidate_full_reduction_series(times, target_flux, comp_flux, airmass,
                                                       jd_times=None, adaptive_summary=None):
    result = {
        'applied': False,
        'failure_reason': "the raw comparison-candidate photometry did not yield a usable light curve.",
        'filter_diagnostics': [],
        'debug_times': np.array([], dtype=float),
        'debug_target_flux': np.array([], dtype=float),
        'debug_comp_flux': np.array([], dtype=float),
        'debug_raw_ratio': np.array([], dtype=float),
        'initial_sigma_keep_mask': np.array([], dtype=bool),
        'time': np.array([], dtype=float),
        'flux': np.array([], dtype=float),
        'unc': np.array([], dtype=float),
        'airmass': np.array([], dtype=float),
        'jd_time': np.array([], dtype=float),
        'target_flux': np.array([], dtype=float),
        'comp_flux': np.array([], dtype=float),
        'source_indices': np.array([], dtype=int),
    }

    prepared = prepare_lightcurve_fit_input_series(
        times,
        target_flux,
        comp_flux,
        airmass,
        jd_times=jd_times,
    )
    result['filter_diagnostics'] = prepared.get('filter_diagnostics', [])
    for key in ('debug_times', 'debug_target_flux', 'debug_comp_flux', 'debug_raw_ratio', 'initial_sigma_keep_mask'):
        if key in prepared:
            result[key] = prepared[key]
    if not prepared.get('applied'):
        result['failure_reason'] = prepared.get(
            'failure_reason',
            "the raw comparison-candidate photometry did not yield a usable light curve.",
        )
        return result

    good_times = np.asarray(prepared['time'], dtype=float)
    good_flux = np.asarray(prepared['flux'], dtype=float)
    good_unc = np.asarray(prepared['unc'], dtype=float)
    good_airmass = np.asarray(prepared['airmass'], dtype=float)
    good_jd_times = np.asarray(prepared['jd_time'], dtype=float)
    good_target_flux = np.asarray(prepared['target_flux'], dtype=float)
    good_comp_flux = np.asarray(prepared['comp_flux'], dtype=float)
    source_indices = np.asarray(prepared['source_indices'], dtype=int)

    adaptive_clip_mask = np.zeros(good_times.shape[0], dtype=bool)
    if adaptive_summary is not None:
        aperture_series = np.asarray(adaptive_summary.get('aperture_series', []), dtype=float)
        annulus_series = np.asarray(adaptive_summary.get('annulus_series', []), dtype=float)
        if aperture_series.ndim == 1 and annulus_series.ndim == 1:
            try:
                selected_apertures = aperture_series[source_indices]
                selected_annuli = annulus_series[source_indices]
            except IndexError:
                selected_apertures = None
                selected_annuli = None
            if (
                selected_apertures is not None
                and selected_apertures.shape == good_times.shape
                and selected_annuli.shape == good_times.shape
            ):
                adaptive_clip_mask = adaptive_aperture_outlier_mask(
                    selected_apertures,
                    selected_annuli,
                )

    if np.count_nonzero(~adaptive_clip_mask) < LIGHTCURVE_MIN_VALID_POINTS:
        result['failure_reason'] = (
            "adaptive-aperture filtering left too few points for a stable comparison-candidate reduction."
        )
        return result

    if np.any(adaptive_clip_mask):
        good_times = good_times[~adaptive_clip_mask]
        good_flux = good_flux[~adaptive_clip_mask]
        good_unc = good_unc[~adaptive_clip_mask]
        good_airmass = good_airmass[~adaptive_clip_mask]
        good_jd_times = good_jd_times[~adaptive_clip_mask]
        good_target_flux = good_target_flux[~adaptive_clip_mask]
        good_comp_flux = good_comp_flux[~adaptive_clip_mask]
        source_indices = source_indices[~adaptive_clip_mask]

    relative_flux_mask = relative_flux_filter_mask(good_flux)
    if np.count_nonzero(relative_flux_mask) < LIGHTCURVE_MIN_VALID_POINTS:
        result['failure_reason'] = (
            "the raw comparison-candidate light curve failed the relative-flux filter before full reduction."
        )
        return result

    result.update({
        'applied': True,
        'failure_reason': None,
        'time': good_times[relative_flux_mask],
        'flux': good_flux[relative_flux_mask],
        'unc': good_unc[relative_flux_mask],
        'airmass': good_airmass[relative_flux_mask],
        'jd_time': good_jd_times[relative_flux_mask],
        'target_flux': good_target_flux[relative_flux_mask],
        'comp_flux': good_comp_flux[relative_flux_mask],
        'source_indices': source_indices[relative_flux_mask],
    })
    return result


def comparison_candidate_coverage_priority(assessment):
    if not isinstance(assessment, dict) or not assessment.get('valid'):
        return 4

    pre_points = int(assessment.get('pre_ingress_points', 0) or 0)
    post_points = int(assessment.get('post_egress_points', 0) or 0)
    transit_fraction = coerce_finite_transit_qc_scalar(
        assessment.get('transit_fraction_observed', np.nan)
    )
    has_two_sided_oot = pre_points > 0 and post_points > 0
    covers_full_window = (
        bool(assessment.get('covers_ingress', False))
        and bool(assessment.get('covers_mid_transit', False))
        and bool(assessment.get('covers_egress', False))
    )

    if has_two_sided_oot and covers_full_window:
        return 0
    if has_two_sided_oot:
        return 1
    if np.isfinite(transit_fraction) and transit_fraction >= 0.75 and assessment.get('covers_mid_transit', False):
        return 2
    if np.isfinite(transit_fraction) and transit_fraction > 0:
        return 3
    return 4


def is_low_one_sided_expected_transit_coverage(assessment):
    if not isinstance(assessment, dict) or not assessment.get('valid'):
        return False

    pre_points = int(assessment.get('pre_ingress_points', 0) or 0)
    post_points = int(assessment.get('post_egress_points', 0) or 0)
    success_label = str(assessment.get('success_label', '')).strip().lower()
    expected_successful = bool(assessment.get('expected_successful', False))
    return (pre_points == 0 or post_points == 0) and (
        success_label in ('very low', 'low') or not expected_successful
    )


def partial_transit_geometry_retry_limits(assessment):
    active = is_low_one_sided_expected_transit_coverage(assessment)
    note = None
    if active:
        note = (
            "Skipped; pre-UltraNest coverage is one-sided/LOW, so EXOTIC does not expand this "
            "geometry posterior range while the transit shape is baseline-degenerate."
        )
    return {
        'active': active,
        'note': note,
        'max_retries': {
            'rprs': PARTIAL_COVERAGE_RPRS_POSTERIOR_MAX_RETRIES,
            'ars': PARTIAL_COVERAGE_ARS_POSTERIOR_MAX_RETRIES,
            'b': PARTIAL_COVERAGE_IMPACT_PARAMETER_POSTERIOR_MAX_RETRIES,
        },
    }


def score_comparison_candidate_lightcurve_scout(prepared_series, eebls_summary, prior):
    result = {
        'score': np.nan,
        'scatter': np.nan,
        'scatter_score': np.nan,
        'depth_score': np.nan,
        'eebls_score': np.nan,
        'eebls_snr': np.nan,
        'eebls_depth': np.nan,
        'expected_depth': np.nan,
    }
    if not isinstance(prepared_series, dict) or not prepared_series.get('applied'):
        return result

    flux = np.asarray(prepared_series.get('flux', []), dtype=float)
    finite_flux = flux[np.isfinite(flux) & (flux > 0)]
    if finite_flux.size < LIGHTCURVE_MIN_VALID_POINTS:
        return result

    baseline = bn.nanmedian(finite_flux)
    if not np.isfinite(baseline) or baseline <= 0:
        return result

    normalized_flux = finite_flux / baseline
    scatter = robust_scatter(normalized_flux - bn.nanmedian(normalized_flux))
    rprs = coerce_finite_transit_qc_scalar(prior.get('rprs', np.nan) if isinstance(prior, dict) else np.nan)
    expected_depth = rprs ** 2 if np.isfinite(rprs) and rprs >= 0 else np.nan
    if np.isfinite(scatter):
        result['scatter'] = float(scatter)
    if np.isfinite(expected_depth):
        result['expected_depth'] = float(expected_depth)

    depth = coerce_finite_transit_qc_scalar((eebls_summary or {}).get('depth', np.nan))
    depth_snr = coerce_finite_transit_qc_scalar((eebls_summary or {}).get('depth_snr', np.nan))
    if np.isfinite(depth):
        result['eebls_depth'] = float(depth)
    if np.isfinite(depth_snr):
        result['eebls_snr'] = float(depth_snr)

    scatter_reference = expected_depth if np.isfinite(expected_depth) and expected_depth > 0 else 0.005
    if np.isfinite(scatter) and scatter >= 0:
        result['scatter_score'] = float(np.clip(1.0 / (1.0 + scatter / max(scatter_reference, 1e-6)), 0.0, 1.0))

    if np.isfinite(depth) and depth > 0 and np.isfinite(expected_depth) and expected_depth > 0:
        depth_ratio = depth / expected_depth
        if np.isfinite(depth_ratio) and depth_ratio > 0:
            result['depth_score'] = float(np.clip(np.exp(-abs(np.log(depth_ratio)) / np.log(2.0)), 0.0, 1.0))

    if np.isfinite(depth_snr) and depth_snr > 0:
        result['eebls_score'] = float(np.clip(depth_snr / 8.0, 0.0, 1.0))

    components = [
        (0.45, result['scatter_score']),
        (0.35, result['depth_score']),
        (0.20, result['eebls_score']),
    ]
    available = [(weight, value) for weight, value in components if np.isfinite(value)]
    if available:
        weight_sum = sum(weight for weight, _ in available)
        result['score'] = float(sum(weight * value for weight, value in available) / weight_sum)
    return result


def build_comparison_candidate_preflight(times, jd_times, airmass, ld, p_dict, target_flux, comp_flux,
                                         adaptive_summary=None, use_eebls_to_initialize_tmid_and_bounds=True):
    prepared = prepare_comparison_candidate_full_reduction_series(
        times,
        target_flux,
        comp_flux,
        airmass,
        jd_times=jd_times,
        adaptive_summary=adaptive_summary,
    )
    try:
        prior = build_comparison_candidate_transit_prior(p_dict, ld)
    except (KeyError, IndexError, TypeError, ValueError):
        return {
            'prepared_series': prepared,
            'coverage_assessment': None,
            'coverage_priority': 4,
            'eebls_summary': None,
            'tmid_search_summary': None,
            'duration_prior': None,
            'scout': {'score': np.nan},
        }

    duration_prior = build_single_transit_duration_prior(p_dict)
    eebls_summary = None
    tmid_search_summary = None
    coverage_assessment = None
    scout = {'score': np.nan}

    if prepared.get('applied'):
        good_times = np.asarray(prepared['time'], dtype=float)
        good_flux = np.asarray(prepared['flux'], dtype=float)
        good_unc = np.asarray(prepared['unc'], dtype=float)
        expected_duration = estimate_transit_duration_from_prior_geometry(prior)
        tmid_search_summary = estimate_ephemeris_tmid_and_bounds(
            good_times,
            p_dict.get('midT', prior.get('tmid', np.nan)),
            prior['per'],
            p_dict.get('midTUnc', 0.01),
            p_dict.get('pPerUnc', 0.0),
            expected_duration=expected_duration,
            sigma_multiplier=35.0,
        )
        prior['tmid'] = tmid_search_summary['tmid']
        lower, upper = tmid_search_summary['bounds']
        if use_eebls_to_initialize_tmid_and_bounds:
            eebls_summary = estimate_tmid_and_bounds_with_eebls(
                good_times,
                good_flux,
                good_unc,
                prior,
                [lower, upper],
            )
        else:
            eebls_summary = {'applied': False, 'depth': np.nan, 'depth_snr': np.nan}
        coverage_assessment = build_expected_transit_coverage_assessment(
            good_times,
            prior,
            flux_values=good_flux,
            flux_errors=good_unc,
            tmid_search_summary=tmid_search_summary,
            duration_prior=duration_prior,
        )
        scout = score_comparison_candidate_lightcurve_scout(prepared, eebls_summary, prior)

    return {
        'prepared_series': prepared,
        'coverage_assessment': coverage_assessment,
        'coverage_priority': comparison_candidate_coverage_priority(coverage_assessment),
        'eebls_summary': eebls_summary,
        'tmid_search_summary': tmid_search_summary,
        'duration_prior': duration_prior,
        'scout': scout,
    }


def comparison_preflight_field_band_limit(plans):
    finite_scores = [
        plan['summary'].get('aggregate_score', np.nan)
        for plan in plans
        if np.isfinite(plan['summary'].get('aggregate_score', np.nan))
    ]
    if not finite_scores:
        return np.inf
    best_score = float(min(finite_scores))
    return best_score + max(
        COMPARISON_PREFLIGHT_FIELD_SCORE_ABSOLUTE_BAND,
        abs(best_score) * COMPARISON_PREFLIGHT_FIELD_SCORE_RELATIVE_BAND,
    )


def rank_comparison_candidate_preflight_plans(plans):
    if not plans:
        return []

    field_band_limit = comparison_preflight_field_band_limit(plans)

    def sort_key(plan):
        preflight = plan.get('preflight') or {}
        scout = preflight.get('scout') or {}
        aggregate_score = plan['summary'].get('aggregate_score', np.inf)
        finite_aggregate = aggregate_score if np.isfinite(aggregate_score) else np.inf
        close_field_band = 0 if finite_aggregate <= field_band_limit else 1
        scout_score = scout.get('score', np.nan)
        scout_sort = -float(scout_score) if np.isfinite(scout_score) else np.inf
        return (
            int(preflight.get('coverage_priority', 4)),
            close_field_band,
            scout_sort,
            finite_aggregate,
            plan.get('field_rank', np.inf),
        )

    return sorted(plans, key=sort_key)


def log_comparison_candidate_preflight_order(plans, ranked_plans):
    if not plans or not ranked_plans:
        return
    original_order = [plan['summary'].get('comp_index') for plan in plans]
    ranked_order = [plan['summary'].get('comp_index') for plan in ranked_plans]
    if original_order == ranked_order:
        return

    log_info(
        "Comparison-star target-fit order adjusted by pre-UltraNest coverage/scout preflight "
        "(Tmid remains free; scout uses coverage, scatter, depth plausibility, and EEBLS SNR)."
    )
    for new_rank, plan in enumerate(ranked_plans, start=1):
        summary = plan['summary']
        preflight = plan.get('preflight') or {}
        coverage = preflight.get('coverage_assessment') or {}
        scout = preflight.get('scout') or {}
        scout_score = scout.get('score', np.nan)
        scout_text = "n/a" if not np.isfinite(scout_score) else f"{scout_score:.3f}"
        scatter = scout.get('scatter', np.nan)
        scatter_text = "n/a" if not np.isfinite(scatter) else f"{100.0 * scatter:.4f}%"
        eebls_snr = scout.get('eebls_snr', np.nan)
        eebls_text = "n/a" if not np.isfinite(eebls_snr) else f"{eebls_snr:.2f}"
        label = summary.get('label', f"Comp {summary.get('comp_index', 0) + 1}")
        log_info(
            f"  Preflight rank {new_rank}: {label} "
            f"(field rank {plan.get('field_rank', 0) + 1}), coverage_priority={preflight.get('coverage_priority', 'n/a')}, "
            f"pre/post={coverage.get('pre_ingress_points', 'n/a')}/{coverage.get('post_egress_points', 'n/a')}, "
            f"scout={scout_text}, scatter={scatter_text}, eebls_snr={eebls_text}."
        )


def match_time_subset_indices(full_times, subset_times, rtol=1e-10, atol=1e-10):
    full_times = np.asarray(full_times, dtype=float).reshape(-1)
    subset_times = np.asarray(subset_times, dtype=float).reshape(-1)
    if subset_times.size == 0:
        return np.array([], dtype=int)
    if full_times.size < subset_times.size:
        return None

    matched_indices = []
    search_start = 0
    for subset_time in subset_times:
        if not np.isfinite(subset_time):
            return None
        remaining = full_times[search_start:]
        matches = np.flatnonzero(np.isclose(remaining, subset_time, rtol=rtol, atol=atol))
        if matches.size == 0:
            return None
        matched_index = search_start + int(matches[0])
        matched_indices.append(matched_index)
        search_start = matched_index + 1

    return np.asarray(matched_indices, dtype=int)


def finalize_comparison_candidate_full_reduction(times, target_flux, comp_flux, airmass, ld, p_dict,
                                                 jd_times=None,
                                                 disable_vertical_flux_normalization=False,
                                                 detrend_on_outoftransit_baseline=True,
                                                 use_impactparameter_rather_than_inclination_to_fit=True,
                                                 use_eebls_to_initialize_tmid_and_bounds=True,
                                                 plot_time_range=None,
                                                 baseline_duration_multiplier=FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT,
                                                 adaptive_summary=None,
                                                 run_fast_ultranest_before_final_run=FAST_ULTRANEST_BEFORE_FINAL_RUN_DEFAULT,
                                                 precomputed_candidate_series=None):
    result = {
        'applied': False,
        'fit': None,
        'good_times': np.array([], dtype=float),
        'good_flux': np.array([], dtype=float),
        'good_unc': np.array([], dtype=float),
        'good_airmass': np.array([], dtype=float),
        'good_jd_times': np.array([], dtype=float),
        'good_target_flux': np.array([], dtype=float),
        'good_comp_flux': np.array([], dtype=float),
        'source_indices': np.array([], dtype=int),
        'data_highres': None,
        'duration_samples': np.array([], dtype=float),
        'failure_reason': "full candidate reduction did not run.",
        'filter_diagnostics': [],
        'note': None,
    }
    if precomputed_candidate_series is None:
        prepared = prepare_comparison_candidate_full_reduction_series(
            times,
            target_flux,
            comp_flux,
            airmass,
            jd_times=jd_times,
            adaptive_summary=adaptive_summary,
        )
    else:
        prepared = precomputed_candidate_series
    result['filter_diagnostics'] = prepared.get('filter_diagnostics', [])
    if not prepared.get('applied'):
        result['failure_reason'] = prepared.get(
            'failure_reason',
            "the raw comparison-candidate photometry did not yield a usable light curve.",
        )
        return result

    good_times = np.asarray(prepared['time'], dtype=float)
    good_flux = np.asarray(prepared['flux'], dtype=float)
    good_unc = np.asarray(prepared['unc'], dtype=float)
    good_airmass = np.asarray(prepared['airmass'], dtype=float)
    good_jd_times = np.asarray(prepared['jd_time'], dtype=float)
    good_target_flux = np.asarray(prepared['target_flux'], dtype=float)
    good_comp_flux = np.asarray(prepared['comp_flux'], dtype=float)
    source_indices = np.asarray(prepared['source_indices'], dtype=int)

    prior = build_comparison_candidate_transit_prior(p_dict, ld)

    expected_duration = estimate_transit_duration_from_prior_geometry(prior)
    tmid_search_summary = estimate_ephemeris_tmid_and_bounds(
        good_times,
        p_dict['midT'],
        prior['per'],
        p_dict['midTUnc'],
        p_dict['pPerUnc'],
        expected_duration=expected_duration,
        sigma_multiplier=35.0,
    )
    prior['tmid'] = tmid_search_summary['tmid']
    lower, upper = tmid_search_summary['bounds']

    eebls_search_summary = None
    if use_eebls_to_initialize_tmid_and_bounds:
        eebls_search_summary = estimate_tmid_and_bounds_with_eebls(
            good_times,
            good_flux,
            good_unc,
            prior,
            [lower, upper],
        )
        if eebls_search_summary.get('applied'):
            prior['tmid'] = eebls_search_summary['tmid']
            lower, upper = eebls_search_summary['bounds']

    skip_final_airmass_fit = False
    airmass_skip_note = None
    final_airmass_span = airmass_span(good_airmass)
    if should_skip_airmass_fit(good_airmass):
        skip_final_airmass_fit = True
        airmass_skip_note = (
            f"Skipped (airmass span {final_airmass_span:.4f} <= {AIRMASS_FLAT_RANGE_THRESHOLD:.2f}); "
            "no airmass correction applied."
        )

    bounds = build_initial_transit_bounds(
        prior,
        [lower, upper],
        ars_unc=p_dict.get('aRsUnc'),
    )
    apply_vertical_flux_normalization_bound(
        prior,
        bounds,
        good_flux,
        disable_vertical_flux_normalization,
    )
    if not skip_final_airmass_fit:
        bounds['a2'] = [-3, 3]
    ensure_pre_final_ultranest_baseline_bounds(prior, bounds, good_flux, fit_a2=True)

    debug_phase_clip_keep_mask = None
    prefit = lc_fitter(
        good_times,
        good_flux,
        good_unc,
        good_airmass,
        prior,
        bounds,
        jd_times=good_jd_times,
        mode='lm',
        use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
    )
    if (
        prefit is not None
        and hasattr(prefit, 'residuals')
        and hasattr(prefit, 'phase')
        and np.shape(prefit.residuals) == np.shape(good_times)
        and np.shape(prefit.phase) == np.shape(good_times)
    ):
        phase_clip_mask = phase_bin_sigma_clip(prefit.residuals, prefit.phase, sigma=3, bins=10)
        min_required_points = max(len(bounds) + 1, LIGHTCURVE_MIN_VALID_POINTS)
        if np.any(phase_clip_mask) and np.count_nonzero(~phase_clip_mask) >= min_required_points:
            debug_phase_clip_keep_mask = np.asarray(~phase_clip_mask, dtype=bool).copy()
            result['filter_diagnostics'].append(build_time_rejection_diagnostic(
                "Final-fit phase residual clip",
                good_times,
                ~phase_clip_mask,
                note="Dropped phase-binned residual outliers before the comparison-candidate ultranest fit.",
            ))
            good_times = good_times[~phase_clip_mask]
            good_flux = good_flux[~phase_clip_mask]
            good_unc = good_unc[~phase_clip_mask]
            good_airmass = good_airmass[~phase_clip_mask]
            good_jd_times = good_jd_times[~phase_clip_mask]
            good_target_flux = good_target_flux[~phase_clip_mask]
            good_comp_flux = good_comp_flux[~phase_clip_mask]
            source_indices = source_indices[~phase_clip_mask]

    full_good_times = np.asarray(good_times, dtype=float)
    full_good_flux = np.asarray(good_flux, dtype=float)
    full_good_unc = np.asarray(good_unc, dtype=float)
    full_good_airmass = np.asarray(good_airmass, dtype=float)
    full_good_jd_times = np.asarray(good_jd_times, dtype=float)
    full_good_target_flux = np.asarray(good_target_flux, dtype=float)
    full_good_comp_flux = np.asarray(good_comp_flux, dtype=float)
    full_source_indices = np.asarray(source_indices, dtype=int)

    fast_binning = {'applied': False, 'note': None}
    fit_times = full_good_times
    fit_flux = full_good_flux
    fit_unc = full_good_unc
    fit_airmass = full_good_airmass
    fit_jd_times = full_good_jd_times
    if run_fast_ultranest_before_final_run:
        fast_binning = build_fast_ultranest_lightcurve_series(
            full_good_times,
            full_good_flux,
            full_good_unc,
            full_good_airmass,
            jd_times=full_good_jd_times,
        )
        if fast_binning.get('applied'):
            log_info(fast_binning['note'])
            fit_times = fast_binning['time']
            fit_flux = fast_binning['flux']
            fit_unc = fast_binning['unc']
            fit_airmass = fast_binning['airmass']
            fit_jd_times = fast_binning['jd_times']

    fit_prior = dict(prior)
    fit_bounds = clone_lightcurve_bounds(bounds)
    ensure_pre_final_ultranest_baseline_bounds(fit_prior, fit_bounds, fit_flux, fit_a2=True)
    pre_ultranest_coverage_assessment = build_expected_transit_coverage_assessment(
        full_good_times,
        prior,
        flux_values=full_good_flux,
        flux_errors=full_good_unc,
        tmid_search_summary=tmid_search_summary,
        duration_prior=build_single_transit_duration_prior(p_dict),
    )

    final_fit, fitted_flux, fitted_unc = fit_final_lightcurve_with_oot_baseline_detrending(
        fit_times,
        fit_flux,
        fit_unc,
        fit_airmass,
        fit_prior,
        fit_bounds,
        jd_times=fit_jd_times,
        skip_airmass_fit=skip_final_airmass_fit,
        airmass_skip_note=airmass_skip_note,
        disable_vertical_flux_normalization=disable_vertical_flux_normalization,
        detrend_on_outoftransit_baseline=detrend_on_outoftransit_baseline,
        use_impactparameter_rather_than_inclination_to_fit=
        use_impactparameter_rather_than_inclination_to_fit,
        plot_time_range=plot_time_range,
        baseline_duration_multiplier=baseline_duration_multiplier,
        expected_planet_dict=p_dict,
        expected_tmid_search_summary=tmid_search_summary,
        eebls_search_summary=eebls_search_summary,
        extend_sparse_posterior_live_points=False,
        keep_ultranest_sampler_for_deferred_extension=not bool(fast_binning.get('applied')),
        fix_baseline_terms_for_final=not bool(fast_binning.get('applied')),
        pre_ultranest_coverage_assessment=pre_ultranest_coverage_assessment,
    )
    annotate_fast_ultranest_binning(final_fit, fast_binning)
    if final_fit is None:
        result['failure_reason'] = "the full comparison-candidate reduction did not converge."
        return result

    final_fit_times = np.asarray(getattr(final_fit, 'time', good_times), dtype=float)
    if not fast_binning.get('applied') and (
        final_fit_times.shape != good_times.shape
        or not np.allclose(final_fit_times, good_times, rtol=1e-10, atol=1e-10)
    ):
        final_time_indices = match_time_subset_indices(good_times, final_fit_times)
        if final_time_indices is not None:
            good_times = good_times[final_time_indices]
            good_airmass = good_airmass[final_time_indices]
            good_jd_times = good_jd_times[final_time_indices]
            good_target_flux = good_target_flux[final_time_indices]
            good_comp_flux = good_comp_flux[final_time_indices]
            source_indices = source_indices[final_time_indices]

    annotate_lightcurve_filter_diagnostics(final_fit, result['filter_diagnostics'])
    annotate_selected_photometry_debug(
        final_fit,
        prepared['debug_times'],
        prepared['debug_target_flux'],
        prepared['debug_comp_flux'],
        prepared['debug_raw_ratio'],
        prepared['initial_sigma_keep_mask'],
        phase_clip_keep_mask_on_sigma_filtered=debug_phase_clip_keep_mask,
    )

    data_highres, duration_samples = estimate_transit_duration_samples_from_fit(final_fit)
    result.update({
        'applied': True,
        'fit': final_fit,
        'good_times': full_good_times if fast_binning.get('applied') else np.asarray(good_times, dtype=float),
        'good_flux': full_good_flux if fast_binning.get('applied') else np.asarray(good_flux, dtype=float),
        'good_unc': full_good_unc if fast_binning.get('applied') else np.asarray(good_unc, dtype=float),
        'good_airmass': full_good_airmass if fast_binning.get('applied') else np.asarray(good_airmass, dtype=float),
        'good_jd_times': full_good_jd_times if fast_binning.get('applied') else np.asarray(good_jd_times, dtype=float),
        'good_target_flux': full_good_target_flux if fast_binning.get('applied') else np.asarray(good_target_flux, dtype=float),
        'good_comp_flux': full_good_comp_flux if fast_binning.get('applied') else np.asarray(good_comp_flux, dtype=float),
        'source_indices': full_source_indices if fast_binning.get('applied') else np.asarray(source_indices, dtype=int),
        'fast_ultranest_binning': fast_binning,
        'fast_fit_good_times': np.asarray(fit_times, dtype=float),
        'fast_fit_good_flux': np.asarray(fitted_flux, dtype=float),
        'fast_fit_good_unc': np.asarray(fitted_unc, dtype=float),
        'fast_fit_good_airmass': np.asarray(fit_airmass, dtype=float),
        'fast_fit_good_jd_times': None if fit_jd_times is None else np.asarray(fit_jd_times, dtype=float),
        'fast_fit_prior': fit_prior,
        'fast_fit_bounds': fit_bounds,
        'skip_airmass_fit': skip_final_airmass_fit,
        'airmass_skip_note': airmass_skip_note,
        'data_highres': data_highres,
        'duration_samples': duration_samples,
        'failure_reason': None,
        'note': 'completed the full comparison-candidate reduction directly from the raw target/reference light curve.',
    })
    return result


def selected_final_live_point_target(enabled=None):
    if enabled is None:
        enabled = should_use_sparse_posterior_live_point_retry(
            os.environ.get(
                SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_ENV,
                SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_DEFAULT,
            )
        )
    base_live_points = get_configured_ultranest_min_num_live_points()
    if not enabled:
        return base_live_points, None
    extension_factor = int(max(1, SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT))
    target_live_points = int(max(
        base_live_points + extension_factor * base_live_points,
        base_live_points + 1,
    ))
    return base_live_points, target_live_points


def baseline_fixed_errors_from_fit(fit):
    errors = getattr(fit, 'errors', {}) if fit is not None else {}
    fixed_errors = {}
    if isinstance(errors, dict):
        for key in ('a0', 'a1', 'a2'):
            value = errors.get(key)
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(value) and value >= 0:
                fixed_errors[key] = value
    if 'a0' in fixed_errors and 'a1' not in fixed_errors:
        fixed_errors['a1'] = fixed_errors['a0']
    return fixed_errors


def build_full_resolution_final_prior_from_previous_fit(previous_fit, p_dict):
    previous_parameters = getattr(previous_fit, 'parameters', {})
    if not isinstance(previous_parameters, dict):
        previous_parameters = {}

    prior = {
        'rprs': p_dict.get('rprs', previous_parameters.get('rprs')),
        'ars': p_dict.get('aRs', previous_parameters.get('ars')),
        'per': p_dict.get('pPer', previous_parameters.get('per')),
        'inc': p_dict.get('inc', previous_parameters.get('inc')),
        'u0': previous_parameters.get('u0', 0.0),
        'u1': previous_parameters.get('u1', 0.0),
        'u2': previous_parameters.get('u2', 0.0),
        'u3': previous_parameters.get('u3', 0.0),
        'ecc': p_dict.get('ecc', previous_parameters.get('ecc', 0.0)),
        'omega': p_dict.get('omega', previous_parameters.get('omega', 0.0)),
        'tmid': p_dict.get('midT', previous_parameters.get('tmid')),
        'a2': previous_parameters.get('a2', 0.0),
        'a0': previous_parameters.get('a0', previous_parameters.get('a1', 1.0)),
    }
    prior['a1'] = previous_parameters.get('a1', prior['a0'])
    prior.update(previous_parameters)
    if 'a0' not in prior and 'a1' in prior:
        prior['a0'] = prior['a1']
    if 'a1' not in prior and 'a0' in prior:
        prior['a1'] = prior['a0']
    return prior


def refit_selected_fast_comparison_on_full_lightcurve(
    selected_result,
    p_dict,
    skip_airmass_fit=False,
    airmass_skip_note=None,
    detrend_on_outoftransit_baseline=True,
    use_impactparameter_rather_than_inclination_to_fit=True,
    plot_time_range=None,
    duration_prior=None,
    sparse_live_point_extension_enabled=None,
):
    previous_fit = selected_result.get('fit') if isinstance(selected_result, dict) else None
    if previous_fit is None or not getattr(previous_fit, 'fast_ultranest_binning_applied', False):
        return None

    times = np.asarray(selected_result.get('good_times'), dtype=float)
    flux_values = np.asarray(selected_result.get('good_flux'), dtype=float)
    flux_errors = np.asarray(selected_result.get('good_unc'), dtype=float)
    airmass = np.asarray(selected_result.get('good_airmass'), dtype=float)
    jd_times = selected_result.get('good_jd_times')
    jd_times = None if jd_times is None else np.asarray(jd_times, dtype=float)
    if not (times.shape == flux_values.shape == flux_errors.shape == airmass.shape):
        log_info(
            "Warning: Could not run the full-resolution selected comparison-star final fit "
            "because the saved fast-fit light-curve arrays were not aligned.",
            warn=True,
        )
        return None
    if jd_times is not None and jd_times.shape != times.shape:
        jd_times = None

    prior = build_full_resolution_final_prior_from_previous_fit(previous_fit, p_dict)
    fallback_bounds = selected_result.get('fast_fit_bounds')
    if not isinstance(fallback_bounds, dict):
        fallback_bounds = getattr(previous_fit, 'bounds', {})
    bounds = get_posterior_refit_final_bounds(previous_fit, fallback_bounds)
    bounds = clone_lightcurve_bounds(bounds)
    for key in ('a0', 'a1', 'a2'):
        bounds.pop(key, None)

    for key in ('rprs', 'tmid', 'ars', 'inc'):
        if key not in bounds:
            if key == 'rprs':
                bounds[key] = build_initial_rprs_bounds(prior.get('rprs', p_dict.get('rprs', 0.1)))
            elif key == 'tmid':
                tmid = prior.get('tmid', p_dict.get('midT', np.nan))
                tmid_unc = p_dict.get('midTUnc', 0.01)
                try:
                    half_width = max(float(tmid_unc) * 3.0, np.finfo(float).eps)
                except (TypeError, ValueError):
                    half_width = 0.01
                bounds[key] = [float(tmid) - half_width, float(tmid) + half_width]
            elif key == 'ars':
                bounds[key] = build_initial_ars_bounds(prior.get('ars', p_dict.get('aRs')), p_dict.get('aRsUnc'))
            elif key == 'inc':
                inc = float(prior.get('inc', p_dict.get('inc', 89.0)))
                bounds[key] = [inc - 5.0, min(90.0, inc + 5.0)]

    fit_flux = flux_values
    fit_unc = flux_errors
    detrend_result = {'applied': False, 'note': 'Disabled; using the full-resolution light curve directly.'}
    if detrend_on_outoftransit_baseline:
        detrend_result = detrend_flux_on_out_of_transit_baseline(
            times,
            flux_values,
            flux_errors,
            previous_fit,
            prior=prior,
        )
        if detrend_result.get('applied'):
            fit_flux = np.asarray(detrend_result['flux'], dtype=float)
            fit_unc = np.asarray(detrend_result['unc'], dtype=float)

    fixed_errors = baseline_fixed_errors_from_fit(previous_fit)
    base_live_points, target_live_points = selected_final_live_point_target(
        sparse_live_point_extension_enabled,
    )
    min_live_points = target_live_points if target_live_points is not None else base_live_points
    coverage_duration_prior = (
        duration_prior if isinstance(duration_prior, dict) else build_single_transit_duration_prior(p_dict)
    )
    pre_ultranest_coverage_assessment = build_expected_transit_coverage_assessment(
        times,
        prior,
        flux_values=fit_flux,
        flux_errors=fit_unc,
        tmid_search_summary=build_ephemeris_tmid_search_summary_for_coverage(
            times,
            p_dict,
            prior=prior,
            duration_prior=coverage_duration_prior,
            sigma_multiplier=35.0,
        ),
        duration_prior=coverage_duration_prior,
    )
    log_expected_transit_coverage_assessment(pre_ultranest_coverage_assessment)

    log_info(
        "Running the selected comparison-star final UltraNest fit on the full-resolution light curve "
        f"with fixed a0/a2 from the previous fast UltraNest fit at {min_live_points} minimum live points."
    )
    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times,
        fit_flux,
        fit_unc,
        airmass,
        prior,
        bounds,
        jd_times=jd_times,
        use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
        duration_prior=duration_prior,
        keep_ultranest_sampler=False,
        fixed_parameter_errors=fixed_errors,
        fixed_flux_baseline=True,
        ultranest_min_num_live_points=min_live_points,
        pre_ultranest_coverage_assessment=pre_ultranest_coverage_assessment,
    )
    annotate_pre_ultranest_transit_coverage(fit, pre_ultranest_coverage_assessment)
    fit = apply_plot_time_range(fit, times if plot_time_range is None else plot_time_range)
    annotate_airmass_fit(fit, airmass, skip_airmass_fit, note=airmass_skip_note)
    annotate_out_of_transit_baseline_parameter_fit(
        fit,
        True,
        note="Used a0 and a2 from the previous fast UltraNest fit for the full-resolution final run.",
        pre_points=0,
        post_points=0,
        a0=prior.get('a0'),
        a0_error=fixed_errors.get('a0'),
        a2=prior.get('a2'),
        a2_error=fixed_errors.get('a2'),
    )
    annotate_out_of_transit_baseline_detrending(
        fit,
        bool(detrend_result.get('applied')),
        note=detrend_result.get('note'),
        slope=detrend_result.get('slope'),
        intercept=detrend_result.get('intercept'),
        pre_points=detrend_result.get('pre_points', 0),
        post_points=detrend_result.get('post_points', 0),
    )
    annotate_fast_ultranest_binning(
        fit,
        {
            'applied': False,
            'original_point_count': int(times.shape[0]),
            'binned_point_count': int(times.shape[0]),
            'note': 'Full-resolution selected comparison-star final run; fast binning was not applied.',
        },
    )
    if target_live_points is not None:
        diagnostics = evaluate_sparse_posterior_sample_support(fit, base_live_points=base_live_points)
        annotate_sparse_posterior_live_point_extension(
            fit,
            True,
            True,
            note=(
                "Applied full-resolution selected comparison-star final UltraNest run "
                f"({base_live_points}->{target_live_points} minimum live points) using fixed a0/a2 "
                "from the previous fast UltraNest fit."
            ),
            diagnostics=diagnostics,
            post_extension_diagnostics=diagnostics,
            base_live_points=base_live_points,
            target_live_points=target_live_points,
            extension_factor=SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT,
        )
    else:
        annotate_sparse_posterior_live_point_extension(fit, False, False)
    annotate_transit_detection_qc(fit)
    clear_fit_ultranest_resume_state(fit)
    return fit, fit_flux, fit_unc


def save_comparison_candidate_full_reduction_outputs(save_dir, provisional_fit, final_fit,
                                                     p_dict, observation_date, comp_index,
                                                     comp_coords=None, min_aperture=None, min_annulus=None,
                                                     adaptive_summary=None, method_label=None,
                                                     selection_summary=None, duration_samples=None,
                                                     data_highres=None):
    if save_dir is None or final_fit is None or observation_date is None or comp_index is None:
        return None

    candidate_dir = comparison_candidate_output_dir(save_dir, comp_index)
    temp_dir = candidate_dir / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    archive_errors = []
    debug_series_path = None
    bestfit_plot_path = None
    triangle_plot_path = None

    debug_fit = provisional_fit if provisional_fit is not None else final_fit
    if debug_fit is not None:
        try:
            debug_series_path = save_selected_photometry_debug_series(
                candidate_dir,
                p_dict['pName'],
                observation_date,
                debug_fit,
            )
        except Exception as exc:
            archive_errors.append(archive_exception_payload(
                "Could not save the selected raw target/reference ratio diagnostics",
                exc,
            ))

    plotter = getattr(final_fit, 'plot_bestfit', None)
    if callable(plotter):
        try:
            plot_kwargs = {}
            if callable_accepts_keyword(plotter, 'show_flux_baseline_label'):
                plot_kwargs['show_flux_baseline_label'] = False
            fig, _ = plotter(**plot_kwargs)
            bestfit_plot_path = temp_dir / safe_output_filename(
                "BestFit",
                p_dict['pName'],
                filename_date_token(observation_date),
                extension="png",
            )
            fig.savefig(bestfit_plot_path)
            plt.close(fig)
        except Exception as exc:
            archive_errors.append(archive_exception_payload(
                "Could not save the final best-fit plot",
                exc,
            ))

    if callable(getattr(final_fit, 'plot_triangle', None)):
        try:
            fig = _plot_triangle_for_output(
                final_fit,
                plot_title=f"Comparison candidate #{int(comp_index) + 1} fit",
            )
            triangle_plot_path = comparison_candidate_triangle_plot_output_path(
                candidate_dir,
                p_dict['pName'],
                observation_date,
                comp_index,
            )
            if fig is not None:
                fig.savefig(triangle_plot_path)
                plt.close(fig)
        except Exception as exc:
            archive_errors.append(archive_exception_payload(
                "Could not save the triangle plot",
                exc,
            ))

    duration_samples = np.asarray([] if duration_samples is None else duration_samples, dtype=float)
    if duration_samples.size == 0:
        measured_duration = getattr(final_fit, 'duration_measured', np.nan)
        if np.isfinite(measured_duration) and measured_duration > 0:
            duration_samples = np.asarray([measured_duration], dtype=float)

    candidate_info_dict = {
        'save': str(candidate_dir),
        'date': observation_date,
    }
    try:
        if data_highres is None:
            data_highres, _ = estimate_transit_duration_samples_from_fit(final_fit, sample_count=1)
        if data_highres is not None:
            plot_final_lightcurve(final_fit, data_highres, p_dict['pName'], candidate_info_dict['save'], observation_date)
    except Exception as exc:
        archive_errors.append(archive_exception_payload(
            "Could not save the final lightcurve plot",
            exc,
        ))

    output_files = OutputFiles(final_fit, p_dict, candidate_info_dict, duration_samples)
    try:
        phase = get_phase(final_fit.time, p_dict['pPer'], final_fit.parameters['tmid'])
        output_files.final_lightcurve(phase)
    except Exception as exc:
        archive_errors.append(archive_exception_payload(
            "Could not save FinalLightCurve CSV",
            exc,
        ))

    try:
        output_files.final_planetary_params(
            phot_opt=True,
            vsp_params=[],
            comp_star=int(comp_index + 1),
            comp_coords=comp_coords,
            min_aper=0 if min_aperture is None else np.round(min_aperture, 2),
            min_annul=(None if min_annulus is None else np.round(min_annulus, 2)),
            adaptive_summary=adaptive_summary,
        )
    except Exception as exc:
        archive_errors.append(archive_exception_payload(
            "Could not save FinalParams JSON",
            exc,
        ))

    summary_path = temp_dir / safe_output_filename(
        "ComparisonCandidateSummary",
        p_dict['pName'],
        filename_date_token(observation_date),
        extension="json",
    )
    summary_payload = {
        'planet_name': p_dict['pName'],
        'observation_date': observation_date,
        'comparison_star': int(comp_index + 1),
        'comparison_position': comp_coords,
        'method_label': method_label,
        'selection_summary': selection_summary or {},
        'parameter_summary': summarize_lightcurve_fit_parameters(final_fit),
        'transit_qc': getattr(final_fit, 'transit_qc', None),
        'saved_debug_series': None if debug_series_path is None else str(debug_series_path),
        'saved_bestfit_plot': None if bestfit_plot_path is None else str(bestfit_plot_path),
        'saved_triangle_plot': None if triangle_plot_path is None else str(triangle_plot_path),
        'archive_errors': archive_errors,
    }
    with summary_path.open('w', encoding='utf-8') as handle:
        json.dump(make_json_safe(summary_payload), handle, indent=4)

    return candidate_dir


def archive_failed_comparison_fit(save_dir, planet_name, observation_date, attempt, method_label=None):
    if save_dir is None or planet_name is None or observation_date is None or not attempt:
        return None

    comp_index = attempt.get('comp_index')
    if comp_index is None:
        return None

    archive_dir = failed_comparison_archive_dir(save_dir, comp_index)
    temp_dir = archive_dir / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    fit = attempt.get('fit')
    archive_errors = []
    debug_series_path = None
    bestfit_plot_path = None

    if fit is not None:
        try:
            debug_series_path = save_selected_photometry_debug_series(
                archive_dir,
                planet_name,
                observation_date,
                fit,
            )
        except Exception as exc:
            archive_errors.append(archive_exception_payload(
                "Could not save the selected raw target/reference ratio diagnostics",
                exc,
            ))

        plotter = getattr(fit, 'plot_bestfit', None)
        if callable(plotter):
            try:
                plot_kwargs = {}
                if callable_accepts_keyword(plotter, 'show_flux_baseline_label'):
                    plot_kwargs['show_flux_baseline_label'] = False
                fig, _ = plotter(**plot_kwargs)
                bestfit_plot_path = temp_dir / safe_output_filename(
                    "BestFit",
                    planet_name,
                    filename_date_token(observation_date),
                    extension="png",
                )
                fig.savefig(bestfit_plot_path)
                plt.close(fig)
            except Exception as exc:
                archive_errors.append(archive_exception_payload(
                    "Could not save the provisional best-fit plot",
                    exc,
                ))

    summary_path = temp_dir / safe_output_filename(
        "FailedFitSummary",
        planet_name,
        filename_date_token(observation_date),
        extension="json",
    )
    summary_payload = {
        'planet_name': planet_name,
        'observation_date': observation_date,
        'comparison_star': None if comp_index is None else int(comp_index + 1),
        'comparison_label': attempt.get('label', f"Comp {comp_index + 1}"),
        'comparison_position': attempt.get('position'),
        'method_label': method_label,
        'failure_reason': attempt.get('failure_reason'),
        'fit_diagnostics': attempt.get('fit_diagnostics') or {},
        'parameter_summary': attempt.get('parameter_summary'),
        'fit_point_count': attempt.get('fit_point_count'),
        'eebls_snr': attempt.get('eebls_snr'),
        'transit_delta_bic': attempt.get('transit_delta_bic'),
        'residual_scatter': attempt.get('residual_scatter'),
        'ktmf_metric': attempt.get('ktmf_metric'),
        'ktmf_contributions': attempt.get('ktmf_contributions') or [],
        'transit_qc': getattr(fit, 'transit_qc', None) if fit is not None else None,
        'saved_debug_series': None if debug_series_path is None else str(debug_series_path),
        'saved_bestfit_plot': None if bestfit_plot_path is None else str(bestfit_plot_path),
        'archive_errors': archive_errors,
    }
    with summary_path.open('w', encoding='utf-8') as handle:
        json.dump(make_json_safe(summary_payload), handle, indent=4)

    return archive_dir


def save_selected_photometry_debug_series(save_dir, planet_name, observation_date, fit):
    if fit is None:
        return None

    debug = getattr(fit, 'selected_photometry_debug', None)
    if not debug:
        return None

    times = np.asarray(debug.get('times'), dtype=float)
    target_flux = np.asarray(debug.get('target_flux'), dtype=float)
    comp_flux = np.asarray(debug.get('comp_flux'), dtype=float)
    raw_ratio = np.asarray(debug.get('raw_ratio'), dtype=float)
    initial_sigma_keep_mask = np.asarray(debug.get('initial_sigma_keep_mask'), dtype=bool)
    phase_clip_keep_mask = np.asarray(
        debug.get('phase_clip_keep_mask_on_sigma_filtered', np.ones(np.count_nonzero(initial_sigma_keep_mask))),
        dtype=bool,
    )

    if not (
        times.shape == target_flux.shape == comp_flux.shape == raw_ratio.shape == initial_sigma_keep_mask.shape
    ):
        return None

    phase_keep_full = np.zeros(times.shape[0], dtype=bool)
    sigma_kept_indices = np.flatnonzero(initial_sigma_keep_mask)
    if sigma_kept_indices.size:
        if phase_clip_keep_mask.shape[0] != sigma_kept_indices.size:
            phase_clip_keep_mask = np.ones(sigma_kept_indices.size, dtype=bool)
        phase_keep_full[sigma_kept_indices] = phase_clip_keep_mask

    output_dir = Path(save_dir) / "temp"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / safe_output_filename(
        "SelectedPhotometryRawRatio",
        planet_name,
        filename_date_token(observation_date),
        extension="csv",
    )

    output_rows = np.column_stack(
        [
            times,
            target_flux,
            comp_flux,
            raw_ratio,
            initial_sigma_keep_mask.astype(int),
            phase_keep_full.astype(int),
        ]
    )
    np.savetxt(
        output_path,
        output_rows,
        delimiter=",",
        header=(
            "BJD_TDB,Target Flux,Comp Flux,Raw Ratio,"
            "Kept After Initial Sigma Clip,Kept After Phase Residual Clip"
        ),
        comments="",
        fmt=["%.8f", "%.8f", "%.8f", "%.8f", "%d", "%d"],
    )
    return output_path


def annotate_rprs_posterior_refit(fit, applied, note=None, history=None):
    annotate_parameter_posterior_refit(fit, 'rprs', applied, note=note, history=history)


def annotate_ars_posterior_refit(fit, applied, note=None, history=None):
    annotate_parameter_posterior_refit(fit, 'ars', applied, note=note, history=history)


def annotate_impact_parameter_posterior_refit(fit, applied, note=None, history=None):
    annotate_parameter_posterior_refit(fit, 'b', applied, note=note, history=history)


def annotate_parameter_posterior_refit(fit, parameter_key, applied, note=None, history=None):
    if fit is None:
        return

    history = [] if history is None else list(history)
    attr_prefix = f"{parameter_key}_posterior_refit"
    setattr(fit, f"{attr_prefix}_applied", bool(applied))
    setattr(fit, f"{attr_prefix}_note", note)
    setattr(fit, f"{attr_prefix}_count", len(history))
    setattr(fit, f"{attr_prefix}_history", history)
    if history:
        latest = history[-1]
        setattr(fit, f"{attr_prefix}_edge", latest.get('edge'))
        setattr(fit, f"{attr_prefix}_mode", latest.get('mode'))
        setattr(fit, f"{attr_prefix}_std", latest.get('std'))
        setattr(fit, f"{attr_prefix}_original_bounds", latest.get('original_bounds'))
        setattr(fit, f"{attr_prefix}_bounds", latest.get('new_bounds'))
    else:
        setattr(fit, f"{attr_prefix}_edge", None)
        setattr(fit, f"{attr_prefix}_mode", None)
        setattr(fit, f"{attr_prefix}_std", None)
        setattr(fit, f"{attr_prefix}_original_bounds", None)
        setattr(fit, f"{attr_prefix}_bounds", None)


def annotate_sparse_posterior_live_point_extension(
    fit,
    enabled,
    applied,
    note=None,
    diagnostics=None,
    post_extension_diagnostics=None,
    base_live_points=None,
    target_live_points=None,
    extension_factor=SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT,
):
    if fit is None:
        return

    fit.sparse_posterior_live_point_extension_enabled = bool(enabled)
    fit.sparse_posterior_live_point_extension_applied = bool(applied)
    fit.sparse_posterior_live_point_extension_note = note
    fit.sparse_posterior_live_point_extension_diagnostics = diagnostics
    fit.sparse_posterior_live_point_extension_post_diagnostics = post_extension_diagnostics
    fit.sparse_posterior_live_point_extension_base_live_points = base_live_points
    fit.sparse_posterior_live_point_extension_target_live_points = target_live_points
    fit.sparse_posterior_live_point_extension_factor = extension_factor


def clear_fit_ultranest_resume_state(fit):
    clear_resume_state = getattr(fit, 'clear_ultranest_resume_state', None)
    if callable(clear_resume_state):
        clear_resume_state()


def callable_accepts_keyword(callable_obj, keyword):
    try:
        signature = inspect.signature(callable_obj)
    except (TypeError, ValueError):
        return False

    if keyword in signature.parameters:
        return True
    return any(
        parameter.kind == inspect.Parameter.VAR_KEYWORD
        for parameter in signature.parameters.values()
    )


def get_configured_ultranest_min_num_live_points():
    return parse_ultranest_min_num_live_points(
        os.environ.get(
            ULTRANEST_MIN_NUM_LIVE_POINTS_ENV,
            ULTRANEST_MIN_NUM_LIVE_POINTS_DEFAULT,
        )
    )


def _effective_sample_count(weights, fallback_count):
    if weights is None:
        return float(fallback_count)

    weights = np.asarray(weights, dtype=float)
    finite_weights = weights[np.isfinite(weights) & (weights > 0)]
    if finite_weights.size == 0:
        return float(fallback_count)

    weight_sum = float(np.sum(finite_weights))
    weight_square_sum = float(np.sum(finite_weights ** 2))
    if not np.isfinite(weight_sum) or not np.isfinite(weight_square_sum) or weight_square_sum <= 0:
        return float(fallback_count)
    return float((weight_sum ** 2) / weight_square_sum)


def _fit_posterior_sample_matrix(fit, parameter_keys):
    parameter_keys = list(parameter_keys)
    if fit is None or not parameter_keys:
        return np.empty((0, 0), dtype=float), None

    sample_points = None
    sample_weights = None
    try:
        sample_points, _, sample_weights = fit._get_triangle_plot_samples()
    except Exception:
        sample_points = None

    if sample_points is not None:
        sample_points = np.asarray(sample_points, dtype=float)
        if sample_points.ndim == 2 and sample_points.shape[0] > 0:
            sampled_keys = list(getattr(fit, 'sampled_keys', []))
            bounds = getattr(fit, 'bounds', {})
            bound_keys = list(bounds.keys()) if isinstance(bounds, dict) else []
            physical_getter = getattr(fit, '_physical_values_from_sample_point', None)
            columns = []
            for key in parameter_keys:
                if key in sampled_keys:
                    key_index = sampled_keys.index(key)
                    if key_index >= sample_points.shape[1]:
                        return np.empty((0, len(parameter_keys)), dtype=float), None
                    columns.append(np.asarray(sample_points[:, key_index], dtype=float))
                elif callable(physical_getter) and bound_keys:
                    columns.append(np.asarray([
                        physical_getter(point, bound_keys, sampled_keys).get(key, np.nan)
                        for point in sample_points
                    ], dtype=float))
                else:
                    break
            else:
                weights = None
                if sample_weights is not None:
                    sample_weights = np.asarray(sample_weights, dtype=float)
                    if sample_weights.ndim == 1 and sample_weights.shape[0] == sample_points.shape[0]:
                        weights = sample_weights
                return np.column_stack(columns), weights

    sample_getter = getattr(fit, 'get_parameter_posterior_samples', None)
    if not callable(sample_getter):
        return np.empty((0, len(parameter_keys)), dtype=float), None

    columns = []
    min_size = None
    for key in parameter_keys:
        values = np.asarray(sample_getter(key), dtype=float).reshape(-1)
        columns.append(values)
        min_size = values.size if min_size is None else min(min_size, values.size)

    if min_size is None or min_size == 0:
        return np.empty((0, len(parameter_keys)), dtype=float), None

    return np.column_stack([values[:min_size] for values in columns]), None


def evaluate_sparse_posterior_sample_support(
    fit,
    parameter_keys=SPARSE_POSTERIOR_RETRY_PARAMETER_KEYS,
    base_live_points=None,
    minimum_effective_samples=None,
    minimum_occupied_bins=SPARSE_POSTERIOR_MIN_OCCUPIED_BINS,
    minimum_occupied_bin_fraction=SPARSE_POSTERIOR_MIN_OCCUPIED_BIN_FRACTION,
    minimum_effective_samples_per_occupied_bin=SPARSE_POSTERIOR_MIN_EFFECTIVE_SAMPLES_PER_OCCUPIED_BIN,
):
    if base_live_points is None:
        base_live_points = get_configured_ultranest_min_num_live_points()

    if minimum_effective_samples is None:
        minimum_effective_samples = max(
            SPARSE_POSTERIOR_MIN_EFFECTIVE_SAMPLES_FLOOR,
            int(np.ceil(SPARSE_POSTERIOR_MIN_EFFECTIVE_SAMPLES_PER_LIVE_POINT * float(base_live_points))),
        )
    minimum_effective_samples = int(max(1, minimum_effective_samples))
    minimum_occupied_bins = int(max(1, minimum_occupied_bins))
    minimum_occupied_bin_fraction = float(np.clip(minimum_occupied_bin_fraction, 0.0, 1.0))
    minimum_effective_samples_per_occupied_bin = float(
        max(0.0, minimum_effective_samples_per_occupied_bin)
    )

    sample_matrix, sample_weights = _fit_posterior_sample_matrix(fit, parameter_keys)
    diagnostics = {
        'sparse': False,
        'reason': None,
        'parameter_keys': list(parameter_keys),
        'base_live_points': int(base_live_points),
        'minimum_effective_samples': minimum_effective_samples,
        'minimum_occupied_bins': minimum_occupied_bins,
        'minimum_occupied_bin_fraction': minimum_occupied_bin_fraction,
        'minimum_effective_samples_per_occupied_bin': minimum_effective_samples_per_occupied_bin,
        'parameters': {},
    }

    if sample_matrix.size == 0 or sample_matrix.shape[0] == 0:
        diagnostics['sparse'] = True
        diagnostics['reason'] = "posterior samples are unavailable for Rp/R*, Tmid, and a/Rs."
        return diagnostics

    sparse_reasons = []
    for column_index, key in enumerate(parameter_keys):
        if column_index >= sample_matrix.shape[1]:
            sample_values = np.array([], dtype=float)
        else:
            sample_values = np.asarray(sample_matrix[:, column_index], dtype=float)
        finite_mask = np.isfinite(sample_values)
        finite_values = sample_values[finite_mask]
        parameter_weights = sample_weights[finite_mask] if sample_weights is not None else None
        sample_count = int(finite_values.size)
        effective_count = _effective_sample_count(parameter_weights, sample_count)

        occupied_bins = 0
        central_count = 0
        bin_count = 0
        occupied_bin_fraction = 0.0
        effective_samples_per_occupied_bin = 0.0
        if sample_count >= 2:
            q05, q95 = np.nanpercentile(finite_values, [5, 95])
            central_mask = (finite_values >= q05) & (finite_values <= q95)
            central_values = finite_values[central_mask]
            central_count = int(central_values.size)
            if np.isfinite(q05) and np.isfinite(q95) and q05 < q95 and central_count > 0:
                bin_count = int(np.clip(np.sqrt(sample_count), 10, 40))
                hist_counts, _ = np.histogram(central_values, bins=bin_count, range=(q05, q95))
                occupied_bins = int(np.count_nonzero(hist_counts > 0))
                occupied_bin_fraction = (
                    float(occupied_bins) / float(bin_count)
                    if bin_count > 0
                    else 0.0
                )
                if occupied_bins > 0:
                    effective_samples_per_occupied_bin = float(effective_count) / float(occupied_bins)

        parameter_diagnostic = {
            'sample_count': sample_count,
            'effective_sample_count': float(effective_count),
            'central_sample_count': central_count,
            'central_bin_count': bin_count,
            'occupied_bins': occupied_bins,
            'occupied_bin_fraction': occupied_bin_fraction,
            'effective_samples_per_occupied_bin': effective_samples_per_occupied_bin,
            'sparse': False,
            'reason': None,
        }

        if effective_count < minimum_effective_samples:
            parameter_diagnostic['sparse'] = True
            parameter_diagnostic['reason'] = (
                f"effective samples {effective_count:.0f} < {minimum_effective_samples}"
            )
        elif occupied_bins and occupied_bins < minimum_occupied_bins:
            parameter_diagnostic['sparse'] = True
            parameter_diagnostic['reason'] = (
                f"central posterior occupies {occupied_bins} histogram bins < {minimum_occupied_bins}"
            )
        elif bin_count and occupied_bin_fraction < minimum_occupied_bin_fraction:
            parameter_diagnostic['sparse'] = True
            parameter_diagnostic['reason'] = (
                f"central posterior occupies {occupied_bin_fraction:.2f} of histogram bins "
                f"< {minimum_occupied_bin_fraction:.2f}"
            )
        elif (
            occupied_bins
            and minimum_effective_samples_per_occupied_bin > 0
            and effective_samples_per_occupied_bin < minimum_effective_samples_per_occupied_bin
        ):
            parameter_diagnostic['sparse'] = True
            parameter_diagnostic['reason'] = (
                f"effective samples per occupied bin {effective_samples_per_occupied_bin:.1f} "
                f"< {minimum_effective_samples_per_occupied_bin:.1f}"
            )

        if parameter_diagnostic['sparse']:
            sparse_reasons.append(f"{key}: {parameter_diagnostic['reason']}")
        diagnostics['parameters'][key] = parameter_diagnostic

    if sparse_reasons:
        diagnostics['sparse'] = True
        diagnostics['reason'] = "; ".join(sparse_reasons)
    else:
        diagnostics['reason'] = "posterior sample support is sufficient for Rp/R*, Tmid, and a/Rs."

    return diagnostics


def sparse_posterior_diagnostics_summary(diagnostics):
    if not isinstance(diagnostics, dict):
        return "posterior sample support diagnostics are unavailable"
    reason = diagnostics.get('reason')
    if reason:
        return str(reason)
    return "posterior sample support diagnostics are unavailable"


def extend_sparse_posterior_live_points_if_needed(
    fit,
    enabled=None,
    extension_factor=SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT,
    require_sparse=True,
    extension_label="sparse-posterior",
):
    if enabled is None:
        enabled = should_use_sparse_posterior_live_point_retry(
            os.environ.get(
                SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_ENV,
                SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_DEFAULT,
            )
        )

    if not enabled:
        annotate_sparse_posterior_live_point_extension(fit, False, False)
        clear_fit_ultranest_resume_state(fit)
        return fit

    base_live_points = get_configured_ultranest_min_num_live_points()
    diagnostics = evaluate_sparse_posterior_sample_support(fit, base_live_points=base_live_points)
    if not diagnostics.get('sparse'):
        if require_sparse:
            annotate_sparse_posterior_live_point_extension(
                fit,
                True,
                False,
                note=f"Not needed; {sparse_posterior_diagnostics_summary(diagnostics)}",
                diagnostics=diagnostics,
                base_live_points=base_live_points,
                extension_factor=extension_factor,
            )
            clear_fit_ultranest_resume_state(fit)
            return fit

    extender = getattr(fit, 'extend_ultranest_fit', None)
    if not callable(extender):
        note = (
            "Skipped; the retained UltraNest sampler state is unavailable for an additive "
            f"{extension_label} live-point extension."
        )
        log_info(f"Warning: {note}", warn=True)
        annotate_sparse_posterior_live_point_extension(
            fit,
            True,
            False,
            note=note,
            diagnostics=diagnostics,
            base_live_points=base_live_points,
            extension_factor=extension_factor,
        )
        clear_fit_ultranest_resume_state(fit)
        return fit

    extension_factor = int(max(1, extension_factor))
    target_live_points = int(max(
        base_live_points + extension_factor * base_live_points,
        base_live_points + 1,
    ))
    try:
        current_max_ncalls = int(float(getattr(fit, 'max_ncalls', 2e5)))
    except (TypeError, ValueError):
        current_max_ncalls = int(2e5)
    target_max_ncalls = int(max(current_max_ncalls, current_max_ncalls * (extension_factor + 1)))
    if diagnostics.get('sparse'):
        log_info(
            "Posterior samples for Rp/R*, Tmid, and a/Rs are sparse "
            f"({sparse_posterior_diagnostics_summary(diagnostics)}); continuing UltraNest "
            f"from {base_live_points} to {target_live_points} minimum live points "
            "using the retained final-pass sampler bounds."
        )
    else:
        log_info(
            f"Continuing the {extension_label} UltraNest fit from {base_live_points} "
            f"to {target_live_points} minimum live points using the retained final-pass "
            f"sampler bounds ({sparse_posterior_diagnostics_summary(diagnostics)})."
        )
    applied = bool(extender(min_num_live_points=target_live_points, max_ncalls=target_max_ncalls))
    post_diagnostics = evaluate_sparse_posterior_sample_support(fit, base_live_points=base_live_points)

    if applied and post_diagnostics.get('sparse'):
        note = (
            f"Applied additive {extension_label} UltraNest extension "
            f"({base_live_points}->{target_live_points} minimum live points), but "
            f"{sparse_posterior_diagnostics_summary(post_diagnostics)}"
        )
        log_info(
            "Warning: sparse posterior support remains after the additive UltraNest extension; "
            "please inspect the triangle plot carefully.",
            warn=True,
        )
    elif applied:
        note = (
            f"Applied additive {extension_label} UltraNest extension "
            f"({base_live_points}->{target_live_points} minimum live points)."
        )
    else:
        note = (
            f"Skipped; UltraNest did not continue the additive {extension_label} extension "
            "from the retained sampler state."
        )

    annotate_sparse_posterior_live_point_extension(
        fit,
        True,
        applied,
        note=note,
        diagnostics=diagnostics,
        post_extension_diagnostics=post_diagnostics,
        base_live_points=base_live_points,
        target_live_points=target_live_points,
        extension_factor=extension_factor,
    )
    clear_fit_ultranest_resume_state(fit)
    return fit


def extend_selected_comparison_live_points_if_needed(fit, enabled=None):
    return extend_sparse_posterior_live_points_if_needed(
        fit,
        enabled=enabled,
        extension_factor=SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT,
        require_sparse=False,
        extension_label="selected comparison-star final",
    )


def build_initial_rprs_bounds(
    rprs,
    lower_scale=INITIAL_RPRS_BOUND_LOWER_SCALE,
    upper_scale=INITIAL_RPRS_BOUND_UPPER_SCALE,
):
    try:
        rprs = float(rprs)
        lower_scale = float(lower_scale)
        upper_scale = float(upper_scale)
    except (TypeError, ValueError):
        return [RPRS_SEARCH_BOUND_MIN, RPRS_SEARCH_BOUND_MAX]

    if not np.isfinite(rprs) or rprs <= 0:
        return [RPRS_SEARCH_BOUND_MIN, RPRS_SEARCH_BOUND_MAX]
    if rprs >= RPRS_SEARCH_BOUND_MAX:
        return [RPRS_SEARCH_BOUND_MIN, RPRS_SEARCH_BOUND_MAX]

    lower_bound = max(RPRS_SEARCH_BOUND_MIN, lower_scale * rprs)
    upper_bound = min(RPRS_SEARCH_BOUND_MAX, upper_scale * rprs)
    if not np.isfinite(upper_bound) or upper_bound <= lower_bound:
        lower_bound = RPRS_SEARCH_BOUND_MIN
        upper_bound = RPRS_SEARCH_BOUND_MAX

    return [float(lower_bound), float(upper_bound)]


def build_initial_ars_bounds(
    ars,
    ars_unc=None,
    sigma_multiplier=INITIAL_ARS_BOUND_SIGMA_MULTIPLIER,
    fallback_relative_half_width=INITIAL_ARS_BOUND_FALLBACK_RELATIVE_HALF_WIDTH,
):
    try:
        ars = float(ars)
    except (TypeError, ValueError):
        ars = np.nan

    try:
        ars_unc = float(ars_unc)
    except (TypeError, ValueError):
        ars_unc = np.nan

    if not np.isfinite(ars) or ars <= ARS_SEARCH_BOUND_MIN:
        return [float(ARS_SEARCH_BOUND_MIN), float(ARS_SEARCH_BOUND_FALLBACK_MAX)]

    if np.isfinite(ars_unc) and ars_unc > 0:
        half_width = float(max(ARS_SEARCH_BOUND_MIN, sigma_multiplier * ars_unc))
    else:
        half_width = float(max(ARS_SEARCH_BOUND_MIN, fallback_relative_half_width * ars))

    lower_bound = max(float(ARS_SEARCH_BOUND_MIN), float(ars - half_width))
    upper_bound = float(ars + half_width)
    if not np.isfinite(upper_bound) or upper_bound <= lower_bound:
        upper_bound = float(lower_bound + max(np.finfo(float).eps, ARS_SEARCH_BOUND_MIN))

    return [float(lower_bound), float(upper_bound)]


def build_initial_transit_bounds(prior, tmid_bounds, ars_unc=None, inclination_half_width=5.0):
    lower, upper = [float(value) for value in np.asarray(tmid_bounds, dtype=float).reshape(-1)[:2]]
    # Keep ars ahead of inc so the internal impact-parameter parameterization
    # uses the sampled ars value when converting inclination to b.
    return {
        'rprs': build_initial_rprs_bounds(prior['rprs']),
        'tmid': [lower, upper],
        'ars': build_initial_ars_bounds(prior['ars'], ars_unc=ars_unc),
        'inc': [prior['inc'] - inclination_half_width, min(90, prior['inc'] + inclination_half_width)],
    }


def clone_lightcurve_bounds(bounds):
    return {
        key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
        for key, value in bounds.items()
    }


def annotate_posterior_refit_final_bounds(fit, bounds):
    if fit is None:
        return
    fit.posterior_refit_final_bounds = clone_lightcurve_bounds(bounds)


def get_posterior_refit_final_bounds(fit, fallback_bounds):
    effective_bounds = clone_lightcurve_bounds(fallback_bounds)
    fit_bounds = getattr(fit, 'posterior_refit_final_bounds', None)
    if not isinstance(fit_bounds, dict):
        fit_bounds = {}
        for key in ('rprs', 'ars', 'inc'):
            refit_bounds = getattr(fit, f'{key}_posterior_refit_bounds', None)
            if refit_bounds is not None:
                fit_bounds[key] = refit_bounds

    for key, value in fit_bounds.items():
        if not isinstance(value, (list, tuple, np.ndarray)):
            continue
        try:
            lower_bound, upper_bound = [
                float(bound) for bound in np.asarray(value, dtype=float).reshape(-1)[:2]
            ]
        except (TypeError, ValueError, IndexError):
            continue
        if (
            np.isfinite(lower_bound)
            and np.isfinite(upper_bound)
            and lower_bound < upper_bound
        ):
            effective_bounds[key] = [lower_bound, upper_bound]
    return effective_bounds


def sanitize_parameter_search_bounds(bounds, key, minimum_bound, maximum_bound=None, fallback_maximum=None):
    sanitized = clone_lightcurve_bounds(bounds)
    if key not in sanitized:
        return sanitized

    try:
        lower_bound, upper_bound = [
            float(value) for value in np.asarray(sanitized[key], dtype=float).reshape(-1)[:2]
        ]
    except (TypeError, ValueError, IndexError):
        sanitized[key] = [
            float(minimum_bound),
            float(maximum_bound if maximum_bound is not None else fallback_maximum),
        ]
        return sanitized

    if not np.isfinite(lower_bound) or not np.isfinite(upper_bound):
        sanitized[key] = [
            float(minimum_bound),
            float(maximum_bound if maximum_bound is not None else fallback_maximum),
        ]
        return sanitized

    lower_bound = max(float(minimum_bound), float(lower_bound))
    if maximum_bound is not None:
        upper_bound = min(float(maximum_bound), float(upper_bound))
    if lower_bound >= upper_bound:
        sanitized[key] = [
            float(minimum_bound),
            float(maximum_bound if maximum_bound is not None else fallback_maximum),
        ]
    else:
        sanitized[key] = [lower_bound, float(upper_bound)]
    return sanitized


def sanitize_rprs_search_bounds(bounds):
    return sanitize_parameter_search_bounds(
        bounds,
        'rprs',
        RPRS_SEARCH_BOUND_MIN,
        maximum_bound=RPRS_SEARCH_BOUND_MAX,
        fallback_maximum=RPRS_SEARCH_BOUND_MAX,
    )


def sanitize_ars_search_bounds(bounds):
    return sanitize_parameter_search_bounds(
        bounds,
        'ars',
        ARS_SEARCH_BOUND_MIN,
        fallback_maximum=ARS_SEARCH_BOUND_FALLBACK_MAX,
    )


def sanitize_inclination_search_bounds(bounds):
    return sanitize_parameter_search_bounds(
        bounds,
        'inc',
        INCLINATION_SEARCH_BOUND_MIN,
        maximum_bound=INCLINATION_SEARCH_BOUND_MAX,
        fallback_maximum=INCLINATION_SEARCH_BOUND_MAX,
    )


def sanitize_retry_search_bounds(bounds):
    return sanitize_inclination_search_bounds(sanitize_ars_search_bounds(sanitize_rprs_search_bounds(bounds)))


def impact_parameter_scale_for_retry(values):
    try:
        ars = float(values['ars'])
    except (KeyError, TypeError, ValueError):
        return np.nan

    try:
        ecc = float(values.get('ecc', 0.0))
    except (TypeError, ValueError):
        ecc = 0.0

    try:
        omega = np.deg2rad(float(values.get('omega', 0.0)))
    except (TypeError, ValueError):
        omega = 0.0

    denom = 1.0 + ecc * np.sin(omega)
    if np.isclose(denom, 0.0):
        denom = np.finfo(float).eps
    scale = ars * (1.0 - ecc ** 2) / denom
    return float(scale) if np.isfinite(scale) and scale > 0 else np.nan


def impact_parameter_scale_range_for_retry(prior, bounds):
    values = dict(prior)
    ars_candidates = []
    if 'ars' in bounds:
        try:
            ars_candidates.extend(
                float(value) for value in np.asarray(bounds['ars'], dtype=float).reshape(-1)[:2]
            )
        except (TypeError, ValueError, IndexError):
            pass
    if 'ars' in values:
        try:
            ars_candidates.append(float(values['ars']))
        except (TypeError, ValueError):
            pass

    scales = []
    for ars_value in ars_candidates:
        candidate_values = dict(values)
        candidate_values['ars'] = ars_value
        scale = impact_parameter_scale_for_retry(candidate_values)
        if np.isfinite(scale) and scale > 0:
            scales.append(scale)

    if not scales:
        scale = impact_parameter_scale_for_retry(values)
        if np.isfinite(scale) and scale > 0:
            scales.append(scale)

    if not scales:
        return np.nan, np.nan
    return float(np.nanmin(scales)), float(np.nanmax(scales))


def inclination_from_impact_parameter_for_retry(impact_parameter, scale):
    if not np.isfinite(scale) or scale <= 0:
        return np.nan
    try:
        impact_parameter = float(impact_parameter)
    except (TypeError, ValueError):
        return np.nan
    cosi = np.clip(impact_parameter / scale, -1.0, 1.0)
    return float(np.rad2deg(np.arccos(cosi)))


def impact_parameter_retry_proposed_inclination_bounds(diagnostics, current_prior, current_bounds):
    if not diagnostics:
        return None
    if 'inc' not in current_bounds:
        return None

    try:
        previous_lower, previous_upper = [
            float(value) for value in np.asarray(current_bounds['inc'], dtype=float).reshape(-1)[:2]
        ]
        proposed_b_lower, proposed_b_upper = [
            float(value) for value in np.asarray(diagnostics.get('bounds'), dtype=float).reshape(-1)[:2]
        ]
    except (TypeError, ValueError, IndexError):
        return None

    if (
        not np.isfinite(previous_lower)
        or not np.isfinite(previous_upper)
        or previous_lower >= previous_upper
        or not np.isfinite(proposed_b_lower)
        or not np.isfinite(proposed_b_upper)
        or proposed_b_lower >= proposed_b_upper
    ):
        return None

    min_scale, max_scale = impact_parameter_scale_range_for_retry(current_prior, current_bounds)
    if not np.isfinite(min_scale) or not np.isfinite(max_scale):
        return None

    new_lower = previous_lower
    new_upper = previous_upper
    clipped_edge = diagnostics.get('edge')

    if clipped_edge in ('upper', None):
        inc_for_upper_b = inclination_from_impact_parameter_for_retry(proposed_b_upper, max_scale)
        if np.isfinite(inc_for_upper_b):
            new_lower = min(new_lower, inc_for_upper_b)

    if clipped_edge in ('lower', None):
        inc_for_lower_b = inclination_from_impact_parameter_for_retry(max(0.0, proposed_b_lower), min_scale)
        if np.isfinite(inc_for_lower_b):
            new_upper = max(new_upper, inc_for_lower_b)

    return [float(new_lower), float(new_upper)]


def clamp_parameter_prior_to_bounds(prior, bounds, key):
    clamped = dict(prior)
    if key not in clamped or key not in bounds:
        return clamped

    try:
        parameter_value = float(clamped[key])
        lower_bound, upper_bound = [
            float(value) for value in np.asarray(bounds[key], dtype=float).reshape(-1)[:2]
        ]
    except (TypeError, ValueError, IndexError):
        return clamped

    if (
        np.isfinite(parameter_value) and np.isfinite(lower_bound) and np.isfinite(upper_bound)
        and lower_bound < upper_bound
    ):
        clamped[key] = float(np.clip(parameter_value, lower_bound, upper_bound))
    return clamped


def clamp_rprs_prior_to_bounds(prior, bounds):
    return clamp_parameter_prior_to_bounds(prior, bounds, 'rprs')


def clamp_ars_prior_to_bounds(prior, bounds):
    return clamp_parameter_prior_to_bounds(prior, bounds, 'ars')


def clamp_inclination_prior_to_bounds(prior, bounds):
    return clamp_parameter_prior_to_bounds(prior, bounds, 'inc')


def clamp_retry_priors_to_bounds(prior, bounds):
    return clamp_inclination_prior_to_bounds(
        clamp_ars_prior_to_bounds(clamp_rprs_prior_to_bounds(prior, bounds), bounds),
        bounds,
    )


def enforce_minimum_parameter_retry_half_width(
    mode,
    bounds,
    min_half_width,
    minimum_bound,
    maximum_bound=None,
):
    try:
        lower_bound, upper_bound = [
            float(value) for value in np.asarray(bounds, dtype=float).reshape(-1)[:2]
        ]
    except (TypeError, ValueError, IndexError):
        return bounds

    if not np.isfinite(lower_bound) or not np.isfinite(upper_bound) or lower_bound >= upper_bound:
        return bounds

    center = float(mode) if np.isfinite(mode) else float(0.5 * (lower_bound + upper_bound))
    half_width = max(float(min_half_width), 0.0)
    expanded_lower = min(lower_bound, center - half_width)
    expanded_upper = max(upper_bound, center + half_width)

    if expanded_lower < minimum_bound:
        if maximum_bound is None:
            expanded_upper = expanded_upper + (minimum_bound - expanded_lower)
        else:
            expanded_upper = min(
                maximum_bound,
                expanded_upper + (minimum_bound - expanded_lower),
            )
        expanded_lower = minimum_bound
    if maximum_bound is not None and expanded_upper > maximum_bound:
        expanded_lower = max(
            minimum_bound,
            expanded_lower - (expanded_upper - maximum_bound),
        )
        expanded_upper = maximum_bound

    return [float(expanded_lower), float(expanded_upper)]


def enforce_minimum_rprs_retry_half_width(mode, bounds, min_half_width=RPRS_RETRY_MIN_HALF_WIDTH):
    return enforce_minimum_parameter_retry_half_width(
        mode,
        bounds,
        min_half_width,
        RPRS_SEARCH_BOUND_MIN,
        maximum_bound=RPRS_SEARCH_BOUND_MAX,
    )


def enforce_minimum_ars_retry_half_width(mode, bounds, min_half_width=ARS_RETRY_MIN_HALF_WIDTH):
    return enforce_minimum_parameter_retry_half_width(
        mode,
        bounds,
        min_half_width,
        ARS_SEARCH_BOUND_MIN,
    )


def run_nested_lightcurve_fit_with_rprs_posterior_retry(
    times,
    flux_values,
    flux_errors,
    airmass,
    prior,
    bounds,
    jd_times=None,
    use_impactparameter_rather_than_inclination_to_fit=True,
    max_rprs_retries=RPRS_POSTERIOR_MAX_RETRIES_DEFAULT,
    duration_prior=None,
    max_ars_retries=ARS_POSTERIOR_MAX_RETRIES_DEFAULT,
    max_impact_parameter_retries=IMPACT_PARAMETER_POSTERIOR_MAX_RETRIES_DEFAULT,
    keep_ultranest_sampler=False,
    baseline_fit_mask=None,
    fixed_parameter_errors=None,
    fixed_flux_baseline=False,
    ultranest_min_num_live_points=None,
    pre_ultranest_coverage_assessment=None,
):
    def impact_parameter_retry_available(fit, local_bounds):
        if not use_impactparameter_rather_than_inclination_to_fit or 'inc' not in local_bounds:
            return False
        if getattr(fit, 'impact_parameter_sampled_directly', False):
            return False
        sampled_keys = getattr(fit, 'sampled_keys', []) or []
        sample_bounds = getattr(fit, 'sample_bounds', {})
        return 'b' in sampled_keys or (isinstance(sample_bounds, dict) and 'b' in sample_bounds)

    def identity_retry_bounds(diagnostics, local_prior, local_bounds, config):
        return diagnostics.get('bounds') if diagnostics else None

    def impact_parameter_retry_bounds(diagnostics, local_prior, local_bounds, config):
        return impact_parameter_retry_proposed_inclination_bounds(
            diagnostics,
            local_prior,
            local_bounds,
        )

    def normal_retry_expands(previous_bounds, new_bounds, clipped_edge, config):
        previous_lower, previous_upper = [
            float(value) for value in np.asarray(previous_bounds, dtype=float).reshape(-1)[:2]
        ]
        new_lower, new_upper = [
            float(value) for value in np.asarray(new_bounds, dtype=float).reshape(-1)[:2]
        ]
        if clipped_edge == 'upper':
            return new_upper > previous_upper + 1e-12
        if clipped_edge == 'lower':
            return new_lower < previous_lower - 1e-12
        return new_lower < previous_lower - 1e-12 or new_upper > previous_upper + 1e-12

    def impact_parameter_retry_expands(previous_bounds, new_bounds, clipped_edge, config):
        previous_lower, previous_upper = [
            float(value) for value in np.asarray(previous_bounds, dtype=float).reshape(-1)[:2]
        ]
        new_lower, new_upper = [
            float(value) for value in np.asarray(new_bounds, dtype=float).reshape(-1)[:2]
        ]
        if clipped_edge == 'upper':
            return new_lower < previous_lower - 1e-12
        if clipped_edge == 'lower':
            return new_upper > previous_upper + 1e-12
        return new_lower < previous_lower - 1e-12 or new_upper > previous_upper + 1e-12

    partial_retry_limits = partial_transit_geometry_retry_limits(pre_ultranest_coverage_assessment)
    retry_configs = [
        {
            'key': 'rprs',
            'diagnostic_key': 'rprs',
            'bounds_key': 'rprs',
            'label': 'Rp/R*',
            'sanitize_bounds': sanitize_rprs_search_bounds,
            'enforce_half_width': enforce_minimum_rprs_retry_half_width,
            'propose_bounds': identity_retry_bounds,
            'expands_bounds': normal_retry_expands,
            'max_retries': min(
                max_rprs_retries,
                partial_retry_limits['max_retries']['rprs'],
            ) if partial_retry_limits['active'] else max_rprs_retries,
            'min_bound': RPRS_SEARCH_BOUND_MIN,
            'max_bound': RPRS_SEARCH_BOUND_MAX,
            'prior_mode_key': 'rprs',
            'annotate': annotate_rprs_posterior_refit,
        },
        {
            'key': 'ars',
            'diagnostic_key': 'ars',
            'bounds_key': 'ars',
            'label': 'a/Rs',
            'sanitize_bounds': sanitize_ars_search_bounds,
            'enforce_half_width': enforce_minimum_ars_retry_half_width,
            'propose_bounds': identity_retry_bounds,
            'expands_bounds': normal_retry_expands,
            'max_retries': min(
                max_ars_retries,
                partial_retry_limits['max_retries']['ars'],
            ) if partial_retry_limits['active'] else max_ars_retries,
            'min_bound': ARS_SEARCH_BOUND_MIN,
            'max_bound': None,
            'prior_mode_key': 'ars',
            'annotate': annotate_ars_posterior_refit,
        },
        {
            'key': 'b',
            'diagnostic_key': 'b',
            'bounds_key': 'inc',
            'label': 'impact parameter',
            'sanitize_bounds': sanitize_inclination_search_bounds,
            'enforce_half_width': lambda mode, bounds: bounds,
            'propose_bounds': impact_parameter_retry_bounds,
            'expands_bounds': impact_parameter_retry_expands,
            'max_retries': min(
                max_impact_parameter_retries,
                partial_retry_limits['max_retries']['b'],
            ) if partial_retry_limits['active'] else max_impact_parameter_retries,
            'min_bound': INCLINATION_SEARCH_BOUND_MIN,
            'max_bound': INCLINATION_SEARCH_BOUND_MAX,
            'prior_mode_key': None,
            'available': impact_parameter_retry_available,
            'annotate': annotate_impact_parameter_posterior_refit,
        },
    ]

    def build_fit(local_prior, local_bounds):
        local_bounds = sanitize_retry_search_bounds(local_bounds)
        local_prior = clamp_retry_priors_to_bounds(local_prior, local_bounds)
        fit_kwargs = {
            'jd_times': jd_times,
            'mode': 'ns',
            'use_impactparameter_rather_than_inclination_to_fit':
            use_impactparameter_rather_than_inclination_to_fit,
        }
        if isinstance(duration_prior, dict) and duration_prior.get('applied'):
            fit_kwargs['duration_prior'] = duration_prior
        if keep_ultranest_sampler and callable_accepts_keyword(lc_fitter, 'keep_ultranest_sampler'):
            fit_kwargs['keep_ultranest_sampler'] = True
        if baseline_fit_mask is not None and callable_accepts_keyword(lc_fitter, 'baseline_fit_mask'):
            fit_kwargs['baseline_fit_mask'] = baseline_fit_mask
        if fixed_parameter_errors and callable_accepts_keyword(lc_fitter, 'fixed_parameter_errors'):
            fit_kwargs['fixed_parameter_errors'] = fixed_parameter_errors
        if fixed_flux_baseline and callable_accepts_keyword(lc_fitter, 'fixed_flux_baseline'):
            fit_kwargs['fixed_flux_baseline'] = True
        if (
            ultranest_min_num_live_points is not None
            and callable_accepts_keyword(lc_fitter, 'ultranest_min_num_live_points')
        ):
            fit_kwargs['ultranest_min_num_live_points'] = ultranest_min_num_live_points
        fit = lc_fitter(
            times,
            flux_values,
            flux_errors,
            airmass,
            local_prior,
            local_bounds,
            **fit_kwargs,
        )
        annotate_duration_prior(fit, duration_prior)
        return fit

    current_bounds = sanitize_retry_search_bounds(bounds)
    current_prior = clamp_retry_priors_to_bounds(prior, current_bounds)
    retry_histories = {config['key']: [] for config in retry_configs}
    retry_notes = {config['key']: None for config in retry_configs}
    latest_diagnostics = {config['key']: None for config in retry_configs}
    blocked_retry_keys = set()
    fit = build_fit(current_prior, current_bounds)

    while True:
        diagnostics_getter = getattr(fit, "get_parameter_posterior_recenter_diagnostics", None)
        if not callable(diagnostics_getter):
            latest_diagnostics = {config['key']: None for config in retry_configs}
            break

        retry_config = None
        diagnostics = None
        for config in retry_configs:
            key = config['key']
            diagnostic_key = config.get('diagnostic_key', key)
            bounds_key = config.get('bounds_key', key)
            if bounds_key not in current_bounds:
                continue

            available = config.get('available')
            if callable(available) and not available(fit, current_bounds):
                continue

            parameter_diagnostics = diagnostics_getter(diagnostic_key)
            latest_diagnostics[key] = parameter_diagnostics
            if key in blocked_retry_keys:
                continue

            new_bounds = parameter_diagnostics.get('bounds') if parameter_diagnostics else None
            if len(retry_histories[key]) >= int(max(0, config['max_retries'])):
                if (
                    partial_retry_limits['active']
                    and parameter_diagnostics
                    and parameter_diagnostics.get('clipped')
                    and retry_notes[key] is None
                ):
                    retry_notes[key] = partial_retry_limits['note']
                continue

            if parameter_diagnostics and parameter_diagnostics.get('clipped') and new_bounds is not None:
                retry_config = config
                diagnostics = parameter_diagnostics
                break

        if retry_config is None:
            break

        key = retry_config['key']
        bounds_key = retry_config.get('bounds_key', key)
        label = retry_config['label']
        new_bounds = retry_config.get('propose_bounds', identity_retry_bounds)(
            diagnostics,
            current_prior,
            current_bounds,
            retry_config,
        )
        try:
            new_lower, new_upper = [float(value) for value in new_bounds]
        except (TypeError, ValueError):
            retry_notes[key] = f"Skipped; the automatic {label} retry proposed malformed bounds."
            blocked_retry_keys.add(key)
            continue
        if not np.isfinite(new_lower) or not np.isfinite(new_upper) or new_lower >= new_upper:
            retry_notes[key] = f"Skipped; the automatic {label} retry proposed invalid bounds."
            blocked_retry_keys.add(key)
            continue

        previous_bounds = current_bounds.get(bounds_key)
        clamped_bounds = retry_config['sanitize_bounds']({bounds_key: [new_lower, new_upper]}).get(
            bounds_key,
            [new_lower, new_upper],
        )
        clamped_bounds = retry_config['enforce_half_width'](
            diagnostics.get('mode', np.nan),
            clamped_bounds,
        )
        clamped_bounds = retry_config['sanitize_bounds']({bounds_key: clamped_bounds}).get(bounds_key, clamped_bounds)
        new_lower, new_upper = [float(value) for value in clamped_bounds]
        if previous_bounds is not None:
            previous_lower, previous_upper = [float(value) for value in np.asarray(previous_bounds, dtype=float).reshape(-1)[:2]]
            clipped_edge = diagnostics.get('edge')
            expands_sampled_range = retry_config.get('expands_bounds', normal_retry_expands)(
                previous_bounds,
                [new_lower, new_upper],
                clipped_edge,
                retry_config,
            )

            if not expands_sampled_range:
                maximum_bound = retry_config['max_bound']
                if (
                    maximum_bound is not None and
                    previous_lower <= retry_config['min_bound'] + 1e-12 and
                    previous_upper >= maximum_bound - 1e-12
                ):
                    retry_notes[key] = (
                        f"Skipped; the automatic {label} retry reached the maximum exoplanet "
                        f"search range [{retry_config['min_bound']:.6f}, {maximum_bound:.6f}]."
                    )
                else:
                    retry_notes[key] = f"Skipped; the automatic {label} retry did not expand the sampled range."
                blocked_retry_keys.add(key)
                continue

        retry_histories[key].append({
            'attempt': len(retry_histories[key]) + 1,
            'edge': diagnostics.get('edge'),
            'mode': float(diagnostics.get('mode', np.nan)),
            'std': float(diagnostics.get('std', np.nan)),
            'original_bounds': None if previous_bounds is None else [float(previous_bounds[0]), float(previous_bounds[1])],
            'new_bounds': [new_lower, new_upper],
        })
        log_info(
            f"{label} posterior is truncated against the "
            f"{diagnostics.get('edge', 'active')} search bound; retrying nested fit "
            f"with {label} centered at {diagnostics.get('mode', np.nan):.6f} "
            f"and sigma {diagnostics.get('std', np.nan):.6f} "
            f"over [{new_lower:.6f}, {new_upper:.6f}]."
        )

        updated_bounds = clone_lightcurve_bounds(current_bounds)
        updated_bounds[bounds_key] = [new_lower, new_upper]
        updated_bounds = sanitize_retry_search_bounds(updated_bounds)

        updated_prior = dict(current_prior)
        fit_parameters = getattr(fit, 'parameters', {})
        if isinstance(fit_parameters, dict):
            for bound_key in updated_bounds:
                if bound_key in fit_parameters:
                    updated_prior[bound_key] = fit_parameters[bound_key]
        prior_mode_key = retry_config.get('prior_mode_key', key)
        if prior_mode_key is not None and np.isfinite(diagnostics.get('mode', np.nan)):
            updated_prior[prior_mode_key] = float(diagnostics['mode'])
        updated_prior = clamp_retry_priors_to_bounds(updated_prior, updated_bounds)

        current_prior = updated_prior
        current_bounds = updated_bounds
        fit = build_fit(current_prior, current_bounds)

    final_diagnostics_getter = getattr(fit, "get_parameter_posterior_recenter_diagnostics", None)
    annotate_posterior_refit_final_bounds(fit, current_bounds)
    for config in retry_configs:
        key = config['key']
        diagnostic_key = config.get('diagnostic_key', key)
        bounds_key = config.get('bounds_key', key)
        label = config['label']
        history = retry_histories[key]
        final_diagnostics = None
        available = config.get('available')
        config_available = not callable(available) or available(fit, current_bounds)
        if callable(final_diagnostics_getter) and bounds_key in current_bounds and config_available:
            final_diagnostics = final_diagnostics_getter(diagnostic_key)
        elif latest_diagnostics.get(key) is not None:
            final_diagnostics = latest_diagnostics[key]

        if history:
            note = f"Applied {len(history)} automatic {label} posterior range refit(s)."
            if final_diagnostics and final_diagnostics.get('clipped'):
                retry_label = "retry" if len(history) == 1 else "retries"
                note = (
                    f"{note} The posterior still hugs the {final_diagnostics.get('edge')} bound after "
                    f"{len(history)} {retry_label}."
                )
                log_info(
                    f"Warning: {label} posterior still appears truncated after the automatic retries; "
                    "please inspect the triangle plot carefully.",
                    warn=True,
                )
        elif retry_notes[key] is not None:
            note = retry_notes[key]
        elif final_diagnostics is not None and final_diagnostics.get('reason'):
            note = f"Not needed; {final_diagnostics['reason']}"
        else:
            note = "Not evaluated; posterior diagnostics are unavailable for this fit."

        config['annotate'](fit, bool(history), note=note, history=history)
    return fit


def log_info(string, warn=False, error=False):
    if error:
        print(f"\033[31m {string}\033[0m", flush=True)
    elif warn:
        print(f"\033[34m {string}\033[0m", flush=True)
    else:
        print(string, flush=True)
    log.debug(string)
    _reset_runtime_traceback_watchdog()
    return True


def _find_runtime_handler(handler_name):
    for handler in log.handlers:
        if getattr(handler, "_exotic_runtime_handler_name", None) == handler_name:
            return handler
    return None


def _runtime_traceback_watchdog_seconds():
    try:
        return float(os.environ.get(
            _RUNTIME_TRACEBACK_WATCHDOG_SECONDS_ENV,
            _RUNTIME_TRACEBACK_WATCHDOG_DEFAULT_SECONDS,
        ))
    except (TypeError, ValueError):
        return _RUNTIME_TRACEBACK_WATCHDOG_DEFAULT_SECONDS


def _reset_runtime_traceback_watchdog():
    global _RUNTIME_TRACEBACK_WATCHDOG_ACTIVE

    if not _RUNTIME_LOGGING_CONFIGURED:
        return

    timeout = _runtime_traceback_watchdog_seconds()
    if timeout <= 0:
        cancel_runtime_traceback_watchdog()
        return

    try:
        faulthandler.cancel_dump_traceback_later()
    except Exception:
        pass

    try:
        faulthandler.dump_traceback_later(timeout, repeat=False, file=sys.stdout)
        _RUNTIME_TRACEBACK_WATCHDOG_ACTIVE = True
    except Exception:
        _RUNTIME_TRACEBACK_WATCHDOG_ACTIVE = False


def cancel_runtime_traceback_watchdog():
    global _RUNTIME_TRACEBACK_WATCHDOG_ACTIVE

    if not _RUNTIME_TRACEBACK_WATCHDOG_ACTIVE:
        return

    try:
        faulthandler.cancel_dump_traceback_later()
    except Exception:
        pass
    _RUNTIME_TRACEBACK_WATCHDOG_ACTIVE = False


def configure_runtime_logging():
    global _RUNTIME_LOGGING_CONFIGURED

    logging.root.setLevel(logging.DEBUG)
    log.setLevel(logging.DEBUG)

    if _find_runtime_handler(_RUNTIME_FILE_HANDLER_NAME) is None:
        try:
            file_handler = TimedRotatingFileHandler(filename="exotic.log", when="midnight", backupCount=2)
        except Exception as exc:
            print(f"Warning: Could not initialize exotic.log ({exc}).")
        else:
            file_handler._exotic_runtime_handler_name = _RUNTIME_FILE_HANDLER_NAME
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(
                logging.Formatter(
                    "%(asctime)s.%(msecs)03d [%(threadName)-12.12s] %(levelname)-5.5s  "
                    "%(funcName)s:%(lineno)d - %(message)s",
                    "%Y-%m-%dT%H:%M:%S",
                )
            )
            log.addHandler(file_handler)

    console_handler = _find_runtime_handler(_RUNTIME_CONSOLE_HANDLER_NAME)
    if console_handler is None:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler._exotic_runtime_handler_name = _RUNTIME_CONSOLE_HANDLER_NAME
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(logging.Formatter("%(message)s"))
        log.addHandler(console_handler)
    else:
        try:
            console_handler.setStream(sys.stdout)
        except Exception:
            console_handler.stream = sys.stdout

    try:
        faulthandler.enable(file=sys.stdout, all_threads=True)
    except Exception:
        pass

    _RUNTIME_LOGGING_CONFIGURED = True
    _reset_runtime_traceback_watchdog()


def _logger_has_current_stdout_handler(logger):
    current_stdout = sys.stdout
    active_logger = logger
    while active_logger:
        for handler in active_logger.handlers:
            if getattr(handler, "stream", None) is current_stdout:
                return True
        if not getattr(active_logger, "propagate", False):
            break
        active_logger = active_logger.parent
    return False


def _write_exception_traceback_to_stdout(message, exc_type, exc_value, exc_traceback):
    traceback_text = ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))
    try:
        print(f"\n{message}", file=sys.stdout, flush=True)
        print(traceback_text, file=sys.stdout, end="", flush=True)
    except Exception:
        try:
            print(f"\n{message}", file=sys.__stdout__, flush=True)
            print(traceback_text, file=sys.__stdout__, end="", flush=True)
        except Exception:
            pass


def _log_exception_with_fallback(message, exc_type, exc_value, exc_traceback):
    wrote_to_logger = False
    try:
        log.error(message, exc_info=(exc_type, exc_value, exc_traceback))
        wrote_to_logger = True
    except Exception:
        pass

    if not wrote_to_logger or not _logger_has_current_stdout_handler(log):
        _write_exception_traceback_to_stdout(message, exc_type, exc_value, exc_traceback)


def _handle_unhandled_exception(exc_type, exc_value, exc_traceback):
    global _UNHANDLED_EXCEPTION_LOGGED

    if exc_type is not None and issubclass(exc_type, KeyboardInterrupt):
        return

    if _UNHANDLED_EXCEPTION_LOGGED:
        return

    _UNHANDLED_EXCEPTION_LOGGED = True
    _log_exception_with_fallback("Unhandled exception during EXOTIC run", exc_type, exc_value, exc_traceback)


def _handle_thread_exception(args):
    if args.exc_type is not None and issubclass(args.exc_type, KeyboardInterrupt):
        return

    thread_name = args.thread.name if args.thread is not None else "unknown"
    _log_exception_with_fallback(
        f"Unhandled exception in thread '{thread_name}'",
        args.exc_type,
        args.exc_value,
        args.exc_traceback,
    )


def install_exception_hooks():
    global _EXCEPTION_HOOKS_INSTALLED

    if _EXCEPTION_HOOKS_INSTALLED:
        return

    sys.excepthook = _handle_unhandled_exception
    threading.excepthook = _handle_thread_exception
    _EXCEPTION_HOOKS_INSTALLED = True


def should_log_plate_solution_path(wcs_file):
    if not wcs_file:
        return False

    normalized_path = os.fspath(wcs_file).replace("\\", "/")
    return normalized_path != "/tmp" and not normalized_path.startswith("/tmp/")


def log_mid_transit_range_warning_once(array_times, tmid_prior):
    global _mid_transit_warning_reported
    if _mid_transit_warning_reported:
        return

    _mid_transit_warning_reported = True
    # Keep this warning in plain black text and show it once per run.
    log_info("\nWarning:")
    log_info(" Estimated mid-transit time is not within the observations")
    log_info(" Check Period & Mid-transit time in inits.json. Make sure the uncertainties are not 0 or Nan.")
    log_info(f"  obs start:{array_times.min()}")
    log_info(f"    obs end:{array_times.max()}")
    log_info(f" tmid prior:{tmid_prior}\n")


def relative_flux_filter_mask(relative_flux, max_relative_flux=RELATIVE_FLUX_MAX):
    relative_flux = np.asarray(relative_flux, dtype=float)
    return (
        np.isfinite(relative_flux)
        & np.greater(relative_flux, 0)
    )


def valid_flux_ratio_mask(relative_flux):
    relative_flux = np.asarray(relative_flux, dtype=float)
    return np.isfinite(relative_flux) & np.greater(relative_flux, 0)


def valid_comparison_frame_mask(flux_values):
    flux_values = np.asarray(flux_values, dtype=float)
    return np.isfinite(flux_values) & (flux_values > 0)


def robust_flux_floor_mask(
    flux_values,
    min_fraction_of_median=ROBUST_FLUX_MIN_FRACTION_OF_MEDIAN,
    min_points=ROBUST_FLUX_MIN_POINTS,
):
    flux_values = np.asarray(flux_values, dtype=float)
    valid = np.isfinite(flux_values) & (flux_values > 0)
    if np.count_nonzero(valid) < max(LIGHTCURVE_MIN_VALID_POINTS, int(min_points)):
        return valid

    center, _ = sigma_clipped_nanmedian(flux_values[valid], sigma=4.0, max_iters=3)
    if not np.isfinite(center) or center <= 0:
        center = bn.nanmedian(flux_values[valid])
    if not np.isfinite(center) or center <= 0:
        return valid

    floor = float(min_fraction_of_median) * float(center)
    if not np.isfinite(floor) or floor <= 0:
        return valid

    return valid & np.greater_equal(flux_values, floor)


def robust_target_reference_flux_mask(target_flux, reference_flux):
    target_mask = robust_flux_floor_mask(target_flux)
    if reference_flux is None:
        return target_mask

    reference_mask = robust_flux_floor_mask(reference_flux)
    return target_mask & reference_mask


def is_fast_aperture_mask_enabled(config_value):
    if config_value is None:
        return False
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'fast', 'center', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'exact', 'off'):
            return False

    log_info("Warning: Invalid 'Fast Aperture Mask (y/n)' value; using exact mode.", warn=True)
    return False


def is_comp_star_required(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off'):
            return False

    log_info("Warning: Invalid 'require_comp_star' value; requiring a comparison star.", warn=True)
    return True


def is_target_driven_comp_selection_enabled(config_value):
    if config_value is None:
        return False
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off'):
            return False

    log_info("Warning: Invalid target-driven comparison selection value; using comp-driven selection.", warn=True)
    return False


def should_skip_low_comparison_coverage_rejection(config_value):
    if config_value is None:
        return False
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info("Warning: Invalid 'skip_low_comparison_coverage_rejection' value; keeping coverage rejection enabled.",
             warn=True)
    return False


def should_fit_lightcurve_to_every_comparison_candidate(config_value):
    if config_value is None:
        return False
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info("Warning: Invalid 'fit_lightcurve_to_every_comparison_candidate' value; defaulting to disabled.",
             warn=True)
    return False


def should_use_sparse_posterior_live_point_retry(config_value):
    if config_value is None:
        return SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_DEFAULT
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'use_sparse_posterior_live_point_retry' value; "
        "defaulting to enabled.",
        warn=True,
    )
    return SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_DEFAULT


def should_run_fast_ultranest_before_final_run(config_value):
    if config_value is None:
        return FAST_ULTRANEST_BEFORE_FINAL_RUN_DEFAULT
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'run fast ultranest before final run' value; "
        "defaulting to enabled.",
        warn=True,
    )
    return FAST_ULTRANEST_BEFORE_FINAL_RUN_DEFAULT


def configure_sparse_posterior_live_point_retry(config_value):
    enabled = should_use_sparse_posterior_live_point_retry(config_value)
    os.environ[SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_ENV] = "1" if enabled else "0"
    return enabled


def should_pick_comparison_by_eebls_snr(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'pick_comparison_by_eebls_snr' value; "
        "keeping EEBLS SNR comparison selection enabled.",
        warn=True,
    )
    return True


def should_use_deviation_from_expected_transit_in_qc(config_value):
    if config_value is None:
        return TRANSIT_QC_USE_DEVIATION_FROM_EXPECTED_DEFAULT
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'use_deviation_from_expected_transit_in_qc' value; "
        "keeping expected-value deviation QC enabled.",
        warn=True,
    )
    return TRANSIT_QC_USE_DEVIATION_FROM_EXPECTED_DEFAULT


def parse_deviation_from_expected_transit_in_qc_sigma(config_value):
    if config_value is None:
        return TRANSIT_QC_DEVIATION_SIGMA_DEFAULT

    try:
        sigma_value = float(config_value)
    except (TypeError, ValueError):
        log_info(
            "Warning: Invalid 'deviation_from_expected_transit_in_qc_sigma' value; using the default 5 sigma.",
            warn=True,
        )
        return TRANSIT_QC_DEVIATION_SIGMA_DEFAULT

    if not np.isfinite(sigma_value) or sigma_value <= 0:
        log_info(
            "Warning: Non-positive 'deviation_from_expected_transit_in_qc_sigma' value; using the default 5 sigma.",
            warn=True,
        )
        return TRANSIT_QC_DEVIATION_SIGMA_DEFAULT
    return float(sigma_value)


def should_assess_all_comparisons_before_selecting_best(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'assess_all_comparisons_before_selecting_best' value; "
        "defaulting to assess all comparisons.",
        warn=True,
    )
    return True


def should_exit_at_first_qc_pass_solution(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'exit_at_first_qc_pass_solution' value; "
        "defaulting to exit at the first QC PASS solution.",
        warn=True,
    )
    return True


def parse_ultranest_min_num_live_points(config_value):
    if config_value is None:
        return ULTRANEST_MIN_NUM_LIVE_POINTS_DEFAULT

    if isinstance(config_value, str) and config_value.strip() == "":
        return ULTRANEST_MIN_NUM_LIVE_POINTS_DEFAULT

    try:
        live_points = int(float(str(config_value).strip()))
    except (TypeError, ValueError):
        log_info(
            "Warning: Invalid 'minimum number of live points for ultranest' value; "
            f"defaulting to {ULTRANEST_MIN_NUM_LIVE_POINTS_DEFAULT}.",
            warn=True,
        )
        return ULTRANEST_MIN_NUM_LIVE_POINTS_DEFAULT

    if live_points <= 0:
        log_info(
            "Warning: 'minimum number of live points for ultranest' must be positive; "
            f"defaulting to {ULTRANEST_MIN_NUM_LIVE_POINTS_DEFAULT}.",
            warn=True,
        )
        return ULTRANEST_MIN_NUM_LIVE_POINTS_DEFAULT

    return live_points


def configure_ultranest_min_num_live_points(config_value):
    live_points = parse_ultranest_min_num_live_points(config_value)
    os.environ[ULTRANEST_MIN_NUM_LIVE_POINTS_ENV] = str(live_points)
    return live_points


def parse_rprs_search_bound_max(config_value):
    if config_value is None:
        return RPRS_SEARCH_BOUND_MAX_DEFAULT

    if isinstance(config_value, str) and config_value.strip() == "":
        return RPRS_SEARCH_BOUND_MAX_DEFAULT

    try:
        max_bound = float(str(config_value).strip())
    except (TypeError, ValueError):
        log_info(
            "Warning: Invalid 'rprs_search_bound_max' value; "
            f"defaulting to {RPRS_SEARCH_BOUND_MAX_DEFAULT:.3f}.",
            warn=True,
        )
        return RPRS_SEARCH_BOUND_MAX_DEFAULT

    if not np.isfinite(max_bound) or max_bound <= RPRS_SEARCH_BOUND_MIN:
        log_info(
            "Warning: 'rprs_search_bound_max' must be finite and positive; "
            f"defaulting to {RPRS_SEARCH_BOUND_MAX_DEFAULT:.3f}.",
            warn=True,
        )
        return RPRS_SEARCH_BOUND_MAX_DEFAULT

    if max_bound > RPRS_SEARCH_BOUND_ABSOLUTE_MAX:
        log_info(
            "Warning: 'rprs_search_bound_max' exceeds the absolute safety ceiling "
            f"of {RPRS_SEARCH_BOUND_ABSOLUTE_MAX:.3f}; clamping to that ceiling.",
            warn=True,
        )
        return RPRS_SEARCH_BOUND_ABSOLUTE_MAX

    return float(max_bound)


def configure_rprs_search_bound_max(config_value):
    global RPRS_SEARCH_BOUND_MAX
    RPRS_SEARCH_BOUND_MAX = parse_rprs_search_bound_max(config_value)
    return RPRS_SEARCH_BOUND_MAX


def log_ultranest_mpi_status():
    status = get_mpi_status()
    size = int(status.get("size") or 1)
    rank = int(status.get("rank") or 0)
    if size <= 1 or rank != 0:
        return status

    if status.get("available"):
        log_info(f"UltraNest MPI mode detected: {size} process(es).")
    else:
        log_info(
            "Warning: MPI launch detected, but mpi4py is unavailable; "
            "UltraNest cannot coordinate MPI workers until mpi4py is installed.",
            warn=True,
        )
    return status


def validate_ultranest_mpi_runtime():
    status = get_mpi_status()
    size = int(status.get("size") or 1)
    if size <= 1:
        return status

    message = (
        "EXOTIC was launched under MPI, which duplicates the full reduction on every rank. "
        "Start EXOTIC once and set EXOTIC_ULTRANEST_WORKERS to control UltraNest CPU parallelism."
    )
    if int(status.get("rank") or 0) == 0:
        log_info(f"Error: {message}", error=True)
    raise RuntimeError(message)


def configure_windows_multiprocessing_main_spec():
    if sys.platform != "win32":
        return False

    configured = False
    spawn_executable = _windows_python_spawn_executable()
    if spawn_executable:
        multiprocessing.set_executable(spawn_executable)
        if getattr(sys, "frozen", False):
            sys.frozen = False
        configured = True

    main_module = sys.modules.get("__main__")
    if main_module is None:
        return configured

    main_file = getattr(main_module, "__file__", None)
    if not main_file or os.path.basename(os.fspath(main_file)).lower() not in {"exotic.exe", "exotic-script.py"}:
        return configured

    if getattr(main_module, "__spec__", None) is not None:
        main_module.__spec__ = None
    main_module.__file__ = None
    if getattr(main_module, "__package__", None) is not None:
        main_module.__package__ = None

    return True


def ProcessPoolExecutor(*args, **kwargs):
    if sys.platform == "win32":
        return ThreadPoolExecutor(*args, **kwargs)
    return _ProcessPoolExecutor(*args, **kwargs)


def _windows_python_spawn_executable():
    candidates = [
        getattr(sys, "_base_executable", None),
        sys.executable,
        os.path.join(sys.exec_prefix, "python.exe"),
        os.path.join(getattr(sys, "base_exec_prefix", sys.exec_prefix), "python.exe"),
    ]
    for candidate in candidates:
        if not candidate:
            continue
        executable = os.fspath(candidate)
        if os.path.basename(executable).lower() in {"python.exe", "pythonw.exe"}:
            return executable
    return None


def should_use_psf_photometry(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info("Warning: Invalid 'use_psf_photometry' value; keeping PSF photometry enabled.", warn=True)
    return True


def should_use_aperture_photometry(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info("Warning: Invalid 'use_aperture_photometry' value; keeping aperture photometry enabled.", warn=True)
    return True


def should_use_eebls_to_initialize_tmid_and_bounds(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'use_eebls_to_initialize_tmid_and_bounds' value; "
        "keeping the EEBLS transit initializer enabled.",
        warn=True,
    )
    return True


def should_detect_bad_pixels_before_photometry(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'detect_bad_pixels_before_photometry' value; keeping bad-pixel precheck enabled.",
        warn=True,
    )
    return True


def get_multiprocess_bad_pixel_precheck_processes(config_value):
    if config_value is None:
        return None
    if isinstance(config_value, bool):
        if not config_value:
            return None
        return os.cpu_count() or 1
    if isinstance(config_value, (int, float)):
        if np.isfinite(config_value) and int(config_value) > 0:
            return int(config_value)
        return None
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('', 'n', 'no', 'false', '0', 'off'):
            return None
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return os.cpu_count() or 1
        try:
            parsed = float(normalized)
        except ValueError:
            parsed = np.nan
        if np.isfinite(parsed) and int(parsed) > 0:
            return int(parsed)

    log_info(
        "Warning: Invalid 'multiprocess_bad_pixel_precheck' value; keeping bad-pixel precheck multiprocessing disabled.",
        warn=True,
    )
    return None


def is_adaptive_aperture_mode_enabled(config_value):
    if config_value is None:
        return False
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info("Warning: Invalid 'use_adaptive_apertures' value; using fixed apertures.", warn=True)
    return False


def should_ignore_header_wcs(config_value):
    if config_value is None:
        return False
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info("Warning: Invalid 'Ignore WCS in Header and Do Manual Alignment? (y/n)' value; "
             "using header WCS when available.", warn=True)
    return False


def get_bad_wcs_threshold_fraction(config_value):
    default_fraction = SPARSE_MISSING_WCS_DROP_THRESHOLD
    default_percent = default_fraction * 100.0
    if config_value is None:
        return default_fraction

    if isinstance(config_value, str):
        normalized = config_value.strip()
        if normalized == "":
            return default_fraction
        if normalized.endswith('%'):
            normalized = normalized[:-1].strip()
    else:
        normalized = config_value

    try:
        threshold_percent = float(normalized)
    except (TypeError, ValueError):
        log_info(f"Warning: Invalid 'bad_wcs_threshold_percent' value; using default {default_percent:g}%.",
                 warn=True)
        return default_fraction

    if not np.isfinite(threshold_percent) or threshold_percent < 0 or threshold_percent > 100:
        log_info(f"Warning: Invalid 'bad_wcs_threshold_percent' value; using default {default_percent:g}%.",
                 warn=True)
        return default_fraction

    return threshold_percent / 100.0


def get_pointing_rejection_sigma(config_value):
    default_sigma = 4.0
    if config_value is None:
        return default_sigma

    if isinstance(config_value, str):
        normalized = config_value.strip()
        if normalized == "":
            return default_sigma
    else:
        normalized = config_value

    try:
        sigma = float(normalized)
    except (TypeError, ValueError):
        log_info(
            f"Warning: Invalid 'pointing_rejection_sigma' value; using default {default_sigma:g}.",
            warn=True,
        )
        return default_sigma

    if not np.isfinite(sigma):
        log_info(
            f"Warning: Invalid 'pointing_rejection_sigma' value; using default {default_sigma:g}.",
            warn=True,
        )
        return default_sigma

    if sigma < 0:
        log_info(
            f"Warning: Invalid 'pointing_rejection_sigma' value; using default {default_sigma:g}.",
            warn=True,
        )
        return default_sigma

    if sigma == 0:
        return None

    return sigma


def is_vertical_flux_normalization_disabled(config_value):
    if config_value is None:
        return False
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info("Warning: Invalid 'disable vertical flux normalization' value; using default enabled normalization.", warn=True)
    return False


def is_out_of_transit_baseline_detrending_enabled(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'detrend_on_outoftransit_baseline' value; using default enabled setting.",
        warn=True,
    )
    return True


def get_final_fit_baseline_duration_multiplier(config_value):
    if config_value is None:
        return FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT
    if isinstance(config_value, (int, float)) and np.isfinite(config_value) and config_value >= 0:
        return float(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip()
        if normalized == "":
            return FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT
        try:
            parsed = float(normalized)
        except ValueError:
            parsed = np.nan
        if np.isfinite(parsed) and parsed >= 0:
            return float(parsed)

    log_info(
        "Warning: Invalid 'final_fit_baseline_duration_multiplier' value; "
        f"using default {FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT:.1f}.",
        warn=True,
    )
    return FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT


def estimate_transit_duration_from_prior_geometry(prior):
    try:
        period = float(prior['per'])
        rprs = float(prior['rprs'])
        ars = float(prior['ars'])
        inc = float(prior['inc'])
    except (KeyError, TypeError, ValueError):
        return np.nan

    if (
        not np.isfinite(period) or period <= 0
        or not np.isfinite(rprs) or rprs < 0
        or not np.isfinite(ars) or ars <= 0
        or not np.isfinite(inc)
    ):
        return np.nan

    ecc = prior.get('ecc', 0.0)
    omega = np.deg2rad(prior.get('omega', 0.0))
    sin_inc = np.sin(np.deg2rad(inc))
    if not np.isfinite(sin_inc) or sin_inc <= 0:
        return np.nan

    impact_scale = ars * (1.0 - ecc ** 2) / max(np.finfo(float).eps, 1.0 + ecc * np.sin(omega))
    impact_parameter = impact_scale * np.cos(np.deg2rad(inc))
    chord_sq = (1.0 + rprs) ** 2 - impact_parameter ** 2
    if not np.isfinite(chord_sq) or chord_sq <= 0 or not np.isfinite(impact_scale) or impact_scale <= 0:
        return np.nan

    argument = np.sqrt(chord_sq) / (impact_scale * sin_inc)
    argument = float(np.clip(argument, -1.0, 1.0))
    duration = (period / np.pi) * np.arcsin(argument)
    return float(duration) if np.isfinite(duration) and duration > 0 else np.nan


def _cacheable_duration_prior_scalar(value):
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return None
    return None if not np.isfinite(numeric_value) else float(numeric_value)


def _duration_prior_scalar_from_cache(value):
    return np.nan if value is None else float(value)


@lru_cache(maxsize=128)
def _cached_single_transit_duration_prior(
    period,
    rprs,
    ars,
    inc,
    ecc,
    omega,
    period_unc,
    rprs_unc,
    ars_unc,
    inc_unc,
):
    prior = {
        'per': _duration_prior_scalar_from_cache(period),
        'rprs': _duration_prior_scalar_from_cache(rprs),
        'ars': _duration_prior_scalar_from_cache(ars),
        'inc': _duration_prior_scalar_from_cache(inc),
        'ecc': _duration_prior_scalar_from_cache(ecc),
        'omega': _duration_prior_scalar_from_cache(omega),
    }
    expected_duration = estimate_transit_duration_from_prior_geometry(prior)
    if not np.isfinite(expected_duration) or expected_duration <= 0:
        return {
            'applied': False,
            'expected_duration': np.nan,
            'sigma_log_duration': np.nan,
            'relative_sigma': np.nan,
            'source': 'unavailable',
            'sample_count': 0,
            'note': (
                "Not applied; could not estimate a physical transit duration from the published single-transit priors."
            ),
        }

    period_unc = _duration_prior_scalar_from_cache(period_unc)
    rprs_unc = _duration_prior_scalar_from_cache(rprs_unc)
    ars_unc = _duration_prior_scalar_from_cache(ars_unc)
    inc_unc = _duration_prior_scalar_from_cache(inc_unc)

    fallback_sigma_log = float(np.log1p(DURATION_PRIOR_FALLBACK_RELATIVE_SIGMA))
    min_sigma_log = float(np.log1p(DURATION_PRIOR_MIN_RELATIVE_SIGMA))
    sigma_log_duration = fallback_sigma_log
    sample_count = 0
    source = "fallback relative width"

    if any(
        value is not None and value > 0
        for value in (period_unc, rprs_unc, ars_unc, inc_unc)
    ):
        rng = np.random.default_rng(0)
        sample_draw_count = int(max(DURATION_PRIOR_MONTE_CARLO_SAMPLES, 1))
        period_draws = np.full(sample_draw_count, prior['per'], dtype=float)
        rprs_draws = np.full(sample_draw_count, prior['rprs'], dtype=float)
        ars_draws = np.full(sample_draw_count, prior['ars'], dtype=float)
        inc_draws = np.full(sample_draw_count, prior['inc'], dtype=float)

        if period_unc is not None and period_unc > 0:
            period_draws = rng.normal(prior['per'], period_unc, sample_draw_count)
        if rprs_unc is not None and rprs_unc > 0:
            rprs_draws = rng.normal(prior['rprs'], rprs_unc, sample_draw_count)
        if ars_unc is not None and ars_unc > 0:
            ars_draws = rng.normal(prior['ars'], ars_unc, sample_draw_count)
        if inc_unc is not None and inc_unc > 0:
            inc_draws = rng.normal(prior['inc'], inc_unc, sample_draw_count)

        inc_draws = np.clip(inc_draws, 1e-6, 89.999999)
        durations = np.full(sample_draw_count, np.nan, dtype=float)
        for index in range(sample_draw_count):
            durations[index] = estimate_transit_duration_from_prior_geometry({
                'per': period_draws[index],
                'rprs': rprs_draws[index],
                'ars': ars_draws[index],
                'inc': inc_draws[index],
                'ecc': prior['ecc'],
                'omega': prior['omega'],
            })

        valid_durations = durations[np.isfinite(durations) & (durations > 0)]
        sample_count = int(valid_durations.size)
        if valid_durations.size >= DURATION_PRIOR_MIN_VALID_MONTE_CARLO_SAMPLES:
            log_offsets = np.log(valid_durations / expected_duration)
            lower_offset, upper_offset = np.nanpercentile(log_offsets, [16, 84])
            sigma_log_duration = float(max(0.5 * (upper_offset - lower_offset), min_sigma_log))
            source = "published geometry uncertainties"
        else:
            source = "fallback relative width (insufficient valid uncertainty samples)"
    else:
        source = "fallback relative width (missing published geometry uncertainties)"

    sigma_log_duration = float(max(sigma_log_duration, min_sigma_log))
    relative_sigma = float(np.expm1(sigma_log_duration))
    note = (
        f"Applied; expected duration={expected_duration:.6f} day(s) with an approximate 1-sigma width of "
        f"{relative_sigma * 100.0:.1f}% from {source}"
    )
    if sample_count > 0:
        note += f" ({sample_count} propagated sample(s))."
    else:
        note += "."

    return {
        'applied': True,
        'expected_duration': float(expected_duration),
        'sigma_log_duration': float(sigma_log_duration),
        'relative_sigma': relative_sigma,
        'source': source,
        'sample_count': sample_count,
        'note': note,
    }


def build_single_transit_duration_prior(planet_dict):
    if not isinstance(planet_dict, dict):
        return {
            'applied': False,
            'expected_duration': np.nan,
            'sigma_log_duration': np.nan,
            'relative_sigma': np.nan,
            'source': 'unavailable',
            'sample_count': 0,
            'note': "Not applied; missing published single-transit planet metadata.",
        }

    return dict(
        _cached_single_transit_duration_prior(
            _cacheable_duration_prior_scalar(planet_dict.get('pPer', np.nan)),
            _cacheable_duration_prior_scalar(planet_dict.get('rprs', np.nan)),
            _cacheable_duration_prior_scalar(planet_dict.get('aRs', np.nan)),
            _cacheable_duration_prior_scalar(planet_dict.get('inc', np.nan)),
            _cacheable_duration_prior_scalar(planet_dict.get('ecc', 0.0)),
            _cacheable_duration_prior_scalar(planet_dict.get('omega', 0.0)),
            _cacheable_duration_prior_scalar(planet_dict.get('pPerUnc', np.nan)),
            _cacheable_duration_prior_scalar(planet_dict.get('rprsUnc', np.nan)),
            _cacheable_duration_prior_scalar(planet_dict.get('aRsUnc', np.nan)),
            _cacheable_duration_prior_scalar(planet_dict.get('incUnc', np.nan)),
        )
    )


def estimate_ephemeris_tmid_and_bounds(
    times,
    prior_tmid,
    period,
    midt_unc,
    per_unc,
    expected_duration=np.nan,
    sigma_multiplier=25.0,
):
    summary = {
        'method': 'ephemeris',
        'applied': False,
        'tmid': float(prior_tmid) if np.isfinite(prior_tmid) else np.nan,
        'uncertainty': np.nan,
        'bounds': [np.nan, np.nan],
        'cycle_index': np.nan,
        'propagated_half_width': np.nan,
        'half_width': np.nan,
        'observations_bracket_expected_transit': False,
        'duration_capped': False,
        'observed_window_capped': False,
        'note': 'Using ephemeris-derived Tmid bounds.',
    }

    try:
        prior_tmid = float(prior_tmid)
        period = float(period)
        midt_unc = float(midt_unc)
        per_unc = float(per_unc)
        sigma_multiplier = float(sigma_multiplier)
    except (TypeError, ValueError):
        summary['note'] = 'Using ephemeris-derived Tmid bounds with invalid prior metadata.'
        return summary

    if not np.isfinite(prior_tmid) or not np.isfinite(period) or period <= 0:
        summary['note'] = 'Using ephemeris-derived Tmid bounds with invalid Tmid/period metadata.'
        return summary

    valid_times = np.asarray(times, dtype=float)
    valid_times = valid_times[np.isfinite(valid_times)]
    if valid_times.size == 0:
        summary['bounds'] = [prior_tmid, prior_tmid]
        summary['note'] = 'Using ephemeris-derived Tmid bounds with no finite observation times.'
        return summary

    phases = (valid_times - prior_tmid) / period
    cycle_index = float(np.floor(phases).max())
    tmid = float(prior_tmid + cycle_index * period)
    propagated_uncertainty = np.sqrt(midt_unc ** 2 + (cycle_index * per_unc) ** 2)

    propagated_half_width = np.abs(sigma_multiplier * midt_unc + cycle_index * sigma_multiplier * per_unc)
    max_half_width = 0.25 * period
    if not np.isfinite(propagated_half_width) or propagated_half_width <= 0:
        half_width = max_half_width
        propagated_half_width = np.nan
    else:
        half_width = min(float(propagated_half_width), max_half_width)

    cadence = np.nan
    if valid_times.size > 1:
        cadence = np.nanmedian(np.diff(np.sort(valid_times)))

    lower = float(tmid - half_width)
    upper = float(tmid + half_width)
    if np.isfinite(expected_duration) and expected_duration > 0:
        coverage_margin = 0.5 * float(expected_duration)
        if np.isfinite(cadence) and cadence > 0:
            coverage_margin = max(coverage_margin, 3.0 * cadence)

        pre_points = int(np.count_nonzero(valid_times < tmid - coverage_margin))
        post_points = int(np.count_nonzero(valid_times > tmid + coverage_margin))
        bracketed = pre_points > 0 and post_points > 0
        summary['observations_bracket_expected_transit'] = bracketed

        cadence_floor = 0.0
        if np.isfinite(cadence) and cadence > 0:
            cadence_floor = 5.0 * cadence
        duration_cap = max(
            EPHEMERIS_BRACKETED_TMID_HALF_WIDTH_DURATION_MULTIPLIER * float(expected_duration),
            cadence_floor,
        )
        if bracketed and np.isfinite(duration_cap) and duration_cap > 0 and duration_cap < half_width:
            half_width = float(duration_cap)
            summary['duration_capped'] = True
            lower = float(tmid - half_width)
            upper = float(tmid + half_width)

        if bracketed:
            observed_lower = float(np.nanmin(valid_times) + 0.5 * float(expected_duration))
            observed_upper = float(np.nanmax(valid_times) - 0.5 * float(expected_duration))
            if (
                np.isfinite(observed_lower)
                and np.isfinite(observed_upper)
                and observed_upper > observed_lower
            ):
                tightened_lower = max(lower, observed_lower)
                tightened_upper = min(upper, observed_upper)
                if tightened_upper > tightened_lower and (
                    tightened_lower > lower + 1e-12 or tightened_upper < upper - 1e-12
                ):
                    lower = float(tightened_lower)
                    upper = float(tightened_upper)
                    summary['observed_window_capped'] = True

    half_width = max(float(tmid - lower), float(upper - tmid))
    summary.update({
        'tmid': tmid,
        'uncertainty': propagated_uncertainty,
        'bounds': [lower, upper],
        'cycle_index': cycle_index,
        'propagated_half_width': propagated_half_width,
        'half_width': half_width,
        'applied': summary['duration_capped'] or summary['observed_window_capped'],
    })
    if summary['observed_window_capped']:
        summary['note'] = (
            "Ephemeris-derived Tmid bounds were intersected with the observed time span needed to contain the "
            f"full expected transit; using bounds=[{lower:.6f}, {upper:.6f}] instead of the wider propagated "
            f"half-width {float(propagated_half_width):.6f} day(s)."
        )
    elif summary['duration_capped']:
        summary['note'] = (
            "Ephemeris-derived Tmid bounds were narrowed to the expected-transit timescale because the "
            f"observations bracket the expected transit; using bounds=[{lower:.6f}, {upper:.6f}] "
            f"instead of the wider propagated half-width {float(propagated_half_width):.6f} day(s)."
        )
    else:
        summary['note'] = f"Using ephemeris-derived Tmid bounds [{lower:.6f}, {upper:.6f}]."

    return summary


def estimate_tmid_and_bounds_with_eebls(times, flux_values, flux_errors, prior, fallback_bounds):
    summary = {
        'method': 'ephemeris',
        'applied': False,
        'tmid': float(prior.get('tmid', np.nan)),
        'bounds': [float(fallback_bounds[0]), float(fallback_bounds[1])],
        'duration': np.nan,
        'depth': np.nan,
        'depth_snr': np.nan,
        'note': 'EEBLS transit initializer did not run.',
    }

    times = np.asarray(times, dtype=float)
    flux_values = np.asarray(flux_values, dtype=float)
    flux_errors = np.asarray(flux_errors, dtype=float)
    period = float(prior.get('per', np.nan))
    if not np.isfinite(period) or period <= 0:
        summary['note'] = 'EEBLS transit initializer skipped: invalid orbital period.'
        return summary

    valid = np.isfinite(times) & np.isfinite(flux_values) & (flux_values > 0)
    if flux_errors.shape == flux_values.shape:
        valid &= np.isfinite(flux_errors) & (flux_errors > 0)
    else:
        flux_errors = np.full_like(flux_values, np.nan, dtype=float)

    if np.count_nonzero(valid) < max(LIGHTCURVE_MIN_VALID_POINTS, EEBLS_MIN_VALID_POINTS):
        summary['note'] = 'EEBLS transit initializer skipped: not enough valid points.'
        return summary

    valid_times = np.asarray(times[valid], dtype=float)
    valid_flux = np.asarray(flux_values[valid], dtype=float)
    valid_errors = np.asarray(flux_errors[valid], dtype=float)
    sort_index = np.argsort(valid_times)
    fit_times = valid_times[sort_index]
    fit_flux = valid_flux[sort_index]
    fit_errors = valid_errors[sort_index]
    cadence = np.nanmedian(np.diff(fit_times))
    if not np.isfinite(cadence) or cadence <= 0:
        cadence = max(np.finfo(float).eps, 0.005 * period)

    x = fit_times - np.nanmedian(fit_times)
    baseline = np.ones_like(fit_flux, dtype=float)
    design = np.column_stack((np.ones_like(x), x))
    if np.count_nonzero(np.isfinite(x)) >= 2:
        weights = np.ones_like(fit_flux, dtype=float)
        finite_error_mask = np.isfinite(fit_errors) & (fit_errors > 0)
        if np.any(finite_error_mask):
            weights[finite_error_mask] = 1.0 / (fit_errors[finite_error_mask] ** 2)
            weights[~finite_error_mask] = 0.0
            if not np.any(weights > 0):
                weights = np.ones_like(fit_flux, dtype=float)
        sqrt_weights = np.sqrt(weights)
        try:
            coeffs, _, _, _ = np.linalg.lstsq(design * sqrt_weights[:, None], fit_flux * sqrt_weights, rcond=None)
            baseline = coeffs[0] + coeffs[1] * x
            if not np.all(np.isfinite(baseline)) or np.any(baseline <= 0):
                baseline = np.ones_like(fit_flux, dtype=float)
        except np.linalg.LinAlgError:
            baseline = np.ones_like(fit_flux, dtype=float)

    detrended_flux = fit_flux / baseline
    detrended_flux /= np.nanmedian(detrended_flux)
    detrended_errors = fit_errors / baseline
    if not np.all(np.isfinite(detrended_errors)) or np.any(detrended_errors <= 0):
        detrended_errors = None

    expected_duration = estimate_transit_duration_from_prior_geometry(prior)
    if not np.isfinite(expected_duration) or expected_duration <= 0:
        expected_duration = 0.05 * period

    min_duration = max(3.0 * cadence, EEBLS_DURATION_MIN_FRACTION * expected_duration)
    max_duration = min(0.25 * period, max(min_duration * 1.5, EEBLS_DURATION_MAX_FRACTION * expected_duration))
    if not np.isfinite(min_duration) or not np.isfinite(max_duration) or max_duration <= 0 or min_duration > max_duration:
        summary['note'] = 'EEBLS transit initializer skipped: invalid duration search grid.'
        return summary

    durations = np.linspace(min_duration, max_duration, EEBLS_DURATION_GRID_SIZE)
    durations = np.unique(durations[np.isfinite(durations) & (durations > 0)])
    if durations.size == 0:
        summary['note'] = 'EEBLS transit initializer skipped: empty duration search grid.'
        return summary

    try:
        bls = BoxLeastSquares(fit_times, detrended_flux, dy=detrended_errors)
        results = bls.power(period, durations, objective='snr')
    except Exception as exc:
        summary['note'] = f'EEBLS transit initializer failed: {type(exc).__name__}: {exc}'
        return summary

    power = np.asarray(results.power, dtype=float)
    if power.size == 0 or not np.any(np.isfinite(power)):
        summary['note'] = 'EEBLS transit initializer skipped: no finite search power values were returned.'
        return summary

    best_index = int(np.nanargmax(power))
    tmid = float(np.asarray(results.transit_time, dtype=float)[best_index])
    duration = float(np.asarray(results.duration, dtype=float)[best_index])
    depth = float(np.asarray(results.depth, dtype=float)[best_index])
    depth_snr = float(np.asarray(results.depth_snr, dtype=float)[best_index])
    if (
        not np.isfinite(tmid)
        or not np.isfinite(duration) or duration <= 0
        or not np.isfinite(depth) or depth <= 0
        or not np.isfinite(depth_snr) or depth_snr <= 0
    ):
        summary['note'] = 'EEBLS transit initializer skipped: the best-fitting transit candidate was not physical.'
        return summary

    coverage_margin = max(0.5 * duration, 3.0 * cadence)
    pre_points = int(np.count_nonzero(np.isfinite(fit_times) & (fit_times < tmid - coverage_margin)))
    post_points = int(np.count_nonzero(np.isfinite(fit_times) & (fit_times > tmid + coverage_margin)))
    if pre_points == 0 or post_points == 0:
        summary.update({
            'method': 'eebls',
            'applied': False,
            'tmid': tmid,
            'duration': duration,
            'depth': depth,
            'depth_snr': depth_snr,
        })
        summary['note'] = (
            "EEBLS transit initializer found a box-like event, but it is not bracketed by data on both sides "
            f"({pre_points} pre-point(s), {post_points} post-point(s)); keeping the EEBLS depth SNR only and "
            "falling back to the non-EEBLS Tmid bounds."
        )
        return summary

    duration_for_bounds = duration
    if np.isfinite(expected_duration) and expected_duration > 0:
        duration_for_bounds = max(duration_for_bounds, 0.75 * expected_duration)
    half_width = min(
        0.25 * period,
        max(EEBLS_TMID_HALF_WIDTH_DURATION_MULTIPLIER * duration_for_bounds, 5.0 * cadence),
    )
    if not np.isfinite(half_width) or half_width <= 0:
        summary['note'] = 'EEBLS transit initializer skipped: invalid Tmid search half-width.'
        return summary

    summary.update({
        'method': 'eebls',
        'applied': True,
        'tmid': tmid,
        'bounds': [float(tmid - half_width), float(tmid + half_width)],
        'duration': duration,
        'depth': depth,
        'depth_snr': depth_snr,
        'note': (
            "EEBLS transit initializer found a box-like transit candidate at "
            f"Tmid={tmid:.6f} day(s) with duration={duration:.6f} day(s), depth={depth:.5f}, "
            f"depth_snr={depth_snr:.2f}, and bounds=[{tmid - half_width:.6f}, {tmid + half_width:.6f}]."
        ),
    })
    return summary


def annotate_lightcurve_tmid_search(fit, summary):
    if fit is None:
        return

    fit.initial_tmid_search_method = summary.get('method')
    fit.initial_tmid_search_applied = bool(summary.get('applied'))
    fit.initial_tmid_search_tmid = summary.get('tmid')
    fit.initial_tmid_search_uncertainty = summary.get('uncertainty')
    fit.initial_tmid_search_bounds = summary.get('bounds')
    fit.initial_tmid_search_duration = summary.get('duration')
    fit.initial_tmid_search_depth = summary.get('depth')
    fit.initial_tmid_search_depth_snr = summary.get('depth_snr')
    fit.initial_tmid_search_note = summary.get('note')


def annotate_lightcurve_eebls_diagnostic(fit, summary):
    if fit is None:
        return

    summary = {} if summary is None else dict(summary)
    fit.eebls_diagnostic_computed = bool(summary)
    fit.eebls_diagnostic_method = summary.get('method')
    fit.eebls_diagnostic_applied = bool(summary.get('applied'))
    fit.eebls_diagnostic_tmid = summary.get('tmid')
    fit.eebls_diagnostic_bounds = summary.get('bounds')
    fit.eebls_diagnostic_duration = summary.get('duration')
    fit.eebls_diagnostic_depth = summary.get('depth')
    fit.eebls_diagnostic_depth_snr = summary.get('depth_snr')
    fit.eebls_diagnostic_note = summary.get('note')


def extract_lightcurve_fit_eebls_snr(fit):
    if fit is None:
        return np.nan

    for attr_name in ('eebls_diagnostic_depth_snr', 'initial_tmid_search_depth_snr'):
        value = getattr(fit, attr_name, np.nan)
        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(numeric_value):
            return numeric_value

    return np.nan


def ensure_lightcurve_fit_eebls_diagnostic(fit):
    if fit is None:
        return None

    existing_snr = extract_lightcurve_fit_eebls_snr(fit)
    if np.isfinite(existing_snr):
        return {
            'applied': True,
            'depth_snr': existing_snr,
            'note': 'Existing EEBLS diagnostic reused.',
        }

    times = np.asarray(getattr(fit, 'time', np.array([])), dtype=float)
    flux_values = np.asarray(getattr(fit, 'data', np.array([])), dtype=float)
    if times.ndim != 1 or times.size < LIGHTCURVE_MIN_VALID_POINTS or flux_values.shape != times.shape:
        return None

    dataerr_obj = getattr(fit, 'dataerr', None)
    flux_errors = None if dataerr_obj is None else np.asarray(dataerr_obj, dtype=float)
    if flux_errors is None or flux_errors.shape != times.shape:
        flux_errors = np.full(times.shape, 1.0, dtype=float)
    else:
        finite_positive = np.isfinite(flux_errors) & (flux_errors > 0)
        if not np.any(finite_positive):
            flux_errors = np.full(times.shape, 1.0, dtype=float)
        elif not np.all(finite_positive):
            replacement = float(np.nanmedian(flux_errors[finite_positive]))
            flux_errors = np.where(finite_positive, flux_errors, replacement)

    fit_prior = getattr(fit, 'prior', {}) or {}
    parameters = getattr(fit, 'parameters', {}) or {}
    period = fit_prior.get('per', parameters.get('per', np.nan))
    tmid = fit_prior.get('tmid', parameters.get('tmid', np.nan))
    try:
        period = float(period)
        tmid = float(tmid)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(period) or period <= 0 or not np.isfinite(tmid):
        return None

    fallback_bounds = getattr(fit, 'initial_tmid_search_bounds', None)
    if fallback_bounds is None:
        fallback_bounds = getattr(fit, 'eebls_diagnostic_bounds', None)
    if fallback_bounds is None:
        bounds = getattr(fit, 'bounds', None)
        fallback_bounds = bounds.get('tmid') if isinstance(bounds, dict) else None

    try:
        lower, upper = [float(value) for value in np.asarray(fallback_bounds, dtype=float).reshape(-1)[:2]]
    except (TypeError, ValueError, IndexError):
        lower = float(tmid - 0.25 * period)
        upper = float(tmid + 0.25 * period)
    if not np.isfinite(lower) or not np.isfinite(upper) or upper <= lower:
        lower = float(tmid - 0.25 * period)
        upper = float(tmid + 0.25 * period)

    eebls_summary = estimate_tmid_and_bounds_with_eebls(
        times,
        flux_values,
        flux_errors,
        {'per': period, 'tmid': tmid},
        [lower, upper],
    )
    annotate_lightcurve_eebls_diagnostic(fit, eebls_summary)
    return eebls_summary


def should_use_impactparameter_rather_than_inclination_to_fit(config_value):
    if config_value is None:
        return True
    if isinstance(config_value, bool):
        return config_value
    if isinstance(config_value, (int, float)):
        return bool(config_value)
    if isinstance(config_value, str):
        normalized = config_value.strip().lower()
        if normalized in ('y', 'yes', 'true', '1', 'on'):
            return True
        if normalized in ('n', 'no', 'false', '0', 'off', ''):
            return False

    log_info(
        "Warning: Invalid 'use_impactparameter_rather_than_inclination_to_fit' value; "
        "using impact parameter for nested fitting.",
        warn=True,
    )
    return True


def apply_vertical_flux_normalization_bound(prior, bounds, flux_values, disabled):
    finite_flux = np.asarray(flux_values, dtype=float)
    finite_flux = finite_flux[np.isfinite(finite_flux) & (finite_flux > 0)]
    baseline_guess = 1.0 if finite_flux.size == 0 else float(np.nanmedian(finite_flux))
    if not np.isfinite(baseline_guess) or baseline_guess <= 0:
        baseline_guess = 1.0

    prior['a0'] = baseline_guess
    prior['a1'] = baseline_guess

    if not disabled:
        # Some paths deliver an approximately unity-normalized light curve, while
        # others still carry an arbitrary positive baseline. Keep the legacy
        # near-unity bound only when the working series is already close to 1.
        if 0.95 <= baseline_guess <= 1.05:
            bounds['a0'] = [0.95, 1.05]
        else:
            lower = max(np.finfo(float).eps, baseline_guess * 0.75)
            upper = baseline_guess * 1.25
            bounds['a0'] = [lower, upper]


def ensure_pre_final_ultranest_baseline_bounds(prior, bounds, flux_values, fit_a2=True):
    finite_flux = np.asarray(flux_values, dtype=float)
    finite_flux = finite_flux[np.isfinite(finite_flux) & (finite_flux > 0)]
    baseline_guess = prior.get('a0', prior.get('a1', np.nan))
    try:
        baseline_guess = float(baseline_guess)
    except (TypeError, ValueError):
        baseline_guess = np.nan
    if not np.isfinite(baseline_guess) or baseline_guess <= 0:
        baseline_guess = 1.0 if finite_flux.size == 0 else float(np.nanmedian(finite_flux))
    if not np.isfinite(baseline_guess) or baseline_guess <= 0:
        baseline_guess = 1.0

    prior['a0'] = baseline_guess
    prior['a1'] = baseline_guess
    if 'a0' not in bounds and 'a1' not in bounds:
        if 0.95 <= baseline_guess <= 1.05:
            bounds['a0'] = [0.95, 1.05]
        else:
            lower = max(np.finfo(float).eps, baseline_guess * 0.75)
            upper = baseline_guess * 1.25
            bounds['a0'] = [lower, upper]

    if fit_a2:
        prior['a2'] = prior.get('a2', 0.0)
        if 'a2' not in bounds:
            bounds['a2'] = list(TRANSIT_QC_DEFAULT_A2_BOUNDS)


def _weighted_mean_with_fallback(values, weights=None):
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if not np.any(finite):
        return np.nan

    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        valid_weights = finite & np.isfinite(weights) & (weights > 0)
        if np.any(valid_weights):
            return float(np.sum(values[valid_weights] * weights[valid_weights]) / np.sum(weights[valid_weights]))

    return float(np.nanmean(values[finite]))


def build_fast_ultranest_lightcurve_series(
    times,
    flux_values,
    flux_errors,
    airmass,
    jd_times=None,
    max_points=FAST_ULTRANEST_MAX_BINNED_POINTS,
    min_points_to_bin=FAST_ULTRANEST_MIN_POINTS_TO_BIN,
):
    times = np.asarray(times, dtype=float)
    flux_values = np.asarray(flux_values, dtype=float)
    flux_errors = np.asarray(flux_errors, dtype=float)
    airmass = np.asarray(airmass, dtype=float)
    jd_array = None if jd_times is None else np.asarray(jd_times, dtype=float)

    base_result = {
        'applied': False,
        'note': None,
        'time': times,
        'flux': flux_values,
        'unc': flux_errors,
        'airmass': airmass,
        'jd_times': jd_array,
        'original_point_count': int(times.shape[0]),
        'binned_point_count': int(times.shape[0]),
        'bin_indices': None,
    }

    if not (times.shape == flux_values.shape == flux_errors.shape == airmass.shape):
        base_result['note'] = 'Skipped; light-curve arrays were not aligned for fast UltraNest binning.'
        return base_result
    if jd_array is not None and jd_array.shape != times.shape:
        base_result['note'] = 'Skipped; JD timestamps were not aligned for fast UltraNest binning.'
        return base_result

    point_count = int(times.shape[0])
    if point_count <= int(min_points_to_bin):
        base_result['note'] = (
            f"Skipped; {point_count} point(s) did not exceed the fast UltraNest "
            f"binning threshold of {int(min_points_to_bin)}."
        )
        return base_result

    max_points = int(max(1, max_points))
    target_points = min(max_points, point_count)
    valid = (
        np.isfinite(times)
        & np.isfinite(flux_values)
        & np.isfinite(flux_errors)
        & (flux_errors > 0)
        & np.isfinite(airmass)
    )
    if jd_array is not None:
        valid &= np.isfinite(jd_array)
    if np.count_nonzero(valid) <= target_points:
        base_result['note'] = 'Skipped; too few finite points remained for fast UltraNest binning.'
        return base_result

    ordered_indices = np.flatnonzero(valid)[np.argsort(times[valid])]
    chunks = [chunk for chunk in np.array_split(ordered_indices, target_points) if chunk.size > 0]
    if len(chunks) >= point_count or not chunks:
        base_result['note'] = 'Skipped; fast UltraNest binning would not reduce the light curve.'
        return base_result

    binned_time = []
    binned_flux = []
    binned_unc = []
    binned_airmass = []
    binned_jd = [] if jd_array is not None else None
    for chunk in chunks:
        chunk_unc = flux_errors[chunk]
        weights = np.zeros(chunk_unc.shape, dtype=float)
        valid_unc = np.isfinite(chunk_unc) & (chunk_unc > 0)
        weights[valid_unc] = 1.0 / (chunk_unc[valid_unc] ** 2)
        binned_time.append(_weighted_mean_with_fallback(times[chunk], weights))
        binned_flux.append(_weighted_mean_with_fallback(flux_values[chunk], weights))
        if np.any(weights > 0):
            binned_unc.append(float(np.sqrt(1.0 / np.sum(weights[weights > 0]))))
        else:
            scatter = float(np.nanstd(flux_values[chunk]))
            binned_unc.append(scatter / np.sqrt(max(chunk.size, 1)) if np.isfinite(scatter) else np.nan)
        binned_airmass.append(_weighted_mean_with_fallback(airmass[chunk], weights))
        if jd_array is not None:
            binned_jd.append(_weighted_mean_with_fallback(jd_array[chunk], weights))

    binned_time = np.asarray(binned_time, dtype=float)
    binned_flux = np.asarray(binned_flux, dtype=float)
    binned_unc = np.asarray(binned_unc, dtype=float)
    binned_airmass = np.asarray(binned_airmass, dtype=float)
    finite_binned = (
        np.isfinite(binned_time)
        & np.isfinite(binned_flux)
        & np.isfinite(binned_unc)
        & (binned_unc > 0)
        & np.isfinite(binned_airmass)
    )
    if binned_jd is not None:
        binned_jd = np.asarray(binned_jd, dtype=float)
        finite_binned &= np.isfinite(binned_jd)

    if np.count_nonzero(finite_binned) < LIGHTCURVE_MIN_VALID_POINTS:
        base_result['note'] = 'Skipped; fast UltraNest binning produced too few finite bins.'
        return base_result

    result = dict(base_result)
    result.update({
        'applied': True,
        'time': binned_time[finite_binned],
        'flux': binned_flux[finite_binned],
        'unc': binned_unc[finite_binned],
        'airmass': binned_airmass[finite_binned],
        'jd_times': None if binned_jd is None else binned_jd[finite_binned],
        'binned_point_count': int(np.count_nonzero(finite_binned)),
        'bin_indices': [chunk.tolist() for i, chunk in enumerate(chunks) if finite_binned[i]],
        'note': (
            f"Using fast UltraNest binning for pre-final runs: "
            f"{point_count} point(s) -> {int(np.count_nonzero(finite_binned))} binned point(s)."
        ),
    })
    return result


def annotate_fast_ultranest_binning(fit, binning_result):
    if fit is None or not isinstance(binning_result, dict):
        return
    fit.fast_ultranest_binning_applied = bool(binning_result.get('applied', False))
    fit.fast_ultranest_original_point_count = int(binning_result.get('original_point_count', 0))
    fit.fast_ultranest_binned_point_count = int(binning_result.get('binned_point_count', 0))
    fit.fast_ultranest_binning_note = binning_result.get('note')


def summarize_initial_fit_transit_coverage(
    times,
    fit,
    flux_values=None,
    flux_errors=None,
    depth_fraction=OUT_OF_TRANSIT_BASELINE_DEPTH_FRACTION,
):
    times = np.asarray(times, dtype=float)
    transit_model = np.asarray(getattr(fit, 'transit', []), dtype=float)
    if flux_values is None:
        flux_values = np.ones_like(times, dtype=float)
    else:
        flux_values = np.asarray(flux_values, dtype=float)

    if transit_model.shape != times.shape or flux_values.shape != times.shape:
        return {
            'valid': False,
            'note': 'initial fit did not provide a transit model aligned with the light curve.',
        }

    valid = (
        np.isfinite(times)
        & np.isfinite(flux_values)
        & (flux_values > 0)
        & np.isfinite(transit_model)
    )
    if flux_errors is not None:
        flux_errors = np.asarray(flux_errors, dtype=float)
        if flux_errors.shape == flux_values.shape:
            valid &= np.isfinite(flux_errors) & (flux_errors > 0)

    if np.count_nonzero(valid) < 3:
        return {
            'valid': False,
            'note': 'not enough finite flux points remain to isolate the modeled transit window.',
        }

    depth = np.clip(1.0 - transit_model, 0.0, None)
    max_depth = np.nanmax(depth[valid])
    if not np.isfinite(max_depth) or max_depth <= 0:
        return {
            'valid': False,
            'note': 'initial fit did not produce a measurable transit depth for baseline isolation.',
        }

    threshold = max(1e-6, depth_fraction * max_depth)
    in_transit = valid & (depth > threshold)
    if not np.any(in_transit):
        return {
            'valid': False,
            'note': 'could not isolate ingress and egress from the initial fit.',
        }

    ingress_time = float(np.nanmin(times[in_transit]))
    egress_time = float(np.nanmax(times[in_transit]))
    oot_mask = valid & ((times < ingress_time) | (times > egress_time))

    mid_transit = float(getattr(fit, 'parameters', {}).get('tmid', np.nanmedian(times[valid])))
    pre_mask = oot_mask & (times < mid_transit)
    post_mask = oot_mask & (times > mid_transit)
    pre_points = int(np.count_nonzero(pre_mask))
    post_points = int(np.count_nonzero(post_mask))
    has_two_sided_oot = pre_points > 0 and post_points > 0

    summary = {
        'valid': True,
        'note': None,
        'mid_transit': mid_transit,
        'ingress_time': ingress_time,
        'egress_time': egress_time,
        'in_transit_mask': in_transit,
        'oot_mask': oot_mask,
        'pre_mask': pre_mask,
        'post_mask': post_mask,
        'pre_points': pre_points,
        'post_points': post_points,
        'has_two_sided_oot': has_two_sided_oot,
    }
    if not has_two_sided_oot:
        summary['note'] = 'need out-of-transit coverage on both sides of transit to fit a linear baseline.'
    return summary


def summarize_prior_transit_coverage(
    times,
    prior,
    flux_values=None,
    flux_errors=None,
):
    times = np.asarray(times, dtype=float)
    if flux_values is None:
        flux_values = np.ones_like(times, dtype=float)
    else:
        flux_values = np.asarray(flux_values, dtype=float)

    if times.shape != flux_values.shape:
        return {
            'valid': False,
            'note': 'prior-based transit coverage could not be aligned with the light curve.',
        }

    valid = np.isfinite(times) & np.isfinite(flux_values) & (flux_values > 0)
    if flux_errors is not None:
        flux_errors = np.asarray(flux_errors, dtype=float)
        if flux_errors.shape == flux_values.shape:
            valid &= np.isfinite(flux_errors) & (flux_errors > 0)

    if np.count_nonzero(valid) < 3:
        return {
            'valid': False,
            'note': 'not enough finite flux points remain to evaluate prior-based transit coverage.',
        }

    try:
        mid_transit = float(prior.get('tmid', np.nan))
    except (AttributeError, TypeError, ValueError):
        mid_transit = np.nan
    if not np.isfinite(mid_transit):
        return {
            'valid': False,
            'note': 'prior-based transit coverage skipped: no finite ephemeris-centered Tmid was available.',
        }

    duration = estimate_transit_duration_from_prior_geometry(prior)
    if not np.isfinite(duration) or duration <= 0:
        return {
            'valid': False,
            'note': 'prior-based transit coverage skipped: could not estimate a physical transit duration from the priors.',
        }

    ingress_time = float(mid_transit - 0.5 * duration)
    egress_time = float(mid_transit + 0.5 * duration)
    in_transit = valid & (times >= ingress_time) & (times <= egress_time)
    if not np.any(in_transit):
        return {
            'valid': False,
            'note': 'prior-based transit coverage skipped: the ephemeris-centered transit window does not overlap the observations.',
        }

    oot_mask = valid & ((times < ingress_time) | (times > egress_time))
    pre_mask = oot_mask & (times < mid_transit)
    post_mask = oot_mask & (times > mid_transit)
    pre_points = int(np.count_nonzero(pre_mask))
    post_points = int(np.count_nonzero(post_mask))
    has_two_sided_oot = pre_points > 0 and post_points > 0

    summary = {
        'valid': True,
        'note': None,
        'mid_transit': mid_transit,
        'ingress_time': ingress_time,
        'egress_time': egress_time,
        'in_transit_mask': in_transit,
        'oot_mask': oot_mask,
        'pre_mask': pre_mask,
        'post_mask': post_mask,
        'pre_points': pre_points,
        'post_points': post_points,
        'has_two_sided_oot': has_two_sided_oot,
        'used_prior_ephemeris': True,
        'duration': duration,
    }
    if not has_two_sided_oot:
        summary['note'] = (
            'need out-of-transit coverage on both sides of the ephemeris-centered transit window '
            'to fit a linear baseline.'
        )
    return summary


def _coverage_duration_from_context(prior, duration_prior=None):
    if isinstance(duration_prior, dict):
        duration = coerce_finite_transit_qc_scalar(duration_prior.get('expected_duration', np.nan))
        if np.isfinite(duration) and duration > 0:
            return float(duration)
    return estimate_transit_duration_from_prior_geometry(prior)


def _coverage_tmid_from_context(prior, tmid_search_summary=None):
    if isinstance(tmid_search_summary, dict):
        tmid = coerce_finite_transit_qc_scalar(tmid_search_summary.get('tmid', np.nan))
        if np.isfinite(tmid):
            return float(tmid)
    try:
        return float(prior.get('tmid', np.nan))
    except (AttributeError, TypeError, ValueError):
        return np.nan


def build_ephemeris_tmid_search_summary_for_coverage(
    times,
    planet_dict,
    prior=None,
    duration_prior=None,
    sigma_multiplier=35.0,
):
    if not isinstance(planet_dict, dict):
        return None
    prior = prior if isinstance(prior, dict) else {}

    prior_tmid = coerce_finite_transit_qc_scalar(
        planet_dict.get('midT', prior.get('tmid', np.nan))
    )
    period = coerce_finite_transit_qc_scalar(
        planet_dict.get('pPer', prior.get('per', np.nan))
    )
    midt_unc = coerce_finite_transit_qc_scalar(planet_dict.get('midTUnc', 0.0))
    per_unc = coerce_finite_transit_qc_scalar(planet_dict.get('pPerUnc', 0.0))
    if not np.isfinite(midt_unc):
        midt_unc = 0.0
    if not np.isfinite(per_unc):
        per_unc = 0.0
    if not np.isfinite(prior_tmid) or not np.isfinite(period) or period <= 0:
        return None

    coverage_prior = dict(prior)
    coverage_prior.setdefault('tmid', prior_tmid)
    coverage_prior.setdefault('per', period)
    coverage_prior.setdefault('rprs', planet_dict.get('rprs', np.nan))
    coverage_prior.setdefault('ars', planet_dict.get('aRs', np.nan))
    coverage_prior.setdefault('inc', planet_dict.get('inc', np.nan))
    coverage_prior.setdefault('ecc', planet_dict.get('ecc', 0.0))
    coverage_prior.setdefault('omega', planet_dict.get('omega', 0.0))
    expected_duration = _coverage_duration_from_context(
        coverage_prior,
        duration_prior=duration_prior,
    )

    return estimate_ephemeris_tmid_and_bounds(
        times,
        prior_tmid,
        period,
        midt_unc,
        per_unc,
        expected_duration=expected_duration,
        sigma_multiplier=sigma_multiplier,
    )


def expected_transit_observed_segment(
    observed_start,
    observed_end,
    ingress_time,
    mid_transit,
    egress_time,
):
    if observed_end < ingress_time:
        return "pre-transit baseline only"
    if observed_start > egress_time:
        return "post-transit baseline only"

    pieces = []
    if observed_start < ingress_time:
        pieces.append("pre-ingress baseline")
    if observed_start <= ingress_time <= observed_end:
        pieces.append("ingress")
    if observed_start <= mid_transit <= observed_end:
        pieces.append("mid-transit")
    if observed_start <= egress_time <= observed_end:
        pieces.append("egress")
    if observed_end > egress_time:
        pieces.append("post-egress baseline")
    if not pieces:
        if observed_end < mid_transit:
            return "inside the first half of transit"
        if observed_start > mid_transit:
            return "inside the second half of transit"
        return "inside the expected transit"
    return " plus ".join(pieces)


def score_expected_transit_model_success(
    transit_fraction_observed,
    covers_ingress,
    covers_mid_transit,
    covers_egress,
    pre_points,
    post_points,
):
    has_two_sided_baseline = pre_points > 0 and post_points > 0
    if transit_fraction_observed <= 0:
        return "very low", 0.05
    if transit_fraction_observed < 0.25:
        return "very low", 0.15
    if not has_two_sided_baseline:
        if transit_fraction_observed >= 0.9 and covers_ingress and covers_egress:
            return "moderate", 0.50
        if transit_fraction_observed >= 0.5 and covers_mid_transit:
            return "low", 0.35
        return "low", 0.25
    if transit_fraction_observed >= 0.9 and covers_ingress and covers_egress:
        return "high", 0.85
    if transit_fraction_observed >= 0.65 and covers_mid_transit and (covers_ingress or covers_egress):
        return "moderate", 0.65
    if transit_fraction_observed >= 0.4:
        return "low", 0.40
    return "low", 0.25


def build_expected_transit_coverage_assessment(
    times,
    prior,
    flux_values=None,
    flux_errors=None,
    tmid_search_summary=None,
    duration_prior=None,
):
    times = np.asarray(times, dtype=float)
    if flux_values is None:
        flux_values = np.ones_like(times, dtype=float)
    else:
        flux_values = np.asarray(flux_values, dtype=float)

    base = {
        'valid': False,
        'point_count': 0,
        'observed_start': np.nan,
        'observed_end': np.nan,
        'observed_span': np.nan,
        'expected_tmid': np.nan,
        'expected_duration': np.nan,
        'expected_ingress_time': np.nan,
        'expected_egress_time': np.nan,
        'overlap_duration': 0.0,
        'transit_fraction_observed': 0.0,
        'pre_ingress_points': 0,
        'in_transit_points': 0,
        'post_egress_points': 0,
        'covers_ingress': False,
        'covers_mid_transit': False,
        'covers_egress': False,
        'observed_segment': 'unknown',
        'success_label': 'unknown',
        'success_chance': np.nan,
        'expected_successful': False,
        'note': 'Could not evaluate expected transit coverage before UltraNest.',
    }

    if times.shape != flux_values.shape:
        base['note'] = 'Could not evaluate expected transit coverage because time and flux arrays were misaligned.'
        return base

    valid = np.isfinite(times) & np.isfinite(flux_values) & (flux_values > 0)
    if flux_errors is not None:
        flux_errors = np.asarray(flux_errors, dtype=float)
        if flux_errors.shape == flux_values.shape:
            valid &= np.isfinite(flux_errors) & (flux_errors > 0)

    if np.count_nonzero(valid) < 3:
        base['note'] = 'Could not evaluate expected transit coverage because too few finite light-curve points remain.'
        return base

    finite_times = np.sort(times[valid])
    observed_start = float(finite_times[0])
    observed_end = float(finite_times[-1])
    observed_span = float(observed_end - observed_start)
    mid_transit = _coverage_tmid_from_context(prior, tmid_search_summary=tmid_search_summary)
    duration = _coverage_duration_from_context(prior, duration_prior=duration_prior)
    base.update({
        'point_count': int(finite_times.size),
        'observed_start': observed_start,
        'observed_end': observed_end,
        'observed_span': observed_span,
        'expected_tmid': mid_transit,
        'expected_duration': duration,
    })

    if not np.isfinite(mid_transit):
        base['note'] = 'Could not evaluate expected transit coverage because no finite ephemeris Tmid was available.'
        return base
    if not np.isfinite(duration) or duration <= 0:
        base['note'] = 'Could not evaluate expected transit coverage because the expected transit duration is unavailable.'
        return base

    ingress_time = float(mid_transit - 0.5 * duration)
    egress_time = float(mid_transit + 0.5 * duration)
    in_transit_mask = valid & (times >= ingress_time) & (times <= egress_time)
    pre_mask = valid & (times < ingress_time)
    post_mask = valid & (times > egress_time)
    overlap_start = max(observed_start, ingress_time)
    overlap_end = min(observed_end, egress_time)
    overlap_duration = max(0.0, float(overlap_end - overlap_start))
    transit_fraction_observed = float(np.clip(overlap_duration / duration, 0.0, 1.0))
    covers_ingress = observed_start <= ingress_time <= observed_end
    covers_mid_transit = observed_start <= mid_transit <= observed_end
    covers_egress = observed_start <= egress_time <= observed_end
    observed_segment = expected_transit_observed_segment(
        observed_start,
        observed_end,
        ingress_time,
        mid_transit,
        egress_time,
    )
    success_label, success_chance = score_expected_transit_model_success(
        transit_fraction_observed,
        covers_ingress,
        covers_mid_transit,
        covers_egress,
        int(np.count_nonzero(pre_mask)),
        int(np.count_nonzero(post_mask)),
    )
    expected_successful = success_chance >= 0.5

    if expected_successful:
        note = (
            "The observed timestamps appear to contain enough of the expected transit window "
            "for a constrained nested fit."
        )
    elif transit_fraction_observed <= 0:
        note = (
            "The observed timestamps do not overlap the expected transit window; "
            "UltraNest is unlikely to recover a constrained transit solution."
        )
    elif int(np.count_nonzero(pre_mask)) == 0 or int(np.count_nonzero(post_mask)) == 0:
        note = (
            "The expected transit is not bracketed by out-of-transit data on both sides; "
            "UltraNest may chase partial-transit or baseline-degenerate solutions."
        )
    else:
        note = (
            "The expected transit is only partially observed; UltraNest may return broad or "
            "edge-hugging posteriors."
        )

    base.update({
        'valid': True,
        'expected_ingress_time': ingress_time,
        'expected_egress_time': egress_time,
        'overlap_duration': overlap_duration,
        'transit_fraction_observed': transit_fraction_observed,
        'pre_ingress_points': int(np.count_nonzero(pre_mask)),
        'in_transit_points': int(np.count_nonzero(in_transit_mask)),
        'post_egress_points': int(np.count_nonzero(post_mask)),
        'covers_ingress': bool(covers_ingress),
        'covers_mid_transit': bool(covers_mid_transit),
        'covers_egress': bool(covers_egress),
        'observed_segment': observed_segment,
        'success_label': success_label,
        'success_chance': float(success_chance),
        'expected_successful': bool(expected_successful),
        'note': note,
    })
    return base


def _format_minutes_from_days(days):
    try:
        value = float(days) * 24.0 * 60.0
    except (TypeError, ValueError):
        return "n/a"
    return "n/a" if not np.isfinite(value) else f"{value:.1f} min"


def log_expected_transit_coverage_assessment(assessment, indent="  "):
    if not isinstance(assessment, dict):
        return

    if not assessment.get('valid'):
        log_info(
            f"{indent}Warning: pre-UltraNest transit coverage assessment unavailable: "
            f"{assessment.get('note', 'unknown reason')}",
            warn=True,
        )
        return

    success_label = str(assessment.get('success_label', 'unknown')).upper()
    success_chance = coerce_finite_transit_qc_scalar(assessment.get('success_chance', np.nan))
    success_text = success_label
    if np.isfinite(success_chance):
        success_text = f"{success_label} (~{100.0 * float(success_chance):.0f}%)"

    warn = not bool(assessment.get('expected_successful', False))
    log_info(f"{indent}Pre-UltraNest transit coverage assessment:", warn=warn)
    log_info(
        f"{indent}  Data time range: {assessment['observed_start']:.8f} to "
        f"{assessment['observed_end']:.8f} BJD_TDB "
        f"({_format_minutes_from_days(assessment.get('observed_span'))}, "
        f"{assessment.get('point_count', 0)} point(s)).",
        warn=warn,
    )
    log_info(
        f"{indent}  Expected transit window: ingress {assessment['expected_ingress_time']:.8f}, "
        f"mid {assessment['expected_tmid']:.8f}, egress {assessment['expected_egress_time']:.8f} "
        f"BJD_TDB (duration {_format_minutes_from_days(assessment.get('expected_duration'))}).",
        warn=warn,
    )
    log_info(
        f"{indent}  Observed coverage: {assessment.get('observed_segment', 'unknown')}; "
        f"{100.0 * assessment.get('transit_fraction_observed', 0.0):.1f}% of the expected transit "
        f"window with {assessment.get('pre_ingress_points', 0)} pre-ingress, "
        f"{assessment.get('in_transit_points', 0)} in-transit, and "
        f"{assessment.get('post_egress_points', 0)} post-egress point(s).",
        warn=warn,
    )
    log_info(
        f"{indent}  Estimated fit success: {success_text}. {assessment.get('note', '')}",
        warn=warn,
    )


def extract_baseline_corrected_lightcurve_arrays(fit):
    times = np.asarray(getattr(fit, 'time', []), dtype=float)
    if times.ndim != 1 or times.size == 0:
        return None, None, None

    flux_values = None
    flux_source = None

    detrended = np.asarray(getattr(fit, 'detrended', []), dtype=float)
    if detrended.shape == times.shape:
        valid_detrended = np.isfinite(detrended) & (detrended > 0)
        if np.any(valid_detrended):
            flux_values = detrended.copy()
            flux_source = "current detrended light curve"

    if flux_values is None:
        data = np.asarray(getattr(fit, 'data', []), dtype=float)
        airmass_model = np.asarray(getattr(fit, 'airmass_model', []), dtype=float)
        if data.shape == times.shape and airmass_model.shape == times.shape:
            with np.errstate(divide='ignore', invalid='ignore'):
                flux_values = np.divide(data, airmass_model)
            flux_source = "current flux ratio divided by the fitted airmass/baseline model"
        elif data.shape == times.shape:
            flux_values = data.copy()
            flux_source = "current raw flux ratio"
        else:
            return None, None, None

    flux_errors = None
    detrended_errors = np.asarray(getattr(fit, 'detrendederr', []), dtype=float)
    if detrended_errors.shape == times.shape:
        valid_detrended_errors = np.isfinite(detrended_errors) & (detrended_errors > 0)
        if np.any(valid_detrended_errors):
            flux_errors = detrended_errors.copy()

    if flux_errors is None:
        data_errors = np.asarray(getattr(fit, 'dataerr', []), dtype=float)
        airmass_model = np.asarray(getattr(fit, 'airmass_model', []), dtype=float)
        if data_errors.shape == times.shape and airmass_model.shape == times.shape:
            with np.errstate(divide='ignore', invalid='ignore'):
                flux_errors = np.divide(data_errors, airmass_model)
        elif data_errors.shape == times.shape:
            flux_errors = data_errors.copy()
        else:
            flux_errors = np.full(times.shape, np.nan, dtype=float)

    return flux_values, flux_errors, flux_source


def prepare_final_fit_lightcurve_series(
    fit,
    depth_fraction=OUT_OF_TRANSIT_BASELINE_DEPTH_FRACTION,
):
    times = np.asarray(getattr(fit, 'time', []), dtype=float)
    if times.ndim != 1 or times.size == 0:
        return {
            'applied': False,
            'note': 'could not prepare a final-fit light curve because the current fit had no time samples.',
        }

    flux_values, flux_errors, flux_source = extract_baseline_corrected_lightcurve_arrays(fit)
    if flux_values is None:
        return {
            'applied': False,
            'note': 'could not derive a baseline-corrected light curve for final-fit preparation from the current fit.',
        }

    flux_values = np.asarray(flux_values, dtype=float)
    flux_errors = np.asarray(flux_errors, dtype=float)

    valid_flux = np.isfinite(flux_values) & (flux_values > 0)
    valid_errors = np.isfinite(flux_errors) & (flux_errors > 0)

    if np.count_nonzero(valid_flux) < LIGHTCURVE_MIN_VALID_POINTS:
        return {
            'applied': False,
            'note': 'not enough finite baseline-corrected flux points remained for final-fit preparation.',
        }

    coverage_summary = summarize_initial_fit_transit_coverage(
        times,
        fit,
        flux_values=flux_values,
        flux_errors=flux_errors if np.any(valid_errors) else None,
        depth_fraction=depth_fraction,
    )

    baseline_mask = valid_flux
    used_two_sided_oot = False
    pre_points = coverage_summary.get('pre_points', 0)
    post_points = coverage_summary.get('post_points', 0)

    if coverage_summary.get('valid') and coverage_summary.get('has_two_sided_oot'):
        candidate_baseline_mask = coverage_summary['oot_mask'] & valid_flux
        if np.count_nonzero(candidate_baseline_mask) >= LIGHTCURVE_MIN_VALID_POINTS:
            baseline_mask = candidate_baseline_mask
            used_two_sided_oot = True

    baseline_level, baseline_scatter = sigma_clipped_nanmedian(flux_values[baseline_mask])
    if not np.isfinite(baseline_level) or baseline_level <= 0:
        fallback_mask = valid_flux
        baseline_level, baseline_scatter = sigma_clipped_nanmedian(flux_values[fallback_mask])
        baseline_mask = fallback_mask

    if not np.isfinite(baseline_level) or baseline_level <= 0:
        return {
            'applied': False,
            'note': 'could not determine a positive baseline level for final-fit preparation.',
        }

    normalized_flux = flux_values / baseline_level
    normalized_unc = flux_errors / baseline_level

    if used_two_sided_oot:
        uncertainty_mask = baseline_mask & np.isfinite(normalized_unc) & (normalized_unc > 0)
        observed_scatter = np.nanstd(normalized_flux[baseline_mask])
        predicted_unc = np.nanmedian(normalized_unc[uncertainty_mask]) if np.any(uncertainty_mask) else np.nan
        if np.isfinite(observed_scatter) and observed_scatter > 0 and np.isfinite(predicted_unc) and predicted_unc > 0:
            normalized_unc *= observed_scatter / predicted_unc

    valid_normalized_unc = np.isfinite(normalized_unc) & (normalized_unc > 0)
    if not np.any(valid_normalized_unc):
        fallback_unc = baseline_scatter / baseline_level
        if not np.isfinite(fallback_unc) or fallback_unc <= 0:
            fallback_unc = np.nanstd(normalized_flux[baseline_mask])
        if not np.isfinite(fallback_unc) or fallback_unc <= 0:
            fallback_unc = np.finfo(float).eps
        normalized_unc = np.full(times.shape, fallback_unc, dtype=float)
    else:
        fallback_unc = np.nanmedian(normalized_unc[valid_normalized_unc])
        if not np.isfinite(fallback_unc) or fallback_unc <= 0:
            fallback_unc = np.nanstd(normalized_flux[baseline_mask])
        if not np.isfinite(fallback_unc) or fallback_unc <= 0:
            fallback_unc = np.finfo(float).eps
        normalized_unc[~valid_normalized_unc] = fallback_unc

    if used_two_sided_oot:
        note = (
            f"Prepared the final-fit input light curve from the {flux_source} and normalized it with "
            f"{pre_points} pre-ingress and {post_points} post-egress modeled out-of-transit point(s)."
        )
    else:
        note = (
            f"Prepared the final-fit input light curve from the {flux_source} and normalized it with a "
            f"sigma-clipped full-series baseline because the current fit only bracketed one side of transit "
            f"({pre_points} pre-ingress and {post_points} post-egress modeled out-of-transit point(s))."
        )

    return {
        'applied': True,
        'flux': normalized_flux,
        'unc': normalized_unc,
        'note': note,
        'source': flux_source,
        'coverage_summary': coverage_summary,
        'used_two_sided_oot': used_two_sided_oot,
        'baseline_level': float(baseline_level),
    }


def fit_airmass_baseline_parameters_on_out_of_transit(
    times,
    flux_values,
    flux_errors,
    airmass,
    fit,
    prior=None,
    bounds=None,
    depth_fraction=OUT_OF_TRANSIT_BASELINE_DEPTH_FRACTION,
):
    times = np.asarray(times, dtype=float)
    flux_values = np.asarray(flux_values, dtype=float)
    flux_errors = np.asarray(flux_errors, dtype=float)
    airmass = np.asarray(airmass, dtype=float)
    prior = {} if prior is None else dict(prior)
    bounds = {} if bounds is None else dict(bounds)

    base_result = {
        'applied': False,
        'note': 'out-of-transit baseline parameter fitting did not run.',
        'oot_mask': None,
        'pre_points': 0,
        'post_points': 0,
        'a0': np.nan,
        'a0_error': np.nan,
        'a2': prior.get('a2', 0.0),
        'a2_error': np.nan,
        'used_prior_ephemeris': False,
    }

    if not (times.shape == flux_values.shape == flux_errors.shape == airmass.shape):
        base_result['note'] = 'light-curve arrays could not be aligned for out-of-transit baseline fitting.'
        return base_result

    coverage_summary = summarize_initial_fit_transit_coverage(
        times,
        fit,
        flux_values=flux_values,
        flux_errors=flux_errors,
        depth_fraction=depth_fraction,
    )
    if not coverage_summary.get('valid') and prior:
        prior_coverage = summarize_prior_transit_coverage(
            times,
            prior,
            flux_values=flux_values,
            flux_errors=flux_errors,
        )
        if prior_coverage.get('valid'):
            coverage_summary = prior_coverage
            base_result['used_prior_ephemeris'] = True

    if not coverage_summary.get('valid'):
        base_result['note'] = coverage_summary.get(
            'note',
            'could not isolate out-of-transit points for baseline parameter fitting.',
        )
        return base_result

    oot_mask = np.asarray(coverage_summary.get('oot_mask'), dtype=bool)
    finite_mask = (
        oot_mask
        & np.isfinite(times)
        & np.isfinite(flux_values)
        & (flux_values > 0)
        & np.isfinite(flux_errors)
        & (flux_errors > 0)
        & np.isfinite(airmass)
    )
    point_count = int(np.count_nonzero(finite_mask))
    fit_a2 = 'a2' in bounds
    min_points = 3 if fit_a2 else 2
    base_result['pre_points'] = coverage_summary.get('pre_points', 0)
    base_result['post_points'] = coverage_summary.get('post_points', 0)

    if point_count < min_points:
        base_result['note'] = (
            f"only {point_count} finite out-of-transit point(s) were available; "
            f"need at least {min_points} to fit baseline parameters."
        )
        return base_result

    reference_airmass = transit_qc_airmass_reference(airmass)
    x = airmass[finite_mask] - reference_airmass
    y = flux_values[finite_mask]
    yerr = flux_errors[finite_mask]

    a0_bounds = bounds.get('a0') or bounds.get('a1') or [0.5, 1.5]
    try:
        a0_lower, a0_upper = np.asarray(a0_bounds, dtype=float).reshape(-1)[:2]
    except (TypeError, ValueError, IndexError):
        a0_lower, a0_upper = 0.5, 1.5
    if not np.isfinite(a0_lower) or a0_lower <= 0:
        a0_lower = max(np.nanmedian(y) * 0.5, np.finfo(float).eps)
    if not np.isfinite(a0_upper) or a0_upper <= a0_lower:
        a0_upper = max(np.nanmedian(y) * 1.5, a0_lower * 1.01)

    if fit_a2:
        try:
            a2_lower, a2_upper = np.asarray(bounds.get('a2'), dtype=float).reshape(-1)[:2]
        except (TypeError, ValueError, IndexError):
            a2_lower, a2_upper = -3.0, 3.0
        if not np.isfinite(a2_lower) or not np.isfinite(a2_upper) or a2_lower >= a2_upper:
            a2_lower, a2_upper = -3.0, 3.0
    else:
        a2_lower = a2_upper = float(prior.get('a2', 0.0) or 0.0)

    initial_a0 = float(np.clip(np.nanmedian(y), a0_lower, a0_upper))
    initial_a2 = float(prior.get('a2', 0.0) or 0.0)
    if fit_a2:
        initial_a2 = float(np.clip(initial_a2, a2_lower, a2_upper))

    if fit_a2:
        initial = np.array([np.log(initial_a0), initial_a2], dtype=float)
        lower_bounds = np.array([np.log(a0_lower), a2_lower], dtype=float)
        upper_bounds = np.array([np.log(a0_upper), a2_upper], dtype=float)
    else:
        initial = np.array([np.log(initial_a0)], dtype=float)
        lower_bounds = np.array([np.log(a0_lower)], dtype=float)
        upper_bounds = np.array([np.log(a0_upper)], dtype=float)

    def residuals(params):
        log_a0 = params[0]
        a2_value = params[1] if fit_a2 else initial_a2
        model = np.exp(log_a0) * np.exp(a2_value * x)
        return (y - model) / yerr

    try:
        result = least_squares(
            residuals,
            x0=initial,
            bounds=(lower_bounds, upper_bounds),
            jac='3-point',
            loss='linear',
        )
    except (ValueError, np.linalg.LinAlgError):
        base_result['note'] = 'weighted out-of-transit baseline parameter fit failed.'
        return base_result

    if not getattr(result, 'success', False) or not np.all(np.isfinite(result.x)):
        base_result['note'] = 'weighted out-of-transit baseline parameter fit did not converge.'
        return base_result

    log_a0 = float(result.x[0])
    a0 = float(np.exp(log_a0))
    a2 = float(result.x[1] if fit_a2 else initial_a2)
    jacobian = np.asarray(result.jac, dtype=float)
    residual_vector = np.asarray(result.fun, dtype=float)
    dof = max(1, residual_vector.size - result.x.size)
    reduced_chi2 = np.sum(residual_vector ** 2) / dof
    covariance = None
    if jacobian.ndim == 2 and jacobian.shape[0] >= jacobian.shape[1]:
        try:
            covariance = np.linalg.pinv(jacobian.T @ jacobian)
            covariance *= max(float(reduced_chi2), 1.0)
        except np.linalg.LinAlgError:
            covariance = None

    if covariance is not None and covariance.shape[0] >= 1:
        log_a0_error = float(np.sqrt(max(covariance[0, 0], 0.0)))
        a0_error = abs(a0) * log_a0_error
    else:
        a0_error = np.nan
    if covariance is not None and fit_a2 and covariance.shape[0] >= 2:
        a2_error = float(np.sqrt(max(covariance[1, 1], 0.0)))
    else:
        a2_error = 0.0 if not fit_a2 else np.nan

    if not np.isfinite(a0_error) or a0_error <= 0:
        a0_error = float(np.nanmedian(yerr))
    if fit_a2 and (not np.isfinite(a2_error) or a2_error <= 0):
        airmass_span_value = np.nanmax(x) - np.nanmin(x)
        if np.isfinite(airmass_span_value) and airmass_span_value > 0:
            a2_error = float(np.nanmedian(yerr / np.maximum(y, np.finfo(float).eps)) / airmass_span_value)
        else:
            a2_error = 0.0

    side_note = (
        f"{coverage_summary.get('pre_points', 0)} pre-ingress and "
        f"{coverage_summary.get('post_points', 0)} post-egress out-of-transit point(s)"
    )
    if base_result['used_prior_ephemeris']:
        side_note += " from the ephemeris-centered transit window"

    return {
        'applied': True,
        'note': (
            "Fitted a0"
            + (" and a2" if fit_a2 else "")
            + f" using only {side_note}; these baseline terms are fixed/profiled in the final transit fit."
        ),
        'oot_mask': finite_mask,
        'pre_points': coverage_summary.get('pre_points', 0),
        'post_points': coverage_summary.get('post_points', 0),
        'a0': a0,
        'a0_error': float(a0_error),
        'a2': a2,
        'a2_error': float(a2_error),
        'used_prior_ephemeris': base_result['used_prior_ephemeris'],
    }


def detrend_flux_on_out_of_transit_baseline(
    times,
    flux_values,
    flux_errors,
    fit,
    prior=None,
    depth_fraction=OUT_OF_TRANSIT_BASELINE_DEPTH_FRACTION,
):
    times = np.asarray(times, dtype=float)
    flux_values = np.asarray(flux_values, dtype=float)
    flux_errors = np.asarray(flux_errors, dtype=float)
    if flux_errors.shape != flux_values.shape:
        flux_errors = np.ones_like(flux_values, dtype=float)

    coverage_summary = summarize_initial_fit_transit_coverage(
        times,
        fit,
        flux_values=flux_values,
        flux_errors=flux_errors,
        depth_fraction=depth_fraction,
    )
    used_prior_coverage = False
    if prior is not None and (
        (not coverage_summary.get('valid'))
        or (not coverage_summary.get('has_two_sided_oot'))
    ):
        prior_coverage_summary = summarize_prior_transit_coverage(
            times,
            prior,
            flux_values=flux_values,
            flux_errors=flux_errors,
        )
        if prior_coverage_summary.get('valid') and prior_coverage_summary.get('has_two_sided_oot'):
            coverage_summary = prior_coverage_summary
            used_prior_coverage = True

    if not coverage_summary.get('valid'):
        return {
            'applied': False,
            'note': coverage_summary.get('note', 'could not isolate a transit window for baseline fitting.'),
        }

    if not coverage_summary.get('has_two_sided_oot'):
        return {
            'applied': False,
            'note': coverage_summary.get(
                'note',
                'need out-of-transit coverage on both sides of transit to fit a linear baseline.',
            ),
            'pre_points': coverage_summary.get('pre_points', 0),
            'post_points': coverage_summary.get('post_points', 0),
        }

    oot_mask = coverage_summary['oot_mask']
    mid_transit = coverage_summary['mid_transit']
    ingress_time = coverage_summary['ingress_time']
    egress_time = coverage_summary['egress_time']
    pre_points = coverage_summary['pre_points']
    post_points = coverage_summary['post_points']

    x = times[oot_mask] - mid_transit
    if np.allclose(x, x[0]):
        return {
            'applied': False,
            'note': 'out-of-transit timestamps do not span enough time to fit a line.',
            'pre_points': pre_points,
            'post_points': post_points,
        }

    design = np.column_stack((np.ones_like(x), x))
    oot_errors = flux_errors[oot_mask]
    weights = np.ones_like(x, dtype=float)
    valid_weights = np.isfinite(oot_errors) & (oot_errors > 0)
    if np.any(valid_weights):
        weights = np.zeros_like(x, dtype=float)
        weights[valid_weights] = 1.0 / (oot_errors[valid_weights] ** 2)
        if not np.any(weights > 0):
            weights = np.ones_like(x, dtype=float)

    sqrt_weights = np.sqrt(weights)
    try:
        coeffs, _, _, _ = np.linalg.lstsq(design * sqrt_weights[:, None], flux_values[oot_mask] * sqrt_weights, rcond=None)
    except np.linalg.LinAlgError:
        return {
            'applied': False,
            'note': 'linear out-of-transit baseline fit failed.',
            'pre_points': pre_points,
            'post_points': post_points,
        }

    intercept, slope = coeffs
    baseline = intercept + slope * (times - mid_transit)
    if not np.all(np.isfinite(baseline)) or np.any(baseline <= 0):
        return {
            'applied': False,
            'note': 'linear baseline prediction became non-physical for part of the light curve.',
            'pre_points': pre_points,
            'post_points': post_points,
        }

    return {
        'applied': True,
        'note': (
            (
                "Applied weighted linear out-of-transit baseline detrending using the "
                "ephemeris-centered transit window from the priors because the fitted transit "
                "window was one-sided. "
                if used_prior_coverage else
                "Applied weighted linear out-of-transit baseline detrending using "
            )
            + f"{pre_points} pre-ingress and {post_points} post-egress points."
        ),
        'flux': flux_values / baseline,
        'unc': flux_errors / baseline,
        'baseline': baseline,
        'slope': float(slope),
        'intercept': float(intercept),
        'pre_points': pre_points,
        'post_points': post_points,
        'ingress_time': ingress_time,
        'egress_time': egress_time,
        'used_prior_ephemeris': used_prior_coverage,
    }


def estimate_transit_duration_from_fit(fit):
    if fit is None:
        return np.nan

    for attribute_name in ('duration_expected', 'duration_measured'):
        duration = getattr(fit, attribute_name, np.nan)
        if np.isfinite(duration) and duration > 0:
            return float(duration)

    times = np.asarray(getattr(fit, 'time', []), dtype=float)
    transit_model = np.asarray(getattr(fit, 'transit', []), dtype=float)
    if times.shape != transit_model.shape or times.size == 0:
        return np.nan

    in_transit = np.isfinite(times) & np.isfinite(transit_model) & (transit_model < 1)
    if not np.any(in_transit):
        return np.nan

    transit_times = np.sort(times[in_transit])
    if transit_times.size == 1:
        sorted_times = np.sort(times[np.isfinite(times)])
        if sorted_times.size < 2:
            return np.nan
        cadence = np.nanmedian(np.diff(sorted_times))
        return float(cadence) if np.isfinite(cadence) and cadence > 0 else np.nan

    cadence = np.nanmedian(np.diff(np.sort(times[np.isfinite(times)])))
    if not np.isfinite(cadence) or cadence <= 0:
        cadence = 0.0
    duration = (transit_times[-1] - transit_times[0]) + cadence
    return float(duration) if np.isfinite(duration) and duration > 0 else np.nan


def build_nested_tmid_refinement_from_initial_fit(times, flux_values, flux_errors, prior, bounds, fit):
    original_tmid_bounds = clone_lightcurve_bounds(bounds).get('tmid')
    base_plan = {
        'applied': False,
        'note': 'Not needed; using the original nested-sampling Tmid bounds.',
        'prior': dict(prior),
        'bounds': clone_lightcurve_bounds(bounds),
        'original_tmid_bounds': original_tmid_bounds,
        'refined_tmid_bounds': original_tmid_bounds,
    }

    if fit is None or not hasattr(fit, 'parameters') or not isinstance(fit.parameters, dict):
        base_plan['note'] = "Skipped; the initial LM fit did not provide fitted parameters for nested-sampling refinement."
        return base_plan

    tmid = fit.parameters.get('tmid', np.nan)
    if not np.isfinite(tmid):
        base_plan['note'] = "Skipped; the initial LM fit did not return a finite Tmid."
        return base_plan

    coverage_summary = summarize_initial_fit_transit_coverage(
        times,
        fit,
        flux_values=flux_values,
        flux_errors=flux_errors,
    )
    if not coverage_summary.get('valid'):
        base_plan['note'] = (
            "Skipped; the initial LM fit did not provide a usable modeled transit window for nested-sampling refinement."
        )
        return base_plan

    if not coverage_summary.get('has_two_sided_oot'):
        pre_points = coverage_summary.get('pre_points', 0)
        post_points = coverage_summary.get('post_points', 0)
        base_plan['note'] = (
            "Skipped; the initial LM fit only captured one side of the modeled transit, "
            "so tightening the nested-sampling Tmid bounds would lock onto a partial-transit solution "
            f"({pre_points} pre-ingress and {post_points} post-egress out-of-transit point(s))."
        )
        return base_plan

    duration = estimate_transit_duration_from_fit(fit)
    period = fit.parameters.get('per', prior.get('per', np.nan))
    if not np.isfinite(duration) or duration <= 0:
        base_plan['note'] = "Skipped; the initial LM fit did not produce a measurable transit duration."
        return base_plan
    if np.isfinite(period) and duration >= period:
        base_plan['note'] = "Skipped; the initial LM fit returned a non-physical transit duration."
        return base_plan

    valid_times = np.asarray(times, dtype=float)
    valid_times = valid_times[np.isfinite(valid_times)]
    cadence = np.nan
    if valid_times.size > 1:
        cadence = np.nanmedian(np.diff(np.sort(valid_times)))

    half_width = FINAL_FIT_TMID_HALF_DURATION_MULTIPLIER * duration
    if np.isfinite(cadence) and cadence > 0:
        half_width = max(half_width, 3.0 * cadence)
    if not np.isfinite(half_width) or half_width <= 0:
        base_plan['note'] = "Skipped; the nested-sampling Tmid refinement half-width was not physical."
        return base_plan

    refined_lower = float(tmid - half_width)
    refined_upper = float(tmid + half_width)
    if original_tmid_bounds is not None and len(original_tmid_bounds) == 2:
        original_lower = float(original_tmid_bounds[0])
        original_upper = float(original_tmid_bounds[1])
        refined_lower = max(original_lower, refined_lower)
        refined_upper = min(original_upper, refined_upper)

    if not np.isfinite(refined_lower) or not np.isfinite(refined_upper) or refined_upper <= refined_lower:
        base_plan['note'] = "Skipped; the refined nested-sampling Tmid bounds collapsed to an invalid range."
        return base_plan

    refined_tmid_bounds = [refined_lower, refined_upper]
    base_plan['refined_tmid_bounds'] = refined_tmid_bounds
    if original_tmid_bounds is not None and np.allclose(
        np.asarray(original_tmid_bounds, dtype=float),
        np.asarray(refined_tmid_bounds, dtype=float),
        atol=1e-12,
        rtol=0.0,
    ):
        base_plan['note'] = (
            "Not needed; the initial LM fit already sat inside the original nested-sampling Tmid bounds."
        )
        return base_plan

    refined_prior = dict(prior)
    for key in ('rprs', 'ars', 'tmid', 'inc', 'a2'):
        if key in refined_prior and key in fit.parameters:
            refined_prior[key] = fit.parameters[key]

    refined_bounds = clone_lightcurve_bounds(bounds)
    refined_bounds['tmid'] = refined_tmid_bounds

    base_plan.update({
        'applied': True,
        'note': (
            "Using the initial LM fit to recenter nested-sampling Tmid bounds to "
            f"[{refined_lower:.6f}, {refined_upper:.6f}] around Tmid={tmid:.6f}."
        ),
        'prior': refined_prior,
        'bounds': refined_bounds,
    })
    return base_plan


def build_final_fit_prefit_refinement_plan(
    times,
    flux_values,
    flux_errors,
    airmass,
    prior,
    bounds,
    fit,
    jd_times=None,
    baseline_duration_multiplier=FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT,
):
    times = np.asarray(times, dtype=float)
    flux_values = np.asarray(flux_values, dtype=float)
    flux_errors = np.asarray(flux_errors, dtype=float)
    airmass = np.asarray(airmass, dtype=float)
    jd_array = None if jd_times is None else np.asarray(jd_times, dtype=float)

    original_tmid_bounds = clone_lightcurve_bounds(bounds).get('tmid')
    base_plan = {
        'applied': False,
        'note': None,
        'baseline_duration_multiplier': float(baseline_duration_multiplier),
        'duration': np.nan,
        'trimmed_pre_points': 0,
        'trimmed_post_points': 0,
        'original_point_count': int(times.shape[0]),
        'refined_point_count': int(times.shape[0]),
        'original_tmid_bounds': original_tmid_bounds,
        'refined_tmid_bounds': original_tmid_bounds,
        'times': times,
        'flux': flux_values,
        'unc': flux_errors,
        'airmass': airmass,
        'jd_times': jd_array,
        'prior': dict(prior),
        'bounds': clone_lightcurve_bounds(bounds),
    }

    if fit is None:
        base_plan['note'] = "Skipped; the initial nested fit did not return a solution."
        return base_plan

    duration = estimate_transit_duration_from_fit(fit)
    base_plan['duration'] = duration
    if not np.isfinite(duration) or duration <= 0:
        base_plan['note'] = "Skipped; the initial nested fit did not produce a measurable transit duration."
        return base_plan

    fit_parameters = getattr(fit, 'parameters', {})
    tmid = fit_parameters.get('tmid', prior.get('tmid', np.nan))
    if not np.isfinite(tmid):
        base_plan['note'] = "Skipped; the initial nested fit did not return a finite Tmid."
        return base_plan

    coverage_summary = summarize_initial_fit_transit_coverage(
        times,
        fit,
        flux_values=flux_values,
        flux_errors=flux_errors,
    )
    if coverage_summary.get('valid') and not coverage_summary.get('has_two_sided_oot'):
        pre_points = coverage_summary.get('pre_points', 0)
        post_points = coverage_summary.get('post_points', 0)
        base_plan['note'] = (
            "Skipped; the initial nested fit only captured one side of the modeled transit, "
            "so tightening the second-pass Tmid bounds would lock onto a partial-transit solution "
            f"({pre_points} pre-ingress and {post_points} post-egress out-of-transit point(s))."
        )
        return base_plan

    period = fit_parameters.get('per', prior.get('per', np.nan))
    if np.isfinite(period) and duration >= period:
        base_plan['note'] = "Skipped; the initial nested fit returned a non-physical transit duration."
        return base_plan

    half_duration = FINAL_FIT_TMID_HALF_DURATION_MULTIPLIER * duration
    if not np.isfinite(half_duration) or half_duration <= 0:
        base_plan['note'] = "Skipped; the initial nested fit returned an invalid transit duration."
        return base_plan

    keep_half_width = duration * (0.5 + baseline_duration_multiplier)
    lower_window = tmid - keep_half_width
    upper_window = tmid + keep_half_width
    keep_mask = np.isfinite(times) & (times >= lower_window) & (times <= upper_window)

    trimmed_pre_points = int(np.count_nonzero(np.isfinite(times) & (times < lower_window)))
    trimmed_post_points = int(np.count_nonzero(np.isfinite(times) & (times > upper_window)))
    kept_points = int(np.count_nonzero(keep_mask))
    min_required_points = max(LIGHTCURVE_MIN_VALID_POINTS, len(bounds) + 1)
    if kept_points < min_required_points:
        keep_mask = np.ones(times.shape[0], dtype=bool)
        kept_points = int(keep_mask.sum())
        trimmed_pre_points = 0
        trimmed_post_points = 0
        trim_note = (
            "kept the full light curve because trimming would leave too few points for a stable refit"
        )
    elif trimmed_pre_points or trimmed_post_points:
        trim_note = (
            f"trimmed {trimmed_pre_points} pre-ingress and {trimmed_post_points} post-egress point(s)"
        )
    else:
        trim_note = "kept the full light curve because no extra baseline points fell outside the target window"

    refined_lower = float(tmid - half_duration)
    refined_upper = float(tmid + half_duration)
    if original_tmid_bounds is not None and len(original_tmid_bounds) == 2:
        original_lower = float(original_tmid_bounds[0])
        original_upper = float(original_tmid_bounds[1])
        refined_lower = max(original_lower, refined_lower)
        refined_upper = min(original_upper, refined_upper)

    if not np.isfinite(refined_lower) or not np.isfinite(refined_upper) or refined_upper <= refined_lower:
        base_plan['note'] = "Skipped; the refined final-fit Tmid bounds collapsed to an invalid range."
        return base_plan

    refined_tmid_bounds = [refined_lower, refined_upper]
    base_plan['refined_tmid_bounds'] = refined_tmid_bounds

    bounds_changed = False
    if original_tmid_bounds is None:
        bounds_changed = True
    else:
        bounds_changed = not np.allclose(
            np.asarray(original_tmid_bounds, dtype=float),
            np.asarray(refined_tmid_bounds, dtype=float),
            atol=1e-12,
            rtol=0.0,
        )

    if not np.any(keep_mask):
        base_plan['note'] = "Skipped; no valid points remained inside the requested prefit window."
        return base_plan

    refined_times = times[keep_mask]
    refined_flux = flux_values[keep_mask]
    refined_unc = flux_errors[keep_mask]
    refined_airmass = airmass[keep_mask]
    refined_jd_times = None if jd_array is None else jd_array[keep_mask]

    refined_prior = dict(prior)
    if isinstance(fit_parameters, dict):
        for key in ('rprs', 'ars', 'tmid', 'inc', 'a2'):
            if key in refined_prior and key in fit_parameters:
                refined_prior[key] = fit_parameters[key]

    refined_bounds = clone_lightcurve_bounds(bounds)
    refined_bounds['tmid'] = refined_tmid_bounds

    base_plan.update({
        'applied': bool(trimmed_pre_points or trimmed_post_points or bounds_changed),
        'note': (
            "Using an initial nested fit to estimate a transit duration of "
            f"{duration:.6f} day(s), {trim_note}, and setting the second-pass "
            f"Tmid bounds to [{refined_tmid_bounds[0]:.6f}, {refined_tmid_bounds[1]:.6f}]."
        ),
        'trimmed_pre_points': trimmed_pre_points,
        'trimmed_post_points': trimmed_post_points,
        'refined_point_count': kept_points,
        'times': refined_times,
        'flux': refined_flux,
        'unc': refined_unc,
        'airmass': refined_airmass,
        'jd_times': refined_jd_times,
        'prior': refined_prior,
        'bounds': refined_bounds,
    })

    if not base_plan['applied']:
        base_plan['note'] = (
            "Not needed; the initial final-fit solution already used the desired baseline window "
            "and the Tmid bounds already matched the modeled transit duration."
        )

    return base_plan


def fit_final_lightcurve_with_oot_baseline_detrending(
    times,
    flux_values,
    flux_errors,
    airmass,
    prior,
    bounds,
    jd_times=None,
    skip_airmass_fit=False,
    airmass_skip_note=None,
    disable_vertical_flux_normalization=False,
    detrend_on_outoftransit_baseline=True,
    use_impactparameter_rather_than_inclination_to_fit=True,
    plot_time_range=None,
    baseline_duration_multiplier=FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT,
    expected_planet_dict=None,
    expected_tmid_search_summary=None,
    eebls_search_summary=None,
    duration_prior=None,
    extend_sparse_posterior_live_points=True,
    keep_ultranest_sampler_for_deferred_extension=False,
    fix_baseline_terms_for_final=True,
    pre_ultranest_coverage_assessment=None,
):
    if duration_prior is None and expected_planet_dict is not None:
        duration_prior = build_single_transit_duration_prior(expected_planet_dict)
    if pre_ultranest_coverage_assessment is None:
        pre_ultranest_coverage_assessment = build_expected_transit_coverage_assessment(
            times,
            prior,
            flux_values=flux_values,
            flux_errors=flux_errors,
            tmid_search_summary=expected_tmid_search_summary,
            duration_prior=duration_prior,
        )
    log_expected_transit_coverage_assessment(pre_ultranest_coverage_assessment)
    sparse_posterior_live_point_extension_enabled = should_use_sparse_posterior_live_point_retry(
        os.environ.get(
            SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_ENV,
            SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_DEFAULT,
        )
    )
    keep_ultranest_for_sparse_extension = (
        sparse_posterior_live_point_extension_enabled
        and (
            extend_sparse_posterior_live_points
            or keep_ultranest_sampler_for_deferred_extension
        )
    )

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times,
        flux_values,
        flux_errors,
        airmass,
        prior,
        bounds,
        jd_times=jd_times,
        use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
        duration_prior=duration_prior,
        keep_ultranest_sampler=keep_ultranest_for_sparse_extension,
        pre_ultranest_coverage_assessment=pre_ultranest_coverage_assessment,
    )
    fit = apply_plot_time_range(fit, times if plot_time_range is None else plot_time_range)
    annotate_airmass_fit(fit, airmass, skip_airmass_fit, note=airmass_skip_note)
    annotate_transit_qc_fit_context(
        fit,
        planet_dict=expected_planet_dict,
        tmid_search_summary=expected_tmid_search_summary,
        eebls_search_summary=eebls_search_summary,
    )
    annotate_pre_ultranest_transit_coverage(fit, pre_ultranest_coverage_assessment)

    effective_bounds = get_posterior_refit_final_bounds(fit, bounds)
    prefit_plan = build_final_fit_prefit_refinement_plan(
        times,
        flux_values,
        flux_errors,
        airmass,
        prior,
        effective_bounds,
        fit,
        jd_times=jd_times,
        baseline_duration_multiplier=baseline_duration_multiplier,
    )
    working_times = prefit_plan['times']
    working_flux = prefit_plan['flux']
    working_unc = prefit_plan['unc']
    working_airmass = prefit_plan['airmass']
    working_jd_times = prefit_plan['jd_times']
    working_prior = prefit_plan['prior']
    working_bounds = prefit_plan['bounds']

    if prefit_plan.get('applied'):
        log_info("Applying final-fit prefit refinement before the optional baseline detrending pass.")
        log_info(prefit_plan['note'])
        apply_vertical_flux_normalization_bound(
            working_prior,
            working_bounds,
            working_flux,
            disable_vertical_flux_normalization,
        )
        fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
            working_times,
            working_flux,
            working_unc,
            working_airmass,
            working_prior,
            working_bounds,
            jd_times=working_jd_times,
            use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
            duration_prior=duration_prior,
            keep_ultranest_sampler=keep_ultranest_for_sparse_extension,
            pre_ultranest_coverage_assessment=pre_ultranest_coverage_assessment,
        )
        fit = apply_plot_time_range(fit, working_times if plot_time_range is None else plot_time_range)
        annotate_airmass_fit(fit, working_airmass, skip_airmass_fit, note=airmass_skip_note)
        annotate_transit_qc_fit_context(
            fit,
            planet_dict=expected_planet_dict,
            tmid_search_summary=expected_tmid_search_summary,
            eebls_search_summary=eebls_search_summary,
        )
        annotate_pre_ultranest_transit_coverage(fit, pre_ultranest_coverage_assessment)

    working_bounds = get_posterior_refit_final_bounds(fit, working_bounds)

    annotate_final_fit_prefit_refinement(
        fit,
        prefit_plan.get('applied', False),
        note=prefit_plan.get('note'),
        baseline_duration_multiplier=baseline_duration_multiplier,
        duration=prefit_plan.get('duration'),
        original_point_count=prefit_plan.get('original_point_count'),
        refined_point_count=prefit_plan.get('refined_point_count'),
        trimmed_pre_points=prefit_plan.get('trimmed_pre_points', 0),
        trimmed_post_points=prefit_plan.get('trimmed_post_points', 0),
        original_tmid_bounds=prefit_plan.get('original_tmid_bounds'),
        refined_tmid_bounds=prefit_plan.get('refined_tmid_bounds'),
    )

    if detrend_on_outoftransit_baseline:
        baseline_parameter_result = fit_airmass_baseline_parameters_on_out_of_transit(
            working_times,
            working_flux,
            working_unc,
            working_airmass,
            fit,
            prior=working_prior,
            bounds=working_bounds,
        )
    else:
        baseline_parameter_result = {
            'applied': False,
            'note': 'Disabled with out-of-transit baseline detrending.',
            'pre_points': 0,
            'post_points': 0,
        }
    if baseline_parameter_result.get('applied') and not fix_baseline_terms_for_final:
        baseline_parameter_result = {
            'applied': False,
            'note': 'Deferred; pre-final UltraNest runs keep a0 and a2 as simultaneous fitted parameters.',
            'pre_points': baseline_parameter_result.get('pre_points', 0),
            'post_points': baseline_parameter_result.get('post_points', 0),
        }
    baseline_fit_mask = None
    baseline_fixed_errors = {}
    baseline_constrained_prior = dict(working_prior)
    baseline_constrained_bounds = clone_lightcurve_bounds(working_bounds)
    if baseline_parameter_result.get('applied'):
        log_info("Prepared out-of-transit airmass/baseline parameter constraints for the final transit refit.")
        log_info(baseline_parameter_result['note'])
        baseline_fit_mask = np.asarray(baseline_parameter_result['oot_mask'], dtype=bool)
        baseline_fixed_errors = {
            'a2': baseline_parameter_result.get('a2_error', 0.0),
        }
        baseline_constrained_prior['a0'] = baseline_parameter_result['a0']
        baseline_constrained_prior['a1'] = baseline_parameter_result['a0']
        baseline_constrained_prior['a2'] = baseline_parameter_result['a2']
        for key in ('rprs', 'ars', 'tmid', 'inc'):
            if key in baseline_constrained_prior and key in getattr(fit, 'parameters', {}):
                baseline_constrained_prior[key] = fit.parameters[key]
        baseline_constrained_bounds.pop('a0', None)
        baseline_constrained_bounds.pop('a1', None)
        baseline_constrained_bounds.pop('a2', None)
    else:
        annotate_out_of_transit_baseline_parameter_fit(
            fit,
            False,
            note=baseline_parameter_result.get('note'),
            pre_points=baseline_parameter_result.get('pre_points', 0),
            post_points=baseline_parameter_result.get('post_points', 0),
        )

    def run_oot_baseline_parameter_refit_if_needed(current_fit):
        if not baseline_parameter_result.get('applied'):
            return current_fit

        refit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
            working_times,
            working_flux,
            working_unc,
            working_airmass,
            baseline_constrained_prior,
            baseline_constrained_bounds,
            jd_times=working_jd_times,
            use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
            duration_prior=duration_prior,
            keep_ultranest_sampler=keep_ultranest_for_sparse_extension,
            baseline_fit_mask=baseline_fit_mask,
            fixed_parameter_errors=baseline_fixed_errors,
            fixed_flux_baseline=True,
            pre_ultranest_coverage_assessment=pre_ultranest_coverage_assessment,
        )
        refit = apply_plot_time_range(refit, working_times if plot_time_range is None else plot_time_range)
        annotate_airmass_fit(refit, working_airmass, skip_airmass_fit, note=airmass_skip_note)
        annotate_transit_qc_fit_context(
            refit,
            planet_dict=expected_planet_dict,
            tmid_search_summary=expected_tmid_search_summary,
            eebls_search_summary=eebls_search_summary,
        )
        annotate_pre_ultranest_transit_coverage(refit, pre_ultranest_coverage_assessment)
        annotate_final_fit_prefit_refinement(
            refit,
            prefit_plan.get('applied', False),
            note=prefit_plan.get('note'),
            baseline_duration_multiplier=baseline_duration_multiplier,
            duration=prefit_plan.get('duration'),
            original_point_count=prefit_plan.get('original_point_count'),
            refined_point_count=prefit_plan.get('refined_point_count'),
            trimmed_pre_points=prefit_plan.get('trimmed_pre_points', 0),
            trimmed_post_points=prefit_plan.get('trimmed_post_points', 0),
            original_tmid_bounds=prefit_plan.get('original_tmid_bounds'),
            refined_tmid_bounds=prefit_plan.get('refined_tmid_bounds'),
        )
        annotate_out_of_transit_baseline_parameter_fit(
            refit,
            True,
            note=baseline_parameter_result.get('note'),
            pre_points=baseline_parameter_result.get('pre_points', 0),
            post_points=baseline_parameter_result.get('post_points', 0),
            a0=baseline_parameter_result.get('a0'),
            a0_error=baseline_parameter_result.get('a0_error'),
            a2=baseline_parameter_result.get('a2'),
            a2_error=baseline_parameter_result.get('a2_error'),
        )
        return refit

    if not detrend_on_outoftransit_baseline:
        fit = run_oot_baseline_parameter_refit_if_needed(fit)
        annotate_out_of_transit_baseline_detrending(
            fit,
            False,
            note="Disabled; using the direct nested-sampling fit.",
        )
        annotate_transit_detection_qc(fit)
        if extend_sparse_posterior_live_points:
            fit = extend_sparse_posterior_live_points_if_needed(
                fit,
                enabled=sparse_posterior_live_point_extension_enabled,
            )
        annotate_transit_detection_qc(fit)
        return fit, working_flux, working_unc

    detrend_result = detrend_flux_on_out_of_transit_baseline(
        working_times,
        working_flux,
        working_unc,
        fit,
        prior=working_prior,
    )
    if not detrend_result.get('applied'):
        note = f"Skipped; {detrend_result.get('note', 'unable to fit an out-of-transit baseline.')}"
        log_info(f"Optional out-of-transit baseline detrending skipped: {detrend_result.get('note', 'unknown reason')}")
        fit = run_oot_baseline_parameter_refit_if_needed(fit)
        annotate_out_of_transit_baseline_detrending(
            fit,
            False,
            note=note,
            pre_points=detrend_result.get('pre_points', 0),
            post_points=detrend_result.get('post_points', 0),
        )
        annotate_transit_detection_qc(fit)
        if extend_sparse_posterior_live_points:
            fit = extend_sparse_posterior_live_points_if_needed(
                fit,
                enabled=sparse_posterior_live_point_extension_enabled,
            )
        annotate_transit_detection_qc(fit)
        return fit, working_flux, working_unc

    log_info("Applying optional out-of-transit linear baseline detrending and refitting final light curve.")
    log_info(detrend_result['note'])

    refit_prior = dict(working_prior)
    for key in ('rprs', 'ars', 'tmid', 'inc', 'a2'):
        if key in refit_prior and key in fit.parameters:
            refit_prior[key] = fit.parameters[key]

    refit_bounds = clone_lightcurve_bounds(working_bounds)
    apply_vertical_flux_normalization_bound(
        refit_prior,
        refit_bounds,
        detrend_result['flux'],
        disable_vertical_flux_normalization,
    )
    if baseline_parameter_result.get('applied'):
        refit_prior['a0'] = baseline_parameter_result['a0']
        refit_prior['a1'] = baseline_parameter_result['a0']
        refit_prior['a2'] = baseline_parameter_result['a2']
        refit_bounds.pop('a0', None)
        refit_bounds.pop('a1', None)
        refit_bounds.pop('a2', None)

    refit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        working_times,
        detrend_result['flux'],
        detrend_result['unc'],
        working_airmass,
        refit_prior,
        refit_bounds,
        jd_times=working_jd_times,
        use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
        duration_prior=duration_prior,
        keep_ultranest_sampler=keep_ultranest_for_sparse_extension,
        baseline_fit_mask=baseline_fit_mask,
        fixed_parameter_errors=baseline_fixed_errors,
        fixed_flux_baseline=bool(baseline_parameter_result.get('applied')),
        pre_ultranest_coverage_assessment=pre_ultranest_coverage_assessment,
    )
    refit = apply_plot_time_range(refit, working_times if plot_time_range is None else plot_time_range)
    annotate_airmass_fit(refit, working_airmass, skip_airmass_fit, note=airmass_skip_note)
    annotate_transit_qc_fit_context(
        refit,
        planet_dict=expected_planet_dict,
        tmid_search_summary=expected_tmid_search_summary,
        eebls_search_summary=eebls_search_summary,
    )
    annotate_pre_ultranest_transit_coverage(refit, pre_ultranest_coverage_assessment)
    annotate_final_fit_prefit_refinement(
        refit,
        prefit_plan.get('applied', False),
        note=prefit_plan.get('note'),
        baseline_duration_multiplier=baseline_duration_multiplier,
        duration=prefit_plan.get('duration'),
        original_point_count=prefit_plan.get('original_point_count'),
        refined_point_count=prefit_plan.get('refined_point_count'),
        trimmed_pre_points=prefit_plan.get('trimmed_pre_points', 0),
        trimmed_post_points=prefit_plan.get('trimmed_post_points', 0),
        original_tmid_bounds=prefit_plan.get('original_tmid_bounds'),
        refined_tmid_bounds=prefit_plan.get('refined_tmid_bounds'),
    )
    annotate_out_of_transit_baseline_detrending(
        refit,
        True,
        note=detrend_result['note'],
        slope=detrend_result['slope'],
        intercept=detrend_result['intercept'],
        pre_points=detrend_result['pre_points'],
        post_points=detrend_result['post_points'],
    )
    annotate_out_of_transit_baseline_parameter_fit(
        refit,
        bool(baseline_parameter_result.get('applied')),
        note=baseline_parameter_result.get('note'),
        pre_points=baseline_parameter_result.get('pre_points', 0),
        post_points=baseline_parameter_result.get('post_points', 0),
        a0=baseline_parameter_result.get('a0'),
        a0_error=baseline_parameter_result.get('a0_error'),
        a2=baseline_parameter_result.get('a2'),
        a2_error=baseline_parameter_result.get('a2_error'),
    )
    annotate_transit_detection_qc(refit)
    if extend_sparse_posterior_live_points:
        refit = extend_sparse_posterior_live_points_if_needed(
            refit,
            enabled=sparse_posterior_live_point_extension_enabled,
        )
    annotate_transit_detection_qc(refit)
    return refit, detrend_result['flux'], detrend_result['unc']


def psf_sigma_from_fit(psf_row, fallback_sigma=np.nan):
    try:
        sigx = float(psf_row[3])
        sigy = float(psf_row[4])
        sigma = 0.5 * (sigx + sigy)
    except (IndexError, TypeError, ValueError):
        sigma = np.nan

    if np.isfinite(sigma) and sigma > 0:
        return float(sigma)

    if np.isfinite(fallback_sigma) and fallback_sigma > 0:
        return float(fallback_sigma)

    return np.nan


def representative_psf_sigma(psf_rows, fallback_sigma=np.nan):
    try:
        sigmas = np.asarray(psf_rows[:, 3], dtype=float) + np.asarray(psf_rows[:, 4], dtype=float)
    except (IndexError, TypeError, ValueError):
        sigmas = np.array([], dtype=float)

    if sigmas.size:
        sigmas *= 0.5
        sigmas[~np.isfinite(sigmas) | (sigmas <= 0)] = np.nan
        center, _ = sigma_clipped_nanmedian(sigmas)
        if np.isfinite(center) and center > 0:
            return float(center)

    if np.isfinite(fallback_sigma) and fallback_sigma > 0:
        return float(fallback_sigma)

    return np.nan


def summarize_adaptive_aperture_usage(psf_rows, aperture_scale, annulus_scale, fallback_sigma=np.nan):
    try:
        aperture_scale = float(aperture_scale)
        annulus_scale = float(annulus_scale)
    except (TypeError, ValueError):
        return None

    if not np.isfinite(aperture_scale) or not np.isfinite(annulus_scale):
        return None

    rows = np.asarray(psf_rows)
    if rows.ndim != 2 or rows.shape[0] == 0:
        return None

    frame_sigma = np.array(
        [psf_sigma_from_fit(row, fallback_sigma=fallback_sigma) for row in rows],
        dtype=float,
    )
    frame_sigma[~np.isfinite(frame_sigma) | (frame_sigma <= 0)] = np.nan

    aperture_series = aperture_scale * frame_sigma
    annulus_series = annulus_scale * frame_sigma
    fwhm_series = GAUSSIAN_SIGMA_TO_FWHM * frame_sigma

    geometry_rows = [
        resolve_sky_annulus_geometry(aperture_radius, annulus_width, psf_sigma=sigma)
        if np.isfinite(aperture_radius) and np.isfinite(annulus_width) and np.isfinite(sigma)
        else None
        for aperture_radius, annulus_width, sigma in zip(aperture_series, annulus_series, frame_sigma)
    ]
    sky_inner_series = np.array(
        [np.nan if geometry is None else geometry['inner_radius'] for geometry in geometry_rows],
        dtype=float,
    )
    sky_outer_series = np.array(
        [np.nan if geometry is None else geometry['outer_radius'] for geometry in geometry_rows],
        dtype=float,
    )
    sky_pixel_series = np.array(
        [np.nan if geometry is None else geometry['effective_sky_pixels'] for geometry in geometry_rows],
        dtype=float,
    )

    if not np.any(np.isfinite(aperture_series)) or not np.any(np.isfinite(annulus_series)):
        return None

    return {
        'aperture_sigma': aperture_scale,
        'annulus_sigma': annulus_scale,
        'frame_sigma': frame_sigma,
        'fwhm_series': fwhm_series,
        'aperture_series': aperture_series,
        'annulus_series': annulus_series,
        'sky_inner_series': sky_inner_series,
        'sky_outer_series': sky_outer_series,
        'sky_pixel_series': sky_pixel_series,
        'aperture_median': float(np.nanmedian(aperture_series)),
        'aperture_std': float(np.nanstd(aperture_series)),
        'aperture_min': float(np.nanmin(aperture_series)),
        'aperture_max': float(np.nanmax(aperture_series)),
        'annulus_median': float(np.nanmedian(annulus_series)),
        'annulus_std': float(np.nanstd(annulus_series)),
        'annulus_min': float(np.nanmin(annulus_series)),
        'annulus_max': float(np.nanmax(annulus_series)),
    }


def update_photometry_adaptive_summary(photometry_info, use_adaptive_apertures, aperture_values, annulus_values,
                                       psf_rows, fallback_sigma=np.nan):
    photometry_info['adaptive_summary'] = None

    if (not use_adaptive_apertures) or photometry_info.get('min_aperture') in (None, 0):
        return None

    a_idx = photometry_info.get('aperture_index')
    an_idx = photometry_info.get('annulus_index')
    if a_idx is None or an_idx is None or aperture_values is None or annulus_values is None:
        return None

    aperture_grid = np.asarray(aperture_values, dtype=float)
    annulus_grid = np.asarray(annulus_values, dtype=float)
    if a_idx >= aperture_grid.size or an_idx >= annulus_grid.size:
        return None

    photometry_info['adaptive_summary'] = summarize_adaptive_aperture_usage(
        psf_rows,
        aperture_grid[a_idx],
        annulus_grid[an_idx],
        fallback_sigma=fallback_sigma,
    )
    return photometry_info['adaptive_summary']


def reported_photometry_aperture_radii(photometry_info):
    adaptive_summary = photometry_info.get('adaptive_summary')
    if adaptive_summary is None:
        return photometry_info.get('min_aperture'), photometry_info.get('min_annulus')

    aperture = adaptive_summary['aperture_median']
    if photometry_info.get('min_aperture') is not None and photometry_info['min_aperture'] < 0:
        aperture = -aperture
    return aperture, adaptive_summary['annulus_median']


def build_observing_background_series(psf_data, aper_data, photometry_info, comp_star_count):
    use_aperture_background = photometry_info.get('min_aperture') != 0
    a_idx = photometry_info.get('aperture_index')
    an_idx = photometry_info.get('annulus_index')

    if use_aperture_background and aper_data is not None and a_idx is not None and an_idx is not None:
        background_series = {
            'target': np.asarray(aper_data['target_bg'][:, a_idx, an_idx], dtype=float),
        }
        for comp_idx in range(comp_star_count):
            ckey = f"comp{comp_idx + 1}"
            bg_key = f"{ckey}_bg"
            if bg_key in aper_data:
                background_series[ckey] = np.asarray(aper_data[bg_key][:, a_idx, an_idx], dtype=float)
        return background_series

    background_series = {
        'target': np.asarray(psf_data['target'][:, 6], dtype=float),
    }
    for comp_idx in range(comp_star_count):
        ckey = f"comp{comp_idx + 1}"
        if ckey in psf_data:
            background_series[ckey] = np.asarray(psf_data[ckey][:, 6], dtype=float)
    return background_series


def resolve_frame_aperture_radii(apertures, annuli, adaptive_apertures=False, frame_sigma=np.nan,
                                 fallback_sigma=np.nan):
    aperture_values = np.asarray(apertures, dtype=float)
    annulus_values = np.asarray(annuli, dtype=float)

    if not adaptive_apertures:
        return aperture_values, annulus_values

    sigma_to_use = float(frame_sigma) if np.isfinite(frame_sigma) and frame_sigma > 0 else np.nan
    if (not np.isfinite(sigma_to_use) or sigma_to_use <= 0) and np.isfinite(fallback_sigma) and fallback_sigma > 0:
        sigma_to_use = float(fallback_sigma)
    if not np.isfinite(sigma_to_use) or sigma_to_use <= 0:
        sigma_to_use = 1.0

    return aperture_values * sigma_to_use, annulus_values * sigma_to_use


# Initialze plate status log
plateStatus = PlateStatus(log_info)

def sigma_clip(ogdata, sigma=3, dt=21, po=2, times=None):
    values = np.asarray(ogdata, dtype=float)
    nanmask = np.isnan(values)
    valid_mask = ~nanmask
    valid_indices = np.flatnonzero(valid_mask)

    if not (po < dt <= valid_indices.size):
        return nanmask

    segment_ranges = []
    if times is not None:
        time_values = np.asarray(times, dtype=float)
        if time_values.shape == values.shape:
            valid_times = time_values[valid_mask]
            finite_valid_times = np.isfinite(valid_times)
            if np.all(finite_valid_times):
                cadence = np.nanmedian(np.diff(valid_times)) if valid_times.size > 1 else np.nan
                if np.isfinite(cadence) and cadence > 0:
                    gap_threshold = max(
                        5.0 * cadence,
                        0.25 * int(dt) * cadence,
                    )
                    local_start = 0
                    for local_index, gap in enumerate(np.diff(valid_times), start=1):
                        if gap > gap_threshold:
                            segment_ranges.append((local_start, local_index))
                            local_start = local_index
                    segment_ranges.append((local_start, valid_times.size))

    if not segment_ranges:
        segment_ranges = [(0, valid_indices.size)]

    clipped_mask = np.zeros(valid_indices.size, dtype=bool)
    for start, stop in segment_ranges:
        local_values = values[valid_indices[start:stop]]
        if not (po < dt <= local_values.size):
            continue

        mdata = savgol_filter(local_values, window_length=dt, polyorder=po)
        # mdata = median_filter(local_values, dt)
        res = local_values - mdata
        if res.size == 0:
            continue
        # Vectorized bootstrap estimate avoids Python-loop overhead in tight runs.
        sample_size = min(25, res.size)
        bootstrap_samples = np.random.choice(res, size=(100, sample_size), replace=True)
        std = bn.nanmedian(bn.nanstd(bootstrap_samples, axis=1))
        # std = np.nanstd(res) # biased from large outliers
        if not np.isfinite(std) or std <= 0:
            continue
        sigmask = np.abs(res) > sigma * std
        clipped_mask[start:stop] = sigmask

    nanmask[valid_indices] = clipped_mask
    return nanmask


def adaptive_aperture_outlier_mask(aperture_series, annulus_series=None, sigma=4.5, window=15, polyorder=2):
    aperture_series = np.asarray(aperture_series, dtype=float)
    combined_mask = _adaptive_series_outlier_mask(
        aperture_series,
        sigma=sigma,
        window=window,
        polyorder=polyorder,
    )

    if annulus_series is None:
        return combined_mask

    annulus_series = np.asarray(annulus_series, dtype=float)
    annulus_mask = _adaptive_series_outlier_mask(
        annulus_series,
        sigma=sigma,
        window=window,
        polyorder=polyorder,
    )
    return combined_mask | annulus_mask


def _adaptive_series_outlier_mask(series, sigma=4.5, window=15, polyorder=2):
    values = np.asarray(series, dtype=float)
    nanmask = ~np.isfinite(values)
    valid_indices = np.flatnonzero(~nanmask)
    if valid_indices.size < max(polyorder + 3, 7):
        return nanmask

    valid_values = values[valid_indices]
    window_length = min(int(window), valid_values.size)
    if window_length % 2 == 0:
        window_length -= 1

    if window_length >= polyorder + 2:
        trend = savgol_filter(valid_values, window_length=window_length, polyorder=polyorder, mode='interp')
        residuals = valid_values - trend
        scatter = robust_scatter(residuals)
        center = trend
    else:
        scatter = np.nan
        center = np.full(valid_values.shape, bn.nanmedian(valid_values))

    if not np.isfinite(scatter) or scatter <= 0:
        center = np.full(valid_values.shape, bn.nanmedian(valid_values))
        residuals = valid_values - center
        scatter = robust_scatter(residuals)
        if not np.isfinite(scatter) or scatter <= 0:
            return nanmask

    local_mask = np.abs(valid_values - center) > sigma * scatter
    outlier_mask = nanmask.copy()
    outlier_mask[valid_indices] = local_mask
    return outlier_mask


def robust_scatter(data):
    values = np.asarray(data, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size < 2:
        return np.nan

    center = bn.nanmedian(finite)
    mad = bn.nanmedian(np.abs(finite - center))
    if np.isfinite(mad) and mad > 0:
        return 1.4826 * mad

    scatter = bn.nanstd(finite)
    if np.isfinite(scatter) and scatter > 0:
        return scatter

    return np.nan


def phase_bin_sigma_clip(values, phase, sigma=3, bins=10, min_points=5, max_iters=3):
    values = np.asarray(values, dtype=float)
    phase = np.asarray(phase, dtype=float)
    nanmask = ~np.isfinite(values) | ~np.isfinite(phase)
    valid_indices = np.flatnonzero(~nanmask)

    if valid_indices.size < max(min_points, 3):
        return nanmask

    phase_valid = phase[valid_indices]
    min_phase = np.nanmin(phase_valid)
    max_phase = np.nanmax(phase_valid)
    if not np.isfinite(min_phase) or not np.isfinite(max_phase) or min_phase == max_phase:
        return nanmask

    bin_count = max(1, int(bins))
    edges = np.linspace(min_phase, max_phase, bin_count + 1)
    bin_ids = np.searchsorted(edges[1:-1], phase_valid, side='right')
    keep_mask = np.ones(valid_indices.size, dtype=bool)
    values_valid = values[valid_indices]

    for bin_id in range(bin_count):
        local_positions = np.flatnonzero(bin_ids == bin_id)
        if local_positions.size < min_points:
            continue

        local_keep = np.ones(local_positions.size, dtype=bool)
        for _ in range(max_iters):
            candidate_values = values_valid[local_positions][local_keep]
            if candidate_values.size < min_points:
                break

            center = bn.nanmedian(candidate_values)
            scatter = robust_scatter(candidate_values)
            if not np.isfinite(scatter) or scatter <= 0:
                break

            within_limits = np.abs(candidate_values - center) <= sigma * scatter
            if np.all(within_limits):
                break

            local_keep[np.flatnonzero(local_keep)[~within_limits]] = False

        keep_mask[local_positions] &= local_keep

    nanmask[valid_indices] = ~keep_mask
    return nanmask


def apply_lightcurve_mask(lightcurve, mask, sort_index=None):
    if lightcurve is None:
        return

    mask = np.asarray(mask, dtype=bool)
    target_length = mask.shape[0]
    if sort_index is not None:
        sort_index = np.asarray(sort_index)
        target_length = sort_index.shape[0]

    array_attrs = (
        'time',
        'data',
        'airmass',
        'transit',
        'jd_times',
        'phase',
        'residuals',
        'model',
        'detrended',
        'detrendederr',
        'dataerr',
        'airmass_model',
        'wf',
    )

    for attr in array_attrs:
        if not hasattr(lightcurve, attr):
            continue

        values = getattr(lightcurve, attr)
        if values is None:
            continue

        array_values = np.asarray(values)
        if array_values.ndim == 0 or array_values.shape[0] != target_length:
            continue

        if sort_index is not None:
            array_values = array_values[sort_index]
        setattr(lightcurve, attr, array_values[mask])


def apply_plot_time_range(lightcurve, time_values):
    if lightcurve is None:
        return lightcurve

    values = np.asarray(time_values, dtype=float).reshape(-1)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return lightcurve

    lightcurve.plot_time_range = (float(np.min(finite)), float(np.max(finite)))
    updater = getattr(lightcurve, "_update_plot_geometry", None)
    if callable(updater):
        updater()

    return lightcurve


def build_time_rejection_diagnostic(stage, times, keep_mask, note=None):
    times = np.asarray(times, dtype=float).reshape(-1)
    keep_mask = np.asarray(keep_mask, dtype=bool).reshape(-1)
    if times.shape[0] != keep_mask.shape[0]:
        return None

    finite_mask = np.isfinite(times)
    input_point_count = int(np.count_nonzero(finite_mask))
    dropped_times = np.asarray(times[finite_mask & ~keep_mask], dtype=float)
    kept_point_count = int(np.count_nonzero(finite_mask & keep_mask))
    dropped_point_count = int(dropped_times.size)

    cadence = np.nan
    finite_times = np.sort(times[finite_mask])
    if finite_times.size > 1:
        cadence = np.nanmedian(np.diff(finite_times))

    dropped_ranges = []
    if dropped_times.size:
        dropped_times = np.sort(dropped_times)
        gap_threshold = np.inf
        if np.isfinite(cadence) and cadence > 0:
            gap_threshold = max(
                TIME_REJECTION_GROUP_GAP_CADENCE_MULTIPLIER * cadence,
                np.finfo(float).eps,
            )

        range_start = float(dropped_times[0])
        range_end = float(dropped_times[0])
        range_count = 1
        for current_time in dropped_times[1:]:
            current_time = float(current_time)
            if np.isfinite(gap_threshold) and (current_time - range_end) <= gap_threshold:
                range_end = current_time
                range_count += 1
                continue

            dropped_ranges.append({
                'start': range_start,
                'end': range_end,
                'count': int(range_count),
            })
            range_start = current_time
            range_end = current_time
            range_count = 1

        dropped_ranges.append({
            'start': range_start,
            'end': range_end,
            'count': int(range_count),
        })

    return {
        'stage': stage,
        'note': note,
        'input_point_count': input_point_count,
        'kept_point_count': kept_point_count,
        'dropped_point_count': dropped_point_count,
        'cadence': float(cadence) if np.isfinite(cadence) else np.nan,
        'dropped_ranges': dropped_ranges,
        'first_dropped_time': (float(dropped_times[0]) if dropped_times.size else np.nan),
        'last_dropped_time': (float(dropped_times[-1]) if dropped_times.size else np.nan),
    }


def format_time_rejection_diagnostic(diagnostic, max_ranges=TIME_REJECTION_RANGE_DISPLAY_LIMIT):
    if not diagnostic:
        return None

    dropped_ranges = diagnostic.get('dropped_ranges') or []
    if not dropped_ranges:
        return (
            f"{diagnostic.get('stage', 'frame filter')}: removed 0/"
            f"{diagnostic.get('input_point_count', 0)} frame(s)."
        )

    display_ranges = dropped_ranges[:max_ranges]
    range_parts = []
    for range_summary in display_ranges:
        start = float(range_summary['start'])
        end = float(range_summary['end'])
        count = int(range_summary['count'])
        if count <= 1 or np.isclose(start, end):
            range_parts.append(f"{start:.8f} ({count} frame)")
        else:
            frame_label = "frame" if count == 1 else "frames"
            range_parts.append(f"{start:.8f} to {end:.8f} ({count} {frame_label})")

    if len(dropped_ranges) > len(display_ranges):
        remaining = len(dropped_ranges) - len(display_ranges)
        range_parts.append(f"... {remaining} more range(s)")

    message = (
        f"{diagnostic.get('stage', 'frame filter')}: removed "
        f"{diagnostic.get('dropped_point_count', 0)}/{diagnostic.get('input_point_count', 0)} frame(s); "
        f"BJD range(s): {'; '.join(range_parts)}."
    )
    note = diagnostic.get('note')
    if note:
        message += f" {note}"
    return message


def log_lightcurve_filter_diagnostics(diagnostics, header="Lightcurve frame rejection diagnostics", only_removed=True):
    normalized = [diagnostic for diagnostic in (diagnostics or []) if diagnostic]
    if only_removed:
        normalized = [diagnostic for diagnostic in normalized if diagnostic.get('dropped_point_count', 0) > 0]
    if not normalized:
        return

    log_info(f"\n{header}:")
    for diagnostic in normalized:
        diagnostic_text = format_time_rejection_diagnostic(diagnostic)
        if diagnostic_text:
            log_info(f"  {diagnostic_text}")


def exp_offset(hdr, time_unit, exp):
    """Returns exposure offset (in days) of more than 0 if headers reveals
    the time was estimated at the start of the exposure rather than the middle
    """
    if 'start' in hdr.comments[time_unit]:
        return exp / (2.0 * 60.0 * 60.0 * 24.0)
    return 0.0


def ut_date(hdr, time_unit, exp):
    """Converts the Gregorian Date to Julian Date from the header and returns it
    along with the exposure offset
    """
    if time_unit == 'DATE-OBS':
        greg_date = hdr[time_unit] if 'T' in hdr[time_unit] else f"{hdr[time_unit]}T{hdr['TIME-OBS']}"
    else:
        greg_date = hdr[time_unit]

    dt = dup.parse(greg_date)
    atime = Time(dt)

    julian_time = atime.jd
    offset = exp_offset(hdr, time_unit, exp)

    return julian_time + offset


def julian_date(hdr, time_unit, exp):
    """Returns Julian Date from the header along with the exposure offset.
    If the image is taken from MicroObservatory (MJD-OBS),
    add a timing offset (2400000.5) due to being less precise
    """
    time_offset = 2400000.5 if time_unit == 'MJD-OBS' else 0.0

    julian_time = float(hdr[time_unit]) + time_offset
    offset = exp_offset(hdr, time_unit, exp)

    return julian_time + offset

def get_exp_time(hdr):
    exp_list = ["EXPTIME", "EXPOSURE", "EXP"]
    exp_time = next((exptime for exptime in exp_list if exptime in hdr), None)
    return hdr[exp_time] if exp_time is not None else 0.0

def img_time_jd(hdr):
    """Converts time from the header file to the Julian Date (JD, if needed)
    and adds an exposure offset (if needed)

    Parameters
    ----------
    hdr : astropy.io.fits.header.Header
        A header file that includes the time from when the image was taken
    Returns
    -------
    float
        Time of when the image was taken in the JD with exposure offset
    """
    time_list = ['UT-OBS', 'JULIAN', 'MJD-OBS', 'DATE-OBS']

    exp = get_exp_time(hdr)
    hdr_time = next((time_unit for time_unit in time_list if time_unit in hdr), None)

    if hdr_time == 'MJD_OBS':
        hdr_time = hdr_time if "epoch" not in hdr.comments[hdr_time] else 'DATE-OBS'

    if hdr_time in ['UT-OBS', 'DATE-OBS']:
        return ut_date(hdr, hdr_time, exp)
    return julian_date(hdr, hdr_time, exp)


def img_time_bjd_tdb(hdr, p_dict, info_dict):
    """Converts time from the header file to BJD-TDB time (if needed)
    and adds an exposure offset (if needed)

    Parameters
    ----------
    hdr : astropy.io.fits.header.Header
        A header file that includes the time from when the image was taken
    p_dict: planetary settings dictionary
    info_dict: observatory settings dictionary

    Returns
    -------
    float
        Time of when the image was taken in BJD-TDB with exposure offset
    """
    # Check for BJD time first (preference)
    time_list = ['BJD_TDB', 'BJD_TBD', 'BJD']
    exp = get_exp_time(hdr)

    hdr_time = next((time for time in time_list if time in hdr), None)
    # Not found, get julian date

    if hdr_time is None:
        time_list = ['UT-OBS', 'JULIAN', 'MJD-OBS', 'DATE-OBS']
        hdr_time = next((time for time in time_list if time in hdr), None)
        if hdr_time == 'MJD_OBS':
            hdr_time = hdr_time if "epoch" not in hdr.comments[hdr_time] else 'DATE-OBS'
        if hdr_time in ['UT-OBS', 'DATE-OBS']:
            jd_time = ut_date(hdr, hdr_time, exp)
        else:
            jd_time = julian_date(hdr, hdr_time, exp)
        # And convert to BJD_TDB
        bjd_time = convert_jd_to_bjd([jd_time], p_dict, info_dict)[0]
    else:   # Else, already BJD - convert and adjust for exposure
        bjd_time = julian_date(hdr, hdr_time, exp)
    return bjd_time

def air_mass(hdr, ra, dec, lat, long, elevation, time):
    """Scrapes or calculates the airmass at the time of when the image was taken.
    Airmass(X): X = sec(z), z = secant of the zenith angle (angle between zenith and star)

    Parameters
    ----------
    hdr : astropy.io.fits.header.Header
        A header file that may include the airmass or altitude from when the image was taken
    ra : float
        Right Ascension
    dec : float
        Declination
    lat : float
        Latitude
    long : float
        Longitude
    elevation : float
        Elevation/Altitude

    Returns
    -------
    float
        Airmass value
    """
    if 'AIRMASS' in hdr:
        am = float(hdr['AIRMASS'])
    elif 'TELALT' in hdr:
        alt = float(hdr['TELALT'])
        cos_am = np.cos((np.pi / 180) * (90.0 - alt))
        am = 1 / cos_am
    else:
        pointing = SkyCoord(f"{ra} {dec}", unit=(u.deg, u.deg), frame='icrs')

        location = EarthLocation.from_geodetic(lat=lat * u.deg, lon=long * u.deg, height=elevation)
        time = Time(time, format='jd', scale='utc', location=location)
        point_altaz = pointing.transform_to(AltAz(obstime=time, location=location))
        am = float(point_altaz.secz)
    return am


def flux_conversion(fluxes, errors, flux_format):
    """Converting differential magnitudes to fluxes and calculating its errors
    """
    conv = 1000.0 if flux_format == 'millimagnitude' else 1.0

    pos_err = 10.0 ** (-0.4 * ((fluxes + errors) / conv))
    neg_err = 10.0 ** (-0.4 * ((fluxes - errors) / conv))
    fluxes = 10.0 ** (-0.4 * (fluxes / conv))

    pos_err_dist = abs(pos_err - fluxes)
    neg_err_dist = abs(neg_err - fluxes)
    mean_errors = (pos_err_dist * neg_err_dist) ** 0.5

    return fluxes, mean_errors


# Check for difference between NEA and initialization file
def check_parameters(init_parameters, parameters):
    different = False
    uncert = 1 / 36

    for key, value in parameters.items():
        if key in ['ra', 'dec'] and init_parameters[key]:
            if not parameters[key] - uncert <= init_parameters[key] <= parameters[key] + uncert:
                different = True
                break
            continue
        if value != init_parameters[key]:
            different = True
            break

    if different:
        log_info("\nDifference(s) found between initialization file parameters and "
                 "those scraped by EXOTIC from the NASA Exoplanet Archive."
                 "\nWould you like:"
                 "\n  (1) EXOTIC to adopt of all of your defined parameters or"
                 "\n  (2) to review the ones scraped from the Archive that differ?")
        opt = user_input("Enter 1 or 2: ", type_=int, values=[1, 2])

        if opt == 2:
            return True
        else:
            return False


# --------PLANETARY PARAMETERS UI------------------------------------------
# Get the user's confirmation of values that will later be used in lightcurve fit
def get_planetary_parameters(candplanetbool, userpdict, pdict=None):
    log_info("*******************************************")
    log_info("Planetary Parameters for Lightcurve Fitting")

    # The order of planet_params list must match the pDict that is declared when scraping the NASA Exoplanet Archive
    planet_params = ["Target Star RA in the form: HH:MM:SS (ignore the decimal values)",
                     "Target Star DEC in form: <sign>DD:MM:SS (ignore the decimal values and don't forget the '+' or '-' sign!)",
                     "Planet's Name",
                     "Host Star's Name",
                     "Orbital Period (days)",
                     "Orbital Period Uncertainty (days) \n(Keep in mind that 1.2e-34 is the same as 1.2 x 10^-34)",
                     "Published Mid-Transit Time (BJD_UTC)",
                     "Mid-Transit Time Uncertainty (BJD-UTC)",
                     "Ratio of Planet to Stellar Radius (Rp/Rs)",
                     "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty",
                     "Ratio of Distance to Stellar Radius (a/Rs)",
                     "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty",
                     "Orbital Inclination (deg)",
                     "Orbital Inclination (deg) Uncertainty",
                     "Argument of Periastron (deg)",
                     "Orbital Eccentricity (0 if null)",
                     "Star Effective Temperature (K)",
                     "Star Effective Temperature Positive Uncertainty (K)",
                     "Star Effective Temperature Negative Uncertainty (K)",
                     "Star Metallicity ([FE/H])",
                     "Star Metallicity Positive Uncertainty ([FE/H])",
                     "Star Metallicity Negative Uncertainty ([FE/H])",
                     "Star Surface Gravity (log(g))",
                     "Star Surface Gravity Positive Uncertainty (log(g))",
                     "Star Surface Gravity Negative Uncertainty (log(g))",
                     "Star Distance (pc)",
                     "Star Proper Motion RA (mas/yr)",
                     "Star Proper Motion DEC (mas/yr)"]

    # Conversion between hours to degrees if user entered ra and dec
    if userpdict['ra'] is None:
        userpdict['ra'] = user_input(f"\nEnter the {planet_params[0]}: ", type_=str)
    if userpdict['dec'] is None:
        userpdict['dec'] = user_input(f"\nEnter the {planet_params[1]}: ", type_=str)
    if type(userpdict['ra']) and type(userpdict['dec']) is str:
        userpdict['ra'], userpdict['dec'] = radec_hours_to_degree(userpdict['ra'], userpdict['dec'])

    radeclist = ['ra', 'dec']
    if not candplanetbool:
        for idx, item in enumerate(radeclist):
            uncert = 20 / 3600
            if pdict[item] - uncert <= userpdict[item] <= pdict[item] + uncert:
                continue
            else:
                log_info(f"\n\nWarning: {pdict['pName']} initialization file's {planet_params[idx]} does not match "
                         "the value scraped by EXOTIC from the NASA Exoplanet Archive.\n", warn=True)
                log_info(f"\tNASA Exoplanet Archive value (degrees): {pdict[item]}", warn=True)
                log_info(f"\tInitialization file value (degrees): {userpdict[item]}", warn=True)
                log_info("\nWould you like to:"
                         "\n  (1) use NASA Exoplanet Archive value, "
                         "\n  (2) use initialization file value, or "
                         "\n  (3) enter in a new value.", warn=True)
                option = user_input("Which option do you choose? (1/2/3): ", type_=int, values=[1, 2, 3])

                if option == 1:
                    userpdict[item] = pdict[item]
                elif option == 2:
                    continue
                else:
                    userpdict['ra'] = user_input(f"Enter the {planet_params[0]}: ", type_=str)
                    userpdict['dec'] = user_input(f"Enter the {planet_params[1]}: ", type_=str)
                    break

    if type(userpdict['ra']) and type(userpdict['dec']) is str:
        userpdict['ra'], userpdict['dec'] = radec_hours_to_degree(userpdict['ra'], userpdict['dec'])

    # Exoplanet confirmed in NASA Exoplanet Archive
    if not candplanetbool:
        log_info(f"*** Here are the values scraped from the NASA Exoplanet Archive for {pdict['pName']} that were not "
                 "set (or set to null) in your initialization file. ***")

        for i, key in enumerate(userpdict):
            if key in ('ra', 'dec'):
                continue
            if key in ('pName', 'sName'):
                userpdict[key] = pdict[key]
            # Initialization planetary parameters match NEA
            if pdict[key] == userpdict[key]:
                continue
            # Initialization planetary parameters don't match NASA Exoplanet Archive
            if userpdict[key] is not None:
                log_info(f"\n\nWarning: {pdict['pName']} initialization file's {planet_params[i]} does not match "
                         "the value scraped by EXOTIC from the NASA Exoplanet Archive.\n", warn=True)
                log_info(f"\tNASA Exoplanet Archive value: {pdict[key]}", warn=True)
                log_info(f"\tInitialization file value: {userpdict[key]}", warn=True)
                log_info("\nWould you like to: "
                         "\n  (1) use NASA Exoplanet Archive value, "
                         "\n  (2) use initialization file value, or "
                         "\n  (3) enter in a new value.", warn=True)
                option = user_input("Which option do you choose? (1/2/3): ", type_=int, values=[1, 2, 3])
                if option == 1:
                    userpdict[key] = pdict[key]
                elif option == 2:
                    continue
                else:
                    userpdict[key] = user_input(f"Enter the {planet_params[i]}: ", type_=type(userpdict[key]))
            # Did not use initialization file or null
            else:
                log_info(f"\n {pdict['pName']} {planet_params[i]}: {pdict[key]}")
                agreement = user_input("Do you agree? (y/n): ", type_=str, values=['y', 'n'])
                if agreement == 'y':
                    userpdict[key] = pdict[key]
                else:
                    userpdict[key] = user_input(f"Enter the {planet_params[i]}: ", type_=type(pdict[key]))

    # Exoplanet not confirmed in NASA Exoplanet Archive
    else:
        for i, key in enumerate(userpdict):
            if key in ('ra', 'dec'):
                continue
            # Used initialization file and is not empty
            if userpdict[key] is not None:
                agreement = user_input(f"{planet_params[i]}: {userpdict[key]} \nDo you agree? (y/n): ",
                                       type_=str, values=['y', 'n'])
                if agreement == 'y':
                    continue
                else:
                    userpdict[key] = user_input(f"Enter the {planet_params[i]}: ", type_=type(userpdict[key]))
            # Did not use initialization file
            else:
                if key in ('pName', 'sName'):
                    userpdict[key] = user_input(f"\nEnter the {planet_params[i]}: ", type_=str)
                else:
                    userpdict[key] = user_input(f"Enter the {planet_params[i]}: ", type_=float)
    return userpdict


# Conversion of Right Ascension and Declination: hours -> degrees
def radec_hours_to_degree(ra, dec):
    while True:
        try:
            ra_value = str(ra).strip()
            dec_value = str(dec).strip()

            # Accept either sexagesimal RA strings (HH:MM:SS) or decimal RA degrees.
            # A decimal-like value with no separators should be treated as degrees.
            ra_unit = u.hourangle if any(sep in ra_value for sep in (':', ' ')) else u.deg

            # Declination can be provided as either sexagesimal or decimal degrees.
            dec_unit = u.deg
            if any(sep in dec_value for sep in (':', ' ')):
                dec_value = dec_value.replace(':', ' ')

            if ra_unit is u.hourangle:
                ra_value = ra_value.replace(':', ' ')

            c = SkyCoord(ra=ra_value, dec=dec_value, unit=(ra_unit, dec_unit))
            return c.ra.degree, c.dec.degree
        except ValueError:
            log_info("Error: The format entered for Right Ascension and/or Declination is not correct, "
                     "please try again.", error=True)
            ra = input("Input the Right Ascension of target (HH:MM:SS): ")
            dec = input("Input the Declination of target (<sign>DD:MM:SS): ")


def check_all_standard_filters(ld, observed_filter):
    if ld.check_standard(observed_filter):
        return True
    elif observed_filter['filter']:
        filter_name = observed_filter['filter'].lower().replace(' ', '')
        filter_name = re.sub(ld_re_punct_p, '', filter_name)
        filter_abbreviation = next((filter_abbr for filter_abbr in LimbDarkening.fwhm_names_nonspecific.keys()
                                    if filter_name == filter_abbr.lower()), None)
        filter_desc = next((filter_desc for filter_desc in LimbDarkening.fwhm_names_nonspecific.values()
                            if filter_name == re.sub(ld_re_punct_p, '', filter_desc.lower().replace(' ', ''))),
                           None)

        if filter_abbreviation:
            observed_filter['filter'] = LimbDarkening.fwhm_names_nonspecific.get(filter_abbreviation)
            observed_filter['name'] = filter_abbreviation
            custom_range(ld, observed_filter)
            return True

        if filter_desc:
            observed_filter['filter'] = filter_desc
            observed_filter['name'] = next((k for k, v in LimbDarkening.fwhm_names_nonspecific.items() if v == filter_desc))
            custom_range(ld, observed_filter)
            return True

    return False


def custom_range(ld, observed_filter):
    while True:
        if ld.check_fwhm(observed_filter):
            ld.set_filter(observed_filter['name'], observed_filter['filter'],
                          float(observed_filter['wl_min']), float(observed_filter['wl_max']))
            return
        else:
            observed_filter['wl_min'] = user_input(f"FWHM minimum wavelength (nm):", type_=str)
            observed_filter['wl_max'] = user_input(f"FWHM maximum wavelength (nm):", type_=str)


def standard_filter(ld, observed_filter):
    LimbDarkening.standard_list()

    while True:
        if not observed_filter['filter']:
            observed_filter['filter'] = user_input("\nPlease enter in the Filter Name or Abbreviation "
                                           "(EX: Johnson V, V, STB, RJ): ", type_=str)

        if check_all_standard_filters(ld, observed_filter):
            return
        else:
            log_info("\nError: The entered filter is not in the provided list of standard filters.", warn=True)
            observed_filter['filter'] = None


def user_entered_ld(ld, observed_filter):
    order = ['first', 'second', 'third', 'fourth']

    input_list = [(f"\nEnter in your {order[i]} nonlinear term:",
                   f"\nEnter in your {order[i]} nonlinear term uncertainty:") for i in range(len(order))]
    ld_ = [(user_input(input_[0], type_=float), user_input(input_[1], type_=float)) for input_ in input_list]

    custom_range(ld, observed_filter)
    ld.set_ld(ld_[0], ld_[1], ld_[2], ld_[3])


def nonlinear_ld(ld, info_dict):
    user_entered = False
    observed_filter = {
        'filter': info_dict['filter'],
        'name': None,
        'wl_min': info_dict['wl_min'],
        'wl_max': info_dict['wl_max']
    }
    ld.check_fwhm(observed_filter)

    if not check_all_standard_filters(ld, observed_filter):
        if observed_filter['wl_min'] and observed_filter['wl_max']:
            custom_range(ld, observed_filter)
            ld.set_filter('N/A', "Custom", float(observed_filter['wl_min']), float(observed_filter['wl_max']))
        else:
            opt = info_dict.get('ld_uncertainties')

            if isinstance(opt, str):
                opt = opt.lower().strip()

            if opt not in ('y', 'n'):
                opt = user_input("\nWould you like EXOTIC to calculate your limb darkening parameters "
                                 "with uncertainties? (y/n):", type_=str, values=['y', 'n'])

            if opt == 'y':
                opt = user_input("Please enter 1 to use a standard filter or 2 for a customized filter:",
                                 type_=int, values=[1, 2])
                if opt == 1:
                    observed_filter['filter'] = None
                    standard_filter(ld, observed_filter)
                elif opt == 2:
                    custom_range(ld, observed_filter)
                    ld.set_filter('N/A', "Custom", float(observed_filter['wl_min']), float(observed_filter['wl_max']))
            else:
                user_entered_ld(ld, observed_filter)
                user_entered = True

    if not user_entered:
        ld.calculate_ld()

    info_dict['filter'] = ld.filter_name
    info_dict['filter_desc'] = ld.filter_desc
    info_dict['wl_min'] = ld.wl_min
    info_dict['wl_max'] = ld.wl_max


def get_ld_values(planet_dict, info_dict):
    ld_obj = LimbDarkening(planet_dict)
    nonlinear_ld(ld_obj, info_dict)

    ld0 = ld_obj.ld0
    ld1 = ld_obj.ld1
    ld2 = ld_obj.ld2
    ld3 = ld_obj.ld3
    ld = [ld0[0], ld1[0], ld2[0], ld3[0]]

    return ld, ld0, ld1, ld2, ld3


def corruption_check(files):
    valid_files = []
    for file in files:
        plateStatus.setCurrentFilename(file)
        try:
            with fits.open(name=file, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True) as hdu1:
                valid_files.append(file)
        except OSError as e:
            # Since google collab can have problems with initial load of big data sets from google
            # drive, lets pause and retry this once when we fail: if the file was corrupted the first time,
            # nothing will get better...
            log_info(f"Warning: retrying verify\n\t-File: {file}\n\t-Reason: {e}", warn=True)
            sleep(5)
            try: 
                with fits.open(name=file, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True) as hdu1:
                    valid_files.append(file)
            except OSError as e:
                log.debug(f"Warning: corrupted file found and removed from reduction\n\t-File: {file}\n\t-Reason: {e}")
                plateStatus.fitsFormatError(e)
    return valid_files

def check_wcs(fits_file, save_directory, plate_opt, rt=False, use_nextastro_astrometry=False,
              ra=None, dec=None, pixel_scale=None, ignore_header_wcs=False):
    wcs_file = None

    if not ignore_header_wcs and search_wcs(fits_file).is_celestial:
        if plate_opt == 'y' and not rt:
            log_info("Your FITS files already have WCS (World Coordinate System) information in their headers. "
                     "EXOTIC will use the existing header WCS and skip external plate solving.")
        else:
            log_info("Your FITS files have WCS (World Coordinate System) information in their headers. "
                     "EXOTIC will proceed to use these. "
                     "NOTE: If you do not trust your WCS coordinates, "
                     "please restart EXOTIC after enabling plate solutions via astrometry.net.")
        return fits_file

    if plate_opt == 'y' and not rt:
        wcs_file = get_wcs(fits_file, save_directory, use_nextastro_astrometry=use_nextastro_astrometry, ra=ra, dec=dec, pixel_scale=pixel_scale)
    if ignore_header_wcs:
        if wcs_file:
            log_info("Ignoring FITS header WCS for alignment and using the legacy image-to-image alignment path.")
        else:
            log_info("Ignoring FITS header WCS and using the legacy image-to-image alignment path.")
        return wcs_file

    return wcs_file


def search_wcs(file):
    header = get_first_image_header(file)
    return search_wcs_from_header(header)


def search_wcs_from_header(header):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=FITSFixedWarning)
        return WCS(header)
        # return WCS(fits.open(file)[('SCI', 1)].header)


def get_first_image_header(file_name):
    extension = 0
    header = fits.getheader(filename=file_name, ext=extension)
    while header.get('NAXIS', 0) == 0:
        extension += 1
        header = fits.getheader(filename=file_name, ext=extension)
    return header


def collect_celestial_wcs_coverage(inputfiles):
    has_celestial_wcs = []
    missing_wcs_files = []
    for file_name in inputfiles:
        file_has_celestial_wcs = False
        try:
            image_header = get_first_image_header(file_name)
            file_has_celestial_wcs = search_wcs_from_header(image_header).is_celestial
        except Exception:
            file_has_celestial_wcs = False

        has_celestial_wcs.append(file_has_celestial_wcs)
        if not file_has_celestial_wcs:
            missing_wcs_files.append(str(file_name))

    return np.array(has_celestial_wcs, dtype=bool), missing_wcs_files


def evaluate_celestial_wcs_coverage(inputfiles):
    has_celestial_wcs, missing_wcs_files = collect_celestial_wcs_coverage(inputfiles)
    total_files = len(inputfiles)
    all_have_celestial_wcs = total_files > 0 and bool(has_celestial_wcs.all())
    return all_have_celestial_wcs, missing_wcs_files


def log_file_preview(file_names, label):
    if not file_names:
        return

    preview = ", ".join([_display_filename(file_name) for file_name in file_names[:3]])
    remainder = len(file_names) - 3
    if remainder > 0:
        preview = f"{preview}, ... (+{remainder} more)"
    log.debug(f"{label}: {preview}")


def log_missing_celestial_wcs_preview(missing_wcs_files):
    log_file_preview(missing_wcs_files, "Files without usable celestial WCS")


def format_file_preview_for_user(file_names, limit=6):
    if not file_names:
        return ""

    display_names = [_display_filename(file_name) for file_name in file_names[:limit]]
    remainder = len(file_names) - len(display_names)
    if remainder > 0:
        display_names.append(f"... (+{remainder} more)")
    return ", ".join(display_names)


def leading_rejected_reference_prefix(ordered_inputfiles, dropped_files):
    if ordered_inputfiles is None:
        return [], None

    ordered_inputfiles = [str(file_name) for file_name in ordered_inputfiles]
    dropped_lookup = {str(file_name) for file_name in (dropped_files or [])}
    leading_rejected = []
    next_candidate = None

    for file_name in ordered_inputfiles:
        if file_name in dropped_lookup:
            leading_rejected.append(file_name)
            continue
        next_candidate = file_name
        break

    return leading_rejected, next_candidate


def abort_if_reference_frame_rejected(reference_file, dropped_files, ordered_inputfiles=None,
                                      rejection_label="Pointing precheck"):
    if reference_file is None or not dropped_files:
        return False

    reference_file = str(reference_file)
    dropped_files = [str(file_name) for file_name in dropped_files]
    if reference_file not in dropped_files:
        return False

    other_dropped_files = [file_name for file_name in dropped_files if file_name != reference_file]
    leading_rejected_files, next_reference_candidate = leading_rejected_reference_prefix(
        ordered_inputfiles,
        dropped_files,
    )
    if not leading_rejected_files:
        leading_rejected_files = [reference_file]

    log_info(
        f"Error: {rejection_label} rejected the first usable image "
        f"({_display_filename(reference_file)}). EXOTIC uses that frame as the reference image for "
        "the supplied target and comparison-star pixel coordinates, so it is not safe to continue "
        "with a different reference image.",
        error=True,
    )
    if other_dropped_files:
        log_info(
            f"{rejection_label} also rejected {len(other_dropped_files)} other frame(s): "
            f"{format_file_preview_for_user(other_dropped_files)}",
            error=True,
        )

    leading_preview = format_file_preview_for_user(leading_rejected_files)
    if len(leading_rejected_files) == 1:
        removal_instruction = (
            f"Please remove or move this rejected frame and run again: {leading_preview}."
        )
    else:
        removal_instruction = (
            f"Please remove or move these leading rejected frames and run again: {leading_preview}."
        )

    if next_reference_candidate is not None:
        removal_instruction += (
            f" The next remaining frame would be "
            f"{_display_filename(next_reference_candidate)}."
        )
    else:
        removal_instruction += (
            " No non-rejected frame remains after that prefix, so this dataset still would not "
            "have a usable reference image."
        )

    log_info(
        removal_instruction,
        error=True,
    )
    log_info(
        "If you need to keep those frames, reorder the dataset so a good reference image comes first, "
        "or set optional_info 'pointing_rejection_sigma' to 0 to disable this precheck.",
        error=True,
    )
    return True


def collect_wcs_frame_center_pointings(inputfiles):
    positions = np.full((len(inputfiles), 2), np.nan, dtype=float)
    usable_mask = np.zeros(len(inputfiles), dtype=bool)
    usable_indices = []
    ra_values = []
    dec_values = []

    for index, file_name in enumerate(inputfiles):
        try:
            image_header = get_first_image_header(file_name)
            wcs = search_wcs_from_header(image_header)
            if not wcs.is_celestial:
                continue

            width = int(image_header.get("NAXIS1", 0))
            height = int(image_header.get("NAXIS2", 0))
            if width <= 0 or height <= 0:
                continue

            center_x = (width - 1) / 2.0
            center_y = (height - 1) / 2.0
            ra_deg, dec_deg = wcs.pixel_to_world_values(center_x, center_y)
            if not np.isfinite(ra_deg) or not np.isfinite(dec_deg):
                continue

            usable_indices.append(index)
            ra_values.append(float(ra_deg))
            dec_values.append(float(dec_deg))
        except Exception:
            continue

    if not usable_indices:
        return positions, usable_mask

    coords = SkyCoord(ra=np.asarray(ra_values) * u.deg, dec=np.asarray(dec_values) * u.deg, frame='icrs')
    reference_coord = coords[0]
    delta_lon, delta_lat = reference_coord.spherical_offsets_to(coords)
    offsets = np.column_stack((delta_lon.to_value(u.arcsec), delta_lat.to_value(u.arcsec)))

    for index, offset in zip(usable_indices, offsets):
        positions[index] = offset
        usable_mask[index] = True

    return positions, usable_mask


def log_pointing_precheck_alignment_progress(i, total_files, file_name):
    log_info(
        f"Pointing precheck alignment progress: file {i + 1} of {total_files} : "
        f"{_display_filename(file_name)}"
    )


def _pointing_precheck_return(positions, usable_mask, alignment_transforms, return_transforms):
    if return_transforms:
        return positions, usable_mask, alignment_transforms
    return positions, usable_mask


def _filter_alignment_transform_cache(alignment_transforms, retained_files):
    if not alignment_transforms:
        return {}

    retained_keys = {str(file_name) for file_name in retained_files}
    return {
        file_key: tform
        for file_key, tform in alignment_transforms.items()
        if file_key in retained_keys
    }


def collect_transform_frame_pointings(inputfiles, frame_loader=None, return_transforms=False,
                                      multiprocess_transformations=None,
                                      generalDark=None, generalBias=None, generalFlat=None,
                                      demosaic_fmt=None, demosaic_out=None, demosaic_mult=None):
    positions = np.full((len(inputfiles), 2), np.nan, dtype=float)
    usable_mask = np.zeros(len(inputfiles), dtype=bool)
    alignment_transforms = {}

    if len(inputfiles) == 0:
        return _pointing_precheck_return(positions, usable_mask, alignment_transforms, return_transforms)

    if frame_loader is None:
        frame_loader = lambda file_name: load_calibrated_reduction_image(
            file_name,
            generalDark,
            generalBias,
            generalFlat,
            demosaic_fmt,
            demosaic_out,
            demosaic_mult,
        )

    total_files = len(inputfiles)
    log_pointing_precheck_alignment_progress(0, total_files, inputfiles[0])
    try:
        reference_image = frame_loader(inputfiles[0])
    except Exception as exc:
        log_info(
            f"Warning: pointing precheck alignment fallback could not load the reference frame "
            f"{_display_filename(inputfiles[0])} ({exc}).",
            warn=True,
        )
        return _pointing_precheck_return(positions, usable_mask, alignment_transforms, return_transforms)

    if getattr(reference_image, "ndim", 0) != 2:
        log_info("Warning: pointing precheck alignment fallback requires 2-D images; skipping.", warn=True)
        return _pointing_precheck_return(positions, usable_mask, alignment_transforms, return_transforms)

    height, width = reference_image.shape
    reference_anchor = np.array([[(width - 1) / 2.0, (height - 1) / 2.0]], dtype=float)
    positions[0] = reference_anchor[0]
    usable_mask[0] = True
    alignment_transforms[str(inputfiles[0])] = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])

    if multiprocess_transformations is not None and multiprocess_transformations > 0 and len(inputfiles) > 1:
        try:
            positions, usable_mask, alignment_transforms = build_multiprocess_pointing_precheck_transforms(
                inputfiles,
                multiprocess_transformations,
                reference_anchor,
                generalDark=generalDark,
                generalBias=generalBias,
                generalFlat=generalFlat,
                demosaic_fmt=demosaic_fmt,
                demosaic_out=demosaic_out,
                demosaic_mult=demosaic_mult,
            )
            return _pointing_precheck_return(positions, usable_mask, alignment_transforms, return_transforms)
        except Exception as exc:
            log_info(
                "Warning: pointing precheck multiprocessing failed; falling back to serial alignment "
                f"({exc}).",
                warn=True,
            )

    for index, file_name in enumerate(inputfiles[1:], start=1):
        log_pointing_precheck_alignment_progress(index, total_files, file_name)
        try:
            image_data = frame_loader(file_name)
            if getattr(image_data, "ndim", 0) != 2:
                continue

            tform = transformation(
                image_data,
                file_name,
                report_failure=False,
                reference_image=reference_image,
            )
            mapped_anchor = np.asarray(tform(reference_anchor), dtype=float).reshape(-1, 2)[0]
            if np.all(np.isfinite(mapped_anchor)):
                positions[index] = mapped_anchor
                usable_mask[index] = True
                alignment_transforms[str(file_name)] = tform
        except Exception:
            continue

    return _pointing_precheck_return(positions, usable_mask, alignment_transforms, return_transforms)


def sigma_clip_pointing_positions(positions, sigma=3.0, max_iters=5):
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError("positions must be an Nx2 array")

    finite_mask = np.all(np.isfinite(positions), axis=1)
    keep_mask = finite_mask.copy()
    if np.count_nonzero(keep_mask) < POINTING_REJECTION_MIN_FRAMES:
        return keep_mask

    sigma = float(sigma)
    for _ in range(max_iters):
        candidate_positions = positions[keep_mask]
        if candidate_positions.shape[0] < POINTING_REJECTION_MIN_FRAMES:
            break

        center = np.nanmedian(candidate_positions, axis=0)
        deltas = candidate_positions - center
        radial_offsets = np.hypot(deltas[:, 0], deltas[:, 1])

        scatter_x = robust_scatter(deltas[:, 0])
        scatter_y = robust_scatter(deltas[:, 1])
        radial_scatter = robust_scatter(radial_offsets)

        if not np.isfinite(scatter_x) or scatter_x <= 0:
            scatter_x = radial_scatter
        if not np.isfinite(scatter_y) or scatter_y <= 0:
            scatter_y = radial_scatter

        if (not np.isfinite(scatter_x) or scatter_x <= 0
                or not np.isfinite(scatter_y) or scatter_y <= 0):
            break

        normalized_distance = np.sqrt((deltas[:, 0] / scatter_x) ** 2 + (deltas[:, 1] / scatter_y) ** 2)
        current_keep = normalized_distance <= sigma
        if np.all(current_keep):
            break

        updated_keep = keep_mask.copy()
        updated_keep[np.flatnonzero(keep_mask)] = current_keep
        if np.array_equal(updated_keep, keep_mask):
            break
        keep_mask = updated_keep

    return keep_mask


def filter_pointing_outlier_frames(inputfiles, pointing_rejection_sigma=None, ignore_header_wcs=False,
                                   frame_loader=None, return_alignment_transforms=False,
                                   multiprocess_transformations=None,
                                   generalDark=None, generalBias=None, generalFlat=None,
                                   demosaic_fmt=None, demosaic_out=None, demosaic_mult=None):
    inputfiles = np.array(inputfiles)
    keep_mask = np.ones(len(inputfiles), dtype=bool)
    alignment_transforms = {}

    def format_result(result_inputfiles, result_keep_mask, dropped_files):
        if return_alignment_transforms:
            return (
                result_inputfiles,
                result_keep_mask,
                dropped_files,
                _filter_alignment_transform_cache(alignment_transforms, result_inputfiles),
            )
        return result_inputfiles, result_keep_mask, dropped_files

    if len(inputfiles) == 0 or pointing_rejection_sigma is None:
        return format_result(inputfiles, keep_mask, [])

    if len(inputfiles) < POINTING_REJECTION_MIN_FRAMES:
        log_info(
            f"Pointing precheck skipped: only {len(inputfiles)} frame(s); "
            f"need at least {POINTING_REJECTION_MIN_FRAMES}.",
        )
        return format_result(inputfiles, keep_mask, [])

    positions = None
    usable_mask = None
    mode_label = None

    if not ignore_header_wcs:
        wcs_positions, wcs_usable_mask = collect_wcs_frame_center_pointings(inputfiles)
        usable_wcs_count = int(np.count_nonzero(wcs_usable_mask))
        if usable_wcs_count == len(inputfiles):
            positions = wcs_positions
            usable_mask = wcs_usable_mask
            mode_label = "WCS"
        elif usable_wcs_count > 0:
            log_info(
                f"Pointing precheck: usable WCS-derived pointing centers found for "
                f"{usable_wcs_count}/{len(inputfiles)} frame(s); falling back to alignment-derived positions."
            )
        else:
            log_info("Pointing precheck: no usable WCS-derived pointing centers found; using alignment-derived positions.")

    if positions is None:
        positions, usable_mask, alignment_transforms = collect_transform_frame_pointings(
            inputfiles,
            frame_loader=frame_loader,
            return_transforms=True,
            multiprocess_transformations=multiprocess_transformations,
            generalDark=generalDark,
            generalBias=generalBias,
            generalFlat=generalFlat,
            demosaic_fmt=demosaic_fmt,
            demosaic_out=demosaic_out,
            demosaic_mult=demosaic_mult,
        )
        mode_label = "alignment"

    usable_count = int(np.count_nonzero(usable_mask))
    if usable_count < POINTING_REJECTION_MIN_FRAMES:
        log_info(
            f"Pointing precheck skipped: only {usable_count} usable {mode_label}-derived pointing estimate(s); "
            f"need at least {POINTING_REJECTION_MIN_FRAMES}.",
        )
        return format_result(inputfiles, keep_mask, [])

    keep_mask[np.flatnonzero(usable_mask)] = sigma_clip_pointing_positions(
        positions[usable_mask],
        sigma=pointing_rejection_sigma,
        max_iters=POINTING_REJECTION_MAX_ITERS,
    )

    dropped_files = inputfiles[~keep_mask].tolist()
    if not dropped_files:
        log_info(
            f"Pointing precheck ({mode_label}): no frames exceeded the "
            f"{float(pointing_rejection_sigma):g}-sigma pointing threshold."
        )
        return format_result(inputfiles, keep_mask, [])

    retained_files = inputfiles[keep_mask]
    log_info(
        f"Pointing precheck ({mode_label}): {len(retained_files)}/{len(inputfiles)} frame(s) remain after "
        f"dropping {len(dropped_files)} file(s) beyond {float(pointing_rejection_sigma):g} sigma from the "
        "median pointing."
    )
    log_file_preview(dropped_files, "Pointing precheck dropped files")
    return format_result(retained_files, keep_mask, dropped_files)


def filter_sparse_missing_wcs_frames(inputfiles, ignore_header_wcs=False, max_missing_fraction=None):
    inputfiles = np.array(inputfiles)
    keep_mask = np.ones(len(inputfiles), dtype=bool)
    if ignore_header_wcs or len(inputfiles) == 0:
        return inputfiles, keep_mask, []

    if max_missing_fraction is None:
        max_missing_fraction = SPARSE_MISSING_WCS_DROP_THRESHOLD

    keep_mask, missing_wcs_files = collect_celestial_wcs_coverage(inputfiles)
    missing_count = len(missing_wcs_files)
    total_files = len(inputfiles)
    if missing_count == 0:
        return inputfiles, keep_mask, []

    missing_fraction = missing_count / total_files
    if missing_count < total_files and missing_fraction < max_missing_fraction:
        retained_files = inputfiles[keep_mask]
        threshold_percent = max_missing_fraction * 100.0
        log_info(
            f"WCS precheck: {len(retained_files)}/{total_files} files have celestial WCS. "
            f"Dropping {missing_count} file(s) without celestial WCS because they are below the "
            f"{threshold_percent:g}% threshold."
        )
        log_missing_celestial_wcs_preview(missing_wcs_files)
        return retained_files, keep_mask, missing_wcs_files

    return inputfiles, np.ones(total_files, dtype=bool), []


def should_use_multiprocess_transform_precompute(inputfiles, requested_processes, ignore_header_wcs=False):
    if requested_processes is None or requested_processes <= 0:
        return False

    if ignore_header_wcs:
        log_info("Header WCS ignore override enabled. Keeping multiprocessing transformation precompute.")
        return True

    all_have_celestial_wcs, missing_wcs_files = evaluate_celestial_wcs_coverage(inputfiles)
    if all_have_celestial_wcs:
        log_info("All input FITS files have celestial WCS in their headers. "
                 "Skipping multiprocessing transformation precompute.")
        return False

    total_files = len(inputfiles)
    missing_count = len(missing_wcs_files)
    log_info(f"WCS precheck: {total_files - missing_count}/{total_files} files have celestial WCS. "
             "Keeping multiprocessing transformation precompute for fallback alignment.")
    log_missing_celestial_wcs_preview(missing_wcs_files)

    return True


def get_wcs(file, directory="", use_nextastro_astrometry=False, ra=None, dec=None, pixel_scale=None):
    astrometry_service = 'NextAstro astrometry server (https://astrometry.nextastro.org/)' if use_nextastro_astrometry else 'nova.astrometry.net'
    log_info("\nGetting the plate solution for your imaging file to translate pixel coordinates on the sky. "
             f"\nUsing astrometry service: {astrometry_service}."
             "\nPlease wait....")

    if use_nextastro_astrometry:
        print("Contacting NextAstro Astrometry Server")
        nextastro_solver = NextAstroPlateSolution(
            file=file,
            directory=directory,
            ra=ra,
            dec=dec,
            pixel_scale=pixel_scale,
            suppress_fail_warning=True,
            message_logger=log_info
        )
        wcs_file = nextastro_solver.plate_solution()
        if wcs_file:
            return wcs_file

        nextastro_bad_gateway = nextastro_solver.last_http_status == 502
        if nextastro_bad_gateway:
            log_info("NextAstro Server not responding. Will try nova.astrometry.net")
        else:
            log_info("NextAstro astrometry server did not return a solution; falling back to nova.astrometry.net.")
        print("Communication with nova.astrometry.net")
        nova_solver = PlateSolution(file=file, directory=directory, ra=ra, dec=dec,
                                    pixel_scale=pixel_scale, suppress_fail_warning=True)
        wcs_file = nova_solver.plate_solution()
        if wcs_file:
            return wcs_file
        if nextastro_bad_gateway:
            log_info("NextAstro Server not responding. Both astrometry methods trialed, pushing forward without astrometry solution")
            return False
        return PlateSolution.fail(nova_solver.last_error_type or 'plate solution lookup')

    animate_toggle(True)
    nova_solver = PlateSolution(file=file, directory=directory, ra=ra, dec=dec,
                                pixel_scale=pixel_scale, suppress_fail_warning=True)
    wcs_file = nova_solver.plate_solution()
    if wcs_file:
        animate_toggle()
        return wcs_file

    log_info("nova.astrometry.net did not return a solution; falling back to NextAstro astrometry server.")
    print("Contacting NextAstro Astrometry Server")
    nextastro_solver = NextAstroPlateSolution(
        file=file,
        directory=directory,
        ra=ra,
        dec=dec,
        pixel_scale=pixel_scale,
        suppress_fail_warning=True,
        message_logger=log_info
    )
    wcs_file = nextastro_solver.plate_solution()
    animate_toggle()
    if wcs_file:
        return wcs_file
    if nextastro_solver.last_http_status == 502:
        log_info("NextAstro Server not responding. Both astrometry methods trialed, pushing forward without astrometry solution")
        return False
    return PlateSolution.fail(nextastro_solver.last_error_type or 'plate solution lookup',
                              service_name=f'NextAstro ({nextastro_solver.api_url})')


# Getting the right ascension and declination for every pixel in imaging file if there is a plate solution
def _resolve_wcs_image_dimensions(header, image_shape=None):
    width = header.get('NAXIS1', header.get('ZNAXIS1'))
    height = header.get('NAXIS2', header.get('ZNAXIS2'))
    if width is not None and height is not None:
        return int(width), int(height)

    if image_shape is not None and len(image_shape) >= 2:
        return int(image_shape[-1]), int(image_shape[-2])

    wcs_header = WCS(header)
    if wcs_header.pixel_shape is not None and len(wcs_header.pixel_shape) >= 2:
        return int(wcs_header.pixel_shape[0]), int(wcs_header.pixel_shape[1])

    if wcs_header.array_shape is not None and len(wcs_header.array_shape) >= 2:
        return int(wcs_header.array_shape[1]), int(wcs_header.array_shape[0])

    raise KeyError("Keyword 'NAXIS1' not found.")


def get_ra_dec(header, image_shape=None):
    wcs_header = WCS(header)
    width, height = _resolve_wcs_image_dimensions(header, image_shape=image_shape)
    xaxis = np.arange(width)
    yaxis = np.arange(height)
    x, y = np.meshgrid(xaxis, yaxis)
    return wcs_header.all_pix2world(x, y, 1)


def deg_to_pix(exp_ra, exp_dec, ra_list, dec_list):
    dist = (ra_list - exp_ra) ** 2 + (dec_list - exp_dec) ** 2
    return np.unravel_index(dist.argmin(), dist.shape)


def project_target_pixel_wcs(exp_ra, exp_dec, ra_list, dec_list, wcs_header=None):
    if wcs_header is not None:
        try:
            x_pixel, y_pixel = WCS(wcs_header).all_world2pix(exp_ra, exp_dec, 0)
            x_pixel = float(np.asarray(x_pixel).reshape(-1)[0])
            y_pixel = float(np.asarray(y_pixel).reshape(-1)[0])
            if np.isfinite(x_pixel) and np.isfinite(y_pixel):
                return x_pixel, y_pixel
        except Exception as exc:
            log.debug(f"Direct WCS pixel projection failed; falling back to grid search: {exc}")

    calculated_y_pixel, calculated_x_pixel = deg_to_pix(exp_ra, exp_dec, ra_list, dec_list)
    return float(calculated_x_pixel), float(calculated_y_pixel)


def pixel_within_image(x_pixel, y_pixel, image_shape, margin=0.0):
    height, width = image_shape[:2]
    return (
        np.isfinite(x_pixel)
        and np.isfinite(y_pixel)
        and margin <= x_pixel < (width - margin)
        and margin <= y_pixel < (height - margin)
    )


def any_projected_coord_out_of_frame(coords, image_shape):
    for x_pixel, y_pixel in np.asarray(coords, dtype=float):
        if not pixel_within_image(x_pixel, y_pixel, image_shape):
            return True
    return False


def check_target_pixel_wcs(input_x_pixel, input_y_pixel, info_dict, ra_list, dec_list, image_data, obs_time,
                           non_interactive_run=False, wcs_header=None):
    """
    Verify the provided pixel coordinates match the target's right ascension and declination.
    """
    updated_ra, updated_dec = update_coordinates_with_proper_motion(info_dict, obs_time)

    calculated_x_pixel, calculated_y_pixel = project_target_pixel_wcs(
        updated_ra, updated_dec, ra_list, dec_list, wcs_header=wcs_header
    )

    if not pixel_within_image(calculated_x_pixel, calculated_y_pixel, image_data.shape):
        log_info("Warning: WCS-derived target pixel coordinates fall outside the image; "
                 "keeping the input target coordinates.", warn=True)
        return input_x_pixel, input_y_pixel

    centroid_margin = 7.5
    if not pixel_within_image(calculated_x_pixel, calculated_y_pixel, image_data.shape, margin=centroid_margin):
        log_info("Warning: WCS-derived target pixel coordinates are too close to the image edge for "
                 "centroid fitting; keeping the input target coordinates.", warn=True)
        return input_x_pixel, input_y_pixel

    centroid_x, centroid_y, sigma_x, sigma_y = get_psf_parameters(image_data, calculated_x_pixel, calculated_y_pixel)

    return check_coordinates(input_x_pixel, input_y_pixel, centroid_x, centroid_y, sigma_x, sigma_y,
                             calculated_x_pixel, calculated_y_pixel, non_interactive_run=non_interactive_run)


def get_psf_parameters(image_data, x_pixel, y_pixel):
    try:
        psf_data = fit_centroid(image_data, [x_pixel, y_pixel], 0)
    except Exception as exc:
        log.debug(f"Centroid fit failed while validating WCS target coordinates: {exc}")
        return np.nan, np.nan, np.nan, np.nan
    return psf_data[0], psf_data[1], psf_data[3], psf_data[4]


def check_coordinates(input_x_pixel, input_y_pixel, centroid_x, centroid_y, sigma_x, sigma_y,
                      calculated_x_pixel, calculated_y_pixel, non_interactive_run=False):
    while True:
        try:
            validate_pixel_coordinates(input_x_pixel, input_y_pixel, centroid_x, centroid_y, sigma_x, sigma_y)
            return input_x_pixel, input_y_pixel
        except ValueError:
            if non_interactive_run:
                if np.isfinite(centroid_x) and np.isfinite(centroid_y):
                    log_info("Proceeding with WCS-derived centroided target coordinates due to "
                             "--non-interactive-run.", warn=True)
                    return centroid_x, centroid_y
                log_info("Proceeding with WCS-derived target pixel coordinates due to "
                         "--non-interactive-run (centroid unavailable).", warn=True)
                return calculated_x_pixel, calculated_y_pixel
            new_x_pixel, new_y_pixel = prompt_user_for_coordinates(input_x_pixel, input_y_pixel,
                                                                   calculated_x_pixel, calculated_y_pixel)
            if new_x_pixel == input_x_pixel and new_y_pixel == input_y_pixel:
                return input_x_pixel, input_y_pixel
            else:
                input_x_pixel, input_y_pixel = new_x_pixel, new_y_pixel


def validate_pixel_coordinates(input_x_pixel, input_y_pixel, centroid_x, centroid_y, sigma_x, sigma_y):
    """
    Validating the provided pixel coordinates are within 5 PSF of the expected coordinates.
    """
    x_min = centroid_x - (sigma_x * 5)
    x_max = centroid_x + (sigma_x * 5)
    y_min = centroid_y - (sigma_y * 5)
    y_max = centroid_y + (sigma_y * 5)

    if not (x_min <= input_x_pixel <= x_max):
        log_info("\nWarning: The X Pixel Coordinate entered does not match the target's Right Ascension.", warn=True)
        raise ValueError
    if not (y_min <= input_y_pixel <= y_max):
        log_info("\nWarning: The Y Pixel Coordinate entered does not match the target's Declination.", warn=True)
        raise ValueError


def prompt_user_for_coordinates(input_x_pixel, input_y_pixel, calculated_x_pixel, calculated_y_pixel):
    log_info(f"Your input pixel coordinates: [{input_x_pixel}, {input_y_pixel}]")
    log_info(f"EXOTIC's calculated pixel coordinates: [{calculated_x_pixel}, {calculated_y_pixel}]")
    opt = user_input("Would you like to re-enter the pixel coordinates? (y/n): ", type_=str, values=['y', 'n'])

    if opt == 'y':
        use_suggested = user_input(
            f"Here are the suggested pixel coordinates:"
            f"  X Pixel: {calculated_x_pixel}"
            f"  Y Pixel: {calculated_y_pixel}"
            "\nWould you like to use these? (y/n): ",
            type_=str, values=['y', 'n']
        )

        if use_suggested == 'y':
            return calculated_x_pixel, calculated_y_pixel
        else:
            input_x_pixel = user_input("Please re-enter the target star's X Pixel Coordinate: ", type_=int)
            input_y_pixel = user_input("Please re-enter the target star's Y Pixel Coordinate: ", type_=int)

    return input_x_pixel, input_y_pixel


# Checks if comparison star is variable via querying SIMBAD
def query_variable_star_apis(ra, dec):
    # Convert comparison star coordinates from pixel to WCS
    sample = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    return vsx_variable(sample.ra.deg, sample.dec.deg)
    #
    # # Query SIMBAD and search identifier result table to determine if comparison star is variable in any form
    # # This is a secondary check if GAIA query returns inconclusive results
    # star_name = simbad_query(sample)
    # if not star_name:
    #     log_info("Warning: Your comparison star cannot be resolved in the SIMBAD star database; "
    #              "EXOTIC cannot check if it is variable or not. "
    #              "\nEXOTIC will still include this star in the reduction. "
    #              "\nPlease proceed with caution as we cannot check for stellar variability.\n", warn=True)
    #     return False
    # else:
    #     identifiers = Simbad.query_objectids(star_name)
    #
    #     for currName in identifiers:
    #         if "V*" in currName[0]:
    #             return True
    #     return False


@retry(stop=stop_after_delay(30))
def vsx_auid(ra, dec, radius=0.01, maglimit=14):
    try:
        url = f"https://www.aavso.org/vsx/index.php?view=api.list&ra={ra}&dec={dec}&radius={radius}&tomag={maglimit}&format=json"
        result = requests.get(url, timeout=30)
        result.raise_for_status()
        vsx_objects = extract_vsx_objects(result.json())
        if not vsx_objects:
            return False
        return vsx_objects[0].get('AUID', False) or False
    except Exception:
        log.info("\nThe target star does not have an AUID.")
        return False


def extract_vsx_objects(payload):
    if not isinstance(payload, dict):
        return []

    vsx_objects = payload.get('VSXObjects', [])
    if isinstance(vsx_objects, dict):
        vsx_object = vsx_objects.get('VSXObject', [])
        if isinstance(vsx_object, dict):
            return [vsx_object]
        if isinstance(vsx_object, list):
            return vsx_object
        return []

    if isinstance(vsx_objects, list):
        return vsx_objects

    return []


def describe_retry_exception(err):
    if isinstance(err, RetryError):
        last_attempt = getattr(err, 'last_attempt', None)
        attempt_number = getattr(last_attempt, 'attempt_number', None)
        try:
            root_cause = last_attempt.exception() if last_attempt is not None else None
        except Exception:
            root_cause = None

        if root_cause is not None:
            attempt_text = f" after {attempt_number} attempts" if attempt_number else ""
            return (f"{err.__class__.__name__}{attempt_text} "
                    f"({root_cause.__class__.__name__}: {root_cause})")

    return f"{err.__class__.__name__}: {err}"


def extract_http_status_code_from_error(err):
    if err is None:
        return None

    response = getattr(err, 'response', None)
    status_code = getattr(response, 'status_code', None)
    if status_code is not None:
        try:
            return int(status_code)
        except (TypeError, ValueError):
            return None

    status_match = re.search(r"\bHTTP\s+(\d{3})\b", str(err), flags=re.IGNORECASE)
    if status_match is None:
        return None

    try:
        return int(status_match.group(1))
    except (TypeError, ValueError):
        return None


def should_retry_nextastro_variability_error(err):
    if isinstance(err, requests.exceptions.RequestException):
        status_code = extract_http_status_code_from_error(err)
        return status_code is None or status_code in NEXTASTRO_VARIABILITY_RETRYABLE_HTTP_STATUS_CODES

    if isinstance(err, RuntimeError):
        status_code = extract_http_status_code_from_error(err)
        return status_code in NEXTASTRO_VARIABILITY_RETRYABLE_HTTP_STATUS_CODES

    return False


def submit_nextastro_variability_request(api_url, payload, content_encoding=None):
    request_body, headers, content_encoding, raw_size, compressed_size = build_compressed_json_request(
        payload,
        content_encoding=content_encoding,
    )
    log_info(
        "NextAstro variability request compression: "
        f"{content_encoding} ({compressed_size} bytes sent; {raw_size} bytes raw)"
    )
    return requests.post(api_url, data=request_body, headers=headers, timeout=30), content_encoding


@retry(stop=stop_after_delay(30))
def vsx_variable(ra, dec, radius=0.01, maglimit=14):
    default_vsx_error = None
    try:
        url = f"https://www.aavso.org/vsx/index.php?view=api.list&ra={ra}&dec={dec}&radius={radius}&tomag={maglimit}&format=json"
        result = requests.get(url, timeout=30)
        result.raise_for_status()
        vsx_objects = extract_vsx_objects(result.json())
        if not vsx_objects:
            return False

        first_vsx_object = vsx_objects[0]
        var = first_vsx_object.get('Category', '')

        if isinstance(var, str) and var.lower() == "variable":
            vname = first_vsx_object.get('Name')
            vdec = first_vsx_object.get('Declination2000')
            vra = first_vsx_object.get('RA2000')
            log_info(f"\nVSX variable check found {vname} at RA {vra}, DEC {vdec}\n"
                     f"and will be removed from reduction.", warn=True)
            return True
        return False
    except Exception as err:
        default_vsx_error = err

    try:
        log_info(f"\nDefault VSX request failed ({default_vsx_error}); falling back to NextAstro VSX server.", warn=True)
        fallback_result = nextastro_variability_test([(ra, dec)])
        is_variable = bool(fallback_result[0])
        if is_variable:
            log_info("\nNextAstro VSX fallback flagged this star as variable and it will be removed from reduction.", warn=True)
        return is_variable
    except Exception:
        return False

def build_comp_ra_dec(ra_wcs, dec_wcs, comp_stars):
    comp_ra_dec = []
    for _, comp_star in enumerate(comp_stars[:]):
        comp_ra_dec.append([ra_wcs[int(comp_star[1])][int(comp_star[0])],
            dec_wcs[int(comp_star[1])][int(comp_star[0])]])
    return comp_ra_dec


@retry(
    stop=stop_after_attempt(NEXTASTRO_VARIABILITY_MAX_RETRY_ATTEMPTS),
    wait=wait_fixed(NEXTASTRO_VARIABILITY_RETRY_WAIT_SECONDS),
    retry=retry_if_exception(should_retry_nextastro_variability_error),
)
def nextastro_variability_test(comp_ra_dec):
    api_url = 'https://photometry.nextastro.org/variability_test'

    payload = [{'ra': float(ra), 'dec': float(dec)} for ra, dec in comp_ra_dec]
    log_info(f"NextAstro variability request JSON: {json.dumps(payload)}")
    result, content_encoding = submit_nextastro_variability_request(api_url, payload)
    if result.status_code == 415 and content_encoding == 'zstd':
        log_info(
            "NextAstro variability server rejected zstd-compressed request (HTTP 415); "
            "retrying this request once with gzip.",
            warn=True,
        )
        result, content_encoding = submit_nextastro_variability_request(
            api_url,
            payload,
            content_encoding='gzip',
        )
    if result.status_code != 200:
        raise RuntimeError(f"NextAstro variability server returned HTTP {result.status_code}.")

    body = result.json()
    log_info(f"NextAstro variability response JSON: {json.dumps(body)}")
    if not isinstance(body, list) or len(body) != len(payload):
        raise RuntimeError("NextAstro variability server returned an unexpected response format.")

    variability_flags = []
    for index, star in enumerate(body):
        is_in_vsx = int(star.get('is_in_vsx', 0))
        if is_in_vsx not in [0, 1]:
            raise RuntimeError(f"Unexpected is_in_vsx value ({is_in_vsx}) for star index {index}.")
        variability_flags.append(bool(is_in_vsx))

    return variability_flags


def _finite_float(value, default=None):
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return parsed if np.isfinite(parsed) else default


def normalize_nextastro_filter_key(obs_filter):
    return re.sub(r"[^a-z0-9]", "", str(obs_filter or "").lower())


def nextastro_photometry_band_candidates(obs_filter):
    filter_key = normalize_nextastro_filter_key(obs_filter)
    direct_map = {
        'u': [('umag', 'err_umag', 'u')],
        'johnsonu': [('umag', 'err_umag', 'u')],
        'su': [('umag', 'err_umag', 'u')],
        'up': [('umag', 'err_umag', 'u')],
        'b': [('Bmag', 'err_Bmag', 'B')],
        'johnsonb': [('Bmag', 'err_Bmag', 'B')],
        'photographicb': [('Bmag', 'err_Bmag', 'B')],
        'bb': [('Bmag', 'err_Bmag', 'B')],
        'pb': [('Bmag', 'err_Bmag', 'B')],
        'v': [('Vmag', 'err_Vmag', 'V')],
        'johnsonv': [('Vmag', 'err_Vmag', 'V')],
        'bv': [('Vmag', 'err_Vmag', 'V')],
        'cv': [('Vmag', 'err_Vmag', 'V')],
        'clearv': [('Vmag', 'err_Vmag', 'V')],
        'clearunfilteredreducedtovsequence': [('Vmag', 'err_Vmag', 'V')],
        'mobscv': [('Vmag', 'err_Vmag', 'V')],
        'c': [('Vmag', 'err_Vmag', 'V')],
        'clear': [('Vmag', 'err_Vmag', 'V')],
        'lum': [('Vmag', 'err_Vmag', 'V')],
        'luminance': [('Vmag', 'err_Vmag', 'V')],
        'sg': [('g', 'dg', 'g')],
        'sloang': [('g', 'dg', 'g')],
        'sdssg': [('g', 'dg', 'g')],
        'photographicg': [('g', 'dg', 'g')],
        'gp': [('g', 'dg', 'g')],
        'g': [('g', 'dg', 'g')],
        'pg': [('g', 'dg', 'g')],
        'tg': [('g', 'dg', 'g')],
        'sr': [('r', 'dr', 'r')],
        'sloanr': [('r', 'dr', 'r')],
        'sdssr': [('r', 'dr', 'r')],
        'johnsonr': [('r', 'dr', 'r')],
        'cousinsr': [('r', 'dr', 'r')],
        'clearunfilteredreducedtorsequence': [('r', 'dr', 'r')],
        'photographicr': [('r', 'dr', 'r')],
        'rp': [('r', 'dr', 'r')],
        'r': [('r', 'dr', 'r')],
        'rc': [('r', 'dr', 'r')],
        'rj': [('r', 'dr', 'r')],
        'pr': [('r', 'dr', 'r')],
        'tr': [('r', 'dr', 'r')],
        'cr': [('r', 'dr', 'r')],
        'si': [('i', 'di', 'i')],
        'sloani': [('i', 'di', 'i')],
        'sdssi': [('i', 'di', 'i')],
        'johnsoni': [('i', 'di', 'i')],
        'cousinsi': [('i', 'di', 'i')],
        'ip': [('i', 'di', 'i')],
        'i': [('i', 'di', 'i')],
        'ic': [('i', 'di', 'i')],
        'ij': [('i', 'di', 'i')],
        'sz': [('z', 'dz', 'z')],
        'sloanz': [('z', 'dz', 'z')],
        'sdssz': [('z', 'dz', 'z')],
        'panstarrszshort': [('z', 'dz', 'z')],
        'zp': [('z', 'dz', 'z')],
        'z': [('z', 'dz', 'z')],
        'zs': [('z', 'dz', 'z')],
    }

    fallback = [
        ('Vmag', 'err_Vmag', 'V'),
        ('g', 'dg', 'g'),
        ('r', 'dr', 'r'),
        ('i', 'di', 'i'),
        ('Bmag', 'err_Bmag', 'B'),
        ('z', 'dz', 'z'),
        ('umag', 'err_umag', 'u'),
    ]
    candidates = list(direct_map.get(filter_key, []))
    candidates.extend(candidate for candidate in fallback if candidate not in candidates)
    return candidates


def nextastro_catalog_rows(catalog_response):
    if not isinstance(catalog_response, dict):
        return []
    rows = catalog_response.get('rows', [])
    if not isinstance(rows, list):
        return []
    columns = catalog_response.get('columns', [])
    if catalog_response.get('row_format') == 'arrays':
        return [
            {column: row[index] if index < len(row) else None for index, column in enumerate(columns)}
            for row in rows
            if isinstance(row, list)
        ]
    return [row for row in rows if isinstance(row, dict)]


def row_nextastro_magnitude(row, band_candidates):
    for priority, (mag_column, error_column, band_label) in enumerate(band_candidates):
        magnitude = _finite_float(row.get(mag_column))
        magnitude_error = _finite_float(row.get(error_column))
        if magnitude is None or magnitude_error is None:
            continue
        return {
            'priority': priority,
            'mag': magnitude,
            'error': abs(magnitude_error),
            'mag_band': band_label,
            'mag_column': mag_column,
            'mag_error_column': error_column,
        }
    return None


def sky_separation_arcsec(ra_a, dec_a, ra_b, dec_b):
    first = SkyCoord(float(ra_a) * u.deg, float(dec_a) * u.deg, frame='fk5')
    second = SkyCoord(float(ra_b) * u.deg, float(dec_b) * u.deg, frame='fk5')
    return float(first.separation(second).arcsec)


def nextastro_photometry_catalog_match(catalog_response, ra, dec, obs_filter,
                                       max_separation_arcsec=NEXTASTRO_PHOTOMETRY_MATCH_RADIUS_ARCSEC):
    band_candidates = nextastro_photometry_band_candidates(obs_filter)
    matches = []
    for row in nextastro_catalog_rows(catalog_response):
        row_ra = _finite_float(row.get('ra'))
        row_dec = _finite_float(row.get('dec'))
        if row_ra is None or row_dec is None:
            continue
        magnitude = row_nextastro_magnitude(row, band_candidates)
        if magnitude is None:
            continue
        separation = sky_separation_arcsec(ra, dec, row_ra, row_dec)
        if separation > max_separation_arcsec:
            continue
        matches.append({
            **magnitude,
            'catalog_ra': row_ra,
            'catalog_dec': row_dec,
            'source_id': row.get('source_id'),
            'id': row.get('id'),
            'separation_arcsec': separation,
        })

    if not matches:
        return None
    matches.sort(key=lambda match: (match['priority'], match['separation_arcsec']))
    return matches[0]


@retry(
    stop=stop_after_attempt(NEXTASTRO_VARIABILITY_MAX_RETRY_ATTEMPTS),
    wait=wait_fixed(NEXTASTRO_VARIABILITY_RETRY_WAIT_SECONDS),
    retry=retry_if_exception(should_retry_nextastro_variability_error),
)
def nextastro_photometry_cone_query(ra, dec, radius_arcsec, columns=None):
    api_url = f'{NEXTASTRO_PHOTOMETRY_API_URL}/cone_query'
    payload = {
        'columns': list(columns or NEXTASTRO_PHOTOMETRY_COLUMNS),
        'ra': float(ra),
        'dec': float(dec),
        'radius_arcsec': float(radius_arcsec),
    }
    log_info(f"NextAstro photometry catalog request JSON: {json.dumps(payload)}")
    result = requests.post(api_url, json=payload, timeout=30)
    if result.status_code != 200:
        raise RuntimeError(f"NextAstro photometry catalog returned HTTP {result.status_code}.")

    body = result.json()
    if not isinstance(body, dict) or not isinstance(body.get('rows'), list):
        raise RuntimeError("NextAstro photometry catalog returned an unexpected response format.")
    log_info(
        "NextAstro photometry catalog response JSON: "
        f"{json.dumps({'count': body.get('count'), 'columns': body.get('columns')})}"
    )
    return body


def nextastro_photometry_catalog_for_wcs(wcs_file, axis, img_scale, obs_filter):
    if not wcs_file or img_scale is None:
        return None
    image_width, image_height = float(axis[0]), float(axis[1])
    if not (np.isfinite(image_width) and np.isfinite(image_height) and np.isfinite(float(img_scale))):
        return None

    wcs_hdr = search_wcs(wcs_file)
    center_ra, center_dec = wcs_hdr.pixel_to_world_values(image_width / 2.0, image_height / 2.0)
    radius_arcsec = 0.5 * float(img_scale) * float(np.hypot(image_width, image_height))
    radius_arcsec += NEXTASTRO_PHOTOMETRY_FIELD_PADDING_ARCSEC
    log_info(
        "\nQuerying NextAstro photometry catalog for the full reduced field "
        f"(radius={radius_arcsec:.1f} arcsec)."
    )
    return nextastro_photometry_cone_query(center_ra, center_dec, radius_arcsec)


def nextastro_photometry_for_coordinate(ra, dec, obs_filter,
                                        radius_arcsec=NEXTASTRO_PHOTOMETRY_MATCH_RADIUS_ARCSEC):
    catalog_response = nextastro_photometry_cone_query(ra, dec, radius_arcsec)
    return nextastro_photometry_catalog_match(
        catalog_response,
        ra,
        dec,
        obs_filter,
        max_separation_arcsec=radius_arcsec,
    )


def nextastro_calibration_label(match):
    source_id = match.get('source_id') or match.get('id')
    if source_id not in (None, ''):
        return f"NextAstro-{source_id}"
    return f"RA{match['ra']:.6f}_DEC{match['dec']:.6f}"


def merge_nextastro_calibration_stars(comp_stars, comp_ra_dec, obs_filter, existing_comp_stars=None,
                                      field_catalog=None):
    calibration_stars = dict(existing_comp_stars or {})
    existing_positions = {
        tuple(value.get('pos', []))
        for value in calibration_stars.values()
        if isinstance(value, dict)
    }

    added_count = 0
    for index, (comp_pos, comp_radec) in enumerate(zip(comp_stars, comp_ra_dec)):
        if tuple(comp_pos) in existing_positions:
            continue
        comp_ra = _finite_float(comp_radec[0])
        comp_dec = _finite_float(comp_radec[1])
        if comp_ra is None or comp_dec is None:
            continue

        match = None
        if field_catalog is not None:
            match = nextastro_photometry_catalog_match(field_catalog, comp_ra, comp_dec, obs_filter)
        if match is None:
            try:
                match = nextastro_photometry_for_coordinate(comp_ra, comp_dec, obs_filter)
            except Exception as exc:
                log_info(
                    f"Warning: NextAstro photometry catalog lookup failed for comparison star #{index + 1} "
                    f"({describe_retry_exception(exc)}).",
                    warn=True,
                )
                continue
        if match is None:
            log_info(
                f"Warning: NextAstro photometry catalog did not find a usable magnitude for "
                f"comparison star #{index + 1}.",
                warn=True,
            )
            continue

        match.update({
            'ra': comp_ra,
            'dec': comp_dec,
            'pos': list(comp_pos),
            'catalog_source': 'NextAstro photometry catalog',
            'is_aavso_vsp': False,
            'observed_filter': obs_filter,
        })
        label = nextastro_calibration_label(match)
        unique_label = label
        duplicate_index = 2
        while unique_label in calibration_stars:
            unique_label = f"{label}-{duplicate_index}"
            duplicate_index += 1
        calibration_stars[unique_label] = match
        existing_positions.add(tuple(comp_pos))
        added_count += 1
        log_info(
            f"NextAstro photometry calibration for comparison star #{index + 1}: "
            f"{match['mag_band']}={match['mag']:.5f} +/- {match['error']:.5f}, "
            f"RA={comp_ra:.7f}, Dec={comp_dec:.7f}, "
            f"catalog separation={match['separation_arcsec']:.2f} arcsec."
        )

    if added_count:
        log_info(f"Added {added_count} NextAstro photometry catalog comparison star calibration(s).")
    return calibration_stars


def nextastro_prereduced_calibration_star(phot_comp_star, obs_filter):
    if not isinstance(phot_comp_star, dict):
        return None, None

    comp_ra = _finite_float(phot_comp_star.get('ra'))
    comp_dec = _finite_float(phot_comp_star.get('dec'))
    if comp_ra is None or comp_dec is None:
        return None, None

    match = nextastro_photometry_for_coordinate(comp_ra, comp_dec, obs_filter)
    if match is None:
        return None, None

    match.update({
        'ra': comp_ra,
        'dec': comp_dec,
        'pos': [
            phot_comp_star.get('x', ''),
            phot_comp_star.get('y', ''),
        ],
        'catalog_source': 'NextAstro photometry catalog',
        'is_aavso_vsp': False,
        'observed_filter': obs_filter,
    })
    label = nextastro_calibration_label(match)
    log_info(
        "NextAstro photometry calibration for pre-reduced comparison star: "
        f"{match['mag_band']}={match['mag']:.5f} +/- {match['error']:.5f}, "
        f"RA={comp_ra:.7f}, Dec={comp_dec:.7f}, "
        f"catalog separation={match['separation_arcsec']:.2f} arcsec."
    )
    return label, match


def check_for_variable_stars(ra_wcs, dec_wcs, comp_stars, use_nextastro_variability_server=False):
    if use_nextastro_variability_server and comp_stars:
        try:
            log_info("\nChecking for variability using NextAstro variability server.")
            comp_ra_dec = build_comp_ra_dec(ra_wcs, dec_wcs, comp_stars)
            variability_flags = nextastro_variability_test(comp_ra_dec)

            for i, (comp_star, is_variable) in enumerate(zip(comp_stars[:], variability_flags)):
                log_info(f"\nChecking for variability in Comparison Star #{i + 1}:"
                         f"\n\tPixel X: {comp_star[0]} Pixel Y: {comp_star[1]}"
                         f"\n\tNextAstro flagged variable: {is_variable}")
                if is_variable:
                    comp_stars.remove(comp_star)
            return
        except Exception as e:
            log_info(f"\nWarning: NextAstro variability server check failed ({describe_retry_exception(e)}). "
                     "Falling back to individual VSX variability checks.", warn=True)

    for i, comp_star in enumerate(comp_stars[:]):
        ra = ra_wcs[int(comp_star[1])][int(comp_star[0])]
        dec = dec_wcs[int(comp_star[1])][int(comp_star[0])]

        log_info(f"\nChecking for variability in Comparison Star #{i + 1}:"
                 f"\n\tPixel X: {comp_star[0]} Pixel Y: {comp_star[1]}")
        if query_variable_star_apis(ra, dec):
            comp_stars.remove(comp_star)

# Apply calibrations if applicable
def apply_cals(image_data, gen_dark, gen_bias, gen_flat, i):
    if gen_dark is not None and gen_dark.size != 0:
        if i == 0:
            log_info("Dark subtracting images.")
        image_data = image_data - gen_dark
    elif gen_bias is not None and gen_bias.size != 0:  # if a dark is not available, then at least subtract off the pedestal via the bias
        if i == 0:
            log_info("Bias-correcting images.")
        image_data = image_data - gen_bias
    else:
        pass

    if gen_flat is not None and gen_flat.size != 0:
        if i == 0:
            log_info("Flattening images.")
        gen_flat[gen_flat == 0] = 1
        image_data = image_data / gen_flat
    return image_data

def calculate_demosaic_mult(demosaic_out): 
    if not demosaic_out:
        return None       
    # Build vector to convert RBG pixels to single output
    if isinstance(demosaic_out, list):
        demosaic_mult = np.array(demosaic_out)
    elif demosaic_out == 'red':
        demosaic_mult = np.array([ 1.0, 0.0, 0.0 ])
    elif demosaic_out == 'green':
        demosaic_mult = np.array([ 0.0, 1.0, 0.0 ])
    elif demosaic_out == 'blue':
        demosaic_mult = np.array([ 0.0, 0.0, 1.0 ])
    elif demosaic_out == 'gray':
        demosaic_mult = np.array([ 0.299, 0.587, 0.114 ])   # Same as rbg2gray
    elif demosaic_out == 'blueblock':
        demosaic_mult = np.array([ 0.299, 0.587, 0.0 ]) # drop blue, same mix of red, green as gray
    else:   # Green default
        demosaic_mult = np.array([ 0.0, 1.0, 0.0 ])
    # Normalize
    demosaic_mult = demosaic_mult / (demosaic_mult[0]+demosaic_mult[1]+demosaic_mult[2])
    return demosaic_mult

# If demosaic requested, process
def demosaic_img(image_data, demosaic_fmt, demosaic_out, demosaic_mult, i):
    if demosaic_fmt:
        if i == 0:
            log_info(f"Demosaicing images (mapping {demosaic_fmt} to {demosaic_out})")
        img_dtype = image_data.dtype    # Save data type
        new_image_data = demosaicing_CFA_Bayer_bilinear(image_data, demosaic_fmt)
        image_data = (new_image_data @ demosaic_mult).astype(img_dtype)
    return image_data

def vsp_query(file, axis, obs_filter, img_scale, maglimit=14, user_comp_stars=None, user_targ_star=None):
    if user_comp_stars is None:
        user_comp_stars = []

    vsp_comp_stars_info = {}
    vsp_star_count = 0
    observed_filter = obs_filter

    # Build combined list for comps and target - there are known cases when AAVsO comps have planets (XO-2 N)
    # Plus, we don't want comp too close to target
    targ_and_comp_stars = user_comp_stars[:]
    if user_targ_star is not None:
        targ_and_comp_stars.append(user_targ_star)

    wcs_hdr = search_wcs(file)
    fov = (img_scale * max(axis)) / 60
    ra, dec = wcs_hdr.pixel_to_world_values(axis[0] // 2, axis[1] // 2)
    # Respect limits from AAVSO API (as reported by API error messages)
    if fov > 180 and maglimit > 12:
        maglimit = 12

    url = f"https://www.aavso.org/apps/vsp/api/chart/?format=json&ra={ra:5f}&dec={dec:5f}&fov={fov}&maglimit={maglimit}"
    result = requests.get(url)
    data = result.json()
    chart_id = data['chartid']

    if obs_filter == "CV":
        obs_filter = "V"
    elif obs_filter == "R":
        obs_filter = "Rc"

    if data['photometry']:
        for star in data['photometry']:
            ra_deg, dec_deg = radec_hours_to_degree(star['ra'], star['dec'])
            ra_pix, dec_pix = wcs_hdr.world_to_pixel_values(ra_deg, dec_deg)

            if (ra_pix < axis[0] and dec_pix < axis[1]) and (ra_pix > 1 and dec_pix > 1):
                vsp_star = [int(ra_pix.min()), int(dec_pix.min())]
                exist, vsp_star = check_comp_star_exists(targ_and_comp_stars, vsp_star)

                if obs_filter in [band['band'] for band in star['bands']]:
                    star_info = next(band for band in star['bands'] if band['band'] == obs_filter)

                    vsp_comp_stars_info[star['auid']] = {
                        'pos': vsp_star,
                        'mag': star_info['mag'],
                        'error': star_info['error'],
                        'ra': ra_deg,
                        'dec': dec_deg,
                        'catalog_ra': ra_deg,
                        'catalog_dec': dec_deg,
                        'mag_band': obs_filter,
                        'observed_filter': observed_filter,
                        'catalog_source': 'AAVSO VSP',
                        'is_aavso_vsp': True,
                    }

                    if not exist:
                        vsp_star_count = add_vsp_star(vsp_star_count, user_comp_stars, vsp_star)

            if len(vsp_comp_stars_info) > 1:
                break

    if not vsp_star_count:
        log_info("\nNo comparison stars were gathered from AAVSO.\n")

    return vsp_comp_stars_info, chart_id


def add_vsp_star(vsp_star_count, user_comp_stars, vsp_star):
    user_comp_stars.append(vsp_star)
    log_info(f"\nAdded Comparison Star #{len(user_comp_stars)}, coordinates {vsp_star} from AAVSO")

    return vsp_star_count + 1


def check_comp_star_exists(user_stars, vsp_star, tol=10):
    """Checks if a comparison star from VSP exists in the user-entered
    comparison star list

    Parameters
    ----------
    user_stars : list
        A header file that may include the airmass or altitude from when the image was taken
    vsp_star : list
        Right Ascension
    tol : float
        Declination

    Returns
    -------
    bool
        True if VSP star exists in user entered stars, otherwise False
    list
        Pixel coordinate of either the user entered star (exists), otherwise pixel coordinates
        of VSP
    """
    for user_star in user_stars:
        pixel_distance = [abs(star1 - star2) for star1, star2 in zip(user_star, vsp_star)]

        if all(i <= tol for i in pixel_distance):
            return True, user_star
    return False, vsp_star



TRANSFORM_TIMING_STAGES = [
    'astroalign_direct',
    'fft_translation',
    'astroalign_filtered',
    'astroalign_mask',
    'imreg_dft',
]

_TRANSFORM_TIMING_STATS = {
    stage: {'count': 0, 'success': 0, 'total_s': 0.0} for stage in TRANSFORM_TIMING_STAGES
}
_TRANSFORM_TIMING_STATS['mask_loops_skipped'] = 0

PHOTOMETRY_TIMING_STAGES = ['fit_centroid', 'aperPhot']
_PHOTOMETRY_TIMING_STATS = {
    stage: {'count': 0, 'total_s': 0.0} for stage in PHOTOMETRY_TIMING_STAGES
}


def reset_transform_timing_stats():
    for stage in TRANSFORM_TIMING_STAGES:
        _TRANSFORM_TIMING_STATS[stage] = {'count': 0, 'success': 0, 'total_s': 0.0}
    _TRANSFORM_TIMING_STATS['mask_loops_skipped'] = 0


def reset_photometry_timing_stats():
    for stage in PHOTOMETRY_TIMING_STAGES:
        _PHOTOMETRY_TIMING_STATS[stage] = {'count': 0, 'total_s': 0.0}


def _record_transform_stage_timing(stage, elapsed_s, success):
    stage_stats = _TRANSFORM_TIMING_STATS[stage]
    stage_stats['count'] += 1
    stage_stats['total_s'] += elapsed_s
    if success:
        stage_stats['success'] += 1


def _record_photometry_stage_timing(stage, elapsed_s):
    stage_stats = _PHOTOMETRY_TIMING_STATS[stage]
    stage_stats['count'] += 1
    stage_stats['total_s'] += elapsed_s


def log_transform_timing_stats(prefix='Transformation timing summary'):
    logged_any = False
    lines = []

    for stage in TRANSFORM_TIMING_STAGES:
        stage_stats = _TRANSFORM_TIMING_STATS[stage]
        if stage_stats['count'] == 0:
            continue

        avg_ms = 1000.0 * stage_stats['total_s'] / stage_stats['count']
        lines.append(
            f"{stage}: calls={stage_stats['count']}, success={stage_stats['success']}, "
            f"avg_ms={avg_ms:.2f}, total_s={stage_stats['total_s']:.2f}"
        )
        logged_any = True

    if _TRANSFORM_TIMING_STATS['mask_loops_skipped']:
        lines.append(f"astroalign_mask_loops_skipped={_TRANSFORM_TIMING_STATS['mask_loops_skipped']}")
        logged_any = True

    if logged_any:
        log_info(f"{prefix}: " + " | ".join(lines))


def log_photometry_timing_stats(prefix='Photometry timing summary'):
    logged_any = False
    lines = []

    for stage in PHOTOMETRY_TIMING_STAGES:
        stage_stats = _PHOTOMETRY_TIMING_STATS[stage]
        if stage_stats['count'] == 0:
            continue

        avg_ms = 1000.0 * stage_stats['total_s'] / stage_stats['count']
        lines.append(
            f"{stage}: calls={stage_stats['count']}, avg_ms={avg_ms:.2f}, total_s={stage_stats['total_s']:.2f}"
        )
        logged_any = True

    if logged_any:
        log_info(f"{prefix}: " + " | ".join(lines))


def log_reduction_timing_overview(prefix='Reduction timing overview'):
    transform_total_s = sum(_TRANSFORM_TIMING_STATS[stage]['total_s'] for stage in TRANSFORM_TIMING_STAGES)
    photometry_total_s = sum(_PHOTOMETRY_TIMING_STATS[stage]['total_s'] for stage in PHOTOMETRY_TIMING_STAGES)
    combined_total_s = transform_total_s + photometry_total_s
    if combined_total_s <= 0:
        return

    dominant_bucket = 'transform'
    if photometry_total_s > transform_total_s:
        dominant_bucket = 'photometry'

    transform_pct = 100.0 * transform_total_s / combined_total_s
    photometry_pct = 100.0 * photometry_total_s / combined_total_s
    log_info(
        f"{prefix}: transform_total_s={transform_total_s:.2f} ({transform_pct:.1f}%), "
        f"photometry_total_s={photometry_total_s:.2f} ({photometry_pct:.1f}%), "
        f"dominant={dominant_bucket}"
    )


def _display_filename(file_name):
    return str(file_name).replace("\\", "/").rsplit("/", 1)[-1]


def format_plate_solution_reference(wcs_file):
    return f"Here is the filename where we got the WCS from: {_display_filename(wcs_file)}"


# Aligns imaging data from .fits file to easily track the host and comparison star's positions
def transformation(image_data, file_name, roi=1, report_failure=True, reference_image=None):
    start_time = perf_counter()
    display_file_name = _display_filename(file_name)

    if report_failure:
        plateStatus.setCurrentFilename(file_name)

    # crop image to ROI
    if reference_image is None:
        current_image = image_data[0]
        reference_image = image_data[1]
    else:
        current_image = image_data

    reference_cache = _get_reference_transform_cache(reference_image, roi)
    roix = reference_cache['roix']
    roiy = reference_cache['roiy']
    roi_reference = reference_cache['roi_reference']
    roi_current = current_image[roiy, roix]

    if roi_reference.shape != roi_current.shape or roi_reference.size == 0:
        log.debug(
            f"Warning: Following image failed pre-alignment checks in "
            f"{perf_counter() - start_time:.2f}s - {display_file_name}"
        )
        if report_failure:
            plateStatus.alignmentError()
        return SimilarityTransform(scale=1, rotation=0, translation=[0, 0])

    fft_tform = None

    # Fast FFT translation estimate before more expensive fallback stages.
    # Most cadence images are dominated by small translations, so this stage
    # can often solve alignment without invoking significantly slower
    # feature-matching methods.
    stage_start = perf_counter()
    try:
        shift, error, _ = phase_cross_correlation(roi_current, roi_reference, upsample_factor=4)
        if np.all(np.isfinite(shift)) and np.isfinite(error):
            max_shift = max(abs(shift[0]), abs(shift[1]))
            if max_shift <= max(roi_current.shape):
                fft_tform = SimilarityTransform(scale=1, rotation=0, translation=[-shift[1], -shift[0]])
                fft_high_confidence = error <= 0.1 and max_shift <= max(roi_current.shape) * 0.25
                _record_transform_stage_timing('fft_translation', perf_counter() - stage_start, True)
                if fft_high_confidence:
                    log.debug(
                        f"Transformation solved via high-confidence FFT in "
                        f"{perf_counter() - start_time:.2f}s for {display_file_name}"
                    )
                    return fft_tform
            else:
                _record_transform_stage_timing('fft_translation', perf_counter() - stage_start, False)
        else:
            _record_transform_stage_timing('fft_translation', perf_counter() - stage_start, False)
    except Exception:
        _record_transform_stage_timing('fft_translation', perf_counter() - stage_start, False)

    # Find transformation from .FITS files and catch exceptions if not able to.
    stage_start = perf_counter()
    try:
        results = aa.find_transform(roi_reference, roi_current)
        _record_transform_stage_timing('astroalign_direct', perf_counter() - stage_start, True)
        log.debug(
            f"Transformation solved via astroalign direct pass in "
            f"{perf_counter() - start_time:.2f}s for {display_file_name}"
        )
        return results[0]
    except Exception:
        _record_transform_stage_timing('astroalign_direct', perf_counter() - stage_start, False)

    # One cheap filtered pass to suppress noise and retry astroalign.
    filtered_current = gaussian_filter(roi_current, sigma=1.0)
    filtered_reference = reference_cache['filtered_reference']

    stage_start = perf_counter()
    try:
        results = aa.find_transform(filtered_reference, filtered_current)
        _record_transform_stage_timing('astroalign_filtered', perf_counter() - stage_start, True)
        log.debug(
            f"Transformation solved via filtered astroalign in "
            f"{perf_counter() - start_time:.2f}s for {display_file_name}"
        )
        return results[0]
    except Exception:
        _record_transform_stage_timing('astroalign_filtered', perf_counter() - stage_start, False)

    for p in [99, 98, 95, 90]:
        base_mask1 = reference_cache['reference_masks'][p]
        p_cur = np.percentile(roi_current, p)
        base_mask0 = roi_current > p_cur

        for it in [2, 1, 0]:
            # create binary mask to align image
            mask1 = base_mask1
            mask0 = base_mask0

            if it > 0:
                mask1 = binary_erosion(mask1, iterations=it)
                mask0 = binary_erosion(mask0, iterations=it)

            stage_start = perf_counter()
            try:
                results = aa.find_transform(mask1, mask0)
                _record_transform_stage_timing('astroalign_mask', perf_counter() - stage_start, True)
                log.debug(
                    f"Transformation solved via mask astroalign (p={p}, erode={it}) in "
                    f"{perf_counter() - start_time:.2f}s for {display_file_name}"
                )
                return results[0]
            except Exception:
                _record_transform_stage_timing('astroalign_mask', perf_counter() - stage_start, False)

    stage_start = perf_counter()
    try:
        result1 = ird.similarity(roi_reference, roi_current, numiter=3)
        _record_transform_stage_timing('imreg_dft', perf_counter() - stage_start, True)
        log.debug(
            f"Transformation solved via imreg_dft fallback in "
            f"{perf_counter() - start_time:.2f}s for {display_file_name}"
        )
        return SimilarityTransform(scale=result1['scale'], rotation=np.radians(result1['angle']),
                                   translation=[-1 * result1['tvec'][1], -1 * result1['tvec'][0]])
    except Exception:
        _record_transform_stage_timing('imreg_dft', perf_counter() - stage_start, False)

    if fft_tform is not None:
        log.debug(
            f"Transformation fell back to FFT translation in "
            f"{perf_counter() - start_time:.2f}s for {display_file_name}"
        )
        return fft_tform

    log.debug(
        f"Warning: Following image failed to align in "
        f"{perf_counter() - start_time:.2f}s - {display_file_name}"
    )
    if report_failure:
        plateStatus.alignmentError()
    return SimilarityTransform(scale=1, rotation=0, translation=[0, 0])

def load_image_data(file_name):
    hdul = fits.open(name=file_name, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
    extension = 0
    image_header = hdul[extension].header
    while image_header["NAXIS"] == 0:
        extension += 1
        image_header = hdul[extension].header

    image_data = hdul[extension].data
    hdul.close()
    return image_data


def persistent_bad_pixel_count_threshold(frame_count, minimum_fraction=BAD_PIXEL_DETECTION_FRACTION):
    if frame_count <= 0:
        return 1
    return max(1, int(np.floor(float(minimum_fraction) * frame_count)) + 1)


def detect_frame_bad_pixels(image_data,
                            outlier_sigma=BAD_PIXEL_OUTLIER_SIGMA,
                            isolation_sigma=BAD_PIXEL_ISOLATION_SIGMA,
                            isolation_ratio=BAD_PIXEL_ISOLATION_RATIO):
    values = np.asarray(image_data, dtype=float)
    if values.ndim != 2 or values.size == 0:
        return np.zeros(values.shape[:2], dtype=bool)

    finite_mask = np.isfinite(values)
    if np.count_nonzero(finite_mask) < BAD_PIXEL_NEIGHBOR_FOOTPRINT.sum():
        return np.zeros(values.shape, dtype=bool)

    working = np.array(values, copy=True)
    frame_median = bn.nanmedian(working[finite_mask])
    if not np.isfinite(frame_median):
        frame_median = 0.0
    working[~finite_mask] = frame_median

    neighbor_median = median_filter(working, footprint=BAD_PIXEL_NEIGHBOR_FOOTPRINT, mode='mirror')
    neighbor_max = maximum_filter(working, footprint=BAD_PIXEL_NEIGHBOR_FOOTPRINT, mode='mirror')
    residual = working - neighbor_median

    global_scatter = robust_scatter(residual[finite_mask])
    if not np.isfinite(global_scatter) or global_scatter <= 0:
        global_scatter = robust_scatter(working[finite_mask])
    if not np.isfinite(global_scatter) or global_scatter <= 0:
        return np.zeros(values.shape, dtype=bool)

    local_scatter = 1.4826 * median_filter(
        np.abs(residual),
        footprint=BAD_PIXEL_NEIGHBOR_FOOTPRINT,
        mode='mirror',
    )
    diff_threshold = np.maximum(outlier_sigma * local_scatter, BAD_PIXEL_GLOBAL_SIGMA * global_scatter)
    isolation_threshold = max(isolation_sigma * global_scatter, 1.0)
    neighbor_scale = np.maximum(np.abs(neighbor_max), 1.0)

    with np.errstate(divide='ignore', invalid='ignore'):
        isolation_ratio_values = np.divide(np.abs(working), neighbor_scale)

    return (
        finite_mask
        & (residual > diff_threshold)
        & ((working - neighbor_max) > isolation_threshold)
        & (isolation_ratio_values >= isolation_ratio)
    )


_BAD_PIXEL_PRECHECK_POOL_CONTEXT = {}


def _bad_pixel_precheck_pool_initializer(generalDark, generalBias, generalFlat,
                                         demosaic_fmt, demosaic_out, demosaic_mult):
    global _BAD_PIXEL_PRECHECK_POOL_CONTEXT
    suppress_inherited_tk_cleanup_in_worker()
    _BAD_PIXEL_PRECHECK_POOL_CONTEXT = {
        'generalDark': generalDark,
        'generalBias': generalBias,
        'generalFlat': generalFlat,
        'demosaic_fmt': demosaic_fmt,
        'demosaic_out': demosaic_out,
        'demosaic_mult': demosaic_mult,
    }


def _load_bad_pixel_precheck_worker_frame(file_name):
    context = _BAD_PIXEL_PRECHECK_POOL_CONTEXT
    hdul = fits.open(name=file_name, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
    extension = 0
    image_header = hdul[extension].header
    while image_header["NAXIS"] == 0:
        extension += 1
        image_header = hdul[extension].header

    image_data = hdul[extension].data
    hdul.close()

    image_data = apply_cals(
        image_data,
        context.get('generalDark'),
        context.get('generalBias'),
        context.get('generalFlat'),
        1,
    )
    image_data = demosaic_img(
        image_data,
        context.get('demosaic_fmt'),
        context.get('demosaic_out'),
        context.get('demosaic_mult'),
        1,
    )
    return image_data


def _bad_pixel_precheck_task(task):
    index, file_name = task
    try:
        frame_data = _load_bad_pixel_precheck_worker_frame(file_name)
    except Exception as exc:
        return {
            'index': index,
            'file_name': file_name,
            'usable': False,
            'error': str(exc),
        }

    frame_mask = detect_frame_bad_pixels(frame_data)
    if frame_mask.ndim != 2:
        return {
            'index': index,
            'file_name': file_name,
            'usable': False,
            'not_2d': True,
        }

    return {
        'index': index,
        'file_name': file_name,
        'usable': True,
        'mask': frame_mask,
    }


def _merge_bad_pixel_precheck_mask(detection_counts, frame_mask, file_name):
    frame_mask = np.asarray(frame_mask, dtype=bool)
    if frame_mask.ndim != 2:
        log_info(
            f"Warning: skipping bad-pixel precheck for {_display_filename(file_name)} because the frame is not 2-D.",
            warn=True,
        )
        return detection_counts, False

    if detection_counts is None:
        detection_counts = np.zeros(frame_mask.shape, dtype=np.uint32)
    elif detection_counts.shape != frame_mask.shape:
        log_info(
            "Warning: skipping bad-pixel precheck for "
            f"{_display_filename(file_name)} because its shape {frame_mask.shape} does not match "
            f"the reference frame shape {detection_counts.shape}.",
            warn=True,
        )
        return detection_counts, False

    detection_counts += frame_mask.astype(np.uint32)
    return detection_counts, True


def _scan_bad_pixel_precheck_frames_serial(inputfiles, frame_loader):
    total_files = len(inputfiles)
    detection_counts = None
    scanned_files = 0

    for index, file_name in enumerate(inputfiles):
        plateStatus.setCurrentFilename(file_name)
        try:
            frame_data = frame_loader(file_name)
        except Exception as exc:
            log_info(
                f"Warning: skipping bad-pixel precheck for {_display_filename(file_name)} ({exc}).",
                warn=True,
            )
            continue

        frame_mask = detect_frame_bad_pixels(frame_data)
        detection_counts, usable = _merge_bad_pixel_precheck_mask(detection_counts, frame_mask, file_name)
        if usable:
            scanned_files += 1

        completed = index + 1
        if completed == total_files or completed % BAD_PIXEL_PROGRESS_LOG_INTERVAL == 0:
            log_info(f"Bad-pixel precheck progress: {completed}/{total_files}")

    return detection_counts, scanned_files


def _scan_bad_pixel_precheck_frames_multiprocess(inputfiles, max_processes,
                                                 generalDark=None, generalBias=None, generalFlat=None,
                                                 demosaic_fmt=None, demosaic_out=None, demosaic_mult=None):
    total_files = len(inputfiles)
    max_workers = min(max_processes, os.cpu_count() or 1, total_files, MAX_MULTIPROCESS_BAD_PIXEL_WORKERS)
    detection_counts = None
    scanned_files = 0

    log_info(
        "Using multiprocessing for bad-pixel precheck "
        f"with {max_workers} worker(s) across {total_files} image(s)."
    )

    tasks = [(index, str(file_name)) for index, file_name in enumerate(inputfiles)]
    with suppress_tk_cleanup_during_process_pool():
        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_bad_pixel_precheck_pool_initializer,
            initargs=(
                generalDark,
                generalBias,
                generalFlat,
                demosaic_fmt,
                demosaic_out,
                demosaic_mult,
            ),
        ) as executor:
            futures = [executor.submit(_bad_pixel_precheck_task, task) for task in tasks]
            completed = 0
            for future in as_completed(futures):
                result = future.result()
                file_name = result.get('file_name')
                if result.get('usable'):
                    detection_counts, usable = _merge_bad_pixel_precheck_mask(
                        detection_counts,
                        result.get('mask'),
                        file_name,
                    )
                    if usable:
                        scanned_files += 1
                elif result.get('not_2d'):
                    log_info(
                        f"Warning: skipping bad-pixel precheck for {_display_filename(file_name)} "
                        "because the frame is not 2-D.",
                        warn=True,
                    )
                else:
                    log_info(
                        f"Warning: skipping bad-pixel precheck for {_display_filename(file_name)} "
                        f"({result.get('error')}).",
                        warn=True,
                    )

                completed += 1
                if completed == total_files or completed % BAD_PIXEL_PROGRESS_LOG_INTERVAL == 0:
                    log_info(f"Bad-pixel precheck progress: {completed}/{total_files}")

    return detection_counts, scanned_files


def build_persistent_bad_pixel_map(inputfiles, frame_loader, save_directory=None,
                                   minimum_fraction=BAD_PIXEL_DETECTION_FRACTION,
                                   minimum_frames=BAD_PIXEL_PRECHECK_MIN_FRAMES,
                                   max_processes=None, generalDark=None, generalBias=None,
                                   generalFlat=None, demosaic_fmt=None, demosaic_out=None,
                                   demosaic_mult=None):
    inputfiles = list(inputfiles)
    total_files = len(inputfiles)
    if total_files < minimum_frames:
        log_info(
            f"Bad-pixel precheck skipped: only {total_files} frame(s); need at least {minimum_frames} frames.",
        )
        return None

    try:
        max_processes = int(max_processes) if max_processes is not None else None
    except (TypeError, ValueError):
        max_processes = None

    if max_processes is not None and max_processes > 1 and total_files > 1:
        try:
            detection_counts, scanned_files = _scan_bad_pixel_precheck_frames_multiprocess(
                inputfiles,
                max_processes,
                generalDark=generalDark,
                generalBias=generalBias,
                generalFlat=generalFlat,
                demosaic_fmt=demosaic_fmt,
                demosaic_out=demosaic_out,
                demosaic_mult=demosaic_mult,
            )
        except Exception as exc:
            log_info(
                f"Warning: bad-pixel precheck multiprocessing failed ({exc}); falling back to serial scanning.",
                warn=True,
            )
            detection_counts, scanned_files = _scan_bad_pixel_precheck_frames_serial(inputfiles, frame_loader)
    else:
        detection_counts, scanned_files = _scan_bad_pixel_precheck_frames_serial(inputfiles, frame_loader)

    if detection_counts is None or scanned_files < minimum_frames:
        log_info(
            f"Bad-pixel precheck skipped: only {scanned_files} usable frame(s); need at least {minimum_frames}.",
            warn=True,
        )
        return None

    required_count = persistent_bad_pixel_count_threshold(scanned_files, minimum_fraction)
    bad_pixel_mask = detection_counts >= required_count
    coord_y, coord_x = np.nonzero(bad_pixel_mask)

    counts_path = None
    mask_path = None
    if save_directory is not None:
        temp_dir = Path(save_directory) / "temp"
        temp_dir.mkdir(parents=True, exist_ok=True)
        counts_path = temp_dir / BAD_PIXEL_COUNTS_FILENAME
        mask_path = temp_dir / BAD_PIXEL_MASK_FILENAME
        fits.writeto(counts_path, detection_counts.astype(np.int32), overwrite=True)
        fits.writeto(mask_path, bad_pixel_mask.astype(np.uint8), overwrite=True)

    threshold_percent = minimum_fraction * 100.0
    summary = (
        f"Bad-pixel precheck: identified {int(np.count_nonzero(bad_pixel_mask))} persistent bad pixel(s) "
        f"after scanning {scanned_files}/{total_files} frame(s) with a >{threshold_percent:g}% recurrence threshold "
        f"({required_count}+ detections)."
    )
    if counts_path is not None and mask_path is not None:
        summary += f" Saved {counts_path.name} and {mask_path.name} to temp/."
    log_info(summary)

    return {
        'count_image': detection_counts,
        'mask': bad_pixel_mask,
        'coord_y': coord_y.astype(int),
        'coord_x': coord_x.astype(int),
        'required_count': required_count,
        'minimum_fraction': float(minimum_fraction),
        'frame_count': scanned_files,
        'counts_path': counts_path,
        'mask_path': mask_path,
    }


def repair_bad_pixels_in_frame(image_data, bad_pixel_reference):
    if bad_pixel_reference is None:
        return image_data

    coord_y = bad_pixel_reference.get('coord_y')
    coord_x = bad_pixel_reference.get('coord_x')
    if coord_y is None or coord_x is None:
        mask = np.asarray(bad_pixel_reference.get('mask'), dtype=bool)
        if mask.size == 0:
            return image_data
        coord_y, coord_x = np.nonzero(mask)

    coord_y = np.asarray(coord_y, dtype=int)
    coord_x = np.asarray(coord_x, dtype=int)
    if coord_y.size == 0 or coord_x.size == 0:
        return image_data

    repaired = np.array(image_data, dtype=float, copy=True)
    valid_coords = (
        (coord_y >= 0) & (coord_y < repaired.shape[0])
        & (coord_x >= 0) & (coord_x < repaired.shape[1])
    )
    if not np.any(valid_coords):
        return repaired

    coord_y = coord_y[valid_coords]
    coord_x = coord_x[valid_coords]
    repaired[coord_y, coord_x] = np.nan

    padded = np.pad(repaired, 1, mode='edge')
    yp = coord_y + 1
    xp = coord_x + 1
    neighbors = np.stack([
        padded[yp - 1, xp - 1],
        padded[yp - 1, xp],
        padded[yp - 1, xp + 1],
        padded[yp, xp - 1],
        padded[yp, xp + 1],
        padded[yp + 1, xp - 1],
        padded[yp + 1, xp],
        padded[yp + 1, xp + 1],
    ], axis=0)

    fill_values = np.nanmedian(neighbors, axis=0)
    if np.any(~np.isfinite(fill_values)):
        frame_median = bn.nanmedian(repaired)
        if not np.isfinite(frame_median):
            frame_median = 0.0
        fill_values[~np.isfinite(fill_values)] = frame_median

    repaired[coord_y, coord_x] = fill_values
    return repaired


def transformation_task(i, file_name, reference_file):
    if i == 0:
        return i, SimilarityTransform(scale=1, rotation=0, translation=[0, 0])

    image_data = load_image_data(file_name)
    reference_image = load_image_data(reference_file)
    # Multiprocess pre-computation should not emit plate-status warnings; the
    # serial reduction path decides whether the fallback transform is needed.
    return i, transformation(image_data, file_name, report_failure=False, reference_image=reference_image)


_TRANSFORM_REFERENCE_IMAGE = None
_TRANSFORM_REFERENCE_CACHE = None


def _build_reference_transform_cache(reference_image, roi):
    height = reference_image.shape[0]
    width = reference_image.shape[1]
    roix = slice(int(width * (0.5 - roi / 2)), int(width * (0.5 + roi / 2)))
    roiy = slice(int(height * (0.5 - roi / 2)), int(height * (0.5 + roi / 2)))

    roi_reference = reference_image[roiy, roix]

    cache = {
        'ref_id': id(reference_image),
        'shape': reference_image.shape,
        'roi': roi,
        'roix': roix,
        'roiy': roiy,
        'roi_reference': roi_reference,
        'filtered_reference': gaussian_filter(roi_reference, sigma=1.0),
    }

    reference_masks = {}
    for p in [99, 98, 95, 90]:
        p_ref = np.percentile(roi_reference, p)
        reference_masks[p] = roi_reference > p_ref
    cache['reference_masks'] = reference_masks

    return cache


def _get_reference_transform_cache(reference_image, roi):
    global _TRANSFORM_REFERENCE_CACHE

    if (_TRANSFORM_REFERENCE_CACHE is None
            or _TRANSFORM_REFERENCE_CACHE['ref_id'] != id(reference_image)
            or _TRANSFORM_REFERENCE_CACHE['shape'] != reference_image.shape
            or _TRANSFORM_REFERENCE_CACHE['roi'] != roi):
        _TRANSFORM_REFERENCE_CACHE = _build_reference_transform_cache(reference_image, roi)

    return _TRANSFORM_REFERENCE_CACHE


def _transformation_pool_initializer(reference_file):
    global _TRANSFORM_REFERENCE_IMAGE, _TRANSFORM_REFERENCE_CACHE
    suppress_inherited_tk_cleanup_in_worker()
    _TRANSFORM_REFERENCE_IMAGE = load_image_data(reference_file)
    _TRANSFORM_REFERENCE_CACHE = None


def transformation_task_with_cached_reference(i, file_name):
    image_data = load_image_data(file_name)
    return i, transformation(image_data, file_name, report_failure=False, reference_image=_TRANSFORM_REFERENCE_IMAGE)


class _ParallelPlateStatusRecorder:
    def __init__(self):
        self.warnings = []

    def setCurrentFilename(self, filename):
        return self

    def outOfFrameWarning(self, starIndex):
        self.warnings.append(('out_of_frame', int(starIndex), np.nan, np.nan))

    def lowFluxAmplitudeWarning(self, starIndex, xc, yc):
        self.warnings.append(('low_flux', int(starIndex), float(xc), float(yc)))

    def alignmentError(self):
        self.warnings.append(('alignment_error', -1, np.nan, np.nan))


_PLATE_STATUS_SWAP_LOCK = threading.RLock()
_ALIGNMENT_POOL_CONTEXT = {}


def _alignment_pool_initializer(reference_file, generalDark, generalBias, generalFlat,
                                demosaic_fmt, demosaic_out, demosaic_mult, bad_pixel_reference):
    global _ALIGNMENT_POOL_CONTEXT, _TRANSFORM_REFERENCE_IMAGE, _TRANSFORM_REFERENCE_CACHE
    suppress_inherited_tk_cleanup_in_worker()
    _ALIGNMENT_POOL_CONTEXT = {
        'generalDark': generalDark,
        'generalBias': generalBias,
        'generalFlat': generalFlat,
        'demosaic_fmt': demosaic_fmt,
        'demosaic_out': demosaic_out,
        'demosaic_mult': demosaic_mult,
        'bad_pixel_reference': bad_pixel_reference,
    }
    _TRANSFORM_REFERENCE_IMAGE = load_calibrated_reduction_image(
        reference_file,
        generalDark,
        generalBias,
        generalFlat,
        demosaic_fmt,
        demosaic_out,
        demosaic_mult,
        bad_pixel_reference=bad_pixel_reference,
    )
    _TRANSFORM_REFERENCE_CACHE = None


def _load_alignment_worker_frame(file_name):
    context = _ALIGNMENT_POOL_CONTEXT
    hdul = fits.open(name=file_name, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
    extension = 0
    image_header = hdul[extension].header
    while image_header["NAXIS"] == 0:
        extension += 1
        image_header = hdul[extension].header

    image_data = hdul[extension].data
    hdul.close()

    image_data = apply_cals(
        image_data,
        context.get('generalDark'),
        context.get('generalBias'),
        context.get('generalFlat'),
        1,
    )
    image_data = demosaic_img(
        image_data,
        context.get('demosaic_fmt'),
        context.get('demosaic_out'),
        context.get('demosaic_mult'),
        1,
    )
    image_data = repair_bad_pixels_in_frame(image_data, context.get('bad_pixel_reference'))
    return image_header, image_data


def _pointing_precheck_alignment_task(task):
    i, file_name, reference_anchor = task
    try:
        _, image_data = _load_alignment_worker_frame(file_name)
        if getattr(image_data, "ndim", 0) != 2:
            return {
                'index': i,
                'file_name': file_name,
                'usable': False,
                'position': np.array([np.nan, np.nan], dtype=float),
                'transform': None,
            }

        tform = transformation(
            image_data,
            file_name,
            report_failure=False,
            reference_image=_TRANSFORM_REFERENCE_IMAGE,
        )
        mapped_anchor = np.asarray(tform(reference_anchor), dtype=float).reshape(-1, 2)[0]
        usable = bool(np.all(np.isfinite(mapped_anchor)))
        return {
            'index': i,
            'file_name': file_name,
            'usable': usable,
            'position': mapped_anchor,
            'transform': tform if usable else None,
        }
    except Exception as exc:
        return {
            'index': i,
            'file_name': file_name,
            'usable': False,
            'position': np.array([np.nan, np.nan], dtype=float),
            'transform': None,
            'error': str(exc),
        }


def build_multiprocess_pointing_precheck_transforms(inputfiles, max_processes, reference_anchor,
                                                   generalDark=None, generalBias=None, generalFlat=None,
                                                   demosaic_fmt=None, demosaic_out=None, demosaic_mult=None):
    total_jobs = len(inputfiles)
    positions = np.full((total_jobs, 2), np.nan, dtype=float)
    usable_mask = np.zeros(total_jobs, dtype=bool)
    alignment_transforms = {}
    if total_jobs == 0:
        return positions, usable_mask, alignment_transforms

    reference_anchor = np.asarray(reference_anchor, dtype=float).reshape(-1, 2)
    positions[0] = reference_anchor[0]
    usable_mask[0] = True
    alignment_transforms[str(inputfiles[0])] = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])
    if total_jobs == 1:
        return positions, usable_mask, alignment_transforms

    max_workers = min(max_processes, os.cpu_count() or 1, total_jobs - 1, MAX_MULTIPROCESS_TRANSFORM_WORKERS)
    log_info(
        "Using multiprocessing for pointing precheck alignment "
        f"with {max_workers} worker(s) across {total_jobs} image(s)."
    )

    tasks = [
        (i, str(file_name), reference_anchor)
        for i, file_name in enumerate(inputfiles)
        if i != 0
    ]

    with suppress_tk_cleanup_during_process_pool():
        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_alignment_pool_initializer,
            initargs=(
                str(inputfiles[0]),
                generalDark,
                generalBias,
                generalFlat,
                demosaic_fmt,
                demosaic_out,
                demosaic_mult,
                None,
            ),
        ) as executor:
            futures = [executor.submit(_pointing_precheck_alignment_task, task) for task in tasks]
            completed = 1
            for future in as_completed(futures):
                result = future.result()
                index = result['index']
                if result.get('usable'):
                    positions[index] = result['position']
                    usable_mask[index] = True
                    alignment_transforms[result['file_name']] = result['transform']
                completed += 1
                if completed == total_jobs or completed % 10 == 0:
                    log_info(f"Pointing precheck alignment progress: {completed}/{total_jobs}")

    return positions, usable_mask, alignment_transforms


def _fit_alignment_candidate_psfs(image_data, predicted_coords, target_fast_centroid, frame_fast_centroid):
    global plateStatus
    predicted_coords = np.asarray(predicted_coords, dtype=float)
    with _PLATE_STATUS_SWAP_LOCK:
        original_plate_status = plateStatus
        recorder = _ParallelPlateStatusRecorder()
        plateStatus = recorder
        try:
            psf_rows = {
                'target': fit_centroid_or_warn_out_of_frame(
                    image_data,
                    choose_centroid_seed_position(predicted_coords[0], None),
                    0,
                    fast_mode=target_fast_centroid,
                )
            }
            for comp_idx in range(max(0, predicted_coords.shape[0] - 1)):
                psf_rows[f"comp{comp_idx + 1}"] = fit_centroid_or_warn_out_of_frame(
                    image_data,
                    choose_centroid_seed_position(predicted_coords[comp_idx + 1], None),
                    comp_idx + 1,
                    fast_mode=frame_fast_centroid,
                )
            return {
                'coords': predicted_coords,
                'psf_rows': psf_rows,
                'warnings': list(recorder.warnings),
            }
        finally:
            plateStatus = original_plate_status


def _parallel_alignment_task(task):
    (
        i,
        file_name,
        target_and_comp_pixels,
        target_and_comp_radec,
        ignore_header_wcs,
        target_fast_centroid,
        frame_fast_centroid,
        compute_fallback_transform,
        first_frame_uses_input_comp_pixels,
        precomputed_fallback_transform,
    ) = task

    target_and_comp_pixels = np.asarray(target_and_comp_pixels, dtype=float)
    if target_and_comp_radec is not None:
        target_and_comp_radec = np.asarray(target_and_comp_radec, dtype=float)

    image_header, image_data = _load_alignment_worker_frame(file_name)
    result = {
        'index': i,
        'file_name': file_name,
        'wcs': None,
        'fallback': None,
    }

    if not ignore_header_wcs and target_and_comp_radec is not None:
        try:
            wcs_hdr = search_wcs_from_header(image_header)
            if wcs_hdr.is_celestial:
                pix_x, pix_y = wcs_hdr.world_to_pixel_values(
                    target_and_comp_radec[:, 0],
                    target_and_comp_radec[:, 1],
                )
                pix_x = np.asarray(pix_x, dtype=float).reshape(-1)
                pix_y = np.asarray(pix_y, dtype=float).reshape(-1)
                projected_coords = np.column_stack((pix_x, pix_y))
                if i == 0:
                    projected_coords[0] = target_and_comp_pixels[0]
                    if first_frame_uses_input_comp_pixels:
                        projected_coords = np.array(target_and_comp_pixels, dtype=float, copy=True)

                wcs_candidate = _fit_alignment_candidate_psfs(
                    image_data,
                    projected_coords,
                    target_fast_centroid,
                    frame_fast_centroid,
                )
                wcs_candidate['projected_off_frame'] = any_projected_coord_out_of_frame(
                    projected_coords,
                    image_data.shape,
                )
                result['wcs'] = wcs_candidate
        except Exception as exc:
            result['wcs_error'] = str(exc)

    if precomputed_fallback_transform is not None or compute_fallback_transform or result['wcs'] is None:
        if precomputed_fallback_transform is not None:
            tform = precomputed_fallback_transform
        elif i == 0:
            tform = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])
        else:
            tform = transformation(
                image_data,
                file_name,
                report_failure=False,
                reference_image=_TRANSFORM_REFERENCE_IMAGE,
            )
        transformed_coords = np.asarray(tform(target_and_comp_pixels), dtype=float)
        result['fallback'] = _fit_alignment_candidate_psfs(
            image_data,
            transformed_coords,
            target_fast_centroid,
            frame_fast_centroid,
        )

    return result


def _replay_parallel_alignment_warnings(file_name, warnings):
    if not warnings:
        return

    plateStatus.setCurrentFilename(file_name)
    for warning_type, star_index, xc, yc in warnings:
        if warning_type == 'out_of_frame':
            plateStatus.outOfFrameWarning(star_index)
        elif warning_type == 'low_flux':
            plateStatus.lowFluxAmplitudeWarning(star_index, xc, yc)
        elif warning_type == 'alignment_error':
            plateStatus.alignmentError()


def _store_alignment_candidate_psfs(candidate, frame_index, psf_data, comp_keys):
    psf_data['target'][frame_index] = candidate['psf_rows']['target']
    for comp_idx, comp_key in enumerate(comp_keys):
        psf_data[comp_key][frame_index] = candidate['psf_rows'].get(
            f"comp{comp_idx + 1}",
            _nan_psf_result(),
        )


def _update_reference_comp_offsets(psf_data, tar_comp_dist, comp_keys):
    target_row = psf_data['target'][0]
    if not centroid_position_is_finite(target_row):
        return

    for comp_key in comp_keys:
        comp_row = psf_data[comp_key][0]
        if not centroid_position_is_finite(comp_row):
            continue
        tar_comp_dist[comp_key][0] = abs(int(comp_row[0]) - int(target_row[0]))
        tar_comp_dist[comp_key][1] = abs(int(comp_row[1]) - int(target_row[1]))


def apply_parallel_alignment_result(result, frame_index, psf_data, tar_comp_dist, comp_keys):
    wcs_candidate = result.get('wcs')
    selected_candidate = None
    selected_source = 'fallback'

    if wcs_candidate is not None:
        comp_psf_rows = {
            comp_key: wcs_candidate['psf_rows'].get(f"comp{comp_idx + 1}", _nan_psf_result())
            for comp_idx, comp_key in enumerate(comp_keys)
        }
        previous_comp_psf_rows = {}
        if frame_index != 0:
            previous_comp_psf_rows = {comp_key: psf_data[comp_key][frame_index - 1] for comp_key in comp_keys}

        wcs_alignment_decision = should_keep_header_wcs_alignment(
            wcs_candidate.get('projected_off_frame', False),
            frame_index,
            wcs_candidate['psf_rows']['target'],
            previous_target_psf_row=None if frame_index == 0 else psf_data['target'][frame_index - 1],
            comp_psf_rows=comp_psf_rows,
            previous_comp_psf_rows=previous_comp_psf_rows,
            expected_offsets=tar_comp_dist,
        )
        if wcs_alignment_decision['use_wcs_alignment']:
            selected_candidate = wcs_candidate
            selected_source = 'wcs'

    if selected_candidate is None:
        selected_candidate = result.get('fallback') or wcs_candidate

    if selected_candidate is None:
        selected_candidate = {
            'psf_rows': {'target': _nan_psf_result()},
            'warnings': [('alignment_error', -1, np.nan, np.nan)],
        }

    _store_alignment_candidate_psfs(selected_candidate, frame_index, psf_data, comp_keys)
    _replay_parallel_alignment_warnings(result.get('file_name'), selected_candidate.get('warnings'))
    if frame_index == 0:
        _update_reference_comp_offsets(psf_data, tar_comp_dist, comp_keys)

    return selected_source


def build_multiprocess_alignment_results(inputfiles, max_processes, target_and_comp_pixels,
                                         target_and_comp_radec=None, ignore_header_wcs=False,
                                         generalDark=None, generalBias=None, generalFlat=None,
                                         demosaic_fmt=None, demosaic_out=None, demosaic_mult=None,
                                         bad_pixel_reference=None, use_fast_centroid_cadence=False,
                                         use_adaptive_apertures=False, compute_fallback_transform=True,
                                         first_frame_uses_input_comp_pixels=False,
                                         precomputed_fallback_transforms=None):
    total_jobs = len(inputfiles)
    if total_jobs == 0:
        return []

    max_workers = min(max_processes, os.cpu_count() or 1, total_jobs, MAX_MULTIPROCESS_TRANSFORM_WORKERS)
    results = [None] * total_jobs

    log_info(
        "Using multiprocessing for alignment "
        f"with {max_workers} worker(s) across {total_jobs} image(s)."
    )

    tasks = []
    for i, file_name in enumerate(inputfiles):
        precomputed_fallback_transform = None
        if precomputed_fallback_transforms:
            precomputed_fallback_transform = precomputed_fallback_transforms.get(str(file_name))
        frame_fast_centroid = should_use_fast_centroid(i) if use_fast_centroid_cadence else False
        target_fast_centroid = (
            should_use_fast_target_centroid(i, adaptive_apertures=use_adaptive_apertures)
            if use_fast_centroid_cadence else False
        )
        tasks.append((
            i,
            str(file_name),
            target_and_comp_pixels,
            target_and_comp_radec,
            ignore_header_wcs,
            target_fast_centroid,
            frame_fast_centroid,
            compute_fallback_transform,
            first_frame_uses_input_comp_pixels,
            precomputed_fallback_transform,
        ))

    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=_alignment_pool_initializer,
        initargs=(
            str(inputfiles[0]),
            generalDark,
            generalBias,
            generalFlat,
            demosaic_fmt,
            demosaic_out,
            demosaic_mult,
            bad_pixel_reference,
        ),
    ) as executor:
        futures = [executor.submit(_parallel_alignment_task, task) for task in tasks]
        completed = 0
        for future in as_completed(futures):
            result = future.result()
            results[result['index']] = result
            completed += 1
            if completed == total_jobs or completed % 10 == 0:
                log_info(f"Multiprocessing alignment progress: {completed}/{total_jobs}")

    return results


MAX_MULTIPROCESS_TRANSFORM_WORKERS = 8
SPARSE_MISSING_WCS_DROP_THRESHOLD = 0.03
POINTING_REJECTION_MIN_FRAMES = 5
POINTING_REJECTION_MAX_ITERS = 5

# Automatic aperture-grid tuning constants (in PSF sigma units)
APERTURE_SIGMA_MIN = 1.5
APERTURE_SIGMA_MAX = 6.0
ANNULUS_SIGMA_MIN = 6.0
ANNULUS_SIGMA_MAX = 15.0
GAUSSIAN_SIGMA_TO_FWHM = 2.355
SKY_ANNULUS_MIN_GAP_PIXELS = 2.0
SKY_ANNULUS_MIN_FWHM_MULTIPLIER = 2.0
SKY_ANNULUS_MIN_EFFECTIVE_PIXELS = 250.0
SKY_BACKGROUND_SIGMA_CLIP = 3.0
SKY_BACKGROUND_SIGMA_CLIP_MAX_ITERS = 3
APERTURE_AUTOTUNE_COARSE_APER_POINTS = 5
APERTURE_AUTOTUNE_COARSE_ANNULUS_POINTS = 4
APERTURE_AUTOTUNE_REFINED_APER_POINTS = 6
APERTURE_AUTOTUNE_REFINED_ANNULUS_POINTS = 6
APERTURE_AUTOTUNE_APER_HALF_WIDTH_SIGMA = 0.9
APERTURE_AUTOTUNE_ANNULUS_HALF_WIDTH_SIGMA = 2.0
APERTURE_AUTOTUNE_MIN_FRAMES = 8
APERTURE_AUTOTUNE_MAX_FRAMES = 12

# Refit full PSF moments periodically; use a faster moment estimator for most frames.
CENTROID_FULL_FIT_CADENCE = 6


def build_multiprocess_transformations(inputfiles, max_processes):
    reference_file = str(inputfiles[0])
    max_workers = min(max_processes, os.cpu_count() or 1, len(inputfiles), MAX_MULTIPROCESS_TRANSFORM_WORKERS)
    transforms = {}
    total_jobs = len(inputfiles)

    log_info(
        "Using multiprocessing for transformations "
        f"with {max_workers} worker(s) across {total_jobs} image(s)."
    )

    transforms[0] = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])

    with suppress_tk_cleanup_during_process_pool():
        with ProcessPoolExecutor(max_workers=max_workers, initializer=_transformation_pool_initializer,
                                 initargs=(reference_file,)) as executor:
            futures = [executor.submit(transformation_task_with_cached_reference, i, str(file_name))
                       for i, file_name in enumerate(inputfiles) if i != 0]

            completed = 1
            for future in as_completed(futures):
                i, tform = future.result()
                transforms[i] = tform
                completed += 1

                if completed == total_jobs or completed % 10 == 0:
                    log_info(f"Multiprocessing transformations progress: {completed}/{total_jobs}")

    return transforms


def log_alignment_progress(i, total_jobs, file_name, use_multiprocess_progress):
    if use_multiprocess_progress:
        completed = i + 1
        if completed == total_jobs or completed % 10 == 0:
            log_info(f"Multiprocessing alignment progress: {completed}/{total_jobs}")
        return

    display_file_name = _display_filename(file_name)
    sys.stdout.write(f"Aligning frame {i + 1} of {total_jobs} : {display_file_name}\n")
    log.debug(f"Aligning frame {i + 1} of {total_jobs} : {display_file_name}\n")
    sys.stdout.flush()


def get_img_scale(hdr, wcs_file, pixel_init):
    if wcs_file:
        wcs_hdr = get_first_image_header(wcs_file)
        astrometry_scale = [key.value.split(' ') for key in wcs_hdr._cards if 'scale:' in str(key.value)]

        if astrometry_scale:
            img_scale_num = astrometry_scale[0][1]
            img_scale_units = astrometry_scale[0][2]
        else:
            wcs = WCS(wcs_hdr).proj_plane_pixel_scales()
            img_scale_num = (wcs[0].value + wcs[1].value) / 2 * 3600  # Convert to arcsec/pixel
            img_scale_units = "arcsec/pixel"
    elif 'IM_SCALE' in hdr:
        img_scale_num = hdr['IM_SCALE']
        img_scale_units = hdr.comments['IM_SCALE']
    elif 'PIXSCALE' in hdr:
        img_scale_num = hdr['PIXSCALE']
        img_scale_units = hdr.comments['PIXSCALE']
    elif pixel_init:
        img_scale_num = pixel_init
        img_scale_units = "arcsec/pixel"
    else:
        log_info("Not able to find Image Scale in the Image Header.")
        img_scale_num = user_input("Please enter Image Scale (arcsec/pixel): ", type_=float)
        img_scale_units = "arcsec/pixel"

    img_scale = f"Image scale in {img_scale_units}: {round_to_2(float(img_scale_num))}"

    return img_scale, float(img_scale_num)


def exp_time_med(exptimes):
    # exposure time
    consistent_et = False
    if len(exptimes) > 0:
        consistent_et = all(elem == exptimes[0] for elem in exptimes)

    exptimes = np.array(exptimes)

    if consistent_et:
        return exptimes[0]
    else:
        return np.median(exptimes)


def update_coordinates_with_proper_motion(info_dict, time_obs):
    parameter_names = {
        'dist': 'Distance (pc)',
        'pm_ra': 'Proper Motion RA (mas/yr)',
        'pm_dec': 'Proper Motion DEC (mas/yr)'
    }

    numeric_values = {}
    missing_values = []

    for key in ['dist', 'pm_ra', 'pm_dec']:
        raw_value = info_dict.get(key, 0.0)

        try:
            parsed_value = float(raw_value)
        except (TypeError, ValueError):
            parsed_value = 0.0

        numeric_values[key] = parsed_value
        if parsed_value == 0.0:
            missing_values.append(parameter_names[key])

    if missing_values:
        missing_values = ", ".join(missing_values)
        log_info("Warning: Cannot account for proper motion due to missing values in: "
                 f"\n{missing_values}. If you find your target or comparisons are not detected well, please "
                 f"re-run and fill in values in the initialization file to account for proper motion", warn=True)

        return info_dict['ra'], info_dict['dec']
    else:
        time_j2000 = Time(2000.0, format='jyear')
        time_obs = Time(time_obs, format='jd')

        coord = SkyCoord(
            ra=info_dict['ra'] * u.deg,
            dec=info_dict['dec'] * u.deg,
            distance=numeric_values['dist'] * u.pc,
            pm_ra_cosdec=numeric_values['pm_ra'] * u.mas / u.yr,
            pm_dec=numeric_values['pm_dec'] * u.mas / u.yr,
            frame="icrs",
            obstime=time_j2000
        )

        updated_coord = coord.apply_space_motion(new_obstime=time_obs)
        return updated_coord.ra.deg, updated_coord.dec.deg


def gaussian_psf(x, y, x0, y0, a, sigx, sigy, rot, b):
    rx = (x - x0) * np.cos(rot) - (y - y0) * np.sin(rot)
    ry = (x - x0) * np.sin(rot) + (y - y0) * np.cos(rot)
    gausx = np.exp(-rx ** 2 / (2 * sigx ** 2))
    gausy = np.exp(-ry ** 2 / (2 * sigy ** 2))
    return a * gausx * gausy + b


def mesh_box(pos, box, maxx=0, maxy=0):
    pos = [int(np.round(pos[0])), int(np.round(pos[1]))]
    if maxx:
        x = np.arange(max(0,pos[0] - box), min(maxx, pos[0] + box + 1))
    else:
        x = np.arange(max(0,pos[0] - box), pos[0] + box + 1)
    if maxy:
        y = np.arange(max(0,pos[1] - box), min(maxy, pos[1] + box + 1))
    else:
        y = np.arange(max(0,pos[1] - box), pos[1] + box + 1)
    xv, yv = np.meshgrid(x, y)
    return xv.astype(int), yv.astype(int)


def should_use_fast_centroid(frame_index):
    return frame_index % CENTROID_FULL_FIT_CADENCE != 0


def should_use_fast_target_centroid(frame_index, adaptive_apertures=False):
    return should_use_fast_centroid(frame_index) and not adaptive_apertures


def _fit_centroid_moments(subarray, xv, yv, pos, box):
    background = bn.nanmedian(subarray)
    weights = subarray - background
    weights = np.where(np.isfinite(weights) & (weights > 0), weights, 0.0)
    wsum = np.sum(weights)

    if not np.isfinite(wsum) or wsum <= 0:
        floor = np.nanmin(subarray)
        weights = subarray - floor
        weights = np.where(np.isfinite(weights) & (weights > 0), weights, 0.0)
        wsum = np.sum(weights)

    if not np.isfinite(wsum) or wsum <= 0:
        return np.empty(7) * np.nan

    wx = float(np.sum(xv * weights) / wsum)
    wy = float(np.sum(yv * weights) / wsum)

    # Keep centroid near the expected star location in crowded fields.
    wx = float(np.clip(wx, pos[0] - box * 0.5, pos[0] + box * 0.5))
    wy = float(np.clip(wy, pos[1] - box * 0.5, pos[1] + box * 0.5))

    dx = xv - wx
    dy = yv - wy
    var_x = float(np.sum(weights * dx * dx) / wsum)
    var_y = float(np.sum(weights * dy * dy) / wsum)
    cov_xy = float(np.sum(weights * dx * dy) / wsum)

    sigx = float(np.clip(np.sqrt(max(var_x, 0.25)), 0.5, 20.0))
    sigy = float(np.clip(np.sqrt(max(var_y, 0.25)), 0.5, 20.0))
    rot = float(0.5 * np.arctan2(2.0 * cov_xy, var_x - var_y)) if np.isfinite(cov_xy) else 0.0
    amp = float(max(np.nanmax(subarray) - background, 0.0))

    return np.array([wx, wy, amp, sigx, sigy, rot, float(background)], dtype=float)


def _has_usable_centroid_signal(subarray, amplitude, min_snr=5.0):
    if not np.isfinite(amplitude) or amplitude <= 0:
        return False

    scatter = float(bn.nanstd(subarray))
    if not np.isfinite(scatter) or scatter <= 0:
        return True

    return amplitude >= (min_snr * scatter)


def _nan_psf_result():
    return np.full(7, np.nan, dtype=float)


def fit_centroid_or_warn_out_of_frame(data, pos, starIndex, **kwargs):
    if not pixel_within_image(pos[0], pos[1], data.shape):
        plateStatus.outOfFrameWarning(starIndex)
        return _nan_psf_result()
    return fit_centroid(data, pos, starIndex, **kwargs)


def fractional_flux_change_within_limit(current_amplitude, previous_amplitude, limit=0.5):
    if (not np.isfinite(current_amplitude)
            or not np.isfinite(previous_amplitude)
            or previous_amplitude == 0):
        return False

    return np.abs((current_amplitude - previous_amplitude) / previous_amplitude) <= limit


def centroid_position_is_finite(psf_row):
    try:
        coords = np.asarray(psf_row[:2], dtype=float)
    except (TypeError, ValueError, IndexError):
        return False

    return bool(np.all(np.isfinite(coords)))


def choose_centroid_seed_position(predicted_pos, previous_psf_row=None, max_offset_pixels=5.0):
    predicted = np.asarray(predicted_pos, dtype=float).reshape(-1)
    if predicted.size < 2 or not np.all(np.isfinite(predicted[:2])):
        return np.array([np.nan, np.nan], dtype=float)

    if not centroid_position_is_finite(previous_psf_row):
        return np.array(predicted[:2], dtype=float)

    previous = np.asarray(previous_psf_row[:2], dtype=float)
    if np.hypot(*(previous - predicted[:2])) > float(max_offset_pixels):
        return np.array(predicted[:2], dtype=float)

    return previous.astype(float, copy=True)


def centroid_offset_matches_reference(psf_a, psf_b, expected_dx, expected_dy,
                                      tolerance=WCS_REFERENCE_GEOMETRY_TOLERANCE_PIXELS):
    if not centroid_position_is_finite(psf_a) or not centroid_position_is_finite(psf_b):
        return False
    if not np.isfinite(expected_dx) or not np.isfinite(expected_dy):
        return False

    dx = float(abs(float(psf_a[0]) - float(psf_b[0])))
    dy = float(abs(float(psf_a[1]) - float(psf_b[1])))
    tolerance = float(tolerance)

    return (
        abs(dx - float(expected_dx)) <= tolerance
        and abs(dy - float(expected_dy)) <= tolerance
    )


def should_keep_header_wcs_alignment(
    projected_off_frame,
    frame_index,
    target_psf_row,
    previous_target_psf_row=None,
    comp_psf_rows=None,
    previous_comp_psf_rows=None,
    expected_offsets=None,
    tolerance=WCS_REFERENCE_GEOMETRY_TOLERANCE_PIXELS,
    min_geometry_match_fraction=WCS_MIN_GEOMETRY_MATCH_FRACTION,
):
    decision = {
        'use_wcs_alignment': False,
        'reason': 'missing_target_centroid',
        'target_flux_change_ok': True,
        'comp_flux_change_ok': True,
        'geometry_match_count': 0,
        'geometry_test_count': 0,
    }

    if projected_off_frame:
        decision.update(use_wcs_alignment=True, reason='projected_off_frame')
        return decision

    if not centroid_position_is_finite(target_psf_row):
        return decision

    if frame_index == 0:
        decision.update(use_wcs_alignment=True, reason='first_frame')
        return decision

    if previous_target_psf_row is not None:
        decision['target_flux_change_ok'] = fractional_flux_change_within_limit(
            target_psf_row[2],
            previous_target_psf_row[2],
        )

    comp_psf_rows = {} if comp_psf_rows is None else dict(comp_psf_rows)
    previous_comp_psf_rows = {} if previous_comp_psf_rows is None else dict(previous_comp_psf_rows)
    expected_offsets = {} if expected_offsets is None else dict(expected_offsets)

    for key, comp_row in comp_psf_rows.items():
        prev_comp_row = previous_comp_psf_rows.get(key)
        if prev_comp_row is not None:
            decision['comp_flux_change_ok'] = (
                decision['comp_flux_change_ok']
                and fractional_flux_change_within_limit(comp_row[2], prev_comp_row[2])
            )

        if not centroid_position_is_finite(comp_row):
            continue

        expected_offset = expected_offsets.get(key)
        if expected_offset is None:
            continue

        try:
            expected_dx = float(expected_offset[0])
            expected_dy = float(expected_offset[1])
        except (TypeError, ValueError, IndexError):
            continue

        decision['geometry_test_count'] += 1
        if centroid_offset_matches_reference(
            comp_row,
            target_psf_row,
            expected_dx,
            expected_dy,
            tolerance=tolerance,
        ):
            decision['geometry_match_count'] += 1

    geometry_test_count = decision['geometry_test_count']
    if geometry_test_count == 0:
        decision.update(use_wcs_alignment=True, reason='finite_target_only')
        return decision

    minimum_matches = max(1, int(np.ceil(float(min_geometry_match_fraction) * geometry_test_count)))
    if decision['geometry_match_count'] >= minimum_matches:
        decision.update(use_wcs_alignment=True, reason='geometry_match')
    else:
        decision['reason'] = 'geometry_mismatch'

    return decision


# Method fits a 2D gaussian function that matches the star_psf to the star image and returns its pixel coordinates
def fit_centroid(data, pos, starIndex, psf_function=gaussian_psf, box=15, weightedcenter=False, fast_mode=False):
    stage_start = perf_counter()
    # get sub field in image
    try:
        xv, yv = mesh_box(pos, box, maxx=data.shape[1], maxy=data.shape[0])
        subarray = data[yv, xv]
        try:
            init = [np.nanmax(subarray) - np.nanmin(subarray), 1, 1, 0, np.nanmin(subarray)]
        except ValueError as ve:
            # Handle null subfield - cannot solve
            plateStatus.outOfFrameWarning(starIndex)
            log.debug(f"Warning: empty subfield for fit_centroid at {np.round(pos, 2)}")
            return _nan_psf_result()

        moment_fit = _fit_centroid_moments(subarray, xv, yv, pos, box)
        if np.isfinite(moment_fit[0]):
            wx, wy = moment_fit[0], moment_fit[1]
            init = [moment_fit[2], moment_fit[3], moment_fit[4], moment_fit[5], moment_fit[6]]
            if fast_mode:
                if _has_usable_centroid_signal(subarray, init[0]):
                    return moment_fit

                plateStatus.lowFluxAmplitudeWarning(starIndex, pos[0], pos[1])
                log.debug(
                    f"Warning: Measured fast centroid amplitude is really low---"
                    f"are you sure there is a star at {np.round(pos, 2)}?"
                )
                return _nan_psf_result()
        else:
            # compute flux weighted centroid in x and y
            wx = np.sum(xv[0] * subarray.sum(0)) / subarray.sum(0).sum()
            wy = np.sum(yv[:, 0] * subarray.sum(1)) / subarray.sum(1).sum()

        # lower bound: [xc, yc, amp, sigx, sigy, rotation,  bg]
        lo = [pos[0] - box * 0.5, pos[1] - box * 0.5, 0, 0.5, 0.5, -np.pi / 4, np.nanmin(subarray) - 1]
        up = [pos[0] + box * 0.5, pos[1] + box * 0.5, 1e7, 20, 20, np.pi / 4, np.nanmax(subarray) + 1]
        x0 = np.array([*pos, *init], dtype=float)
        lo_arr = np.array(lo, dtype=float)
        up_arr = np.array(up, dtype=float)
        if np.all(np.isfinite(x0)):
            x0 = np.clip(x0, lo_arr + 1e-6, up_arr - 1e-6)
        has_usable_signal = _has_usable_centroid_signal(subarray, init[0])

        def fcn2min(pars):
            model = psf_function(xv, yv, *pars)
            return (subarray - model).flatten()

        try:
            res = least_squares(fcn2min, x0=x0, bounds=[lo, up], jac='2-point', xtol=None, method='trf')
        except Exception as exc:
            if has_usable_signal and np.isfinite(moment_fit[0]):
                log.debug(f"Centroid PSF fit failed at {np.round(pos, 2)}; using moment centroid instead: {exc}")
                return moment_fit

            if not has_usable_signal:
                plateStatus.lowFluxAmplitudeWarning(starIndex, pos[0], pos[1])
                log.debug(f"Warning: Measured flux amplitude is really low---are you sure there is a star at {np.round(pos, 2)}?")
            else:
                log.debug(f"Centroid PSF fit failed at {np.round(pos, 2)}; attempting LM fallback: {exc}")

            try:
                res = least_squares(fcn2min, x0=x0, jac='2-point', xtol=1e-12, method='lm')
            except Exception as lm_exc:
                log.debug(f"Centroid LM fallback failed at {np.round(pos, 2)}: {lm_exc}")
                return _nan_psf_result()

        # Preserve the solved PSF center for subpixel tracking by default.
        # The weighted-center override remains available as an explicit legacy option.
        if weightedcenter:
            res.x[0] = wx
            res.x[1] = wy
        if np.isfinite(moment_fit[6]):
            res.x[6] = moment_fit[6]

        return res.x
    finally:
        _record_photometry_stage_timing('fit_centroid', perf_counter() - stage_start)


def sigma_clipped_nanmedian(data, sigma=3.0, max_iters=3):
    clipped = np.array(data, dtype=float, copy=True)
    if clipped.size == 0:
        return np.nan, np.nan

    clipped[~np.isfinite(clipped)] = np.nan
    nan_count_prev = np.count_nonzero(np.isnan(clipped))

    for _ in range(max_iters):
        center = bn.nanmedian(clipped)
        scatter = bn.nanstd(clipped)
        if not np.isfinite(center):
            return np.nan, np.nan
        if not np.isfinite(scatter) or scatter <= 0:
            break

        clipped[np.abs(clipped - center) > sigma * scatter] = np.nan
        nan_count = np.count_nonzero(np.isnan(clipped))
        if nan_count == nan_count_prev:
            break
        nan_count_prev = nan_count

    return bn.nanmedian(clipped), bn.nanstd(clipped)


def normalize_flux_series_to_approximate_unity(
    flux_values,
    flux_errors=None,
    sigma=3.0,
    max_iters=3,
    min_points=LIGHTCURVE_MIN_VALID_POINTS,
):
    flux_values = np.asarray(flux_values, dtype=float)
    normalized_flux = np.array(flux_values, dtype=float, copy=True)
    normalized_unc = None if flux_errors is None else np.array(flux_errors, dtype=float, copy=True)

    valid_flux = np.isfinite(flux_values) & (flux_values > 0)
    if np.count_nonzero(valid_flux) < max(1, int(min_points)):
        return normalized_flux, normalized_unc, np.nan

    baseline_level, _ = sigma_clipped_nanmedian(
        flux_values[valid_flux],
        sigma=sigma,
        max_iters=max_iters,
    )
    if not np.isfinite(baseline_level) or baseline_level <= 0:
        baseline_level = bn.nanmedian(flux_values[valid_flux])
    if not np.isfinite(baseline_level) or baseline_level <= 0:
        return normalized_flux, normalized_unc, np.nan

    normalized_flux[valid_flux] = flux_values[valid_flux] / baseline_level
    if normalized_unc is not None and normalized_unc.shape == flux_values.shape:
        finite_unc = np.isfinite(normalized_unc)
        normalized_unc[finite_unc] = normalized_unc[finite_unc] / baseline_level

    return normalized_flux, normalized_unc, float(baseline_level)


def weighted_nanpercentile(values, weights, percentile):
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(valid):
        return np.nan

    values = values[valid]
    weights = weights[valid]
    sort_index = np.argsort(values, kind='mergesort')
    values = values[sort_index]
    weights = weights[sort_index]

    total_weight = float(np.sum(weights))
    if not np.isfinite(total_weight) or total_weight <= 0:
        return np.nan

    if values.size == 1:
        return float(values[0])

    cumulative = (np.cumsum(weights) - 0.5 * weights) / total_weight
    target = float(np.clip(percentile, 0.0, 100.0)) / 100.0
    return float(np.interp(target, cumulative, values, left=values[0], right=values[-1]))


def weighted_nanstd(values, weights):
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(valid):
        return np.nan

    values = values[valid]
    weights = weights[valid]
    total_weight = float(np.sum(weights))
    if not np.isfinite(total_weight) or total_weight <= 0:
        return np.nan

    mean = float(np.sum(weights * values) / total_weight)
    variance = float(np.sum(weights * (values - mean) ** 2) / total_weight)
    return float(np.sqrt(max(variance, 0.0)))


def sigma_clipped_weighted_median(values, weights, sigma=3.0, max_iters=3, high_only=False):
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(keep):
        return np.nan, np.nan

    for _ in range(max_iters):
        center = weighted_nanpercentile(values[keep], weights[keep], 50.0)
        scatter = weighted_nanstd(values[keep], weights[keep])
        if not np.isfinite(center):
            return np.nan, np.nan
        if not np.isfinite(scatter) or scatter <= 0:
            break

        if high_only:
            updated_keep = keep & (values <= center + sigma * scatter)
        else:
            updated_keep = keep & (np.abs(values - center) <= sigma * scatter)
        if np.array_equal(updated_keep, keep):
            break
        keep = updated_keep

        if not np.any(keep):
            return np.nan, np.nan

    return weighted_nanpercentile(values[keep], weights[keep], 50.0), weighted_nanstd(values[keep], weights[keep])


def psf_fwhm_from_sigma(sigma):
    try:
        sigma = float(sigma)
    except (TypeError, ValueError):
        return np.nan

    if not np.isfinite(sigma) or sigma <= 0:
        return np.nan

    return float(GAUSSIAN_SIGMA_TO_FWHM * sigma)


def resolve_sky_annulus_geometry(
    aperture_radius,
    annulus_width,
    psf_sigma=np.nan,
    minimum_gap_pixels=SKY_ANNULUS_MIN_GAP_PIXELS,
    minimum_fwhm_multiplier=SKY_ANNULUS_MIN_FWHM_MULTIPLIER,
    minimum_sky_pixels=SKY_ANNULUS_MIN_EFFECTIVE_PIXELS,
):
    aperture_radius = abs(float(aperture_radius))
    annulus_width = max(float(annulus_width), 0.0)

    inner_radius = aperture_radius + float(minimum_gap_pixels)
    fwhm = psf_fwhm_from_sigma(psf_sigma)
    if np.isfinite(fwhm):
        inner_radius = max(inner_radius, float(minimum_fwhm_multiplier) * fwhm)

    outer_radius = inner_radius + annulus_width
    effective_sky_pixels = np.pi * max(outer_radius ** 2 - inner_radius ** 2, 0.0)

    if minimum_sky_pixels is not None and np.isfinite(minimum_sky_pixels) and minimum_sky_pixels > 0:
        minimum_outer_radius = float(np.sqrt(inner_radius ** 2 + float(minimum_sky_pixels) / np.pi))
        if minimum_outer_radius > outer_radius:
            outer_radius = minimum_outer_radius
            effective_sky_pixels = np.pi * max(outer_radius ** 2 - inner_radius ** 2, 0.0)

    return {
        'inner_radius': float(inner_radius),
        'outer_radius': float(outer_radius),
        'annulus_width': float(max(outer_radius - inner_radius, 0.0)),
        'effective_sky_pixels': float(effective_sky_pixels),
        'fwhm': float(fwhm) if np.isfinite(fwhm) else np.nan,
    }


# Method calculates the flux of the star (uses the skybg_phot method to do background sub)
def aperPhot(data, starIndex, xc, yc, r=5, dr=5, fast_mode=False, sigma_hint=np.nan):
    stage_start = perf_counter()
    try:
        # Check for invalid coordinates
        if np.isnan(xc) or np.isnan(yc):
            return 0, 0

        # Calculate background if dr > 0
        if dr > 0:
            sky_geometry = resolve_sky_annulus_geometry(r, dr, psf_sigma=sigma_hint)
            bgflux, sigmabg, Nbg = skybg_phot(
                data,
                starIndex,
                xc,
                yc,
                sky_geometry['inner_radius'],
                sky_geometry['annulus_width'],
                fast_mode=fast_mode,
            )
            if not np.isfinite(bgflux):
                return np.nan, bgflux
        else:
            bgflux, sigmabg, Nbg = 0, 0, 0

        # Create aperture and mask
        aperture = CircularAperture(positions=[(xc, yc)], r=r)
        mask_method = 'center' if fast_mode else 'exact'
        mask = aperture.to_mask(method=mask_method)[0]
        data_cutout = mask.cutout(data)

        # Check if aperture is valid
        if data_cutout is None:
            # Aperture is partially or fully outside the image
            return 0, bgflux    # Return zero flux but valid background

        # Calculate and return aperture sum
        aperture_sum = (mask.data * (data_cutout - bgflux)).sum()
        return aperture_sum, bgflux
    finally:
        _record_photometry_stage_timing('aperPhot', perf_counter() - stage_start)


def skybg_phot(data, starIndex, xc, yc, r=10, dr=5, ptol=99, debug=False, fast_mode=False):
    # The sky annulus uses an inner radius r and an outer radius r + dr.
    # Callers are responsible for choosing r and dr from the aperture radius and PSF size.
    annulus = CircularAnnulus(positions=[(xc, yc)], r_in=float(r), r_out=float(r + dr))
    mask_method = 'center' if fast_mode else 'exact'
    annulus_mask = annulus.to_mask(method=mask_method)[0]
    annulus_cutout = annulus_mask.cutout(data, fill_value=np.nan)

    if annulus_cutout is None:
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: empty sky background annulus for {xc:.1f}, {yc:.1f}."
                 f"\nCheck if star is present or close to border.")
        return np.nan, np.nan, 0

    annulus_cutout = np.asarray(annulus_cutout, dtype=float)
    annulus_weights = np.asarray(annulus_mask.data, dtype=float)
    valid_mask = np.isfinite(annulus_cutout) & np.isfinite(annulus_weights) & (annulus_weights > 0)
    if not np.any(valid_mask):
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: no valid sky background pixels for {xc:.1f}, {yc:.1f}."
                 f"\nCheck if star is present or close to border.")
        return np.nan, np.nan, 0

    annulus_pixels = annulus_cutout[valid_mask]
    annulus_pixel_weights = annulus_weights[valid_mask]

    try:
        cutoff = weighted_nanpercentile(annulus_pixels, annulus_pixel_weights, ptol)
    except (IndexError, ValueError):
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: IndexError, problem computing sky bg for {xc:.1f}, {yc:.1f}."
                 f"\nCheck if star is present or close to border.")
        return np.nan, np.nan, 0

    if not np.isfinite(cutoff):
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: invalid cutoff while computing sky bg for {xc:.1f}, {yc:.1f}.")
        return np.nan, np.nan, 0

    clipped_keep = annulus_pixels <= cutoff
    clipped_pixels = annulus_pixels[clipped_keep]
    clipped_weights = annulus_pixel_weights[clipped_keep]
    if clipped_pixels.size == 0:
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: percentile clipping removed all sky background pixels for {xc:.1f}, {yc:.1f}.")
        return np.nan, np.nan, 0

    dat = np.full_like(annulus_cutout, np.nan, dtype=float)
    dat[valid_mask] = annulus_cutout[valid_mask]
    dat[valid_mask & (annulus_cutout > cutoff)] = np.nan

    if debug:
        minb = float(np.nanmin(annulus_pixels))
        maxb = float(np.nanmean(annulus_pixels) + 3 * np.nanstd(annulus_pixels))
        bgsky = np.full_like(annulus_cutout, np.nan, dtype=float)
        bgsky[valid_mask] = annulus_cutout[valid_mask]
        cmed, _ = sigma_clipped_weighted_median(
            clipped_pixels,
            clipped_weights,
            sigma=SKY_BACKGROUND_SIGMA_CLIP,
            max_iters=SKY_BACKGROUND_SIGMA_CLIP_MAX_ITERS,
            high_only=True,
        )
        amed, _ = sigma_clipped_weighted_median(
            annulus_pixels,
            annulus_pixel_weights,
            sigma=SKY_BACKGROUND_SIGMA_CLIP,
            max_iters=SKY_BACKGROUND_SIGMA_CLIP_MAX_ITERS,
            high_only=True,
        )

        fig, ax = plt.subplots(2, 2, figsize=(9, 9))
        im = ax[0, 0].imshow(annulus_cutout, vmin=minb, vmax=maxb, cmap='inferno')
        ax[0, 0].set_title("Original Data")
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax[0, 0])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        ax[1, 0].hist(annulus_pixels, label=f'Sky Annulus ({np.nanmedian(annulus_pixels):.1f}, {amed:.1f})',
                      alpha=0.5, bins=np.arange(minb, maxb))
        ax[1, 0].hist(clipped_pixels, label=f'Clipped ({np.nanmedian(clipped_pixels):.1f}, {cmed:.1f})', alpha=0.5,
                      bins=np.arange(minb, maxb))
        ax[1, 0].legend(loc='best')
        ax[1, 0].set_title("Sky Background")
        ax[1, 0].set_xlabel("Pixel Value")

        ax[1, 1].imshow(dat, vmin=minb, vmax=maxb, cmap='inferno')
        ax[1, 1].set_title("Clipped Sky Background")

        ax[0, 1].imshow(bgsky, vmin=minb, vmax=maxb, cmap='inferno')
        ax[0, 1].set_title("Sky Annulus")
        plt.tight_layout()
        plt.show()
    sky_median, sky_sigma = sigma_clipped_weighted_median(
        clipped_pixels,
        clipped_weights,
        sigma=SKY_BACKGROUND_SIGMA_CLIP,
        max_iters=SKY_BACKGROUND_SIGMA_CLIP_MAX_ITERS,
        high_only=True,
    )
    return sky_median, sky_sigma, float(np.sum(annulus_pixel_weights))

def process_dark_frames(dark_files):
    """Process dark frames and return the master dark."""
    if not dark_files:
        return None
    # Dark files whose median is much higher than the overall dark files median will be filtered
    # e.g. to discard saturated dark files that may negatively affect the master dark used to calibrate the science frames
    # First pass: collect all dark frame medians
    darks_medians = [(dark_file, np.nanmedian(fits.getdata(dark_file))) for dark_file in dark_files]

    d_median = np.median([median for _, median in darks_medians])
    threshold = 1.7  # 70% higher than overall median

    # Second pass: collect valid dark frames
    darks_img_list = []
    for dark_file, dark_median in darks_medians:
        median_ratio = dark_median / d_median
        if median_ratio > threshold:
            log_info(
                f"\nWarning: Skipping suspicious dark frame {dark_file}: "
                f"median/overall_median = {median_ratio:.2f}\n",
                warn=True
            )
            continue
        dark_data = fits.getdata(dark_file)
        darks_img_list.append(dark_data)
            
    return np.median(darks_img_list, axis=0) if darks_img_list else None

def process_bias_frames(bias_files):
    """Process bias frames and return the master bias."""
    if not bias_files:
        return None
        
    biases_img_list = [fits.getdata(bias_file) for bias_file in bias_files]  
    return np.median(biases_img_list, axis=0) if biases_img_list else None

def process_flat_frames(flat_files, master_bias=None):
    """Process flat frames and return the normalized master flat."""
    if not flat_files:
        return None
        
    flats_img_list = [fits.getdata(flat_file) for flat_file in flat_files]      
    master_flat = np.median(flats_img_list, axis=0)
    # Bias subtract after creating master flat
    if master_bias is not None:
        master_flat = master_flat - master_bias
    # Normalize
    medi = np.median(master_flat)
    return master_flat / medi

def convert_jd_to_bjd(non_bjd, p_dict, info_dict):
    global _BJD_FALLBACK_WARNING_LOGGED

    try:
        goodTimes = JDUTC_to_BJDTDB(non_bjd, ra=p_dict['ra'], dec=p_dict['dec'], lat=info_dict['lat'],
                                    longi=info_dict['long'], alt=info_dict['elev'])[0]
    except Exception as exc:
        if not _BJD_FALLBACK_WARNING_LOGGED:
            _BJD_FALLBACK_WARNING_LOGGED = True
            try:
                log.warning(
                    "barycorrpy JDUTC_to_BJDTDB conversion failed; falling back to astropy light-travel-time "
                    "conversion for this run.",
                    exc_info=True,
                )
            except Exception:
                traceback.print_exception(type(exc), exc, exc.__traceback__, file=sys.stdout)
        targetloc = SkyCoord(p_dict['ra'], p_dict['dec'], unit=(u.deg, u.deg), frame='icrs')
        obsloc = EarthLocation(lat=info_dict['lat'], lon=info_dict['long'], height=info_dict['elev'])
        timesToConvert = Time(non_bjd, format='jd', scale='utc', location=obsloc)
        ltt_bary = timesToConvert.light_travel_time(targetloc)
        time_barycentre = timesToConvert.tdb  + ltt_bary
        goodTimes = time_barycentre.value

    return goodTimes


def calculate_variablility(fit_lc_ref, fit_lc_best):
    info_ref = None

    mask_oot_ref = (fit_lc_ref.transit == 1)
    mask_oot_best = (fit_lc_best.transit == 1)

    intx_times = np.intersect1d(fit_lc_best.jd_times[mask_oot_best], fit_lc_ref.jd_times[mask_oot_ref])

    if intx_times.any():
        mask_ref = np.isin(fit_lc_ref.jd_times, intx_times)
        mask_best = np.isin(fit_lc_best.jd_times, intx_times)

        norm_flux_ref = (fit_lc_ref.data / np.nanmedian(fit_lc_ref.data[mask_ref]))[mask_ref]
        norm_flux_best = (fit_lc_best.data / np.nanmedian(fit_lc_best.data[mask_best]))[mask_best]

        info_ref = {
            'fit_lc': fit_lc_ref,
            'mask_ref': mask_ref,
            'res': norm_flux_best - norm_flux_ref,
        }

    return info_ref


def choose_comp_star_variability(fit_lc_refs, fit_lc_best, ref_comp, comp_stars, vsp_comp_stars, save):
    colors = ["firebrick", "darkorange", "olivedrab", "lightseagreen", "steelblue", "rebeccapurple", "mediumvioletred"]
    markers = ['.', 'v', 's', 'D', '^']
    k = 0

    labels = {tuple(value['pos']): key for key, value in vsp_comp_stars.items()}

    for i, ckey in enumerate(fit_lc_refs.keys()):
        if i >= len(colors):
            i = 0
        if k >= len(markers):
            k = 0
        ref_comp[ckey] = calculate_variablility(fit_lc_refs[ckey]['myfit'], fit_lc_best)

        if ref_comp[ckey]:
            plt.errorbar(ref_comp[ckey]['fit_lc'].jd_times[ref_comp[ckey]['mask_ref']], ref_comp[ckey]['res'],
                         fmt=markers[k], color=colors[i], label=f"{labels[tuple(fit_lc_refs[ckey]['pos'])]}")
        k += 1

    plot_variable_residuals(save)

    std_devs = {key: np.std(value['res']) for key, value in ref_comp.items() if value}
    if not std_devs:
        raise RuntimeError("No usable comparison-star residuals were available for stellar variability calibration.")
    min_std_dev = min(std_devs, key=lambda y: abs(std_devs[y]))

    return comp_stars[min_std_dev]


def stellar_variability_label(comp_label, comp_star):
    if comp_star.get('is_aavso_vsp', True):
        return comp_label
    comp_ra = _finite_float(comp_star.get('ra'))
    comp_dec = _finite_float(comp_star.get('dec'))
    if comp_ra is not None and comp_dec is not None:
        return f"RA={comp_ra:.7f} Dec={comp_dec:.7f}"
    return comp_label


def build_stellar_variability_params_from_fit(lc_fit, comp_star, comp_pos, comp_label, save, s_name,
                                              observed_filter=None):
    comp_mag = _finite_float(comp_star.get('mag'))
    comp_mag_error = _finite_float(comp_star.get('error'))
    if comp_mag is None or comp_mag_error is None:
        raise RuntimeError("Comparison-star magnitude or magnitude uncertainty is unavailable.")
    observed_filter = observed_filter or comp_star.get('observed_filter')

    fit_data = np.asarray(getattr(lc_fit, 'data', []), dtype=float)
    fit_airmass_model = np.asarray(
        getattr(lc_fit, 'airmass_model', np.ones_like(fit_data)),
        dtype=float,
    )
    fit_airmass = np.asarray(getattr(lc_fit, 'airmass', np.ones_like(fit_data)), dtype=float)
    fit_times = np.asarray(getattr(lc_fit, 'jd_times', getattr(lc_fit, 'time', [])), dtype=float)
    transit_model = np.asarray(getattr(lc_fit, 'transit', np.ones_like(fit_data)), dtype=float)

    if not (fit_data.shape == fit_airmass_model.shape == fit_airmass.shape == fit_times.shape):
        raise RuntimeError("Lightcurve arrays have inconsistent shapes for stellar variability output.")

    if transit_model.shape == fit_data.shape:
        mask_ref = transit_model == 1
    else:
        mask_ref = np.ones_like(fit_data, dtype=bool)

    with np.errstate(divide='ignore', invalid='ignore'):
        detrended_all = np.divide(fit_data, fit_airmass_model)
    if np.count_nonzero(mask_ref) == 0:
        mask_ref = np.isfinite(detrended_all)

    detrended = detrended_all[mask_ref]
    selected_airmass_model = fit_airmass_model[mask_ref]
    selected_data = fit_data[mask_ref]
    selected_times = fit_times[mask_ref]
    selected_airmass = fit_airmass[mask_ref]

    oot_scatter = np.nanstd(detrended)
    median_data = np.nanmedian(selected_data)
    with np.errstate(divide='ignore', invalid='ignore'):
        norm_flux_unc = oot_scatter * selected_airmass_model / median_data
        target_mag = comp_mag - (2.5 * np.log10(detrended))
        target_mag_error = (
            comp_mag_error ** 2
            + (-2.5 * norm_flux_unc / (detrended * np.log(10))) ** 2
        ) ** 0.5

    valid = (
        np.isfinite(selected_times)
        & np.isfinite(selected_airmass)
        & np.isfinite(target_mag)
        & np.isfinite(target_mag_error)
    )
    if np.count_nonzero(valid) == 0:
        raise RuntimeError("No finite stellar variability magnitude points were produced.")

    display_label = stellar_variability_label(comp_label, comp_star)
    vsp_params = []
    for time_value, airmass_value, mag_value, mag_error_value in zip(
        selected_times[valid],
        selected_airmass[valid],
        target_mag[valid],
        target_mag_error[valid],
    ):
        vsp_params.append({
            'time': time_value,
            'airmass': airmass_value,
            'mag': mag_value,
            'mag_err': mag_error_value,
            'cname': display_label,
            'cmag': comp_mag,
            'cmag_err': comp_mag_error,
            'pos': comp_pos,
            'comp_ra': comp_star.get('ra'),
            'comp_dec': comp_star.get('dec'),
            'catalog_ra': comp_star.get('catalog_ra'),
            'catalog_dec': comp_star.get('catalog_dec'),
            'catalog_source': comp_star.get('catalog_source', 'AAVSO VSP'),
            'is_aavso_vsp': bool(comp_star.get('is_aavso_vsp', True)),
            'mag_band': comp_star.get('mag_band', 'V'),
            'observed_filter': observed_filter,
            'source_id': comp_star.get('source_id'),
            'catalog_id': comp_star.get('id'),
            'separation_arcsec': comp_star.get('separation_arcsec'),
        })

    plot_stellar_variability(vsp_params, save, s_name, display_label)
    return vsp_params


def stellar_variability(fit_lc_refs, fit_lc_best, comp_stars, vsp_comp_stars, vsp_ind, best_comp, save, s_name,
                        observed_filter=None):
    info_comps = {}

    try:
        if best_comp is None or (best_comp not in vsp_ind):
            comp_pos = choose_comp_star_variability(fit_lc_refs, fit_lc_best, info_comps, comp_stars, vsp_comp_stars,
                                                    save)
        else:
            comp_pos = comp_stars[best_comp]
            info_comps[best_comp] = calculate_variablility(fit_lc_refs[best_comp]['myfit'], fit_lc_best)
    except Exception as e:
        log_info(f"Error selecting or calculating variability for comparison star: {e}", warn=True)
        return []

    try:
        comp_star = next(vsp_comp_stars[ckey] for ckey in vsp_comp_stars.keys() if comp_pos == vsp_comp_stars[ckey]['pos'])
        vsp_auid_comp = next(key for key, value in vsp_comp_stars.items() if value['pos'] == comp_pos)
    except StopIteration:
        log_info("Comparison star or VSP AUID not found.", warn=True)
        return []

    try:
        info_comp = info_comps[comp_stars.index(comp_pos)]
        return build_stellar_variability_params_from_fit(
            info_comp['fit_lc'],
            comp_star,
            comp_pos,
            vsp_auid_comp,
            save,
            s_name,
            observed_filter=observed_filter,
        )
    except KeyError as e:
        log_info(f"Key error in processing stellar variability: {e}", warn=True)
        return []
    except Exception as e:
        log_info(f"Error in processing stellar variability: {e}", warn=True)
        return []


# Mid-Transit Time Prior Helper Functions
def numberOfTransitsAway(timeData, period, originalT):
    return int((np.nanmin(timeData) - originalT) / period) + 1


def nearestTransitTime(timeData, period, originalT):
    nearT = ((numberOfTransitsAway(timeData, period, originalT) * period) + originalT)
    return nearT


def save_comp_ra_dec(wcs_file, ra_file, dec_file, comp_coords):
    comp_ra, comp_dec = None, None

    if wcs_file:
        comp_ra = ra_file[int(comp_coords[1])][int(comp_coords[0])]
        comp_dec = dec_file[int(comp_coords[1])][int(comp_coords[0])]

    comp_star = {
        'ra': str(comp_ra) if comp_ra else comp_ra,
        'dec': str(comp_dec) if comp_dec else comp_dec,
        'x': str(comp_coords[0]) if comp_coords[0] else comp_coords[0],
        'y': str(comp_coords[1]) if comp_coords[1] else comp_coords[1]
    }

    return comp_star


def realTimeReduce(i, target_name, p_dict, info_dict, ax, use_nextastro_astrometry=False, multiprocess_transformations=None):
    timeList, airMassList, exptimes, norm_flux = [], [], [], []
    ignore_header_wcs = should_ignore_header_wcs(info_dict.get('ignore_header_wcs'))
    bad_wcs_threshold_fraction = get_bad_wcs_threshold_fraction(info_dict.get('bad_wcs_threshold_percent'))
    pointing_rejection_sigma = get_pointing_rejection_sigma(info_dict.get('pointing_rejection_sigma'))
    detect_bad_pixels_before_photometry = should_detect_bad_pixels_before_photometry(
        info_dict.get('detect_bad_pixels_before_photometry', 'y')
    )
    multiprocess_bad_pixel_precheck = get_multiprocess_bad_pixel_precheck_processes(
        info_dict.get('multiprocess_bad_pixel_precheck', 'n')
    )

    plateStatus.initializeFilenames(info_dict['images'])
    inputfiles = corruption_check(info_dict['images'])
    # time sort images
    times = []
    for ifile in inputfiles:
        plateStatus.setCurrentFilename(ifile)
        extension = 0
        header = fits.getheader(filename=ifile, ext=extension)
        while header['NAXIS'] == 0:
            extension += 1
            header = fits.getheader(filename=ifile, ext=extension)
        obsTime = img_time_bjd_tdb(header, p_dict, info_dict)
        times.append(obsTime)
        plateStatus.setObsTime(obsTime)

    si = np.argsort(times)
    inputfiles = np.array(inputfiles)[si]
    inputfiles, _, dropped_wcs_files = filter_sparse_missing_wcs_frames(
        inputfiles,
        ignore_header_wcs=ignore_header_wcs,
        max_missing_fraction=bad_wcs_threshold_fraction,
    )
    if dropped_wcs_files:
        plateStatus.initializeFilenames(list(inputfiles))
    pointing_precheck_inputfiles = np.array(inputfiles, copy=True)
    pointing_reference_file = inputfiles[0] if len(inputfiles) else None
    inputfiles, _, dropped_pointing_files, pointing_alignment_transforms = filter_pointing_outlier_frames(
        inputfiles,
        pointing_rejection_sigma=pointing_rejection_sigma,
        ignore_header_wcs=ignore_header_wcs,
        return_alignment_transforms=True,
        multiprocess_transformations=multiprocess_transformations,
    )
    if dropped_pointing_files:
        if abort_if_reference_frame_rejected(
            pointing_reference_file,
            dropped_pointing_files,
            ordered_inputfiles=pointing_precheck_inputfiles,
        ):
            ax.clear()
            ax.set_title(target_name)
            ax.set_ylabel('Normalized Flux')
            ax.set_xlabel('Time (JD)')
            ax.text(
                0.5,
                0.5,
                "Reference image rejected by pointing precheck.\nSee log for details.",
                transform=ax.transAxes,
                ha='center',
                va='center',
            )
            plt.close(ax.figure)
            return
        plateStatus.initializeFilenames(list(inputfiles))

    bad_pixel_reference = None
    if detect_bad_pixels_before_photometry:
        log_info(
            "Bad-pixel precheck enabled: scanning frames for persistent isolated high-count outliers before photometry."
        )
        bad_pixel_reference = build_persistent_bad_pixel_map(
            inputfiles,
            load_image_data,
            save_directory=info_dict['save'],
            max_processes=multiprocess_bad_pixel_precheck,
        )
    else:
        log_info("Bad-pixel precheck disabled per optional_info setting.")

    exotic_UIprevTPX = info_dict['tar_coords'][0]
    exotic_UIprevTPY = info_dict['tar_coords'][1]

    plateStatus.setCurrentFilename(inputfiles[0])
    wcs_file = check_wcs(inputfiles[0], info_dict['save'], info_dict['plate_opt'], rt=True,
                         use_nextastro_astrometry=use_nextastro_astrometry,
                         ra=p_dict.get('ra'), dec=p_dict.get('dec'), pixel_scale=info_dict.get('pixel_scale'),
                         ignore_header_wcs=ignore_header_wcs)
    comp_star = info_dict['comp_stars']
    tar_radec, comp_radec = None, []
    first_image = fits.getdata(inputfiles[0])

    if wcs_file:
        wcs_header = get_first_image_header(wcs_file)

        ra_file, dec_file = get_ra_dec(wcs_header, image_shape=first_image.shape)
        tar_radec = (ra_file[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)],
                     dec_file[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)])

        ra = ra_file[int(comp_star[1])][int(comp_star[0])]
        dec = dec_file[int(comp_star[1])][int(comp_star[0])]

        comp_radec.append((ra, dec))

    target_and_comp_radec = None
    if tar_radec is not None and comp_radec:
        target_and_comp_radec = np.array([tar_radec, comp_radec[0]], dtype=float)
    target_and_comp_pixels = np.array(
        [[exotic_UIprevTPX, exotic_UIprevTPY], comp_star],
        dtype=float,
    )

    centroid_reference_image = load_image_data(inputfiles[0])
    centroid_reference_image = repair_bad_pixels_in_frame(centroid_reference_image, bad_pixel_reference)
    targ_sig_xy = fit_centroid(centroid_reference_image, [exotic_UIprevTPX, exotic_UIprevTPY], 0)[3:5]
    del centroid_reference_image

    # aperture and annulus scale factors in PSF sigma units
    aper_sigma = 3 * max(targ_sig_xy)
    annulus_sigma = 10
    fast_aperture_mask = is_fast_aperture_mask_enabled(info_dict.get('fast_aperture_mask'))
    use_adaptive_apertures = is_adaptive_aperture_mode_enabled(info_dict.get('use_adaptive_apertures'))
    aper = np.nan
    annulus = np.nan
    sigma = np.nan
    if use_adaptive_apertures:
        log_info("Adaptive aperture scaling enabled for realtime photometry.")

    # alloc psf fitting param
    psf_data = {
        # x-cent, y-cent, amplitude, sigma-x, sigma-y, rotation, offset
        'target': np.zeros((len(inputfiles), 7)),  # PSF fit
        'comp': np.zeros((len(inputfiles), 7))
    }
    tar_comp_dist = {
        'comp': np.zeros(2, dtype=int)
    }

    # open files, calibrate, align, photometry
    reset_transform_timing_stats()
    reset_photometry_timing_stats()
    multiprocess_alignment_results = None
    use_multiprocess_alignment = multiprocess_transformations is not None and multiprocess_transformations > 0
    if use_multiprocess_alignment:
        multiprocess_alignment_results = build_multiprocess_alignment_results(
            inputfiles,
            multiprocess_transformations,
            target_and_comp_pixels,
            target_and_comp_radec=target_and_comp_radec,
            ignore_header_wcs=ignore_header_wcs,
            bad_pixel_reference=bad_pixel_reference,
            use_fast_centroid_cadence=True,
            use_adaptive_apertures=use_adaptive_apertures,
            compute_fallback_transform=True,
            first_frame_uses_input_comp_pixels=True,
            precomputed_fallback_transforms=pointing_alignment_transforms,
        )
    use_multiprocess_transform_precompute = False
    fallback_transforms = pointing_alignment_transforms
    for i, fileName in enumerate(inputfiles):
        plateStatus.setCurrentFilename(fileName)
        hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False,
                         ignore_missing_end=True)
        frame_fast_centroid = should_use_fast_centroid(i)
        target_fast_centroid = should_use_fast_target_centroid(i, adaptive_apertures=use_adaptive_apertures)

        extension = 0
        image_header = hdul[extension].header
        while image_header["NAXIS"] == 0:
            extension += 1
            image_header = hdul[extension].header

        # TIME
        timeVal = img_time_bjd_tdb(image_header, p_dict, info_dict)
        timeList.append(timeVal)

        # IMAGES
        imageData = hdul[extension].data
        imageData = repair_bad_pixels_in_frame(imageData, bad_pixel_reference)

        if i == 0:
            firstImage = np.copy(imageData)

        if multiprocess_alignment_results is not None:
            apply_parallel_alignment_result(
                multiprocess_alignment_results[i],
                i,
                psf_data,
                tar_comp_dist,
                ['comp'],
            )
        else:
            use_wcs_alignment = False
            if not ignore_header_wcs:
                try:
                    wcs_hdr = search_wcs_from_header(image_header)
                    use_wcs_alignment = wcs_hdr.is_celestial
                except Exception:
                    use_wcs_alignment = False

            if use_wcs_alignment:
                try:
                    if i == 0:
                        tx, ty = exotic_UIprevTPX, exotic_UIprevTPY
                        cx, cy = comp_star
                    else:
                        pix_x, pix_y = wcs_hdr.world_to_pixel_values(
                            target_and_comp_radec[:, 0],
                            target_and_comp_radec[:, 1],
                        )
                        pix_x = np.asarray(pix_x, dtype=float).reshape(-1)
                        pix_y = np.asarray(pix_y, dtype=float).reshape(-1)
                        tx, ty = pix_x[0], pix_y[0]
                        cx, cy = pix_x[1], pix_y[1]

                    projected_coords = np.array([[tx, ty], [cx, cy]], dtype=float)
                    projected_off_frame = any_projected_coord_out_of_frame(projected_coords, imageData.shape)
                    target_seed = choose_centroid_seed_position(
                        [tx, ty],
                        None if i == 0 else psf_data['target'][i - 1],
                    )
                    comp_seed = choose_centroid_seed_position(
                        [cx, cy],
                        None if i == 0 else psf_data['comp'][i - 1],
                    )

                    psf_data['target'][i] = fit_centroid_or_warn_out_of_frame(
                        imageData,
                        target_seed,
                        0,
                        fast_mode=target_fast_centroid,
                    )
                    psf_data['comp'][i] = fit_centroid_or_warn_out_of_frame(
                        imageData,
                        comp_seed,
                        1,
                        fast_mode=frame_fast_centroid,
                    )

                    if i == 0:
                        tar_comp_dist['comp'][0] = abs(int(psf_data['comp'][0][0]) - int(psf_data['target'][0][0]))
                        tar_comp_dist['comp'][1] = abs(int(psf_data['comp'][0][1]) - int(psf_data['target'][0][1]))
                    wcs_alignment_decision = should_keep_header_wcs_alignment(
                        projected_off_frame,
                        i,
                        psf_data['target'][i],
                        previous_target_psf_row=None if i == 0 else psf_data['target'][i - 1],
                        comp_psf_rows={'comp': psf_data['comp'][i]},
                        previous_comp_psf_rows={} if i == 0 else {'comp': psf_data['comp'][i - 1]},
                        expected_offsets={'comp': tar_comp_dist['comp']},
                    )
                    use_wcs_alignment = wcs_alignment_decision['use_wcs_alignment']
                except Exception:
                    use_wcs_alignment = False

            log_alignment_progress(
                i,
                len(inputfiles),
                fileName,
                use_multiprocess_transform_precompute,
            )

            if not use_wcs_alignment:
                cached_tform = fallback_transforms.get(str(fileName)) if fallback_transforms else None
                if cached_tform is not None:
                    tform = cached_tform
                elif i == 0:
                    tform = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])
                else:
                    tform = transformation(imageData, fileName, reference_image=firstImage)

                transformed_coords = np.asarray(tform(target_and_comp_pixels), dtype=float)
                tx, ty = transformed_coords[0]
                target_seed = choose_centroid_seed_position(
                    [tx, ty],
                    None if i == 0 else psf_data['target'][i - 1],
                )
                psf_data['target'][i] = fit_centroid_or_warn_out_of_frame(
                    imageData,
                    target_seed,
                    0,
                    fast_mode=target_fast_centroid,
                )

                cx, cy = transformed_coords[1]
                comp_seed = choose_centroid_seed_position(
                    [cx, cy],
                    None if i == 0 else psf_data['comp'][i - 1],
                )
                psf_data['comp'][i] = fit_centroid_or_warn_out_of_frame(
                    imageData,
                    comp_seed,
                    1,
                    fast_mode=frame_fast_centroid,
                )

                if i == 0:
                    tar_comp_dist['comp'][0] = abs(int(psf_data['comp'][0][0]) - int(psf_data['target'][0][0]))
                    tar_comp_dist['comp'][1] = abs(int(psf_data['comp'][0][1]) - int(psf_data['target'][0][1]))

        # aperture photometry
        frame_sigma = psf_sigma_from_fit(psf_data['target'][i], fallback_sigma=sigma)
        if i == 0:
            sigma = frame_sigma
            if not np.isfinite(sigma) or sigma <= 0:
                log_info("Warning: Initial PSF sigma is invalid; using sigma=1.0 for aperture photometry.", warn=True)
                sigma = 1.0

        if use_adaptive_apertures:
            frame_aper, frame_annulus = resolve_frame_aperture_radii(
                [aper_sigma],
                [annulus_sigma],
                adaptive_apertures=True,
                frame_sigma=frame_sigma,
                fallback_sigma=sigma,
            )
            aper = float(frame_aper[0])
            annulus = float(frame_annulus[0])
        elif i == 0:
            aper, annulus = resolve_frame_aperture_radii(
                [aper_sigma],
                [annulus_sigma],
                adaptive_apertures=True,
                frame_sigma=sigma,
                fallback_sigma=sigma,
            )
            aper = float(aper[0])
            annulus = float(annulus[0])

        comp_frame_sigma = psf_sigma_from_fit(psf_data['comp'][i], fallback_sigma=frame_sigma)
        tFlux = aperPhot(
            imageData,
            0,
            psf_data['target'][i, 0],
            psf_data['target'][i, 1],
            aper,
            annulus,
            fast_mode=fast_aperture_mask,
            sigma_hint=frame_sigma,
        )[0]
        cFlux = aperPhot(
            imageData,
            1,
            psf_data['comp'][i, 0],
            psf_data['comp'][i, 1],
            aper,
            annulus,
            fast_mode=fast_aperture_mask,
            sigma_hint=comp_frame_sigma,
        )[0]
        norm_flux.append(tFlux / cFlux)

        # close file + delete from memory
        hdul.close()
        del hdul
        # Replaced each loop, so clean up
        del imageData

    log_transform_timing_stats('Transformation timing summary (real-time reduce)')
    log_photometry_timing_stats('Photometry timing summary (real-time reduce)')
    log_reduction_timing_overview('Reduction timing overview (real-time reduce)')

    ax.clear()
    ax.set_title(target_name)
    ax.set_ylabel('Normalized Flux')
    ax.set_xlabel('Time (JD)')
    ax.plot(timeList, norm_flux, 'bo')


def fit_lightcurve(times, tFlux, cFlux, airmass, ld, pDict, jd_times=None,
                   allow_mid_transit_range_warning=True, disable_vertical_flux_normalization=False,
                   final_fit_mode='lm',
                   use_impactparameter_rather_than_inclination_to_fit=True,
                   plot_time_range=None,
                   use_eebls_to_initialize_tmid_and_bounds=True,
                   compute_eebls_diagnostics=False):
    plot_time_range = np.asarray(times if plot_time_range is None else plot_time_range, dtype=float)
    prepared = prepare_lightcurve_fit_input_series(
        times,
        tFlux,
        cFlux,
        airmass,
        jd_times=jd_times,
    )
    if not prepared.get('applied'):
        return None, None, None

    filter_diagnostics = prepared['filter_diagnostics']
    debug_times = prepared['debug_times']
    debug_target_flux = prepared['debug_target_flux']
    debug_comp_flux = prepared['debug_comp_flux']
    debug_raw_ratio = prepared['debug_raw_ratio']
    debug_initial_sigma_keep_mask = prepared['initial_sigma_keep_mask']
    arrayFinalFlux = prepared['flux']
    f1 = prepared['target_flux']
    f2 = prepared['comp_flux']
    arrayNormUnc = prepared['unc']
    arrayTimes = prepared['time']
    arrayJDTimes = prepared['jd_time']
    arrayAirmass = prepared['airmass']
    skip_airmass_fit = prepared['skip_airmass_fit']


    # -----LM LIGHTCURVE FIT--------------------------------------
    prior = {
        'rprs': pDict['rprs'],  # Rp/Rs
        'ars': pDict['aRs'],  # a/Rs
        'per': pDict['pPer'],  # Period [day]
        'inc': pDict['inc'],  # Inclination [deg]
        'u0': ld[0], 'u1': ld[1], 'u2': ld[2], 'u3': ld[3],  # limb darkening (nonlinear)
        'ecc': pDict['ecc'],  # Eccentricity
        'omega': pDict['omega'],  # Arg of periastron
        'tmid': pDict['midT'],  # time of mid transit [day]
        'a2': 0,  # Flux lower bound
    }

    arrayPhases = (arrayTimes - pDict['midT']) / prior['per']
    expected_duration = estimate_transit_duration_from_prior_geometry(prior)
    tmid_search_summary = estimate_ephemeris_tmid_and_bounds(
        arrayTimes,
        pDict['midT'],
        prior['per'],
        pDict['midTUnc'],
        pDict['pPerUnc'],
        expected_duration=expected_duration,
        sigma_multiplier=25.0,
    )
    prior['tmid'] = tmid_search_summary['tmid']
    lower, upper = tmid_search_summary['bounds']
    ephemeris_tmid_search_summary = tmid_search_summary

    if (
        allow_mid_transit_range_warning
        and np.floor(arrayPhases).max() - np.floor(arrayPhases).min() == 0
    ):
        log_mid_transit_range_warning_once(arrayTimes, prior['tmid'])

    if tmid_search_summary.get('duration_capped'):
        log_info(tmid_search_summary['note'])
    eebls_search_summary = None
    if use_eebls_to_initialize_tmid_and_bounds or compute_eebls_diagnostics:
        eebls_search_summary = estimate_tmid_and_bounds_with_eebls(
            arrayTimes,
            arrayFinalFlux,
            arrayNormUnc,
            prior,
            [lower, upper],
        )
        if use_eebls_to_initialize_tmid_and_bounds and eebls_search_summary.get('applied'):
            tmid_search_summary = eebls_search_summary
            prior['tmid'] = tmid_search_summary['tmid']
            lower, upper = tmid_search_summary['bounds']

    mybounds = build_initial_transit_bounds(
        prior,
        [lower, upper],
        ars_unc=pDict.get('aRsUnc'),
    )
    apply_vertical_flux_normalization_bound(
        prior,
        mybounds,
        arrayFinalFlux,
        disable_vertical_flux_normalization,
    )
    if not skip_airmass_fit:
        mybounds['a2'] = [-1, 1]

    if arrayTimes.shape[0] < LIGHTCURVE_MIN_VALID_POINTS:
        return None, None, None

    if np.isnan(arrayTimes).any() or np.isnan(arrayFinalFlux).any() or np.isnan(arrayNormUnc).any():
        log_info("\nWarning: NANs in time, flux or error", warn=True)

    myfit = lc_fitter(
        arrayTimes,
        arrayFinalFlux,
        arrayNormUnc,
        arrayAirmass,
        prior,
        mybounds,
        jd_times=arrayJDTimes,
        mode='lm',
        use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
    )
    myfit = apply_plot_time_range(myfit, plot_time_range)
    annotate_airmass_fit(myfit, arrayAirmass, skip_airmass_fit)
    annotate_lightcurve_filter_diagnostics(myfit, filter_diagnostics)
    annotate_lightcurve_tmid_search(myfit, tmid_search_summary)
    annotate_lightcurve_eebls_diagnostic(myfit, eebls_search_summary)

    if (
        myfit is not None
        and hasattr(myfit, 'residuals')
        and hasattr(myfit, 'phase')
        and np.shape(myfit.residuals) == np.shape(arrayTimes)
        and np.shape(myfit.phase) == np.shape(arrayTimes)
    ):
        phase_clip_mask = phase_bin_sigma_clip(myfit.residuals, myfit.phase, sigma=3, bins=10)
        min_required_points = max(len(mybounds) + 1, 5)
        if np.count_nonzero(~phase_clip_mask) >= min_required_points and np.any(phase_clip_mask):
            filter_diagnostics.append(build_time_rejection_diagnostic(
                "Phase-binned residual clip",
                arrayTimes,
                ~phase_clip_mask,
                note="Dropped phase-binned residual outliers after the initial LM fit before refitting.",
            ))
            arrayFinalFlux = arrayFinalFlux[~phase_clip_mask]
            arrayNormUnc = arrayNormUnc[~phase_clip_mask]
            arrayTimes = arrayTimes[~phase_clip_mask]
            arrayJDTimes = arrayJDTimes[~phase_clip_mask]
            arrayAirmass = arrayAirmass[~phase_clip_mask]
            f1 = f1[~phase_clip_mask]
            f2 = f2[~phase_clip_mask]

            myfit = lc_fitter(
                arrayTimes,
                arrayFinalFlux,
                arrayNormUnc,
                arrayAirmass,
                prior,
                mybounds,
                jd_times=arrayJDTimes,
                mode='lm',
                use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
            )
            myfit = apply_plot_time_range(myfit, plot_time_range)
            annotate_airmass_fit(myfit, arrayAirmass, skip_airmass_fit)
            annotate_lightcurve_filter_diagnostics(myfit, filter_diagnostics)
            annotate_lightcurve_tmid_search(myfit, tmid_search_summary)
            annotate_lightcurve_eebls_diagnostic(myfit, eebls_search_summary)

    debug_phase_clip_keep_mask = np.ones(np.count_nonzero(debug_initial_sigma_keep_mask), dtype=bool)
    if final_fit_mode == 'ns' and myfit is not None:
        duration_prior = build_single_transit_duration_prior(pDict)
        pre_ultranest_coverage_assessment = build_expected_transit_coverage_assessment(
            arrayTimes,
            prior,
            flux_values=arrayFinalFlux,
            flux_errors=arrayNormUnc,
            tmid_search_summary=ephemeris_tmid_search_summary,
            duration_prior=duration_prior,
        )
        log_expected_transit_coverage_assessment(pre_ultranest_coverage_assessment)
        nested_refinement = build_nested_tmid_refinement_from_initial_fit(
            arrayTimes,
            arrayFinalFlux,
            arrayNormUnc,
            prior,
            mybounds,
            myfit,
        )
        myfit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
            arrayTimes,
            arrayFinalFlux,
            arrayNormUnc,
            arrayAirmass,
            nested_refinement['prior'],
            nested_refinement['bounds'],
            jd_times=arrayJDTimes,
            use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
            duration_prior=duration_prior,
        )
        annotate_pre_ultranest_transit_coverage(myfit, pre_ultranest_coverage_assessment)
        myfit = apply_plot_time_range(myfit, plot_time_range)
        annotate_airmass_fit(myfit, arrayAirmass, skip_airmass_fit)
        annotate_lightcurve_filter_diagnostics(myfit, filter_diagnostics)
        if myfit is not None:
            annotate_lightcurve_tmid_search(myfit, tmid_search_summary)
            annotate_lightcurve_eebls_diagnostic(myfit, eebls_search_summary)
            annotate_nested_tmid_refinement(
                myfit,
                nested_refinement.get('applied', False),
                note=nested_refinement.get('note'),
                original_tmid_bounds=nested_refinement.get('original_tmid_bounds'),
                refined_tmid_bounds=nested_refinement.get('refined_tmid_bounds'),
            )

    if myfit is not None:
        if 'phase_clip_mask' in locals():
            debug_phase_clip_keep_mask = np.asarray(~phase_clip_mask, dtype=bool).copy()
        annotate_selected_photometry_debug(
            myfit,
            debug_times,
            debug_target_flux,
            debug_comp_flux,
            debug_raw_ratio,
            debug_initial_sigma_keep_mask,
            phase_clip_keep_mask_on_sigma_filtered=debug_phase_clip_keep_mask,
        )
        annotate_transit_qc_expected_values(myfit, pDict)
        annotate_transit_detection_qc(myfit)

    return myfit, f1, f2


def diagnose_lightcurve_fit_inputs(times, tflux, cflux, airmass, enforce_relative_flux_max=True):
    times = np.asarray(times, dtype=float)
    tflux = np.asarray(tflux, dtype=float)
    cflux = np.asarray(cflux, dtype=float)
    airmass = np.asarray(airmass, dtype=float)

    diagnostics = {
        'input_point_count': int(times.shape[0]),
        'has_reference_flux': False,
        'relative_flux_point_count': 0,
        'sigma_clip_point_count': 0,
        'usable_point_count': 0,
        'failed_stage': None,
        'failure_reason': None,
    }

    if diagnostics['input_point_count'] <= 1:
        diagnostics.update({
            'relative_flux_point_count': diagnostics['input_point_count'],
            'failed_stage': 'coverage',
            'failure_reason': (
                f"only {diagnostics['input_point_count']} frame(s) remained after masking invalid "
                "comparison flux; need at least 2 to fit."
            ),
        })
        return diagnostics

    si = np.argsort(times)
    times_sorted = times[si]
    tflux_sorted = tflux[si]
    cflux_sorted = cflux[si]
    with np.errstate(divide='ignore', invalid='ignore'):
        flux_ratio_sorted = np.divide(tflux_sorted, cflux_sorted)

    has_reference_flux = not np.allclose(cflux_sorted, 1.0)
    diagnostics['has_reference_flux'] = bool(has_reference_flux)
    diagnostics['relative_flux_point_count'] = int(times_sorted.shape[0])

    if has_reference_flux:
        finite_ratio_mask = np.isfinite(flux_ratio_sorted)
        nonpositive_ratio_mask = finite_ratio_mask & np.less_equal(flux_ratio_sorted, 0)
        finite_ratio_values = flux_ratio_sorted[finite_ratio_mask]
        ratio_range_text = "finite ratio range=n/a"
        if finite_ratio_values.size:
            ratio_range_text = (
                f"finite ratio range={np.nanmin(finite_ratio_values):.4f} to "
                f"{np.nanmax(finite_ratio_values):.4f}"
            )
        nonfinite_ratio_count = int(np.count_nonzero(~finite_ratio_mask))
        nonpositive_ratio_count = int(np.count_nonzero(nonpositive_ratio_mask))
        rejected_ratio_count = nonfinite_ratio_count + nonpositive_ratio_count
        relative_flux_mask = valid_flux_ratio_mask(flux_ratio_sorted)
        rejection_detail = (
            f"(non-finite={nonfinite_ratio_count}, non-positive={nonpositive_ratio_count}, "
            f"{ratio_range_text})."
        )
        rejection_context = "invalid target/reference ratio screening"
        diagnostics['relative_flux_point_count'] = int(np.count_nonzero(relative_flux_mask))
        times_sorted = times_sorted[relative_flux_mask]
        tflux_sorted = tflux_sorted[relative_flux_mask]
        cflux_sorted = cflux_sorted[relative_flux_mask]
        flux_ratio_sorted = flux_ratio_sorted[relative_flux_mask]
        airmass_sorted = airmass[si][relative_flux_mask]
        if diagnostics['relative_flux_point_count'] <= 1:
            diagnostics.update({
                'failed_stage': 'relative_flux_filter',
                'failure_reason': (
                    "relative-flux filtering left "
                    f"{diagnostics['relative_flux_point_count']} usable point(s); "
                    f"rejected {rejected_ratio_count}/{diagnostics['input_point_count']} frame(s) "
                    f"during {rejection_context} {rejection_detail}"
                ),
            })
            return diagnostics
    else:
        airmass_sorted = airmass[si]

    dt = np.mean(np.diff(times_sorted))
    if np.isfinite(dt) and dt > 0:
        ndt = int(25. / 24. / 60. / dt) * 2 + 1
    else:
        ndt = 5
    if ndt > len(times_sorted):
        ndt = int(len(times_sorted) / 4) * 2 + 1
    filtered_data = sigma_clip(flux_ratio_sorted, sigma=3, dt=max(5, ndt), times=times_sorted)
    valid_mask = ~filtered_data
    diagnostics['sigma_clip_point_count'] = int(np.count_nonzero(valid_mask))
    if diagnostics['sigma_clip_point_count'] <= 1:
        diagnostics.update({
            'failed_stage': 'sigma_clip',
            'failure_reason': (
                "sigma clipping left "
                f"{diagnostics['sigma_clip_point_count']} usable point(s); not enough data remained "
                "for a lightcurve fit."
            ),
        })
        return diagnostics

    arrayFinalFlux = flux_ratio_sorted[valid_mask]
    f1 = tflux_sorted[valid_mask]
    sigf1 = f1 ** 0.5
    f2 = cflux_sorted[valid_mask]
    sigf2 = f2 ** 0.5
    if np.sum(cflux) == len(cflux):
        arrayNormUnc = sigf1
    else:
        arrayNormUnc = np.sqrt((sigf1 / f2) ** 2 + (sigf2 * f1 / f2 ** 2) ** 2)
    arrayTimes = times_sorted[valid_mask]
    arrayAirmass = airmass_sorted[valid_mask]

    nanmask = np.isnan(arrayFinalFlux) | np.isnan(arrayNormUnc) | np.isnan(arrayTimes) | np.isnan(arrayAirmass)
    nanmask = nanmask | np.less_equal(arrayFinalFlux, 0) | np.less_equal(arrayNormUnc, 0)
    nanmask = nanmask | np.isinf(arrayFinalFlux) | np.isinf(arrayNormUnc) | np.isinf(arrayTimes) | np.isinf(
        arrayAirmass
    )
    diagnostics['usable_point_count'] = int(np.count_nonzero(~nanmask))
    if diagnostics['usable_point_count'] <= 1:
        diagnostics.update({
            'failed_stage': 'nan_filter',
            'failure_reason': (
                "filtering non-finite or non-positive flux/uncertainty values left "
                f"{diagnostics['usable_point_count']} usable point(s); need at least 2."
            ),
        })
    elif diagnostics['usable_point_count'] < LIGHTCURVE_MIN_VALID_POINTS:
        diagnostics.update({
            'failed_stage': 'minimum_points',
            'failure_reason': (
                f"only {diagnostics['usable_point_count']} usable point(s) remained after filtering; "
                f"need at least {LIGHTCURVE_MIN_VALID_POINTS} for a lightcurve fit."
            ),
        })

    return diagnostics


def ensure_lightcurve_fit_failure_reason(diagnostics, fit_result, failed_stage, failure_reason):
    diagnostics = {} if diagnostics is None else dict(diagnostics)
    if fit_result is None and diagnostics.get('failure_reason') is None:
        diagnostics.update({
            'failed_stage': failed_stage,
            'failure_reason': failure_reason,
        })
    return diagnostics


def prepare_lightcurve_fit_input_series(
    times,
    target_flux,
    comp_flux,
    airmass,
    jd_times=None,
):
    times = np.asarray(times, dtype=float)
    target_flux = np.asarray(target_flux, dtype=float)
    comp_flux = np.asarray(comp_flux, dtype=float)
    airmass = np.asarray(airmass, dtype=float)
    jd_times_array = None if jd_times is None else np.asarray(jd_times, dtype=float)

    prepared = {
        'applied': False,
        'failure_reason': None,
        'filter_diagnostics': [],
        'debug_times': np.array([], dtype=float),
        'debug_target_flux': np.array([], dtype=float),
        'debug_comp_flux': np.array([], dtype=float),
        'debug_raw_ratio': np.array([], dtype=float),
        'initial_sigma_keep_mask': np.array([], dtype=bool),
        'time': np.array([], dtype=float),
        'flux': np.array([], dtype=float),
        'unc': np.array([], dtype=float),
        'jd_time': None,
        'airmass': np.array([], dtype=float),
        'target_flux': np.array([], dtype=float),
        'comp_flux': np.array([], dtype=float),
        'source_indices': np.array([], dtype=int),
        'skip_airmass_fit': False,
        'approximate_baseline_level': np.nan,
    }

    if not (
        times.ndim == target_flux.ndim == comp_flux.ndim == airmass.ndim == 1
        and times.shape == target_flux.shape == comp_flux.shape == airmass.shape
    ):
        prepared['failure_reason'] = (
            "lightcurve inputs must be 1D arrays with matching lengths before fitting."
        )
        return prepared

    if jd_times_array is not None and jd_times_array.shape != times.shape:
        jd_times_array = None

    plot_indices = np.argsort(times)
    times_sorted = times[plot_indices]
    target_flux_sorted = target_flux[plot_indices]
    comp_flux_sorted = comp_flux[plot_indices]
    source_indices = np.asarray(plot_indices, dtype=int)
    with np.errstate(divide='ignore', invalid='ignore'):
        flux_ratio_sorted = np.divide(target_flux_sorted, comp_flux_sorted)

    filter_diagnostics = []
    has_reference_flux = not np.allclose(comp_flux_sorted, 1.0)
    if has_reference_flux:
        flux_ratio_mask = valid_flux_ratio_mask(flux_ratio_sorted)
        filter_diagnostics.append(build_time_rejection_diagnostic(
            "Target/reference ratio filter",
            times_sorted,
            flux_ratio_mask,
            note="Dropped non-finite or non-positive target/reference ratios before fitting.",
        ))
        times_sorted = times_sorted[flux_ratio_mask]
        target_flux_sorted = target_flux_sorted[flux_ratio_mask]
        comp_flux_sorted = comp_flux_sorted[flux_ratio_mask]
        flux_ratio_sorted = flux_ratio_sorted[flux_ratio_mask]
        source_indices = source_indices[flux_ratio_mask]
        if jd_times_array is None:
            jd_times_sorted = times_sorted.copy()
        else:
            jd_times_sorted = jd_times_array[plot_indices][flux_ratio_mask]
        airmass_sorted = airmass[plot_indices][flux_ratio_mask]
    else:
        jd_times_sorted = times_sorted.copy() if jd_times_array is None else jd_times_array[plot_indices]
        airmass_sorted = airmass[plot_indices]

    if len(times_sorted) <= 1:
        prepared['failure_reason'] = "too few valid points remained after the target/reference ratio filter."
        prepared['filter_diagnostics'] = filter_diagnostics
        return prepared

    debug_times = np.asarray(times_sorted, dtype=float).copy()
    debug_target_flux = np.asarray(target_flux_sorted, dtype=float).copy()
    debug_comp_flux = np.asarray(comp_flux_sorted, dtype=float).copy()
    debug_raw_ratio = np.asarray(flux_ratio_sorted, dtype=float).copy()

    dt = np.mean(np.diff(times_sorted))
    ndt = int(25. / 24. / 60. / dt) * 2 + 1
    if ndt > len(times_sorted):
        ndt = int(len(times_sorted)/4) * 2 + 1
    filtered_data = sigma_clip(flux_ratio_sorted, sigma=3, dt=max(5, ndt), times=times_sorted)
    valid_mask = ~filtered_data
    initial_sigma_keep_mask = np.asarray(valid_mask, dtype=bool).copy()
    filter_diagnostics.append(build_time_rejection_diagnostic(
        "Initial sigma clip",
        times_sorted,
        valid_mask,
        note="Dropped 3-sigma target/reference-ratio outliers before the first lightcurve fit.",
    ))

    flux = flux_ratio_sorted[valid_mask]
    filtered_target_flux = target_flux_sorted[valid_mask]
    filtered_comp_flux = comp_flux_sorted[valid_mask]
    if np.sum(comp_flux) == len(comp_flux):
        unc = filtered_target_flux ** 0.5
    else:
        sigf1 = filtered_target_flux ** 0.5
        sigf2 = filtered_comp_flux ** 0.5
        unc = np.sqrt((sigf1 / filtered_comp_flux) ** 2 + (sigf2 * filtered_target_flux / filtered_comp_flux ** 2) ** 2)
    fit_times = times_sorted[valid_mask]
    fit_jd_times = jd_times_sorted[valid_mask]
    fit_airmass = airmass_sorted[valid_mask]
    source_indices = source_indices[valid_mask]

    nanmask = np.isnan(flux) | np.isnan(unc) | np.isnan(fit_times) | np.isnan(fit_airmass) | np.less_equal(flux, 0) | np.less_equal(unc, 0)
    nanmask = nanmask | np.isinf(flux) | np.isinf(unc) | np.isinf(fit_times) | np.isinf(fit_airmass)
    filter_diagnostics.append(build_time_rejection_diagnostic(
        "Finite/positive photometry filter",
        fit_times,
        ~nanmask,
        note="Dropped non-finite or non-positive flux, uncertainty, time, or airmass values.",
    ))

    if np.sum(~nanmask) <= 1:
        prepared['failure_reason'] = "too few valid points remained after removing non-finite or non-positive photometry."
        prepared.update({
            'filter_diagnostics': filter_diagnostics,
            'debug_times': debug_times,
            'debug_target_flux': debug_target_flux,
            'debug_comp_flux': debug_comp_flux,
            'debug_raw_ratio': debug_raw_ratio,
            'initial_sigma_keep_mask': initial_sigma_keep_mask,
        })
        return prepared

    normalized_flux, normalized_unc, approximate_baseline_level = normalize_flux_series_to_approximate_unity(
        flux[~nanmask],
        unc[~nanmask],
    )

    prepared.update({
        'applied': True,
        'filter_diagnostics': filter_diagnostics,
        'debug_times': debug_times,
        'debug_target_flux': debug_target_flux,
        'debug_comp_flux': debug_comp_flux,
        'debug_raw_ratio': debug_raw_ratio,
        'initial_sigma_keep_mask': initial_sigma_keep_mask,
        'time': fit_times[~nanmask],
        'flux': normalized_flux,
        'unc': normalized_unc,
        'jd_time': fit_jd_times[~nanmask],
        'airmass': fit_airmass[~nanmask],
        'target_flux': filtered_target_flux[~nanmask],
        'comp_flux': filtered_comp_flux[~nanmask],
        'source_indices': source_indices[~nanmask],
        'skip_airmass_fit': should_skip_airmass_fit(fit_airmass[~nanmask]),
        'approximate_baseline_level': approximate_baseline_level,
    })
    return prepared


def cheap_lightcurve_prescore(tFlux, cFlux, airmass, enforce_relative_flux_max=True):
    with np.errstate(divide='ignore', invalid='ignore'):
        flux_ratio = np.divide(tFlux, cFlux)

    finite_mask = valid_flux_ratio_mask(flux_ratio) & np.isfinite(airmass)
    if np.count_nonzero(finite_mask) < 5:
        return np.inf

    x_vals = airmass[finite_mask]
    y_vals = flux_ratio[finite_mask]

    if should_skip_airmass_fit(x_vals):
        detrended = y_vals / bn.nanmedian(y_vals)
    else:
        slope, intercept = np.polyfit(x_vals, y_vals, 1)
        trend = slope * x_vals + intercept
        with np.errstate(divide='ignore', invalid='ignore'):
            detrended = np.divide(y_vals, trend)

    return bn.nanstd(detrended)


def evaluate_lightcurve_candidate(task):
    (
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times,
        plot_time_range,
        disable_vertical_flux_normalization,
        use_impactparameter_rather_than_inclination_to_fit,
        use_eebls_to_initialize_tmid_and_bounds,
        compute_eebls_diagnostics,
    ) = task
    fit_diagnostics = diagnose_lightcurve_fit_inputs(
        times,
        tflux,
        cflux,
        airmass,
        enforce_relative_flux_max=False,
    )
    myfit, tflux_fit, cflux_fit = fit_lightcurve(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times,
        allow_mid_transit_range_warning=False,
        disable_vertical_flux_normalization=disable_vertical_flux_normalization,
        final_fit_mode='ns',
        use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
        plot_time_range=plot_time_range,
        use_eebls_to_initialize_tmid_and_bounds=use_eebls_to_initialize_tmid_and_bounds,
        compute_eebls_diagnostics=compute_eebls_diagnostics,
    )
    fit_diagnostics = ensure_lightcurve_fit_failure_reason(
        fit_diagnostics,
        myfit,
        failed_stage='lightcurve_fit',
        failure_reason="the lightcurve fitter did not converge to a usable solution.",
    )
    transit_qc_failure_reason = lightcurve_fit_transit_qc_failure_reason(myfit)
    if transit_qc_failure_reason is not None:
        fit_diagnostics = dict(fit_diagnostics)
        fit_diagnostics.update({
            'failed_stage': 'transit_qc',
            'failure_reason': transit_qc_failure_reason,
        })

    return {
        'myfit': myfit,
        'accepted': myfit is not None and transit_qc_failure_reason is None,
        'eebls_snr': extract_lightcurve_fit_eebls_snr(myfit),
        'transit_delta_bic': extract_lightcurve_fit_transit_delta_bic(myfit),
        'residual_scatter': extract_lightcurve_fit_residual_scatter(myfit),
        'ktmf_metric': extract_lightcurve_fit_ktmf_metric(myfit),
        'ktmf_contributions': extract_lightcurve_fit_ktmf_contributions(myfit),
        'transit_qc_status': getattr(myfit, 'transit_qc_status', None),
        'transit_qc_summary': getattr(myfit, 'transit_qc_summary', None),
        'rejected_by_transit_qc': transit_qc_failure_reason is not None,
        'fit_diagnostics': fit_diagnostics,
        'failure_reason': fit_diagnostics.get('failure_reason'),
        'fit_point_count': 0 if tflux_fit is None else int(len(tflux_fit)),
    }, tflux_fit, cflux_fit


def build_target_fit_candidate_jobs(psf_data, aper_data, apers, annuli, airmass, comp_stars, sigma,
                                    require_comp_star=True,
                                    skip_low_comparison_coverage_rejection=False,
                                    use_psf_photometry=True,
                                    use_aperture_photometry=True):
    candidate_jobs = []
    comp_star_count = len(comp_stars)

    if use_psf_photometry and comp_star_count > 0:
        target_flux = 2 * np.pi * psf_data['target'][:, 2] * psf_data['target'][:, 3] * psf_data['target'][:, 4]
        target_flux_mask = robust_flux_floor_mask(target_flux)
        psf_comp_flux_map = {
            f"comp{comp_idx + 1}": 2 * np.pi * psf_data[f"comp{comp_idx + 1}"][:, 2]
            * psf_data[f"comp{comp_idx + 1}"][:, 3]
            * psf_data[f"comp{comp_idx + 1}"][:, 4]
            for comp_idx in range(comp_star_count)
        }
        psf_comp_coverage = comparison_star_coverage_summary(
            psf_comp_flux_map,
            skip_rejection=skip_low_comparison_coverage_rejection,
            validity_mask_func=robust_flux_floor_mask,
        )
        for comp_idx in range(comp_star_count):
            ckey = f"comp{comp_idx + 1}"
            if psf_comp_coverage[ckey]['coverage_rejected']:
                continue

            comp_flux = psf_comp_flux_map[ckey]
            psf_mask = target_flux_mask & robust_flux_floor_mask(comp_flux)
            candidate_jobs.append({
                'method': 'psf',
                'a': None,
                'an': None,
                'aper': 0.0,
                'annulus': float(15 * sigma),
                'comp_index': comp_idx,
                'ckey': ckey,
                'mask': psf_mask,
                'coverage_count': psf_comp_coverage[ckey]['coverage_count'],
                'coverage_total_frame_count': psf_comp_coverage[ckey]['coverage_total_frame_count'],
                'coverage_reference_count': psf_comp_coverage[ckey]['coverage_reference_count'],
                'coverage_min_required_count': psf_comp_coverage[ckey]['coverage_min_required_count'],
                'coverage_rejected': psf_comp_coverage[ckey]['coverage_rejected'],
                'prescore': cheap_lightcurve_prescore(
                    target_flux[psf_mask],
                    comp_flux[psf_mask],
                    airmass[psf_mask],
                    enforce_relative_flux_max=False,
                ),
            })

    if use_aperture_photometry and aper_data is not None and apers is not None and annuli is not None:
        for a, aper in enumerate(apers):
            for an, annulus in enumerate(annuli):
                target_flux = aper_data['target'][:, a, an]
                aperture_comp_flux_map = {
                    f"comp{comp_idx + 1}": aper_data[f"comp{comp_idx + 1}"][:, a, an]
                    for comp_idx in range(comp_star_count)
                }
                aperture_comp_coverage = comparison_star_coverage_summary(
                    aperture_comp_flux_map,
                    skip_rejection=skip_low_comparison_coverage_rejection,
                )

                if not require_comp_star:
                    candidate_jobs.append({
                        'method': 'aperture',
                        'a': a,
                        'an': an,
                        'aper': float(aper),
                        'annulus': float(annulus),
                        'comp_index': None,
                        'ckey': None,
                        'mask': np.ones(target_flux.shape[0], dtype=bool),
                        'coverage_count': int(target_flux.shape[0]),
                        'coverage_total_frame_count': int(target_flux.shape[0]),
                        'coverage_reference_count': float(target_flux.shape[0]),
                        'coverage_min_required_count': LIGHTCURVE_MIN_VALID_POINTS,
                        'coverage_rejected': False,
                        'prescore': cheap_lightcurve_prescore(
                            target_flux,
                            np.ones(target_flux.shape[0]),
                            airmass,
                            enforce_relative_flux_max=False,
                        ),
                    })

                for comp_idx in range(comp_star_count):
                    ckey = f"comp{comp_idx + 1}"
                    if aperture_comp_coverage[ckey]['coverage_rejected']:
                        continue

                    comp_series = aperture_comp_flux_map[ckey]
                    aper_mask = valid_comparison_frame_mask(comp_series)
                    candidate_jobs.append({
                        'method': 'aperture',
                        'a': a,
                        'an': an,
                        'aper': float(aper),
                        'annulus': float(annulus),
                        'comp_index': comp_idx,
                        'ckey': ckey,
                        'mask': aper_mask,
                        'coverage_count': aperture_comp_coverage[ckey]['coverage_count'],
                        'coverage_total_frame_count': aperture_comp_coverage[ckey]['coverage_total_frame_count'],
                        'coverage_reference_count': aperture_comp_coverage[ckey]['coverage_reference_count'],
                        'coverage_min_required_count': aperture_comp_coverage[ckey]['coverage_min_required_count'],
                        'coverage_rejected': aperture_comp_coverage[ckey]['coverage_rejected'],
                        'prescore': cheap_lightcurve_prescore(
                            target_flux[aper_mask],
                            comp_series[aper_mask],
                            airmass[aper_mask],
                            enforce_relative_flux_max=False,
                        ),
                    })

    return candidate_jobs


def target_fit_candidate_task(candidate, times, jd_times, airmass, ld, p_dict, psf_data, aper_data,
                              plot_time_range=None,
                              disable_vertical_flux_normalization=False,
                              use_impactparameter_rather_than_inclination_to_fit=True,
                              use_eebls_to_initialize_tmid_and_bounds=True,
                              compute_eebls_diagnostics=True):
    candidate_mask = np.asarray(candidate['mask'], dtype=bool)

    if candidate['method'] == 'psf':
        target_flux = 2 * np.pi * psf_data['target'][:, 2] * psf_data['target'][:, 3] * psf_data['target'][:, 4]
        if candidate['ckey'] is None:
            comp_flux = np.ones(target_flux.shape[0], dtype=float)
        else:
            comp_flux = (
                2 * np.pi * psf_data[candidate['ckey']][:, 2]
                * psf_data[candidate['ckey']][:, 3]
                * psf_data[candidate['ckey']][:, 4]
            )
    else:
        target_flux = aper_data['target'][:, candidate['a'], candidate['an']]
        if candidate['ckey'] is None:
            comp_flux = np.ones(target_flux.shape[0], dtype=float)
        else:
            comp_flux = aper_data[candidate['ckey']][:, candidate['a'], candidate['an']]

    return (
        times[candidate_mask],
        target_flux[candidate_mask],
        comp_flux[candidate_mask],
        airmass[candidate_mask],
        ld,
        p_dict,
        jd_times[candidate_mask],
        plot_time_range,
        disable_vertical_flux_normalization,
        use_impactparameter_rather_than_inclination_to_fit,
        use_eebls_to_initialize_tmid_and_bounds,
        compute_eebls_diagnostics,
    )


def run_target_driven_photometry_search(times, jd_times, airmass, ld, p_dict, comp_stars, psf_data, aper_data,
                                        apers, annuli, sigma,
                                        require_comp_star=True,
                                        plot_time_range=None,
                                        disable_vertical_flux_normalization=False,
                                        skip_low_comparison_coverage_rejection=False,
                                        use_psf_photometry=True,
                                        use_aperture_photometry=True,
                                        multiprocess_lightcurve_fits=None,
                                        use_impactparameter_rather_than_inclination_to_fit=True,
                                        use_eebls_to_initialize_tmid_and_bounds=True,
                                        pick_comparison_by_eebls_snr=True):
    candidate_jobs = build_target_fit_candidate_jobs(
        psf_data,
        aper_data,
        apers,
        annuli,
        airmass,
        comp_stars,
        sigma,
        require_comp_star=require_comp_star,
        skip_low_comparison_coverage_rejection=skip_low_comparison_coverage_rejection,
        use_psf_photometry=use_psf_photometry,
        use_aperture_photometry=use_aperture_photometry,
    )
    evaluated_candidates = list(candidate_jobs)
    for candidate_order, candidate in enumerate(evaluated_candidates):
        candidate.setdefault('candidate_order', candidate_order)
    if not evaluated_candidates:
        return {
            'candidate_jobs': candidate_jobs,
            'evaluated_candidates': evaluated_candidates,
            'shortlist': evaluated_candidates,
            'candidate_summaries': [],
            'best_candidate': None,
            'best_fit_lc': None,
            'selected_ktmf_metric': np.nan,
            'selected_transit_delta_bic': np.nan,
            'selection_metric': 'ktmf',
            'selected_eebls_snr': np.nan,
            'flux_tar': None,
            'flux_ref': None,
        }

    fit_tasks = [
        target_fit_candidate_task(
            candidate,
            times,
            jd_times,
            airmass,
            ld,
            p_dict,
            psf_data,
            aper_data,
            plot_time_range=plot_time_range,
            disable_vertical_flux_normalization=disable_vertical_flux_normalization,
            use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
            use_eebls_to_initialize_tmid_and_bounds=use_eebls_to_initialize_tmid_and_bounds,
            compute_eebls_diagnostics=True,
        )
        for candidate in evaluated_candidates
    ]

    if multiprocess_lightcurve_fits is not None and multiprocess_lightcurve_fits > 0:
        log_info(f"Using multiprocessing for candidate lightcurve fits ({multiprocess_lightcurve_fits} processes).")
        with suppress_tk_cleanup_during_process_pool():
            with ProcessPoolExecutor(
                max_workers=multiprocess_lightcurve_fits,
                initializer=suppress_inherited_tk_cleanup_in_worker,
            ) as executor:
                fit_results = list(executor.map(evaluate_lightcurve_candidate, fit_tasks))
    else:
        fit_results = [evaluate_lightcurve_candidate(task) for task in fit_tasks]

    candidate_summaries = []
    successful_candidates = []
    for candidate, result in zip(evaluated_candidates, fit_results):
        fit_meta, tflux_fit, cflux_fit = result
        summary = summarize_target_fit_candidate(candidate, fit_meta, comp_stars)
        candidate_summaries.append(summary)

        if fit_meta is None or fit_meta.get('myfit') is None or not fit_meta.get('accepted', True):
            continue
        successful_candidates.append((summary, candidate, fit_meta, tflux_fit, cflux_fit))

    best_candidate = None
    best_fit_lc = None
    best_ktmf_metric = np.nan
    best_transit_delta_bic = np.nan
    best_tflux = None
    best_cflux = None
    selection_metric = 'ktmf'
    selected_eebls_snr = np.nan
    if successful_candidates:
        has_ktmf = any(np.isfinite(item[0].get('ktmf_metric', np.nan)) for item in successful_candidates)
        has_eebls = any(np.isfinite(item[0].get('eebls_snr', np.nan)) for item in successful_candidates)
        has_delta_bic = any(np.isfinite(item[0].get('transit_delta_bic', np.nan)) for item in successful_candidates)

        if has_ktmf:
            selection_metric = 'ktmf'
            if pick_comparison_by_eebls_snr:
                selected_entry = min(
                    successful_candidates,
                    key=lambda item: (
                        0 if np.isfinite(item[0].get('ktmf_metric', np.nan)) else 1,
                        -item[0].get('ktmf_metric', np.nan) if np.isfinite(item[0].get('ktmf_metric', np.nan)) else np.inf,
                        0 if np.isfinite(item[0].get('eebls_snr', np.nan)) else 1,
                        -item[0].get('eebls_snr', np.nan) if np.isfinite(item[0].get('eebls_snr', np.nan)) else np.inf,
                        0 if np.isfinite(item[0].get('transit_delta_bic', np.nan)) else 1,
                        -item[0].get('transit_delta_bic', np.nan) if np.isfinite(item[0].get('transit_delta_bic', np.nan)) else np.inf,
                        item[0].get('candidate_order', np.inf),
                    ),
                )
            else:
                selected_entry = min(
                    successful_candidates,
                    key=lambda item: (
                        0 if np.isfinite(item[0].get('ktmf_metric', np.nan)) else 1,
                        -item[0].get('ktmf_metric', np.nan) if np.isfinite(item[0].get('ktmf_metric', np.nan)) else np.inf,
                        0 if np.isfinite(item[0].get('transit_delta_bic', np.nan)) else 1,
                        -item[0].get('transit_delta_bic', np.nan) if np.isfinite(item[0].get('transit_delta_bic', np.nan)) else np.inf,
                        0 if np.isfinite(item[0].get('eebls_snr', np.nan)) else 1,
                        -item[0].get('eebls_snr', np.nan) if np.isfinite(item[0].get('eebls_snr', np.nan)) else np.inf,
                        item[0].get('candidate_order', np.inf),
                    ),
                )
        elif pick_comparison_by_eebls_snr and has_eebls:
            selection_metric = 'eebls_snr'
            selected_entry = min(
                successful_candidates,
                key=lambda item: (
                    0 if np.isfinite(item[0].get('eebls_snr', np.nan)) else 1,
                    -item[0].get('eebls_snr', np.nan) if np.isfinite(item[0].get('eebls_snr', np.nan)) else np.inf,
                    0 if np.isfinite(item[0].get('transit_delta_bic', np.nan)) else 1,
                    -item[0].get('transit_delta_bic', np.nan) if np.isfinite(item[0].get('transit_delta_bic', np.nan)) else np.inf,
                    item[0].get('candidate_order', np.inf),
                ),
            )
        elif has_delta_bic:
            selection_metric = 'transit_delta_bic'
            selected_entry = min(
                successful_candidates,
                key=lambda item: (
                    0 if np.isfinite(item[0].get('transit_delta_bic', np.nan)) else 1,
                    -item[0].get('transit_delta_bic', np.nan) if np.isfinite(item[0].get('transit_delta_bic', np.nan)) else np.inf,
                    0 if np.isfinite(item[0].get('eebls_snr', np.nan)) else 1,
                    -item[0].get('eebls_snr', np.nan) if np.isfinite(item[0].get('eebls_snr', np.nan)) else np.inf,
                    item[0].get('candidate_order', np.inf),
                ),
            )
        elif has_eebls:
            selection_metric = 'eebls_snr'
            selected_entry = min(
                successful_candidates,
                key=lambda item: (
                    0 if np.isfinite(item[0].get('eebls_snr', np.nan)) else 1,
                    -item[0].get('eebls_snr', np.nan) if np.isfinite(item[0].get('eebls_snr', np.nan)) else np.inf,
                    item[0].get('candidate_order', np.inf),
                ),
            )
        else:
            selection_metric = 'comparison_field_rank'
            selected_entry = min(
                successful_candidates,
                key=lambda item: item[0].get('candidate_order', np.inf),
            )

        selected_summary, best_candidate, fit_meta, best_tflux, best_cflux = selected_entry
        best_fit_lc = fit_meta['myfit']
        best_ktmf_metric = selected_summary.get('ktmf_metric', np.nan)
        best_transit_delta_bic = selected_summary.get('transit_delta_bic', np.nan)
        selected_eebls_snr = selected_summary.get('eebls_snr', np.nan)
        best_identity = target_fit_candidate_identity(best_candidate)
        for summary in candidate_summaries:
            summary['selected'] = target_fit_candidate_identity(summary) == best_identity

    return {
        'candidate_jobs': candidate_jobs,
        'evaluated_candidates': evaluated_candidates,
        'shortlist': evaluated_candidates,
        'candidate_summaries': candidate_summaries,
        'best_candidate': best_candidate,
        'best_fit_lc': best_fit_lc,
        'selected_ktmf_metric': best_ktmf_metric,
        'selected_transit_delta_bic': best_transit_delta_bic,
        'selection_metric': selection_metric,
        'selected_eebls_snr': selected_eebls_snr,
        'flux_tar': best_tflux,
        'flux_ref': best_cflux,
    }


def selected_photometry_method_label(photometry_info):
    min_aperture = photometry_info.get('min_aperture')
    min_annulus = photometry_info.get('min_annulus')

    if min_aperture == 0:
        return "PSF photometry"
    if min_aperture is None:
        return "Photometry"

    aper_text = abs(float(min_aperture))
    if min_annulus is None or not np.isfinite(min_annulus):
        return f"Aperture photometry (aper={aper_text:.2f}px)"
    return f"Aperture photometry (aper={aper_text:.2f}px, annulus={float(min_annulus):.2f}px)"


def format_comp_star_position(position):
    if position is None:
        return "x=n/a, y=n/a"

    try:
        x_pos, y_pos = position
        return f"x={float(x_pos):.1f}, y={float(y_pos):.1f}"
    except (TypeError, ValueError):
        return f"coords={position}"


def deduplicate_comparison_star_coords(comp_stars, min_separation_pixels=COMPARISON_STAR_DUPLICATE_DISTANCE_PIXELS):
    """Normalize comparison-star coordinates without merging nearby stars.

    The function name is historical. User-provided comparison-star selections
    are intentional inputs, so nearby stars must remain distinct candidates.
    """
    if comp_stars is None:
        return [], []

    normalized_coords = []
    for coord in comp_stars:
        try:
            x_pos, y_pos = float(coord[0]), float(coord[1])
        except (TypeError, ValueError, IndexError):
            continue

        normalized_coords.append([x_pos, y_pos])

    return normalized_coords, []


def format_comp_star_coverage_text(summary):
    coverage_text = (
        f"{summary['coverage_count']} valid frame(s)"
        f" out of {summary.get('coverage_total_frame_count', 'n/a')} total"
        f"; min_required={summary.get('coverage_min_required_count', 0)}"
    )
    coverage_median = summary.get('coverage_reference_count', np.nan)
    if np.isfinite(coverage_median):
        coverage_text += f"; peer_median={coverage_median:.1f}"
    return coverage_text


def format_eebls_snr(value):
    if value is None:
        return "n/a"

    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return "n/a"

    return "n/a" if not np.isfinite(numeric_value) else f"{numeric_value:.2f}"


def extract_lightcurve_fit_transit_delta_bic(fit):
    if fit is None:
        return np.nan

    transit_qc = getattr(fit, 'transit_qc', None)
    if isinstance(transit_qc, dict):
        value = transit_qc.get('delta_bic', np.nan)
        if np.isfinite(value):
            return float(value)

    value = getattr(fit, 'transit_qc_delta_bic', np.nan)
    return np.nan if not np.isfinite(value) else float(value)


def extract_lightcurve_fit_residual_scatter(fit):
    if fit is None:
        return np.nan

    transit_qc = getattr(fit, 'transit_qc', None)
    if isinstance(transit_qc, dict):
        value = transit_qc.get('residual_scatter', np.nan)
        if np.isfinite(value):
            return float(value)

    value = getattr(fit, 'transit_qc_residual_scatter', np.nan)
    if np.isfinite(value):
        return float(value)

    value = getattr(fit, 'res_stdev', np.nan)
    if np.isfinite(value):
        return float(value)

    residuals = np.asarray(getattr(fit, 'residuals', np.array([])), dtype=float)
    data = np.asarray(getattr(fit, 'data', np.array([])), dtype=float)
    if residuals.shape != data.shape or residuals.size == 0:
        return np.nan

    median_flux = np.nanmedian(data)
    if not np.isfinite(median_flux) or median_flux == 0:
        return np.nan
    return float(np.std(residuals) / median_flux)


def extract_lightcurve_fit_ktmf_metric(fit):
    if fit is None:
        return np.nan

    transit_qc = getattr(fit, 'transit_qc', None)
    if isinstance(transit_qc, dict):
        value = transit_qc.get('ktmf_metric', np.nan)
        if np.isfinite(value):
            return float(value)

    value = getattr(fit, 'transit_qc_ktmf_metric', np.nan)
    return np.nan if not np.isfinite(value) else float(value)


def extract_lightcurve_fit_ktmf_contributions(fit):
    if fit is None:
        return []

    transit_qc = getattr(fit, 'transit_qc', None)
    if isinstance(transit_qc, dict):
        contributions = transit_qc.get('ktmf_contributions')
        if isinstance(contributions, list):
            return contributions

    contributions = getattr(fit, 'transit_qc_ktmf_contributions', None)
    return contributions if isinstance(contributions, list) else []


def compact_comparison_attempt_for_output(attempt):
    if not isinstance(attempt, dict):
        return {}

    return {
        'rank': attempt.get('rank'),
        'comp_index': attempt.get('comp_index'),
        'label': attempt.get('label'),
        'selected': attempt.get('selected'),
        'selection_reason': attempt.get('selection_reason'),
        'ktmf_metric': attempt.get('ktmf_metric'),
        'ktmf_contributions': attempt.get('ktmf_contributions') or [],
        'transit_delta_bic': attempt.get('transit_delta_bic'),
        'eebls_snr': attempt.get('eebls_snr'),
        'residual_scatter': attempt.get('residual_scatter'),
        'fit_point_count': attempt.get('fit_point_count'),
        'transit_qc_status': attempt.get('transit_qc_status'),
        'transit_qc_summary': attempt.get('transit_qc_summary'),
        'rejected_by_transit_qc': attempt.get('rejected_by_transit_qc'),
        'failure_reason': attempt.get('failure_reason'),
    }


def format_transit_delta_bic(value):
    if value is None:
        return "n/a"

    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return "n/a"

    return "n/a" if not np.isfinite(numeric_value) else f"{numeric_value:.2f}"


def format_residual_scatter(value):
    if value is None:
        return "n/a"

    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return "n/a"

    return "n/a" if not np.isfinite(numeric_value) else f"{numeric_value * 100.0:.4f}%"


def format_ktmf_metric(value):
    if value is None:
        return "n/a"

    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return "n/a"

    return "n/a" if not np.isfinite(numeric_value) else f"{numeric_value:.2f}/5.00"


def format_ktmf_contribution(contribution):
    if not isinstance(contribution, dict):
        return "KTMF contribution: unavailable"

    label = contribution.get('label', 'Unknown component')
    detail = contribution.get('detail') or 'n/a'
    available = bool(contribution.get('available'))
    points = float(contribution.get('points', 0.0) or 0.0)
    max_points = float(contribution.get('max_points', 0.0) or 0.0)
    score = contribution.get('score', np.nan)

    if available and np.isfinite(score):
        return (
            f"KTMF contribution: {label} +{points:.2f}/{max_points:.2f} "
            f"(score={score:.2f}; {detail})"
        )

    return f"KTMF contribution: {label} +0.00/0.00 (unavailable; {detail})"


def summarize_lightcurve_fit_assessment(fit):
    if fit is None:
        return None

    fit_method = getattr(fit, 'ns_type', None) or getattr(fit, 'fit_method', None) or 'lm'
    try:
        rprs_retry_count = int(getattr(fit, 'rprs_posterior_refit_count', 0) or 0)
    except (TypeError, ValueError):
        rprs_retry_count = 0
    try:
        ars_retry_count = int(getattr(fit, 'ars_posterior_refit_count', 0) or 0)
    except (TypeError, ValueError):
        ars_retry_count = 0
    try:
        b_retry_count = int(getattr(fit, 'b_posterior_refit_count', 0) or 0)
    except (TypeError, ValueError):
        b_retry_count = 0

    return {
        'fit_method': fit_method,
        'duration_prior_applied': bool(getattr(fit, 'duration_prior_applied', False)),
        'duration_prior_note': getattr(fit, 'duration_prior_note', None),
        'pre_ultranest_transit_coverage_valid': bool(
            getattr(fit, 'pre_ultranest_transit_coverage_valid', False)
        ),
        'pre_ultranest_transit_coverage_status': getattr(
            fit,
            'pre_ultranest_transit_coverage_status',
            None,
        ),
        'pre_ultranest_transit_coverage_chance': getattr(
            fit,
            'pre_ultranest_transit_coverage_chance',
            np.nan,
        ),
        'pre_ultranest_transit_coverage_note': getattr(
            fit,
            'pre_ultranest_transit_coverage_note',
            None,
        ),
        'rprs_posterior_refit_applied': bool(getattr(fit, 'rprs_posterior_refit_applied', False)),
        'rprs_posterior_refit_count': rprs_retry_count,
        'rprs_posterior_refit_note': getattr(fit, 'rprs_posterior_refit_note', None),
        'ars_posterior_refit_applied': bool(getattr(fit, 'ars_posterior_refit_applied', False)),
        'ars_posterior_refit_count': ars_retry_count,
        'ars_posterior_refit_note': getattr(fit, 'ars_posterior_refit_note', None),
        'b_posterior_refit_applied': bool(getattr(fit, 'b_posterior_refit_applied', False)),
        'b_posterior_refit_count': b_retry_count,
        'b_posterior_refit_note': getattr(fit, 'b_posterior_refit_note', None),
        'sparse_posterior_live_point_extension_applied': bool(
            getattr(fit, 'sparse_posterior_live_point_extension_applied', False)
        ),
        'sparse_posterior_live_point_extension_note': getattr(
            fit,
            'sparse_posterior_live_point_extension_note',
            None,
        ),
        'prefit_refinement_applied': bool(getattr(fit, 'prefit_refinement_applied', False)),
        'prefit_refinement_note': getattr(fit, 'prefit_refinement_note', None),
        'oot_baseline_parameter_fit_applied': bool(
            getattr(fit, 'oot_baseline_parameter_fit_applied', False)
        ),
        'oot_baseline_parameter_fit_note': getattr(fit, 'oot_baseline_parameter_fit_note', None),
        'oot_baseline_detrending_applied': bool(getattr(fit, 'oot_baseline_detrending_applied', False)),
        'oot_baseline_detrending_note': getattr(fit, 'oot_baseline_detrending_note', None),
        'airmass_fit_skipped': bool(getattr(fit, 'airmass_fit_skipped', False)),
        'airmass_correction_note': getattr(fit, 'airmass_correction_note', None),
        'nested_tmid_refinement_applied': bool(getattr(fit, 'nested_tmid_refinement_applied', False)),
        'nested_tmid_refinement_note': getattr(fit, 'nested_tmid_refinement_note', None),
        'ultranest_error_fallbacks': getattr(fit, 'ultranest_error_fallbacks', {}) or {},
    }


def best_available_attempt_fit(attempt):
    if not isinstance(attempt, dict):
        return None
    return (
        attempt.get('full_reduction_fit')
        or attempt.get('fit')
        or attempt.get('provisional_fit')
    )


def log_lightcurve_fit_assessment_lines(fit, indent="    "):
    assessment = summarize_lightcurve_fit_assessment(fit)
    if not assessment:
        return

    rprs_retry_status = "not applied"
    if assessment['rprs_posterior_refit_applied']:
        retry_count = assessment['rprs_posterior_refit_count']
        rprs_retry_status = (
            f"applied ({retry_count} refit(s))"
            if retry_count > 0 else
            "applied"
        )
    ars_retry_status = "not applied"
    if assessment['ars_posterior_refit_applied']:
        retry_count = assessment['ars_posterior_refit_count']
        ars_retry_status = (
            f"applied ({retry_count} refit(s))"
            if retry_count > 0 else
            "applied"
        )
    b_retry_status = "not applied"
    if assessment['b_posterior_refit_applied']:
        retry_count = assessment['b_posterior_refit_count']
        b_retry_status = (
            f"applied ({retry_count} refit(s))"
            if retry_count > 0 else
            "applied"
        )
    prefit_status = "applied" if assessment['prefit_refinement_applied'] else "not applied"
    oot_parameter_status = "applied" if assessment['oot_baseline_parameter_fit_applied'] else "not applied"
    oot_status = "applied" if assessment['oot_baseline_detrending_applied'] else "not applied"
    duration_prior_status = "applied" if assessment['duration_prior_applied'] else "not applied"
    sparse_extension_status = (
        "applied"
        if assessment['sparse_posterior_live_point_extension_applied']
        else "not applied"
    )

    log_info(
        f"{indent}fit assessment: fit_method={assessment['fit_method']}, "
        f"duration_prior={duration_prior_status}, "
        f"Rp/R* posterior retry={rprs_retry_status}, "
        f"a/Rs posterior retry={ars_retry_status}, "
        f"impact parameter posterior retry={b_retry_status}, "
        f"sparse posterior extension={sparse_extension_status}, "
        f"prefit_refinement={prefit_status}, "
        f"oot_baseline_parameter_fit={oot_parameter_status}, "
        f"oot_baseline_detrending={oot_status}"
    )
    if assessment.get('duration_prior_note'):
        log_info(f"{indent}Duration prior note: {assessment['duration_prior_note']}")
    if assessment.get('pre_ultranest_transit_coverage_note'):
        status = assessment.get('pre_ultranest_transit_coverage_status') or 'unknown'
        chance = coerce_finite_transit_qc_scalar(
            assessment.get('pre_ultranest_transit_coverage_chance', np.nan)
        )
        chance_text = f", chance~{100.0 * float(chance):.0f}%" if np.isfinite(chance) else ""
        log_info(
            f"{indent}Pre-UltraNest coverage note: status={str(status).upper()}{chance_text}; "
            f"{assessment['pre_ultranest_transit_coverage_note']}"
        )
    if assessment.get('rprs_posterior_refit_note'):
        log_info(f"{indent}Rp/R* posterior retry note: {assessment['rprs_posterior_refit_note']}")
    if assessment.get('ars_posterior_refit_note'):
        log_info(f"{indent}a/Rs posterior retry note: {assessment['ars_posterior_refit_note']}")
    if assessment.get('b_posterior_refit_note'):
        log_info(f"{indent}Impact parameter posterior retry note: {assessment['b_posterior_refit_note']}")
    if assessment.get('sparse_posterior_live_point_extension_note'):
        log_info(
            f"{indent}Sparse posterior live-point extension note: "
            f"{assessment['sparse_posterior_live_point_extension_note']}"
        )
    if assessment.get('prefit_refinement_note'):
        log_info(f"{indent}Prefit refinement note: {assessment['prefit_refinement_note']}")
    if assessment.get('oot_baseline_parameter_fit_note'):
        log_info(f"{indent}OOT baseline parameter-fit note: {assessment['oot_baseline_parameter_fit_note']}")
    if assessment.get('oot_baseline_detrending_note'):
        log_info(f"{indent}OOT baseline detrending note: {assessment['oot_baseline_detrending_note']}")
    if assessment.get('airmass_correction_note'):
        log_info(f"{indent}Airmass correction note: {assessment['airmass_correction_note']}")
    if assessment.get('nested_tmid_refinement_note'):
        log_info(f"{indent}Nested Tmid refinement note: {assessment['nested_tmid_refinement_note']}")
    if assessment.get('ultranest_error_fallbacks'):
        fallback_keys = ", ".join(sorted(assessment['ultranest_error_fallbacks']))
        log_info(
            f"{indent}UltraNest uncertainty fallback note: replaced posterior summary "
            f"error(s) for {fallback_keys} using the sampled log-likelihood neighborhood."
        )


def log_comparison_candidate_evaluation_start(comp_summary, rank, ranked_count, method_label, fit_diagnostics):
    if comp_summary is None:
        return

    label = comp_summary.get('label', f"Comp {comp_summary.get('comp_index', 0) + 1}")
    position_text = format_comp_star_position(comp_summary.get('position'))
    coverage_text = format_comp_star_coverage_text({
        'coverage_count': comp_summary.get('coverage_count', 0),
        'coverage_total_frame_count': comp_summary.get('coverage_total_frame_count', 0),
        'coverage_reference_count': comp_summary.get('coverage_reference_count', np.nan),
        'coverage_min_required_count': comp_summary.get('coverage_min_required_count', 0),
    })
    suitability_score = comp_summary.get('aggregate_score', np.nan)
    suitability_text = "n/a" if not np.isfinite(suitability_score) else f"{suitability_score * 100.0:.4f}%"
    usable_point_count = 0 if fit_diagnostics is None else fit_diagnostics.get('usable_point_count', 0)

    log_info(
        f"\nStarting comparison-star target-fit evaluation for {label} ({position_text}) "
        f"[rank {rank + 1}/{ranked_count}] with {method_label}."
    )
    log_info(
        f"  Candidate inputs: suitability={suitability_text}, coverage={coverage_text}, "
        f"usable_after_filters={usable_point_count}."
    )
    log_info("  Preparing comparison-candidate light curve for the full reduction.")


def log_comparison_candidate_evaluation_result(attempt):
    if not attempt:
        return

    fit = best_available_attempt_fit(attempt)
    fit_method = "n/a"
    assessment = summarize_lightcurve_fit_assessment(fit)
    if assessment:
        fit_method = assessment.get('fit_method', 'n/a')

    qc_status = attempt.get('transit_qc_status')
    qc_text = "n/a" if not qc_status else str(qc_status).upper()
    if attempt.get('fit') is None:
        status_text = "FAILED"
    elif attempt.get('rejected_by_transit_qc', False):
        status_text = f"REJECTED ({qc_text})"
    elif attempt.get('full_reduction_applied', False):
        status_text = f"COMPLETE ({qc_text})"
    else:
        status_text = "PROVISIONAL ONLY"

    reason_text = (
        attempt.get('transit_qc_summary')
        or attempt.get('failure_reason')
        or attempt.get('selection_reason')
        or "completed comparison-star target-fit evaluation."
    )
    log_info(
        f"Completed comparison-star target-fit evaluation for {attempt.get('label', 'comparison candidate')}: "
        f"status={status_text}, fit_method={fit_method}, fit_points={attempt.get('fit_point_count', 0)}, "
        f"transit_qc={qc_text}, transit_delta_bic={format_transit_delta_bic(attempt.get('transit_delta_bic', np.nan))}, "
        f"residual_scatter={format_residual_scatter(attempt.get('residual_scatter', np.nan))}, "
        f"ktmf={format_ktmf_metric(attempt.get('ktmf_metric', np.nan))}, reason={reason_text}"
    )
    parameter_summary = attempt.get('parameter_summary')
    if parameter_summary:
        log_info(f"  parameters: {parameter_summary}")
    log_lightcurve_fit_assessment_lines(fit, indent="  ")
    if attempt.get('final_output_dir'):
        log_info(f"  outputs: {attempt['final_output_dir']}")
    for contribution in attempt.get('ktmf_contributions', []):
        log_info(f"  {format_ktmf_contribution(contribution)}")


def comparison_selection_metric_label(selection_metric):
    if selection_metric == 'first_qc_pass':
        return "First QC PASS"
    if selection_metric == 'promising_partial':
        return "Promising Partial"
    if selection_metric == 'comparison_field_rank':
        return "Comparison-Field Rank"
    if selection_metric == 'ktmf':
        return "KTMF"
    if selection_metric == 'eebls_snr':
        return "EEBLS SNR"
    return "transit-vs-flat Delta BIC"


def should_stop_after_promising_partial_comparison_attempt(attempt):
    if not isinstance(attempt, dict):
        return False
    if attempt.get('fit') is None or not attempt.get('full_reduction_applied', False):
        return False
    if attempt.get('rejected_by_transit_qc', False):
        return False

    status = str(attempt.get('transit_qc_status') or '').strip().lower()
    if status != 'marginal':
        return False

    try:
        coverage_priority = int(attempt.get('preflight_coverage_priority', 4))
    except (TypeError, ValueError):
        coverage_priority = 4
    if coverage_priority > 2:
        return False

    ktmf_metric = coerce_finite_transit_qc_scalar(attempt.get('ktmf_metric', np.nan))
    delta_bic = coerce_finite_transit_qc_scalar(attempt.get('transit_delta_bic', np.nan))
    return (
        np.isfinite(ktmf_metric)
        and ktmf_metric >= PROMISING_PARTIAL_COMPARISON_KTMF_MIN
        and np.isfinite(delta_bic)
        and delta_bic >= TRANSIT_QC_DELTA_BIC_PASS_THRESHOLD
    )


def select_preferred_comparison_attempt(attempts, pick_comparison_by_eebls_snr=True):
    selected_result = None
    selection_metric = 'ktmf'
    if not attempts:
        return selected_result, selection_metric

    has_ktmf = any(np.isfinite(attempt.get('ktmf_metric', np.nan)) for attempt in attempts)
    has_eebls = any(np.isfinite(attempt.get('eebls_snr', np.nan)) for attempt in attempts)
    has_delta_bic = any(np.isfinite(attempt.get('transit_delta_bic', np.nan)) for attempt in attempts)

    if has_ktmf:
        selection_metric = 'ktmf'
        if pick_comparison_by_eebls_snr:
            selected_result = min(
                attempts,
                key=lambda attempt: (
                    0 if np.isfinite(attempt.get('ktmf_metric', np.nan)) else 1,
                    -attempt.get('ktmf_metric', np.nan) if np.isfinite(attempt.get('ktmf_metric', np.nan)) else np.inf,
                    0 if np.isfinite(attempt.get('eebls_snr', np.nan)) else 1,
                    -attempt.get('eebls_snr', np.nan) if np.isfinite(attempt.get('eebls_snr', np.nan)) else np.inf,
                    0 if np.isfinite(attempt.get('transit_delta_bic', np.nan)) else 1,
                    -attempt.get('transit_delta_bic', np.nan) if np.isfinite(attempt.get('transit_delta_bic', np.nan)) else np.inf,
                    attempt.get('rank', np.inf),
                ),
            )
        else:
            selected_result = min(
                attempts,
                key=lambda attempt: (
                    0 if np.isfinite(attempt.get('ktmf_metric', np.nan)) else 1,
                    -attempt.get('ktmf_metric', np.nan) if np.isfinite(attempt.get('ktmf_metric', np.nan)) else np.inf,
                    0 if np.isfinite(attempt.get('transit_delta_bic', np.nan)) else 1,
                    -attempt.get('transit_delta_bic', np.nan) if np.isfinite(attempt.get('transit_delta_bic', np.nan)) else np.inf,
                    0 if np.isfinite(attempt.get('eebls_snr', np.nan)) else 1,
                    -attempt.get('eebls_snr', np.nan) if np.isfinite(attempt.get('eebls_snr', np.nan)) else np.inf,
                    attempt.get('rank', np.inf),
                ),
            )
    elif pick_comparison_by_eebls_snr and has_eebls:
        selection_metric = 'eebls_snr'
        selected_result = min(
            attempts,
            key=lambda attempt: (
                0 if np.isfinite(attempt.get('eebls_snr', np.nan)) else 1,
                -attempt.get('eebls_snr', np.nan) if np.isfinite(attempt.get('eebls_snr', np.nan)) else np.inf,
                0 if np.isfinite(attempt.get('transit_delta_bic', np.nan)) else 1,
                -attempt.get('transit_delta_bic', np.nan) if np.isfinite(attempt.get('transit_delta_bic', np.nan)) else np.inf,
                attempt.get('rank', np.inf),
            ),
        )
    elif has_delta_bic:
        selection_metric = 'transit_delta_bic'
        selected_result = min(
            attempts,
            key=lambda attempt: (
                0 if np.isfinite(attempt.get('transit_delta_bic', np.nan)) else 1,
                -attempt.get('transit_delta_bic', np.nan) if np.isfinite(attempt.get('transit_delta_bic', np.nan)) else np.inf,
                0 if np.isfinite(attempt.get('eebls_snr', np.nan)) else 1,
                -attempt.get('eebls_snr', np.nan) if np.isfinite(attempt.get('eebls_snr', np.nan)) else np.inf,
                attempt.get('rank', np.inf),
            ),
        )
    elif has_eebls:
        selection_metric = 'eebls_snr'
        selected_result = min(
            attempts,
            key=lambda attempt: (
                0 if np.isfinite(attempt.get('eebls_snr', np.nan)) else 1,
                -attempt.get('eebls_snr', np.nan) if np.isfinite(attempt.get('eebls_snr', np.nan)) else np.inf,
                attempt.get('rank', np.inf),
            ),
        )
    else:
        selected_result = min(attempts, key=lambda attempt: attempt.get('rank', np.inf))

    return selected_result, selection_metric


def target_fit_candidate_identity(candidate):
    return (
        candidate.get('method'),
        candidate.get('a'),
        candidate.get('an'),
        candidate.get('comp_index'),
    )


def summarize_target_fit_candidate(candidate, fit_meta, comp_stars):
    comp_index = candidate.get('comp_index')
    fit_meta = {} if fit_meta is None else dict(fit_meta)
    fit_result = fit_meta.get('myfit')
    return {
        'label': "Target-only" if comp_index is None else f"Comp {comp_index + 1}",
        'position': None if comp_index is None else comp_stars[comp_index],
        'selected': False,
        'method_label': comparison_method_label(candidate),
        'prescore': candidate.get('prescore', np.inf),
        'fit': fit_result,
        'coverage_count': candidate.get('coverage_count', 0),
        'coverage_total_frame_count': candidate.get('coverage_total_frame_count', 0),
        'coverage_reference_count': candidate.get('coverage_reference_count', np.nan),
        'coverage_min_required_count': candidate.get('coverage_min_required_count', 0),
        'coverage_rejected': candidate.get('coverage_rejected', False),
        'fit_point_count': fit_meta.get('fit_point_count', 0),
        'eebls_snr': fit_meta.get('eebls_snr', np.nan),
        'transit_delta_bic': fit_meta.get('transit_delta_bic', np.nan),
        'residual_scatter': fit_meta.get('residual_scatter', np.nan),
        'ktmf_metric': fit_meta.get('ktmf_metric', np.nan),
        'ktmf_contributions': fit_meta.get('ktmf_contributions') or [],
        'fit_diagnostics': fit_meta.get('fit_diagnostics') or {},
        'failure_reason': fit_meta.get('failure_reason'),
        'parameter_summary': summarize_lightcurve_fit_parameters(fit_result),
        'comp_index': comp_index,
        'a': candidate.get('a'),
        'an': candidate.get('an'),
        'method': candidate.get('method'),
        'candidate_order': candidate.get('candidate_order'),
        'transit_qc_status': fit_meta.get('transit_qc_status'),
        'transit_qc_summary': fit_meta.get('transit_qc_summary'),
        'rejected_by_transit_qc': fit_meta.get('rejected_by_transit_qc', False),
    }


def log_comparison_calibration_fit_attempt_summaries(attempts, method_label):
    if not attempts:
        return

    log_info("\nComparison-star calibration target-fit diagnostics:")
    log_info(f"Photometry method: {method_label}")

    for attempt in attempts:
        selected_label = " [selected]" if attempt.get('selected') else ""
        position_text = format_comp_star_position(attempt.get('position'))
        diagnostics = attempt.get('fit_diagnostics') or {}
        usable_point_count = diagnostics.get('usable_point_count', 0)
        coverage_text = format_comp_star_coverage_text(attempt)
        suitability_score = attempt.get('aggregate_score', np.inf)
        suitability_text = "n/a" if not np.isfinite(suitability_score) else f"{suitability_score * 100.0:.4f}%"
        eebls_text = format_eebls_snr(attempt.get('eebls_snr', np.nan))
        transit_delta_bic_text = format_transit_delta_bic(attempt.get('transit_delta_bic', np.nan))
        residual_text = format_residual_scatter(attempt.get('residual_scatter', np.nan))
        ktmf_text = format_ktmf_metric(attempt.get('ktmf_metric', np.nan))
        reason_text = attempt.get('selection_reason') or attempt.get(
            'failure_reason',
            "selected: strongest KTMF among the evaluated comparison stars",
        )
        if attempt.get('failed_run_dir'):
            reason_text += f"; archived={attempt['failed_run_dir']}"
        log_info(
            f"  {attempt['label']}{selected_label} ({position_text}): "
            f"suitability={suitability_text}, coverage={coverage_text}, "
            f"usable_after_filters={usable_point_count}, fit_points={attempt.get('fit_point_count', 0)}, "
            f"eebls_snr={eebls_text}, transit_delta_bic={transit_delta_bic_text}, "
            f"residual_scatter={residual_text}, ktmf={ktmf_text}, reason={reason_text}"
        )
        parameter_summary = attempt.get('parameter_summary')
        if parameter_summary:
            log_info(f"    parameters: {parameter_summary}")
        log_lightcurve_fit_assessment_lines(best_available_attempt_fit(attempt), indent="    ")
        for contribution in attempt.get('ktmf_contributions', []):
            log_info(f"    {format_ktmf_contribution(contribution)}")


def log_target_fit_candidate_summaries(candidate_summaries, max_entries=10):
    if not candidate_summaries:
        return

    log_info("\nTarget-fit candidate diagnostics:")

    displayed_summaries = candidate_summaries[:max_entries]
    for summary in displayed_summaries:
        selected_label = " [selected]" if summary.get('selected') else ""
        position_text = format_comp_star_position(summary.get('position'))
        diagnostics = summary.get('fit_diagnostics') or {}
        usable_point_count = diagnostics.get('usable_point_count', 0)
        coverage_text = format_comp_star_coverage_text(summary)
        prescore = summary.get('prescore', np.inf)
        prescore_text = "n/a" if not np.isfinite(prescore) else f"{prescore * 100.0:.4f}%"
        eebls_text = format_eebls_snr(summary.get('eebls_snr', np.nan))
        transit_delta_bic_text = format_transit_delta_bic(summary.get('transit_delta_bic', np.nan))
        residual_text = format_residual_scatter(summary.get('residual_scatter', np.nan))
        ktmf_text = format_ktmf_metric(summary.get('ktmf_metric', np.nan))
        reason_text = summary.get(
            'failure_reason',
            "selected: strongest KTMF in the evaluated candidate set",
        )
        log_info(
            f"  {summary['label']}{selected_label} ({position_text}) with {summary['method_label']}: "
            f"prescore={prescore_text}, coverage={coverage_text}, "
            f"usable_after_filters={usable_point_count}, fit_points={summary.get('fit_point_count', 0)}, "
            f"eebls_snr={eebls_text}, transit_delta_bic={transit_delta_bic_text}, "
            f"residual_scatter={residual_text}, ktmf={ktmf_text}, reason={reason_text}"
        )
        parameter_summary = summary.get('parameter_summary')
        if parameter_summary:
            log_info(f"    parameters: {parameter_summary}")
        for contribution in summary.get('ktmf_contributions', []):
            log_info(f"    {format_ktmf_contribution(contribution)}")

    if len(candidate_summaries) > len(displayed_summaries):
        log_info(
            f"  ... omitted {len(candidate_summaries) - len(displayed_summaries)} additional "
            "target-fit candidate(s); consider increasing the log limit if you need the full list."
        )


def comparison_calibration_selection_reason(summary, best_comp_score):
    if summary.get('selected'):
        return "selected: lowest suitability score among coverage-qualified, sigma-clip-qualified comparison stars for this method"

    if summary.get('coverage_rejected'):
        return (
            "not selected: low coverage "
            f"({summary['coverage_count']} < {summary['coverage_min_required_count']} valid frames)"
        )

    if summary.get('suitability_outlier_rejected'):
        threshold = summary.get('suitability_high_threshold', np.nan)
        if np.isfinite(threshold):
            return (
                "not selected: suitability score was rejected by high-side sigma clipping "
                f"({summary['aggregate_score'] * 100.0:.4f}% > {threshold * 100.0:.4f}%)"
            )
        return "not selected: suitability score was rejected by high-side sigma clipping"

    aggregate_score = summary.get('aggregate_score', np.inf)
    if not np.isfinite(aggregate_score):
        return "not selected: no usable ensemble or pairwise calibration score"

    if np.isfinite(best_comp_score):
        score_gap = aggregate_score - best_comp_score
        if np.isfinite(score_gap) and score_gap > 0:
            return (
                "not selected: suitability score was "
                f"{score_gap * 100.0:.4f}% above the selected comparison star"
            )

    return "not selected: another comparison star ranked better for this photometry method"


def comparison_candidate_fit_selection_reason(summary, photometry_info):
    if summary.get('failure_reason'):
        return summary['failure_reason']

    selection_basis = photometry_info.get('selection_basis', 'target_fit')
    selection_metric = photometry_info.get('selection_metric', 'ktmf')
    selected_comp_num = photometry_info.get('comp_star_num')
    selected_eebls_snr = photometry_info.get('comparison_eebls_snr', np.nan)
    selected_transit_delta_bic = photometry_info.get('comparison_transit_delta_bic', np.nan)
    selected_ktmf_metric = photometry_info.get('comparison_ktmf_metric', np.nan)
    candidate_eebls_snr = summary.get('eebls_snr', np.nan)
    candidate_transit_delta_bic = summary.get('transit_delta_bic', np.nan)
    candidate_ktmf_metric = summary.get('ktmf_metric', np.nan)

    if summary.get('selected'):
        if selection_basis == 'comparison_field':
            return "selected: comparison-field calibration ranked this star best for the chosen photometry method"
        if selection_basis == 'comparison_field_retry':
            return (
                "selected: comparison-field calibration fell back to this star "
                "after better-ranked candidates failed target fitting"
            )
        if selection_basis == 'comparison_field_qc_fallback':
            return (
                "selected: best available comparison-star fit after all completed candidates were rejected "
                "by transit QC"
            )
        if selection_metric == 'first_qc_pass':
            return "selected: first completed comparison-star candidate with PASS transit QC"
        if selection_metric == 'ktmf' and np.isfinite(candidate_ktmf_metric):
            return "selected: highest KTMF in the chosen search"
        if selection_metric == 'eebls_snr' and np.isfinite(candidate_eebls_snr):
            return "selected: highest EEBLS SNR in the chosen search"
        if np.isfinite(candidate_transit_delta_bic):
            return "selected: strongest transit-vs-flat Delta BIC in the chosen search"
        return "selected: strongest transit evidence in the chosen search"

    if selection_basis == 'comparison_field':
        if selected_comp_num is None:
            return "not selected: comparison-field calibration chose a different candidate"
        return f"not selected: comparison-field calibration chose Comp {selected_comp_num}"
    if selection_basis == 'comparison_field_retry':
        if selected_comp_num is None:
            return "not selected: comparison-field fallback chose a different candidate"
        return (
            "not selected: comparison-field fallback chose "
            f"Comp {selected_comp_num} after better-ranked candidate(s) failed target fitting"
        )
    if selection_basis == 'comparison_field_qc_fallback':
        if selected_comp_num is None:
            return "not selected: comparison-field QC fallback chose a different candidate"
        return (
            "not selected: comparison-field QC fallback chose "
            f"Comp {selected_comp_num} as the best available rejected fit"
        )
    if selection_metric == 'first_qc_pass':
        if selected_comp_num is None:
            return "not selected: search stopped after another candidate reached PASS transit QC"
        return (
            "not selected: search stopped after "
            f"Comp {selected_comp_num} reached PASS transit QC"
        )

    if selection_metric == 'eebls_snr' and np.isfinite(selected_eebls_snr):
        if not np.isfinite(candidate_eebls_snr):
            return "not selected: no finite EEBLS SNR was available for this candidate"
        if candidate_eebls_snr < selected_eebls_snr - 1e-12:
            return (
                "not selected: EEBLS SNR was "
                f"{candidate_eebls_snr:.2f} vs {selected_eebls_snr:.2f} for the selected fit"
            )
        if candidate_eebls_snr > selected_eebls_snr + 1e-12:
            return (
                "not selected: this post-selection diagnostic fit has a stronger "
                "EEBLS box signal than the selected fit; the earlier search did not choose it"
            )

    if np.isfinite(candidate_ktmf_metric) and np.isfinite(selected_ktmf_metric):
        if candidate_ktmf_metric < selected_ktmf_metric - 1e-12:
            return (
                "not selected: KTMF was "
                f"{candidate_ktmf_metric:.2f} vs {selected_ktmf_metric:.2f} for the selected fit"
            )
        if candidate_ktmf_metric > selected_ktmf_metric + 1e-12:
            return (
                "not selected: this post-selection diagnostic fit has a higher KTMF than the selected fit; "
                "the earlier target-fit search did not choose it"
            )

    if np.isfinite(candidate_transit_delta_bic) and np.isfinite(selected_transit_delta_bic):
        if candidate_transit_delta_bic < selected_transit_delta_bic - 1e-12:
            return (
                "not selected: transit-vs-flat Delta BIC was "
                f"{candidate_transit_delta_bic:.2f} vs {selected_transit_delta_bic:.2f} for the selected fit"
            )
        if candidate_transit_delta_bic > selected_transit_delta_bic + 1e-12:
            return (
                "not selected: this post-selection diagnostic fit has stronger transit evidence "
                "than the selected fit; the earlier target-fit search did not choose it"
            )

    if selected_comp_num is None:
        return "not selected: another candidate remained preferred in the target-fit search"
    return f"not selected: Comp {selected_comp_num} remained preferred in the target-fit search"


def format_fit_parameter_with_uncertainty(value, error=None, scale=1.0, suffix=""):
    if value is None or not np.isfinite(value):
        return "n/a"

    scaled_value = float(value) * scale
    if error is None or not np.isfinite(error) or error < 0:
        return f"{round_to_2(scaled_value)}{suffix}"

    scaled_error = float(error) * abs(scale)
    return f"{round_to_2(scaled_value, scaled_error)} +/- {round_to_2(scaled_error)}{suffix}"


def summarize_lightcurve_fit_parameters(fit):
    if fit is None or not hasattr(fit, 'parameters'):
        return None

    parameters = getattr(fit, 'parameters', {}) or {}
    errors = getattr(fit, 'errors', {}) or {}
    fit_method = getattr(fit, 'ns_type', 'lm')
    summary_parts = [
        f"fit_method={fit_method}",
        f"Tmid={format_fit_parameter_with_uncertainty(parameters.get('tmid'), errors.get('tmid'))}",
        f"Rp/R*={format_fit_parameter_with_uncertainty(parameters.get('rprs'), errors.get('rprs'))}",
    ]

    rprs = parameters.get('rprs')
    rprs_err = errors.get('rprs')
    depth = None if rprs is None else 100.0 * float(rprs) ** 2
    depth_err = None
    if rprs is not None and rprs_err is not None and np.isfinite(rprs) and np.isfinite(rprs_err):
        depth_err = 200.0 * float(rprs) * float(rprs_err)
    summary_parts.append(f"depth={format_fit_parameter_with_uncertainty(depth, depth_err, suffix='%')}")
    summary_parts.append(f"inc={format_fit_parameter_with_uncertainty(parameters.get('inc'), errors.get('inc'))}")

    if getattr(fit, 'airmass_fit_skipped', False):
        summary_parts.append("airmass=skipped")
    else:
        baseline_key = 'a0' if 'a0' in parameters else 'a1'
        summary_parts.append(
            f"{baseline_key}={format_fit_parameter_with_uncertainty(parameters.get(baseline_key), errors.get(baseline_key))}"
        )
        summary_parts.append(
            f"a2={format_fit_parameter_with_uncertainty(parameters.get('a2'), errors.get('a2'))}"
        )

    return ", ".join(summary_parts)


def log_comparison_candidate_fit_summaries(candidate_fit_summaries, photometry_info):
    if not candidate_fit_summaries:
        return

    selection_basis = photometry_info.get('selection_basis', 'target_fit').replace('_', '-')
    log_info("\nComparison-star lightcurve fit diagnostics:")
    log_info(f"Selection basis: {selection_basis}")
    log_info(
        "Selection metric: "
        f"{comparison_selection_metric_label(photometry_info.get('selection_metric', 'ktmf'))}"
    )

    for summary in candidate_fit_summaries:
        selected_label = " [selected]" if summary.get('selected') else ""
        position_text = format_comp_star_position(summary.get('position'))
        diagnostics = summary.get('fit_diagnostics') or {}
        usable_point_count = diagnostics.get('usable_point_count', 0)
        coverage_text = format_comp_star_coverage_text(summary)
        eebls_text = format_eebls_snr(summary.get('eebls_snr', np.nan))
        transit_delta_bic_text = format_transit_delta_bic(summary.get('transit_delta_bic', np.nan))
        residual_text = format_residual_scatter(summary.get('residual_scatter', np.nan))
        ktmf_text = format_ktmf_metric(summary.get('ktmf_metric', np.nan))
        reason_text = comparison_candidate_fit_selection_reason(summary, photometry_info)
        log_info(
            f"  {summary['label']}{selected_label} ({position_text}): "
            f"coverage={coverage_text}, "
            f"usable_after_filters={usable_point_count}, fit_points={summary['fit_point_count']}, "
            f"eebls_snr={eebls_text}, transit_delta_bic={transit_delta_bic_text}, "
            f"residual_scatter={residual_text}, ktmf={ktmf_text}, reason={reason_text}"
        )
        parameter_summary = summary.get('parameter_summary')
        if parameter_summary:
            log_info(f"    parameters: {parameter_summary}")
        for contribution in summary.get('ktmf_contributions', []):
            log_info(f"    {format_ktmf_contribution(contribution)}")


def fit_lightcurve_to_every_comparison_candidate(times, jd_times, airmass, ld, p_dict, comp_stars,
                                                 psf_data, aper_data, photometry_info,
                                                 plot_time_range=None,
                                                 disable_vertical_flux_normalization=False,
                                                 skip_low_comparison_coverage_rejection=False,
                                                 use_impactparameter_rather_than_inclination_to_fit=True,
                                                 use_eebls_to_initialize_tmid_and_bounds=True):
    if photometry_info.get('best_fit_lc') is None or not comp_stars:
        return []

    use_psf_photometry = photometry_info.get('min_aperture') == 0
    if use_psf_photometry:
        target_flux = 2 * np.pi * psf_data['target'][:, 2] * psf_data['target'][:, 3] * psf_data['target'][:, 4]
        comp_flux_map = {
            f"comp{comp_index + 1}": 2 * np.pi * psf_data[f"comp{comp_index + 1}"][:, 2]
            * psf_data[f"comp{comp_index + 1}"][:, 3]
            * psf_data[f"comp{comp_index + 1}"][:, 4]
            for comp_index in range(len(comp_stars))
        }
    else:
        aperture_index = photometry_info.get('aperture_index')
        annulus_index = photometry_info.get('annulus_index')
        if aperture_index is None or annulus_index is None:
            return []
        target_flux = aper_data['target'][:, aperture_index, annulus_index]
        comp_flux_map = {
            f"comp{comp_index + 1}": aper_data[f"comp{comp_index + 1}"][:, aperture_index, annulus_index]
            for comp_index in range(len(comp_stars))
        }

    candidate_fit_summaries = []
    selected_comp_star_num = photometry_info.get('comp_star_num')
    coverage_summary = comparison_star_coverage_summary(
        comp_flux_map,
        skip_rejection=skip_low_comparison_coverage_rejection,
        validity_mask_func=(
            robust_flux_floor_mask if use_psf_photometry else valid_comparison_frame_mask
        ),
    )

    for comp_index, position in enumerate(comp_stars):
        label = f"Comp {comp_index + 1}"
        ckey = f"comp{comp_index + 1}"
        comp_flux_series = comp_flux_map[ckey]

        if use_psf_photometry:
            fit_mask = robust_target_reference_flux_mask(target_flux, comp_flux_series)
        else:
            fit_mask = valid_comparison_frame_mask(comp_flux_series)
        coverage_count = coverage_summary[ckey]['coverage_count']
        coverage_total_frame_count = coverage_summary[ckey]['coverage_total_frame_count']
        coverage_reference_count = coverage_summary[ckey]['coverage_reference_count']
        coverage_min_required_count = coverage_summary[ckey]['coverage_min_required_count']
        coverage_rejected = coverage_summary[ckey]['coverage_rejected']
        fit_result, target_fit_flux, comp_fit_flux = None, None, None
        fit_diagnostics = {
            'input_point_count': int(times.shape[0]),
            'has_reference_flux': True,
            'relative_flux_point_count': 0,
            'sigma_clip_point_count': 0,
            'usable_point_count': 0,
            'failed_stage': 'coverage',
            'failure_reason': (
                f"only {coverage_count} frame(s) had finite positive comparison flux; need at least 2 to fit."
            ),
        }
        if coverage_rejected:
            fit_diagnostics['failure_reason'] = (
                "comparison candidate rejected after iterative low-coverage clipping "
                f"({coverage_count} < {coverage_min_required_count} valid frame(s); "
                f"peer median={coverage_reference_count:.1f})."
            )
        elif coverage_count > 1:
            fit_diagnostics = diagnose_lightcurve_fit_inputs(
                times[fit_mask],
                target_flux[fit_mask],
                comp_flux_series[fit_mask],
                airmass[fit_mask],
                enforce_relative_flux_max=False,
            )
        if not coverage_rejected and coverage_count > 1 and fit_diagnostics['failure_reason'] is None:
            fit_result, target_fit_flux, comp_fit_flux = fit_lightcurve(
                times[fit_mask],
                target_flux[fit_mask],
                comp_flux_series[fit_mask],
                airmass[fit_mask],
                ld,
                p_dict,
                jd_times[fit_mask],
                allow_mid_transit_range_warning=False,
                disable_vertical_flux_normalization=disable_vertical_flux_normalization,
                final_fit_mode='ns',
                use_impactparameter_rather_than_inclination_to_fit=use_impactparameter_rather_than_inclination_to_fit,
                plot_time_range=plot_time_range,
                use_eebls_to_initialize_tmid_and_bounds=use_eebls_to_initialize_tmid_and_bounds,
                compute_eebls_diagnostics=True,
            )
            fit_diagnostics = ensure_lightcurve_fit_failure_reason(
                fit_diagnostics,
                fit_result,
                failed_stage='nested_fit',
                failure_reason="the nested lightcurve fitter did not converge to a usable solution.",
            )

        fit_point_count = 0 if target_fit_flux is None else int(len(target_fit_flux))
        parameter_summary = summarize_lightcurve_fit_parameters(fit_result)

        candidate_fit_summaries.append({
            'comp_index': comp_index,
            'label': label,
            'position': position,
            'selected': selected_comp_star_num == comp_index + 1,
            'fit': fit_result,
            'eebls_snr': extract_lightcurve_fit_eebls_snr(fit_result),
            'transit_delta_bic': extract_lightcurve_fit_transit_delta_bic(fit_result),
            'residual_scatter': extract_lightcurve_fit_residual_scatter(fit_result),
            'ktmf_metric': extract_lightcurve_fit_ktmf_metric(fit_result),
            'ktmf_contributions': extract_lightcurve_fit_ktmf_contributions(fit_result),
            'coverage_count': coverage_count,
            'coverage_total_frame_count': coverage_total_frame_count,
            'coverage_reference_count': coverage_reference_count,
            'coverage_min_required_count': coverage_min_required_count,
            'coverage_rejected': coverage_rejected,
            'fit_point_count': fit_point_count,
            'fit_diagnostics': fit_diagnostics,
            'failure_reason': fit_diagnostics.get('failure_reason'),
            'fit_method': None if fit_result is None else getattr(fit_result, 'ns_type', 'lm'),
            'parameter_summary': parameter_summary,
        })

    return candidate_fit_summaries


def normalize_flux_series(flux_values, validity_mask_func=valid_comparison_frame_mask):
    flux_values = np.asarray(flux_values, dtype=float)
    normalized = np.full(flux_values.shape, np.nan, dtype=float)
    finite_mask = validity_mask_func(flux_values)
    if np.count_nonzero(finite_mask) < 5:
        return normalized

    flux_median = bn.nanmedian(flux_values[finite_mask])
    if not np.isfinite(flux_median) or flux_median <= 0:
        return normalized

    normalized[finite_mask] = flux_values[finite_mask] / flux_median
    return normalized


def normalized_ratio_series(numerator_flux, denominator_flux):
    numerator_flux = np.asarray(numerator_flux, dtype=float)
    denominator_flux = np.asarray(denominator_flux, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.divide(numerator_flux, denominator_flux)
    ratio[~np.isfinite(ratio)] = np.nan
    ratio[ratio <= 0] = np.nan
    return ratio


def build_normalized_comp_ensemble(normalized_flux_map, exclude_key):
    ensemble_members = [flux for key, flux in normalized_flux_map.items() if key != exclude_key]
    if not ensemble_members:
        return None

    ensemble_stack = np.vstack(ensemble_members)
    valid_mask = np.any(np.isfinite(ensemble_stack), axis=0)
    if not np.any(valid_mask):
        return None

    ensemble = np.full(ensemble_stack.shape[1], np.nan, dtype=float)
    ensemble[valid_mask] = np.nanmedian(ensemble_stack[:, valid_mask], axis=0)
    return ensemble


def comparison_star_coverage_summary(comp_flux_map,
                                     min_fraction=COMPARISON_STAR_MIN_COVERAGE_FRACTION,
                                     min_points=COMPARISON_STAR_MIN_VALID_FRAMES,
                                     skip_rejection=False,
                                     validity_mask_func=valid_comparison_frame_mask):
    comp_keys = list(comp_flux_map.keys())
    if not comp_keys:
        return {}

    coverage_counts = {
        key: int(np.count_nonzero(validity_mask_func(comp_flux_map[key])))
        for key in comp_keys
    }
    total_frame_count = max(np.asarray(comp_flux_map[key]).shape[0] for key in comp_keys)
    effective_min_points = int(min_points) if total_frame_count >= int(min_points) else 0
    active_keys = list(comp_keys)
    coverage_reference_count = float(np.nanmedian([coverage_counts[key] for key in active_keys]))
    coverage_min_required_count = max(effective_min_points, 0)
    coverage_scatter = np.nan

    for _ in range(COMPARISON_STAR_COVERAGE_MAX_ITERS):
        active_counts = np.asarray([coverage_counts[key] for key in active_keys], dtype=float)
        if active_counts.size == 0:
            break

        coverage_reference_count = float(np.nanmedian(active_counts))
        coverage_scatter = robust_scatter(active_counts)
        threshold_candidates = [
            effective_min_points,
            int(np.ceil(float(min_fraction) * coverage_reference_count)),
        ]
        if np.isfinite(coverage_scatter) and coverage_scatter > 0:
            threshold_candidates.append(
                int(np.ceil(coverage_reference_count - COMPARISON_STAR_COVERAGE_SIGMA * coverage_scatter))
            )
        coverage_min_required_count = max(threshold_candidates)

        kept_keys = [key for key in active_keys if coverage_counts[key] >= coverage_min_required_count]
        if len(kept_keys) == len(active_keys):
            break
        active_keys = kept_keys

    coverage_summary = {}
    active_key_set = set(active_keys)
    for key in comp_keys:
        coverage_summary[key] = {
            'coverage_count': coverage_counts[key],
            'coverage_total_frame_count': total_frame_count,
            'coverage_reference_count': coverage_reference_count,
            'coverage_median_count': coverage_reference_count,
            'coverage_scatter': coverage_scatter,
            'coverage_min_required_count': coverage_min_required_count,
            'coverage_rejected': False if skip_rejection else key not in active_key_set,
        }

    return coverage_summary


def apply_comparison_star_suitability_outlier_rejection(
    comp_summaries,
    sigma=COMPARISON_STAR_SUITABILITY_OUTLIER_SIGMA,
    min_candidates=COMPARISON_STAR_SUITABILITY_MIN_CANDIDATES,
    eligible_indices=None,
):
    if eligible_indices is None:
        eligible_indices = [
            index
            for index, summary in enumerate(comp_summaries)
            if (
                not summary.get('coverage_rejected')
                and np.isfinite(summary.get('aggregate_score', np.inf))
            )
        ]
    else:
        eligible_indices = [
            int(index)
            for index in eligible_indices
            if (
                0 <= int(index) < len(comp_summaries)
                and not comp_summaries[int(index)].get('coverage_rejected')
                and np.isfinite(comp_summaries[int(index)].get('aggregate_score', np.inf))
            )
        ]
    clipping_candidate_floor = max(3, int(min_candidates))
    reference_score = np.nan
    scatter = np.nan
    high_threshold = np.nan
    kept_indices = list(eligible_indices)
    rejected_index_set = set()

    if len(eligible_indices) >= clipping_candidate_floor:
        eligible_scores = np.asarray(
            [comp_summaries[index]['aggregate_score'] for index in eligible_indices],
            dtype=float,
        )
        if eligible_scores.size and np.any(np.isfinite(eligible_scores)):
            reference_score = float(np.nanmedian(eligible_scores))
            scatter = robust_scatter(eligible_scores)
            if np.isfinite(scatter) and scatter > 0:
                high_threshold = reference_score + float(sigma) * scatter
                kept_indices = [
                    index for index in eligible_indices
                    if comp_summaries[index]['aggregate_score'] <= high_threshold
                ]
                rejected_index_set = set(eligible_indices) - set(kept_indices)

    for index, summary in enumerate(comp_summaries):
        summary['suitability_outlier_rejected'] = index in rejected_index_set
        summary['suitability_reference_score'] = reference_score
        summary['suitability_scatter'] = scatter
        summary['suitability_high_threshold'] = high_threshold

    return {
        'eligible_indices': eligible_indices,
        'active_indices': kept_indices,
        'rejected_indices': sorted(rejected_index_set),
        'reference_score': reference_score,
        'scatter': scatter,
        'high_threshold': high_threshold,
    }


def comparison_star_image_outlier_summary(
    normalized_flux_map,
    active_keys,
    sigma=COMPARISON_IMAGE_OUTLIER_SIGMA,
    min_active_stars=COMPARISON_IMAGE_OUTLIER_MIN_ACTIVE_STARS,
    min_valid_pairs=COMPARISON_IMAGE_OUTLIER_MIN_VALID_PAIRS,
):
    active_keys = [key for key in active_keys if key in normalized_flux_map]
    series_length = 0
    for key in active_keys:
        flux_values = np.asarray(normalized_flux_map.get(key), dtype=float)
        if flux_values.ndim == 1:
            series_length = flux_values.shape[0]
            break

    keep_mask = np.ones(series_length, dtype=bool)
    summary = {
        'frame_keep_mask': keep_mask,
        'rejected_frame_indices': [],
        'rejected_frame_count': 0,
        'valid_pair_counts': np.zeros(series_length, dtype=int),
        'outlier_pair_counts': np.zeros(series_length, dtype=int),
        'available_pair_count': 0,
        'required_valid_pair_count': 0,
        'sigma': float(sigma),
    }

    if series_length == 0 or len(active_keys) < max(2, int(min_active_stars)):
        return summary

    pairwise_valid_flags = []
    pairwise_outlier_flags = []
    for index, key in enumerate(active_keys):
        numerator_flux = np.asarray(normalized_flux_map[key], dtype=float)
        if numerator_flux.ndim != 1 or numerator_flux.shape[0] != series_length:
            continue

        for other_key in active_keys[index + 1:]:
            denominator_flux = np.asarray(normalized_flux_map[other_key], dtype=float)
            if denominator_flux.ndim != 1 or denominator_flux.shape[0] != series_length:
                continue

            ratio = normalized_ratio_series(numerator_flux, denominator_flux)
            valid_mask = np.isfinite(ratio)
            if np.count_nonzero(valid_mask) < LIGHTCURVE_MIN_VALID_POINTS:
                continue

            center, scatter = sigma_clipped_nanmedian(ratio[valid_mask], sigma=4.0, max_iters=3)
            if not np.isfinite(center):
                center = float(bn.nanmedian(ratio[valid_mask]))
            robust_pair_scatter = robust_scatter(ratio[valid_mask] - center)
            if np.isfinite(robust_pair_scatter) and robust_pair_scatter > 0:
                scatter = robust_pair_scatter
            if not np.isfinite(scatter) or scatter <= 0:
                scatter = robust_scatter(ratio[valid_mask] - center)
            if not np.isfinite(scatter) or scatter <= 0:
                scatter = COMPARISON_IMAGE_OUTLIER_MIN_SCATTER

            outlier_mask = valid_mask & np.greater(np.abs(ratio - center), float(sigma) * scatter)
            pairwise_valid_flags.append(valid_mask)
            pairwise_outlier_flags.append(outlier_mask)

    available_pair_count = len(pairwise_valid_flags)
    required_valid_pair_count = max(int(min_valid_pairs), len(active_keys) - 1)
    summary['available_pair_count'] = available_pair_count
    summary['required_valid_pair_count'] = required_valid_pair_count
    if available_pair_count < required_valid_pair_count:
        return summary

    valid_pair_counts = np.sum(np.vstack(pairwise_valid_flags), axis=0).astype(int)
    outlier_pair_counts = np.sum(np.vstack(pairwise_outlier_flags), axis=0).astype(int)
    rejected_mask = (
        (valid_pair_counts >= required_valid_pair_count)
        & (outlier_pair_counts == valid_pair_counts)
        & (valid_pair_counts > 0)
    )
    keep_mask = ~rejected_mask

    summary.update({
        'frame_keep_mask': keep_mask,
        'rejected_frame_indices': np.flatnonzero(rejected_mask).astype(int).tolist(),
        'rejected_frame_count': int(np.count_nonzero(rejected_mask)),
        'valid_pair_counts': valid_pair_counts,
        'outlier_pair_counts': outlier_pair_counts,
    })
    return summary


def comparison_star_stability_summary(comp_flux_map, airmass, skip_low_coverage_rejection=False,
                                      validity_mask_func=valid_comparison_frame_mask):
    if not comp_flux_map:
        return {
            'pairwise_matrix': np.empty((0, 0), dtype=float),
            'comp_summaries': [],
            'field_score': np.inf,
            'best_comp_index': None,
            'best_comp_score': np.inf,
            'suitability_outlier_rejected_count': 0,
            'suitability_high_threshold': np.nan,
            'suitability_reference_score': np.nan,
            'suitability_scatter': np.nan,
            'field_image_keep_mask': np.ones(0, dtype=bool),
            'image_outlier_rejected_count': 0,
            'image_outlier_sigma': COMPARISON_IMAGE_OUTLIER_SIGMA,
            'image_outlier_required_valid_pairs': 0,
            'image_outlier_available_pairs': 0,
            'image_outlier_valid_pair_counts': np.zeros(0, dtype=int),
            'image_outlier_outlier_pair_counts': np.zeros(0, dtype=int),
        }

    comp_keys = list(comp_flux_map.keys())
    normalized_flux_map = {
        key: normalize_flux_series(comp_flux_map[key], validity_mask_func=validity_mask_func)
        for key in comp_keys
    }
    coverage_summary = comparison_star_coverage_summary(
        comp_flux_map,
        skip_rejection=skip_low_coverage_rejection,
        validity_mask_func=validity_mask_func,
    )
    coverage_qualified_keys = [
        key for key in comp_keys
        if not coverage_summary[key]['coverage_rejected']
    ]

    def build_stability_iteration(active_keys, frame_keep_mask=None):
        active_key_set = set(active_keys)
        if comp_keys:
            reference_shape = normalized_flux_map[comp_keys[0]].shape
        else:
            reference_shape = ()
        if frame_keep_mask is None:
            working_flux_map = normalized_flux_map
        else:
            frame_keep_mask = np.asarray(frame_keep_mask, dtype=bool)
            if frame_keep_mask.shape != reference_shape:
                working_flux_map = normalized_flux_map
            else:
                working_flux_map = {
                    key: np.where(frame_keep_mask, normalized_flux_map[key], np.nan)
                    for key in comp_keys
                }
        active_flux_map = {
            eligible_key: working_flux_map[eligible_key]
            for eligible_key in active_keys
        }
        pairwise_matrix = np.full((len(comp_keys), len(comp_keys)), np.nan, dtype=float)
        comp_summaries = []

        for i, key in enumerate(comp_keys):
            normalized_flux = working_flux_map[key]
            self_score = cheap_lightcurve_prescore(normalized_flux, np.ones(normalized_flux.shape[0]), airmass)
            pairwise_scores = []
            pairwise_series = {}

            for j, other_key in enumerate(comp_keys):
                if i == j or other_key not in active_key_set:
                    continue
                other_flux = working_flux_map[other_key]
                score = cheap_lightcurve_prescore(normalized_flux, other_flux, airmass)
                pairwise_matrix[i, j] = score
                pairwise_series[f"vs {j + 1}"] = normalized_ratio_series(normalized_flux, other_flux)
                if np.isfinite(score):
                    pairwise_scores.append(float(score))

            ensemble_flux = build_normalized_comp_ensemble(active_flux_map, key)
            ensemble_score = np.inf
            ensemble_ratio_series = np.full(normalized_flux.shape, np.nan, dtype=float)
            if ensemble_flux is not None:
                ensemble_score = cheap_lightcurve_prescore(normalized_flux, ensemble_flux, airmass)
                ensemble_ratio_series = normalized_ratio_series(normalized_flux, ensemble_flux)

            if pairwise_scores:
                pairwise_median = float(np.nanmedian(pairwise_scores))
                pairwise_max = float(np.nanmax(pairwise_scores))
                pairwise_upper = float(np.nanpercentile(pairwise_scores, 75))
            else:
                pairwise_median = np.inf
                pairwise_max = np.inf
                pairwise_upper = np.inf

            aggregate_inputs = [score for score in (ensemble_score, pairwise_upper) if np.isfinite(score)]
            aggregate_score = max(aggregate_inputs) if aggregate_inputs else self_score
            if coverage_summary[key]['coverage_rejected']:
                aggregate_score = np.inf

            comp_summaries.append({
                'comp_index': i,
                'key': key,
                'label': f"Comp {i + 1}",
                'pairwise_median_score': pairwise_median,
                'pairwise_max_score': pairwise_max,
                'ensemble_score': float(ensemble_score) if np.isfinite(ensemble_score) else np.inf,
                'self_score': float(self_score) if np.isfinite(self_score) else np.inf,
                'aggregate_score': float(aggregate_score) if np.isfinite(aggregate_score) else np.inf,
                'valid_pair_count': len(pairwise_scores),
                'pairwise_ratio_series': pairwise_series,
                'ensemble_ratio_series': ensemble_ratio_series,
                'coverage_count': coverage_summary[key]['coverage_count'],
                'coverage_total_frame_count': coverage_summary[key]['coverage_total_frame_count'],
                'coverage_reference_count': coverage_summary[key]['coverage_reference_count'],
                'coverage_min_required_count': coverage_summary[key]['coverage_min_required_count'],
                'coverage_rejected': coverage_summary[key]['coverage_rejected'],
                'suitability_outlier_rejected': False,
                'suitability_reference_score': np.nan,
                'suitability_scatter': np.nan,
                'suitability_high_threshold': np.nan,
            })

        return pairwise_matrix, comp_summaries

    active_keys = list(coverage_qualified_keys)
    rejected_outlier_keys = set()
    rejection_metadata = {}
    for _ in range(COMPARISON_STAR_SUITABILITY_MAX_ITERS):
        _, iteration_summaries = build_stability_iteration(active_keys)
        active_indices = [comp_keys.index(key) for key in active_keys]
        outlier_summary = apply_comparison_star_suitability_outlier_rejection(
            iteration_summaries,
            eligible_indices=active_indices,
        )
        newly_rejected_indices = outlier_summary['rejected_indices']
        if not newly_rejected_indices:
            break

        newly_rejected_keys = [comp_keys[index] for index in newly_rejected_indices]
        if len(active_keys) - len(newly_rejected_keys) < 3:
            break

        for index in newly_rejected_indices:
            key = comp_keys[index]
            if key in rejection_metadata:
                continue
            rejection_metadata[key] = {
                'reference_score': outlier_summary['reference_score'],
                'scatter': outlier_summary['scatter'],
                'high_threshold': outlier_summary['high_threshold'],
            }
        rejected_outlier_keys.update(newly_rejected_keys)
        active_keys = [key for key in active_keys if key not in rejected_outlier_keys]

    image_outlier_summary = comparison_star_image_outlier_summary(
        normalized_flux_map,
        active_keys,
    )
    field_image_keep_mask = image_outlier_summary['frame_keep_mask']
    pairwise_matrix, comp_summaries = build_stability_iteration(
        active_keys,
        frame_keep_mask=field_image_keep_mask,
    )
    final_active_scores = np.asarray(
        [
            summary['aggregate_score']
            for summary in comp_summaries
            if summary['key'] in set(active_keys) and np.isfinite(summary['aggregate_score'])
        ],
        dtype=float,
    )
    final_reference_score = np.nan
    final_scatter = np.nan
    final_high_threshold = np.nan
    if final_active_scores.size and np.any(np.isfinite(final_active_scores)):
        final_reference_score = float(np.nanmedian(final_active_scores))
        final_scatter = robust_scatter(final_active_scores)
        if np.isfinite(final_scatter) and final_scatter > 0:
            final_high_threshold = (
                final_reference_score + COMPARISON_STAR_SUITABILITY_OUTLIER_SIGMA * final_scatter
            )

    for summary in comp_summaries:
        rejection_info = rejection_metadata.get(summary['key'])
        if rejection_info is not None:
            summary['suitability_outlier_rejected'] = True
            summary['suitability_reference_score'] = rejection_info['reference_score']
            summary['suitability_scatter'] = rejection_info['scatter']
            summary['suitability_high_threshold'] = rejection_info['high_threshold']
        else:
            summary['suitability_outlier_rejected'] = False
            summary['suitability_reference_score'] = final_reference_score
            summary['suitability_scatter'] = final_scatter
            summary['suitability_high_threshold'] = final_high_threshold

    finite_comp_scores = [
        summary['aggregate_score']
        for summary in comp_summaries
        if (
            np.isfinite(summary['aggregate_score'])
            and not summary.get('suitability_outlier_rejected')
        )
    ]
    field_score = float(np.nanmedian(finite_comp_scores)) if finite_comp_scores else np.inf
    best_comp_index = None
    best_comp_score = np.inf
    for summary in comp_summaries:
        if summary.get('suitability_outlier_rejected'):
            continue
        if summary['aggregate_score'] < best_comp_score:
            best_comp_score = summary['aggregate_score']
            best_comp_index = summary['comp_index']

    return {
        'pairwise_matrix': pairwise_matrix,
        'comp_summaries': comp_summaries,
        'field_score': field_score,
        'best_comp_index': best_comp_index,
        'best_comp_score': best_comp_score,
        'suitability_outlier_rejected_count': len(rejected_outlier_keys),
        'suitability_high_threshold': final_high_threshold,
        'suitability_reference_score': final_reference_score,
        'suitability_scatter': final_scatter,
        'field_image_keep_mask': field_image_keep_mask,
        'image_outlier_rejected_count': image_outlier_summary['rejected_frame_count'],
        'image_outlier_sigma': image_outlier_summary['sigma'],
        'image_outlier_required_valid_pairs': image_outlier_summary['required_valid_pair_count'],
        'image_outlier_available_pairs': image_outlier_summary['available_pair_count'],
        'image_outlier_valid_pair_counts': image_outlier_summary['valid_pair_counts'],
        'image_outlier_outlier_pair_counts': image_outlier_summary['outlier_pair_counts'],
    }


def comparison_field_sort_key(summary):
    return (
        np.inf if summary.get('field_score') is None else summary['field_score'],
        np.inf if summary.get('best_comp_score') is None else summary['best_comp_score'],
    )


def initialize_aperture_data_store(frame_count, aperture_count, annulus_count, comp_star_count):
    aper_shape = (frame_count, aperture_count, annulus_count)
    aper_data = {
        'target': np.full(aper_shape, np.nan, dtype=float),
        'target_bg': np.full(aper_shape, np.nan, dtype=float),
    }

    for comp_idx in range(comp_star_count):
        ckey = f"comp{comp_idx + 1}"
        aper_data[ckey] = np.full(aper_shape, np.nan, dtype=float)
        aper_data[f"{ckey}_bg"] = np.full(aper_shape, np.nan, dtype=float)

    return aper_data


def compute_star_aperture_grid(data, star_index, xc, yc, apertures, annuli, fast_mode=False, sigma_hint=np.nan):
    flux_grid = np.full((len(apertures), len(annuli)), np.nan, dtype=float)
    bg_grid = np.full((len(apertures), len(annuli)), np.nan, dtype=float)

    if np.isnan(xc) or np.isnan(yc):
        return flux_grid, bg_grid

    mask_method = 'center' if fast_mode else 'exact'

    for a_idx, aperture_radius in enumerate(apertures):
        aperture = CircularAperture(positions=[(xc, yc)], r=float(aperture_radius))
        mask = aperture.to_mask(method=mask_method)[0]
        data_cutout = mask.cutout(data)

        mask_area = None
        raw_aperture_sum = None
        if data_cutout is not None:
            mask_area = np.sum(mask.data)
            raw_aperture_sum = (mask.data * data_cutout).sum()

        for an_idx, annulus_width in enumerate(annuli):
            stage_start = perf_counter()
            try:
                if annulus_width > 0:
                    sky_geometry = resolve_sky_annulus_geometry(
                        aperture_radius=float(aperture_radius),
                        annulus_width=float(annulus_width),
                        psf_sigma=sigma_hint,
                    )
                    bgflux, _, _ = skybg_phot(
                        data,
                        star_index,
                        xc,
                        yc,
                        sky_geometry['inner_radius'],
                        sky_geometry['annulus_width'],
                        fast_mode=fast_mode,
                    )
                else:
                    bgflux = 0

                bg_grid[a_idx, an_idx] = bgflux

                if data_cutout is None:
                    flux_grid[a_idx, an_idx] = 0
                else:
                    flux_grid[a_idx, an_idx] = raw_aperture_sum - bgflux * mask_area
            finally:
                _record_photometry_stage_timing('aperPhot', perf_counter() - stage_start)

    return flux_grid, bg_grid


def populate_aperture_data_for_frame(image_data, frame_index, psf_data, comp_star_count, aper_data, apertures, annuli,
                                     fast_aperture_mask, adaptive_apertures=False, fallback_sigma=np.nan):
    frame_sigma = psf_sigma_from_fit(psf_data['target'][frame_index], fallback_sigma=fallback_sigma)
    frame_apertures, frame_annuli = resolve_frame_aperture_radii(
        apertures,
        annuli,
        adaptive_apertures=adaptive_apertures,
        frame_sigma=frame_sigma,
        fallback_sigma=fallback_sigma,
    )

    target_flux, target_bg = compute_star_aperture_grid(
        image_data,
        0,
        psf_data['target'][frame_index, 0],
        psf_data['target'][frame_index, 1],
        frame_apertures,
        frame_annuli,
        fast_mode=fast_aperture_mask,
        sigma_hint=frame_sigma,
    )
    aper_data['target'][frame_index] = target_flux
    aper_data['target_bg'][frame_index] = target_bg

    for comp_idx in range(comp_star_count):
        ckey = f"comp{comp_idx + 1}"
        comp_sigma = psf_sigma_from_fit(psf_data[ckey][frame_index], fallback_sigma=frame_sigma)
        comp_flux, comp_bg = compute_star_aperture_grid(
            image_data,
            comp_idx + 1,
            psf_data[ckey][frame_index, 0],
            psf_data[ckey][frame_index, 1],
            frame_apertures,
            frame_annuli,
            fast_mode=fast_aperture_mask,
            sigma_hint=comp_sigma,
        )
        aper_data[ckey][frame_index] = comp_flux
        aper_data[f"{ckey}_bg"][frame_index] = comp_bg


def load_calibrated_reduction_image(file_name, generalDark, generalBias, generalFlat,
                                    demosaic_fmt, demosaic_out, demosaic_mult,
                                    bad_pixel_reference=None):
    hdul = fits.open(name=file_name, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
    extension = 0
    image_header = hdul[extension].header
    while image_header["NAXIS"] == 0:
        extension += 1
        image_header = hdul[extension].header

    image_data = hdul[extension].data
    hdul.close()

    image_data = apply_cals(image_data, generalDark, generalBias, generalFlat, 1)
    image_data = demosaic_img(image_data, demosaic_fmt, demosaic_out, demosaic_mult, 1)
    image_data = repair_bad_pixels_in_frame(image_data, bad_pixel_reference)
    return image_data


def _refined_sigma_grid(center, lower_bound, upper_bound, half_width, points):
    low = max(lower_bound, center - half_width)
    high = min(upper_bound, center + half_width)
    if high <= low:
        low, high = lower_bound, upper_bound
    return np.linspace(low, high, points)


def auto_tune_aperture_sigma_grid(coarse_apertures_sigma, coarse_annuli_sigma, coarse_aper_data, comp_star_count,
                                  subset_airmass, require_comp_star=True,
                                  skip_low_comparison_coverage_rejection=False):
    best_candidate = None
    best_score = np.inf

    for a_idx, aperture_sigma in enumerate(coarse_apertures_sigma):
        for an_idx, annulus_sigma in enumerate(coarse_annuli_sigma):
            comp_flux_map = {
                f"comp{comp_idx + 1}": coarse_aper_data[f"comp{comp_idx + 1}"][:, a_idx, an_idx]
                for comp_idx in range(comp_star_count)
            }
            field_summary = comparison_star_stability_summary(
                comp_flux_map,
                subset_airmass,
                skip_low_coverage_rejection=skip_low_comparison_coverage_rejection,
            )
            field_score = field_summary['field_score']
            if np.isfinite(field_score) and comparison_field_sort_key(field_summary) < (best_score, np.inf):
                best_score = field_score
                best_candidate = {
                    'aper_sigma': float(aperture_sigma),
                    'annulus_sigma': float(annulus_sigma),
                    'comp_index': field_summary['best_comp_index'],
                }

    if best_candidate is None:
        center_aper_sigma = float(np.median(coarse_apertures_sigma))
        center_annulus_sigma = float(np.median(coarse_annuli_sigma))
        fallback_comp_index = None
        if comp_star_count > 0:
            fallback_comp_index = 0
        elif require_comp_star:
            fallback_comp_index = None
        best_candidate = {
            'aper_sigma': center_aper_sigma,
            'annulus_sigma': center_annulus_sigma,
            'comp_index': fallback_comp_index,
        }
    else:
        center_aper_sigma = best_candidate['aper_sigma']
        center_annulus_sigma = best_candidate['annulus_sigma']

    refined_apertures_sigma = _refined_sigma_grid(
        center_aper_sigma,
        APERTURE_SIGMA_MIN,
        APERTURE_SIGMA_MAX,
        APERTURE_AUTOTUNE_APER_HALF_WIDTH_SIGMA,
        APERTURE_AUTOTUNE_REFINED_APER_POINTS,
    )
    refined_annuli_sigma = _refined_sigma_grid(
        center_annulus_sigma,
        ANNULUS_SIGMA_MIN,
        ANNULUS_SIGMA_MAX,
        APERTURE_AUTOTUNE_ANNULUS_HALF_WIDTH_SIGMA,
        APERTURE_AUTOTUNE_REFINED_ANNULUS_POINTS,
    )

    return refined_apertures_sigma, refined_annuli_sigma, best_candidate, best_score


def comparison_method_label(candidate):
    if candidate.get('method') == 'psf':
        return "PSF photometry"
    return f"Aperture photometry (aper={candidate['aper']:.2f}px, annulus={candidate['annulus']:.2f}px)"


def select_comparison_calibrated_photometry(psf_data, aper_data, apers, annuli, airmass, comp_stars, sigma,
                                           skip_low_comparison_coverage_rejection=False,
                                           use_psf_photometry=True,
                                           use_aperture_photometry=True):
    candidate_summaries = []
    comp_star_count = len(comp_stars)

    if comp_star_count == 0:
        return None

    if use_psf_photometry:
        psf_flux_map = {
            f"comp{comp_idx + 1}": 2 * np.pi * psf_data[f"comp{comp_idx + 1}"][:, 2]
            * psf_data[f"comp{comp_idx + 1}"][:, 3]
            * psf_data[f"comp{comp_idx + 1}"][:, 4]
            for comp_idx in range(comp_star_count)
        }
        psf_summary = comparison_star_stability_summary(
            psf_flux_map,
            airmass,
            skip_low_coverage_rejection=skip_low_comparison_coverage_rejection,
            validity_mask_func=robust_flux_floor_mask,
        )
        psf_summary.update({
            'method': 'psf',
            'a': None,
            'an': None,
            'aper': 0.0,
            'annulus': float(15 * sigma),
        })
        candidate_summaries.append(psf_summary)

    if use_aperture_photometry and aper_data is not None and apers is not None and annuli is not None:
        for a_idx, aperture in enumerate(apers):
            for an_idx, annulus in enumerate(annuli):
                comp_flux_map = {
                    f"comp{comp_idx + 1}": aper_data[f"comp{comp_idx + 1}"][:, a_idx, an_idx]
                    for comp_idx in range(comp_star_count)
                }
                candidate_summary = comparison_star_stability_summary(
                    comp_flux_map,
                    airmass,
                    skip_low_coverage_rejection=skip_low_comparison_coverage_rejection,
                )
                candidate_summary.update({
                    'method': 'aperture',
                    'a': a_idx,
                    'an': an_idx,
                    'aper': float(aperture),
                    'annulus': float(annulus),
                })
                candidate_summaries.append(candidate_summary)

    finite_candidates = [
        candidate for candidate in candidate_summaries
        if np.isfinite(candidate['field_score']) and candidate['best_comp_index'] is not None
    ]
    if not finite_candidates:
        return None

    finite_candidates.sort(key=comparison_field_sort_key)
    best_candidate = finite_candidates[0]
    best_comp_index = best_candidate['best_comp_index']
    method_label = comparison_method_label(best_candidate)
    comp_summaries = []
    for summary in best_candidate['comp_summaries']:
        comp_summary = dict(summary)
        comp_summary['position'] = comp_stars[comp_summary['comp_index']]
        comp_summary['selected'] = comp_summary['comp_index'] == best_comp_index
        comp_summary['selection_reason'] = comparison_calibration_selection_reason(
            comp_summary,
            best_candidate['best_comp_score'],
        )
        comp_summaries.append(comp_summary)

    best_candidate['comp_summaries'] = comp_summaries
    best_candidate['method_label'] = method_label
    return best_candidate


def ranked_comparison_calibration_summaries(comparison_calibration):
    if comparison_calibration is None:
        return []

    ranked_summaries = []
    for summary in comparison_calibration.get('comp_summaries', []):
        aggregate_score = summary.get('aggregate_score', np.inf)
        if summary.get('coverage_rejected'):
            continue
        if summary.get('suitability_outlier_rejected'):
            continue
        if not np.isfinite(aggregate_score):
            continue
        ranked_summaries.append(summary)

    ranked_summaries.sort(
        key=lambda summary: (
            summary.get('aggregate_score', np.inf),
            summary.get('comp_index', np.inf),
        )
    )
    return ranked_summaries


def fit_ranked_comparison_calibration_candidates(times, jd_times, airmass, ld, p_dict, comparison_calibration,
                                                 psf_data, aper_data, target_psf_flux,
                                                 plot_time_range=None,
                                                 disable_vertical_flux_normalization=False,
                                                 detrend_on_outoftransit_baseline=True,
                                                 use_impactparameter_rather_than_inclination_to_fit=True,
                                                 use_eebls_to_initialize_tmid_and_bounds=True,
                                                 pick_comparison_by_eebls_snr=True,
                                                 assess_all_comparisons_before_selecting_best=True,
                                                 exit_at_first_qc_pass_solution=True,
                                                 final_fit_baseline_duration_multiplier=
                                                 FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT,
                                                 use_adaptive_apertures=False,
                                                 adaptive_aperture_values=None,
                                                 adaptive_annulus_values=None,
                                                 fallback_sigma=np.nan,
                                                 run_fast_ultranest_before_final_run=
                                                 FAST_ULTRANEST_BEFORE_FINAL_RUN_DEFAULT,
                                                 save_dir=None,
                                                 planet_name=None,
                                                 observation_date=None):
    ranked_summaries = ranked_comparison_calibration_summaries(comparison_calibration)
    if not ranked_summaries:
        return {
            'ranked_summaries': [],
            'attempts': [],
            'selected_result': None,
        }

    method = comparison_calibration['method']
    method_label = comparison_calibration.get('method_label', method)
    aperture_index = comparison_calibration.get('a')
    annulus_index = comparison_calibration.get('an')
    if method == 'psf':
        target_flux = target_psf_flux
    else:
        target_flux = aper_data['target'][:, aperture_index, annulus_index]

    adaptive_summary = build_comparison_candidate_adaptive_summary(
        comparison_calibration,
        psf_data,
        use_adaptive_apertures=use_adaptive_apertures,
        adaptive_aperture_values=adaptive_aperture_values,
        adaptive_annulus_values=adaptive_annulus_values,
        fallback_sigma=fallback_sigma,
    )
    field_image_keep_mask = np.asarray(
        comparison_calibration.get('field_image_keep_mask', np.ones(times.shape[0], dtype=bool)),
        dtype=bool,
    )
    if field_image_keep_mask.shape != times.shape:
        field_image_keep_mask = np.ones(times.shape[0], dtype=bool)
    field_image_clip_diagnostic = None
    if np.any(~field_image_keep_mask):
        required_pairs = comparison_calibration.get('image_outlier_required_valid_pairs', 0)
        sigma_threshold = comparison_calibration.get('image_outlier_sigma', COMPARISON_IMAGE_OUTLIER_SIGMA)
        field_image_clip_diagnostic = build_time_rejection_diagnostic(
            "Comparison-field image clip",
            times,
            field_image_keep_mask,
            note=(
                "Dropped frames flagged after comparison-star suitability clipping because every "
                f"valid pairwise comparison was more than {sigma_threshold:.2f} sigma from its flat-line median "
                f"(min valid pair count={required_pairs})."
            ),
        )

    preflight_plans = []
    for field_rank, comp_summary in enumerate(ranked_summaries):
        comp_index = comp_summary['comp_index']
        ckey = comp_summary.get('key', f"comp{comp_index + 1}")
        if method == 'psf':
            comp_flux = (
                2 * np.pi * psf_data[ckey][:, 2]
                * psf_data[ckey][:, 3]
                * psf_data[ckey][:, 4]
            )
        else:
            comp_flux = aper_data[ckey][:, aperture_index, annulus_index]

        fit_mask = field_image_keep_mask.copy()
        if method == 'psf':
            fit_mask &= robust_target_reference_flux_mask(target_flux, comp_flux)

        fit_diagnostics = diagnose_lightcurve_fit_inputs(
            times[fit_mask],
            target_flux[fit_mask],
            comp_flux[fit_mask],
            airmass[fit_mask],
            enforce_relative_flux_max=False,
        )
        preflight = build_comparison_candidate_preflight(
            times[fit_mask],
            jd_times[fit_mask],
            airmass[fit_mask],
            ld,
            p_dict,
            target_flux[fit_mask],
            comp_flux[fit_mask],
            adaptive_summary=adaptive_summary,
            use_eebls_to_initialize_tmid_and_bounds=use_eebls_to_initialize_tmid_and_bounds,
        )
        preflight_plans.append({
            'field_rank': field_rank,
            'summary': comp_summary,
            'ckey': ckey,
            'comp_flux': comp_flux,
            'fit_mask': fit_mask,
            'fit_diagnostics': fit_diagnostics,
            'preflight': preflight,
        })

    ranked_preflight_plans = rank_comparison_candidate_preflight_plans(preflight_plans)
    log_comparison_candidate_preflight_order(preflight_plans, ranked_preflight_plans)

    attempts = []
    stopped_after_first_qc_pass = False
    stopped_after_promising_partial = False
    for rank, plan in enumerate(ranked_preflight_plans):
        comp_summary = plan['summary']
        comp_index = comp_summary['comp_index']
        ckey = plan['ckey']
        comp_flux = plan['comp_flux']
        fit_mask = plan['fit_mask']
        fit_diagnostics = plan['fit_diagnostics']
        preflight = plan.get('preflight') or {}
        log_comparison_candidate_evaluation_start(
            comp_summary,
            rank,
            len(ranked_preflight_plans),
            method_label,
            fit_diagnostics,
        )
        log_info(
            "  Full reduction starting. Optional out-of-transit baseline detrending is "
            f"{'enabled' if detrend_on_outoftransit_baseline else 'disabled'}."
        )
        final_reduction = finalize_comparison_candidate_full_reduction(
            times[fit_mask],
            target_flux[fit_mask],
            comp_flux[fit_mask],
            airmass[fit_mask],
            ld,
            p_dict,
            jd_times=jd_times[fit_mask],
            disable_vertical_flux_normalization=disable_vertical_flux_normalization,
            detrend_on_outoftransit_baseline=detrend_on_outoftransit_baseline,
            use_impactparameter_rather_than_inclination_to_fit=
            use_impactparameter_rather_than_inclination_to_fit,
            use_eebls_to_initialize_tmid_and_bounds=use_eebls_to_initialize_tmid_and_bounds,
            plot_time_range=plot_time_range,
            baseline_duration_multiplier=final_fit_baseline_duration_multiplier,
            adaptive_summary=adaptive_summary,
            run_fast_ultranest_before_final_run=run_fast_ultranest_before_final_run,
            precomputed_candidate_series=preflight.get('prepared_series'),
        )
        fit_result = final_reduction.get('fit') if final_reduction.get('applied') else None
        tflux_fit = final_reduction.get('good_target_flux')
        cflux_fit = final_reduction.get('good_comp_flux')
        fit_diagnostics = ensure_lightcurve_fit_failure_reason(
            fit_diagnostics,
            fit_result,
            failed_stage='full_candidate_reduction',
            failure_reason=final_reduction.get(
                'failure_reason',
                "the raw comparison-candidate photometry did not converge to a usable fully reduced solution.",
            ),
        )
        if fit_result is None:
            log_info(
                f"  {comp_summary.get('label', f'Comp {comp_index + 1}')}: "
                "the raw comparison-candidate light curve did not converge to a usable fully reduced fit."
            )
        selection_fit = fit_result
        transit_qc_failure_reason = lightcurve_fit_transit_qc_failure_reason(selection_fit)
        if transit_qc_failure_reason is not None:
            fit_diagnostics = dict(fit_diagnostics)
            fit_diagnostics.update({
                'failed_stage': 'transit_qc',
                'failure_reason': transit_qc_failure_reason,
            })
        if field_image_clip_diagnostic is not None:
            attached_fit_ids = set()
            for fit_candidate in (fit_result, final_reduction.get('fit'), selection_fit):
                if fit_candidate is None or id(fit_candidate) in attached_fit_ids:
                    continue
                attached_fit_ids.add(id(fit_candidate))
                prepend_lightcurve_filter_diagnostic(fit_candidate, field_image_clip_diagnostic)

        attempt = {
            'rank': rank,
            'field_rank': plan.get('field_rank'),
            'comp_index': comp_index,
            'ckey': ckey,
            'label': comp_summary.get('label', f"Comp {comp_index + 1}"),
            'position': comp_summary.get('position'),
            'aggregate_score': comp_summary.get('aggregate_score', np.inf),
            'coverage_count': comp_summary.get('coverage_count', 0),
            'coverage_total_frame_count': comp_summary.get('coverage_total_frame_count', 0),
            'coverage_reference_count': comp_summary.get('coverage_reference_count', np.nan),
            'coverage_min_required_count': comp_summary.get('coverage_min_required_count', 0),
            'coverage_rejected': comp_summary.get('coverage_rejected', False),
            'fit': selection_fit,
            'provisional_fit': None,
            'full_reduction_fit': final_reduction.get('fit'),
            'good_times': final_reduction.get('good_times'),
            'good_flux': final_reduction.get('good_flux'),
            'good_unc': final_reduction.get('good_unc'),
            'good_airmass': final_reduction.get('good_airmass'),
            'good_jd_times': final_reduction.get('good_jd_times'),
            'tflux_fit': tflux_fit,
            'cflux_fit': cflux_fit,
            'source_indices': final_reduction.get('source_indices'),
            'duration_samples': final_reduction.get('duration_samples'),
            'data_highres': final_reduction.get('data_highres'),
            'fit_diagnostics': fit_diagnostics,
            'eebls_snr': extract_lightcurve_fit_eebls_snr(selection_fit),
            'transit_delta_bic': extract_lightcurve_fit_transit_delta_bic(selection_fit),
            'residual_scatter': extract_lightcurve_fit_residual_scatter(selection_fit),
            'ktmf_metric': extract_lightcurve_fit_ktmf_metric(selection_fit),
            'ktmf_contributions': extract_lightcurve_fit_ktmf_contributions(selection_fit),
            'fit_point_count': 0 if tflux_fit is None else int(len(np.asarray(tflux_fit, dtype=float))),
            'failure_reason': fit_diagnostics.get('failure_reason'),
            'parameter_summary': summarize_lightcurve_fit_parameters(selection_fit),
            'transit_qc_status': getattr(selection_fit, 'transit_qc_status', None),
            'transit_qc_summary': getattr(selection_fit, 'transit_qc_summary', None),
            'rejected_by_transit_qc': final_reduction.get('applied', False) and transit_qc_failure_reason is not None,
            'selected': False,
            'selection_reason': None,
            'search_stopped_after_qc_pass': False,
            'search_stopped_after_promising_partial': False,
            'failed_run_dir': None,
            'final_output_dir': None,
            'full_reduction_applied': final_reduction.get('applied', False),
            'full_reduction_note': final_reduction.get('note'),
            'fast_ultranest_binning': final_reduction.get('fast_ultranest_binning'),
            'skip_airmass_fit': final_reduction.get('skip_airmass_fit', False),
            'airmass_skip_note': final_reduction.get('airmass_skip_note'),
            'preflight_coverage_priority': preflight.get('coverage_priority'),
            'preflight_scout_score': (preflight.get('scout') or {}).get('score', np.nan),
        }
        if final_reduction.get('applied') and selection_fit is not None and save_dir is not None:
            final_output_dir = save_comparison_candidate_full_reduction_outputs(
                save_dir,
                None,
                selection_fit,
                p_dict,
                observation_date,
                comp_index,
                comp_coords=comp_summary.get('position'),
                min_aperture=(0 if comparison_calibration['method'] == 'psf' else comparison_calibration.get('aper')),
                min_annulus=comparison_calibration.get('annulus'),
                adaptive_summary=adaptive_summary,
                method_label=comparison_calibration.get('method_label'),
                selection_summary={
                    'ktmf_metric': extract_lightcurve_fit_ktmf_metric(selection_fit),
                    'transit_delta_bic': extract_lightcurve_fit_transit_delta_bic(selection_fit),
                    'eebls_snr': extract_lightcurve_fit_eebls_snr(selection_fit),
                    'transit_qc_status': getattr(selection_fit, 'transit_qc_status', None),
                    'transit_qc_summary': getattr(selection_fit, 'transit_qc_summary', None),
                },
                duration_samples=final_reduction.get('duration_samples'),
                data_highres=final_reduction.get('data_highres'),
            )
            if final_output_dir is not None:
                attempt['final_output_dir'] = str(final_output_dir)
        if transit_qc_failure_reason is not None:
            attempt['selection_reason'] = (
                "rejected: transit detection QC flagged this comparison as a poor transit candidate"
            )
            archive_dir = archive_failed_comparison_fit(
                save_dir,
                planet_name,
                observation_date,
                attempt,
                method_label=comparison_calibration.get('method_label'),
            )
            if archive_dir is not None:
                attempt['failed_run_dir'] = str(archive_dir)
        log_comparison_candidate_evaluation_result(attempt)
        attempts.append(attempt)
        if (
            exit_at_first_qc_pass_solution
            and attempt.get('fit') is not None
            and attempt.get('full_reduction_applied', False)
            and lightcurve_fit_transit_qc_passed(selection_fit)
        ):
            attempt['search_stopped_after_qc_pass'] = True
            stopped_after_first_qc_pass = True
            log_info(
                "Stopping comparison-star candidate search after the first transit-QC PASS fit "
                f"({attempt['label']})."
            )
            break
        if (
            exit_at_first_qc_pass_solution
            and should_stop_after_promising_partial_comparison_attempt(attempt)
        ):
            attempt['search_stopped_after_promising_partial'] = True
            stopped_after_promising_partial = True
            log_info(
                "Stopping comparison-star candidate search after a promising partial-coverage "
                f"MARGINAL fit ({attempt['label']}); proceeding to selected full-resolution confirmation."
            )
            break

    selected_result = None
    completed_attempts = [
        attempt
        for attempt in attempts
        if (
            attempt.get('fit') is not None
            and attempt.get('full_reduction_applied', False)
        )
    ]
    successful_attempts = [
        attempt
        for attempt in completed_attempts
        if not attempt.get('rejected_by_transit_qc', False)
    ]
    first_qc_pass_attempt = next(
        (
            attempt for attempt in attempts
            if attempt.get('search_stopped_after_qc_pass', False)
        ),
        None,
    )
    first_promising_partial_attempt = next(
        (
            attempt for attempt in attempts
            if attempt.get('search_stopped_after_promising_partial', False)
        ),
        None,
    )
    selection_metric = 'ktmf'
    fallback_to_qc_rejected = False
    if first_qc_pass_attempt is not None:
        selected_result = first_qc_pass_attempt
        selection_metric = 'first_qc_pass'
    elif first_promising_partial_attempt is not None:
        selected_result = first_promising_partial_attempt
        selection_metric = 'promising_partial'
    elif successful_attempts:
        selected_result, selection_metric = select_preferred_comparison_attempt(
            successful_attempts,
            pick_comparison_by_eebls_snr=pick_comparison_by_eebls_snr,
        )
    else:
        qc_rejected_attempts = [
            attempt for attempt in completed_attempts
            if attempt.get('rejected_by_transit_qc', False)
        ]
        selected_result, selection_metric = select_preferred_comparison_attempt(
            qc_rejected_attempts,
            pick_comparison_by_eebls_snr=pick_comparison_by_eebls_snr,
        )
        fallback_to_qc_rejected = selected_result is not None

    if selected_result is not None:
        selected_result['selected'] = True
        selected_result['selected_despite_transit_qc'] = fallback_to_qc_rejected
        selected_ktmf_metric = selected_result.get('ktmf_metric', np.nan)
        selected_transit_delta_bic = selected_result.get('transit_delta_bic', np.nan)
        selected_eebls_snr = selected_result.get('eebls_snr', np.nan)

        for attempt in attempts:
            if attempt is selected_result:
                if fallback_to_qc_rejected:
                    attempt['selection_reason'] = (
                        "selected as best available fallback: all completed comparison-star "
                        "target fits were rejected by transit QC"
                    )
                    if selection_metric == 'ktmf' and np.isfinite(selected_ktmf_metric):
                        attempt['selection_reason'] += (
                            f"; this fit had the highest KTMF ({format_ktmf_metric(selected_ktmf_metric)})"
                        )
                    elif selection_metric == 'eebls_snr' and np.isfinite(selected_eebls_snr):
                        attempt['selection_reason'] += (
                            f"; this fit had the highest EEBLS SNR ({selected_eebls_snr:.2f})"
                        )
                    elif np.isfinite(selected_transit_delta_bic):
                        attempt['selection_reason'] += (
                            "; this fit had the strongest transit-vs-flat Delta BIC "
                            f"({format_transit_delta_bic(selected_transit_delta_bic)})"
                        )
                elif attempt.get('search_stopped_after_qc_pass', False):
                    attempt['selection_reason'] = (
                        "selected: first completed comparison-star candidate with PASS transit QC"
                    )
                elif attempt.get('search_stopped_after_promising_partial', False):
                    attempt['selection_reason'] = (
                        "selected: first partial-coverage comparison-star candidate with promising "
                        "MARGINAL transit diagnostics"
                    )
                elif selection_metric == 'ktmf' and np.isfinite(selected_ktmf_metric):
                    attempt['selection_reason'] = (
                        "selected: highest KTMF among the evaluated "
                        "comparison-star calibration candidates"
                    )
                elif selection_metric == 'eebls_snr' and np.isfinite(selected_eebls_snr):
                    attempt['selection_reason'] = (
                        "selected: highest EEBLS SNR among the evaluated "
                        "comparison-star calibration candidates"
                    )
                else:
                    attempt['selection_reason'] = (
                        "selected: strongest transit-vs-flat Delta BIC among the evaluated "
                        "comparison-star calibration candidates"
                    )
                continue
            if (
                attempt.get('fit') is not None
                and attempt.get('full_reduction_applied', False)
                and (not attempt.get('rejected_by_transit_qc', False) or fallback_to_qc_rejected)
            ):
                if selection_metric == 'ktmf' and np.isfinite(selected_ktmf_metric):
                    attempt['selection_reason'] = (
                        "not selected: KTMF "
                        f"{format_ktmf_metric(attempt.get('ktmf_metric', np.nan))} was lower than the selected "
                        f"{format_ktmf_metric(selected_ktmf_metric)}"
                    )
                elif selection_metric == 'first_qc_pass':
                    attempt['selection_reason'] = (
                        "not selected: search stopped after the first comparison-star candidate "
                        "with PASS transit QC"
                    )
                elif selection_metric == 'promising_partial':
                    attempt['selection_reason'] = (
                        "not selected: search stopped after the first partial-coverage comparison-star "
                        "candidate with promising MARGINAL transit diagnostics"
                    )
                elif selection_metric == 'eebls_snr' and np.isfinite(selected_eebls_snr):
                    if np.isfinite(attempt.get('eebls_snr', np.nan)):
                        attempt['selection_reason'] = (
                            "not selected: EEBLS SNR "
                            f"{attempt['eebls_snr']:.2f} was lower than the selected "
                            f"{selected_eebls_snr:.2f}"
                        )
                    else:
                        attempt['selection_reason'] = (
                            "not selected: no finite EEBLS SNR was available for this candidate"
                        )
                else:
                    attempt['selection_reason'] = (
                        "not selected: transit-vs-flat Delta BIC "
                        f"{format_transit_delta_bic(attempt.get('transit_delta_bic', np.nan))} was lower than the selected "
                        f"{format_transit_delta_bic(selected_transit_delta_bic)}"
                    )

        selected_fit = selected_result.get('fit')
        if selected_fit is not None:
            full_resolution_refit_applied = False
            full_resolution_refit = refit_selected_fast_comparison_on_full_lightcurve(
                selected_result,
                p_dict,
                skip_airmass_fit=bool(selected_result.get('skip_airmass_fit', False)),
                airmass_skip_note=selected_result.get('airmass_skip_note'),
                detrend_on_outoftransit_baseline=detrend_on_outoftransit_baseline,
                use_impactparameter_rather_than_inclination_to_fit=
                use_impactparameter_rather_than_inclination_to_fit,
                plot_time_range=plot_time_range,
                duration_prior=build_single_transit_duration_prior(p_dict),
            )
            if full_resolution_refit is not None:
                full_resolution_refit_applied = True
                selected_result['fit'], selected_result['good_flux'], selected_result['good_unc'] = full_resolution_refit
                selected_result['full_reduction_note'] = (
                    "selected candidate rerun on the full-resolution light curve after fast UltraNest search."
                )
            else:
                selected_result['fit'] = extend_selected_comparison_live_points_if_needed(selected_fit)
            selected_result['full_reduction_fit'] = selected_result['fit']
            if (
                full_resolution_refit_applied
                or getattr(selected_result['fit'], 'sparse_posterior_live_point_extension_applied', False)
            ):
                annotate_transit_detection_qc(selected_result['fit'])
                selected_result['eebls_snr'] = extract_lightcurve_fit_eebls_snr(selected_result['fit'])
                selected_result['transit_delta_bic'] = extract_lightcurve_fit_transit_delta_bic(selected_result['fit'])
                selected_result['residual_scatter'] = extract_lightcurve_fit_residual_scatter(selected_result['fit'])
                selected_result['ktmf_metric'] = extract_lightcurve_fit_ktmf_metric(selected_result['fit'])
                selected_result['ktmf_contributions'] = extract_lightcurve_fit_ktmf_contributions(selected_result['fit'])
                selected_result['parameter_summary'] = summarize_lightcurve_fit_parameters(selected_result['fit'])
                selected_result['transit_qc_status'] = getattr(selected_result['fit'], 'transit_qc_status', None)
                selected_result['transit_qc_summary'] = getattr(selected_result['fit'], 'transit_qc_summary', None)
                data_highres, duration_samples = estimate_transit_duration_samples_from_fit(selected_result['fit'])
                selected_result['data_highres'] = data_highres
                selected_result['duration_samples'] = duration_samples
                if save_dir is not None:
                    final_output_dir = save_comparison_candidate_full_reduction_outputs(
                        save_dir,
                        None,
                        selected_result['fit'],
                        p_dict,
                        observation_date,
                        selected_result['comp_index'],
                        comp_coords=selected_result.get('position'),
                        min_aperture=(0 if comparison_calibration['method'] == 'psf' else comparison_calibration.get('aper')),
                        min_annulus=comparison_calibration.get('annulus'),
                        adaptive_summary=adaptive_summary,
                        method_label=comparison_calibration.get('method_label'),
                        selection_summary={
                            'ktmf_metric': selected_result.get('ktmf_metric', np.nan),
                            'transit_delta_bic': selected_result.get('transit_delta_bic', np.nan),
                            'eebls_snr': selected_result.get('eebls_snr', np.nan),
                            'transit_qc_status': selected_result.get('transit_qc_status'),
                            'transit_qc_summary': selected_result.get('transit_qc_summary'),
                        },
                        duration_samples=selected_result.get('duration_samples'),
                        data_highres=selected_result.get('data_highres'),
                    )
                    if final_output_dir is not None:
                        selected_result['final_output_dir'] = str(final_output_dir)

    for attempt in attempts:
        if selected_result is not None and attempt is selected_result:
            continue
        clear_fit_ultranest_resume_state(attempt.get('fit'))
        full_reduction_fit = attempt.get('full_reduction_fit')
        if full_reduction_fit is not attempt.get('fit'):
            clear_fit_ultranest_resume_state(full_reduction_fit)

    return {
        'ranked_summaries': ranked_summaries,
        'attempts': attempts,
        'selected_result': selected_result,
        'selection_metric': selection_metric,
        'stopped_after_first_qc_pass': stopped_after_first_qc_pass,
        'stopped_after_promising_partial': stopped_after_promising_partial,
    }


def parse_args():
    parser = argparse.ArgumentParser(description="Using a JSON initialization file to bypass user inputs for EXOTIC.")
    parser.add_argument('-rt', '--realtime',
                        nargs='?', default=None, type=str, const='',
                        help="Plots transit in real-time while observing with a telescope. "
                             "An initialization file (e.g., inits.json) is optional to use with this command.")
    parser.add_argument('-red', '--reduce',
                        nargs='?', default=None, type=str, const='',
                        help="Performs aperture photometry on FITS files and a reduction on dataset. "
                             "An initialization file (e.g., inits.json) is optional to use with this command.")
    parser.add_argument('-pre', '--prereduced',
                        nargs='?', default=None, type=str, const='',
                        help="Performs a reduction on dataset using the nested sampler only. "
                             "An initialization file (e.g., inits.json) is optional to use with this command.")
    parser.add_argument('-phot', '--photometry',
                        nargs='?', default=None, type=str, const='',
                        help="Performs only aperture photometry on FITS files. "
                             "An initialization file (e.g., inits.json) is optional to use with this command.")
    parser.add_argument('-ov', '--override',
                        action='store_true',
                        help="Adopts all JSON planetary parameters, which will override the NASA Exoplanet Archive. "
                             "Can be used as an additional argument with -rt (--realtime), -red (--reduce), "
                             "-pre (--prereduced), and -phot (--photometry)."
                             "Do not combine with the -nea, --nasaexoarch argument.")
    parser.add_argument('-nea', '--nasaexoarch',
                        action='store_true',
                        help="Adopts all the NASA Exoplanet Archive planetary parameters from "
                             "https://exoplanetarchive.ipac.caltech.edu. "
                             "Can be used as an additional argument with -rt (--realtime), -red (--reduce), "
                             "-pre (--prereduced), and -phot (--photometry)."
                             "Do not combine with the -ov, --override argument.")
    parser.add_argument('--use-nextastro-astrometry',
                        action='store_true',
                        help="Use NextAstro's astrometry service (https://astrometry.nextastro.org/) instead of nova.astrometry.net for plate solving.")
    parser.add_argument('--use-nextastro-variability-server',
                        action='store_true',
                        help="Use NextAstro's variability server for a batch comparison-star variability check. "
                             "If the service returns an error, EXOTIC falls back to individual VSX checks.")
    parser.add_argument('--non-interactive-run',
                        action='store_true',
                        help="Run without interactive prompts for target pixel-coordinate mismatch checks. "
                             "If a mismatch is detected, EXOTIC logs a warning and proceeds with the "
                             "user-provided coordinates.")
    parser.add_argument('--multiprocess-transformations',
                        type=int,
                        default=None,
                        help="Use multiprocessing for frame alignment and fallback image transformations. "
                             "Provide an integer number of processes to use.")
    parser.add_argument('--multiprocess-lightcurve-fits',
                        type=int,
                        default=None,
                        help="Use multiprocessing while evaluating candidate lightcurve fits. "
                             "Provide an integer number of processes to use.")
    return parser.parse_args()


def _main_impl():
    # command line args
    args = parse_args()
    if args.multiprocess_transformations is not None and args.multiprocess_transformations < 1:
        raise ValueError("--multiprocess-transformations requires an integer greater than 0.")
    if args.multiprocess_lightcurve_fits is not None and args.multiprocess_lightcurve_fits < 1:
        raise ValueError("--multiprocess-lightcurve-fits requires an integer greater than 0.")
    configure_windows_multiprocessing_main_spec()
    validate_ultranest_mpi_runtime()

    log.debug("*************************")
    log.debug("EXOTIC reduction log file")
    log.debug("*************************")
    log.debug("Starting ...")
    log.debug("")
    log.debug(f"Python Version: {sys.version}")

    log_info("\n*************************************************************")
    log_info("Welcome to the EXOplanet Transit Interpretation Code (EXOTIC)")
    log_info(f"Version {__version__}")
    log_info("*************************************************************\n")

    # ---INITIALIZATION-------------------------------------------------------

    epw_md5 = None

    userpDict = {'ra': None, 'dec': None, 'pName': None, 'sName': None, 'pPer': None, 'pPerUnc': None,
                 'midT': None, 'midTUnc': None, 'rprs': None, 'rprsUnc': None, 'aRs': None, 'aRsUnc': None,
                 'inc': None, 'incUnc': None, 'omega': None, 'ecc': None, 'teff': None,
                 'teffUncPos': None, 'teffUncNeg': None, 'met': None, 'metUncPos': None, 'metUncNeg': None,
                 'logg': None, 'loggUncPos': None, 'loggUncNeg': None, 'dist': None, 'pm_ra': None, 'pm_dec': None}

    # ---USER INPUTS--------------------------------------------------------------------------
    if isinstance(args.realtime, str):
        reduction_opt = 1
    elif isinstance(args.reduce, str) or isinstance(args.prereduced, str) or isinstance(args.photometry, str):
        reduction_opt = 2
    else:
        reduction_opt = user_input("\nPlease select Reduction method:"
                                   "\n\t1: Real Time Reduction (for analyzing your data while observing)"
                                   "\n\t2: Complete Reduction (for analyzing your data after an observing run)"
                                   "\nEnter 1 or 2: ", type_=int, values=[1, 2])

    if not (args.reduce or args.prereduced or args.realtime or args.photometry):
        file_cmd_opt = user_input("\nPlease select how to input your initial parameters:"
                                  "\n\t1: Command Line"
                                  "\n\t2: Input File (inits.json)"
                                  "\nEnter 1 or 2: ", type_=int, values=[1, 2])
    else:
        file_cmd_opt = 2

    if reduction_opt == 1:
        log_info("\n**************************************************************")
        log_info("Real Time Reduction ('Control + C'  or close the plot to quit)")
        log_info("**************************************************************\n")

        if file_cmd_opt == 2:
            init_opt = 'y'
        else:
            init_opt = 'n'

        inputs_obj = Inputs(init_opt=init_opt)

        if init_opt == 'y':
            init_path, userpDict = inputs_obj.search_init(args.realtime, userpDict)

        exotic_infoDict, userpDict['pName'] = inputs_obj.real_time(userpDict['pName'])

        while True:
            carry_on = user_input(f"\nType continue after the first image has been taken and saved: ", type_=str)
            if carry_on.lower().strip() == 'continue':
                break

        log_info("Real Time Plotting ('Control + C' or close the plot to quit)")
        log_info("\nPlease be patient. It will take at least 15 seconds for the first image to get plotted.")

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(userpDict['pName'])
        ax.set_ylabel('Normalized Flux')
        ax.set_xlabel('Time (JD)')

        anim = FuncAnimation(
            fig,
            realTimeReduce,
            fargs=(userpDict['pName'], userpDict, exotic_infoDict, ax, args.use_nextastro_astrometry, args.multiprocess_transformations),
            interval=15000
        )
        plt.show()

    # ----USER INPUTS----------------------------------------------------------
    else:
        log_info("\n**************************")
        log_info("Complete Reduction Routine")
        log_info("**************************")

        init_path, wcs_file, wcs_header, ra_wcs, dec_wcs, vsp_params, auid = None, None, None, None, None, None, None
        generalDark, generalBias, generalFlat = np.empty(shape=(0, 0)), np.empty(shape=(0, 0)), np.empty(shape=(0, 0))
        demosaic_fmt = None
        demosaic_out = None
        precheck_inputfile_count = None
        post_wcs_inputfile_count = None
        post_pointing_inputfile_count = None
        dropped_wcs_files = []
        dropped_pointing_files = []
        ignore_header_wcs = False
        bad_wcs_threshold_fraction = np.nan
        pointing_rejection_sigma = np.nan
        detect_bad_pixels_before_photometry = None
        bad_pixel_reference = None

        if isinstance(args.reduce, str):
            fitsortext = 1
            init_path = args.reduce
        elif isinstance(args.prereduced, str):
            fitsortext = 2
            init_path = args.prereduced
        elif isinstance(args.photometry, str):
            fitsortext = 1
            init_path = args.photometry
        else:
            fitsortext = user_input("\nPlease select method:"
                                    "\n\t1: Perform Aperture Photometry on FITS files"
                                    "\n\t2: Fit lightcurve for Pre-reduced Data in a .txt format"
                                    "\nEnter 1 or 2: ", type_=int, values=[1, 2])

        if file_cmd_opt == 2:
            init_opt = 'y'
        else:
            init_opt = 'n'

        inputs_obj = Inputs(init_opt=init_opt)

        if init_opt == 'y':
            init_path, userpDict = inputs_obj.search_init(init_path, userpDict)

        if fitsortext == 1:
            exotic_infoDict, userpDict['pName'] = inputs_obj.complete_red(userpDict['pName'])
        else:
            exotic_infoDict, userpDict['pName'] = inputs_obj.prereduced(userpDict['pName'])
            for motion_key in ('dist', 'pm_ra', 'pm_dec'):
                current_motion_value = userpDict.get(motion_key)
                if current_motion_value is None or (isinstance(current_motion_value, str) and not current_motion_value.strip()):
                    header_motion_value = exotic_infoDict.get(motion_key)
                    if header_motion_value is not None:
                        userpDict[motion_key] = header_motion_value
        disable_vertical_flux_normalization = is_vertical_flux_normalization_disabled(
            exotic_infoDict.get('disable_vertical_flux_normalization', False)
        )
        detrend_on_outoftransit_baseline = is_out_of_transit_baseline_detrending_enabled(
            exotic_infoDict.get('detrend_on_outoftransit_baseline', True)
        )
        final_fit_baseline_duration_multiplier = get_final_fit_baseline_duration_multiplier(
            exotic_infoDict.get(
                'final_fit_baseline_duration_multiplier',
                FINAL_FIT_BASELINE_DURATION_MULTIPLIER_DEFAULT,
            )
        )
        use_eebls_tmid_initializer = should_use_eebls_to_initialize_tmid_and_bounds(
            exotic_infoDict.get('use_eebls_to_initialize_tmid_and_bounds', 'y')
        )
        pick_comparison_by_eebls_snr = should_pick_comparison_by_eebls_snr(
            exotic_infoDict.get('pick_comparison_by_eebls_snr', 'y')
        )
        use_impactparameter_rather_than_inclination_to_fit = (
            should_use_impactparameter_rather_than_inclination_to_fit(
                exotic_infoDict.get('use_impactparameter_rather_than_inclination_to_fit', 'y')
            )
        )
        run_fast_ultranest_before_final_run = should_run_fast_ultranest_before_final_run(
            exotic_infoDict.get(
                'run_fast_ultranest_before_final_run',
                FAST_ULTRANEST_BEFORE_FINAL_RUN_DEFAULT,
            )
        )
        ultranest_min_num_live_points = configure_ultranest_min_num_live_points(
            exotic_infoDict.get(
                'ultranest_min_num_live_points',
                ULTRANEST_MIN_NUM_LIVE_POINTS_DEFAULT,
            )
        )
        rprs_search_bound_max = configure_rprs_search_bound_max(
            exotic_infoDict.get(
                'rprs_search_bound_max',
                RPRS_SEARCH_BOUND_MAX_DEFAULT,
            )
        )
        log_info(f"UltraNest minimum live points: {ultranest_min_num_live_points}.")
        log_info(f"Rp/R* maximum search bound: {rprs_search_bound_max:.3f}.")
        if run_fast_ultranest_before_final_run:
            log_info(
                "Fast pre-final UltraNest enabled: comparison-candidate UltraNest search runs "
                f"with at most {FAST_ULTRANEST_MAX_BINNED_POINTS} binned light-curve point(s) "
                f"when more than {FAST_ULTRANEST_MIN_POINTS_TO_BIN} points are available."
            )
        else:
            log_info("Fast pre-final UltraNest disabled per optional_info setting.")
        use_sparse_posterior_live_point_retry = configure_sparse_posterior_live_point_retry(
            exotic_infoDict.get(
                'use_sparse_posterior_live_point_retry',
                SPARSE_POSTERIOR_LIVE_POINT_RETRY_ENABLED_DEFAULT,
            )
        )
        if use_sparse_posterior_live_point_retry:
            if run_fast_ultranest_before_final_run:
                log_info(
                    "Selected comparison-star live-point extension enabled: comparison candidates "
                    "are ranked with fast pre-final UltraNest fits, then the chosen final comparison "
                    "fit reruns on the full-resolution light curve with "
                    f"{SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT}x additional minimum live points."
                )
            else:
                log_info(
                    "Selected comparison-star live-point extension enabled: comparison candidates "
                    "are ranked at the configured UltraNest live-point count, then the chosen final "
                    f"comparison fit continues with {SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT}x "
                    "additional minimum live points using its retained final-pass bounds."
                )
        log_ultranest_mpi_status()

        # Make a temp directory of helpful files
        Path(Path(exotic_infoDict['save']) / "temp").mkdir(exist_ok=True)

        if not args.override:
            nea_obj = NASAExoplanetArchive(planet=userpDict['pName'])
            userpDict['pName'], CandidatePlanetBool, pDict = nea_obj.planet_info()
        else:
            pDict = userpDict
            CandidatePlanetBool = False
        # Seed random number generator (for run to run consistency)
        if exotic_infoDict['random_seed']:
            log_info(f"Setting random number seed to {exotic_infoDict['random_seed']}")
        else:
            exotic_infoDict['random_seed'] = int.from_bytes(hashlib.sha256(f"{pDict['pName']}:{pDict['midT']}".encode()).digest()[0:4], byteorder='little')
            log_info(f"Generated random number seed {exotic_infoDict['random_seed']}")
        np.random.seed(exotic_infoDict['random_seed'])

        if fitsortext == 1:
            # Only do the dark correction if user selects this option
            generalDark = process_dark_frames(exotic_infoDict['darks'])
            generalBias = process_bias_frames(exotic_infoDict['biases'])
            generalFlat = process_flat_frames(exotic_infoDict['flats'], generalBias)

            if exotic_infoDict['demosaic_fmt']:
                demosaic_fmt = exotic_infoDict['demosaic_fmt'].upper()
            if exotic_infoDict['demosaic_out']:
                demosaic_out = exotic_infoDict['demosaic_out']     
            demosaic_mult = calculate_demosaic_mult(demosaic_out)   

        if file_cmd_opt == 2:
            if args.nasaexoarch:
                pass
            elif args.override:
                if type(pDict['ra']) and type(pDict['dec']) is str:
                    pDict['ra'], pDict['dec'] = radec_hours_to_degree(pDict['ra'], pDict['dec'])
            else:
                diff = False

                if type(userpDict['ra']) and type(userpDict['dec']) is str:
                    userpDict['ra'], userpDict['dec'] = radec_hours_to_degree(userpDict['ra'], userpDict['dec'])

                if not CandidatePlanetBool:
                    diff = check_parameters(userpDict, pDict)
                if diff:
                    pDict = get_planetary_parameters(CandidatePlanetBool, userpDict, pdict=pDict)
                else:
                    pDict = userpDict
        else:
            pDict = get_planetary_parameters(CandidatePlanetBool, userpDict, pdict=pDict)

        # check for Nans + Zeros
        for k in pDict:
            if k == 'rprs' and (pDict[k] == 0 or np.isnan(pDict[k])):
                log_info(f"Error: {k} value is 0 or NaN. Please use a non-zero value in inits.json", error=True)
                pDict[k] = 0.8 # instead of 1 since priors on RpRs are 0 to RpRs*1.25
                log_info("EXOTIC will override the Rp/Rs value.")
            if "Unc" in k:
                if not pDict[k]:
                    log_info(f"Warning: {k} uncertainty is 0. Please use a non-zero value in inits.json", warn=True)
                    pDict[k] = 1
                elif pDict[k] == 0 or np.isnan(pDict[k]):
                    log_info(f"Warning: {k} uncertainty is 0. Please use a non-zero value in inits.json", warn=True)
                    pDict[k] = 1
            elif pDict[k] is None:
                log_info(f"Warning: {k} is None. Please use a numeric value in inits.json", warn=True)
                pDict[k] = 0

        if fitsortext == 1:
            log_info("\n**************************"
                     "\nStarting Reduction Process"
                     "\n**************************\n")

            #########################################
            # FLUX DATA EXTRACTION AND MANIPULATION
            #########################################

            airMassList, exptimes = [], []

            plateStatus.initializeFilenames(exotic_infoDict['images'])
            inputfiles = corruption_check(exotic_infoDict['images'])
            # time sort images
            times, jd_times = [], []
            log_info(f"Reading FITS timestamps and converting to BJD_TDB for {len(inputfiles)} frame(s).")
            for file_index, file in enumerate(inputfiles):
                extension = 0
                plateStatus.setCurrentFilename(file)
                header = fits.getheader(filename=file, ext=extension)
                while header['NAXIS'] == 0:
                    extension += 1
                    header = fits.getheader(filename=file, ext=extension)
                obsTime = img_time_bjd_tdb(header, pDict, exotic_infoDict)
                times.append(obsTime)
                plateStatus.setObsTime(obsTime)
                jd_times.append(img_time_jd(header))
                completed = file_index + 1
                if completed == len(inputfiles) or completed % 25 == 0:
                    log_info(f"Timestamp conversion progress: {completed}/{len(inputfiles)}")

            extension = 0
            plateStatus.setCurrentFilename(inputfiles[0])
            header = fits.getheader(filename=inputfiles[0], ext=extension)
            while header['NAXIS'] == 0:
                extension += 1
                header = fits.getheader(filename=inputfiles[0], ext=extension)

            # checks for MOBS data
            if 'CREATOR' in header:
                if 'MicroObservatory' in header['CREATOR'] and 'MOBS' not in exotic_infoDict['second_obs'].upper():
                    if exotic_infoDict['second_obs'].upper() != "":
                        exotic_infoDict['second_obs'] += ",MOBS"
                    else:
                        exotic_infoDict['second_obs'] = "MOBS"
                    exotic_infoDict['filter'] = "MObs CV"
                    exotic_infoDict['elev'] = 1268
                    exotic_infoDict['lat'] = 31.675467
                    exotic_infoDict['long'] = -110.951376
                    exotic_infoDict['pixel_bin'] = "2x2"

            exotic_infoDict.setdefault('observed_filter', exotic_infoDict.get('filter'))
            log_info("Calculating limb-darkening coefficients.")
            ld, ld0, ld1, ld2, ld3 = get_ld_values(pDict, exotic_infoDict)
            log_info("Limb-darkening coefficients ready.")

            # check for EPW_MD5 checksum
            if 'EPW_MD5' in header:
                epw_md5 = header['EPW_MD5']

            si = np.argsort(times)
            times = np.array(times)[si]
            jd_times = np.array(jd_times)[si]
            inputfiles = np.array(inputfiles)[si]
            precheck_inputfile_count = int(len(inputfiles))
            finite_plot_times = times[np.isfinite(times)]
            full_plot_time_range = None
            if finite_plot_times.size:
                full_plot_time_range = (float(np.min(finite_plot_times)), float(np.max(finite_plot_times)))
            ignore_header_wcs = should_ignore_header_wcs(exotic_infoDict.get('ignore_header_wcs'))
            bad_wcs_threshold_fraction = get_bad_wcs_threshold_fraction(
                exotic_infoDict.get('bad_wcs_threshold_percent')
            )
            pointing_rejection_sigma = get_pointing_rejection_sigma(
                exotic_infoDict.get('pointing_rejection_sigma')
            )
            detect_bad_pixels_before_photometry = should_detect_bad_pixels_before_photometry(
                exotic_infoDict.get('detect_bad_pixels_before_photometry', 'y')
            )
            multiprocess_bad_pixel_precheck = get_multiprocess_bad_pixel_precheck_processes(
                exotic_infoDict.get('multiprocess_bad_pixel_precheck', 'n')
            )
            inputfiles, wcs_keep_mask, dropped_wcs_files = filter_sparse_missing_wcs_frames(
                inputfiles,
                ignore_header_wcs=ignore_header_wcs,
                max_missing_fraction=bad_wcs_threshold_fraction,
            )
            if dropped_wcs_files:
                times = times[wcs_keep_mask]
                jd_times = jd_times[wcs_keep_mask]
                plateStatus.initializeFilenames(list(inputfiles))
            post_wcs_inputfile_count = int(len(inputfiles))
            pointing_precheck_inputfiles = np.array(inputfiles, copy=True)
            pointing_reference_file = inputfiles[0] if len(inputfiles) else None
            inputfiles, pointing_keep_mask, dropped_pointing_files, pointing_alignment_transforms = filter_pointing_outlier_frames(
                inputfiles,
                pointing_rejection_sigma=pointing_rejection_sigma,
                ignore_header_wcs=ignore_header_wcs,
                frame_loader=lambda file_name: load_calibrated_reduction_image(
                    file_name,
                    generalDark,
                    generalBias,
                    generalFlat,
                    demosaic_fmt,
                    demosaic_out,
                    demosaic_mult,
                ),
                return_alignment_transforms=True,
                multiprocess_transformations=args.multiprocess_transformations,
                generalDark=generalDark,
                generalBias=generalBias,
                generalFlat=generalFlat,
                demosaic_fmt=demosaic_fmt,
                demosaic_out=demosaic_out,
                demosaic_mult=demosaic_mult,
            )
            if dropped_pointing_files:
                if abort_if_reference_frame_rejected(
                    pointing_reference_file,
                    dropped_pointing_files,
                    ordered_inputfiles=pointing_precheck_inputfiles,
                ):
                    return
                times = times[pointing_keep_mask]
                jd_times = jd_times[pointing_keep_mask]
                finite_plot_times = times[np.isfinite(times)]
                full_plot_time_range = None
                if finite_plot_times.size:
                    full_plot_time_range = (float(np.min(finite_plot_times)), float(np.max(finite_plot_times)))
                plateStatus.initializeFilenames(list(inputfiles))
            post_pointing_inputfile_count = int(len(inputfiles))

            bad_pixel_reference = None
            if detect_bad_pixels_before_photometry:
                log_info(
                    "Bad-pixel precheck enabled: scanning calibrated frames for persistent isolated "
                    "high-count outliers before plate-solve checks and photometry."
                )
                bad_pixel_reference = build_persistent_bad_pixel_map(
                    inputfiles,
                    lambda file_name: load_calibrated_reduction_image(
                        file_name,
                        generalDark,
                        generalBias,
                        generalFlat,
                        demosaic_fmt,
                        demosaic_out,
                        demosaic_mult,
                    ),
                    save_directory=exotic_infoDict['save'],
                    max_processes=multiprocess_bad_pixel_precheck,
                    generalDark=generalDark,
                    generalBias=generalBias,
                    generalFlat=generalFlat,
                    demosaic_fmt=demosaic_fmt,
                    demosaic_out=demosaic_out,
                    demosaic_mult=demosaic_mult,
                )
            else:
                log_info("Bad-pixel precheck disabled per optional_info setting.")

            exotic_UIprevTPX = exotic_infoDict['tar_coords'][0]
            exotic_UIprevTPY = exotic_infoDict['tar_coords'][1]

            # fit target in the first image and use it to determine aperture and annulus range
            inc = 0
            for ifile in inputfiles:
                plateStatus.setCurrentFilename(ifile)
                if bad_pixel_reference is not None:
                    first_image = load_calibrated_reduction_image(
                        ifile,
                        generalDark,
                        generalBias,
                        generalFlat,
                        demosaic_fmt,
                        demosaic_out,
                        demosaic_mult,
                        bad_pixel_reference=bad_pixel_reference,
                    )
                else:
                    first_image = fits.getdata(ifile)
                try:
                    initial_centroid = fit_centroid(first_image, [exotic_UIprevTPX, exotic_UIprevTPY], 0)
                    if np.isnan(initial_centroid[0]):
                        inc += 1
                    else:
                        break
                except Exception:
                    inc += 1
                finally:
                    del first_image

            if inc > 0:
                log_info(f"Skipping first {inc} files - Target star not found")
                inputfiles = inputfiles[inc:]
                times = times[inc:]
                jd_times = jd_times[inc:]
                pointing_alignment_transforms = {}
            plateStatus.setCurrentFilename(inputfiles[0])
            header = get_first_image_header(inputfiles[0])

            # For astrometry hints, prioritize coordinates explicitly provided by the user
            # (from inits.json / CLI) over values scraped from NASA Exoplanet Archive.
            hint_ra = userpDict.get('ra', pDict.get('ra'))
            hint_dec = userpDict.get('dec', pDict.get('dec'))

            wcs_file = check_wcs(inputfiles[0], exotic_infoDict['save'], exotic_infoDict['plate_opt'],
                                 use_nextastro_astrometry=args.use_nextastro_astrometry,
                                 ra=hint_ra, dec=hint_dec, pixel_scale=exotic_infoDict.get('pixel_scale'),
                                 ignore_header_wcs=ignore_header_wcs)
            img_scale_str, img_scale = get_img_scale(header, wcs_file, exotic_infoDict['pixel_scale'])
            plateStatus.initializeComparisonStarCount(len(exotic_infoDict['comp_stars']))
            ra_dec_tar, ra_dec_wcs = None, []
            chart_id, vsp_comp_stars, vsp_list = None, {}, []

            if wcs_file:
                if should_log_plate_solution_path(wcs_file):
                    log_info(f"\n{format_plate_solution_reference(wcs_file)}")
                reference_image = fits.getdata(inputfiles[0])
                wcs_header = get_first_image_header(wcs_file)
                ra_wcs, dec_wcs = get_ra_dec(wcs_header, image_shape=reference_image.shape)

                exotic_UIprevTPX, exotic_UIprevTPY = check_target_pixel_wcs(exotic_UIprevTPX, exotic_UIprevTPY,
                                                                            pDict, ra_wcs, dec_wcs,
                                                                            reference_image,
                                                                            jd_times[0],
                                                                            non_interactive_run=args.non_interactive_run,
                                                                            wcs_header=wcs_header)
                ra_dec_tar = (ra_wcs[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)],
                             dec_wcs[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)])

                auid = vsx_auid(ra_dec_tar[0], ra_dec_tar[1])

                check_for_variable_stars(ra_wcs, dec_wcs, exotic_infoDict['comp_stars'],
                                         use_nextastro_variability_server=args.use_nextastro_variability_server)

                if exotic_infoDict['aavso_comp'] == 'y':
                    vsp_comp_stars, chart_id = vsp_query(wcs_file,[header['NAXIS1'], header['NAXIS2']],
                                                         exotic_infoDict['filter'], img_scale,
                                                         user_comp_stars=exotic_infoDict['comp_stars'],
                                                         user_targ_star = [ exotic_UIprevTPX, exotic_UIprevTPY ])
                    vsp_list = [vsp_star['pos'] for vsp_star in vsp_comp_stars.values()]

                while not exotic_infoDict['comp_stars']:
                    log_info("\nThere are no comparison stars left as all of them were indicated as variable stars."
                             "\nPlease reenter new comparison star coordinates.")
                    exotic_infoDict['comp_stars'] = comparison_star_coords(exotic_infoDict['comp_stars'], False)
                    check_for_variable_stars(ra_wcs, dec_wcs, exotic_infoDict['comp_stars'],
                                             use_nextastro_variability_server=args.use_nextastro_variability_server)

                exotic_infoDict['comp_stars'], duplicate_comp_messages = deduplicate_comparison_star_coords(
                    exotic_infoDict['comp_stars']
                )
                for duplicate_message in duplicate_comp_messages:
                    log_info(duplicate_message)

                # Build RA/Dec for comp after list is finalized (avoid off by one issues, etc
                ra_dec_wcs = build_comp_ra_dec(ra_wcs, dec_wcs, exotic_infoDict['comp_stars'])
                nextastro_field_catalog = None
                try:
                    nextastro_field_catalog = nextastro_photometry_catalog_for_wcs(
                        wcs_file,
                        [header['NAXIS1'], header['NAXIS2']],
                        img_scale,
                        exotic_infoDict['filter'],
                    )
                except Exception as exc:
                    log_info(
                        "\nWarning: NextAstro full-field photometry catalog lookup failed "
                        f"({describe_retry_exception(exc)}). Will try per-comparison catalog lookups.",
                        warn=True,
                    )
                vsp_comp_stars = merge_nextastro_calibration_stars(
                    exotic_infoDict['comp_stars'],
                    ra_dec_wcs,
                    exotic_infoDict['filter'],
                    existing_comp_stars=vsp_comp_stars,
                    field_catalog=nextastro_field_catalog,
                )
                vsp_list = [vsp_star['pos'] for vsp_star in vsp_comp_stars.values()]
                plateStatus.initializeComparisonStarCount(len(exotic_infoDict['comp_stars']))
            else:
                exotic_infoDict['comp_stars'], duplicate_comp_messages = deduplicate_comparison_star_coords(
                    exotic_infoDict['comp_stars']
                )
                for duplicate_message in duplicate_comp_messages:
                    log_info(duplicate_message)
                plateStatus.initializeComparisonStarCount(len(exotic_infoDict['comp_stars']))

            # alloc psf fitting param
            psf_data = {
                # x-cent, y-cent, amplitude, sigma-x, sigma-y, rotation, offset
                'target': np.zeros((len(inputfiles), 7)),  # PSF fit
            }
            tar_comp_dist = {}
            vsp_num = []
            comp_star_count = len(exotic_infoDict['comp_stars'])
            require_comp_star = is_comp_star_required(exotic_infoDict.get('require_comp_star', 'y'))
            target_driven_comp_selection = is_target_driven_comp_selection_enabled(
                exotic_infoDict.get('target_driven_comp_selection', 'n')
            )
            skip_low_comp_coverage_rejection = should_skip_low_comparison_coverage_rejection(
                exotic_infoDict.get('skip_low_comparison_coverage_rejection', 'n')
            )
            if skip_low_comp_coverage_rejection:
                log_info("Skipping low-coverage comparison-star rejection per optional_info setting.")
            fit_every_comparison_candidate = should_fit_lightcurve_to_every_comparison_candidate(
                exotic_infoDict.get('fit_lightcurve_to_every_comparison_candidate', 'n')
            )
            use_deviation_from_expected_transit_in_qc = should_use_deviation_from_expected_transit_in_qc(
                exotic_infoDict.get('use_deviation_from_expected_transit_in_qc', True)
            )
            deviation_from_expected_transit_in_qc_sigma = parse_deviation_from_expected_transit_in_qc_sigma(
                exotic_infoDict.get('deviation_from_expected_transit_in_qc_sigma', 5.0)
            )
            assess_all_comparisons_before_selecting_best = should_assess_all_comparisons_before_selecting_best(
                exotic_infoDict.get('assess_all_comparisons_before_selecting_best', 'y')
            )
            exit_at_first_qc_pass_solution = should_exit_at_first_qc_pass_solution(
                exotic_infoDict.get('exit_at_first_qc_pass_solution', 'y')
            )
            use_psf_photometry = should_use_psf_photometry(
                exotic_infoDict.get('use_psf_photometry', 'y')
            )
            use_aperture_photometry = should_use_aperture_photometry(
                exotic_infoDict.get('use_aperture_photometry', 'y')
            )
            use_adaptive_apertures = is_adaptive_aperture_mode_enabled(
                exotic_infoDict.get('use_adaptive_apertures', False)
            )
            if not use_psf_photometry and not use_aperture_photometry:
                log_info("Error: both PSF and aperture photometry are disabled in optional_info.", error=True)
                return
            if not use_psf_photometry:
                log_info("PSF photometry disabled per optional_info setting.")
            if not use_aperture_photometry:
                log_info("Aperture photometry disabled per optional_info setting.")
            if not use_eebls_tmid_initializer:
                log_info("EEBLS transit initializer disabled per optional_info setting.")
            if not pick_comparison_by_eebls_snr:
                log_info("Comparison-star selection by EEBLS SNR disabled per optional_info setting.")
            if not use_deviation_from_expected_transit_in_qc:
                log_info("Expected-value transit QC deviation checks disabled per optional_info setting.")
            if target_driven_comp_selection:
                log_info(
                    "Warning: target-driven comparison selection is no longer used; "
                    "EXOTIC will run comparison-star calibration followed by full candidate reductions.",
                    warn=True,
                )
            if not assess_all_comparisons_before_selecting_best:
                log_info(
                    "Warning: 'assess_all_comparisons_before_selecting_best' is now ignored; "
                    "comparison-star target-fit search is controlled by 'exit_at_first_qc_pass_solution'.",
                    warn=True,
                )
            if not exit_at_first_qc_pass_solution:
                log_info(
                    "Comparison-star candidate search will evaluate all ranked candidates before selection "
                    "because 'exit_at_first_qc_pass_solution' is disabled."
                )

            pDict['use_deviation_from_expected_transit_in_qc'] = use_deviation_from_expected_transit_in_qc
            pDict['deviation_from_expected_transit_in_qc_sigma'] = deviation_from_expected_transit_in_qc_sigma

            for i, coord in enumerate(exotic_infoDict['comp_stars']):
                ckey = f"comp{i + 1}"
                if coord in vsp_list:
                    vsp_num.append(i)
                psf_data[ckey] = np.zeros((len(inputfiles), 7))
                tar_comp_dist[ckey] = np.zeros(2)

            coarse_tune_frames = 0
            coarse_apertures_sigma = None
            coarse_annuli_sigma = None
            if use_aperture_photometry:
                coarse_tune_frames = min(len(inputfiles), APERTURE_AUTOTUNE_MAX_FRAMES)
                if len(inputfiles) >= APERTURE_AUTOTUNE_MIN_FRAMES:
                    coarse_tune_frames = max(APERTURE_AUTOTUNE_MIN_FRAMES, coarse_tune_frames)
                coarse_apertures_sigma = np.linspace(
                    APERTURE_SIGMA_MIN,
                    APERTURE_SIGMA_MAX,
                    APERTURE_AUTOTUNE_COARSE_APER_POINTS,
                )
                coarse_annuli_sigma = np.linspace(
                    ANNULUS_SIGMA_MIN,
                    ANNULUS_SIGMA_MAX,
                    APERTURE_AUTOTUNE_COARSE_ANNULUS_POINTS,
                )
                log_info(
                    "Automatic aperture tuning enabled: "
                    f"coarse_grid={len(coarse_apertures_sigma)}x{len(coarse_annuli_sigma)}, "
                    f"coarse_frames={coarse_tune_frames}."
                )

            sigma = np.nan
            coarse_aperture_values = None
            coarse_annulus_values = None
            aperture_values = None
            annulus_values = None
            apers = None
            annuli = None
            aperture_grid_tuned = False
            coarse_aper_data = None
            if use_aperture_photometry:
                coarse_aper_data = initialize_aperture_data_store(
                    coarse_tune_frames,
                    len(coarse_apertures_sigma),
                    len(coarse_annuli_sigma),
                    comp_star_count,
                )
            aper_data = None
            coarse_frame_cache = [None] * coarse_tune_frames if use_aperture_photometry else []

            target_and_comp_radec = None
            if ra_dec_tar is not None and ra_dec_wcs:
                target_and_comp_radec = np.array([ra_dec_tar, *ra_dec_wcs], dtype=float)
            target_and_comp_pixels = np.array(
                [[exotic_UIprevTPX, exotic_UIprevTPY], *exotic_infoDict['comp_stars']],
                dtype=float,
            )
            fast_aperture_mask = is_fast_aperture_mask_enabled(exotic_infoDict.get('fast_aperture_mask'))
            if use_aperture_photometry and use_adaptive_apertures:
                log_info("Adaptive aperture scaling enabled: evaluating aperture candidates in PSF sigma units per frame.")

            # open files, calibrate, align, photometry
            reset_transform_timing_stats()
            reset_photometry_timing_stats()
            multiprocess_alignment_results = None
            use_multiprocess_alignment = (
                args.multiprocess_transformations is not None and args.multiprocess_transformations > 0
            )
            comp_alignment_keys = [f"comp{j + 1}" for j in range(comp_star_count)]
            if use_multiprocess_alignment:
                multiprocess_alignment_results = build_multiprocess_alignment_results(
                    inputfiles,
                    args.multiprocess_transformations,
                    target_and_comp_pixels,
                    target_and_comp_radec=target_and_comp_radec,
                    ignore_header_wcs=ignore_header_wcs,
                    generalDark=generalDark,
                    generalBias=generalBias,
                    generalFlat=generalFlat,
                    demosaic_fmt=demosaic_fmt,
                    demosaic_out=demosaic_out,
                    demosaic_mult=demosaic_mult,
                    bad_pixel_reference=bad_pixel_reference,
                    use_fast_centroid_cadence=False,
                    use_adaptive_apertures=use_adaptive_apertures,
                    compute_fallback_transform=True,
                    precomputed_fallback_transforms=pointing_alignment_transforms,
                )
            use_multiprocess_transform_precompute = False
            fallback_transforms = pointing_alignment_transforms
            for i, fileName in enumerate(inputfiles):
                plateStatus.setCurrentFilename(fileName)
                hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False,
                                 ignore_missing_end=True)
                # Final reductions should always use the full centroid fit so the
                # centroid series does not inherit the fast moment-estimator cadence.
                frame_fast_centroid = False
                target_fast_centroid = False

                extension = 0
                image_header = hdul[extension].header
                while image_header["NAXIS"] == 0:
                    extension += 1
                    image_header = hdul[extension].header

                airMassList.append(air_mass(image_header, pDict['ra'], pDict['dec'], exotic_infoDict['lat'], exotic_infoDict['long'],
                                            exotic_infoDict['elev'], jd_times[i]))

                exptimes.append(get_exp_time(image_header))

                # IMAGES
                imageData = hdul[extension].data

                # CALS
                imageData = apply_cals(imageData, generalDark, generalBias, generalFlat, i)
                # Demosaic, if needed
                imageData = demosaic_img(imageData, demosaic_fmt, demosaic_out, demosaic_mult, i)
                imageData = repair_bad_pixels_in_frame(imageData, bad_pixel_reference)

                if i == 0:
                    firstImage = np.copy(imageData)

                if multiprocess_alignment_results is not None:
                    apply_parallel_alignment_result(
                        multiprocess_alignment_results[i],
                        i,
                        psf_data,
                        tar_comp_dist,
                        comp_alignment_keys,
                    )
                else:
                    use_wcs_alignment = False
                    if not ignore_header_wcs:
                        try:
                            wcs_hdr = search_wcs_from_header(image_header)
                            use_wcs_alignment = wcs_hdr.is_celestial
                        except Exception:
                            use_wcs_alignment = False

                    if use_wcs_alignment:
                        try:
                            pix_x = pix_y = None
                            if target_and_comp_radec is not None:
                                pix_x, pix_y = wcs_hdr.world_to_pixel_values(
                                    target_and_comp_radec[:, 0],
                                    target_and_comp_radec[:, 1],
                                )
                                pix_x = np.asarray(pix_x, dtype=float).reshape(-1)
                                pix_y = np.asarray(pix_y, dtype=float).reshape(-1)

                            if i == 0:
                                tx, ty = exotic_UIprevTPX, exotic_UIprevTPY
                            else:
                                tx, ty = pix_x[0], pix_y[0]

                            projected_coords = np.array(
                                [[tx, ty], *np.column_stack((pix_x[1:], pix_y[1:]))] if pix_x is not None else [[tx, ty]],
                                dtype=float,
                            )
                            projected_off_frame = any_projected_coord_out_of_frame(projected_coords, imageData.shape)
                            target_seed = choose_centroid_seed_position(
                                [tx, ty],
                                None if i == 0 else psf_data['target'][i - 1],
                            )

                            psf_data['target'][i] = fit_centroid_or_warn_out_of_frame(
                                imageData,
                                target_seed,
                                0,
                                fast_mode=target_fast_centroid,
                            )

                            # TODO: Add check for flux on target/comp stars relative to others in the field
                            # in case of cloudy data, large changes, etc.
                            current_comp_psf_rows = {}
                            previous_comp_psf_rows = {}
                            for j in range(len(exotic_infoDict['comp_stars'])):
                                ckey = f"comp{j + 1}"

                                cx, cy = pix_x[j + 1], pix_y[j + 1]
                                comp_seed = choose_centroid_seed_position(
                                    [cx, cy],
                                    None if i == 0 else psf_data[ckey][i - 1],
                                )
                                psf_data[ckey][i] = fit_centroid_or_warn_out_of_frame(
                                    imageData,
                                    comp_seed,
                                    j + 1,
                                    fast_mode=frame_fast_centroid,
                                )

                                current_comp_psf_rows[ckey] = psf_data[ckey][i]
                                if i != 0:
                                    previous_comp_psf_rows[ckey] = psf_data[ckey][i - 1]
                                else:
                                    tar_comp_dist[ckey][0] = abs(int(psf_data[ckey][0][0]) - int(psf_data['target'][0][0]))
                                    tar_comp_dist[ckey][1] = abs(int(psf_data[ckey][0][1]) - int(psf_data['target'][0][1]))

                            wcs_alignment_decision = should_keep_header_wcs_alignment(
                                projected_off_frame,
                                i,
                                psf_data['target'][i],
                                previous_target_psf_row=None if i == 0 else psf_data['target'][i - 1],
                                comp_psf_rows=current_comp_psf_rows,
                                previous_comp_psf_rows=previous_comp_psf_rows,
                                expected_offsets=tar_comp_dist,
                            )
                            use_wcs_alignment = wcs_alignment_decision['use_wcs_alignment']
                        except Exception:
                            use_wcs_alignment = False

                    log_alignment_progress(
                        i,
                        len(inputfiles),
                        fileName,
                        use_multiprocess_transform_precompute,
                    )

                    if not use_wcs_alignment:
                        cached_tform = fallback_transforms.get(str(fileName)) if fallback_transforms else None
                        if cached_tform is not None:
                            tform = cached_tform
                        elif i == 0:
                            tform = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])
                        else:
                            tform = transformation(imageData, fileName, reference_image=firstImage)

                        transformed_coords = np.asarray(tform(target_and_comp_pixels), dtype=float)
                        tx, ty = transformed_coords[0]
                        target_seed = choose_centroid_seed_position(
                            [tx, ty],
                            None if i == 0 else psf_data['target'][i - 1],
                        )
                        psf_data['target'][i] = fit_centroid_or_warn_out_of_frame(
                            imageData,
                            target_seed,
                            0,
                            fast_mode=target_fast_centroid,
                        )

                        for j, coord in enumerate(exotic_infoDict['comp_stars']):
                            ckey = f"comp{j + 1}"

                            cx, cy = transformed_coords[j + 1]
                            comp_seed = choose_centroid_seed_position(
                                [cx, cy],
                                None if i == 0 else psf_data[ckey][i - 1],
                            )
                            psf_data[ckey][i] = fit_centroid_or_warn_out_of_frame(
                                imageData,
                                comp_seed,
                                j + 1,
                                fast_mode=frame_fast_centroid,
                            )

                            if i == 0:
                                tar_comp_dist[ckey][0] = abs(int(psf_data[ckey][0][0]) - int(psf_data['target'][0][0]))
                                tar_comp_dist[ckey][1] = abs(int(psf_data[ckey][0][1]) - int(psf_data['target'][0][1]))

                # aperture photometry
                if use_aperture_photometry and i == 0:
                    sigma = psf_sigma_from_fit(psf_data['target'][0])
                    if not np.isfinite(sigma) or sigma <= 0:
                        log_info("Warning: Initial PSF sigma is invalid; using sigma=1.0 for automatic aperture tuning.", warn=True)
                        sigma = 1.0
                    if use_adaptive_apertures:
                        coarse_aperture_values = coarse_apertures_sigma
                        coarse_annulus_values = coarse_annuli_sigma
                    else:
                        coarse_aperture_values = coarse_apertures_sigma * sigma
                        coarse_annulus_values = coarse_annuli_sigma * sigma

                if use_aperture_photometry and i < coarse_tune_frames:
                    coarse_frame_cache[i] = np.array(imageData, copy=True)
                    populate_aperture_data_for_frame(
                        imageData,
                        i,
                        psf_data,
                        comp_star_count,
                        coarse_aper_data,
                        coarse_aperture_values,
                        coarse_annulus_values,
                        fast_aperture_mask,
                        adaptive_apertures=use_adaptive_apertures,
                        fallback_sigma=sigma,
                    )

                    if i == coarse_tune_frames - 1:
                        subset_airmass = np.asarray(airMassList[:coarse_tune_frames], dtype=float)
                        refined_apertures_sigma, refined_annuli_sigma, best_coarse_candidate, best_coarse_score = auto_tune_aperture_sigma_grid(
                            coarse_apertures_sigma,
                            coarse_annuli_sigma,
                            coarse_aper_data,
                            comp_star_count,
                            subset_airmass,
                            require_comp_star=require_comp_star,
                            skip_low_comparison_coverage_rejection=skip_low_comp_coverage_rejection,
                        )
                        if use_adaptive_apertures:
                            aperture_values = refined_apertures_sigma
                            annulus_values = refined_annuli_sigma
                        else:
                            aperture_values = refined_apertures_sigma * sigma
                            annulus_values = refined_annuli_sigma * sigma
                        apers = refined_apertures_sigma * sigma
                        annuli = refined_annuli_sigma * sigma
                        aper_data = initialize_aperture_data_store(len(inputfiles), len(apers), len(annuli), comp_star_count)
                        aperture_grid_tuned = True

                        best_comp_label = "none"
                        if best_coarse_candidate['comp_index'] is not None:
                            best_comp_label = str(best_coarse_candidate['comp_index'] + 1)
                        score_text = "n/a" if not np.isfinite(best_coarse_score) else f"{best_coarse_score:.5f}"
                        log_info(
                            "Auto-tuned aperture grid: "
                            f"coarse_best=(aper={best_coarse_candidate['aper_sigma']:.2f} sigma, "
                            f"annulus={best_coarse_candidate['annulus_sigma']:.2f} sigma, comp={best_comp_label}, score={score_text}), "
                            f"refined_grid={len(refined_apertures_sigma)}x{len(refined_annuli_sigma)}."
                        )

                        log_info(f"Backfilling refined aperture photometry for the first {coarse_tune_frames} frame(s).")
                        for backfill_idx in range(coarse_tune_frames):
                            backfill_image = coarse_frame_cache[backfill_idx]
                            loaded_from_disk = False
                            if backfill_image is None:
                                backfill_image = load_calibrated_reduction_image(
                                    inputfiles[backfill_idx],
                                    generalDark,
                                    generalBias,
                                    generalFlat,
                                    demosaic_fmt,
                                    demosaic_out,
                                    demosaic_mult,
                                    bad_pixel_reference=bad_pixel_reference,
                                )
                                loaded_from_disk = True
                            try:
                                populate_aperture_data_for_frame(
                                    backfill_image,
                                    backfill_idx,
                                    psf_data,
                                    comp_star_count,
                                    aper_data,
                                    aperture_values,
                                    annulus_values,
                                    fast_aperture_mask,
                                    adaptive_apertures=use_adaptive_apertures,
                                    fallback_sigma=sigma,
                                )
                            finally:
                                if loaded_from_disk:
                                    del backfill_image
                                coarse_frame_cache[backfill_idx] = None
                elif use_aperture_photometry:
                    if not aperture_grid_tuned:
                        # Defensive fallback for unexpected control flow.
                        aperture_values = coarse_aperture_values
                        annulus_values = coarse_annulus_values
                        if use_adaptive_apertures:
                            apers = coarse_apertures_sigma * sigma
                            annuli = coarse_annuli_sigma * sigma
                        else:
                            apers = coarse_aperture_values
                            annuli = coarse_annulus_values
                        aper_data = initialize_aperture_data_store(len(inputfiles), len(apers), len(annuli), comp_star_count)
                        aperture_grid_tuned = True

                    populate_aperture_data_for_frame(
                        imageData,
                        i,
                        psf_data,
                        comp_star_count,
                        aper_data,
                        aperture_values,
                        annulus_values,
                        fast_aperture_mask,
                        adaptive_apertures=use_adaptive_apertures,
                        fallback_sigma=sigma,
                    )

                # close file + delete from memory
                hdul.close()
                del hdul
                del imageData

            log_transform_timing_stats('Transformation timing summary (full reduction)')
            log_photometry_timing_stats('Photometry timing summary (full reduction)')
            log_reduction_timing_overview('Reduction timing overview (full reduction)')

            # filter bad images
            badmask = np.isnan(psf_data["target"][:, 0]) | (psf_data["target"][:, 0] == 0)
            if aper_data is not None:
                badmask = badmask | (aper_data["target"][:, 0, 0] == 0) | np.isnan(aper_data["target"][:, 0, 0])
            goodmask = ~badmask
            global_frame_filter_diagnostic = build_time_rejection_diagnostic(
                "Target centroid/aperture validity filter",
                times,
                goodmask,
                note="Dropped frames before photometry selection because the target centroid or target aperture photometry was invalid.",
            )
            if global_frame_filter_diagnostic is not None and global_frame_filter_diagnostic['dropped_point_count'] > 0:
                log_lightcurve_filter_diagnostics(
                    [global_frame_filter_diagnostic],
                    header="Global reduction frame rejections before photometry selection",
                )
            if np.sum(goodmask) == 0:
                log_info("No images to fit...check reference image for alignment (first image of sequence)")

            # convert to numpy arrays - strip all bad data
            times = times[goodmask]
            jd_times = jd_times[goodmask]
            airmass = np.array(airMassList)[goodmask]
            psf_data["target"] = psf_data["target"][goodmask]
            if aper_data is not None:
                aper_data["target"] = aper_data["target"][goodmask]
                aper_data["target_bg"] = aper_data["target_bg"][goodmask]
            for j in range(len(exotic_infoDict['comp_stars'])):
                ckey = f"comp{j + 1}"
                psf_data[ckey] = psf_data[ckey][goodmask]
                if aper_data is not None:
                    aper_data[ckey] = aper_data[ckey][goodmask]
                    aper_data[f"{ckey}_bg"] = aper_data[f"{ckey}_bg"][goodmask]

            sigma_display = representative_psf_sigma(psf_data['target'], fallback_sigma=sigma)
            if not np.isfinite(sigma_display) or sigma_display <= 0:
                sigma_display = 1.0

            if aperture_values is not None and annulus_values is not None:
                if use_adaptive_apertures:
                    apers = np.asarray(aperture_values, dtype=float) * sigma_display
                    annuli = np.asarray(annulus_values, dtype=float) * sigma_display
                else:
                    apers = np.asarray(aperture_values, dtype=float)
                    annuli = np.asarray(annulus_values, dtype=float)

            exotic_infoDict['exposure'] = exp_time_med(exptimes)

            # save PSF data to disk using savetxt
            np.savetxt(Path(exotic_infoDict['save']) / "temp" / "psf_data_target.txt", psf_data["target"],
                          header="#x_centroid, y_centroid, amplitude, sigma_x, sigma_y, rotation offset",
                          fmt="%.6f")
                        # x-cent, y-cent, amplitude, sigma-x, sigma-y, rotation, offset

            # PSF flux
            tFlux = 2 * np.pi * psf_data['target'][:, 2] * psf_data['target'][:, 3] * psf_data['target'][:, 4]

            ref_flux = {}
            if vsp_list:
                ref_flux = {i: None for i in vsp_num}

            flux_values = {
                'flux_tar': None,
                'flux_ref': None,
                'flux_unc_tar': None,
                'flux_unc_ref': None
            }

            centroid_positions = {
                'x_targ': None,
                'y_targ': None,
                'x_ref': None,
                'y_ref': None
            }

            photometry_info = {
                'best_fit_lc': None,
                'comp_star_num': None,
                'comp_star_coords': None,
                'min_aperture': None,
                'min_annulus': None,
                'aperture_index': None,
                'annulus_index': None,
                'adaptive_summary': None,
                'calibration_field_score': np.inf,
                'selection_basis': 'target_fit',
                'selection_metric': 'ktmf',
                'comparison_ktmf_metric': np.nan,
                'comparison_eebls_snr': np.nan,
                'comparison_transit_delta_bic': np.nan,
            }

            comparison_calibration = None
            comparison_calibration = select_comparison_calibrated_photometry(
                psf_data,
                aper_data,
                apers,
                annuli,
                airmass,
                exotic_infoDict['comp_stars'],
                sigma_display,
                skip_low_comparison_coverage_rejection=skip_low_comp_coverage_rejection,
                use_psf_photometry=use_psf_photometry,
                use_aperture_photometry=use_aperture_photometry,
            )

            if comparison_calibration is not None:
                log_info("\nCalibrating comparison stars before target fitting. Please wait.")
                log_info(f"Comparison-star field method: {comparison_calibration['method_label']}")
                log_info(f"Comparison-star field score: {comparison_calibration['field_score'] * 100.0:.4f}%")
                if comparison_calibration.get('suitability_outlier_rejected_count', 0) > 0:
                    threshold = comparison_calibration.get('suitability_high_threshold', np.nan)
                    if np.isfinite(threshold):
                        log_info(
                            "Comparison-star field sigma clipping rejected "
                            f"{comparison_calibration['suitability_outlier_rejected_count']} high-suitability "
                            f"outlier(s) above {threshold * 100.0:.4f}% before target-fit evaluation."
                        )
                    else:
                        log_info(
                            "Comparison-star field sigma clipping rejected "
                            f"{comparison_calibration['suitability_outlier_rejected_count']} high-suitability "
                            "outlier(s) before target-fit evaluation."
                        )
                if comparison_calibration.get('image_outlier_rejected_count', 0) > 0:
                    required_pairs = comparison_calibration.get('image_outlier_required_valid_pairs', 0)
                    sigma_threshold = comparison_calibration.get('image_outlier_sigma', COMPARISON_IMAGE_OUTLIER_SIGMA)
                    log_info(
                        "Comparison-star field image clipping rejected "
                        f"{comparison_calibration['image_outlier_rejected_count']} frame(s) after suitability clipping "
                        f"because every valid pairwise comparison was more than {sigma_threshold:.2f} sigma from "
                        f"its flat-line median (min valid pair count={required_pairs})."
                    )
                for summary in comparison_calibration['comp_summaries']:
                    aggregate_text = "n/a" if not np.isfinite(summary['aggregate_score']) else f"{summary['aggregate_score'] * 100.0:.4f}%"
                    ensemble_text = "n/a" if not np.isfinite(summary['ensemble_score']) else f"{summary['ensemble_score'] * 100.0:.4f}%"
                    pairwise_text = "n/a" if not np.isfinite(summary['pairwise_median_score']) else f"{summary['pairwise_median_score'] * 100.0:.4f}%"
                    selected_label = " [selected]" if summary['selected'] else ""
                    position_text = format_comp_star_position(summary['position'])
                    coverage_text = f"coverage={format_comp_star_coverage_text(summary)}"
                    if summary['coverage_rejected']:
                        coverage_text += " [rejected: low coverage]"
                    if summary.get('suitability_outlier_rejected'):
                        coverage_text += " [rejected: high suitability outlier]"
                    log_info(
                        f"  {summary['label']}{selected_label} ({position_text}): suitability={aggregate_text}, "
                        f"ensemble={ensemble_text}, pairwise_median={pairwise_text}, "
                        f"valid_pairs={summary['valid_pair_count']}, {coverage_text}, "
                        f"reason={summary['selection_reason']}"
                    )

                try:
                    plot_comp_star_pairwise_matrix(
                        comparison_calibration['pairwise_matrix'],
                        comparison_calibration['best_comp_index'],
                        pDict['pName'],
                        exotic_infoDict['save'],
                        exotic_infoDict['date'],
                        comparison_calibration['method_label'],
                    )
                    plot_comp_star_calibration_series(
                        times,
                        comparison_calibration['comp_summaries'],
                        pDict['pName'],
                        exotic_infoDict['save'],
                        exotic_infoDict['date'],
                        comparison_calibration['method_label'],
                    )
                    plot_individual_comp_star_calibration_series(
                        times,
                        comparison_calibration['comp_summaries'],
                        pDict['pName'],
                        exotic_infoDict['save'],
                        exotic_infoDict['date'],
                        comparison_calibration['method_label'],
                    )
                    plot_comp_star_suitability(
                        comparison_calibration['comp_summaries'],
                        pDict['pName'],
                        exotic_infoDict['save'],
                        exotic_infoDict['date'],
                        comparison_calibration['method_label'],
                    )
                    save_comp_star_calibration_summary(
                        exotic_infoDict['save'],
                        pDict['pName'],
                        exotic_infoDict['date'],
                        comparison_calibration['method_label'],
                        comparison_calibration['field_score'],
                        comparison_calibration['comp_summaries'],
                        comparison_calibration['best_comp_index'],
                    )
                except Exception as e:
                    log_info(f"Warning: Could not save comparison-star calibration outputs ({e}).", warn=True)

                comparison_fit_search = fit_ranked_comparison_calibration_candidates(
                    times,
                    jd_times,
                    airmass,
                    ld,
                    pDict,
                    comparison_calibration,
                    psf_data,
                    aper_data,
                    tFlux,
                    plot_time_range=full_plot_time_range,
                    disable_vertical_flux_normalization=disable_vertical_flux_normalization,
                    detrend_on_outoftransit_baseline=detrend_on_outoftransit_baseline,
                    use_impactparameter_rather_than_inclination_to_fit=
                    use_impactparameter_rather_than_inclination_to_fit,
                    use_eebls_to_initialize_tmid_and_bounds=use_eebls_tmid_initializer,
                    pick_comparison_by_eebls_snr=pick_comparison_by_eebls_snr,
                    assess_all_comparisons_before_selecting_best=assess_all_comparisons_before_selecting_best,
                    exit_at_first_qc_pass_solution=exit_at_first_qc_pass_solution,
                    final_fit_baseline_duration_multiplier=final_fit_baseline_duration_multiplier,
                    use_adaptive_apertures=use_adaptive_apertures,
                    adaptive_aperture_values=aperture_values,
                    adaptive_annulus_values=annulus_values,
                    fallback_sigma=sigma_display,
                    run_fast_ultranest_before_final_run=run_fast_ultranest_before_final_run,
                    save_dir=exotic_infoDict['save'],
                    planet_name=pDict['pName'],
                    observation_date=exotic_infoDict['date'],
                )
                comparison_calibration['ranked_fit_comp_indices'] = [
                    summary['comp_index'] for summary in comparison_fit_search['ranked_summaries']
                ]
                comparison_calibration['fit_attempt_summaries'] = comparison_fit_search['attempts']

                selected_attempt = comparison_fit_search['selected_result']
                fit_attempts = comparison_fit_search['attempts']
                if fit_attempts:
                    comparison_calibration['selected_fit_diagnostics'] = (
                        selected_attempt['fit_diagnostics']
                        if selected_attempt is not None
                        else fit_attempts[-1]['fit_diagnostics']
                    )
                    if selected_attempt is None or len(fit_attempts) > 1:
                        log_comparison_calibration_fit_attempt_summaries(
                            fit_attempts,
                            comparison_calibration['method_label'],
                        )

                if selected_attempt is not None:
                    selected_comp_index = selected_attempt['comp_index']
                    selected_ckey = selected_attempt['ckey']
                    selected_comp_coords = exotic_infoDict['comp_stars'][selected_comp_index]
                    selected_min_aperture = 0 if comparison_calibration['method'] == 'psf' else comparison_calibration['aper']
                    selected_min_annulus = comparison_calibration['annulus']
                    selected_a = None if comparison_calibration['method'] == 'psf' else comparison_calibration['a']
                    selected_an = None if comparison_calibration['method'] == 'psf' else comparison_calibration['an']
                    myfit = selected_attempt['fit']
                    tFlux1 = selected_attempt['tflux_fit']
                    cFlux1 = selected_attempt['cflux_fit']
                    selected_source_indices = np.asarray(
                        selected_attempt.get('source_indices', np.arange(len(tFlux1), dtype=int)),
                        dtype=int,
                    )
                    if selected_attempt.get('search_stopped_after_qc_pass', False):
                        selection_basis = 'first_qc_pass'
                    elif selected_attempt.get('search_stopped_after_promising_partial', False):
                        selection_basis = 'promising_partial'
                    elif selected_attempt.get('selected_despite_transit_qc', False):
                        selection_basis = 'comparison_field_qc_fallback'
                    elif selected_comp_index == comparison_calibration['best_comp_index']:
                        selection_basis = 'comparison_field'
                    else:
                        selection_basis = 'comparison_field_retry'
                    if selection_basis == 'first_qc_pass':
                        log_info(
                            "Comparison-star calibration target-fit selection chose "
                            f"Comp {selected_comp_index + 1} with {comparison_calibration['method_label']} "
                            "because it was the first candidate to pass transit QC."
                        )
                    elif selection_basis == 'promising_partial':
                        log_info(
                            "Comparison-star calibration target-fit selection chose "
                            f"Comp {selected_comp_index + 1} with {comparison_calibration['method_label']} "
                            "because pre-UltraNest preflight and the candidate fit indicated a promising "
                            "partial-coverage MARGINAL solution."
                        )
                    elif selection_basis == 'comparison_field_qc_fallback':
                        fallback_selection_metric = comparison_fit_search.get('selection_metric', 'ktmf')
                        if fallback_selection_metric == 'ktmf':
                            fallback_metric_value = format_ktmf_metric(
                                selected_attempt.get('ktmf_metric', np.nan)
                            )
                        elif fallback_selection_metric == 'eebls_snr':
                            fallback_metric_value = format_eebls_snr(
                                selected_attempt.get('eebls_snr', np.nan)
                            )
                        else:
                            fallback_metric_value = format_transit_delta_bic(
                                selected_attempt.get('transit_delta_bic', np.nan)
                            )
                        log_info(
                            "Warning: all completed comparison-star target fits were rejected by transit QC; "
                            "continuing with the best available fit "
                            f"(Comp {selected_comp_index + 1}, "
                            f"{comparison_selection_metric_label(fallback_selection_metric)}="
                            f"{fallback_metric_value}) so final outputs are still produced.",
                            warn=True,
                        )
                    elif selection_basis == 'comparison_field_retry':
                        retry_count = selected_attempt['rank']
                        log_info(
                            "Comparison-star calibration target-fit selection chose "
                            f"Comp {selected_comp_index + 1} with {comparison_calibration['method_label']} "
                            f"after evaluating {retry_count} better-ranked field-stability candidate(s); "
                            f"it delivered the best {comparison_selection_metric_label(comparison_fit_search['selection_metric'])} "
                            "among successful fits."
                        )

                    photometry_info.update(best_fit_lc=myfit,
                                           comp_star_num=selected_comp_index + 1,
                                           comp_star_coords=selected_comp_coords,
                                           min_aperture=selected_min_aperture,
                                           min_annulus=selected_min_annulus,
                                           aperture_index=selected_a,
                                           annulus_index=selected_an,
                                           reuse_selected_full_reduction_fit=bool(
                                               selected_attempt.get('full_reduction_applied', False)
                                               and selected_attempt.get('fit') is not None
                                           ),
                                           selected_source_indices=selected_source_indices,
                                           selected_fit_good_times=selected_attempt.get('good_times'),
                                           selected_fit_good_flux=selected_attempt.get('good_flux'),
                                           selected_fit_good_unc=selected_attempt.get('good_unc'),
                                           selected_fit_good_airmass=selected_attempt.get('good_airmass'),
                                           selected_fit_duration_samples=selected_attempt.get('duration_samples'),
                                           selected_fit_data_highres=selected_attempt.get('data_highres'),
                                           selected_fit_final_output_dir=selected_attempt.get('final_output_dir'),
                                           calibration_field_score=comparison_calibration['field_score'],
                                           selection_basis=selection_basis,
                                           selection_metric=comparison_fit_search.get('selection_metric', 'ktmf'),
                                           selected_comparison_selection_reason=selected_attempt.get('selection_reason'),
                                           selected_comparison_attempt=compact_comparison_attempt_for_output(selected_attempt),
                                           comparison_fit_attempt_summaries=[
                                               compact_comparison_attempt_for_output(attempt)
                                               for attempt in comparison_fit_search.get('attempts', [])
                                           ],
                                           comparison_ktmf_metric=selected_attempt.get('ktmf_metric', np.nan),
                                           selected_comparison_ktmf_contributions=selected_attempt.get('ktmf_contributions') or [],
                                           comparison_eebls_snr=selected_attempt.get('eebls_snr', np.nan),
                                           comparison_transit_delta_bic=selected_attempt.get('transit_delta_bic', np.nan),
                                           selected_comparison_fit_point_count=selected_attempt.get('fit_point_count'),
                                           selected_comparison_transit_qc_status=selected_attempt.get('transit_qc_status'),
                                           selected_comparison_transit_qc_summary=selected_attempt.get('transit_qc_summary'))

                    flux_values.update(flux_tar=tFlux1, flux_ref=cFlux1,
                                       flux_unc_tar=tFlux1 ** 0.5, flux_unc_ref=cFlux1 ** 0.5)

                    centroid_positions.update(x_targ=psf_data["target"][selected_source_indices, 0],
                                              y_targ=psf_data["target"][selected_source_indices, 1],
                                              x_ref=psf_data[selected_ckey][selected_source_indices, 0],
                                              y_ref=psf_data[selected_ckey][selected_source_indices, 1])

                    if selected_comp_index in vsp_num:
                        ref_flux[selected_comp_index] = {
                            'myfit': myfit,
                            'pos': exotic_infoDict['comp_stars'][selected_comp_index]
                        }

                    if vsp_num:
                        if comparison_calibration['method'] == 'psf':
                            for j in vsp_num:
                                ckey = f"comp{j + 1}"
                                cFlux = 2 * np.pi * psf_data[ckey][:, 2] * psf_data[ckey][:, 3] * psf_data[ckey][:, 4]
                                vsp_fit, _, _ = fit_lightcurve(
                                    times, tFlux, cFlux, airmass, ld, pDict, jd_times,
                                    disable_vertical_flux_normalization=disable_vertical_flux_normalization,
                                    use_impactparameter_rather_than_inclination_to_fit=
                                    use_impactparameter_rather_than_inclination_to_fit,
                                    plot_time_range=full_plot_time_range,
                                    use_eebls_to_initialize_tmid_and_bounds=use_eebls_tmid_initializer,
                                )
                                ref_flux[j] = {
                                    'myfit': vsp_fit,
                                    'pos': exotic_infoDict['comp_stars'][j]
                                }
                        else:
                            best_a = comparison_calibration['a']
                            best_an = comparison_calibration['an']
                            best_target_flux = aper_data['target'][:, best_a, best_an]
                            for j in vsp_num:
                                ckey = f"comp{j + 1}"
                                aper_mask = np.isfinite(aper_data[ckey][:, best_a, best_an])
                                cFlux = aper_data[ckey][aper_mask][:, best_a, best_an]
                                vsp_fit, _, _ = fit_lightcurve(
                                    times[aper_mask], best_target_flux[aper_mask], cFlux,
                                    airmass[aper_mask], ld, pDict, jd_times[aper_mask],
                                    disable_vertical_flux_normalization=disable_vertical_flux_normalization,
                                    use_impactparameter_rather_than_inclination_to_fit=
                                    use_impactparameter_rather_than_inclination_to_fit,
                                    plot_time_range=full_plot_time_range,
                                    use_eebls_to_initialize_tmid_and_bounds=use_eebls_tmid_initializer,
                                )
                                ref_flux[j] = {
                                    'myfit': vsp_fit,
                                    'pos': exotic_infoDict['comp_stars'][j]
                                }
                else:
                    if fit_attempts:
                        failed_attempt = fit_attempts[-1]
                        failed_comp_index = failed_attempt['comp_index']
                        failure_reason = failed_attempt['fit_diagnostics'].get(
                            'failure_reason',
                            "the full comparison-star candidate reduction did not converge to a usable solution.",
                        )
                        attempted_count = len(fit_attempts)
                        ranked_count = len(comparison_fit_search['ranked_summaries'])
                        log_info(
                            "Error: Comparison-star calibration exhausted "
                            f"{attempted_count}/{ranked_count} ranked comparison star(s) for "
                            f"{comparison_calibration['method_label']} without a usable fully reduced target fit "
                            f"(last attempt: Comp {failed_comp_index + 1}; reason: {failure_reason}).",
                            error=True,
                        )
                    else:
                        log_info(
                            "Error: Comparison-star calibration did not produce any coverage-qualified "
                            "comparison stars to fully reduce against the target fit.",
                            error=True,
                        )
                    return

            update_photometry_adaptive_summary(
                photometry_info,
                use_adaptive_apertures,
                aperture_values,
                annulus_values,
                psf_data['target'],
                fallback_sigma=sigma,
            )

            if require_comp_star and photometry_info['comp_star_num'] is None:
                log_info(
                    "Error: require_comp_star is enabled, but every evaluated comparison-star candidate "
                    "was rejected or failed to complete a usable full reduction. See the comparison-star "
                    "calibration fit diagnostics above for the per-candidate failure reasons.",
                    error=True,
                )
                return

            log_info("\n\n*********************************************")
            if np.isfinite(photometry_info['calibration_field_score']):
                log_info(f"Comparison-Star Field Score: {round(photometry_info['calibration_field_score'] * 100, 4)}%")
            summary_min_aperture = photometry_info.get('min_aperture')
            if photometry_info.get('comp_star_num') is not None or (
                summary_min_aperture is not None and summary_min_aperture < 0
            ):
                log_info(
                    "Comparison Selection Metric: "
                    f"{comparison_selection_metric_label(photometry_info.get('selection_metric', 'ktmf'))}"
                )
            if np.isfinite(photometry_info.get('comparison_ktmf_metric', np.nan)):
                log_info(f"Selected Comparison KTMF: {photometry_info['comparison_ktmf_metric']:.2f} / 5.00")
            if np.isfinite(photometry_info.get('comparison_eebls_snr', np.nan)):
                log_info(f"Selected Comparison EEBLS SNR: {photometry_info['comparison_eebls_snr']:.2f}")
            if np.isfinite(photometry_info.get('comparison_transit_delta_bic', np.nan)):
                log_info(
                    "Selected Comparison Transit Delta BIC: "
                    f"{photometry_info['comparison_transit_delta_bic']:.2f}"
                )
            selected_method_label = selected_photometry_method_label(photometry_info)
            display_aperture, display_annulus = reported_photometry_aperture_radii(photometry_info)
            adaptive_summary = photometry_info.get('adaptive_summary')
            if photometry_info['min_aperture'] == 0:  # psf
                log_info(f"Best Comparison Star: #{photometry_info['comp_star_num']}")
                log_info("Optimal Method: PSF photometry")
            elif photometry_info['min_aperture'] < 0:  # no comp star
                log_info("Best Comparison Star: None")
                if adaptive_summary is not None:
                    log_info(f"Optimal Aperture: {abs(display_aperture):.2f} +/- {adaptive_summary['aperture_std']:.2f} px")
                    log_info(f"Optimal Annulus: {display_annulus:.2f} +/- {adaptive_summary['annulus_std']:.2f} px")
                    log_info(f"Adaptive Aperture Scale: {adaptive_summary['aperture_sigma']:.2f} sigma")
                    log_info(f"Adaptive Annulus Scale: {adaptive_summary['annulus_sigma']:.2f} sigma")
                    log_info(f"Aperture Range: {adaptive_summary['aperture_min']:.2f} to {adaptive_summary['aperture_max']:.2f} px")
                    log_info(f"Annulus Range: {adaptive_summary['annulus_min']:.2f} to {adaptive_summary['annulus_max']:.2f} px")
                else:
                    log_info(f"Optimal Aperture: {abs(np.round(display_aperture, 2))}")
                    log_info(f"Optimal Annulus: {np.round(display_annulus, 2)}")
            else:
                log_info(f"Best Comparison Star: #{photometry_info['comp_star_num']}")
                if adaptive_summary is not None:
                    log_info(f"Optimal Aperture: {display_aperture:.2f} +/- {adaptive_summary['aperture_std']:.2f} px")
                    log_info(f"Optimal Annulus: {display_annulus:.2f} +/- {adaptive_summary['annulus_std']:.2f} px")
                    log_info(f"Adaptive Aperture Scale: {adaptive_summary['aperture_sigma']:.2f} sigma")
                    log_info(f"Adaptive Annulus Scale: {adaptive_summary['annulus_sigma']:.2f} sigma")
                    log_info(f"Aperture Range: {adaptive_summary['aperture_min']:.2f} to {adaptive_summary['aperture_max']:.2f} px")
                    log_info(f"Annulus Range: {adaptive_summary['annulus_min']:.2f} to {adaptive_summary['annulus_max']:.2f} px")
                else:
                    log_info(f"Optimal Aperture: {np.round(display_aperture, 2)}")
                    log_info(f"Optimal Annulus: {np.round(display_annulus, 2)}")
            log_info("*********************************************\n")

            best_fit_lc = photometry_info['best_fit_lc']
            bestCompStar = photometry_info['comp_star_num']
            comp_coords = photometry_info['comp_star_coords']
            log_lightcurve_filter_diagnostics(
                getattr(best_fit_lc, 'frame_filter_diagnostics', []),
                header="Selected lightcurve frame rejections during target fitting",
            )
            try:
                selected_photometry_debug_path = save_selected_photometry_debug_series(
                    exotic_infoDict['save'],
                    pDict['pName'],
                    exotic_infoDict['date'],
                    best_fit_lc,
                )
                if selected_photometry_debug_path is not None:
                    log_info(
                        f"Saved selected raw target/reference ratio diagnostics to "
                        f"{selected_photometry_debug_path}."
                    )
            except Exception as e:
                log_info(
                    f"Warning: Could not save selected raw target/reference ratio diagnostics ({e}).",
                    warn=True,
                )

            if fit_every_comparison_candidate and exotic_infoDict['comp_stars']:
                candidate_fit_summaries = fit_lightcurve_to_every_comparison_candidate(
                    times,
                    jd_times,
                    airmass,
                    ld,
                    pDict,
                    exotic_infoDict['comp_stars'],
                    psf_data,
                    aper_data,
                    photometry_info,
                    plot_time_range=full_plot_time_range,
                    disable_vertical_flux_normalization=disable_vertical_flux_normalization,
                    skip_low_comparison_coverage_rejection=skip_low_comp_coverage_rejection,
                    use_impactparameter_rather_than_inclination_to_fit=
                    use_impactparameter_rather_than_inclination_to_fit,
                    use_eebls_to_initialize_tmid_and_bounds=use_eebls_tmid_initializer,
                )
                saved_candidate_fit_count = sum(1 for summary in candidate_fit_summaries if summary['fit'] is not None)
                failed_candidate_fit_count = len(candidate_fit_summaries) - saved_candidate_fit_count
                if candidate_fit_summaries:
                    log_comparison_candidate_fit_summaries(candidate_fit_summaries, photometry_info)
                    try:
                        plot_comp_star_candidate_lightcurve_fits(
                            candidate_fit_summaries,
                            pDict['pName'],
                            exotic_infoDict['save'],
                            exotic_infoDict['date'],
                            selected_method_label,
                        )
                        log_info(
                            f"Saved {saved_candidate_fit_count} comparison-candidate lightcurve fit plot(s) to temp/."
                        )
                        if failed_candidate_fit_count:
                            log_info(
                                f"Skipped {failed_candidate_fit_count} comparison candidate(s) that did not yield a usable lightcurve fit."
                            )
                    except Exception as e:
                        log_info(f"Warning: Could not save comparison-candidate lightcurve plots ({e}).", warn=True)

            # save psf_data to disk for best comparison star
            if bestCompStar:
                np.savetxt(Path(exotic_infoDict['save']) / "temp" / "psf_data_comp.txt", psf_data[f"comp{bestCompStar}"],
                            header="#x_centroid, y_centroid, amplitude, sigma_x, sigma_y, rotation offset",
                            fmt="%.6f")

            reuse_selected_full_reduction_fit = bool(
                photometry_info.get('reuse_selected_full_reduction_fit', False)
                and best_fit_lc is not None
            )
            psf_selection_indices = np.asarray(
                photometry_info.get('selected_source_indices', np.arange(len(best_fit_lc.time))),
                dtype=int,
            )
            if psf_selection_indices.shape[0] != len(best_fit_lc.time):
                psf_selection_indices = np.arange(len(best_fit_lc.time), dtype=int)

            # sigma clip
            if reuse_selected_full_reduction_fit:
                log_info(
                    "Reusing the selected comparison-star full-reduction ultranest fit; "
                    "skipping duplicate selected-only final-fit clipping."
                )
                si = np.arange(len(best_fit_lc.time), dtype=int)
                time_clip_mask = np.zeros(len(best_fit_lc.time), dtype=bool)
                phase_clip_mask = np.zeros_like(time_clip_mask, dtype=bool)
                adaptive_clip_mask = np.zeros_like(time_clip_mask, dtype=bool)
            else:
                si = np.argsort(best_fit_lc.time)
                dt = np.mean(np.diff(np.sort(best_fit_lc.time)))
                ndt = int(30. / 24. / 60. / dt) * 2 + 1  # ~30 minutes
                time_clip_mask = sigma_clip(best_fit_lc.data[si], sigma=3, dt=ndt, times=best_fit_lc.time[si])
                phase_clip_mask = np.zeros_like(time_clip_mask, dtype=bool)
                if hasattr(best_fit_lc, 'residuals') and hasattr(best_fit_lc, 'phase'):
                    phase_clip_mask = phase_bin_sigma_clip(best_fit_lc.residuals[si], best_fit_lc.phase[si], sigma=3, bins=10)
                adaptive_clip_mask = np.zeros_like(time_clip_mask, dtype=bool)
                if use_adaptive_apertures and adaptive_summary is not None:
                    sorted_apertures = np.asarray(adaptive_summary['aperture_series'], dtype=float)[si]
                    sorted_annuli = np.asarray(adaptive_summary['annulus_series'], dtype=float)[si]
                    retained_mask = adaptive_aperture_outlier_mask(sorted_apertures[~time_clip_mask & ~phase_clip_mask],
                                                                   sorted_annuli[~time_clip_mask & ~phase_clip_mask])
                    adaptive_clip_mask[~time_clip_mask & ~phase_clip_mask] = retained_mask
            gi = ~(time_clip_mask | phase_clip_mask | adaptive_clip_mask)  # good indexs
            prefinal_filter_diagnostics = [
                build_time_rejection_diagnostic(
                    "Final-fit time sigma clip",
                    np.asarray(best_fit_lc.time, dtype=float)[si],
                    ~time_clip_mask,
                    note="Dropped time-series outliers before the final fit.",
                ),
                build_time_rejection_diagnostic(
                    "Final-fit phase residual clip",
                    np.asarray(best_fit_lc.time, dtype=float)[si],
                    ~phase_clip_mask,
                    note="Dropped phase-binned residual outliers before the final fit.",
                ),
                build_time_rejection_diagnostic(
                    "Final-fit adaptive-aperture clip",
                    np.asarray(best_fit_lc.time, dtype=float)[si],
                    ~adaptive_clip_mask,
                    note="Dropped adaptive-aperture radius outliers before the final fit.",
                ),
            ]
            phase_clip_removed = np.count_nonzero(phase_clip_mask & ~time_clip_mask)
            if phase_clip_removed:
                log_info(f"Removed {phase_clip_removed} phase-binned residual outlier(s) before final fit.")
            adaptive_clip_removed = np.count_nonzero(adaptive_clip_mask)
            if adaptive_clip_removed:
                log_info(f"Removed {adaptive_clip_removed} adaptive-aperture radius outlier(s) before final fit.")

            if np.isnan(best_fit_lc.data).all():
                log_info("Error: No valid photometry data found.", error=True)
                return

            apply_lightcurve_mask(best_fit_lc, gi, sort_index=si)

            goodTimes = best_fit_lc.time
            goodAirmasses = best_fit_lc.airmass

            if reuse_selected_full_reduction_fit:
                selected_good_flux = photometry_info.get('selected_fit_good_flux')
                selected_good_unc = photometry_info.get('selected_fit_good_unc')
                if selected_good_flux is not None and np.shape(selected_good_flux) == np.shape(goodTimes):
                    goodFluxes = np.asarray(selected_good_flux, dtype=float)
                else:
                    goodFluxes = np.asarray(best_fit_lc.detrended, dtype=float)
                if selected_good_unc is not None and np.shape(selected_good_unc) == np.shape(goodTimes):
                    goodNormUnc = np.asarray(selected_good_unc, dtype=float)
                else:
                    goodNormUnc = np.asarray(best_fit_lc.detrendederr, dtype=float)
            else:
                final_fit_series = prepare_final_fit_lightcurve_series(best_fit_lc)
                if not final_fit_series.get('applied'):
                    log_info(
                        f"Warning: {final_fit_series.get('note', 'could not prepare the final-fit light curve from the selected fit.')} "
                        "Falling back to the current detrended light curve arrays.",
                        warn=True,
                    )
                    goodFluxes = np.asarray(best_fit_lc.detrended, dtype=float)
                    goodNormUnc = np.asarray(best_fit_lc.detrendederr, dtype=float)
                else:
                    log_info(final_fit_series['note'])
                    goodFluxes = np.asarray(final_fit_series['flux'], dtype=float)
                    goodNormUnc = np.asarray(final_fit_series['unc'], dtype=float)

            centroid_positions.update(x_targ=centroid_positions['x_targ'][si][gi],
                                      y_targ=centroid_positions['y_targ'][si][gi],
                                      x_ref=centroid_positions['x_ref'][si][gi],
                                      y_ref=centroid_positions['y_ref'][si][gi])

            flux_values.update(flux_tar=flux_values['flux_tar'][si][gi],
                               flux_ref=flux_values['flux_ref'][si][gi],
                               flux_unc_tar=flux_values['flux_unc_tar'][si][gi],
                               flux_unc_ref=flux_values['flux_unc_ref'][si][gi])

            relative_flux_mask = relative_flux_filter_mask(goodFluxes)
            prefinal_filter_diagnostics.append(build_time_rejection_diagnostic(
                "Final relative-flux validity filter",
                goodTimes,
                relative_flux_mask,
                note="Dropped non-finite or non-positive normalized flux values before the final fit.",
            ))
            if np.count_nonzero(relative_flux_mask) == 0:
                log_info(
                    "Error: No valid photometry data found after removing non-finite or non-positive relative flux values.",
                    error=True,
                )
                return

            log_lightcurve_filter_diagnostics(
                prefinal_filter_diagnostics,
                header="Selected lightcurve frame rejections before the final fit",
            )
            apply_lightcurve_mask(best_fit_lc, relative_flux_mask)

            goodTimes = goodTimes[relative_flux_mask]
            goodFluxes = goodFluxes[relative_flux_mask]
            goodNormUnc = goodNormUnc[relative_flux_mask]
            goodAirmasses = goodAirmasses[relative_flux_mask]

            centroid_positions.update(x_targ=centroid_positions['x_targ'][relative_flux_mask],
                                      y_targ=centroid_positions['y_targ'][relative_flux_mask],
                                      x_ref=centroid_positions['x_ref'][relative_flux_mask],
                                      y_ref=centroid_positions['y_ref'][relative_flux_mask])

            flux_values.update(flux_tar=flux_values['flux_tar'][relative_flux_mask],
                               flux_ref=flux_values['flux_ref'][relative_flux_mask],
                               flux_unc_tar=flux_values['flux_unc_tar'][relative_flux_mask],
                               flux_unc_ref=flux_values['flux_unc_ref'][relative_flux_mask])

            psf_selection_indices = psf_selection_indices[si][gi][relative_flux_mask]
            obs_stats_sort_index = psf_selection_indices
            obs_stats_keep_mask = np.ones(psf_selection_indices.shape[0], dtype=bool)

            update_photometry_adaptive_summary(
                photometry_info,
                use_adaptive_apertures,
                aperture_values,
                annulus_values,
                psf_data['target'][psf_selection_indices],
                fallback_sigma=sigma,
            )
            display_aperture, display_annulus = reported_photometry_aperture_radii(photometry_info)


            if photometry_info['min_aperture'] == 0:
                opt_method = "PSF"
                # Calculate min_aper and min_annulus using the stdev
                stdev_fov = (psf_data['target'][:, 3] + psf_data['target'][:, 4]) * 0.5
                min_aper_fov = float(5 * stdev_fov.mean())
                min_annulus_fov = float(15 * stdev_fov.mean())
            else:
                opt_method = "Aperture"
                min_aper_fov = float(display_aperture)
                min_annulus_fov = float(display_annulus)

            fov_aperture = min_aper_fov if opt_method == "PSF" else float(display_aperture)
            fov_annulus = min_annulus_fov if opt_method == "PSF" else float(display_annulus)
            fov_sky_geometry = resolve_sky_annulus_geometry(
                fov_aperture,
                fov_annulus,
                psf_sigma=sigma_display,
            )

            plot_fov(fov_aperture, fov_annulus, sigma_display,
                     centroid_positions['x_targ'][0], centroid_positions['y_targ'][0],
                     centroid_positions['x_ref'][0], centroid_positions['y_ref'][0],
                     firstImage, img_scale_str, pDict['pName'], exotic_infoDict['save'],
                     exotic_infoDict['date'], opt_method, min_aper_fov, min_annulus_fov,
                     sky_inner_radius=fov_sky_geometry['inner_radius'],
                     sky_outer_radius=fov_sky_geometry['outer_radius'])

            plot_centroids(centroid_positions['x_targ'], centroid_positions['y_targ'],
                           centroid_positions['x_ref'], centroid_positions['y_ref'],
                           goodTimes, pDict['pName'], exotic_infoDict['save'], exotic_infoDict['date'])

            plot_flux(goodTimes, flux_values['flux_tar'], flux_values['flux_unc_tar'],
                      flux_values['flux_ref'], flux_values['flux_unc_ref'],
                      goodFluxes, goodNormUnc, goodAirmasses, pDict['pName'], exotic_infoDict['save'],
                      exotic_infoDict['date'])

            adaptive_summary = photometry_info.get('adaptive_summary')
            if adaptive_summary is not None:
                plot_adaptive_aperture_diagnostics(
                    goodTimes,
                    adaptive_summary['aperture_series'],
                    adaptive_summary['annulus_series'],
                    adaptive_summary['fwhm_series'],
                    goodAirmasses,
                    pDict['pName'],
                    exotic_infoDict['save'],
                    exotic_infoDict['date'],
                    adaptive_summary['aperture_sigma'],
                    adaptive_summary['annulus_sigma'],
                )

            # TODO: convert the exoplanet archive mid transit time to bjd - need to take into account observatory location listed in Exoplanet Archive
            # tMidtoC = astropy.time.Time(timeMidTransit, format='jd', scale='utc')
            # forPhaseResult = JDUTC_to_BJDTDB(tMidtoC, ra=raDeg, dec=decDeg, lat=lati, longi=longit, alt=2000)
            # bjdMidTOld = float(forPhaseResult[0])
            # bjdMidTOld = pDict['midT']

            # goodPhasesList = []
            # convert all the phases based on the updated bjd times
            # for convertedTime in goodTimes:
            #     bjdPhase = getPhase(float(convertedTime), pDict['pPer'], bjdMidTOld)
            #     goodPhasesList.append(bjdPhase)
            # goodPhases = np.array(goodPhasesList)

            # Calculate the standard deviation of the normalized flux values
            # standardDev1 = np.std(goodFluxes)

            if vsp_comp_stars:
                if not bestCompStar:
                    vsp_params = stellar_variability(ref_flux, best_fit_lc, exotic_infoDict['comp_stars'],
                                                      vsp_comp_stars, vsp_num, None, exotic_infoDict['save'],
                                                      pDict['sName'],
                                                      observed_filter=exotic_infoDict.get('observed_filter',
                                                                                          exotic_infoDict.get('filter')))
                else:
                    vsp_params = stellar_variability(ref_flux, best_fit_lc, exotic_infoDict['comp_stars'],
                                                      vsp_comp_stars, vsp_num, bestCompStar - 1, exotic_infoDict['save'],
                                                      pDict['sName'],
                                                      observed_filter=exotic_infoDict.get('observed_filter',
                                                                                          exotic_infoDict.get('filter')))

            log_info("\n\nOutput File Saved")
        else:
            goodTimes, goodFluxes, goodNormUnc, goodAirmasses = [], [], [], []
            bestCompStar, comp_coords = None, None
            exotic_infoDict.setdefault('observed_filter', exotic_infoDict.get('filter'))
            ld, ld0, ld1, ld2, ld3 = get_ld_values(pDict, exotic_infoDict)

            with exotic_infoDict['prered_file'].open('r') as f:
                for processed_data in f:
                    try:
                        processed_data = processed_data.split(',')
                        goodTimes.append(float(processed_data[0]))
                        goodFluxes.append(float(processed_data[1]))
                        goodNormUnc.append(float(processed_data[2]))
                        goodAirmasses.append(float(processed_data[3]))
                    except ValueError:
                        continue

            goodTimes = np.array(goodTimes)
            goodFluxes = np.array(goodFluxes)
            goodNormUnc = np.array(goodNormUnc)
            goodAirmasses = np.array(goodAirmasses)

            if exotic_infoDict['file_time'] != 'BJD_TDB':
                missing_location = [
                    label for key, label in (('long', 'longitude'), ('lat', 'latitude'), ('elev', 'elevation'))
                    if exotic_infoDict.get(key) is None
                ]
                if missing_location:
                    log_info("Error: Longitude, latitude, and elevation are required to convert "
                             f"pre-reduced {exotic_infoDict['file_time']} timestamps to BJD_TDB.", error=True)
                    return
                time_offset = 2400000.5 if exotic_infoDict['file_time'] == 'MJD_UTC' else 0.0
                goodTimes = convert_jd_to_bjd([time_ + time_offset for time_ in goodTimes], pDict, exotic_infoDict)

            if exotic_infoDict['file_units'] != 'flux':
                print("check flux convert")
                goodFluxes, goodNormUnc = flux_conversion(goodFluxes, goodNormUnc, exotic_infoDict['file_units'])

            relative_flux_mask = relative_flux_filter_mask(goodFluxes)
            if np.count_nonzero(relative_flux_mask) == 0:
                log_info(
                    "Error: No valid photometry data found after removing non-finite or non-positive relative flux values.",
                    error=True,
                )
                return

            goodTimes = goodTimes[relative_flux_mask]
            goodFluxes = goodFluxes[relative_flux_mask]
            goodNormUnc = goodNormUnc[relative_flux_mask]
            goodAirmasses = goodAirmasses[relative_flux_mask]
            goodFluxes, goodNormUnc, _ = normalize_flux_series_to_approximate_unity(
                goodFluxes,
                goodNormUnc,
            )
            finite_plot_times = goodTimes[np.isfinite(goodTimes)]
            full_plot_time_range = None
            if finite_plot_times.size:
                full_plot_time_range = (float(np.min(finite_plot_times)), float(np.max(finite_plot_times)))

        # for k in myfit.bounds.keys():
        #     print(f"{myfit.parameters[k]:.6f} +- {myfit.errors[k]}")

        if args.photometry:
            log_info("\nPhotometric Extraction Complete.")
            return

        log_info("\n")
        log_info("****************************************")
        log_info("Fitting a Light Curve Model to Your Data")
        log_info("****************************************\n")

        reuse_selected_final_model = bool(
            fitsortext == 1
            and photometry_info.get('reuse_selected_full_reduction_fit', False)
            and photometry_info.get('best_fit_lc') is not None
        )

        ##########################
        # NESTED SAMPLING FITTING
        ##########################

        prior = {
            'rprs': pDict['rprs'],  # Rp/Rs
            'ars': pDict['aRs'],  # a/Rs
            'per': pDict['pPer'],  # Period [day]
            'inc': pDict['inc'],  # Inclination [deg]
            'u0': ld0[0], 'u1': ld1[0], 'u2': ld2[0], 'u3': ld3[0],  # limb darkening (nonlinear)
            'ecc': pDict['ecc'],  # Eccentricity
            'omega': pDict['omega'],  # Arg of periastron
            'tmid': pDict['midT'],  # time of mid transit [day]
            'a2': 0,  # Flux lower bound
        }

        phase = (goodTimes - prior['tmid']) / prior['per']
        expected_duration = estimate_transit_duration_from_prior_geometry(prior)
        ephemeris_tmid_search_summary = estimate_ephemeris_tmid_and_bounds(
            goodTimes,
            pDict['midT'],
            prior['per'],
            pDict['midTUnc'],
            pDict['pPerUnc'],
            expected_duration=expected_duration,
            sigma_multiplier=35.0,
        )
        prior['tmid'] = ephemeris_tmid_search_summary['tmid']
        lower, upper = ephemeris_tmid_search_summary['bounds']

        if np.floor(phase).max() - np.floor(phase).min() == 0:
            log_info("Error: Estimated mid-transit not in observation range (check priors or observation time)", error=True)
            log_info(f"start:{np.min(goodTimes)}", error=True)
            log_info(f"  end:{np.max(goodTimes)}", error=True)
            log_info(f"prior:{prior['tmid']}", error=True)

        if ephemeris_tmid_search_summary.get('duration_capped'):
            log_info(ephemeris_tmid_search_summary['note'])
        eebls_tmid_search_summary = None
        if use_eebls_tmid_initializer:
            eebls_tmid_search_summary = estimate_tmid_and_bounds_with_eebls(
                goodTimes,
                goodFluxes,
                goodNormUnc,
                prior,
                [lower, upper],
            )
            log_info(eebls_tmid_search_summary['note'])
            if eebls_tmid_search_summary.get('applied'):
                prior['tmid'] = eebls_tmid_search_summary['tmid']
                lower, upper = eebls_tmid_search_summary['bounds']

        final_airmass_span = airmass_span(goodAirmasses)
        airmass_skip_note = None
        skip_final_airmass_fit = bool(exotic_infoDict.get('airmass_already_corrected'))
        if skip_final_airmass_fit:
            airmass_skip_note = (
                "Skipped (input AAVSO file already reports AIRMASS, AIRMASS CORRECTION FUNCTION); "
                "no airmass correction applied."
            )
            log_info(
                "Input AAVSO file reports AIRMASS, AIRMASS CORRECTION FUNCTION; "
                "skipping airmass fitting and applying no airmass correction."
            )
        elif should_skip_airmass_fit(goodAirmasses):
            skip_final_airmass_fit = True
            log_info(
                f"Airmass span {final_airmass_span:.4f} <= {AIRMASS_FLAT_RANGE_THRESHOLD:.2f}; "
                "skipping airmass fitting and applying no airmass correction."
            )

        mybounds = build_initial_transit_bounds(
            prior,
            [lower, upper],
            ars_unc=pDict.get('aRsUnc'),
        )
        apply_vertical_flux_normalization_bound(
            prior,
            mybounds,
            goodFluxes,
            disable_vertical_flux_normalization,
        )
        if not skip_final_airmass_fit:
            mybounds['a2'] = [-3, 3]

        if np.isnan(goodFluxes).all():
            log_info("Error: No valid photometry data found.", error=True)
            return

        if reuse_selected_final_model:
            myfit = photometry_info['best_fit_lc']
            goodTimes = np.asarray(getattr(myfit, 'time', goodTimes), dtype=float)
            goodAirmasses = np.asarray(getattr(myfit, 'airmass', goodAirmasses), dtype=float)
            reused_flux = photometry_info.get('selected_fit_good_flux')
            reused_unc = photometry_info.get('selected_fit_good_unc')
            if reused_flux is not None and np.shape(reused_flux) == np.shape(goodTimes):
                goodFluxes = np.asarray(reused_flux, dtype=float)
            else:
                goodFluxes = np.asarray(getattr(myfit, 'detrended', goodFluxes), dtype=float)
            if reused_unc is not None and np.shape(reused_unc) == np.shape(goodTimes):
                goodNormUnc = np.asarray(reused_unc, dtype=float)
            else:
                goodNormUnc = np.asarray(getattr(myfit, 'detrendederr', goodNormUnc), dtype=float)
            log_info(
                "Using the selected comparison-star full-reduction ultranest fit for final outputs; "
                "no additional final nested-sampling fit is being run."
            )
        else:
            # final light curve fit
            myfit, goodFluxes, goodNormUnc = fit_final_lightcurve_with_oot_baseline_detrending(
                goodTimes,
                goodFluxes,
                goodNormUnc,
                goodAirmasses,
                prior,
                mybounds,
                skip_airmass_fit=skip_final_airmass_fit,
                airmass_skip_note=airmass_skip_note,
                disable_vertical_flux_normalization=disable_vertical_flux_normalization,
                detrend_on_outoftransit_baseline=detrend_on_outoftransit_baseline,
                use_impactparameter_rather_than_inclination_to_fit=
                use_impactparameter_rather_than_inclination_to_fit,
                plot_time_range=full_plot_time_range,
                baseline_duration_multiplier=final_fit_baseline_duration_multiplier,
                expected_planet_dict=pDict,
                expected_tmid_search_summary=ephemeris_tmid_search_summary,
                eebls_search_summary=eebls_tmid_search_summary,
            )
        if (
            reuse_selected_final_model
            and getattr(myfit, 'sparse_posterior_live_point_extension_note', None) is None
        ):
            myfit = extend_selected_comparison_live_points_if_needed(myfit)
            annotate_transit_detection_qc(myfit)
        # myfit.dataerr *= np.sqrt(myfit.chi2 / myfit.data.shape[0])  # scale errorbars by sqrt(rchi2)
        # myfit.detrendederr *= np.sqrt(myfit.chi2 / myfit.data.shape[0])

        if fitsortext != 1 and not vsp_params:
            try:
                calibration_label, calibration_star = nextastro_prereduced_calibration_star(
                    exotic_infoDict.get('phot_comp_star'),
                    exotic_infoDict.get('filter'),
                )
                if calibration_star:
                    vsp_params = build_stellar_variability_params_from_fit(
                        myfit,
                        calibration_star,
                        calibration_star.get('pos'),
                        calibration_label,
                        exotic_infoDict['save'],
                        pDict['sName'],
                        observed_filter=exotic_infoDict.get('observed_filter', exotic_infoDict.get('filter')),
                    )
                    if not auid:
                        auid = vsx_auid(pDict['ra'], pDict['dec'])
                else:
                    log_info(
                        "\nWarning: Could not create pre-reduced stellar variability output because "
                        "no comparison-star RA/Dec with a usable NextAstro catalog magnitude was available.",
                        warn=True,
                    )
            except Exception as exc:
                log_info(
                    f"\nWarning: Could not create pre-reduced stellar variability output "
                    f"({describe_retry_exception(exc)}).",
                    warn=True,
                )

        # estimate transit duration
        pars = dict(**myfit.parameters)
        times = np.linspace(np.min(myfit.time), np.max(myfit.time), 1000)
        data_highres = transit(times, pars)
        dt = np.diff(times).mean()
        durs = []
        for r in range(1000):
            # randomize parameters
            for k in myfit.errors:
                pars[k] = np.random.normal(myfit.parameters[k], myfit.errors[k])

            data = transit(times, pars)
            tmask = data < 1
            durs.append(tmask.sum() * dt)

        plot_final_lightcurve(myfit, data_highres, pDict['pName'], exotic_infoDict['save'], exotic_infoDict['date'])

        if fitsortext == 1:
            observing_background_series = build_observing_background_series(
                psf_data,
                aper_data,
                photometry_info,
                len(exotic_infoDict['comp_stars']),
            )
            plot_obs_stats(myfit, exotic_infoDict['comp_stars'], psf_data, obs_stats_sort_index,
                           obs_stats_keep_mask, pDict['pName'],
                           exotic_infoDict['save'], exotic_infoDict['date'],
                           relative_flux_mask=None,
                           background_series=observing_background_series)

        #######################################################################
        # print final extracted planetary parameters
        #######################################################################

        log_info("\n*********************************************************")
        log_info("FINAL PLANETARY PARAMETERS\n")
        log_info(f"          Mid-Transit Time [BJD_TDB]: {round_to_2(myfit.parameters['tmid'], myfit.errors['tmid'])} +/- {round_to_2(myfit.errors['tmid'])}")
        log_info(f"  Radius Ratio (Planet/Star) [Rp/R*]: {round_to_2(myfit.parameters['rprs'], myfit.errors['rprs'])} +/- {round_to_2(myfit.errors['rprs'])}")
        log_info(f"           Transit depth [(Rp/R*)^2]: {round_to_2(100. * (myfit.parameters['rprs'] ** 2.))} +/- {round_to_2(100. * 2. * myfit.parameters['rprs'] * myfit.errors['rprs'])} [%]")
        log_info(f"           Orbital Inclination [inc]: {round_to_2(myfit.parameters['inc'], myfit.errors['inc'])} +/- {round_to_2(myfit.errors['inc'])}")
        ars_text = format_parameter_with_error(myfit.parameters.get('ars'), myfit.errors.get('ars'))
        if ars_text is not None:
            log_info(f" Ratio of Distance to Stellar Radius [a/Rs]: {ars_text}")
        impact_parameter, impact_error = fit_impact_parameter_value_error(myfit)
        impact_text = format_parameter_with_error(impact_parameter, impact_error)
        if impact_text is not None:
            log_info(f"                 Impact Parameter [b]: {impact_text}")
        if getattr(myfit, 'airmass_fit_skipped', False):
            log_info(f"                 Airmass correction: {myfit.airmass_correction_note}")
        else:
            log_info(f"               Airmass coefficient 1: {round_to_2(myfit.parameters['a1'], myfit.errors['a1'])} +/- {round_to_2(myfit.errors['a1'])}")
            log_info(f"               Airmass coefficient 2: {round_to_2(myfit.parameters['a2'], myfit.errors['a2'])} +/- {round_to_2(myfit.errors['a2'])}")
        transit_qc = getattr(myfit, 'transit_qc', None)
        if transit_qc:
            residual_scatter = transit_qc.get('residual_scatter', np.nan)
            if np.isfinite(residual_scatter):
                log_info(f"Residual scatter around full model fit: {residual_scatter * 100.0:.4f}%")
            qc_status = str(transit_qc.get('status', 'unknown')).upper()
            qc_summary = transit_qc.get('summary')
            if qc_summary:
                log_info(f"                Transit detection QC: {qc_status} - {qc_summary}")
            else:
                log_info(f"                Transit detection QC: {qc_status}")
            if np.isfinite(transit_qc.get('deviation_from_expected_value', np.nan)):
                log_info(
                    f"      Deviation From Expected Value: {transit_qc['deviation_from_expected_value']:.2f} / 1.00"
                )
            if np.isfinite(transit_qc.get('tmid_deviation_sigma', np.nan)):
                log_info(
                    f"          Expected-value Tmid sigma: {transit_qc['tmid_deviation_sigma']:.2f}"
                )
            if np.isfinite(transit_qc.get('tmid_deviation_minutes', np.nan)):
                log_info(
                    f"         Expected-value Tmid offset: {transit_qc['tmid_deviation_minutes']:.2f} minutes"
                )
            if np.isfinite(transit_qc.get('tmid_deviation_threshold_minutes', np.nan)):
                log_info(
                    "      Expected-value Tmid QC window: "
                    f"{transit_qc['tmid_deviation_threshold_minutes']:.2f} minutes"
                )
            if np.isfinite(transit_qc.get('rprs_deviation_sigma', np.nan)):
                log_info(
                    f"         Expected-value Rp/R* sigma: {transit_qc['rprs_deviation_sigma']:.2f}"
                )
            if np.isfinite(transit_qc.get('ktmf_metric', np.nan)):
                log_info(f"                             KTMF: {transit_qc['ktmf_metric']:.2f} / 5.00")
            for contribution in transit_qc.get('ktmf_contributions', []):
                log_info(f"      {format_ktmf_contribution(contribution)}")
        if fitsortext == 1:
            if np.isfinite(photometry_info.get('calibration_field_score', np.inf)):
                log_info(f"        Comparison-Star Field Score: {round_to_2(100. * photometry_info['calibration_field_score'])} %")
            display_aperture, display_annulus = reported_photometry_aperture_radii(photometry_info)
            adaptive_summary = photometry_info.get('adaptive_summary')
            if photometry_info['min_aperture'] >= 0:
                log_info(f"                Best Comparison Star: #{bestCompStar} - {comp_coords}")
            else:
                log_info("                 Best Comparison Star: None")
            if photometry_info['min_aperture'] == 0:
                log_info("                       Optimal Method: PSF photometry")
            else:
                if adaptive_summary is not None:
                    log_info(f"                    Optimal Aperture: {abs(display_aperture):.2f} +/- {adaptive_summary['aperture_std']:.2f} px")
                    log_info(f"                     Optimal Annulus: {display_annulus:.2f} +/- {adaptive_summary['annulus_std']:.2f} px")
                    log_info(f"              Adaptive Aperture Scale: {adaptive_summary['aperture_sigma']:.2f} sigma")
                    log_info(f"               Adaptive Annulus Scale: {adaptive_summary['annulus_sigma']:.2f} sigma")
                    log_info(f"                     Aperture Range: {adaptive_summary['aperture_min']:.2f} to {adaptive_summary['aperture_max']:.2f} px")
                    log_info(f"                      Annulus Range: {adaptive_summary['annulus_min']:.2f} to {adaptive_summary['annulus_max']:.2f} px")
                else:
                    log_info(f"                    Optimal Aperture: {abs(np.round(display_aperture, 2))}")
                    log_info(f"                     Optimal Annulus: {np.round(display_annulus, 2)}")
        log_info(f"              Transit Duration [day]: {round_to_2(np.mean(durs), np.std(durs))} +/- {round_to_2(np.std(durs))}")
        log_info("*********************************************************")

        ##########
        # SAVE DATA
        ##########

        selected_triangle_source_dir = (
            photometry_info.get('selected_fit_final_output_dir')
            if reuse_selected_final_model
            else None
        )
        save_final_triangle_plot(
            myfit,
            exotic_infoDict['save'],
            pDict['pName'],
            exotic_infoDict['date'],
            source_dir=selected_triangle_source_dir,
        )

        if vsp_params:
            AIDoutput_files = AIDOutputFiles(myfit, pDict, exotic_infoDict, auid, chart_id, vsp_params)
        output_files = OutputFiles(myfit, pDict, exotic_infoDict, durs)
        error_txt = "\n\tPlease report this issue on the Exoplanet Watch Slack Channel in #data-reductions."

        try:
            phase = get_phase(myfit.time, pDict['pPer'], myfit.parameters['tmid'])
            output_files.final_lightcurve(phase)
        except Exception as e:
            log_info(f"\nError: Could not create FinalLightCurve.csv. {error_txt}\n\t{e}", error=True)
        try:
            if fitsortext == 1:
                display_aperture, display_annulus = reported_photometry_aperture_radii(photometry_info)
                output_files.final_planetary_params(phot_opt=True, vsp_params=vsp_params,
                                                    comp_star=bestCompStar, comp_coords=comp_coords,
                                                    min_aper=np.round(display_aperture, 2),
                                                    min_annul=np.round(display_annulus, 2),
                                                    adaptive_summary=photometry_info.get('adaptive_summary'),
                                                    photometry_info=photometry_info,
                                                    publish_to_root=True)
            else:
                output_files.final_planetary_params(
                    phot_opt=False,
                    vsp_params=vsp_params,
                    publish_to_root=True,
                )
        except Exception as e:
            log_info(f"\nError: Could not create FinalParams.json. {error_txt}\n\t{e}", error=True)
        try:
            if bestCompStar:
                exotic_infoDict['phot_comp_star'] = save_comp_ra_dec(wcs_file, ra_wcs, dec_wcs, comp_coords)
            aavso_photometry_info = photometry_info if fitsortext == 1 else None
            aavso_frame_filtering_info = None
            aavso_astrometry_info = None
            aavso_bad_pixel_info = None
            if fitsortext == 1:
                aavso_frame_filtering_info = {
                    'initial_frame_count': precheck_inputfile_count,
                    'after_missing_wcs_filter_frame_count': post_wcs_inputfile_count,
                    'final_prephotometry_frame_count': post_pointing_inputfile_count,
                    'ignore_header_wcs': ignore_header_wcs,
                    'bad_wcs_threshold_percent': (
                        100.0 * bad_wcs_threshold_fraction
                        if np.isfinite(bad_wcs_threshold_fraction)
                        else np.nan
                    ),
                    'pointing_rejection_sigma': pointing_rejection_sigma,
                    'dropped_missing_wcs_files': dropped_wcs_files,
                    'dropped_pointing_files': dropped_pointing_files,
                }
                aavso_astrometry_info = {
                    'wcs_file': str(wcs_file) if wcs_file else None,
                    'coordinate_source': 'wcs' if wcs_file else 'input_pixels',
                    'ignore_header_wcs': ignore_header_wcs,
                    'plate_solution_option': exotic_infoDict.get('plate_opt'),
                    'target_input_pixel': exotic_infoDict.get('tar_coords'),
                    'comparison_input_pixels': exotic_infoDict.get('comp_stars'),
                    'target_ra_dec_deg': ra_dec_tar,
                    'comparison_ra_dec_deg': ra_dec_wcs,
                    'catalog_ra_dec_deg': [pDict.get('ra'), pDict.get('dec')],
                    'gaia_distance_pc': pDict.get('dist'),
                    'proper_motion_ra_mas_yr': pDict.get('pm_ra'),
                    'proper_motion_dec_mas_yr': pDict.get('pm_dec'),
                }
                aavso_bad_pixel_info = {
                    'enabled': bool(detect_bad_pixels_before_photometry),
                    'detected': bad_pixel_reference is not None,
                }
                if bad_pixel_reference is not None:
                    bad_pixel_mask = np.asarray(bad_pixel_reference.get('mask'), dtype=bool)
                    aavso_bad_pixel_info.update({
                        'bad_pixel_count': int(np.count_nonzero(bad_pixel_mask)),
                        'frame_count': bad_pixel_reference.get('frame_count'),
                        'required_count': bad_pixel_reference.get('required_count'),
                        'minimum_fraction': bad_pixel_reference.get('minimum_fraction'),
                        'counts_path': bad_pixel_reference.get('counts_path'),
                        'mask_path': bad_pixel_reference.get('mask_path'),
                    })
            output_files.aavso(
                exotic_infoDict['phot_comp_star'],
                goodAirmasses,
                ld0,
                ld1,
                ld2,
                ld3,
                epw_md5,
                photometry_info=aavso_photometry_info,
                astrometry_info=aavso_astrometry_info,
                frame_filtering_info=aavso_frame_filtering_info,
                bad_pixel_info=aavso_bad_pixel_info,
            )
        except Exception as e:
            log_info(f"\nError: Could not create AAVSO.txt. {error_txt}\n\t{e}", error=True)
        try:
            if vsp_params:
                AIDoutput_files.aavso()
        except Exception as e:
            log_info(f"\nError: Could not create AID_AAVSO.txt. {error_txt}\n\t{e}", error=True)
        try:
            output_files.plate_status(plateStatus)
        except Exception as e:
            log_info(f"\nError: Could not create plate_status.csv. {error_txt}\n\t{e}", error=True)

        log_info("Output Files Saved")

        log_info("\n************************")
        log_info("End of Reduction Process")
        log_info("************************")

        log_info("\n\n************************")
        log_info("EXOTIC has successfully run!!!")
        log_info("It is now safe to close this window.")
        log_info("************************")

        log.debug("Stopped ...")


def main():
    global _UNHANDLED_EXCEPTION_LOGGED

    _UNHANDLED_EXCEPTION_LOGGED = False
    configure_runtime_logging()
    install_exception_hooks()

    try:
        return _main_impl()
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as exc:
        _handle_unhandled_exception(type(exc), exc, exc.__traceback__)
        raise
    finally:
        cancel_runtime_traceback_watchdog()


def cli():
    global _UNHANDLED_EXCEPTION_LOGGED

    _UNHANDLED_EXCEPTION_LOGGED = False
    configure_runtime_logging()
    install_exception_hooks()

    try:
        return main()
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as exc:
        _handle_unhandled_exception(type(exc), exc, exc.__traceback__)
        raise


if __name__ == "__main__":
    raise SystemExit(cli())
