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
import json
import hashlib
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
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
from photutils.aperture import CircularAperture
import re
import requests
# scipy imports
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
from scipy.ndimage import binary_erosion, gaussian_filter
from skimage.registration import phase_cross_correlation
from skimage.transform import SimilarityTransform
# error handling for scraper
from tenacity import retry, stop_after_delay
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
try:  # nea
    from .api.nea import NASAExoplanetArchive
except ImportError:  # package import
    from api.nea import NASAExoplanetArchive
try:  # output files
    from output_files import OutputFiles, AIDOutputFiles, save_comp_star_calibration_summary
except ImportError:  # package import
    from .output_files import OutputFiles, AIDOutputFiles, save_comp_star_calibration_summary
try:
    from plate_status import PlateStatus
except ImportError:
    from .plate_status import PlateStatus
try:  # plots
    from plots import plot_fov, plot_centroids, plot_obs_stats, plot_final_lightcurve, plot_flux, \
        plot_stellar_variability, plot_variable_residuals, plot_comp_star_pairwise_matrix, \
        plot_comp_star_calibration_series, plot_comp_star_suitability
except ImportError:  # package import
    from .plots import plot_fov, plot_centroids, plot_obs_stats, plot_final_lightcurve, plot_flux, \
        plot_stellar_variability, plot_variable_residuals, plot_comp_star_pairwise_matrix, \
        plot_comp_star_calibration_series, plot_comp_star_suitability
try:  # tools
    from utils import round_to_2, user_input
except ImportError: # package import
    from .utils import round_to_2, user_input
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
_mid_transit_warning_reported = False
RELATIVE_FLUX_MAX = 2.0
AIRMASS_FLAT_RANGE_THRESHOLD = 0.05


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


def log_info(string, warn=False, error=False):
    if error:
        print(f"\033[31m {string}\033[0m")
    elif warn:
        print(f"\033[34m {string}\033[0m")
    else:
        print(string)
    log.debug(string)
    return True


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
    return np.isfinite(relative_flux) & np.less_equal(relative_flux, max_relative_flux)


def is_fast_aperture_mask_enabled(config_value):
    if config_value is None:
        return True
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

    log_info("Warning: Invalid 'Fast Aperture Mask (y/n)' value; using fast mode.", warn=True)
    return True


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


def apply_vertical_flux_normalization_bound(prior, bounds, flux_values, disabled):
    finite_flux = np.asarray(flux_values, dtype=float)
    finite_flux = finite_flux[np.isfinite(finite_flux) & (finite_flux > 0)]
    baseline_guess = 1.0 if finite_flux.size == 0 else float(np.nanmedian(finite_flux))
    baseline_guess = float(np.clip(baseline_guess, 0.95, 1.05))

    prior['a0'] = baseline_guess
    prior['a1'] = baseline_guess

    if not disabled:
        bounds['a0'] = [0.95, 1.05]


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

def sigma_clip(ogdata, sigma=3, dt=21, po=2):
    nanmask = np.isnan(ogdata)

    if po < dt <= len(ogdata[~nanmask]):
        mdata = savgol_filter(ogdata[~nanmask], window_length=dt, polyorder=po)
        # mdata = median_filter(ogdata[~nanmask], dt)
        res = ogdata[~nanmask] - mdata
        # Vectorized bootstrap estimate avoids Python-loop overhead in tight runs.
        sample_size = min(25, res.size)
        bootstrap_samples = np.random.choice(res, size=(100, sample_size), replace=True)
        std = bn.nanmedian(bn.nanstd(bootstrap_samples, axis=1))
        # std = np.nanstd(res) # biased from large outliers
        sigmask = np.abs(res) > sigma * std
        nanmask[~nanmask] = sigmask

    return nanmask


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

    if plate_opt == 'y' and not rt:
        wcs_file = get_wcs(fits_file, save_directory, use_nextastro_astrometry=use_nextastro_astrometry, ra=ra, dec=dec, pixel_scale=pixel_scale)
    if ignore_header_wcs:
        if wcs_file:
            log_info("Ignoring FITS header WCS for alignment and using the legacy image-to-image alignment path.")
        else:
            log_info("Ignoring FITS header WCS and using the legacy image-to-image alignment path.")
        return wcs_file
    if not wcs_file:
        if search_wcs(fits_file).is_celestial:
            log_info("Your FITS files have WCS (World Coordinate System) information in their headers. "
                     "EXOTIC will proceed to use these. "
                     "NOTE: If you do not trust your WCS coordinates, "
                     "please restart EXOTIC after enabling plate solutions via astrometry.net.")
            wcs_file = fits_file

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


def evaluate_celestial_wcs_coverage(inputfiles):
    missing_wcs_files = []
    for file_name in inputfiles:
        try:
            image_header = get_first_image_header(file_name)
            if not search_wcs_from_header(image_header).is_celestial:
                missing_wcs_files.append(str(file_name))
        except Exception:
            missing_wcs_files.append(str(file_name))

    total_files = len(inputfiles)
    all_have_celestial_wcs = total_files > 0 and len(missing_wcs_files) == 0
    return all_have_celestial_wcs, missing_wcs_files


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
    if missing_count > 0:
        preview = ", ".join(missing_wcs_files[:3])
        remainder = missing_count - 3
        if remainder > 0:
            preview = f"{preview}, ... (+{remainder} more)"
        log.debug(f"Files without usable celestial WCS: {preview}")

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
            suppress_fail_warning=True
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
        suppress_fail_warning=True
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
def get_ra_dec(header):
    wcs_header = WCS(header)
    xaxis = np.arange(header['NAXIS1'])
    yaxis = np.arange(header['NAXIS2'])
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


@retry(stop=stop_after_delay(30))
def nextastro_variability_test(comp_ra_dec):
    api_url = 'https://photometry.nextastro.org/variability_test'

    payload = [{'ra': float(ra), 'dec': float(dec)} for ra, dec in comp_ra_dec]
    log_info(f"NextAstro variability request JSON: {json.dumps(payload)}")
    result = requests.post(api_url, json=payload, timeout=30)
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
            log_info(f"\nWarning: NextAstro variability server check failed ({e}). "
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
                        'error': star_info['error']
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
    _TRANSFORM_REFERENCE_IMAGE = load_image_data(reference_file)
    _TRANSFORM_REFERENCE_CACHE = None


def transformation_task_with_cached_reference(i, file_name):
    image_data = load_image_data(file_name)
    return i, transformation(image_data, file_name, report_failure=False, reference_image=_TRANSFORM_REFERENCE_IMAGE)


MAX_MULTIPROCESS_TRANSFORM_WORKERS = 8

# Automatic aperture-grid tuning constants (in PSF sigma units)
APERTURE_SIGMA_MIN = 1.5
APERTURE_SIGMA_MAX = 6.0
ANNULUS_SIGMA_MIN = 6.0
ANNULUS_SIGMA_MAX = 15.0
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


def log_finding_transformation_progress(i, total_jobs, file_name, use_multiprocess_progress):
    if use_multiprocess_progress:
        completed = i + 1
        if completed == total_jobs or completed % 10 == 0:
            log_info(f"Multiprocessing finding transformations progress: {completed}/{total_jobs}")
        return

    display_file_name = _display_filename(file_name)
    sys.stdout.write(f"Finding transformation {i + 1} of {total_jobs} : {display_file_name}\n")
    log.debug(f"Finding transformation {i + 1} of {total_jobs} : {display_file_name}\n")
    sys.stdout.flush()


def get_img_scale(hdr, wcs_file, pixel_init):
    if wcs_file:
        wcs_hdr = fits.getheader(wcs_file)
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


def centroid_offset_matches_reference(psf_a, psf_b, expected_dx, expected_dy, tolerance=1):
    x_values = [psf_a[0], psf_b[0]]
    y_values = [psf_a[1], psf_b[1]]
    if not np.all(np.isfinite(x_values + y_values)):
        return False

    return (
        expected_dx - tolerance <= abs(int(psf_a[0]) - int(psf_b[0])) <= expected_dx + tolerance
        and expected_dy - tolerance <= abs(int(psf_a[1]) - int(psf_b[1])) <= expected_dy + tolerance
    )


# Method fits a 2D gaussian function that matches the star_psf to the star image and returns its pixel coordinates
def fit_centroid(data, pos, starIndex, psf_function=gaussian_psf, box=15, weightedcenter=True, fast_mode=False):
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
                return moment_fit
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

        # override psf fit results with weighted centroid
        if weightedcenter:
            res.x[0] = wx
            res.x[1] = wy

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


# Method calculates the flux of the star (uses the skybg_phot method to do background sub)
def aperPhot(data, starIndex, xc, yc, r=5, dr=5, fast_mode=True):
    stage_start = perf_counter()
    try:
        # Check for invalid coordinates
        if np.isnan(xc) or np.isnan(yc):
            return 0, 0

        # Calculate background if dr > 0
        if dr > 0:
            bgflux, sigmabg, Nbg = skybg_phot(data, starIndex, xc, yc, r + 2, dr)
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


def skybg_phot(data, starIndex, xc, yc, r=10, dr=5, ptol=99, debug=False):
    # create a crude annulus to mask out bright background pixels
    # the box will not extend beyond the borders of the image
    image_height, image_width = data.shape
    xv, yv = mesh_box([xc, yc], np.round(r + dr), maxx=image_width, maxy=image_height)
    if xv.size == 0 or yv.size == 0:
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: empty sky background box for {xc:.1f}, {yc:.1f}."
                 f"\nCheck if star is present or close to border.")
        return np.nan, np.nan, 0

    r_inner2 = float(r) ** 2
    r_outer2 = float(r + dr) ** 2
    rv2 = (xv - xc) ** 2 + (yv - yc) ** 2
    mask = (rv2 > r_inner2) & (rv2 < r_outer2)
    if not np.any(mask):
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: empty sky background annulus for {xc:.1f}, {yc:.1f}."
                 f"\nCheck if star is present or close to border.")
        return np.nan, np.nan, 0

    annulus_pixels = np.asarray(data[yv, xv][mask], dtype=float)
    if annulus_pixels.size == 0:
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: no valid sky background pixels for {xc:.1f}, {yc:.1f}."
                 f"\nCheck if star is present or close to border.")
        return np.nan, np.nan, 0

    try:
        cutoff = np.nanpercentile(annulus_pixels, ptol)
    except (IndexError, ValueError):
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: IndexError, problem computing sky bg for {xc:.1f}, {yc:.1f}."
                 f"\nCheck if star is present or close to border.")
        return np.nan, np.nan, 0

    dat = np.array(data[yv, xv], dtype=float)
    dat[dat > cutoff] = np.nan  # ignore pixels brighter than percentile

    if debug:
        minb = data[yv, xv][mask].min()
        maxb = data[yv, xv][mask].mean() + 3 * data[yv, xv][mask].std()
        nanmask = np.nan * np.zeros(mask.shape)
        nanmask[mask] = 1
        bgsky = data[yv, xv] * nanmask
        cmed, _ = sigma_clipped_nanmedian(dat.flatten(), sigma=3.0, max_iters=3)
        amed, _ = sigma_clipped_nanmedian(bgsky.flatten(), sigma=3.0, max_iters=3)

        fig, ax = plt.subplots(2, 2, figsize=(9, 9))
        im = ax[0, 0].imshow(data[yv, xv], vmin=minb, vmax=maxb, cmap='inferno')
        ax[0, 0].set_title("Original Data")
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax[0, 0])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        ax[1, 0].hist(bgsky.flatten(), label=f'Sky Annulus ({np.nanmedian(bgsky):.1f}, {amed:.1f})',
                      alpha=0.5, bins=np.arange(minb, maxb))
        ax[1, 0].hist(dat.flatten(), label=f'Clipped ({np.nanmedian(dat):.1f}, {cmed:.1f})', alpha=0.5,
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
    dat_flat = dat.ravel()
    sky_median, sky_sigma = sigma_clipped_nanmedian(dat_flat, sigma=3.0, max_iters=3)
    return sky_median, sky_sigma, np.sum(mask)

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
    try:
        goodTimes = JDUTC_to_BJDTDB(non_bjd, ra=p_dict['ra'], dec=p_dict['dec'], lat=info_dict['lat'],
                                    longi=info_dict['long'], alt=info_dict['elev'])[0]
    except:
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
    min_std_dev = min(std_devs, key=lambda y: abs(std_devs[y]))

    return comp_stars[min_std_dev]


def stellar_variability(fit_lc_refs, fit_lc_best, comp_stars, vsp_comp_stars, vsp_ind, best_comp, save, s_name):
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
        Mc, Mc_err = comp_star['mag'], comp_star['error']

        info_comp = info_comps[comp_stars.index(comp_pos)]
        lc_fit = info_comp['fit_lc']
        mask_ref = info_comp['mask_ref']

        oot_scatter = np.std((lc_fit.data / lc_fit.airmass_model)[mask_ref])
        norm_flux_unc = oot_scatter * lc_fit.airmass_model[mask_ref]
        norm_flux_unc /= np.nanmedian(lc_fit.data[mask_ref])

        model = np.exp(lc_fit.parameters['a2'] * lc_fit.airmass_model[mask_ref])
        flux = lc_fit.data[mask_ref]
        detrended = flux / model

        Mt = Mc - (2.5 * np.log10(detrended))
        Mt_err = (Mc_err ** 2 + (-2.5 * norm_flux_unc / (detrended * np.log(10))) ** 2) ** 0.5
    except KeyError as e:
        log_info(f"Key error in processing stellar variability: {e}", warn=True)
        return []
    except Exception as e:
        log_info(f"Error in processing stellar variability: {e}", warn=True)
        return []

    try:
        vsp_params = [{
            'time': lc_fit.jd_times[mask_ref][i],
            'airmass': lc_fit.airmass[mask_ref][i],
            'mag': mt,
            'mag_err': Mt_err[i],
            'cname': vsp_auid_comp,
            'cmag': Mc,
            'pos': comp_pos
        } for i, mt in enumerate(Mt)]

        plot_stellar_variability(vsp_params, save, s_name, vsp_auid_comp)
    except Exception as e:
        log_info(f"Error in plotting or finalizing stellar variability data: {e}", warn=True)
        return []

    return vsp_params


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

    use_multiprocess_transform_precompute = should_use_multiprocess_transform_precompute(
        inputfiles, multiprocess_transformations, ignore_header_wcs=ignore_header_wcs
    )
    fallback_transforms = {}
    if use_multiprocess_transform_precompute:
        fallback_transforms = build_multiprocess_transformations(inputfiles, multiprocess_transformations)

    exotic_UIprevTPX = info_dict['tar_coords'][0]
    exotic_UIprevTPY = info_dict['tar_coords'][1]

    plateStatus.setCurrentFilename(inputfiles[0])
    wcs_file = check_wcs(inputfiles[0], info_dict['save'], info_dict['plate_opt'], rt=True,
                         use_nextastro_astrometry=use_nextastro_astrometry,
                         ra=p_dict.get('ra'), dec=p_dict.get('dec'), pixel_scale=info_dict.get('pixel_scale'),
                         ignore_header_wcs=ignore_header_wcs)
    comp_star = info_dict['comp_stars']
    tar_radec, comp_radec = None, []

    if wcs_file:
        wcs_header = fits.getheader(filename=wcs_file)

        ra_file, dec_file = get_ra_dec(wcs_header)
        tar_radec = (ra_file[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)],
                     dec_file[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)])

        ra = ra_file[int(comp_star[1])][int(comp_star[0])]
        dec = dec_file[int(comp_star[1])][int(comp_star[0])]

        comp_radec.append((ra, dec))

    target_and_comp_radec = None
    if tar_radec is not None and comp_radec:
        target_and_comp_radec = np.array([tar_radec, comp_radec[0]], dtype=float)

    first_image = fits.getdata(inputfiles[0])
    targ_sig_xy = fit_centroid(first_image, [exotic_UIprevTPX, exotic_UIprevTPY], 0)[3:5]

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
    for i, fileName in enumerate(inputfiles):
        plateStatus.setCurrentFilename(fileName)
        hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False,
                         ignore_missing_end=True)
        frame_fast_centroid = should_use_fast_centroid(i)

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

        if i == 0:
            firstImage = np.copy(imageData)

        log_finding_transformation_progress(
            i,
            len(inputfiles),
            fileName,
            use_multiprocess_transform_precompute,
        )

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

                psf_data['target'][i] = fit_centroid_or_warn_out_of_frame(
                    imageData,
                    [tx, ty],
                    0,
                    fast_mode=frame_fast_centroid,
                )
                psf_data['comp'][i] = fit_centroid_or_warn_out_of_frame(
                    imageData,
                    [cx, cy],
                    1,
                    fast_mode=frame_fast_centroid,
                )

                target_flux_change_ok = True
                comp_valid = True
                if projected_off_frame:
                    use_wcs_alignment = True
                elif i != 0:
                    target_flux_change_ok = fractional_flux_change_within_limit(
                        psf_data['target'][i][2],
                        psf_data['target'][i - 1][2],
                    )
                    comp_valid = (
                        centroid_offset_matches_reference(
                            psf_data['comp'][i],
                            psf_data['target'][i],
                            tar_comp_dist['comp'][0],
                            tar_comp_dist['comp'][1],
                        )
                        and fractional_flux_change_within_limit(
                            psf_data['comp'][i][2],
                            psf_data['comp'][i - 1][2],
                        )
                    )
                    use_wcs_alignment = target_flux_change_ok and comp_valid
                else:
                    tar_comp_dist['comp'][0] = abs(int(psf_data['comp'][0][0]) - int(psf_data['target'][0][0]))
                    tar_comp_dist['comp'][1] = abs(int(psf_data['comp'][0][1]) - int(psf_data['target'][0][1]))
                    use_wcs_alignment = True
            except Exception:
                use_wcs_alignment = False

        if not use_wcs_alignment:
            if i == 0:
                tform = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])
            else:
                tform = fallback_transforms[i] if i in fallback_transforms else transformation(imageData, fileName, reference_image=firstImage)

            transformed_coords = np.asarray(
                tform(np.array([[exotic_UIprevTPX, exotic_UIprevTPY], comp_star], dtype=float)),
                dtype=float,
            )
            tx, ty = transformed_coords[0]
            psf_data['target'][i] = fit_centroid_or_warn_out_of_frame(
                imageData,
                [tx, ty],
                0,
                fast_mode=frame_fast_centroid,
            )

            cx, cy = transformed_coords[1]
            psf_data['comp'][i] = fit_centroid_or_warn_out_of_frame(
                imageData,
                [cx, cy],
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

        tFlux = aperPhot(imageData, 0, psf_data['target'][i, 0], psf_data['target'][i, 1], aper, annulus,
                         fast_mode=fast_aperture_mask)[0]
        cFlux = aperPhot(imageData, 1, psf_data['comp'][i, 0], psf_data['comp'][i, 1], aper, annulus,
                         fast_mode=fast_aperture_mask)[0]
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
                   allow_mid_transit_range_warning=True, disable_vertical_flux_normalization=False):
    # remove outliers
    si = np.argsort(times)
    times_sorted = times[si]
    tflux_sorted = tFlux[si]
    cflux_sorted = cFlux[si]
    with np.errstate(divide='ignore', invalid='ignore'):
        flux_ratio_sorted = np.divide(tflux_sorted, cflux_sorted)

    has_reference_flux = not np.allclose(cflux_sorted, 1.0)
    if has_reference_flux:
        relative_flux_mask = relative_flux_filter_mask(flux_ratio_sorted)
        times_sorted = times_sorted[relative_flux_mask]
        tflux_sorted = tflux_sorted[relative_flux_mask]
        cflux_sorted = cflux_sorted[relative_flux_mask]
        flux_ratio_sorted = flux_ratio_sorted[relative_flux_mask]
        jd_times_sorted = jd_times[si][relative_flux_mask]
        airmass_sorted = airmass[si][relative_flux_mask]
    else:
        jd_times_sorted = jd_times[si]
        airmass_sorted = airmass[si]

    if len(times_sorted) <= 1:
        return None, None, None

    dt = np.mean(np.diff(times_sorted))
    ndt = int(25. / 24. / 60. / dt) * 2 + 1
    if ndt > len(times_sorted):
        ndt = int(len(times_sorted)/4) * 2 + 1
    filtered_data = sigma_clip(flux_ratio_sorted, sigma=3, dt=max(5, ndt))
    valid_mask = ~filtered_data

    arrayFinalFlux = flux_ratio_sorted[valid_mask]
    f1 = tflux_sorted[valid_mask]
    sigf1 = f1 ** 0.5
    f2 = cflux_sorted[valid_mask]
    sigf2 = f2 ** 0.5
    if np.sum(cFlux) == len(cFlux):
        arrayNormUnc = sigf1
    else:
        arrayNormUnc = np.sqrt((sigf1 / f2) ** 2 + (sigf2 * f1 / f2 ** 2) ** 2)
    arrayTimes = times_sorted[valid_mask]
    arrayJDTimes = jd_times_sorted[valid_mask]
    arrayAirmass = airmass_sorted[valid_mask]

    # remove nans
    nanmask = np.isnan(arrayFinalFlux) | np.isnan(arrayNormUnc) | np.isnan(arrayTimes) | np.isnan(
        arrayAirmass) | np.less_equal(arrayFinalFlux, 0) | np.less_equal(arrayNormUnc, 0)
    nanmask = nanmask | np.isinf(arrayFinalFlux) | np.isinf(arrayNormUnc) | np.isinf(arrayTimes) | np.isinf(
        arrayAirmass)

    if np.sum(~nanmask) <= 1:
        return None, None, None
    else:
        arrayFinalFlux = arrayFinalFlux[~nanmask]
        arrayNormUnc = arrayNormUnc[~nanmask]
        arrayTimes = arrayTimes[~nanmask]
        arrayJDTimes = arrayJDTimes[~nanmask]
        arrayAirmass = arrayAirmass[~nanmask]
        f1 = f1[~nanmask]
        f2 = f2[~nanmask]

    skip_airmass_fit = should_skip_airmass_fit(arrayAirmass)


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
    prior['tmid'] = pDict['midT'] + np.floor(arrayPhases).max() * prior['per']

    upper = prior['tmid'] + np.abs(25 * pDict['midTUnc'] + np.floor(arrayPhases).max() * 25 * pDict['pPerUnc'])
    lower = prior['tmid'] - np.abs(25 * pDict['midTUnc'] + np.floor(arrayPhases).max() * 25 * pDict['pPerUnc'])

    if upper > prior['tmid'] + 0.25 * prior['per']:
        upper = prior['tmid'] + 0.25 * prior['per']
    if lower < prior['tmid'] - 0.25 * prior['per']:
        lower = prior['tmid'] - 0.25 * prior['per']

    if (
        allow_mid_transit_range_warning
        and np.floor(arrayPhases).max() - np.floor(arrayPhases).min() == 0
    ):
        log_mid_transit_range_warning_once(arrayTimes, prior['tmid'])

    mybounds = {
        'rprs': [0, prior['rprs'] * 1.25],
        'tmid': [lower, upper],
        'inc': [prior['inc'] - 5, min(90, prior['inc'] + 5)],
    }
    apply_vertical_flux_normalization_bound(
        prior,
        mybounds,
        arrayFinalFlux,
        disable_vertical_flux_normalization,
    )
    if not skip_airmass_fit:
        mybounds['a2'] = [-1, 1]

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
        mode='lm'
    )
    annotate_airmass_fit(myfit, arrayAirmass, skip_airmass_fit)

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
                mode='lm'
            )
            annotate_airmass_fit(myfit, arrayAirmass, skip_airmass_fit)

    return myfit, f1, f2


def cheap_lightcurve_prescore(tFlux, cFlux, airmass):
    with np.errstate(divide='ignore', invalid='ignore'):
        flux_ratio = np.divide(tFlux, cFlux)

    finite_mask = np.isfinite(flux_ratio) & np.isfinite(airmass) & (flux_ratio > 0)
    if not np.allclose(cFlux, 1.0):
        finite_mask &= relative_flux_filter_mask(flux_ratio)
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
    times, tflux, cflux, airmass, ld, p_dict, jd_times, disable_vertical_flux_normalization = task
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
    )
    if myfit is None:
        return None, tflux_fit, cflux_fit

    res_std = myfit.residuals.std() / np.median(myfit.data)
    return {
        'myfit': myfit,
        'res_std': res_std,
    }, tflux_fit, cflux_fit


def normalize_flux_series(flux_values):
    flux_values = np.asarray(flux_values, dtype=float)
    normalized = np.full(flux_values.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(flux_values) & (flux_values > 0)
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
    ratio[ratio > RELATIVE_FLUX_MAX] = np.nan
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


def comparison_star_stability_summary(comp_flux_map, airmass):
    if not comp_flux_map:
        return {
            'pairwise_matrix': np.empty((0, 0), dtype=float),
            'comp_summaries': [],
            'field_score': np.inf,
            'best_comp_index': None,
            'best_comp_score': np.inf,
        }

    comp_keys = list(comp_flux_map.keys())
    normalized_flux_map = {key: normalize_flux_series(comp_flux_map[key]) for key in comp_keys}
    pairwise_matrix = np.full((len(comp_keys), len(comp_keys)), np.nan, dtype=float)
    comp_summaries = []

    for i, key in enumerate(comp_keys):
        normalized_flux = normalized_flux_map[key]
        self_score = cheap_lightcurve_prescore(normalized_flux, np.ones(normalized_flux.shape[0]), airmass)
        pairwise_scores = []
        pairwise_series = {}

        for j, other_key in enumerate(comp_keys):
            if i == j:
                continue
            other_flux = normalized_flux_map[other_key]
            score = cheap_lightcurve_prescore(normalized_flux, other_flux, airmass)
            pairwise_matrix[i, j] = score
            pairwise_series[f"vs {j + 1}"] = normalized_ratio_series(normalized_flux, other_flux)
            if np.isfinite(score):
                pairwise_scores.append(float(score))

        ensemble_flux = build_normalized_comp_ensemble(normalized_flux_map, key)
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
        })

    finite_comp_scores = [summary['aggregate_score'] for summary in comp_summaries if np.isfinite(summary['aggregate_score'])]
    field_score = float(np.nanmedian(finite_comp_scores)) if finite_comp_scores else np.inf
    best_comp_index = None
    best_comp_score = np.inf
    for summary in comp_summaries:
        if summary['aggregate_score'] < best_comp_score:
            best_comp_score = summary['aggregate_score']
            best_comp_index = summary['comp_index']

    return {
        'pairwise_matrix': pairwise_matrix,
        'comp_summaries': comp_summaries,
        'field_score': field_score,
        'best_comp_index': best_comp_index,
        'best_comp_score': best_comp_score,
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


def compute_star_aperture_grid(data, star_index, xc, yc, apertures, annuli, fast_mode=True):
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
                    bgflux, _, _ = skybg_phot(data, star_index, xc, yc, float(aperture_radius) + 2, float(annulus_width))
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
    )
    aper_data['target'][frame_index] = target_flux
    aper_data['target_bg'][frame_index] = target_bg

    for comp_idx in range(comp_star_count):
        ckey = f"comp{comp_idx + 1}"
        comp_flux, comp_bg = compute_star_aperture_grid(
            image_data,
            comp_idx + 1,
            psf_data[ckey][frame_index, 0],
            psf_data[ckey][frame_index, 1],
            frame_apertures,
            frame_annuli,
            fast_mode=fast_aperture_mask,
        )
        aper_data[ckey][frame_index] = comp_flux
        aper_data[f"{ckey}_bg"][frame_index] = comp_bg


def load_calibrated_reduction_image(file_name, generalDark, generalBias, generalFlat,
                                    demosaic_fmt, demosaic_out, demosaic_mult):
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
    return image_data


def _refined_sigma_grid(center, lower_bound, upper_bound, half_width, points):
    low = max(lower_bound, center - half_width)
    high = min(upper_bound, center + half_width)
    if high <= low:
        low, high = lower_bound, upper_bound
    return np.linspace(low, high, points)


def auto_tune_aperture_sigma_grid(coarse_apertures_sigma, coarse_annuli_sigma, coarse_aper_data, comp_star_count,
                                  subset_airmass, require_comp_star=True):
    best_candidate = None
    best_score = np.inf

    for a_idx, aperture_sigma in enumerate(coarse_apertures_sigma):
        for an_idx, annulus_sigma in enumerate(coarse_annuli_sigma):
            comp_flux_map = {
                f"comp{comp_idx + 1}": coarse_aper_data[f"comp{comp_idx + 1}"][:, a_idx, an_idx]
                for comp_idx in range(comp_star_count)
            }
            field_summary = comparison_star_stability_summary(comp_flux_map, subset_airmass)
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


def select_comparison_calibrated_photometry(psf_data, aper_data, apers, annuli, airmass, comp_stars, sigma):
    candidate_summaries = []
    comp_star_count = len(comp_stars)

    if comp_star_count == 0:
        return None

    psf_flux_map = {
        f"comp{comp_idx + 1}": 2 * np.pi * psf_data[f"comp{comp_idx + 1}"][:, 2]
        * psf_data[f"comp{comp_idx + 1}"][:, 3]
        * psf_data[f"comp{comp_idx + 1}"][:, 4]
        for comp_idx in range(comp_star_count)
    }
    psf_summary = comparison_star_stability_summary(psf_flux_map, airmass)
    psf_summary.update({
        'method': 'psf',
        'a': None,
        'an': None,
        'aper': 0.0,
        'annulus': float(15 * sigma),
    })
    candidate_summaries.append(psf_summary)

    for a_idx, aperture in enumerate(apers):
        for an_idx, annulus in enumerate(annuli):
            comp_flux_map = {
                f"comp{comp_idx + 1}": aper_data[f"comp{comp_idx + 1}"][:, a_idx, an_idx]
                for comp_idx in range(comp_star_count)
            }
            candidate_summary = comparison_star_stability_summary(comp_flux_map, airmass)
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
        comp_summaries.append(comp_summary)

    best_candidate['comp_summaries'] = comp_summaries
    best_candidate['method_label'] = method_label
    return best_candidate


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
                        help="Use multiprocessing when finding image transformations. "
                             "Provide an integer number of processes to use.")
    parser.add_argument('--multiprocess-lightcurve-fits',
                        type=int,
                        default=None,
                        help="Use multiprocessing while evaluating candidate lightcurve fits. "
                             "Provide an integer number of processes to use.")
    return parser.parse_args()


def main():
    # command line args
    args = parse_args()
    if args.multiprocess_transformations is not None and args.multiprocess_transformations < 1:
        raise ValueError("--multiprocess-transformations requires an integer greater than 0.")
    if args.multiprocess_lightcurve_fits is not None and args.multiprocess_lightcurve_fits < 1:
        raise ValueError("--multiprocess-lightcurve-fits requires an integer greater than 0.")

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
            for file in inputfiles:
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

            ld, ld0, ld1, ld2, ld3 = get_ld_values(pDict, exotic_infoDict)

            # check for EPW_MD5 checksum
            if 'EPW_MD5' in header:
                epw_md5 = header['EPW_MD5']

            si = np.argsort(times)
            times = np.array(times)[si]
            jd_times = np.array(jd_times)[si]
            inputfiles = np.array(inputfiles)[si]
            
            exotic_UIprevTPX = exotic_infoDict['tar_coords'][0]
            exotic_UIprevTPY = exotic_infoDict['tar_coords'][1]

            # fit target in the first image and use it to determine aperture and annulus range
            inc = 0
            for ifile in inputfiles:
                plateStatus.setCurrentFilename(ifile)
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
            plateStatus.setCurrentFilename(inputfiles[0])

            # For astrometry hints, prioritize coordinates explicitly provided by the user
            # (from inits.json / CLI) over values scraped from NASA Exoplanet Archive.
            hint_ra = userpDict.get('ra', pDict.get('ra'))
            hint_dec = userpDict.get('dec', pDict.get('dec'))
            ignore_header_wcs = should_ignore_header_wcs(exotic_infoDict.get('ignore_header_wcs'))

            wcs_file = check_wcs(inputfiles[0], exotic_infoDict['save'], exotic_infoDict['plate_opt'],
                                 use_nextastro_astrometry=args.use_nextastro_astrometry,
                                 ra=hint_ra, dec=hint_dec, pixel_scale=exotic_infoDict.get('pixel_scale'),
                                 ignore_header_wcs=ignore_header_wcs)
            img_scale_str, img_scale = get_img_scale(header, wcs_file, exotic_infoDict['pixel_scale'])
            plateStatus.initializeComparisonStarCount(len(exotic_infoDict['comp_stars']))
            ra_dec_tar, ra_dec_wcs = None, []
            chart_id, vsp_comp_stars, vsp_list = None, None, []

            if wcs_file:
                if should_log_plate_solution_path(wcs_file):
                    log_info(f"\nHere is the path to your plate solution: {wcs_file}")
                wcs_header = fits.getheader(filename=wcs_file)
                ra_wcs, dec_wcs = get_ra_dec(wcs_header)

                exotic_UIprevTPX, exotic_UIprevTPY = check_target_pixel_wcs(exotic_UIprevTPX, exotic_UIprevTPY,
                                                                            pDict, ra_wcs, dec_wcs,
                                                                            fits.getdata(inputfiles[0]),
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
                # Build RA/Dec for comp after list is finalized (avoid off by one issues, etc
                ra_dec_wcs = build_comp_ra_dec(ra_wcs, dec_wcs, exotic_infoDict['comp_stars'])
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
            use_adaptive_apertures = is_adaptive_aperture_mode_enabled(
                exotic_infoDict.get('use_adaptive_apertures', False)
            )

            for i, coord in enumerate(exotic_infoDict['comp_stars']):
                ckey = f"comp{i + 1}"
                if coord in vsp_list:
                    vsp_num.append(i)
                psf_data[ckey] = np.zeros((len(inputfiles), 7))
                tar_comp_dist[ckey] = np.zeros(2)

            coarse_tune_frames = min(len(inputfiles), APERTURE_AUTOTUNE_MAX_FRAMES)
            if len(inputfiles) >= APERTURE_AUTOTUNE_MIN_FRAMES:
                coarse_tune_frames = max(APERTURE_AUTOTUNE_MIN_FRAMES, coarse_tune_frames)
            coarse_apertures_sigma = np.linspace(APERTURE_SIGMA_MIN, APERTURE_SIGMA_MAX, APERTURE_AUTOTUNE_COARSE_APER_POINTS)
            coarse_annuli_sigma = np.linspace(ANNULUS_SIGMA_MIN, ANNULUS_SIGMA_MAX, APERTURE_AUTOTUNE_COARSE_ANNULUS_POINTS)
            log_info(
                "Automatic aperture tuning enabled: "
                f"coarse_grid={len(coarse_apertures_sigma)}x{len(coarse_annuli_sigma)}, "
                f"coarse_frames={coarse_tune_frames}."
            )

            sigma = None
            coarse_aperture_values = None
            coarse_annulus_values = None
            aperture_values = None
            annulus_values = None
            apers = None
            annuli = None
            aperture_grid_tuned = False
            coarse_aper_data = initialize_aperture_data_store(
                coarse_tune_frames,
                len(coarse_apertures_sigma),
                len(coarse_annuli_sigma),
                comp_star_count,
            )
            aper_data = None
            coarse_frame_cache = [None] * coarse_tune_frames

            use_multiprocess_transform_precompute = should_use_multiprocess_transform_precompute(
                inputfiles, args.multiprocess_transformations, ignore_header_wcs=ignore_header_wcs
            )
            fallback_transforms = {}
            if use_multiprocess_transform_precompute:
                fallback_transforms = build_multiprocess_transformations(inputfiles, args.multiprocess_transformations)

            target_and_comp_radec = None
            if ra_dec_tar is not None and ra_dec_wcs:
                target_and_comp_radec = np.array([ra_dec_tar, *ra_dec_wcs], dtype=float)
            target_and_comp_pixels = np.array(
                [[exotic_UIprevTPX, exotic_UIprevTPY], *exotic_infoDict['comp_stars']],
                dtype=float,
            )
            fast_aperture_mask = is_fast_aperture_mask_enabled(exotic_infoDict.get('fast_aperture_mask'))
            if use_adaptive_apertures:
                log_info("Adaptive aperture scaling enabled: evaluating aperture candidates in PSF sigma units per frame.")

            # open files, calibrate, align, photometry
            reset_transform_timing_stats()
            reset_photometry_timing_stats()
            for i, fileName in enumerate(inputfiles):
                plateStatus.setCurrentFilename(fileName)
                hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False,
                                 ignore_missing_end=True)
                frame_fast_centroid = should_use_fast_centroid(i)

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

                if i == 0:
                    firstImage = np.copy(imageData)

                log_finding_transformation_progress(
                    i,
                    len(inputfiles),
                    fileName,
                    use_multiprocess_transform_precompute,
                )

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

                        psf_data['target'][i] = fit_centroid_or_warn_out_of_frame(
                            imageData,
                            [tx, ty],
                            0,
                            fast_mode=frame_fast_centroid,
                        )

                        # TODO: Add check for flux on target/comp stars relative to others in the field
                        # in case of cloudy data, large changes, etc.
                        target_flux_change_ok = True
                        if not projected_off_frame and i != 0:
                            target_flux_change_ok = fractional_flux_change_within_limit(
                                psf_data['target'][i][2],
                                psf_data['target'][i - 1][2],
                            )

                        comp_valid = True
                        for j in range(len(exotic_infoDict['comp_stars'])):
                            ckey = f"comp{j + 1}"

                            cx, cy = pix_x[j + 1], pix_y[j + 1]
                            psf_data[ckey][i] = fit_centroid_or_warn_out_of_frame(
                                imageData,
                                [cx, cy],
                                j + 1,
                                fast_mode=frame_fast_centroid,
                            )

                            if projected_off_frame:
                                continue
                            if i != 0:
                                comp_valid = comp_valid and (
                                    centroid_offset_matches_reference(
                                        psf_data[ckey][i],
                                        psf_data['target'][i],
                                        tar_comp_dist[ckey][0],
                                        tar_comp_dist[ckey][1],
                                    )
                                    and fractional_flux_change_within_limit(
                                        psf_data[ckey][i][2],
                                        psf_data[ckey][i - 1][2],
                                    )
                                )
                            else:
                                tar_comp_dist[ckey][0] = abs(int(psf_data[ckey][0][0]) - int(psf_data['target'][0][0]))
                                tar_comp_dist[ckey][1] = abs(int(psf_data[ckey][0][1]) - int(psf_data['target'][0][1]))

                        use_wcs_alignment = projected_off_frame or (target_flux_change_ok and comp_valid)
                    except Exception:
                        use_wcs_alignment = False

                if not use_wcs_alignment:
                    if i == 0:
                        tform = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])
                    else:
                        tform = fallback_transforms[i] if i in fallback_transforms else transformation(imageData, fileName, reference_image=firstImage)

                    transformed_coords = np.asarray(tform(target_and_comp_pixels), dtype=float)
                    tx, ty = transformed_coords[0]
                    psf_data['target'][i] = fit_centroid_or_warn_out_of_frame(
                        imageData,
                        [tx, ty],
                        0,
                        fast_mode=frame_fast_centroid,
                    )

                    for j, coord in enumerate(exotic_infoDict['comp_stars']):
                        ckey = f"comp{j + 1}"

                        cx, cy = transformed_coords[j + 1]
                        psf_data[ckey][i] = fit_centroid_or_warn_out_of_frame(
                            imageData,
                            [cx, cy],
                            j + 1,
                            fast_mode=frame_fast_centroid,
                        )

                        if i == 0:
                            tar_comp_dist[ckey][0] = abs(int(psf_data[ckey][0][0]) - int(psf_data['target'][0][0]))
                            tar_comp_dist[ckey][1] = abs(int(psf_data[ckey][0][1]) - int(psf_data['target'][0][1]))

                # aperture photometry
                if i == 0:
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

                if i < coarse_tune_frames:
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
                else:
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
            badmask = np.isnan(psf_data["target"][:, 0]) | (psf_data["target"][:, 0] == 0) | (aper_data["target"][:, 0, 0] == 0) | np.isnan(
                aper_data["target"][:, 0, 0])
            goodmask = ~badmask
            if np.sum(goodmask) == 0:
                log_info("No images to fit...check reference image for alignment (first image of sequence)")

            # convert to numpy arrays - strip all bad data
            times = times[goodmask]
            jd_times = jd_times[goodmask]
            airmass = np.array(airMassList)[goodmask]
            psf_data["target"] = psf_data["target"][goodmask]
            aper_data["target"] = aper_data["target"][goodmask]
            aper_data["target_bg"] = aper_data["target_bg"][goodmask]
            for j in range(len(exotic_infoDict['comp_stars'])):
                ckey = f"comp{j + 1}"
                psf_data[ckey] = psf_data[ckey][goodmask]
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
                'min_std': 100000,
                'min_aperture': None,
                'min_annulus': None,
                'calibration_field_score': np.inf,
                'selection_basis': 'target_fit',
            }

            comparison_calibration = None
            if target_driven_comp_selection:
                log_info("\nUsing target-driven comparison-star selection per optional_info setting.")
            else:
                comparison_calibration = select_comparison_calibrated_photometry(
                    psf_data,
                    aper_data,
                    apers,
                    annuli,
                    airmass,
                    exotic_infoDict['comp_stars'],
                    sigma_display,
                )

            if comparison_calibration is not None:
                log_info("\nCalibrating comparison stars before target fitting. Please wait.")
                log_info(f"Comparison-star field method: {comparison_calibration['method_label']}")
                log_info(f"Comparison-star field score: {comparison_calibration['field_score'] * 100.0:.4f}%")
                for summary in comparison_calibration['comp_summaries']:
                    aggregate_text = "n/a" if not np.isfinite(summary['aggregate_score']) else f"{summary['aggregate_score'] * 100.0:.4f}%"
                    ensemble_text = "n/a" if not np.isfinite(summary['ensemble_score']) else f"{summary['ensemble_score'] * 100.0:.4f}%"
                    pairwise_text = "n/a" if not np.isfinite(summary['pairwise_median_score']) else f"{summary['pairwise_median_score'] * 100.0:.4f}%"
                    selected_label = " [selected]" if summary['selected'] else ""
                    log_info(
                        f"  {summary['label']}{selected_label}: suitability={aggregate_text}, "
                        f"ensemble={ensemble_text}, pairwise_median={pairwise_text}, "
                        f"valid_pairs={summary['valid_pair_count']}"
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

                selected_comp_index = comparison_calibration['best_comp_index']
                selected_ckey = f"comp{selected_comp_index + 1}"
                selected_comp_coords = exotic_infoDict['comp_stars'][selected_comp_index]
                selected_min_aperture = 0 if comparison_calibration['method'] == 'psf' else comparison_calibration['aper']
                selected_min_annulus = comparison_calibration['annulus']

                if comparison_calibration['method'] == 'psf':
                    selected_target_flux = tFlux
                    selected_comp_flux = (
                        2 * np.pi * psf_data[selected_ckey][:, 2] * psf_data[selected_ckey][:, 3] * psf_data[selected_ckey][:, 4]
                    )
                else:
                    best_a = comparison_calibration['a']
                    best_an = comparison_calibration['an']
                    selected_target_flux = aper_data['target'][:, best_a, best_an]
                    selected_comp_flux = aper_data[selected_ckey][:, best_a, best_an]

                myfit, tFlux1, cFlux1 = fit_lightcurve(
                    times, selected_target_flux, selected_comp_flux, airmass, ld, pDict, jd_times,
                    disable_vertical_flux_normalization=disable_vertical_flux_normalization,
                )
                if myfit is not None:
                    res_std = myfit.residuals.std() / np.median(myfit.data)
                    photometry_info.update(best_fit_lc=myfit,
                                           comp_star_num=selected_comp_index + 1,
                                           comp_star_coords=selected_comp_coords,
                                           min_std=res_std,
                                           min_aperture=selected_min_aperture,
                                           min_annulus=selected_min_annulus,
                                           calibration_field_score=comparison_calibration['field_score'],
                                           selection_basis='comparison_field')

                    flux_values.update(flux_tar=tFlux1, flux_ref=cFlux1,
                                       flux_unc_tar=tFlux1 ** 0.5, flux_unc_ref=cFlux1 ** 0.5)

                    centroid_positions.update(x_targ=psf_data["target"][:, 0], y_targ=psf_data["target"][:, 1],
                                              x_ref=psf_data[selected_ckey][:, 0], y_ref=psf_data[selected_ckey][:, 1])

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
                                )
                                ref_flux[j] = {
                                    'myfit': vsp_fit,
                                    'pos': exotic_infoDict['comp_stars'][j]
                                }
                else:
                    log_info("Warning: Comparison-star calibration selected a photometry setup that failed target fitting."
                             " Falling back to target-driven photometry selection.", warn=True)

            if photometry_info['best_fit_lc'] is None:
                # Legacy fallback when comparison-star-only calibration cannot determine a usable setup.
                for j in range(len(exotic_infoDict['comp_stars'])):
                    ckey = f"comp{j + 1}"

                    cFlux = 2 * np.pi * psf_data[ckey][:, 2] * psf_data[ckey][:, 3] * psf_data[ckey][:, 4]
                    myfit, tFlux1, cFlux1 = fit_lightcurve(
                        times, tFlux, cFlux, airmass, ld, pDict, jd_times,
                        disable_vertical_flux_normalization=disable_vertical_flux_normalization,
                    )
                    res_std = np.inf

                    if myfit is not None:
                        for k in myfit.bounds.keys():
                            log.debug(f"  {k}: {myfit.parameters[k]:.6f}")

                        log.debug("The Residual Standard Deviation is: "
                                  f"{round(100 * myfit.residuals.std() / np.median(myfit.data), 6)}%")
                        log.debug(f"The Mean Squared Error is: {round(np.sum(myfit.residuals ** 2), 6)}\n")

                        res_std = myfit.residuals.std() / np.median(myfit.data)

                    if photometry_info['min_std'] > res_std and myfit is not None:
                        photometry_info.update(best_fit_lc=myfit,
                                               comp_star_num=j + 1, comp_star_coords=exotic_infoDict['comp_stars'][j],
                                               min_std=res_std, min_aperture=0, min_annulus=15 * sigma_display,
                                               selection_basis='target_fit')

                        flux_values.update(flux_tar=tFlux1, flux_ref=cFlux1,
                                           flux_unc_tar=tFlux1 ** 0.5, flux_unc_ref=cFlux1 ** 0.5)

                        centroid_positions.update(x_targ=psf_data["target"][:, 0], y_targ=psf_data["target"][:, 1],
                                                  x_ref=psf_data[ckey][:, 0], y_ref=psf_data[ckey][:, 1])

                    if j in vsp_num:
                        ref_flux[j] = {
                            'myfit': myfit,
                            'pos': exotic_infoDict['comp_stars'][j]
                        }

                log_info("\nComputing best comparison star, aperture, and sky annulus from the target lightcurve. Please wait.")

                candidate_jobs = []
                for a, aper in enumerate(apers):
                    for an, annulus in enumerate(annuli):
                        target_flux = aper_data['target'][:, a, an]

                        if not require_comp_star:
                            candidate_jobs.append({
                                'a': a,
                                'an': an,
                                'aper': aper,
                                'annulus': annulus,
                                'comp_index': None,
                                'ckey': None,
                                'mask': np.ones(target_flux.shape[0], dtype=bool),
                                'prescore': cheap_lightcurve_prescore(target_flux, np.ones(target_flux.shape[0]), airmass),
                            })

                        for j in range(len(exotic_infoDict['comp_stars'])):
                            ckey = f"comp{j + 1}"
                            comp_series = aper_data[ckey][:, a, an]
                            aper_mask = np.isfinite(comp_series)
                            comp_flux = comp_series[aper_mask]
                            candidate_jobs.append({
                                'a': a,
                                'an': an,
                                'aper': aper,
                                'annulus': annulus,
                                'comp_index': j,
                                'ckey': ckey,
                                'mask': aper_mask,
                                'prescore': cheap_lightcurve_prescore(target_flux[aper_mask], comp_flux, airmass[aper_mask]),
                            })

                finite_candidates = [c for c in candidate_jobs if np.isfinite(c['prescore'])]
                if finite_candidates:
                    finite_candidates.sort(key=lambda candidate: candidate['prescore'])
                    shortlist_count = max(50, int(0.35 * len(finite_candidates)))
                    shortlist = finite_candidates[:min(len(finite_candidates), shortlist_count)]
                else:
                    shortlist = []

                if not shortlist:
                    shortlist = candidate_jobs

                fit_tasks = []
                for candidate in shortlist:
                    candidate_mask = candidate['mask']
                    target_flux = aper_data['target'][:, candidate['a'], candidate['an']][candidate_mask]
                    if candidate['comp_index'] is None:
                        comp_flux = np.ones(target_flux.shape[0])
                    else:
                        comp_flux = aper_data[candidate['ckey']][:, candidate['a'], candidate['an']][candidate_mask]

                    fit_tasks.append((
                        times[candidate_mask],
                        target_flux,
                        comp_flux,
                        airmass[candidate_mask],
                        ld,
                        pDict,
                        jd_times[candidate_mask],
                        disable_vertical_flux_normalization,
                    ))

                fit_results = []
                if args.multiprocess_lightcurve_fits is not None and args.multiprocess_lightcurve_fits > 0:
                    log_info(f"Using multiprocessing for candidate lightcurve fits ({args.multiprocess_lightcurve_fits} processes).")
                    with ProcessPoolExecutor(max_workers=args.multiprocess_lightcurve_fits) as executor:
                        fit_results = list(executor.map(evaluate_lightcurve_candidate, fit_tasks))
                else:
                    fit_results = [evaluate_lightcurve_candidate(task) for task in fit_tasks]

                best_candidate = None
                for candidate, result in zip(shortlist, fit_results):
                    fit_meta, tFlux1, cFlux1 = result
                    if fit_meta is None:
                        continue

                    myfit = fit_meta['myfit']
                    res_std = fit_meta['res_std']

                    if photometry_info['min_std'] > res_std:
                        best_candidate = candidate
                        photometry_info.update(best_fit_lc=myfit,
                                               comp_star_num=(None if candidate['comp_index'] is None else candidate['comp_index'] + 1),
                                               comp_star_coords=(None if candidate['comp_index'] is None else exotic_infoDict['comp_stars'][candidate['comp_index']]),
                                               min_std=res_std,
                                               min_aperture=(-candidate['aper'] if candidate['comp_index'] is None else candidate['aper']),
                                               min_annulus=candidate['annulus'],
                                               selection_basis='target_fit')

                        flux_values.update(flux_tar=tFlux1, flux_ref=cFlux1,
                                           flux_unc_tar=tFlux1 ** 0.5, flux_unc_ref=cFlux1 ** 0.5)

                        x_ref_data = psf_data['target'][:, 0]
                        y_ref_data = psf_data['target'][:, 1]
                        if candidate['ckey'] is not None:
                            x_ref_data = psf_data[candidate['ckey']][:, 0]
                            y_ref_data = psf_data[candidate['ckey']][:, 1]

                        centroid_positions.update(x_targ=psf_data["target"][:, 0], y_targ=psf_data["target"][:, 1],
                                                  x_ref=x_ref_data, y_ref=y_ref_data)

                if best_candidate is not None and vsp_num:
                    best_a = best_candidate['a']
                    best_an = best_candidate['an']
                    best_target_flux = aper_data['target'][:, best_a, best_an]
                    for j in vsp_num:
                        ckey = f"comp{j + 1}"
                        aper_mask = np.isfinite(aper_data[ckey][:, best_a, best_an])
                        cFlux = aper_data[ckey][aper_mask][:, best_a, best_an]
                        vsp_fit, _, _ = fit_lightcurve(
                            times[aper_mask], best_target_flux[aper_mask], cFlux,
                            airmass[aper_mask], ld, pDict, jd_times[aper_mask],
                            disable_vertical_flux_normalization=disable_vertical_flux_normalization,
                        )
                        ref_flux[j] = {
                            'myfit': vsp_fit,
                            'pos': exotic_infoDict['comp_stars'][j]
                        }

            if require_comp_star and photometry_info['comp_star_num'] is None:
                log_info("Error: require_comp_star is enabled, but no valid comparison star could be selected.", error=True)
                return

            log_info("\n\n*********************************************")
            if np.isfinite(photometry_info['calibration_field_score']):
                log_info(f"Comparison-Star Field Score: {round(photometry_info['calibration_field_score'] * 100, 4)}%")
            if photometry_info['min_aperture'] == 0:  # psf
                log_info(f"Best Comparison Star: #{photometry_info['comp_star_num']}")
                log_info(f"Target-Fit Residual Scatter: {round(photometry_info['min_std'] * 100, 4)}%")
                log_info("Optimal Method: PSF photometry")
            elif photometry_info['min_aperture'] < 0:  # no comp star
                log_info("Best Comparison Star: None")
                log_info(f"Target-Fit Residual Scatter: {round(photometry_info['min_std'] * 100, 4)}%")
                log_info(f"Optimal Aperture: {abs(np.round(photometry_info['min_aperture'], 2))}")
                log_info(f"Optimal Annulus: {np.round(photometry_info['min_annulus'], 2)}")
            else:
                log_info(f"Best Comparison Star: #{photometry_info['comp_star_num']}")
                log_info(f"Target-Fit Residual Scatter: {round(photometry_info['min_std'] * 100, 4)}%")
                log_info(f"Optimal Aperture: {np.round(photometry_info['min_aperture'], 2)}")
                log_info(f"Optimal Annulus: {np.round(photometry_info['min_annulus'], 2)}")
            log_info("*********************************************\n")

            best_fit_lc = photometry_info['best_fit_lc']
            bestCompStar = photometry_info['comp_star_num']
            comp_coords = photometry_info['comp_star_coords']

            # save psf_data to disk for best comparison star
            if bestCompStar:
                np.savetxt(Path(exotic_infoDict['save']) / "temp" / "psf_data_comp.txt", psf_data[f"comp{bestCompStar}"],
                            header="#x_centroid, y_centroid, amplitude, sigma_x, sigma_y, rotation offset",
                            fmt="%.6f")

            # sigma clip
            si = np.argsort(best_fit_lc.time)
            dt = np.mean(np.diff(np.sort(best_fit_lc.time)))
            ndt = int(30. / 24. / 60. / dt) * 2 + 1  # ~30 minutes
            time_clip_mask = sigma_clip(best_fit_lc.data[si], sigma=3, dt=ndt)
            phase_clip_mask = np.zeros_like(time_clip_mask, dtype=bool)
            if hasattr(best_fit_lc, 'residuals') and hasattr(best_fit_lc, 'phase'):
                phase_clip_mask = phase_bin_sigma_clip(best_fit_lc.residuals[si], best_fit_lc.phase[si], sigma=3, bins=10)
            gi = ~(time_clip_mask | phase_clip_mask)  # good indexs
            phase_clip_removed = np.count_nonzero(phase_clip_mask & ~time_clip_mask)
            if phase_clip_removed:
                log_info(f"Removed {phase_clip_removed} phase-binned residual outlier(s) before final fit.")

            # Calculate the proper timeseries uncertainties from the residuals of the out-of-transit data
            OOT = (best_fit_lc.transit == 1)  # find out-of-transit portion of the lightcurve

            if sum(OOT) <= 1:
                OOTscatter = np.std(best_fit_lc.residuals)
                goodNormUnc = OOTscatter * best_fit_lc.airmass_model
                goodNormUnc = goodNormUnc / np.nanmedian(best_fit_lc.data)
                goodFluxes = best_fit_lc.data / np.nanmedian(best_fit_lc.data)
            else:
                OOTscatter = np.std((best_fit_lc.data / best_fit_lc.airmass_model)[OOT])  # calculate the scatter in the data
                goodNormUnc = OOTscatter * best_fit_lc.airmass_model  # scale this scatter back up by the airmass model and then adopt these as the uncertainties
                goodNormUnc = goodNormUnc / np.nanmedian(best_fit_lc.data[OOT])
                goodFluxes = best_fit_lc.data / np.nanmedian(best_fit_lc.data[OOT])

            if np.isnan(best_fit_lc.data).all():
                log_info("Error: No valid photometry data found.", error=True)
                return

            apply_lightcurve_mask(best_fit_lc, gi, sort_index=si)

            goodTimes = best_fit_lc.time
            goodFluxes = goodFluxes[si][gi]
            goodNormUnc = goodNormUnc[si][gi]
            goodAirmasses = best_fit_lc.airmass

            centroid_positions.update(x_targ=centroid_positions['x_targ'][si][gi],
                                      y_targ=centroid_positions['y_targ'][si][gi],
                                      x_ref=centroid_positions['x_ref'][si][gi],
                                      y_ref=centroid_positions['y_ref'][si][gi])

            flux_values.update(flux_tar=flux_values['flux_tar'][si][gi],
                               flux_ref=flux_values['flux_ref'][si][gi],
                               flux_unc_tar=flux_values['flux_unc_tar'][si][gi],
                               flux_unc_ref=flux_values['flux_unc_ref'][si][gi])

            relative_flux_mask = relative_flux_filter_mask(goodFluxes)
            if np.count_nonzero(relative_flux_mask) == 0:
                log_info("Error: No valid photometry data found after removing relative flux values above 2.", error=True)
                return

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


            if photometry_info['min_aperture'] == 0:
                opt_method = "PSF"
                # Calculate min_aper and min_annulus using the stdev
                stdev_fov = (psf_data['target'][:, 3] + psf_data['target'][:, 4]) * 0.5
                min_aper_fov = float(5 * stdev_fov.mean())
                min_annulus_fov = float(15 * stdev_fov.mean())
            else:
                opt_method = "Aperture"
                min_aper_fov = float(photometry_info['min_aperture'])
                min_annulus_fov = float(photometry_info['min_annulus'])
            
            plot_fov(photometry_info['min_aperture'], photometry_info['min_annulus'], sigma_display,
                     centroid_positions['x_targ'][0], centroid_positions['y_targ'][0],
                     centroid_positions['x_ref'][0], centroid_positions['y_ref'][0],
                     firstImage, img_scale_str, pDict['pName'], exotic_infoDict['save'], exotic_infoDict['date'], opt_method, min_aper_fov, min_annulus_fov)

            plot_centroids(centroid_positions['x_targ'], centroid_positions['y_targ'],
                           centroid_positions['x_ref'], centroid_positions['y_ref'],
                           goodTimes, pDict['pName'], exotic_infoDict['save'], exotic_infoDict['date'])

            plot_flux(goodTimes, flux_values['flux_tar'], flux_values['flux_unc_tar'],
                      flux_values['flux_ref'], flux_values['flux_unc_ref'],
                      goodFluxes, goodNormUnc, goodAirmasses, pDict['pName'], exotic_infoDict['save'],
                      exotic_infoDict['date'])

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
                                                     pDict['sName'])
                else:
                    vsp_params = stellar_variability(ref_flux, best_fit_lc, exotic_infoDict['comp_stars'],
                                                     vsp_comp_stars, vsp_num, bestCompStar - 1, exotic_infoDict['save'],
                                                     pDict['sName'])

            log_info("\n\nOutput File Saved")
        else:
            goodTimes, goodFluxes, goodNormUnc, goodAirmasses = [], [], [], []
            bestCompStar, comp_coords = None, None
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
                log_info("Error: No valid photometry data found after removing relative flux values above 2.", error=True)
                return

            goodTimes = goodTimes[relative_flux_mask]
            goodFluxes = goodFluxes[relative_flux_mask]
            goodNormUnc = goodNormUnc[relative_flux_mask]
            goodAirmasses = goodAirmasses[relative_flux_mask]

        # for k in myfit.bounds.keys():
        #     print(f"{myfit.parameters[k]:.6f} +- {myfit.errors[k]}")

        if args.photometry:
            log_info("\nPhotometric Extraction Complete.")
            return

        log_info("\n")
        log_info("****************************************")
        log_info("Fitting a Light Curve Model to Your Data")
        log_info("****************************************\n")

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
        prior['tmid'] = pDict['midT'] + np.floor(phase).max() * prior['per']
        upper = pDict['midT'] + 35 * pDict['midTUnc'] + np.floor(phase).max() * (pDict['pPer'] + 35 * pDict['pPerUnc'])
        lower = pDict['midT'] - 35 * pDict['midTUnc'] + np.floor(phase).max() * (pDict['pPer'] - 35 * pDict['pPerUnc'])

        # clip bounds so they're within 1 orbit
        if upper > prior['tmid'] + 0.25*prior['per']:
            upper = prior['tmid'] + 0.25*prior['per']
        if lower < prior['tmid'] - 0.25*prior['per']:
            lower = prior['tmid'] - 0.25*prior['per']

        if np.floor(phase).max() - np.floor(phase).min() == 0:
            log_info("Error: Estimated mid-transit not in observation range (check priors or observation time)", error=True)
            log_info(f"start:{np.min(goodTimes)}", error=True)
            log_info(f"  end:{np.max(goodTimes)}", error=True)
            log_info(f"prior:{prior['tmid']}", error=True)

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

        mybounds = {
            'rprs': [0, prior['rprs'] * 1.25],
            'tmid': [lower, upper],
            'inc': [prior['inc'] - 5, min(90, prior['inc'] + 5)],
        }
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

        # final light curve fit
        myfit = lc_fitter(goodTimes, goodFluxes, goodNormUnc, goodAirmasses, prior, mybounds, mode='ns')
        annotate_airmass_fit(myfit, goodAirmasses, skip_final_airmass_fit, note=airmass_skip_note)
        # myfit.dataerr *= np.sqrt(myfit.chi2 / myfit.data.shape[0])  # scale errorbars by sqrt(rchi2)
        # myfit.detrendederr *= np.sqrt(myfit.chi2 / myfit.data.shape[0])

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
            plot_obs_stats(myfit, exotic_infoDict['comp_stars'], psf_data, si, gi, pDict['pName'],
                           exotic_infoDict['save'], exotic_infoDict['date'],
                           relative_flux_mask=relative_flux_mask)

        #######################################################################
        # print final extracted planetary parameters
        #######################################################################

        log_info("\n*********************************************************")
        log_info("FINAL PLANETARY PARAMETERS\n")
        log_info(f"          Mid-Transit Time [BJD_TDB]: {round_to_2(myfit.parameters['tmid'], myfit.errors['tmid'])} +/- {round_to_2(myfit.errors['tmid'])}")
        log_info(f"  Radius Ratio (Planet/Star) [Rp/R*]: {round_to_2(myfit.parameters['rprs'], myfit.errors['rprs'])} +/- {round_to_2(myfit.errors['rprs'])}")
        log_info(f"           Transit depth [(Rp/R*)^2]: {round_to_2(100. * (myfit.parameters['rprs'] ** 2.))} +/- {round_to_2(100. * 2. * myfit.parameters['rprs'] * myfit.errors['rprs'])} [%]")
        log_info(f"           Orbital Inclination [inc]: {round_to_2(myfit.parameters['inc'], myfit.errors['inc'])} +/- {round_to_2(myfit.errors['inc'])}")
        if getattr(myfit, 'airmass_fit_skipped', False):
            log_info(f"                 Airmass correction: {myfit.airmass_correction_note}")
        else:
            log_info(f"               Airmass coefficient 1: {round_to_2(myfit.parameters['a1'], myfit.errors['a1'])} +/- {round_to_2(myfit.errors['a1'])}")
            log_info(f"               Airmass coefficient 2: {round_to_2(myfit.parameters['a2'], myfit.errors['a2'])} +/- {round_to_2(myfit.errors['a2'])}")
        log_info(f"                    Residual scatter: {round_to_2(100. * np.std(myfit.residuals / np.median(myfit.data)))} %")
        if fitsortext == 1:
            if np.isfinite(photometry_info.get('calibration_field_score', np.inf)):
                log_info(f"        Comparison-Star Field Score: {round_to_2(100. * photometry_info['calibration_field_score'])} %")
            if photometry_info['min_aperture'] >= 0:
                log_info(f"                Best Comparison Star: #{bestCompStar} - {comp_coords}")
            else:
                log_info("                 Best Comparison Star: None")
            if photometry_info['min_aperture'] == 0:
                log_info("                       Optimal Method: PSF photometry")
            else:
                log_info(f"                    Optimal Aperture: {abs(np.round(photometry_info['min_aperture'], 2))}")
                log_info(f"                     Optimal Annulus: {np.round(photometry_info['min_annulus'], 2)}")
        log_info(f"              Transit Duration [day]: {round_to_2(np.mean(durs), np.std(durs))} +/- {round_to_2(np.std(durs))}")
        log_info("*********************************************************")

        ##########
        # SAVE DATA
        ##########

        fig = myfit.plot_triangle()
        fig.savefig(Path(exotic_infoDict['save']) / "temp" /
                    f"Triangle_{pDict['pName']}_{exotic_infoDict['date']}.png")

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
                output_files.final_planetary_params(phot_opt=True, vsp_params=vsp_params,
                                                    comp_star=bestCompStar, comp_coords=comp_coords,
                                                    min_aper=np.round(photometry_info['min_aperture'], 2),
                                                    min_annul=np.round(photometry_info['min_annulus'], 2))
            else:
                output_files.final_planetary_params(phot_opt=False, vsp_params=vsp_params)
        except Exception as e:
            log_info(f"\nError: Could not create FinalParams.json. {error_txt}\n\t{e}", error=True)
        try:
            if bestCompStar:
                exotic_infoDict['phot_comp_star'] = save_comp_ra_dec(wcs_file, ra_wcs, dec_wcs, comp_coords)
            output_files.aavso(exotic_infoDict['phot_comp_star'], goodAirmasses, ld0, ld1, ld2, ld3, epw_md5)
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


if __name__ == "__main__":
    # configure logger for standalone execution
    logging.root.setLevel(logging.DEBUG)
    fileFormatter = logging.Formatter("%(asctime)s.%(msecs)03d [%(threadName)-12.12s] %(levelname)-5.5s  "
                                      "%(funcName)s:%(lineno)d - %(message)s", f"%Y-%m-%dT%H:%M:%S")
    fileHandler = TimedRotatingFileHandler(filename="exotic.log", when="midnight", backupCount=2)
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(fileFormatter)
    consoleFormatter = logging.Formatter("%(message)s")
    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setFormatter(consoleFormatter)
    consoleHandler.setLevel(logging.INFO)
    log.addHandler(fileHandler)
    log.addHandler(consoleHandler)
    main()
