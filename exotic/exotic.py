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
import hashlib
from time import sleep
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
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
# UTC to BJD converter import
from barycorrpy.utc_tdb import JDUTC_to_BJDTDB
import copy
# julian conversion imports
import dateutil.parser as dup
import imreg_dft as ird
from pathlib import Path
import pyvo as vo
import logging
from logging.handlers import TimedRotatingFileHandler
from matplotlib.animation import FuncAnimation
# Pyplot imports
import matplotlib.pyplot as plt
import numpy as np
# photometry
from photutils.aperture import CircularAperture
import pandas as pd
import re
import requests
# scipy imports
from scipy.optimize import least_squares
from scipy.stats import mode
from scipy.signal import savgol_filter
from scipy.ndimage import binary_erosion
from skimage.util import view_as_windows
from skimage.transform import SimilarityTransform
from skimage.color import rgb2gray
# error handling for scraper
from tenacity import retry, stop_after_delay
# color, color_demosaicing
from colour_demosaicing import demosaicing_CFA_Bayer_bilinear
# ########## EXOTIC imports ##########
try:  # light curve numerics
    from .api.elca import lc_fitter, binner, transit, get_phase
except ImportError:  # package import
    from api.elca import lc_fitter, binner, transit, get_phase
try:  # output files
    from inputs import Inputs, comparison_star_coords
except ImportError:  # package import
    from .inputs import Inputs, comparison_star_coords
try:  # ld
    from .api.ld import LimbDarkening, ld_re_punct_p
except ImportError:  # package import
    from api.ld import LimbDarkening, ld_re_punct_p
try:  # plate solution
    from .api.plate_solution import PlateSolution
except ImportError:  # package import
    from api.plate_solution import PlateSolution
try:  # nea
    from .api.nea import NASAExoplanetArchive
except ImportError:  # package import
    from api.nea import NASAExoplanetArchive
try:  # output files
    from output_files import OutputFiles, AIDOutputFiles
except ImportError:  # package import
    from .output_files import OutputFiles, AIDOutputFiles
try:
    from plate_status import PlateStatus
except ImportError:
    from .plate_status import PlateStatus
try:  # plots
    from plots import plot_fov, plot_centroids, plot_obs_stats, plot_final_lightcurve, plot_flux, \
        plot_stellar_variability, plot_variable_residuals
except ImportError:  # package import
    from .plots import plot_fov, plot_centroids, plot_obs_stats, plot_final_lightcurve, plot_flux, \
        plot_stellar_variability, plot_variable_residuals
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

def log_info(string, warn=False, error=False):
    if error:
        print(f"\033[31m {string}\033[0m")
    elif warn:
        print(f"\033[33m {string}\033[0m")
    else:
        print(string)
    log.debug(string)
    return True

# Initialze plate status log
plateStatus = PlateStatus(log_info)

def sigma_clip(ogdata, sigma=3, dt=21, po=2):
    nanmask = np.isnan(ogdata)

    if po < dt <= len(ogdata[~nanmask]):
        mdata = savgol_filter(ogdata[~nanmask], window_length=dt, polyorder=po)
        # mdata = median_filter(ogdata[~nanmask], dt)
        res = ogdata[~nanmask] - mdata
        std = np.nanmedian([np.nanstd(np.random.choice(res, 25)) for i in range(100)])
        # std = np.nanstd(res) # biased from large outliers
        sigmask = np.abs(res) > sigma * std
        nanmask[~nanmask] = sigmask

    return nanmask


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
            ra = ra.replace(' ', '').replace(':', ' ')
            dec = dec.replace(' ', '').replace(':', ' ')
            c = SkyCoord(ra + ' ' + dec, unit=(u.hourangle, u.deg))
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

    if check_all_standard_filters(ld, observed_filter):
        pass
    elif observed_filter['wl_min'] and observed_filter['wl_max']:
        custom_range(ld, observed_filter)
        ld.set_filter('N/A', "Custom", float(observed_filter['wl_min']), float(observed_filter['wl_max']))
    else:
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

def check_wcs(fits_file, save_directory, plate_opt, rt=False):
    wcs_file = None

    if plate_opt == 'y' and not rt:
        wcs_file = get_wcs(fits_file, save_directory)
    if not wcs_file:
        if search_wcs(fits_file).is_celestial:
            log_info("Your FITS files have WCS (World Coordinate System) information in their headers. "
                     "EXOTIC will proceed to use these. "
                     "NOTE: If you do not trust your WCS coordinates, "
                     "please restart EXOTIC after enabling plate solutions via astrometry.net.")
            wcs_file = fits_file

    return wcs_file


def search_wcs(file):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=FITSFixedWarning)
        header = fits.getheader(filename=file)
        return WCS(header)
        # return WCS(fits.open(file)[('SCI', 1)].header)


def get_wcs(file, directory=""):
    log_info("\nGetting the plate solution for your imaging file to translate pixel coordinates on the sky. "
             "\nPlease wait....")
    animate_toggle(True)
    wcs_obj = PlateSolution(file=file, directory=directory)
    wcs_file = wcs_obj.plate_solution()
    animate_toggle()
    return wcs_file


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


def check_target_pixel_wcs(input_x_pixel, input_y_pixel, info_dict, ra_list, dec_list, image_data, obs_time):
    """
    Verify the provided pixel coordinates match the target's right ascension and declination.
    """
    updated_ra, updated_dec = update_coordinates_with_proper_motion(info_dict, obs_time)

    calculated_y_pixel, calculated_x_pixel = deg_to_pix(updated_ra, updated_dec, ra_list, dec_list)

    centroid_x, centroid_y, sigma_x, sigma_y = get_psf_parameters(image_data, calculated_x_pixel, calculated_y_pixel)

    return check_coordinates(input_x_pixel, input_y_pixel, centroid_x, centroid_y, sigma_x, sigma_y,
                             calculated_x_pixel, calculated_y_pixel)


def get_psf_parameters(image_data, x_pixel, y_pixel):
    psf_data = fit_centroid(image_data, [x_pixel, y_pixel], 0)
    return psf_data[0], psf_data[1], psf_data[3], psf_data[4]


def check_coordinates(input_x_pixel, input_y_pixel, centroid_x, centroid_y, sigma_x, sigma_y,
                      calculated_x_pixel, calculated_y_pixel):
    while True:
        try:
            validate_pixel_coordinates(input_x_pixel, input_y_pixel, centroid_x, centroid_y, sigma_x, sigma_y)
            return input_x_pixel, input_y_pixel
        except ValueError:
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
    # radius = u.Quantity(20.0, u.arcsec)
    # # Query GAIA first to check for variability using the phot_variable_flag trait
    # gaia_result = gaia_query(sample, radius)
    # if not gaia_result:
    #     log_info("Warning: Your comparison star cannot be resolved in the Gaia star database; "
    #              "EXOTIC cannot check if it is variable or not. "
    #              "\nEXOTIC will still include this star in the reduction. "
    #              "\nPlease proceed with caution as we cannot check for stellar variability.\n", warn=True)
    # else:
    #     # Individually go through the phot_variable_flag indicator for each star to see if variable or not
    #     variableFlagList = gaia_result.columns["phot_variable_flag"]
    #     constantCounter = 0
    #     for currFlag in variableFlagList:
    #         if currFlag == "VARIABLE":
    #             return True
    #         elif currFlag == "NOT_AVAILABLE":
    #             continue
    #         elif currFlag == "CONSTANT":
    #             constantCounter += 1
    #     if constantCounter == len(variableFlagList):
    #         return False
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
        url = f"https://vsx.aavso.org/index.php?view=api.list&ra={ra}&dec={dec}&radius={radius}&tomag={maglimit}&format=json"
        result = requests.get(url)
        return result.json()['VSXObjects']['VSXObject'][0]['AUID']
    except Exception:
        log.info("\nThe target star does not have an AUID.")
        return False


@retry(stop=stop_after_delay(30))
def vsx_variable(ra, dec, radius=0.01, maglimit=14):
    try:
        url = f"https://vsx.aavso.org/index.php?view=api.list&ra={ra}&dec={dec}&radius={radius}&tomag={maglimit}&format=json"
        result = requests.get(url)
        var = result.json()['VSXObjects']['VSXObject'][0]['Category']

        if var.lower() == "variable":
            vname = result.json()['VSXObjects']['VSXObject'][0]['Name']
            vdec = result.json()['VSXObjects']['VSXObject'][0]['Declination2000']
            vra = result.json()['VSXObjects']['VSXObject'][0]['RA2000']
            log_info(f"\nVSX variable check found {vname} at RA {vra}, DEC {vdec}\n"
                     f"and will be removed from reduction.", warn=True)
            return True
        return False
    except Exception:
        return False

def build_comp_ra_dec(ra_wcs, dec_wcs, comp_stars):
    comp_ra_dec = []
    for _, comp_star in enumerate(comp_stars[:]):
        comp_ra_dec.append([ra_wcs[int(comp_star[1])][int(comp_star[0])],
            dec_wcs[int(comp_star[1])][int(comp_star[0])]])
    return comp_ra_dec

def check_for_variable_stars(ra_wcs, dec_wcs, comp_stars):
    for i, comp_star in enumerate(comp_stars[:]):
        ra = ra_wcs[int(comp_star[1])][int(comp_star[0])]
        dec = dec_wcs[int(comp_star[1])][int(comp_star[0])]

        log_info(f"\nChecking for variability in Comparison Star #{i + 1}:"
                 f"\n\tPixel X: {comp_star[0]} Pixel Y: {comp_star[1]}")
        if query_variable_star_apis(ra, dec):
            comp_stars.remove(comp_star)


@retry(stop=stop_after_delay(30))
def gaia_query(sample, radius):
    try:
        gaia_query = Gaia.cone_search(sample, radius)
        return gaia_query.get_results()
    except Exception:
        return False


@retry(stop=stop_after_delay(30))
def simbad_query(sample):
    try:
        simbad_result = Simbad.query_region(sample, radius=20 * u.arcsec)
        return simbad_result['MAIN_ID'][0].decode("utf-8")
    except Exception:
        return False


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

    url = f"https://apps.aavso.org/vsp/api/chart/?format=json&ra={ra:5f}&dec={dec:5f}&fov={fov}&maglimit={maglimit}"
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


# Aligns imaging data from .fits file to easily track the host and comparison star's positions
def transformation(image_data, file_name, roi=1):
    # crop image to ROI
    height = image_data.shape[1]
    width = image_data.shape[2]
    roix = slice(int(width * (0.5 - roi / 2)), int(width * (0.5 + roi / 2)))
    roiy = slice(int(height * (0.5 - roi / 2)), int(height * (0.5 + roi / 2)))

    # Find transformation from .FITS files and catch exceptions if not able to.
    try:
        results = aa.find_transform(image_data[1][roiy, roix], image_data[0][roiy, roix])
        return results[0]
    except Exception:
        ws = 5
        # smooth image and try to align again
        windows = view_as_windows(image_data[0], (ws,ws), step=1)
        medimg = np.median(windows, axis=(2,3))

        windows = view_as_windows(image_data[1], (ws,ws), step=1)
        medimg1 = np.median(windows, axis=(2,3))

        try:
            results = aa.find_transform(medimg1[roiy, roix], medimg[roiy, roix])
            return results[0]
        except Exception:
            pass

        for p in [99, 98, 95, 90]:
            for it in [2, 1, 0]:

                # create binary mask to align image
                mask1 = image_data[1][roiy, roix] > np.percentile(image_data[1][roiy, roix], p)
                mask1 = binary_erosion(mask1, iterations=it)

                mask0 = image_data[0][roiy, roix] > np.percentile(image_data[0][roiy, roix], p)
                mask0 = binary_erosion(mask0, iterations=it)

                try:
                    results = aa.find_transform(mask1, mask0)
                    return results[0]
                except Exception:
                    try:
                        result1 = ird.similarity(image_data[1][roiy, roix], image_data[0][roiy, roix], numiter=3)
                        return SimilarityTransform(scale=result1['scale'], rotation=np.radians(result1['angle']),
                                                   translation=[-1 * result1['tvec'][1], -1 * result1['tvec'][0]])
                    except Exception:
                        pass

    log.debug(f"Warning: Following image failed to align - {file_name}")
    plateStatus.alignmentError()
    return SimilarityTransform(scale=1, rotation=0, translation=[0, 0])


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

    missing_values = [parameter_names[key] for key in ['dist', 'pm_ra', 'pm_dec'] if info_dict.get(key, 0.0) == 0.0]

    if missing_values:
        missing_values = ", ".join(missing_values)
        log_info("Warning: Cannot account for proper motion due to missing values in: "
                 f"\n{missing_values}. Please re-run and fill in values in the initialization file to account "
                 f"for proper motion", warn=True)

        return info_dict['ra'], info_dict['dec']
    else:
        time_j2000 = Time(2000.0, format='jyear')
        time_obs = Time(time_obs, format='jd')

        coord = SkyCoord(
            ra=info_dict['ra'] * u.deg,
            dec=info_dict['dec'] * u.deg,
            distance=info_dict['dist'] * u.pc,
            pm_ra_cosdec=info_dict['pm_ra'] * u.mas / u.yr,
            pm_dec=info_dict['pm_dec'] * u.mas / u.yr,
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


# Method fits a 2D gaussian function that matches the star_psf to the star image and returns its pixel coordinates
def fit_centroid(data, pos, starIndex, psf_function=gaussian_psf, box=15, weightedcenter=True):
    # get sub field in image
    xv, yv = mesh_box(pos, box, maxx=data.shape[1], maxy=data.shape[0])
    subarray = data[yv, xv]
    try:
        init = [np.nanmax(subarray) - np.nanmin(subarray), 1, 1, 0, np.nanmin(subarray)]
    except ValueError as ve:
        # Handle null subfield - cannot solve
        plateStatus.outOfFrameWarning(starIndex)
        log.debug(f"Warning: empty subfield for fit_centroid at {np.round(pos, 2)}")
        return np.empty(7) * np.nan
    # compute flux weighted centroid in x and y
    wx = np.sum(xv[0]*subarray.sum(0))/subarray.sum(0).sum()
    wy = np.sum(yv[:,0]*subarray.sum(1))/subarray.sum(1).sum()

    # lower bound: [xc, yc, amp, sigx, sigy, rotation,  bg]
    lo = [pos[0] - box * 0.5, pos[1] - box * 0.5, 0, 0.5, 0.5, -np.pi / 4, np.nanmin(subarray) - 1]
    up = [pos[0] + box * 0.5, pos[1] + box * 0.5, 1e7, 20, 20, np.pi / 4, np.nanmax(subarray) + 1]

    def fcn2min(pars):
        model = psf_function(xv, yv, *pars)
        return (subarray - model).flatten()

    try:
        res = least_squares(fcn2min, x0=[*pos, *init], bounds=[lo, up], jac='3-point', xtol=None, method='trf')
    except:
        # Report low flux warning 
        plateStatus.lowFluxAmplitudeWarning(starIndex, pos[0], pos[1])
        log.debug(f"Warning: Measured flux amplitude is really low---are you sure there is a star at {np.round(pos, 2)}?")

        res = least_squares(fcn2min, x0=[*pos, *init], jac='3-point', xtol=None, method='lm')

    # override psf fit results with weighted centroid
    if weightedcenter:
        res.x[0] = wx
        res.x[1] = wy

    return res.x


# Method calculates the flux of the star (uses the skybg_phot method to do background sub)
def aperPhot(data, starIndex, xc, yc, r=5, dr=5):
    # Check for invalid coordinates
    if np.isnan(xc) or np.isnan(yc):
        return 0, 0
    
    # Calculate background if dr > 0
    if dr > 0:
        bgflux, sigmabg, Nbg = skybg_phot(data, starIndex, xc, yc, r + 2, dr)
    else:
        bgflux, sigmabg, Nbg = 0, 0, 0
    
    # Create aperture and mask
    aperture = CircularAperture(positions=[(xc, yc)], r=r)
    mask = aperture.to_mask(method='exact')[0]
    data_cutout = mask.cutout(data)
    
    # Check if aperture is valid
    if data_cutout is None:
        # Aperture is partially or fully outside the image
        return 0, bgflux    # Return zero flux but valid background
    
    # Calculate and return aperture sum
    aperture_sum = (mask.data * (data_cutout - bgflux)).sum()
    return aperture_sum, bgflux


def skybg_phot(data, starIndex, xc, yc, r=10, dr=5, ptol=99, debug=False):
    # create a crude annulus to mask out bright background pixels
    # the box will not extend beyond the borders of the image
    image_height, image_width = data.shape
    xv, yv = mesh_box([xc, yc], np.round(r + dr), maxx=image_width, maxy=image_height)
    rv = ((xv - xc) ** 2 + (yv - yc) ** 2) ** 0.5
    mask = (rv > r) & (rv < (r + dr))
    try:
        cutoff = np.nanpercentile(data[yv, xv][mask], ptol)
    except IndexError:
        plateStatus.skyBackgroundWarning(starIndex, xc, yc)
        log.debug(f"Warning: IndexError, problem computing sky bg for {xc:.1f}, {yc:.1f}."
                 f"\nCheck if star is present or close to border.")

        # create pixel wise mask on entire image
        x = np.arange(data.shape[1])
        y = np.arange(data.shape[0])
        xv, yv = np.meshgrid(x, y)
        rv = ((xv - xc) ** 2 + (yv - yc) ** 2) ** 0.5
        mask = (rv > r) & (rv < (r + dr))
        cutoff = np.nanpercentile(data[yv, xv][mask], ptol)

    dat = np.array(data[yv, xv], dtype=float)
    dat[dat > cutoff] = np.nan  # ignore pixels brighter than percentile

    if debug:
        minb = data[yv, xv][mask].min()
        maxb = data[yv, xv][mask].mean() + 3 * data[yv, xv][mask].std()
        nanmask = np.nan * np.zeros(mask.shape)
        nanmask[mask] = 1
        bgsky = data[yv, xv] * nanmask
        cmode = mode(dat.flatten(), nan_policy='omit', keepdims=True).mode[0]
        amode = mode(bgsky.flatten(), nan_policy='omit', keepdims=True).mode[0]

        fig, ax = plt.subplots(2, 2, figsize=(9, 9))
        im = ax[0, 0].imshow(data[yv, xv], vmin=minb, vmax=maxb, cmap='inferno')
        ax[0, 0].set_title("Original Data")
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax[0, 0])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        ax[1, 0].hist(bgsky.flatten(), label=f'Sky Annulus ({np.nanmedian(bgsky):.1f}, {amode:.1f})',
                      alpha=0.5, bins=np.arange(minb, maxb))
        ax[1, 0].hist(dat.flatten(), label=f'Clipped ({np.nanmedian(dat):.1f}, {cmode:.1f})', alpha=0.5,
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
    return mode(dat.flatten(), nan_policy='omit', keepdims=True).mode[0], np.nanstd(dat.flatten()), np.sum(mask)

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


def realTimeReduce(i, target_name, p_dict, info_dict, ax):
    timeList, airMassList, exptimes, norm_flux = [], [], [], []

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
    exotic_UIprevTPX = info_dict['tar_coords'][0]
    exotic_UIprevTPY = info_dict['tar_coords'][1]

    plateStatus.setCurrentFilename(inputfiles[0])
    wcs_file = check_wcs(inputfiles[0], info_dict['save'], info_dict['plate_opt'], rt=True)
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

    first_image = fits.getdata(inputfiles[0])
    targ_sig_xy = fit_centroid(first_image, [exotic_UIprevTPX, exotic_UIprevTPY], 0)[3:5]

    # aperture size in stdev (sigma) of PSF
    aper = 3 * max(targ_sig_xy)
    annulus = 10

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
    for i, fileName in enumerate(inputfiles):
        plateStatus.setCurrentFilename(fileName)
        hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False,
                         ignore_missing_end=True)

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

        sys.stdout.write(f"Finding transformation {i + 1} of {len(inputfiles)} : {fileName}\n")
        log.debug(f"Finding transformation {i + 1} of {len(inputfiles)} : {fileName}\n")
        sys.stdout.flush()

        try:
            wcs_hdr = search_wcs(fileName)
            if not wcs_hdr.is_celestial:
                raise Exception

            if i == 0:
                tx, ty = exotic_UIprevTPX, exotic_UIprevTPY
            else:
                pix_coords = wcs_hdr.world_to_pixel_values(tar_radec[0], tar_radec[1])
                tx, ty = pix_coords[0].take(0), pix_coords[1].take(0)

            psf_data['target'][i] = fit_centroid(imageData, [tx, ty], 0)

            if i != 0 and np.abs((psf_data['target'][i][2] - psf_data['target'][i - 1][2])
                                 / psf_data['target'][i - 1][2]) > 0.5:
                raise Exception

            pix_coords = wcs_hdr.world_to_pixel_values(comp_radec[0][0], comp_radec[0][1])
            cx, cy = pix_coords[0].take(0), pix_coords[1].take(0)
            psf_data['comp'][i] = fit_centroid(imageData, [cx, cy], 1)

            if i != 0:
                if not (tar_comp_dist['comp'][0] - 1 <= abs(int(psf_data['comp'][0][0]) - int(psf_data['target'][i][0])) <= tar_comp_dist['comp'][0] + 1 and
                        tar_comp_dist['comp'][1] - 1 <= abs(int(psf_data['comp'][0][1]) - int(psf_data['target'][i][1])) <= tar_comp_dist['comp'][1] + 1) or \
                        np.abs((psf_data['comp'][i][2] - psf_data['comp'][i - 1][2]) / psf_data['comp'][i - 1][2]) > 0.5:
                    raise Exception
            else:
                tar_comp_dist['comp'][0] = abs(int(psf_data['comp'][0][0]) - int(psf_data['target'][0][0]))
                tar_comp_dist['comp'][1] = abs(int(psf_data['comp'][0][1]) - int(psf_data['target'][0][1]))
        except Exception:
            if i == 0:
                tform = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])
            else:
                tform = transformation(np.array([imageData, firstImage]), fileName)

            tx, ty = tform([exotic_UIprevTPX, exotic_UIprevTPY])[0]
            psf_data['target'][i] = fit_centroid(imageData, [tx, ty], 0)

            cx, cy = tform(comp_star)[0]
            psf_data['comp'][i] = fit_centroid(imageData, [cx, cy], 1)

            if i == 0:
                tar_comp_dist['comp'][0] = abs(int(psf_data['comp'][0][0]) - int(psf_data['target'][0][0]))
                tar_comp_dist['comp'][1] = abs(int(psf_data['comp'][0][1]) - int(psf_data['target'][0][1]))

        # aperture photometry
        if i == 0:
            sigma = float((psf_data['target'][0][3] + psf_data['target'][0][4]) * 0.5)
            aper *= sigma
            annulus *= sigma

        tFlux = aperPhot(imageData, 0, psf_data['target'][i, 0], psf_data['target'][i, 1], aper, annulus)[0]
        cFlux = aperPhot(imageData, 1, psf_data['comp'][i, 0], psf_data['comp'][i, 1], aper, annulus)[0]
        norm_flux.append(tFlux / cFlux)

        # close file + delete from memory
        hdul.close()
        del hdul
        # Replaced each loop, so clean up
        del imageData

    ax.clear()
    ax.set_title(target_name)
    ax.set_ylabel('Normalized Flux')
    ax.set_xlabel('Time (JD)')
    ax.plot(timeList, norm_flux, 'bo')


def fit_lightcurve(times, tFlux, cFlux, airmass, ld, pDict, jd_times=None):
    # remove outliers
    si = np.argsort(times)
    dt = np.mean(np.diff(np.sort(times)))
    ndt = int(25. / 24. / 60. / dt) * 2 + 1
    if ndt > len(times):
        ndt = int(len(times)/4) * 2 + 1
    filtered_data = sigma_clip((tFlux / cFlux)[si], sigma=3, dt=max(5,ndt))
    arrayFinalFlux = (tFlux / cFlux)[si][~filtered_data]
    f1 = tFlux[si][~filtered_data]
    sigf1 = f1 ** 0.5
    f2 = cFlux[si][~filtered_data]
    sigf2 = f2 ** 0.5
    if np.sum(cFlux) == len(cFlux):
        arrayNormUnc = sigf1
    else:
        arrayNormUnc = np.sqrt((sigf1 / f2) ** 2 + (sigf2 * f1 / f2 ** 2) ** 2)
    arrayTimes = times[si][~filtered_data]
    arrayJDTimes = jd_times[si][~filtered_data]
    arrayAirmass = airmass[si][~filtered_data]

    # remove nans
    nanmask = np.isnan(arrayFinalFlux) | np.isnan(arrayNormUnc) | np.isnan(arrayTimes) | np.isnan(
        arrayAirmass) | np.less_equal(arrayFinalFlux, 0) | np.less_equal(arrayNormUnc, 0)
    nanmask = nanmask | np.isinf(arrayFinalFlux) | np.isinf(arrayNormUnc) | np.isinf(arrayTimes) | np.isinf(
        arrayAirmass)

    if np.sum(~nanmask) <= 1:
        log_info('No data left after filtering', warn=True)
        return None, None, None
    else:
        arrayFinalFlux = arrayFinalFlux[~nanmask]
        arrayNormUnc = arrayNormUnc[~nanmask]
        arrayTimes = arrayTimes[~nanmask]
        arrayJDTimes = arrayJDTimes[~nanmask]
        arrayAirmass = arrayAirmass[~nanmask]


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
        'a1': arrayFinalFlux.mean(),  # max() - arrayFinalFlux.min(), #mid Flux
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

    if np.floor(arrayPhases).max() - np.floor(arrayPhases).min() == 0:
        log_info("\nWarning:", warn=True)
        log_info(" Estimated mid-transit time is not within the observations", warn=True)
        log_info(" Check Period & Mid-transit time in inits.json. Make sure the uncertainties are not 0 or Nan.", warn=True)
        log_info(f"  obs start:{arrayTimes.min()}", warn=True)
        log_info(f"    obs end:{arrayTimes.max()}", warn=True)
        log_info(f" tmid prior:{prior['tmid']}\n", warn=True)

    mybounds = {
        'rprs': [0, prior['rprs'] * 1.25],
        'tmid': [lower, upper],
        'inc': [prior['inc'] - 5, min(90, prior['inc'] + 5)],
        'a1': [0.5 * min(arrayFinalFlux), 2 * max(arrayFinalFlux)],
        'a2': [-1, 1]
    }

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

    return myfit, f1, f2


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
    return parser.parse_args()


def main():
    # command line args
    args = parse_args()

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

        anim = FuncAnimation(fig, realTimeReduce, fargs=(userpDict['pName'], userpDict, exotic_infoDict, ax), interval=15000)
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

            wcs_file = check_wcs(inputfiles[0], exotic_infoDict['save'], exotic_infoDict['plate_opt'])
            img_scale_str, img_scale = get_img_scale(header, wcs_file, exotic_infoDict['pixel_scale'])
            plateStatus.initializeComparisonStarCount(len(exotic_infoDict['comp_stars']))
            ra_dec_tar, ra_dec_wcs = None, []
            chart_id, vsp_comp_stars, vsp_list = None, None, []

            if wcs_file:
                log_info(f"\nHere is the path to your plate solution: {wcs_file}")
                wcs_header = fits.getheader(filename=wcs_file)
                ra_wcs, dec_wcs = get_ra_dec(wcs_header)

                exotic_UIprevTPX, exotic_UIprevTPY = check_target_pixel_wcs(exotic_UIprevTPX, exotic_UIprevTPY,
                                                                            pDict, ra_wcs, dec_wcs,
                                                                            fits.getdata(inputfiles[0]),
                                                                            jd_times[0])
                ra_dec_tar = (ra_wcs[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)],
                             dec_wcs[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)])

                auid = vsx_auid(ra_dec_tar[0], ra_dec_tar[1])

                check_for_variable_stars(ra_wcs, dec_wcs, exotic_infoDict['comp_stars'])

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
                    check_for_variable_stars(ra_wcs, dec_wcs, exotic_infoDict['comp_stars'])
                # Build RA/Dec for comp after list is finalized (avoid off by one issues, etc
                ra_dec_wcs = build_comp_ra_dec(ra_wcs, dec_wcs, exotic_infoDict['comp_stars'])
                plateStatus.initializeComparisonStarCount(len(exotic_infoDict['comp_stars']))

            # aperture sizes in stdev (sigma) of PSF
            apers = np.linspace(1.5, 6, 20)
            annuli = np.linspace(6, 15, 19)

            # alloc psf fitting param
            psf_data = {
                # x-cent, y-cent, amplitude, sigma-x, sigma-y, rotation, offset
                'target': np.zeros((len(inputfiles), 7)),  # PSF fit
            }
            aper_data = {
                'target': np.zeros((len(inputfiles), len(apers), len(annuli))),
                'target_bg': np.zeros((len(inputfiles), len(apers), len(annuli)))
            }
            tar_comp_dist = {}
            vsp_num = []

            for i, coord in enumerate(exotic_infoDict['comp_stars']):
                ckey = f"comp{i + 1}"
                if coord in vsp_list:
                    vsp_num.append(i)
                psf_data[ckey] = np.zeros((len(inputfiles), 7))
                aper_data[ckey] = np.zeros((len(inputfiles), len(apers), len(annuli)))
                aper_data[f"{ckey}_bg"] = np.zeros((len(inputfiles), len(apers), len(annuli)))
                tar_comp_dist[ckey] = np.zeros(2)

            # open files, calibrate, align, photometry
            for i, fileName in enumerate(inputfiles):
                plateStatus.setCurrentFilename(fileName)
                hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False,
                                 ignore_missing_end=True)

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

                sys.stdout.write(f"Finding transformation {i + 1} of {len(inputfiles)} : {fileName}\n")
                log.debug(f"Finding transformation {i + 1} of {len(inputfiles)} : {fileName}\n")
                sys.stdout.flush()

                try:
                    wcs_hdr = search_wcs(fileName)
                    if not wcs_hdr.is_celestial:
                        raise Exception

                    if i == 0:
                        tx, ty = exotic_UIprevTPX, exotic_UIprevTPY
                    else:
                        pix_coords = wcs_hdr.world_to_pixel_values(ra_dec_tar[0], ra_dec_tar[1])
                        tx, ty = pix_coords[0].take(0), pix_coords[1].take(0)

                    psf_data['target'][i] = fit_centroid(imageData, [tx, ty], 0)

                    # TODO: Add check for flux on target/comp stars relative to others in the field
                    # in case of cloudy data, large changes, etc.
                    if i != 0 and np.abs((psf_data['target'][i][2] - psf_data['target'][i - 1][2])
                                         / psf_data['target'][i - 1][2]) > 0.5:
                        raise Exception

                    for j in range(len(exotic_infoDict['comp_stars'])):
                        ckey = f"comp{j + 1}"

                        pix_coords = wcs_hdr.world_to_pixel_values(ra_dec_wcs[j][0], ra_dec_wcs[j][1])
                        cx, cy = pix_coords[0].take(0), pix_coords[1].take(0)
                        psf_data[ckey][i] = fit_centroid(imageData, [cx, cy], j+1)

                        if i != 0:
                            if not (tar_comp_dist[ckey][0] - 1 <= abs(int(psf_data[ckey][i][0]) - int(psf_data['target'][i][0])) <= tar_comp_dist[ckey][0] + 1 and
                                    tar_comp_dist[ckey][1] - 1 <= abs(int(psf_data[ckey][i][1]) - int(psf_data['target'][i][1])) <= tar_comp_dist[ckey][1] + 1) or \
                                    np.abs((psf_data[ckey][i][2] - psf_data[ckey][i - 1][2]) / psf_data[ckey][i - 1][2]) > 0.5:
                                raise Exception
                        else:
                            tar_comp_dist[ckey][0] = abs(int(psf_data[ckey][0][0]) - int(psf_data['target'][0][0]))
                            tar_comp_dist[ckey][1] = abs(int(psf_data[ckey][0][1]) - int(psf_data['target'][0][1]))
                except Exception:
                    if i == 0:
                        tform = SimilarityTransform(scale=1, rotation=0, translation=[0, 0])
                    else:
                        tform = transformation(np.array([imageData, firstImage]), fileName)

                    tx, ty = tform([exotic_UIprevTPX, exotic_UIprevTPY])[0]
                    psf_data['target'][i] = fit_centroid(imageData, [tx, ty], 0)

                    for j, coord in enumerate(exotic_infoDict['comp_stars']):
                        ckey = f"comp{j + 1}"

                        cx, cy = tform(coord)[0]
                        psf_data[ckey][i] = fit_centroid(imageData, [cx, cy], j+1)

                        if i == 0:
                            tar_comp_dist[ckey][0] = abs(int(psf_data[ckey][0][0]) - int(psf_data['target'][0][0]))
                            tar_comp_dist[ckey][1] = abs(int(psf_data[ckey][0][1]) - int(psf_data['target'][0][1]))

                # aperture photometry
                if i == 0:
                    sigma = float((psf_data['target'][0][3] + psf_data['target'][0][4]) * 0.5)
                    apers *= sigma
                    annuli *= sigma

                for a, aper in enumerate(apers):
                    for an, annulus in enumerate(annuli):
                        if not np.isnan(psf_data['target'][i, 0]):
                            aper_data["target"][i][a][an], aper_data["target_bg"][i][a][an] = aperPhot(imageData, 0,
                                                                                                        psf_data['target'][i, 0],
                                                                                                        psf_data['target'][i, 1],
                                                                                                        aper, annulus)
                        else:
                            aper_data["target"][i][a][an] = np.nan
                            aper_data["target_bg"][i][a][an] = np.nan
                        # loop through comp stars
                        for j in range(len(exotic_infoDict['comp_stars'])):
                            ckey = f"comp{j + 1}"
                            if not np.isnan(psf_data[ckey][i][0]):
                                aper_data[ckey][i][a][an], \
                                aper_data[f"{ckey}_bg"][i][a][an] = aperPhot(imageData, j + 1, psf_data[ckey][i, 0],
                                                                            psf_data[ckey][i, 1], aper, annulus)
                            else:
                                aper_data[ckey][i][a][an] = np.nan
                                aper_data[f"{ckey}_bg"][i][a][an] = np.nan

                # close file + delete from memory
                hdul.close()
                del hdul
                del imageData

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
                'min_annulus': None
            }

            # loop over comp stars
            for j in range(len(exotic_infoDict['comp_stars'])):
                ckey = f"comp{j + 1}"

                cFlux = 2 * np.pi * psf_data[ckey][:, 2] * psf_data[ckey][:, 3] * psf_data[ckey][:, 4]
                myfit, tFlux1, cFlux1 = fit_lightcurve(times, tFlux, cFlux, airmass, ld, pDict, jd_times)

                if myfit is not None:
                    for k in myfit.bounds.keys():
                        log.debug(f"  {k}: {myfit.parameters[k]:.6f}")

                    log.debug("The Residual Standard Deviation is: "
                            f"{round(100 * myfit.residuals.std() / np.median(myfit.data), 6)}%")
                    log.debug(f"The Mean Squared Error is: {round(np.sum(myfit.residuals ** 2), 6)}\n")

                    res_std = myfit.residuals.std() / np.median(myfit.data)

                if photometry_info['min_std'] > res_std and myfit is not None:
                    photometry_info.update(best_fit_lc=copy.deepcopy(myfit),
                                           comp_star_num=j + 1, comp_star_coords=exotic_infoDict['comp_stars'][j],
                                           min_std=res_std, min_aperture=0, min_annulus=15 * sigma)

                    flux_values.update(flux_tar=tFlux1, flux_ref=cFlux1,
                                       flux_unc_tar=tFlux1 ** 0.5, flux_unc_ref=cFlux1 ** 0.5)

                    centroid_positions.update(x_targ=psf_data["target"][:, 0], y_targ=psf_data["target"][:, 1],
                                              x_ref=psf_data[ckey][:, 0], y_ref=psf_data[ckey][:, 1])

                if j in vsp_num:
                    ref_flux[j] = {
                        'myfit': copy.deepcopy(myfit),
                        'pos': exotic_infoDict['comp_stars'][j]
                    }

            log_info("\nComputing best comparison star, aperture, and sky annulus. Please wait.")

            # Aperture Photometry
            for a, aper in enumerate(apers):
                for an, annulus in enumerate(annuli):
                    tFlux = aper_data['target'][:, a, an]
                    ref_flux_opt, ref_flux_opt2, backtrack = False, False, True
                    temp_ref_flux = {i: None for i in vsp_num}

                    # fit without a comparison star
                    myfit, tFlux1, cFlux1 = fit_lightcurve(times, tFlux, np.ones(tFlux.shape[0]), airmass, ld, pDict, jd_times)

                    if myfit is not None:
                        for k in myfit.bounds.keys():
                            log.debug(f"  {k}: {myfit.parameters[k]:.6f}")

                        log.debug("The Residual Standard Deviation is: "
                                f"{round(100 * myfit.residuals.std() / np.median(myfit.data), 6)}%")
                        log.debug(f"The Mean Squared Error is: {round(np.sum(myfit.residuals ** 2), 6)}\n")

                        res_std = myfit.residuals.std() / np.median(myfit.data)
                    if photometry_info['min_std'] > res_std and myfit is not None:
                        ref_flux_opt = True

                        photometry_info.update(best_fit_lc=copy.deepcopy(myfit),
                                               comp_star_num=None, comp_star_coords=None,
                                               min_std=res_std, min_aperture=-aper, min_annulus=annulus)

                        flux_values.update(flux_tar=tFlux1, flux_ref=cFlux1,
                                           flux_unc_tar=tFlux1 ** 0.5, flux_unc_ref=cFlux1 ** 0.5)

                        centroid_positions.update(x_targ=psf_data["target"][:, 0], y_targ=psf_data["target"][:, 1],
                                                  x_ref=psf_data[ckey][:, 0], y_ref=psf_data[ckey][:, 1])

                    # try to fit data with comp star
                    for j in range(len(exotic_infoDict['comp_stars'])):
                        ckey = f"comp{j + 1}"
                        aper_mask = np.isfinite(aper_data[ckey][:, a, an])
                        cFlux = aper_data[ckey][aper_mask][:, a, an]

                        myfit, tFlux1, cFlux1 = fit_lightcurve(times[aper_mask], tFlux[aper_mask], cFlux, airmass[aper_mask], ld, pDict, jd_times[aper_mask])

                        if myfit is not None:
                            if j in vsp_num:
                                temp_ref_flux[j] = {
                                    'myfit': copy.deepcopy(myfit),
                                    'pos': exotic_infoDict['comp_stars'][j]
                                }

                            for k in myfit.bounds.keys():
                                log.debug(f"  {k}: {myfit.parameters[k]:.6f}")

                            log.debug("The Residual Standard Deviation is: "
                                    f"{round(100 * myfit.residuals.std() / np.median(myfit.data), 6)}%")
                            log.debug(f"The Mean Squared Error is: {round(np.sum(myfit.residuals ** 2), 6)}\n")

                            res_std = myfit.residuals.std() / np.median(myfit.data)
                        if photometry_info['min_std'] > res_std and myfit is not None:  # If the standard deviation is less than the previous min
                            ref_flux_opt2 = True

                            photometry_info.update(best_fit_lc=copy.deepcopy(myfit),
                                                   comp_star_num=j + 1,
                                                   comp_star_coords=exotic_infoDict['comp_stars'][j],
                                                   min_std=res_std, min_aperture=aper, min_annulus=annulus)

                            flux_values.update(flux_tar=tFlux1, flux_ref=cFlux1,
                                               flux_unc_tar=tFlux1 ** 0.5, flux_unc_ref=cFlux1 ** 0.5)

                            centroid_positions.update(x_targ=psf_data["target"][:, 0], y_targ=psf_data["target"][:, 1],
                                                      x_ref=psf_data[ckey][:, 0], y_ref=psf_data[ckey][:, 1])

                        if ref_flux_opt or ref_flux_opt2:
                            if j in vsp_num:
                                ref_flux[j] = {
                                    'myfit': copy.deepcopy(myfit),
                                    'pos': exotic_infoDict['comp_stars'][j]
                                }

                            if backtrack:
                                for i, value in enumerate(temp_ref_flux.values()):
                                    if value is not None and i != j:
                                        ref_flux[i] = value
                                backtrack = False

            log_info("\n\n*********************************************")
            if photometry_info['min_aperture'] == 0:  # psf
                log_info(f"Best Comparison Star: #{photometry_info['comp_star_num']}")
                log_info(f"Minimum Residual Scatter: {round(photometry_info['min_std'] * 100, 4)}%")
                log_info("Optimal Method: PSF photometry")
            elif photometry_info['min_aperture'] < 0:  # no comp star
                log_info("Best Comparison Star: None")
                log_info(f"Minimum Residual Scatter: {round(photometry_info['min_std'] * 100, 4)}%")
                log_info(f"Optimal Aperture: {abs(np.round(photometry_info['min_aperture'], 2))}")
                log_info(f"Optimal Annulus: {np.round(photometry_info['min_annulus'], 2)}")
            else:
                log_info(f"Best Comparison Star: #{photometry_info['comp_star_num']}")
                log_info(f"Minimum Residual Scatter: {round(photometry_info['min_std'] * 100, 4)}%")
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
            gi = ~sigma_clip(best_fit_lc.data[si], sigma=3, dt=ndt)  # good indexs

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

            best_fit_lc.time = best_fit_lc.time[si][gi]
            best_fit_lc.data = best_fit_lc.data[si][gi]
            best_fit_lc.airmass = best_fit_lc.airmass[si][gi]
            best_fit_lc.transit = best_fit_lc.transit[si][gi]
            best_fit_lc.jd_times = best_fit_lc.jd_times[si][gi]

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
            
            plot_fov(photometry_info['min_aperture'], photometry_info['min_annulus'], sigma,
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
                time_offset = 2400000.5 if exotic_infoDict['file_time'] == 'MJD_UTC' else 0.0
                goodTimes = convert_jd_to_bjd([time_ + time_offset for time_ in goodTimes], pDict, exotic_infoDict)

            if exotic_infoDict['file_units'] != 'flux':
                print("check flux convert")
                goodFluxes, goodNormUnc = flux_conversion(goodFluxes, goodNormUnc, exotic_infoDict['file_units'])

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

        mybounds = {
            'rprs': [0, prior['rprs'] * 1.25],
            'tmid': [lower, upper],
            'inc': [prior['inc'] - 5, min(90, prior['inc'] + 5)],
            'a2': [-3, 3],
        }

        if np.isnan(goodFluxes).all():
            log_info("Error: No valid photometry data found.", error=True)
            return

        # final light curve fit
        myfit = lc_fitter(goodTimes, goodFluxes, goodNormUnc, goodAirmasses, prior, mybounds, mode='ns')
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
                           exotic_infoDict['save'], exotic_infoDict['date'])

        #######################################################################
        # print final extracted planetary parameters
        #######################################################################

        log_info("\n*********************************************************")
        log_info("FINAL PLANETARY PARAMETERS\n")
        log_info(f"          Mid-Transit Time [BJD_TDB]: {round_to_2(myfit.parameters['tmid'], myfit.errors['tmid'])} +/- {round_to_2(myfit.errors['tmid'])}")
        log_info(f"  Radius Ratio (Planet/Star) [Rp/R*]: {round_to_2(myfit.parameters['rprs'], myfit.errors['rprs'])} +/- {round_to_2(myfit.errors['rprs'])}")
        log_info(f"           Transit depth [(Rp/R*)^2]: {round_to_2(100. * (myfit.parameters['rprs'] ** 2.))} +/- {round_to_2(100. * 2. * myfit.parameters['rprs'] * myfit.errors['rprs'])} [%]")
        log_info(f"           Orbital Inclination [inc]: {round_to_2(myfit.parameters['inc'], myfit.errors['inc'])} +/- {round_to_2(myfit.errors['inc'])}")
        log_info(f"               Airmass coefficient 1: {round_to_2(myfit.parameters['a1'], myfit.errors['a1'])} +/- {round_to_2(myfit.errors['a1'])}")
        log_info(f"               Airmass coefficient 2: {round_to_2(myfit.parameters['a2'], myfit.errors['a2'])} +/- {round_to_2(myfit.errors['a2'])}")
        log_info(f"                    Residual scatter: {round_to_2(100. * np.std(myfit.residuals / np.median(myfit.data)))} %")
        if fitsortext == 1:
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
