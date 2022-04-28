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
# Image alignment import
import astroalign as aa
aa.PIXEL_TOL = 1
# aa.NUM_NEAREST_NEIGHBORS=10
# astropy imports
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.io import fits
import astropy.time
from astropy.visualization import astropy_mpl_style
from astropy.wcs import WCS, FITSFixedWarning
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
# UTC to BJD converter import
from barycorrpy import utc_tdb
# julian conversion imports
import dateutil.parser as dup
from pathlib import Path
import pyvo as vo
import logging
from logging.handlers import TimedRotatingFileHandler
from matplotlib.animation import FuncAnimation
# Pyplot imports
import matplotlib.pyplot as plt
# from numba import njit
import numpy as np
# photometry
from photutils import CircularAperture
# scipy imports
from scipy.optimize import least_squares
from scipy.stats import mode
from scipy.signal import savgol_filter
from scipy.ndimage import binary_erosion
from skimage.util import view_as_windows
from skimage.transform import SimilarityTransform
# error handling for scraper
from tenacity import retry, stop_after_delay

# ########## EXOTIC imports ##########
try:  # light curve numerics
    from .api.elca import lc_fitter, binner, transit, get_phase
except ImportError:  # package import
    from api.elca import lc_fitter, binner, transit, get_phase
try:  # output files
    from inputs import Inputs
except ImportError:  # package import
    from .inputs import Inputs
try:  # plate solution
    from .api.plate_solution import PlateSolution
except ImportError:  # package import
    from api.plate_solution import PlateSolution
try:  # nonlinear limb darkening numerics
    from .api.gael_ld import createldgrid
except ImportError:  # package import
    from api.gael_ld import createldgrid
try:  # filters
    from .api.filters import fwhm
except ImportError:  # package import
    from api.filters import fwhm
try:  # nea
    from .api.nea import NASAExoplanetArchive
except ImportError:  # package import
    from api.nea import NASAExoplanetArchive
try:  # output files
    from output_files import OutputFiles
except ImportError:  # package import
    from .output_files import OutputFiles
try:  # plots
    from plots import plot_fov, plot_centroids, plot_obs_stats, plot_final_lightcurve, plot_flux
except ImportError:  # package import
    from .plots import plot_fov, plot_centroids, plot_obs_stats, plot_final_lightcurve, plot_flux
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


def sigma_clip(ogdata, sigma=3, dt=21, po=2):
    nanmask = np.isnan(ogdata)

    if po < dt <= len(ogdata):
        mdata = savgol_filter(ogdata[~nanmask], window_length=dt, polyorder=po)
        # mdata = median_filter(ogdata[~nanmask], dt)
        res = ogdata[~nanmask] - mdata
        std = np.nanmedian([np.nanstd(np.random.choice(res, 25)) for i in range(100)])
        # std = np.nanstd(res) # biased from large outliers
        sigmask = np.abs(res) > sigma * std
        nanmask[~nanmask] = sigmask

    return nanmask


# Get Julian time, don't need to divide by 2 since assume mid-EXPOSURE
# Find separate funciton in code that does julian conversion to BJD_TDB
# Method that gets and returns the julian time of the observation
def getJulianTime(header):
    exptime_offset = 0

    exp = header.get('EXPTIME')  # checking for variation in .fits header format
    if exp:
        exp = exp
    else:
        exp = header.get('EXPOSURE')

    # Grab the BJD first
    if 'BJD_TDB' in header:
        julianTime = float(header['BJD_TDB'])
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in header.comments['BJD_TDB']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    elif "BJD_TBD" in header:  # looks like TheSky misspells BJD_TDB as BJD_TBD -> hardcoding this common misspelling
        julianTime = float(header['BJD_TBD'])
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in header.comments['BJD_TBD']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    elif 'BJD' in header:
        julianTime = float(header['BJD'])
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in header.comments['BJD']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    # then the DATE-OBS
    elif "UT-OBS" in header:
        gDateTime = header['UT-OBS']  # gets the gregorian date and time from the fits file header
        dt = dup.parse(gDateTime)
        atime = astropy.time.Time(dt)
        julianTime = atime.jd
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in header.comments['UT-OBS']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    # Then Julian Date
    elif 'JULIAN' in header:
        julianTime = float(header['JULIAN'])
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in header.comments['JULIAN']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    # Then MJD-OBS last, as in the MicroObservatory headers, it has less precision
    elif ("MJD-OBS" in header) and ("epoch" not in header.comments['MJD-OBS']):
        julianTime = float(header["MJD-OBS"]) + 2400000.5
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in header.comments['MJD-OBS']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    else:
        if 'T' in header['DATE-OBS']:
            gDateTime = header['DATE-OBS']  # gets the gregorian date and time from the fits file header
        else:
            gDateTime = f"{header['DATE-OBS']}T{header['TIME-OBS']}"

        dt = dup.parse(gDateTime)
        atime = astropy.time.Time(dt)
        julianTime = atime.jd
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in header.comments['DATE-OBS']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now

    # If the mid-exposure time is given in the fits header, then no offset is needed to calculate the mid-exposure time
    return julianTime + exptime_offset


# Method that gets and returns the airmass from the fits file (Really the Altitude)
def getAirMass(image_header, ra, dec, lati, longit, elevation):
    # Grab airmass from image header; if not listed, calculate it from TELALT; if that isn't listed, then calculate it the hard way
    if 'AIRMASS' in image_header:
        am = float(image_header['AIRMASS'])
    elif 'TELALT' in image_header:
        alt = float(image_header['TELALT'])  # gets the airmass from the fits file header in (sec(z)) (Secant of the zenith angle)
        cosam = np.cos((np.pi / 180) * (90.0 - alt))
        am = 1 / cosam
    else:
        # pointing = SkyCoord(str(astropy.coordinates.Angle(raStr+" hours").deg)+" "+str(astropy.coordinates.Angle(decStr+" degrees").deg ), unit=(u.deg, u.deg), frame='icrs')
        pointing = SkyCoord(str(ra) + " " + str(dec), unit=(u.deg, u.deg), frame='icrs')

        location = EarthLocation.from_geodetic(lat=lati * u.deg, lon=longit * u.deg, height=elevation)
        atime = astropy.time.Time(getJulianTime(image_header), format='jd', scale='utc', location=location)
        pointingAltAz = pointing.transform_to(AltAz(obstime=atime, location=location))
        am = float(pointingAltAz.secz)
    return am


# Convert time units to BJD_TDB if pre-reduced file not in proper units
def timeConvert(timeList, timeFormat, pDict, info_dict):
    # Perform appropriate conversion for each time format if needed
    if timeFormat == 'JD_UTC':
        convertedTimes = utc_tdb.JDUTC_to_BJDTDB(timeList, ra=pDict['ra'], dec=pDict['dec'],
                                                 lat=info_dict['lat'], longi=info_dict['long'], alt=info_dict['elev'])
    elif timeFormat == 'MJD_UTC':
        convertedTimes = utc_tdb.JDUTC_to_BJDTDB(timeList + 2400000.5, ra=pDict['ra'], dec=pDict['dec'],
                                                 lat=info_dict['lat'], longi=info_dict['long'], alt=info_dict['elev'])
    timeList = convertedTimes[0]

    return timeList


# Convert magnitude units to flux if pre-reduced file not in flux already
def fluxConvert(fluxList, errorList, fluxFormat):
    # If units already in flux, do nothing, perform appropriate conversions to flux otherwise
    if fluxFormat == 'magnitude':
        convertedPositiveErrors = 10. ** ((-1. * (fluxList + errorList)) / 2.5)
        convertedNegativeErrors = 10. ** ((-1. * (fluxList - errorList)) / 2.5)
        fluxList = 10. ** ((-1. * fluxList) / 2.5)
    elif fluxFormat == 'millimagnitude':
        convertedPositiveErrors = 10. ** ((-1. * ((fluxList + errorList) / 1000.)) / 2.5)
        convertedNegativeErrors = 10. ** ((-1. * ((fluxList - errorList) / 1000.)) / 2.5)
        fluxList = 10. ** (-1. * (fluxList / 1000.) / 2.5)

    # Use distance from mean of upper/lower error bounds to calculate new sigma values
    positiveErrorDistance = abs(convertedPositiveErrors - fluxList)
    negativeErrorDistance = abs(convertedNegativeErrors - fluxList)
    meanErrorList = (positiveErrorDistance * negativeErrorDistance) ** 0.5

    return fluxList, meanErrorList


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
                     "Orbital Eccentricity (0 if null)",
                     "Star Effective Temperature (K)",
                     "Star Effective Temperature Positive Uncertainty (K)",
                     "Star Effective Temperature Negative Uncertainty (K)",
                     "Star Metallicity ([FE/H])",
                     "Star Metallicity Positive Uncertainty ([FE/H])",
                     "Star Metallicity Negative Uncertainty ([FE/H])",
                     "Star Surface Gravity (log(g))",
                     "Star Surface Gravity Positive Uncertainty (log(g))",
                     "Star Surface Gravity Negative Uncertainty (log(g))"]

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


class LimbDarkening:

    def __init__(self, teff=None, teffpos=None, teffneg=None, met=None, metpos=None, metneg=None,
                 logg=None, loggpos=None, loggneg=None, wl_min=None, wl_max=None, filter_type=None):
        self.priors = {'T*': teff, 'T*_uperr': teffpos, 'T*_lowerr': teffneg,
                       'FEH*': met, 'FEH*_uperr': metpos, 'FEH*_lowerr': metneg,
                       'LOGG*': logg, 'LOGG*_uperr': loggpos, 'LOGG*_lowerr': loggneg}
        self.filter_type = filter_type
        self.wl_min = wl_min
        self.wl_max = wl_max
        self.fwhm = fwhm
        self.ld0 = self.ld1 = self.ld2 = self.ld3 = None

    def nonlinear_ld(self):
        if self.filter_type and not (self.wl_min or self.wl_max):
            self._standard()
        elif self.wl_min or self.wl_max:
            self._custom()
        else:
            opt = user_input("\nWould you like EXOTIC to calculate your limb darkening parameters "
                             "with uncertainties? (y/n):", type_=str, values=['y', 'n'])

            if opt == 'y':
                opt = user_input("Please enter 1 to use a standard filter or 2 for a customized filter:",
                                 type_=int, values=[1, 2])
                if opt == 1:
                    self._standard()
                elif opt == 2:
                    self._custom()
            else:
                self._user_entered()
        return self.ld0, self.ld1, self.ld2, self.ld3, self.filter_type, self.wl_min * 1000, self.wl_max * 1000

    def _standard_list(self):
        log_info("\n\n***************************")
        log_info("Limb Darkening Coefficients")
        log_info("***************************")
        log_info("\nThe standard bands that are available for limb darkening parameters (https://www.aavso.org/filters)"
                 "\nas well as filters for MObs and LCO (0.4m telescope) datasets:\n")
        for key, value in self.fwhm.items():
            log_info(f"\t{key[1]}: {key[0]} - ({value[0]:.2f}-{value[1]:.2f}) nm")

    def _standard(self):
        while True:
            try:
                if not self.filter_type:
                    self._standard_list()
                    self.filter_type = user_input("\nPlease enter in the filter type (EX: Johnson V, V, STB, RJ):",
                                                  type_=str)
                for key, value in self.fwhm.items():
                    if self.filter_type in (key[0], key[1]) and self.filter_type != 'N/A':
                        self.filter_type = (key[0], key[1])
                        break
                else:
                    raise KeyError
                break
            except KeyError:
                log_info("\nError: The entered filter is not in the provided list of standard filters.", error=True)
                self.filter_type = None

        self.wl_min = self.fwhm[self.filter_type][0]
        self.wl_max = self.fwhm[self.filter_type][1]
        self.filter_type = self.filter_type[1]
        self._calculate_ld()

    def _custom(self):
        self.filter_type = 'N/A'
        if not self.wl_min:
            self.wl_min = user_input("FWHM Minimum wavelength (nm):", type_=float)
        if not self.wl_max:
            self.wl_max = user_input("FWHM Maximum wavelength (nm):", type_=float)
        self._calculate_ld()

    def _user_entered(self):
        self.filter_type = user_input("\nEnter in your filter name:", type_=str)
        ld_0 = user_input("\nEnter in your first nonlinear term:", type_=float)
        ld0_unc = user_input("Enter in your first nonlinear term uncertainty:", type_=float)
        ld_1 = user_input("\nEnter in your second nonlinear term:", type_=float)
        ld1_unc = user_input("Enter in your second nonlinear term uncertainty:", type_=float)
        ld_2 = user_input("\nEnter in your third nonlinear term:", type_=float)
        ld2_unc = user_input("Enter in your third nonlinear term uncertainty:", type_=float)
        ld_3 = user_input("\nEnter in your fourth nonlinear term:", type_=float)
        ld3_unc = user_input("Enter in your fourth nonlinear term uncertainty:", type_=float)
        self.ld0, self.ld1, self.ld2, self.ld3 = (ld_0, ld0_unc), (ld_1, ld1_unc), (ld_2, ld2_unc), (ld_3, ld3_unc)

        log.debug(f"Filter name: {self.filter_type}")
        log.debug(f"User-defined nonlinear limb-darkening coefficients: {ld_0}+/-{ld0_unc}, {ld_1}+/-{ld1_unc}, "
                  f"{ld_2}+/-{ld2_unc}, {ld_3}+/-{ld3_unc}")

    def _calculate_ld(self):
        self.wl_min = self.wl_min / 1000
        self.wl_max = self.wl_max / 1000
        ld_params = createldgrid(np.array([self.wl_min]), np.array([self.wl_max]), self.priors)
        self.ld0 = ld_params['LD'][0][0], ld_params['ERR'][0][0]
        self.ld1 = ld_params['LD'][1][0], ld_params['ERR'][1][0]
        self.ld2 = ld_params['LD'][2][0], ld_params['ERR'][2][0]
        self.ld3 = ld_params['LD'][3][0], ld_params['ERR'][3][0]

        log.debug("EXOTIC-calculated nonlinear limb-darkening coefficients: ")
        log.debug(f"{ld_params['LD'][0][0]} +/- + {ld_params['ERR'][0][0]}")
        log.debug(f"{ld_params['LD'][1][0]} +/- + {ld_params['ERR'][1][0]}")
        log.debug(f"{ld_params['LD'][2][0]} +/- + {ld_params['ERR'][2][0]}")
        log.debug(f"{ld_params['LD'][3][0]} +/- + {ld_params['ERR'][3][0]}")


def corruption_check(files):
    valid_files = []
    for file in files:
        try:
            hdul = fits.open(name=file, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
            valid_files.append(file)
        except OSError as e:
            log_info(f"Warning: corrupted file found and removed from reduction\n\t-File: {file}\n\t-Reason: {e}", warn=True)
        finally:
            if getattr(hdul, "close", None) and callable(hdul.close):
                hdul.close()
            del hdul

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


def get_wcs(file, directory=""):
    log_info("\nGetting the plate solution for your imaging file to translate pixel coordinates on the sky. "
             "\nPlease wait....")
    animate_toggle(True)
    wcs_obj = PlateSolution(file=file, directory=directory)
    wcs_file = wcs_obj.plate_solution()
    animate_toggle()
    return wcs_file


# Getting the right ascension and declination for every pixel in imaging file if there is a plate solution
def get_radec(header):
    wcs_header = WCS(header)
    xaxis = np.arange(header['NAXIS1'])
    yaxis = np.arange(header['NAXIS2'])
    x, y = np.meshgrid(xaxis, yaxis)
    return wcs_header.all_pix2world(x, y, 1)


def deg_to_pix(exp_ra, exp_dec, ra_list, dec_list):
    dist = (ra_list - exp_ra) ** 2 + (dec_list - exp_dec) ** 2
    return np.unravel_index(dist.argmin(), dist.shape)


# Check the ra and dec against the plate solution to see if the user entered in the correct values
def check_targetpixelwcs(pixx, pixy, expra, expdec, ralist, declist):
    while True:
        try:
            uncert = 1 / 36
            # Margins are within 100 arcseconds
            if not (expra - uncert <= ralist[int(pixy)][int(pixx)] <= expra + uncert):
                log_info("\nWarning: The X Pixel Coordinate entered does not match the target's Right Ascension.", warn=True)
                raise ValueError
            if not (expdec - uncert <= declist[int(pixy)][int(pixx)] <= expdec + uncert):
                log_info("\nWarning: The Y Pixel Coordinate entered does not match the target's Declination.", warn=True)
                raise ValueError
            return pixx, pixy
        except ValueError:
            log_info(f"Your input pixel coordinates: [{pixx}, {pixy}]")
            exotic_pixy, exotic_pixx = deg_to_pix(expra, expdec, ralist, declist)
            log_info(f"EXOTIC's calculated pixel coordinates: [{exotic_pixx}, {exotic_pixy}]")

            opt = user_input("Would you like to re-enter the pixel coordinates? (y/n): ", type_=str, values=['y', 'n'])

            # User wants to change their coordinates
            if opt == 'y':
                searchopt = user_input(f"Here are the suggested pixel coordinates:"
                                       f"  X Pixel: {exotic_pixx}"
                                       f"  Y Pixel: {exotic_pixy}"
                                       "\nWould you like to use these? (y/n): ",
                                       type_=str, values=['y', 'n'])
                # Use the coordinates found by code
                if searchopt == 'y':
                    return exotic_pixx, exotic_pixy
                # User enters their own coordinates to be re-checked
                else:
                    pixx = user_input("Please re-enter the target star's X Pixel Coordinate: ", type_=int)
                    pixy = user_input("Please re-enter the target star's Y Pixel Coordinate: ", type_=int)
            else:
                # User does not want to change coordinates even though they don't match the expected ra and dec
                return pixx, pixy


# Checks if comparison star is variable via querying SIMBAD
def variableStarCheck(ra, dec):
    # Convert comparison star coordinates from pixel to WCS
    sample = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    radius = u.Quantity(20.0, u.arcsec)

    # Query GAIA first to check for variability using the phot_variable_flag trait
    gaia_result = gaia_query(sample, radius)
    if not gaia_result:
        log_info("Warning: Your comparison star cannot be resolved in the Gaia star database; "
                 "EXOTIC cannot check if it is variable or not. "
                 "\nEXOTIC will still include this star in the reduction. "
                 "\nPlease proceed with caution as we cannot check for stellar variability.\n", warn=True)
    else:
        # Individually go through the phot_variable_flag indicator for each star to see if variable or not
        variableFlagList = gaia_result.columns["phot_variable_flag"]
        constantCounter = 0
        for currFlag in variableFlagList:
            if currFlag == "VARIABLE":
                return True
            elif currFlag == "NOT_AVAILABLE":
                continue
            elif currFlag == "CONSTANT":
                constantCounter += 1
        if constantCounter == len(variableFlagList):
            return False

    # Query SIMBAD and search identifier result table to determine if comparison star is variable in any form
    # This is a secondary check if GAIA query returns inconclusive results
    star_name = simbad_query(sample)
    if not star_name:
        log_info("Warning: Your comparison star cannot be resolved in the SIMBAD star database; "
                 "EXOTIC cannot check if it is variable or not. "
                 "\nEXOTIC will still include this star in the reduction. "
                 "\nPlease proceed with caution as we cannot check for stellar variability.\n", warn=True)
        return False
    else:
        identifiers = Simbad.query_objectids(star_name)

        for currName in identifiers:
            if "V*" in currName[0]:
                return True
        return False


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
    if gen_dark.size != 0:
        if i == 0:
            log_info("Dark subtracting images.")
        image_data = image_data - gen_dark
    elif gen_bias.size != 0:  # if a dark is not available, then at least subtract off the pedestal via the bias
        if i == 0:
            log_info("Bias-correcting images.")
        image_data = image_data - gen_bias
    else:
        pass

    if gen_flat.size != 0:
        if i == 0:
            log_info("Flattening images.")
        gen_flat[gen_flat == 0] = 1
        image_data = image_data / gen_flat

    return image_data


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
    except Exception as ee:
        log_info(ee)

        ws = 5
        # smooth image and try to align again
        windows = view_as_windows(image_data[0], (ws,ws), step=1)
        medimg = np.median(windows, axis=(2,3))

        windows = view_as_windows(image_data[1], (ws,ws), step=1)
        medimg1 = np.median(windows, axis=(2,3))

        try:
            results = aa.find_transform(medimg1[roiy, roix], medimg[roiy, roix])
            return results[0]
        except Exception as ee:
            log_info(ee)

        log_info(file_name)

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
                except Exception as ee:
                    log_info(ee)

    log_info(f"Warning: Alignment failed - {file_name}", warn=True)
    return SimilarityTransform(scale=1, rotation=0, translation=[0, 0])


def get_pixel_scale(wcs_header, header, pixel_init):
    astrometry_scale = None

    if wcs_header:
        astrometry_scale = [key.value.split(' ') for key in wcs_header._cards if 'scale:' in str(key.value)]

    if astrometry_scale:
        image_scale_num = astrometry_scale[0][1]
        image_scale_units = astrometry_scale[0][2]
        image_scale = f"Image scale in {image_scale_units}: {image_scale_num}"
    elif 'IM_SCALE' in header:
        image_scale_num = header['IM_SCALE']
        image_scale_units = header.comments['IM_SCALE']
        image_scale = f"Image scale in {image_scale_units}: {image_scale_num}"
    elif 'PIXSCALE' in header:
        image_scale_num = header['PIXSCALE']
        image_scale_units = header.comments['PIXSCALE']
        image_scale = f"Image scale in {image_scale_units}: {image_scale_num}"
    elif pixel_init:
        image_scale = f"Image scale in arcsecs/pixel: {pixel_init}"
    else:
        log_info("Not able to find Image Scale in the Image Header.")
        image_scale_num = user_input("Please enter Image Scale (arcsec/pixel): ", type_=float)
        image_scale = f"Image scale in arcsecs/pixel: {image_scale_num}"
    return image_scale


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


# finds target in WCS image after applying proper motion correction from SIMBAD
def find_target(target, hdufile, verbose=False):
    # query simbad to get proper motions
    service = vo.dal.TAPService("http://simbad.u-strasbg.fr/simbad/sim-tap")
    # http://simbad.u-strasbg.fr/simbad/tap/tapsearch.html
    query = '''
    SELECT basic.OID, ra, dec, main_id, pmra, pmdec
    FROM basic JOIN ident ON oidref = oid
    WHERE id = '{}';
    '''.format(target)

    result = service.search(query)
    # TODO check that simbad returned a value

    # set up astropy object
    coord = SkyCoord(
        ra=result['ra'][0] * u.deg,
        dec=result['dec'][0] * u.deg,
        distance=1 * u.pc,
        pm_ra_cosdec=result['pmra'][0] * u.mas / u.yr,
        pm_dec=result['pmdec'][0] * u.mas / u.yr,
        frame="icrs",
        obstime=astropy.time.Time("2000-1-1T00:00:00")
    )

    hdu = fits.open(hdufile)[0]

    try:
        dateobs = hdu.header['DATE_OBS']
    except:
        dateobs = hdu.header['DATE']

    # ignore timezone
    if len(dateobs.split('-')) == 4:
        dateobs = '-'.join(dateobs.split('-')[:-1])

    t = astropy.time.Time(dateobs, format='isot', scale='utc')
    coordpm = coord.apply_space_motion(new_obstime=t)

    # wcs coordinate translation
    wcs = WCS(hdu.header)

    pixcoord = wcs.wcs_world2pix([[coordpm.ra.value, coordpm.dec.value]], 0)

    if verbose:
        print("Simbad:", result)
        print("\nObs Date:", t)
        print("NEW:", coordpm.ra, coordpm.dec)
        print("\nTarget Location:", np.round(pixcoord[0], 2))

    return pixcoord[0]


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
def fit_centroid(data, pos, psf_function=gaussian_psf, box=15, weightedcenter=True):
    # get sub field in image
    xv, yv = mesh_box(pos, box, maxx=data.shape[1], maxy=data.shape[0])
    subarray = data[yv, xv]
    init = [np.nanmax(subarray) - np.nanmin(subarray), 1, 1, 0, np.nanmin(subarray)]
    
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
        log_info(f"Warning: Measured flux amplitude is really low---are you sure there is a star at {np.round(pos, 2)}?", warn=True)

        res = least_squares(fcn2min, x0=[*pos, *init], jac='3-point', xtol=None, method='lm')

    # override psf fit results with weighted centroid
    if weightedcenter:
        res.x[0] = wx
        res.x[1] = wy

    return res.x


# Method calculates the flux of the star (uses the skybg_phot method to do background sub)
def aperPhot(data, xc, yc, r=5, dr=5):
    if dr > 0:
        bgflux, sigmabg, Nbg = skybg_phot(data, xc, yc, r + 2, dr)
    else:
        bgflux, sigmabg, Nbg = 0, 0, 0

    aperture = CircularAperture(positions=[(xc, yc)], r=r)
    mask = aperture.to_mask(method='exact')[0]
    data_cutout = mask.cutout(data)
    aperture_sum = (mask.data * (data_cutout - bgflux)).sum()

    return aperture_sum, bgflux


def skybg_phot(data, xc, yc, r=10, dr=5, ptol=99, debug=False):
    # create a crude annulus to mask out bright background pixels
    xv, yv = mesh_box([xc, yc], np.round(r + dr))
    rv = ((xv - xc) ** 2 + (yv - yc) ** 2) ** 0.5
    mask = (rv > r) & (rv < (r + dr))
    try:
        cutoff = np.nanpercentile(data[yv, xv][mask], ptol)
    except IndexError:
        log_info(f"Warning: IndexError, problem computing sky bg for {xc:.1f}, {yc:.1f}."
                 f"\nCheck if star is present or close to border.", warn=True)

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
        cmode = mode(dat.flatten(), nan_policy='omit').mode[0]
        amode = mode(bgsky.flatten(), nan_policy='omit').mode[0]

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
    return mode(dat.flatten(), nan_policy='omit').mode[0], np.nanstd(dat.flatten()), np.sum(mask)


# Mid-Transit Time Prior Helper Functions
def numberOfTransitsAway(timeData, period, originalT):
    return int((np.nanmin(timeData) - originalT) / period) + 1


def nearestTransitTime(timeData, period, originalT):
    nearT = ((numberOfTransitsAway(timeData, period, originalT) * period) + originalT)
    return nearT


def save_comp_radec(wcs_file, ra_file, dec_file, comp_coords):
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


def realTimeReduce(i, target_name, info_dict, ax):
    timeList, airMassList, exptimes, norm_flux = [], [], [], []

    inputfiles = corruption_check(info_dict['images'])

    # time sort images
    times = []
    for ifile in inputfiles:
        extension = 0
        header = fits.getheader(filename=ifile, ext=extension)
        while header['NAXIS'] == 0:
            extension += 1
            header = fits.getheader(filename=ifile, ext=extension)
        times.append(getJulianTime(header))

    si = np.argsort(times)
    inputfiles = np.array(inputfiles)[si]
    exotic_UIprevTPX = info_dict['tar_coords'][0]
    exotic_UIprevTPY = info_dict['tar_coords'][1]

    wcs_file = check_wcs(inputfiles[0], info_dict['save'], info_dict['plate_opt'], rt=True)
    comp_star = info_dict['comp_stars']
    tar_radec, comp_radec = None, []

    if wcs_file:
        wcs_header = fits.getheader(filename=wcs_file)

        ra_file, dec_file = get_radec(wcs_header)
        tar_radec = (ra_file[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)],
                     dec_file[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)])

        ra = ra_file[int(comp_star[1])][int(comp_star[0])]
        dec = dec_file[int(comp_star[1])][int(comp_star[0])]

        comp_radec.append((ra, dec))

    first_image = fits.getdata(inputfiles[0])
    targ_sig_xy = fit_centroid(first_image, [exotic_UIprevTPX, exotic_UIprevTPY])[3:5]

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
        hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False,
                         ignore_missing_end=True)

        extension = 0
        image_header = hdul[extension].header
        while image_header["NAXIS"] == 0:
            extension += 1
            image_header = hdul[extension].header

        # TIME
        timeVal = getJulianTime(image_header)
        timeList.append(timeVal)

        # IMAGES
        imageData = hdul[extension].data

        if i == 0:
            firstImage = np.copy(imageData)

        sys.stdout.write(f"Finding transformation {i + 1} of {len(inputfiles)}\r")
        log.debug(f"Finding transformation {i + 1} of {len(inputfiles)}\r")
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

            psf_data['target'][i] = fit_centroid(imageData, [tx, ty])

            if i != 0 and np.abs((psf_data['target'][i][2] - psf_data['target'][i - 1][2])
                                 / psf_data['target'][i - 1][2]) > 0.5:
                raise Exception

            pix_coords = wcs_hdr.world_to_pixel_values(comp_radec[0][0], comp_radec[0][1])
            cx, cy = pix_coords[0].take(0), pix_coords[1].take(0)
            psf_data['comp'][i] = fit_centroid(imageData, [cx, cy])

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
            psf_data['target'][i] = fit_centroid(imageData, [tx, ty])

            cx, cy = tform(comp_star)[0]
            psf_data['comp'][i] = fit_centroid(imageData, [cx, cy])

            if i == 0:
                tar_comp_dist['comp'][0] = abs(int(psf_data['comp'][0][0]) - int(psf_data['target'][0][0]))
                tar_comp_dist['comp'][1] = abs(int(psf_data['comp'][0][1]) - int(psf_data['target'][0][1]))

        # aperture photometry
        if i == 0:
            sigma = float((psf_data['target'][0][3] + psf_data['target'][0][4]) * 0.5)
            aper *= sigma
            annulus *= sigma

        tFlux = aperPhot(imageData, psf_data['target'][i, 0], psf_data['target'][i, 1], aper, annulus)[0]
        cFlux = aperPhot(imageData, psf_data['comp'][i, 0], psf_data['comp'][i, 1], aper, annulus)[0]
        norm_flux.append(tFlux / cFlux)

        # close file + delete from memory
        hdul.close()
        del hdul
    del imageData

    ax.clear()
    ax.set_title(target_name)
    ax.set_ylabel('Normalized Flux')
    ax.set_xlabel('Time (JD)')
    ax.plot(timeList, norm_flux, 'bo')


def fit_lightcurve(times, tFlux, cFlux, airmass, ld, pDict):
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
    arrayAirmass = airmass[si][~filtered_data]

    # remove nans
    nanmask = np.isnan(arrayFinalFlux) | np.isnan(arrayNormUnc) | np.isnan(arrayTimes) | np.isnan(
        arrayAirmass) | np.less_equal(arrayFinalFlux, 0) | np.less_equal(arrayNormUnc, 0)
    nanmask = nanmask | np.isinf(arrayFinalFlux) | np.isinf(arrayNormUnc) | np.isinf(arrayTimes) | np.isinf(
        arrayAirmass)

    if np.sum(~nanmask) <= 1:
        log_info('No data left after filtering', warn=True)
    else:
        arrayFinalFlux = arrayFinalFlux[~nanmask]
        arrayNormUnc = arrayNormUnc[~nanmask]
        arrayTimes = arrayTimes[~nanmask]
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

    if np.floor(arrayPhases).max() - np.floor(arrayPhases).min() == 0:
        log_info("\nWarning:", warn=True)
        log_info(" Estimated mid-transit time is not within the observations", warn=True)
        log_info(" Check Period & Mid-transit time in inits.json. Make sure the uncertainties are not 0 or Nan.", warn=True)
        log_info(f"  obs start:{arrayTimes.min()}", warn=True)
        log_info(f"    obs end:{arrayTimes.max()}", warn=True)
        log_info(f" tmid prior:{prior['tmid']}\n", warn=True)

    mybounds = {
        'rprs': [0, pDict['rprs'] * 1.25],
        'tmid': [lower, upper],
        'ars': [max(pDict['aRs'] - 15 * pDict['aRsUnc'], 0.), pDict['aRs'] + 15 * pDict['aRsUnc']],
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
        mode='lm'
    )

    return myfit


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

    xTargCent, yTargCent, xRefCent, yRefCent, finXTargCent, finYTargCent, finXRefCent, finYRefCent = ([] for m in range(8))

    minSTD = 100000  # sets the initial minimum standard deviation absurdly high so it can be replaced immediately
    # minChi2 = 100000

    userpDict = {'ra': None, 'dec': None, 'pName': None, 'sName': None, 'pPer': None, 'pPerUnc': None,
                 'midT': None, 'midTUnc': None, 'rprs': None, 'rprsUnc': None, 'aRs': None, 'aRsUnc': None,
                 'inc': None, 'incUnc': None, 'ecc': None, 'teff': None,
                 'teffUncPos': None, 'teffUncNeg': None, 'met': None, 'metUncPos': None, 'metUncNeg': None,
                 'logg': None, 'loggUncPos': None, 'loggUncNeg': None}

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

        anim = FuncAnimation(fig, realTimeReduce, fargs=(userpDict['pName'], exotic_infoDict, ax), interval=15000)
        plt.show()

    # ----USER INPUTS----------------------------------------------------------
    else:
        log_info("\n**************************")
        log_info("Complete Reduction Routine")
        log_info("**************************")

        init_path, wcs_file, wcs_header, ra_file, dec_file = None, None, None, None, None
        generalDark, generalBias, generalFlat = np.empty(shape=(0, 0)), np.empty(shape=(0, 0)), np.empty(shape=(0, 0))

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

        if fitsortext == 1:
            # Only do the dark correction if user selects this option
            if exotic_infoDict['darks']:
                darksImgList = []
                for darkFile in exotic_infoDict['darks']:
                    darkData = fits.getdata(darkFile)
                    darksImgList.append(darkData)
                generalDark = np.median(darksImgList, axis=0)

            if exotic_infoDict['biases']:
                biasesImgList = []
                for biasFile in exotic_infoDict['biases']:
                    biasData = fits.getdata(biasFile)
                    biasesImgList.append(biasData)
                generalBias = np.median(biasesImgList, axis=0)

            if exotic_infoDict['flats']:
                flatsImgList = []
                for flatFile in exotic_infoDict['flats']:
                    flatData = fits.getdata(flatFile)
                    flatsImgList.append(flatData)
                notNormFlat = np.median(flatsImgList, axis=0)

                # if the bias exists, bias subtract the flatfield
                if exotic_infoDict['biases']:
                    notNormFlat = notNormFlat - generalBias

                # NORMALIZE
                medi = np.median(notNormFlat)
                generalFlat = notNormFlat / medi

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

        ld_obj = LimbDarkening(teff=pDict['teff'], teffpos=pDict['teffUncPos'], teffneg=pDict['teffUncNeg'],
                               met=pDict['met'], metpos=pDict['metUncPos'], metneg=pDict['metUncNeg'],
                               logg=pDict['logg'], loggpos=pDict['loggUncPos'], loggneg=pDict['loggUncNeg'],
                               wl_min=exotic_infoDict['wl_min'], wl_max=exotic_infoDict['wl_max'],
                               filter_type=exotic_infoDict['filter'])
        ld0, ld1, ld2, ld3, exotic_infoDict['filter'], exotic_infoDict['wl_min'], exotic_infoDict['wl_max'] = \
            ld_obj.nonlinear_ld()
        ld = [ld0[0], ld1[0], ld2[0], ld3[0]]

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

            timeList, airMassList, exptimes = [], [], []

            inputfiles = corruption_check(exotic_infoDict['images'])

            # time sort images
            times = []
            for ifile in inputfiles:
                extension = 0
                header = fits.getheader(filename=ifile, ext=extension)
                while header['NAXIS'] == 0:
                    extension += 1
                    header = fits.getheader(filename=ifile, ext=extension)
                times.append(getJulianTime(header))

            # checks for MOBS data
            mobs_header = fits.getheader(filename=inputfiles[0], ext=0)
            if 'CREATOR' in mobs_header:
                if 'MicroObservatory' in mobs_header['CREATOR'] and 'MOBS' not in exotic_infoDict['second_obs'].upper():
                    if exotic_infoDict['second_obs'].upper() != "N/A":
                        exotic_infoDict['second_obs'] += ",MOBS"
                    else:
                        exotic_infoDict['second_obs'] = "MOBS"

            si = np.argsort(times)
            inputfiles = np.array(inputfiles)[si]
            exotic_UIprevTPX = exotic_infoDict['tar_coords'][0]
            exotic_UIprevTPY = exotic_infoDict['tar_coords'][1]

            # fit target in the first image and use it to determine aperture and annulus range
            inc = 0
            for ifile in inputfiles:
                first_image = fits.getdata(ifile)
                try:
                    get_first = fit_centroid(first_image, [exotic_UIprevTPX, exotic_UIprevTPY])
                    break
                except Exception:
                    inc += 1
                finally:
                    del first_image

            inputfiles = inputfiles[inc:]
            wcs_file = check_wcs(inputfiles[0], exotic_infoDict['save'], exotic_infoDict['plate_opt'])
            compStarList = exotic_infoDict['comp_stars']
            tar_radec, comp_radec = None, []

            if wcs_file:
                log_info(f"\nHere is the path to your plate solution: {wcs_file}")
                wcs_header = fits.getheader(filename=wcs_file)
                ra_file, dec_file = get_radec(wcs_header)

                # Checking pixel coordinates against plate solution
                exotic_UIprevTPX, exotic_UIprevTPY = check_targetpixelwcs(exotic_UIprevTPX, exotic_UIprevTPY,
                                                                          pDict['ra'], pDict['dec'], ra_file, dec_file)
                tar_radec = (ra_file[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)],
                             dec_file[int(exotic_UIprevTPY)][int(exotic_UIprevTPX)])

                for compn, comp in enumerate(exotic_infoDict['comp_stars'][:]):
                    ra = ra_file[int(comp[1])][int(comp[0])]
                    dec = dec_file[int(comp[1])][int(comp[0])]
                    comp_radec.append((ra, dec))

                    log_info(f"\nChecking for variability in Comparison Star #{compn+1}:"
                             f"\n\tPixel X: {comp[0]} Pixel Y: {comp[1]}")
                    if variableStarCheck(ra_file[int(comp[1])][int(comp[0])], dec_file[int(comp[1])][int(comp[0])]):
                        log_info("\nCurrent comparison star is variable, proceeding to next star.")
                        exotic_infoDict['comp_stars'].remove(comp)
                compStarList = exotic_infoDict['comp_stars']

            # aperture sizes in stdev (sigma) of PSF
            apers = np.linspace(3, 6, 7)
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

            for i, coord in enumerate(compStarList):
                ckey = f"comp{i + 1}"
                psf_data[ckey] = np.zeros((len(inputfiles), 7))
                aper_data[ckey] = np.zeros((len(inputfiles), len(apers), len(annuli)))
                aper_data[f"{ckey}_bg"] = np.zeros((len(inputfiles), len(apers), len(annuli)))
                tar_comp_dist[ckey] = np.zeros(2)

            # open files, calibrate, align, photometry
            for i, fileName in enumerate(inputfiles):
                hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False,
                                 ignore_missing_end=True)

                extension = 0
                image_header = hdul[extension].header
                while image_header["NAXIS"] == 0:
                    extension += 1
                    image_header = hdul[extension].header

                # TIME
                timeVal = getJulianTime(image_header)
                timeList.append(timeVal)

                # AIRMASS
                airMass = getAirMass(image_header, pDict['ra'], pDict['dec'], exotic_infoDict['lat'], exotic_infoDict['long'],
                                     exotic_infoDict['elev'])  # gets the airmass at the time the image was taken
                airMassList.append(airMass)  # adds that airmass value to the list of airmasses

                # EXPOSURE_TIME
                exp = image_header.get('EXPTIME')  # checking for variation in .fits header format
                if exp:
                    exptimes.append(image_header['EXPTIME'])
                else:
                    exptimes.append(image_header['EXPOSURE'])

                # IMAGES
                imageData = hdul[extension].data

                # CALS
                imageData = apply_cals(imageData, generalDark, generalBias, generalFlat, i)

                if i == 0:
                    image_scale = get_pixel_scale(wcs_header, image_header, exotic_infoDict['pixel_scale'])

                    firstImage = np.copy(imageData)

                sys.stdout.write(f"Finding transformation {i + 1} of {len(inputfiles)}\r")
                log.debug(f"Finding transformation {i + 1} of {len(inputfiles)}\r")
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

                    psf_data['target'][i] = fit_centroid(imageData, [tx, ty])

                    # TODO: Add check for flux on target/comp stars relative to others in the field
                    # in case of cloudy data, large changes, etc.
                    if i != 0 and np.abs((psf_data['target'][i][2] - psf_data['target'][i - 1][2])
                                         / psf_data['target'][i - 1][2]) > 0.5:
                        raise Exception

                    for j in range(len(compStarList)):
                        ckey = f"comp{j + 1}"

                        pix_coords = wcs_hdr.world_to_pixel_values(comp_radec[j][0], comp_radec[j][1])
                        cx, cy = pix_coords[0].take(0), pix_coords[1].take(0)
                        psf_data[ckey][i] = fit_centroid(imageData, [cx, cy])

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
                    psf_data['target'][i] = fit_centroid(imageData, [tx, ty])

                    for j, coord in enumerate(compStarList):
                        ckey = f"comp{j + 1}"

                        cx, cy = tform(coord)[0]
                        psf_data[ckey][i] = fit_centroid(imageData, [cx, cy])

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
                        aper_data["target"][i][a][an], aper_data["target_bg"][i][a][an] = aperPhot(imageData,
                                                                                                   psf_data['target'][i, 0],
                                                                                                   psf_data['target'][i, 1],
                                                                                                   aper, annulus)

                        # loop through comp stars
                        for j in range(len(compStarList)):
                            ckey = f"comp{j + 1}"
                            aper_data[ckey][i][a][an], aper_data[f"{ckey}_bg"][i][a][an] = aperPhot(imageData,
                                                                                                    psf_data[ckey][i, 0],
                                                                                                    psf_data[ckey][i, 1],
                                                                                                    aper, annulus)

                # close file + delete from memory
                hdul.close()
                del hdul
            del imageData

            # filter bad images
            badmask = (psf_data["target"][:, 0] == 0) | (aper_data["target"][:, 0, 0] == 0) | np.isnan(
                aper_data["target"][:, 0, 0])
            if np.sum(~badmask) == 0:
                log_info("No images to fit...check reference image for alignment (first image of sequence)")

            # convert to numpy arrays
            times = np.array(timeList)[~badmask]
            airmass = np.array(airMassList)[~badmask]
            psf_data["target"] = psf_data["target"][~badmask]
            si = np.argsort(times)

            exotic_infoDict['exposure'] = exp_time_med(exptimes)

            # PSF flux
            tFlux = 2 * np.pi * psf_data['target'][:, 2] * psf_data['target'][:, 3] * psf_data['target'][:, 4]

            # loop over comp stars
            for j in range(len(compStarList)):
                ckey = f"comp{j + 1}"
                psf_data[ckey] = psf_data[ckey][~badmask]

                cFlux = 2 * np.pi * psf_data[ckey][:, 2] * psf_data[ckey][:, 3] * psf_data[ckey][:, 4]
                myfit = fit_lightcurve(times, tFlux, cFlux, airmass, ld, pDict)

                for k in myfit.bounds.keys():
                    log.debug(f"  {k}: {myfit.parameters[k]:.6f}")

                log.debug("The Residual Standard Deviation is: "
                          f"{round(100 * myfit.residuals.std() / np.median(myfit.data), 6)}%")
                log.debug(f"The Mean Squared Error is: {round(np.sum(myfit.residuals ** 2), 6)}\n")

                resstd = myfit.residuals.std() / np.median(myfit.data)
                if minSTD > resstd:  # If the standard deviation is less than the previous min
                    bestCompStar = j + 1
                    comp_coords = compStarList[j]
                    minSTD = resstd
                    minAperture = 0
                    minAnnulus = 15 * sigma
                    # arrayNormUnc = myfit.dataerr

                    # sets the lists we want to print to correspond to the optimal aperature
                    goodFluxes = np.copy(myfit.data)
                    # goodNormUnc = np.copy(myfit.dataerr)
                    nonBJDTimes = np.copy(myfit.time)
                    goodAirmasses = np.copy(myfit.airmass)
                    goodTargets = tFlux
                    goodReferences = cFlux
                    goodTUnc = tFlux ** 0.5
                    goodRUnc = cFlux ** 0.5
                    bestlmfit = myfit

                    finXTargCent = psf_data["target"][:, 0]
                    finYTargCent = psf_data["target"][:, 1]
                    finXRefCent = psf_data[ckey][:, 0]
                    finYRefCent = psf_data[ckey][:, 1]

            log_info("\nComputing best comparison star, aperture, and sky annulus. Please wait.")

            # Aperture Photometry
            for a, aper in enumerate(apers):
                for an, annulus in enumerate(annuli):
                    tFlux = aper_data['target'][:, a, an]

                    # fit without a comparison star
                    myfit = fit_lightcurve(times, tFlux, np.ones(tFlux.shape[0]), airmass, ld, pDict)

                    for k in myfit.bounds.keys():
                        log.debug(f"  {k}: {myfit.parameters[k]:.6f}")

                    log.debug("The Residual Standard Deviation is: "
                              f"{round(100 * myfit.residuals.std() / np.median(myfit.data), 6)}%")
                    log.debug(f"The Mean Squared Error is: {round(np.sum(myfit.residuals ** 2), 6)}\n")

                    resstd = myfit.residuals.std() / np.median(myfit.data)
                    if minSTD > resstd:  # If the standard deviation is less than the previous min
                        minSTD = resstd
                        minAperture = -aper
                        minAnnulus = annulus
                        # arrayNormUnc = tFlux ** 0.5

                        # sets the lists we want to print to correspond to the optimal aperature
                        goodFluxes = np.copy(myfit.data)
                        # goodNormUnc = np.copy(myfit.dataerr)
                        nonBJDTimes = np.copy(myfit.time)
                        # nonBJDPhases = np.copy(myfit.phase)
                        goodAirmasses = np.copy(myfit.airmass)
                        goodTargets = tFlux
                        goodReferences = cFlux
                        goodTUnc = tFlux ** 0.5
                        goodRUnc = cFlux ** 0.5
                        # goodResids = myfit.residuals
                        bestlmfit = myfit

                        finXTargCent = psf_data["target"][:, 0]
                        finYTargCent = psf_data["target"][:, 1]
                        finXRefCent = psf_data[ckey][:, 0]
                        finYRefCent = psf_data[ckey][:, 1]

                    # try to fit data with comp star
                    for j in range(len(compStarList)):
                        ckey = f"comp{j + 1}"
                        cFlux = aper_data[ckey][:, a, an]

                        myfit = fit_lightcurve(times, tFlux, cFlux, airmass, ld, pDict)

                        for k in myfit.bounds.keys():
                            log.debug(f"  {k}: {myfit.parameters[k]:.6f}")

                        log.debug("The Residual Standard Deviation is: "
                                  f"{round(100 * myfit.residuals.std() / np.median(myfit.data), 6)}%")
                        log.debug(f"The Mean Squared Error is: {round(np.sum(myfit.residuals ** 2), 6)}\n")

                        resstd = myfit.residuals.std() / np.median(myfit.data)
                        if minSTD > resstd:  # If the standard deviation is less than the previous min
                            bestCompStar = j + 1
                            comp_coords = compStarList[j]
                            minSTD = resstd
                            minAperture = aper
                            minAnnulus = annulus
                            # arrayNormUnc = arrayNormUnc

                            # sets the lists we want to print to correspond to the optimal aperature
                            goodFluxes = np.copy(myfit.data)
                            # goodNormUnc = np.copy(myfit.dataerr)
                            nonBJDTimes = np.copy(myfit.time)
                            # nonBJDPhases = np.copy(myfit.phase)
                            goodAirmasses = np.copy(myfit.airmass)
                            goodTargets = tFlux
                            goodReferences = cFlux
                            goodTUnc = tFlux ** 0.5
                            goodRUnc = cFlux ** 0.5
                            # goodResids = myfit.residuals
                            bestlmfit = myfit

                            finXTargCent = psf_data["target"][:, 0]
                            finYTargCent = psf_data["target"][:, 1]
                            finXRefCent = psf_data[ckey][:, 0]
                            finYRefCent = psf_data[ckey][:, 1]

            # log best fit
            log_info("\n\n*********************************************")
            if minAperture == 0:  # psf
                log_info(f"Best Comparison Star: #{bestCompStar}")
                log_info(f"Minimum Residual Scatter: {round(minSTD * 100, 4)}%")
                log_info("Optimal Method: PSF photometry")
                bestaperture = "PSF photometry"
            elif minAperture < 0:  # no comp star
                log_info("Best Comparison Star: None")
                log_info(f"Minimum Residual Scatter: {round(minSTD * 100, 4)}%")
                log_info(f"Optimal Aperture: {abs(np.round(minAperture, 2))}")
                log_info(f"Optimal Annulus: {np.round(minAnnulus, 2)}")
                bestCompStar, comp_coords = None, None
                bestaperture = str(abs(np.round(minAperture, 2)))
            else:
                log_info(f"Best Comparison Star: #{bestCompStar}")
                log_info(f"Minimum Residual Scatter: {round(minSTD * 100, 4)}%")
                log_info(f"Optimal Aperture: {np.round(minAperture, 2)}")
                log_info(f"Optimal Annulus: {np.round(minAnnulus, 2)}")
                bestaperture = str(abs(np.round(minAperture, 2)))
            log_info("*********************************************\n")

            # Take the BJD times from the image headers
            if "BJD_TDB" in image_header or "BJD" in image_header or "BJD_TBD" in image_header:
                goodTimes = nonBJDTimes
            # If not in there, then convert all the final times into BJD - using astropy alone
            else:
                log_info("No Barycentric Julian Dates (BJDs) in Image Headers for standardizing time format. "
                         "Converting all JDs to BJD_TDBs.")
                log_info("Please be patient- this step can take a few minutes.")

                animate_toggle(True)
                try:
                    resultos = utc_tdb.JDUTC_to_BJDTDB(nonBJDTimes, ra=pDict['ra'], dec=pDict['dec'],
                                                       lat=exotic_infoDict['lat'], longi=exotic_infoDict['long'],
                                                       alt=exotic_infoDict['elev'])
                    goodTimes = resultos[0]
                except:
                    targetloc = astropy.coordinates.SkyCoord(pDict['ra'], pDict['dec'], unit=(u.deg, u.deg),
                                                             frame='icrs')
                    obsloc = astropy.coordinates.EarthLocation(lat=exotic_infoDict['lat'], lon=exotic_infoDict['long'],
                                                               height=exotic_infoDict['elev'])
                    timesToConvert = astropy.time.Time(nonBJDTimes, format='jd', scale='utc', location=obsloc)
                    ltt_bary = timesToConvert.light_travel_time(targetloc)
                    time_barycentre = timesToConvert.tdb + ltt_bary
                    goodTimes = time_barycentre.value

                animate_toggle()

            # sigma clip
            si = np.argsort(goodTimes)
            dt = np.mean(np.diff(np.sort(goodTimes)))
            ndt = int(30. / 24. / 60. / dt) * 2 + 1  # ~30 minutes
            gi = ~sigma_clip(goodFluxes[si], sigma=3, dt=ndt)  # good indexs

            # Calculate the proper timeseries uncertainties from the residuals of the out-of-transit data
            OOT = (bestlmfit.transit == 1)  # find out-of-transit portion of the lightcurve

            if sum(OOT) <= 1:
                OOTscatter = np.std(bestlmfit.residuals)
                goodNormUnc = OOTscatter * bestlmfit.airmass_model
                goodNormUnc = goodNormUnc / np.nanmedian(goodFluxes)
                goodFluxes = goodFluxes / np.nanmedian(goodFluxes)
            else:
                OOTscatter = np.std((bestlmfit.data / bestlmfit.airmass_model)[OOT])  # calculate the scatter in the data
                goodNormUnc = OOTscatter * bestlmfit.airmass_model  # scale this scatter back up by the airmass model and then adopt these as the uncertainties
                goodNormUnc = goodNormUnc / np.nanmedian(goodFluxes[OOT])
                goodFluxes = goodFluxes / np.nanmedian(goodFluxes[OOT])

            if np.isnan(goodFluxes).all():
                log_info("Error: No valid photometry data found.", error=True)
                return

            goodTimes = goodTimes[si][gi]
            goodFluxes = goodFluxes[si][gi]
            goodNormUnc = goodNormUnc[si][gi]
            goodAirmasses = goodAirmasses[si][gi]

            plot_fov(minAperture, minAnnulus, sigma, finXTargCent[0], finYTargCent[0], finXRefCent[0], finYRefCent[0],
                     firstImage, image_scale, pDict['pName'], exotic_infoDict['save'], exotic_infoDict['date'])

            # Centroid position plots
            plot_centroids(finXTargCent[si][gi], finYTargCent[si][gi], finXRefCent[si][gi], finYRefCent[si][gi],
                           goodTimes, pDict['pName'], exotic_infoDict['save'], exotic_infoDict['date'])

            # TODO: convert the exoplanet archive mid transit time to bjd - need to take into account observatory location listed in Exoplanet Archive
            # tMidtoC = astropy.time.Time(timeMidTransit, format='jd', scale='utc')
            # forPhaseResult = utc_tdb.JDUTC_to_BJDTDB(tMidtoC, ra=raDeg, dec=decDeg, lat=lati, longi=longit, alt=2000)
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

            plot_flux(goodTimes, goodTargets[si][gi], goodTUnc[si][gi], goodReferences[si][gi], goodRUnc[si][gi],
                      goodFluxes, goodNormUnc, goodAirmasses, pDict['pName'], exotic_infoDict['save'],
                      exotic_infoDict['date'])

            log_info("\n\nOutput File Saved")
        else:
            goodTimes, goodFluxes, goodNormUnc, goodAirmasses = [], [], [], []
            bestCompStar, comp_coords = None, None

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
                goodTimes = timeConvert(goodTimes, exotic_infoDict['file_time'], pDict, exotic_infoDict)

            if exotic_infoDict['file_units'] != 'flux':
                print("check flux convert")
                import pdb; pdb.set_trace()
                goodFluxes, goodNormUnc = fluxConvert(goodFluxes, goodNormUnc, exotic_infoDict['file_units'])

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
        if upper > prior['tmid'] + 0.5*prior['per']:
            upper = prior['tmid'] + 0.5*prior['per']
        if lower < prior['tmid'] - 0.5*prior['per']:
            lower = prior['tmid'] - 0.5*prior['per']

        if np.floor(phase).max() - np.floor(phase).min() == 0:
            log_info("Error: Estimated mid-transit not in observation range (check priors or observation time)", error=True)
            log_info(f"start:{np.min(goodTimes)}", error=True)
            log_info(f"  end:{np.max(goodTimes)}", error=True)
            log_info(f"prior:{prior['tmid']}", error=True)

        mybounds = {
            'rprs': [0, pDict['rprs'] * 1.25],
            'tmid': [lower, upper],
            'ars': [max(pDict['aRs'] - 15 * pDict['aRsUnc'], 0.), pDict['aRs'] + 15 * pDict['aRsUnc']],
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
            plot_obs_stats(myfit, compStarList, psf_data, si, gi, pDict['pName'],
                           exotic_infoDict['save'], exotic_infoDict['date'])

        #######################################################################
        # print final extracted planetary parameters
        #######################################################################

        log_info("\n*********************************************************")
        log_info("FINAL PLANETARY PARAMETERS\n")
        log_info(f"          Mid-Transit Time [BJD_TDB]: {round_to_2(myfit.parameters['tmid'], myfit.errors['tmid'])} +/- {round_to_2(myfit.errors['tmid'])}")
        log_info(f"  Radius Ratio (Planet/Star) [Rp/Rs]: {round_to_2(myfit.parameters['rprs'], myfit.errors['rprs'])} +/- {round_to_2(myfit.errors['rprs'])}")
        log_info(f"           Transit depth [(Rp/Rs)^2]: {round_to_2(100. * (myfit.parameters['rprs'] ** 2.))} +/- {round_to_2(100. * 2. * myfit.parameters['rprs'] * myfit.errors['rprs'])} [%]")
        log_info(f" Semi Major Axis/ Star Radius [a/Rs]: {round_to_2(myfit.parameters['ars'], myfit.errors['ars'])} +/- {round_to_2(myfit.errors['ars'])}")
        log_info(f"               Airmass coefficient 1: {round_to_2(myfit.parameters['a1'], myfit.errors['a1'])} +/- {round_to_2(myfit.errors['a1'])}")
        log_info(f"               Airmass coefficient 2: {round_to_2(myfit.parameters['a2'], myfit.errors['a2'])} +/- {round_to_2(myfit.errors['a2'])}")
        log_info(f"                    Residual scatter: {round_to_2(100. * np.std(myfit.residuals / np.median(myfit.data)))} %")
        if fitsortext == 1:
            if minAperture >= 0:
                log_info(f"                Best Comparison Star: #{bestCompStar} - {comp_coords}")
            else:
                log_info("                 Best Comparison Star: None")
            if minAperture == 0:
                log_info("                       Optimal Method: PSF photometry")
            else:
                log_info(f"                    Optimal Aperture: {abs(np.round(minAperture, 2))}")
                log_info(f"                     Optimal Annulus: {np.round(minAnnulus, 2)}")
        log_info(f"              Transit Duration [day]: {round_to_2(np.mean(durs), np.std(durs))} +/- {round_to_2(np.std(durs))}")
        log_info("*********************************************************")

        ##########
        # SAVE DATA
        ##########

        fig = myfit.plot_triangle()
        fig.savefig(Path(exotic_infoDict['save']) / "temp" /
                    f"Triangle_{pDict['pName']}_{exotic_infoDict['date']}.png")

        output_files = OutputFiles(myfit, pDict, exotic_infoDict, durs)
        error_txt = "\n\tPlease report this issue on the Exoplanet Watch Slack Channel in #data-reductions."

        try:
            phase = get_phase(myfit.time, pDict['pPer'], myfit.parameters['tmid'])
            output_files.final_lightcurve(phase)
        except Exception as e:
            log_info(f"\nError: Could not create FinalLightCurve.csv. {error_txt}\n\t{e}", error=True)
        try:
            if fitsortext == 1:
                output_files.final_planetary_params(phot_opt=True, comp_star=bestCompStar, comp_coords=comp_coords,
                                                    min_aper=np.round(minAperture, 2),
                                                    min_annul=np.round(minAnnulus, 2))
            else:
                output_files.final_planetary_params(phot_opt=False)
        except Exception as e:
            log_info(f"\nError: Could not create FinalParams.json. {error_txt}\n\t{e}", error=True)
        try:
            if bestCompStar:
                exotic_infoDict['phot_comp_star'] = save_comp_radec(wcs_file, ra_file, dec_file, comp_coords)
            output_files.aavso(exotic_infoDict['phot_comp_star'], goodAirmasses, ld0, ld1, ld2, ld3)
        except Exception as e:
            log_info(f"\nError: Could not create AAVSO.txt. {error_txt}\n\t{e}", error=True)

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
