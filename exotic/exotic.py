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
done_flag_animate_exotic = None
th_animate_exotic = None

import threading
import time
import sys


def animate():
    import itertools
    for c in itertools.cycle(['|', '/', '-', '\\']):
        if done_flag_animate_exotic:
            sys.stdout.write('\rThinking ... DONE!')
            sys.stdout.write('\n')
            break
        sys.stdout.write('\rThinking ' + c + ' ... ')
        sys.stdout.flush()
        time.sleep(0.11)
    sys.stdout.flush()


def animate_toggle(enabled=False):
    """
    Console feature intended to be run in a single synchronous block.
    """
    global done_flag_animate_exotic
    global th_animate_exotic
    counter_wait = 0
    if enabled:
        done_flag_animate_exotic = False
        th_animate_exotic = threading.Thread(target=animate, daemon=True)
        th_animate_exotic.start()
    else:
        done_flag_animate_exotic = True
        if th_animate_exotic is not None:  # kill animate thread before proceeding
            while th_animate_exotic.is_alive() and counter_wait <= 30:
                time.sleep(.37)
                counter_wait += 1
            th_animate_exotic.join()  # lightest solution
        th_animate_exotic = None
    return not done_flag_animate_exotic


if __name__ == "__main__":
    print("Importing Python Packages - please wait.")
    animate_toggle(True)

# ########## IMPORTS -- PRELOAD ANIMATION END   ##########

import argparse
from datetime import datetime
# Image alignment import
import astroalign as aa

aa.PIXEL_TOL = 1
# aa.NUM_NEAREST_NEIGHBORS=10
# astropy imports
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.io import fits
import astropy.time
from astropy.visualization import astropy_mpl_style, ZScaleInterval, ImageNormalize
from astropy.visualization.stretch import LinearStretch, SquaredStretch, SqrtStretch, LogStretch
from astropy.wcs import WCS, FITSFixedWarning
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astroscrappy import detect_cosmics
# UTC to BJD converter import
from barycorrpy import utc_tdb
# julian conversion imports
import dateutil.parser as dup
# Nested Sampling imports
import dynesty
import glob as g
from pathlib import Path
from io import StringIO
import pyvo as vo
import json
import logging
from logging.handlers import TimedRotatingFileHandler
from matplotlib.animation import FuncAnimation
# Pyplot imports
import matplotlib.pyplot as plt

import numpy as np
import os
# data processing
import pandas
# photometry
from photutils import CircularAperture
from photutils import aperture_photometry
import requests
# scipy imports
from scipy.ndimage import median_filter, generic_filter
from scipy.optimize import least_squares
from scipy.stats import mode
from scipy.signal import savgol_filter
from scipy.ndimage import binary_dilation, label
# cross correlation imports
from skimage.registration import phase_cross_correlation
from skimage.transform import rescale
# error handling for scraper
from tenacity import retry, retry_if_exception_type, retry_if_result, \
    stop_after_attempt, wait_exponential
import warnings

# ########## EXOTIC imports ##########
try:  # science constants
    from constants import *
except ImportError:
    from .constants import *
try:  # light curve numerics
    from .api.elca import lc_fitter, binner, transit
except ImportError:  # package import
    from api.elca import lc_fitter, binner, transit
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


def sigma_clip(ogdata, sigma=3, dt=21):
    nanmask = np.isnan(ogdata)
    mdata = savgol_filter(ogdata[~nanmask], dt, 2)
    # mdata = median_filter(ogdata[~nanmask], dt)
    res = ogdata[~nanmask] - mdata
    std = np.nanmedian([np.nanstd(np.random.choice(res, 25)) for i in range(100)])
    # std = np.nanstd(res) # biased from large outliers
    sigmask = np.abs(res) > sigma * std
    nanmask[~nanmask] = sigmask
    return nanmask


def dms_to_dd(dms_in):
    """
    Quick helper method to convert long/lat values in degree-minute-second (dms) form
    (using ':' separators) to decimal (dd) form
    :param dms_in: DMS long/lat value, colon separated
    :return float: Properly signed long/lat value in decimal float form
    """
    if dms_in is None or isinstance(dms_in, str) is False or str(dms_in).count(":") != 2:
        raise ValueError("Invalid DMS input provided for calculations. ...")
    # clean string of errant leading/trailing/internal spaces
    dms = str(dms_in).strip().replace(" ", "")
    degrees, minutes, seconds = dms.split(":")
    dec = abs(float(degrees)) + float(minutes) / 60. + float(seconds) / 3600.
    if float(degrees) < 0.:
        dec = dec * -1.
    return dec


# ################### START ARCHIVE SCRAPER (PRIORS) ##########################
class NASAExoplanetArchive:

    def __init__(self, planet=None, candidate=False):
        self.planet = planet
        # self.candidate = candidate
        self.pl_dict = None

        # CONFIGURATIONS
        self.requests_timeout = 16, 512  # connection timeout, response timeout in secs.

    def planet_info(self, fancy=False):
        log.info(f"\nLooking up {self.planet}. Please wait. ...")
        self.planet, candidate = self._new_scrape(filename="eaConf.json", target=self.planet)

        if not candidate:
            with open("eaConf.json", "r") as confirmed:
                data = json.load(confirmed)
                planets = [data[i]['pl_name'] for i in range(len(data))]
                idx = planets.index(self.planet)
                self._get_params(data[idx])
                log.info(f"Successfully found {self.planet} in the NASA Exoplanet Archive!")

        # fancy keys matching inits fil
        if fancy:
            coord = SkyCoord(ra=self.pl_dict['ra'] * u.degree, dec=self.pl_dict['dec'] * u.degree)
            rastr = coord.to_string('hmsdms', sep=':').split(' ')[0]
            decstr = coord.to_string('hmsdms', sep=':').split(' ')[1]

            flabels = {
                "Target Star RA": rastr,
                "Target Star Dec": decstr,
                "Planet Name": self.pl_dict['pName'],
                "Host Star Name": self.pl_dict['sName'],
                "Orbital Period (days)": self.pl_dict['pPer'],
                "Orbital Period Uncertainty": self.pl_dict['pPerUnc'],
                "Published Mid-Transit Time (BJD-UTC)": self.pl_dict['midT'],
                "Mid-Transit Time Uncertainty": self.pl_dict['midTUnc'],
                "Ratio of Planet to Stellar Radius (Rp/Rs)": self.pl_dict['rprs'],
                "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": self.pl_dict['rprsUnc'],
                "Ratio of Distance to Stellar Radius (a/Rs)": self.pl_dict['aRs'],
                "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": self.pl_dict['aRsUnc'],
                "Orbital Inclination (deg)": self.pl_dict['inc'],
                "Orbital Inclination (deg) Uncertainity": self.pl_dict['incUnc'],
                "Orbital Eccentricity (0 if null)": self.pl_dict['ecc'],
                "Star Effective Temperature (K)": self.pl_dict['teff'],
                "Star Effective Temperature (+) Uncertainty": self.pl_dict['teffUncPos'],
                "Star Effective Temperature (-) Uncertainty": self.pl_dict['teffUncNeg'],
                "Star Metallicity ([FE/H])": self.pl_dict['met'],
                "Star Metallicity (+) Uncertainty": self.pl_dict['metUncPos'],
                "Star Metallicity (-) Uncertainty": self.pl_dict['metUncNeg'],
                "Star Surface Gravity (log(g))": self.pl_dict['logg'],
                "Star Surface Gravity (+) Uncertainty": self.pl_dict['loggUncPos'],
                "Star Surface Gravity (-) Uncertainty": self.pl_dict['loggUncNeg']
            }

            return json.dumps(flabels, indent=4)

        return self.planet, candidate, self.pl_dict

    @staticmethod
    def dataframe_to_jsonfile(dataframe, filename):
        jsondata = json.loads(dataframe.to_json(orient='table', index=False))
        with open(filename, "w") as f:
            f.write(json.dumps(jsondata['data'], indent=4))

    def _tap_query(self, base_url, query, dataframe=True):
        # table access protocol query

        # build url
        uri_full = base_url
        for k in query:
            if k != "format":
                uri_full += f"{k} {query[k]} "

        uri_full = f"{uri_full[:-1]} &format={query.get('format', 'csv')}"
        uri_full = uri_full.replace(' ', '+')
        log.info(uri_full)

        response = requests.get(uri_full, timeout=self.requests_timeout)
        # TODO check status_code?

        if dataframe:
            return pandas.read_csv(StringIO(response.text))
        else:
            return response.text

    def resolve_name(self):
        uri_ipac_base = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
        uri_ipac_query = {
            # Table columns: https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
            "select": "pl_name,hostname",
            "from": "ps",  # Table name
            "where": "tran_flag = 1",
            "order by": "pl_pubdate desc",
            "format": "csv"
        }

        if self.planet:
            uri_ipac_query["where"] += " and pl_name = '{}'".format(self.planet)

        default = self._tap_query(uri_ipac_base, uri_ipac_query)

        if len(default) == 0:
            return False
        else:
            return True

    @retry(stop=stop_after_attempt(3),
           wait=wait_exponential(multiplier=1, min=17, max=1024),
           retry=(retry_if_exception_type(requests.exceptions.RequestException) |
                  retry_if_exception_type(ConnectionError)))
    def _new_scrape(self, filename="eaConf.json", target=None):

        # scrape_new()
        uri_ipac_base = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
        uri_ipac_query = {
            # Table columns: https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
            "select": "pl_name,hostname,tran_flag,pl_massj,pl_radj,pl_radjerr1,"
                      "pl_ratdor,pl_ratdorerr1,pl_ratdorerr2,pl_orbincl,pl_orbinclerr1,pl_orbinclerr2,"
                      "pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbeccen,"
                      "pl_orblper,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,"
                      "pl_trandep,pl_trandeperr1,pl_trandeperr2,"
                      "pl_ratror,pl_ratrorerr1,pl_ratrorerr2,"
                      "st_teff,st_tefferr1,st_tefferr2,st_met,st_meterr1,st_meterr2,"
                      "st_logg,st_loggerr1,st_loggerr2,st_mass,st_rad,st_raderr1,ra,dec,pl_pubdate",
            "from": "ps",  # Table name
            "where": "tran_flag = 1 and default_flag = 1",
            "order by": "pl_pubdate desc",
            "format": "csv"
        }

        if target:
            uri_ipac_query["where"] += f" and pl_name = '{target}'"

        default = self._tap_query(uri_ipac_base, uri_ipac_query)

        # fill in missing columns
        uri_ipac_query['where'] = 'tran_flag=1'

        if target:
            uri_ipac_query["where"] += f" and pl_name = '{target}'"

        extra = self._tap_query(uri_ipac_base, uri_ipac_query)

        if len(default) == 0:
            target = input(
                f"Cannot find target ({target}) in NASA Exoplanet Archive. Check case sensitivity and spacing and"
                "\nre-enter the planet's name or type candidate if this is a planet candidate: ")
            if target.strip().lower() == 'candidate':
                target = user_input("\nPlease enter candidate planet's name: ", type_=str)
                return target, True
            else:
                return self._new_scrape(filename="eaConf.json", target=target)
        else:
            # replaces NEA default with most recent publication
            default.iloc[0] = extra.iloc[0]

            # for each planet
            for i in default.pl_name:

                # extract rows for each planet
                ddata = default.loc[default.pl_name == i]
                edata = extra.loc[extra.pl_name == i]

                # for each nan column in default
                nans = ddata.isna()
                for k in ddata.keys():
                    if nans[k].bool():  # if col value is nan
                        if not edata[k].isna().all():  # if replacement data exists
                            # replace with first index
                            default.loc[default.pl_name == i, k] = edata[k][edata[k].notna()].values[0]
                            # TODO could use mean for some variables (not mid-transit)
                            # log.info(i,k,edata[k][edata[k].notna()].values[0])
                        else:
                            # permanent nans - require manual entry
                            if k == 'pl_orblper':  # omega
                                default.loc[default.pl_name == i, k] = 0
                            elif k == 'pl_ratdor':  # a/R*
                                # Kepler's 3rd law
                                semi = SA(ddata.st_mass.values[0], ddata.pl_orbper.values[0])
                                default.loc[default.pl_name == i, k] = semi * AU / (
                                        ddata.st_rad.values[0] * R_SUN)
                            elif k == 'pl_orbincl':  # inclination
                                default.loc[default.pl_name == i, k] = 90
                            elif k == "pl_orbeccen":  # eccentricity
                                default.loc[default.pl_name == i, k] = 0
                            elif k == "st_met":  # [Fe/H]
                                default.loc[default.pl_name == i, k] = 0
                            else:
                                default.loc[default.pl_name == i, k] = 0

            NASAExoplanetArchive.dataframe_to_jsonfile(default, filename)
            return target, False

    def _get_params(self, data):
        # translate data from Archive keys to Ethan Keys
        try:
            rprs = np.sqrt(data['pl_trandep'] / 100.)
            rprserr = np.sqrt(np.abs((data['pl_trandeperr1'] / 100.) * (data['pl_trandeperr2'] / 100.))) / (2. * rprs)
        except (KeyError, TypeError):
            try:
                rprs = data['pl_ratror']
                rprserr = np.sqrt(np.abs(data['pl_ratrorerr1'] * data['pl_ratrorerr2']))
            except (KeyError, TypeError):
                rp = data['pl_radj'] * R_JUP
                rperr = np.sqrt(np.abs(data['pl_radjerr1'] * data['pl_radjerr2'])) * R_JUP
                rs = data['st_rad'] * R_SUN
                rserr = np.sqrt(np.abs(data['st_raderr1'] * data['st_raderr2'])) * R_SUN
                rprserr = ((rperr / rs) ** 2 + (-rp * rserr / rs ** 2) ** 2) ** 0.5
                rprs = rp / rs

        self.pl_dict = {
            'ra': data['ra'],
            'dec': data['dec'],
            'pName': data['pl_name'],
            'sName': data['hostname'],
            'pPer': data['pl_orbper'],
            'pPerUnc': np.sqrt(np.abs(data['pl_orbpererr1'] * data['pl_orbpererr2'])),

            'midT': data['pl_tranmid'],
            'midTUnc': np.sqrt(np.abs(data['pl_tranmiderr1'] * data['pl_tranmiderr2'])),
            'rprs': rprs,
            'rprsUnc': rprserr,
            'aRs': data['pl_ratdor'],
            'aRsUnc': np.sqrt(np.abs(data.get('pl_ratdorerr1', 1) * data['pl_ratdorerr2'])),
            'inc': data['pl_orbincl'],
            'incUnc': np.sqrt(np.abs(data['pl_orbinclerr1'] * data['pl_orbinclerr2'])),

            'ecc': data.get('pl_orbeccen', 0),
            'teff': data['st_teff'],
            'teffUncPos': data['st_tefferr1'],
            'teffUncNeg': data['st_tefferr2'],
            'met': data['st_met'],
            'metUncPos': max(0.01, data['st_meterr1']),
            'metUncNeg': min(-0.01, data['st_meterr2']),
            'logg': data['st_logg'],
            'loggUncPos': data['st_loggerr1'],
            'loggUncNeg': data['st_loggerr2']
        }

        if self.pl_dict['aRsUnc'] == 0:
            self.pl_dict['aRsUnc'] = 0.1

        if self.pl_dict['incUnc'] == 0:
            self.pl_dict['incUnc'] = 0.1


# ################### END ARCHIVE SCRAPER (PRIORS) ############################


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
            gDateTime = "{}T{}".format(header['DATE-OBS'], header['TIME-OBS'])

        dt = dup.parse(gDateTime)
        atime = astropy.time.Time(dt)
        julianTime = atime.jd
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in header.comments['DATE-OBS']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now

    # If the mid-exposure time is given in the fits header, then no offset is needed to calculate the mid-exposure time
    return julianTime + exptime_offset


# Method that gets and returns the current phase of the target
def getPhase(curTime, pPeriod, tMid):
    phase = (curTime - tMid - 0.5 * pPeriod) / pPeriod % 1
    return phase - 0.5


# Method that gets and returns the airmass from the fits file (Really the Altitude)
def getAirMass(hdul, ra, dec, lati, longit, elevation):
    # Determine the correct image header extension to use
    extension = 0
    image_header = hdul[extension].header
    while image_header["NAXIS"] == 0:
        extension += 1
        image_header = hdul[extension].header

    # Grab airmass from image header; if not listed, calculate it from TELALT; if that isn't listed, then calculate it the hard way
    if 'AIRMASS' in hdul[extension].header:
        am = float(hdul[extension].header['AIRMASS'])
    elif 'TELALT' in hdul[extension].header:
        alt = float(hdul[extension].header[
                        'TELALT'])  # gets the airmass from the fits file header in (sec(z)) (Secant of the zenith angle)
        cosam = np.cos((np.pi / 180) * (90.0 - alt))
        am = 1 / cosam
    else:
        # pointing = SkyCoord(str(astropy.coordinates.Angle(raStr+" hours").deg)+" "+str(astropy.coordinates.Angle(decStr+" degrees").deg ), unit=(u.deg, u.deg), frame='icrs')
        pointing = SkyCoord(str(ra) + " " + str(dec), unit=(u.deg, u.deg), frame='icrs')

        location = EarthLocation.from_geodetic(lat=lati * u.deg, lon=longit * u.deg, height=elevation)
        atime = astropy.time.Time(getJulianTime(hdul[extension].header), format='jd', scale='utc', location=location)
        pointingAltAz = pointing.transform_to(AltAz(obstime=atime, location=location))
        am = float(pointingAltAz.secz)
    return am


def user_input(prompt, type_, val1=None, val2=None, val3=None):
    while True:
        try:
            result = type_(input(prompt))
            log.debug(f"{prompt}{result}")
        except ValueError:
            log.info("Sorry, not a valid datatype.")
            continue
        if type_ == str and val1 and val2 and val3:
            result = result.lower().strip()
            if result not in (val1, val2, val3):
                log.info("Sorry, your response was not valid.")
            else:
                return result
        elif type_ == int and val1 and val2 and val3:
            if result not in (val1, val2, val3):
                log.info("Sorry, your response was not valid.")
            else:
                return result
        elif type_ == int and val1 and val2 and val3:
            if result not in (val1, val2, val3):
                log.info("Sorry, your response was not valid.")
            else:
                return result
        else:
            return result


def prereduced_file():
    while True:
        try:
            file = user_input("Enter the path and file name of your data file: ", type_=str)
            if file == "ok":
                file = "/Users/rzellem/Documents/EXOTIC/sample-data/NormalizedFluxHAT-P-32 bDecember 17, 2017.txt"
                # file = "/Users/rzellem/Downloads/fluxorama.csv
                log.info("Hello, Rob.")

            file = Path(file)

            if file.is_file():
                return file
            else:
                raise FileNotFoundError
        except FileNotFoundError:
            log.info("Error: Data file not found. Please try again.")


def save_directory(directory):
    while True:
        try:
            if not directory:
                directory = user_input("Enter the directory to save the results and plots into "
                                       "or type new to create one: ", type_=str)
            if directory == 'new':
                directory = create_directory()
            else:
                if not Path(directory).is_dir():
                    raise NotADirectoryError
            return directory
        except (NotADirectoryError, OSError):
            log.info("Error: the directory entered does not exist. Please try again. Make sure to follow this "
                     "\nformatting (using whichever directory you choose): /sample-data/results")
            directory = None


# Create a save directory within the current working directory
def create_directory():
    save_path = Path.cwd()
    while True:
        directory = user_input("Enter the name for your new directory: ", type_=str)
        try:
            save_path = save_path / directory
            Path(save_path).mkdir()
        except OSError:
            log.info(f"Creation of the directory {save_path}/{directory} failed.")
        else:
            log.info(f"Successfully created the directory {save_path}.")
            return save_path


# Check user's inits.json for user information and planetary parameters
def check_init_file(init, dict_info, dict_params):
    with init.open('r') as json_file:
        data = json.load(json_file)

    comparison_info = {'fitsdir': 'Directory with FITS files',
                       'saveplot': 'Directory to Save Plots',
                       'flatsdir': 'Directory of Flats',
                       'darksdir': 'Directory of Darks',
                       'biasesdir': 'Directory of Biases',
                       'aavsonum': 'AAVSO Observer Code (N/A if none)',
                       'secondobs': 'Secondary Observer Codes (N/A if none)',
                       'date': 'Observation date',
                       'lat': 'Obs. Latitude',
                       'long': 'Obs. Longitude',
                       'elev': 'Obs. Elevation (meters)',
                       'ctype': 'Camera Type (CCD or DSLR)',
                       'pixelbin': 'Pixel Binning',
                       'filter': 'Filter Name (aavso.org/filters)',
                       'wl_min': 'Filter Minimum Wavelength (nm)',
                       'wl_max': 'Filter Maximum Wavelength (nm)',
                       'notes': 'Observing Notes',
                       'tarcoords': 'Target Star X & Y Pixel',
                       'compstars': 'Comparison Star(s) X & Y Pixel',
                       'plate_opt': 'Plate Solution? (y/n)',
                       'pixel_scale': 'Pixel Scale (Ex: 5.21 arcsecs/pixel)',
                       'image_align': 'Align Images? (y/n)'}

    comparison_parameters = {'ra': 'Target Star RA',
                             'dec': 'Target Star Dec',
                             'pName': "Planet Name",
                             'sName': "Host Star Name",
                             'pPer': 'Orbital Period (days)',
                             'pPerUnc': 'Orbital Period Uncertainty',
                             'midT': 'Published Mid-Transit Time (BJD-UTC)',
                             'midTUnc': 'Mid-Transit Time Uncertainty',
                             'rprs': 'Ratio of Planet to Stellar Radius (Rp/Rs)',
                             'rprsUnc': 'Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty',
                             'aRs': 'Ratio of Distance to Stellar Radius (a/Rs)',
                             'aRsUnc': 'Ratio of Distance to Stellar Radius (a/Rs) Uncertainty',
                             'inc': 'Orbital Inclination (deg)',
                             'incUnc': 'Orbital Inclination (deg) Uncertainity',
                             'ecc': 'Orbital Eccentricity (0 if null)',
                             'teff': 'Star Effective Temperature (K)',
                             'teffUncPos': 'Star Effective Temperature (+) Uncertainty',
                             'teffUncNeg': 'Star Effective Temperature (-) Uncertainty',
                             'met': 'Star Metallicity ([FE/H])',
                             'metUncPos': 'Star Metallicity (+) Uncertainty',
                             'metUncNeg': 'Star Metallicity (-) Uncertainty',
                             'logg': 'Star Surface Gravity (log(g))',
                             'loggUncPos': 'Star Surface Gravity (+) Uncertainty',
                             'loggUncNeg': 'Star Surface Gravity (-) Uncertainty'}

    dict_info = get_init_params(comparison_info, dict_info, data['user_info'])
    dict_info = get_init_params(comparison_info, dict_info, data['optional_info'])
    dict_params = get_init_params(comparison_parameters, dict_params, data['planetary_parameters'])

    return dict_info, dict_params


def get_init_params(comp, dict1, dict2):
    for key, value in comp.items():
        try:
            dict1[key] = dict2[value]
        except KeyError:
            pass
    return dict1


# Get inits.json file from user input
def get_initialization_file(info_dict, user_pdict, args_init):
    cwd = Path.cwd()
    log.info(f"\nYour current working directory is: {cwd}")
    log.info(f"Potential initialization files I've found in {cwd} are: ")
    [log.info(f"\t{file}") for file in cwd.glob('*.json') if file.is_file()]

    while True:
        try:
            if not args_init:
                init_file = user_input("\nPlease enter the Directory and Filename of your Initialization File: ",
                                       type_=str)
            else:
                init_file = args_init
            if init_file == 'ok':
                init_file = '/Users/rzellem/Documents/EXOTIC/inits.json'
            return check_init_file(Path(init_file), info_dict, user_pdict)
        except FileNotFoundError:
            log.info("Error: Initialization file not found. Please try again.")
        except IsADirectoryError:
            log.info("Error: Entered a directory. Please try again.")
        finally:
            args_init = None


class InitializationFile:

    def __init__(self, info, planet_name=None):
        self.info = info
        self.planet_name = planet_name

    def get_info(self):
        if self.info['fitsdir'] is None:
            self.image_directory()
        if self.info['saveplot'] is None:
            self.save_directory()
        if self.info['aavsonum'] or self.info['secondobs'] or self.info['ctype'] or self.info['pixelbin'] \
                or self.info['filter'] or self.info['notes'] is None:
            self.initial()
        if self.planet_name is None:
            self.planet()
        if self.info['lat'] is None:
            self.latitude()
        if self.info['long'] is None:
            self.longitude()
        if self.info['tarcoords'] is None:
            self.target_star_coords()
        if self.info['compstars'] is None:
            self.comparison_star_coords()
        return self.info, self.planet_name

    def image_directory(self):
        self.info['fitsdir'] = user_input("Please enter the Directory of Imaging Files: ", type_=str)

    def save_directory(self):
        self.info['saveplot'] = user_input("Please enter the directory to save the results and plots into "
                                           "or type new to create one: ", type_=str)

    def initial(self):
        notes = ['Please enter your AAVSO Observer Account Number (type N/A if you do not currently have an account): ',
                 'Please enter your comma-separated secondary observer codes (type N/A if only 1 observer code): ',
                 'Please enter the camera type (CCD or DSLR): ',
                 'Please enter the pixel binning: ',
                 'Please enter the filter name: ',
                 'Please enter any observing notes (seeing, weather, etc.): ']
        i = 0

        for key, value in self.info.items():
            if key in ('aavsonum', 'secondobs', 'ctype', 'pixelbin', 'filter', 'notes') and value is None:
                self.info[key] = input(notes[i])
                i += 1

    def latitude(self):
        self.info['latitude'] = user_input("Please enter the latitude of where you observed (deg) "
                                           "(Don't forget the sign where North is '+' and South is '-'): ", type_=str)

    def longitude(self):
        self.info['longitude'] = user_input("Please enter the longitude of where you observed (deg) "
                                            "(Don't forget the sign where East is '+' and West is '-'): ", type_=str)

    def target_star_coords(self, pname):
        x_pix = user_input(f"\n{pname} X Pixel Coordinate: ", type_=int)
        y_pix = user_input(f"\n{pname} Y Pixel Coordinate: ", type_=int)

        self.info['tarcoords'] = [x_pix, y_pix]

    def comparison_star_coords(self):
        num_comp_stars = user_input("How many comparison stars would you like to use? (1-10): ", type_=int)
        comp_stars = []

        for num in range(num_comp_stars):
            x_pix = user_input(f"Comparison Star {num + 1} X Pixel Coordinate: ", type_=int)
            y_pix = user_input(f"Comparison Star {num + 1} Y Pixel Coordinate: ", type_=int)
            comp_stars.append((x_pix, y_pix))

        self.info['compstars'] = comp_stars

    def exposure(self):
        self.info['exposure'] = user_input("Please enter your exposure time (seconds): ", type_=int)

    def pixel_scale(self):
        self.info['flatsdir'] = user_input("Please enter the size of your pixel (e.g., 5 arcsec/pixel): ", type_=str)

    def planet(self):
        while True:
            self.planet_name = user_input("\nPlease enter the Planet Name. Make sure it matches the case sensitive "
                                          "name and spacing used on Exoplanet Archive "
                                          "(https://exoplanetarchive.ipac.caltech.edu/index.html): ", type_=str)

            if not self.planet_name[-2].isspace():
                log.info("The convention on the NASA Exoplanet Archive "
                         "(https://exoplanetarchive.ipac.caltech.edu/index.html) \nis to have a space between the star "
                         "name and the planet letter. \nPlease confirm that you have properly input the planet's name."
                         "\nPlease confirm:"
                         f"\n  (1) {self.planet_name} is correct."
                         "\n  (2) The planet name needs to be changed.")
                planetnameconfirm = user_input("\nPlease select 1 or 2: ", type_=int, val1=1, val2=2)
            else:
                break
            if planetnameconfirm == 1:
                break


def planet_name(planet):
    while True:
        if not planet:
            planet = user_input("\nEnter the Planet Name. Make sure it matches the case sensitive name "
                                "and spacing used on Exoplanet Archive "
                                "(https://exoplanetarchive.ipac.caltech.edu/index.html): ", type_=str)

        if not planet[-2].isspace():
            log.info("The convention on the NASA Exoplanet Archive "
                     "(https://exoplanetarchive.ipac.caltech.edu/index.html) is to have a space between the "
                     "star name and the planet letter. Please confirm that you have properly "
                     "input the planet's name."
                     "\nPlease confirm:"
                     f"\n  (1) {planet} is correct."
                     "\n  (2) The planet name needs to be changed.")

            opt = user_input("\nPlease select 1 or 2: ", type_=int, val1=1, val2=2)
        else:
            opt = 1

        if opt == 1:
            return planet
        planet = None


def obs_date(date):
    while True:
        try:
            if not date:
                date = user_input("\nPlease enter the Observation Date (DD-Month-YYYY): ", type_=str)
            if date != datetime.strptime(date, '%d-%B-%Y').strftime('%d-%B-%Y'):
                raise ValueError
            return date
        except ValueError:
            log.info('\nThe entered Observation Date format is incorrect.')
            date = None


def latitude(lat):
    while True:
        try:
            if not lat:
                lat = user_input("Enter the latitude (in degrees) of where you observed. "
                                 "Don't forget the sign where North is '+' and South is '-'! "
                                 "(Example: +50.4): ", type_=str)

            lat = lat.replace(' ', '')
            if lat[0] != '+' and lat[0] != '-':
                raise ValueError("You forgot the sign for the latitude! North is '+' and South is '-'. "
                                 "Please try again.")

            # Convert to float if latitude in decimal. If latitude is in +/-HH:MM:SS format, convert to a float.
            try:
                lat = float(lat)
            except ValueError:
                lat = float(dms_to_dd(lat))

            if lat <= -90.00 or lat >= 90.00:
                raise ValueError("Your latitude is out of range. Please enter a latitude between -90 and +90 (deg)")
            return lat
        except ValueError as err:
            log.info(err.args)
            lat = None


def longitude(long):
    while True:
        try:
            if not long:
                long = user_input("Enter the longitude (in degrees) of where you observed. "
                                  "(Don't forget the sign where East is '+' and West is '-')! "
                                  "(Example: -32.12): ", type_=str)

            long = long.replace(' ', '')
            if long[0] != '+' and long[0] != '-':
                raise ValueError("You forgot the sign for the longitude! East is '+' and West is '-'. "
                                 "Please try again.")

            # Convert to float if longitude in decimal. If longitude is in +/-HH:MM:SS format, convert to a float.
            try:
                long = float(long)
            except ValueError:
                long = float(dms_to_dd(long))

            if long <= -180.00 or long >= 180.00:
                raise ValueError("Your longitude is out of range. Please enter a longitude between -180 and +180 (deg)")
            return long
        except ValueError as err:
            log.info(err.args)
            long = None


def elevation(elev, lat, long):
    while True:
        try:
            if not elev:
                elev = open_elevation(lat, long)
                if not elev:
                    log.info("EXOTIC could not retrieve elevation.")
                    elev = user_input("Enter the elevation (in meters) of where you observed: ", type_=float)
            return float(elev)
        except ValueError:
            log.info("The entered elevation is incorrect.")
            elev = None


def is_false(value):
    return value is False


def result_if_max_retry_count(retry_state):
    pass


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
       retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
       retry_error_callback=result_if_max_retry_count)
def open_elevation(lat, long):
    query = f"https://api.open-elevation.com/api/v1/lookup?locations={lat},{long}"

    try:
        r = requests.get(query).json()
        return r['results'][0]['elevation']
    except requests.exceptions.RequestException:
        return False


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
    uncert = 20 / 3600

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
        log.info("\nDifference(s) found between initialization file parameters and "
                 "those scraped by EXOTIC from the NASA Exoplanet Archive."
                 "\n Would you like:"
                 "\n  (1) EXOTIC to adopt of all of your defined parameters or"
                 "\n  (2) to review the ones scraped from the Archive that differ?")
        opt = user_input("\nPlease enter 1 or 2: ", type_=str, val1='1', val2='2')

        if opt == '2':
            return True
        else:
            return False


# --------PLANETARY PARAMETERS UI------------------------------------------
# Get the user's confirmation of values that will later be used in lightcurve fit
def get_planetary_parameters(candplanetbool, userpdict, pdict=None):
    log.info("*******************************************")
    log.info("Planetary Parameters for Lightcurve Fitting")

    # The order of planet_params list must match the pDict that is declared when scraping the NASA Exoplanet Archive
    planet_params = ["Target Star RA in the form: HH:MM:SS (ignore the decimal values)",
                     "Target Star DEC in form: <sign>DD:MM:SS (ignore the decimal values and don't forget the '+' or '-' sign!)",
                     "Planet's Name",
                     "Host Star's Name",
                     "Orbital Period (days)",
                     "Orbital Period Uncertainty (days) \n(Keep in mind that 1.2e-34 is the same as 1.2 x 10^-34)",
                     "Published Mid-Transit Time (BJD_UTC)",
                     "Mid-Transit Time Uncertainty (JD)",
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
                log.info(f"\n\n*** WARNING: {pdict['pName']} initialization file's {planet_params[idx]} does not match "
                         "the value scraped by EXOTIC from the NASA Exoplanet Archive. ***\n")
                log.info(f"\tNASA Exoplanet Archive value (degrees): {pdict[item]}")
                log.info(f"\tInitialization file value (degrees): {userpdict[item]}")
                log.info("\nWould you like to:"
                         "\n  (1) use NASA Exoplanet Archive value, "
                         "\n  (2) use initialization file value, or "
                         "\n  (3) enter in a new value.")
                option = user_input("Which option do you choose? (1/2/3): ", type_=int, val1=1, val2=2, val3=3)

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
        log.info(f"*** Here are the values scraped from the NASA Exoplanet Archive for {pdict['pName']} that were not "
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
                log.info(f"\n\n*** WARNING: {pdict['pName']} initialization file's {planet_params[i]} does not match "
                         "the value scraped by EXOTIC from the NASA Exoplanet Archive. ***\n")
                log.info(f"\tNASA Exoplanet Archive value: {pdict[key]}")
                log.info(f"\tInitialization file value: {userpdict[key]}")
                log.info("\nWould you like to: "
                         "\n  (1) use NASA Exoplanet Archive value, "
                         "\n  (2) use initialization file value, or "
                         "\n  (3) enter in a new value.")
                option = user_input("Which option do you choose? (1/2/3): ", type_=int, val1=1, val2=2, val3=3)
                if option == 1:
                    userpdict[key] = pdict[key]
                elif option == 2:
                    continue
                else:
                    userpdict[key] = user_input(f"Enter the {planet_params[i]}: ", type_=type(userpdict[key]))
            # Did not use initialization file or null
            else:
                log.info(f"\n {pdict['pName']} {planet_params[i]}: {pdict[key]}")
                agreement = user_input("Do you agree? (y/n): ", type_=str, val1='y', val2='n')
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
                                       type_=str, val1='y', val2='n')
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
            log.info("Error: The format is not correct, please try again.")
            ra = input("Input the right ascension of target (HH:MM:SS): ")
            dec = input("Input the declination of target (<sign>DD:MM:SS): ")


def round_to_2(*args):
    x = args[0]
    if len(args) == 1:
        y = args[0]
    else:
        y = args[1]
    if np.floor(y) >= 1. or y == 0.0:
        roundval = 2
    else:
        roundval = -int(np.floor(np.log10(abs(y)))) + 1
    return round(x, roundval)


# Check if user's directory contains imaging FITS files that are able to be reduced
def check_imaging_files(directory, filename):
    file_extensions = ['.fits', '.fit', '.fts', '.fz']
    input_files = []

    while True:
        try:
            if Path(directory).is_dir():
                directory = Path(directory)
                for ext in file_extensions:
                    for file in directory.iterdir():
                        if file.is_file() and file.name.lower().endswith(ext.lower()) and file.name[0:2] not in (
                        'ref', 'wcs'):
                            input_files.append(str(file))
                    if input_files:
                        return directory, input_files
                if not input_files:
                    raise FileNotFoundError
            else:
                raise NotADirectoryError
        except FileNotFoundError:
            opt = user_input(f"\nError: {filename} files not found with .fits, .fit or .fts extensions in {directory}."
                             "\nWould you like to enter in an alternate image extension in addition to .FITS? (y/n): ",
                             type_=str, val1='y', val2='n')
            if opt == 'y':
                add_ext = user_input("Please enter the extension you want to add (EX: .FITS): ", type_=str)
                file_extensions.append(add_ext)
            else:
                directory = user_input(f"Enter the directory path where {filename} files are located: ", type_=str)
        except (NotADirectoryError, OSError):
            log.info("\nError: No such directory exists when searching for FITS files. Please try again.")
            directory = user_input(f"Enter the directory path where {filename} files are located: ", type_=str)


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
        self._standard_list()

        if self.filter_type and not (self.wl_min or self.wl_max):
            self._standard()
        elif self.wl_min or self.wl_max:
            self._custom()
        else:
            opt = user_input("\nWould you like EXOTIC to calculate your limb darkening parameters "
                             "with uncertainties? (y/n):", type_=str, val1='y', val2='n')

            if opt == 'y':
                opt = user_input("Please enter 1 to use a standard filter or 2 for a customized filter:",
                                 type_=int, val1=1, val2=2)
                if opt == 1:
                    self._standard()
                elif opt == 2:
                    self._custom()
            else:
                self._user_entered()
        return self.ld0, self.ld1, self.ld2, self.ld3, self.filter_type, self.wl_min * 1000, self.wl_max * 1000

    def _standard_list(self):
        log.info("\n\n***************************")
        log.info("Limb Darkening Coefficients")
        log.info("***************************")
        log.info("\nThe standard bands that are available for limb darkening parameters (https://www.aavso.org/filters)"
                 "\nas well as filters for MObs and LCO (0.4m telescope) datasets:\n")
        for key, value in self.fwhm.items():
            log.info(f"\t{key[1]}: {key[0]} - ({value[0]:.2f}-{value[1]:.2f}) nm")

    def _standard(self):
        while True:
            try:
                if not self.filter_type:
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
                log.info("\nError: The entered filter is not in the provided list of standard filters.")
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
            log.info(f"Found corrupted file and removing from reduction: {file}, ({e})")
        finally:
            if getattr(hdul, "close", None) and callable(hdul.close):
                hdul.close()
            del hdul

    return valid_files


def check_wcs(fits_file, save_directory, plate_opt):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=FITSFixedWarning)
        header = fits.getheader(filename=fits_file)
        wcs_header = WCS(header)
        wcs_exists = wcs_header.is_celestial

    if plate_opt not in ('y', 'n'):
        plate_opt = user_input("\nWould you like to upload the your image for a plate solution?"
                               "\nDISCLAIMER: One of your imaging files will be publicly viewable "
                               "on nova.astrometry.net. (y/n): ", type_=str, val1='y', val2='n')

    if plate_opt == 'y':
        return get_wcs(fits_file, save_directory)
    elif plate_opt == 'n':
        if wcs_exists:
            log.info("Your FITS files have WCS information in their headers. EXOTIC will proceed to use these. "
                     "NOTE: If you do not trust your WCS coordinates, "
                     "please restart EXOTIC after enabling plate solutions via astrometry.net.")
            return fits_file
        else:
            return False


def get_wcs(file, directory=""):
    log.info("\nGetting the plate solution for your imaging file. Please wait. ...")
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


# Check the ra and dec against the plate solution to see if the user entered in the correct values
def check_targetpixelwcs(pixx, pixy, expra, expdec, ralist, declist):
    while True:
        try:
            uncert = 20 / 3600
            # Margins are within 20 arcseconds
            if expra - uncert >= ralist[pixy][pixx] or ralist[pixy][pixx] >= expra + uncert:
                log.info("Error: The X Pixel Coordinate entered does not match the target's right ascension.")
                raise ValueError
            if expdec - uncert >= declist[pixy][pixx] or declist[pixy][pixx] >= expdec + uncert:
                log.info("Error: The Y Pixel Coordinate entered does not match the target's declination.")
                raise ValueError
            return pixx, pixy
        except ValueError:
            opt = user_input("Would you like to re-enter the pixel coordinates? (y/n): ", type_=str, val1='y', val2='n')

            # User wants to change their coordinates
            if opt == 'y':
                # Checks for the closest pixel location in ralist and declist for expected ra and dec
                dist = (ralist - expra) ** 2 + (declist - expdec) ** 2
                pixy, pixx = np.unravel_index(dist.argmin(), dist.shape)
                searchopt = user_input(f"Here are the suggested pixel coordinates:"
                                       f"  X Pixel: {pixx}"
                                       f"  Y Pixel: {pixy}"
                                       "\nWould you like to use these? (y/n): ",
                                       type_=str, val1='y', val2='n')
                # Use the coordinates found by code
                if searchopt == 'y':
                    return pixx, pixy
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

    # Query GAIA first to check for variability using the phot_variable_flag trait
    radius = u.Quantity(20.0, u.arcsec)
    try:
        gaiaQuery = Gaia.cone_search_async(sample, radius)
        gaiaResult = gaiaQuery.get_results()
    except Exception:
        log.info("Not able to query information from Simbad.")
        return False

    # Individually go through the phot_variable_flag indicator for each star to see if variable or not
    variableFlagList = gaiaResult.columns["phot_variable_flag"]
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

    # query SIMBAD and search identifier result table to determine if comparison star is variable in any form
    # This is a secondary check if GAIA query returns inconclusive results
    simbad_result = Simbad.query_region(sample, radius=20 * u.arcsec)
    try:
        starName = simbad_result['MAIN_ID'][0].decode("utf-8")
    except:
        log.info("Your star cannot be resolved in SIMBAD. Proceed with caution.")
        return False
    identifiers = Simbad.query_objectids(starName)
    for currName in identifiers:
        if "V*" in currName[0]:
            return True
    return False


# Aligns imaging data from .fits file to easily track the host and comparison star's positions
def transformation(image_data, num_images, file_name, count, roi=1):
    pos = np.zeros((1, 2))

    # crop image to ROI
    height = image_data.shape[1]
    width = image_data.shape[2]
    roix = slice(int(width * (0.5 - roi / 2)), int(width * (0.5 + roi / 2)))
    roiy = slice(int(height * (0.5 - roi / 2)), int(height * (0.5 + roi / 2)))

    # scale data
    npix = np.product(image_data[0].shape)
    if npix > 16e6: # 16 Mega Pixel
        sf = 0.25
    elif npix > 7.9e6: # ~8 MP
        sf = 0.33
    elif npix > 1e6: # more than 1 Mega Pixel
        sf = 0.5
    else:
        sf = 1

    if sf != 1:
        image_data = np.array([
            rescale(image_data[0], sf),
            rescale(image_data[1], sf),
        ])

    # Find transformation from .FITS files and catch exceptions if not able to.
    try:
        sys.stdout.write(f"Finding transformation {count + 1} of {num_images}\r")
        log.debug(f"Finding transformation {count + 1} of {num_images}\r")
        sys.stdout.flush()

        results = aa.find_transform(image_data[1][roiy, roix], image_data[0][roiy, roix])
        rot = results[0].rotation
        pos[0] = results[0].translation/sf

    except Exception as ee:
        log.info(ee)
        log.info(f"Image {count + 1} of {num_images} failed to align, passing on image: {file_name}")

        rot = 0
        pos = np.zeros((1, 2))

    return pos, rot


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
        image_scale = f"Image scale in arc-secs/pixel: {pixel_init}"
    else:
        log.info("Cannot find the pixel scale in the image header.")
        image_scale_num = user_input("Please enter the size of your pixel (e.g., 5 arc-sec/pixel): ", type_=float)
        image_scale = f"Image scale in arc-secs/pixel: {image_scale_num}"
    return image_scale


# Will remove later from code as these are older metadata formatting replaced by -XC. Kept for compatibility
def previous_data_format(pdict, ld_0, ld_1, ld_2, ld_3, my_fit):
    return (f"#FILTER={exotic_infoDict['filter']}\n"
            f"#PRIORS=Period={round_to_2(pdict['pPer'], pdict['pPerUnc'])} +/- {round_to_2(pdict['pPerUnc'])},a/R*={round_to_2(pdict['aRs'], pdict['aRsUnc'])} +/- {round_to_2(pdict['aRsUnc'])}"
            f",inc={round_to_2(pdict['inc'], pdict['incUnc'])} +/- {round_to_2(pdict['incUnc'])},ecc={round_to_2(pdict['ecc'])}"
            f",u0={round_to_2(ld_0[0], ld_0[1])} +/- {round_to_2(ld_0[1])},u1={round_to_2(ld_1[0], ld_1[1])} +/- {round_to_2(ld_1[1])},u2={round_to_2(ld_2[0], ld_2[1])} +/- {round_to_2(ld_2[1])}"
            f",u3={round_to_2(ld_3[0], ld_3[1])} +/- {round_to_2(ld_3[1])}\n"
            f"#RESULTS=Tc={round_to_2(my_fit.parameters['tmid'], my_fit.errors['tmid'])} +/- {round_to_2(my_fit.errors['tmid'])}"
            f",Rp/R*={round_to_2(my_fit.parameters['rprs'], my_fit.errors['rprs'])} +/- {round_to_2(my_fit.errors['rprs'])}"
            f",Am1={round_to_2(my_fit.parameters['a1'], my_fit.errors['a1'])} +/- {round_to_2(my_fit.errors['a1'])}"
            f",Am2={round_to_2(my_fit.parameters['a2'], my_fit.errors['a2'])} +/- {round_to_2(my_fit.errors['a2'])}\n")


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
        print("")
        print("Target Location:", np.round(pixcoord[0], 2))

    return pixcoord[0]


# defines the star point spread function as a 2D Gaussian
def star_psf(x, y, x0, y0, a, sigx, sigy, b):
    gaus = a * np.exp(-(x - x0) ** 2 / (2 * sigx ** 2)) * np.exp(-(y - y0) ** 2 / (2 * sigy ** 2)) + b
    return gaus


# Method uses the Full Width Half Max to estimate the standard deviation of the star's psf
def estimate_sigma(x, maxidx=-1):
    if maxidx == -1:
        maxidx = np.argmax(x)
    lower = np.abs(x - 0.5 * np.nanmax(x))[:maxidx].argmin()
    upper = np.abs(x - 0.5 * np.nanmax(x))[maxidx:].argmin() + maxidx
    FWHM = upper - lower
    return FWHM / (2 * np.sqrt(2 * np.log(2)))


def gaussian_psf(x, y, x0, y0, a, sigx, sigy, rot, b):
    rx = (x - x0) * np.cos(rot) - (y - y0) * np.sin(rot)
    ry = (x - x0) * np.sin(rot) + (y - y0) * np.cos(rot)
    gausx = np.exp(-rx ** 2 / (2 * sigx ** 2))
    gausy = np.exp(-ry ** 2 / (2 * sigy ** 2))
    return a * gausx * gausy + b


def fit_psf(data, pos, init, lo, up, psf_function=gaussian_psf, lossfn='linear', method='trf', box=10):
    xv, yv = mesh_box(pos, box)

    def fcn2min(pars):
        model = psf_function(xv, yv, *pars)
        return (data[yv, xv] - model).flatten()

    if method == 'trf':
        res = least_squares(fcn2min, x0=[*pos, *init], bounds=[lo, up], loss=lossfn, jac='3-point', method='dogbox',
                            xtol=None, ftol=1e-3, tr_options='exact')
    else:
        res = least_squares(fcn2min, x0=[*pos, *init], loss=lossfn, jac='3-point', method=method)
    return res.x


def mesh_box(pos, box):
    pos = [int(np.round(pos[0])), int(np.round(pos[1]))]
    x = np.arange(pos[0] - box, pos[0] + box + 1)
    y = np.arange(pos[1] - box, pos[1] + box + 1)
    xv, yv = np.meshgrid(x, y)
    return xv.astype(int), yv.astype(int)


# Method fits a 2D gaussian function that matches the star_psf to the star image and returns its pixel coordinates
def fit_centroid(data, pos, init=[], box=10, debug=False):
    # get sub field in image
    xv, yv = mesh_box(pos, box)

    # weighted flux centroid
    wfx = np.sum(np.unique(xv) * (data[yv, xv].sum(0) - data[yv, xv].sum(0).min())) / np.sum(
        (data[yv, xv].sum(0) - data[yv, xv].sum(0).min()))
    wfy = np.sum(np.unique(yv) * (data[yv, xv].sum(1) - data[yv, xv].sum(1).min())) / np.sum(
        (data[yv, xv].sum(1) - data[yv, xv].sum(1).min()))

    if len(init) == 5:
        pass
    else:
        init = [np.nanmax(data[yv, xv]) - np.nanmin(data[yv, xv]), 1, 1, 0, np.nanmin(data[yv, xv])]

    try:
        # fit gaussian PSF
        pars = fit_psf(
            data,
            [wfx, wfy],  # position estimate
            init,  # initial guess: [amp, sigx, sigy, rotation, bg]
            [wfx - box * 0.5, wfy - box * 0.5, 0, 0.5, 0.5, -np.pi / 4, np.nanmin(data) - 1],
            # lower bound: [xc, yc, amp, sigx, sigy, rotation,  bg]
            [wfx + box * 0.5, wfy + box * 0.5, 1e7, 20, 20, np.pi / 4, np.nanmax(data[yv, xv]) + 1],  # upper bound
            psf_function=gaussian_psf, method='trf',
            box=box  # only fit a subregion +/- 5 px from centroid
        )
    except:
        log.info(f"WARNING: trouble fitting Gaussian PSF to star at {wfx},{wfy}")
        log.info("  check location of comparison star in the first few images")
        log.info("  fitting parameters are out of bounds")
        log.info(f"  init: {init}")
        log.info(f" lower: {[wfx - 5, wfy - 5, 0, 0, 0, -np.pi / 4, np.nanmin(data) - 1]}")
        log.info(f" upper: {[wfx + 5, wfy + 5, 1e7, 20, 20, np.pi / 4, np.nanmax(data[yv, xv]) + 1]}")

        # use LM in unbounded optimization
        pars = fit_psf(
            data, [wfx, wfy], init,
            [wfx - 5, wfy - 5, 0, 0, 0, -PI / 4, np.nanmin(data) - 1],
            [wfx + 5, wfy + 5, 1e7, 20, 20, PI / 4, np.nanmax(data[yv, xv]) + 1],
            psf_function=gaussian_psf,
            box=box, method='lm'
        )

    if pars[2] <= 10:
        log.info(
            f"CAUTION: Measured flux amplitude is really low---are you sure there is a star at {np.round(pos, 2)}?")

    return pars


# Method calculates the flux of the star (uses the skybg_phot method to do backgorund sub)
def aperPhot(data, xc, yc, r=5, dr=5):
    if dr > 0:
        bgflux, sigmabg, Nbg = skybg_phot(data, xc, yc, r + 2, dr)
    else:
        bgflux, sigmabg, Nbg = 0, 0
    positions = [(xc, yc)]
    bdata = data - bgflux

    apertures = CircularAperture(positions, r=r)
    phot_table = aperture_photometry(bdata, apertures, method='exact')

    return float(phot_table['aperture_sum']), bgflux


def skybg_phot(data, xc, yc, r=10, dr=5, ptol=99, debug=False):
    # create a crude annulus to mask out bright background pixels
    xv, yv = mesh_box([xc, yc], np.round(r + dr))
    rv = ((xv - xc) ** 2 + (yv - yc) ** 2) ** 0.5
    mask = (rv > r) & (rv < (r + dr))
    try:
        cutoff = np.nanpercentile(data[yv, xv][mask], ptol)
    except IndexError:
        log.info(f"IndexError, problem computing sky bg for {xc:.1f}, {yc:.1f}. Check if star is present or close to border.")
        cutoff = np.nanpercentile(data[yv, xv], ptol)

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

        ax[1, 0].hist(bgsky.flatten(), label='Sky Annulus ({:.1f}, {:.1f})'.format(np.nanmedian(bgsky), amode),
                      alpha=0.5, bins=np.arange(minb, maxb))
        ax[1, 0].hist(dat.flatten(), label='Clipped ({:.1f}, {:.1f})'.format(np.nanmedian(dat), cmode), alpha=0.5,
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


# make plots of the centroid positions as a function of time
def plotCentroids(xTarg, yTarg, xRef, yRef, times, targetname, date):
    times = np.array(times)
    # X TARGET
    plt.figure()
    plt.plot(times - np.nanmin(times), xTarg, '-bo')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('X Pixel Position')
    plt.title(targetname + ' X Centroid Position ' + date)
    plt.savefig(Path(exotic_infoDict['saveplot']) / "temp" / f"XCentroidPos_{targetname}_{date}.png")
    plt.close()

    # Y TARGET
    plt.figure()
    plt.plot(times - np.nanmin(times), yTarg, '-bo')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('Y Pixel Position')
    plt.title(targetname + ' Y Centroid Position ' + date)
    plt.savefig(Path(exotic_infoDict['saveplot']) / "temp" / f"YCentroidPos_{targetname}_{date}.png")
    plt.close()

    # X COMP
    plt.figure()
    plt.plot(times - np.nanmin(times), xRef, '-ro')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('X Pixel Position')
    plt.title('Comp Star X Centroid Position ' + date)
    plt.savefig(Path(exotic_infoDict['saveplot']) / "temp" / f"CompStarXCentroidPos_{targetname}_{date}.png")
    plt.close()

    # Y COMP
    plt.figure()
    plt.plot(times - np.nanmin(times), yRef, '-ro')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('Y Pixel Position')
    plt.title('Comp Star Y Centroid Position ' + date)
    plt.savefig(Path(exotic_infoDict['saveplot']) / "temp" / f"CompStarYCentroidPos_{targetname}_{date}.png")
    plt.close()

    # X DISTANCE BETWEEN TARGET AND COMP
    plt.figure()
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('X Pixel Distance')
    for e in range(0, len(xTarg)):
        plt.plot(times[e] - np.nanmin(times), abs(int(xTarg[e]) - int(xRef[e])), 'bo')
    plt.title('Distance between Target and Comparison X position')
    plt.savefig(Path(exotic_infoDict['saveplot']) / "temp" / f"XCentroidDistance_{targetname}_{date}.png")
    plt.close()

    # Y DISTANCE BETWEEN TARGET AND COMP
    plt.figure()
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('Y Pixel Difference')

    for d in range(0, len(yTarg)):
        plt.plot(times[d] - np.nanmin(times), abs(int(yTarg[d]) - int(yRef[d])), 'bo')
    plt.title('Difference between Target and Comparison Y position')
    plt.savefig(Path(exotic_infoDict['saveplot']) / "temp" / f"YCentroidDistance_{targetname}_{date}.png")
    plt.close()


def psf_format(data, pos, init=[]):
    target = fit_centroid(data, pos, init=init, box=10)

    return {
        'x': target[0],
        'y': target[1],
        'amp': target[2],
        'sig_x': target[3],
        'sig_y': target[4],
        'rot': target[5],
        'off': target[6]
    }


def realTimeReduce(i, target_name, ax, distFC, real_time_imgs, UIprevTPX, UIprevTPY, UIprevRPX, UIprevRPY):
    targetFluxVals = []
    referenceFluxVals = []
    normalizedFluxVals = []
    fileNameList = []
    timeList = []

    for file_name in real_time_imgs:
        extension = 0
        header = fits.getheader(filename=file_name, ext=extension)
        while header['NAXIS'] == 0:
            extension += 1
            header = fits.getheader(filename=file_name, ext=extension)
        timeVal = getJulianTime(header)
        timeList.append(timeVal)
        fileNameList.append(file_name)

    # Time sorts the file names based on the fits file header
    timeSortedNames = [x for _, x in sorted(zip(timeList, fileNameList))]

    # Time sorts file
    timeList = np.array(timeList)
    timeList = timeList[np.argsort(np.array(timeList))]

    # Extracts data from the image file and puts it in a 2D numpy array: firstImageData
    firstImageData = fits.getdata(timeSortedNames[0], ext=0)
    firstImageHeader = fits.getheader(timeSortedNames[0], ext=0)

    # fit first image
    targ = psf_format(firstImageData, [UIprevTPX, UIprevTPY])
    ref = psf_format(firstImageData, [UIprevRPX, UIprevRPY])

    # just use one aperture and annulus
    apertureR = 3 * max(targ['sig_x'], targ['sig_y'])
    annulusR = 10

    for i, imageFile in enumerate(timeSortedNames):

        hdul = fits.open(name=imageFile, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)

        # Extracts data from the image file and puts it in a 2D numpy array: imageData
        imageData = fits.getdata(imageFile, ext=0)

        # Find the target star in the image and get its pixel coordinates if it is the first file
        if i == 0:
            # Initializing the star location guess as the user inputted pixel coordinates
            prevTPX, prevTPY, prevRPX, prevRPY = UIprevTPX, UIprevTPY, UIprevRPX, UIprevRPY
            prevTSigX, prevTSigY, prevRSigX, prevRSigY = targ['sig_x'], targ['sig_y'], ref['sig_x'], ref['sig_y']

            prevImageData = imageData  # no shift should be registered

        # ---FLUX CALCULATION WITH BACKGROUND SUBTRACTION---------------------------------

        # corrects for any image shifts that result from a tracking slip
        shift, error, diffphase = phase_cross_correlation(prevImageData, imageData)
        xShift = shift[1]
        yShift = shift[0]

        prevTPX = prevTPX - xShift
        prevTPY = prevTPY - yShift
        prevRPX = prevRPX - xShift
        prevRPY = prevRPY - yShift

        # --------GAUSSIAN FIT AND CENTROIDING----------------------------------------------

        txmin = int(prevTPX) - distFC  # left
        txmax = int(prevTPX) + distFC  # right
        tymin = int(prevTPY) - distFC  # top
        tymax = int(prevTPY) + distFC  # bottom

        targSearchA = imageData[tymin:tymax, txmin:txmax]

        # Set reference search area
        rxmin = int(prevRPX) - distFC  # left
        rxmax = int(prevRPX) + distFC  # right
        rymin = int(prevRPY) - distFC  # top
        rymax = int(prevRPY) + distFC  # bottom

        refSearchA = imageData[rymin:rymax, rxmin:rxmax]

        # Guess at Gaussian Parameters and feed them in to help gaussian fitter

        tGuessAmp = targSearchA.max() - targSearchA.min()

        # Fits Centroid for Target
        myPriors = [tGuessAmp, prevTSigX, prevTSigY, 0, targSearchA.min()]
        targ = psf_format(imageData, [prevTPX, prevTPY], init=myPriors)
        tpsfFlux = 2 * PI * targ['amp'] * targ['sig_x'] * targ['sig_y']
        currTPX = targ['x']
        currTPY = targ['y']

        # Fits Centroid for Reference
        rGuessAmp = refSearchA.max() - refSearchA.min()
        myRefPriors = [rGuessAmp, prevRSigX, prevRSigY, 0, refSearchA.min()]
        ref = psf_format(imageData, [prevRPX, prevRPY], init=myRefPriors)
        rpsfFlux = 2 * PI * ref['amp'] * ref['sig_x'] * ref['sig_y']
        currRPX = ref['x']
        currRPY = ref['y']

        # gets the flux value of the target star and
        tFluxVal, tTotCts = aperPhot(imageData, currTPX, currTPY, apertureR, annulusR)
        targetFluxVals.append(tFluxVal)  # adds tFluxVal to the total list of flux values of target star

        # gets the flux value of the reference star and subtracts the background light
        rFluxVal, rTotCts = aperPhot(imageData, currRPX, currRPY, apertureR, annulusR)
        referenceFluxVals.append(rFluxVal)  # adds rFluxVal to the total list of flux values of reference star

        normalizedFluxVals.append((tFluxVal / rFluxVal))

        # UPDATE PIXEL COORDINATES and SIGMAS
        # target
        prevTPX = currTPX
        prevTPY = currTPY
        prevTSigX = targ['sig_x']
        prevTSigY = targ['sig_y']
        # reference
        prevRPX = currRPX
        prevRPY = currRPY
        prevRSigX = ref['sig_x']
        prevRSigY = ref['sig_y']

        # UPDATE FILE COUNT
        prevImageData = imageData

        hdul.close()
        del hdul

    ax.clear()
    ax.set_title(target_name)
    ax.set_ylabel('Normalized Flux')
    ax.set_xlabel('Time (jd)')
    ax.plot(timeList, normalizedFluxVals, 'bo')


def fit_lightcurve(times, tFlux, cFlux, airmass, ld, pDict):
    # remove outliers
    si = np.argsort(times)
    dt = np.mean(np.diff(np.sort(times)))
    ndt = int(25. / 24. / 60. / dt) * 2 + 1
    filtered_data = sigma_clip((tFlux / cFlux)[si], sigma=3, dt=ndt)
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
        'omega': 0,  # Arg of periastron
        'tmid': pDict['midT'],  # time of mid transit [day]
        'a1': arrayFinalFlux.mean(),  # max() - arrayFinalFlux.min(), #mid Flux
        'a2': 0,  # Flux lower bound
    }

    arrayPhases = (arrayTimes - pDict['midT']) / prior['per']
    prior['tmid'] = pDict['midT'] + np.floor(arrayPhases).max() * prior['per']
    upper = prior['tmid'] + np.abs(25 * pDict['midTUnc'] + np.floor(arrayPhases).max() * 25 * pDict['pPerUnc'])
    lower = prior['tmid'] - np.abs(25 * pDict['midTUnc'] + np.floor(arrayPhases).max() * 25 * pDict['pPerUnc'])

    if np.floor(arrayPhases).max() - np.floor(arrayPhases).min() == 0:
        log.info("\nWARNING!")
        log.info(" Estimated mid-transit time is not within the observations")
        log.info(" Check Period & Mid-transit time in inits.json. Make sure the uncertainties are not 0 or Nan.")
        log.info(f"  obs start:{arrayTimes.min()}")
        log.info(f"    obs end:{arrayTimes.max()}")
        log.info(f" tmid prior:{prior['tmid']}\n")

    mybounds = {
        'rprs': [0, pDict['rprs'] * 1.25],
        'tmid': [lower, upper],
        'ars': [pDict['aRs'] - 5 * pDict['aRsUnc'], pDict['aRs'] + 5 * pDict['aRsUnc']],

        'a1': [0.5 * min(arrayFinalFlux), 2 * max(arrayFinalFlux)],
        'a2': [-1, 1]
    }

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
                        default='', type=str,
                        help="Plots transit in real-time while observing with a telescope. "
                             "An initialization file (e.g., inits.json) is required to use this command.")
    parser.add_argument('-red', '--reduce',
                        default='', type=str,
                        help="Performs aperture photometry on FITS files and a reduction on dataset. "
                             "An initialization file (e.g., inits.json) is required to use this command.")
    parser.add_argument('-pre', '--prereduced',
                        default='', type=str,
                        help="Performs a reduction on dataset using the nested sampler only. "
                             "An initialization file (e.g., inits.json) is required to use this command.")
    parser.add_argument('-ov', '--override',
                        action='store_true',
                        help="Adopts all JSON planetary parameters, which will override the NASA Exoplanet Archive. "
                             "Can be used as an additional argument with -red, --reduce and -pre, --prereduced. "
                             "Do not combine with the -nea, --nasaexoarch argument.")
    parser.add_argument('-nea', '--nasaexoarch',
                        action='store_true',
                        help="Adopts all the NASA Exoplanet Archive planetary parameters from "
                             "https://exoplanetarchive.ipac.caltech.edu. Can be used as an additional argument with "
                             "-red, --reduce and -pre, --prereduced. "
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

    log.info("\n*************************************************************")
    log.info("Welcome to the EXOplanet Transit Interpretation Code (EXOTIC)")
    log.info(f"Version {__version__}")
    log.info("*************************************************************\n")

    # ---INITIALIZATION-------------------------------------------------------
    global exotic_infoDict

    exotic_UIprevTPX, exotic_UIprevTPY, exotic_UIprevRPX, exotic_UIprevRPY = 0, 0, 0, 0

    fileNameList, timeSortedNames, xTargCent, yTargCent, xRefCent, yRefCent, finXTargCent, finYTargCent, finXRefCent, finYRefCent = (
        [] for m in range(10))

    minSTD = 100000  # sets the initial minimum standard deviation absurdly high so it can be replaced immediately
    minChi2 = 100000
    distFC = 10  # gaussian search area

    # ---USER INPUTS--------------------------------------------------------------------------
    if args.realtime:
        reduction_opt = 1
    elif args.reduce or args.prereduced:
        reduction_opt = 2
    else:
        reduction_opt = user_input("\nEnter '1' for Real Time Reduction or '2' for for Complete Reduction: ",
                                   type_=int, val1=1, val2=2)

    #############################
    # Real Time Reduction Routine
    #############################

    if reduction_opt == 1:
        log.info("\n**************************************************************")
        log.info("Real Time Reduction ('Control + C'  or close the plot to quit)")
        log.info("**************************************************************\n")

        if not args.realtime:
            real_time_dir = user_input("Enter the Directory Path of imaging files: ", type_=str)
        else:
            real_time = Path(args.realtime)
            with real_time.open('r') as json_file:
                info = json.load(json_file)
                real_time_dir = info['user_info']['Directory with FITS files']
                targetName = info['planetary_parameters']['Planet Name']
                exotic_UIprevTPX = info['user_info']['Target Star X & Y Pixel'][0]
                exotic_UIprevTPY = info['user_info']['Target Star X & Y Pixel'][1]
                exotic_UIprevRPX = info['user_info']['Comparison Star(s) X & Y Pixel'][0][0]
                exotic_UIprevRPY = info['user_info']['Comparison Star(s) X & Y Pixel'][0][1]

        real_time_dir, real_time_imgs = check_imaging_files(real_time_dir, 'imaging')

        targetName = planet_name(targetName)

        while True:
            carry_on = user_input(f"Type continue after the first image has been taken and saved: ", type_=str)
            if carry_on != 'continue':
                continue
            break

        if not args.realtime:
            exotic_UIprevTPX = user_input(f"{targetName} X Pixel Coordinate: ", type_=int)
            exotic_UIprevTPY = user_input(f"{targetName} Y Pixel Coordinate: ", type_=int)
            exotic_UIprevRPX = user_input("Comp Star X Pixel Coordinate: ", type_=int)
            exotic_UIprevRPY = user_input("Comp Star Y Pixel Coordinate: ", type_=int)

        log.info("Real Time Plotting ('Control + C' or close the plot to quit)")
        log.info("\nPlease be patient. It will take at least 15 seconds for the first image to get plotted.")

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(targetName)
        ax.set_ylabel('Normalized Flux')
        ax.set_xlabel('Time (jd)')

        anim = FuncAnimation(fig, realTimeReduce,
                             fargs=(targetName, ax, distFC, real_time_imgs, exotic_UIprevTPX, exotic_UIprevTPY,
                                    exotic_UIprevRPX, exotic_UIprevRPY), interval=15000)  # refresh every 15 seconds
        plt.show()

    ###########################
    # Complete Reduction Routine
    ###########################

    # ----USER INPUTS----------------------------------------------------------
    else:
        log.info("\n**************************")
        log.info("Complete Reduction Routine")
        log.info("**************************")

        init_path, wcs_file, wcs_header = None, None, None
        compStarList = []

        exotic_infoDict = {'fitsdir': None, 'saveplot': None, 'flatsdir': None, 'darksdir': None, 'biasesdir': None,
                           'aavsonum': None, 'secondobs': None, 'date': None, 'lat': None, 'long': None, 'elev': None,
                           'ctype': None, 'pixelbin': None, 'filter': None, 'wl_min': None, 'wl_max': None,
                           'notes': None,
                           'tarcoords': None, 'compstars': None, 'plate_opt': None, 'pixel_scale': None, 
                           'image_align':None }

        userpDict = {'ra': None, 'dec': None, 'pName': None, 'sName': None, 'pPer': None, 'pPerUnc': None,
                     'midT': None, 'midTUnc': None, 'rprs': None, 'rprsUnc': None, 'aRs': None, 'aRsUnc': None,
                     'inc': None, 'incUnc': None, 'ecc': None, 'teff': None,
                     'teffUncPos': None, 'teffUncNeg': None, 'met': None, 'metUncPos': None, 'metUncNeg': None,
                     'logg': None, 'loggUncPos': None, 'loggUncNeg': None}

        if args.reduce:
            fitsortext = 1
            init_path = args.reduce
        elif args.prereduced:
            fitsortext = 2
            init_path = args.prereduced
        else:
            fitsortext = user_input("Enter '1' to perform aperture photometry on fits files or '2' to start with "
                                    "pre-reduced data in a .txt format: ", type_=int, val1=1, val2=2)

        if not args.reduce and not args.prereduced:
            fileorcommandline = user_input("\nHow would you like to input your initial parameters? "
                                           "Enter '1' to use the Command Line or '2' to use an input file: ",
                                           type_=int, val1=1, val2=2)
        else:
            fileorcommandline = 2

        # Read in input file rather than using the command line
        if fileorcommandline == 2:
            exotic_infoDict, userpDict = get_initialization_file(exotic_infoDict, userpDict, init_path)
            init_obj = InitializationFile(exotic_infoDict, userpDict['pName'])
            exotic_infoDict, userpDict['pName'] = init_obj.get_info()

            if exotic_infoDict['flatsdir'] is None:
                flatsBool = False
            else:
                flatsBool = True

            if exotic_infoDict['darksdir'] is None:
                darksBool = False
            else:
                darksBool = True

            if exotic_infoDict['biasesdir'] is None:
                biasesBool = False
            else:
                biasesBool = True

            if flatsBool + darksBool + biasesBool:
                cals = 'y'
            else:
                cals = 'n'

            # Initial position of target star
            exotic_UIprevTPX = exotic_infoDict['tarcoords'][0]
            exotic_UIprevTPY = exotic_infoDict['tarcoords'][1]

            # Read in locations of comp stars
            for i in exotic_infoDict['compstars']:
                compStarList.append((i[0], i[1]))

        if fitsortext == 1:
            # File directory name and initial guess at target and comp star locations on image.
            if fileorcommandline == 1:
                exotic_infoDict['fitsdir'] = user_input("\nEnter the directory path where imaging files are located. "
                                                        "(Example using the sample data: "
                                                        "sample-data/HatP32Dec202017): ", type_=str)

            exotic_infoDict['fitsdir'], inputfiles = check_imaging_files(exotic_infoDict['fitsdir'], 'imaging')
        else:
            data_file = prereduced_file()
            exotic_infoDict['exposure'] = user_input("Please enter your image exposure time (seconds): ", type_=int)

        exotic_infoDict['saveplot'] = save_directory(exotic_infoDict['saveplot'])

        # Make a temp directory of helpful files
        Path(Path(exotic_infoDict['saveplot']) / "temp").mkdir(exist_ok=True)

        userpDict['pName'] = planet_name(userpDict['pName'])

        if not args.override:
            nea_obj = NASAExoplanetArchive(planet=userpDict['pName'])
            userpDict['pName'], CandidatePlanetBool, pDict = nea_obj.planet_info()
        else:
            pDict = userpDict
            CandidatePlanetBool = False

        exotic_infoDict['date'] = obs_date(exotic_infoDict['date'])

        if fitsortext == 1:
            exotic_infoDict['lat'] = latitude(exotic_infoDict['lat'])
            exotic_infoDict['long'] = longitude(exotic_infoDict['long'])
            exotic_infoDict['elev'] = elevation(exotic_infoDict['elev'], exotic_infoDict['lat'], exotic_infoDict['long'])

            # TARGET STAR
            if fileorcommandline == 1:
                exotic_UIprevTPX = user_input(f"\n{userpDict['pName']} X Pixel Coordinate: ", type_=int)
                exotic_UIprevTPY = user_input(f"\n{userpDict['pName']} Y Pixel Coordinate: ", type_=int)
                numCompStars = user_input("How many comparison stars would you like to use? (1-10): ", type_=int)

                for num in range(numCompStars):
                    xpix = user_input(f"Comparison Star {num + 1} X Pixel Coordinate: ", type_=int)
                    ypix = user_input(f"Comparison Star {num + 1} Y Pixel Coordinate: ", type_=int)
                    compStarList.append((xpix, ypix))

            # ---HANDLE CALIBRATION IMAGES------------------------------------------------
            if fileorcommandline == 1:
                cals = user_input("\nDo you have any calibration images (flats, darks or biases)? (y/n): ",
                                  type_=str, val1='y', val2='n')

            # if they have cals, handle them by calculating the median flat, dark or bias
            if cals == 'y':

                # darks
                if fileorcommandline == 1:
                    darks = user_input("\nDo you have darks? (y/n): ", type_=str, val1='y', val2='n')
                    if darks == 'y':
                        darksBool = True
                        exotic_infoDict['darksdir'] = user_input("Enter the directory path to your darks "
                                                                 "(must be in their own separate folder): ",
                                                                 type_=str)  # +"/*.FITS"
                    else:
                        darksBool = False

                # Only do the dark correction if user selects this option
                if darksBool:
                    exotic_infoDict['darksdir'], inputdarks = check_imaging_files(exotic_infoDict['darksdir'], 'darks')
                    darksImgList = []
                    for darkFile in inputdarks:
                        darkData = fits.getdata(darkFile, ext=0)
                        darksImgList.append(darkData)
                    generalDark = np.median(darksImgList, axis=0)

                # biases
                if fileorcommandline == 1:
                    biases = user_input("\nDo you have biases? (y/n): ", type_=str, val1='y', val2='n')
                    if biases == 'y':
                        biasesBool = True
                        exotic_infoDict['biasesdir'] = user_input("Enter the directory path to your biases "
                                                                  "(must be in their own separate folder): ",
                                                                  type_=str)  # +"/*.FITS"
                    else:
                        biasesBool = False

                if biasesBool:
                    # Add / to end of directory if user does not input it
                    exotic_infoDict['biasesdir'], inputbiases = check_imaging_files(exotic_infoDict['biasesdir'],
                                                                                    'biases')
                    biasesImgList = []
                    for biasFile in inputbiases:
                        biasData = fits.getdata(biasFile, ext=0)
                        biasesImgList.append(biasData)
                    generalBias = np.median(biasesImgList, axis=0)

                # flats
                if fileorcommandline == 1:
                    flats = user_input("\nDo you have flats? (y/n): ", type_=str, val1='y', val2='n')

                    if flats == 'y':
                        flatsBool = True
                        exotic_infoDict['flatsdir'] = user_input("Enter the directory path to your flats "
                                                                 "(must be in their own separate folder): ",
                                                                 type_=str)  # +"/*.FITS"
                    else:
                        flatsBool = False

                if flatsBool:
                    exotic_infoDict['flatsdir'], inputflats = check_imaging_files(exotic_infoDict['flatsdir'], 'flats')
                    flatsImgList = []
                    for flatFile in inputflats:
                        flatData = fits.getdata(flatFile, ext=0)
                        flatsImgList.append(flatData)
                    notNormFlat = np.median(flatsImgList, axis=0)

                    # if the bias exists, bias subtract the flatfield
                    if biasesBool:
                        notNormFlat = notNormFlat - generalBias

                    # NORMALIZE
                    medi = np.median(notNormFlat)
                    generalFlat = notNormFlat / medi
            else:
                flatsBool = False
                darksBool = False
                biasesBool = False

        log.info("***************************************\n")

        # Handle AAVSO Formatting
        if fileorcommandline == 1:
            exotic_infoDict['aavsonum'] = user_input("Please enter your AAVSO Observer Account Number "
                                                     "(type N/A if you do not currently have an account): ", type_=str)
            exotic_infoDict['secondobs'] = user_input("Please enter your comma-separated secondary observer codes "
                                                      "(or type N/A if only 1 observer code): ", type_=str)
            exotic_infoDict['ctype'] = user_input("Please enter your camera type (CCD or DSLR): ", type_=str)
            exotic_infoDict['pixelbin'] = user_input("Please enter your pixel binning: ", type_=str)
            # exotic_infoDict['exposure'] = user_input("Please enter your exposure time (seconds): ", type_=int)
            exotic_infoDict['filter'] = user_input("Please enter your filter name from the options at "
                                                   "http://astroutils.astronomy.ohio-state.edu/exofast/limbdark.shtml: ",
                                                   type_=str)
            exotic_infoDict['notes'] = user_input("Please enter any observing notes (seeing, weather, etc.)."
                                                  "If none, leave blank and press enter: ", type_=str)

        if not exotic_infoDict['notes'].replace(' ', ''):
            exotic_infoDict['notes'] = "na"

        if fileorcommandline == 2:
            if args.nasaexoarch:
                pass
            elif args.override:
                pDict['ra'], pDict['dec'] = radec_hours_to_degree(pDict['ra'], pDict['dec'])
            else:
                diff = False

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
        ld0, ld1, ld2, ld3, exotic_infoDict['filter'], exotic_infoDict['wl_min'], exotic_infoDict[
            'wl_max'] = ld_obj.nonlinear_ld()
        ld = [ld0[0], ld1[0], ld2[0], ld3[0]]

        # check for Nans + Zeros
        for k in pDict:
            if "Unc" in k:
                if not pDict[k]:
                    log.info(f"WARNING! {k} uncertainty is 0. Please use a non-zero value in inits.json")
                    pDict[k] = 1
                elif pDict[k] == 0 or np.isnan(pDict[k]):
                    log.info(f"WARNING! {k} uncertainty is 0. Please use a non-zero value in inits.json")
                    pDict[k] = 1
            elif pDict[k] is None:
                log.info(f"WARNING! {k} is None. Please use a numeric value in inits.json")
                pDict[k] = 0

        if fitsortext == 1:
            log.info("\n**************************"
                     "\nStarting Reduction Process"
                     "\n**************************\n")

            #########################################
            # FLUX DATA EXTRACTION AND MANIPULATION
            #########################################

            allImageData, timeList, airMassList, exptimes = [], [], [], []

            # TODO filter input files to get good reference for image alignment

            inputfiles = corruption_check(inputfiles)

            # time sort images
            times = []
            for file in inputfiles:
                extension = 0
                header = fits.getheader(filename=file, ext=extension)
                while header['NAXIS'] == 0:
                    extension += 1
                    header = fits.getheader(filename=file, ext=extension)
                times.append(getJulianTime(header))

            si = np.argsort(times)
            inputfiles = np.array(inputfiles)[si]

            # fit target in the first image and use it to determine aperture and annulus range
            inc = 0
            for file in inputfiles:
                first_image = fits.getdata(file, ext=0)
                try:
                    args = fit_centroid(first_image, [exotic_UIprevTPX, exotic_UIprevTPY], box=10)
                    break
                except Exception:
                    inc += 1
                finally:
                    del first_image, args

            inputfiles = inputfiles[inc:]
            wcs_file = check_wcs(inputfiles[0], exotic_infoDict['saveplot'], exotic_infoDict['plate_opt'])

            if wcs_file:
                log.info(f"\nHere is the path to your plate solution: {wcs_file}")
                wcs_header = fits.getheader(filename=wcs_file)
                rafile, decfile = get_radec(wcs_header)

                # Checking pixel coordinates against plate solution
                exotic_UIprevTPX, exotic_UIprevTPY = check_targetpixelwcs(exotic_UIprevTPX, exotic_UIprevTPY,
                                                                          pDict['ra'], pDict['dec'], rafile, decfile)

                for comp in compStarList:
                    log.info("\nChecking for variability in Comparison Star: \n"
                             f"Pixel X: {comp[0]} Pixel Y: {comp[1]}")
                    if variableStarCheck(rafile[comp[1]][comp[0]], decfile[comp[1]][comp[0]]):
                        log.info("\nCurrent comparison star is variable, proceeding to next star.")
                        compStarList.remove(comp)

            # alloc psf fitting param
            psf_data = {
                # x-cent, y-cent, amplitude, sigma-x, sigma-y, rotation, offset
                'target': np.zeros((len(inputfiles), 7)),  # PSF fit
                'target_align': np.zeros((len(inputfiles), 2)),  # image alignment estimate
            }

            # aperture sizes in stdev (sigma) of PSF
            apers = np.linspace(2, 6, 10)
            annuli = np.linspace(6, 15, 10)

            aper_data = {
                'target': np.zeros((len(inputfiles), len(apers), len(annuli))),
                'target_err': np.zeros((len(inputfiles), len(apers), len(annuli))),
                'target_bg': np.zeros((len(inputfiles), len(apers), len(annuli)))
            }

            for i, coord in enumerate(compStarList):
                ckey = "comp{}".format(i + 1)
                psf_data[ckey] = np.zeros((len(inputfiles), 7))
                psf_data[ckey + "_align"] = np.zeros((len(inputfiles), 2))
                aper_data[ckey] = np.zeros((len(inputfiles), len(apers), len(annuli)))
                aper_data[ckey + "_bg"] = np.zeros((len(inputfiles), len(apers), len(annuli)))

            alignmentBool = exotic_infoDict.get('image_align', True)
            if alignmentBool == 'y':
                alignmentBool = True
            elif alignmentBool == 'n':
                alignmentBool = False

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
                airMass = getAirMass(hdul, pDict['ra'], pDict['dec'], exotic_infoDict['lat'], exotic_infoDict['long'],
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

                # apply cals if applicable
                if darksBool:
                    if i == 0:
                        log.info("Dark subtracting images.")
                    imageData = imageData - generalDark
                elif biasesBool: # if a dark is not available, then at least subtract off the pedestal via the bias
                    if i == 0:
                        log.info("Bias-correcting images.")
                    imageData = imageData - generalBias
                else:
                    pass

                if flatsBool:
                    if i == 0:
                        log.info("Flattening images.")
                    generalFlat[generalFlat == 0] = 1
                    imageData = imageData / generalFlat

                if i == 0:
                    image_scale = get_pixel_scale(wcs_header, image_header, exotic_infoDict['pixel_scale'])

                    log.info(f"Reference Image for Alignment: {fileName}")
                    firstImage = np.copy(imageData)

                    #log.info("\nAligning your images from FITS files. Please wait.")

                # Image Alignment
                if alignmentBool:
                    apos, arot = transformation(np.array([firstImage, imageData]), len(inputfiles), fileName, i)
                else:
                    apos = np.array([[0,0]])
                    arot = 0

                # Fit PSF for target star
                if (np.pi - 0.1) <= np.abs(arot) <= (np.pi + 0.1):
                    xrot = exotic_UIprevTPX * np.cos(arot) - exotic_UIprevTPY * np.sin(arot) + apos[0][0]
                    yrot = exotic_UIprevTPX * np.sin(arot) + exotic_UIprevTPY * np.cos(arot) + apos[0][1]
                else:
                    xrot = exotic_UIprevTPX * np.cos(arot) - exotic_UIprevTPY * np.sin(arot) - apos[0][0]
                    yrot = exotic_UIprevTPX * np.sin(arot) + exotic_UIprevTPY * np.cos(arot) - apos[0][1]

                psf_data["target_align"][i] = [xrot,yrot]
                if i == 0:
                    psf_data["target"][i] = fit_centroid(imageData, [xrot, yrot], box=10)
                else:
                    if alignmentBool:
                        psf_data["target"][i] = fit_centroid(
                            imageData,
                            [xrot, yrot],
                            psf_data["target"][0][2:],  # reference psf in first image
                            box=10)
                    else:
                        # use previous PSF as prior
                        psf_data["target"][i] = fit_centroid(
                            imageData,
                            psf_data["target"][i-1][:2],
                            psf_data["target"][i-1][2:],  # reference psf in first image
                            box=10)

                        # check for change in amplitude of PSF
                        if np.abs( (psf_data["target"][i][2]-psf_data["target"][i-1][2])/psf_data["target"][i-1][2]) > 0.5:
                            log.info("Can't find target. Trying to align...")

                            apos, arot = transformation(np.array([firstImage, imageData]), len(inputfiles), fileName, i)

                            # Fit PSF for target star
                            if 3.0 <= np.abs(arot) <= 3.3:
                                xrot = exotic_UIprevTPX * np.cos(arot) - exotic_UIprevTPY * np.sin(arot) + apos[0][0]
                                yrot = exotic_UIprevTPX * np.sin(arot) + exotic_UIprevTPY * np.cos(arot) + apos[0][1]
                            else:
                                xrot = exotic_UIprevTPX * np.cos(arot) - exotic_UIprevTPY * np.sin(arot) - apos[0][0]
                                yrot = exotic_UIprevTPX * np.sin(arot) + exotic_UIprevTPY * np.cos(arot) - apos[0][1]

                            psf_data["target"][i] = fit_centroid( imageData, [xrot, yrot], psf_data["target"][0][2:], box=10)

                # fit for the centroids in all images
                for j,coord in enumerate(compStarList):
                    ckey = "comp{}".format(j+1)
                    # apply transformation
                    if (np.pi - 0.1) <= np.abs(arot) <= (np.pi + 0.1):
                        xrot = coord[0] * np.cos(arot) - coord[1] * np.sin(arot) + apos[0][0]
                        yrot = coord[0] * np.sin(arot) + coord[1] * np.cos(arot) + apos[0][1]
                    else:
                        xrot = coord[0] * np.cos(arot) - coord[1] * np.sin(arot) - apos[0][0]
                        yrot = coord[0] * np.sin(arot) + coord[1] * np.cos(arot) - apos[0][1]

                    psf_data[ckey+"_align"][i] = [xrot,yrot]
                    if i == 0:
                        psf_data[ckey][i] = fit_centroid(imageData, [xrot, yrot], box=10)
                    else:
                        if alignmentBool:
                            psf_data[ckey][i] = fit_centroid(
                                imageData,
                                [xrot, yrot],
                                psf_data[ckey][0][2:],  # initialize with psf in first image
                                box=10)
                        else:
                            # use previous PSF as prior
                            psf_data[ckey][i] = fit_centroid(
                                imageData,
                                psf_data[ckey][i-1][:2],
                                psf_data[ckey][i-1][2:],
                                box=10)

                            # check for change in amplitude of PSF
                            if np.abs( (psf_data[ckey][i][2]-psf_data[ckey][i-1][2])/psf_data[ckey][i-1][2]) > 0.5:
                                log.info("Can't find target. Trying to align...")

                                apos, arot = transformation(np.array([firstImage, imageData]), len(inputfiles), fileName, i)

                                # Fit PSF for target star
                                if 3.0 <= np.abs(arot) <= 3.3:
                                    xrot = exotic_UIprevTPX * np.cos(arot) - exotic_UIprevTPY * np.sin(arot) + apos[0][0]
                                    yrot = exotic_UIprevTPX * np.sin(arot) + exotic_UIprevTPY * np.cos(arot) + apos[0][1]
                                else:
                                    xrot = exotic_UIprevTPX * np.cos(arot) - exotic_UIprevTPY * np.sin(arot) - apos[0][0]
                                    yrot = exotic_UIprevTPX * np.sin(arot) + exotic_UIprevTPY * np.cos(arot) - apos[0][1]

                                psf_data[ckey][i] = fit_centroid( imageData, [xrot, yrot], psf_data[ckey][0][2:], box=10)

                # aperture photometry
                if i == 0:
                    sigma = float((psf_data["target"][0][3] + psf_data["target"][0][4]) * 0.5)
                    apers *= sigma
                    annuli *= sigma

                for a, aper in enumerate(apers):
                    for an, annulus in enumerate(annuli):
                        aper_data["target"][i][a][an], aper_data["target_bg"][i][a][an] = aperPhot(imageData,
                                                                                                   psf_data["target"][
                                                                                                       i, 0],
                                                                                                   psf_data["target"][
                                                                                                       i, 1],
                                                                                                   aper, annulus)

                        # loop through comp stars
                        for j, coord in enumerate(compStarList):
                            ckey = "comp{}".format(j + 1)
                            aper_data[ckey][i][a][an], aper_data[ckey + "_bg"][i][a][an] = aperPhot(imageData,
                                                                                                    psf_data[ckey][
                                                                                                        i, 0],
                                                                                                    psf_data[ckey][
                                                                                                        i, 1],
                                                                                                    aper, annulus)

                # close file + delete from memory
                hdul.close()
                del hdul

            # filter bad images
            badmask = (psf_data["target"][:, 0] == 0) | (aper_data["target"][:, 0, 0] == 0) | np.isnan(
                aper_data["target"][:, 0, 0])
            if np.sum(~badmask) == 0:
                log.error("No images to fit...check reference image for alignment (first image of sequence)")

            # convert to numpy arrays
            times = np.array(timeList)[~badmask]
            airmass = np.array(airMassList)[~badmask]
            psf_data["target"] = psf_data["target"][~badmask]
            psf_data["target_align"] = psf_data["target_align"][~badmask]
            si = np.argsort(times)

            # exposure time
            consistent_et = False
            if len(exptimes) > 0:
                consistent_et = all(elem == exptimes[0] for elem in exptimes)

            exptimes = np.array(exptimes)

            if consistent_et:
                exotic_infoDict['exposure'] = exptimes[0]
            else:
                exotic_infoDict['exposure'] = np.median(exptimes)

            # PSF flux
            tFlux = 2 * np.pi * psf_data['target'][:, 2] * psf_data['target'][:, 3] * psf_data['target'][:, 4]

            # loop over comp stars
            for j, coord in enumerate(compStarList):
                ckey = "comp{}".format(j + 1)
                psf_data[ckey] = psf_data[ckey][~badmask]
                psf_data[ckey + "_align"] = psf_data[ckey + "_align"][~badmask]

                cFlux = 2 * np.pi * psf_data[ckey][:, 2] * psf_data[ckey][:, 3] * psf_data[ckey][:, 4]
                myfit = fit_lightcurve(times, tFlux, cFlux, airmass, ld, pDict)
                for k in myfit.bounds.keys():
                    log.debug("  {}: {:.6f}".format(k, myfit.parameters[k]))

                log.debug('The Residual Standard Deviation is: ' + str(
                    round(100 * myfit.residuals.std() / np.median(myfit.data), 6)) + "%")
                log.debug('The Mean Squared Error is: ' + str(round(np.sum(myfit.residuals ** 2), 6)) + '\n')

                resstd = myfit.residuals.std() / np.median(myfit.data)
                if minSTD > resstd:  # If the standard deviation is less than the previous min
                    bestCompStar = j + 1
                    minSTD = resstd
                    minAperture = 0
                    minAnnulus = 15 * sigma
                    arrayNormUnc = myfit.dataerr

                    # sets the lists we want to print to correspond to the optimal aperature
                    goodFluxes = np.copy(myfit.data)
                    goodNormUnc = np.copy(myfit.dataerr)
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

            log.info("Computing best aperture...")

            # Aperture Photometry
            for a, aper in enumerate(apers):

                for an, annulus in enumerate(annuli):
                    tFlux = aper_data['target'][:, a, an]
                    tFlux_err = aper_data['target_err'][:, a, an]

                    # fit without a comparison star
                    myfit = fit_lightcurve(times, tFlux, np.ones(tFlux.shape), airmass, ld, pDict)

                    for k in myfit.bounds.keys():
                        log.debug("  {}: {:.6f}".format(k, myfit.parameters[k]))

                    log.debug('The Residual Standard Deviation is: ' + str(
                        round(100 * myfit.residuals.std() / np.median(myfit.data), 6)) + "%")
                    log.debug('The Mean Squared Error is: ' + str(round(np.sum(myfit.residuals ** 2), 6)) + '\n')

                    resstd = myfit.residuals.std() / np.median(myfit.data)
                    if minSTD > resstd:  # If the standard deviation is less than the previous min
                        minSTD = resstd
                        minAperture = -aper
                        minAnnulus = annulus
                        arrayNormUnc = arrayNormUnc

                        # sets the lists we want to print to correspond to the optimal aperature
                        goodFluxes = np.copy(myfit.data)
                        goodNormUnc = np.copy(myfit.dataerr)
                        nonBJDTimes = np.copy(myfit.time)
                        nonBJDPhases = np.copy(myfit.phase)
                        goodAirmasses = np.copy(myfit.airmass)
                        goodTargets = tFlux
                        goodReferences = cFlux
                        goodTUnc = tFlux ** 0.5
                        goodRUnc = cFlux ** 0.5
                        goodResids = myfit.residuals
                        bestlmfit = myfit

                        finXTargCent = psf_data["target"][:, 0]
                        finYTargCent = psf_data["target"][:, 1]
                        finXRefCent = psf_data[ckey][:, 0]
                        finYRefCent = psf_data[ckey][:, 1]

                    # try to fit data with comp star
                    for j, coord in enumerate(compStarList):
                        ckey = "comp{}".format(j + 1)
                        cFlux = aper_data[ckey][:, a, an]

                        myfit = fit_lightcurve(times, tFlux, cFlux, airmass, ld, pDict)

                        for k in myfit.bounds.keys():
                            log.debug("  {}: {:.6f}".format(k, myfit.parameters[k]))

                        log.debug('The Residual Standard Deviation is: ' + str(
                            round(100 * myfit.residuals.std() / np.median(myfit.data), 6)) + "%")
                        log.debug('The Mean Squared Error is: ' + str(round(np.sum(myfit.residuals ** 2), 6)) + '\n')

                        resstd = myfit.residuals.std() / np.median(myfit.data)
                        if minSTD > resstd:  # If the standard deviation is less than the previous min
                            bestCompStar = j + 1
                            comp_coords = coord
                            minSTD = resstd
                            minAperture = aper
                            minAnnulus = annulus
                            arrayNormUnc = arrayNormUnc

                            # sets the lists we want to print to correspond to the optimal aperature
                            goodFluxes = np.copy(myfit.data)
                            goodNormUnc = np.copy(myfit.dataerr)
                            nonBJDTimes = np.copy(myfit.time)
                            nonBJDPhases = np.copy(myfit.phase)
                            goodAirmasses = np.copy(myfit.airmass)
                            goodTargets = tFlux
                            goodReferences = cFlux
                            goodTUnc = tFlux ** 0.5
                            goodRUnc = cFlux ** 0.5
                            goodResids = myfit.residuals
                            bestlmfit = myfit

                            finXTargCent = psf_data["target"][:, 0]
                            finYTargCent = psf_data["target"][:, 1]
                            finXRefCent = psf_data[ckey][:, 0]
                            finYRefCent = psf_data[ckey][:, 1]

            # log best fit
            if minAperture == 0:  # psf
                log.info('*********************************************')
                log.info('Best Comparison Star: #' + str(bestCompStar))
                log.info('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                log.info('Optimal Method: PSF photometry')
                log.info('********************************************\n')

            elif minAperture < 0:  # no comp star
                log.info('*********************************************')
                log.info('Best Comparison Star: None')
                log.info('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                log.info('Optimal Aperture: ' + str(abs(np.round(minAperture, 2))))
                log.info('Optimal Annulus: ' + str(np.round(minAnnulus, 2)))
                log.info('********************************************\n')
                bestCompStar = None

            else:
                log.info('*********************************************')
                log.info('Best Comparison Star: #' + str(bestCompStar))
                log.info('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                log.info('Optimal Aperture: ' + str(np.round(minAperture, 2)))
                log.info('Optimal Annulus: ' + str(np.round(minAnnulus, 2)))
                log.info('********************************************\n')

            # Take the BJD times from the image headers
            if "BJD_TDB" in image_header or "BJD" in image_header or "BJD_TBD" in image_header:
                goodTimes = nonBJDTimes
            # If not in there, then convert all the final times into BJD - using astropy alone
            else:
                log.info("No BJDs in Image Headers. Converting all JDs to BJD_TDBs.")
                log.info("Please be patient- this step can take a few minutes.")

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
            if len(OOT) == 0: # if user does not get any out-of-transit data, normalize by the max data instead
                OOT = (bestlmfit.transit == np.nanmax(bestlmfit.transit))
            OOTscatter = np.std((bestlmfit.data / bestlmfit.airmass_model)[OOT])  # calculate the scatter in the data
            goodNormUnc = OOTscatter * bestlmfit.airmass_model  # scale this scatter back up by the airmass model and then adopt these as the uncertainties

            # Normalize by OOT per AAVSO Database upload requirements
            goodNormUnc = goodNormUnc / np.nanmedian(goodFluxes[OOT])
            goodFluxes = goodFluxes / np.nanmedian(goodFluxes[OOT])

            goodTimes = goodTimes[si][gi]
            goodFluxes = goodFluxes[si][gi]
            goodNormUnc = goodNormUnc[si][gi]
            goodAirmasses = goodAirmasses[si][gi]

            picframe = 10 * (minAperture + 15 * sigma)
            pltx = [max([0, min([finXTargCent[0], finXRefCent[0]]) - picframe]),
                    min([np.shape(firstImage)[0], max([finXTargCent[0], finXRefCent[0]]) + picframe])]
            plty = [max([0, min([finYTargCent[0], finYRefCent[0]]) - picframe]),
                    min([np.shape(firstImage)[1], max([finYTargCent[0], finYRefCent[0]]) + picframe])]
            plt.close()

            for stretch in [LinearStretch(), SquaredStretch(), SqrtStretch(), LogStretch()]:
                fig, ax = plt.subplots()
                # Draw apertures and sky annuli
                target_circle = plt.Circle((finXTargCent[0], finYTargCent[0]), minAperture, color='lime', fill=False,
                                           ls='-', label='Target')
                target_circle_sky = plt.Circle((finXTargCent[0], finYTargCent[0]), minAperture + minAnnulus,
                                               color='lime', fill=False, ls='--', lw=.5)
                if minAperture >= 0:
                    ref_circle = plt.Circle((finXRefCent[0], finYRefCent[0]), minAperture, color='r', fill=False,
                                            ls='-.', label='Comp')
                    ref_circle_sky = plt.Circle((finXRefCent[0], finYRefCent[0]), minAperture + minAnnulus, color='r',
                                                fill=False, ls='--', lw=.5)

                med_img = median_filter(firstImage, (4, 4))[int(pltx[0]):round(int(pltx[1])),
                          int(plty[0]):round(int(plty[1]))]
                norm = ImageNormalize(firstImage, interval=ZScaleInterval(), stretch=stretch)
                plt.imshow(firstImage, norm=norm, origin='lower', cmap='Greys_r', interpolation=None,
                           vmin=np.percentile(med_img, 5), vmax=np.percentile(med_img, 99))
                plt.plot(finXTargCent[0], finYTargCent[0], marker='+', color='lime')
                ax.add_artist(target_circle)
                ax.add_artist(target_circle_sky)
                if minAperture >= 0:
                    ax.add_artist(ref_circle)
                    ax.add_artist(ref_circle_sky)
                    plt.plot(finXRefCent[0], finYRefCent[0], '+r')
                plt.xlabel("x-axis [pixel]")
                plt.ylabel("y-axis [pixel]")
                plt.title(f"FOV for {pDict['pName']}\n({image_scale} arcsec/pix)")
                plt.xlim(pltx[0], pltx[1])
                plt.ylim(plty[0], plty[1])
                ax.grid(False)
                plt.plot(0, 0, color='lime', ls='-', label='Target')
                if minAperture >= 0:
                    plt.plot(0, 0, color='r', ls='-.', label='Comp')
                l = plt.legend(framealpha=0.75)
                for text in l.get_texts():
                    text.set_color("white")
                apos = '\''
                Path(exotic_infoDict['saveplot']).mkdir(parents=True, exist_ok=True)
                plt.savefig(Path(exotic_infoDict['saveplot']) /
                            f"FOV_{pDict['pName']}_{exotic_infoDict['date']}_"
                            f"{str(stretch.__class__).split('.')[-1].split(apos)[0]}.pdf", bbox_inches='tight')
                plt.savefig(Path(exotic_infoDict['saveplot']) /
                            f"FOV_{pDict['pName']}_{exotic_infoDict['date']}_"
                            f"{str(stretch.__class__).split('.')[-1].split(apos)[0]}.png", bbox_inches='tight')
                plt.close()

            log.info(f"\nFOV file saved as: {exotic_infoDict['saveplot']}/FOV_{pDict['pName']}_"
                     f"{exotic_infoDict['date']}_{str(stretch.__class__).split('.')[-1].split(apos)[0]}.pdf")

            # Centroid position plots
            plotCentroids(finXTargCent[si][gi], finYTargCent[si][gi], finXRefCent[si][gi], finYRefCent[si][gi],
                          goodTimes, pDict['pName'], exotic_infoDict['date'])

            # TODO: convert the exoplanet archive mid transit time to bjd - need to take into account observatory location listed in Exoplanet Archive
            # tMidtoC = astropy.time.Time(timeMidTransit, format='jd', scale='utc')
            # forPhaseResult = utc_tdb.JDUTC_to_BJDTDB(tMidtoC, ra=raDeg, dec=decDeg, lat=lati, longi=longit, alt=2000)
            # bjdMidTOld = float(forPhaseResult[0])
            bjdMidTOld = pDict['midT']

            goodPhasesList = []
            # convert all the phases based on the updated bjd times
            for convertedTime in goodTimes:
                bjdPhase = getPhase(float(convertedTime), pDict['pPer'], bjdMidTOld)
                goodPhasesList.append(bjdPhase)
            goodPhases = np.array(goodPhasesList)

            # Calculate the standard deviation of the normalized flux values
            standardDev1 = np.std(goodFluxes)

            ######################################
            # PLOTS ROUND 1
            ####################################
            # Make plots of raw target and reference values
            plt.figure()
            plt.errorbar(goodTimes, goodTargets[si][gi], yerr=goodTUnc[si][gi], linestyle='None', fmt='-o')
            plt.xlabel('Time (BJD)')
            plt.ylabel('Total Flux')
            # plt.rc('grid', linestyle="-", color='black')
            # plt.grid(True)
            plt.title(pDict['pName'] + ' Raw Flux Values ' + exotic_infoDict['date'])
            plt.savefig(Path(exotic_infoDict['saveplot']) / "temp" /
                        f"TargetRawFlux_{pDict['pName']}_{exotic_infoDict['date']}.png")
            plt.close()

            plt.figure()
            plt.errorbar(goodTimes, goodReferences[si][gi], yerr=goodRUnc[si][gi], linestyle='None', fmt='-o')
            plt.xlabel('Time (BJD)')
            plt.ylabel('Total Flux')
            # plt.rc('grid', linestyle="-", color='black')
            # plt.grid(True)
            plt.title('Comparison Star Raw Flux Values ' + exotic_infoDict['date'])
            plt.savefig(Path(exotic_infoDict['saveplot']) / "temp" /
                        f"CompRawFlux_{pDict['pName']}_{exotic_infoDict['date']}.png")
            plt.close()

            # Plots final reduced light curve (after the 3 sigma clip)
            plt.figure()
            plt.errorbar(goodTimes, goodFluxes, yerr=goodNormUnc, linestyle='None', fmt='-bo')
            plt.xlabel('Time (BJD)')
            plt.ylabel('Normalized Flux')
            # plt.rc('grid', linestyle="-", color='black')
            # plt.grid(True)
            plt.title(pDict['pName'] + ' Normalized Flux vs. Time ' + exotic_infoDict['date'])
            plt.savefig(Path(exotic_infoDict['saveplot']) / "temp" /
                        f"NormalizedFluxTime_{pDict['pName']}_{exotic_infoDict['date']}.png")
            plt.close()

            # Save normalized flux to text file prior to MCMC
            params_file = Path(
                exotic_infoDict['saveplot']) / f"NormalizedFlux_{pDict['pName']}_{exotic_infoDict['date']}.txt"
            with params_file.open('w') as f:
                f.write("BJD,Norm Flux,Norm Err,AM\n")

                for ti, fi, erri, ami in zip(goodTimes, goodFluxes, goodNormUnc, goodAirmasses):
                    f.write(f"{round(ti, 8)},{round(fi, 7)},{round(erri, 6)},{round(ami, 2)}\n")

            log.info("\nOutput File Saved")
        else:
            goodTimes, goodFluxes, goodNormUnc, goodAirmasses = [], [], [], []
            bestCompStar = None

            with data_file.open('r') as f:
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

            # Ask user for time format and convert it if not in BJD_TDB
            log.info("NOTE: If your file is not in one of the following formats, "
                     "please re-reduce your data into one of the time formats recognized by EXOTIC.")

            while True:
                time_format = user_input("Which of the following time formats is your data file stored in? "
                                         "\nBJD_TDB / JD_UTC / MJD_UTC: ", type_=str)
                time_format = time_format.upper().strip()

                if time_format not in ['BJD_TDB', 'JD_UTC', 'MJD_UTC']:
                    log.info("Invalid entry; please try again.")
                else:
                    break

            if time_format != 'BJD_TDB':
                goodTimes = timeConvert(goodTimes, time_format, pDict, exotic_infoDict)

            # Ask user for flux units and convert to flux if in magnitude/millimagnitude
            log.info("NOTE: If your file is not in one of the following formats, "
                     "please re-reduce your data into one of the time formats recognized by EXOTIC.")

            while True:
                flux_format = user_input("Which of the following units of flux is your data file stored in? "
                                         "\nflux / magnitude / millimagnitude: ", type_=str)
                flux_format = flux_format.lower().strip()

                if flux_format not in ['flux', 'magnitude', 'millimagnitude']:
                    log.info("Invalid entry; please try again.")
                else:
                    break

            if flux_format != 'flux':
                goodFluxes, goodNormUnc = fluxConvert(goodFluxes, goodNormUnc, flux_format)

        # for k in myfit.bounds.keys():
        #     print("{:.6f} +- {}".format( myfit.parameters[k], myfit.errors[k]))

        log.info("\n")
        log.info("****************************************")
        log.info("Fitting a Light Curve Model to Your Data")
        log.info("****************************************\n")

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
            'omega': 0,  # Arg of periastron
            'tmid': pDict['midT'],  # time of mid transit [day]
            'a1': goodFluxes.mean(),  # max() - arrayFinalFlux.min(), #mid Flux
            'a2': 0,  # Flux lower bound
        }

        phase = (goodTimes - prior['tmid']) / prior['per']
        prior['tmid'] = pDict['midT'] + np.floor(phase).max() * prior['per']
        upper = pDict['midT'] + 35 * pDict['midTUnc'] + np.floor(phase).max() * (pDict['pPer'] + 35 * pDict['pPerUnc'])
        lower = pDict['midT'] - 35 * pDict['midTUnc'] + np.floor(phase).max() * (pDict['pPer'] - 35 * pDict['pPerUnc'])

        if np.floor(phase).max() - np.floor(phase).min() == 0:
            log.info("ERROR: Estimated mid-transit not in observation range (check priors or observation time)")
            log.info(f"start:{goodTimes.min()}")
            log.info(f"  end:{goodTimes.max()}")
            log.info(f"prior:{prior['tmid']}")

        mybounds = {
            'rprs': [0, pDict['rprs'] * 1.25],
            'tmid': [lower, upper],
            'ars': [pDict['aRs'] - 1, pDict['aRs'] + 1],
            'a1': [min(0, np.nanmin(goodFluxes)), 3 * np.nanmax(goodFluxes)],
            'a2': [-3, 3],
        }

        # final light curve fit
        myfit = lc_fitter(goodTimes, goodFluxes, goodNormUnc, goodAirmasses, prior, mybounds, mode='ns')
        # myfit.dataerr *= np.sqrt(myfit.chi2 / myfit.data.shape[0])  # scale errorbars by sqrt(rchi2)
        # myfit.detrendederr *= np.sqrt(myfit.chi2 / myfit.data.shape[0])

        # estimate transit duration
        pars = dict(**myfit.parameters)
        times = np.linspace(np.min(myfit.time), np.max(myfit.time), 1000)
        dt = np.diff(times).mean()
        durs = []
        for r in range(1000):
            # randomize parameters
            for k in myfit.errors:
                pars[k] = np.random.normal(myfit.parameters[k], myfit.errors[k])

            data = transit(times, pars)
            tmask = data < 1
            durs.append(tmask.sum()*dt)

        ########################
        # PLOT FINAL LIGHT CURVE
        ########################
        f, (ax_lc, ax_res) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

        ax_lc.set_title(pDict['pName'])
        ax_res.set_xlabel('Phase')

        # clip plot to get rid of white space
        ax_res.set_xlim([min(myfit.phase), max(myfit.phase)])
        ax_lc.set_xlim([min(myfit.phase), max(myfit.phase)])

        # making borders and tick labels black
        ax_lc.spines['bottom'].set_color('black')
        ax_lc.spines['top'].set_color('black')
        ax_lc.spines['right'].set_color('black')
        ax_lc.spines['left'].set_color('black')
        ax_lc.tick_params(axis='x', colors='black')
        ax_lc.tick_params(axis='y', colors='black')

        ax_res.spines['bottom'].set_color('black')
        ax_res.spines['top'].set_color('black')
        ax_res.spines['right'].set_color('black')
        ax_res.spines['left'].set_color('black')
        ax_res.tick_params(axis='x', colors='black')
        ax_res.tick_params(axis='y', colors='black')

        # residual plot
        ax_res.errorbar(myfit.phase, myfit.residuals / np.median(myfit.data), yerr=myfit.detrendederr, color='gray',
                        marker='o', markersize=5, linestyle='None', mec='None', alpha=0.75)
        ax_res.plot(myfit.phase, np.zeros(len(myfit.phase)), 'r-', lw=2, alpha=1, zorder=100)
        ax_res.set_ylabel('Residuals')
        ax_res.set_ylim([-3 * np.nanstd(myfit.residuals / np.median(myfit.data)),
                         3 * np.nanstd(myfit.residuals / np.median(myfit.data))])

        correctedSTD = np.std(myfit.residuals / np.median(myfit.data))
        ax_lc.errorbar(myfit.phase, myfit.detrended, yerr=myfit.detrendederr, ls='none',
                       marker='o', color='gray', markersize=5, mec='None', alpha=0.75)
        ax_lc.plot(myfit.phase, myfit.transit, 'r', zorder=1000, lw=2)

        ax_lc.set_ylabel('Relative Flux')
        ax_lc.get_xaxis().set_visible(False)

        ax_res.errorbar(binner(myfit.phase, len(myfit.residuals) // 10),
                        binner(myfit.residuals / np.median(myfit.data), len(myfit.residuals) // 10),
                        yerr=
                        binner(myfit.residuals / np.median(myfit.data), len(myfit.residuals) // 10, myfit.detrendederr)[
                            1],
                        fmt='s', ms=5, mfc='b', mec='None', ecolor='b', zorder=10)
        ax_lc.errorbar(binner(myfit.phase, len(myfit.phase) // 10),
                       binner(myfit.detrended, len(myfit.detrended) // 10),
                       yerr=
                       binner(myfit.residuals / np.median(myfit.data), len(myfit.residuals) // 10, myfit.detrendederr)[
                           1],
                       fmt='s', ms=5, mfc='b', mec='None', ecolor='b', zorder=10)

        # remove vertical whitespace
        f.subplots_adjust(hspace=0)

        # For some reason, saving as a pdf crashed on Rob's laptop...so adding in a try statement to save it as a pdf if it can, otherwise, png
        Path(exotic_infoDict['saveplot']).mkdir(parents=True, exist_ok=True)
        try:
            f.savefig(Path(exotic_infoDict['saveplot']) /
                      f"FinalLightCurve_{pDict['pName']}_{exotic_infoDict['date']}.pdf", bbox_inches="tight")
            f.savefig(Path(exotic_infoDict['saveplot']) /
                      f"FinalLightCurve_{pDict['pName']}_{exotic_infoDict['date']}.png", bbox_inches="tight")
        except:
            f.savefig(Path(exotic_infoDict['saveplot']) /
                      f"FinalLightCurve_{pDict['pName']}_{exotic_infoDict['date']}.png", bbox_inches="tight")
        plt.close()

        ###################################################################################

        # triangle plot
        fig, axs = dynesty.plotting.cornerplot(myfit.results, labels=list(mybounds.keys()), quantiles_2d=[0.4, 0.85],
                                               smooth=0.015, show_titles=True, use_math_text=True, title_fmt='.2e',
                                               hist2d_kwargs={'alpha': 1, 'zorder': 2, 'fill_contours': False})
        dynesty.plotting.cornerpoints(myfit.results, labels=list(mybounds.keys()),
                                      fig=[fig, axs[1:, :-1]], plot_kwargs={'alpha': 0.1, 'zorder': 1, })
        fig.savefig(Path(exotic_infoDict['saveplot']) / "temp" /
                    f"Triangle_{pDict['pName']}_{exotic_infoDict['date']}.png")
        plt.close()

        # write output to text file
        params_file = Path(
            exotic_infoDict['saveplot']) / f"FinalLightCurve_{pDict['pName']}_{exotic_infoDict['date']}.csv"
        with params_file.open('w') as f:
            f.write(f"# FINAL TIMESERIES OF {pDict['pName']}\n")
            f.write("# BJD_TDB,Orbital Phase,Flux,Uncertainty,Model,Airmass\n")
            phase = getPhase(myfit.time, pDict['pPer'], myfit.parameters['tmid'])

            for bjdi, phasei, fluxi, fluxerri, modeli, ami in zip(myfit.time, phase, myfit.detrended,
                                                                  myfit.dataerr / myfit.airmass_model, myfit.transit,
                                                                  myfit.airmass_model):
                f.write(f"{bjdi}, {phasei}, {fluxi}, {fluxerri}, {modeli}, {ami}\n")

        ##########
        # PSF data
        ##########

        fig, ax = plt.subplots(3,2, figsize=(12,10))
        fig.suptitle(f"Observing Statistics - Target - {exotic_infoDict['date']}")
        ax[0,0].plot(myfit.time, psf_data['target'][si,0][gi], 'k.')
        ax[0,0].set_ylabel("X-Centroid [px]")
        ax[0,1].plot(myfit.time, psf_data['target'][si,1][gi], 'k.')
        ax[0,1].set_ylabel("Y-Centroid [px]")
        ax[1,0].plot(myfit.time, 2.355*0.5*(psf_data['target'][si,3][gi] + psf_data['target'][si,4][gi]), 'k.')
        ax[1,0].set_ylabel("Seeing [px]")
        ax[1,1].plot(myfit.time, myfit.airmass, 'k.')
        ax[1,1].set_ylabel("Airmass")
        ax[2,0].plot(myfit.time, psf_data['target'][si,2][gi], 'k.')
        ax[2,1].plot(myfit.time, psf_data['target'][si,6][gi], 'k.')
        ax[2,0].set_ylabel("Amplitude [ADU]")
        ax[2,1].set_ylabel("Background [ADU]")
        ax[0,0].set_xlabel("Time [BJD]")
        ax[0,1].set_xlabel("Time [BJD]")
        ax[1,0].set_xlabel("Time [BJD]")
        ax[1,1].set_xlabel("Time [BJD]")
        ax[2,0].set_xlabel("Time [BJD]")
        ax[2,1].set_xlabel("Time [BJD]")
        plt.tight_layout()

        try:
            fig.savefig(Path(exotic_infoDict['saveplot']) /
                        f"Observing_Statistics_target_{exotic_infoDict['date']}.png", bbox_inches="tight")
        except:
            pass
        fig.savefig(Path(exotic_infoDict['saveplot']) /
                    f"Observing_Statistics_target_{exotic_infoDict['date']}.pdf", bbox_inches="tight")
        plt.close()

        # PSF DATA for COMP STARS
        for j,coord in enumerate(compStarList):
            ctitle = "Comp Star {}".format(j+1)
            ckey = "comp{}".format(j+1)

            fig, ax = plt.subplots(3,2, figsize=(12,10))
            fig.suptitle(f"Observing Statistics - {ctitle} - {exotic_infoDict['date']}")
            ax[0,0].plot(myfit.time, psf_data[ckey][si,0][gi], 'k.')
            ax[0,0].set_ylabel("X-Centroid [px]")
            ax[0,1].plot(myfit.time, psf_data[ckey][si,1][gi], 'k.')
            ax[0,1].set_ylabel("Y-Centroid [px]")
            ax[1,0].plot(myfit.time, 2.355*0.5*(psf_data[ckey][si,3][gi] + psf_data[ckey][si,4][gi]), 'k.')
            ax[1,0].set_ylabel("Seeing [px]")
            ax[1,1].plot(myfit.time, myfit.airmass, 'k.')
            ax[1,1].set_ylabel("Airmass")
            ax[2,0].plot(myfit.time, psf_data[ckey][si,2][gi], 'k.')
            ax[2,1].plot(myfit.time, psf_data[ckey][si,6][gi], 'k.')
            ax[2,0].set_ylabel("Amplitude [ADU]")
            ax[2,1].set_ylabel("Background [ADU]")
            ax[0,0].set_xlabel("Time [BJD_TBD]")
            ax[0,1].set_xlabel("Time [BJD_TBD]")
            ax[1,0].set_xlabel("Time [BJD_TBD]")
            ax[1,1].set_xlabel("Time [BJD_TBD]")
            ax[2,0].set_xlabel("Time [BJD_TBD]")
            ax[2,1].set_xlabel("Time [BJD_TBD]")
            plt.tight_layout()

            try:
                fig.savefig(Path(exotic_infoDict['saveplot']) /
                            f"Observing_Statistics_{ckey}_{exotic_infoDict['date']}.pdf", bbox_inches="tight")
            except:
                pass
            fig.savefig(Path(exotic_infoDict['saveplot']) /
                        f"Observing_Statistics_{ckey}_{exotic_infoDict['date']}.png", bbox_inches="tight")
            plt.close()


        #######################################################################
        # print final extracted planetary parameters
        #######################################################################

        log.info("\n*********************************************************")
        log.info("FINAL PLANETARY PARAMETERS\n")
        log.info(
            f"              Mid-Transit Time [BJD_TDB]: {round_to_2(myfit.parameters['tmid'], myfit.errors['tmid'])} +/- {round_to_2(myfit.errors['tmid'])}")
        log.info(
            f"  Radius Ratio (Planet/Star) [Rp/Rs]: {round_to_2(myfit.parameters['rprs'], myfit.errors['rprs'])} +/- {round_to_2(myfit.errors['rprs'])}")
        log.info(
            f" Semi Major Axis/ Star Radius [a/Rs]: {round_to_2(myfit.parameters['ars'], myfit.errors['ars'])} +/- {round_to_2(myfit.errors['ars'])}")
        log.info(
            f"               Airmass coefficient 1: {round_to_2(myfit.parameters['a1'], myfit.errors['a1'])} +/- {round_to_2(myfit.errors['a1'])}")
        log.info(
            f"               Airmass coefficient 2: {round_to_2(myfit.parameters['a2'], myfit.errors['a2'])} +/- {round_to_2(myfit.errors['a2'])}")
        log.info(
            f"                    Residual scatter: {round_to_2(100. * np.std(myfit.residuals / np.median(myfit.data)))} %")
        log.info(
            f"              Transit Duration [day]: {round_to_2(np.mean(durs))} +/- {round_to_2(np.std(durs))}")
        log.info("*********************************************************")

        ##########
        # SAVE DATA
        ##########

        params_file = Path(exotic_infoDict['saveplot']) / f"FinalParams_{pDict['pName']}_{exotic_infoDict['date']}.json"
        params_num = {
            "Mid-Transit Time (Tmid)": f"{round_to_2(myfit.parameters['tmid'], myfit.errors['tmid'])} +/- {round_to_2(myfit.errors['tmid'])} BJD_TDB",
            "Ratio of Planet to Stellar Radius (Rp/Rs)": f"{round_to_2(myfit.parameters['rprs'], myfit.errors['rprs'])} +/- {round_to_2(myfit.errors['rprs'])}",
            "Transit depth (Rp/Rs)^2": f"{round_to_2(100. * (myfit.parameters['rprs'] ** 2.))} +/- {round_to_2(100. * 2. * myfit.parameters['rprs'] * myfit.errors['rprs'])} [%]",
            "Semi Major Axis/Star Radius (a/Rs)": f"{round_to_2(myfit.parameters['ars'], myfit.errors['ars'])} +/- {round_to_2(myfit.errors['ars'])} ",
            "Airmass coefficient 1 (a1)": f"{round_to_2(myfit.parameters['a1'], myfit.errors['a1'])} +/- {round_to_2(myfit.errors['a1'])}",
            "Airmass coefficient 2 (a2)": f"{round_to_2(myfit.parameters['a2'], myfit.errors['a2'])} +/- {round_to_2(myfit.errors['a2'])}",
            "Scatter in the residuals of the lightcurve fit is": f"{round_to_2(100. * np.std(myfit.residuals / np.median(myfit.data)))} %",
            "Transit Duration (day)":f"{round_to_2(np.mean(durs))} +/- {round_to_2(np.std(durs))}"
        }
        final_params = {'FINAL PLANETARY PARAMETERS': params_num}

        # write output to json file
        with params_file.open('w') as f:
            json.dump(final_params, f, indent=4)

        log.info(f"\nFinal Planetary Parameters have been saved in {exotic_infoDict['saveplot']} as "
                 f"{pDict['pName']}_{exotic_infoDict['date']}.json\n")

        if bestCompStar:
            comp_ra = None
            comp_dec = None

            if wcs_file:
                comp_ra = rafile[comp_coords[1]][comp_coords[0]]
                comp_dec = decfile[comp_coords[1]][comp_coords[0]]

            comp_star = [{'ra': str(comp_ra) if comp_ra else comp_ra,
                          'dec': str(comp_dec) if comp_dec else comp_dec,
                          'x': str(comp_coords[0]) if comp_coords[0] else comp_coords[0],
                          'y': str(comp_coords[1]) if comp_coords[1] else comp_coords[1]}]
        else:
            comp_star = []

        filter_dict = {'name': exotic_infoDict['filter'],
                       'fwhm': [
                           str(exotic_infoDict['wl_min']) if exotic_infoDict['wl_min'] else exotic_infoDict['wl_min'],
                           str(exotic_infoDict['wl_max']) if exotic_infoDict['wl_max'] else exotic_infoDict['wl_max']]}

        priors_dict = {'Period': {'value': str(round_to_2(pDict['pPer'], pDict['pPerUnc'])),
                                  'uncertainty': str(round_to_2(pDict['pPerUnc'])) if pDict['pPerUnc'] else pDict[
                                      'pPerUnc']},
                       'a/R*': {'value': str(round_to_2(pDict['aRs'], pDict['aRsUnc'])),
                                'uncertainty': str(round_to_2(pDict['aRsUnc'])) if pDict['aRsUnc'] else pDict[
                                    'aRsUnc']},
                       'inc': {'value': str(round_to_2(pDict['inc'], pDict['incUnc'])),
                               'uncertainty': str(round_to_2(pDict['incUnc'])) if pDict['incUnc'] else pDict['incUnc']},
                       'ecc': {'value': str(round_to_2(pDict['ecc'])), 'uncertainty': None},
                       'u0': {'value': str(round_to_2(ld0[0], ld0[1])), 'uncertainty': str(round_to_2(ld0[1]))},
                       'u1': {'value': str(round_to_2(ld1[0], ld1[1])), 'uncertainty': str(round_to_2(ld1[1]))},
                       'u2': {'value': str(round_to_2(ld2[0], ld2[1])), 'uncertainty': str(round_to_2(ld2[1]))},
                       'u3': {'value': str(round_to_2(ld3[0], ld3[1])), 'uncertainty': str(round_to_2(ld3[1]))}}

        round_to_2(myfit.parameters['a1'], myfit.errors['a1'])

        results_dict = {'Tc': {'value': str(round_to_2(myfit.parameters['tmid'], myfit.errors['tmid'])),
                               'uncertainty': str(round_to_2(myfit.errors['tmid']))},
                        'Rp/R*': {'value': str(round_to_2(myfit.parameters['rprs'], myfit.errors['rprs'])),
                                  'uncertainty': str(round_to_2(myfit.errors['rprs']))},
                        'Am1': {'value': str(round_to_2(myfit.parameters['a1'], myfit.errors['a1'])),
                                'uncertainty': str(round_to_2(myfit.errors['a1']))},
                        'Am2': {'value': str(round_to_2(myfit.parameters['a2'], myfit.errors['a2'])),
                                'uncertainty': str(round_to_2(myfit.errors['a2']))},
                        'Duration':{'value':str(round_to_2(np.mean(durs))),
                                    'uncertainty':str(round_to_2(np.std(durs)))}}

        params_file = Path(exotic_infoDict['saveplot']) / f"AAVSO_{pDict['pName']}_{exotic_infoDict['date']}.txt"
        with params_file.open('w') as f:
            f.write("#TYPE=EXOPLANET\n"  # fixed
                    f"#OBSCODE={exotic_infoDict['aavsonum']}\n"  # UI
                    f"#SECONDARY_OBSCODE={exotic_infoDict['secondobs']}\n"  # UI
                    f"#SOFTWARE=EXOTIC v{__version__}\n"  # fixed
                    "#DELIM=,\n"  # fixed
                    "#DATE_TYPE=BJD_TDB\n"  # fixed
                    f"#OBSTYPE={exotic_infoDict['ctype']}\n"
                    f"#STAR_NAME={pDict['sName']}\n"  # code yields
                    f"#EXOPLANET_NAME={pDict['pName']}\n"  # code yields
                    f"#BINNING={exotic_infoDict['pixelbin']}\n"  # user input
                    f"#EXPOSURE_TIME={exotic_infoDict.get('exposure', -1)}\n"  # UI
                    f"#FILTER-XC={json.dumps(filter_dict)}\n"
                    f"#COMP_STAR-XC={json.dumps(comp_star)}\n"
                    f"#NOTES={exotic_infoDict['notes']}\n"
                    "#DETREND_PARAMETERS=AIRMASS, AIRMASS CORRECTION FUNCTION\n"  # fixed
                    "#MEASUREMENT_TYPE=Rnflux\n"  # fixed
                    f"#PRIORS-XC={json.dumps(priors_dict)}\n"  # code yields
                    f"#RESULTS-XC={json.dumps(results_dict)}\n")  # code yields

            # Older formatting, will remove later
            f.write(previous_data_format(pDict, ld0, ld1, ld2, ld3, myfit))

            f.write("# EXOTIC is developed by Exoplanet Watch (exoplanets.nasa.gov/exoplanet-watch/), a citizen science project managed by NASA’s Jet Propulsion Laboratory on behalf of NASA’s Universe of Learning. This work is supported by NASA under award number NNX16AC65A to the Space Telescope Science Institute.\n"
                    "# Use of this data is governed by the AAVSO Data Usage Guidelines: aavso.org/data-usage-guidelines\n")

            f.write("#DATE,FLUX,MERR,DETREND_1,DETREND_2\n")
            for aavsoC in range(0, len(myfit.time)):
                # f.write(f"{round(myfit.time[aavsoC], 8)},{round(myfit.data[aavsoC] / myfit.parameters['a1'], 7)},"
                #         f"{round(myfit.dataerr[aavsoC] / myfit.parameters['a1'], 7)},{round(goodAirmasses[aavsoC], 7)},"
                #         f"{round(myfit.airmass_model[aavsoC] / myfit.parameters['a1'], 7)}\n")
                f.write(f"{round(myfit.time[aavsoC], 8)},{round(myfit.data[aavsoC], 7)},"
                        f"{round(myfit.dataerr[aavsoC], 7)},{round(goodAirmasses[aavsoC], 7)},"
                        f"{round(myfit.airmass_model[aavsoC], 7)}\n")

        log.info("Output File Saved")

        log.info("\n************************")
        log.info("End of Reduction Process")
        log.info("************************")

        log.debug("Stopped ...")


if __name__ == "__main__":
    # configure logger for standalone execution
    logging.root.setLevel(logging.NOTSET)
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
