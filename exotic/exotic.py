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
# from astroscrappy import detect_cosmics
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
import matplotlib.patheffects as PathEffects

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
from scipy.ndimage import binary_dilation, label, binary_erosion
# cross correlation imports
from skimage.transform import SimilarityTransform
# error handling for scraper
from tenacity import retry, retry_if_exception_type, stop_after_attempt, \
    stop_after_delay, wait_exponential
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
try: # output files
    from output_files import OutputFiles
except ImportError:  # package import
    from .output_files import OutputFiles
try: # tools
    from util import round_to_2, user_input
except ImportError: # package import
    from .util import round_to_2, user_input
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


def log_info(string):
    print(string)
    log.debug(string)
    return True


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


# ################### START ARCHIVE SCRAPER (PRIORS) ##########################
class NASAExoplanetArchive:

    def __init__(self, planet=None, candidate=False):
        self.planet = planet
        # self.candidate = candidate
        self.pl_dict = None

        # CONFIGURATIONS
        self.requests_timeout = 16, 512  # connection timeout, response timeout in secs.

    def planet_info(self, fancy=False):
        log_info(f"\nLooking up {self.planet} on the NASA Exoplanet Archive. Please wait....")
        self.planet, candidate = self._new_scrape(filename="eaConf.json", target=self.planet)

        if not candidate:
            with open("eaConf.json", "r") as confirmed:
                data = json.load(confirmed)
                planets = [data[i]['pl_name'] for i in range(len(data))]
                idx = planets.index(self.planet)
                self._get_params(data[idx])
                log_info(f"Successfully found {self.planet} in the NASA Exoplanet Archive!")

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
                "Orbital Inclination (deg) Uncertainty": self.pl_dict['incUnc'],
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
        # log_info(uri_full)
        log.debug(uri_full)

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
                            # log_info(i,k,edata[k][edata[k].notna()].values[0])
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


def planet_name(planet):
    while True:
        if not planet:
            planet = user_input("\nEnter the Planet Name. Make sure it matches the case sensitive name "
                                "and spacing used on Exoplanet Archive "
                                "(https://exoplanetarchive.ipac.caltech.edu/index.html): ", type_=str)

        if not planet[-2].isspace():
            log_info("The convention on the NASA Exoplanet Archive "
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
                log_info(f"\n\n*** WARNING: {pdict['pName']} initialization file's {planet_params[idx]} does not match "
                         "the value scraped by EXOTIC from the NASA Exoplanet Archive. ***\n")
                log_info(f"\tNASA Exoplanet Archive value (degrees): {pdict[item]}")
                log_info(f"\tInitialization file value (degrees): {userpdict[item]}")
                log_info("\nWould you like to:"
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
                log_info(f"\n\n*** WARNING: {pdict['pName']} initialization file's {planet_params[i]} does not match "
                         "the value scraped by EXOTIC from the NASA Exoplanet Archive. ***\n")
                log_info(f"\tNASA Exoplanet Archive value: {pdict[key]}")
                log_info(f"\tInitialization file value: {userpdict[key]}")
                log_info("\nWould you like to: "
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
                log_info(f"\n {pdict['pName']} {planet_params[i]}: {pdict[key]}")
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
            log_info("Error: The format entered for Right Ascension and/or Declination is not correct, "
                     "please try again.")
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
                log_info("\nError: The entered filter is not in the provided list of standard filters.")
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
            log_info(f"Found corrupted file and removing from reduction: {file}, ({e})")
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

    if plate_opt == 'y':
        return get_wcs(fits_file, save_directory)
    elif plate_opt == 'n':
        if wcs_exists:
            log_info("Your FITS files have WCS (World Coordinate System) information in their headers. "
                     "EXOTIC will proceed to use these. "
                     "NOTE: If you do not trust your WCS coordinates, "
                     "please restart EXOTIC after enabling plate solutions via astrometry.net.")
            return fits_file
        else:
            return False


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


# Check the ra and dec against the plate solution to see if the user entered in the correct values
def check_targetpixelwcs(pixx, pixy, expra, expdec, ralist, declist):
    while True:
        try:
            uncert = 1 / 36
            # Margins are within 100 arcseconds
            if not (expra - uncert <= ralist[int(pixy)][int(pixx)] <= expra + uncert):
                log_info("\n*** Warning: The X Pixel Coordinate entered does not match the target's Right Ascension. ***")
                raise ValueError
            if not (expdec - uncert <= declist[int(pixy)][int(pixx)] <= expdec + uncert):
                log_info("\n*** Warning: The Y Pixel Coordinate entered does not match the target's Declination. ***")
                raise ValueError
            return pixx, pixy
        except ValueError:
            log_info(f"Your input pixel coordinates: [{pixx}, {pixy}]")

            dist = (ralist - expra) ** 2 + (declist - expdec) ** 2
            exotic_pixy, exotic_pixx = np.unravel_index(dist.argmin(), dist.shape)
            log_info(f"EXOTIC's calculated pixel coordinates: [{exotic_pixx}, {exotic_pixy}]")

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
    radius = u.Quantity(20.0, u.arcsec)

    # Query GAIA first to check for variability using the phot_variable_flag trait
    gaia_result = gaia_query(sample, radius)
    if not gaia_result:
        log_info("*** WARNING: Your comparison star cannot be resolved in the Gaia star database; "
                 "EXOTIC cannot check if it is variable or not. "
                 "***\nEXOTIC will still include this star in the reduction. "
                 "\nPlease proceed with caution as we cannot check for stellar variability.\n")
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
        log_info("*** WARNING: Your comparison star cannot be resolved in the SIMBAD star database; "
                 "EXOTIC cannot check if it is variable or not. "
                 "***\nEXOTIC will still include this star in the reduction. "
                 "\nPlease proceed with caution as we cannot check for stellar variability.\n")
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


# Aligns imaging data from .fits file to easily track the host and comparison star's positions
def transformation(image_data, num_images, file_name, count, roi=1):
    pos = np.zeros((1, 2))

    # crop image to ROI
    height = image_data.shape[1]
    width = image_data.shape[2]
    roix = slice(int(width * (0.5 - roi / 2)), int(width * (0.5 + roi / 2)))
    roiy = slice(int(height * (0.5 - roi / 2)), int(height * (0.5 + roi / 2)))

    # Find transformation from .FITS files and catch exceptions if not able to.
    try:
        sys.stdout.write(f"Finding transformation {count + 1} of {num_images}\r")
        log.debug(f"Finding transformation {count + 1} of {num_images}\r")
        sys.stdout.flush()

        results = aa.find_transform(image_data[1][roiy, roix], image_data[0][roiy, roix])
        return results[0]
    except Exception as ee:
        log_info(ee)

        for p in [99, 98, 95, 90]:
            for it in [2,1,0]:

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
    
    log_info(f"alignment failed: {file_name}")
    return SimilarityTransform(scale=1, rotation=0, translation=[0,0])


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
        log_info("Not able to find Pixel Scale in the Image Header.")
        image_scale_num = user_input("Please enter the size of your pixel (e.g., 5 arc-sec/pixel): ", type_=float)
        image_scale = f"Image scale in arc-secs/pixel: {image_scale_num}"
    return image_scale


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


def gaussian_psf(x, y, x0, y0, a, sigx, sigy, rot, b):
    rx = (x - x0) * np.cos(rot) - (y - y0) * np.sin(rot)
    ry = (x - x0) * np.sin(rot) + (y - y0) * np.cos(rot)
    gausx = np.exp(-rx ** 2 / (2 * sigx ** 2))
    gausy = np.exp(-ry ** 2 / (2 * sigy ** 2))
    return a * gausx * gausy + b


def mesh_box(pos, box):
    pos = [int(np.round(pos[0])), int(np.round(pos[1]))]
    x = np.arange(pos[0] - box, pos[0] + box + 1)
    y = np.arange(pos[1] - box, pos[1] + box + 1)
    xv, yv = np.meshgrid(x, y)
    return xv.astype(int), yv.astype(int)


# Method fits a 2D gaussian function that matches the star_psf to the star image and returns its pixel coordinates
def fit_centroid(data, pos, init=[], psf_function=gaussian_psf, box=10):
    # get sub field in image
    xv, yv = mesh_box(pos, box)

    init = [np.nanmax(data[yv, xv]) - np.nanmin(data[yv, xv]), 1, 1, 0, np.nanmin(data[yv, xv])]

    # lower bound: [xc, yc, amp, sigx, sigy, rotation,  bg]
    lo = [pos[0] - box * 0.5, pos[1] - box * 0.5, 0, 0.5, 0.5, -np.pi / 4, np.nanmin(data[yv,xv]) - 1]
    up = [pos[0] + box * 0.5, pos[1] + box * 0.5, 1e7, 20, 20, np.pi / 4, np.nanmax(data[yv,xv]) + 1]

    def fcn2min(pars):
        model = psf_function(xv, yv, *pars)
        return (data[yv, xv] - model).flatten()

    try:
        res = least_squares(fcn2min, x0=[*pos, *init], bounds=[lo, up], jac='3-point', xtol=None, method='trf')
    except:
        log_info(f"CAUTION: Measured flux amplitude is really low---are you sure there is a star at {np.round(pos, 2)}?")

        res = least_squares(fcn2min, x0=[*pos, *init], jac='3-point', xtol=None, method='lm')

    return res.x


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
        log_info(f"IndexError, problem computing sky bg for {xc:.1f}, {yc:.1f}. "
                 f"Check if star is present or close to border.")

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
    plt.savefig(Path(exotic_infoDict['save']) / "temp" / f"XCentroidPos_{targetname}_{date}.png")
    plt.close()

    # Y TARGET
    plt.figure()
    plt.plot(times - np.nanmin(times), yTarg, '-bo')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('Y Pixel Position')
    plt.title(targetname + ' Y Centroid Position ' + date)
    plt.savefig(Path(exotic_infoDict['save']) / "temp" / f"YCentroidPos_{targetname}_{date}.png")
    plt.close()

    # X COMP
    plt.figure()
    plt.plot(times - np.nanmin(times), xRef, '-ro')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('X Pixel Position')
    plt.title('Comp Star X Centroid Position ' + date)
    plt.savefig(Path(exotic_infoDict['save']) / "temp" / f"CompStarXCentroidPos_{targetname}_{date}.png")
    plt.close()

    # Y COMP
    plt.figure()
    plt.plot(times - np.nanmin(times), yRef, '-ro')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('Y Pixel Position')
    plt.title('Comp Star Y Centroid Position ' + date)
    plt.savefig(Path(exotic_infoDict['save']) / "temp" / f"CompStarYCentroidPos_{targetname}_{date}.png")
    plt.close()

    # X DISTANCE BETWEEN TARGET AND COMP
    plt.figure()
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('X Pixel Distance')
    for e in range(0, len(xTarg)):
        plt.plot(times[e] - np.nanmin(times), abs(int(xTarg[e]) - int(xRef[e])), 'bo')
    plt.title('Distance between Target and Comparison X position')
    plt.savefig(Path(exotic_infoDict['save']) / "temp" / f"XCentroidDistance_{targetname}_{date}.png")
    plt.close()

    # Y DISTANCE BETWEEN TARGET AND COMP
    plt.figure()
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('Y Pixel Difference')

    for d in range(0, len(yTarg)):
        plt.plot(times[d] - np.nanmin(times), abs(int(yTarg[d]) - int(yRef[d])), 'bo')
    plt.title('Difference between Target and Comparison Y position')
    plt.savefig(Path(exotic_infoDict['save']) / "temp" / f"YCentroidDistance_{targetname}_{date}.png")
    plt.close()


def psf_format(data, pos, init=[]):
    target = fit_centroid(data, pos, init=init)

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
        tform = transformation(np.array([imageData, firstImageData]), len(timeSortedNames), imageFile, i)

        # apply transform
        tx, ty = tform([UIprevTPX, UIprevTPY])[0]
        rx, ry = tform([UIprevRPX, UIprevRPY])[0]

        print(f"Coordinates of target star: ({tx},{ty})")
        print(f"Coordinates of comp star: ({rx},{ry})")

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
    ax.set_xlabel('Time (JD)')
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
        log_info("\nWARNING!")
        log_info(" Estimated mid-transit time is not within the observations")
        log_info(" Check Period & Mid-transit time in inits.json. Make sure the uncertainties are not 0 or Nan.")
        log_info(f"  obs start:{arrayTimes.min()}")
        log_info(f"    obs end:{arrayTimes.max()}")
        log_info(f" tmid prior:{prior['tmid']}\n")

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

    log_info("\n*************************************************************")
    log_info("Welcome to the EXOplanet Transit Interpretation Code (EXOTIC)")
    log_info(f"Version {__version__}")
    log_info("*************************************************************\n")

    # ---INITIALIZATION-------------------------------------------------------
    global exotic_infoDict

    fileNameList, timeSortedNames, xTargCent, yTargCent, xRefCent, yRefCent, finXTargCent, finYTargCent, finXRefCent, finYRefCent = (
        [] for m in range(10))

    minSTD = 100000  # sets the initial minimum standard deviation absurdly high so it can be replaced immediately
    minChi2 = 100000
    distFC = 10  # gaussian search area

    userpDict = {'ra': None, 'dec': None, 'pName': None, 'sName': None, 'pPer': None, 'pPerUnc': None,
                 'midT': None, 'midTUnc': None, 'rprs': None, 'rprsUnc': None, 'aRs': None, 'aRsUnc': None,
                 'inc': None, 'incUnc': None, 'ecc': None, 'teff': None,
                 'teffUncPos': None, 'teffUncNeg': None, 'met': None, 'metUncPos': None, 'metUncNeg': None,
                 'logg': None, 'loggUncPos': None, 'loggUncNeg': None}

    # ---USER INPUTS--------------------------------------------------------------------------
    if args.realtime:
        reduction_opt = 1
    elif args.reduce or args.prereduced:
        reduction_opt = 2
    else:
        reduction_opt = user_input("\nPlease select: \n\t1: for Real Time Reduction (for analyzing your data while observing) \n\t2: for for Complete Reduction (for analyzing your data after an observing run). \nPlease enter 1 or 2: ",
                                   type_=int, val1=1, val2=2)

    if reduction_opt == 1:
        log_info("\n**************************************************************")
        log_info("Real Time Reduction ('Control + C'  or close the plot to quit)")
        log_info("**************************************************************\n")

        if args.realtime:
            init_opt = 'y'
        else:
            init_opt = 'n'

        inputs_obj = Inputs(init_opt=init_opt)

        if init_opt == 'y':
            init_path, userpDict = inputs_obj.search_init(args.realtime, userpDict)

        userpDict['pName'] = planet_name(userpDict['pName'])

        exotic_infoDict = inputs_obj.real_time(userpDict['pName'])
        exotic_UIprevTPX = exotic_infoDict['tar_coords'][0]
        exotic_UIprevTPY = exotic_infoDict['tar_coords'][1]
        exotic_UIprevRPX = exotic_infoDict['comp_stars'][0]
        exotic_UIprevRPY = exotic_infoDict['comp_stars'][1]

        while True:
            carry_on = user_input(f"\nType continue after the first image has been taken and saved: ", type_=str)
            if carry_on != 'continue':
                continue
            break

        log_info("Real Time Plotting ('Control + C' or close the plot to quit)")
        log_info("\nPlease be patient. It will take at least 15 seconds for the first image to get plotted.")

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(userpDict['pName'])
        ax.set_ylabel('Normalized Flux')
        ax.set_xlabel('Time (JD)')

        anim = FuncAnimation(fig, realTimeReduce,
                             fargs=(userpDict['pName'], ax, distFC, exotic_infoDict['images'], exotic_UIprevTPX, exotic_UIprevTPY,
                                    exotic_UIprevRPX, exotic_UIprevRPY), interval=15000)  # refresh every 15 seconds
        plt.show()

    # ----USER INPUTS----------------------------------------------------------
    else:
        log_info("\n**************************")
        log_info("Complete Reduction Routine")
        log_info("**************************")

        init_path, wcs_file, wcs_header = None, None, None

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

        if fileorcommandline == 2:
            init_opt = 'y'
        else:
            init_opt = 'n'

        inputs_obj = Inputs(init_opt=init_opt)

        if init_opt == 'y':
            init_path, userpDict = inputs_obj.search_init(init_path, userpDict)

        userpDict['pName'] = planet_name(userpDict['pName'])

        if fitsortext == 1:
            exotic_infoDict = inputs_obj.complete_red(userpDict['pName'])
        else:
            exotic_infoDict = inputs_obj.prereduced()

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
                    darkData = fits.getdata(darkFile, ext=0)
                    darksImgList.append(darkData)
                generalDark = np.median(darksImgList, axis=0)

            if exotic_infoDict['biases']:
                biasesImgList = []
                for biasFile in exotic_infoDict['biases']:
                    biasData = fits.getdata(biasFile, ext=0)
                    biasesImgList.append(biasData)
                generalBias = np.median(biasesImgList, axis=0)

            if exotic_infoDict['flats']:
                flatsImgList = []
                for flatFile in exotic_infoDict['flats']:
                    flatData = fits.getdata(flatFile, ext=0)
                    flatsImgList.append(flatData)
                notNormFlat = np.median(flatsImgList, axis=0)

                # if the bias exists, bias subtract the flatfield
                if exotic_infoDict['biases']:
                    notNormFlat = notNormFlat - generalBias

                # NORMALIZE
                medi = np.median(notNormFlat)
                generalFlat = notNormFlat / medi

        # log_info("***************************************\n")

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
        ld0, ld1, ld2, ld3, exotic_infoDict['filter'], exotic_infoDict['wl_min'], exotic_infoDict['wl_max'] = \
            ld_obj.nonlinear_ld()
        ld = [ld0[0], ld1[0], ld2[0], ld3[0]]

        # check for Nans + Zeros
        for k in pDict:
            if k=='rprs' and (pDict[k] == 0 or np.isnan(pDict[k])):
                log_info(f"ERROR! {k} value is 0 or NaN. Please use a non-zero value in inits.json")
                pDict[k] = 0.8 # instead of 1 since priors on RpRs are 0 to RpRs*1.25
                log_info("EXOTIC will override the Rp/Rs value.")
            if k=='aRs' and (pDict[k] < 1 or np.isnan(pDict[k])):
                log_info(f"WARNING! {k} value is <1 or NaN. Please use a non-zero value in inits.json")
                pDict[k] = user_input("\nPlease enter candidate planet's name: ", type_=float)
            if "Unc" in k:
                if not pDict[k]:
                    log_info(f"WARNING! {k} uncertainty is 0. Please use a non-zero value in inits.json")
                    pDict[k] = 1
                elif pDict[k] == 0 or np.isnan(pDict[k]):
                    log_info(f"WARNING! {k} uncertainty is 0. Please use a non-zero value in inits.json")
                    pDict[k] = 1
            elif pDict[k] is None:
                log_info(f"WARNING! {k} is None. Please use a numeric value in inits.json")
                pDict[k] = 0

        if fitsortext == 1:
            log_info("\n**************************"
                     "\nStarting Reduction Process"
                     "\n**************************\n")

            #########################################
            # FLUX DATA EXTRACTION AND MANIPULATION
            #########################################

            allImageData, timeList, airMassList, exptimes = [], [], [], []

            # TODO filter input files to get good reference for image alignment

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

            si = np.argsort(times)
            inputfiles = np.array(inputfiles)[si]
            exotic_UIprevTPX = exotic_infoDict['tar_coords'][0]
            exotic_UIprevTPY = exotic_infoDict['tar_coords'][1]

            # fit target in the first image and use it to determine aperture and annulus range
            inc = 0
            for ifile in inputfiles:
                first_image = fits.getdata(ifile, ext=0)
                try:
                    args = fit_centroid(first_image, [exotic_UIprevTPX, exotic_UIprevTPY])
                    break
                except Exception:
                    inc += 1
                finally:
                    del first_image

            inputfiles = inputfiles[inc:]
            wcs_file = check_wcs(inputfiles[0], exotic_infoDict['save'], exotic_infoDict['plate_opt'])
            compStarList = exotic_infoDict['comp_stars']

            if wcs_file:
                log_info(f"\nHere is the path to your plate solution: {wcs_file}")
                wcs_header = fits.getheader(filename=wcs_file)
                rafile, decfile = get_radec(wcs_header)

                # Checking pixel coordinates against plate solution
                exotic_UIprevTPX, exotic_UIprevTPY = check_targetpixelwcs(exotic_UIprevTPX, exotic_UIprevTPY,
                                                                          pDict['ra'], pDict['dec'], rafile, decfile)

                for compn, comp in enumerate(exotic_infoDict['comp_stars'][:]):
                    log_info("\nChecking for variability in Comparison Star #"+str(compn+1)+" : \n"
                             f"Pixel X: {comp[0]} Pixel Y: {comp[1]}")
                    if variableStarCheck(rafile[int(comp[1])][int(comp[0])], decfile[int(comp[1])][int(comp[0])]):
                            log_info("\nCurrent comparison star is variable, proceeding to next star.")
                            exotic_infoDict['comp_stars'].remove(comp)
                compStarList = exotic_infoDict['comp_stars']

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
                'target_bg': np.zeros((len(inputfiles), len(apers), len(annuli)))
            }

            for i, coord in enumerate(compStarList):
                ckey = "comp{}".format(i + 1)
                psf_data[ckey] = np.zeros((len(inputfiles), 7))
                psf_data[ckey + "_align"] = np.zeros((len(inputfiles), 2))
                aper_data[ckey] = np.zeros((len(inputfiles), len(apers), len(annuli)))
                aper_data[ckey + "_bg"] = np.zeros((len(inputfiles), len(apers), len(annuli)))

            if exotic_infoDict['img_align_opt'] == 'n':
                alignmentBool = False
            else:
                alignmentBool = True

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
                if exotic_infoDict['darks']:
                    if i == 0:
                        log_info("Dark subtracting images.")
                    imageData = imageData - generalDark
                elif exotic_infoDict['biases']: # if a dark is not available, then at least subtract off the pedestal via the bias
                    if i == 0:
                        log_info("Bias-correcting images.")
                    imageData = imageData - generalBias
                else:
                    pass

                if exotic_infoDict['flats']:
                    if i == 0:
                        log_info("Flattening images.")
                    generalFlat[generalFlat == 0] = 1
                    imageData = imageData / generalFlat

                if i == 0:
                    image_scale = get_pixel_scale(wcs_header, image_header, exotic_infoDict['pixel_scale'])

                    firstImage = np.copy(imageData)

                    if alignmentBool:
                        log_info("\n\nAligning your images. Please wait.")

                # Image Alignment
                if alignmentBool:
                    # log_info(f"\nReference Image for Alignment: {fileName}\n")
                    tform = transformation(np.array([imageData, firstImage]), len(inputfiles), fileName, i)
                else:
                    tform = SimilarityTransform(scale=1, rotation=0, translation=[0,0])

                # apply transform
                xrot, yrot = tform([exotic_UIprevTPX, exotic_UIprevTPY])[0]

                psf_data["target_align"][i] = [xrot, yrot]
                if i == 0:
                    psf_data["target"][i] = fit_centroid(imageData, [xrot, yrot])
                else:
                    if alignmentBool:
                        psf_data["target"][i] = fit_centroid(
                            imageData,
                            [xrot, yrot],
                            psf_data["target"][0][2:])
                    else:
                        # use previous PSF as prior
                        psf_data["target"][i] = fit_centroid(
                            imageData,
                            psf_data["target"][i-1][:2],
                            psf_data["target"][i-1][2:])

                        # check for change in amplitude of PSF
                        if np.abs( (psf_data["target"][i][2]-psf_data["target"][i-1][2])/psf_data["target"][i-1][2]) > 0.5:
                            log_info("Can't find target. Trying to align...")

                            tform = transformation(np.array([imageData, firstImage]), len(inputfiles), fileName, i)
                            
                            xrot, yrot = tform([exotic_UIprevTPX, exotic_UIprevTPY])[0]

                            psf_data["target"][i] = fit_centroid( imageData, [xrot, yrot], psf_data["target"][0][2:])

                # fit for the centroids in all images
                for j,coord in enumerate(compStarList):
                    ckey = "comp{}".format(j+1)

                    # apply transformation
                    xrot, yrot = tform(coord)[0]

                    psf_data[ckey+"_align"][i] = [xrot,yrot]
                    if i == 0:
                        psf_data[ckey][i] = fit_centroid(imageData, [xrot, yrot])
                    else:
                        if alignmentBool:
                            psf_data[ckey][i] = fit_centroid(
                                imageData,
                                [xrot, yrot],
                                psf_data[ckey][0][2:])
                        else:
                            # use previous PSF as prior
                            psf_data[ckey][i] = fit_centroid(
                                imageData,
                                psf_data[ckey][i-1][:2],
                                psf_data[ckey][i-1][2:])

                            # check for change in amplitude of PSF
                            if np.abs( (psf_data[ckey][i][2]-psf_data[ckey][i-1][2])/psf_data[ckey][i-1][2]) > 0.5:
                                #log_info("Can't find comp. Trying to align...")

                                tform = transformation(np.array([imageData, firstImage]), len(inputfiles), fileName, i)

                                xrot, yrot = tform(coord)[0]

                                psf_data[ckey][i] = fit_centroid( imageData, [xrot, yrot], psf_data[ckey][0][2:])

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
                log_info("No images to fit...check reference image for alignment (first image of sequence)")

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
                    comp_coords = coord
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

            log_info("Computing best comparison star, aperture, and sky annulus. Please wait.")

            # Aperture Photometry
            for a, aper in enumerate(apers):

                for an, annulus in enumerate(annuli):
                    tFlux = aper_data['target'][:, a, an]

                    # fit without a comparison star
                    myfit = fit_lightcurve(times, tFlux, np.ones(tFlux.shape[0]), airmass, ld, pDict)

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
                        arrayNormUnc = tFlux**0.5

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
            log_info("\n\n*********************************************")
            if minAperture == 0:  # psf
                log_info(f"Best Comparison Star: #{bestCompStar}")
                log_info(f"Minimum Residual Scatter: {round(minSTD * 100, 4)}%")
                log_info("Optimal Method: PSF photometry")
            elif minAperture < 0:  # no comp star
                log_info("Best Comparison Star: None")
                log_info(f"Minimum Residual Scatter: {round(minSTD * 100, 4)}%")
                log_info(f"Optimal Aperture: {abs(np.round(minAperture, 2))}")
                log_info(f"Optimal Annulus: {np.round(minAnnulus, 2)}")
                bestCompStar = None
            else:
                log_info(f"Best Comparison Star: #{bestCompStar}")
                log_info(f"Minimum Residual Scatter: {round(minSTD * 100, 4)}%")
                log_info(f"Optimal Aperture: {np.round(minAperture, 2)}")
                log_info(f"Optimal Annulus: {np.round(minAnnulus, 2)}")
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
                    min([np.shape(firstImage)[1], max([finXTargCent[0], finXRefCent[0]]) + picframe])]
            plty = [max([0, min([finYTargCent[0], finYRefCent[0]]) - picframe]),
                    min([np.shape(firstImage)[0], max([finYTargCent[0], finYRefCent[0]]) + picframe])]
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
                    text.set_color("k")
                    text.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='white')])
                apos = '\''
                Path(exotic_infoDict['save']).mkdir(parents=True, exist_ok=True)
                plt.savefig(Path(exotic_infoDict['save']) /
                            f"FOV_{pDict['pName']}_{exotic_infoDict['date']}_"
                            f"{str(stretch.__class__).split('.')[-1].split(apos)[0]}.pdf", bbox_inches='tight')
                plt.savefig(Path(exotic_infoDict['save']) /
                            f"FOV_{pDict['pName']}_{exotic_infoDict['date']}_"
                            f"{str(stretch.__class__).split('.')[-1].split(apos)[0]}.png", bbox_inches='tight')
                plt.close()

            log_info(f"\nFOV file saved as: {exotic_infoDict['save']}/FOV_{pDict['pName']}_"
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
            plt.savefig(Path(exotic_infoDict['save']) / "temp" /
                        f"TargetRawFlux_{pDict['pName']}_{exotic_infoDict['date']}.png")
            plt.close()

            plt.figure()
            plt.errorbar(goodTimes, goodReferences[si][gi], yerr=goodRUnc[si][gi], linestyle='None', fmt='-o')
            plt.xlabel('Time (BJD)')
            plt.ylabel('Total Flux')
            # plt.rc('grid', linestyle="-", color='black')
            # plt.grid(True)
            plt.title('Comparison Star Raw Flux Values ' + exotic_infoDict['date'])
            plt.savefig(Path(exotic_infoDict['save']) / "temp" /
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
            plt.savefig(Path(exotic_infoDict['save']) / "temp" /
                        f"NormalizedFluxTime_{pDict['pName']}_{exotic_infoDict['date']}.png")
            plt.close()

            # Save normalized flux to text file prior to MCMC
            params_file = Path(
                exotic_infoDict['save']) / f"NormalizedFlux_{pDict['pName']}_{exotic_infoDict['date']}.txt"
            with params_file.open('w') as f:
                f.write("BJD,Norm Flux,Norm Err,AM\n")

                for ti, fi, erri, ami in zip(goodTimes, goodFluxes, goodNormUnc, goodAirmasses):
                    f.write(f"{round(ti, 8)},{round(fi, 7)},{round(erri, 6)},{round(ami, 2)}\n")

            log_info("\n\nOutput File Saved")
        else:
            goodTimes, goodFluxes, goodNormUnc, goodAirmasses = [], [], [], []
            bestCompStar = None

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
                goodFluxes, goodNormUnc = fluxConvert(goodFluxes, goodNormUnc, exotic_infoDict['file_units'])

        # for k in myfit.bounds.keys():
        #     print("{:.6f} +- {}".format( myfit.parameters[k], myfit.errors[k]))

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
            log_info("ERROR: Estimated mid-transit not in observation range (check priors or observation time)")
            log_info(f"start:{goodTimes.min()}")
            log_info(f"  end:{goodTimes.max()}")
            log_info(f"prior:{prior['tmid']}")

        mybounds = {
            'rprs': [0, pDict['rprs'] * 1.25],
            'tmid': [lower, upper],
            'ars': [pDict['aRs'] - 5 * pDict['aRsUnc'], pDict['aRs'] + 5 * pDict['aRsUnc']],
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
        data_highres = transit(times, pars)
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
        ax_lc.plot(np.linspace(np.nanmin(myfit.phase), np.nanmax(myfit.phase), 1000), data_highres, 'r', zorder=1000, lw=2)

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
        Path(exotic_infoDict['save']).mkdir(parents=True, exist_ok=True)
        try:
            f.savefig(Path(exotic_infoDict['save']) /
                      f"FinalLightCurve_{pDict['pName']}_{exotic_infoDict['date']}.pdf", bbox_inches="tight")
            f.savefig(Path(exotic_infoDict['save']) /
                      f"FinalLightCurve_{pDict['pName']}_{exotic_infoDict['date']}.png", bbox_inches="tight")
        except:
            f.savefig(Path(exotic_infoDict['save']) /
                      f"FinalLightCurve_{pDict['pName']}_{exotic_infoDict['date']}.png", bbox_inches="tight")
        plt.close()


        ###################################################################################

        # triangle plot
        fig, axs = dynesty.plotting.cornerplot(myfit.results, labels=list(mybounds.keys()), quantiles_2d=[0.4, 0.85],
                                               smooth=0.015, show_titles=True, use_math_text=True, title_fmt='.2e',
                                               hist2d_kwargs={'alpha': 1, 'zorder': 2, 'fill_contours': False})
        dynesty.plotting.cornerpoints(myfit.results, labels=list(mybounds.keys()),
                                      fig=[fig, axs[1:, :-1]], plot_kwargs={'alpha': 0.1, 'zorder': 1, })
        fig.savefig(Path(exotic_infoDict['save']) / "temp" /
                    f"Triangle_{pDict['pName']}_{exotic_infoDict['date']}.png")
        plt.close()

        phase = getPhase(myfit.time, pDict['pPer'], myfit.parameters['tmid'])

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
            fig.savefig(Path(exotic_infoDict['save']) /
                        f"Observing_Statistics_target_{exotic_infoDict['date']}.png", bbox_inches="tight")
        except:
            pass
        fig.savefig(Path(exotic_infoDict['save']) /
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
                fig.savefig(Path(exotic_infoDict['save']) /
                            f"Observing_Statistics_{ckey}_{exotic_infoDict['date']}.pdf", bbox_inches="tight")
            except:
                pass
            fig.savefig(Path(exotic_infoDict['save']) /
                        f"Observing_Statistics_{ckey}_{exotic_infoDict['date']}.png", bbox_inches="tight")
            plt.close()


        #######################################################################
        # print final extracted planetary parameters
        #######################################################################

        log_info("\n*********************************************************")
        log_info("FINAL PLANETARY PARAMETERS\n")
        log_info(f"              Mid-Transit Time [BJD_TDB]: {round_to_2(myfit.parameters['tmid'], myfit.errors['tmid'])} +/- {round_to_2(myfit.errors['tmid'])}")
        log_info(f"  Radius Ratio (Planet/Star) [Rp/Rs]: {round_to_2(myfit.parameters['rprs'], myfit.errors['rprs'])} +/- {round_to_2(myfit.errors['rprs'])}")
        log_info(f" Semi Major Axis/ Star Radius [a/Rs]: {round_to_2(myfit.parameters['ars'], myfit.errors['ars'])} +/- {round_to_2(myfit.errors['ars'])}")
        log_info(f"               Airmass coefficient 1: {round_to_2(myfit.parameters['a1'], myfit.errors['a1'])} +/- {round_to_2(myfit.errors['a1'])}")
        log_info(f"               Airmass coefficient 2: {round_to_2(myfit.parameters['a2'], myfit.errors['a2'])} +/- {round_to_2(myfit.errors['a2'])}")
        log_info(f"                    Residual scatter: {round_to_2(100. * np.std(myfit.residuals / np.median(myfit.data)))} %")
        log_info(f"              Transit Duration [day]: {round_to_2(np.mean(durs))} +/- {round_to_2(np.std(durs))}")
        log_info("*********************************************************")

        comp_star = []

        if bestCompStar:
            comp_ra = None
            comp_dec = None

            if wcs_file:
                comp_ra = rafile[int(comp_coords[1])][int(comp_coords[0])]
                comp_dec = decfile[int(comp_coords[1])][int(comp_coords[0])]

            comp_star.append({
                'ra': str(comp_ra) if comp_ra else comp_ra,
                'dec': str(comp_dec) if comp_dec else comp_dec,
                'x': str(comp_coords[0]) if comp_coords[0] else comp_coords[0],
                'y': str(comp_coords[1]) if comp_coords[1] else comp_coords[1]
            })

        ##########
        # SAVE DATA
        ##########

        output_files = OutputFiles(myfit, pDict, exotic_infoDict, durs)
        error_txt = "\nPlease report this issue on the Exoplanet Watch Slack Channel in #data-reductions."

        try:
            output_files.final_lightcurve(phase)
        except Exception as e:
            log_info(f"Error: Could not create FinalLightCurve.csv. {error_txt}\n{e}")
        try:
            output_files.final_planetary_params()
        except Exception as e:
            log_info(f"Error: Could not create FinalParams.json. {error_txt}\n{e}")
        try:
            output_files.aavso(comp_star, goodAirmasses, ld0, ld1, ld2, ld3)
        except Exception as e:
            log_info(f"Error: Could not create AAVSO.txt. {error_txt}\n{e}")

        log_info("Output File Saved")

        log_info("\n************************")
        log_info("End of Reduction Process")
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
