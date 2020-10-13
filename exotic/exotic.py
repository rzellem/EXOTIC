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
# Image alignment import
import astroalign as aa
# astropy imports
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
from io import StringIO
import json
import logging
from logging.handlers import TimedRotatingFileHandler
from matplotlib.animation import FuncAnimation
# Pyplot imports
import matplotlib.pyplot as plt_exotic
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
# cross correlation imports
from skimage.registration import phase_cross_correlation
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
    from .api.elca import lc_fitter, binner
except ImportError:  # package import
    from api.elca import lc_fitter, binner
try:  # nonlinear limb darkening numerics
    from .api.gaelLDNL import createldgrid
except ImportError:  # package import
    from api.gaelLDNL import createldgrid
try:  # simple version
    from .version import __version__
except ImportError:  # package import
    from version import __version__

animate_toggle()  # CLOSE PRELOAD ANIMATION
# -- IMPORTS END --------------------------------------------------------------

# ################### START PROPERTIES/SETTINGS ############################# #
# GLOBALS (set in main before method calls)
exotic_infoDict = dict()
exotic_UIprevTPX, exotic_UIprevTPY, exotic_UIprevRPX, exotic_UIprevRPY = 0, 0, 0, 0
exotic_distFC = 0
exotic_ax1 = plt_exotic.figure()  # placeholder

# SETTINGS
plt_exotic.style.use(astropy_mpl_style)
# To increase memory allocation for EXOTIC; allows for more fits files
# resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
# ################### END PROPERTIES/SETTINGS ############################### #

# logging -- https://docs.python.org/3/library/logging.html
log = logging.getLogger(__name__)


def sigma_clip(ogdata, sigma=3, dt=21):
    nanmask = np.isnan(ogdata)
    mdata = savgol_filter(ogdata[~nanmask], dt, 2)
    res = ogdata[~nanmask] - mdata
    std = np.nanmedian([np.nanstd(np.random.choice(res, 50)) for i in range(100)])
    # std = np.nanstd(res) # biased from large outliers
    sigmask = np.abs(res) > sigma*std
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

    def planet_info(self):
        log.info(f"Looking up {self.planet}. Please wait. ...")
        self.planet, candidate = self._new_scrape(filename="eaConf.json", target=self.planet)

        if not candidate:
            with open("eaConf.json", "r") as confirmed:
                data = json.load(confirmed)
                planets = [data[i]['pl_name'] for i in range(len(data))]
                idx = planets.index(self.planet)
                self._get_params(data[idx])
                log.info('Successfully found {self.planet} in the NASA Exoplanet Archive!')

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

        uri_full = uri_full[:-1] + "&format={}".format(query.get("format", "csv"))
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
            target = input(f"Cannot find target ({target}) in NASA Exoplanet Archive. Check case sensitivity and spacing and"
                           "\nre-enter the planet's name or type candidate if this is a planet candidate: ")
            if target.strip().lower() == 'candidate':
                target = input("\nPlease enter candidate planet's name: ")
                log.debug("Candidate Planet's Name: %s " % target)
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
                rperr = np.sqrt(np.abs(data['pl_radjerr1']*data['pl_radjerr2'])) * R_JUP
                rs = data['st_rad'] * R_SUN
                rserr = np.sqrt(np.abs(data['st_raderr1']*data['st_raderr2'])) * R_SUN
                rprserr = ((rperr / rs) ** 2 + (-rp * rserr / rs ** 2) ** 2) ** 0.5
                rprs = rp / rs
        self.pl_dict = {
            'ra': data['ra'],
            'dec': data['dec'],
            'pName': data['pl_name'],
            'sName': data['hostname'],
            'pPer': data['pl_orbper'],
            'pPerUnc': np.sqrt(np.abs(data['pl_orbpererr1']*data['pl_orbpererr2'])),

            'midT': data['pl_tranmid'],
            'midTUnc': np.sqrt(np.abs(data['pl_tranmiderr1']*data['pl_tranmiderr2'])),
            'rprs': rprs,
            'rprsUnc': rprserr,
            'aRs': data['pl_ratdor'],
            'aRsUnc': np.sqrt(np.abs(data['pl_ratdorerr1']*data['pl_ratdorerr2'])),
            'inc': data['pl_orbincl'],
            'incUnc': np.sqrt(np.abs(data['pl_orbinclerr1']*data['pl_orbinclerr2'])),

            'ecc': data.get('pl_orbeccen', 0),
            'teff': data['st_teff'],
            'teffUncPos': data['st_tefferr1'],
            'teffUncNeg': data['st_tefferr2'],
            'met': data['st_met'],
            'metUncPos': data['st_meterr1'],
            'metUncNeg': data['st_meterr2'],
            'logg': data['st_logg'],
            'loggUncPos': data['st_loggerr1'],
            'loggUncNeg': data['st_loggerr2']
        }
# ################### END ARCHIVE SCRAPER (PRIORS) ############################


# Get Julian time, don't need to divide by 2 since assume mid-EXPOSURE
# Find separate funciton in code that does julian conversion to BJD_TDB
# Method that gets and returns the julian time of the observation
def getJulianTime(hdul):
    exptime_offset = 0
    imageheader = hdul[0].header

    exp = imageheader.get('EXPTIME')  # checking for variation in .fits header format
    if exp:
        exp = exp
    else:
        exp = imageheader.get('EXPOSURE')

    # Grab the BJD first
    if 'BJD_TDB' in hdul[0].header:
        julianTime = float(hdul[0].header['BJD_TDB'])
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in hdul[0].header.comments['BJD_TDB']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    elif 'BJD' in hdul[0].header:
        julianTime = float(hdul[0].header['BJD'])
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in hdul[0].header.comments['BJD']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    # then the DATE-OBS
    elif "UT-OBS" in hdul[0].header:
        gDateTime = hdul[0].header['UT-OBS']  # gets the gregorian date and time from the fits file header
        dt = dup.parse(gDateTime)
        atime = astropy.time.Time(dt)
        julianTime = atime.jd
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in hdul[0].header.comments['UT-OBS']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    # Then Julian Date
    elif 'JULIAN' in hdul[0].header:
        julianTime = float(hdul[0].header['JULIAN'])
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in hdul[0].header.comments['JULIAN']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    # Then MJD-OBS last, as in the MicroObservatory headers, it has less precision
    elif "MJD-OBS" in hdul[0].header:
        julianTime = float(hdul[0].header["MJD-OBS"]) + 2400000.5
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in hdul[0].header.comments['MJD-OBS']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now
    else:
        gDateTime = hdul[0].header['DATE-OBS']  # gets the gregorian date and time from the fits file header
        dt = dup.parse(gDateTime)
        atime = astropy.time.Time(dt)
        julianTime = atime.jd
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in hdul[0].header.comments['DATE-OBS']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now

    # If the mid-exposure time is given in the fits header, then no offset is needed to calculate the mid-exposure time
    return julianTime + exptime_offset


# Method that gets and returns the current phase of the target
def getPhase(curTime, pPeriod, tMid):
    phase = (curTime - tMid - 0.5*pPeriod) / pPeriod % 1
    return phase - 0.5


# Method that gets and returns the airmass from the fits file (Really the Altitude)
def getAirMass(hdul, ra, dec, lati, longit, elevation):
    # Grab airmass from image header; if not listed, calculate it from TELALT; if that isn't listed, then calculate it the hard way
    if 'AIRMASS' in hdul[0].header:
        am = float(hdul[0].header['AIRMASS'])
    elif 'TELALT' in hdul[0].header:
        alt = float(hdul[0].header['TELALT'])  # gets the airmass from the fits file header in (sec(z)) (Secant of the zenith angle)
        cosam = np.cos((PI / 180) * (90.0 - alt))
        am = 1 / cosam
    else:
        # pointing = SkyCoord(str(astropy.coordinates.Angle(raStr+" hours").deg)+" "+str(astropy.coordinates.Angle(decStr+" degrees").deg ), unit=(u.deg, u.deg), frame='icrs')
        pointing = SkyCoord(str(ra)+" "+str(dec), unit=(u.deg, u.deg), frame='icrs')

        location = EarthLocation.from_geodetic(lat=lati*u.deg, lon=longit*u.deg, height=elevation)
        atime = astropy.time.Time(getJulianTime(hdul), format='jd', scale='utc', location=location)
        pointingAltAz = pointing.transform_to(AltAz(obstime=atime, location=location))
        am = float(pointingAltAz.secz)
    return am


# Validate user input
def user_input(prompt, type_, val1=None, val2=None, val3=None, val4=None):
    while True:
        try:
            option = type_(input(prompt))
        except ValueError:
            print('Sorry, not a valid data type.')
            continue
        if type_ == str and val1 and val2:
            option = option.lower().replace(' ', '')
            if option not in (val1, val2):
                print("Sorry, your response was not valid.")
            else:
                return option
        elif type_ == int and val1 and val2 and val3:
            if option not in (val1, val2, val3):
                print("Sorry, your response was not valid.")
            else:
                return option
        elif type_ == int and val1 and val2 and val3 and val4:
            if option not in (val1, val2, val3, val4):
                print("Sorry, your response was not valid.")
            else:
                return option
        elif type_ == int and val1 and val2:
            if option not in (val1, val2):
                print("Sorry, your response was not valid.")
            else:
                return option
        else:
            return option


def get_save_directory(save_directory):
    while True:
        try:
            if save_directory == 'new':
                save_directory = create_directory()
            else:
                save_directory = os.path.join(save_directory, '')
                if not os.path.isdir(save_directory):
                    raise OSError
            return save_directory
        except OSError:
            print('Error: the directory entered does not exist. Please try again. Make sure to follow this formatting (using whichever directory you choose): /sample-data/results')
            save_directory = input("Enter the directory to save the results and plots into or type new to create one: ")

            log.debug(
                '\nError: the directory entered does not exist. Please try again. Make sure to follow this formatting (using whichever directory you choose): /sample-data/results')
            log.debug("Enter the directory to save the results and plots into or type new to create one: "+save_directory)


# Create a save directory within the current working directory
def create_directory():
    save_path = ""
    while True:
        try:
            directory_name = input('Enter the name for your new directory: ')
            log.debug('Enter the name for your new directory: '+directory_name)
            save_path = os.path.join(os.getcwd() + directory_name, '')
            os.mkdir(save_path)
        except OSError:
            log.info(f'Creation of the directory {save_path} failed.')
        else:
            log.info(f'Successfully created the directory {save_path}.')
            return save_path


# Check user's inits.json for user information and planetary parameters
def inits_file(inits_path, dict_info, dict_params):
    with open(inits_path) as json_file:
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
                       'pixel_scale': 'Pixel Scale (Ex: 5.21 arcsecs/pixel)'}

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
def get_initialization_file(infodict, userpdict):
    print(f"\nYour current working directory is: {os.getcwd()}")
    print(f"\nPotential initialization files I've found in {os.getcwd()} are: ")
    [print(i) for i in g.glob(os.getcwd() + "/*.json")]

    log.debug("Your current working directory is: " + os.getcwd())
    log.debug("Potential initialization files I've found in {} are: ".format(os.getcwd()))
    [log.debug("\t"+i) for i in g.glob(os.getcwd() + "/*.json")]

    while True:
        try:
            initfilename = str(input("\nPlease enter the Directory and Filename of your Initialization File: "))
            log.debug("Please enter the Directory and Filename of your Initialization File: "+initfilename)
            if initfilename == 'ok':
                initfilename = "/Users/rzellem/Documents/EXOTIC/inits.json"
            return inits_file(initfilename, infodict, userpdict)
        except FileNotFoundError:
            print("Error: Initialization file not found. Please try again.")
        except IsADirectoryError:
            print('Error: Entered a directory. Please try again.')


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
        if self.info['elev'] is None:
            self.elevation()
        if self.info['tarcoords'] is None:
            self.target_star_coords()
        if self.info['compstars'] is None:
            self.comparison_star_coords()
        return self.info, self.planet_name

    def image_directory(self):
        self.info['fitsdir'] = input('Please enter the Directory of Imaging Files: ')
        log.debug('Please enter the Directory of Imaging Files: '+self.info['fitsdir'])

    def save_directory(self):
        self.info['saveplot'] = input('Please enter the directory to save the results and plots into or type new to create one: ')
        log.debug(
            '\nPlease enter the directory to save the results and plots into or type new to create one: '+self.info['saveplot'])

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
        self.info['latitude'] = input("Please enter the longitude of where you observed (deg) "
                                      "(Don't forget the sign where East is '+' and West is '-'): ")
        log.debug("Please enter the longitude of where you observed (deg) "+self.info['latitude'])

    def longitude(self):
        self.info['longitude'] = input("Please enter the longitude of where you observed (deg) "
                                       "(Don't forget the sign where East is '+' and West is '-'): ")
        log.debug("Please enter the longitude of where you observed (deg) "+self.info['longitude'])

    def elevation(self):
        self.info['elev'] = user_input('Please enter the elevation (in meters) of where you observed: ', type_=float)
        log.debug('Please enter the elevation (in meters) of where you observed: '+str(self.info['elev']))

    def target_star_coords(self, pname):
        x_pix = user_input('\n{} X Pixel Coordinate: '.format(pname), type_=int)
        y_pix = user_input('\n{} Y Pixel Coordinate: '.format(pname), type_=int)

        log.debug(f'{pname} X Pixel Coordinate: {x_pix}')
        log.debug(f'{pname} Y Pixel Coordinate: {y_pix}')
        self.info['tarcoords'] = [x_pix, y_pix]

    def comparison_star_coords(self):
        num_comp_stars = user_input('How many comparison stars would you like to use? (1-10) ', type_=int)
        log.debug('How many comparison stars would you like to use? (1-10) '+str(num_comp_stars))
        comp_stars = []

        for num in range(num_comp_stars):
            x_pix = user_input('Comparison Star {} X Pixel Coordinate: '.format(num + 1), type_=int)
            y_pix = user_input('Comparison Star {} Y Pixel Coordinate: '.format(num + 1), type_=int)

            log.debug(f'Comparison Star {str((num + 1))} X Pixel Coordinate: {str(x_pix)}')
            log.debug(f'Comparison Star {str((num + 1))} Y Pixel Coordinate: {str(y_pix)}')
            comp_stars.append((x_pix, y_pix))
        self.info['compstars'] = comp_stars

    def exposure(self):
        self.info['exposure'] = user_input('Please enter your exposure time (seconds): ', type_=int)
        log.debug('Please enter your exposure time (seconds): '+str(self.info['exposure']))

    def pixel_scale(self):
        self.info['flatsdir'] = input('Please enter the size of your pixel (Ex: 5 arcsec/pixel): ')
        log.debug('Please enter the size of your pixel (Ex: 5 arcsec/pixel): '+str(self.info['flatsdir']))

    def planet(self):
        # self.planet_name = input('\nPlease enter the Planet Name. Make sure it matches the case sensitive name used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): ')
        # log.debug(
        #     '\nPlease enter the Planet Name. Make sure it matches the case sensitive name used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): '+self.planet_name)

        while True:
            self.planet_name = str(input(
                "\nPlease enter the Planet Name. Make sure it matches the case sensitive name and spacing used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): "))
            log.debug(
                "Please enter the Planet Name. Make sure it matches the case sensitive name and spacing used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): " +
                userpDict['pName'])

            if not userpDict['pName'][-2].isspace():
                print(
                    "The convention on the NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html) is to have a space between the star name and the planet letter. Please confirm that you have properly input the planet's name.")
                log.debug(
                    "The convention on the NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html) is to have a space between the star name and the planet letter. Please confirm that you have properly input the planet's name.")
                planetnameconfirm = user_input('\nPlease confirm:\n  (1) ' + userpDict[
                    'pName'] + ' is correct.\n  (2) The planet name needs to be changed.\nPlease select 1 or 2: ',
                                               type_=int, val1=1, val2=2)
                log.debug('\nPlease confirm:\n  (1) ' + userpDict[
                    'pName'] + ' is correct.\n  (2) The planet name needs to be changed.\nPlease select 1 or 2: ' + str(
                    planetnameconfirm))
            else:
                break
            if planetnameconfirm == 1:
                break


# Convert time units to BJD_TDB if pre-reduced file not in proper units
def timeConvert(timeList, timeFormat, pDict, info_dict):
    # if timeFormat is already in BJD_TDB, just do nothing
    # Perform appropriate conversion for each time format if needed
    if timeFormat == "JD_UTC":
        convertedTimes = utc_tdb.JDUTC_to_BJDTDB(timeList, ra=pDict['ra'], dec=pDict['dec'], lat=info_dict['lat'], longi=info_dict['long'], alt=info_dict['elev'])
        timeList = convertedTimes[0]
    elif timeFormat == "MJD_UTC":
        convertedTimes = utc_tdb.JDUTC_to_BJDTDB(timeList + 2400000.5, ra=pDict['ra'], dec=pDict['dec'], lat=info_dict['lat'], longi=info_dict['long'], alt=info_dict['elev'])
        timeList = convertedTimes[0]
    return timeList


# Convert magnitude units to flux if pre-reduced file not in flux already
def fluxConvert(fluxList, errorList, fluxFormat):
    # If units already in flux, do nothing, perform appropriate conversions to flux otherwise
    if fluxFormat == "magnitude":
        convertedPositiveErrors = 10. ** ((-1. * (fluxList + errorList)) / 2.5)
        convertedNegativeErrors = 10. ** ((-1. * (fluxList - errorList)) / 2.5)
        fluxList = 10. ** ((-1. * fluxList) / 2.5)
    if fluxFormat == "millimagnitude":
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
    for key, value in parameters.items():
        if value != init_parameters[key]:
            different = True
            break

    if different:
        opt = user_input('\nDifference(s) found between initialization file parameters and those scraped by EXOTIC from the NASA Exoplanet Archive. '
                         '\n Would you like:\n  (1) EXOTIC to adopt of all of your defined parameters or\n  (2) to review the ones scraped from the Archive that differ? \nPlease enter 1 or 2: ', type_=str, val1='1', val2='2')
        log.debug('Difference(s) found between initialization file parameters and those scraped by EXOTIC from the NASA Exoplanet Archive. '
                  '\n Would you like:\n  (1) EXOTIC to adopt of all of your defined parameters or\n  (2) to review the ones scraped from the Archive that differ? \nPlease enter 1 or 2: '+opt)

        if opt == '2':
            return True
        else:
            return False


# --------PLANETARY PARAMETERS UI------------------------------------------
# Get the user's confirmation of values that will later be used in lightcurve fit
def get_planetary_parameters(candplanetbool, userpdict, pdict=None):
    print('\n*******************************************')
    print("Planetary Parameters for Lightcurve Fitting\n")

    log.debug('*******************************************')
    log.debug("Planetary Parameters for Lightcurve Fitting")

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
        userpdict['ra'] = input('\nEnter the %s: ' % planet_params[0])
    if userpdict['dec'] is None:
        userpdict['dec'] = input('\nEnter the %s: ' % planet_params[1])
    userpdict['ra'], userpdict['dec'] = radec_hours_to_degree(userpdict['ra'], userpdict['dec'])

    radeclist = ['ra', 'dec']
    if not candplanetbool:
        for idx, item in enumerate(radeclist):
            uncert = 20 / 3600
            if pdict[item] - uncert <= userpdict[item] <= pdict[item] + uncert:
                continue
            else:
                print("\n\n*** WARNING: %s initialization file's %s does not match the value scraped by EXOTIC from the NASA Exoplanet Archive. ***\n" % (pdict['pName'], planet_params[idx]))
                print("\tNASA Exoplanet Archive value: %s" % pdict[item])
                print("\tInitialization file value: %s" % userpdict[item])
                print("\nWould you like to:\n  (1) use NASA Exoplanet Archive value, \n  (2) use initialization file value, or \n  (3) enter in a new value.")
                option = user_input('Which option do you choose? (1/2/3): ', type_=int, val1=1, val2=2, val3=3)

                log.debug("*** WARNING: %s initialization file's %s does not match the value scraped by EXOTIC from the NASA Exoplanet Archive. ***\n" % (pdict['pName'], planet_params[idx]))
                log.debug("\tNASA Exoplanet Archive value: %s" % pdict[item])
                log.debug("\tInitialization file value: %s" % userpdict[item])
                log.debug("Would you like to: \n  (1) use NASA Exoplanet Archive value, \n  (2) use initialization file value, or \n  (3) enter in a new value.")
                log.debug('Which option do you choose? (1/2/3): '+str(option))

                if option == 1:
                    userpdict[item] = pdict[item]
                elif option == 2:
                    continue
                else:
                    userpdict['ra'] = user_input('Enter the ' + planet_params[0] + ': ', type_=str)
                    userpdict['dec'] = user_input('Enter the ' + planet_params[1] + ': ', type_=str)
                    break

    if type(userpdict['ra']) and type(userpdict['dec']) is str:
        userpdict['ra'], userpdict['dec'] = radec_hours_to_degree(userpdict['ra'], userpdict['dec'])

    # Exoplanet confirmed in NASA Exoplanet Archive
    if not candplanetbool:
        log.debug('*** Here are the values scraped from the NASA Exoplanet Archive for %s that were not set (or set to null) in your initialization file. ***' %
                  pdict['pName'])
        # print('For each planetary parameter, enter "y" if you agree and "n" if you disagree.')
        # print('enter "1" to use NASA Exoplanet Archive value, "2" to use initialization file value, or "3" to enter a new value if you ')
        # print('decided to use an initialization file.')

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
                print("\n\n*** WARNING: %s initialization file's %s does not match the value scraped by EXOTIC from the NASA Exoplanet Archive. ***\n" % (pdict['pName'], planet_params[i]))
                print("\tNASA Exoplanet Archive value: %s" % pdict[key])
                print("\tInitialization file value: %s" % userpdict[key])
                print("\nWould you like to: \n  (1) use NASA Exoplanet Archive value, \n  (2) use initialization file value, or \n  (3) enter in a new value.")
                option = user_input('Which option do you choose? (1/2/3): ', type_=int, val1=1, val2=2, val3=3)

                log.debug("*** WARNING: %s initialization file's %s does not match the value scraped by EXOTIC from the NASA Exoplanet Archive. ***\n" % (pdict['pName'], planet_params[i]))
                log.debug("\tNASA Exoplanet Archive value: %s" % pdict[key])
                log.debug("\tInitialization file value: %s" % userpdict[key])
                log.debug("Would you like to: \n  (1) use NASA Exoplanet Archive value, \n  (2) use initialization file value, or \n  (3) enter in a new value.")
                log.debug('Which option do you choose? (1/2/3): '+str(option))
                if option == 1:
                    userpdict[key] = pdict[key]
                elif option == 2:
                    continue
                else:
                    userpdict[key] = user_input('Enter the ' + planet_params[i] + ': ', type_=type(userpdict[key]))
                    log.debug('Enter the ' + planet_params[i] + ': ' + str(userpdict[key]))
            # Did not use initialization file or null
            else:
                print('\n' + pdict['pName'] + ' ' + planet_params[i] + ': ' + str(pdict[key]))
                agreement = user_input('Do you agree? (y/n): ', type_=str, val1='y', val2='n')
                log.debug('Do you agree? (y/n): '+str(agreement))
                if agreement == 'y':
                    userpdict[key] = pdict[key]
                else:
                    userpdict[key] = user_input('Enter the ' + planet_params[i] + ': ', type_=type(pdict[key]))
                    log.debug('Enter the ' + planet_params[i] + ': '+str(userpdict[key]))

    # Exoplanet not confirmed in NASA Exoplanet Archive
    else:
        for i, key in enumerate(userpdict):
            if key in ('ra', 'dec'):
                continue
            # Used initialization file and is not empty
            if userpdict[key] is not None:
                agreement = user_input('%s: %s \nDo you agree? (y/n): '
                                       % (planet_params[i], userpdict[key]), type_=str, val1='y', val2='n')
                log.debug('%s: %s \nDo you agree? (y/n): '+agreement)
                if agreement == 'y':
                    continue
                else:
                    userpdict[key] = user_input('Enter the ' + planet_params[i] + ': ', type_=type(userpdict[key]))
                    log.debug('Enter the ' + planet_params[i] + ': '+userpdict[key])
            # Did not use initialization file
            else:
                if key in ('pName', 'sName'):
                    userpdict[key] = input('\nEnter the ' + planet_params[i] + ': ')
                else:
                    userpdict[key] = user_input('Enter the ' + planet_params[i] + ': ', type_=float)
                log.debug('Enter the ' + planet_params[i] + ': ' + str(userpdict[key]))
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
            print('Error: The format is not correct, please try again.')
            ra = input('Input the right ascension of target (HH:MM:SS): ')
            dec = input('Input the declination of target (<sign>DD:MM:SS): ')


def round_to_2(*args):
    x = args[0]
    if len(args) == 1:
        y = args[0]
    else:
        y = args[1]
    if np.floor(y) >= 1.:
        roundval = 2
    else:
        roundval = -int(np.floor(np.log10(abs(y))))+1
    return round(x, roundval)


# Check if user's directory contains imaging FITS files that are able to be reduced
def check_imaging_files(directory, filename):
    file_extensions = ['.fits', '.fit', '.fts']
    inputfiles = []

    while True:
        try:
            directory = os.path.join(directory, '')
            if os.path.isdir(directory):
                for exti in file_extensions:
                    for file in os.listdir(directory):
                        if file.lower().endswith(exti.lower()) and file[0:2] not in ('ref', 'wcs'):
                            inputfiles.append(os.path.join(directory, file))
                    if inputfiles:
                        return directory, inputfiles
                if not inputfiles:
                    raise FileNotFoundError
            else:
                raise OSError
        except FileNotFoundError:
            extaddoption = user_input("\nError: " + filename + " files not found with .fits, .fit or .fts extensions in " + directory +
                                      ".\nWould you like to enter in an alternate image extension in addition to .FITS? (y/n): ", type_=str, val1='y', val2='n')
            log.debug("Error: " + filename + " files not found with .fits, .fit or .fts extensions in " + directory +
                      ".\nWould you like to enter in an alternate image extension in addition to .FITS? (y/n): "+extaddoption)
            if extaddoption == 'y':
                addit_extension = input('Please enter the extension you want to add (EX: .FITS): ')
                file_extensions.append(addit_extension)
                log.debug('Please enter the extension you want to add (EX: .FITS): '+addit_extension)
            else:
                directory = input("Enter the directory path where " + filename + " files are located: ")
                log.debug("Enter the directory path where " + filename + " files are located: "+directory)
        except OSError:
            print("\nError: No such directory exists when searching for FITS files. Please try again.")
            directory = input("Enter the directory path where " + filename + " files are located: ")
            log.debug("Error: No such directory exists when searching for FITS files. Please try again.")
            log.debug("Enter the directory path where " + filename + " files are located: " + directory)


class LimbDarkening:

    def __init__(self, teff=None, teffpos=None, teffneg=None, met=None, metpos=None, metneg=None,
                 logg=None, loggpos=None, loggneg=None, wl_min=None, wl_max=None, filter_type=None):
        self.priors = {'T*': teff, 'T*_uperr': teffpos, 'T*_lowerr': teffneg,
                       'FEH*': met, 'FEH*_uperr': metpos, 'FEH*_lowerr': metneg,
                       'LOGG*': logg, 'LOGG*_uperr': loggpos, 'LOGG*_lowerr': loggneg}
        self.filter_type = filter_type
        self.wl_min = wl_min
        self.wl_max = wl_max
        self.ld0 = self.ld1 = self.ld2 = self.ld3 = None

        # Source for FWHM band wavelengths (units: nm): https://www.aavso.org/filters
        # Near-Infrared
        self.fwhm = {('J NIR 1.2micron', 'J'): (1040.00, 1360.00), ('H NIR 1.6micron', 'H'): (1420.00, 1780.00),
                     ('K NIR 2.2micron', 'K'): (2015.00, 2385.00),

                     # Sloan
                     ('Sloan u', 'SU'): (321.80, 386.80), ('Sloan g', 'SG'): (402.50, 551.50),
                     ('Sloan r', 'SR'): (553.10, 693.10), ('Sloan i', 'SI'): (697.50, 827.50),
                     ('Sloan z', 'SZ'): (841.20, 978.20),

                     # Stromgren
                     ('Stromgren b', 'STB'): (459.55, 478.05), ('Stromgren y', 'STY'): (536.70, 559.30),
                     ('Stromgren Hbw', 'STHBW'): (481.50, 496.50), ('Stromgren Hbn', 'STHBN'): (487.50, 484.50),

                     # Johnson
                     ('Johnson U', 'U'): (333.80, 398.80), ('Johnson B', 'B'): (391.60, 480.60),
                     ('Johnson V', 'V'): (502.80, 586.80), ('Johnson R', 'RJ'): (590.00, 810.00),
                     ('Johnson I', 'IJ'): (780.00, 1020.00),

                     # Cousins
                     ('Cousins R', 'R'): (561.70, 719.70), ('Cousins I', 'I'): (721.00, 875.00),

                     # MObs Clear Filter, Source: Martin Fowler
                     ('MObs CV', 'CV'): (350.00, 850.00),

                     # Astrodon CBB: George Silvis: https://astrodon.com/products/astrodon-exo-planet-filter/
                     ('Astrodon ExoPlanet-BB', 'CBB'): (500.00, 1000.00),

                     # LCO, Source: Kalee Tock & Michael Fitzgerald, https://lco.global/observatory/instruments/filters/
                     ('LCO Bessell B', 'N/A'): (391.60, 480.60), ('LCO Bessell V', 'N/A'): (502.80, 586.80),
                     ('LCO Pan-STARRS w', 'N/A'): (404.20, 845.80), ('LCO Pan-STARRS w', 'N/A'): (404.20, 845.80),
                     ('LCO Pan-STARRS zs', 'N/A'): (818.00, 922.00), ("LCO SDSS g'", 'N/A'): (402.00, 552.00),
                     ("LCO SDSS r'", 'N/A'): (552.00, 691.00), ("LCO SDSS i'", 'N/A'): (690.00, 819.00)}

    def nonlinear_ld(self):
        self._standard_list()

        if self.filter_type and not (self.wl_min or self.wl_max):
            self._standard()
        elif self.wl_min or self.wl_max:
            self._custom()
        else:
            option = user_input('\nWould you like EXOTIC to calculate your limb darkening parameters '
                                'with uncertainties? (y/n): ', type_=str, val1='y', val2='n')
            log.debug('Would you like EXOTIC to calculate your limb darkening parameters '
                      'with uncertainties? (y/n): '+str(option))

            if option == 'y':
                opt = user_input('Please enter 1 to use a standard filter or 2 for a customized filter: ',
                                 type_=int, val1=1, val2=2)
                log.debug('Please enter 1 to use a standard filter or 2 for a customized filter: '+str(opt))
                if opt == 1:
                    self._standard()
                elif opt == 2:
                    self._custom()
            else:
                self._user_entered()
        return self.ld0, self.ld1, self.ld2, self.ld3, self.filter_type

    def _standard_list(self):
        print('\n\n***************************')
        print('Limb Darkening Coefficients')
        print('***************************')
        log.debug('***************************'
                  '\nLimb Darkening Coefficients'
                  '\n***************************')
        print('\nStandard bands available to filter for limb darkening parameters (https://www.aavso.org/filters)'
              '\nas well as filters for MObs and LCO (0.4m telescope) datasets:\n')
        for key, value in self.fwhm.items():
            print('\t{}: {} - ({:.2f}-{:.2f}) nm'.format(key[1], key[0], value[0], value[1]))

    def _standard(self):
        while True:
            try:
                if not self.filter_type:
                    self.filter_type = input('\nPlease enter in the filter type (EX: Johnson V, V, STB, RJ): ')
                    log.debug('Please enter in the filter type (EX: Johnson V, V, STB, RJ): '+str(self.filter_type))
                for key, value in self.fwhm.items():
                    if self.filter_type in (key[0], key[1]) and self.filter_type != 'N/A':
                        self.filter_type = (key[0], key[1])
                        break
                else:
                    raise KeyError
                break
            except KeyError:
                print('\nError: The entered filter is not in the provided list of standard filters.')
                log.debug('Error: The entered filter is not in the provided list of standard filters.')
                self.filter_type = None

        self.wl_min = self.fwhm[self.filter_type][0]
        self.wl_max = self.fwhm[self.filter_type][1]
        self.filter_type = self.filter_type[1]
        self._calculate_ld()

    def _custom(self):
        self.filter_type = 'N/A'
        if not self.wl_min:
            self.wl_min = float(input('FWHM Minimum wavelength (nm): '))
            log.debug('FWHM Minimum wavelength (nm): '+str(self.wl_min))
        if not self.wl_max:
            self.wl_max = float(input('FWHM Maximum wavelength (nm): '))
            log.debug('FWHM Maximum wavelength (nm): ' + str(self.wl_max))
        self._calculate_ld()

    def _user_entered(self):
        self.filter_type = input('\nEnter in your filter name: ')
        ld_0 = user_input('\nEnter in your first nonlinear term: ', type_=float)
        ld0_unc = user_input('Enter in your first nonlinear term uncertainty: ', type_=float)
        ld_1 = user_input('\nEnter in your second nonlinear term: ', type_=float)
        ld1_unc = user_input('Enter in your second nonlinear term uncertainty: ', type_=float)
        ld_2 = user_input('\nEnter in your third nonlinear term: ', type_=float)
        ld2_unc = user_input('Enter in your third nonlinear term uncertainty: ', type_=float)
        ld_3 = user_input('\nEenter in your fourth nonlinear term: ', type_=float)
        ld3_unc = user_input('Enter in your fourth nonlinear term uncertainty: ', type_=float)
        self.ld0, self.ld1, self.ld2, self.ld3 = (ld_0, ld0_unc), (ld_1, ld1_unc), (ld_2, ld2_unc), (ld_3, ld3_unc)

        log.debug('Filter name: ' + str(self.filter_type))
        log.debug(f"User-defined nonlinear limb-darkening coefficients: {ld_0}+/-{ld0_unc}, {ld_1}+/-{ld1_unc}, {ld_2}+/-{ld2_unc}, {ld_3}+/-{ld3_unc}")

    def _calculate_ld(self):
        self.wl_min = [self.wl_min / 1000]
        self.wl_max = [self.wl_max / 1000]
        ld_params = createldgrid(np.array(self.wl_min), np.array(self.wl_max), self.priors)
        self.ld0 = ld_params['LD'][0][0], ld_params['ERR'][0][0]
        self.ld1 = ld_params['LD'][1][0], ld_params['ERR'][1][0]
        self.ld2 = ld_params['LD'][2][0], ld_params['ERR'][2][0]
        self.ld3 = ld_params['LD'][3][0], ld_params['ERR'][3][0]

        log.debug("EXOTIC-calculated nonlinear limb-darkening coefficients: ")
        log.debug(""+str(ld_params['LD'][0][0]) + " +/- " + str(ld_params['ERR'][0][0]))
        log.debug(""+str(ld_params['LD'][1][0]) + " +/- " + str(ld_params['ERR'][1][0]))
        log.debug(""+str(ld_params['LD'][2][0]) + " +/- " + str(ld_params['ERR'][2][0]))
        log.debug(""+str(ld_params['LD'][3][0]) + " +/- " + str(ld_params['ERR'][3][0]))


def check_wcs(fits_file, save_directory, plate_opt):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=FITSFixedWarning)
        hdulist = fits.open(name=fits_file, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
        header = hdulist[0].header
        wcsheader = WCS(header)
        wcs_exists = wcsheader.is_celestial
        hdulist.close()
        del hdulist

    if plate_opt == 'y':
        return get_wcs(fits_file, save_directory)
    elif plate_opt == 'n':
        if wcs_exists:
            trust_wcs = user_input('The imaging data from your file has WCS information. Do you trust this? (y/n): ',
                                   type_=str, val1='y', val2='n')
            log.debug('The imaging data from your file has WCS information. Do you trust this? (y/n): '+trust_wcs)

            if trust_wcs == 'y':
                return fits_file
            else:
                return False
    else:
        if wcs_exists:
            trust_wcs = user_input('The imaging data from your file has WCS information. Do you trust this? (y/n): ',
                                   type_=str, val1='y', val2='n')
            log.debug(
                '\nThe imaging data from your file has WCS information. Do you trust this? (y/n): ' + trust_wcs)
        else:
            trust_wcs = 'n'

        if trust_wcs == 'n':
            opt = user_input("\nWould you like to upload the your image for a plate solution?"
                             "\nDISCLAIMER: One of your imaging files will be publicly viewable on nova.astrometry.net."
                             " (y/n): ", type_=str, val1='y', val2='n')
            log.debug("Would you like to upload the your image for a plate solution?"
                      "\nDISCLAIMER: One of your imaging files will be publicly viewable on nova.astrometry.net."
                      " (y/n): "+opt)
            if opt == 'y':
                return get_wcs(fits_file, save_directory)
            else:
                return False
        else:
            return fits_file


def get_wcs(file, directory):
    log.info("")
    log.info("Getting the plate solution for your imaging file. Please wait. ...")
    animate_toggle(True)
    wcs_obj = PlateSolution(file=file, directory=directory)
    wcs_file = wcs_obj.plate_solution()
    animate_toggle()
    return wcs_file


def is_false(value):
    return value is False


def result_if_max_retry_count(retry_state):
    pass


class PlateSolution:
    default_url = 'http://nova.astrometry.net/api/'
    default_apikey = {'apikey': 'vfsyxlmdxfryhprq'}

    def __init__(self, apiurl=default_url, apikey=default_apikey, file=None, directory=None):
        self.apiurl = apiurl
        self.apikey = apikey
        self.file = file
        self.directory = directory

    def plate_solution(self):
        session = self._login()
        if not session:
            return PlateSolution.fail('Login')

        sub_id = self._upload(session)
        if not sub_id:
            return PlateSolution.fail('Upload')

        sub_url = self._get_url('submissions/%s' % sub_id)
        job_id = self._sub_status(sub_url)
        if not job_id:
            return PlateSolution.fail('Submission ID')

        job_url = self._get_url('jobs/%s' % job_id)
        download_url = self.apiurl.replace('/api/', '/new_fits_file/%s/' % job_id)
        wcs_file = os.path.join(self.directory, 'wcs_image.fits')
        wcs_file = self._job_status(job_url, wcs_file, download_url)
        if not wcs_file:
            return PlateSolution.fail('Job Status')
        else:
            print('\n\nSuccess. ')
            log.debug("WCS file creation successful.")
            return wcs_file

    def _get_url(self, service):
        return self.apiurl + service

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(ConnectionError) |
                  retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _login(self):
        r = requests.post(self._get_url('login'), data={'request-json': json.dumps(self.apikey)})
        if r.json()['status'] == 'success':
            return r.json()['session']
        return False

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(ConnectionError) |
                  retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _upload(self, session):
        files = {'file': open(self.file, 'rb')}
        headers = {'request-json': json.dumps({"session": session}), 'allow_commercial_use': 'n',
                   'allow_modifications': 'n', 'publicly_visible': 'n'}

        r = requests.post(self.apiurl + 'upload', files=files, data=headers)

        if r.json()['status'] == 'success':
            return r.json()['subid']
        return False

    @retry(stop=stop_after_attempt(45), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(ConnectionError) |
                  retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _sub_status(self, sub_url):
        r = requests.get(sub_url)
        if r.json()['job_calibrations']:
            return r.json()['jobs'][0]
        return False

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(ConnectionError) |
                  retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _job_status(self, job_url, wcs_file, download_url):
        r = requests.get(job_url)
        if r.json()['status'] == 'success':
            r = requests.get(download_url)
            with open(wcs_file, 'wb') as f:
                f.write(r.content)
            return wcs_file
        return False

    @staticmethod
    def fail(error_type):
        log.info("")
        log.info('WARNING: After multiple attempts, EXOTIC could not retrieve a plate solution from nova.astrometry.net'
                 f' due to {error_type}. EXOTIC will continue reducing data without a plate solution.')
        return False


# Getting the right ascension and declination for every pixel in imaging file if there is a plate solution
def get_radec(hdulWCSrd):
    wcsheader = WCS(hdulWCSrd[0].header)
    xaxis = np.arange(hdulWCSrd[0].header['NAXIS1'])
    yaxis = np.arange(hdulWCSrd[0].header['NAXIS2'])
    x, y = np.meshgrid(xaxis, yaxis)
    return wcsheader.all_pix2world(x, y, 0)


# Check the ra and dec against the plate solution to see if the user entered in the correct values
def check_targetpixelwcs(pixx, pixy, expra, expdec, ralist, declist):
    while True:
        try:
            uncert = 20 / 3600
            # Margins are within 20 arcseconds
            if expra - uncert >= ralist[pixy][pixx] or ralist[pixy][pixx] >= expra + uncert:
                print('\nError: The X Pixel Coordinate entered does not match the right ascension.')
                log.debug('Error: The X Pixel Coordinate entered does not match the right ascension.')
                raise ValueError
            if expdec - uncert >= declist[pixy][pixx] or declist[pixy][pixx] >= expdec + uncert:
                print('\nError: The Y Pixel Coordinate entered does not match the declination.')
                log.debug('Error: The Y Pixel Coordinate entered does not match the declination.')
                raise ValueError
            return pixx, pixy
        except ValueError:
            repixopt = user_input('Would you like to re-enter the pixel coordinates? (y/n): ', type_=str, val1='y', val2='n')
            log.debug('Would you like to re-enter the pixel coordinates? (y/n): '+repixopt)

            # User wants to change their coordinates
            if repixopt == 'y':
                # Checks for the closest pixel location in ralist and declist for expected ra and dec
                dist = (ralist - expra) ** 2 + (declist - expdec) ** 2
                pixy, pixx = np.unravel_index(dist.argmin(), dist.shape)
                searchopt = user_input('Here are the suggested pixel coordinates: X Pixel: %s Y Pixel: %s'
                                       '\nWould you like to use these? (y/n): ' % (pixx, pixy), type_=str, val1='y', val2='n')
                log.debug('Here are the suggested pixel coordinates: X Pixel: %s Y Pixel: %s'
                          '\nWould you like to use these? (y/n): %s' % (pixx, pixy, searchopt))
                # Use the coordinates found by code
                if searchopt == 'y':
                    return pixx, pixy
                # User enters their own coordinates to be re-checked
                else:
                    pixx = user_input("Please re-enter the target star's X Pixel Coordinate: ", type_=int)
                    pixy = user_input("Please re-enter the target star's Y Pixel Coordinate: ", type_=int)

                    log.debug("Please re-enter the target star's X Pixel Coordinate: %s" % pixx)
                    log.debug("Please re-enter the target star's Y Pixel Coordinate: %s" % pixy)
            else:
                # User does not want to change coordinates even though they don't match the expected ra and dec
                return pixx, pixy


# Checks if comparison star is variable via querying SIMBAD
def variableStarCheck(refx, refy, hdulWCS):
    # Read in WCS data from plate solution file and convert comparison star coordinates from pixel to WCS
    w = WCS(hdulWCS[0].header)
    world = w.wcs_pix2world(np.array([[refx, refy]], dtype=np.float64), 1)
    ra = world[0][0]
    dec = world[0][1]
    sample = SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')

    # Query GAIA first to check for variability using the phot_variable_flag trait
    radius = u.Quantity(20.0, u.arcsec)
    gaiaQuery = Gaia.cone_search_async(sample, radius)
    gaiaResult = gaiaQuery.get_results()

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
    simbad_result = Simbad.query_region(sample, radius=20*u.arcsec)
    try:
        starName = simbad_result['MAIN_ID'][0].decode("utf-8")
    except:
        print("Your star cannot be resolved in SIMBAD. Proceed with caution.")
        log.debug("Your star cannot be resolved in SIMBAD. Proceed with caution.")
        return False
    identifiers = Simbad.query_objectids(starName)
    for currName in identifiers:
        if "V*" in currName[0]:
            return True
    return False


# Aligns imaging data from .fits file to easily track the host and comparison star's positions
def image_alignment(imagedata):
    print("\nAligning your images from FITS files. Please wait.")
    boollist = []
    notaligned = 0

    # Align images from .FITS files and catch exceptions if images can't be aligned.
    # boollist for discarded images to delete .FITS data from airmass and times.
    for i, image_file in enumerate(imagedata):
        try:
            sys.stdout.write('Aligning Image %s of %s\r' % (str(i+1), str(len(imagedata))))
            sys.stdout.flush()
            imagedata[i], footprint = aa.register(image_file, imagedata[0])
            boollist.append(True)
        except:
            print('Image %s of %s failed to align' % (str(i+1), str(len(imagedata))))
            notaligned += 1
            boollist.append(False)

    imagedata = imagedata[boollist]

    if notaligned > 0:
        print('\n\n*********************************************************************')
        print('WARNING: From the given imaging files: %s of %s were not aligned.'
              % (str(notaligned), str(len(imagedata) + notaligned)))
        print('*********************************************************************')

        log.debug('*********************************************************************')
        log.debug('WARNING: From the given imaging files: %s of %s were not aligned.' %
                  (str(notaligned), str(len(imagedata) + notaligned)))
        log.debug('*********************************************************************')
        time.sleep(5)

    return imagedata, boollist


def get_pixel_scale(hdul, file, header, pixel_init):
    if file:
        for i in range(len(hdul[0].header)):
            if hdul[0].header['COMMENT'][i].split(' ')[0] == 'scale:':
                imscalen = float(hdul[0].header['COMMENT'][i].split(' ')[1])
                break
        imagescale = 'Image scale in arc-secs/pixel: {:.2f}'.format(imscalen)
    elif "IM_SCALE" in header:
        imscalen = header['IM_SCALE']
        imscaleunits = header.comments['IM_SCALE']
        imagescale = imscaleunits + ": " + str(imscalen)
    elif "PIXSCALE" in header:
        imscalen = header['PIXSCALE']
        imscaleunits = header.comments['PIXSCALE']
        imagescale = imscaleunits + ": " + str(imscalen)
    elif pixel_init:
        imagescale = "Image scale: " + str(pixel_init)
    else:
        print("Cannot find the pixel scale in the image header.")
        imscalen = input("Please enter the size of your pixel (e.g., 5 arc-sec/pixel). ")
        imagescale = "Image scale: " + str(imscalen)

        log.debug("Cannot find the pixel scale in the image header.")
        log.debug("Please enter the size of your pixel (e.g., 5 arc-sec/pixel). "+imscalen)
    return imagescale


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
    rx = (x-x0)*np.cos(rot) - (y-y0)*np.sin(rot)
    ry = (x-x0)*np.sin(rot) + (y-y0)*np.cos(rot)
    gausx = np.exp(-rx**2 / (2*sigx**2))
    gausy = np.exp(-ry**2 / (2*sigy**2))
    return a*gausx*gausy + b


def fit_psf(data, pos, init, lo, up, psf_function=gaussian_psf, lossfn='linear', method='trf', box=15):
    xv, yv = mesh_box(pos, box)

    def fcn2min(pars):
        model = psf_function(xv, yv, *pars)
        return (data[yv, xv]-model).flatten()

    if method == 'trf':
        res = least_squares(fcn2min, x0=[*pos, *init], bounds=[lo, up], loss=lossfn, jac='3-point', method=method)
    else:
        res = least_squares(fcn2min, x0=[*pos, *init], loss=lossfn, jac='3-point', method=method)
    return res.x


def mesh_box(pos, box):
    pos = [int(np.round(pos[0])), int(np.round(pos[1]))]
    x = np.arange(pos[0]-box, pos[0]+box+1)
    y = np.arange(pos[1]-box, pos[1]+box+1)
    xv, yv = np.meshgrid(x, y)
    return xv.astype(int), yv.astype(int)


# Method fits a 2D gaussian function that matches the star_psf to the star image and returns its pixel coordinates
def fit_centroid(data, pos, init=None, box=10):

    # get sub field in image
    xv, yv = mesh_box(pos, box)
    wx, wy = pos  # could take flux weighted centroid if not crowded

    if init:
        pass
    else:
        init = [np.nanmax(data[yv, xv])-np.nanmin(data[yv, xv]), 1, 1, 0, np.nanmin(data[yv, xv])]

    try:
        # fit gaussian PSF
        pars = fit_psf(
            data,
            [wx, wy],  # position estimate
            init,    # initial guess: [amp, sigx, sigy, rotation, bg]
            [wx-5, wy-5, 0, 0, 0, -PI/4, np.nanmin(data)-1],  # lower bound: [xc, yc, amp, sigx, sigy, rotation,  bg]
            [wx+5, wy+5, 1e7, 20, 20, PI/4, np.nanmax(data[yv, xv])+1],  # upper bound
            psf_function=gaussian_psf, method='trf',
            box=box  # only fit a subregion +/- 5 px from centroid
        )
    except:
        print("WARNING trouble fitting Gaussian PSF to star at {},{}".format(wx, wy))
        print("  check location of comparison star in the first few images")
        print("  fitting parameters are out of bounds")
        print("  init:", init)
        print(" lower:", [wx-5, wy-5, 0, 0, 0, -PI/4, np.nanmin(data)-1])
        print(" upper:", [wx+5, wy+5, 1e7, 20, 20, PI/4, np.nanmax(data[yv, xv])+1])

        log.debug("WARNING trouble fitting Gaussian PSF to star at {},{}".format(wx, wy))
        log.debug("  check location of comparison star in the first few images")
        log.debug("  fitting parameters are out of bounds")
        log.debug("  init:", init)
        log.debug(" lower:", [wx-5, wy-5, 0, 0, 0, -np.pi/4, np.nanmin(data)-1])
        log.debug(" upper:", [wx+5, wy+5, 1e7, 20, 20, np.pi/4, np.nanmax(data[yv, xv])+1])

        # use LM in unbounded optimization
        pars = fit_psf(
            data, [wx, wy], init,
            [wx-5, wy-5, 0, 0, 0, -PI/4, np.nanmin(data)-1],
            [wx+5, wy+5, 1e7, 20, 20, PI/4, np.nanmax(data[yv, xv])+1],
            psf_function=gaussian_psf,
            box=box, method='lm'
        )
    return pars


# Method calculates the flux of the star (uses the skybg_phot method to do backgorund sub)
def getFlux(data, xc, yc, r=5, dr=5):

    if dr > 0:
        bgflux = skybg_phot(data, xc, yc, r+2, dr)
    else:
        bgflux = 0
    positions = [(xc, yc)]
    bdata = data-bgflux

    apertures = CircularAperture(positions, r=r)
    phot_table = aperture_photometry(bdata, apertures, method='exact')

    return float(phot_table['aperture_sum']), bgflux


def skybg_phot(data, xc, yc, r=10, dr=5, ptol=85, debug=False):
    # create a crude annulus to mask out bright background pixels
    xv, yv = mesh_box([xc, yc], np.round(r+dr))
    rv = ((xv-xc)**2 + (yv-yc)**2)**0.5
    mask = (rv > r) & (rv < (r+dr))
    cutoff = np.nanpercentile(data[yv, xv][mask], ptol)
    dat = np.array(data[yv, xv], dtype=float)
    dat[dat > cutoff] = np.nan  # ignore pixels brighter than percentile

    if debug:
        minb = data[yv, xv][mask].min()
        maxb = data[yv, xv][mask].mean() + 3 * data[yv, xv][mask].std()
        nanmask = np.nan*np.zeros(mask.shape)
        nanmask[mask] = 1
        bgsky = data[yv, xv]*nanmask
        cmode = mode(dat.flatten(), nan_policy='omit').mode[0]
        amode = mode(bgsky.flatten(), nan_policy='omit').mode[0]

        fig, ax = plt_exotic.subplots(2, 2, figsize=(9, 9))
        im = ax[0, 0].imshow(data[yv, xv], vmin=minb, vmax=maxb, cmap='inferno')
        ax[0, 0].set_title("Original Data")
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax[0, 0])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        ax[1, 0].hist(bgsky.flatten(), label='Sky Annulus ({:.1f}, {:.1f})'.format(np.nanmedian(bgsky), amode), alpha=0.5, bins=np.arange(minb, maxb))
        ax[1, 0].hist(dat.flatten(), label='Clipped ({:.1f}, {:.1f})'.format(np.nanmedian(dat), cmode), alpha=0.5, bins=np.arange(minb, maxb))
        ax[1, 0].legend(loc='best')
        ax[1, 0].set_title("Sky Background")
        ax[1, 0].set_xlabel("Pixel Value")

        ax[1, 1].imshow(dat, vmin=minb, vmax=maxb, cmap='inferno')
        ax[1, 1].set_title("Clipped Sky Background")

        ax[0, 1].imshow(bgsky, vmin=minb, vmax=maxb, cmap='inferno')
        ax[0, 1].set_title("Sky Annulus")
        plt_exotic.tight_layout()
        plt_exotic.show()
    return min(np.nanmedian(dat), mode(dat.flatten(), nan_policy='omit').mode[0])


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
    plt_exotic.figure()
    plt_exotic.plot(times - np.nanmin(times), xTarg, '-bo')
    plt_exotic.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt_exotic.ylabel('X Pixel Position')
    plt_exotic.title(targetname + ' X Centroid Position ' + date)
    plt_exotic.savefig(exotic_infoDict['saveplot'] + 'temp/XCentroidPosition' + targetname + date + '.png')
    plt_exotic.close()

    # Y TARGET
    plt_exotic.figure()
    plt_exotic.plot(times - np.nanmin(times), yTarg, '-bo')
    plt_exotic.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt_exotic.ylabel('Y Pixel Position')
    plt_exotic.title(targetname + ' Y Centroid Position ' + date)
    plt_exotic.savefig(exotic_infoDict['saveplot'] + 'temp/YCentroidPos' + targetname + date + '.png')
    plt_exotic.close()

    # X COMP
    plt_exotic.figure()
    plt_exotic.plot(times - np.nanmin(times), xRef, '-ro')
    plt_exotic.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt_exotic.ylabel('X Pixel Position')
    plt_exotic.title('Comp Star X Centroid Position ' + date)
    plt_exotic.savefig(exotic_infoDict['saveplot'] + 'temp/CompStarXCentroidPos' + date + '.png')
    plt_exotic.close()

    # Y COMP
    plt_exotic.figure()
    plt_exotic.plot(times - np.nanmin(times), yRef, '-ro')
    plt_exotic.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt_exotic.ylabel('Y Pixel Position')
    plt_exotic.title('Comp Star Y Centroid Position ' + date)
    plt_exotic.savefig(exotic_infoDict['saveplot'] + 'temp/CompStarYCentroidPos' + date + '.png')
    plt_exotic.close()

    # X DISTANCE BETWEEN TARGET AND COMP
    plt_exotic.figure()
    plt_exotic.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt_exotic.ylabel('X Pixel Distance')
    for e in range(0, len(xTarg)):
        plt_exotic.plot(times[e] - np.nanmin(times), abs(int(xTarg[e]) - int(xRef[e])), 'bo')
    plt_exotic.title('Distance between Target and Comparison X position')
    plt_exotic.savefig(exotic_infoDict['saveplot'] + 'temp/XCentroidDistance' + targetname + date + '.png')
    plt_exotic.close()

    # Y DISTANCE BETWEEN TARGET AND COMP
    plt_exotic.figure()
    plt_exotic.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt_exotic.ylabel('Y Pixel Difference')
    d = 0
    for d in range(0, len(yTarg)):
        plt_exotic.plot(times[d] - np.nanmin(times), abs(int(yTarg[d]) - int(yRef[d])), 'bo')
    plt_exotic.title('Difference between Target and Comparison Y position')
    plt_exotic.savefig(exotic_infoDict['saveplot'] + 'temp/YCentroidDistance' + targetname + date + '.png')
    plt_exotic.close()


def realTimeReduce(i, target_name):
    targetFluxVals = []
    referenceFluxVals = []
    normalizedFluxVals = []
    fileNameList = []
    timeSortedNames = []
    timeList = []
    timesListed = []

    # -------TIME SORT THE FILES--------------------------------------------------------------------------------
    directoryP = ""
    directToWatch = str(input("Enter the Directory Path where .FITS or .FTS Image Files are located: "))
    # Add / to end of directory if user does not input it
    if directToWatch[-1] != "/":
        directToWatch += "/"
    directoryP = directToWatch

    while len(g.glob(directoryP)) == 0:
        print("Error: .FITS files not found in " + directoryP)
        directToWatch = str(input("Enter the Directory Path where FITS Image Files are located: "))
        # Add / to end of directory if user does not input it
        if directToWatch[-1] != "/":
            directToWatch += "/"
        directoryP = directToWatch

    fileNumber = 1
    for fileName in g.glob(directoryP):  # Loop through all the fits files and time sorts
        fitsHead = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
        # TIME
        timeVal = getJulianTime(fitsHead)  # gets the julian time registered in the fits header
        timeList.append(timeVal)  # adds to time value list
        fileNameList.append(fileName)
        fitsHead.close()  # close stream
        del fitsHead

    # Time sorts the file names based on the fits file header
    timeSortedNames = [x for _, x in sorted(zip(timeList, fileNameList))]

    # sorts the times for later plotting use
    # sortedTimeList = sorted(timeList)

    # hdul = fits.open(name=timeSortedNames[0], memmap=False, cache=False, lazy_load_hdus=False)  # opens the fits file
    # Extracts data from the image file and puts it in a 2D numpy array: firstImageData
    firstImageData = fits.getdata(timeSortedNames[0], ext=0)

    # fit first image
    targx, targy, targamplitude, targsigX, targsigY, targrot, targoff = fit_centroid(firstImageData, [exotic_UIprevTPX, exotic_UIprevTPY], box=10)
    refx, refy, refamplitude, refsigX, refsigY, refrot, refoff = fit_centroid(firstImageData, [exotic_UIprevRPX, exotic_UIprevRPY], box=10)

    # just use one aperture and annulus
    apertureR = 3 * max(targsigX, targsigY)
    annulusR = 10

    for imageFile in timeSortedNames:

        hdul = fits.open(name=imageFile, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
        # Extracts data from the image file and puts it in a 2D numpy array: imageData
        currTime = getJulianTime(hdul)
        hdul[0].data.scale('float32')
        imageData = hdul['ext', 0].data  # fits.getdata(imageFile, ext=0)
        header = hdul[0].header  # fits.getheader(imageFile)

        hdul.close()  # close the stream
        del hdul

        # Find the target star in the image and get its pixel coordinates if it is the first file
        if fileNumber == 1:
            # Initializing the star location guess as the user inputted pixel coordinates
            prevTPX, prevTPY, prevRPX, prevRPY = exotic_UIprevTPX, exotic_UIprevTPY, exotic_UIprevRPX, exotic_UIprevRPY
            prevTSigX, prevTSigY, prevRSigX, prevRSigY = targsigX, targsigY, refsigX, refsigY

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

        txmin = int(prevTPX) - exotic_distFC  # left
        txmax = int(prevTPX) + exotic_distFC  # right
        tymin = int(prevTPY) - exotic_distFC  # top
        tymax = int(prevTPY) + exotic_distFC  # bottom

        targSearchA = imageData[tymin:tymax, txmin:txmax]

        # Set reference search area
        rxmin = int(prevRPX) - exotic_distFC  # left
        rxmax = int(prevRPX) + exotic_distFC  # right
        rymin = int(prevRPY) - exotic_distFC  # top
        rymax = int(prevRPY) + exotic_distFC  # bottom

        refSearchA = imageData[rymin:rymax, rxmin:rxmax]

        # Guess at Gaussian Parameters and feed them in to help gaussian fitter

        tGuessAmp = targSearchA.max() - targSearchA.min()

        # Fits Centroid for Target
        myPriors = [tGuessAmp, prevTSigX, prevTSigY, 0, targSearchA.min()]
        tx, ty, tamplitude, tsigX, tsigY, trot, toff = fit_centroid(imageData, [prevTPX, prevTPY], init=myPriors, box=10)
        tpsfFlux = 2*PI*tamplitude*tsigX*tsigY
        currTPX = tx
        currTPY = ty

        # Fits Centroid for Reference
        rGuessAmp = refSearchA.max() - refSearchA.min()
        myRefPriors = [rGuessAmp, prevRSigX, prevRSigY, 0, refSearchA.min()]
        rx, ry, ramplitude, rsigX, rsigY, rrot, roff = fit_centroid(imageData, [prevRPX, prevRPY], init=myRefPriors, box=10)
        rpsfFlux = 2*PI*ramplitude*rsigX*rsigY
        currRPX = rx
        currRPY = ry

        # gets the flux value of the target star and
        tFluxVal, tTotCts = getFlux(imageData, currTPX, currTPY, apertureR, annulusR)
        targetFluxVals.append(tFluxVal)  # adds tFluxVal to the total list of flux values of target star

        # gets the flux value of the reference star and subracts the background light
        rFluxVal, rTotCts = getFlux(imageData, currRPX, currRPY, apertureR, annulusR)
        referenceFluxVals.append(rFluxVal)  # adds rFluxVal to the total list of flux values of reference star

        normalizedFluxVals.append((tFluxVal / rFluxVal))

        # TIME
        timesListed.append(currTime)

        # UPDATE PIXEL COORDINATES and SIGMAS
        # target
        prevTPX = currTPX
        prevTPY = currTPY
        prevTSigX = tsigX
        prevTSigY = tsigY
        # reference
        prevRPX = currRPX
        prevRPY = currRPY
        prevRSigX = rsigX
        prevTSigY = rsigY

        # UPDATE FILE COUNT
        prevImageData = imageData
        fileNumber = fileNumber + 1

    # EXIT THE FILE LOOP

    exotic_ax1.clear()
    exotic_ax1.set_title(target_name)
    exotic_ax1.set_ylabel('Normalized Flux')
    exotic_ax1.set_xlabel('Time (jd)')
    exotic_ax1.plot(timesListed, normalizedFluxVals, 'bo')


def parse_args():
    parser = argparse.ArgumentParser()

    help_ = "choose an inits file"
    parser.add_argument("-i", "--init", help=help_, type=str, default="")

    return parser.parse_args()


def main():

    # command line args
    args = parse_args()

    log.debug("*************************")
    log.debug("EXOTIC reduction log file")
    log.debug("*************************")
    log.debug("Starting ...")
    log.info("")
    log.info(f"Python Version: {str(sys.version)}")

    log.info("*************************************************************")
    log.info("Welcome to the EXOplanet Transit Interpretation Code (EXOTIC)")
    log.info(f"Version {__version__}")
    log.info("*************************************************************")

    # ---INITIALIZATION-------------------------------------------------------
    global exotic_infoDict
    global exotic_UIprevTPX, exotic_UIprevTPY, exotic_UIprevRPX, exotic_UIprevRPY
    global exotic_distFC
    global exotic_ax1

    targetFluxVals, referenceFluxVals, normalizedFluxVals, targUncertanties, refUncertanties, timeList, phasesList, airMassList = (
        [] for l in range(8))

    fileNameList, timeSortedNames, xTargCent, yTargCent, xRefCent, yRefCent, finXTargCent, finYTargCent, finXRefCent, finYRefCent = (
        [] for m in range(10))

    timesListed = []  # sorted times of observation
    fileNumber = 1  # initializes file number to one
    minSTD = 100000  # sets the initial minimum standard deviation absurdly high so it can be replaced immediately
    minChi2 = 100000
    exotic_distFC = 10  # gaussian search area
    context = {}

    # ---USER INPUTS--------------------------------------------------------------------------
    realTimeAns = user_input('Enter "1" for Real Time Reduction or "2" for for Complete Reduction: ', type_=int, val1=1, val2=2)

    #############################
    # Real Time Reduction Routine
    #############################

    if realTimeAns == 1:
        log.info('**************************************************************')
        log.info('Real Time Reduction ("Control + C"  or close the plot to quit)')
        log.info('**************************************************************')

        directToWatch = str(input("Enter the Directory Path of imaging files: "))
        log.debug("Enter the Directory Path of imaging files: "+directToWatch)
        directoryP = ""
        directoryP = directToWatch
        directToWatch, inputfiles = check_imaging_files(directToWatch, 'imaging')

        # targetName = str(input("Please Enter the Planet Name. Make sure it matches the case sensitive name used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): "))
        # log.debug("Enter the Planet Name. Make sure it matches the case sensitive name used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): "+targetName)

        while True:
            targetName = str(input(
                "\nPlease enter the Planet Name. Make sure it matches the case sensitive name and spacing used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): "))
            log.debug(
                "Please enter the Planet Name. Make sure it matches the case sensitive name and spacing used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): " +
                userpDict['pName'])

            if not userpDict['pName'][-2].isspace():
                print(
                    "The convention on the NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html) is to have a space between the star name and the planet letter. Please confirm that you have properly input the planet's name.")
                log.debug(
                    "The convention on the NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html) is to have a space between the star name and the planet letter. Please confirm that you have properly input the planet's name.")
                planetnameconfirm = user_input('\nPlease confirm:\n  (1) ' + userpDict[
                    'pName'] + ' is correct.\n  (2) The planet name needs to be changed.\nPlease select 1 or 2: ',
                                               type_=int, val1=1, val2=2)
                log.debug('\nPlease confirm:\n  (1) ' + userpDict[
                    'pName'] + ' is correct.\n  (2) The planet name needs to be changed.\nPlease select 1 or 2: ' + str(
                    planetnameconfirm))
            else:
                break
            if planetnameconfirm == 1:
                break

        while True:
            try:
                carryOn = input('Type continue after the first image has been taken and saved: ')
                if carryOn != 'continue':
                    raise ValueError
                break
            except ValueError:
                continue

        exotic_UIprevTPX = user_input(targetName + " X Pixel Coordinate: ", type_=int)
        exotic_UIprevTPY = user_input(targetName + " Y Pixel Coordinate: ", type_=int)
        exotic_UIprevRPX = user_input("Comp Star X Pixel Coordinate: ", type_=int)
        exotic_UIprevRPY = user_input("Comp Star Y Pixel Coordinate: ", type_=int)

        log.debug(targetName + " X Pixel Coordinate: " + str(exotic_UIprevTPX))
        log.debug(targetName + " Y Pixel Coordinate: " + str(exotic_UIprevTPY))
        log.debug("Comp Star X Pixel Coordinate: " + str(exotic_UIprevRPX))
        log.debug("Comp Star Y Pixel Coordinate: " + str(exotic_UIprevRPY))

        print('Real Time Plotting ("Control + C" or close the plot to quit)')
        print('\nPlease be patient. It will take at least 15 seconds for the first image to get plotted.')

        fig = plt_exotic.figure()
        exotic_ax1 = fig.add_subplot(1, 1, 1)
        exotic_ax1.set_title(targetName)
        exotic_ax1.set_ylabel('Normalized Flux')
        exotic_ax1.set_xlabel('Time (jd)')

        anim = FuncAnimation(fig, realTimeReduce, fargs=targetName, interval=15000)  # refresh every 15 seconds
        plt_exotic.show()

    ###########################
    # Complete Reduction Routine
    ###########################

    # ----USER INPUTS----------------------------------------------------------
    else:
        log.info('**************************')
        log.info('Complete Reduction Routine')
        log.info('**************************')

        directoryP = ""
        compStarList = []

        exotic_infoDict = {'fitsdir': None, 'saveplot': None, 'flatsdir': None, 'darksdir': None, 'biasesdir': None,
                    'aavsonum': None, 'secondobs': None, 'date': None, 'lat': None, 'long': None, 'elev': None,
                    'ctype': None, 'pixelbin': None, 'filter': None, 'wl_min': None, 'wl_max': None, 'notes': None,
                    'tarcoords': None, 'compstars': None, 'plate_opt': None, 'pixel_scale': None}

        userpDict = {'ra': None, 'dec': None, 'pName': None, 'sName': None, 'pPer': None, 'pPerUnc': None,
                     'midT': None, 'midTUnc': None, 'rprs': None, 'rprsUnc': None, 'aRs': None, 'aRsUnc': None,
                     'inc': None, 'incUnc': None, 'ecc': None, 'teff': None,
                     'teffUncPos': None, 'teffUncNeg': None, 'met': None, 'metUncPos': None, 'metUncNeg': None,
                     'logg': None, 'loggUncPos': None, 'loggUncNeg': None}

        fitsortext = user_input('Enter "1" to perform aperture photometry on fits files or "2" to start '
                                    'with pre-reduced data in a .txt format: ', type_=int, val1=1, val2=2)

        log.debug('Enter "1" to perform aperture photometry on fits files or "2" to start '
                                '\nwith pre-reduced data in a .txt format: '+str(fitsortext))

        fileorcommandline = user_input('How would you like to input your initial parameters? '
                                        'Enter "1" to use the Command Line or "2" to use an input file: ',
                                        type_=int, val1=1, val2=2)

        log.debug('How would you like to input your initial parameters? '
                  '\nEnter "1" to use the Command Line or "2" to use an input file: '+str(fileorcommandline))

        # os.path.join(os.path.split(os.getcwd())[0], '')

        # Read in input file rather than using the command line
        if fileorcommandline == 2:
            exotic_infoDict, userpDict = get_initialization_file(exotic_infoDict, userpDict)
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
                exotic_infoDict['fitsdir'] = str(input("\nEnter the directory path where imaging files are located. (Example using the sample data: sample-data/HatP32Dec202017): "))
                log.debug("Enter the directory path where imaging files are located. (Example using the sample data: sample-data/HatP32Dec202017): " + str(exotic_infoDict['fitsdir']))

            exotic_infoDict['fitsdir'], inputfiles = check_imaging_files(exotic_infoDict['fitsdir'], 'imaging')
        else:
            datafile = str(input("Enter the path and filename of your data file: "))
            log.debug("Enter the path and filename of your data file: "+str(datafile))
            if datafile == 'ok':
                datafile = "/Users/rzellem/Documents/EXOTIC/sample-data/NormalizedFluxHAT-P-32 bDecember 17, 2017.txt"
                log.debug("Hello, Rob.")
                # datafile = "/Users/rzellem/Downloads/fluxorama.csv
            try:
                initf = open(datafile, 'r')
            except FileNotFoundError:
                print("ERROR: Data file not found. Please try again.")
                log.debug("ERROR: Data file not found. Please try again.")
                sys.exit()

            exotic_infoDict['exposure'] = user_input("Please enter your image exposure time in seconds: ", type_=int)
            log.debug("Please enter your image exposure time in seconds: " + str(exotic_infoDict['exposure']))

            processeddata = initf.readlines()

        if fileorcommandline == 1:
            exotic_infoDict['saveplot'] = input("Enter the directory to save the results and plots into or type new to create one: ")
            log.debug("Enter the directory to save the results and plots into or type new to create one: " + str(exotic_infoDict['saveplot']))

        exotic_infoDict['saveplot'] = get_save_directory(exotic_infoDict['saveplot'])

        # Make a temp directory of helpful files
        try:
            os.makedirs(exotic_infoDict['saveplot'] + "temp/")
        except FileExistsError:
            # directory already exists
            pass

        if fileorcommandline == 1:
            while True:
                userpDict['pName'] = str(input("\nEnter the Planet Name. Make sure it matches the case sensitive name and spacing used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): "))
                log.debug("Enter the Planet Name. Make sure it matches the case sensitive name and spacing used on Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html): "+userpDict['pName'])

                if not userpDict['pName'][-2].isspace():
                    print("The convention on the NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html) is to have a space between the star name and the planet letter. Please confirm that you have properly input the planet's name.")
                    log.debug(
                        "The convention on the NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/index.html) is to have a space between the star name and the planet letter. Please confirm that you have properly input the planet's name.")
                    planetnameconfirm = user_input('\nPlease confirm:\n  (1) '+userpDict['pName']+' is correct.\n  (2) The planet name needs to be changed.\nPlease select 1 or 2: ', type_=int, val1=1, val2=2)
                    log.debug('\nPlease confirm:\n  (1) ' + userpDict[
                        'pName'] + ' is correct.\n  (2) The planet name needs to be changed.\nPlease select 1 or 2: '+str(planetnameconfirm))
                else:
                    break
                if planetnameconfirm==1:
                    break

        nea_obj = NASAExoplanetArchive(planet=userpDict['pName'])
        userpDict['pName'], CandidatePlanetBool, pDict = nea_obj.planet_info()

        # observation date
        if fileorcommandline == 1:
            exotic_infoDict['date'] = str(input("\nEnter the Observation Date (MM-DD-YYYY): "))
            log.debug("Enter the Observation Date (MM-DD-YYYY): " + str(exotic_infoDict['date']))

        # Using a / in your date can screw up the file paths- this will check user's date
        while "/" in exotic_infoDict['date']:
            print("Do not use / in your date. Please try again.")
            exotic_infoDict['date'] = str(input("\nEnter the Observation Date (MM-DD-YYYY): "))
            log.debug("Do not use / in your date. Please try again.")
            log.debug("Enter the Observation Date (MM-DD-YYYY): " + str(exotic_infoDict['date']))

        if fitsortext == 1:
            if fileorcommandline == 1:
                exotic_infoDict['lat'] = input("Enter the latitude (in degrees) of where you observed. "
                                        "Don't forget the sign where North is '+' and South is '-'! "
                                        "(Example: +50.4): ")
                log.debug("Enter the latitude of where you observed (deg) " + str(exotic_infoDict['lat']))
            # Latitude
            while True:
                try:
                    exotic_infoDict['lat'] = exotic_infoDict['lat'].replace(' ', '')
                    if exotic_infoDict['lat'][0] != '+' and exotic_infoDict['lat'][0] != '-':
                        raise ValueError("You forgot the sign for the latitude! North is '+' and South is '-'. Please try again.")
                    lati = float(exotic_infoDict['lat'])
                    if lati <= -90.00 or lati >= 90.00:
                        raise ValueError('Your latitude is out of range. Please enter a latitude between -90 and +90 (deg)')
                    break
                # check to make sure they have a sign
                except ValueError as err:
                    print(err.args)
                    exotic_infoDict['lat'] = input("Enter the latitude (in degrees) of where you observed. "
                                            "Don't forget the sign where North is '+' and South is '-'! "
                                            "(Example: +50.4): ")
                    log.debug("Enter the latitude (in degrees) of where you observed. " + str(exotic_infoDict['lat']))

            if fileorcommandline == 1:
                exotic_infoDict['long'] = input("Enter the longitude (in degrees) of where you observed. "
                                         "(Don't forget the sign where East is '+' and West is '-')! "
                                         "(Example: -32.12): ")
                log.debug("Enter the longitude (in degrees) of where you observed. " + str(exotic_infoDict['lat']))
            # Longitude
            while True:
                try:
                    exotic_infoDict['long'] = exotic_infoDict['long'].replace(' ', '')
                    if exotic_infoDict['long'][0] != '+' and exotic_infoDict['long'][0] != '-':
                        raise ValueError("You forgot the sign for the longitude! East is '+' and West is '-'. Please try again.")
                    longit = float(exotic_infoDict['long'])
                    if longit <= -180.00 or longit >= 180.00:
                        raise ValueError('Your longitude is out of range. Please enter a longitude between -180 and +180 (deg)')
                    break
                # check to make sure they have a sign
                except ValueError as err:
                    print(err.args)
                    exotic_infoDict['long'] = input("Enter the longitude (in degrees) of where you observed. "
                                             "(Don't forget the sign where East is '+' and West is '-')! "
                                             "(Example: -32.12): ")
                    log.debug("Enter the longitude (in degrees) of where you observed. " + str(exotic_infoDict['long']))

            if fileorcommandline == 1:
                exotic_infoDict['elev'] = user_input("Enter the elevation (in meters) of where you observed: ", type_=float)
                log.debug("Enter the elevation (in meters) of where you observed (deg) " + str(exotic_infoDict['elev']))

            # TARGET STAR
            if fileorcommandline == 1:
                exotic_UIprevTPX = user_input('\n' + userpDict['pName'] + " X Pixel Coordinate: ", type_=int)
                exotic_UIprevTPY = user_input(userpDict['pName'] + " Y Pixel Coordinate: ", type_=int)
                numCompStars = user_input("How many comparison stars would you like to use? (1-10) ", type_=int)

                log.debug('' + userpDict['pName'] + " X Pixel Coordinate: " + str(exotic_UIprevTPX))
                log.debug('' + userpDict['pName'] + " Y Pixel Coordinate: " + str(exotic_UIprevTPY))
                log.debug("How many comparison stars would you like to use? (1-10) "+str(numCompStars))

                for num in range(numCompStars):
                    xpix = user_input("Comparison Star %s X Pixel Coordinate: " % str(num+1), type_=int)
                    ypix = user_input("Comparison Star %s Y Pixel Coordinate: " % str(num+1), type_=int)
                    log.debug(f"Comparison Star {str(num+1)} X Pixel Coordinate: {xpix} ")
                    log.debug(f"Comparison Star {str(num+1)} Y Pixel Coordinate: {ypix} ")
                    compStarList.append((xpix, ypix))

            # ---HANDLE CALIBRATION IMAGES------------------------------------------------
            if fileorcommandline == 1:
                cals = user_input('\nDo you have any calibration images (flats, darks or biases)? (y/n): ', type_=str, val1='y', val2='n')
                log.debug('Do you have any calibration images (flats, darks or biases)? (y/n): '+cals)

            # if they have cals, handle them by calculating the median flat, dark or bias
            if cals == 'y':

                # flats
                # THIS DOES NOT ACCOUNT FOR CALIBRATING THE FLATS, WHICH COULD BE TAKEN AT A DIFFERENT EXPOSURE TIME
                if fileorcommandline == 1:
                    flats = user_input('\nDo you have flats? (y/n): ', type_=str, val1='y', val2='n')
                    log.debug('Do you have flats? (y/n): '+flats)
                    if flats == 'y':
                        flatsBool = True
                        exotic_infoDict['flatsdir'] = str(input('Enter the directory path to your flats (must be in their own separate folder): '))  # +"/*.FITS"
                        log.debug('Enter the directory path to your flats (must be in their own separate folder): ' + exotic_infoDict['flatsdir'])
                    else:
                        flatsBool = False

                    if flatsBool:
                        exotic_infoDict['flatsdir'], inputflats = check_imaging_files(exotic_infoDict['flatsdir'], 'flats')
                        flatsImgList = []
                        for flatFile in inputflats:
                            flatData = fits.getdata(flatFile, ext=0)
                            flatsImgList.append(flatData)
                        notNormFlat = np.median(flatsImgList, axis=0)

                        # NORMALIZE
                        medi = np.median(notNormFlat)
                        generalFlat = notNormFlat / medi
                else:
                    flatsBool = False

                # darks
                if fileorcommandline == 1:
                    darks = user_input('\nDo you have darks? (y/n): ', type_=str, val1='y', val2='n')
                    log.debug('Do you have darks? (y/n): '+darks)
                    if darks == 'y':
                        darksBool = True
                        exotic_infoDict['darksdir'] = str(input('Enter the directory path to your darks (must be in their own separate folder): '))  # +"/*.FITS"
                        log.debug('Enter the directory path to your darks (must be in their own separate folder): ' + exotic_infoDict['darksdir'])
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
                    biases = user_input('\nDo you have biases? (y/n): ', type_=str, val1='y', val2='n')
                    log.debug('Do you have biases? (y/n): '+biases)
                    if biases == 'y':
                        biasesBool = True
                        exotic_infoDict['biasesdir'] = str(input('Enter the directory path to your biases (must be in their own separate folder): '))  # +"/*.FITS"
                        log.debug('Enter the directory path to your biases (must be in their own separate folder): ' + exotic_infoDict['biasesdir'])
                    else:
                        biasesBool = False

                if biasesBool:
                    # Add / to end of directory if user does not input it
                    exotic_infoDict['biasesdir'], inputbiases = check_imaging_files(exotic_infoDict['biasesdir'], 'biases')
                    biasesImgList = []
                    for biasFile in inputbiases:
                        biasData = fits.getdata(biasFile, ext=0)
                        biasesImgList.append(biasData)
                    generalBias = np.median(biasesImgList, axis=0)
            else:
                flatsBool = False
                darksBool = False
                biasesBool = False

        print("\n***************************************")
        log.debug("***************************************")

        # Handle AAVSO Formatting
        if fileorcommandline == 1:
            # userNameEmails = str(input('Please enter your name(s) and email address(es) in the format: Your Name (youremail@example.com), Next Name (nextemail@example.com), etc.  '))
            exotic_infoDict['aavsonum'] = str(input('Please enter your AAVSO Observer Account Number (type N/A if you do not currently have an account): '))
            exotic_infoDict['secondobs'] = str(input('Please enter your comma-separated secondary observer codes (or type N/A if only 1 observer code): '))
            exotic_infoDict['ctype'] = str(input("Please enter your camera type (CCD or DSLR): "))
            exotic_infoDict['pixelbin'] = str(input('Please enter your pixel binning: '))
            # infoDict['exposure'] = user_input('Please enter your exposure time (seconds): ', type_=int)
            exotic_infoDict['filter'] = str(input('Please enter your filter name from the options at '
                                           'http://astroutils.astronomy.ohio-state.edu/exofast/limbdark.shtml: '))
            exotic_infoDict['notes'] = str(input('Please enter any observing notes (seeing, weather, etc.): '))

            log.debug('Please enter your AAVSO Observer Account Number (type N/A if you do not currently have an account): ' + exotic_infoDict['aavsonum'])
            log.debug('Please enter your comma-separated secondary observer codes (or type N/A if only 1 observer code): ' + exotic_infoDict['secondobs'])
            log.debug("Please enter your camera type (CCD or DSLR): " + exotic_infoDict['ctype'])
            log.debug('Please enter your pixel binning: ' + exotic_infoDict['pixelbin'])
            log.debug('Please enter your filter name from the options at '
                                           'http://astroutils.astronomy.ohio-state.edu/exofast/limbdark.shtml: ' + exotic_infoDict['filter'])
            log.debug('Please enter any observing notes (seeing, weather, etc.): ' + exotic_infoDict['notes'])

        if fileorcommandline == 2:
            diff = check_parameters(userpDict, pDict)
            if diff:
                pDict = get_planetary_parameters(CandidatePlanetBool, userpDict, pdict=pDict)
        else:
            pDict = get_planetary_parameters(CandidatePlanetBool, userpDict, pdict=pDict)

        ld_obj = LimbDarkening(teff=pDict['teff'], teffpos=pDict['teffUncPos'], teffneg=pDict['teffUncNeg'],
                               met=pDict['met'], metpos=pDict['metUncPos'], metneg=pDict['metUncNeg'],
                               logg=pDict['logg'], loggpos=pDict['loggUncPos'], loggneg=pDict['loggUncNeg'],
                               wl_min=exotic_infoDict['wl_min'], wl_max=exotic_infoDict['wl_max'], filter_type=exotic_infoDict['filter'])
        ld0, ld1, ld2, ld3, exotic_infoDict['filter'] = ld_obj.nonlinear_ld()

        if fitsortext == 1:
            print('\n**************************')
            print('Starting Reduction Process')
            print('**************************')

            log.debug('**************************'
                                '\nStarting Reduction Process'
                                '\n**************************')

            #########################################
            # FLUX DATA EXTRACTION AND MANIPULATION
            #########################################

            # Loop placed to check user-entered x and y target coordinates against WCS.
            while True:
                fileNumber = 1
                allImageData, timeList, fileNameList, timesListed, airMassList, fileNameStr, exptimes = [], [], [], [], [], [], []

                # ----TIME SORT THE FILES-------------------------------------------------------------
                for fileName in inputfiles:  # Loop through all the fits files in the directory and executes data reduction

                    # fitsHead = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False)  # opens the file

                    # FOR 61'' DATA ONLY: ONLY REDUCE DATA FROM B FILTER
                    # if fitsHead[0].header ['FILTER']== 'Harris-B':
                    #     #TIME
                    #     timeVal = getJulianTime(fitsHead) #gets the julian time registered in the fits header
                    #     timeList.append(timeVal) #adds to time value list
                    #     fileNameList.append (fileName)
                    # fitsHead.close()  # close stream
                    # del fitsHead

                    # Keeps a list of file names
                    fileNameStr.append(fileName)

                    try:
                        hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
                    except OSError as e:
                        log.info(f'Found corrupted file and removing from reduction: {fileName}, ({e})')
                        fileNameStr.remove(fileName)
                        if getattr(hdul, "close", None) and callable(hdul.close):
                            hdul.close()
                        del hdul
                        continue

                    imageheader = hdul[0].header
                    # TIME
                    timeVal = getJulianTime(hdul)  # gets the julian time registered in the fits header
                    timeList.append(timeVal)  # adds to time value list
                    fileNameList.append(fileName)

                    # TIME
                    currTime = getJulianTime(hdul)
                    timesListed.append(currTime)

                    # AIRMASS
                    airMass = getAirMass(hdul, pDict['ra'], pDict['dec'], lati, longit, exotic_infoDict['elev'])  # gets the airmass at the time the image was taken
                    airMassList.append(airMass)  # adds that airmass value to the list of airmasses

                    # IMAGES
                    hdul[0].scale('float32')
                    allImageData.append(hdul[0].data)

                    # EXPOSURE_TIME
                    exp = imageheader.get('EXPTIME')  # checking for variation in .fits header format
                    if exp:
                        exptimes.append(imageheader['EXPTIME'])
                    else:
                        exptimes.append(imageheader['EXPOSURE'])

                    hdul.close()  # closes the file to avoid using up all of computer's resources
                    del hdul

                consistent_et = False
                if len(exptimes) > 0:
                    consistent_et = all(elem == exptimes[0] for elem in exptimes)

                exptimes = np.array(exptimes)

                if consistent_et:
                    # print("All Elements in List are Equal")
                    exotic_infoDict['exposure'] = exptimes[0]
                else:
                    # print("All Elements in List are Not Equal")
                    exotic_infoDict['exposure'] = np.median(exptimes)
                    # print(infoDict['exposure'])

                # Recast list as numpy arrays
                allImageData = np.array(allImageData)
                timesListed = np.array(timesListed)
                airMassList = np.array(airMassList)

                # TODO: Is this dead code? The vars pointing and location are undefined.
                # TODO: comment out conditional block?
                # If all of the airmasses == 1, then you need to calculate the airmass for the user
                if set(airMassList) == 1:
                    # pointingAltAz = pointing.transform_to(AltAz(obstime=t, location=location))
                    pass

                # # Time sorts the file names based on the fits file header
                # timeSortedNames = [x for _, x in sorted(zip(timeList, fileNameList))]
                # tsnCopy = timeSortedNames

                # sorts the times for later plotting use
                sortedallImageData = allImageData[np.argsort(timeList)]
                timesListed = timesListed[np.argsort(timeList)]
                airMassList = airMassList[np.argsort(timeList)]
                # sortedTimeList = sorted(timeList)

                # print("\nEXOTIC now has the option to filter the raw images for cosmic rays. Typically, images do not need this filter. However, if you run into an error while running EXOTIC, give this a try. As a heads up, this can take a few minutes.")
                # cosmicrayfilter = user_input("\nDo you want to filter the raw images for cosmic rays? (y/n): ")
                # if cosmicrayfilter.lower() == "yes" or cosmicrayfilter.lower() == "y":
                #     cosmicrayfilter_bool = True
                # else:
                #     cosmicrayfilter_bool = False

                # TODO add option to inits file
                cosmicrayfilter_bool = True
                if cosmicrayfilter_bool:
                    log.info("Filtering your data for cosmic rays.")
                    targx, targy, targamplitude, targsigX, targsigY, targrot, targoff = fit_centroid(sortedallImageData[0], [exotic_UIprevTPX, exotic_UIprevTPY], box=10)
                    psffwhm = 2.355*(targsigX+targsigY)/2
                    log.debug("FWHM in 1st image: {:.2f} px".format(np.round(psffwhm, 2)))
                    log.debug("STDEV before: {:.2f}".format(np.std(sortedallImageData, 0).mean()))

                    # # -------COSMIC RAY FILTERING-----------------------------------------------------------------------

                    animate_toggle(True)

                    for ii in range(len(sortedallImageData)):

                        # remove nans                        
                        nanmask = np.isnan(sortedallImageData[ii])
                        if nanmask.sum() > 0:
                            bg = generic_filter(sortedallImageData[ii], np.nanmedian, (3, 3))
                            sortedallImageData[ii][nanmask] = bg[nanmask]
                        
                        mask, clean = detect_cosmics(sortedallImageData[ii], psfmodel='gauss',  psffwhm=psffwhm, psfsize=2*round(psffwhm)+1,  sepmed=False, sigclip = 4.25, niter=2, objlim=10, cleantype='idw', verbose=False)                        
                        sortedallImageData[ii] = clean
                        
                        # TODO move to function
                        # if ii == 0:
                        # f,ax = plt_exotic.subplots(1,3,figsize=(18,8))
                        # ax[0].imshow(np.log10(ogdata),vmin=np.percentile(np.log10(bg),5), vmax=np.percentile(np.log10(bg),99))
                        # ax[0].set_title("Original Data")
                        # ax[1].imshow(mask,cmap='binary_r')
                        # ax[1].set_title("Cosmic Ray Mask")
                        # ax[2].imshow(np.log10(sortedallImageData[ii]),vmin=np.percentile(np.log10(bg),5), vmax=np.percentile(np.log10(bg),99))
                        # ax[2].set_title("Corrected Image")
                        # plt_exotic.tight_layout()
                        # plt_exotic.show()

                    animate_toggle()

                log.debug("STDEV after: {:.2f}".format(np.std(sortedallImageData, 0).mean()))

                # -------OPTIMAL COMP STAR, APERTURE, AND ANNULUS CALCULATION----------------------------------------

                # Loops through all of the possible aperture and annulus radius
                # guess at optimal aperture by doing a gaussian fit and going out 3 sigma as an estimate

                # hdul = fits.open(name=timeSortedNames[0], memmap=False, cache=False, lazy_load_hdus=False)  # opens the fits file
                # firstImageData = hdul['ext', 0].data  # fits.getdata(timeSortedNames[0], ext=0)
                firstimagecounter = 0
                firstImageData = sortedallImageData[firstimagecounter]

                # Sometimes the first image is a bad one...in that case, we iterate until we do not fail
                while True:
                    # fit Target in the first image and use it to determine aperture and annulus range
                    try:
                        targx, targy, targamplitude, targsigX, targsigY, targrot, targoff = fit_centroid(firstImageData, [exotic_UIprevTPX, exotic_UIprevTPY],
                                                                                                         box=10)
                        break
                    # If the first image is a bad one, then move on to the next image
                    except Exception:
                        firstimagecounter += 1
                        firstImageData = sortedallImageData[firstimagecounter]

                # Filter the other data as well
                sortedallImageData = sortedallImageData[firstimagecounter:]
                timesListed = timesListed[firstimagecounter:]
                airMassList = airMassList[firstimagecounter:]
                # sortedTimeList = sortedTimeList[firstimagecounter:]

                # apply cals correction if applicable
                if darksBool:
                    print("Dark subtracting images.")
                    log.debug("Dark subtracting images.")
                    sortedallImageData = sortedallImageData - generalDark
                elif biasesBool:
                    print("Bias-correcting images.")
                    log.debug("Bias-correcting images.")
                    sortedallImageData = sortedallImageData - generalBias
                else:
                    pass

                if flatsBool:
                    print("Flattening images.")
                    log.debug("Flattening images.")
                    sortedallImageData = sortedallImageData / generalFlat

                # Reference File
                refFile = os.path.join(exotic_infoDict['saveplot'], 'ref_file_%s_%s'
                                       % (firstimagecounter, os.path.split(fileNameStr[firstimagecounter])[-1]))

                # Removes existing file of reference file
                try:
                    os.remove(refFile)
                except OSError:
                    pass
                convertToFITS = fits.PrimaryHDU(data=sortedallImageData[0])
                convertToFITS.writeto(refFile)
                print('\nHere is the path to the reference imaging file used by EXOTIC: \n' + refFile)
                log.debug('Here is the path to the reference imaging file used by EXOTIC: \n' + refFile)
                wcsFile = check_wcs(refFile, exotic_infoDict['saveplot'], exotic_infoDict['plate_opt'])
                hdulWCS = None

                # Check pixel coordinates by converting to WCS. If not correct, loop over again
                if wcsFile:
                    log.info(f'Here is the path to your plate solution: {str(wcsFile)}')
                    hdulWCS = fits.open(name=wcsFile, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
                    rafile, decfile = get_radec(hdulWCS)

                    # Save previously entered x and y pixel coordinates before checking against plate solution
                    saveUIprevTPX, saveUIprevTPY = exotic_UIprevTPX, exotic_UIprevTPY
                    exotic_UIprevTPX, exotic_UIprevTPY = check_targetpixelwcs(exotic_UIprevTPX, exotic_UIprevTPY, pDict['ra'],
                                                                              pDict['dec'], rafile, decfile)
                    # If the coordinates were not changed, do not loop over again
                    if exotic_UIprevTPX == saveUIprevTPX and exotic_UIprevTPY == saveUIprevTPY:
                        break
                else:
                    break

            # Image Alignment
            sortedallImageData, boollist = image_alignment(sortedallImageData)

            timesListed = timesListed[boollist]
            airMassList = airMassList[boollist]

            minAperture = max(1, int(2 * max(targsigX, targsigY)))
            maxAperture = int(5 * max(targsigX, targsigY) + 1)
            minAnnulus = 2
            maxAnnulus = 5
            # exit()
            # fit centroids for first image to determine priors to be used later
            for compCounter in range(0, len(compStarList)):
                print('\n\n***************************************************************')
                print('Determining Optimal Aperture and Annulus Size for Comp Star #' + str(compCounter + 1))
                print('***************************************************************')

                log.debug('***************************************************************')
                log.debug('Determining Optimal Aperture and Annulus Size for Comp Star #' + str(compCounter + 1))
                log.debug('***************************************************************\n')

                # #just in case comp star drifted off and timeSortedNames had to be altered, reset it for the new comp star
                # timeSortedNames = tsnCopy

                exotic_UIprevRPX, exotic_UIprevRPY = compStarList[compCounter]
                print('Target X: ' + str(round(targx)) + ' Target Y: ' + str(round(targy)))
                log.debug('Target X: ' + str(round(targx)) + ' Target Y: ' + str(round(targy)))
                refx, refy, refamplitude, refsigX, refsigY, retrot, refoff = fit_centroid(firstImageData, [exotic_UIprevRPX, exotic_UIprevRPY],
                                                                                          box=10)
                print('Comparison X: ' + str(round(refx)) + ' Comparison Y: ' + str(round(refy)) + '\n')
                log.debug('Comparison X: ' + str(round(refx)) + ' Comparison Y: ' + str(round(refy)) + '\n')

                # If plate solution was generated, use it to check if the comparison stars selected are variable
                # If yes, skip determining optimal aperture and annulus for that comparison star
                if wcsFile:
                    print("Checking for variability in current comparison star... ")
                    log.debug("Checking for variability in current comparison star... ")
                    if variableStarCheck(refx, refy, hdulWCS):
                        print("Current comparison star is variable, proceeding to next star.")
                        log.debug("Current comparison star is variable, proceeding to next star.")
                        continue

                # determines the aperture and annulus combinations to iterate through based on the sigmas of the LM fit
                aperture_min = int(3 * np.nanmax([targsigX, targsigY]))
                aperture_max = int(5 * np.nanmax([targsigX, targsigY]))

                # run through apertures based on PSF shape
                if aperture_min <= 1:
                    aperture_sizes = np.arange(1, 10, 2)
                else:
                    aperture_sizes = np.round(np.linspace(aperture_min, aperture_max, 10), 2)

                aperture_sizes = np.append(aperture_sizes, -1*aperture_sizes)  # no comparison star
                aperture_sizes = np.append(aperture_sizes, 0)  # PSF fit

                # single annulus size
                annulus_sizes = [10, 12, 15]

                target_fits = {}
                ref_fits = {}
                reg_trans = {}

                for apertureR in aperture_sizes:  # aperture loop
                    for annulusR in annulus_sizes:  # annulus loop # no need
                        # don't reprocess
                        if apertureR < 0 < compCounter:
                            continue

                        # only do PSF fit in first annulus for loop
                        if apertureR == 0 and annulusR != annulus_sizes[0]:
                            continue

                        if apertureR == 0:
                            print('Testing Comparison Star #' + str(compCounter+1) + ' with a PSF photometry.')
                            log.debug('Testing Comparison Star #' + str(compCounter + 1) + ' with a PSF photometry.')
                        elif apertureR < 0 and compCounter == 0:
                            print('Testing NO Comparison Star with a '+str(abs(apertureR))+' pixel aperture and a '+str(abs(annulusR))+' pixel annulus.')
                            log.debug('Testing NO Comparison Star with a ' + str(
                                abs(apertureR)) + ' pixel aperture and a ' + str(abs(annulusR)) + ' pixel annulus.')
                        else:
                            print('Testing Comparison Star #' + str(compCounter+1) + ' with a '+str(apertureR)+' pixel aperture and a '+str(annulusR)+' pixel annulus.')
                            log.debug('Testing Comparison Star #' + str(compCounter + 1) + ' with a ' + str(
                                apertureR) + ' pixel aperture and a ' + str(annulusR) + ' pixel annulus.')

                        for fileNumber, imageData in enumerate(sortedallImageData):
                            if apertureR == 0: # psf fit
                                continue

                            # Find the target star in the image and get its pixel coordinates if it is the first file
                            if fileNumber == 0:
                                # Initializing the star location guess as the user inputted pixel coordinates
                                prevTPX, prevTPY, prevRPX, prevRPY = exotic_UIprevTPX, exotic_UIprevTPY, exotic_UIprevRPX, exotic_UIprevRPY  # 398, 275, 419, 203
                                prevTSigX, prevTSigY, prevRSigX, prevRSigY = targsigX, targsigY, refsigX, refsigY
                                prevImageData = imageData  # no shift should be registered

                            # ------ CENTROID FITTING ----------------------------------------

                            # corrects for any image shifts that result from a tracking slip
                            # shift, error, diffphase = phase_cross_correlation(prevImageData, imageData)
                            if fileNumber in reg_trans.keys():
                                shift, error, diffphase = reg_trans[fileNumber]
                            else:
                                shift, error, diffphase = phase_cross_correlation(prevImageData, imageData)
                                reg_trans[fileNumber] = [shift, error, diffphase]

                            xShift = shift[1]
                            yShift = shift[0]

                            prevTPX = prevTPX - xShift
                            prevTPY = prevTPY - yShift
                            prevRPX = prevRPX - xShift
                            prevRPY = prevRPY - yShift

                            # set target search area
                            txmin = int(prevTPX) - exotic_distFC  # left
                            txmax = int(prevTPX) + exotic_distFC  # right
                            tymin = int(prevTPY) - exotic_distFC  # top
                            tymax = int(prevTPY) + exotic_distFC  # bottom

                            # boolean that represents if either the target or comp star gets too close to the detector
                            driftBool = False

                            # check if your target star is too close to the edge of the detector
                            if txmin <= 0 or tymin <= 0 or txmax >= len(imageData) or tymax >= len(imageData[0]):
                                print('*************************************************************************************')
                                print('WARNING: In image '+str(fileNumber)+', your target star has drifted too close to the edge of the detector.')
                                # tooClose = int(input('Enter "1" to pick a new comparison star or enter "2" to continue using the same comp star, with the images with all the remaining images ignored \n'))
                                print('All the remaining images after image #'+str(fileNumber-1)+' will be ignored')

                                log.debug('*************************************************************************************'
                                          '\nWARNING: In image '+str(fileNumber)+', your target star has drifted too close to the edge of the detector.'
                                          '\nAll the remaining images after image #'+str(fileNumber-1)+' will be ignored\n')

                                driftBool = True

                                # mask off the rest of timeSortedNames and then ignore the rest of the procedure until

                            # Set reference search area
                            rxmin = int(prevRPX) - exotic_distFC  # left
                            rxmax = int(prevRPX) + exotic_distFC  # right
                            rymin = int(prevRPY) - exotic_distFC  # top
                            rymax = int(prevRPY) + exotic_distFC  # bottom

                            # check if the reference is too close to the edge of the detector
                            if rxmin <= 0 or rymin <= 0 or rxmax >= len(imageData[0]) or rymax >= len(imageData):
                                print('*************************************************************************************')
                                print('WARNING: In image '+str(fileNumber)+', your reference star has drifted too close to the edge of the detector.')
                                # tooClose = int(input('Enter "1" to pick a new comparison star or enter "2" to continue using the same comp star, with the images with all the remaining images ignored \n'))
                                print('All the remaining images after image #'+str(fileNumber-1)+' will be ignored for this comparison star')
                                print('*************************************************************************************')

                                log.debug('*************************************************************************************'
                                          '\nWARNING: In image '+str(fileNumber)+', your target star has drifted too close to the edge of the detector.'
                                          '\nAll the remaining images after image #'+str(fileNumber-1)+' will be ignored for this comparison star\n')

                                driftBool = True

                            # if the star isn't too close, then proceed as normal
                            if not driftBool:
                                targSearchA = imageData[tymin:tymax, txmin:txmax]
                                refSearchA = imageData[rymin:rymax, rxmin:rxmax]

                                targPos = [prevTPX, prevTPY]

                                # get minimum background value bigger than 0
                                targImFlat = np.sort(np.array(targSearchA).ravel())

                                # Initialize the variable
                                tGuessBkg = 0
                                for el in targImFlat:
                                    if el > 0:
                                        tGuessBkg = el
                                        break

                                refImFlat = np.sort(np.array(refSearchA).ravel())
                                for rel in refImFlat:
                                    if rel > 0:
                                        rGuessBkg = rel
                                        break

                                # Guess at Gaussian Parameters and feed them in to help gaussian fitter

                                tGuessAmp = targSearchA.max() - tGuessBkg
                                if tGuessAmp < 0:
                                    print('ERROR: the Darks have a higher pixel counts than the image itself')
                                    log.debug('ERROR: the Darks have a higher pixel counts than the image itself')
                                myPriors = [tGuessAmp, prevTSigX, prevTSigY, 0, tGuessBkg]

                                # tx, ty, tamplitude, tsigX, tsigY, toff = fit_centroid(imageData, targPos,
                                #                                                     init=myPriors, box=distFC)
                                if fileNumber in target_fits.keys():
                                    tx, ty, tamplitude, tsigX, tsigY, toff = target_fits[fileNumber]
                                else:
                                    tx, ty, tamplitude, tsigX, tsigY, trot, toff = fit_centroid(imageData, targPos,
                                                                                                init=myPriors, box=exotic_distFC)
                                    target_fits[fileNumber] = [tx, ty, tamplitude, tsigX, tsigY, toff]

                                currTPX = tx
                                currTPY = ty

                                # append to list of target centroid positions for later plotting
                                xTargCent.append(currTPX)
                                yTargCent.append(currTPY)

                                rGuessAmp = refSearchA.max() - rGuessBkg
                                myRefPriors = [rGuessAmp, prevRSigX, prevRSigY, 0, rGuessBkg]
                                # rx, ry, ramplitude, rsigX, rsigY, roff = fit_centroid(imageData, [prevRPX, prevRPY],
                                # init=myRefPriors, box=distFC)
                                if fileNumber in ref_fits.keys():
                                    rx, ry, ramplitude, rsigX, rsigY, roff = ref_fits[fileNumber]
                                else:
                                    rx, ry, ramplitude, rsigX, rsigY, rrot, roff = fit_centroid(imageData, [prevRPX, prevRPY],
                                                                                                init=myRefPriors, box=exotic_distFC)
                                    ref_fits[fileNumber] = [rx, ry, ramplitude, rsigX, rsigY, roff]
                                currRPX = rx
                                currRPY = ry

                                # append to list of reference centroid positions for later plotting
                                xRefCent.append(currRPX)
                                yRefCent.append(currRPY)

                                if tamplitude < 0 or tsigX < 0 or tsigY < 0:  # gets rid of negative amplitude values that indicate it couldn't fit gaussian
                                    print('Could not fit 2D gaussian to Target for File Number' + str(fileNumber))
                                    log.debug('Could not fit 2D gaussian to Target for File Number' + str(fileNumber))

                                elif ramplitude < 0 or rsigX < 0 or rsigY < 0:  # gets rid of negative amplitude values that indicate it couldn't fit gaussian
                                    print('Could not fit 2D gaussian to Comparison Star for File Number' + str(fileNumber))
                                    log.debug('Could not fit 2D gaussian to Comparison Star for File Number' + str(
                                        fileNumber))

                                else:
                                    # ------FLUX CALCULATION WITH BACKGROUND SUBTRACTION----------------------------------
                                    
                                    # gets the flux value of the target star and subtracts the background light
                                    tFluxVal, tTotCts = getFlux(imageData, currTPX, currTPY, abs(apertureR), annulusR)
                                    targetFluxVals.append(tFluxVal)  # adds tFluxVal to the total list of flux values of target star
                                    targUncertanties.append(np.sqrt(tFluxVal))  # uncertanty on each point is the sqrt of the total counts

                                    # gets the flux value of the reference star and subracts the background light
                                    rFluxVal, rTotCts = getFlux(imageData, currRPX, currRPY, abs(apertureR), annulusR)
                                    referenceFluxVals.append(rFluxVal)  # adds rFluxVal to the total list of flux values of reference star
                                    refUncertanties.append(np.sqrt(rFluxVal))

                                    # # TIME
                                    # currTime = getJulianTime(hDul)
                                    # timesListed.append(currTime)

                                    # ORBITAL PHASE
                                    currentPhase = getPhase(currTime, pDict['pPer'], pDict['midT'])
                                    phasesList.append(currentPhase)  # adds to list of phases

                                    # # AIRMASS
                                    # airMass = getAirMass(hDul)  # gets the airmass at the time the image was taken
                                    # airMassList.append(airMass)  # adds that airmass value to the list of airmasses

                                    # UPDATE PIXEL COORDINATES and SIGMAS
                                    # target
                                    prevTPX = currTPX
                                    prevTPY = currTPY
                                    prevTSigX = tsigX
                                    prevTSigY = tsigY
                                    # reference
                                    prevRPX = currRPX
                                    prevRPY = currRPY
                                    prevRSigX = rsigX
                                    prevTSigY = rsigY

                                # UPDATE FILE COUNT
                                prevImageData = imageData
                                # fileNumber = fileNumber + 1
                                # hDul.close()  # close the stream

                            # otherwise, mask off the rest of the files from time sorted names including the current one
                            else:
                                print("\nFiltering data to account for drifting target.")
                                log.debug("Filtering data to account for drifting target.")
                                # timeSortedNames = timeSortedNames[:fileNumber]

                                # TIME
                                timesListed = timesListed[:fileNumber]

                                # AIRMASS
                                airMassList = airMassList[:fileNumber]

                                # ALL IMAGES
                                sortedallImageData = sortedallImageData[:fileNumber]

                                # boollist = boollist[:fileNumber]

                                break

                        # EXIT THE FILE LOOP

                        # Convert Everything to numpy Arrays
                        arrayTimes = np.array(timesListed)
                        arrayPhases = np.array(phasesList)
                        arrayAirmass = np.array(airMassList)

                        if apertureR == 0:  # psf fit
                            tpsfflux = []
                            rpsfflux = []
                            for k in target_fits.keys():
                                xc, yc, amp, sigx, sigy, off = target_fits[k]
                                tpsfflux.append(2*PI*sigx*sigy*amp)
                                xc, yc, amp, sigx, sigy, off = ref_fits[k]
                                rpsfflux.append(2*PI*sigx*sigy*amp)
                            arrayReferences = np.array(rpsfflux)
                            arrayTUnc = arrayFinalFlux**0.5
                            arrayRUnc = arrayReferences**0.5

                            arrayFinalFlux = np.array(tpsfflux) / arrayReferences
                            arrayNormUnc = arrayTUnc / arrayReferences
                        elif apertureR < 0:  # no comp star
                            arrayReferences = np.array(referenceFluxVals)
                            arrayTUnc = np.array(targUncertanties)
                            arrayRUnc = np.array(refUncertanties)

                            arrayFinalFlux = np.array(targetFluxVals)
                            arrayNormUnc = np.array(targUncertanties)
                        else:
                            # aperture phot
                            arrayReferences = np.array(referenceFluxVals)
                            arrayTUnc = np.array(targUncertanties)
                            arrayRUnc = np.array(refUncertanties)

                            arrayFinalFlux = np.array(targetFluxVals) / arrayReferences
                            arrayNormUnc = arrayTUnc / arrayReferences

                        # Execute sigma_clip
                        try:
                            filtered_data = sigma_clip(arrayFinalFlux, sigma=3)
                        except TypeError:
                            filtered_data = sigma_clip(arrayFinalFlux, sigma=3)

                        # -----LM LIGHTCURVE FIT--------------------------------------
                        prior = {
                            'rprs': pDict['rprs'],    # Rp/Rs
                            'ars': pDict['aRs'],      # a/Rs
                            'per': pDict['pPer'],     # Period [day]
                            'inc': pDict['inc'],      # Inclination [deg]
                            'u0': ld0[0], 'u1': ld1[0], 'u2': ld2[0], 'u3': ld3[0],  # limb darkening (nonlinear)
                            'ecc': pDict['ecc'],     # Eccentricity
                            'omega': 0,          # Arg of periastron
                            'tmid': pDict['midT'],    # time of mid transit [day]
                            'a1': arrayFinalFlux.mean(),  # max() - arrayFinalFlux.min(), #mid Flux
                            'a2': 0,             # Flux lower bound
                        }

                        phase = (arrayTimes[~filtered_data]-pDict['midT'])/prior['per']
                        prior['tmid'] = pDict['midT'] + np.floor(phase).max()*prior['per']
                        upper = prior['tmid'] + np.abs(25*pDict['midTUnc'] + np.floor(phase).max()*25*pDict['pPerUnc'])
                        lower = prior['tmid'] - np.abs(25*pDict['midTUnc'] + np.floor(phase).max()*25*pDict['pPerUnc'])

                        if np.floor(phase).max()-np.floor(phase).min() == 0:
                            print("WARNING!")
                            print(" Estimated mid-transit time is not within the observations")
                            print(" Check Period & Mid-transit time in inits.json. Make sure the uncertainties are not 0 or Nan.")
                            print('  obs start:', arrayTimes[~filtered_data].min())
                            print('    obs end:', arrayTimes[~filtered_data].max())
                            print(' tmid prior:', prior['tmid'])

                            log.debug("WARNING!")
                            log.debug(" Estimated mid-transit time is not within the observations")
                            log.debug(" Check Period & Mid-transit time in inits.json. Make sure the uncertainties are not 0 or Nan.")
                            log.debug('  obs start:', arrayTimes[~filtered_data].min())
                            log.debug('    obs end:', arrayTimes[~filtered_data].max())
                            log.debug(' tmid prior:', prior['tmid'])

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

                        mybounds = {
                            'rprs': [0, pDict['rprs']+3*pDict['rprsUnc']],
                            'tmid': [max(lower, arrayTimes[~filtered_data].min()), min(arrayTimes[~filtered_data].max(), upper)],
                            'ars': [pDict['aRs']-5*pDict['aRsUnc'], pDict['aRs']+5*pDict['aRsUnc']],

                            'a1': [min(0,min(arrayFinalFlux[~filtered_data])), 3*max(arrayFinalFlux[~filtered_data])],
                            'a2': [-3, 3],
                            # 'a3':[0, max(arrayFinalFlux[~filtered_data])]
                        }

                        myfit = lc_fitter(
                            arrayTimes[~filtered_data],
                            arrayFinalFlux[~filtered_data],
                            arrayNormUnc[~filtered_data],
                            arrayAirmass[~filtered_data],
                            prior,
                            mybounds,
                            mode='lm'
                        )

                        for k in myfit.bounds.keys():
                            print("  {}: {:.6f}".format(k, myfit.parameters[k]))  # , myfit.errors[k]))

                        print('The Residual Standard Deviation is: ' + str(round(100*myfit.residuals.std()/np.median(myfit.data), 6))+"%")
                        print('The Mean Squared Error is: ' + str(round(np.sum(myfit.residuals**2), 6)) + '\n')

                        log.debug('The Residual Standard Deviation is: ' + str(round(100*myfit.residuals.std()/np.median(myfit.data), 6))+"%")
                        log.debug('The Mean Squared Error is: ' + str(round(np.sum(myfit.residuals**2), 6)) + '\n')

                        resstd = myfit.residuals.std()/np.median(myfit.data)
                        if minSTD > resstd:  # If the standard deviation is less than the previous min
                            bestCompStar = compCounter + 1
                            minSTD = resstd  # set the minimum standard deviation to that

                            arrayNormUnc = arrayNormUnc * np.sqrt(myfit.chi2 / myfit.data.shape[0])  # scale errorbars by sqrt(rchi2)
                            minAnnulus = annulusR  # then set min aperature and annulus to those values
                            minAperture = apertureR

                            if len(xTargCent) > 0:
                                # gets the centroid trace plots to ensure tracking is working
                                finXTargCentArray = np.array(xTargCent)
                                finYTargCentArray = np.array(yTargCent)
                                finXRefCentArray = np.array(xRefCent)
                                finYRefCentArray = np.array(yRefCent)
                            else:
                                # PSF fit usually
                                finXTargCentArray = np.array([target_fits[k][0] for k in target_fits])
                                finYTargCentArray = np.array([target_fits[k][1] for k in target_fits])
                                finYRefCentArray = np.array([ref_fits[k][1] for k in ref_fits])
                                finXRefCentArray = np.array([ref_fits[k][0] for k in ref_fits])
                                arrayPhases = getPhase(arrayTimes, pDict['pPer'], myfit.parameters['tmid'])

                            # APPLY DATA FILTER
                            # apply data filter sets the lists we want to print to correspond to the optimal aperature
                            finXTargCent = finXTargCentArray[~filtered_data]
                            finYTargCent = finYTargCentArray[~filtered_data]
                            finXRefCent = finXRefCentArray[~filtered_data]
                            finYRefCent = finYRefCentArray[~filtered_data]

                            # sets the lists we want to print to correspond to the optimal aperature
                            goodFluxes = arrayFinalFlux[~filtered_data]
                            nonBJDTimes = arrayTimes[~filtered_data]
                            nonBJDPhases = arrayPhases[~filtered_data]
                            goodAirmasses = arrayAirmass[~filtered_data]
                            goodTargets = arrayFinalFlux[~filtered_data]
                            goodReferences = arrayReferences[~filtered_data]
                            goodTUnc = arrayTUnc[~filtered_data]
                            goodRUnc = arrayRUnc[~filtered_data]
                            goodNormUnc = arrayNormUnc[~filtered_data]
                            goodResids = myfit.residuals
                            bestlmfit = myfit

                        # Reinitialize the the arrays to be empty
                        # airMassList = []
                        phasesList = []
                        # timesListed = []
                        targetFluxVals = []
                        referenceFluxVals = []
                        targUncertanties = []
                        refUncertanties = []
                        normUncertainties = []
                        xTargCent = []
                        yTargCent = []
                        xRefCent = []
                        yRefCent = []

                    # Exit aperture loop
                # Exit annulus loop
            # Exit the Comp Stars Loop

            if minAperture == 0:  # psf
                print('\n*********************************************')
                print('Best Comparison Star: #' + str(bestCompStar))
                print('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                print('Optimal Method: PSF photometry')
                print('********************************************\n')

                log.debug('*********************************************')
                log.debug('Best Comparison Star: #' + str(bestCompStar))
                log.debug('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                log.debug('Optimal Method: PSF photometry')
                log.debug('********************************************\n')

            elif minAperture < 0:  # no comp star
                print('\n*********************************************')
                print('Best Comparison Star: None')
                print('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                print('Optimal Aperture: ' + str(abs(minAperture)))
                print('Optimal Annulus: ' + str(minAnnulus))
                print('********************************************\n')

                log.debug('*********************************************')
                log.debug('Best Comparison Star: None')
                log.debug('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                log.debug('Optimal Aperture: ' + str(abs(minAperture)))
                log.debug('Optimal Annulus: ' + str(minAnnulus))
                log.debug('********************************************\n')

            else:
                print('\n*********************************************')
                print('Best Comparison Star: #' + str(bestCompStar))
                print('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                print('Optimal Aperture: ' + str(minAperture))
                print('Optimal Annulus: ' + str(minAnnulus))
                print('********************************************\n')

                log.debug('*********************************************')
                log.debug('Best Comparison Star: #' + str(bestCompStar))
                log.debug('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                log.debug('Optimal Aperture: ' + str(minAperture))
                log.debug('Optimal Annulus: ' + str(minAnnulus))
                log.debug('********************************************\n')

            # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # # Save a text file of the RA and DEC of the target and comp
            # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # outParamsFile = open(saveDirectory + targetName + date + '.radec', 'w+')
            # outParamsFile.write('#RA, Dec, Target = 0 / Ref Star = 1, Centroid [pix]\n')
            # outParamsFile.write(raStr+","+decStr+",0,"+str(minAperture)+"\n")
            # outParamsFile.close()
 
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Save an image of the FOV
            # (for now, take the first image; later will sum all of the images up)
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            imscale = get_pixel_scale(hdulWCS, wcsFile, imageheader, exotic_infoDict['pixel_scale'])
            if hdulWCS:
                hdulWCS.close()  # close stream
                del hdulWCS
            # imwidth = np.shape(sortedallImageData[0])[1]
            # imheight = np.shape(sortedallImageData[0])[0]
            picframe = 10*(minAperture+minAnnulus)
            pltx = [min([finXTargCent[0], finXRefCent[0]])-picframe, max([finXTargCent[0], finXRefCent[0]])+picframe]
            # FORwidth = pltx[1]-pltx[0]
            plty = [min([finYTargCent[0], finYRefCent[0]])-picframe, max([finYTargCent[0], finYRefCent[0]])+picframe]
            # FORheight = plty[1]-plty[0]

            fig, ax = plt_exotic.subplots()
            target_circle = plt_exotic.Circle((finXTargCent[0], finYTargCent[0]), minAperture, color='lime', fill=False, ls='-', label='Target')
            target_circle_sky = plt_exotic.Circle((finXTargCent[0], finYTargCent[0]), minAperture + minAnnulus, color='lime', fill=False, ls='--', lw=.5)
            if minAperture >= 0:
                ref_circle = plt_exotic.Circle((finXRefCent[0], finYRefCent[0]), minAperture, color='r', fill=False, ls='-.', label='Comp')
                ref_circle_sky = plt_exotic.Circle((finXRefCent[0], finYRefCent[0]), minAperture + minAnnulus, color='r', fill=False, ls='--', lw=.5)

            plt_exotic.close()

            for stretch in [LinearStretch(), SquaredStretch(), SqrtStretch(), LogStretch()]:
                med_img = median_filter(sortedallImageData[0], (4, 4))[int(pltx[0]):round(int(pltx[1])),int(plty[0]):round(int(plty[1]))]
                norm = ImageNormalize(sortedallImageData[0], interval=ZScaleInterval(), stretch=stretch)
                plt_exotic.imshow(sortedallImageData[0], norm=norm, origin='lower', cmap='Greys_r', interpolation=None, vmin=np.percentile(med_img, 5), vmax=np.percentile(med_img, 99))
                plt_exotic.plot(finXTargCent[0], finYTargCent[0], marker='+', color='lime')
                ax.add_artist(target_circle)
                ax.add_artist(target_circle_sky)
                if minAperture >= 0:
                    ax.add_artist(ref_circle)
                    ax.add_artist(ref_circle_sky)
                    plt_exotic.plot(finXRefCent[0], finYRefCent[0], '+r')
                plt_exotic.xlabel("x-axis [pixel]")
                plt_exotic.ylabel("y-axis [pixel]")
                plt_exotic.title("FOV for " + pDict['pName'] + "\n(" + imscale + ")")
                plt_exotic.xlim(pltx[0], pltx[1])
                plt_exotic.ylim(plty[0], plty[1])
                ax.grid(False)
                plt_exotic.plot(0, 0, color='lime', ls='-', label='Target')
                if minAperture >= 0:
                    plt_exotic.plot(0, 0, color='r', ls='-.', label='Comp')
                l = plt_exotic.legend(frameon=None, framealpha=0.1)
                for text in l.get_texts():
                    text.set_color("white")
                plt_exotic.savefig(exotic_infoDict['saveplot'] + "FOV" + pDict['pName'] + exotic_infoDict['date'] +"_"+ str(stretch.__class__).split(".")[-1].split("'")[0] + ".pdf", bbox_inches='tight')
                plt_exotic.close()

            print("\nFOV file saved as: " + exotic_infoDict['saveplot'] + "FOV" + pDict['pName'] + exotic_infoDict['date'] + ".pdf\n")
            log.debug("FOV file saved as: " + exotic_infoDict['saveplot'] + "FOV" + pDict['pName'] + exotic_infoDict['date'] + ".pdf\n")

            # Take the BJD times from the image headers
            if "BJD_TDB" in imageheader:
                goodTimes = nonBJDTimes
            # If not in there, then convert all the final times into BJD - using astropy alone
            else:
                print("No BJDs in Image Headers. Converting all JDs to BJD_TDBs.")
                print("Please be patient- this step can take a few minutes.")

                log.debug("No BJDs in Image Headers. Converting all JDs to BJD_TDBs.")
                # targetloc = astropy.coordinates.SkyCoord(raStr, decStr, unit=(astropy.units.deg,astropy.units.deg), frame='icrs')
                # obsloc = astropy.coordinates.EarthLocation(lat=lati, lon=longit)
                # timesToConvert = astropy.time.Time(nonBJDTimes, format='jd', scale='utc', location=obsloc)
                # ltt_bary = timesToConvert.light_travel_time(targetloc)
                # time_barycentre = timesToConvert.tdb + ltt_bary
                # resultos = time_barycentre.value
                # goodTimes = resultos
                animate_toggle(True)
                resultos = utc_tdb.JDUTC_to_BJDTDB(nonBJDTimes, ra=pDict['ra'], dec=pDict['dec'], lat=lati, longi=longit, alt=exotic_infoDict['elev'])
                goodTimes = resultos[0]
                animate_toggle()

            # Centroid position plots
            plotCentroids(finXTargCent, finYTargCent, finXRefCent, finYRefCent, goodTimes, pDict['pName'], exotic_infoDict['date'])

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

            # another 3 sigma clip based on residuals of the best LM fit
            try:
                interFilter = sigma_clip(goodResids, sigma=3)
            except TypeError:
                interFilter = sigma_clip(goodResids, sigma=3)

            goodFluxes = goodFluxes[~interFilter]
            goodTimes = goodTimes[~interFilter]
            goodPhases = goodPhases[~interFilter]
            goodAirmasses = goodAirmasses[~interFilter]
            goodTargets = goodTargets[~interFilter]
            goodReferences = goodReferences[~interFilter]
            goodTUnc = goodTUnc[~interFilter]
            goodRUnc = goodRUnc[~interFilter]
            goodNormUnc = goodNormUnc[~interFilter]

            # Calculate the standard deviation of the normalized flux values
            standardDev1 = np.std(goodFluxes)

            ######################################
            # PLOTS ROUND 1
            ####################################
            # Make plots of raw target and reference values
            plt_exotic.figure()
            plt_exotic.errorbar(goodTimes, goodTargets, yerr=goodTUnc, linestyle='None', fmt='-o')
            plt_exotic.xlabel('Time (BJD)')
            plt_exotic.ylabel('Total Flux')
            # plt_exotic.rc('grid', linestyle="-", color='black')
            # plt_exotic.grid(True)
            plt_exotic.title(pDict['pName'] + ' Raw Flux Values ' + exotic_infoDict['date'])
            plt_exotic.savefig(exotic_infoDict['saveplot'] + 'temp/TargetRawFlux' + pDict['pName'] + exotic_infoDict['date'] + '.png')
            plt_exotic.close()

            plt_exotic.figure()
            plt_exotic.errorbar(goodTimes, goodReferences, yerr=goodRUnc, linestyle='None', fmt='-o')
            plt_exotic.xlabel('Time (BJD)')
            plt_exotic.ylabel('Total Flux')
            # plt_exotic.rc('grid', linestyle="-", color='black')
            # plt_exotic.grid(True)
            plt_exotic.title('Comparison Star Raw Flux Values ' + exotic_infoDict['date'])
            plt_exotic.savefig(exotic_infoDict['saveplot'] + 'temp/CompRawFlux' + pDict['pName'] + exotic_infoDict['date'] + '.png')
            plt_exotic.close()

            # Plots final reduced light curve (after the 3 sigma clip)
            plt_exotic.figure()
            plt_exotic.errorbar(goodTimes, goodFluxes, yerr=goodNormUnc, linestyle='None', fmt='-bo')
            plt_exotic.xlabel('Time (BJD)')
            plt_exotic.ylabel('Normalized Flux')
            # plt_exotic.rc('grid', linestyle="-", color='black')
            # plt_exotic.grid(True)
            plt_exotic.title(pDict['pName'] + ' Normalized Flux vs. Time ' + exotic_infoDict['date'])
            plt_exotic.savefig(exotic_infoDict['saveplot'] + 'NormalizedFluxTime' + pDict['pName'] + exotic_infoDict['date'] + '.png')
            plt_exotic.close()

            # Save normalized flux to text file prior to MCMC
            outParamsFile = open(exotic_infoDict['saveplot'] + 'NormalizedFlux' + pDict['pName'] + exotic_infoDict['date'] + '.txt', 'w+')
            outParamsFile.write(str("BJD") + ',' + str("Norm Flux") + ',' + str("Norm Err") + ',' + str("AM") + '\n')
            for ti, fi, erri, ami in zip(goodTimes, goodFluxes, goodNormUnc, goodAirmasses):
                outParamsFile.write(str(round(ti, 8)) + ',' + str(round(fi, 7)) + ',' + str(round(erri, 6)) + ',' + str(round(ami, 2)) + '\n')
            # CODE YIELDED DATA IN PREV LINE FORMAT
            outParamsFile.close()
            print('\nOutput File Saved')
            log.debug('Output File Saved')
        else:
            goodTimes, goodFluxes, goodNormUnc, goodAirmasses = [], [], [], []
            for i in processeddata:
                try:
                    goodTimes.append(float(i.split(",")[0]))
                    goodFluxes.append(float(i.split(",")[1]))
                    goodNormUnc.append(float(i.split(",")[2]))
                    goodAirmasses.append(float(i.split(",")[3]))
                except ValueError:
                    continue

            goodTimes = np.array(goodTimes)
            goodFluxes = np.array(goodFluxes)
            goodNormUnc = np.array(goodNormUnc)
            goodAirmasses = np.array(goodAirmasses)

            # Ask user for time format and convert it if not in BJD_TDB
            validTimeFormats = ['BJD_TDB', "MJD_UTC", "JD_UTC"]
            formatEntered = False
            print("\nNOTE: If your file is not in one of the following formats, please re-reduce your data into one of the time formats recognized by EXOTIC.")
            log.debug("NOTE: If your file is not in one of the following formats, please re-reduce your data into one of the time formats recognized by EXOTIC.")
            while not formatEntered:
                print("Which of the following time formats is your data file stored in? (Type q to quit)")
                timeFormat = str(input("BJD_TDB / JD_UTC / MJD_UTC: "))
                log.debug("Which of the following time formats is your data file stored in? (Type q to quit) "+timeFormat)
                if (timeFormat.upper()).strip() == 'Q':
                    sys.exit()
                elif (timeFormat.upper()).strip() not in validTimeFormats:
                    print("\nInvalid entry; please try again.")
                    log.debug("Invalid entry; please try again.")
                else:
                    formatEntered = True
            timeFormat = (timeFormat.upper()).strip()
            goodTimes = timeConvert(goodTimes, timeFormat, pDict, exotic_infoDict)

            #Ask user for flux units and convert to flux if in magnitude/millimagnitude
            validFluxFormats = ['flux', "magnitude", "millimagnitude"]
            formatEntered = False
            print("\nNOTE: If your file is not in one of the following formats, please rereduce your data into one of the time formats recognized by EXOTIC.")
            log.debug("NOTE: If your file is not in one of the following formats, please re-reduce your data into one of the time formats recognized by EXOTIC.")
            while not formatEntered:
                print("Which of the following units of flux is your data file stored in? (Type q to quit)")
                fluxFormat = str(input("flux / magnitude / millimagnitude: "))
                log.debug("Which of the following units of flux is your data file stored in? (Type q to quit)"+fluxFormat)
                if (fluxFormat.upper()).strip() == 'Q':
                    sys.exit()
                elif (fluxFormat.lower()).strip() not in validFluxFormats:
                    print("\nInvalid entry; please try again.")
                    log.debug("Invalid entry; please try again.")
                else:
                    formatEntered = True
            fluxFormat = (fluxFormat.lower()).strip()
            if fluxFormat != "flux":
                goodFluxes, goodNormUnc = fluxConvert(goodFluxes, goodNormUnc, fluxFormat)

            bjdMidTOld = goodTimes[0]
            standardDev1 = np.std(goodFluxes)

        print('\n****************************************')
        print('Fitting a Light Curve Model to Your Data')
        print('****************************************\n')

        log.debug('****************************************')
        log.debug('Fitting a Light Curve Model to Your Data')
        log.debug('****************************************\n')

        ##########################
        # NESTED SAMPLING FITTING
        ##########################

        try:
            prior = {
                'rprs': pDict['rprs'],    # Rp/Rs
                'ars': pDict['aRs'],      # a/Rs
                'per': pDict['pPer'],     # Period [day]
                'inc': pDict['inc'],      # Inclination [deg]
                'u0': ld0[0], 'u1': ld1[0], 'u2': ld2[0], 'u3': ld3[0],  # limb darkening (nonlinear)
                'ecc': pDict['ecc'],     # Eccentricity
                'omega': 0,          # Arg of periastron
                'tmid': pDict['midT'],    # time of mid transit [day]
                'a1': bestlmfit.parameters['a1'],  # mid Flux
                'a2': bestlmfit.parameters['a2'],  # Flux lower bound
            }
        except:
            prior = {
                'rprs': pDict['rprs'],    # Rp/Rs
                'ars': pDict['aRs'],      # a/Rs
                'per': pDict['pPer'],     # Period [day]
                'inc': pDict['inc'],      # Inclination [deg]
                'u0': ld0[0], 'u1': ld1[0], 'u2': ld2[0], 'u3': ld3[0],  # limb darkening (nonlinear)
                'ecc': pDict['ecc'],     # Eccentricity
                'omega': 0,          # Arg of periastron
                'tmid': pDict['midT'],    # time of mid transit [day]
                'a1': goodFluxes.mean(),  # max() - arrayFinalFlux.min(), #mid Flux
                'a2': 0,             # Flux lower bound
            }

        phase = (goodTimes-prior['tmid'])/prior['per']
        prior['tmid'] = pDict['midT'] + np.floor(phase).max()*prior['per']
        upper = pDict['midT'] + 25*pDict['midTUnc'] + np.floor(phase).max()*(pDict['pPer']+25*pDict['pPerUnc'])
        lower = pDict['midT'] - 25*pDict['midTUnc'] + np.floor(phase).max()*(pDict['pPer']-25*pDict['pPerUnc'])

        if np.floor(phase).max()-np.floor(phase).min() == 0:
            print('ERROR: Estimated mid-transit not in observation range (check priors or observation time)')
            print('start:', goodTimes.min())
            print('  end:', goodTimes.max())
            print('prior:', prior['tmid'])

            log.debug('ERROR: Estimated mid-transit not in observation range (check priors or observation time)')
            log.debug('start:', goodTimes.min())
            log.debug('  end:', goodTimes.max())
            log.debug('prior:', prior['tmid'])

        try:
            mybounds = {
                'rprs': [pDict['rprs']-3*pDict['rprsUnc'], pDict['rprs']+3*pDict['rprsUnc']],
                'tmid': [max(lower, goodTimes.min()), min(goodTimes.max(), upper)],
                'ars': [pDict['aRs']-5*pDict['aRsUnc'], pDict['aRs']+5*pDict['aRsUnc']],

                'a1': [bestlmfit.parameters['a1']*0.75, bestlmfit.parameters['a1']*1.25],
                'a2': [bestlmfit.parameters['a2']-0.25, bestlmfit.parameters['a2']+0.25],
            }
        except:
            mybounds = {
                'rprs': [pDict['rprs']-3*pDict['rprsUnc'], pDict['rprs']+3*pDict['rprsUnc']],
                'tmid': [max(lower, goodTimes.min()), min(goodTimes.max(), upper)],
                'ars': [pDict['aRs']-5*pDict['aRsUnc'], pDict['aRs']+5*pDict['aRsUnc']],
                'a1': [ min(0,np.nanmin(goodFluxes)), 3*np.nanmax(goodFluxes)],
                'a2': [-3, 3],
            }

        # fitting method in elca.py
        myfit = lc_fitter(goodTimes, goodFluxes, goodNormUnc, goodAirmasses, prior, mybounds, mode='ns')

        # for k in myfit.bounds.keys():
        #     print("{:.6f} +- {}".format( myfit.parameters[k], myfit.errors[k]))

        ########################
        # PLOT FINAL LIGHT CURVE
        ########################
        f, (ax_lc, ax_res) = plt_exotic.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

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
        ax_res.errorbar(myfit.phase, myfit.residuals/np.median(myfit.data), yerr=myfit.detrendederr, color='gray', marker='o', markersize=5, linestyle='None', mec='None', alpha=0.75)
        ax_res.plot(myfit.phase, np.zeros(len(myfit.phase)), 'r-', lw=2, alpha=1, zorder=100)
        ax_res.set_ylabel('Residuals')
        ax_res.set_ylim([-3 * np.nanstd(myfit.residuals/np.median(myfit.data)), 3 * np.nanstd(myfit.residuals/np.median(myfit.data))])

        correctedSTD = np.std(myfit.residuals/np.median(myfit.data))
        ax_lc.errorbar(myfit.phase, myfit.detrended, yerr=myfit.detrendederr, ls='none',
                       marker='o', color='gray', markersize=5, mec='None', alpha=0.75)
        ax_lc.plot(myfit.phase, myfit.transit, 'r', zorder=1000, lw=2)

        ax_lc.set_ylabel('Relative Flux')
        ax_lc.get_xaxis().set_visible(False)

        ax_res.errorbar(binner(myfit.phase, len(myfit.residuals) // 10), binner(myfit.residuals/np.median(myfit.data), len(myfit.residuals) // 10),
                        yerr=binner(myfit.residuals/np.median(myfit.data), len(myfit.residuals) // 10, myfit.detrendederr)[1],
                        fmt='s', ms=5, mfc='b', mec='None', ecolor='b', zorder=10)
        ax_lc.errorbar(binner(myfit.phase, len(myfit.phase) // 10),
                       binner(myfit.detrended, len(myfit.detrended) // 10),
                       yerr=binner(myfit.residuals/np.median(myfit.data), len(myfit.residuals) // 10, myfit.detrendederr)[1],
                       fmt='s', ms=5, mfc='b', mec='None', ecolor='b', zorder=10)

        # remove vertical whitespace
        f.subplots_adjust(hspace=0)

        # For some reason, saving as a pdf crashed on Rob's laptop...so adding in a try statement to save it as a pdf if it can, otherwise, png
        try:
            f.savefig(exotic_infoDict['saveplot'] + 'FinalLightCurve' + pDict['pName'] + exotic_infoDict['date'] + ".pdf", bbox_inches="tight")
        except AttributeError:
            f.savefig(exotic_infoDict['saveplot'] + 'FinalLightCurve' + pDict['pName'] + exotic_infoDict['date'] + ".png", bbox_inches="tight")
        plt_exotic.close()

        ###################################################################################

        # triangle plot
        fig, axs = dynesty.plotting.cornerplot(myfit.results, labels=list(mybounds.keys()), quantiles_2d=[0.4, 0.85], smooth=0.015, show_titles=True, use_math_text=True, title_fmt='.2e', hist2d_kwargs={'alpha': 1, 'zorder': 2, 'fill_contours': False})
        dynesty.plotting.cornerpoints(myfit.results, labels=list(mybounds.keys()), fig=[fig, axs[1:, :-1]], plot_kwargs={'alpha': 0.1, 'zorder': 1, })
        fig.savefig(exotic_infoDict['saveplot'] + 'temp/Triangle_{}_{}.png'.format(pDict['pName'], exotic_infoDict['date']))
        plt_exotic.close()

        # write output to text file
        outParamsFile = open(exotic_infoDict['saveplot'] + 'FinalLightCurve' + pDict['pName'] + exotic_infoDict['date'] + '.csv', 'w+')
        outParamsFile.write('# FINAL TIMESERIES OF ' + pDict['pName'] + '\n')
        outParamsFile.write('# BJD_TDB,Orbital Phase,Flux,Uncertainty,Model,Airmass\n')

        phase = getPhase(myfit.time, pDict['pPer'], myfit.parameters['tmid'])

        for bjdi, phasei, fluxi, fluxerri, modeli, ami in zip(myfit.time, phase, myfit.detrended, myfit.dataerr/myfit.airmass_model, myfit.transit, myfit.airmass_model):

            outParamsFile.write("{}, {}, {}, {}, {}, {}\n".format(bjdi, phasei, fluxi, fluxerri, modeli, ami))

        outParamsFile.close()

        #######################################################################
        # print final extracted planetary parameters
        #######################################################################

        log.info('*********************************************************')
        log.info('FINAL PLANETARY PARAMETERS\n')
        log.info('              Mid-Transit Time [BJD]: {} +- {} '.format(round_to_2(myfit.parameters['tmid'], myfit.errors['tmid']), round_to_2(myfit.errors['tmid'])))
        log.info('  Radius Ratio (Planet/Star) [Rp/Rs]: {} +- {} '.format(round_to_2(myfit.parameters['rprs'], myfit.errors['rprs']), round_to_2(myfit.errors['rprs'])))
        log.info(' Semi Major Axis/ Star Radius [a/Rs]: {} +- {} '.format(round_to_2(myfit.parameters['ars'], myfit.errors['ars']), round_to_2(myfit.errors['ars'])))
        log.info('               Airmass coefficient 1: {} +- {} '.format(round_to_2(myfit.parameters['a1'], myfit.errors['a1']), round_to_2(myfit.errors['a1'])))
        log.info('               Airmass coefficient 2: {} +- {} '.format(round_to_2(myfit.parameters['a2'], myfit.errors['a2']), round_to_2(myfit.errors['a2'])))
        log.info('The scatter in the residuals of the lightcurve fit is: {} %'.format(round_to_2(100. * np.std(myfit.residuals/np.median(myfit.data)))))
        log.info('*********************************************************')

        ##########
        # SAVE DATA
        ##########

        # TODO write as json
        # write output to text file
        outParamsFile = open(exotic_infoDict['saveplot'] + 'FinalParams' + pDict['pName'] + exotic_infoDict['date'] + '.txt', 'w+')
        outParamsFile.write('FINAL PLANETARY PARAMETERS\n')
        outParamsFile.write('')
        outParamsFile.write(' Mid-Transit Time: ' + str(round_to_2(myfit.parameters['tmid'], myfit.errors['tmid'])) + ' +/- ' + str(round_to_2(myfit.errors['tmid'])) + ' (BJD)\n')
        outParamsFile.write(' Ratio of Planet to Stellar Radius: ' + str(round_to_2(myfit.parameters['rprs'], myfit.errors['rprs'])) + ' +/- ' + str(
            round_to_2(myfit.errors['rprs'])) + ' (Rp/Rs)\n')
        outParamsFile.write(' transit depth uncertainty: ' + str(round_to_2(100. * 2. * myfit.parameters['rprs'] * myfit.errors['rprs'])) + ' (%)\n')
        outParamsFile.write(' airmass coefficient 1: ' + str(round_to_2(myfit.parameters['a1'], myfit.errors['a1'])) + ' +/- ' + str(round_to_2(myfit.errors['a1'])) + '\n')
        outParamsFile.write(' airmass coefficient 2: ' + str(round_to_2(myfit.parameters['a2'], myfit.errors['a2'])) + ' +/- ' + str(round_to_2(myfit.errors['a2'])) + '\n')
        outParamsFile.write(' scatter in the residuals of the lightcurve fit is: ' + str(round_to_2(100. * np.std(myfit.residuals/np.median(myfit.data)))) + '%\n')
        outParamsFile.close()
        print('\nFinal Planetary Parameters have been saved in ' + exotic_infoDict['saveplot'] + ' as '
              + pDict['pName'] + exotic_infoDict['date'] + '.txt' + '\n')

        # AAVSO Format
        userCode = exotic_infoDict['aavsonum']
        secuserCode = exotic_infoDict['secondobs']
        # else:
        outParamsFile = open(exotic_infoDict['saveplot'] + 'AAVSO' + pDict['pName'] + exotic_infoDict['date'] + '.txt', 'w+')
        outParamsFile.write('#TYPE=EXOPLANET\n')  # fixed
        outParamsFile.write('#OBSCODE=' + exotic_infoDict['aavsonum'] + '\n')  # UI
        outParamsFile.write('#SECONDARY_OBSCODE=' + exotic_infoDict['secondobs'] + '\n')  # UI
        outParamsFile.write('#SOFTWARE=EXOTIC v' + __version__ + '\n')  # fixed
        outParamsFile.write('#DELIM=,\n')  # fixed
        outParamsFile.write('#DATE_TYPE=BJD_TDB\n')  # fixed
        outParamsFile.write('#OBSTYPE=' + exotic_infoDict['ctype'] + '\n')
        outParamsFile.write('#STAR_NAME=' + pDict['sName'] + '\n')  # code yields
        outParamsFile.write('#EXOPLANET_NAME=' + pDict['pName'] + '\n')  # code yields
        outParamsFile.write('#BINNING=' + exotic_infoDict['pixelbin'] + '\n')  # user input
        outParamsFile.write('#EXPOSURE_TIME=' + str(exotic_infoDict['exposure']) + '\n')  # UI
        outParamsFile.write('#FILTER=' + exotic_infoDict['filter'] + '\n')
        outParamsFile.write('#NOTES=' + exotic_infoDict['notes'] + '\n')
        outParamsFile.write('#DETREND_PARAMETERS=AIRMASS, AIRMASS CORRECTION FUNCTION\n')  # fixed
        outParamsFile.write('#MEASUREMENT_TYPE=Rnflux\n')  # fixed
        # outParamsFile.write('#PRIORS=Period=' + str(planetPeriod) + ' +/- ' + str(ogPeriodErr) + ',a/R*=' + str(
        #     semi) + ',Tc=' + str(round(bjdMidTranCur, 8)) + ' +/- ' + str(round(propMidTUnct, 8)) + ',T0=' + str(
        #     round(bjdMidTOld, 8)) + ' +/- ' + str(round(ogMidTErr, 8)) + ',inc=' + str(inc) + ',ecc=' + str(
        #     eccent) + ',u1=' + str(linearLimb) + ',u2=' + str(quadLimb) + '\n')  # code yields
        outParamsFile.write('#PRIORS=Period=' + str(pDict['pPer']) + ' +/- ' + str(pDict['pPerUnc']) + ',a/R*='
                            + str(pDict['aRs']) + ',inc=' + str(pDict['inc']) + ',ecc=' + str(pDict['ecc']) + ',u0='
                            + str(ld0[0]) + ' +/- ' + str(ld0[1]) + ',u1=' + str(ld1[0]) + ' +/- ' + str(ld1[1])
                            + ',u2=' + str(ld2[0]) + ' +/- ' + str(ld2[1]) + ',u3=' + str(ld3[0]) + ' +/- '
                            + str(ld3[1]) + '\n')
        # code yields
        outParamsFile.write(
            '#RESULTS=Tc=' + str(round(myfit.parameters['tmid'], 8)) + ' +/- ' + str(round(myfit.errors['tmid'], 8)) + ',Rp/R*=' + str(
                round(myfit.parameters['rprs'], 6)) + ' +/- ' + str(round(myfit.errors['rprs'], 6)) + ',Am1=' + str(
                round(myfit.parameters['a1'], 5)) + ' +/- ' + str(round(myfit.errors['a1'], 5)) + ',Am2=' + str(
                round(myfit.parameters['a2'], 5)) + ' +/- ' + str(round(myfit.errors['a2'], 5)) + '\n')  # code yields
        # outParamsFile.write('#NOTES= ' + userNameEmails + '\n')
        outParamsFile.write('#DATE,FLUX,MERR,DETREND_1,DETREND_2\n')
        for aavsoC in range(0, len(myfit.time)):
            outParamsFile.write(
                str(round(myfit.time[aavsoC], 8)) + ',' + str(round(myfit.data[aavsoC], 7)) + ',' + str(
                    round(myfit.dataerr[aavsoC], 7)) + ',' + str(round(goodAirmasses[aavsoC], 7)) + ',' + str(
                    round(myfit.airmass_model[aavsoC], 7)) + '\n')

        # CODE YIELDED DATA IN PREV LINE FORMAT
        outParamsFile.close()
        print('Output File Saved')

        log.info('************************')
        log.info('End of Reduction Process')
        log.info('************************')

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
