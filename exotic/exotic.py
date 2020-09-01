# Copyright (c) 2019-2020, California Institute of Technology.
# All rights reserved.  Based on Government Sponsored Research under contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#    1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#    3. Neither the name of the California Institute of Technology (Caltech), its operating division the Jet Propulsion Laboratory (JPL), the National Aeronautics and Space Administration (NASA), nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE CALIFORNIA INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


####################################################################
# EXOplanet Transit Interpretation Code (EXOTIC)
#
# Authors: Ethan Blaser, Rob Zellem, Kyle Pearson, Tamim Fatahi, Marlee Smith, Aaron Tran, John Engelke, Sujay Nair, Jon Varghese, Michael Fitzgerald
# Supplemental Code: Kyle Pearson, Gael Roudier, and Jason Eastman
####################################################################

# --IMPORTS -----------------------------------------------------------
print("Importing Python Packages - please wait.")

# EXOTIC version number
# Now adhering to the Semantic Versioning 2.0.0
# Given a version number MAJOR.MINOR.PATCH, increment the:
# MAJOR version when you make incompatible API changes,
# MINOR version when you add functionality in a backwards compatible manner, and
# PATCH version when you make backwards compatible bug fixes.
# Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.
# https://semver.org, e.g. __version__ = "0.14.4" from the version import
try:  # module import
    from .version import __version__
except ImportError:  # package import
    from version import __version__

import itertools
import threading
import time
import sys

print('Python Version: %s' % sys.version)

# To increase memory allocation for EXOTIC; allows for more fits files
# import resource
# resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))


# here is the animation
def animate():
    for c in itertools.cycle(['|', '/', '-', '\\']):
        if done:
            break
        sys.stdout.write('\rThinking ' + c)
        sys.stdout.flush()
        time.sleep(0.1)
    # sys.stdout.write('\nDone!     \n')
    # print('\nDone!     \n')


done = False
t = threading.Thread(target=animate, daemon=True)
t.start()

import os
import json
import logging
import platform
import warnings
import argparse
import glob as g
from io import StringIO

# data processing
import pandas
import requests
import numpy as np

# julian conversion imports
import dateutil.parser as dup

# UTC to BJD converter import
from barycorrpy import utc_tdb

# Curve fitting imports
from scipy.optimize import least_squares
from scipy.ndimage import median_filter
from scipy.stats import mode

# Pyplot imports
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from matplotlib.animation import FuncAnimation

plt.style.use(astropy_mpl_style)

# Nested Sampling imports
import dynesty
try:  # module import
    from .api.elca import lc_fitter, binner
except ImportError:  # package import
    from api.elca import lc_fitter, binner

# astropy imports
import astropy.time
import astropy.coordinates
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.wcs import WCS
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astropy.utils.exceptions import AstropyWarning

# Image alignment import
import astroalign as aa

# Nonlinear Limb Darkening Calculations import
try:  # module import
    from .api.gaelLDNL import createldgrid
except ImportError:  # package import
    from api.gaelLDNL import createldgrid

# photometry
from photutils import CircularAperture
from photutils import aperture_photometry

# cross corrolation imports
from skimage.registration import phase_cross_correlation

# error handling for scraper
from tenacity import retry, retry_if_exception_type, stop_after_attempt, \
    wait_exponential, wait_random

# long process here
# time.sleep(10)
done = True

# ################### START PROPERTIES ########################################
# GLOBALS (set in main before method calls)
infoDict = dict()
UIprevTPX, UIprevTPY, UIprevRPX, UIprevRPY = 0, 0, 0, 0
distFC = 0
ax1 = plt.figure()  # placeholder

# CONFIGURATIONS
requests_timeout = 16, 512  # connection timeout, response timeout in secs.

# SHARED CONSTANTS
pi = 3.14159
au = 1.496e11  # m
rsun = 6.955e8  # m
rjup = 7.1492e7  # m
G = 0.00029591220828559104  # day, AU, Msun

# SHARED LAMBDAS
# keplerian semi-major axis (au)
sa = lambda m, P: (G*m*P**2/(4*pi**2))**(1./3)
# ################### END PROPERTIES ##########################################


def sigma_clip(ogdata, sigma=3, dt=20):
    mdata = median_filter(ogdata, dt)
    res = ogdata - mdata
    std = np.nanmedian([np.nanstd(np.random.choice(res,50)) for i in range(100)])
    #std = np.nanstd(res) # biased from large outliers
    return np.abs(res) > sigma*std


# ################### START ARCHIVE SCRAPER (PRIORS) ##########################
def dataframe_to_jsonfile(dataframe, filename):
    jsondata = json.loads(dataframe.to_json(orient='table', index=False))
    with open(filename, "w") as f:
        f.write(json.dumps(jsondata['data'], indent=4))


def tap_query(base_url, query, dataframe=True):
    # table access protocol query

    # build url
    uri_full = base_url
    for k in query:
        if k != "format":
            uri_full += "{} {} ".format(k, query[k])

    uri_full = uri_full[:-1] + "&format={}".format(query.get("format", "csv"))
    uri_full = uri_full.replace(' ', '+')
    print(uri_full)

    response = requests.get(uri_full, timeout=requests_timeout)
    # TODO check status_code?

    if dataframe:
        return pandas.read_csv(StringIO(response.text))
    else:
        return response.text

@retry(stop=stop_after_attempt(3),
       wait=wait_exponential(multiplier=1, min=17, max=1024),
       retry=(retry_if_exception_type(requests.exceptions.RequestException) |
       retry_if_exception_type(ConnectionError)))
def new_scrape(filename="eaConf.json", target=None):

    # scrape_new()
    uri_ipac_base = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
    uri_ipac_query = {
        # Table columns: https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
        "select"   : "pl_name,hostname,tran_flag,pl_massj,pl_radj,pl_radjerr1,"
                     "pl_ratdor,pl_ratdorerr1,pl_orbincl,pl_orbinclerr1,"
                     "pl_orbper,pl_orbpererr1,pl_orbeccen,"
                     "pl_orblper,pl_tranmid,pl_tranmiderr1,"
                     "pl_trandep,pl_trandeperr1,pl_trandeperr2,"
                     "pl_ratror,pl_ratrorerr1,pl_ratrorerr2,"
                     "st_teff,st_tefferr1,st_tefferr2,st_met,st_meterr1,st_meterr2,"
                     "st_logg,st_loggerr1,st_loggerr2,st_mass,st_rad,st_raderr1,ra,dec,pl_pubdate",
        "from"     : "ps",  # Table name
        "where"    : "tran_flag = 1 and default_flag = 1",
        "order by" : "pl_pubdate desc",
        "format"   : "csv"
    }

    if target:
        uri_ipac_query["where"] += " and pl_name = '{}'".format(target)

    default = tap_query(uri_ipac_base, uri_ipac_query)

    # fill in missing columns
    uri_ipac_query['where'] = 'tran_flag=1'

    if target:
        uri_ipac_query["where"] += " and pl_name = '{}'".format(target)

    extra = tap_query(uri_ipac_base, uri_ipac_query)

    if len(default) == 0:
        target = input("Cannot find target ({}) in NASA Exoplanet Archive. Check case sensitivity and re-enter the"
                       "\nplanet's name or type candidate if this is a planet candidate: ".format(target))
        if target.strip().lower() == 'candidate':
            target = input("\nPlease enter candidate planet's name: ")
            return target, True
        else:
            return new_scrape(filename="eaConf.json", target=target)
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
                        # print(i,k,edata[k][edata[k].notna()].values[0])
                    else:
                        # permanent nans - require manual entry
                        if k == 'pl_orblper':  # omega
                            default.loc[default.pl_name == i, k] = 0
                        elif k == 'pl_ratdor':  # a/R*
                            # Kepler's 3rd law
                            semi = sa(ddata.st_mass.values[0], ddata.pl_orbper.values[0])
                            default.loc[default.pl_name == i, k] = semi*au / (ddata.st_rad.values[0]*rsun)
                        elif k == 'pl_orbincl':  # inclination
                            default.loc[default.pl_name == i, k] = 90
                        elif k == "pl_orbeccen":  # eccentricity
                            default.loc[default.pl_name == i, k] = 0
                        elif k == "st_met":  # [Fe/H]
                            default.loc[default.pl_name == i, k] = 0

        dataframe_to_jsonfile(default, filename)
        return target, False


def new_getParams(data):
    # translate data from Archive keys to Ethan Keys
    try:
        rprs = np.sqrt(data['pl_trandep']/100.)
        rprserr = np.sqrt(np.abs((data['pl_trandeperr1']/100.)*(data['pl_trandeperr2']/100.)))/(2.*rprs)
    except (KeyError, TypeError):
        try:
            rprs = data['pl_ratror']
            rprserr = np.sqrt(np.abs(data['pl_ratrorerr1']*data['pl_ratrorerr2']))
        except (KeyError, TypeError):
            rp = data['pl_radj']*rjup
            rperr = data['pl_radjerr1']*rjup
            rs = data['st_rad']*rsun
            rserr = data['st_raderr1']*rsun
            rprserr = ((rperr/rs)**2 + (-rp*rserr/rs**2)**2 )**0.5
            rprs = rp/rs
    planetDictionary = {
        'ra': data['ra'],
        'dec': data['dec'],
        'pName': data['pl_name'],
        'sName': data['hostname'],
        'pPer': data['pl_orbper'],
        'pPerUnc': data['pl_orbpererr1'],

        'midT': data['pl_tranmid'],
        'midTUnc': data['pl_tranmiderr1'],
        'rprs': rprs,
        'rprsUnc': rprserr,
        'aRs': data['pl_ratdor'],
        'aRsUnc': data['pl_ratdorerr1'],
        'inc': data['pl_orbincl'],
        'incUnc': data['pl_orbinclerr1'],

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

    return planetDictionary
# ################### END ARCHIVE SCRAPER (PRIORS) ############################

#Get Julian time, don't need to divide by 2 since assume mid-EXPOSURE
#Find separate funciton in code that does julian conversion to BJD_TDB

# Method that gets and returns the julian time of the observation
def getJulianTime(hdul):
    exptime_offset = 0
    imageheader = hdul[0].header

    exp = imageheader.get('EXPTIME')  #checking for variation in .fits header format
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
        time = astropy.time.Time(dt)
        julianTime = time.jd
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
        time = astropy.time.Time(dt)
        julianTime = time.jd
        # If the time is from the beginning of the observation, then need to calculate mid-exposure time
        if "start" in hdul[0].header.comments['DATE-OBS']:
            exptime_offset = exp / 2. / 60. / 60. / 24.  # assume exptime is in seconds for now

    # If the mid-exposure time is given in the fits header, then no offset is needed to calculate the mid-exposure time
    return julianTime + exptime_offset


# Method that gets and returns the current phase of the target
def getPhase(curTime, pPeriod, tMid):
    return ((curTime - tMid + 0.5*pPeriod) / pPeriod) % 1


# Method that gets and returns the airmass from the fits file (Really the Altitude)
def getAirMass(hdul, ra, dec, lati, longit, elevation):
    # Grab airmass from image header; if not listed, calculate it from TELALT; if that isn't listed, then calculate it the hard way
    if 'AIRMASS' in hdul[0].header:
        am = float(hdul[0].header['AIRMASS'])
    elif 'TELALT' in hdul[0].header:
        alt = float(hdul[0].header['TELALT'])  # gets the airmass from the fits file header in (sec(z)) (Secant of the zenith angle)
        cosam = np.cos((np.pi / 180) * (90.0 - alt))
        am = 1 / cosam
    else:
        # pointing = SkyCoord(str(astropy.coordinates.Angle(raStr+" hours").deg)+" "+str(astropy.coordinates.Angle(decStr+" degrees").deg ), unit=(u.deg, u.deg), frame='icrs')
        pointing = SkyCoord(str(ra)+" "+str(dec), unit=(u.deg, u.deg), frame='icrs')

        location = EarthLocation.from_geodetic(lat=lati*u.deg, lon=longit*u.deg, height=elevation)
        time = astropy.time.Time(getJulianTime(hdul), format='jd', scale='utc', location=location)
        pointingAltAz = pointing.transform_to(AltAz(obstime=time, location=location))
        am = float(pointingAltAz.secz)
    return am


# Validate user input
def user_input(prompt, type_, val1=None, val2=None, val3=None):
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
        elif type_ == int and val1 and val2:
            if option not in (val1, val2):
                print("Sorry, your response was not valid.")
            else:
                return option
        else:
            return option


# Create a save directory within the current working directory
def create_directory(directory):
    while True:
        try:
            directoryname = input('Enter the name for your new directory: ')
            savedirectory = os.path.join(directory + directoryname, '')
            os.mkdir(savedirectory)
        except OSError:
            print("Creation of the directory %s failed" % savedirectory)
        else:
            print("Successfully created the directory %s " % savedirectory)
            return savedirectory


# Check user's inits.json for user information and planetary parameters
def inits_file(initspath, dictinfo, dictplanet):
    with open(initspath) as json_file:
        data = json.load(json_file)

    initdictinfo = data['user_info']
    initdictplanet = data['planetary_parameters']

    for key, key2 in zip(dictinfo, initdictinfo[0]):
        dictinfo[key] = initdictinfo[0][key2]

    for key, key2 in zip(dictplanet, initdictplanet[0]):
        dictplanet[key] = initdictplanet[0][key2]

    return dictinfo, dictplanet


# Get planet parameters from NEA if not a candidate
def get_planet_name(targetname):

    planetdict = None
    print("\nLooking up {}- please wait.".format(targetname))
    # check to make sure the target can be found in the exoplanet archive right after they enter its name
    targetname, candidateplanet = new_scrape(filename="eaConf.json", target=targetname)
    if not candidateplanet:
        with open("eaConf.json", "r") as confirmedFile:
            data = json.load(confirmedFile)
            planets = [data[i]['pl_name'] for i in range(len(data))]
            idx = planets.index(targetname)
            planetdict = new_getParams(data[idx])
            print('\nSuccessfully found {} in the NASA Exoplanet Archive!'.format(targetname))

    return targetname, candidateplanet, planetdict

#Convert time units to BJD_TDB if pre-reduced file not in proper units
def timeConvert(timeList, timeFormat, pDict, info_dict):
    #if timeFormat is already in BJD_TDB, just do nothing
    #Perform appropriate conversion for each time format if needed
    if timeFormat == "JD_UTC":
        convertedTimes = utc_tdb.JDUTC_to_BJDTDB(timeList, ra=pDict['ra'], dec=pDict['dec'], lat=info_dict['lat'], longi=info_dict['long'], alt=info_dict['elev'])
        timeList = convertedTimes[0]
    elif timeFormat == "MJD_UTC":
        convertedTimes = utc_tdb.JDUTC_to_BJDTDB(timeList + 2400000.5, ra=pDict['ra'], dec=pDict['dec'], lat=info_dict['lat'], longi=info_dict['long'], alt=info_dict['elev'])
        timeList = convertedTimes[0]
    return timeList

#Convert magnitude units to flux if pre-reduced file not in flux already
def fluxConvert(fluxList, errorList, fluxFormat):
    #If units already in flux, do nothing, perform appropriate conversions to flux otherwise
    if fluxFormat == "magnitude":
        convertedPositiveErrors = 10. ** ((-1. * (fluxList + errorList)) / 2.5)
        convertedNegativeErrors = 10. ** ((-1. * (fluxList - errorList)) / 2.5)
        fluxList = 10. ** ((-1. * fluxList) / 2.5)
    if fluxFormat == "millimagnitude":
        convertedPositiveErrors = 10. ** ((-1. * ((fluxList + errorList) / 1000.)) / 2.5)
        convertedNegativeErrors = 10. ** ((-1. * ((fluxList - errorList) / 1000.)) / 2.5)
        fluxList = 10. ** ((-1. * (fluxList / 1000.) / 2.5))
    #Use distance from mean of upper/lower error bounds to calculate new sigma values
    positiveErrorDistance = abs(convertedPositiveErrors - fluxList)
    negativeErrorDistance = abs(convertedNegativeErrors - fluxList)
    meanErrorList = (positiveErrorDistance * negativeErrorDistance) ** (0.5)
    return fluxList, meanErrorList

# --------PLANETARY PARAMETERS UI------------------------------------------
# Get the user's confirmation of values that will later be used in lightcurve fit
def get_planetary_parameters(candplanetbool, userpdict, pdict=None):
    print('\n*******************************************')
    print("Planetary Parameters for Lightcurve Fitting\n")

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
                print("\n\n*** WARNING: %s initialization file's %s does not match the NASA Exoplanet Archive. ***\n" % (pdict['pName'], planet_params[idx]))
                print("\tNASA Exoplanet Archive value: %s" % pdict[item])
                print("\tInitialization file value: %s" % userpdict[item])
                print("\nWould you like to: (1) use NASA Exoplanet Archive value, (2) use initialization file value, or (3) enter in a new value.")
                option = user_input('Which option do you choose? (1/2/3): ', type_=int, val1=1, val2=2, val3=3)
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
        print('\n\n*** Here are the values scraped from the NASA Exoplanet Archive for %s that were not set (or set to null) in your initialization file. ***' % pdict['pName'])
        print('For each planetary parameter, enter "y" if you agree and "n" if you disagree.')
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
                print("\n\n*** WARNING: %s initialization file's %s does not match the NASA Exoplanet Archive. ***\n" % (pdict['pName'], planet_params[i]))
                print("\tNASA Exoplanet Archive value: %s" % pdict[key])
                print("\tInitialization file value: %s" % userpdict[key])
                print("\nWould you like to: (1) use NASA Exoplanet Archive value, (2) use initialization file value, or (3) enter in a new value.")
                option = user_input('Which option do you choose? (1/2/3): ', type_=int, val1=1, val2=2, val3=3)
                if option == 1:
                    userpdict[key] = pdict[key]
                elif option == 2:
                    continue
                else:
                    userpdict[key] = user_input('Enter the ' + planet_params[i] + ': ', type_=type(userpdict[key]))
            # Did not use initialization file or null
            else:
                print('\n' + pdict['pName'] + ' ' + planet_params[i] + ': ' + str(pdict[key]))
                agreement = user_input('Do you agree? (y/n): ', type_=str, val1='y', val2='n')
                if agreement == 'y':
                    userpdict[key] = pdict[key]
                else:
                    userpdict[key] = user_input('Enter the ' + planet_params[i] + ': ', type_=type(pdict[key]))

    # Exoplanet not confirmed in NASA Exoplanet Archive
    else:
        for i, key in enumerate(userpdict):
            if key in ('ra', 'dec'):
                continue
            # Used initialization file and is not empty
            if userpdict[key] is not None:
                agreement = user_input('%s: %s \nDo you agree? (y/n): '
                                       % (planet_params[i], userpdict[key]), type_=str, val1='y', val2='n')
                if agreement == 'y':
                    continue
                else:
                    userpdict[key] = user_input('Enter the ' + planet_params[i] + ': ', type_=type(userpdict[key]))
            # Did not use initialization file
            else:
                if key in ('pName', 'sName'):
                    userpdict[key] = input('\nEnter the ' + planet_params[i] + ': ')
                else:
                    userpdict[key] = user_input('Enter the ' + planet_params[i] + ': ', type_=float)
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
                        if file.lower().endswith(exti.lower()) and file[0:2] not in ('ref', 'new'):
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
            if extaddoption == 'y':
                file_extensions.append(input('Please enter the extension you want to add (EX: .FITS): '))
            else:
                directory = input("Enter the directory path where " + filename + " files are located: ")
        except OSError:
            print("\nError: No such directory exists when searching for FITS files. Please try again.")
            directory = input("Enter the directory path where " + filename + " files are located: ")


# Calculating Limb Darkening Parameters using LDTK
def ld_nonlinear(teff, teffpos, teffneg, met, metpos, metneg, logg, loggpos, loggneg):
                     # Source for FWHM band wavelengths (units: nm): https://www.aavso.org/filters
                     # Near-Infrared
    minmaxwavelen = {('J NIR 1.2micron', 'J'): (1040.00, 1360.00), ('H NIR 1.6micron', 'H'): (1420.00, 1780.00),
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

                     # LCO, Source: Kalee Tock & Michael Fitzgerald, https://lco.global/observatory/instruments/filters/
                     ('LCO Bessell B', 'N/A'): (391.60, 480.60), ('LCO Bessell V', 'N/A'): (502.80, 586.80),
                     ('LCO Pan-STARRS w', 'N/A'): (404.20, 845.80), ("LCO SDSS u'", 'N/A'): (325.50, 382.50),
                     ("LCO SDSS g'", 'N/A'): (402.00, 552.00), ("LCO SDSS r'", 'N/A'): (552.00, 691.00),
                     ("LCO SDSS i'", 'N/A'): (690.00, 819.00)}

    print('\n***************************')
    print('Limb Darkening Coefficients')
    print('***************************')

    print('\nStandard bands available to filter for limb darkening parameters (https://www.aavso.org/filters)'
          '\nas well as filters for MObs and LCO (0.4m telescope) datasets:\n')
    for key, value in minmaxwavelen.items():
        print('\t{}: {} - ({:.2f}-{:.2f}) nm'.format(key[1], key[0], value[0], value[1]))

    ldopt = user_input('\nWould you like EXOTIC to calculate your limb darkening parameters with uncertainties? (y/n): ',
                       type_=str, val1='y', val2='n')

    # User decides to allow EXOTIC to calculate limb darkening parameters
    if ldopt == 'y':
        standcustomopt = user_input('Please enter 1 to use a standard filter or 2 for a customized filter: ',
                                    type_=int, val1=1, val2=2)

        # Standard filters calculating limb darkening parameters
        if standcustomopt == 1:
            while True:
                try:
                    filtername = input('\nPlease enter in the filter type (EX: Johnson V, V, STB, RJ): ')
                    for key, value in minmaxwavelen.items():
                        if filtername in (key[0], key[1]) and filtername != 'N/A':
                            filtername = (key[0], key[1])
                            break
                    else:
                        raise KeyError
                    break
                except KeyError:
                    print('Error: The entered filter is not in the provided list of standard filters.')

            wlmin = [minmaxwavelen[filtername][0] / 1000]
            wlmax = [minmaxwavelen[filtername][1] / 1000]
            filtername = filtername[1]

        # Custom filters calculating limb darkening parameters
        else:
            wlmin = [float(input('FWHM Minimum wavelength (nm): ')) / 1000]
            wlmax = [float(input('FWHM Maximum wavelength (nm): ')) / 1000]
            filtername = 'N/A'


        priors = {'T*': teff, 'T*_uperr': teffpos, 'T*_lowerr': teffneg,
                  'FEH*': met, 'FEH*_uperr': metpos, 'FEH*_lowerr': metneg,
                  'LOGG*': logg, 'LOGG*_uperr': loggpos, 'LOGG*_lowerr': loggneg}

        ldparams = createldgrid(np.array(wlmin), np.array(wlmax), priors)

        nlld0 = ldparams['LD'][0][0], ldparams['ERR'][0][0]
        nlld1 = ldparams['LD'][1][0], ldparams['ERR'][1][0]
        nlld2 = ldparams['LD'][2][0], ldparams['ERR'][2][0]
        nlld3 = ldparams['LD'][3][0], ldparams['ERR'][3][0]

    # User enters in their own limb darkening parameters with uncertainties
    else:
        filtername = input('\nEnter in your filter name: ')
        nlld0 = user_input('\nEnter in your first nonlinear term: ', type_=float)
        nlld0unc = user_input('Enter in your first nonlinear term uncertainty: ', type_=float)
        nlld1 = user_input('\nEnter in your second nonlinear term: ', type_=float)
        nlld1unc = user_input('Enter in your second nonlinear term uncertainty: ', type_=float)
        nlld2 = user_input('\nEnter in your third nonlinear term: ', type_=float)
        nlld2unc = user_input('Enter in your third nonlinear term uncertainty: ', type_=float)
        nlld3 = user_input('\nEenter in your fourth nonlinear term: ', type_=float)
        nlld3unc = user_input('Enter in your fourth nonlinear term uncertainty: ', type_=float)
        nlld0, nlld1, nlld2, nlld3 = (nlld0, nlld0unc), (nlld1, nlld1unc), (nlld2, nlld2unc), (nlld3, nlld3unc)

    return nlld0, nlld1, nlld2, nlld3, filtername


# Checks for corrupted FITS files
def check_file_corruption(files):

    with warnings.catch_warnings():
        warnings.simplefilter('error', category=AstropyWarning)
        for file in files:
            try:
                with fits.open(file, checksum=True, ignore_missing_end=True) as hdul:
                    pass
            except (AstropyWarning, OSError) as e:
                print('Found corrupted file and removing from reduction: {}, ({})'.format(file,e))
                files.remove(file)
        return files


# Check for WCS in the user's imaging data and possibly plate solves.
def check_wcs(fits_file, saveDirectory):
    hdulist = fits.open(name=fits_file, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
    header = hdulist[0].header
    hdulist.close()  # close stream
    del hdulist

    # MJD seems sometimes throw off an error. Deleted since not important for plate solving
    try:
        del header['MJD-OBS']
    except:
        pass

    # Gets the WCS of the header and checks to see if it exists
    wcsheader = WCS(header)
    wcsExists = wcsheader.is_celestial

    # If the fits file has WCS info, ask the user if they trust it
    if wcsExists:
        trustWCS = user_input('The imaging data from your file has WCS information. Do you trust this? (y/n): ', type_=str, val1='y', val2='n')
    else:
        trustWCS = 'n'

    if trustWCS == 'n':
        plateSol = user_input("\nWould you like to upload the your image for a plate solution?"
                              "\nDISCLAIMER: One of your imaging files will be publicly viewable on nova.astrometry.net. (y/n): ", type_=str, val1='y', val2='n')
        # Plate solve the fits file
        if plateSol == 'y':
            print("\nGetting the plate solution for your imaging file. Please wait.")
            global done
            done = False
            t = threading.Thread(target=animate, daemon=True)
            t.start()

            # Plate solves the first imaging file
            imagingFile = fits_file
            wcsFile = plate_solution(imagingFile, saveDirectory)
            done = True

            # Return plate solution from nova.astrometry.net
            return wcsFile
        else:
            # User either did not want a plate solution or had one in file header and decided not to use it,
            # therefore return nothing
            return False
    else:
        # User trusted their imaging file header's WCS
        return fits_file


# Gets the WCS of a .fits file for the user from nova.astrometry.net w/ API key
def plate_solution(fits_file, saveDirectory):
    default_url = 'http://nova.astrometry.net/api/'

    # Login to Exoplanet Watch's profile w/ API key. If session fails, allow 5 attempts of
    # rejoining before returning False and informing user of technical failure.
    for i in range(5):
        try:
            r = requests.post(default_url + 'login', data={'request-json': json.dumps({"apikey": "vfsyxlmdxfryhprq"})})
            sess = r.json()['session']
            break
        except Exception:
            if i == 4:
                print('Imaging file could not receive a plate solution due to technical difficulties '
                      'from nova.astrometry.net. Please try again later. Data reduction will continue.')
                return False
            time.sleep(5)
            continue

    # Saves session number to upload imaging file
    files = {'file': open(fits_file, 'rb')}
    headers = {'request-json': json.dumps({"session": sess}), 'allow_commercial_use': 'n',
               'allow_modifications': 'n', 'publicly_visible': 'n'}

    # Uploads the .fits file to nova.astrometry.net
    r = requests.post(default_url + 'upload', files=files, data=headers)
    if r.json()['status'] != 'success':
        print('Imaging file could not receive a plate solution due to technical difficulties '
              'from nova.astrometry.net. Please try again later. Data reduction will continue.')
        return False

    # Saves submission id for checking on the status of image uploaded
    sub_id = r.json()['subid']
    submissions_url = 'http://nova.astrometry.net/api/submissions/%s' % sub_id

    # Once the image has successfully been plate solved, the following loop will break
    while True:
        r = requests.get(submissions_url)
        if r.json()['job_calibrations']:
            break
        time.sleep(5)

    # Checks the job id's status for parameters
    job_id = r.json()['jobs']
    job_url = 'http://nova.astrometry.net/api/jobs/%s' % job_id[0]
    wcs_file = os.path.join(saveDirectory, 'newfits.fits')

    # Checks the job id's status
    while True:
        r = requests.get(job_url)
        if r.json()['status'] == 'success':
            fits_download_url = 'http://nova.astrometry.net/new_fits_file/%s' % job_id[0]
            r = requests.get(fits_download_url)
            with open(wcs_file, 'wb') as f:
                f.write(r.content)
            print('\n\nSuccess. ')
            return wcs_file
        elif r.json()['status'] == 'failure':
            print('\n\n.FITS file has failed to be given WCS.')
            return False
        time.sleep(5)


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
                raise ValueError
            if expdec - uncert >= declist[pixy][pixx] or declist[pixy][pixx] >= expdec + uncert:
                print('\nError: The Y Pixel Coordinate entered does not match the declination.')
                raise ValueError
            return pixx, pixy
        except ValueError:
            repixopt = user_input('Would you like to re-enter the pixel coordinates? (y/n): ', type_=str, val1='y', val2='n')

            # User wants to change their coordinates
            if repixopt == 'y':
                # Checks for the closest pixel location in ralist and declist for expected ra and dec
                dist = (ralist - expra) ** 2 + (declist - expdec) ** 2
                pixy, pixx = np.unravel_index(dist.argmin(), dist.shape)
                searchopt = user_input('Here are the suggested pixel coordinates: X Pixel: %s Y Pixel: %s'
                                       '\nWould you like to use these? (y/n): ' % (pixx, pixy), type_=str, val1='y', val2='n')
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
def variableStarCheck(refx, refy, hdulWCS):
    #Read in WCS data from plate solution file and convert comparison star coordinates from pixel to WCS
    w = WCS(hdulWCS[0].header)
    world = w.wcs_pix2world(np.array([[refx, refy]], dtype = np.float64), 1)
    ra = world[0][0]
    dec = world[0][1]
    sample = SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')

    #Query GAIA first to check for variability using the phot_variable_flag trait
    radius = u.Quantity(20.0, u.arcsec)
    gaiaQuery = Gaia.cone_search_async(sample, radius)
    gaiaResult = gaiaQuery.get_results()

    #Individually go through the phot_variable_flag indicator for each star to see if variable or not
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

    #query SIMBAD and search identifier result table to determine if comparison star is variable in any form
    #This is a secondary check if GAIA query returns inconclusive results
    simbad_result = Simbad.query_region(sample, radius=20*u.arcsec)
    try:
        starName = simbad_result['MAIN_ID'][0].decode("utf-8")
    except:
        print("Your star cannot be resolved in SIMBAD. Proceed with caution.")
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
        time.sleep(5)

    return imagedata, boollist


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


def gaussian_psf(x,y,x0,y0,a,sigx,sigy,rot, b):
    rx = (x-x0)*np.cos(rot) - (y-y0)*np.sin(rot)
    ry = (x-x0)*np.sin(rot) + (y-y0)*np.cos(rot)
    gausx = np.exp(-(rx)**2 / (2*sigx**2) )
    gausy = np.exp(-(ry)**2 / (2*sigy**2) )
    return a*gausx*gausy + b


def fit_psf(data,pos,init,lo,up,psf_function=gaussian_psf,lossfn='linear',method='trf',box=15):
    xv,yv = mesh_box(pos, box)

    def fcn2min(pars):
        model = psf_function(xv,yv,*pars)
        return (data[yv,xv]-model).flatten()

    if method == 'trf':
        res = least_squares(fcn2min,x0=[*pos,*init],bounds=[lo,up],loss=lossfn,jac='3-point',method=method)
    else:
        res = least_squares(fcn2min,x0=[*pos,*init],loss=lossfn,jac='3-point',method=method)
    return res.x


def mesh_box(pos,box):
    pos = [int(np.round(pos[0])), int(np.round(pos[1]))]
    x = np.arange(pos[0]-box, pos[0]+box+1)
    y = np.arange(pos[1]-box, pos[1]+box+1)
    xv, yv = np.meshgrid(x, y)
    return xv.astype(int), yv.astype(int)


# Method fits a 2D gaussian function that matches the star_psf to the star image and returns its pixel coordinates
def fit_centroid(data, pos, init=None, box=10):

    # get sub field in image
    xv, yv = mesh_box(pos, box)
    wx, wy = pos # could take flux weighted centroid if not crowded

    if init:
        pass
    else:
        init = [np.nanmax(data[yv,xv])-np.nanmin(data[yv,xv]), 1, 1, 0, np.nanmin(data[yv,xv])]

    try:
        # fit gaussian PSF
        pars = fit_psf(
            data,
            [wx, wy], # position estimate
            init,    # initial guess: [amp, sigx, sigy, rotation, bg]
            [wx-5, wy-5, 0, 0, 0, -np.pi/4, np.nanmin(data)-1 ], # lower bound: [xc, yc, amp, sigx, sigy, rotation,  bg]
            [wx+5, wy+5, 1e7, 20, 20, np.pi/4, np.nanmax(data[yv,xv])+1 ], # upper bound
            psf_function=gaussian_psf, method='trf',
            box=box  # only fit a subregion +/- 5 px from centroid
        )
    except:
        print("WARNING trouble fitting Gaussian PSF to star at {},{}".format(wx,wy))
        print("  check location of comparison star in the first few images")
        print("  fitting parameters are out of bounds")
        print("  init:",init)
        print(" lower:",[wx-5, wy-5, 0, 0, 0, -np.pi/4, np.nanmin(data)-1 ] )
        print(" upper:",[wx+5, wy+5, 1e7, 20, 20, np.pi/4, np.nanmax(data[yv,xv])+1 ])

        # use LM in unbounded optimization
        pars = fit_psf(
            data, [wx, wy], init,
            [wx-5, wy-5, 0, 0, 0, -np.pi/4, np.nanmin(data)-1 ],
            [wx+5, wy+5, 1e7, 20, 20, np.pi/4, np.nanmax(data[yv,xv])+1 ],
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
    bdata[bdata < 0] = 0

    apertures = CircularAperture(positions, r=r)
    phot_table = aperture_photometry(bdata, apertures, method='exact')

    return float(phot_table['aperture_sum']), bgflux

def skybg_phot(data, xc, yc, r=10, dr=5, ptol=75):
    # create a crude annulus to mask out bright background pixels
    xv, yv = mesh_box([xc, yc], np.round(r+dr))
    rv = ((xv-xc)**2 + (yv-yc)**2)**0.5
    mask = (rv > r) & (rv < (r+dr))
    cutoff = np.percentile(data[yv, xv][mask], ptol)
    dat = np.array(data[yv,xv],dtype=float)
    dat[dat > cutoff] = np.nan # ignore pixels brighter than percntile 
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
    plt.figure()
    plt.plot(times - np.nanmin(times), xTarg, '-bo')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('X Pixel Position')
    plt.title(targetname + ' X Centroid Position ' + date)
    plt.savefig(infoDict['saveplot'] + 'temp/XCentroidPosition' + targetname + date + '.png')
    plt.close()

    # Y TARGET
    plt.figure()
    plt.plot(times - np.nanmin(times), yTarg, '-bo')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('Y Pixel Position')
    plt.title(targetname + ' Y Centroid Position ' + date)
    plt.savefig(infoDict['saveplot'] + 'temp/YCentroidPos' + targetname + date + '.png')
    plt.close()

    # X COMP
    plt.figure()
    plt.plot(times - np.nanmin(times), xRef, '-ro')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('X Pixel Position')
    plt.title('Comp Star X Centroid Position ' + date)
    plt.savefig(infoDict['saveplot'] + 'temp/CompStarXCentroidPos' + date + '.png')
    plt.close()

    # Y COMP
    plt.figure()
    plt.plot(times - np.nanmin(times), yRef, '-ro')
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('Y Pixel Position')
    plt.title('Comp Star Y Centroid Position ' + date)
    plt.savefig(infoDict['saveplot'] + 'temp/CompStarYCentroidPos' + date + '.png')
    plt.close()

    # X DISTANCE BETWEEN TARGET AND COMP
    plt.figure()
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('X Pixel Distance')
    for e in range(0, len(xTarg)):
        plt.plot(times[e] - np.nanmin(times), abs(int(xTarg[e]) - int(xRef[e])), 'bo')
    plt.title('Distance between Target and Comparison X position')
    plt.savefig(infoDict['saveplot'] + 'temp/XCentroidDistance' + targetname + date + '.png')
    plt.close()

    # Y DISTANCE BETWEEN TARGET AND COMP
    plt.figure()
    plt.xlabel('Time (JD-' + str(np.nanmin(times)) + ')')
    plt.ylabel('Y Pixel Difference')
    d = 0
    for d in range(0, len(yTarg)):
        plt.plot(times[d] - np.nanmin(times), abs(int(yTarg[d]) - int(yRef[d])), 'bo')
    plt.title('Difference between Target and Comparison Y position')
    plt.savefig(infoDict['saveplot'] + 'temp/YCentroidDistance' + targetname + date + '.png')
    plt.close()


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
    targx, targy, targamplitude, targsigX, targsigY, targrot, targoff = fit_centroid(firstImageData, [UIprevTPX, UIprevTPY], box=10)
    refx, refy, refamplitude, refsigX, refsigY, refrot, refoff = fit_centroid(firstImageData, [UIprevRPX, UIprevRPY], box=10)

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
            prevTPX, prevTPY, prevRPX, prevRPY = UIprevTPX, UIprevTPY, UIprevRPX, UIprevRPY
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
        tx, ty, tamplitude, tsigX, tsigY, trot, toff = fit_centroid(imageData, [prevTPX, prevTPY], init=myPriors, box=10)
        tpsfFlux = 2*np.pi*tamplitude*tsigX*tsigY
        currTPX = tx
        currTPY = ty

        # Fits Centroid for Reference
        rGuessAmp = refSearchA.max() - refSearchA.min()
        myRefPriors = [rGuessAmp, prevRSigX, prevRSigY, 0, refSearchA.min()]
        rx, ry, ramplitude, rsigX, rsigY, rrot, roff = fit_centroid(imageData, [prevRPX, prevRPY], init=myRefPriors, box=10)
        rpsfFlux = 2*np.pi*ramplitude*rsigX*rsigY
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

    ax1.clear()
    ax1.set_title(target_name)
    ax1.set_ylabel('Normalized Flux')
    ax1.set_xlabel('Time (jd)')
    ax1.plot(timesListed, normalizedFluxVals, 'bo')


def parse_args():
    # TODO
    parser = argparse.ArgumentParser()

    help_ = "Choose a target to process"
    parser.add_argument("-t", "--target", help=help_, type=str, default="all")

    return parser.parse_args()


def main():
    # TODO use text based interface if no command line arguments

    print('\n')
    print('*************************************************************')
    print('Welcome to the EXOplanet Transit Interpretation Code (EXOTIC)')
    print("Version ", __version__)
    print('*************************************************************\n')

    # ---INITIALIZATION-------------------------------------------------------
    global infoDict
    global UIprevTPX, UIprevTPY, UIprevRPX, UIprevRPY
    global distFC
    global ax1

    targetFluxVals, referenceFluxVals, normalizedFluxVals, targUncertanties, refUncertanties, timeList, phasesList, airMassList = (
        [] for l in range(8))

    fileNameList, timeSortedNames, xTargCent, yTargCent, xRefCent, yRefCent, finXTargCent, finYTargCent, finXRefCent, finYRefCent = (
        [] for m in range(10))

    timesListed = []  # sorted times of observation
    fileNumber = 1  # initializes file number to one
    minSTD = 100000  # sets the initial minimum standard deviation absurdly high so it can be replaced immediately
    minChi2 = 100000
    distFC = 10  # gaussian search area
    context = {}

    # ---USER INPUTS--------------------------------------------------------------------------

    realTimeAns = user_input('Enter "1" for Real Time Reduction or "2" for for Complete Reduction: ', type_=int, val1=1, val2=2)

    #############################
    # Real Time Reduction Routine
    #############################

    if realTimeAns == 1:
        print('\n**************************************************************')
        print('Real Time Reduction ("Control + C"  or close the plot to quit)')
        print('**************************************************************\n')

        directToWatch = str(input("Enter the Directory Path of imaging files: "))
        directoryP = ""
        directoryP = directToWatch
        directToWatch, inputfiles = check_imaging_files(directToWatch, 'imaging')

        targetName = str(input("Enter the Planet Name: "))

        while True:
            try:
                carryOn = input('Type continue after the first image has been taken and saved: ')
                if carryOn != 'continue':
                    raise ValueError
                break
            except ValueError:
                continue

        UIprevTPX = user_input(targetName + " X Pixel Coordinate: ", type_=int)
        UIprevTPY = user_input(targetName + " Y Pixel Coordinate: ", type_=int)
        UIprevRPX = user_input("Comp Star X Pixel Coordinate: ", type_=int)
        UIprevRPY = user_input("Comp Star Y Pixel Coordinate: ", type_=int)

        print('Real Time Plotting ("Control + C" or close the plot to quit)')
        print('\nPlease be patient. It will take at least 15 seconds for the first image to get plotted.')

        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.set_title(targetName)
        ax1.set_ylabel('Normalized Flux')
        ax1.set_xlabel('Time (jd)')

        anim = FuncAnimation(fig, realTimeReduce, fargs=(targetName), interval=15000)  # refresh every 15 seconds
        plt.show()

    ###########################
    # Complete Reduction Routine
    ###########################

    # ----USER INPUTS----------------------------------------------------------
    else:
        print('\n**************************')
        print('Complete Reduction Routine')
        print('**************************\n')

        directoryP = ""
        compStarList = []

        infoDict = {'fitsdir': None, 'saveplot': None, 'flatsdir': None, 'darksdir': None, 'biasesdir': None,
                    'aavsonum': None, 'secondobs': None, 'date': None, 'lat': None, 'long': None,'elev': None,
                    'ctype': None, 'pixelbin': None, 'filter': None, 'notes': None, 'tarcoords': None, 'compstars': None}

        userpDict = {'ra': None, 'dec': None, 'pName': None, 'sName': None, 'pPer': None, 'pPerUnc': None,
                     'midT': None, 'midTUnc': None, 'rprs': None, 'rprsUnc': None, 'aRs': None, 'aRsUnc': None,
                     'inc': None, 'incUnc': None, 'ecc': None, 'teff': None,
                     'teffUncPos': None, 'teffUncNeg': None, 'met': None, 'metUncPos': None, 'metUncNeg': None,
                     'logg': None, 'loggUncPos': None, 'loggUncNeg': None}

        fitsortext = user_input('Enter "1" to perform aperture photometry on fits files or "2" to start '
                                'with pre-reduced data in a .txt format: ', type_=int, val1=1, val2=2)

        fileorcommandline = user_input('How would you like to input your initial parameters? '
                                       'Enter "1" to use the Command Line or "2" to use an input file: ', type_=int, val1=1, val2=2)

        cwd = os.getcwd() # os.path.join(os.path.split(os.getcwd())[0], '')

        # Read in input file rather than using the command line
        if fileorcommandline == 2:
            print("\nYour current working directory is: ", cwd)
            print("\nPotential initialization files I've found in {} are: ".format(cwd))
            [print(i) for i in g.glob(cwd + "/*.json")]

            # Parse input file
            while True:
                try:
                    initfilename = str(input("\nPlease enter the Directory and Filename of your Initialization File: "))
                    if initfilename == 'ok':
                        initfilename = "/Users/rzellem/Documents/EXOTIC/inits.json"
                    infoDict, userpDict = inits_file(initfilename, infoDict, userpDict)
                    break
                except FileNotFoundError:
                    print("Error: Initialization file not found. Please try again.")
                except IsADirectoryError:
                    print('Error: Entered a directory. Please try again.')

            if infoDict['flatsdir'] == "n":
                flatsBool = False
            else:
                flatsBool = True

            if infoDict['darksdir'] == "n":
                darksBool = False
            else:
                darksBool = True

            if infoDict['biasesdir'] == "n":
                biasesBool = False
            else:
                biasesBool = True

            if flatsBool + darksBool + biasesBool:
                cals = 'y'
            else:
                cals = 'n'

            # Initial position of target star
            UIprevTPX = infoDict['tarcoords'][0]
            UIprevTPY = infoDict['tarcoords'][1]

            # Read in locations of comp stars
            for i in infoDict['compstars']:
                compStarList.append((i[0], i[1]))

        if fitsortext == 1:
            # File directory name and initial guess at target and comp star locations on image.
            if fileorcommandline == 1:
                infoDict['fitsdir'] = str(input("\nEnter the Directory of imaging files: "))

            infoDict['fitsdir'], inputfiles = check_imaging_files(infoDict['fitsdir'], 'imaging')
        else:
            datafile = str(input("Enter the path and filename of your data file: "))
            if datafile == 'ok':
                datafile = "/Users/rzellem/Documents/EXOTIC/sample-data/NormalizedFluxHAT-P-32 bDecember 17, 2017.txt"
                # datafile = "/Users/rzellem/Downloads/fluxorama.csv
            try:
                initf = open(datafile, 'r')
            except FileNotFoundError:
                print("Data file not found. Please try again.")
                sys.exit()

            infoDict['exposure'] = user_input("Please enter your image exposure time in seconds: ", type_=int)

            processeddata = initf.readlines()

        if fileorcommandline == 1:
            infoDict['saveplot'] = str(input("Enter the directory to save plots into or type new to create one: "))

        # Check to see if the save directory exists
        while True:
            try:
                if infoDict['saveplot'] == 'new':
                    infoDict['saveplot'] = create_directory(cwd)
                    break
                infoDict['saveplot'] = os.path.join(infoDict['saveplot'], '')
                # In case the user forgets the trailing / for the folder
                if os.path.isdir(infoDict['saveplot']):
                    break
                raise OSError
            except OSError:
                print('Error: the entered save directory does not exist. Please try again.')
                infoDict['saveplot'] = input("Enter the directory to save plots into or type new to create one: ")

        # Make a temp directory of helpful files
        try:
            os.makedirs(infoDict['saveplot'] + "temp/")
        except FileExistsError:
            # directory already exists
            pass

        if fileorcommandline == 1:
            userpDict['pName'] = str(input("\nEnter the Planet Name: "))

        userpDict['pName'], CandidatePlanetBool, pDict = get_planet_name(userpDict['pName'])

        # observation date
        if fileorcommandline == 1:
            infoDict['date'] = str(input("\nEnter the Observation Date: "))

        # Using a / in your date can screw up the file paths- this will check user's date
        while "/" in infoDict['date']:
            print("Do not use / in your date. Please try again.")
            infoDict['date'] = str(input("\nEnter the Observation Date: "))

        if fitsortext == 1:
            if fileorcommandline == 1:
                infoDict['lat'] = input("Enter the latitude of where you observed (deg) "
                                        "(Don't forget the sign where North is '+' and South is '-'): ")
            # Latitude
            while True:
                try:
                    infoDict['lat'] = infoDict['lat'].replace(' ', '')
                    if infoDict['lat'][0] != '+' and infoDict['lat'][0] != '-':
                        raise ValueError("You forgot the sign for the latitude! North is '+' and South is '-'. Please try again.")
                    lati = float(infoDict['lat'])
                    if lati <= -90.00 or lati >= 90.00:
                        raise ValueError('Your latitude is out of range. Please enter a latitude between -90 and +90 (deg)')
                    break
                # check to make sure they have a sign
                except ValueError as err:
                    print(err.args)
                    infoDict['lat'] = input("Enter the latitude of where you observed (deg) "
                                            "(Don't forget the sign where North is '+' and South is '-'): ")

            if fileorcommandline == 1:
                infoDict['long'] = input("Enter the longitude of where you observed (deg) "
                                         "(Don't forget the sign where East is '+' and West is '-'): ")
            # Longitude
            while True:
                try:
                    infoDict['long'] = infoDict['long'].replace(' ', '')
                    if infoDict['long'][0] != '+' and infoDict['long'][0] != '-':
                        raise ValueError("You forgot the sign for the longitude! East is '+' and West is '-'. Please try again.")
                    longit = float(infoDict['long'])
                    if longit <= -180.00 or longit >= 180.00:
                        raise ValueError('Your longitude is out of range. Please enter a longitude between -180 and +180 (deg)')
                    break
                # check to make sure they have a sign
                except ValueError as err:
                    print(err.args)
                    infoDict['long'] = input("Enter the longitude of where you observed (deg) "
                                             "(Don't forget the sign where East is '+' and West is '-'): ")

            if fileorcommandline == 1:
                infoDict['elev'] = user_input("Enter the elevation (in meters) of where you observed: ", type_=float)

            # TARGET STAR
            if fileorcommandline == 1:
                UIprevTPX = user_input('\n' + userpDict['pName'] + " X Pixel Coordinate: ", type_=int)
                UIprevTPY = user_input(userpDict['pName'] + " Y Pixel Coordinate: ", type_=int)
                numCompStars = user_input("How many comparison stars would you like to use? (1-10) ", type_=int)

                for num in range(numCompStars):
                    xpix = user_input("Comparison Star %s X Pixel Coordinate: " % str(num+1), type_=int)
                    ypix = user_input("Comparison Star %s Y Pixel Coordinate: " % str(num+1), type_=int)
                    compStarList.append((xpix, ypix))

            # ---HANDLE CALIBRATION IMAGES------------------------------------------------
            if fileorcommandline == 1:
                cals = user_input('\nDo you have any calibration images (flats, darks or biases)? (y/n): ', type_=str, val1='y', val2='n')

            # if they have cals, handle them by calculating the median flat, dark or bias
            if cals == 'y':

                # flats
                # THIS DOES NOT ACCOUNT FOR CALIBRATING THE FLATS, WHICH COULD BE TAKEN AT A DIFFERENT EXPOSURE TIME
                if fileorcommandline == 1:
                    flats = user_input('\nDo you have flats? (y/n): ', type_=str, val1='y', val2='n')
                    if flats == 'y':
                        flatsBool = True
                        infoDict['flatsdir'] = str(input('Enter the directory path to your flats (must be in their own separate folder): '))  # +"/*.FITS"
                    else:
                        flatsBool = False

                    if flatsBool:
                        infoDict['flatsdir'], inputflats = check_imaging_files(infoDict['flatsdir'], 'flats')
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
                    if darks == 'y':
                        darksBool = True
                        infoDict['darksdir'] = str(input('Enter the directory path to your darks (must be in their own separate folder): '))  # +"/*.FITS"
                    else:
                        darksBool = False

                # Only do the dark correction if user selects this option
                if darksBool:
                    infoDict['darksdir'], inputdarks = check_imaging_files(infoDict['darksdir'], 'darks')
                    darksImgList = []
                    for darkFile in inputdarks:
                        darkData = fits.getdata(darkFile, ext=0)
                        darksImgList.append(darkData)
                    generalDark = np.median(darksImgList, axis=0)

                # biases
                if fileorcommandline == 1:
                    biases = user_input('\nDo you have biases? (y/n): ', type_=str, val1='y', val2='n')
                    if biases == 'y':
                        biasesBool = True
                        infoDict['biasesdir'] = str(input('Enter the directory path to your biases (must be in their own separate folder): '))  # +"/*.FITS"
                    else:
                        biasesBool = False

                if biasesBool:
                    # Add / to end of directory if user does not input it
                    infoDict['biasesdir'], inputbiases = check_imaging_files(infoDict['biasesdir'], 'biases')
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

        # Handle AAVSO Formatting
        if fileorcommandline == 1:
            # userNameEmails = str(input('Please enter your name(s) and email address(es) in the format: Your Name (youremail@example.com), Next Name (nextemail@example.com), etc.  '))
            infoDict['aavsonum'] = str(input('Please enter your AAVSO Observer Account Number (type N/A if you do not currently have an account): '))
            infoDict['secondobs'] = str(input('Please enter your comma-separated secondary observer codes (or type N/A if only 1 observer code): '))
            infoDict['ctype'] = str(input("Please enter your camera type (CCD or DSLR): "))
            infoDict['pixelbin'] = str(input('Please enter your pixel binning: '))
            # infoDict['exposure'] = user_input('Please enter your exposure time (seconds): ', type_=int)
            infoDict['filter'] = str(input('Please enter your filter name from the options at '
                                           'http://astroutils.astronomy.ohio-state.edu/exofast/limbdark.shtml: '))
            infoDict['notes'] = str(input('Please enter any observing notes (seeing, weather, etc.): '))


        pDict = get_planetary_parameters(CandidatePlanetBool, userpDict, pdict=pDict)

        ld0, ld1, ld2, ld3, infoDict['filter'] = ld_nonlinear(pDict['teff'], pDict['teffUncPos'], pDict['teffUncNeg'],
                                                              pDict['met'], pDict['metUncNeg'], pDict['metUncPos'],
                                                              pDict['logg'], pDict['loggUncPos'], pDict['loggUncNeg'])

        # If fits files are used, check that they are not corrupted
        if fitsortext == 1: 
            inputfiles = check_file_corruption(inputfiles)

        exptimes = list()

        if fitsortext == 1:
            print('\n**************************')
            print('Starting Reduction Process')
            print('**************************')

            #########################################
            # FLUX DATA EXTRACTION AND MANIPULATION
            #########################################

            # Loop placed to check user-entered x and y target coordinates against WCS.
            while True:
                fileNumber = 1
                allImageData, timeList, fileNameList, timesListed, airMassList, fileNameStr = [], [], [], [], [], []

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

                    hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
                    imageheader = hdul[0].header
                    # TIME
                    timeVal = getJulianTime(hdul)  # gets the julian time registered in the fits header
                    timeList.append(timeVal)  # adds to time value list
                    fileNameList.append(fileName)

                    # TIME
                    currTime = getJulianTime(hdul)
                    timesListed.append(currTime)

                    # AIRMASS
                    airMass = getAirMass(hdul, pDict['ra'], pDict['dec'], lati, longit, infoDict['elev'])  # gets the airmass at the time the image was taken
                    airMassList.append(airMass)  # adds that airmass value to the list of airmasses

                    # IMAGES
                    hdul[0].scale('float32')
                    allImageData.append(hdul[0].data)

                    # EXPOSURE_TIME
                    exp = imageheader.get('EXPTIME')  #checking for variation in .fits header format
                    if exp:
                        exptimes.append(imageheader['EXPTIME'])
                    else:
                        exptimes.append(imageheader['EXPOSURE'])

                    hdul.close()  # closes the file to avoid using up all of computer's resources
                    del hdul

                consistent_et = False
                if len(exptimes) > 0 :
                    consistent_et = all(elem == exptimes[0] for elem in exptimes)

                exptimes = np.array(exptimes)

                if consistent_et :
                    #print("All Elements in List are Equal")
                    infoDict['exposure'] = exptimes[0]
                else:
                    #print("All Elements in List are Not Equal")
                    infoDict['exposure'] = np.median(exptimes)
                    #print(infoDict['exposure'])

                # Recast list as numpy arrays
                allImageData = np.array(allImageData)
                timesListed = np.array(timesListed)
                airMassList = np.array(airMassList)

                # TODO: Is this dead code? The vars pointing and location are undefined.
                # TODO: comment out conditional block?
                # If all of the airmasses == 1, then you need to calculate the airmass for the user
                if set(airMassList) == 1:
                    pointingAltAz = pointing.transform_to(AltAz(obstime=t, location=location))

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

                # The cosmic ray filter isn't really working for now...so let's just turn it off
                cosmicrayfilter_bool = False
                if cosmicrayfilter_bool:
                    print("\nFiltering your data for cosmic rays.")
                    done = False
                    t = threading.Thread(target=animate, daemon=True)
                    t.start()
                    # # -------COSMIC RAY FILTERING-----------------------------------------------------------------------
                    # For now, this is a simple median filter...in the future, should use something more smart later
                    for xi in np.arange(np.shape(sortedallImageData)[-2]):
                        # print("Filtering for cosmic rays in image row: "+str(xi)+"/"+str(np.shape(sortedallImageData)[-2]))
                        for yi in np.arange(np.shape(sortedallImageData)[-1]):
                            # Simple median filter
                            idx = np.abs(sortedallImageData[:, xi, yi]-np.nanmedian(sortedallImageData[:, xi, yi])) > 5*np.nanstd(sortedallImageData[:, xi, yi])
                            sortedallImageData[idx, xi, yi] = np.nanmedian(sortedallImageData[:, xi, yi])
                            # Filter iteratively until no more 5sigma outliers exist - not currently working, so keep commented out for now
                            # while sum(idx) > 0:
                            #     # sortedallImageData[idx,xi,yi] = np.nanmedian(sortedallImageData[:,xi,yi])
                            #     idx = np.abs(sortedallImageData[:,xi,yi]-np.nanmedian(sortedallImageData[:,xi,yi])) >  5*np.nanstd(sortedallImageData[:,xi,yi])
                    done = True

                # if len(sortedTimeList) == 0:
                #     print("Error: .FITS files not found in " + directoryP)
                #     sys.exit()

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
                        targx, targy, targamplitude, targsigX, targsigY, targrot, targoff = fit_centroid(firstImageData, [UIprevTPX, UIprevTPY],
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
                    sortedallImageData = sortedallImageData - generalDark
                elif biasesBool:
                    print("Bias-correcting images.")
                    sortedallImageData = sortedallImageData - generalBias
                else:
                    pass

                if flatsBool:
                    print("Flattening images.")
                    sortedallImageData = sortedallImageData / generalFlat

                # Reference File
                refFile = os.path.join(infoDict['saveplot'], 'ref_file_%s_%s'
                                       % (firstimagecounter, os.path.split(fileNameStr[firstimagecounter])[-1]))

                # Removes existing file of reference file
                try:
                    os.remove(refFile)
                except OSError:
                    pass
                convertToFITS = fits.PrimaryHDU(data=sortedallImageData[0])
                convertToFITS.writeto(refFile)
                print('\nHere is the path to the reference imaging file EXOTIC: \n' + refFile)
                wcsFile = check_wcs(refFile, infoDict['saveplot'])

                # Check pixel coordinates by converting to WCS. If not correct, loop over again
                if wcsFile:
                    print('Here is the path to your plate solution: \n' + wcsFile)
                    hdulWCS = fits.open(name=wcsFile, memmap=False, cache=False, lazy_load_hdus=False, ignore_missing_end=True)
                    rafile, decfile = get_radec(hdulWCS)

                    # Save previously entered x and y pixel coordinates before checking against plate solution
                    saveUIprevTPX, saveUIprevTPY = UIprevTPX, UIprevTPY
                    UIprevTPX, UIprevTPY = check_targetpixelwcs(UIprevTPX, UIprevTPY, pDict['ra'],
                                                                pDict['dec'], rafile, decfile)
                    # If the coordinates were not changed, do not loop over again
                    if UIprevTPX == saveUIprevTPX and UIprevTPY == saveUIprevTPY:
                        break
                else:
                    break

            # Image Alignment
            sortedallImageData, boollist = image_alignment(sortedallImageData)

            timesListed = timesListed[boollist]
            airMassList = airMassList[boollist]

            minAperture = max(1,int(2 * max(targsigX, targsigY)))
            maxAperture = int(5 * max(targsigX, targsigY) + 1)
            minAnnulus = 2
            maxAnnulus = 5
            # exit()
            # fit centroids for first image to determine priors to be used later
            for compCounter in range(0, len(compStarList)):
                print('\n\n***************************************************************')
                print('Determining Optimal Aperture and Annulus Size for Comp Star #' + str(compCounter + 1))
                print('***************************************************************')

                # #just in case comp star drifted off and timeSortedNames had to be altered, reset it for the new comp star
                # timeSortedNames = tsnCopy

                UIprevRPX, UIprevRPY = compStarList[compCounter]
                print('Target X: ' + str(round(targx)) + ' Target Y: ' + str(round(targy)))
                refx, refy, refamplitude, refsigX, refsigY, retrot, refoff = fit_centroid(firstImageData, [UIprevRPX, UIprevRPY],
                                                                                box=10)
                print('Comparison X: ' + str(round(refx)) + ' Comparison Y: ' + str(round(refy)) + '\n')

                #If plate solution was generated, use it to check if the comparison stars selected are variable
                #If yes, skip determining optimal aperture and annulus for that comparison star
                if wcsFile:
                    print("Checking for variability in current comparison star... ")
                    if variableStarCheck(refx, refy, hdulWCS):
                        print("Current comparison star is variable, proceeding to next star.")
                        continue

                # determines the aperture and annulus combinations to iterate through based on the sigmas of the LM fit
                aperture_min = int(3 * np.nanmax([targsigX, targsigY]))
                aperture_max = int(5 * np.nanmax([targsigX, targsigY]))

                # run through apertures based on PSF shape
                if aperture_min <= 1:
                    aperture_sizes = np.arange(1, 10, 2)
                else:
                    aperture_sizes = np.round(np.linspace(aperture_min, aperture_max, 10),2)

                aperture_sizes = np.append(aperture_sizes, -1*aperture_sizes) # no comparison star
                aperture_sizes = np.append(aperture_sizes, 0) # PSF fit
 
                # single annulus size
                annulus_sizes = [10,12,15]

                target_fits = {}
                ref_fits = {}
                reg_trans = {}

                for apertureR in aperture_sizes:  # aperture loop
                    for annulusR in annulus_sizes:  # annulus loop # no need
                        # don't reprocess 
                        if apertureR < 0 and compCounter > 0:
                            continue 
                        
                        # only do PSF fit in first annulus for loop
                        if apertureR == 0 and annulusR != annulus_sizes[0]:
                            continue

                        if apertureR == 0:
                            print('Testing Comparison Star #' + str(compCounter+1) + ' with a PSF photometry.')
                        elif apertureR < 0 and compCounter == 0:
                            print('Testing NO Comparison Star with a '+str(abs(apertureR))+' pixel aperture and a '+str(abs(annulusR))+' pixel annulus.')
                        else:
                            print('Testing Comparison Star #' + str(compCounter+1) + ' with a '+str(apertureR)+' pixel aperture and a '+str(annulusR)+' pixel annulus.')

                        for fileNumber, imageData in enumerate(sortedallImageData):
                            if apertureR == 0: # psf fit
                                continue

                            # Find the target star in the image and get its pixel coordinates if it is the first file
                            if fileNumber == 0:
                                # Initializing the star location guess as the user inputted pixel coordinates
                                prevTPX, prevTPY, prevRPX, prevRPY = UIprevTPX, UIprevTPY, UIprevRPX, UIprevRPY  # 398, 275, 419, 203
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
                            txmin = int(prevTPX) - distFC  # left
                            txmax = int(prevTPX) + distFC  # right
                            tymin = int(prevTPY) - distFC  # top
                            tymax = int(prevTPY) + distFC  # bottom

                            # boolean that represents if either the target or comp star gets too close to the detector
                            driftBool = False

                            #check if your target star is too close to the edge of the detector
                            if txmin <= 0 or tymin <= 0 or txmax >= len(imageData) or tymax >= len(imageData[0]):
                                print('*************************************************************************************')
                                print('WARNING: In image '+str(fileNumber)+', your target star has drifted too close to the edge of the detector.')
                                #tooClose = int(input('Enter "1" to pick a new comparison star or enter "2" to continue using the same comp star, with the images with all the remaining images ignored \n'))
                                print('All the remaining images after image #'+str(fileNumber-1)+' will be ignored')
                                driftBool = True

                                # mask off the rest of timeSortedNames and then ignore the rest of the procedure until

                            # Set reference search area
                            rxmin = int(prevRPX) - distFC  # left
                            rxmax = int(prevRPX) + distFC  # right
                            rymin = int(prevRPY) - distFC  # top
                            rymax = int(prevRPY) + distFC  # bottom

                            # check if the reference is too close to the edge of the detector
                            if (rxmin <= 0 or rymin <= 0 or rxmax >= len(imageData[0]) or rymax >= len(imageData)):
                                print('*************************************************************************************')
                                print('WARNING: In image '+str(fileNumber)+', your reference star has drifted too close to the edge of the detector.')
                                #tooClose = int(input('Enter "1" to pick a new comparison star or enter "2" to continue using the same comp star, with the images with all the remaining images ignored \n'))
                                print('All the remaining images after image #'+str(fileNumber-1)+' will be ignored for this comparison star')
                                print('*************************************************************************************')
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
                                    print('Error: the Darks have a higher pixel counts than the image itself')
                                myPriors = [tGuessAmp, prevTSigX, prevTSigY, 0, tGuessBkg]

                                # tx, ty, tamplitude, tsigX, tsigY, toff = fit_centroid(imageData, targPos,
                                #                                                     init=myPriors, box=distFC)
                                if fileNumber in target_fits.keys():
                                    tx, ty, tamplitude, tsigX, tsigY, toff = target_fits[fileNumber]
                                else:
                                    tx, ty, tamplitude, tsigX, tsigY, trot, toff = fit_centroid(imageData, targPos,
                                                                                          init=myPriors, box=distFC)
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
                                                                                          init=myRefPriors, box=distFC)
                                    ref_fits[fileNumber] = [rx, ry, ramplitude, rsigX, rsigY, roff]
                                currRPX = rx
                                currRPY = ry

                                # append to list of reference centroid positions for later plotting
                                xRefCent.append(currRPX)
                                yRefCent.append(currRPY)

                                if tamplitude < 0 or tsigX < 0 or tsigY < 0:  # gets rid of negative amplitude values that indicate it couldn't fit gaussian
                                    print('Could not fit 2D gaussian to Target for File Number' + str(fileNumber))

                                elif ramplitude < 0 or rsigX < 0 or rsigY < 0:  # gets rid of negative amplitude values that indicate it couldn't fit gaussian
                                    print('Could not fit 2D gaussian to Comparison Star for File Number' + str(fileNumber))

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
                        
                        if apertureR == 0: # psf fit
                            tpsfflux = []
                            rpsfflux = []
                            for k in target_fits.keys():
                                xc,yc,amp,sigx,sigy,off = target_fits[k]
                                tpsfflux.append(2*np.pi*sigx*sigy*amp)
                                xc,yc,amp,sigx,sigy,off = ref_fits[k]
                                rpsfflux.append(2*np.pi*sigx*sigy*amp)
                            arrayReferences = np.array(rpsfflux)
                            arrayTUnc = arrayFinalFlux**0.5
                            arrayRUnc = arrayReferences**0.5

                            arrayFinalFlux = np.array(tpsfflux) / arrayReferences
                            arrayNormUnc = arrayTUnc / arrayReferences
                        elif apertureR < 0: # no comp star
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
                            'rprs':pDict['rprs'],    # Rp/Rs
                            'ars':pDict['aRs'],      # a/Rs
                            'per':pDict['pPer'],     # Period [day]
                            'inc':pDict['inc'],      # Inclination [deg]
                            'u0': ld0[0], 'u1': ld1[0], 'u2': ld2[0], 'u3': ld3[0],  # limb darkening (nonlinear)
                            'ecc': pDict['ecc'],     # Eccentricity
                            'omega':0,          # Arg of periastron
                            'tmid':pDict['midT'],    # time of mid transit [day]
                            'a1': arrayFinalFlux.mean(), #max() - arrayFinalFlux.min(), #mid Flux
                            'a2': 0,             #Flux lower bound
                        }

                        phase = (arrayTimes[~filtered_data]-prior['tmid'])/prior['per']
                        prior['tmid'] = pDict['midT'] + np.floor(phase).max()*prior['per']
                        upper = pDict['midT']+ 25*pDict['midTUnc'] + np.floor(phase).max()*(pDict['pPer']+25*pDict['pPerUnc'])
                        lower = pDict['midT']- 25*pDict['midTUnc'] + np.floor(phase).max()*(pDict['pPer']-25*pDict['pPerUnc'])

                        if np.floor(phase).max()-np.floor(phase).min() == 0:
                            print("WARNING!")
                            print(" Estimated mid-transit time is not within the observations")
                            print(" Check Period & Mid-transit time in inits.json. Make sure the uncertainties are not 0 or Nan.")
                            print('  obs start:', arrayTimes[~filtered_data].min())
                            print('    obs end:', arrayTimes[~filtered_data].max())
                            print(' tmid prior:', prior['tmid'])

                        # check for Nans + Zeros
                        for k in pDict:
                            if "Unc" in k:
                                if not pDict[k]:
                                    print(" WARNING! {} uncertainty is 0. Please use a non-zero value in inits.json".format(k))
                                    pDict[k] = 1
                                elif pDict[k] == 0 or np.isnan(pDict[k]):
                                    print(" WARNING! {} uncertainty is 0. Please use a non-zero value in inits.json".format(k))
                                    pDict[k] = 1
                            elif not pDict[k]:
                                print(" WARNING! {} is None. Please use a numeric value in inits.json".format(k))
                                pDict[k] = 0

                        mybounds = {
                            'rprs':[0, pDict['rprs']+3*pDict['rprsUnc']],
                            'tmid':[max(lower,arrayTimes[~filtered_data].min()),min(arrayTimes[~filtered_data].max(),upper)],
                            'ars':[pDict['aRs']-5*pDict['aRsUnc'], pDict['aRs']+5*pDict['aRsUnc']],

                            'a1':[0, 3*max(arrayFinalFlux[~filtered_data])],
                            'a2':[-3,3],
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
                            print("  {}: {:.6f}".format(k, myfit.parameters[k])) #, myfit.errors[k]))

                        print('The Residual Standard Deviation is: ' + str(round(100*myfit.residuals.std()/np.median(myfit.data), 6))+"%")
                        print('The Mean Squared Error is: ' + str(round( np.sum(myfit.residuals**2), 6)) + '\n')
                        resstd = myfit.residuals.std()/np.median(myfit.data)
                        if minSTD > resstd:  # If the standard deviation is less than the previous min
                            bestCompStar = compCounter + 1
                            minSTD = resstd  # set the minimum standard deviation to that

                            arrayNormUnc = arrayNormUnc * np.sqrt(myfit.chi2 / myfit.data.shape[0])  # scale errorbars by sqrt(rchi2)
                            minAnnulus = annulusR  # then set min aperature and annulus to those values
                            minAperture = apertureR
                            # gets the centroid trace plots to ensure tracking is working
                            finXTargCentArray = np.array(xTargCent)
                            finYTargCentArray = np.array(yTargCent)
                            finXRefCentArray = np.array(xRefCent)
                            finYRefCentArray = np.array(yRefCent)

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

            if minAperture == 0: # psf
                print('\n*********************************************')
                print('Best Comparison Star: #' + str(bestCompStar))
                print('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                print('Optimal Method: PSF photometry')
                print('********************************************\n')
            elif minAperture < 0: # no comp star
                print('\n*********************************************')
                print('Best Comparison Star: None')
                print('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                print('Optimal Aperture: ' + str(abs(minAperture)))
                print('Optimal Annulus: ' + str(minAnnulus))
                print('********************************************\n')
            else:
                print('\n*********************************************')
                print('Best Comparison Star: #' + str(bestCompStar))
                print('Minimum Residual Scatter: ' + str(round(minSTD * 100, 4)) + '%')
                print('Optimal Aperture: ' + str(minAperture))
                print('Optimal Annulus: ' + str(minAnnulus))
                print('********************************************\n')

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
            if wcsFile:
                if hdulWCS[0].header['COMMENT'][135].split(' ')[0] == 'scale:':
                    imscalen = float(hdulWCS[0].header['COMMENT'][135].split(' ')[1])
                    imscaleunits = 'Image scale in arc-secs/pixel'
                    imscale = imscaleunits + ": " + str(round(imscalen, 2))
                else:
                    i = 100
                    while hdulWCS[0].header['COMMENT'][i].split(' ')[0] != 'scale:':
                        i += 1
                    imscalen = float(hdulWCS[0].header['COMMENT'][i].split(' ')[1])
                    imscaleunits = 'Image scale in arc-secs/pixel'
                    imscale = imscaleunits + ": " + str(round(imscalen, 2))
                hdulWCS.close()  # close stream
                del hdulWCS
            elif "IM_SCALE" in imageheader:
                imscalen = imageheader['IM_SCALE']
                imscaleunits = imageheader.comments['IM_SCALE']
                imscale = imscaleunits + ": " + str(imscalen)
            elif "PIXSCALE" in imageheader:
                imscalen = imageheader['PIXSCALE']
                imscaleunits = imageheader.comments['PIXSCALE']
                imscale = imscaleunits + ": " + str(imscalen)
            else:
                print("Cannot find the pixel scale in the image header.")
                # pixscale = input("Do you know the size of your pixels? (y/n) ")
                # if pixscale == 'y' or pixscale == 'Y' or pixscale == 'yes':
                imscalen = input("Please enter the size of your pixel (e.g., 5 arcsec/pixel). ")
                imscale = "Image scale: " + imscalen
            imwidth = np.shape(sortedallImageData[0])[1]
            imheight = np.shape(sortedallImageData[0])[0]
            picframe = 10*(minAperture+minAnnulus)
            pltx = [min([finXTargCent[0], finXRefCent[0]])-picframe, max([finXTargCent[0], finXRefCent[0]])+picframe]
            FORwidth = pltx[1]-pltx[0]
            plty = [min([finYTargCent[0], finYRefCent[0]])-picframe, max([finYTargCent[0], finYRefCent[0]])+picframe]
            FORheight = plty[1]-plty[0]
            fig, ax = plt.subplots()
            target_circle = plt.Circle((finXTargCent[0], finYTargCent[0]), minAperture, color='lime', fill=False, ls='-', label='Target')
            target_circle_sky = plt.Circle((finXTargCent[0], finYTargCent[0]), minAperture+minAnnulus, color='lime', fill=False, ls='--', lw=.5)
            if minAperture >= 0:
                ref_circle = plt.Circle((finXRefCent[0], finYRefCent[0]), minAperture, color='r', fill=False, ls='-.', label='Comp')
                ref_circle_sky = plt.Circle((finXRefCent[0], finYRefCent[0]), minAperture+minAnnulus, color='r', fill=False, ls='--', lw=.5)
            plt.imshow(np.log10(sortedallImageData[0]), origin='lower', cmap='Greys_r', interpolation=None, vmin=np.log10(np.nanmin(sortedallImageData[0][sortedallImageData[0] > 0])), vmax=np.log10(np.nanmax(sortedallImageData[0][sortedallImageData[0] > 0])))  #,vmax=np.nanmax([arrayTargets[0],arrayReferences[0]]))
            plt.plot(finXTargCent[0], finYTargCent[0], marker='+', color='lime')
            ax.add_artist(target_circle)
            ax.add_artist(target_circle_sky)
            if minAperture >= 0:
                ax.add_artist(ref_circle)
                ax.add_artist(ref_circle_sky)
                plt.plot(finXRefCent[0], finYRefCent[0], '+r')
            plt.xlabel("x-axis [pixel]")
            plt.ylabel("y-axis [pixel]")
            plt.title("FOV for " + pDict['pName'] + "\n(" + imscale + ")")
            plt.xlim(pltx[0], pltx[1])
            plt.ylim(plty[0], plty[1])
            ax.grid(False)
            plt.plot(0, 0, color='lime', ls='-', label='Target')
            if minAperture >= 0:
                plt.plot(0, 0, color='r', ls='-.', label='Comp')
            l = plt.legend(frameon=None, framealpha=0)
            for text in l.get_texts():
                text.set_color("white")
            plt.savefig(infoDict['saveplot'] + "FOV" + pDict['pName'] + infoDict['date'] + ".pdf", bbox_inches='tight')
            plt.close()

            print("\nFOV file saved as: "+infoDict['saveplot'] + "FOV" + pDict['pName'] + infoDict['date'] + ".pdf\n")

            # Take the BJD times from the image headers
            if "BJD_TDB" in imageheader:
                goodTimes = nonBJDTimes
            # If not in there, then convert all the final times into BJD - using astropy alone
            else:
                print("No BJDs in Image Headers. Converting all JDs to BJD_TDBs.")
                print("Please be patient- this step can take a few minutes.")
                # targetloc = astropy.coordinates.SkyCoord(raStr, decStr, unit=(astropy.units.deg,astropy.units.deg), frame='icrs')
                # obsloc = astropy.coordinates.EarthLocation(lat=lati, lon=longit)
                # timesToConvert = astropy.time.Time(nonBJDTimes, format='jd', scale='utc', location=obsloc)
                # ltt_bary = timesToConvert.light_travel_time(targetloc)
                # time_barycentre = timesToConvert.tdb + ltt_bary
                # resultos = time_barycentre.value
                # goodTimes = resultos
                done = False
                t = threading.Thread(target=animate, daemon=True)
                t.start()
                resultos = utc_tdb.JDUTC_to_BJDTDB(nonBJDTimes, ra=pDict['ra'], dec=pDict['dec'], lat=lati, longi=longit, alt=infoDict['elev'])
                goodTimes = resultos[0]
                done = True

            # Centroid position plots
            plotCentroids(finXTargCent, finYTargCent, finXRefCent, finYRefCent, goodTimes, pDict['pName'], infoDict['date'])

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
            plt.figure()
            plt.errorbar(goodTimes, goodTargets, yerr=goodTUnc, linestyle='None', fmt='-o')
            plt.xlabel('Time (BJD)')
            plt.ylabel('Total Flux')
            # plt.rc('grid', linestyle="-", color='black')
            # plt.grid(True)
            plt.title(pDict['pName'] + ' Raw Flux Values ' + infoDict['date'])
            plt.savefig(infoDict['saveplot'] + 'temp/TargetRawFlux' + pDict['pName'] + infoDict['date'] + '.png')
            plt.close()

            plt.figure()
            plt.errorbar(goodTimes, goodReferences, yerr=goodRUnc, linestyle='None', fmt='-o')
            plt.xlabel('Time (BJD)')
            plt.ylabel('Total Flux')
            # plt.rc('grid', linestyle="-", color='black')
            # plt.grid(True)
            plt.title('Comparison Star Raw Flux Values ' + infoDict['date'])
            plt.savefig(infoDict['saveplot'] + 'temp/CompRawFlux' + pDict['pName'] + infoDict['date'] + '.png')
            plt.close()

            # Plots final reduced light curve (after the 3 sigma clip)
            plt.figure()
            plt.errorbar(goodPhases, goodFluxes, yerr=goodNormUnc, linestyle='None', fmt='-bo')
            plt.xlabel('Phase')
            plt.ylabel('Normalized Flux')
            # plt.rc('grid', linestyle="-", color='black')
            # plt.grid(True)
            plt.title(pDict['pName'] + ' Normalized Flux vs. Phase ' + infoDict['date'])
            plt.savefig(infoDict['saveplot'] + 'NormalizedFluxPhase' + pDict['pName'] + infoDict['date'] + '.png')
            plt.close()

            # Save normalized flux to text file prior to MCMC
            outParamsFile = open(infoDict['saveplot'] + 'NormalizedFlux' + pDict['pName'] + infoDict['date'] + '.txt', 'w+')
            outParamsFile.write(str("BJD") + ',' + str("Norm Flux") + ',' + str("Norm Err") + ',' + str("AM") + '\n')
            for ti, fi, erri, ami in zip(goodTimes, goodFluxes, goodNormUnc, goodAirmasses):
                outParamsFile.write(str(round(ti, 8)) + ',' + str(round(fi, 7)) + ',' + str(round(erri, 6)) + ',' + str(round(ami, 2)) + '\n')
            # CODE YIELDED DATA IN PREV LINE FORMAT
            outParamsFile.close()
            print('\nOutput File Saved')
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

            #Ask user for time format and convert it if not in BJD_TDB
            validTimeFormats = ['BJD_TDB', "MJD_UTC", "JD_UTC"]
            formatEntered = False
            print("\nNOTE: If your file is not in one of the following formats, please rereduce your data into one of the time formats recognized by EXOTIC.")
            while not formatEntered:
                print("Which of the following time formats is your data file stored in? (Type q to quit)")
                timeFormat = str(input("BJD_TDB / JD_UTC / MJD_UTC: "))
                if (timeFormat.upper()).strip() == 'Q':
                    sys.exit()
                elif (timeFormat.upper()).strip() not in validTimeFormats:
                    print("\nInvalid entry; please try again.")
                else:
                    formatEntered = True
            timeFormat = (timeFormat.upper()).strip()
            goodTimes = timeConvert(goodTimes, timeFormat, pDict, infoDict)

            #Ask user for flux units and convert to flux if in magnitude/millimagnitude
            validFluxFormats = ['flux', "magnitude", "millimagnitude"]
            formatEntered = False
            print("\nNOTE: If your file is not in one of the following formats, please rereduce your data into one of the time formats recognized by EXOTIC.")
            while not formatEntered:
                print("Which of the following units of flux is your data file stored in? (Type q to quit)")
                fluxFormat = str(input("flux / magnitude / millimagnitude: "))
                if (fluxFormat.upper()).strip() == 'Q':
                    sys.exit()
                elif (fluxFormat.lower()).strip() not in validFluxFormats:
                    print("\nInvalid entry; please try again.")
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


        ##########################
        # NESTED SAMPLING FITTING
        ##########################

        try:
            prior = {
                'rprs':pDict['rprs'],    # Rp/Rs
                'ars':pDict['aRs'],      # a/Rs
                'per':pDict['pPer'],     # Period [day]
                'inc':pDict['inc'],      # Inclination [deg]
                'u0': ld0[0], 'u1': ld1[0], 'u2': ld2[0], 'u3': ld3[0],  # limb darkening (nonlinear)
                'ecc': pDict['ecc'],     # Eccentricity
                'omega':0,          # Arg of periastron
                'tmid':pDict['midT'],    # time of mid transit [day]
                'a1': bestlmfit.parameters['a1'], #mid Flux
                'a2': bestlmfit.parameters['a2'], #Flux lower bound
            }
        except:
            prior = {
                'rprs':pDict['rprs'],    # Rp/Rs
                'ars':pDict['aRs'],      # a/Rs
                'per':pDict['pPer'],     # Period [day]
                'inc':pDict['inc'],      # Inclination [deg]
                'u0': ld0[0], 'u1': ld1[0], 'u2': ld2[0], 'u3': ld3[0],  # limb darkening (nonlinear)
                'ecc': pDict['ecc'],     # Eccentricity
                'omega':0,          # Arg of periastron
                'tmid':pDict['midT'],    # time of mid transit [day]
                'a1': goodFluxes.mean(), #max() - arrayFinalFlux.min(), #mid Flux
                'a2': 0,             #Flux lower bound
            }

        phase = (goodTimes-prior['tmid'])/prior['per']
        prior['tmid'] = pDict['midT'] + np.floor(phase).max()*prior['per']
        upper = pDict['midT']+ 25*pDict['midTUnc'] + np.floor(phase).max()*(pDict['pPer']+25*pDict['pPerUnc'])
        lower = pDict['midT']- 25*pDict['midTUnc'] + np.floor(phase).max()*(pDict['pPer']-25*pDict['pPerUnc'])

        if np.floor(phase).max()-np.floor(phase).min() == 0:
            print('Estimated mid-transit not in observation range (check priors or observation time)')
            print('start:', goodTimes.min())
            print('  end:', goodTimes.max())
            print('prior:', prior['tmid'])

        try:
            mybounds = {
                'rprs':[pDict['rprs']-3*pDict['rprsUnc'], pDict['rprs']+3*pDict['rprsUnc']],
                'tmid':[max(lower,goodTimes.min()),min(goodTimes.max(),upper)],
                'ars':[pDict['aRs']-5*pDict['aRsUnc'], pDict['aRs']+5*pDict['aRsUnc']],

                'a1':[bestlmfit.parameters['a1']*0.75, bestlmfit.parameters['a1']*1.25],
                'a2':[bestlmfit.parameters['a2']-0.25, bestlmfit.parameters['a2']+0.25],
            }
        except:
            mybounds = {
                'rprs':[pDict['rprs']-3*pDict['rprsUnc'], pDict['rprs']+3*pDict['rprsUnc']],
                'tmid':[max(lower,goodTimes.min()),min(goodTimes.max(),upper)],
                'ars':[pDict['aRs']-5*pDict['aRsUnc'], pDict['aRs']+5*pDict['aRsUnc']],
                'a1':[0, 3*np.nanmax(goodFluxes)],
                'a2':[-3,3],
            }

        # fitting method in elca.py
        myfit = lc_fitter(goodTimes, goodFluxes, goodNormUnc, goodAirmasses, prior, mybounds, mode='ns')

        # for k in myfit.bounds.keys():
        #     print("{:.6f} +- {}".format( myfit.parameters[k], myfit.errors[k]))

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
        ax_res.errorbar(myfit.phase, myfit.residuals/np.median(myfit.data), yerr=myfit.detrendederr ,color='gray', marker='o', markersize=5, linestyle='None', mec='None', alpha=0.75)
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
            f.savefig(infoDict['saveplot'] + 'FinalLightCurve' + pDict['pName'] + infoDict['date'] + ".pdf", bbox_inches="tight")
        except AttributeError:
            f.savefig(infoDict['saveplot'] + 'FinalLightCurve' + pDict['pName'] + infoDict['date'] + ".png", bbox_inches="tight")
        plt.close()

        ###################################################################################

        # triangle plot
        fig,axs = dynesty.plotting.cornerplot(myfit.results, labels=list(mybounds.keys()), quantiles_2d=[0.4,0.85], smooth=0.015, show_titles=True,use_math_text=True, title_fmt='.2e',hist2d_kwargs={'alpha':1,'zorder':2,'fill_contours':False})
        dynesty.plotting.cornerpoints(myfit.results, labels=list(mybounds.keys()), fig=[fig,axs[1:,:-1]],plot_kwargs={'alpha':0.1,'zorder':1,} )
        fig.savefig(infoDict['saveplot'] + 'temp/Triangle_{}_{}.png'.format(pDict['pName'], infoDict['date']))


        # write output to text file
        outParamsFile = open(infoDict['saveplot'] + 'FinalLightCurve' + pDict['pName'] + infoDict['date'] + '.csv', 'w+')
        outParamsFile.write('# FINAL TIMESERIES OF ' + pDict['pName'] + '\n')
        outParamsFile.write('# BJD_TDB,Orbital Phase,Flux,Uncertainty,Model,Airmass\n')

        phase = (myfit.time - myfit.parameters['tmid'] + 0.5*pDict['pPer'])/pDict['pPer'] % 1

        for bjdi, phasei, fluxi, fluxerri, modeli, ami in zip( myfit.time, phase, myfit.detrended, myfit.dataerr/myfit.airmass_model, myfit.transit, myfit.airmass_model):

            outParamsFile.write("{}, {}, {}, {}, {}, {}\n".format(bjdi, phasei, fluxi, fluxerri, modeli, ami))

        outParamsFile.close()

        #######################################################################
        # print final extracted planetary parameters
        #######################################################################

        print('*********************************************************')
        print('FINAL PLANETARY PARAMETERS\n')
        print('              Mid-Transit Time [BJD]: {} +- {} '.format(round_to_2(myfit.parameters['tmid'],myfit.errors['tmid']), round_to_2(myfit.errors['tmid'])))
        print('  Radius Ratio (Planet/Star) [Rp/Rs]: {} +- {} '.format(round_to_2(myfit.parameters['rprs'],myfit.errors['rprs']), round_to_2(myfit.errors['rprs'])))
        print(' Semi Major Axis/ Star Radius [a/Rs]: {} +- {} '.format(round_to_2(myfit.parameters['ars'],myfit.errors['ars']), round_to_2(myfit.errors['ars'])))
        print('               Airmass coefficient 1: {} +- {} '.format(round_to_2(myfit.parameters['a1'],myfit.errors['a1']), round_to_2(myfit.errors['a1'])))
        print('               Airmass coefficient 2: {} +- {} '.format(round_to_2(myfit.parameters['a2'],myfit.errors['a2']), round_to_2(myfit.errors['a2'])))
        print('The scatter in the residuals of the lightcurve fit is: {} %'.format(round_to_2(100. * np.std(myfit.residuals/np.median(myfit.data)))))
        print('\n*********************************************************')

        ##########
        # SAVE DATA
        ##########

        # TODO write as json
        # write output to text file
        outParamsFile = open(infoDict['saveplot'] + 'FinalParams' + pDict['pName'] + infoDict['date'] + '.txt', 'w+')
        outParamsFile.write('FINAL PLANETARY PARAMETERS\n')
        outParamsFile.write('')
        outParamsFile.write(' Mid-Transit Time: ' + str(round_to_2(myfit.parameters['tmid'],myfit.errors['tmid'])) + ' +/- ' + str(round_to_2(myfit.errors['tmid'])) + ' (BJD)\n')
        outParamsFile.write(' Ratio of Planet to Stellar Radius: ' + str(round_to_2(myfit.parameters['rprs'],myfit.errors['rprs'])) + ' +/- ' + str(
            round_to_2(myfit.errors['rprs'])) + ' (Rp/Rs)\n')
        outParamsFile.write(' transit depth uncertainty: ' + str(round_to_2(100. * 2. * myfit.parameters['rprs'] * myfit.errors['rprs'])) + ' (%)\n')
        outParamsFile.write(' airmass coefficient 1: ' + str(round_to_2(myfit.parameters['a1'],myfit.errors['a1'])) + ' +/- ' + str(round_to_2(myfit.errors['a1'])) + '\n')
        outParamsFile.write(' airmass coefficient 2: ' + str(round_to_2(myfit.parameters['a2'],myfit.errors['a2'])) + ' +/- ' + str(round_to_2(myfit.errors['a2'])) + '\n')
        outParamsFile.write(' scatter in the residuals of the lightcurve fit is: ' + str( round_to_2(100. * np.std(myfit.residuals/np.median(myfit.data)))) + '%\n')
        outParamsFile.close()
        print('\nFinal Planetary Parameters have been saved in ' + infoDict['saveplot'] + ' as '
              + pDict['pName'] + infoDict['date'] + '.txt' + '\n')

        # AAVSO Format
        userCode = infoDict['aavsonum']
        secuserCode = infoDict['secondobs']
        # else:
        outParamsFile = open(infoDict['saveplot'] + 'AAVSO' + pDict['pName'] + infoDict['date'] + '.txt', 'w+')
        outParamsFile.write('#TYPE=EXOPLANET\n')  # fixed
        outParamsFile.write('#OBSCODE=' + infoDict['aavsonum'] + '\n')  # UI
        outParamsFile.write('#SECONDARYOBSCODE=' + infoDict['secondobs'] + '\n')  # UI
        outParamsFile.write('#SOFTWARE=EXOTIC v' + __version__ + '\n')  # fixed
        outParamsFile.write('#DELIM=,\n')  # fixed
        outParamsFile.write('#DATE_TYPE=BJD_TDB\n')  # fixed
        outParamsFile.write('#OBSTYPE=' + infoDict['ctype'] + '\n')
        outParamsFile.write('#STAR_NAME=' + pDict['sName'] + '\n')  # code yields
        outParamsFile.write('#EXOPLANET_NAME=' + pDict['pName'] + '\n')  # code yields
        outParamsFile.write('#BINNING=' + infoDict['pixelbin'] + '\n')  # user input
        outParamsFile.write('#EXPOSURE_TIME=' + str(infoDict['exposure']) + '\n')  # UI
        outParamsFile.write('#FILTER=' + infoDict['filter'] + '\n')
        outParamsFile.write('#NOTES=' + infoDict['notes'] + '\n')
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

        print('\n************************')
        print('End of Reduction Process')
        print('************************')


if __name__ == "__main__":
    main()