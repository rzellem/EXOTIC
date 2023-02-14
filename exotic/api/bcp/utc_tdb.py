from __future__ import division
from __future__ import print_function
import urllib
import os
import datetime
import numpy as np
import sys
import math
import astropy
from astropy.coordinates import EarthLocation
from astropy.coordinates import get_body_barycentric_posvel
from astropy.time import Time

from functools import lru_cache

from . import PINT_erfautils as PINT
from barycorrpy.utils import get_stellar_data, find_hip, CalculatePositionVector
from barycorrpy.PhysicalConstants import *

fpath = os.path.join(os.path.dirname(__file__), 'data')


def staleness_check(file_time, now):
    '''
    Check whether the leap second file needs to be updated.

    Leap second updates generally announced shortly after December 31st or June 30th. Thus assuming a month lag for the list to be updated,
    check if a February 1st or August 1st has passed between when the file was downloaded and when the check is being run.

    INPUT:
        Both inputs are in datetime format
        file_time : Time the leap second file was last downloaded. Stored in leap_log.txt.
        now : UTC time right now.


    Return 1 if need to update, else return 0.

    '''
    if file_time.month > 8:
        year = file_time.year + 1
    else:
        year = file_time.year

    if (file_time < datetime.datetime(year, 2, 1) <= now) or (file_time < datetime.datetime(year, 9, 1) <= now):
        return 1
    else:
        return 0


def leap_download(ls_fpath, log_fpath):
    """
    Download the leap second file and update the log file.
    INPUT:
        ls_fpath : Path to where the file will be saved.
        log_fpath : Path to where the log file is created.

    """

    url = 'http://maia.usno.navy.mil/ser7/tai-utc.dat'

    if sys.version_info.major == 3:
        try:
            request = urllib.request.urlopen(url, timeout=10)  # timeout in 10 seconds if no response from server
            with open(ls_fpath, 'w') as f:
                f.write(request.read())
            flag = 0
            with open(log_fpath, 'w') as f:
                f.write(str(datetime.datetime.utcnow()))  # Write date of download in log file
        except (urllib.error.URLError, IOError):
            flag = 1
    else:
        import urllib2
        try:
            # should this be urllib2.urlretrieve?
            urllib2.urlretrieve(url, ls_fpath)
            flag = 0
            with open(log_fpath, 'w') as f:
                f.write(str(datetime.datetime.utcnow()))  # Write date of download in log file
        except (urllib2.HTTPError, IOError):
            flag = 1

    return flag


def leap_manage(utctime, fpath, leap_update):
    '''

    DEPRECATED:     v0.4.0 onwards
    20210303
    Using Astropy for leap seconds instead of maintaining/updating a separate files.
    --------------------------------------------------------------------------------------

    Calculates the offset between UTC and TAI from the leap second file.
    Also, checks for 'staleness' (how old is the file) of the leap second file.
    When the leap second file gets about 6 months old, it will update the leap second file, to check if there has been
    a new leap second announced.


    INPUT:
        Enter UTC time as Astropy Time Object. In UTC Scale.

        fpath : Path to where the file would be saved.
        leap_update : If True, when the leap second file is more than 6 months old it will attempt to download a new one.
                    If False, then will just give a warning message.

    OUTPUT:
        tai_utc : Leap second offset between TAI and UTC (Scalar)
        warning, error : Warning and Error messages

        Offset in seconds, between TAI and UTC.
    '''

    warning = []
    error = []

    # Location of saved leap second file
    ls_fpath = os.path.join(fpath, 'leap_sec.txt')

    # Location of log file to record how old is the leap second file
    log_fpath = os.path.join(fpath, 'leapsec_log.txt')

    # If log or leap file does not exist, then download the leap second file and create a log file
    if (os.path.isfile(log_fpath) == False or os.path.isfile(ls_fpath) == False):
        if leap_update == True:
            flag = leap_download(ls_fpath=ls_fpath, log_fpath=log_fpath)
            if flag == 0:
                warning += ['WARNING JD = ' + str(
                    utctime.jd) + ' : Downloaded leap second file from http://maia.usno.navy.mil/ser7/tai-utc.dat ']
            else:
                error += ['ERROR : JD = ' + str(
                    utctime.jd) + ' :  Unable to download leap second file. Check internet connection or leap_dir parameter']
                return 0, warning, error
        else:
            error += ['ERROR : JD = ' + str(
                utctime.jd) + ' :  LEAP SECOND FILE / LOG FILE DOES NOT EXIST. Please set leap_update = True to update file. Corrections may not be accurate ']
            return 0, warning, error


    # If file exists then check if need to update
    else:

        file_time = _read_log(log_fpath)
        now = datetime.datetime.utcnow()  # Current time

        if file_time < utctime.datetime:  # Check for file staleness if test time is after the last file update time

            # Check if need to update file
            if staleness_check(file_time, now) == 1:
                warning += ['WARNING : JD = ' + str(
                    utctime.jd) + ' :  Need to update leap second file. Corrections may not be accurate to 1 cm/s with old file']
                if leap_update == True:
                    flag = leap_download(ls_fpath=ls_fpath, log_fpath=log_fpath)
                    if flag == 0:
                        warning += ['JD = ' + str(
                            utctime.jd) + ' : Downloaded leap second file from http://maia.usno.navy.mil/ser7/tai-utc.dat ']
                    else:
                        error += ['ERROR : JD = ' + str(
                            utctime.jd) + ' :  Unable to download leap second file. Check internet connection or leap_dir parameter']
                        return 0, warning, error

                else:
                    error += ['ERROR : JD = ' + str(
                        utctime.jd) + ' :  Leap second file should be updated. Set leap_update = True to download file']
                    return 0, warning, error

    jd, offset = _read_leapfile(ls_fpath)

    if type(utctime) != astropy.time.core.Time:  # Check if it is in Astropy Time Object
        print("Input Time should be as Astropy Time Object")
        raise TypeError("Input Time should be as Astropy Time Object")
    else:
        JDUTC = utctime.jd

    tai_utc = offset[np.max(np.where(JDUTC >= np.array(jd)))]  # Leap second offset value to convert UTC to TAI

    return tai_utc, warning, error


def JDUTC_to_JDTDB(utctime, leap_update=True, fpath=fpath):
    '''
    v0.4.0 onwards
    20210303
    Using Astropy for leap seconds instead of maintaining/updating a separate files.

    DEPRECATED
    Convert JDUTC to JDTDB (Barycentric Dynamical Time)
    INPUT:
        utctime : Enter UTC time as Astropy Time Object. In UTC Scale.
        fpath : Path to where the file would be saved. Default is script directory.
        leap_update : If True, when the leap second file is more than 6 months old it will attempt to download a new one.
                      If False, then will just give a warning message.  Default is True. [Not used with versions >= v0.4.0]

    OUTPUT:
        JDTDB : Julian Date Barycentric Dynamic Time (Astropy Time object)
        JDTT: Julian Date Terrestrial Dynamic time (Astropy Time object)
        warning,error : Warning and Error message if any from the routine




    '''
    warning = []
    error = []

    if astropy.__version__ < '4.0.4':
        error += [
            "ERROR: Using Astropy version < 4.0.4 might not have the latest leap second file. Please update to a newer Astropy version to ensure updated leap seconds. "]

    return utctime.tdb, utctime.tt, warning, error


def JDUTC_to_BJDTDB(JDUTC,
                    starname='', hip_id=None, ra=None, dec=None, epoch=None, pmra=None, pmdec=None, px=None, rv=None,
                    obsname='', lat=0., longi=0., alt=0.,
                    ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__), 'data'), leap_update=True):
    '''
    Time conversion between JDUTC to BJDTDB. See Eastman et al. (2010)
    This code is precise to about 10 ms.

    Calling procedure for JDUTC_to_BJDTDB. Accepts vector time object (i.e., multiple observation JD values).

    INPUT:
        JDUTC : Can enter multiple times in Astropy Time object or as float. Will loop through and find barycentric velocity correction corresponding to those times.
                In UTC Scale. If using float, be careful about format and scale used.
        starname : Name of target. Will query SIMBAD database.
                                OR / AND
        hip_id : Hipparcos Catalog ID. (Integer) . Epoch will be taken to be Catalogue Epoch or J1991.25
                If specified then ra,dec,pmra,pmdec,px, and epoch need not be specified.
                                OR / AND
        ra, dec : RA and Dec of star [degrees].
        epoch : Epoch of coordinates in Julian Date. Default is J2000 or 2451545.
        pmra : Proper motion in RA [mas/year]. Eg. PMRA = d(RA)/dt * cos(dec). Default is 0.
        pmdec : Proper motion in Dec [mas/year]. Default is 0.
        px : Parallax of target [mas]. Default is 0.

        obsname : Name of Observatory as defined in Astropy EarthLocation routine. Can check list by EarthLocation.get_site_names().
                  If obsname is not used, then can enter lat,long,alt.
                                OR
        lat : Latitude of observatory in [degrees]. North (+ve) and South (-ve).
        longi : Longitude of observatory [degrees]. East (+ve) and West (-ve).
        alt : Altitude of observatory [m].

        rv : Radial Velocity of Target [m/s]. Default is 0.
        ephemeris : Name of Ephemeris to be used. List of Ephemeris as queried by jplephem. Default is DE430.
                    For first use Astropy will download the Ephemeris ( for DE430 ~100MB). Options for ephemeris inputs are
                    ['de432s','de430',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de423_for_mercury_and_venus/de423.bsp',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp']
        leap_dir : Directory where leap seconds file will be saved and maintained (STRING). Eg. '/Users/abc/home/savehere/'. Default is script directory. [Not used with versions >= v0.4.0]
        leap_update : If True, when the leap second file is more than 6 months old will attempt to download a new one.
                      If False, then will just give a warning message. Default is True. [Not used with versions >= v0.4.0]

    OUTPUT:
        corr_time : BJDTDB time
        warning : Warning and Error message from the routine.
        status : Status regarding warning and error message. Returns the following -
                0 - No warning or error.
                1 - Warning message.
                2 - Error message.

    Example:
    >>> from astropy.time import Time
    >>> JDUTC = Time(2458000, format='jd', scale='utc')
    >>> utc_tdb.JDUTC_to_BJDTDB(JDUTC,hip_id=8102, lat=-30.169283, longi=-70.806789, alt=2241.9)
    (array([ 2458000.00505211]), [[], []], 0)

    '''

    corr_time = []
    warning = []
    error = []
    status = 0

    # Check for JDUTC type
    if type(JDUTC) != Time:
        warning += [['Warning: Float JDUTC entered. Verify time scale (UTC) and format (JD)']]
        JDUTC = Time(JDUTC, format='jd', scale='utc')

    if JDUTC.isscalar:
        JDUTC = Time([JDUTC])

    star_par = {'ra': ra, 'dec': dec, 'pmra': pmra, 'pmdec': pmdec, 'px': px, 'rv': rv, 'epoch': epoch}
    star_simbad = {'ra': None, 'dec': None, 'pmra': None, 'pmdec': None, 'px': None, 'rv': None, 'epoch': None}
    star_hip = {}
    star_zero = {'ra': 0., 'dec': 0., 'pmra': 0., 'pmdec': 0., 'px': 0., 'rv': 0., 'epoch': 2451545.0}
    star_output = {}

    if starname:
        star_simbad, warning1 = get_stellar_data(starname)
        warning += warning1
    if hip_id:
        if starname:
            warning += ['Warning: Querying SIMBAD and Hipparcos Catalogue']
        star_hip = find_hip(hip_id)

    star_output = star_simbad.copy()
    star_output.update({k: star_hip[k] for k in star_hip if star_hip[k] is not None})
    star_output.update({k: star_par[k] for k in star_par if star_par[k] is not None})
    star_output.update({k: star_zero[k] for k in star_zero if star_output[k] is None})
    warning += ['Following are the stellar positional parameters being used - ', star_output]

    if obsname:
        loc = EarthLocation.of_site(obsname)
        lat = loc.lat.value
        longi = loc.lon.value
        alt = loc.height.value
        warning += [[
                        'Warning: Taking observatory coordinates from Astropy Observatory database. Verify precision. Latitude = %f  Longitude = %f  Altitude = %f' % (
                        lat, longi, alt)]]
    else:
        loc = EarthLocation.from_geodetic(longi, lat, height=alt)

    for jdutc in JDUTC:
        a = _JDUTC_to_BJDTDB(JDUTC=jdutc,
                             loc=loc,
                             ephemeris=ephemeris, **star_output)
        corr_time.append(a[0])
        warning.append(a[1])
        error.append(a[2])

    # Status messages to check for warning or error
    if not all(corr_time): error += ['Check inputs. Error in code']
    if any(error):   status |= 2
    if any(warning): status |= 1
    # Convert corrected from list to numpy array
    corr_time = np.array(corr_time)

    return corr_time, warning + error, status


def _JDUTC_to_BJDTDB(JDUTC,
                     ra=0.0, dec=0.0, epoch=2451545.0, pmra=0.0, pmdec=0.0, px=0.0, rv=0.0,
                     loc=None,
                     ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__), 'data'), leap_update=True):
    '''
    Time conversion between JDUTC to BJDTDB. See Eastman et al. (2010)
    This code is precise to about 10 ms

    See JDUTC_to_BJDTDB() for parameter description.

    '''

    # Convert times to obtain TDB and TT
    JDTDB, JDTT, warning, error = JDUTC_to_JDTDB(JDUTC)
    clock_corr = (JDTDB.jd - JDUTC.jd) * 86400.

    ##### NUTATION, PRECESSION, ETC. #####

    r_pint, v_pint = PINT.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)

    r_eci = r_pint[0]  # [m]
    v_eci = v_pint[0]  # [m/s]

    ##### EPHEMERIDES #####

    earth_geo = get_body_barycentric_posvel('earth', JDTDB, ephemeris=ephemeris)  # [km]
    r_obs = r_eci + earth_geo[0].xyz.value * 1000.  # [m]
    v_geo = earth_geo[1].xyz.value * 1000. / 86400.  # [m/s]

    # Relativistic Addition of Velocities
    v_obs = (v_eci + v_geo) / (1. + v_eci * v_geo / c ** 2)  # [m/s]

    # calculate the Einstein delay relative to the geocenter
    # (TDB accounts for Einstein delay to geocenter)
    einstein_corr = np.sum(r_eci * v_geo) / (c * c)

    ##### Convert Star RA DEC to R0hat vector #####

    r0hat = np.array([math.cos(ra * np.pi / 180.) * math.cos(dec * np.pi / 180.),
                      math.sin(ra * np.pi / 180.) * math.cos(dec * np.pi / 180.),
                      math.sin(dec * np.pi / 180.)])
    # Eq 14 to 17 from Wright and Eastman 2014
    up = [0., 0., 1.]
    east = np.cross(up, r0hat)
    east = east / math.sqrt(sum(east * east))
    north = np.cross(r0hat, east)
    mu = (
                     pmra * east + pmdec * north) / pctoau / 1000  # Divided by 1000 since the Proper motion is in milli-arcseconds.

    ##### Stellar position corrected for motion #####

    epoch0 = 2000. + (epoch - 2451545.) / 365.25
    yearnow = 2000. + (JDTDB.jd - 2451545.) / 365.25

    T = yearnow - epoch0  # [years]
    vpi = rv / 1.e3 * kmstoauyr * (px / 1.e3 / pctoau)  # [rad/yr]
    vel = mu + vpi * r0hat  # [rad/yr] (m in AA)
    r = r0hat + vel * T  # [rad]    (p1 in AA)
    rhat = r / math.sqrt(sum(r * r))

    # Geometric Correction
    geo_corr = np.sum(r_obs * rhat) / c

    delta_t = geo_corr + clock_corr + einstein_corr
    result = JDUTC.jd + delta_t / 86400.

    return result, warning, error


def JDUTC_to_HJDTDB(JDUTC,
                    obsname='', lat=0., longi=0., alt=0.,
                    ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__), 'data'), leap_update=True):
    '''
    Based on Eastman et al. 2010, but converts JDUTC to JDTDB, and then moves the
    clock to the time of emission from the heliocenter.
    This is to serve as a time stamp for the solar RVs as observed by PRV spectrographs (Wright and Kanodia 2020).
    It includes the JDUTC to JDTDB time conversion,  as well as the light travel time from observatory
    to solar center and the Einstein delay.

    Precision has not been explicity tested.

    Calling procedure for JDUTC_to_HJDTDB. Accepts vector time object (i.e., multiple observation JD values).

    INPUT:
        JDUTC : Can enter multiple times in Astropy Time object or as float. Will loop through and find barycentric velocity correction corresponding to those times.
                In UTC Scale. If using float, be careful about format and scale used.


        obsname : Name of Observatory as defined in Astropy EarthLocation routine. Can check list by EarthLocation.get_site_names().
                  If obsname is not used, then can enter lat,long,alt.
                                OR
        lat : Latitude of observatory in [degrees]. North (+ve) and South (-ve).
        longi : Longitude of observatory [degrees]. East (+ve) and West (-ve).
        alt : Altitude of observatory [m].

        ephemeris : Name of Ephemeris to be used. List of Ephemeris as queried by jplephem. Default is DE430.
                    For first use Astropy will download the Ephemeris ( for DE430 ~100MB). Options for ephemeris inputs are
                    ['de432s','de430',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de423_for_mercury_and_venus/de423.bsp',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp']
        leap_dir : Directory where leap seconds file will be saved and maintained (STRING). Eg. '/Users/abc/home/savehere/'. Default is script directory.
        [Not used with versions >= v0.4.0]
        leap_update : If True, when the leap second file is more than 6 months old will attempt to download a new one.
                      If False, then will just give a warning message. Default is True. [Not used with versions >= v0.4.0]

    OUTPUT:
        corr_time : BJDTDB time
        warning : Warning and Error message from the routine.
        status : Status regarding warning and error message. Returns the following -
                0 - No warning or error.
                1 - Warning message.
                2 - Error message.

    Example:
    >>> from astropy.time import Time
    >>> JDUTC = Time(2458000, format='jd', scale='utc')
    >>> utc_tdb.JDUTC_to_HJDTDB(JDUTC,lat=-30.169283, longi=-70.806789, alt=2241.9)
    (array([ 2458000.00662602]), [[], []], 0)

    '''

    corr_time = []
    warning = []
    error = []
    status = 0

    # Check for JDUTC type
    if type(JDUTC) != Time:
        warning += [['Warning: Float JDUTC entered. Verify time scale (UTC) and format (JD)']]
        JDUTC = Time(JDUTC, format='jd', scale='utc')

    if JDUTC.isscalar:
        JDUTC = Time([JDUTC])

    if obsname:
        loc = EarthLocation.of_site(obsname)
        lat = loc.lat.value
        longi = loc.lon.value
        alt = loc.height.value
        warning += [[
                        'Warning: Taking observatory coordinates from Astropy Observatory database. Verify precision. Latitude = %f  Longitude = %f  Altitude = %f' % (
                        lat, longi, alt)]]
    else:
        loc = EarthLocation.from_geodetic(longi, lat, height=alt)

    for jdutc in JDUTC:
        a = _JDUTC_to_HJDTDB(JDUTC=jdutc,
                             loc=loc,
                             ephemeris=ephemeris, leap_dir=leap_dir, leap_update=leap_update)
        corr_time.append(a[0])
        warning.append(a[1])
        error.append(a[2])

    # Status messages to check for warning or error
    if not all(corr_time): error += ['Check inputs. Error in code']
    if any(error):   status |= 2
    if any(warning): status |= 1
    # Convert corrected from list to numpy array
    corr_time = np.array(corr_time)

    return corr_time, warning + error, status


def _JDUTC_to_HJDTDB(JDUTC, loc,
                     ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__), 'data'), leap_update=True):
    '''
    Based on Eastman et al. 2010, but converts JDUTC to JDTDB, and then moves the
    clock to the time of emission from the heliocenter.
    This is to serve as a time stamp for the solar RVs as observed by PRV spectrographs (Wright and Kanodia 2020).
    It includes the JDUTC to JDTDB time conversion,  as well as the light travel time from observatory
    to solar center and the Einstein delay.

    Precision has not been explicity tested.


    See JDUTC_to_HJDTDB() for parameter description.

    '''

    # Convert times to obtain TDB and TT
    JDTDB, JDTT, warning, error = JDUTC_to_JDTDB(JDUTC)
    clock_corr = (JDTDB.jd - JDUTC.jd) * 86400.

    ##### NUTATION, PRECESSION, ETC. #####

    r_pint, v_pint = PINT.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)

    r_eci = r_pint[0]  # [m]
    v_eci = v_pint[0]  # [m/s]

    ##### EPHEMERIDES #####

    earth_geo = get_body_barycentric_posvel('earth', JDTDB, ephemeris=ephemeris)  # [km]
    PosVector_EarthSSB = r_eci + earth_geo[0].xyz.value * 1000.  # [m]
    v_geo = earth_geo[1].xyz.value * 1000. / 86400.  # [m/s]

    # Relativistic Addition of Velocities
    VelVector_EarthSSB = (v_eci + v_geo) / (1. + np.sum(v_eci * v_geo) / c ** 2)  # [m/s]

    ## Ignoring retarded time for Solar vectors

    solar_ephem = get_body_barycentric_posvel('sun', JDTDB, ephemeris=ephemeris)

    PosVector_SolSSB = solar_ephem[0].xyz.value * 1000.  # [m]
    VelVector_SolSSB = solar_ephem[1].xyz.value * 1000. / 86400.  # [m/s]

    PosVector_SolEarth, PosMag_SolEarth, PosHat_SolEarth = CalculatePositionVector(r1=PosVector_SolSSB,
                                                                                   r2=PosVector_EarthSSB)

    geo_corr_SSB = - PosMag_SolEarth / c

    # calculate the Einstein delay relative to the geocenter
    # (TDB accounts for Einstein delay to geocenter)
    einstein_corr = np.sum(r_eci * v_geo) / (c * c)

    delta_t = geo_corr_SSB + clock_corr + einstein_corr
    result = JDUTC.jd + delta_t / 86400.

    return result, warning, error


@lru_cache(maxsize=200)
def _read_log(log_fpath):
    with open(log_fpath, 'r') as f:
        file_time = datetime.datetime.strptime(f.readline(),
                                               '%Y-%m-%d %H:%M:%S.%f')  # Read date of download in log file
    return file_time


@lru_cache(maxsize=200)
def _read_leapfile(ls_fpath):
    """
    Parse in leap second file
    """

    f = open(ls_fpath, 'r')
    jd = []
    offset = []

    for line in f:
        a = (line.split())
        jd.append(float(a[4]))
        offset.append(float(a[6]))
    f.close()

    return jd, offset