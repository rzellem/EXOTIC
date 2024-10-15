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
import astropy.constants as const
from astropy.coordinates import SkyCoord
import astropy.units as u
from io import StringIO
import json
import numpy as np
import os
import pandas
import re
import requests
import time
import urllib.parse
from tenacity import retry, retry_if_exception_type, stop_after_attempt, \
    wait_exponential

# constants
AU = const.au # m
R_SUN = const.R_sun # m
R_JUP = const.R_jup # m

# CALCULATED VALUES
G = const.G.to(AU**3 / (const.M_sun * u.day**2))                    # AU^3 /(msun * day^2)
SA = lambda m, p: (G * m * p ** 2. / (4. * np.pi ** 2.)) ** (1. / 3.)  # Keplerian semi-major axis (au)


def result_if_max_retry_count(retry_state):
    pass


class NASAExoplanetArchive:

    def __init__(self, planet=None, candidate=False):
        self.planet = planet
        # self.candidate = candidate
        self.pl_dict = None

        # CONFIGURATIONS
        self.requests_timeout = 16, 512  # connection timeout, response timeout in secs.

    def planet_info(self, fancy=False):
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
                "Argument of Periastron (deg)": self.pl_dict['omega'],
                "Star Effective Temperature (K)": self.pl_dict['teff'],
                "Star Effective Temperature (+) Uncertainty": self.pl_dict['teffUncPos'],
                "Star Effective Temperature (-) Uncertainty": self.pl_dict['teffUncNeg'],
                "Star Metallicity ([FE/H])": self.pl_dict['met'],
                "Star Metallicity (+) Uncertainty": self.pl_dict['metUncPos'],
                "Star Metallicity (-) Uncertainty": self.pl_dict['metUncNeg'],
                "Star Surface Gravity (log(g))": self.pl_dict['logg'],
                "Star Surface Gravity (+) Uncertainty": self.pl_dict['loggUncPos'],
                "Star Surface Gravity (-) Uncertainty": self.pl_dict['loggUncNeg'],
                "Star Distance (pc)": self.pl_dict['dist'],
                "Star Proper Motion RA (mas/yr)": self.pl_dict['pm_ra'],
                "Star Proper Motion DEC (mas/yr)": self.pl_dict['pm_dec']
            }

            return json.dumps(flabels, indent=4)
        else:
            self.planet, candidate = self._new_scrape(filename="eaConf.json")

            if not candidate:
                with open("eaConf.json", "r") as confirmed:
                    data = json.load(confirmed)
                    planets = [data[i]['pl_name'] for i in range(len(data))]
                    idx = planets.index(self.planet)
                    self._get_params(data[idx])
                    print(f"Successfully found {self.planet} in the NASA Exoplanet Archive!")

            return self.planet, candidate, self.pl_dict

    @staticmethod
    def dataframe_to_jsonfile(dataframe, filename):
        jsondata = json.loads(dataframe.to_json(orient='table', index=False))
        with open(filename, "w") as f:
            f.write(json.dumps(jsondata['data'], indent=4))

    def _tap_query(self, base_url, query, dataframe=True):
        # Build the ADQL query string
        adql_query = ' '.join(f"{k} {v}" for k, v in query.items() if k != "format")
        adql_query = adql_query.strip()  # Remove any trailing space

        # URL-encode the entire ADQL query
        encoded_query = urllib.parse.quote(adql_query)

        # Build the full URL with the encoded query
        # Since base_url already ends with 'query=', we append the encoded query directly
        uri_full = f"{base_url}{encoded_query}&format={query.get('format', 'csv')}"

        # Send the request
        response = requests.get(uri_full, timeout=self.requests_timeout)

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
            uri_ipac_query["where"] += f" and pl_name = '{self.planet}'"

        default = self._tap_query(uri_ipac_base, uri_ipac_query)

        if len(default) == 0:
            return False
        else:
            return True

    @retry(stop=stop_after_attempt(3),
           wait=wait_exponential(multiplier=1, min=10, max=20),
           retry=(retry_if_exception_type(requests.exceptions.RequestException) |
                  retry_if_exception_type(ConnectionError)),
           retry_error_callback=result_if_max_retry_count)
    def planet_names(self, filename="pl_names.json"):
        uri_ipac_base = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
        uri_ipac_query = {
            "select": "pl_name",
            "from": "ps",
            "where": "tran_flag = 1 and default_flag = 1",
            "order by": "pl_pubdate desc",
            "format": "csv"
        }
        default = self._tap_query(uri_ipac_base, uri_ipac_query)

        new_index = [re.sub(r'[^a-zA-Z0-9]', '', planet.lower()) for planet in default.pl_name.values]

        planets = dict(zip(new_index, default.pl_name.values))
        with open(filename, "w") as f:
            f.write(json.dumps(planets, indent=4))

    @retry(stop=stop_after_attempt(3),
           wait=wait_exponential(multiplier=1, min=17, max=1024),
           retry=(retry_if_exception_type(requests.exceptions.RequestException) |
                  retry_if_exception_type(ConnectionError)))
    def _new_scrape(self, filename="eaConf.json"):

        # scrape_new()
        uri_ipac_base = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
        uri_ipac_query = {
            "select": "pl_name,hostname,tran_flag,pl_massj,pl_radj,pl_radjerr1,pl_radjerr2,"
                      "pl_ratdor,pl_ratdorerr1,pl_ratdorerr2,pl_orbincl,pl_orbinclerr1,pl_orbinclerr2,"
                      "pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbeccen,"
                      "pl_orblper,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,"
                      "pl_trandep,pl_trandeperr1,pl_trandeperr2,"
                      "pl_ratror,pl_ratrorerr1,pl_ratrorerr2,"
                      "st_teff,st_tefferr1,st_tefferr2,st_met,st_meterr1,st_meterr2,"
                      "sy_pmra,sy_pmdec,sy_dist,"
                      "st_logg,st_loggerr1,st_loggerr2,st_mass,st_rad,st_raderr1,st_raderr2,ra,dec,pl_pubdate",
            "from": "ps",
            "where": "tran_flag = 1 and default_flag = 1",
            "order by": "pl_pubdate desc",
            "format": "csv"
        }

        if not os.path.exists('pl_names.json') or time.time() - os.path.getmtime('pl_names.json') > 2592000:
            self.planet_names(filename="pl_names.json")
        if os.path.exists('pl_names.json'):
            with open("pl_names.json", "r") as f:
                planets = json.load(f)
                planet_key = re.sub(r'[^a-zA-Z0-9]', '', self.planet.lower())

                planet_exists = planets.get(planet_key, False)

                if planet_exists:
                    self.planet = planet_exists

        print(f"\nLooking up {self.planet} on the NASA Exoplanet Archive. Please wait....")

        if self.planet:
            uri_ipac_query["where"] += f" and pl_name = '{self.planet}'"

        default = self._tap_query(uri_ipac_base, uri_ipac_query)

        # fill in missing columns
        uri_ipac_query['where'] = 'tran_flag=1'

        if self.planet:
            uri_ipac_query["where"] += f" and pl_name = '{self.planet}'"

        extra = self._tap_query(uri_ipac_base, uri_ipac_query)

        if len(default) == 0:
            self.planet = input(f"Cannot find target ({self.planet}) in NASA Exoplanet Archive."
                                f"\nPlease go to https://exoplanetarchive.ipac.caltech.edu to check naming and"
                                "\nre-enter the planet's name or type 'candidate' if this is a planet candidate: ")
            if self.planet.strip().lower() == 'candidate':
                self.planet = user_input("\nPlease enter candidate planet's name: ", type_=str)
                return self.planet, True
            else:
                return self._new_scrape(filename="eaConf.json")
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
                    if nans[k].iloc[0]:  # if col value is nan
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
            return self.planet, False

    def _get_params(self, data):
        # Initialize variables with default values
        rprs = np.nan
        rprserr = np.nan
        rp = np.nan
        rperr = np.nan
        rs = np.nan
        rserr = np.nan

        # compute Rp/Rs
        if 'pl_trandep' in data and data['pl_trandep'] is not None:
            rprs = np.sqrt(data['pl_trandep'] / 100.)
            rprserr = np.sqrt(np.abs((data['pl_trandeperr1'] / 100.) * (data['pl_trandeperr2'] / 100.))) / (2. * rprs)

        elif 'pl_ratror' in data and data['pl_ratror'] is not None:
            rprs = data['pl_ratror']
            rprserr = np.sqrt(np.abs(data['pl_ratrorerr1'] * data['pl_ratrorerr2']))

        # check if rprs is still 0 or nan
        if np.isnan(rprs) or rprs <= 0.:

            # compute with stellar and planetary radius
            if 'pl_radj' in data and data['pl_radj'] is not None and 'st_rad' in data and data['st_rad'] is not None:
                rp = data['pl_radj'] * R_JUP.value
                rperr = np.sqrt(np.abs(data['pl_radjerr1'] * data['pl_radjerr2'])) * R_JUP.value
                rs = data['st_rad'] * R_SUN.value if 'st_rad' in data and data['st_rad'] is not None else np.nan
                rserr = np.sqrt(np.abs(data['st_raderr1'] * data['st_raderr2'])) * R_SUN.value if 'st_raderr1' in data and data['st_raderr1'] is not None else np.nan
                if not np.isnan(rs) and not np.isnan(rp):
                    rprserr = np.sqrt((rperr / rs) ** 2 + (-rp * rserr / rs ** 2) ** 2)
                    rprs = rp / rs

        # compute a/Rs
        if 'pl_ratdor' not in data or data['pl_ratdor'] is None or data['pl_ratdor'] < 1.:
            if 'pl_orbper' in data and data['pl_orbper'] is not None and 'st_rad' in data and data['st_rad'] is not None:
                data['pl_ratdor'] = (data['pl_orbper'] / 365.) ** (2. / 3.) / (data['st_rad'] * R_SUN.to('au')).value
            else:
                print("WARNING: a/Rs could not be calculated due to missing or invalid orbital period or stellar radius.")

        self.pl_dict = {
            'ra': float(data['ra']) if 'ra' in data and data['ra'] is not None else np.nan,
            'dec': float(data['dec']) if 'dec' in data and data['dec'] is not None else np.nan,
            'pName': str(data['pl_name']),
            'sName': str(data['hostname']),
            'pPer': float(data['pl_orbper']) if 'pl_orbper' in data and data['pl_orbper'] is not None else np.nan,
            'pPerUnc': float(np.sqrt(np.abs(data['pl_orbpererr1'] * data['pl_orbpererr2']))) if 'pl_orbpererr1' in data and 'pl_orbpererr2' in data and data['pl_orbpererr1'] is not None and data['pl_orbpererr2'] is not None else np.nan,
            'midT': float(data['pl_tranmid']) if 'pl_tranmid' in data and data['pl_tranmid'] is not None else np.nan,
            'midTUnc': float(np.sqrt(np.abs(data['pl_tranmiderr1'] * data['pl_tranmiderr2']))) if 'pl_tranmiderr1' in data and 'pl_tranmiderr2' in data and data['pl_tranmiderr1'] is not None and data['pl_tranmiderr2'] is not None else np.nan,
            'rprs': float(rprs) if not np.isnan(rprs) else np.nan,
            'rprsUnc': float(rprserr) if not np.isnan(rprserr) else np.nan,
            'aRs': float(data['pl_ratdor']) if 'pl_ratdor' in data and data['pl_ratdor'] is not None else np.nan,
            'aRsUnc': float(np.sqrt(np.abs(data.get('pl_ratdorerr1', 1) * data['pl_ratdorerr2']))) if 'pl_ratdorerr2' in data and data['pl_ratdorerr2'] is not None else 0.1,
            'inc': float(data['pl_orbincl']) if 'pl_orbincl' in data and data['pl_orbincl'] is not None else np.nan,
            'incUnc': float(np.sqrt(np.abs(data['pl_orbinclerr1'] * data['pl_orbinclerr2']))) if 'pl_orbinclerr1' in data and 'pl_orbinclerr2' in data and data['pl_orbinclerr1'] is not None and data['pl_orbinclerr2'] is not None else 0.1,
            'omega': float(data.get('pl_orblper', 0)),
            'ecc': float(data.get('pl_orbeccen', 0)),
            'teff': float(data['st_teff']) if 'st_teff' in data and data['st_teff'] is not None else np.nan,
            'teffUncPos': float(data['st_tefferr1']) if 'st_tefferr1' in data and data['st_tefferr1'] is not None else np.nan,
            'teffUncNeg': float(data['st_tefferr2']) if 'st_tefferr2' in data and data['st_tefferr2'] is not None else np.nan,
            'met': float(data['st_met']) if 'st_met' in data and data['st_met'] is not None else np.nan,
            'metUncPos': float(max(0.01, data['st_meterr1'])) if 'st_meterr1' in data and data['st_meterr1'] is not None else 0.01,
            'metUncNeg': float(min(-0.01, data['st_meterr2'])) if 'st_meterr2' in data and data['st_meterr2'] is not None else -0.01,
            'logg': float(data['st_logg']) if 'st_logg' in data and data['st_logg'] is not None else np.nan,
            'loggUncPos': float(data['st_loggerr1']) if 'st_loggerr1' in data and data['st_loggerr1'] is not None else np.nan,
            'loggUncNeg': float(data['st_loggerr2']) if 'st_loggerr2' in data and data['st_loggerr2'] is not None else np.nan,
            'dist': float(data['sy_dist']) if 'sy_dist' in data and data['sy_dist'] is not None else np.nan,
            'pm_dec': float(data['sy_pmdec']) if 'sy_pmdec' in data and data['sy_pmdec'] is not None else np.nan,
            'pm_ra': float(data['sy_pmra']) if 'sy_pmra' in data and data['sy_pmra'] is not None else np.nan
        }

        if self.pl_dict['aRsUnc'] == 0:
            self.pl_dict['aRsUnc'] = 0.1

        if self.pl_dict['incUnc'] == 0:
            self.pl_dict['incUnc'] = 0.1


# temp, need to remove user_input from NEA
def user_input(prompt, type_, values=None):
    while True:
        try:
            result = type_(input(prompt))
            # log.debug(f"{prompt}{result}")
        except ValueError:
            print("Sorry, not a valid datatype.")
            continue
        if type_ == str and values is not None:
            result = result.lower().strip()
            if result not in values:
                print("Sorry, your response was not valid.")
            else:
                return result
        elif type_ == int and values is not None:
            if result not in values:
                print("Sorry, your response was not valid.")
            else:
                return result
        else:
            return result
