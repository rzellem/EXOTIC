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

from io import StringIO
import json
import numpy as np
import os
import pandas
import requests
import time
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord
from tenacity import retry, retry_if_exception_type, stop_after_attempt, \
    wait_exponential

# constants
AU = const.au                                                       # m
R_SUN = const.R_sun                                                 # m
R_JUP = const.R_jup                                                 # m

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
                "Star Surface Gravity (-) Uncertainty": self.pl_dict['loggUncNeg']
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
        # table access protocol query

        # build url
        uri_full = base_url
        for k in query:
            if k != "format":
                uri_full += f"{k} {query[k]} "

        uri_full = f"{uri_full[:-1]} &format={query.get('format', 'csv')}"
        uri_full = uri_full.replace(' ', '+')
        # log_info(uri_full)
        # log.debug(uri_full)

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

        new_index = [planet.lower().replace(' ', '').replace('-', '') for planet in default.pl_name.values]

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
                for key, value in planets.items():
                    if self.planet.lower().replace(' ', '').replace('-', '') == key:
                        self.planet = value
                        break
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
            return self.planet, False

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
                rp = data['pl_radj'] * R_JUP.value
                rperr = np.sqrt(np.abs(data['pl_radjerr1'] * data['pl_radjerr2'])) * R_JUP.value
                rs = data['st_rad'] * R_SUN.value
                rserr = np.sqrt(np.abs(data['st_raderr1'] * data['st_raderr2'])) * R_SUN.value
                rprserr = ((rperr / rs) ** 2 + (-rp * rserr / rs ** 2) ** 2) ** 0.5
                rprs = rp / rs

        if data['pl_ratdor'] < 1 or np.isnan(data['pl_ratdor']):
            data['pl_ratdor'] = pow((data['pl_orbper'] / 365) ** 2, 1 / 3) / (data['st_rad'] * R_SUN.to('au')).value

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
            'omega': data.get('pl_orblper', 0),
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
