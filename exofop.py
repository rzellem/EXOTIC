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
import pandas as pd
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


class ExoFOPArchive:

    def __init__(self, planet=None, candidate=True):
        self.planet = planet
        self.candidate = candidate
        self.pl_dict = None
        self.requests_timeout = 16, 512  # connection timeout, response timeout in secs.  (CONFIGURATIONS)

    @staticmethod
    def dataframe_to_jsonfile(dataframe, filename):
        jsondata = json.loads(dataframe.to_json(orient='table', index=False))
        with open(filename, "w") as f:
            f.write(json.dumps(jsondata['data'], indent=4))

    # Function to take a host name and return a TIC ID;
    # provides an error if a TIC ID can't be matched
    # (for any confirmed exoplanet, this should return a TIC ID):
    # host star name
    def get_tic(target_name):
        url = 'https://exofop.ipac.caltech.edu/tess/gototicid.php?target=' + target_name + '&json'
        r = requests.get(url)
        assert r.status_code == 200, "gototicid request failed"
        data = r.json()
        assert data['status'] == "OK", "Could not find target: " + target_name
        return (data['TIC'])

    # Function that takes a string of either "stellar" or "planet" to select which table you want,
    # followed by the TIC ID (from above function). Returns all parameters as a pandas dataframe,
    # or False if there are no parameters associated with that TIC:
    def get_parameters(table, target):
        if table == 'stellar':
            url = 'https://exofop.ipac.caltech.edu/tess/download_stellar.php?id=' + target
        elif table == 'planet':
            url = 'https://exofop.ipac.caltech.edu/tess/download_planet.php?id=' + target
        else:
            print("No valid table named", table)
        r = requests.get(url)
        assert r.status_code == 200, "get_parameters request failed"
        query_df = pd.read_csv(StringIO(r.text), sep='|')
        if query_df.empty:
            return False
        return query_df

    # def planet_info(self, fancy=False):
    #     # fancy keys matching inits fil
    #     if fancy:
    #         coord = SkyCoord(ra=self.pl_dict['ra'] * u.degree, dec=self.pl_dict['dec'] * u.degree)
    #         rastr = coord.to_string('hmsdms', sep=':').split(' ')[0]
    #         decstr = coord.to_string('hmsdms', sep=':').split(' ')[1]
    #
    #         flabels = {
    #             "Target Star RA": rastr,
    #             "Target Star Dec": decstr,
    #             "Planet Name": self.pl_dict['pName'],
    #             "Host Star Name": self.pl_dict['sName'],
    #             "Orbital Period (days)": self.pl_dict['pPer'],
    #             "Orbital Period Uncertainty": self.pl_dict['pPerUnc'],
    #             "Published Mid-Transit Time (BJD-UTC)": self.pl_dict['midT'],
    #             "Mid-Transit Time Uncertainty": self.pl_dict['midTUnc'],
    #             "Ratio of Planet to Stellar Radius (Rp/Rs)": self.pl_dict['rprs'],
    #             "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": self.pl_dict['rprsUnc'],
    #             "Ratio of Distance to Stellar Radius (a/Rs)": self.pl_dict['aRs'],
    #             "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": self.pl_dict['aRsUnc'],
    #             "Orbital Inclination (deg)": self.pl_dict['inc'],
    #             "Orbital Inclination (deg) Uncertainty": self.pl_dict['incUnc'],
    #             "Orbital Eccentricity (0 if null)": self.pl_dict['ecc'],
    #             "Argument of Periastron (deg)": self.pl_dict['omega'],
    #             "Star Effective Temperature (K)": self.pl_dict['teff'],
    #             "Star Effective Temperature (+) Uncertainty": self.pl_dict['teffUncPos'],
    #             "Star Effective Temperature (-) Uncertainty": self.pl_dict['teffUncNeg'],
    #             "Star Metallicity ([FE/H])": self.pl_dict['met'],
    #             "Star Metallicity (+) Uncertainty": self.pl_dict['metUncPos'],
    #             "Star Metallicity (-) Uncertainty": self.pl_dict['metUncNeg'],
    #             "Star Surface Gravity (log(g))": self.pl_dict['logg'],
    #             "Star Surface Gravity (+) Uncertainty": self.pl_dict['loggUncPos'],
    #             "Star Surface Gravity (-) Uncertainty": self.pl_dict['loggUncNeg']
    #         }
    #
    #         return json.dumps(flabels, indent=4)
    #     else:
    #         self.planet, candidate = self._new_scrape(filename="eaConf.json")
    #
    #         if not candidate:
    #             with open("eaConf.json", "r") as confirmed:
    #                 data = json.load(confirmed)
    #                 planets = [data[i]['pl_name'] for i in range(len(data))]
    #                 idx = planets.index(self.planet)
    #                 self._get_params(data[idx])
    #                 print(f"Successfully found {self.planet} in ExoFOP!")
    #
    #         return self.planet, candidate, self.pl_dict

    # def _tap_query(self, base_url, query, dataframe=True):
    #     # table access protocol query
    #     # build url
    #     uri_full = base_url
    #     for k in query:
    #         if k != "format":
    #             uri_full += f"{k} {query[k]} "
    #
    #     uri_full = f"{uri_full[:-1]} &format={query.get('format', 'csv')}"
    #     uri_full = uri_full.replace(' ', '+')
    #     # log_info(uri_full)
    #     # log.debug(uri_full)
    #
    #     response = requests.get(uri_full, timeout=self.requests_timeout)
    #     assert response.status_code == 200, "uri_fill request failed"
    #
    #     if dataframe:
    #         return pd.read_csv(StringIO(response.text))
    #     else:
    #         return response.text
    #
    # def resolve_name(self):
    #     uri_ipac_base = "https://exofop.ipac.caltech.edu/tess/"
    #     uri_ipac_query = {
    #         "select": "pl_name,hostname",
    #         "from": "ps",  # Table name
    #         "where": "tran_flag = 1",
    #         "order by": "pl_pubdate desc",
    #         "format": "csv"
    #     }
    #
    #     if self.planet:
    #         uri_ipac_query["where"] += f" and pl_name = '{self.planet}'"
    #
    #     default = self._tap_query(uri_ipac_base, uri_ipac_query)
    #
    #     if len(default) == 0:
    #         return False
    #     else:
    #         return True
    #
    #
    # @retry(stop=stop_after_attempt(3),
    #        wait=wait_exponential(multiplier=1, min=10, max=20),
    #        retry=(retry_if_exception_type(requests.exceptions.RequestException) |
    #               retry_if_exception_type(ConnectionError)),
    #        retry_error_callback=result_if_max_retry_count)
    #
    # def planet_names(self, filename="pl_names.json"):
    #     uri_ipac_base = "https://exofop.ipac.caltech.edu/tess/"
    #     uri_ipac_query = {
    #         "select": "TIC ID",
    #         "from": "ps",
    #         "where": "tran_flag = 1 and default_flag = 1",
    #         "order by": "pl_pubdate desc",
    #         "format": "csv"
    #     }
    #     default = self._tap_query(uri_ipac_base, uri_ipac_query)
    #
    #     new_index = [planet.lower().replace(' ', '').replace('-', '') for planet in default.pl_name.values]
    #
    #     planets = dict(zip(new_index, default.pl_name.values))
    #     with open(filename, "w") as f:
    #         f.write(json.dumps(planets, indent=4))

    @retry(stop=stop_after_attempt(3),
           wait=wait_exponential(multiplier=1, min=17, max=1024),
           retry=(retry_if_exception_type(requests.exceptions.RequestException) |
                  retry_if_exception_type(ConnectionError)))

    # def _new_scrape(self, filename="eaConf.json"):
    #     # scrape_new()
    #     uri_ipac_base = "https://exofop.ipac.caltech.edu/tess/"
    #     uri_ipac_query = {
    #         "select": "TIC ID, TIC, TOI, Telescope, Camera, Filter, Pix Scale (arcsec), PSF (arcsec),"
    #                   "Phot Aperture Rad (pix), Date (UT), Duration (m), # of Obs, Type, Transit Coverage,"
    #                   "Delta Mag, User, Group, Tag, Notes",
    #         "from": "ps",
    #         "where": "tran_flag = 1 and default_flag = 1",
    #         "order by": "pl_pubdate desc",
    #         "format": "csv"
    #     }
    #
    #     if not os.path.exists('pl_names.json') or time.time() - os.path.getmtime('pl_names.json') > 2592000:
    #         self.planet_names(filename="pl_names.json")
    #     if os.path.exists('pl_names.json'):
    #         with open("pl_names.json", "r") as f:
    #             planets = json.load(f)
    #             for key, value in planets.items():
    #                 if self.planet.lower().replace(' ', '').replace('-', '') == key:
    #                     self.planet = value
    #                     break
    #     print(f"\nLooking up {self.planet} on the NASA Exoplanet Archive. Please wait....")
    #
    #     if self.planet:
    #         uri_ipac_query["where"] += f" and pl_name = '{self.planet}'"
    #
    #     default = self._tap_query(uri_ipac_base, uri_ipac_query)
    #
    #     # fill in missing columns
    #     uri_ipac_query['where'] = 'tran_flag=1'
    #
    #     if self.planet:
    #         uri_ipac_query["where"] += f" and pl_name = '{self.planet}'"
    #
    #     extra = self._tap_query(uri_ipac_base, uri_ipac_query)
    #
    #     if len(default) == 0:
    #         self.planet = input(f"Cannot find target ({self.planet}) in NASA Exoplanet Archive."
    #                             f"\nPlease go to https://exoplanetarchive.ipac.caltech.edu to check naming and"
    #                             "\nre-enter the planet's name or type 'candidate' if this is a planet candidate: ")
    #         if self.planet.strip().lower() == 'candidate':
    #             self.planet = user_input("\nPlease enter candidate planet's name: ", type_=str)
    #             return self.planet, True
    #         else:
    #             return self._new_scrape(filename="eaConf.json")
    #     else:
    #         # replaces NEA default with most recent publication
    #         default.iloc[0] = extra.iloc[0]
    #
    #         # for each planet
    #         for i in default.pl_name:
    #
    #             # extract rows for each planet
    #             ddata = default.loc[default.pl_name == i]
    #             edata = extra.loc[extra.pl_name == i]
    #
    #             # for each nan column in default
    #             nans = ddata.isna()
    #             for k in ddata.keys():
    #                 if nans[k].bool():  # if col value is nan
    #                     if not edata[k].isna().all():  # if replacement data exists
    #                         # replace with first index
    #                         default.loc[default.pl_name == i, k] = edata[k][edata[k].notna()].values[0]
    #                         # TODO could use mean for some variables (not mid-transit)
    #                         # log_info(i,k,edata[k][edata[k].notna()].values[0])
    #                     else:
    #                         # permanent nans - require manual entry
    #                         if k == 'pl_orblper':  # omega
    #                             default.loc[default.pl_name == i, k] = 0
    #                         elif k == 'pl_ratdor':  # a/R*
    #                             # Kepler's 3rd law
    #                             semi = SA(ddata.st_mass.values[0], ddata.pl_orbper.values[0])
    #                             default.loc[default.pl_name == i, k] = semi * AU / (
    #                                     ddata.st_rad.values[0] * R_SUN)
    #                         elif k == 'pl_orbincl':  # inclination
    #                             default.loc[default.pl_name == i, k] = 90
    #                         elif k == "pl_orbeccen":  # eccentricity
    #                             default.loc[default.pl_name == i, k] = 0
    #                         elif k == "st_met":  # [Fe/H]
    #                             default.loc[default.pl_name == i, k] = 0
    #                         else:
    #                             default.loc[default.pl_name == i, k] = 0
    #
    #         ExoFOPArchive.dataframe_to_jsonfile(default, filename)
    #         return self.planet, False


#get tbl files for a TOI
def file_get(TIC):
    url = "https://exofop.ipac.caltech.edu/tess/target.php?id="+str(TIC)
    page = requests.get(url)
    data = page.text
    soup = BeautifulSoup(data, 'html.parser')
    x=[[link.get('href'), link.get_text()] for link in soup.find_all('a')]
    df = pd.DataFrame(x, columns = ['url', 'text'])
    df=df.dropna(how='any')
    df=df[df['url'].str.contains("get_file")]
    return df

def file_get_new(TIC):
    TIC_url="https://exofop.ipac.caltech.edu/tess/download_filelist.php?id="+str(TIC)+ '&output=csv'
    # 'https://exofop.ipac.caltech.edu/tess/download_filelist.php?id=' + str(TIC) + '&output=csv'
    test=requests.get(TIC_url)
    assert test.status_code == 200, "download file list failed"
    if TIC == 13883872:
        print(test.text)
    return pd.read_csv(StringIO(test.text))#, sep='|')


def file_download(file_name, file_ID, download_counter):
    print("Downloading ", file_name)
    url="https://exofop.ipac.caltech.edu/tess/get_file.php?id="+str(file_ID)
    r=requests.get(url, allow_redirects=True)
    assert r.status_code == 200, "file download failed"
    destination_file=file_dir+"/"+file_name
    open(destination_file, 'wb').write(r.content)
    return download_counter+1


def bulk_download(obs_type, file_dir=".", TIC, file_search_dict={}, file_ext_dict={}, verbose=False, exotic=False):
    for item in file_search_dict:
        file_search_dict[item]=file_search_dict[item].split('|')
    for item in file_ext_dict:
        file_ext_dict[item]=file_ext_dict[item].split('|')
    print(file_search_dict)
    print(file_ext_dict)
    existing_files=listdir(file_dir)
    page_term=None
    if obs_type.lower()[0]=="i": page_term="imaging"
    if obs_type.lower()[0]=="s": page_term="spect"
    if obs_type.lower()[0]=="t": page_term="tseries"
    if not page_term:
        print("Not a valid category of observations")
        quit()
    source_file="https://exofop.ipac.caltech.edu/tess/download_"+page_term+".php?sort=id&output=pipe"
    print(source_file)
    if exotic == True:
        observed_TIC_list = get_tic(target_name)
    else:
        obs_df=pd.read_csv(source_file, sep='|')
        observed_TOIs=obs_df['TIC ID'].values
        observed_TIC_list=list(dict.fromkeys(observed_TOIs))
    if verbose: print(observed_TIC_list, len(observed_TIC_list))
    print("%.0f TOIs observed" % len(observed_TIC_list))
    already_downloaded=0
    downloads=0
    download_list=[]
    test_df = pd.DataFrame()
    for TIC in observed_TIC_list:
        try:
            test_df = file_get_new(TIC)
            print(TIC)
            for item in file_search_dict:
                data_clean = []
                for value in file_search_dict[item]:
                    condition = test_df[item].str.contains(value)
                    data_clean.append(test_df[test_df[item].str.contains(value)])
                test_df = pd.concat(data_clean, sort=False)
            for item in file_ext_dict:
                data_clean = []
                for value in file_ext_dict[item]:
                    condition = test_df["File Name"].str.endswith(value)
                    data_clean.append(test_df[test_df["File Name"].str.endswith(value)])
                test_df = pd.concat(data_clean, sort=False)

            if verbose: print(test_df)
            # if verbose: print("File extension: ", file_ext)
            # if file_ext:
            #     test_df=test_df[test_df['File Name'].str.endswith(file_ext)]
            # if Type:
            #     test_df = test_df[test_df['Type'] == Type]
            # if verbose: print(test_df)
            for index, row in test_df.iterrows():
                if row['File Name'] in existing_files:
                    print("Have downloaded %s" % row['File Name'])
                    already_downloaded = already_downloaded + 1
                else:
                    print("Don't have %s" % row['File Name'])
                    if verbose: print(row)
                    download_list.append({'file name': row['File Name'], 'file ID': row['File ID']})
        except:
            continue


    print(download_list)
    print("Already downloaded: ", already_downloaded)
    print(len(download_list), " files to download")
    for item in download_list:
        print(item)
        downloads = file_download(item['file name'], item['file ID'], downloads)
    print("%.0f files downloaded" % downloads)
