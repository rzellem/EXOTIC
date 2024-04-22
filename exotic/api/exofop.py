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

# imports
import requests
import json
import numpy as np

# constants

# CALCULATED VALUES

# class

class ExoFOP:
    def __init__(self, tic_code=None):
        self.tic_code = tic_code
        self.base_url = "https://exofop.ipac.caltech.edu/tess/target.php"
        self.data = None

    def query_exofop(self):
        """
        Queries the ExoFOP database for the TIC code and parses the JSON data.
        """
        params = {'id': self.tic_code, 'json': ''}
        response = requests.get(self.base_url, params=params)
        if response.status_code == 200:
            self.data = response.json()
            return self.data
        else:
            raise Exception(f"Failed to retrieve data: Status code {response.status_code}")

    def get_formatted_data(self):
        """
        Formats the JSON data into a dictionary and tries to extract the planet name from 'planet_parameters'.
        """
        if self.data is None:
            print("Data not loaded. Please run 'query_exofop()' first.")
            return None

        formatted_data = {
            'TIC ID': self.tic_code,
            'ra': 'N/A',
            'dec': 'N/A',
            'pName': 'N/A',
            'sName': 'N/A',
            'pPer': 'N/A/',
            'pPerUnc': 'N/A',
            'midT': 'N/A',
            'rprs': 'N/A',
            'rprsUnc': 'N/A',
            'aRs': 'N/A', 
            'aRsUnc': 'N/A',
            'inc': 'N/A',
            'incUnc': 'N/A',
            'omega': 0,
            'ecc': 0,
            'TIC ID': self.tic_code,
            'Discovery Method': self.data.get('basic_info', {}).get('confirmed_planets', 'N/A')
        }

        # Update RA if coordinates are available
        if 'coordinates' in self.data:
            formatted_data['ra'] = self.data['coordinates'].get('ra', 'N/A')

        # Update DEC if coordinates are available
        if 'coordinates' in self.data:
            formatted_data['dec'] = self.data['coordinates'].get('dec', 'N/A')

        planet_params = self.data.get('planet_parameters', [])
        for item in planet_params:
            if 'name' in item and item['name']:
                formatted_data['pName'] = item['name']
            if 'per' in item and item['per']:
                formatted_data['pPer'] = item['per']
                formatted_data['omega'] = item['per']
            if 'per_e' in item and item['per_e']:
                pl_orbpererr1 = float(item['per_e'])
                pl_orbpererr2 = float(item['per_e'])
                pPerUnc = np.sqrt(np.abs(pl_orbpererr1 * (pl_orbpererr2 * -1))) #  or just... abs(pl_orbpererr1)  # Since we're squaring it, sqrt(abs(x^2)) == abs(x)
                formatted_data['pPerUnc'] = pPerUnc
            if 'epoch' in item and item['epoch']:
                formatted_data['midT'] = item['epoch']
            if 'epoch_e' in item and item['epoch_e']:
                pl_tranmiderr1 = float(item['epoch_e'])
                pl_tranmiderr2 = float(item['epoch_e'])
                midTUnc = np.sqrt(np.abs(pl_tranmiderr1 * (pl_tranmiderr2 * -1))) #  or just... abs(pl_tranmiderr1)  # Since we're squaring it, sqrt(abs(x^2)) == abs(x)
                formatted_data['miDTUnc'] = midTUnc
            if 'rr' in item and item['rr']:
                formatted_data['rprs'] = item['rr']
            if 'rr_e' in item and item['rr_e']:
                formatted_data['rprsUnc'] = item['rr_e']
            if 'ar' in item and item['ar']:
                formatted_data['aRs'] = item['ar']
            if 'ar_e' in item and item['ar_e']:
                if item['ar_e'] == 0:
                    formatted_data['aRsUnc'] = 0
                else:
                    pl_ratdorerr1 = float(item['ar_e'])
                    pl_ratdorerr2 = float(item['ar_e'])
                    aRsUnc = np.sqrt(np.abs(pl_ratdorerr1 * (pl_ratdorerr2 * -1))) #  or just... abs(pl_ratdorerr1)  # Since we're squaring it, sqrt(abs(x^2)) == abs(x)
                    formatted_data['aRsUnc'] = aRsUnc
            if 'inc' in item and item['inc']:
                formatted_data['inc'] = item['inc']
            if 'inc_e' in item and item['inc_e']:
                if item['inc_e'] == 0:
                    formatted_data['incUnc'] = 0.1
                else:
                    formatted_data['incUnc'] = item['inc_e']
            if 'ecc' in item and item['ecc']:
                formatted_data['ecc'] = item['ecc']
                break

        # Extract the first star name if available
        star_names = self.data.get('basic_info', {}).get('star_names', '')
        if star_names:
            formatted_data['sName'] = star_names.split(',')[0].strip()

        return formatted_data

# misc

# Fields requesed in the NEA ipac query of NASAExoplanetArchive class (nea.py)

# General Parameters

#     pl_name: Planet Name - The designated name of the exoplanet.
#     hostname: Host Star Name - The name of the star around which the exoplanet orbits.
#     tran_flag: Transit Flag - Indicates whether the planet transits its host star (usually a boolean value).
#     ra: Right Ascension - The right ascension of the star in celestial coordinates, analogous to longitude on Earth.
#     dec: Declination - The declination of the star in celestial coordinates, analogous to latitude on Earth.
#     pl_pubdate: Publication Date - The date when the data was published or last updated.

# Planet-Specific Parameters

#     pl_massj: Planet Mass [Jupiter mass] - The mass of the planet in units of Jupiter masses.
#     pl_radj: Planet Radius [Jupiter radius] - The radius of the planet in units of Jupiter radii.
#     pl_radjerr1, pl_radjerr2: Error in Planet Radius - Positive and negative errors for the planet radius.
#     pl_ratdor: Planet-Star Distance over Star Radius (a/R_star) - This ratio helps in understanding the scale of the planetary system.
#     pl_ratdorerr1, pl_ratdorerr2: Error in Planet-Star Distance over Star Radius.
#     pl_orbincl: Orbital Inclination [degrees] - The angle of the planet's orbit relative to the line of sight from Earth.
#     pl_orbinclerr1, pl_orbinclerr2: Error in Orbital Inclination.
#     pl_orbper: Orbital Period [days] - The time it takes for the planet to complete one orbit around its star.
#     pl_orbpererr1, pl_orbpererr2: Error in Orbital Period.
#     pl_orbeccen: Orbital Eccentricity - A measure of the deviation of the orbit from circularity.
#     pl_orblper: Longitude of Periastron [degrees] - The angular distance where the planet is closest to its star.
#     pl_tranmid: Transit Midpoint [BJD] - The mid-point of when the planet transits across the star.
#     pl_tranmiderr1, pl_tranmiderr2: Error in Transit Midpoint.
#     pl_trandep: Transit Depth [parts per million] - The amount by which the star's light dims during a transit.
#     pl_trandeperr1, pl_trandeperr2: Error in Transit Depth.
#     pl_ratror: Planet to Star Radius Ratio - Ratio of the planet's radius to the star's radius.
#     pl_ratrorerr1, pl_ratrorerr2: Error in Radius Ratio.

# Star-Specific Parameters

#     st_teff: Stellar Effective Temperature [Kelvin] - The surface temperature of the star.
#     st_tefferr1, st_tefferr2: Error in Stellar Effective Temperature.
#     st_met: Stellar Metallicity [dex] - The logarithmic measure of the star's metal content relative to the Sun.
#     st_meterr1, st_meterr2: Error in Stellar Metallicity.
#     st_logg: Stellar Surface Gravity [log10(cm/s^2)] - The logarithm of the star's surface gravity.
#     st_loggerr1, st_loggerr2: Error in Surface Gravity.
#     st_mass: Stellar Mass [Solar mass] - The mass of the star in units of solar masses.
#     st_rad: Stellar Radius [Solar radius] - The radius of the star in units of solar radii.
#     st_raderr1, st_raderr2: Error in Stellar Radius.

# Below is how each item is added to the pl_dict in the NASAExoplanetArchive class (nea.py)
# self.pl_dict = {
# (DONE/)     'ra': data['ra'],
# (DONE/)     'dec': data['dec'],
# (DONE/)     'pName': data['pl_name'],
# (DONE/)     'sName': data['hostname'],
# (DONE/)     'pPer': data['pl_orbper'],
# (DONE/)     'pPerUnc': np.sqrt(np.abs(data['pl_orbpererr1'] * data['pl_orbpererr2'])),

# (DONE/)     'midT': data['pl_tranmid'],
# (DONE/)     'midTUnc': np.sqrt(np.abs(data['pl_tranmiderr1'] * data['pl_tranmiderr2'])),
# (DONE/)     'rprs': rprs,
# (DONE/)     'rprsUnc': rprserr,
# (DONE/)     'aRs': data['pl_ratdor'],
# (DONE/)     'aRsUnc': np.sqrt(np.abs(data.get('pl_ratdorerr1', 1) * data['pl_ratdorerr2'])),
# (DONE/)     'inc': data['pl_orbincl'],
# (DONE/)     'incUnc': np.sqrt(np.abs(data['pl_orbinclerr1'] * data['pl_orbinclerr2'])),
# (DONE/)     'omega': data.get('pl_orblper', 0),
#             'ecc': data.get('pl_orbeccen', 0),
#             'teff': data['st_teff'],
#             'teffUncPos': data['st_tefferr1'],
#             'teffUncNeg': data['st_tefferr2'],
#             'met': data['st_met'],
#             'metUncPos': max(0.01, data['st_meterr1']),
#             'metUncNeg': min(-0.01, data['st_meterr2']),
#             'logg': data['st_logg'],
#             'loggUncPos': data['st_loggerr1'],
#             'loggUncNeg': data['st_loggerr2']
#         }

# (DONE/) if self.pl_dict['aRsUnc'] == 0:
#             self.pl_dict['aRsUnc'] = 0.1

# (DONE/) if self.pl_dict['incUnc'] == 0:
#             self.pl_dict['incUnc'] = 0.1