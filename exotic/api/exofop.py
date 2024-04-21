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

# constants

# CALCULATED VALUES

# class
import requests

import requests
import json

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
            'Planet Name': 'N/A',  # Default if no name is found
            'Host Star Name': self.data.get('basic_info', {}).get('star_names', 'N/A'),
            'Discovery Method': self.data.get('basic_info', {}).get('confirmed_planets', 'N/A')
        }

        # Attempt to extract the planet name from the 'planet_parameters' list
        planet_params = self.data.get('planet_parameters', [])
        for item in planet_params:
            if 'name' in item and item['name']:
                formatted_data['Planet Name'] = item['name']
                break  # Exit after the first valid name is found

        return formatted_data

# misc