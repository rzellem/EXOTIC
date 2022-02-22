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

import requests
from skimage import io
from json import dumps
from tenacity import retry, retry_if_exception_type, retry_if_result, \
    stop_after_attempt, wait_exponential


def is_false(value):
    return value is False


def result_if_max_retry_count(retry_state):
    pass


class vspPlot:

    def __init__(self, host=None, maglimit=14.5, fov=20, dss=False):
        self.host = host
        self.maglimit = maglimit
        self.fov = fov
        self.dss = dss
        self.status = 200
        self.api_url = "http://www.aavso.org/apps/vsp/api/chart/?"

    # Simply takes plot url - not sure if any of this remaining photometry data
    # for stars in the field is useful for our purposes
    def findPlot(self):
        if (self.host is None):
            return None

        plot_data = self._get_json()
        if not plot_data:
            return vspPlot.fail(self.status)

        plot_url = plot_data["image_uri"]
        image = io.imread(plot_url)
        return image

    def _get_url_params(self):
        # "https://app.aavso.org/vsp/api/chart/?star=SS+Cyg&fov=60&maglimit=14.5"
        return self.api_url + "star={}&fov={}&maglimit={}&dss={}&orientation=ccd&format=json" \
            .format(self.host, self.fov, self.maglimit, self.dss)

    def _get_url_charid(self):
        pass

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _get_json(self):
        r = requests.get(self._get_url_params())
        self.status = r.status_code
        if self.status >= 400:
            return False
        else:
            return r.json()

    @staticmethod
    def fail(status):
        # Better ways of handling error cases - leave ambigious if >=400?
        if status == 400:
            error = " due to improper query"
        elif status == 404:
            error = " due to servers not responding"
        elif status == 408:
            error = " due to request timeout"
        else:
            error = ""
        print("WARNING: After multiple attempts, EXOTIC could not retrieve a VSP finder chart from app.aavso.org"
              f"{error}. Please use another method for choosing your comparison stars")
        return None