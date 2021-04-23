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
from json import dumps
from pathlib import Path
from astropy.io.fits import PrimaryHDU, getdata, getheader
from tenacity import retry, retry_if_exception_type, retry_if_result, \
    stop_after_attempt, wait_exponential


def is_false(value):
    return value is False


def result_if_max_retry_count(retry_state):
    pass


class PlateSolution:

    def __init__(self, file=None, directory=None, api_key=None,
                 api_url='http://nova.astrometry.net/api/'):
        if api_key is None:
            api_key = {'apikey': 'vfsyxlmdxfryhprq'}
        self.api_url = api_url
        self.api_key = api_key
        self.file = file
        self.directory = directory

    def plate_solution(self):
        session = self._login()
        if not session:
            return PlateSolution.fail('Login')

        sub_id = self._upload(session)
        if not sub_id:
            return PlateSolution.fail('Upload')

        sub_url = self._get_url(f"submissions/{sub_id}")
        job_id = self._sub_status(sub_url)
        if not job_id:
            return PlateSolution.fail('Submission ID')

        job_url = self._get_url(f"jobs/{job_id}")
        download_url = self.api_url.replace("/api/", f"/wcs_file/{job_id}/")
        wcs_file = Path(self.directory) / "wcs.fits"
        wcs_file = self._job_status(job_url, wcs_file, download_url)
        if not wcs_file:
            return PlateSolution.fail('Job Status')
        else:
            print("WCS file creation successful.")
            return wcs_file

    def _get_url(self, service):
        return self.api_url + service

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _login(self):
        r = requests.post(self._get_url('login'), data={'request-json': dumps(self.api_key)})
        if r.status_code >= 400:
            return False
        elif r.json()['status'] == 'success':
            return r.json()['session']
        return False

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _upload(self, session):
        files = {'file': open(self.file, 'rb')}
        headers = {'request-json': dumps({"session": session}), 'allow_commercial_use': 'n',
                   'allow_modifications': 'n', 'publicly_visible': 'n'}

        r = requests.post(self.api_url + 'upload', files=files, data=headers)

        if r.json()['status'] == 'success':
            return r.json()['subid']
        return False

    @retry(stop=stop_after_attempt(20), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _sub_status(self, sub_url):
        r = requests.get(sub_url)
        if r.json()['job_calibrations']:
            return r.json()['jobs'][0]
        return False

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _job_status(self, job_url, wcs_file, download_url):
        r = requests.get(job_url)
        if r.json()['status'] == 'success':
            r = requests.get(download_url)
            with wcs_file.open('wb') as f:
                f.write(r.content)
            hdu = PrimaryHDU(data=getdata(filename=self.file), header=getheader(filename=wcs_file))
            hdu.writeto(wcs_file, overwrite=True)
            return wcs_file
        return False

    @staticmethod
    def fail(error_type):
        print("WARNING: After multiple attempts, EXOTIC could not retrieve a plate solution from nova.astrometry.net"
              f" due to {error_type}. EXOTIC will continue reducing data without a plate solution.")
        return False
