import logging
import requests
from json import dumps
from pathlib import Path
from tenacity import retry, retry_if_exception_type, retry_if_result, \
    stop_after_attempt, wait_exponential

log = logging.getLogger(__name__)


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

        sub_url = self._get_url(f"submissions/{sub_id}")
        job_id = self._sub_status(sub_url)
        if not job_id:
            return PlateSolution.fail('Submission ID')

        job_url = self._get_url(f"jobs/{job_id}")
        download_url = self.apiurl.replace("/api/", f"/new_fits_file/{job_id}/")
        wcs_file = Path(self.directory) / "wcs_image.fits"
        wcs_file = self._job_status(job_url, wcs_file, download_url)
        if not wcs_file:
            return PlateSolution.fail('Job Status')
        else:
            log.info("WCS file creation successful.")
            return wcs_file

    def _get_url(self, service):
        return self.apiurl + service

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(ConnectionError) |
                  retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _login(self):
        r = requests.post(self._get_url('login'), data={'request-json': dumps(self.apikey)})
        if r.json()['status'] == 'success':
            return r.json()['session']
        return False

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
           retry=(retry_if_result(is_false) | retry_if_exception_type(ConnectionError) |
                  retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _upload(self, session):
        files = {'file': open(self.file, 'rb')}
        headers = {'request-json': dumps({"session": session}), 'allow_commercial_use': 'n',
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
        log.info("WARNING: After multiple attempts, EXOTIC could not retrieve a plate solution from nova.astrometry.net"
                 f" due to {error_type}. EXOTIC will continue reducing data without a plate solution.")
        return False
