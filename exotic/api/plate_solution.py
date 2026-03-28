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
from astropy.io.fits import Header, PrimaryHDU, getdata, getheader
from astropy.stats import sigma_clipped_stats
from json import dumps
from pathlib import Path
import numpy as np
from photutils.detection import DAOStarFinder
import requests
import time
from tenacity import retry, retry_if_exception_type, retry_if_result, \
    stop_after_attempt, wait_exponential

_R_MAX_STOPS_LOW = 7
_R_MAX_STOPS = 10
_R_MAX_SECS = 37
_RQ_TIMEOUT = 16.0
_NEXTASTRO_MAX_SOURCES = 200
_NEXTASTRO_STATUS_MAX_POLLS = 60
_NEXTASTRO_STATUS_POLL_SEC = 2
_NEXTASTRO_IN_PROGRESS_STATUSES = {'queued', 'running'}


def is_false(value):
    return value is False


def result_if_max_retry_count(retry_state):
    pass


class PlateSolution:

    def __init__(self, file=None, directory=None, api_key=None,
                 api_url='http://nova.astrometry.net/api/', ra=None, dec=None,
                 pixel_scale=None, radius=2.0, scale_err=25.0, suppress_fail_warning=False):
        if api_key is None:
            api_key = {'apikey': 'vfsyxlmdxfryhprq'}
        self.api_url = api_url
        self.api_key = api_key
        self.file = file
        self.directory = directory
        self.ra = ra
        self.dec = dec
        self.pixel_scale = pixel_scale
        self.radius = radius
        self.scale_err = scale_err
        self.suppress_fail_warning = suppress_fail_warning
        self.last_error_type = None

    def plate_solution(self):
        self.last_error_type = None
        session = self._login()
        if not session:
            return self._fail('Login')

        sub_id = self._upload(session)
        if not sub_id:
            return self._fail('Upload')

        sub_url = self._get_url(f"submissions/{sub_id}")
        job_id = self._sub_status(sub_url)
        if not job_id:
            return self._fail('Submission ID')

        job_url = self._get_url(f"jobs/{job_id}")
        download_url = self.api_url.replace("/api/", f"/wcs_file/{job_id}/")
        wcs_file = Path(self.directory) / "temp" / "wcs.fits"
        wcs_file = self._job_status(job_url, wcs_file, download_url)
        if not wcs_file:
            return self._fail('Job Status')
        else:
            print("WCS file creation successful.")
            return wcs_file

    def _get_url(self, service):
        return self.api_url + service

    def _fail(self, error_type, service_name='nova.astrometry.net'):
        self.last_error_type = error_type
        if self.suppress_fail_warning:
            return False
        return PlateSolution.fail(error_type, service_name=service_name)

    @retry(stop=stop_after_attempt(_R_MAX_STOPS_LOW), wait=wait_exponential(multiplier=1, min=4, max=_R_MAX_SECS),
           retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _login(self):
        r = requests.post(self._get_url('login'), data={'request-json': dumps(self.api_key)}, timeout=_RQ_TIMEOUT)
        if r.status_code >= 400:
            return False
        elif r.json()['status'] == 'success':
            return r.json()['session']
        return False

    @retry(stop=stop_after_attempt(_R_MAX_STOPS_LOW), wait=wait_exponential(multiplier=1, min=4, max=_R_MAX_SECS),
           retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _upload(self, session):
        request_payload = {"session": session}

        if self.ra is not None and self.dec is not None:
            request_payload.update({
                "center_ra": float(self.ra),
                "center_dec": float(self.dec),
                "radius": float(self.radius)
            })

        if self.pixel_scale not in (None, ""):
            request_payload.update({
                "scale_units": "arcsecperpix",
                "scale_type": "ev",
                "scale_est": float(self.pixel_scale),
                "scale_err": float(self.scale_err)
            })

        headers = {'request-json': dumps(request_payload), 'allow_commercial_use': 'n',
                   'allow_modifications': 'n', 'publicly_visible': 'n'}

        with open(self.file, 'rb') as image_file:
            files = {'file': image_file}
            r = requests.post(self.api_url + 'upload', files=files, data=headers, timeout=_RQ_TIMEOUT)

        if r.json()['status'] == 'success':
            return r.json()['subid']
        return False

    @retry(stop=stop_after_attempt(_R_MAX_STOPS), wait=wait_exponential(multiplier=1, min=4, max=_R_MAX_SECS),
           retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _sub_status(self, sub_url):
        r = requests.get(sub_url, timeout=_RQ_TIMEOUT)
        if r.json()['job_calibrations']:
            return r.json()['jobs'][0]
        return False

    @retry(stop=stop_after_attempt(_R_MAX_STOPS), wait=wait_exponential(multiplier=1, min=4, max=_R_MAX_SECS),
           retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
           retry_error_callback=result_if_max_retry_count)
    def _job_status(self, job_url, wcs_file, download_url):
        r = requests.get(job_url, timeout=_RQ_TIMEOUT)
        if r.json()['status'] == 'success':
            r = requests.get(download_url, timeout=_RQ_TIMEOUT)
            with wcs_file.open('wb') as f:
                f.write(r.content)
            hdu = PrimaryHDU(data=getdata(filename=self.file), header=getheader(filename=wcs_file))
            hdu.writeto(wcs_file, overwrite=True)
            return wcs_file
        return False

    @staticmethod
    def fail(error_type, service_name='nova.astrometry.net'):
        print("WARNING: After multiple attempts, EXOTIC could not retrieve a plate solution from "
              f"{service_name} due to {error_type}. EXOTIC will continue reducing data without a plate solution.")
        return False


class NextAstroPlateSolution:

    def __init__(self, file=None, directory=None, api_url='https://astrometry.nextastro.org/', ra=None, dec=None,
                 pixel_scale=None, suppress_fail_warning=False):
        self.api_url = api_url.rstrip('/')
        self.file = file
        self.directory = directory
        self.ra = ra
        self.dec = dec
        self.pixel_scale = pixel_scale
        self.suppress_fail_warning = suppress_fail_warning
        self.last_error_type = None
        self.last_http_status = None

    def plate_solution(self):
        self.last_error_type = None
        self.last_http_status = None
        self._emit_debug(f"Using NextAstro astrometry server at {self.api_url} for plate solving.")
        source_list = self._generate_source_list()
        if not source_list:
            return self._fail('Source extraction for NextAstro astrometry server')

        request_id = self._submit_solve_request(source_list)
        if not request_id:
            return self._fail('NextAstro solve submission')

        wcs_header = self._poll_for_solution(request_id)
        if not wcs_header:
            return self._fail('NextAstro solve status')

        wcs_file = Path(self.directory) / "temp" / "wcs.fits"
        hdu = PrimaryHDU(data=getdata(filename=self.file), header=wcs_header)
        hdu.writeto(wcs_file, overwrite=True)
        self._emit_debug("WCS file creation successful.")
        return wcs_file

    def _emit_debug(self, message):
        if not self.suppress_fail_warning:
            print(message)

    def _fail(self, error_type):
        self.last_error_type = error_type
        if self.suppress_fail_warning:
            return False
        return PlateSolution.fail(error_type, service_name=f'NextAstro ({self.api_url})')

    def _generate_source_list(self):
        image_data = np.asarray(getdata(filename=self.file), dtype=float)
        if image_data.ndim > 2:
            image_data = image_data.squeeze()

        median, _, std = sigma_clipped_stats(image_data, sigma=3.0)
        if std <= 0:
            std = float(np.nanstd(image_data))
            if std <= 0:
                return None

        finder = DAOStarFinder(fwhm=3.0, threshold=3.5 * std)
        sources = finder(image_data - median)

        if sources is not None and len(sources) > 0:
            bright_sources = self._limit_to_brightest_sources(
                x_coords=sources['xcentroid'],
                y_coords=sources['ycentroid'],
                fluxes=sources['flux']
            )
            if bright_sources is None:
                return None
            return {
                "x": bright_sources["x"],
                "y": bright_sources["y"],
                "flux": bright_sources["flux"],
                "pixel_indexing": "0-based"
            }

        return self._fallback_source_list(image_data, median, std)


    def _fallback_source_list(self, image_data, median, std):
        threshold = median + 3.5 * std
        candidate_indices = np.argwhere(image_data > threshold)
        if candidate_indices.size == 0:
            return None

        candidate_fluxes = image_data[candidate_indices[:, 0], candidate_indices[:, 1]]
        bright_sources = self._limit_to_brightest_sources(
            x_coords=candidate_indices[:, 1],
            y_coords=candidate_indices[:, 0],
            fluxes=candidate_fluxes
        )
        if bright_sources is None:
            return None

        return {
            "x": bright_sources["x"],
            "y": bright_sources["y"],
            "flux": bright_sources["flux"],
            "pixel_indexing": "0-based"
        }

    @staticmethod
    def _limit_to_brightest_sources(x_coords, y_coords, fluxes):
        fluxes = np.asarray(fluxes, dtype=float)
        x_coords = np.asarray(x_coords, dtype=float)
        y_coords = np.asarray(y_coords, dtype=float)

        finite_flux_mask = np.isfinite(fluxes)
        if not np.any(finite_flux_mask):
            return None

        fluxes = fluxes[finite_flux_mask]
        x_coords = x_coords[finite_flux_mask]
        y_coords = y_coords[finite_flux_mask]

        sorted_indices = np.argsort(fluxes)[::-1][:_NEXTASTRO_MAX_SOURCES]

        return {
            "x": x_coords[sorted_indices].tolist(),
            "y": y_coords[sorted_indices].tolist(),
            "flux": fluxes[sorted_indices].tolist()
        }

    @staticmethod
    def _response_body_preview(response, max_chars=240):
        body = getattr(response, 'text', None)
        if body is None:
            content = getattr(response, 'content', b'')
            body = content.decode(errors='replace') if isinstance(content, bytes) else str(content)

        body = " ".join(str(body).split())
        if not body:
            return "<empty response body>"
        if len(body) > max_chars:
            return body[:max_chars - 3] + "..."
        return body

    def _decode_response_json(self, response, context):
        self.last_http_status = getattr(response, 'status_code', None)
        try:
            return response.json()
        except ValueError:
            if response.status_code != 502:
                self._emit_debug(f"[NextAstro] {context} returned non-JSON response "
                                 f"(HTTP {response.status_code}): {self._response_body_preview(response)}")
            return None

    def _submit_solve_request(self, source_list):
        image_data = getdata(filename=self.file)
        payload = {
            "sources": source_list,
            "image": {
                "width": int(image_data.shape[-1]),
                "height": int(image_data.shape[-2]),
            },
            "options": {
                "timeout_sec": 120,
                "max_sources": _NEXTASTRO_MAX_SOURCES
            }
        }

        hints = self._extract_astrometry_hints()
        if hints is not None:
            payload["hints"] = hints

        self._emit_debug(f"[NextAstro] Solve request payload: {payload}")
        response = requests.post(f"{self.api_url}/solve", json=payload, timeout=_RQ_TIMEOUT)
        response_json = self._decode_response_json(response, 'Solve response')
        if response_json is not None and response.status_code != 502:
            self._emit_debug(f"[NextAstro] Solve response: {response_json}")
        if response.status_code >= 400 or response_json is None:
            return False
        if response_json.get('status') in {'queued', 'running'}:
            return response_json.get('request_id')
        return False

    def _extract_astrometry_hints(self):
        hints = {}

        if self.ra is not None and self.dec is not None:
            hints.update({"ra_deg": self.ra, "dec_deg": self.dec})

        if self.pixel_scale not in (None, ""):
            hints.update({"scale_arcsec_per_pix": float(self.pixel_scale), "scale_tolerance_frac": 0.25})

        if not hints:
            return None

        return hints

    def _poll_for_solution(self, request_id):
        latest_status = None
        for _ in range(_NEXTASTRO_STATUS_MAX_POLLS):
            response = requests.get(f"{self.api_url}/status/{request_id}", timeout=_RQ_TIMEOUT)
            response_json = self._decode_response_json(response, 'Status response')
            if response_json is None:
                return False
            if response.status_code >= 400:
                if response.status_code != 502:
                    self._emit_debug(f"[NextAstro] Status response (HTTP {response.status_code}): {response_json}")
                return False

            status = str(response_json.get('status', '')).lower()
            latest_status = response_json.get('status')
            if status == 'solved':
                self._emit_debug(f"[NextAstro] Status response (solved): {response_json}")
                header_dict = response_json.get('solution', {}).get('wcs_header')
                if isinstance(header_dict, dict):
                    return Header(header_dict)
                return False
            if status == 'failed':
                self._emit_debug(f"[NextAstro] Status response (failed): {response_json}")
                return False

            if status not in _NEXTASTRO_IN_PROGRESS_STATUSES:
                self._emit_debug(f"[NextAstro] Status response (unexpected): {response_json}")
                return False

            time.sleep(_NEXTASTRO_STATUS_POLL_SEC)

        self._emit_debug(f"[NextAstro] Polling timed out waiting for terminal status; latest status={latest_status!r}")
        return False
