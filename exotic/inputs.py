import logging
import sys
import json
import math
from pathlib import Path
import requests
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import re

try:
    from utils import user_input, init_params, typecast_check, \
        process_lat_long, find, open_elevation
except ImportError:
    from .utils import user_input, init_params, typecast_check, \
        process_lat_long, find, open_elevation
try:
    from animate import animate_toggle
except ImportError:
    from .animate import animate_toggle
try:
    from api.filters import fwhm as photometric_filters, fwhm_alias as photometric_filter_aliases
except ImportError:
    from .api.filters import fwhm as photometric_filters, fwhm_alias as photometric_filter_aliases


log = logging.getLogger(__name__)
consoleFormatter = logging.Formatter("%(message)s")
consoleHandler = logging.StreamHandler(sys.stdout)
consoleHandler.setFormatter(consoleFormatter)
consoleHandler.setLevel(logging.INFO)
log.addHandler(consoleHandler)

PHOT_COMP_STAR_KEYS = ("ra", "dec", "x", "y")
AAVSO_OBSDATE_HEADER_KEYS = ('OBSDATE',)
AAVSO_LOCATION_HEADER_KEYS = {
    'lat': ('OBSLAT', 'LATITUDE', 'OBS_LATITUDE', 'LAT'),
    'long': ('OBSLON', 'OBSLONG', 'LONGITUDE', 'OBS_LONGITUDE', 'LONG'),
    'elev': ('OBSELEV', 'OBSALT', 'ELEVATION', 'ALTITUDE', 'HEIGHT'),
}
AAVSO_FILTER_HEADER_KEYS = ('FILTER',)
AAVSO_FILTER_XC_HEADER_KEYS = ('FILTER-XC',)
AAVSO_TEXT_HEADER_KEYS = {
    'aavso_num': ('OBSCODE',),
    'second_obs': ('SECONDARY_OBSCODES',),
    'obs_name': ('OBSNAME',),
    'camera': ('OBSTYPE',),
    'pixel_bin': ('BINNING',),
    'notes': ('NOTES',),
    'planet': ('EXOPLANET_NAME',),
    'host_star': ('STAR_NAME',),
}
AAVSO_GAIA_HEADER_KEYS = {
    'dist': ('GAIADIST',),
    'pm_ra': ('GAIAPMRA',),
    'pm_dec': ('GAIAPMDEC', 'GAIADEC'),
}
AAVSO_EXPOSURE_HEADER_KEYS = ('EXPOSURE_TIME', 'EXPTIME', 'EXPOSURE', 'EXP')
AAVSO_TIME_FORMAT_HEADER_KEYS = ('DATE_TYPE',)
AAVSO_MEASUREMENT_TYPE_HEADER_KEYS = ('MEASUREMENT_TYPE',)
AAVSO_DETREND_PARAMETER_HEADER_KEYS = ('DETREND_PARAMETERS',)
AAVSO_ALLOWED_FILE_TIME_FORMATS = {'BJD_TDB', 'JD_UTC', 'MJD_UTC'}
AAVSO_WAVELENGTH_UNIT_FACTORS_TO_NM = {
    'a': 0.1,
    'angstrom': 0.1,
    'angstroms': 0.1,
    'nm': 1.0,
    'nanometer': 1.0,
    'nanometers': 1.0,
    'um': 1000.0,
    'micron': 1000.0,
    'microns': 1000.0,
    'micrometer': 1000.0,
    'micrometers': 1000.0,
    'mum': 1000.0,
}
NEXTASTRO_GAIA_DISTPM_ENDPOINT = 'https://archive.nextastro.org/single_star_gaia_distpm'
NEXTASTRO_REQUEST_TIMEOUT = 30


def is_blank_value(value):
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip().lower() in ('', 'n/a', 'na', 'null', 'none')
    return False


def normalize_aavso_filter_lookup_key(value):
    if is_blank_value(value):
        return None
    return re.sub(r'[\W_]+', '', str(value).strip().lower())


def coerce_finite_float(value):
    if is_blank_value(value):
        return None

    try:
        numeric_value = float(str(value).strip())
    except (TypeError, ValueError):
        return None

    if not math.isfinite(numeric_value):
        return None

    return numeric_value


def radec_to_decimal_degrees(ra, dec):
    if is_blank_value(ra) or is_blank_value(dec):
        return None, None

    ra_value = str(ra).strip()
    dec_value = str(dec).strip()
    ra_unit = u.hourangle if any(separator in ra_value for separator in (':', ' ')) else u.deg

    if ra_unit is u.hourangle:
        ra_value = ra_value.replace(':', ' ')
    if any(separator in dec_value for separator in (':', ' ')):
        dec_value = dec_value.replace(':', ' ')

    try:
        coords = SkyCoord(ra=ra_value, dec=dec_value, unit=(ra_unit, u.deg))
    except ValueError:
        return None, None

    if not math.isfinite(coords.ra.degree) or not math.isfinite(coords.dec.degree):
        return None, None

    return coords.ra.degree, coords.dec.degree


def fetch_nextastro_gaia_distpm(ra_deg, dec_deg):
    response = requests.get(
        NEXTASTRO_GAIA_DISTPM_ENDPOINT,
        params={'ra': ra_deg, 'dec': dec_deg},
        timeout=NEXTASTRO_REQUEST_TIMEOUT,
    )
    response.raise_for_status()

    payload = response.json()
    gaia = payload.get('gaia') if isinstance(payload, dict) else None
    if not isinstance(gaia, dict):
        return {}

    return {
        'dist': coerce_finite_float(gaia.get('distance_pc')),
        'pm_ra': coerce_finite_float(gaia.get('pmra_mas_per_year')),
        'pm_dec': coerce_finite_float(gaia.get('pmdec_mas_per_year')),
    }


def populate_missing_gaia_astrometry(planet_dict):
    missing_keys = [key for key in ('dist', 'pm_ra', 'pm_dec') if is_blank_value(planet_dict.get(key))]
    if not missing_keys:
        return planet_dict

    ra_deg, dec_deg = radec_to_decimal_degrees(planet_dict.get('ra'), planet_dict.get('dec'))
    if ra_deg is None or dec_deg is None:
        return planet_dict

    try:
        gaia_values = fetch_nextastro_gaia_distpm(ra_deg, dec_deg)
    except requests.exceptions.RequestException as exc:
        log_info(f"\nWarning: NextAstro Gaia astrometry lookup failed ({exc}); continuing without missing Gaia values.",
                 warn=True)
        return planet_dict

    filled_keys = []
    for key in missing_keys:
        if gaia_values.get(key) is None:
            continue
        planet_dict[key] = gaia_values[key]
        filled_keys.append(key)

    if filled_keys:
        log_info("\nRetrieved missing Gaia distance/proper motion from NextAstro archive lookup.")

    return planet_dict


AAVSO_FILTER_LOOKUP = {}
for filter_desc, filter_metadata in photometric_filters.items():
    AAVSO_FILTER_LOOKUP[normalize_aavso_filter_lookup_key(filter_desc)] = filter_metadata
    AAVSO_FILTER_LOOKUP[normalize_aavso_filter_lookup_key(filter_metadata.get('name'))] = filter_metadata

for alias, canonical in photometric_filter_aliases.items():
    filter_metadata = photometric_filters.get(canonical)
    if filter_metadata is not None:
        AAVSO_FILTER_LOOKUP[normalize_aavso_filter_lookup_key(alias)] = filter_metadata


class Inputs:

    def __init__(self, init_opt):
        self.init_opt = init_opt
        self.info_dict = {
            'images': None, 'save': None, 'flats': None, 'darks': None, 'biases': None,
            'aavso_num': None, 'second_obs': None, 'obs_name': '', 'date': None, 'lat': None, 'long': None,
            'elev': None, 'camera': None, 'pixel_bin': None, 'filter': None, 'notes': None,
            'plate_opt': None, 'aavso_comp': None, 'tar_coords': None, 'comp_stars': None,
            'prered_file': None, 'file_units': None, 'file_time': None, 'phot_comp_star': None,
            'wl_min': None, 'wl_max': None, 'pixel_scale': None, 'exposure': None,
            'dist': None, 'pm_ra': None, 'pm_dec': None, 'airmass_already_corrected': False,
            'random_seed': None, 'ld_uncertainties': None, "demosaic_fmt": None, "demosaic_out": None,
            'fast_aperture_mask': False, 'require_comp_star': 'y', 'ignore_header_wcs': 'n',
            'prefer_pixel_values_over_wcs_for_target': 'n',
            'target_driven_comp_selection': 'n', 'disable_vertical_flux_normalization': False,
            'stellar_variability_only': False,
            'use_ensemble_photometry_for_stellar_variability': True,
            'photometer_fortuitous_variables': True,
            'use_nextastro_vsx_cache_first': False,
            'detrend_on_outoftransit_baseline': True,
            'final_fit_baseline_duration_multiplier': 1.0,
            'use_eebls_to_initialize_tmid_and_bounds': 'y',
            'pick_comparison_by_eebls_snr': 'y',
            'use_deviation_from_expected_transit_in_qc': True,
            'deviation_from_expected_transit_in_qc_sigma': 5.0,
            'run_final_fit_phase_residual_clip': 'y',
            'exit_at_first_qc_pass_solution': 'y',
            'detect_bad_pixels_before_photometry': 'n',
            'multiprocess_bad_pixel_precheck': 'n',
            'use_impactparameter_rather_than_inclination_to_fit': 'y',
            'use_psf_photometry': 'y', 'use_aperture_photometry': 'y',
            'use_legacy_psf_flux': 'n',
            'psf_seed_track_directory': None,
            'use_adaptive_apertures': False, 'bad_wcs_threshold_percent': 3.0,
            'use_aperture_corrections_and_full_image_fwhm': False,
            'pointing_rejection_sigma': None,
            'reject_overexposed_stars': True,
            'saturation_value': 65535.0,
            'overexposure_threshold_fraction': 0.9,
            'gain_electrons_per_adu': None,
            'read_noise_electrons': None,
            'dark_current_electrons_per_second_per_pixel': None,
            'flat_field_fractional_error': None,
            'telescope_aperture_m': None,
            'scintillation_coefficient': None,
            'skip_low_comparison_coverage_rejection': 'n',
            'fit_lightcurve_to_every_comparison_candidate': 'n',
            'ultranest_min_num_live_points': 200,
            'rprs_search_bound_max': 0.5,
            'restrict_rprs_range': 'y',
            'restrict_rprs_range_percentage': 10.0,
            'use_prior_rprs_when_posterior_pinned': 'y',
            'restrict_ars_range': 'y',
            'restrict_ars_range_percentage': 10.0,
            'use_sparse_posterior_live_point_retry': 'y',
        }
        self.params = {
            'images': imaging_files, 'save': save_directory, 'aavso_num': obs_code, 'second_obs': second_obs_code,
            'date': obs_date, 'lat': latitude, 'long': longitude, 'elev': elevation, 'camera': camera,
            'pixel_bin': pixel_bin, 'notes': obs_notes, 'plate_opt': plate_solution_opt, 'aavso_comp': aavso_comp,
            'tar_coords': target_star_coords, 'comp_stars': comparison_star_coords
        }

    def complete_red(self, planet):
        self.info_dict['images'] = self.params['images'](self.info_dict['images'])

        extension = 0
        hdr = fits.getheader(filename=self.info_dict['images'][0], ext=extension)
        while hdr['NAXIS'] == 0:
            extension += 1
            hdr = fits.getheader(filename=self.info_dict['images'][0], ext=extension)

        for key, value in list(self.params.items()):
            if key == 'elev':
                self.info_dict[key] = self.params[key](self.info_dict[key], self.info_dict['lat'],
                                                       self.info_dict['long'], hdr=hdr)
            elif key == 'tar_coords':
                self.info_dict[key] = self.params[key](self.info_dict[key], planet)
            elif key == 'comp_stars':
                self.info_dict[key] = self.params[key](self.info_dict[key], False)
            elif key == 'images':
                pass
            elif key in ('lat', 'long'):
                self.info_dict[key] = self.params[key](self.info_dict[key], hdr)
            elif key == 'pixel_bin':
                self.info_dict[key] = self.params[key](self.info_dict[key], hdr)
            else:
                self.info_dict[key] = self.params[key](self.info_dict[key])
            if key == 'save':
                self.info_dict['flats'], self.info_dict['darks'], self.info_dict['biases'] = \
                    image_calibrations(self.info_dict['flats'], self.info_dict['darks'],
                                       self.info_dict['biases'], self.init_opt)
                if not planet:
                    planet = planet_name(planet)
                self.info_dict['demosaic_fmt'], self.info_dict['demosaic_out'] = \
                    demosaic_settings(self.info_dict['demosaic_fmt'], self.info_dict['demosaic_out'], self.init_opt)

        return self.info_dict, planet

    def prereduced(self, planet):
        rem_list = ['images', 'plate_opt', 'aavso_comp', 'tar_coords', 'comp_stars']
        [self.params.pop(key) for key in rem_list]
        self.info_dict['aavso_comp'] = 'n'

        self.params.update({'exposure': exposure, 'file_units': data_file_units, 'file_time': data_file_time,
                            'phot_comp_star': phot_comp_star})
        self.info_dict['prered_file'] = prereduced_file(self.info_dict['prered_file'])
        aavso_overrides = parse_aavso_prereduced_overrides(self.info_dict['prered_file'])

        for key in (
            'aavso_num', 'second_obs', 'obs_name', 'lat', 'long', 'elev', 'camera', 'pixel_bin',
            'filter', 'notes', 'wl_min', 'wl_max', 'exposure', 'file_time', 'file_units',
            'dist', 'pm_ra', 'pm_dec'
        ):
            if is_blank_value(self.info_dict.get(key)) and aavso_overrides.get(key) is not None:
                self.info_dict[key] = aavso_overrides[key]
        self.info_dict['airmass_already_corrected'] = bool(aavso_overrides.get('airmass_already_corrected'))

        if not planet and not is_blank_value(aavso_overrides.get('planet')):
            planet = aavso_overrides['planet']
        if not planet:
            planet = planet_name(planet)

        for key, value in list(self.params.items()):
            if key == 'elev':
                self.info_dict[key] = self.params[key](self.info_dict[key], self.info_dict['lat'],
                                                       self.info_dict['long'], required=False)
            elif key == 'lat':
                self.info_dict[key] = self.params[key](self.info_dict[key], required=False)
            elif key == 'long':
                self.info_dict[key] = self.params[key](self.info_dict[key], required=False)
            elif key == 'phot_comp_star':
                self.info_dict[key] = self.params[key](self.info_dict[key], self.info_dict['prered_file'])
            elif key == 'date':
                continue
            else:
                self.info_dict[key] = self.params[key](self.info_dict[key])

        self.info_dict['date'] = prereduced_obs_date(
            self.info_dict.get('date'),
            self.info_dict['prered_file'],
            self.info_dict.get('file_time'),
        )

        return self.info_dict, planet

    def real_time(self, planet):
        rem_list = ['save', 'aavso_num', 'second_obs', 'date', 'lat', 'long', 'elev',
                    'camera', 'pixel_bin', 'filter', 'notes', 'plate_opt']
        [self.params.pop(key) for key in rem_list]

        for key, value in list(self.params.items()):
            if key == 'comp_stars':
                self.info_dict[key] = self.params[key](self.info_dict[key], True)
            elif key == 'tar_coords':
                self.info_dict[key] = self.params[key](self.info_dict[key], planet)
            else:
                self.info_dict[key] = self.params[key](self.info_dict[key])

                if not planet:
                    planet = planet_name(planet)

        return self.info_dict, planet

    def search_init(self, init_file, planet_dict):
        cwd = Path.cwd()

        while True:
            try:
                if not init_file:
                    log_info(f"\nYour current working directory is: {cwd}")
                    log_info(f"Potential initialization files I've found in {cwd} are: ")
                    [log_info(f"\t{file}") for file in cwd.glob('*.json') if file.is_file()]

                    init_file = user_input("\nPlease enter the Directory and Filename of "
                                           "your Initialization File: ", type_=str)
                if init_file == 'ok':
                    init_file = '/Users/rzellem/Documents/EXOTIC/inits.json'
                init_file = Path(init_file)
                planet_params = self.comp_params(init_file, planet_dict)
                return init_file, planet_params
            except (FileNotFoundError, IsADirectoryError) as e:
                log_info(f"Error: Initialization file not found. \n{e}. \nPlease try again.", error=True)

                log_info(f"\nYour current working directory is: {cwd}")
                log_info(f"Potential initialization files I've found in {cwd} are: ")
                [log_info(f"\t{file}") for file in cwd.glob('*.json') if file.is_file()]

                init_file = None
            except ValueError as e:
                log_info(f"\nError: Invalid JSON. Please reformat JSON based on given suggestion:\n\t - {e}",
                         error=True)
                init_file = None

    def comp_params(self, init_file, planet_dict):
        with init_file.open('r') as json_file:
            data = json.load(json_file)

        user_info = {
            'images': 'Directory with FITS files', 'save': 'Directory to Save Plots',
            'flats': 'Directory of Flats', 'darks': 'Directory of Darks', 'biases': 'Directory of Biases',
            'demosaic_fmt': 'Demosaic Format', 'demosaic_out': 'Demosaic Output',
            'aavso_num': ('AAVSO Observer Code (N/A if none)', 'AAVSO Observer Code (blank if none)'),
            'second_obs': ('Secondary Observer Codes (N/A if none)', 'Secondary Observer Codes (blank if none)'),
            'obs_name': 'Observatory Full Title',
            'date': 'Observation date', 'lat': 'Obs. Latitude', 'long': 'Obs. Longitude',
            'elev': ('Obs. Elevation (meters)', 'Obs. Elevation (meters; Note: leave blank if unknown)'),
            'camera': (
                'Camera Type (CCD or DSLR)',
                'Camera Type',
                'Camera Type (e.g., CCD or DSLR)',
                'Camera Type (e.g., CCD or DSLR; Note: if you are using a CMOS, please enter CCD here and then note your actual camera type in "Observing Notes")'
            ),
            'pixel_bin': 'Pixel Binning', 'filter': 'Filter Name (aavso.org/filters)',
            'notes': 'Observing Notes', 'plate_opt': 'Plate Solution? (y/n)',
            'aavso_comp': 'Add Comparison Stars from AAVSO? (y/n)',
            'tar_coords': 'Target Star X & Y Pixel', 'comp_stars': 'Comparison Star(s) X & Y Pixel',
        }
        planet_params = {
            'ra': 'Target Star RA', 'dec': 'Target Star Dec', 'pName': "Planet Name", 'sName': "Host Star Name",
            'pPer': 'Orbital Period (days)', 'pPerUnc': 'Orbital Period Uncertainty',
            'midT': ('Published Mid-Transit Time (BJD-UTC)', 'Published Mid-Transit Time'),
            'midTUnc': 'Mid-Transit Time Uncertainty',
            'rprs': ('Ratio of Planet to Stellar Radius (Rp/Rs)', 'Rp/Rs', 'Rp/R*'),
            'rprsUnc': (
                'Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty',
                'Rp/Rs Uncertainty',
                'Rp/R* Uncertainty',
            ),
            'aRs': ('Ratio of Distance to Stellar Radius (a/Rs)', 'a/Rs', 'a/R*'),
            'aRsUnc': (
                'Ratio of Distance to Stellar Radius (a/Rs) Uncertainty',
                'a/Rs Uncertainty',
                'a/R* Uncertainty',
            ),
            'inc': 'Orbital Inclination (deg)',
            'incUnc': (
                'Orbital Inclination (deg) Uncertainty',
                'Orbital Inclination (deg) Uncertainity',
                'Orbital Inclination Uncertainty',
            ),
            'ecc': ('Orbital Eccentricity (0 if null)', 'Orbital Eccentricity'),
            'teff': 'Star Effective Temperature (K)',
            'omega': 'Argument of Periastron (deg)',
            'teffUncPos': 'Star Effective Temperature (+) Uncertainty',
            'teffUncNeg': 'Star Effective Temperature (-) Uncertainty',
            'met': ('Star Metallicity ([FE/H])', 'Star Metallicity [FE/H]'),
            'metUncPos': 'Star Metallicity (+) Uncertainty',
            'metUncNeg': 'Star Metallicity (-) Uncertainty',
            'logg': 'Star Surface Gravity (log(g))', 'loggUncPos': 'Star Surface Gravity (+) Uncertainty',
            'loggUncNeg': 'Star Surface Gravity (-) Uncertainty',
            'dist': 'Star Distance (pc)',
            'pm_ra': 'Star Proper Motion RA (mas/yr)',
            'pm_dec': 'Star Proper Motion DEC (mas/yr)'
        }
        opt_info = {
            'prered_file': 'Pre-reduced File:', 'file_time': 'Pre-reduced File Time Format (BJD_TDB, JD_UTC, MJD_UTC)',
            'file_units': 'Pre-reduced File Units of Flux (flux, magnitude, millimagnitude)',
            'phot_comp_star': (
                "Comparison Star used in Photometry (leave blank if none)",
                "Comparison Star used in Photometry (blank if none)"
            ),
            'wl_min': 'Filter Minimum Wavelength (nm)', 'wl_max': 'Filter Maximum Wavelength (nm)',
            'ld_uncertainties': 'Calculate Limb Darkening Coefficients with Uncertainties? (y/n)',
            'fast_aperture_mask': ('Fast Aperture Mask (y/n)', 'Use Fast Aperture Mask (y/n)'),
            'require_comp_star': ('require_comp_star', 'Require Comparison Star? (y/n)'),
            'target_driven_comp_selection': (
                'Use target-driven comp selection rather than comp-driven comp selection',
                'target_driven_comp_selection',
            ),
            'ignore_header_wcs': (
                'Ignore WCS in Header and Do Manual Alignment? (y/n)',
                'Ignore WCS in Header and Do Manual Alignment',
                'Ignore WCS in header and do manual alignment',
                'ignore_header_wcs',
            ),
            'prefer_pixel_values_over_wcs_for_target': (
                'prefer_pixel_values_over_wcs_for_target',
                'Prefer Pixel Coordinates to WCS Coordinates if there is a conflict',
                'Prefer Pixel Coordinates to WCS Coordinates if there is a conflict? (y/n)',
            ),
            'disable_vertical_flux_normalization': (
                'disable vertical flux normalization',
                'Disable vertical flux normalization',
            ),
            'stellar_variability_only': (
                'stellar_variability_only',
                'stellar variability only',
                'Stellar Variability Only',
                'Stellar Variability Only? (y/n)',
            ),
            'use_ensemble_photometry_for_stellar_variability': (
                'use_ensemble_photometry_for_stellar_variability',
                'stellar_variability_use_ensemble',
                'Use Ensemble Photometry for Stellar Variability? (y/n)',
            ),
            'photometer_fortuitous_variables': (
                'photometer_fortuitous_variables',
                'Photometer Fortuitous Variables? (y/n)',
            ),
            'use_nextastro_vsx_cache_first': (
                'use_nextastro_vsx_cache_first',
                'Use NextAstro VSX Cache First? (y/n)',
            ),
            'detect_bad_pixels_before_photometry': (
                'detect_bad_pixels_before_photometry',
                'Detect Bad Pixels Before Photometry? (y/n)',
            ),
            'multiprocess_bad_pixel_precheck': (
                'multiprocess_bad_pixel_precheck',
                'Multiprocess Bad-Pixel Precheck? (y/n or process count)',
                'Multiprocess Bad Pixel Precheck? (y/n or process count)',
            ),
            'detrend_on_outoftransit_baseline': (
                'detrend_on_outoftransit_baseline',
                'Detrend on Out-of-Transit Baseline',
                'detrend_on_out_of_transit_baseline',
            ),
            'final_fit_baseline_duration_multiplier': (
                'final_fit_baseline_duration_multiplier',
                'Final Fit Baseline Duration Multiplier',
            ),
            'use_eebls_to_initialize_tmid_and_bounds': (
                'use_eebls_to_initialize_tmid_and_bounds',
                'Use EEBLS to Initialize Tmid and Bounds? (y/n)',
                'Use EEBLS To Initialize Tmid And Bounds? (y/n)',
            ),
            'pick_comparison_by_eebls_snr': (
                'pick_comparison_by_eebls_snr',
                'Pick Comparison by EEBLS SNR? (y/n)',
                'Pick comparison by EEBLS SNR? (y/n)',
            ),
            'use_deviation_from_expected_transit_in_qc': (
                'use_deviation_from_expected_transit_in_qc',
                'Use Deviation From Expected Transit In QC? (y/n)',
            ),
            'deviation_from_expected_transit_in_qc_sigma': (
                'deviation_from_expected_transit_in_qc_sigma',
                'Deviation From Expected Transit In QC Sigma',
            ),
            'run_final_fit_phase_residual_clip': (
                'run_final_fit_phase_residual_clip',
                'Run Final-Fit Phase Residual Clip? (y/n)',
                'Run Final Fit Phase Residual Clip? (y/n)',
            ),
            'exit_at_first_qc_pass_solution': (
                'exit_at_first_qc_pass_solution',
                'exit at first QC PASS solution',
                'Exit at first QC PASS solution',
                'Exit at first QC PASS solution? (y/n)',
                'Exit At First QC PASS Solution? (y/n)',
            ),
            'use_impactparameter_rather_than_inclination_to_fit': (
                'use_impactparameter_rather_than_inclination_to_fit',
                'Use impact parameter rather than inclination to fit? (y/n)',
                'Use Impact Parameter Rather Than Inclination To Fit? (y/n)',
            ),
            'use_psf_photometry': (
                'use_psf_photometry',
                'Use PSF Photometry? (y/n)',
            ),
            'use_legacy_psf_flux': (
                'use_legacy_psf_flux',
                'legacy_psf_flux_mode',
                'Use Legacy PSF Flux? (y/n)',
                'Use Legacy PSF Flux Mode? (y/n)',
            ),
            'psf_seed_track_directory': (
                'psf_seed_track_directory',
                'legacy_psf_seed_track_directory',
                'PSF Seed Track Directory',
                'Legacy PSF Seed Track Directory',
            ),
            'use_aperture_photometry': (
                'use_aperture_photometry',
                'Use Aperture Photometry? (y/n)',
            ),
            'use_adaptive_apertures': (
                'use_adaptive_apertures',
                'Use Adaptive Apertures? (y/n)',
                'Use Adaptive Apertures (y/n)',
            ),
            'use_aperture_corrections_and_full_image_fwhm': (
                'use_aperture_corrections_and_full_image_fwhm',
                'Use Aperture Corrections and Full Image FWHM? (y/n)',
                'Use Aperture Corrections And Full Image FWHM? (y/n)',
            ),
            'reject_overexposed_stars': (
                'reject_overexposed_stars',
                'Reject Overexposed Stars? (y/n)',
                'Reject Overexposed Target and Comparison Stars? (y/n)',
            ),
            'saturation_value': (
                'saturation_value',
                'saturation_value_adu',
                'Saturation Value',
                'SATURATE',
            ),
            'overexposure_threshold_fraction': (
                'overexposure_threshold_fraction',
                'Overexposure Threshold Fraction',
                'Saturation Rejection Threshold Fraction',
            ),
            'skip_low_comparison_coverage_rejection': (
                'skip_low_comparison_coverage_rejection',
                'Skip Low Comparison Coverage Rejection? (y/n)',
            ),
            'fit_lightcurve_to_every_comparison_candidate': (
                'fit_lightcurve_to_every_comparison_candidate',
                'Fit Lightcurve to Every Comparison Candidate? (y/n)',
            ),
            'automatic_optimal_calibration_selector': (
                'automatic_optimal_calibration_selector',
                'Automatic Optimal Calibration Selector? (y/n)',
            ),
            'automatic_optimal_calibration_selector_count': (
                'automatic_optimal_calibration_selector_count',
                'Automatic Optimal Calibration Selector Count',
                'automatic_optimal_calibration_selector_max_stars',
                'Automatic Optimal Calibration Selector Max Stars',
            ),
            'colour_term': (
                'colour_term',
                'color_term',
                'COLTERM',
            ),
            'colour_term_error': (
                'colour_term_error',
                'color_term_error',
                'COLTERR',
            ),
            'colour_term_index': (
                'colour_term_index',
                'color_term_index',
                'COLTIDX',
            ),
            'colour_term_bv': (
                'colour_term_bv',
                'color_term_bv',
                'COLTBV',
            ),
            'colour_term_bv_error': (
                'colour_term_bv_error',
                'color_term_bv_error',
                'COLTBVER',
                'COLTBVERR',
            ),
            'colour_term_bprp': (
                'colour_term_bprp',
                'color_term_bprp',
                'COLTBPRP',
            ),
            'colour_term_bprp_error': (
                'colour_term_bprp_error',
                'color_term_bprp_error',
                'CBPRPERR',
                'COLTBPRPERR',
            ),
            'colour_equation_filter': (
                'colour_equation_filter',
                'color_equation_filter',
                'COLEQFIL',
            ),
            'use_ensemble_photometry_rather_than_single_comp': (
                'use_ensemble_photometry_rather_than_single_comp',
                'Use Ensemble Photometry Rather Than Single Comp? (y/n)',
            ),
            'ultranest_min_num_live_points': (
                'Minimum Number of Live Points for UltraNest',
                'minimum number of live points for ultranest',
                'ultranest_min_num_live_points',
                'ultranest_min_live_points',
                'min_num_live_points',
            ),
            'run_fast_ultranest_before_final_run': (
                'run fast ultranest before final run',
                'Run Fast UltraNest Before Final Run? (y/n)',
                'run_fast_ultranest_before_final_run',
            ),
            'rprs_search_bound_max': (
                'rprs_search_bound_max',
                'max_rprs_search_bound',
                'maximum rprs search bound',
                'maximum Rp/Rs search bound',
                'maximum Rp/R* search bound',
                'Maximum Rp/Rs Search Bound',
                'Maximum Rp/R* Search Bound',
            ),
            'restrict_rprs_range': (
                'restrict_Rp/Rs_range',
                'restrict_Rp/R*_range',
                'restrict_rprs_range',
                'restrict_RpRs_range',
                'Restrict Rp/Rs Range? (y/n)',
                'Restrict Rp/R* Range? (y/n)',
            ),
            'restrict_rprs_range_percentage': (
                'restrict_Rp/Rs_range_percentage',
                'restrict_Rp/R*_range_percentage',
                'restrict_rprs_range_percentage',
                'restrict_RpRs_range_percentage',
                'Restrict Rp/Rs Range Percentage',
                'Restrict Rp/R* Range Percentage',
            ),
            'use_prior_rprs_when_posterior_pinned': (
                'use_prior_Rp/Rs_when_posterior_pinned',
                'use_prior_Rp/R*_when_posterior_pinned',
                'use_prior_rprs_when_posterior_pinned',
                'use_prior_RpRs_when_posterior_pinned',
                'Use Prior Rp/Rs When Posterior Pinned? (y/n)',
                'Use Prior Rp/R* When Posterior Pinned? (y/n)',
            ),
            'restrict_ars_range': (
                'restrict_a/Rs_range',
                'restrict_a/R*_range',
                'restrict_ars_range',
                'restrict_aRs_range',
                'Restrict a/Rs Range? (y/n)',
                'Restrict a/R* Range? (y/n)',
            ),
            'restrict_ars_range_percentage': (
                'restrict_a/Rs_range_percentage',
                'restrict_a/R*_range_percentage',
                'restrict_ars_range_percentage',
                'restrict_aRs_range_percentage',
                'Restrict a/Rs Range Percentage',
                'Restrict a/R* Range Percentage',
            ),
            'use_sparse_posterior_live_point_retry': (
                'use_sparse_posterior_live_point_retry',
                'Use Sparse Posterior Live-Point Retry? (y/n)',
                'Use Sparse Posterior Live Point Retry? (y/n)',
                'Sparse Posterior Live-Point Retry? (y/n)',
            ),
            'bad_wcs_threshold_percent': (
                'bad_wcs_threshold_percent',
                'Bad WCS Threshold Percent',
            ),
            'pointing_rejection_sigma': (
                'pointing_rejection_sigma',
                'Pointing Rejection Sigma',
            ),
            'gain_electrons_per_adu': (
                'gain_electrons_per_adu',
                'gain_e_per_adu',
                'Gain (e-/ADU)',
                'CCD Gain (e-/ADU)',
            ),
            'read_noise_electrons': (
                'read_noise_electrons',
                'read_noise_e',
                'read_noise',
                'Read Noise (e-)',
                'CCD Read Noise (e-)',
            ),
            'dark_current_electrons_per_second_per_pixel': (
                'dark_current_electrons_per_second_per_pixel',
                'dark_current_e_per_s_pix',
                'dark_current',
                'Dark Current (e-/s/pix)',
            ),
            'flat_field_fractional_error': (
                'flat_field_fractional_error',
                'flat_field_fractional_noise',
                'flat_field_error_fraction',
                'Flat Field Fractional Error',
                'Flat-Field Fractional Error',
            ),
            'telescope_aperture_m': (
                'telescope_aperture_m',
                'telescope_aperture_meters',
                'Telescope Aperture (m)',
            ),
            'scintillation_coefficient': (
                'scintillation_coefficient',
                'scintillation_noise_coefficient',
                'Scintillation Coefficient',
            ),
            'pixel_scale': ('Image Scale (Ex: 5.21 arcsecs/pixel)', 'Pixel Scale (Ex: 5.21 arcsecs/pixel)',
                            'Pixel Scale (arsec/pixel)'),
            'exposure': 'Exposure Time (s)',
            'random_seed': 'Random Seed'
        }

        self.info_dict = init_params(user_info, self.info_dict, data['user_info'])
        if self.info_dict['aavso_comp'] is None:
            self.info_dict['aavso_comp'] = 'n'
        self.info_dict = init_params(opt_info, self.info_dict, data['optional_info'])
        planet_dict = init_params(planet_params, planet_dict, data['planetary_parameters'])
        return populate_missing_gaia_astrometry(planet_dict)


def check_imaging_files(directory, img_type):
    file_extensions = ['.fits', '.fit', '.fts', '.fz', '.fits.gz', '.fit.gz', '.fits.fz', 'fit.fz']
    input_files = []

    while True:
        try:
            directory = Path(directory)
            if directory.is_dir() and str(directory).strip():
                for ext in file_extensions:
                    for file in directory.iterdir():
                        if file.is_file() and file.name.lower().endswith(ext.lower()) \
                                and file.name[0:2] not in ('ref', 'wcs'):
                            input_files.append(str(file))
                    if input_files:
                        return input_files
                if not input_files:
                    raise FileNotFoundError
            else:
                raise NotADirectoryError
        except FileNotFoundError:
            log_info(f"\nError: {img_type} files not found with .fits, .fit, .fts, .fz, fit.gz or .fits.gz extensions "
                     f"in {directory}.", error=True)
            opt = user_input("\nWould you like to enter in an alternate image extension in addition to .FITS? (y/n): ",
                             type_=str, values=['y', 'n'])
            if opt == 'y':
                add_ext = user_input("Please enter the extension you want to add (EX: .FITS): ", type_=str)
                file_extensions.append(add_ext)
            else:
                directory = user_input(f"Enter the directory path where {img_type} files are located "
                                       f"(Example using the sample data: sample-data/HatP32Dec202017): ", type_=str)
        except (NotADirectoryError, OSError):
            log_info("\nError: No such directory exists when searching for FITS files. Please try again.", error=True)
            directory = user_input(f"Enter the directory path where {img_type} files are located "
                                   f"(Example using the sample data: sample-data/HatP32Dec202017): ", type_=str)


def imaging_files(directory):
    if not directory:
        directory = user_input("\nEnter the directory path where imaging files are located. "
                               "(Example using the sample data: sample-data/HatP32Dec202017): ", type_=str)
    return check_imaging_files(directory, 'Imaging')


def save_directory(directory):
    while True:
        try:
            if not directory:
                directory = user_input("\nEnter the directory to Save the Results and Plots into "
                                       "or type new to create one: ", type_=str)
            if directory == 'new':
                directory = create_directory()
            else:
                if not Path(directory).is_dir() or directory.replace(' ', '') == '':
                    raise NotADirectoryError
            return directory
        except (NotADirectoryError, OSError):
            log_info("Error: The directory entered does not exist. Please try again. Make sure to follow this "
                     "\nformatting (using whichever directory you choose): /sample-data/results", error=True)
            directory = None


def create_directory():
    save_path = Path.cwd()
    while True:
        directory = user_input("Enter the name for your new directory: ", type_=str)
        try:
            save_path = save_path / directory
            Path(save_path).mkdir()
        except OSError:
            log_info(f"Error: Creation of the directory {save_path}/{directory} failed.", error=True)
        else:
            log_info(f"Successfully created the directory {save_path}.")
            return save_path


def image_calibrations(flats_dir, darks_dir, biases_dir, init):
    opt, flats_list, darks_list, biases_list = None, None, None, None

    if init == 'n':
        opt = user_input("\nDo you have any Calibration Images? (Flats, Darks or Biases)? (y/n): ",
                         type_=str, values=['y', 'n'])

    if opt == 'y' or flats_dir:
        flats_list = check_calibration(flats_dir, 'Flats')
    if opt == 'y' or darks_dir:
        darks_list = check_calibration(darks_dir, 'Darks')
    if opt == 'y' or biases_dir:
        biases_list = check_calibration(biases_dir, 'Biases')

    return flats_list, darks_list, biases_list


def check_calibration(directory, image_type):
    if not directory:
        opt = user_input(f"\nDo you have {image_type}? (y/n): ", type_=str, values=['y', 'n'])
        if opt == 'y':
            directory = user_input(f"Please enter the directory path to your {image_type} "
                                   "(must be in their own separate folder): ", type_=str)
    if directory:
        return check_imaging_files(directory, image_type)
    return None

def demosaic_settings(demosaic_fmt, demosaic_out, init):
    opt = None
    if init == 'n':
        opt = user_input("\nAre images color and require demosaicing? (y/n): ",
            type_=str, values=['y', 'n'])
    if opt == 'y':
        if not demosaic_fmt:
            demosaic_fmt = user_input(f"\nWhat is Bayer pattern for camera? (RGGB, BGGR, GRBG, GBRG): ", type_=str, values=['rggb', 'bggr', 'grbg', 'gbrg'])
            demosaic_fmt = demosaic_fmt.upper()
        if not demosaic_out:
            demosaic_out = user_input(f"\nWhat color channel should be processed? (gray, red, green, blue, blueblock, custom): ", type_=str, values=['gray', 'red', 'green', 'blue', 'blueblock', 'custom'])
        if demosaic_out == 'custom':
            demosaic_red = user_input("\nWhat weight for red channel (0.0-1.0)?", type_=float)
            demosaic_green = user_input("\nWhat weight for green channel (0.0-1.0)?", type_=float)
            demosaic_blue = user_input("\nWhat weight for blue channel (0.0-1.0)?", type_=float)
            demosaic_out = [ demosaic_red, demosaic_green, demosaic_blue ]

    return demosaic_fmt, demosaic_out

def planet_name(planet):
    if not planet:
        planet = user_input("\nPlease enter Planet's name: ", type_=str)
    return planet


def obs_code(code):
    if code is None:
        code = user_input("Please enter your AAVSO Observer Account Number "
                          "(if none, leave blank and press enter): ", type_=str)
    code = code.replace(' ', '')
    if code.lower() == 'n/a':
        code = ""
    return code


def second_obs_code(code):
    if code is None:
        code = user_input("Please enter your comma-separated Secondary Observer Codes "
                          "(if none, leave blank and press enter): ", type_=str)
    code = code.replace(' ', '')
    if code.lower() == 'n/a':
        code = ""
    return code


def obs_date(date):
    while True:
        if not date:
            date = user_input("\nPlease enter the Observation Date: ", type_=str)
        if date:
            break
    if '/' in date:
        date = date.replace('/', '-')
    return date


def normalize_obs_date(date):
    if is_blank_value(date):
        return None

    date = str(date).strip()
    if re.fullmatch(r'\d{8}', date):
        date = f"{date[0:4]}-{date[4:6]}-{date[6:8]}"
    if '/' in date:
        date = date.replace('/', '-')
    return date


def latitude(lat, hdr=None, required=True):
    while True:
        if is_blank_value(lat):
            if hdr:
                lat = find(hdr, ['LATITUDE', 'LAT', 'SITELAT'])
                if lat:
                    return lat
            if not required:
                return None
            lat = user_input("Enter the latitude (in degrees) of where you observed. "
                             "(Don't forget the sign where North is '+' and South is '-')! "
                             "(Example: -32.12): ", type_=str)
        lat = str(lat).strip()

        if lat[0] == '+' or lat[0] == '-':
            # Convert to float if latitude in decimal. If latitude is in +/-HH:MM:SS format, convert to a float.
            try:
                lat = float(lat.strip())
            except ValueError:
                lat = float(process_lat_long(lat, 'latitude'))

            if -90.00 <= lat <= 90.00:
                return lat
            else:
                log_info("Error: Your latitude is out of range. "
                         "Please enter a latitude between -90 and +90 (deg).", error=True)
        else:
            log_info("Error: You forgot the sign for the latitude! North is '+' and South is '-'. Please try again.",
                     error=True)
        lat = None


def longitude(long, hdr=None, required=True):
    while True:
        if is_blank_value(long):
            if hdr:
                long = find(hdr, ['LONGITUD', 'LONG', 'LONGITUDE', 'SITELONG'])
                if long:
                    return long
            if not required:
                return None
            long = user_input("Enter the longitude (in degrees) of where you observed. "
                              "(Don't forget the sign where East is '+' and West is '-')! "
                              "(Example: +152.51): ", type_=str)
        long = str(long).strip()

        if long[0] == '+' or long[0] == '-':
            # Convert to float if longitude in decimal. If longitude is in +/-HH:MM:SS format, convert to a float.
            try:
                long = float(long.strip())
            except ValueError:
                long = float(process_lat_long(long, 'longitude'))

            if -180.00 <= long <= 180.00:
                return long
            else:
                log_info("Error: Your longitude is out of range. "
                         "Please enter a longitude between -180 and +180 (deg).", error=True)
        else:
            log_info("Error: You forgot the sign for the longitude! East is '+' and West is '-'. Please try again.",
                     error=True)
        long = None


def elevation(elev, lat, long, hdr=None, required=True):
    while True:
        try:
            if is_blank_value(elev):
                elev = None
            else:
                elev = typecast_check(type_=float, val=elev)
                if elev is False:
                    raise ValueError

            if elev is None:
                if hdr:
                    elev = find(hdr, ['HEIGHT', 'ELEVATION', 'ELE', 'EL', 'OBSGEO-H', 'ALT-OBS', 'SITEELEV'])
                    if not is_blank_value(elev):
                        return float(elev)
                if not required:
                    return None
                log_info("\nEXOTIC is retrieving elevation based on entered "
                         "latitude and longitude from Open Elevation.")
                animate_toggle(True)
                elev = open_elevation(lat, long)
                animate_toggle()
                if elev is False:
                    log_info("\nWarning: EXOTIC could not retrieve elevation.", warn=True)
                    elev = user_input("Enter the elevation (in meters) of where you observed: ", type_=float)
            return elev
        except ValueError:
            log_info("Error: The entered elevation is incorrect.", error=True)
            elev = None


def camera(c_type):
    if isinstance(c_type, str) and "DSLR" in c_type.strip().upper():
        return "DSLR"
    return "CCD"


def format_fits_binning_axis(value):
    if is_blank_value(value):
        return None

    try:
        binning_value = float(str(value).strip())
    except (TypeError, ValueError):
        return None

    if not math.isfinite(binning_value) or binning_value <= 0:
        return None
    if binning_value.is_integer():
        return str(int(binning_value))
    return str(binning_value)


def fits_header_pixel_bin(hdr):
    if hdr is None:
        return None

    x_binning = format_fits_binning_axis(find(hdr, ['XBINNING']))
    y_binning = format_fits_binning_axis(find(hdr, ['YBINNING']))
    if x_binning is None or y_binning is None:
        return None
    return f"{x_binning}x{y_binning}"


def pixel_bin(pix_bin, hdr=None):
    if is_blank_value(pix_bin):
        pix_bin = fits_header_pixel_bin(hdr)
    if is_blank_value(pix_bin):
        pix_bin = user_input("Please enter the pixel binning: ", type_=str)
    return pix_bin


def obs_notes(notes):
    if not isinstance(notes, str):
        notes = user_input("Please enter any observing notes (seeing, weather, etc.) or leave blank and press enter: ",
                           type_=str)
    if not notes.replace(' ', ''):
        notes = "na"
    return notes


def plate_solution_opt(opt):
    if opt:
        opt = opt.lower().strip()
    if opt not in ('y', 'n'):
        opt = user_input("\nWould you like to upload the your image for a plate solution?"
                         "\nThis will allow EXOTIC to translate your image's pixels into coordinates on the sky."
                         "\nDISCLAIMER: One of your imaging files will be publicly viewable on "
                         "nova.astrometry.net. (y/n): ", type_=str, values=['y', 'n'])
    return opt


def aavso_comp(opt):
    if opt:
        opt = opt.lower().strip()
    if opt not in ('y', 'n'):
        opt = user_input("\nWould you like Comparison Stars added automatically from AAVSO? (y/n): ",
                         type_=str, values=['y', 'n'])
    return opt


def target_star_coords(coords, planet):
    if isinstance(coords, list) and len(coords) == 2:
        pass
    elif isinstance(coords, str) and any(str.isdigit(x) for x in coords):
        coords = re.findall(r"[-+]?(?:\d*\.?\d+)", coords)
        coords = [int(float(coord)) for coord in coords]
    else:
        coords = [user_input(f"\nPlease enter {planet}'s X Pixel Coordinate: ", type_=int),
                  user_input(f"\nPlease enter {planet}'s Y Pixel Coordinate: ", type_=int)]

    return coords


def comparison_star_coords(comp_stars, rt_bool):
    if isinstance(comp_stars, list) and len(comp_stars) >= 1 and \
            all(isinstance(star, list) for star in comp_stars):
        comp_stars = [star for star in comp_stars if star != []]
    elif isinstance(comp_stars, str) and any(str.isdigit(x) for x in comp_stars):
        comp_stars = re.findall(r"[-+]?(?:\d*\.?\d+)", comp_stars)
        comp_stars = [int(float(comp_star)) for comp_star in comp_stars]
        comp_stars = [comp_stars[i:i+2] for i in range(0, len(comp_stars), 2)]
    else:
        comp_stars = []

    if not comp_stars:
        while True:
            if not rt_bool:
                num_comp_stars = user_input("\nHow many Comparison Stars would you like to use? (1 or more): ", type_=int)
                if num_comp_stars >= 1:
                    break
                log_info("\nError: The number of Comparison Stars entered is incorrect.", error=True)
            else:
                num_comp_stars = 1
                break

        for num in range(num_comp_stars):
            x_pix = user_input(f"\nComparison Star {num + 1} X Pixel Coordinate: ", type_=int)
            y_pix = user_input(f"Comparison Star {num + 1} Y Pixel Coordinate: ", type_=int)
            comp_stars.append([x_pix, y_pix])

    if rt_bool and isinstance(comp_stars[0], list):
        comp_stars = comp_stars[0]

    return comp_stars


def exposure(exp):
    exp = typecast_check(type_=float, val=exp)
    if not exp:
        exp = user_input("Please enter your exposure time (seconds): ", type_=float)
    return exp


def prereduced_file(file):
    while True:
        try:
            if not file:
                file = user_input("Enter the path and file name of your data file: ", type_=str)
            if file == "ok":
                file = "/Users/rzellem/Documents/EXOTIC/sample-data/NormalizedFluxHAT-P-32 bDecember 17, 2017.txt"
                # file = "/Users/rzellem/Downloads/fluxorama.csv
                log_info("Hello, Rob.")

            file = Path(file)

            if file.is_file():
                return file
            else:
                raise FileNotFoundError
        except FileNotFoundError:
            log_info("Error: Data file not found. Please try again.", error=True)
            file = None


def blank_phot_comp_star():
    return {key: '' for key in PHOT_COMP_STAR_KEYS}


def normalize_phot_comp_star(comp_star):
    normalized_comp_star = blank_phot_comp_star()

    if not isinstance(comp_star, dict):
        return normalized_comp_star

    for key in PHOT_COMP_STAR_KEYS:
        value = comp_star.get(key, '')
        if value is None:
            continue

        value = str(value).strip()
        normalized_comp_star[key] = '' if value.lower() in ('null', 'none') else value

    return normalized_comp_star


def read_aavso_metadata(prereduced_file_path):
    return dict(parse_aavso_metadata(prereduced_file_path) or [])


def parse_aavso_metadata(prereduced_file_path):
    if not prereduced_file_path:
        return {}

    try:
        with Path(prereduced_file_path).open('r', encoding='utf-8') as file:
            for line in file:
                metadata_line = line.strip()
                if not metadata_line:
                    continue
                if not metadata_line.startswith('#'):
                    break
                if '=' not in metadata_line:
                    continue

                key, value = metadata_line[1:].split('=', 1)
                yield key.strip().upper(), value.strip()
    except (FileNotFoundError, OSError, TypeError):
        return


def first_aavso_metadata_value(metadata, aliases):
    for key in aliases:
        value = metadata.get(key)
        if not is_blank_value(value):
            return value
    return None


def first_aavso_metadata_text(metadata, aliases, allow_blank=False):
    for key in aliases:
        if key not in metadata:
            continue

        value = metadata.get(key)
        if value is None:
            return '' if allow_blank else None

        value = str(value).strip()
        if allow_blank:
            return '' if value.lower() in ('null', 'none') else value
        if not is_blank_value(value):
            return value
    return None


def normalize_aavso_code(value):
    if value is None:
        return None
    value = str(value).strip()
    return '' if value.lower() in ('', 'n/a', 'na', 'null', 'none') else value


def normalize_aavso_blankable_text(value):
    if value is None:
        return None
    value = str(value).strip()
    return '' if value.lower() in ('null', 'none') else value


def normalize_aavso_coordinate_text(value):
    if is_blank_value(value):
        return None

    value = str(value).strip()
    if value[0] in ('+', '-'):
        return value

    try:
        numeric_value = float(value)
    except ValueError:
        return value

    if numeric_value >= 0:
        return f"+{value}"
    return value


def format_aavso_numeric_string(value):
    value = float(value)
    if value.is_integer():
        return f"{value:.1f}"
    return str(value)


def convert_aavso_wavelength_to_nm(value, units='nm'):
    if is_blank_value(value):
        return None

    units_key = 'nm' if units is None else str(units).strip().lower()
    units_key = units_key.replace('µ', 'u').replace('μ', 'u')
    factor = AAVSO_WAVELENGTH_UNIT_FACTORS_TO_NM.get(units_key)
    if factor is None:
        return None

    try:
        return format_aavso_numeric_string(float(str(value).strip()) * factor)
    except ValueError:
        return None


def parse_aavso_json(value):
    if is_blank_value(value):
        return None

    try:
        return json.loads(value)
    except (TypeError, json.JSONDecodeError):
        return None


def parse_aavso_comp_star_from_metadata(metadata):
    comp_star_json = metadata.get('COMP_STAR-XC')
    comp_star = parse_aavso_json(comp_star_json)
    if comp_star is None:
        return blank_phot_comp_star()
    return normalize_phot_comp_star(comp_star)


def lookup_aavso_filter_metadata(*candidates):
    for candidate in candidates:
        lookup_key = normalize_aavso_filter_lookup_key(candidate)
        if lookup_key and lookup_key in AAVSO_FILTER_LOOKUP:
            return AAVSO_FILTER_LOOKUP[lookup_key]
    return None


def parse_aavso_filter_xc_fwhm(filter_metadata):
    if not isinstance(filter_metadata, dict):
        return None, None

    fwhm = filter_metadata.get('fwhm')
    if isinstance(fwhm, dict):
        values = [
            convert_aavso_wavelength_to_nm(fwhm.get('min'), fwhm.get('units', 'nm')),
            convert_aavso_wavelength_to_nm(fwhm.get('max'), fwhm.get('units', 'nm')),
        ]
    elif isinstance(fwhm, (list, tuple)):
        values = []
        for item in fwhm[:2]:
            if isinstance(item, dict):
                values.append(convert_aavso_wavelength_to_nm(item.get('value'), item.get('units', 'nm')))
            else:
                values.append(convert_aavso_wavelength_to_nm(item))
    else:
        values = []

    values = [value for value in values if value is not None]
    if len(values) < 2:
        return None, None

    values = sorted(values[:2], key=float)
    return values[0], values[1]


def parse_aavso_filter_metadata_from_metadata(metadata):
    filter_value = first_aavso_metadata_text(metadata, AAVSO_FILTER_HEADER_KEYS)
    filter_xc = parse_aavso_json(first_aavso_metadata_text(metadata, AAVSO_FILTER_XC_HEADER_KEYS))

    parsed_filter = {
        'filter': filter_value,
        'filter_desc': None,
        'wl_min': None,
        'wl_max': None,
    }

    if isinstance(filter_xc, dict):
        filter_name = normalize_aavso_blankable_text(filter_xc.get('name'))
        filter_desc = normalize_aavso_blankable_text(filter_xc.get('desc'))
        if is_blank_value(parsed_filter['filter']):
            parsed_filter['filter'] = filter_name or filter_desc
        if filter_desc:
            parsed_filter['filter_desc'] = filter_desc

        wl_min, wl_max = parse_aavso_filter_xc_fwhm(filter_xc)
        if wl_min is not None and wl_max is not None:
            parsed_filter['wl_min'] = wl_min
            parsed_filter['wl_max'] = wl_max

    filter_record = lookup_aavso_filter_metadata(
        parsed_filter['filter'],
        parsed_filter['filter_desc'],
    )
    if filter_record is not None:
        if is_blank_value(parsed_filter['filter']):
            parsed_filter['filter'] = filter_record.get('name') or filter_record.get('desc')
        if is_blank_value(parsed_filter['filter_desc']):
            parsed_filter['filter_desc'] = filter_record.get('desc')
        if parsed_filter['wl_min'] is None:
            parsed_filter['wl_min'] = filter_record['fwhm'][0]
        if parsed_filter['wl_max'] is None:
            parsed_filter['wl_max'] = filter_record['fwhm'][1]

    return parsed_filter


def parse_aavso_time_format_from_metadata(metadata):
    value = first_aavso_metadata_text(metadata, AAVSO_TIME_FORMAT_HEADER_KEYS)
    if is_blank_value(value):
        return None

    normalized = value.upper().strip().replace('-', '_').replace(' ', '_')
    if normalized in AAVSO_ALLOWED_FILE_TIME_FORMATS:
        return normalized
    if normalized == 'BJD':
        return 'BJD_TDB'
    if normalized == 'JD':
        return 'JD_UTC'
    if normalized == 'MJD':
        return 'MJD_UTC'
    return None


def parse_aavso_measurement_units_from_metadata(metadata):
    value = first_aavso_metadata_text(metadata, AAVSO_MEASUREMENT_TYPE_HEADER_KEYS)
    if is_blank_value(value):
        return None

    normalized = re.sub(r'[\W_]+', '', value.lower())
    if 'millimag' in normalized or normalized == 'mmag':
        return 'millimagnitude'
    if 'flux' in normalized:
        return 'flux'
    if 'mag' in normalized:
        return 'magnitude'
    return None


def parse_aavso_exposure_from_metadata(metadata):
    value = first_aavso_metadata_text(metadata, AAVSO_EXPOSURE_HEADER_KEYS)
    if is_blank_value(value):
        return None

    try:
        return float(str(value).strip())
    except ValueError:
        return None


def parse_aavso_airmass_detrend_from_metadata(metadata):
    value = first_aavso_metadata_text(metadata, AAVSO_DETREND_PARAMETER_HEADER_KEYS, allow_blank=True)
    if is_blank_value(value):
        return False

    detrend_parameters = {
        re.sub(r'\s+', ' ', item.strip()).upper()
        for item in re.split(r'[;,]', value)
        if item.strip()
    }
    return (
        'AIRMASS' in detrend_parameters
        and 'AIRMASS CORRECTION FUNCTION' in detrend_parameters
    )


def parse_aavso_prereduced_overrides(prereduced_file_path):
    metadata = read_aavso_metadata(prereduced_file_path)
    filter_metadata = parse_aavso_filter_metadata_from_metadata(metadata)

    return {
        'aavso_num': normalize_aavso_code(first_aavso_metadata_text(metadata, AAVSO_TEXT_HEADER_KEYS['aavso_num'], allow_blank=True)),
        'second_obs': normalize_aavso_code(first_aavso_metadata_text(metadata, AAVSO_TEXT_HEADER_KEYS['second_obs'], allow_blank=True)),
        'obs_name': normalize_aavso_blankable_text(first_aavso_metadata_text(metadata, AAVSO_TEXT_HEADER_KEYS['obs_name'], allow_blank=True)),
        'date': normalize_obs_date(first_aavso_metadata_value(metadata, AAVSO_OBSDATE_HEADER_KEYS)),
        'lat': normalize_aavso_coordinate_text(first_aavso_metadata_value(metadata, AAVSO_LOCATION_HEADER_KEYS['lat'])),
        'long': normalize_aavso_coordinate_text(first_aavso_metadata_value(metadata, AAVSO_LOCATION_HEADER_KEYS['long'])),
        'elev': first_aavso_metadata_value(metadata, AAVSO_LOCATION_HEADER_KEYS['elev']),
        'camera': first_aavso_metadata_text(metadata, AAVSO_TEXT_HEADER_KEYS['camera']),
        'pixel_bin': first_aavso_metadata_text(metadata, AAVSO_TEXT_HEADER_KEYS['pixel_bin']),
        'filter': filter_metadata['filter'],
        'filter_desc': filter_metadata['filter_desc'],
        'wl_min': filter_metadata['wl_min'],
        'wl_max': filter_metadata['wl_max'],
        'notes': normalize_aavso_blankable_text(first_aavso_metadata_text(metadata, AAVSO_TEXT_HEADER_KEYS['notes'], allow_blank=True)),
        'file_time': parse_aavso_time_format_from_metadata(metadata),
        'file_units': parse_aavso_measurement_units_from_metadata(metadata),
        'exposure': parse_aavso_exposure_from_metadata(metadata),
        'dist': first_aavso_metadata_text(metadata, AAVSO_GAIA_HEADER_KEYS['dist']),
        'pm_ra': first_aavso_metadata_text(metadata, AAVSO_GAIA_HEADER_KEYS['pm_ra']),
        'pm_dec': first_aavso_metadata_text(metadata, AAVSO_GAIA_HEADER_KEYS['pm_dec']),
        'airmass_already_corrected': parse_aavso_airmass_detrend_from_metadata(metadata),
        'phot_comp_star': parse_aavso_comp_star_from_metadata(metadata),
        'planet': first_aavso_metadata_text(metadata, AAVSO_TEXT_HEADER_KEYS['planet']),
        'host_star': first_aavso_metadata_text(metadata, AAVSO_TEXT_HEADER_KEYS['host_star']),
    }


def parse_aavso_location(prereduced_file_path):
    metadata = read_aavso_metadata(prereduced_file_path)
    parsed_location = {
        key: first_aavso_metadata_value(metadata, aliases)
        for key, aliases in AAVSO_LOCATION_HEADER_KEYS.items()
    }
    parsed_location['lat'] = normalize_aavso_coordinate_text(parsed_location['lat'])
    parsed_location['long'] = normalize_aavso_coordinate_text(parsed_location['long'])
    return parsed_location


def parse_aavso_obsdate(prereduced_file_path):
    metadata = read_aavso_metadata(prereduced_file_path)
    return normalize_obs_date(first_aavso_metadata_value(metadata, AAVSO_OBSDATE_HEADER_KEYS))


def parse_aavso_comp_star(prereduced_file_path):
    metadata = read_aavso_metadata(prereduced_file_path)
    return parse_aavso_comp_star_from_metadata(metadata)


def phot_comp_star(comp_star, prereduced_file_path=None):
    if isinstance(comp_star, dict):
        return normalize_phot_comp_star(comp_star)
    return parse_aavso_comp_star(prereduced_file_path)


def first_prereduced_timestamp(prereduced_file_path):
    if not prereduced_file_path:
        return None

    try:
        with Path(prereduced_file_path).open('r', encoding='utf-8') as file:
            for line in file:
                data_line = line.strip()
                if not data_line or data_line.startswith('#'):
                    continue

                first_column = re.split(r'[\s,]+', data_line, maxsplit=1)[0]
                try:
                    return float(first_column)
                except ValueError:
                    continue
    except (FileNotFoundError, OSError, TypeError):
        return None

    return None


def obs_date_from_first_prereduced_entry(prereduced_file_path, time_format):
    first_timestamp = first_prereduced_timestamp(prereduced_file_path)
    if first_timestamp is None:
        return None

    try:
        if time_format == 'MJD_UTC':
            return Time(first_timestamp, format='mjd', scale='utc').to_value('iso', subfmt='date')
        if time_format == 'BJD_TDB':
            return Time(first_timestamp, format='jd', scale='tdb').to_value('iso', subfmt='date')
        if time_format == 'JD_UTC':
            return Time(first_timestamp, format='jd', scale='utc').to_value('iso', subfmt='date')
    except (TypeError, ValueError):
        return None

    return None


def prereduced_obs_date(date, prereduced_file_path=None, time_format=None):
    aavso_obsdate = parse_aavso_obsdate(prereduced_file_path)
    if aavso_obsdate is not None:
        return aavso_obsdate

    derived_obsdate = obs_date_from_first_prereduced_entry(prereduced_file_path, time_format)
    if derived_obsdate is not None:
        return derived_obsdate

    normalized_date = normalize_obs_date(date)
    if normalized_date is not None:
        return normalized_date

    return ""


def data_file_time(time_format):
    while True:
        if not time_format:
            log_info("\nNOTE: If your file is not in one of the following formats, "
                     "\nplease re-reduce your data into one of the time formats recognized by EXOTIC.")

            time_format = user_input("\nWhich of the following time formats is your data file stored in? "
                                     "\nBJD_TDB / JD_UTC / MJD_UTC: ", type_=str)
        time_format = time_format.upper().strip()

        if time_format not in ['BJD_TDB', 'JD_UTC', 'MJD_UTC']:
            log_info("Warning: Invalid entry; please try again.", warn=True)
            time_format = None
        else:
            return time_format


def data_file_units(units):
    while True:
        if not units:
            log_info("\nNOTE: If your file is not in one of the following units, "
                     "\nplease re-reduce your data into one of the units of flux recognized by EXOTIC.")

            units = user_input("\nWhich of the following units of flux is your data file stored in? "
                               "\nflux / magnitude / millimagnitude: ", type_=str)
        units = units.lower().strip()

        if units not in ['flux', 'magnitude', 'millimagnitude']:
            log_info("Warning: Invalid entry; please try again.", warn=True)
            units = None
        else:
            return units


# temp
def log_info(string, warn=False, error=False):
    if error:
        print(f"\033[91m {string}\033[00m")
    elif warn:
        print(f"\033[34m {string}\033[00m")
    else:
        print(string)
    log.debug(string)
    return True
