import logging
import sys
import json
from pathlib import Path
from astropy.io import fits
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


log = logging.getLogger(__name__)
logging.basicConfig(filename='exotic.log', level=logging.DEBUG)
consoleFormatter = logging.Formatter("%(message)s")
consoleHandler = logging.StreamHandler(sys.stdout)
consoleHandler.setFormatter(consoleFormatter)
consoleHandler.setLevel(logging.INFO)
log.addHandler(consoleHandler)


class Inputs:

    def __init__(self, init_opt):
        self.init_opt = init_opt
        self.info_dict = {
            'images': None, 'save': None, 'flats': None, 'darks': None, 'biases': None,
            'aavso_num': None, 'second_obs': None, 'date': None, 'lat': None, 'long': None,
            'elev': None, 'camera': None, 'pixel_bin': None, 'filter': None, 'notes': None,
            'plate_opt': None, 'aavso_comp': None, 'tar_coords': None, 'comp_stars': None,
            'prered_file': None, 'file_units': None, 'file_time': None, 'phot_comp_star': None,
            'wl_min': None, 'wl_max': None, 'pixel_scale': None, 'exposure': None,
            'random_seed': None, "demosaic_fmt": None, "demosaic_out": None
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
        rem_list = ['images', 'plate_opt', 'tar_coords', 'comp_stars']
        [self.params.pop(key) for key in rem_list]

        self.params.update({'exposure': exposure, 'file_units': data_file_units, 'file_time': data_file_time,
                            'phot_comp_star': phot_comp_star})
        self.info_dict['prered_file'] = prereduced_file(self.info_dict['prered_file'])

        if not planet:
            planet = planet_name(planet)

        for key, value in list(self.params.items()):
            if key == 'elev':
                self.info_dict[key] = self.params[key](self.info_dict[key], self.info_dict['lat'],
                                                       self.info_dict['long'])
            else:
                self.info_dict[key] = self.params[key](self.info_dict[key])

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
            'date': 'Observation date', 'lat': 'Obs. Latitude', 'long': 'Obs. Longitude',
            'elev': ('Obs. Elevation (meters)', 'Obs. Elevation (meters; Note: leave blank if unknown)'),
            'camera': 'Camera Type (CCD or DSLR)',
            'pixel_bin': 'Pixel Binning', 'filter': 'Filter Name (aavso.org/filters)',
            'notes': 'Observing Notes', 'plate_opt': 'Plate Solution? (y/n)',
            'aavso_comp': 'Add Comparison Stars from AAVSO? (y/n)',
            'tar_coords': 'Target Star X & Y Pixel', 'comp_stars': 'Comparison Star(s) X & Y Pixel',
        }
        planet_params = {
            'ra': 'Target Star RA', 'dec': 'Target Star Dec', 'pName': "Planet Name", 'sName': "Host Star Name",
            'pPer': 'Orbital Period (days)', 'pPerUnc': 'Orbital Period Uncertainty',
            'midT': 'Published Mid-Transit Time (BJD-UTC)', 'midTUnc': 'Mid-Transit Time Uncertainty',
            'rprs': 'Ratio of Planet to Stellar Radius (Rp/Rs)',
            'rprsUnc': 'Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty',
            'aRs': 'Ratio of Distance to Stellar Radius (a/Rs)',
            'aRsUnc': 'Ratio of Distance to Stellar Radius (a/Rs) Uncertainty',
            'inc': 'Orbital Inclination (deg)',
            'incUnc': ('Orbital Inclination (deg) Uncertainty', 'Orbital Inclination (deg) Uncertainity'),
            'ecc': 'Orbital Eccentricity (0 if null)', 'teff': 'Star Effective Temperature (K)',
            'omega': 'Argument of Periastron (deg)',
            'teffUncPos': 'Star Effective Temperature (+) Uncertainty',
            'teffUncNeg': 'Star Effective Temperature (-) Uncertainty',
            'met': 'Star Metallicity ([FE/H])', 'metUncPos': 'Star Metallicity (+) Uncertainty',
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
            'phot_comp_star': "Comparison Star used in Photometry (leave blank if none)",
            'wl_min': 'Filter Minimum Wavelength (nm)', 'wl_max': 'Filter Maximum Wavelength (nm)',
            'pixel_scale': ('Image Scale (Ex: 5.21 arcsecs/pixel)', 'Pixel Scale (Ex: 5.21 arcsecs/pixel)'),
            'exposure': 'Exposure Time (s)',
            'random_seed': 'Random Seed'
        }

        self.info_dict = init_params(user_info, self.info_dict, data['user_info'])
        self.info_dict = init_params(opt_info, self.info_dict, data['optional_info'])
        return init_params(planet_params, planet_dict, data['planetary_parameters'])


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


def latitude(lat, hdr=None):
    while True:
        if not lat:
            if hdr:
                lat = find(hdr, ['LATITUDE', 'LAT', 'SITELAT'])
                if lat:
                    return lat
            lat = user_input("Enter the latitude (in degrees) of where you observed. "
                             "(Don't forget the sign where North is '+' and South is '-')! "
                             "(Example: -32.12): ", type_=str)
        lat = lat.strip()

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


def longitude(long, hdr=None):
    while True:
        if not long:
            if hdr:
                long = find(hdr, ['LONGITUD', 'LONG', 'LONGITUDE', 'SITELONG'])
                if long:
                    return long
            long = user_input("Enter the longitude (in degrees) of where you observed. "
                              "(Don't forget the sign where East is '+' and West is '-')! "
                              "(Example: +152.51): ", type_=str)
        long = long.strip()

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


def elevation(elev, lat, long, hdr=None):
    while True:
        try:
            elev = typecast_check(type_=float, val=elev)
            if not elev:
                if hdr:
                    elev = find(hdr, ['HEIGHT', 'ELEVATION', 'ELE', 'EL', 'OBSGEO-H', 'ALT-OBS', 'SITEELEV'])
                    if elev:
                        return int(elev)
                log_info("\nEXOTIC is retrieving elevation based on entered "
                         "latitude and longitude from Open Elevation.")
                animate_toggle(True)
                elev = open_elevation(lat, long)
                animate_toggle()
                if not elev:
                    log_info("\nWarning: EXOTIC could not retrieve elevation.", warn=True)
                    elev = user_input("Enter the elevation (in meters) of where you observed: ", type_=float)
            return elev
        except ValueError:
            log_info("Error: The entered elevation is incorrect.", error=True)
            elev = None


def camera(c_type):
    while True:
        if not c_type:
            c_type = user_input("\nPlease enter the camera type (e.g., CCD or DSLR;\n"
                                "Note: if you are using a CMOS, please enter CCD here and\n"
                                "then note your actual camera type in \"Observing Notes\"): ", type_=str)
        c_type = c_type.strip().upper()
        if c_type not in ["CCD", "DSLR"]:
            c_type = None
        else:
            return c_type


def pixel_bin(pix_bin):
    if not pix_bin:
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
    if isinstance(comp_stars, list) and 1 <= len(comp_stars) <= 10 and \
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
                num_comp_stars = user_input("\nHow many Comparison Stars would you like to use? (1-10): ", type_=int)
                if 1 <= num_comp_stars <= 10:
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


def phot_comp_star(comp_star):
    if not isinstance(comp_star, dict):
        comp_star_opt = user_input("Was a Comparison Star used during Photometry? (y/n): ",
                                   type_=str, values=['y', 'n'])

        comp_star = {
            'ra': user_input("\nEnter Comparison Star RA: ", type_=str) if comp_star_opt == 'y' else '',
            'dec': user_input("Enter Comparison Star DEC: ", type_=str) if comp_star_opt == 'y' else '',
            'x': user_input("\nEnter Comparison Star X Pixel Coordinate: ", type_=str) if comp_star_opt == 'y' else '',
            'y': user_input("Enter Comparison Star Y Pixel Coordinate: ", type_=str) if comp_star_opt == 'y' else ''
        }

    return comp_star


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
        print(f"\033[93m {string}\033[00m")
    else:
        print(string)
    log.debug(string)
    return True
