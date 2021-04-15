import logging
import json
from pathlib import Path
try:
    from util import user_input, dms_to_dd, open_elevation, typecast_check, init_params
except ImportError:
    from .util import user_input, dms_to_dd, open_elevation, typecast_check, init_params


log = logging.getLogger(__name__)


class Inputs:

    def __init__(self, init_opt):
        self.init_opt = init_opt
        self.info_dict = {
            'images': None, 'save': None, 'flats': None, 'darks': None, 'biases': None,
            'aavso_num': None, 'second_obs': None, 'date': None, 'lat': None, 'long': None,
            'elev': None, 'camera': None, 'pixel_bin': None, 'filter': None, 'notes': None,
            'plate_opt': None, 'img_align_opt': None, 'tar_coords': None, 'comp_stars': None,
            'prered_file': None, 'file_units': None, 'file_time': None,
            'wl_min': None, 'wl_max': None, 'pixel_scale': None, 'exposure': None
        }
        self.params = {
            'images': imaging_files, 'save': save_directory, 'aavso_num': obs_code, 'second_obs': second_obs_code,
            'date': obs_date, 'lat': latitude, 'long': longitude, 'elev': elevation, 'camera': camera,
            'pixel_bin': pixel_bin, 'filter': filter_type, 'notes': obs_notes, 'plate_opt': plate_solution_opt,
            'img_align_opt': image_align_opt, 'tar_coords': target_star_coords, 'comp_stars': comparison_star_coords
        }

    def complete_red(self, planet):
        for key, value in list(self.params.items()):
            if key == 'elev':
                self.info_dict[key] = self.params[key](self.info_dict[key], self.info_dict['lat'],
                                                       self.info_dict['long'])
            elif key == 'tar_coords':
                self.info_dict[key] = self.params[key](self.info_dict[key], planet)
            elif key == 'comp_stars':
                self.info_dict[key] = self.params[key](self.info_dict[key], False)
            else:
                self.info_dict[key] = self.params[key](self.info_dict[key])
            if key == 'save':
                self.info_dict['flats'], self.info_dict['darks'], self.info_dict['biases'] = \
                    image_calibrations(self.info_dict['flats'], self.info_dict['darks'],
                                       self.info_dict['biases'], self.init_opt)

        return self.info_dict

    def prereduced(self):
        rem_list = ['images', 'plate_opt', 'img_align_opt', 'tar_coords', 'comp_stars']
        [self.params.pop(key) for key in rem_list]

        self.params.update({'exposure': exposure, 'file_units': data_file_units, 'file_time': data_file_time})
        self.info_dict['prered_file'] = prereduced_file(self.info_dict['prered_file'])

        for key, value in list(self.params.items()):
            if key == 'elev':
                self.info_dict[key] = self.params[key](self.info_dict[key], self.info_dict['lat'],
                                                       self.info_dict['long'])
            else:
                self.info_dict[key] = self.params[key](self.info_dict[key])

        return self.info_dict

    def real_time(self):
        rem_list = ['save', 'aavso_num', 'second_obs', 'date', 'lat', 'long', 'elev',
                    'camera', 'pixel_bin', 'filter', 'notes', 'plate_opt', 'img_align_opt']
        [self.params.pop(key) for key in rem_list]

        for key, value in list(self.params.items()):
            if key == 'comp_stars':
                self.info_dict[key] = self.params[key](self.info_dict[key], True)
            else:
                self.info_dict[key] = self.params[key](self.info_dict[key])

        return self.info_dict

    def search_init(self, init_file, planet_dict):
        cwd = Path.cwd()
        log.info(f"\nYour current working directory is: {cwd}")
        log.info(f"Potential initialization files I've found in {cwd} are: ")
        [log.info(f"\t{file}") for file in cwd.glob('*.json') if file.is_file()]

        while True:
            try:
                if not init_file:
                    init_file = user_input("\nPlease enter the Directory and Filename of "
                                           "your Initialization File: ", type_=str)
                if init_file == 'ok':
                    init_file = '/Users/rzellem/Documents/EXOTIC/inits.json'
                init_file = Path(init_file)
                planet_params = self.comp_params(init_file, planet_dict)
                return init_file, planet_params
            except (FileNotFoundError, IsADirectoryError) as e:
                log.info(f"Error: Initialization file not found. \n{e}. \nPlease try again.")
                init_file = None

    def comp_params(self, init_file, planet_dict):
        with init_file.open('r') as json_file:
            data = json.load(json_file)

        user_info = {
            'images': 'Directory with FITS files', 'save': 'Directory to Save Plots',
            'flats': 'Directory of Flats', 'darks': 'Directory of Darks', 'biases': 'Directory of Biases',
            'aavso_num': 'AAVSO Observer Code (N/A if none)', 'second_obs': 'Secondary Observer Codes (N/A if none)',
            'date': 'Observation date', 'lat': 'Obs. Latitude', 'long': 'Obs. Longitude',
            'elev': 'Obs. Elevation (meters)', 'camera': 'Camera Type (CCD or DSLR)', 'pixel_bin': 'Pixel Binning',
            'filter': 'Filter Name (aavso.org/filters)', 'notes': 'Observing Notes',
            'plate_opt': 'Plate Solution? (y/n)', 'img_align_opt': 'Align Images? (y/n)',
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
            'inc': 'Orbital Inclination (deg)', 'incUnc': 'Orbital Inclination (deg) Uncertainty',
            'ecc': 'Orbital Eccentricity (0 if null)', 'teff': 'Star Effective Temperature (K)',
            'teffUncPos': 'Star Effective Temperature (+) Uncertainty',
            'teffUncNeg': 'Star Effective Temperature (-) Uncertainty',
            'met': 'Star Metallicity ([FE/H])', 'metUncPos': 'Star Metallicity (+) Uncertainty',
            'metUncNeg': 'Star Metallicity (-) Uncertainty',
            'logg': 'Star Surface Gravity (log(g))', 'loggUncPos': 'Star Surface Gravity (+) Uncertainty',
            'loggUncNeg': 'Star Surface Gravity (-) Uncertainty'
        }
        opt_info = {
            'prered_file': 'Pre-reduced File:', 'file_time': 'Pre-reduced File Time Format (BJD_TDB, JD_UTC, MJD_UTC)',
            'file_units': 'Pre-reduced File Units of Flux (flux, magnitude, millimagnitude)',
            'wl_min': 'Filter Minimum Wavelength (nm)', 'wl_max': 'Filter Maximum Wavelength (nm)',
            'pixel_scale': 'Pixel Scale (Ex: 5.21 arcsecs/pixel)', 'exposure': 'Exposure Time (s)',
        }

        self.info_dict = init_params(user_info, self.info_dict, data['user_info'])
        self.info_dict = init_params(opt_info, self.info_dict, data['optional_info'])
        return init_params(planet_params, planet_dict, data['planetary_parameters'])


def check_imaging_files(directory, img_type):
    file_extensions = ['.fits', '.fit', '.fts', '.fz']
    input_files = []

    while True:
        try:
            if Path(directory).is_dir():
                directory = Path(directory)
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
            log.info(f"\nError: {img_type} files not found with .fits, .fit, .fts, or .fz extensions in {directory}.")
            opt = user_input("\nWould you like to enter in an alternate image extension in addition to .FITS? (y/n): ",
                             type_=str, val1='y', val2='n')
            if opt == 'y':
                add_ext = user_input("Please enter the extension you want to add (EX: .FITS): ", type_=str)
                file_extensions.append(add_ext)
            else:
                directory = user_input(f"Enter the directory path where {img_type} files are located "
                                       f"(Example using the sample data: sample-data/HatP32Dec202017): ", type_=str)
        except (NotADirectoryError, OSError):
            log.info("\nError: No such directory exists when searching for FITS files. Please try again.")
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
                directory = user_input("Enter the directory to Save the Results and Plots into "
                                       "or type new to create one: ", type_=str)
            if directory == 'new':
                directory = create_directory()
            else:
                if not Path(directory).is_dir():
                    raise NotADirectoryError
            return directory
        except (NotADirectoryError, OSError):
            log.info("Error: the directory entered does not exist. Please try again. Make sure to follow this "
                     "\nformatting (using whichever directory you choose): /sample-data/results")
            directory = None


def create_directory():
    save_path = Path.cwd()
    while True:
        directory = user_input("Enter the name for your new directory: ", type_=str)
        try:
            save_path = save_path / directory
            Path(save_path).mkdir()
        except OSError:
            log.info(f"Creation of the directory {save_path}/{directory} failed.")
        else:
            log.info(f"Successfully created the directory {save_path}.")
            return save_path


def image_calibrations(flats_dir, darks_dir, biases_dir, init):
    opt, flats_list, darks_list, biases_list = None, None, None, None

    if init == 'n':
        opt = user_input("\nDo you have any Calibration Images? (Flats, Darks or Biases)? (y/n): ",
                         type_=str, val1='y', val2='n')

    if opt == 'y' or flats_dir:
        flats_list = check_calibration(flats_dir, 'Flats')
    if opt == 'y' or darks_dir:
        darks_list = check_calibration(darks_dir, 'Darks')
    if opt == 'y' or biases_dir:
        biases_list = check_calibration(biases_dir, 'Biases')

    return flats_list, darks_list, biases_list


def check_calibration(directory, image_type):
    if not directory:
        opt = user_input(f"\nDo you have {image_type}? (y/n): ", type_=str, val1='y', val2='n')
        if opt == 'y':
            directory = user_input("Please enter the directory path to your Flats "
                                   "(must be in their own separate folder): ", type_=str)
    if directory:
        return check_imaging_files(directory, image_type)
    return None


def obs_code(code):
    if code is None:
        code = user_input("Please enter your AAVSO Observer Account Number "
                          "(if none, leave blank and press enter): ", type_=str)
    if not code.replace(' ', ''):
        code = "N/A"
    return code


def second_obs_code(code):
    if code is None:
        code = user_input("Please enter your comma-separated Secondary Observer Codes "
                          "(if none, leave blank and press enter): ", type_=str)
    if not code.replace(' ', ''):
        code = "N/A"
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


def latitude(lat):
    while True:
        try:
            if not lat:
                lat = user_input("Enter the latitude (in degrees) of where you observed. "
                                 "Don't forget the sign where North is '+' and South is '-'! "
                                 "(Example: +50.4): ", type_=str)

            lat = lat.replace(' ', '')
            if lat[0] != '+' and lat[0] != '-':
                raise ValueError("You forgot the sign for the latitude! North is '+' and South is '-'. "
                                 "Please try again.")

            # Convert to float if latitude in decimal. If latitude is in +/-HH:MM:SS format, convert to a float.
            try:
                lat = float(lat)
            except ValueError:
                lat = float(dms_to_dd(lat))

            if lat <= -90.00 or lat >= 90.00:
                raise ValueError("Your latitude is out of range. Please enter a latitude between -90 and +90 (deg)")
            return lat
        except ValueError as err:
            log.info(err.args)
            lat = None


def longitude(long):
    while True:
        try:
            if not long:
                long = user_input("Enter the longitude (in degrees) of where you observed. "
                                  "(Don't forget the sign where East is '+' and West is '-')! "
                                  "(Example: -32.12): ", type_=str)

            long = long.replace(' ', '')
            if long[0] != '+' and long[0] != '-':
                raise ValueError("You forgot the sign for the longitude! East is '+' and West is '-'. "
                                 "Please try again.")

            # Convert to float if longitude in decimal. If longitude is in +/-HH:MM:SS format, convert to a float.
            try:
                long = float(long)
            except ValueError:
                long = float(dms_to_dd(long))

            if long <= -180.00 or long >= 180.00:
                raise ValueError("Your longitude is out of range. Please enter a longitude between -180 and +180 (deg)")
            return long
        except ValueError as err:
            log.info(err.args)
            long = None


def elevation(elev, lat, long):
    while True:
        try:
            elev = typecast_check(type_=float, val=elev)
            if not elev:
                log.info("\nEXOTIC is retrieving elevation based on entered "
                         "latitude and longitude from Open Elevation.")
                elev = open_elevation(lat, long)
                if not elev:
                    log.info("\nEXOTIC could not retrieve elevation.")
                    elev = user_input("Enter the elevation (in meters) of where you observed: ", type_=float)
            return elev
        except ValueError:
            log.info("The entered elevation is incorrect.")
            elev = None


def camera(c_type):
    if not c_type:
        c_type = user_input("\nPlease enter the camera type (CCD or DSLR): ", type_=str)
    return c_type


def pixel_bin(pix_bin):
    if not pix_bin:
        pix_bin = user_input("Please enter the pixel binning: ", type_=str)
    return pix_bin


def filter_type(f_type):
    if not f_type:
        f_type = user_input("Please enter the filter name: ", type_=str)
    return f_type


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
                         "\nDISCLAIMER: One of your imaging files will be publicly viewable on "
                         "nova.astrometry.net. (y/n): ", type_=str, val1='y', val2='n')
    return opt


def image_align_opt(opt):
    if opt:
        opt = opt.lower().strip()
    if opt not in ('y', 'n'):
        opt = user_input("\nWould you like to align your images (y/n): ", type_=str, val1='y', val2='n')
    return opt


def target_star_coords(coords, planet):
    if isinstance(coords, list) and len(coords) == 2:
        pass
    else:
        coords = [user_input(f"\nPlease enter {planet}'s X Pixel Coordinate: ", type_=float),
                  user_input(f"\nPlease enter {planet}'s Y Pixel Coordinate: ", type_=float)]

    return coords


def comparison_star_coords(comp_stars, rt_bool):
    if isinstance(comp_stars, list) and 1 <= len(comp_stars) <= 10 and \
            all(isinstance(star, list) for star in comp_stars):
        comp_stars = [star for star in comp_stars if star != []]
    else:
        comp_stars = []

    if not comp_stars:
        while True:
            if not rt_bool:
                num_comp_stars = user_input("\nHow many Comparison Stars would you like to use? (1-10): ", type_=int)
                if 1 <= num_comp_stars <= 10:
                    break
                log.info("\nThe number of Comparison Stars entered is incorrect.")
            else:
                num_comp_stars = 1

        for num in range(num_comp_stars):
            x_pix = user_input(f"\nComparison Star {num + 1} X Pixel Coordinate: ", type_=float)
            y_pix = user_input(f"Comparison Star {num + 1} Y Pixel Coordinate: ", type_=float)
            comp_stars.append([x_pix, y_pix])

    if rt_bool and isinstance(comp_stars[1], list):
        comp_stars = comp_stars[1]

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
                log.info("Hello, Rob.")

            file = Path(file)

            if file.is_file():
                return file
            else:
                raise FileNotFoundError
        except FileNotFoundError:
            log.info("Error: Data file not found. Please try again.")
            file = None


def data_file_time(time_format):
    log.info("\nNOTE: If your file is not in one of the following formats, "
             "\nplease re-reduce your data into one of the time formats recognized by EXOTIC.")

    while True:
        if not time_format:
            time_format = user_input("\nWhich of the following time formats is your data file stored in? "
                                     "\nBJD_TDB / JD_UTC / MJD_UTC: ", type_=str)
        time_format = time_format.upper().strip()

        if time_format not in ['BJD_TDB', 'JD_UTC', 'MJD_UTC']:
            log.info("Invalid entry; please try again.")
            time_format = None
        else:
            return time_format


def data_file_units(units):
    log.info("\nNOTE: If your file is not in one of the following units, "
             "\nplease re-reduce your data into one of the units of flux recognized by EXOTIC.")

    while True:
        if not units:
            units = user_input("\nWhich of the following units of flux is your data file stored in? "
                               "\nflux / magnitude / millimagnitude: ", type_=str)
        units = units.lower().strip()

        if units not in ['flux', 'magnitude', 'millimagnitude']:
            log.info("Invalid entry; please try again.")
            units = None
        else:
            return units
