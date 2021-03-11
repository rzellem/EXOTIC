import logging
import json
import requests
from pathlib import Path
from datetime import datetime
from tenacity import retry, retry_if_exception_type, retry_if_result, \
    stop_after_attempt, wait_exponential
try:
    from plate_solution import is_false, result_if_max_retry_count
except ImportError:
    from .plate_solution import is_false, result_if_max_retry_count

log = logging.getLogger(__name__)


def user_input(prompt, type_, val1=None, val2=None, val3=None):
    while True:
        try:
            result = type_(input(prompt))
            log.debug(f"{prompt}{result}")
        except ValueError:
            log.info("Sorry, not a valid datatype.")
            continue
        if type_ == str and val1 and val2 and val3:
            result = result.lower().strip()
            if result not in (val1, val2, val3):
                log.info("Sorry, your response was not valid.")
            else:
                return result
        elif type_ == int and val1 and val2 and val3:
            if result not in (val1, val2, val3):
                log.info("Sorry, your response was not valid.")
            else:
                return result
        elif type_ == int and val1 and val2 and val3:
            if result not in (val1, val2, val3):
                log.info("Sorry, your response was not valid.")
            else:
                return result
        else:
            return result


class UserInputs:

    def __init__(self, init_file, init_opt):
        self.init_file = init_file
        self.init_opt = init_opt
        self.info_dict = {
            'images': None, 'save': None, 'flats': None, 'darks': None, 'biases': None,
            'aavso_num': None, 'second_obs': None, 'date': None, 'lat': None, 'long': None,
            'elev': None, 'camera': None, 'pixel_bin': None, 'filter': None, 'notes': None,
            'plate_opt': None, 'tar_coords': None, 'comp_stars': None, 'pixel_scale': None,
            'exposure': None, 'wl_min': None, 'wl_max': None
        }
        self.planet_dict = {
            'ra': None, 'dec': None, 'p_name': None, 's_name': None, 'per': None, 'per_unc': None,
            'mid_t': None, 'mid_t_unc': None, 'rprs': None, 'rprs_unc': None,
            'ars': None, 'ars_unc': None, 'inc': None, 'inc_unc': None, 'ecc': None,
            'teff': None, 'teff_unc_pos': None, 'teff_unc_neg': None,
            'met': None, 'met_unc_pos': None, 'met_unc_neg': None,
            'logg': None, 'logg_unc_pos': None, 'logg_unc_neg': None
        }

    def complete_red(self):
        params = {
            'images': imaging_files, 'save': save_directory, 'flats': flats, 'darks': darks, 'biases': biases,
            'aavso_num': obs_code, 'second_obs': second_obs_code, 'date': obs_date, 'lat': latitude,
            'long': longitude, 'elev': elevation, 'camera': camera, 'pixel_bin': pixel_bin, 'filter': filter_type,
            'notes': obs_notes, 'plate_opt': plate_solution_opt, 'tar_coords': target_star_coords,
            'comp_stars': comparison_star_coords
        }

        if self.init_opt == 'y':
            self.search_init()

        image_calibrations(self.info_dict['flats'], self.info_dict['darks'], self.info_dict['biases'])

        for key, value in list(params.items()):
            if key == 'elev':
                self.info_dict[key] = params[key](self.info_dict[key], self.info_dict['lat'], self.info_dict['long'])
            elif key == 'tar_coords':
                self.info_dict[key] = params[key](self.info_dict[key], self.planet_dict['p_name'])
            else:
                self.info_dict[key] = params[key](self.info_dict[key])

        return self.info_dict, self.planet_dict

    def search_init(self):
        cwd = Path.cwd()
        log.info(f"\nYour current working directory is: {cwd}")
        log.info(f"Potential initialization files I've found in {cwd} are: ")
        [log.info(f"\t{file}") for file in cwd.glob('*.json') if file.is_file()]

        while True:
            try:
                if not self.init_file:
                    self.init_file = user_input("\nPlease enter the Directory and Filename of "
                                                "your Initialization File: ", type_=str)
                if self.init_file == 'ok':
                    self.init_file = '/Users/rzellem/Documents/EXOTIC/inits.json'
                self.init_file = Path(self.init_file)
                self.comp_params()
                break
            except FileNotFoundError:
                log.info("Error: Initialization file not found. Please try again.")
                self.init_file = None
            except IsADirectoryError:
                log.info("Error: Entered a directory. Please try again.")
                self.init_file = None

    def comp_params(self):
        with self.init_file.open('r') as json_file:
            data = json.load(json_file)

        comp_info = {
            'images': 'Directory with FITS files', 'save': 'Directory to Save Plots', 'flats': 'Directory of Flats',
            'darks': 'Directory of Darks', 'biases': 'Directory of Biases',
            'aavso_num': 'AAVSO Observer Code (N/A if none)', 'second_obs': 'Secondary Observer Codes (N/A if none)',
            'date': 'Observation date', 'lat': 'Obs. Latitude', 'long': 'Obs. Longitude',
            'elev': 'Obs. Elevation (meters)', 'camera': 'Camera Type (CCD or DSLR)', 'pixel_bin': 'Pixel Binning',
            'filter': 'Filter Name (aavso.org/filters)', 'wl_min': 'Filter Minimum Wavelength (nm)',
            'wl_max': 'Filter Maximum Wavelength (nm)', 'notes': 'Observing Notes',
            'plate_opt': 'Plate Solution? (y/n)', 'tar_coords': 'Target Star X & Y Pixel',
            'comp_stars': 'Comparison Star(s) X & Y Pixel', 'pixel_scale': 'Pixel Scale (Ex: 5.21 arcsecs/pixel)'
        }
        comp_planet = {
            'ra': 'Target Star RA', 'dec': 'Target Star Dec', 'p_name': "Planet Name", 's_name': "Host Star Name",
            'per': 'Orbital Period (days)', 'per_unc': 'Orbital Period Uncertainty',
            'mid_t': 'Published Mid-Transit Time (BJD-UTC)', 'mid_t_unc': 'Mid-Transit Time Uncertainty',
            'rprs': 'Ratio of Planet to Stellar Radius (Rp/Rs)',
            'rprs_unc': 'Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty',
            'ars': 'Ratio of Distance to Stellar Radius (a/Rs)',
            'ars_unc': 'Ratio of Distance to Stellar Radius (a/Rs) Uncertainty',
            'inc': 'Orbital Inclination (deg)', 'inc_unc': 'Orbital Inclination (deg) Uncertainity',
            'ecc': 'Orbital Eccentricity (0 if null)', 'teff': 'Star Effective Temperature (K)',
            'teff_unc_pos': 'Star Effective Temperature (+) Uncertainty',
            'teff_unc_neg': 'Star Effective Temperature (-) Uncertainty', 'met': 'Star Metallicity ([FE/H])',
            'met_unc_pos': 'Star Metallicity (+) Uncertainty', 'met_unc_neg': 'Star Metallicity (-) Uncertainty',
            'logg': 'Star Surface Gravity (log(g))',
            'logg_unc_pos': 'Star Surface Gravity (+) Uncertainty',
            'logg_unc_neg': 'Star Surface Gravity (-) Uncertainty'
        }

        self.info_dict = init_params(comp_info, self.info_dict, data['user_info'])
        self.info_dict = init_params(comp_info, self.info_dict, data['optional_info'])
        self.planet_dict = init_params(comp_planet, self.planet_dict, data['planetary_parameters'])


def init_params(comp, dict1, dict2):
    for key, value in comp.items():
        try:
            dict1[key] = dict2[value]
        except KeyError:
            pass
    return dict1


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
                        return directory, input_files
                if not input_files:
                    raise FileNotFoundError
            else:
                raise NotADirectoryError
        except FileNotFoundError:
            opt = user_input(
                f"\nError: {img_type} files not found with .fits, .fit, .fts, or .fz extensions in {directory}."
                "\nWould you like to enter in an alternate image extension in addition to .FITS? (y/n): ",
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


def image_calibrations(flats_dir, darks_dir, biases_dir):
    opt, flats_list, darks_list, biases_list = None, None, None, None

    if [x for x in (flats_dir, darks_dir, biases_dir) if x is None]:
        opt = user_input("\nDo you have any Calibration Images? (Flats, Darks or Biases)? (y/n): ",
                         type_=str, val1='y', val2='n')

    if opt == 'y' or flats_dir:
        flats_dir, flats_list = flats(flats_dir)
    if opt == 'y' or darks_dir:
        darks_dir, darks_list = darks(darks_dir)
    if opt == 'y' or biases_dir:
        biases_dir, biases_list = biases(biases_dir)

    return flats_dir, flats_list, darks_dir, darks_list, biases_dir, biases_list


def flats(directory):
    if directory is None:
        opt = user_input("\nDo you have Flats? (y/n): ", type_=str, val1='y', val2='n')
        if opt == 'y':
            directory = user_input("Please enter the directory path to your Flats "
                                   "(must be in their own separate folder): ", type_=str)
    if directory:
        return check_imaging_files(directory, 'Flats')
    return directory, None


def darks(directory):
    if directory is None:
        opt = user_input("\nDo you have Darks? (y/n): ", type_=str, val1='y', val2='n')
        if opt == 'y':
            directory = user_input("Please enter the directory path to your Darks "
                                   "(must be in their own separate folder): ", type_=str)
    if directory:
        return check_imaging_files(directory, 'Darks')
    return directory, None


def biases(directory):
    if directory is None:
        opt = user_input("\nDo you have Biases? (y/n): ", type_=str, val1='y', val2='n')
        if opt == 'y':
            directory = user_input("Please enter the directory path to your Biases "
                                   "(must be in their own separate folder): ", type_=str)
    if directory:
        return check_imaging_files(directory, 'Biases')
    return directory, None


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
        try:
            if not date:
                date = user_input("\nPlease enter the Observation Date (DD-Month-YYYY): ", type_=str)
            if date != datetime.strptime(date, '%d-%B-%Y').strftime('%d-%B-%Y'):
                raise ValueError
            return date
        except ValueError:
            log.info('\nThe entered Observation Date format is incorrect.')
            date = None


def dms_to_dd(dms_in):
    """
    Quick helper method to convert long/lat values in degree-minute-second (dms) form
    (using ':' separators) to decimal (dd) form
    :param dms_in: DMS long/lat value, colon separated
    :return float: Properly signed long/lat value in decimal float form
    """
    if dms_in is None or isinstance(dms_in, str) is False or str(dms_in).count(":") != 2:
        raise ValueError("Invalid DMS input provided for calculations. ...")
    # clean string of errant leading/trailing/internal spaces
    dms = str(dms_in).strip().replace(" ", "")
    degrees, minutes, seconds = dms.split(":")
    dec = abs(float(degrees)) + float(minutes) / 60. + float(seconds) / 3600.
    if float(degrees) < 0.:
        dec = dec * -1.
    return dec


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
            if not elev:
                log.info("\nEXOTIC is retrieving elevation based on entered "
                         "latitude and longitude from Open Elevation.")
                elev = open_elevation(lat, long)
                if not elev:
                    log.info("\nEXOTIC could not retrieve elevation.")
                    elev = user_input("Enter the elevation (in meters) of where you observed: ", type_=float)
            return float(elev)
        except ValueError:
            log.info("The entered elevation is incorrect.")
            elev = None


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
       retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
       retry_error_callback=result_if_max_retry_count)
def open_elevation(lat, long):
    query = f"https://api.open-elevation.com/api/v1/lookup?locations={lat},{long}"
    try:
        r = requests.get(query).json()
        return r['results'][0]['elevation']
    except requests.exceptions.RequestException:
        return False


def camera(c_type):
    if not c_type:
        c_type = user_input("Please enter the camera type (CCD or DSLR): ", type_=str)
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
        notes = user_input("Please enter any observing notes (seeing, weather, etc.): ", type_=str)
    if not notes.replace(' ', ''):
        notes = "na"
    return notes


def plate_solution_opt(opt):
    if opt not in ('y', 'n'):
        opt = user_input("\nWould you like to upload the your image for a plate solution?"
                         "\nDISCLAIMER: One of your imaging files will be publicly viewable on "
                         "nova.astrometry.net. (y/n): ", type_=str, val1='y', val2='n')
    return opt


def target_star_coords(coords, planet):
    if isinstance(coords, list) and len(coords) == 2:
        pass
    else:
        coords = [user_input(f"\nPlease enter {planet}'s X Pixel Coordinate: ", type_=float),
                  user_input(f"\nPlease enter {planet}'s Y Pixel Coordinate: ", type_=float)]

    return coords


def comparison_star_coords(comp_stars):
    if isinstance(comp_stars, list) and len(comp_stars) >= 1 and all(isinstance(star, list) for star in comp_stars):
        pass
    else:
        while True:
            num_comp_stars = user_input("How many Comparison Stars would you like to use? (1-10): ", type_=int)
            if 1 <= num_comp_stars <= 10:
                break
            log.info("\nThe number of Comparison Stars entered is incorrect.")

        for num in range(num_comp_stars):
            x_pix = user_input(f"\nComparison Star {num + 1} X Pixel Coordinate: ", type_=float)
            y_pix = user_input(f"Comparison Star {num + 1} Y Pixel Coordinate: ", type_=float)
            comp_stars.append((x_pix, y_pix))

    final_comp_stars = [star for star in comp_stars if star != []]

    return final_comp_stars


def exposure(exp):
    if not exp:
        exp = user_input("Please enter your exposure time (seconds): ", type_=int)
    return exp


def prereduced_file():
    while True:
        try:
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


if __name__ == "__main__":
    obj = UserInputs('/Users/abdullahfatahi/Documents/ExoplanetWatch/EXOTIC/inits.json', 'y')
    info, planet = obj.complete_red()
    print(info, planet)
