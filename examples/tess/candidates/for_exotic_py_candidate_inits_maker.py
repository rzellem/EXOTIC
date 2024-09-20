import requests
import os
import time
import json
import csv
from datetime import datetime
from math import degrees, floor
import pandas as pd
import numpy as np
from astropy.constants import G, M_sun, R_sun, au
from astropy import units as u

def download_tic_json(tic_id):
    url = f'https://exofop.ipac.caltech.edu/tess/target.php?id={tic_id}&json'
    response = requests.get(url)
    
    if response.status_code == 200:
        if 'No TIC object found' in response.text:
            print(f'No TIC object found for TIC ID: {tic_id}')
            return False
        else:
            # Create the directory using TIC ID if it doesn't exist
            directory = f'output_inits_files/for_exotic_py_candidate_inits_output_files/{tic_id}_files'
            os.makedirs(directory, exist_ok=True)
            
            # Write to a JSON file
            with open(f'{directory}/{tic_id}.json', 'wb') as file:
                file.write(response.content)
            
            print(f'Downloaded JSON for TIC ID: {tic_id}')
            return True
    else:
        print(f'Failed to download JSON for TIC ID: {tic_id}')
        return False

# Introduction message
print("*****************************************************************************************************")
print("This script will scrape exofop for available exoplanets with their associated parameter entries and make an inits json file for exotic.py. You will be able to choose which ExoFOP entry you would like to make the inits file from. Some parameters are not available for exoplanet entries. You will be given the option to have them estimated using available parameters or to be left as null to fill out later. The process will give you the inits file for entry chosen, original ExoFOP json file and a csv spreadsheet with system parameters for all exoplanet entries.")
print("*****************************************************************************************************")

# Loop until a valid TIC ID is entered and successfully downloaded
while True:
    tic_id = input("Enter the TIC ID: ")

    # Ensure the input is an integer
    try:
        tic_id = int(tic_id)
        if download_tic_json(tic_id):
            break  # Exit the loop if the download was successful
    except ValueError:
        print("Invalid TIC ID. Please enter a correct numeric value.")

# for multiple file version, Delay for a specified number of seconds so we're not a robot. beep boop
time.sleep(1)  

print(f'Download completed and saved to output_inits_files/for_exotic_py_candidate_inits_output_files/{tic_id}_files folder')

def load_json_data(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()
            if not content.strip():  # Check if the file is empty
                print(f"Skipping empty file: {file_path}")
                return None
            return json.loads(content)
    except UnicodeDecodeError:
        with open(file_path, 'r', encoding='latin-1') as file:
            content = file.read()
            if not content.strip():  # Check if the file is empty
                print(f"Skipping empty file: {file_path}")
                return None
            return json.loads(content)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error loading JSON data from {file_path}: {e}")
        return None

def ra_to_hms(ra):
    hours = ra * 24 / 360
    h = floor(hours)
    m = floor((hours - h) * 60)
    s = (hours - h - m/60) * 3600
    return f"{h:02d}:{m:02d}:{s:.2f}"

def dec_to_dms(dec):
    sign = '+' if dec >= 0 else '-'
    dec = abs(dec)
    d = floor(dec)
    m = floor((dec - d) * 60)
    s = (dec - d - m/60) * 3600
    return f"{sign}{d:02d}:{m:02d}:{s:.2f}"

def get_stellar_parameter(stellar_params, param):
    for entry in stellar_params:
        if entry.get(param) is not None:
            return entry.get(param)
    return 'Not available'

def extract_data(data):
    tic_id = data.get('basic_info', {}).get('tic_id', 'Not available')
    ra = data.get('coordinates', {}).get('ra', 'Not available')
    dec = data.get('coordinates', {}).get('dec', 'Not available')
    ra_hms = ra_to_hms(float(ra)) if ra != 'Not available' else 'Not available'
    dec_dms = dec_to_dms(float(dec)) if dec != 'Not available' else 'Not available'
    ra_deg = degrees(float(ra)) if ra != 'Not available' else 'Not available'
    dec_deg = degrees(float(dec)) if dec != 'Not available' else 'Not available'

    planet_data = []
    for planet in data.get('planet_parameters', []):
        name_parts = planet.get('name', '').split('.')
        planet_tic_id = f"{tic_id}.{name_parts[-1]}" if len(name_parts) > 1 else tic_id
        toi = planet.get('toi', '')

        planet_data.append({
            'tic_id': planet_tic_id,
            'toi': toi,
            'epoch': planet.get('epoch', 'Not available'),
            'epoch_e': planet.get('epoch_e', 'Not available'),
            'per': planet.get('per', 'Not available'),
            'per_e': planet.get('per_e', 'Not available'),
            'dur': planet.get('dur', 'Not available'),
            'dur_e': planet.get('dur_e', 'Not available'),
            'dep_p': planet.get('dep_p', 'Not available'),
            'dep_p_e': planet.get('dep_p_e', 'Not available'),
            'rad': planet.get('rad', 'Not available'),
            'rad_e': planet.get('rad_e', 'Not available'),
            'inc': planet.get('inc', 'Not available'),
            'inc_e': planet.get('inc_e', 'Not available'),
            'rr': planet.get('rr', 'Not available'),
            'rr_e': planet.get('rr_e', 'Not available'),
            'ar': planet.get('ar', 'Not available'),
            'ar_e': planet.get('ar_e', 'Not available'),
            'imp': planet.get('imp', 'Not available'),
            'imp_e': planet.get('imp_e', 'Not available'),
            'mass': planet.get('mass', 'Not available'),
            'mass_e': planet.get('mass_e', 'Not available'),
            'sma': planet.get('sma', 'Not available'),
            'sma_e': planet.get('sma_e', 'Not available'),
            'arg': planet.get('arg', 'Not available'),
            'arg_e': planet.get('arg_e', 'Not available'),
            'ins': planet.get('ins', 'Not available'),
            'ins_e': planet.get('ins_e', 'Not available'),
            'eqt': planet.get('eqt', 'Not available'),
            'eqt_e': planet.get('eqt_e', 'Not available'),
            'dens': planet.get('dens', 'Not available'),
            'dens_e': planet.get('dens_e', 'Not available'),
            'time': planet.get('time', 'Not available'),
            'time_e': planet.get('time_e', 'Not available'),
            'vsa': planet.get('vsa', 'Not available'),
            'vsa_e': planet.get('vsa_e', 'Not available'),
            'snr': planet.get('snr', 'Not available'),
            'snr_e': planet.get('snr_e', 'Not available'),
            'pnotes': planet.get('pnotes', 'Not available'),
            'pdate': planet.get('pdate', 'Not available'),
            'puser': planet.get('puser', 'Not available'),
            'pobid': planet.get('pobid', 'Not available'),
            'pgroup': planet.get('pgroup', 'Not available'),
            'ptag': planet.get('ptag', 'Not available')
        })

    stellar_params = data.get('stellar_parameters', [])

    stellar_data = {
        'sdate': get_stellar_parameter(stellar_params, 'sdate'),
        'dist': get_stellar_parameter(stellar_params, 'dist'),
        'teff': get_stellar_parameter(stellar_params, 'teff'),
        'teff_e': get_stellar_parameter(stellar_params, 'teff_e'),
        'mass': get_stellar_parameter(stellar_params, 'mass'),
        'mass_e': get_stellar_parameter(stellar_params, 'mass_e'),
        'logg': get_stellar_parameter(stellar_params, 'logg'),
        'logg_e': get_stellar_parameter(stellar_params, 'logg_e'),
        'srad': get_stellar_parameter(stellar_params, 'srad'),
        'srad_e': get_stellar_parameter(stellar_params, 'srad_e'),
        'logr': get_stellar_parameter(stellar_params, 'logr'),
        'logr_e': get_stellar_parameter(stellar_params, 'logr_e'),
        'sindex': get_stellar_parameter(stellar_params, 'sindex'),
        'sindex_e': get_stellar_parameter(stellar_params, 'sindex_e'),
        'halpha': get_stellar_parameter(stellar_params, 'halpha'),
        'halpha_e': get_stellar_parameter(stellar_params, 'halpha_e'),
        'vsini': get_stellar_parameter(stellar_params, 'vsini'),
        'vsini_e': get_stellar_parameter(stellar_params, 'vsini_e'),
        'rotper': get_stellar_parameter(stellar_params, 'rotper'),
        'rotper_e': get_stellar_parameter(stellar_params, 'rotper_e'),
        'met': get_stellar_parameter(stellar_params, 'met'),
        'met_e': get_stellar_parameter(stellar_params, 'met_e'),
        'dens': get_stellar_parameter(stellar_params, 'dens'),
        'dens_e': get_stellar_parameter(stellar_params, 'dens_e'),
        'lum': get_stellar_parameter(stellar_params, 'lum'),
        'lum_e': get_stellar_parameter(stellar_params, 'lum_e'),
        'rv': get_stellar_parameter(stellar_params, 'rv'),
        'rv_e': get_stellar_parameter(stellar_params, 'rv_e'),
        'age': get_stellar_parameter(stellar_params, 'age'),
        'snr': get_stellar_parameter(stellar_params, 'snr'),
        'snr_e': get_stellar_parameter(stellar_params, 'snr_e'),
        'fitq': get_stellar_parameter(stellar_params, 'fitq'),
        'snotes': get_stellar_parameter(stellar_params, 'snotes'),
        'suser': get_stellar_parameter(stellar_params, 'suser'),
        'sgroup': get_stellar_parameter(stellar_params, 'sgroup'),
        'stag': get_stellar_parameter(stellar_params, 'stag')
    }

    magnitudes = {mag['band']: mag['value'] for mag in data.get('magnitudes', [])}

    imaging_count = len(data.get('imaging', [])) if data.get('imaging') else 0
    spectroscopy_count = len(data.get('spectroscopy', [])) if data.get('spectroscopy') else 0
    time_series_count = len(data.get('time_series', [])) if data.get('time_series') else 0

    return {
        'tic_id': tic_id,
        'ra': ra,
        'dec': dec,
        'ra_hms': ra_hms,
        'dec_dms': dec_dms,
        'ra_deg': ra_deg,
        'dec_deg': dec_deg,
        'planet_data': planet_data,
        'stellar_data': stellar_data,
        'magnitudes': magnitudes,
        'imaging_count': imaging_count,
        'spectroscopy_count': spectroscopy_count,
        'time_series_count': time_series_count
    }

def save_to_csv(data, filename):
    with open(filename, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)

        # Write column descriptions
        writer.writerow(['# TIC ID: TESS Input Catalog identifier'])
        writer.writerow(['# TOI: TESS Object of Interest number'])
        writer.writerow(['# RA (hms): Right Ascension in hh:mm:ss.ss'])
        writer.writerow(['# DEC (dms): Declination in dd:mm:ss.ss '])
        writer.writerow(['# Rp (Earths): Planet radius'])
        writer.writerow(['# Rp_err (Earths): Planet radius error'])
        writer.writerow(['# Per (days): Orbital period'])
        writer.writerow(['# Per_err (days): Orbital period error'])
        writer.writerow(['# Dur (hrs): Transit duration'])
        writer.writerow(['# Dur_err (hrs): Transit duration error'])
        writer.writerow(['# Depth (ppm): Transit depth'])
        writer.writerow(['# Depth_err (ppm): Transit depth error'])
        writer.writerow(['# Inc (deg): Orbital inclination'])
        writer.writerow(['# Inc_err (deg): Orbital inclination error'])
        writer.writerow(['# a/Rstar: Semi-major axis in stellar radii'])
        writer.writerow(['# a/Rstar_err: Semi-major axis error'])
        writer.writerow(['# b: Impact parameter'])
        writer.writerow(['# b_err: Impact parameter error'])
        writer.writerow(['# Mp (Earths): Planet mass'])
        writer.writerow(['# Mp_err (Earths): Planet mass error'])
        writer.writerow(['# a (AU): Semi-major axis'])
        writer.writerow(['# a_err (AU): Semi-major axis error'])
        writer.writerow(['# omega (deg): Argument of periastron'])
        writer.writerow(['# omega_err (deg): Argument of periastron error'])
        writer.writerow(['# Insolation (Earths): Insolation flux'])
        writer.writerow(['# Insolation_err (Earths): Insolation flux error'])
        writer.writerow(['# Teq (K): Equilibrium temperature'])
        writer.writerow(['# Teq_err (K): Equilibrium temperature error'])
        writer.writerow(['# Density (g/cm^3): Planet density'])
        writer.writerow(['# Density_err (g/cm^3): Planet density error'])
        writer.writerow(['# Tp (BJD): Time of periastron passage'])
        writer.writerow(['# Tp_err (BJD): Time of periastron passage error'])
        writer.writerow(['# K (m/s): Radial velocity semi-amplitude'])
        writer.writerow(['# K_err (m/s): Radial velocity semi-amplitude error'])
        writer.writerow(['# SNR: Signal-to-noise ratio'])
        writer.writerow(['# SNR_err: Signal-to-noise ratio error'])
        writer.writerow(['# Notes: Additional notes about the planet'])
        writer.writerow(['# Updated: Date of last update'])
        writer.writerow(['# Updated By: User who last updated the data'])
        writer.writerow(['# Pobid: Planet observation ID'])
        writer.writerow(['# Pgroup: Planet group'])
        writer.writerow(['# Ptag: Planet tag'])
        writer.writerow(['# Rs (Rsun): Stellar radius'])
        writer.writerow(['# Rs_err (Rsun): Stellar radius error'])
        writer.writerow(['# Ms (Msun): Stellar mass'])
        writer.writerow(['# Ms_err (Msun): Stellar mass error'])
        writer.writerow(['# Teff (K): Effective temperature'])
        writer.writerow(['# Teff_err (K): Effective temperature error'])
        writer.writerow(['# logg: Surface gravity'])
        writer.writerow(['# logg_err: Surface gravity error'])
        writer.writerow(['# [Fe/H] (dex): Metallicity'])
        writer.writerow(['# [Fe/H]_err (dex): Metallicity error'])
        writer.writerow(['# rho (g/cm^3): Stellar density'])
        writer.writerow(['# rho_err (g/cm^3): Stellar density error'])
        writer.writerow(['# Luminosity (Lsun): Stellar luminosity'])
        writer.writerow(['# Luminosity_err (Lsun): Stellar luminosity error'])
        writer.writerow(['# Distance (pc): Distance'])
        writer.writerow(['# Age (Gyr): Age'])
        writer.writerow(['# Magnitudes: Apparent magnitudes in various bands'])
        writer.writerow(['# Imaging Observations: Number of images'])
        writer.writerow(['# Spectroscopy Observations: Number of spectra'])
        writer.writerow(['# Time Series Observations: Number of time series observations'])

        # Get the list of all unique magnitude bands from the data
        magnitude_bands = sorted(set(mag for entry in data for mag in entry['magnitudes'].keys()))

        # Write column headers
        headers = ['TIC ID', 'TOI', 'RA (hms)', 'DEC (dms)', 
                   'Epoch (BJD)', 'Epoch_err (BJD)', 'Per (days)', 'Per_err (days)', 'Dur (hrs)',
                   'Dur_err (hrs)', 'Depth (ppm)', 'Depth_err (ppm)', 'Rp (Earths)', 'Rp_err (Earths)', 'Inc (deg)',
                   'Inc_err (deg)', 'Rp/Rs', 'Rp/Rs err' ,'a/Rstar', 'a/Rstar_err', 'b', 'b_err', 'Mp (Earths)', 'Mp_err (Earths)', 'a (AU)',
                   'a_err (AU)', 'omega (deg)', 'omega_err (deg)', 'Insolation (Earths)', 'Insolation_err (Earths)',
                   'Teq (K)', 'Teq_err (K)', 'Density (g/cm^3)', 'Density_err (g/cm^3)', 'Tp (BJD)', 'Tp_err (BJD)',
                   'K (m/s)', 'K_err (m/s)', 'SNR', 'SNR_err', 'Notes', 'Updated', 'Updated By', 'Pobid', 'Pgroup', 'Ptag',
                   'Last Update (UT)', 'Distance (pc)', 'Teff (K)', 'Teff_err (K)', 'Ms (Msun)', 'Ms_err (Msun)',
                   'logg', 'logg_err', 'Rs (Rsun)', 'Rs_err (Rsun)', 'logR', 'logR_err', 'S-index', 'S-index_err',
                   'H-alpha', 'H-alpha_err', 'vsini (km/s)', 'vsini_err (km/s)', 'Prot (days)', 'Prot_err (days)',
                   '[Fe/H] (dex)', '[Fe/H]_err (dex)', 'rho (g/cm^3)', 'rho_err (g/cm^3)', 'Luminosity (Lsun)',
                   'Luminosity_err (Lsun)', 'RV (m/s)', 'RV_err (m/s)', 'Age (Gyr)', 'Stellar_SNR', 'Stellar_SNR_err',
                   'Fit Quality', 'Stellar Notes', 'Updated By', 'Group', 'Tag'] + magnitude_bands + [
                   'Imaging Observations', 'Spectroscopy Observations', 'Time Series Observations']

        writer.writerow(headers)

        for entry in data:
            magnitudes = [entry['magnitudes'].get(mag, 'Not available') for mag in magnitude_bands]
            for planet in entry['planet_data']:
                if '.' in planet['tic_id']:  # Only write rows where tic_id contains a period
                    row = [
                        planet['tic_id'], planet['toi'], entry['ra_hms'], entry['dec_dms'],
                        planet['epoch'], planet['epoch_e'],
                        planet['per'], planet['per_e'], planet['dur'], planet['dur_e'], planet['dep_p'], planet['dep_p_e'],
                        planet['rad'], planet['rad_e'], planet['inc'], planet['inc_e'], planet['rr'], planet['rr_e'],
                        planet['ar'], planet['ar_e'], planet['imp'], planet['imp_e'], 
                        planet['mass'], planet['mass_e'], planet['sma'], planet['sma_e'],
                        planet['arg'], planet['arg_e'], planet['ins'], planet['ins_e'], planet['eqt'], planet['eqt_e'],
                        planet['dens'], planet['dens_e'], planet['time'], planet['time_e'], planet['vsa'], planet['vsa_e'],
                        planet['snr'], planet['snr_e'], planet['pnotes'], planet['pdate'], planet['puser'], planet['pobid'],
                        planet['pgroup'], planet['ptag'],
                        entry['stellar_data'].get('sdate', 'Not available'),
                        entry['stellar_data'].get('dist', 'Not available'),
                        entry['stellar_data'].get('teff', 'Not available'),
                        entry['stellar_data'].get('teff_e', 'Not available'),
                        entry['stellar_data'].get('mass', 'Not available'),
                        entry['stellar_data'].get('mass_e', 'Not available'),
                        entry['stellar_data'].get('logg', 'Not available'),
                        entry['stellar_data'].get('logg_e', 'Not available'),
                        entry['stellar_data'].get('srad', 'Not available'),
                        entry['stellar_data'].get('srad_e', 'Not available'),
                        entry['stellar_data'].get('logr', 'Not available'),
                        entry['stellar_data'].get('logr_e', 'Not available'),
                        entry['stellar_data'].get('sindex', 'Not available'),
                        entry['stellar_data'].get('sindex_e', 'Not available'),
                        entry['stellar_data'].get('halpha', 'Not available'),
                        entry['stellar_data'].get('halpha_e', 'Not available'),
                        entry['stellar_data'].get('vsini', 'Not available'),
                        entry['stellar_data'].get('vsini_e', 'Not available'),
                        entry['stellar_data'].get('rotper', 'Not available'),
                        entry['stellar_data'].get('rotper_e', 'Not available'),
                        entry['stellar_data'].get('met', 'Not available'),
                        entry['stellar_data'].get('met_e', 'Not available'),
                        entry['stellar_data'].get('dens', 'Not available'),
                        entry['stellar_data'].get('dens_e', 'Not available'),
                        entry['stellar_data'].get('lum', 'Not available'),
                        entry['stellar_data'].get('lum_e', 'Not available'),
                        entry['stellar_data'].get('rv', 'Not available'),
                        entry['stellar_data'].get('rv_e', 'Not available'),
                        entry['stellar_data'].get('age', 'Not available'),
                        entry['stellar_data'].get('snr', 'Not available'),
                        entry['stellar_data'].get('snr_e', 'Not available'),
                        entry['stellar_data'].get('fitq', 'Not available'),
                        entry['stellar_data'].get('snotes', 'Not available'),
                        entry['stellar_data'].get('suser', 'Not available'),
                        entry['stellar_data'].get('sgroup', 'Not available'),
                        entry['stellar_data'].get('stag', 'Not available')
                    ] + magnitudes + [
                        entry['imaging_count'],
                        entry['spectroscopy_count'],
                        entry['time_series_count']
                    ]

                    writer.writerow(row)

folder_path = f'output_inits_files/for_exotic_py_candidate_inits_output_files/{tic_id}_files'
extracted_data = []

for file_name in os.listdir(folder_path):
    if not file_name.endswith('.json'):
        print(f"Skipping non-JSON file: {file_name}")
        continue
    file_path = os.path.join(folder_path, file_name)
    print(f"Processing file: {file_path}")
    data = load_json_data(file_path)
    if data:
        extracted_data.append(extract_data(data))

# Save the CSV with the host star TIC ID in the filename
save_to_csv(extracted_data, f'{folder_path}/{tic_id}_exofop_data.csv')
print("CSV file saved")

def display_entries(data):
    relevant_columns = {
        "RA (hms)": "Target Star RA", 
        "DEC (dms)": "Target Star Dec", 
        "TIC ID": ["Planet Name", "Host Star Name"], 
        "Per (days)": "Orbital Period (days)", 
        "Per_err (days)": "Orbital Period Uncertainty", 
        "Epoch (BJD)": "Published Mid-Transit Time (BJD-UTC)", 
        "Epoch_err (BJD)": "Mid-Transit Time Uncertainty", 
        "Rp/Rs": "Ratio of Planet to Stellar Radius (Rp/Rs)", 
        "Rp/Rs err": "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty", 
        "a/Rstar": "Ratio of Distance to Stellar Radius (a/Rs)", 
        "a/Rstar_err": "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty", 
        "Inc (deg)": "Orbital Inclination (deg)", 
        "Inc_err (deg)": "Orbital Inclination (deg) Uncertainty", 
        "omega (deg)": "Argument of Periastron (deg)", 
        "Teff (K)": "Star Effective Temperature (K)", 
        "Teff_err (K)": ["Star Effective Temperature (+) Uncertainty", "Star Effective Temperature (-) Uncertainty"], 
        "[Fe/H] (dex)": "Star Metallicity ([FE/H])", 
        "[Fe/H]_err (dex)": ["Star Metallicity (+) Uncertainty", "Star Metallicity (-) Uncertainty"], 
        "logg": "Star Surface Gravity (log(g))", 
        "logg_err": ["Star Surface Gravity (+) Uncertainty", "Star Surface Gravity (-) Uncertainty"]
    }
    
    for index, row in data.iterrows():
        print(f"Entry {index}:")
        for col, display_name in relevant_columns.items():
            if isinstance(display_name, list):
                for name in display_name:
                    print(f"{name}: {row[col]}")
            else:
                print(f"{display_name}: {row[col]}")
        print("\n")

def get_user_choice(data):
    while True:
        try:
            choice = int(input("Enter the entry number you want to select: "))
            if choice in range(len(data)):
                return choice
            else:
                print("Invalid choice. Please enter a number corresponding to an entry.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def store_parameters(data, choice):
    selected_entry = data.iloc[choice].copy()
    return selected_entry

def print_relevant_parameters(selected_entry):
    relevant_columns = {
        "RA (hms)": "Target Star RA", 
        "DEC (dms)": "Target Star Dec", 
        "TIC ID": ["Planet Name", "Host Star Name"], 
        "Per (days)": "Orbital Period (days)", 
        "Per_err (days)": "Orbital Period Uncertainty", 
        "Epoch (BJD)": "Published Mid-Transit Time (BJD-UTC)", 
        "Epoch_err (BJD)": "Mid-Transit Time Uncertainty", 
        "Rp/Rs": "Ratio of Planet to Stellar Radius (Rp/Rs)", 
        "Rp/Rs err": "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty", 
        "a/Rstar": "Ratio of Distance to Stellar Radius (a/Rs)", 
        "a/Rstar_err": "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty", 
        "Inc (deg)": "Orbital Inclination (deg)", 
        "Inc_err (deg)": "Orbital Inclination (deg) Uncertainty", 
        "omega (deg)": "Argument of Periastron (deg)", 
        "Teff (K)": "Star Effective Temperature (K)", 
        "Teff_err (K)": ["Star Effective Temperature (+) Uncertainty", "Star Effective Temperature (-) Uncertainty"], 
        "[Fe/H] (dex)": "Star Metallicity ([FE/H])", 
        "[Fe/H]_err (dex)": ["Star Metallicity (+) Uncertainty", "Star Metallicity (-) Uncertainty"], 
        "logg": "Star Surface Gravity (log(g))", 
        "logg_err": ["Star Surface Gravity (+) Uncertainty", "Star Surface Gravity (-) Uncertainty"]
    }
    
    for col, display_name in relevant_columns.items():
        if isinstance(display_name, list):
            for name in display_name:
                print(f"{name}: {selected_entry[col]}")
        else:
            print(f"{display_name}: {selected_entry[col]}")

def estimate_missing_parameters(P_days, M_star, R_star, R_p_earth):
    R_p_jupiter = (R_p_earth * u.earthRad).to(u.jupiterRad).value
    P_years = (P_days * u.day).to(u.year)
    M_star_kg = M_star * M_sun

    a = ((G * M_star_kg * P_years**2) / (4 * np.pi**2))**(1/3)
    a_AU = a.to(u.AU).value

    R_star_m = R_star * R_sun
    R_p_m = R_p_jupiter * u.jupiterRad

    pl_ratdor = (a / R_star_m).decompose().value
    pl_ratror = (R_p_m / R_star_m).decompose().value

    return a_AU, pl_ratdor, pl_ratror

def extract_host_star_name(planet_name):
    return str(planet_name).split('.')[0]

def create_inits_file(parameters, file_name):
    Teff_err = parameters.get("Teff_err (K)", [None, None])
    if not isinstance(Teff_err, list):
        Teff_err = [Teff_err, Teff_err]
    
    FeH_err = parameters.get("[Fe/H]_err (dex)", [None, None])
    if not isinstance(FeH_err, list):
        FeH_err = [FeH_err, FeH_err]
    
    logg_err = parameters.get("logg_err", [None, None])
    if not isinstance(logg_err, list):
        logg_err = [logg_err, logg_err]

    parameters = {key: value.item() if isinstance(value, np.generic) else value for key, value in parameters.items()}

    planet_name = parameters["TIC ID"]
    host_star_name = extract_host_star_name(planet_name)

    inits = {
        "planetary_parameters": {
            "Target Star RA": parameters["RA (hms)"],
            "Target Star Dec": parameters["DEC (dms)"],
            "Planet Name": planet_name,
            "Host Star Name": host_star_name,
            "Orbital Period (days)": parameters["Per (days)"],
            "Orbital Period Uncertainty": parameters["Per_err (days)"],
            "Published Mid-Transit Time (BJD-UTC)": parameters["Epoch (BJD)"],
            "Mid-Transit Time Uncertainty": parameters["Epoch_err (BJD)"],
            "Ratio of Planet to Stellar Radius (Rp/Rs)": parameters.get("Rp/Rs"),
            "Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty": parameters.get("Rp/Rs err", 0.1),
            "Ratio of Distance to Stellar Radius (a/Rs)": parameters.get("a/Rstar"),
            "Ratio of Distance to Stellar Radius (a/Rs) Uncertainty": parameters.get("a/Rstar_err", 1.0),
            "Orbital Inclination (deg)": parameters.get("Inc (deg)", 90),
            "Orbital Inclination (deg) Uncertainty": parameters.get("Inc_err (deg)", 8.0),
            "Orbital Eccentricity (0 if null)": parameters.get("Orbital Eccentricity", 0.0),
            "Argument of Periastron (deg)": parameters.get("omega (deg)", 0.0),
            "Star Effective Temperature (K)": parameters["Teff (K)"],
            "Star Effective Temperature (+) Uncertainty": Teff_err[0],
            "Star Effective Temperature (-) Uncertainty": Teff_err[1],
            "Star Metallicity ([FE/H])": parameters.get("[Fe/H] (dex)", 0.0),
            "Star Metallicity (+) Uncertainty": FeH_err[0],
            "Star Metallicity (-) Uncertainty": FeH_err[1],
            "Star Surface Gravity (log(g))": parameters["logg"],
            "Star Surface Gravity (+) Uncertainty": logg_err[0],
            "Star Surface Gravity (-) Uncertainty": logg_err[1]
        },
        "user_info": {
            "Directory with FITS files": parameters.get("Directory with FITS files", None),
            "Directory to Save Plots": parameters.get("Directory to Save Plots", None),
            "Directory of Flats": parameters.get("Directory of Flats", None),
            "Directory of Darks": parameters.get("Directory of Darks", None),
            "Directory of Biases": parameters.get("Directory of Biases", None),
            "AAVSO Observer Code (N/A if none)": parameters.get("AAVSO Observer Code (N/A if none)", "N/A"),
            "Secondary Observer Codes (N/A if none)": parameters.get("Secondary Observer Codes (N/A if none)", "N/A"),
            "Observation date": parameters.get("Observation date", None),
            "Obs. Latitude": parameters.get("Obs. Latitude", None),
            "Obs. Longitude": parameters.get("Obs. Longitude", None),
            "Obs. Elevation (meters)": parameters.get("Obs. Elevation (meters)", None),
            "Camera Type (CCD or DSLR)": parameters.get("Camera Type (CCD or DSLR)", None),
            "Pixel Binning": parameters.get("Pixel Binning", None),
            "Filter Name (aavso.org/filters)": parameters.get("Filter Name (aavso.org/filters)", None),
            "Observing Notes": parameters.get("Observing Notes", "N/A"),
            "Plate Solution? (y/n)": parameters.get("Plate Solution? (y/n)", None),
            "Align Images? (y/n)": parameters.get("Align Images? (y/n)", None),
            "Target Star X & Y Pixel": parameters.get("Target Star X & Y Pixel", None),
            "Comparison Star(s) X & Y Pixel": parameters.get("Comparison Star(s) X & Y Pixel", None)
        },    
        "optional_info": {
            "Pixel Scale (Ex: 5.21 arcsecs/pixel)": parameters.get("Pixel Scale (Ex: 5.21 arcsecs/pixel)", None),
            "Filter Minimum Wavelength (nm)": parameters.get("Filter Minimum Wavelength (nm)", None),
            "Filter Maximum Wavelength (nm)": parameters.get("Filter Maximum Wavelength (nm)", None)
        }
    }
    # Update the filename to include the planet name
    if file_name is None:
        file_name = f'output_inits_files/for_exotic_py_candidate_inits_output_files/{tic_id}_files/{planet_name}_inits.json'

    with open(file_name, 'w') as f:
        json.dump(inits, f, indent=4)
    
    return inits


# Load the CSV file
file_path = f'{folder_path}/{tic_id}_exofop_data.csv'
data = pd.read_csv(file_path, comment='#')

# Display entries for the user to choose from
display_entries(data)

# Get user choice
user_choice = get_user_choice(data)

# Store the selected parameters
stored_parameters = store_parameters(data, user_choice)

# Check for missing parameters and estimate if necessary
P_days = stored_parameters.get("Per (days)")
M_star = stored_parameters.get("Ms (Msun)")
R_star = stored_parameters.get("Rs (Rsun)")
R_p_earth = stored_parameters.get("Rp (Earths)")

missing_parameters = []
if pd.isna(stored_parameters.get("a (AU)")) or pd.isna(stored_parameters.get("Rp/Rs")) or pd.isna(stored_parameters.get("a/Rstar")):
    estimate = input("Do you want to estimate the missing parameters (Rp/Rs, a/Rs, AU)? (y/n): ").strip().lower()
    if estimate == 'y':
        if pd.notna(P_days) and pd.notna(M_star) and pd.notna(R_star) and pd.notna(R_p_earth):
            a_AU, pl_ratdor, pl_ratror = estimate_missing_parameters(P_days, M_star, R_star, R_p_earth)
            if pd.isna(stored_parameters.get("a (AU)")):
                stored_parameters["a (AU)"] = a_AU
            if pd.isna(stored_parameters.get("Rp/Rs")):
                stored_parameters["Rp/Rs"] = pl_ratror
                stored_parameters["Rp/Rs err"] = 0.1
            if pd.isna(stored_parameters.get("a/Rstar")):
                stored_parameters["a/Rstar"] = pl_ratdor
                stored_parameters["a/Rstar_err"] = 1.0
        else:
            missing_parameters = ["P_days", "M_star", "R_star", "R_p_earth"]
            print(f"Missing parameters for estimation: {missing_parameters}")

# Set default values for other parameters if they are missing
if pd.isna(stored_parameters.get("[Fe/H] (dex)")):
    stored_parameters["[Fe/H] (dex)"] = 0.0
    stored_parameters["[Fe/H]_err (dex)"] = 0.1
    print ("Missing metallicity setting to default value 0.0 with error 0.1")

if pd.isna(stored_parameters.get("Inc (deg)")):
    stored_parameters["Inc (deg)"] = 90
    stored_parameters["Inc_err (deg)"] = 8.0
    print ("Missing Inc and Inc_err setting to default value of 90 with error 8.0")

if pd.isna(stored_parameters.get("Orbital Eccentricity")):
    stored_parameters["Orbital Eccentricity"] = 0.0
    print ("Missing Orbital Eccentricity setting to default value 0.0")

if pd.isna(stored_parameters.get("omega (deg)")):
    stored_parameters["omega (deg)"] = 0.0
    print ("Missing omega setting to default value 0.0 deg")

# Prompt for missing parameters if necessary
if "P_days" in missing_parameters and pd.isna(P_days):
    P_days = float(input("Enter Orbital Period (days): "))
if "M_star" in missing_parameters and pd.isna(M_star):
    M_star = float(input("Enter Stellar Mass (Msun): "))
if "R_star" in missing_parameters and pd.isna(R_star):
    R_star = float(input("Enter Stellar Radius (Rsun): "))
if "R_p_earth" in missing_parameters and pd.isna(R_p_earth):
    R_p_earth = float(input("Enter Planet Radius (Earths): "))

# Recalculate missing parameters if all required values are provided
if missing_parameters and pd.notna(P_days) and pd.notna(M_star) and pd.notna(R_star) and pd.notna(R_p_earth):
    a_AU, pl_ratdor, pl_ratror = estimate_missing_parameters(P_days, M_star, R_star, R_p_earth)
    if pd.isna(stored_parameters.get("a (AU)")):
        stored_parameters["a (AU)"] = a_AU
    if pd.isna(stored_parameters.get("Rp/Rs")):
        stored_parameters["Rp/Rs"] = pl_ratror
        stored_parameters["Rp/Rs err"] = 0.1
    if pd.isna(stored_parameters.get("a/Rstar")):
        stored_parameters["a/Rstar"] = pl_ratdor
        stored_parameters["a/Rstar_err"] = 1.0

# Print only the relevant parameters of the chosen entry
print("Relevant Parameters of Chosen Entry:")
print_relevant_parameters(stored_parameters)

# Prompt the user for additional information if required
observer_info = input("Do you want to enter the observer location and telescope information for your targets? If (n) you can fill them in later in the outputted json files or alternatively exotic.py will prompt you for them when you run the target from inits with python exotic.py -ov., (y/n): ").strip().lower()

if observer_info == 'y':
    observer = input("Enter observer name: ").strip()
    fits_dir = input("Enter directory with FITS files: ").strip()
    plots_dir = input("Enter directory to save plots: ").strip()
    flats_dir = input("Enter directory of flats: ").strip()
    darks_dir = input("Enter directory of darks: ").strip()
    biases_dir = input("Enter directory of biases (or leave blank if none): ").strip() or None
    secondary_observers = input("Enter secondary observer codes (or 'N/A' if none): ").strip()
    obs_date = input("Enter observation date: ").strip()
    latitude = input("Enter observation latitude: ").strip()
    longitude = input("Enter observation longitude: ").strip()
    elevation = input("Enter observation elevation (meters): ").strip()
    camera_type = input("Enter camera type (CCD or DSLR): ").strip()
    pixel_binning = input("Enter pixel binning: ").strip()
    filter_name = input("Enter filter name: ").strip()
    observing_notes = input("Enter observing notes (or 'N/A' if none): ").strip()
    plate_solution = input("Plate solution? (y/n): ").strip()
    align_images = input("Align images? (y/n): ").strip()
    target_star_xy = [int(coord) for coord in input("Enter target star X & Y Pixel (comma separated): ").strip().split(',')]
    comparison_stars_xy = [
        [int(coord) for coord in star.strip().split(',')]
        for star in input("Enter comparison stars X & Y Pixel (comma separated for each star): ").strip().split(';')
    ]
    pixel_scale = float(input("Enter pixel scale (Ex: 5.21 arcsecs/pixel): ").strip())
    filter_min_wavelength = int(input("Enter filter minimum wavelength (nm): ").strip())
    filter_max_wavelength = int(input("Enter filter maximum wavelength (nm): ").strip())

    stored_parameters["AAVSO Observer Code (N/A if none)"] = observer
    stored_parameters["Directory with FITS files"] = fits_dir
    stored_parameters["Directory to Save Plots"] = plots_dir
    stored_parameters["Directory of Flats"] = flats_dir
    stored_parameters["Directory of Darks"] = darks_dir
    stored_parameters["Directory of Biases"] = biases_dir
    stored_parameters["Secondary Observer Codes (N/A if none)"] = secondary_observers
    stored_parameters["Observation date"] = obs_date
    stored_parameters["Obs. Latitude"] = latitude
    stored_parameters["Obs. Longitude"] = longitude
    stored_parameters["Obs. Elevation (meters)"] = elevation
    stored_parameters["Camera Type (CCD or DSLR)"] = camera_type
    stored_parameters["Pixel Binning"] = pixel_binning
    stored_parameters["Filter Name (aavso.org/filters)"] = filter_name
    stored_parameters["Observing Notes"] = observing_notes
    stored_parameters["Plate Solution? (y/n)"] = plate_solution
    stored_parameters["Align Images? (y/n)"] = align_images
    stored_parameters["Target Star X & Y Pixel"] = target_star_xy
    stored_parameters["Comparison Star(s) X & Y Pixel"] = comparison_stars_xy
    stored_parameters["Pixel Scale (Ex: 5.21 arcsecs/pixel)"] = pixel_scale
    stored_parameters["Filter Minimum Wavelength (nm)"] = filter_min_wavelength
    stored_parameters["Filter Maximum Wavelength (nm)"] = filter_max_wavelength

# Create the inits file and print the inits dictionary
inits = create_inits_file(stored_parameters, None)

# Print the inits dictionary
print("Inits Dictionary:")
print(json.dumps(inits, indent=4))
print ("All files created")
