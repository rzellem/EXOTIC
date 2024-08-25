# a series of functions for querying radial velocity d
import io
import os
import shutil
import tarfile
import fnmatch
import requests
import numpy as np
import pandas as pd
import urllib.request
from tqdm import tqdm
from astropy.io import fits
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from pyvo.dal import tap
import pyvo as vo

ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
tapobs = vo.dal.TAPService(ESO_TAP_OBS)
tap = vo.dal.TAPService("https://archive.eso.org/tap_obs")

def query_simbad(target):
    pos = SkyCoord.from_name(target)
    target_alias = Simbad.query_objectids(target)
    return target_alias, pos.ra.degree, pos.dec.degree

def query_nasa(target):    
    try:
        result_table = NasaExoplanetArchive.query_criteria(table="pscomppars", where=f"hostname like '{target}%'")
        gamma = float(str(result_table[0]["st_radv"]).split()[0])
    except:
        gamma = 0.0            
    print(f"gamma derived from NasaExoplanetArchive = {gamma}")
    return gamma

def query_dace(target):
    from dace_query.spectroscopy import Spectroscopy

    target_alias, ra, dec = query_simbad(target)
    aliases = [(alias[0].replace(" ", ""), ra, dec) for alias in target_alias]

    for simbad_alias, ra, dec in aliases:
        radial_velocities_table = Spectroscopy.get_timeseries(simbad_alias, sorted_by_instrument=False, output_format="pandas")

        if not radial_velocities_table.empty:
            break
    else:
        print(f"No radial velocity data found for target '{target}'.")
        return pd.DataFrame()

    valid_data = radial_velocities_table[radial_velocities_table['rv_err'] > 0]

    if valid_data.empty:
        print(f"No valid radial velocity data found for target '{target}'.")
        return pd.DataFrame()

    bjd = valid_data['rjd'] + 2400000.0
    rv = valid_data['rv']
    rv_err = valid_data['rv_err']
    ins_name = valid_data['ins_name']

    return pd.DataFrame({
        'bjd': bjd,
        'rv': rv,
        'rv_err': rv_err,
        'ins_name': ins_name
    })



def query_harps(target):
    """
    Query HARPS data for a given target using TapObs and retrieve radial velocity measurements.

    Args:
        target (str): The target name to search for.

    Returns:
        df (pd.DataFrame): A DataFrame containing BJD, radial velocity (RV), RV error,
                           and instrument name ('HARPS') columns.
    """

    # Query SIMBAD for target alias, RA, and Dec
    target_alias, ra, dec = query_simbad(target)

    # Define search radius of 2.5 arcmin in degrees
    sr = 2.5 / 60

    # Create SQL query to fetch observations within the search radius
    query = f"""SELECT *
               FROM ivoa.ObsCore
               WHERE intersects(s_region, circle('', {ra}, {dec}, {sr}))=1"""

    # Search TapObs for observations matching the query
    res = tapobs.search(query=query)

    # Extract product IDs from search results
    product_id_list = tuple(res['dp_id'])

    bjd, rv, rv_err, ins_name = [], [], [], []

    print(f"Processing {len(product_id_list)} product IDs from HARPS...")
    # Loop through product IDs to fetch HARPS data files
    for product_id in tqdm(product_id_list):
        query2 = f"""SELECT archive_id, original_filename, eso_category
                     FROM phase3v2.product_files
                     WHERE eso_category = 'ANCILLARY.HARPSTAR' AND product_id = '{product_id}'"""

        try:
            # Search TapObs for files associated with the current product ID
            res = tapobs.search(query2)
        except Exception as e:
            print(f"Error processing product ID {product_id}: {e}")
            continue

        # Loop through archive IDs to retrieve and extract HARPS data files
        for archive_id in res['archive_id']:
            file_url_table = f'https://dataportal.eso.org/dataportal_new/file/{archive_id}'
            response = requests.get(file_url_table)
            tar_data = io.BytesIO(response.content)

            with tarfile.open(fileobj=tar_data) as tar:
                for member in tar.getmembers():
                    if (member.name.endswith('A.fits') and '_ccf_' in member.name):
                        with fits.open(tar.extractfile(member)) as hdul:
                            rv_value = float(hdul[0].header['HIERARCH ESO DRS CCF RVC']) * 1000  # units = km/s * 1000
                            mjd_value = float(hdul[0].header['MJD-OBS'])  # units = 50000.0000
                            rv_err_value = float(hdul[0].header['HIERARCH ESO DRS DVRMS'])  # estimated error units = m/s

                        bjd_data = mjd_value + 2400000.0
                        if mjd_value > 57357.0:
                            bjd.append(bjd_data)
                            rv.append(rv_value)
                            rv_err.append(rv_err_value)
                            ins_name.append(hdul[0].header['INSTRUME'])

    # Return DataFrame containing BJD, RV, RV error, and instrument name columns
    return pd.DataFrame({
        'bjd': bjd,
        'rv': rv,
        'rv_err': rv_err,
        'ins_name': ins_name
    })


def query_espresso(target):
    target_alias, ra, dec = query_simbad(target)

    sr = 2.5 / 60.0  # search radius of 2.5 arcmin, always expressed in degrees
    query = f"""SELECT *
    FROM ivoa.ObsCore
    WHERE intersects(s_region, circle('', {ra}, {dec}, {sr}))=1
    """
    res = tap.search(query)
    product_id_list = tuple(res['dp_id'])

    bjd = []
    rv = []
    rv_err = []
    ins_name = []

    for product_id in tqdm(product_id_list):
        query = f"""SELECT archive_id, original_filename, eso_category 
        FROM phase3v2.product_files 
        WHERE eso_category = 'ANCILLARY.CCF' AND product_id = '{product_id}'"""
        try:
            res = tap.search(query)
        except:
            continue

        for archive_id in res['archive_id']:
            file_url_table = f'https://dataportal.eso.org/dataportal_new/file/{archive_id}'
            response = requests.get(file_url_table)
            try:
                if response.status_code == 200:
                    with io.BytesIO(response.content) as fits_data:
                        hdul = fits.open(fits_data)
                        ravariation = 36 / 3600

                        if abs(float(hdul[0].header['RA']) - ra) < ravariation and abs(float(hdul[0].header['DEC']) - dec) < ravariation:
                            bjd_data = float(hdul[0].header['HIERARCH ESO QC BJD'])
                            rv_value = float(hdul[0].header['HIERARCH ESO QC CCF RV']) * 1000.0
                            rv_err_value = float(hdul[0].header['HIERARCH ESO QC CCF RV ERROR']) * 1000.0

                            bjd.append(bjd_data)
                            rv.append(rv_value)
                            rv_err.append(rv_err_value)
                            ins_name.append(hdul[0].header['INSTRUME'])

                        hdul.close()
            except:
                print(f"Error during processing {target}, you should consider re-processing this target")

    return pd.DataFrame({
        'bjd': bjd,
        'rv': rv,
        'rv_err': rv_err,
        'ins_name': ins_name
    })

def query_neid(target):
    from pyneid.neid import Neid

    target_alias, ra, dec = query_simbad(target)

    param = {
        'datalevel': 'l2',
        'position': f'circle {ra} {dec} 0.5'
    }

    try:
        Neid.query_criteria(param, format='ipac', outpath='./criteria.tbl')
    except:
        print(f"Possible system issue during attempt to reach NEID. You might try with target {target} later.")
        return pd.DataFrame()

    shutil.rmtree('./dnload_dir', ignore_errors=True)

    try:
        Neid.download('./criteria.tbl', 'l2', 'ipac', './dnload_dir')
        targetfound = True
    except:
        print(f'No data for {target}')
        targetfound = False

    if not targetfound:
        return pd.DataFrame()

    bjd = []
    rv = []
    rv_err = []
    ins_name = []

    for filename in os.listdir('./dnload_dir'):
        if fnmatch.fnmatch(filename, '*.fits'):
            with fits.open('./dnload_dir/' + filename) as hdul:
                try:
                    bjd_data = float(hdul[12].header['ccfjdsum'])
                    rv_value = float(hdul[12].header['ccfrvmod']) * 1000.0 # convert to m/s
                    rv_err_value = float(hdul[12].header['dvrms']) * 1000.0 # convert to m/s
                    bjd.append(bjd_data)
                    rv.append(rv_value)
                    rv_err.append(rv_err_value)
                    ins_name.append('NEID')
                except:
                    pass

    shutil.rmtree('./dnload_dir', ignore_errors=True)

    return pd.DataFrame({
        'bjd': bjd,
        'rv': rv,
        'rv_err': rv_err,
        'ins_name': ins_name
    })

def query_sophie(target):

    filestub = target.replace(" ", "")
    url = f'http://atlas.obs-hp.fr/sophie/sophie.cgi?n=sophiescc&c=o&of=1,leda,simbad&nra=l,simbad,d&d=seq%2C%24link%5Bobjname%5D%2Cdate%2C%24link%5Bslen%5D%2Cmask%2Cccf_offline%2C%24field%5Brv%5D%2Cfwhm%2Cspan%2Cmaxcpp%2Ccontrast%2Clines%2C%24link%5B%27view_head%27%5D%2C%24link%5B%27get_ccf%27%5D%2C%24link%5B%27view_ccf%27%5D%2Cra%2Cdec%2Csseq%2Cdvrms%2Cobstype%2Cexpno%2Cnexp%2Csn26%2Ctargetname%2C%22target_sptype%22%2C%22target_mv%22%2Cmjd%2Cbjd&o={filestub}&a=csv%5B%2C%5D'
    response = urllib.request.urlopen(url)
    csvdata = response.read().decode('utf-8')
    sophielines = csvdata.split("\n")

    header_lines = 42
    bjd_data = []
    rv_data = []
    rv_err_data = []
    ins_name_data = []

    for line in sophielines[header_lines:]:
        sophiefields = line.split(',')

        if len(sophiefields) >= 29 and sophiefields[28].strip():
            bjd_data.append(float(sophiefields[28]))
            rv_value = (float(sophiefields[6])) * 1000.0
            rv_data.append(rv_value)
            rv_err_value = float(sophiefields[19]) * 1000.0
            rv_err_data.append(rv_err_value)
            ins_name_data.append('SOPHIE')

    # remove rv value equal to: 999000.000000
    mask = np.array(rv_data) != 999000.0

    df =  pd.DataFrame({
        'bjd': bjd_data,
        'rv': rv_data,
        'rv_err': rv_err_data,
        'ins_name': ins_name_data
    })

    return df[mask]

def query_harpscat(target):
    from ESOAsg import archive_catalogues
    
    harps_catalogue = archive_catalogues.catalogues_info(all_versions=False, collections='HARPS')
    table_harps = harps_catalogue['table_name'].data[0]
    
    target_alias, ra, dec = query_simbad(target)
    
    sr = 2.5 / 60.  # search radius of 2.5 arcmin, always expressed in degrees
    
    HARPScat = archive_catalogues.get_catalogues(tables=table_harps, type_of_query='async')
    
    bjd = []
    rv = []
    rv_err = []
    ins_name = []
    
    for i in range(len(HARPScat)):
        if abs(HARPScat['tel_targ_alpha'][i] - ra) < sr and abs(HARPScat['tel_targ_delta'][i] - dec) < sr:
            bjd.append(float(HARPScat['drs_bjd'][i]))
            rv.append(float(HARPScat['drs_ccf_rvc'][i]) * 1000.0)
            rv_err.append(float(HARPScat['drs_dvrms'][i]))
            ins_name.append('HARPScat')
    
    rv_mod_df = pd.DataFrame({
        'bjd': bjd,
        'rv': rv,
        'rv_err': rv_err,
        'ins_name': ins_name
    })
    
    return rv_mod_df

def query_hires(target):
    #from hiresprv.auth import login
    from hiresprv.download import Download
    target_alias, ra, dec = query_simbad(target)

    bjd = []
    rv = []
    rv_err = []
    ins_name = []

    print(f"Processing {len(target_alias)} target aliases from Keck HIRES...")

    for alias in target_alias:
        simbad_alias = alias[0].replace(" ", "")
        try:
            data = Download('prv.cookies', './')
        except:
            print('TEMPORARY FATAL ERROR MESSAGE: Keck HIRES archive appears to be unavailable at this time. You may need to generate the file prv.cookies. Check the hiresprv.ipac site for more information.')
            continue

        targetlist = simbad_alias.partition(' ')
        filestub = targetlist[2] if targetlist[0] == "HD" else targetlist[0] + targetlist[2] if len(targetlist) > 0 else simbad_alias

        try:
            rtn = data.rvcurve(filestub)
        except:
            print('TEMPORARY FATAL ERROR MESSAGE: Keck HIRES archive appears to be unavailable at this time. You may need to generate the file prv.cookies. Check the hiresprv.ipac site for more information.')
            continue

        try:
            with open('vst' + filestub + '.csv', 'r') as file:
                next(file)  # Skip the header line
                for line in file:
                    values = line.split(',')
                    bjd_data = float(values[1]) + 2450000.0
                    rv_value = float(values[2]) * 1000.0
                    rv_err_value = float(values[3]) * 1000.0

                    bjd.append(bjd_data)
                    rv.append(rv_value)
                    rv_err.append(rv_err_value)
                    ins_name.append('Keck HIRES')
        except FileNotFoundError:
            continue

        os.remove('vst' + filestub + '.csv')

    return pd.DataFrame({
        'bjd': bjd,
        'rv': rv,
        'rv_err': rv_err,
        'ins_name': ins_name
    })

class RadialVelocityDatabase:
    QUERY_FNS = [
        query_dace,
        query_harps,
        query_espresso,
        query_neid,
        query_sophie
    ]

    def __init__(self, target):
        self.target = target
        self.data = pd.DataFrame()
        self.collect_data()

    def collect_data(self):

        dataframes = []
        for query_fn in self.QUERY_FNS:
            try:
                print(f"Querying data for target '{self.target}' using {query_fn.__name__}...")
                df = query_fn(self.target)
                print(f"Found {len(df)} observations.")

            except Exception as e:
                print(f"Error querying data for target '{self.target}' using {query_fn.__name__}: {e}")
                continue
            
            if not df.empty:
                dataframes.append(df)
        
        if dataframes:
            self.data = pd.concat(dataframes, ignore_index=True)
            self.data.sort_values('bjd', inplace=True)
            self.data.drop_duplicates(subset=['bjd', 'rv', 'rv_err', 'ins_name'], inplace=True)
        else:
            print(f"No data found for target: {self.target}")

    def save_to_csv(self, filename):
        if not self.data.empty:
            self.data.to_csv(filename, index=False)
            print(f"Data saved to {filename}")
        else:
            print("No data to save.")

    def plot_rv_data(self):
        if self.data.empty:
            print("No data to plot.")
            return

        plt.figure(figsize=(12, 6))
        for instrument in self.data['ins_name'].unique():
            instrument_data = self.data[self.data['ins_name'] == instrument]
            plt.errorbar(instrument_data['bjd'], instrument_data['rv'], 
                         yerr=instrument_data['rv_err'], fmt='o', label=instrument)

        plt.xlabel('BJD')
        plt.ylabel('Radial Velocity (m/s)')
        plt.title(f'Radial Velocity Data for {self.target}')
        plt.legend()
        plt.grid(True)
        plt.show()

    def summary_statistics(self):
        if self.data.empty:
            print("No data available for summary statistics.")
            return

        print(f"Summary Statistics for {self.target}:")
        print(f"Total number of observations: {len(self.data)}")
        print(f"Date range: {self.data['bjd'].min():.2f} to {self.data['bjd'].max():.2f}")
        print("\nInstrument-wise statistics:")
        for instrument in self.data['ins_name'].unique():
            instrument_data = self.data[self.data['ins_name'] == instrument]
            print(f"\n{instrument}:")
            print(f"  Number of observations: {len(instrument_data)}")
            print(f"  Mean RV: {instrument_data['rv'].mean():.2f} m/s")
            print(f"  RV standard deviation: {instrument_data['rv'].std():.2f} m/s")
            print(f"  Mean RV error: {instrument_data['rv_err'].mean():.2f} m/s")