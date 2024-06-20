 
# class for RV data input

import os
import fnmatch
from astropy.io import fits
import pandas as pd
import numpy as np
from math import sin,cos,tan,asin,acos,atan,radians,degrees,pi,log,isnan
import math
from statistics import mean
from dace_query.spectroscopy import Spectroscopy
from astroquery.eso import Eso
from astroquery.simbad import Simbad
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from pyneid.neid import Neid
from hiresprv.auth import login
from hiresprv.download import Download
import urllib.request
import shutil
from datetime import datetime
import numbers
import warnings
import requests
import cgi
from pyvo.dal import tap

import matplotlib.pyplot as plt
import pandas

ESO_TAP_OBS = "http://archive.eso.org/tap_obs"

tapobs = tap.TAPService(ESO_TAP_OBS)
import sys
import pyvo as vo
from pyvo.dal import tap

# Verify the version of pyvo 
from pkg_resources import parse_version
if parse_version(vo.__version__) < parse_version('1.4'):
    raise ImportError('pyvo version must be 1.4 or higher')
    
print('pyvo             version {version}'.format(version=vo.__version__))

# Defining the ESO tap service to query for phase 3 products:
tap = vo.dal.TAPService("https://archive.eso.org/tap_obs")

def query_simbad(target):
		try:
			result_table = Simbad.query_object(target)            
			simbad_result_count = len(result_table)
		except:
			warnings.warn('Unable to perform Simbad query. Positional parameters (RA and DEC) are required to load RV data from important archives. This process will terminate.', Error)            
		rasplit = result_table[0][1]
		decsplit = result_table[0][2]
		rasplit = rasplit.split()
		decsplit = decsplit.split()
		rahh = float(rasplit[0])
		ramm = float(rasplit[1])
		rass = float(rasplit[2])
		decd = float(decsplit[0])        
		decm = float(decsplit[1])        
		decs = float(decsplit[2])        
		if decd < 0:
			decm = -decm
			decs = -decs
		radegrees = (rahh*15) + (ramm*.25) + (rass*0.00417)            
		decdegrees = (decd) + (decm/60.0) + (decs/3600.0)
		print("Converted to decimal degrees. RA = ",radegrees, " and DEC = ",decdegrees)
		i=0
		target_alias = Simbad.query_objectids(target)
		return target_alias,radegrees,decdegrees
def query_nasa(target):    
		try:
			result_table = NasaExoplanetArchive.query_criteria(table="pscomppars", where="hostname like '"+target+"%' ")
			cell_value=str(result_table[0]["st_radv"])
			cell_value_split=cell_value.split()
			gamma=float(cell_value_split[0])
		except:
			gamma=0.0            
		print("gamma derived from NasaExoplanetArchive= ",gamma)
		return(gamma)
def convert(list):
		return tuple(list)
def getDispositionFilename( response ):
		"""Get the filename from the Content-Disposition in the response's http header"""
		contentdisposition = response.headers.get('Content-Disposition')
		if contentdisposition == None:
			return None
		value, params = cgi.parse_header(contentdisposition)
		filename = params["filename"]
		return filename
def writeFile( response ):
		"""Write on disk the retrieved file"""
		if response.status_code == 200:
			# The ESO filename can be found in the response header
			filename = getDispositionFilename( response )
			# Let's write on disk the downloaded FITS spectrum using the ESO filename:
			with open(filename, 'wb') as f:
				f.write(response.content)
		return filename                      
def query_dace(target):
		target_alias,ra,dec = query_simbad(target)
		gamma = query_nasa(target)
		rhk = 3.0        
		rv_mod={}
		rv_sum_dace = 0.0
		rv_err_sum_dace = 0.0
		rv_count_dace = 0        
		print("enter DACE")
		print()
		created=False
		alias_count=0
		alias_max=len(target_alias)
		aliasnotfound=True    
		simbad_alias=''
		while alias_count<alias_max and aliasnotfound:
			simbad_alias=target_alias[alias_count][0]                 
			radial_velocities_table = Spectroscopy.get_timeseries(simbad_alias, sorted_by_instrument=False, output_format="pandas")
			print("return length from Spectroscopy for alias= ",len(radial_velocities_table))
			if len(radial_velocities_table) > 0:                
				aliasnotfound=False                
				i=0
				nt = len(radial_velocities_table['rjd'])
				print("nt in normalizer get dace ",nt)
				while i < nt:             
					bjd_data = radial_velocities_table['rjd'][i] + 2400000.0
					if created:
						rv_mod['rjd'] += [radial_velocities_table['rjd'][i]]
						rv_mod['bjd'] += [bjd_data]                    
						rv_mod['BJD'] += [bjd_data]                   
						rv_mod['t'] += [radial_velocities_table['rjd'][i]]                    
						rv_mod['berv'] += [radial_velocities_table['berv'][i]]
						rv_mod['berv_err'] += [-99999]
						rv_mod['bispan'] += [radial_velocities_table['bispan'][i]]
						rv_mod['bispan_err'] += [radial_velocities_table['bispan_err'][i]] 
						rv_mod['drift_noise'] += [radial_velocities_table['drift_noise'][i]]
						rv_mod['drift_noise_err'] += [-99999]
						rv_mod['texp'] += [radial_velocities_table['texp'][i]]
						rv_mod['texp_err'] += [-99999]
						rv_mod['cal_therror'] += [radial_velocities_table['cal_therror'][i]]  
						rv_mod['cal_therror_err'] += [-99999]
						rv_mod['protm08'] += [radial_velocities_table['protm08'][i]] 
						rv_mod['protm08_err'] += [radial_velocities_table['protm08_err'][i]] 
						rv_mod['protn84'] += [radial_velocities_table['protn84'][i]] 
						rv_mod['protn84_err'] += [radial_velocities_table['protn84_err'][i]] 
						rv_mod['fwhm'] += [radial_velocities_table['fwhm'][i]]  
						rv_mod['fwhm_err'] += [radial_velocities_table['fwhm_err'][i]]
						rv_mod['spectroFluxSn20'] += [radial_velocities_table['spectroFluxSn20'][i]]
						rv_mod['spectroFluxSn20_err'] += [-99999]
						rv_mod['cal_thfile'] += [radial_velocities_table['cal_thfile'][i]] 
						rv_mod['rhk_err'] += [radial_velocities_table['rhk_err'][i]]                
						if (math.isnan(float(radial_velocities_table['rhk'][i]))) or (radial_velocities_table['rhk'][i] == -99999.0):     
							rv_mod['rhk'] += [rhk]
						else:
							rv_mod['rhk'] += [radial_velocities_table['rhk'][i]]
						rv_mod['drs_qc'] += [radial_velocities_table['drs_qc'][i]]
						rv_mod['ins_name'] += [radial_velocities_table['ins_name'][i]]
						rv_mod['Telescope'] += [radial_velocities_table['ins_name'][i]]                    
						rv_mod['ins_mode'] += [radial_velocities_table['ins_mode'][i]]
						rv_mod['mask'] += [radial_velocities_table['mask'][i]]
						rv_mod['contrast'] += [radial_velocities_table['contrast'][i]]
						rv_mod['contrast_err'] += [radial_velocities_table['contrast_err'][i]]
						rv_mod['spectroFluxSn50'] += [radial_velocities_table['spectroFluxSn50'][i]]
						rv_mod['spectroFluxSn50_err'] += [-99999]
						rv_mod['naindex'] += [radial_velocities_table['naindex'][i]]
						rv_mod['naindex_err'] += [radial_velocities_table['naindex_err'][i]]
						rv_mod['snca2'] += [radial_velocities_table['snca2'][i]]
						rv_mod['snca2_err'] += [-99999]
						rv_mod['public'] += [radial_velocities_table['public'][i]]
						rv_mod['pub_reference'] += [radial_velocities_table['pub_reference'][i]]
						rv_mod['pub_bibcode'] += [radial_velocities_table['pub_bibcode'][i]]
						rv_mod['sindex'] += [radial_velocities_table['sindex'][i]]
						rv_mod['sindex_err'] += [radial_velocities_table['sindex_err'][i]]
						rv_mod['haindex'] += [radial_velocities_table['haindex'][i]]
						rv_mod['haindex_err'] += [radial_velocities_table['haindex_err'][i]]
						rv_mod['drs_version'] += [radial_velocities_table['drs_version'][i]]  
						rv_mod['caindex'] += [radial_velocities_table['caindex'][i]]
						rv_mod['caindex_err'] += [radial_velocities_table['caindex_err'][i]]
						rv_mod['rjd_err'] += [-99999]
						rv_value = radial_velocities_table['rv'][i]
						rv_mod['rv'] += [rv_value]
						rv_mod['Vel(m/s)'] += [rv_value]
						rv_mod['vel'] += [rv_value]                   
						rv_mod['rv_err'] += [radial_velocities_table['rv_err'][i]]
						rv_mod['ErrVel(m/s)'] += [radial_velocities_table['rv_err'][i]]
						rv_mod['errvel'] += [radial_velocities_table['rv_err'][i]]
						rv_mod['sn26'] += [-99999]                    
						rv_mod['ccf_noise'] += [radial_velocities_table['ccf_noise'][i]]
						rv_mod['ccf_noise_err'] += [-99999]
						rv_mod['drift_used'] += [radial_velocities_table['drift_used'][i]]
						rv_mod['drift_used_err'] += [-99999]
						rv_mod['ccf_asym'] += [radial_velocities_table['ccf_asym'][i]]
						rv_mod['ccf_asym_err'] += [radial_velocities_table['ccf_asym_err'][i]]
						rv_mod['date_night'] += [radial_velocities_table['date_night'][i]]  
						rv_mod['raw_file'] += [radial_velocities_table['raw_file'][i]]
						rv_mod['prog_id'] += ['']
						rv_mod['th_ar'] += [-99999]
						rv_mod['th_ar1'] += [-99999]
						rv_mod['th_ar2'] += [-99999]                        
						rv_mod['target'] += [target]
						rv_mod['RA'] += [ra]
						rv_mod['DEC'] += [dec]                        
					else:
						rv_mod['rjd'] = [radial_velocities_table['rjd'][i]]
						rv_mod['bjd'] = [bjd_data]                    
						rv_mod['BJD'] = [bjd_data]                   
						rv_mod['t'] = [radial_velocities_table['rjd'][i]]                    
						rv_mod['berv'] = [radial_velocities_table['berv'][i]]
						rv_mod['berv_err'] = [-99999]
						rv_mod['bispan'] = [radial_velocities_table['bispan'][i]]
						rv_mod['bispan_err'] = [radial_velocities_table['bispan_err'][i]] 
						rv_mod['drift_noise'] = [radial_velocities_table['drift_noise'][i]]
						rv_mod['drift_noise_err'] = [-99999]
						rv_mod['texp'] = [radial_velocities_table['texp'][i]]
						rv_mod['texp_err'] = [-99999]
						rv_mod['cal_therror'] = [radial_velocities_table['cal_therror'][i]]  
						rv_mod['cal_therror_err'] = [-99999]
						rv_mod['protm08'] = [radial_velocities_table['protm08'][i]] 
						rv_mod['protm08_err'] = [radial_velocities_table['protm08_err'][i]] 
						rv_mod['protn84'] = [radial_velocities_table['protn84'][i]] 
						rv_mod['protn84_err'] = [radial_velocities_table['protn84_err'][i]] 
						rv_mod['fwhm'] = [radial_velocities_table['fwhm'][i]]  
						rv_mod['fwhm_err'] = [radial_velocities_table['fwhm_err'][i]]
						rv_mod['spectroFluxSn20'] = [radial_velocities_table['spectroFluxSn20'][i]]
						rv_mod['spectroFluxSn20_err'] = [-99999]
						rv_mod['cal_thfile'] = [radial_velocities_table['cal_thfile'][i]] 
						rv_mod['rhk_err'] = [radial_velocities_table['rhk_err'][i]]
						if (math.isnan(float(radial_velocities_table['rhk'][i]))) or (radial_velocities_table['rhk'][i] == -99999.0):     
							rv_mod['rhk'] = [rhk]
						else:
							rv_mod['rhk'] = [radial_velocities_table['rhk'][i]]
						rv_mod['drs_qc'] = [radial_velocities_table['drs_qc'][i]]
						rv_mod['ins_name'] = [radial_velocities_table['ins_name'][i]]
						rv_mod['Telescope'] = [radial_velocities_table['ins_name'][i]]                    
						rv_mod['ins_mode'] = [radial_velocities_table['ins_mode'][i]]
						rv_mod['mask'] = [radial_velocities_table['mask'][i]]
						rv_mod['contrast'] = [radial_velocities_table['contrast'][i]]
						rv_mod['contrast_err'] = [radial_velocities_table['contrast_err'][i]]
						rv_mod['spectroFluxSn50'] = [radial_velocities_table['spectroFluxSn50'][i]]
						rv_mod['spectroFluxSn50_err'] = [-99999]
						rv_mod['naindex'] = [radial_velocities_table['naindex'][i]]
						rv_mod['naindex_err'] = [radial_velocities_table['naindex_err'][i]]
						rv_mod['snca2'] = [radial_velocities_table['snca2'][i]]
						rv_mod['snca2_err'] = [-99999]
						rv_mod['public'] = [radial_velocities_table['public'][i]]
						rv_mod['pub_reference'] = [radial_velocities_table['pub_reference'][i]]
						rv_mod['pub_bibcode'] = [radial_velocities_table['pub_bibcode'][i]]
						rv_mod['sindex'] = [radial_velocities_table['sindex'][i]]
						rv_mod['sindex_err'] = [radial_velocities_table['sindex_err'][i]]
						rv_mod['haindex'] = [radial_velocities_table['haindex'][i]]
						rv_mod['haindex_err'] = [radial_velocities_table['haindex_err'][i]]
						rv_mod['drs_version'] = [radial_velocities_table['drs_version'][i]]  
						rv_mod['caindex'] = [radial_velocities_table['caindex'][i]]
						rv_mod['caindex_err'] = [radial_velocities_table['caindex_err'][i]]
						rv_mod['rjd_err'] = [-99999]
						rv_value = radial_velocities_table['rv'][i]                    
						rv_mod['rv'] = [rv_value]
						rv_mod['Vel(m/s)'] = [rv_value]
						rv_mod['vel'] = [rv_value]                   
						rv_mod['rv_err'] = [radial_velocities_table['rv_err'][i]]
						rv_mod['ErrVel(m/s)'] = [radial_velocities_table['rv_err'][i]]
						rv_mod['errvel'] = [radial_velocities_table['rv_err'][i]]
						rv_mod['sn26'] = [-99999]                    
						rv_mod['ccf_noise'] = [radial_velocities_table['ccf_noise'][i]]
						rv_mod['ccf_noise_err'] = [-99999]
						rv_mod['drift_used'] = [radial_velocities_table['drift_used'][i]]
						rv_mod['drift_used_err'] = [-99999]
						rv_mod['ccf_asym'] = [radial_velocities_table['ccf_asym'][i]]
						rv_mod['ccf_asym_err'] = [radial_velocities_table['ccf_asym_err'][i]]
						rv_mod['date_night'] = [radial_velocities_table['date_night'][i]]  
						rv_mod['raw_file'] = [radial_velocities_table['raw_file'][i]]
						rv_mod['prog_id'] = ['']
						rv_mod['th_ar'] = [-99999]
						rv_mod['th_ar1'] = [-99999]
						rv_mod['th_ar2'] = [-99999]                        
						rv_mod['target'] = [target]
						rv_mod['RA'] = [ra]
						rv_mod['DEC'] = [dec]                        
						created=True
					rv_sum_dace = rv_sum_dace + rv_mod['rv'][i]
					rv_err_sum_dace = rv_err_sum_dace + rv_mod['rv_err'][i]
					rv_count_dace += 1
					i+=1
			alias_count+=1
		if rv_count_dace > 0:
			average_rv_dace = rv_sum_dace / rv_count_dace
			average_rv_err_dace = rv_err_sum_dace / rv_count_dace
			average_SNR_dace = rv_sum_dace / rv_err_sum_dace
			format_SNR_dace = "{:.5f}".format(average_SNR_dace)
			print(f'Number of DACE observations is: {rv_count_dace}')
			print(f'Average SNR DACE is: {format_SNR_dace}')
			print()
		rv_mod_df=pd.DataFrame.from_dict(rv_mod)       
		return rv_mod_df
###################################################################################################################           
def query_harps(target):
		target_alias,ra,dec = query_simbad(target)
		gamma = query_nasa(target)
		rhk = 3.0        
		rv_mod={}
		rv_sum_harps = 0.0
		rv_err_sum_harps = 0.0
		rv_count_harps = 0    
		print("enter HARPS15")
		print()
		created=False        
		from astropy.io import fits as pyfits
		hdulist=pyfits.open('ADP.2023-12-04T15_16_53.464.fits')
		HARPS_Table=hdulist[1].data        
		numobs = len(HARPS_Table)
		k=0
		while k<numobs:                
				object_name = HARPS_Table[k][6]
				bjd_data = float(HARPS_Table[k][13])
# On June 3, 2015 the HARPS upgrade to hexagonal fiber was implemented and data prior to that date is unreliable
				ravariation = 2/3600            
				if abs(HARPS_Table[k][35]-ra) < ravariation and abs(HARPS_Table[k][36]-dec) < ravariation and bjd_data >= 2457387.00000:
#					print('target ', object_name)
#					print("RA ",HARPS_Table[k][35], " DEC ",HARPS_Table[k][36]," Program ID ",HARPS_Table[k][2])               
					rjd_data = bjd_data - 2400000.0
					t_data = bjd_data - 2450000.0                    
					if created:
						rv_mod['rjd'] += [rjd_data]
						rv_mod['bjd'] += [float(HARPS_Table[k][13])]
						rv_mod['BJD'] += [float(HARPS_Table[k][13])]
						rv_mod['t'] += [t_data]                       
						berv_value = float(HARPS_Table[k][16])*1000.0
						rv_mod['berv'] += [berv_value]
						rv_mod['berv_err'] += [-99999]
						rv_mod['bispan'] += [-99999]
						rv_mod['bispan_err'] += [-99999] 
						rv_mod['drift_noise'] += [-99999]
						rv_mod['drift_noise_err'] += [-99999]
						rv_mod['texp'] += [-99999]
						rv_mod['texp_err'] += [-99999]
						rv_mod['cal_therror'] += [-99999]  
						rv_mod['cal_therror_err'] += [-99999]
						rv_mod['protm08'] += [-99999] 
						rv_mod['protm08_err'] += [-99999] 
						rv_mod['protn84'] += [-99999] 
						rv_mod['protn84_err'] += [-99999] 
						fwhm_value = float(HARPS_Table[k][23])*1000.0 
						rv_mod['fwhm'] += [fwhm_value]
						rv_mod['fwhm_err'] += [-99999]
						rv_mod['spectroFluxSn20'] += [-99999]
						rv_mod['spectroFluxSn20_err'] += [-99999]
						rv_mod['cal_thfile'] += [-99999] 
						rv_mod['rhk_err'] += [-99999]
						rv_mod['rhk'] += [rhk]
						rv_mod['drs_qc'] += [True]
						rv_mod['ins_name'] += ['HARPS15']
						rv_mod['Telescope'] += ['HARPS15']                         
						rv_mod['ins_mode'] += ['']
						rv_mod['mask'] += ['']
						contrast_value = float(HARPS_Table[k][22])                     
						rv_mod['contrast'] += [contrast_value]
						rv_mod['contrast_err'] += [-99999]
						rv_mod['spectroFluxSn50'] += [-99999]
						rv_mod['spectroFluxSn50_err'] += [-99999]
						rv_mod['naindex'] += [-99999]
						rv_mod['naindex_err'] += [-99999]
						rv_mod['snca2'] += [-99999]
						rv_mod['snca2_err'] += [-99999]
						rv_mod['public'] += [True]
						rv_mod['pub_reference'] += ['none']
						rv_mod['pub_bibcode'] += ['']
						rv_mod['sindex'] += [-99999]
						rv_mod['sindex_err'] += [-99999]
						rv_mod['haindex'] += [-99999]
						rv_mod['haindex_err'] += [-99999]
						rv_mod['drs_version'] += ['Pub']  
						rv_mod['caindex'] += [-99999]
						rv_mod['caindex_err'] += [-99999]
						rv_mod['rjd_err'] += [-99999]
						rv_value = (float(HARPS_Table[k][19])-gamma)*1000.0                        
						rv_mod['rv'] += [rv_value]
						rv_mod['Vel(m/s)'] += [rv_value]
						rv_mod['vel'] += [rv_value]                         
						rv_err_value = float(HARPS_Table[k][24])*1000.0                    
						rv_mod['rv_err'] += [rv_err_value]
						rv_mod['ErrVel(m/s)'] += [rv_err_value]
						rv_mod['errvel'] += [rv_err_value]
						rv_mod['sn26'] += [-99999]                        
						rv_mod['ccf_noise'] += [-99999]
						rv_mod['ccf_noise_err'] += [-99999]
						rv_mod['drift_used'] += [-99999]
						rv_mod['drift_used_err'] += [-99999]
						rv_mod['ccf_asym'] += [-99999]
						rv_mod['ccf_asym_err'] += [-99999]
						rv_mod['date_night'] += ['']  
						rv_mod['raw_file'] += [HARPS_Table[k][1]]
						rv_mod['prog_id'] += [HARPS_Table[k][2]]
						rv_mod['th_ar'] += [-99999]
						rv_mod['th_ar1'] += [-99999]
						rv_mod['th_ar2'] += [-99999]
						rv_mod['target'] += [target]
						rv_mod['RA'] += [HARPS_Table[k][35]]
						rv_mod['DEC'] += [HARPS_Table[k][36]]                        
					else:
						rv_mod['rjd'] = [rjd_data]
						rv_mod['bjd'] = [float(HARPS_Table[k][13])]
						rv_mod['BJD'] = [float(HARPS_Table[k][13])]
						rv_mod['t'] = [t_data]                       
						berv_value = float(HARPS_Table[k][16])*1000.0
						rv_mod['berv'] = [berv_value]
						rv_mod['berv_err'] = [-99999]
						rv_mod['bispan'] = [-99999]
						rv_mod['bispan_err'] = [-99999] 
						rv_mod['drift_noise'] = [-99999]
						rv_mod['drift_noise_err'] = [-99999]
						rv_mod['texp'] = [-99999]
						rv_mod['texp_err'] = [-99999]
						rv_mod['cal_therror'] = [-99999]  
						rv_mod['cal_therror_err'] = [-99999]
						rv_mod['protm08'] = [-99999] 
						rv_mod['protm08_err'] = [-99999] 
						rv_mod['protn84'] = [-99999] 
						rv_mod['protn84_err'] = [-99999] 
						fwhm_value = float(HARPS_Table[k][23])*1000.0 
						rv_mod['fwhm'] = [fwhm_value]
						rv_mod['fwhm_err'] = [-99999]
						rv_mod['spectroFluxSn20'] = [-99999]
						rv_mod['spectroFluxSn20_err'] = [-99999]
						rv_mod['cal_thfile'] = [-99999] 
						rv_mod['rhk_err'] = [-99999]
						rv_mod['rhk'] = [rhk]
						rv_mod['drs_qc'] = [True]
						rv_mod['ins_name'] = ['HARPS15']
						rv_mod['Telescope'] = ['HARPS15']                         
						rv_mod['ins_mode'] = ['']
						rv_mod['mask'] = ['']
						contrast_value = float(HARPS_Table[k][22])                     
						rv_mod['contrast'] = [contrast_value]
						rv_mod['contrast_err'] = [-99999]
						rv_mod['spectroFluxSn50'] = [-99999]
						rv_mod['spectroFluxSn50_err'] = [-99999]
						rv_mod['naindex'] = [-99999]
						rv_mod['naindex_err'] = [-99999]
						rv_mod['snca2'] = [-99999]
						rv_mod['snca2_err'] = [-99999]
						rv_mod['public'] = [True]
						rv_mod['pub_reference'] = ['none']
						rv_mod['pub_bibcode'] = ['']
						rv_mod['sindex'] = [-99999]
						rv_mod['sindex_err'] = [-99999]
						rv_mod['haindex'] = [-99999]
						rv_mod['haindex_err'] = [-99999]
						rv_mod['drs_version'] = ['Pub']  
						rv_mod['caindex'] = [-99999]
						rv_mod['caindex_err'] = [-99999]
						rv_mod['rjd_err'] = [-99999]
						rv_value = (float(HARPS_Table[k][19])-gamma)*1000.0                        
						rv_mod['rv'] = [rv_value]
						rv_mod['Vel(m/s)'] = [rv_value]
						rv_mod['vel'] = [rv_value]                         
						rv_err_value = float(HARPS_Table[k][24])*1000.0                    
						rv_mod['rv_err'] = [rv_err_value]
						rv_mod['ErrVel(m/s)'] = [rv_err_value]
						rv_mod['errvel'] = [rv_err_value]
						rv_mod['sn26'] = [-99999]                        
						rv_mod['ccf_noise'] = [-99999]
						rv_mod['ccf_noise_err'] = [-99999]
						rv_mod['drift_used'] = [-99999]
						rv_mod['drift_used_err'] = [-99999]
						rv_mod['ccf_asym'] = [-99999]
						rv_mod['ccf_asym_err'] = [-99999]
						rv_mod['date_night'] = ['']  
						rv_mod['raw_file'] = [HARPS_Table[k][1]]
						rv_mod['prog_id'] = [HARPS_Table[k][2]]
						rv_mod['th_ar'] = [-99999]
						rv_mod['th_ar1'] = [-99999]
						rv_mod['th_ar2'] = [-99999]
						rv_mod['target'] = [target]
						rv_mod['RA'] = [HARPS_Table[k][35]]
						rv_mod['DEC'] = [HARPS_Table[k][36]]                        
						created=True                        
					rv_sum_harps = rv_sum_harps + rv_value
					rv_err_sum_harps = rv_err_sum_harps + rv_err_value
					rv_count_harps +=1                                       
				k+=1
		if rv_count_harps > 0:
			average_rv_harps = rv_sum_harps / rv_count_harps
			average_rv_err_harps = rv_err_sum_harps / rv_count_harps
			average_SNR_harps = rv_sum_harps / rv_err_sum_harps
			format_SNR_harps = "{:.5f}".format(average_SNR_harps)
			print(f'Number of HARPS observations is: {rv_count_harps}')
			print(f'Average SNR HARPS is: {format_SNR_harps}')
			print()
		rv_mod_df=pd.DataFrame.from_dict(rv_mod)                        
		return rv_mod_df            
###################################################################################################################     
def query_espresso(target):

		target_alias,ra,dec = query_simbad(target)
		gamma = query_nasa(target)
		rhk = 3.0
# Read data from ESO website based upon RA and DEC from Simbad
		sr = 2.5/60. # search radius of 2.5 arcmin, always expressed in degrees
# Cone search: looking for footprints of reduced datasets intersecting a circle of 2.5' around NGC 4666
		query = """SELECT *
		FROM ivoa.ObsCore
		WHERE intersects(s_region, circle('', %f, %f, %f))=1
		""" % (ra , dec, sr)
		res = tap.search(query)
		product_id_list = convert(res['dp_id'])
		product_id_count = len(product_id_list)        
		rv_mod={}
		rv_sum_espresso = 0.0
		rv_err_sum_espresso = 0.0
		rv_count_espresso = 0
		print("enter eso espresso")
		print()
		i=0
		created=False        
		while i<product_id_count:
			query="""SELECT archive_id, original_filename, eso_category 
			FROM phase3v2.product_files 
			WHERE eso_category = 'ANCILLARY.CCF' AND product_id = '%s'""" % (product_id_list[i]) 
			res = tap.search(query)
			max_files = len(res)
			j=0
			while j < max_files:
				file_url_table = ('https://dataportal.eso.org/dataportal_new/file/'+res['archive_id'][j])
				response = requests.get(file_url_table)
				filename = writeFile( response )
				if filename:
					print("Saved file: %s" % (filename))
				else:
					print("Could not get file (status: %d)" % (response.status_code))                    
				hdul = fits.open(str(res['archive_id'][j])+'.fits')
				rjd_data = 0.0
				object_name = hdul[0].header['HIERARCH ESO OBS TARG NAME']
				RA = float(hdul[0].header['RA'])
				DEC = float(hdul[0].header['DEC'])                
				ravariation = 36/3600
				if abs(RA-ra) < ravariation and abs(DEC-dec) < ravariation:    
					bjd_data = float(hdul[0].header['HIERARCH ESO QC BJD'])
					rjd_data = bjd_data - 2400000.0                    
					t_data = bjd_data - 2450000.0                    
					if created:
						rv_mod['rjd'] += [rjd_data]
						rv_mod['bjd'] += [bjd_data]
						rv_mod['BJD'] += [bjd_data]    
						rv_mod['t'] += [t_data]                         
						berv_value = float(hdul[0].header['HIERARCH ESO QC BERV'])*1000.0                       
						rv_mod['berv'] += [berv_value]                       
						rv_mod['berv_err'] += [-99999]
						bis_value = float(hdul[0].header['HIERARCH ESO QC CCF BIS SPAN'])*1000.0
						rv_mod['bispan'] += [bis_value]                        
						bis_err_value = float(hdul[0].header['HIERARCH ESO QC CCF BIS SPAN ERROR'])*1000.0
						rv_mod['bispan_err'] += [bis_err_value]                        
						rv_mod['drift_noise'] += [-99999]
						rv_mod['drift_noise_err'] += [-99999]
						rv_mod['texp'] += [-99999]
						rv_mod['texp_err'] += [-99999]
						rv_mod['cal_therror'] += [-99999] 
						rv_mod['cal_therror_err'] += [-99999]
						rv_mod['protm08'] += [-99999] 
						rv_mod['protm08_err'] += [-99999] 
						rv_mod['protn84'] += [-99999] 
						rv_mod['protn84_err'] += [-99999] 
						fwhm_value = float(hdul[0].header['HIERARCH ESO QC CCF FWHM'])*1000.0
						rv_mod['fwhm'] += [fwhm_value]                        
						fwhm_err_value = float(hdul[0].header['HIERARCH ESO QC CCF FWHM ERROR'])*1000.0
						rv_mod['fwhm_err'] += [fwhm_err_value]                        
						rv_mod['spectroFluxSn20'] += [-99999]
						rv_mod['spectroFluxSn20_err'] += [-99999]
						rv_mod['cal_thfile'] += [-99999] 
						rv_mod['rhk_err'] += [-99999]
						rv_mod['rhk'] += [rhk]
						rv_mod['drs_qc'] += [True]
						rv_mod['ins_name'] += ['ESPRESSO']
						rv_mod['Telescope'] += ['ESPRESSO']                        
						rv_mod['ins_mode'] += ['']
						rv_mod['mask'] += ['']
						contrast_value = float(hdul[0].header['HIERARCH ESO QC CCF CONTRAST'])*1000.0
						rv_mod['contrast'] += [contrast_value]
						contrast_err_value = float(hdul[0].header['HIERARCH ESO QC CCF CONTRAST ERROR'])*1000.0
						rv_mod['contrast_err'] += [contrast_err_value]                        
						rv_mod['spectroFluxSn50'] += [-99999]
						rv_mod['spectroFluxSn50_err'] += [-99999]
						rv_mod['naindex'] += [-99999]
						rv_mod['naindex_err'] += [-99999]
						rv_mod['snca2'] += [-99999]
						rv_mod['snca2_err'] += [-99999]
						rv_mod['public'] += [True]
						rv_mod['pub_reference'] += ['none']
						rv_mod['pub_bibcode'] += ['']
						rv_mod['sindex'] += [-99999]
						rv_mod['sindex_err'] += [-99999]
						rv_mod['haindex'] += [-99999]
						rv_mod['haindex_err'] += [-99999]
						rv_mod['drs_version'] += ['Pub']  
						rv_mod['caindex'] += [-99999]
						rv_mod['caindex_err'] += [-99999]
						rv_mod['rjd_err'] += [-99999]
						rv_value = float(hdul[0].header['HIERARCH ESO QC CCF RV'])*1000.0
						rv_mod['rv'] += [rv_value]
						rv_mod['Vel(m/s)'] += [rv_value]
						rv_mod['vel'] += [rv_value]                         
						rv_err_value = float(hdul[0].header['HIERARCH ESO QC CCF RV ERROR'])*1000.0
						rv_mod['rv_err'] += [rv_err_value]
						rv_mod['ErrVel(m/s)'] += [rv_err_value]
						rv_mod['errvel'] += [rv_err_value]
						rv_mod['sn26'] += [-99999]                        
						rv_mod['ccf_noise'] += [-99999]
						rv_mod['ccf_noise_err'] += [-99999]
						rv_mod['drift_used'] += [-99999]
						rv_mod['drift_used_err'] += [-99999]
						rv_mod['ccf_asym'] += [-99999]
						rv_mod['ccf_asym_err'] += [-99999]
						rv_mod['date_night'] += ['']  
						rv_mod['raw_file'] += ['none']
						rv_mod['prog_id'] += [' ']
						rv_mod['th_ar'] += [-99999]
						rv_mod['th_ar1'] += [-99999]
						rv_mod['th_ar2'] += [-99999]
						rv_mod['target'] += [target]
						rv_mod['RA'] += [RA]
						rv_mod['DEC'] += [DEC]                        
					else:
						rv_mod['rjd'] = [rjd_data]
						rv_mod['bjd'] = [bjd_data]
						rv_mod['BJD'] = [bjd_data]    
						rv_mod['t'] = [t_data]                         
						berv_value = float(hdul[0].header['HIERARCH ESO QC BERV'])*1000.0                       
						rv_mod['berv'] = [berv_value]                       
						rv_mod['berv_err'] = [-99999]
						bis_value = float(hdul[0].header['HIERARCH ESO QC CCF BIS SPAN'])*1000.0
						rv_mod['bispan'] = [bis_value]                        
						bis_err_value = float(hdul[0].header['HIERARCH ESO QC CCF BIS SPAN ERROR'])*1000.0
						rv_mod['bispan_err'] = [bis_err_value]                        
						rv_mod['drift_noise'] = [-99999]
						rv_mod['drift_noise_err'] = [-99999]
						rv_mod['texp'] = [-99999]
						rv_mod['texp_err'] = [-99999]
						rv_mod['cal_therror'] = [-99999] 
						rv_mod['cal_therror_err'] = [-99999]
						rv_mod['protm08'] = [-99999] 
						rv_mod['protm08_err'] = [-99999] 
						rv_mod['protn84'] = [-99999] 
						rv_mod['protn84_err'] = [-99999] 
						fwhm_value = float(hdul[0].header['HIERARCH ESO QC CCF FWHM'])*1000.0
						rv_mod['fwhm'] = [fwhm_value]                        
						fwhm_err_value = float(hdul[0].header['HIERARCH ESO QC CCF FWHM ERROR'])*1000.0
						rv_mod['fwhm_err'] = [fwhm_err_value]                        
						rv_mod['spectroFluxSn20'] = [-99999]
						rv_mod['spectroFluxSn20_err'] = [-99999]
						rv_mod['cal_thfile'] = [-99999] 
						rv_mod['rhk_err'] = [-99999]
						rv_mod['rhk'] = [rhk]
						rv_mod['drs_qc'] = [True]
						rv_mod['ins_name'] = ['ESPRESSO']
						rv_mod['Telescope'] = ['ESPRESSO']                        
						rv_mod['ins_mode'] = ['']
						rv_mod['mask'] = ['']
						contrast_value = float(hdul[0].header['HIERARCH ESO QC CCF CONTRAST'])*1000.0
						rv_mod['contrast'] = [contrast_value]
						contrast_err_value = float(hdul[0].header['HIERARCH ESO QC CCF CONTRAST ERROR'])*1000.0
						rv_mod['contrast_err'] = [contrast_err_value]                        
						rv_mod['spectroFluxSn50'] = [-99999]
						rv_mod['spectroFluxSn50_err'] = [-99999]
						rv_mod['naindex'] = [-99999]
						rv_mod['naindex_err'] = [-99999]
						rv_mod['snca2'] = [-99999]
						rv_mod['snca2_err'] = [-99999]
						rv_mod['public'] = [True]
						rv_mod['pub_reference'] = ['none']
						rv_mod['pub_bibcode'] = ['']
						rv_mod['sindex'] = [-99999]
						rv_mod['sindex_err'] = [-99999]
						rv_mod['haindex'] = [-99999]
						rv_mod['haindex_err'] = [-99999]
						rv_mod['drs_version'] = ['Pub']  
						rv_mod['caindex'] = [-99999]
						rv_mod['caindex_err'] = [-99999]
						rv_mod['rjd_err'] = [-99999]
						rv_value = float(hdul[0].header['HIERARCH ESO QC CCF RV'])*1000.0
						rv_mod['rv'] = [rv_value]
						rv_mod['Vel(m/s)'] = [rv_value]
						rv_mod['vel'] = [rv_value]                         
						rv_err_value = float(hdul[0].header['HIERARCH ESO QC CCF RV ERROR'])*1000.0
						rv_mod['rv_err'] = [rv_err_value]
						rv_mod['ErrVel(m/s)'] = [rv_err_value]
						rv_mod['errvel'] = [rv_err_value]
						rv_mod['sn26'] = [-99999]                        
						rv_mod['ccf_noise'] = [-99999]
						rv_mod['ccf_noise_err'] = [-99999]
						rv_mod['drift_used'] = [-99999]
						rv_mod['drift_used_err'] = [-99999]
						rv_mod['ccf_asym'] = [-99999]
						rv_mod['ccf_asym_err'] = [-99999]
						rv_mod['date_night'] = ['']  
						rv_mod['raw_file'] = ['none']
						rv_mod['prog_id'] = [' ']
						rv_mod['th_ar'] = [-99999]
						rv_mod['th_ar1'] = [-99999]
						rv_mod['th_ar2'] = [-99999]
						rv_mod['target'] = [target]
						rv_mod['RA'] = [RA]
						rv_mod['DEC'] = [DEC]                        
						created=True
					rv_sum_espresso = rv_sum_espresso + rv_value
					rv_err_sum_espresso = rv_err_sum_espresso + rv_err_value
					rv_count_espresso +=1
				hdul.close                    
				j+=1                    
			i+=1                    
		if rv_count_espresso > 0:
			average_rv_espresso = rv_sum_espresso / rv_count_espresso
			average_rv_err_espresso = rv_err_sum_espresso / rv_count_espresso
			average_SNR_espresso = rv_sum_espresso / rv_err_sum_espresso
			format_SNR_espresso = "{:.5f}".format(average_SNR_espresso)
			print(f'Number of ESPRESSO observations is: {rv_count_espresso}')
			print(f'Average ESPRESSO is: {format_SNR_espresso}')
			print()
		rv_mod_df=pd.DataFrame.from_dict(rv_mod)                        
		return rv_mod_df
###################################################################################################################
def query_neid(target):
		target_alias,ra,dec = query_simbad(target)
		gamma = query_nasa(target)
		rhk = 3.0        
		rv_mod={}
		rv_sum_neid = 0.0
		rv_err_sum_neid = 0.0
		rv_count_neid = 0
		print("enter neid")
		print()
		created=False        
		param = dict()
		param['datalevel'] = 'l2'
		param['position'] = 'circle '+str(ra)+' '+str(dec)+' 0.5'
		print("RA from Simbad= ",ra)
		print("QDEC from Simbad= ",dec)
		print("pos parameter= ",param['position'])        
		Neid.query_criteria (param,format='ipac',outpath='./criteria.tbl')
		try:        
			shutil.rmtree('./dnload_dir')
		except:            
			print('NEID dnload_dir not present, no issue to proceed')            
		try:
			Neid.download ('./criteria.tbl','l2','ipac','./dnload_dir')            
			targetfound=True
		except:
			print('no data for ',target)
			print()
			targetfound=False
		if targetfound:            
			for filename in os.listdir('./dnload_dir'):    
				if fnmatch.fnmatch(filename, '*.fits'):
					hdul = fits.open('./dnload_dir/'+filename)                    
					object_name = hdul[0].header['OBJECT']
					try:
						bjd_data = float(hdul[12].header['ccfjdsum'])
						level2header = True
					except:
						level2header = False
					if level2header:            
						if created:
							bjd_data = float(hdul[12].header['ccfjdsum'])
							rjd_data = bjd_data - 2400000.0
							t_data = bjd_data - 2450000.0                
							rv_mod['BJD'] += [bjd_data]
							rv_mod['rjd'] += [rjd_data]
							rv_mod['bjd'] += [bjd_data]                
							rv_mod['t'] += [t_data]                
							rv_mod['berv'] += [-99999]                       
							rv_mod['berv_err'] += [-99999]
							rv_mod['bispan'] += [-99999]                        
							rv_mod['bispan_err'] += [-99999]
							rv_mod['drift_noise'] += [-99999]
							rv_mod['drift_noise_err'] += [-99999]
							rv_mod['texp'] += [-99999]
							rv_mod['texp_err'] += [-99999]
							rv_mod['cal_therror'] += [-99999] 
							rv_mod['cal_therror_err'] += [-99999]
							rv_mod['protm08'] += [-99999] 
							rv_mod['protm08_err'] += [-99999] 
							rv_mod['protn84'] += [-99999] 
							rv_mod['protn84_err'] += [-99999]              
							rv_mod['fwhm'] += [-99999]                        
							rv_mod['fwhm_err'] += [-99999]
							rv_mod['spectroFluxSn20'] += [-99999]
							rv_mod['spectroFluxSn20_err'] += [-99999]
							rv_mod['cal_thfile'] += [-99999] 
							rv_mod['rhk_err'] += [-99999]
							rv_mod['rhk'] += [rhk]
							rv_mod['drs_qc'] += [True]
							rv_mod['ins_name'] += ['NEID']
							rv_mod['Telescope'] += ['NEID']                        
							rv_mod['ins_mode'] += ['']
							rv_mod['mask'] += ['']                
							rv_mod['contrast'] += [-99999]
							rv_mod['contrast_err'] += [-99999]
							rv_mod['spectroFluxSn50'] += [-99999]
							rv_mod['spectroFluxSn50_err'] += [-99999]
							rv_mod['naindex'] += [-99999]
							rv_mod['naindex_err'] += [-99999]
							rv_mod['snca2'] += [-99999]
							rv_mod['snca2_err'] += [-99999]
							rv_mod['public'] += [True]
							rv_mod['pub_reference'] += ['none']
							rv_mod['pub_bibcode'] += ['']
							rv_mod['sindex'] += [-99999]
							rv_mod['sindex_err'] += [-99999]
							rv_mod['haindex'] += [-99999]
							rv_mod['haindex_err'] += [-99999]
							rv_mod['drs_version'] += ['Pub']  
							rv_mod['caindex'] += [-99999]
							rv_mod['caindex_err'] += [-99999]
							rv_mod['rjd_err'] += [-99999]                
							rv_value = float(hdul[12].header['ccfrvmod'])*1000.0
							rv_mod['rv'] += [rv_value]
							rv_mod['Vel(m/s)'] += [rv_value]
							rv_mod['vel'] += [rv_value]                
							rv_err_value = float(hdul[12].header['dvrms'])*1000.0
							rv_mod['rv_err'] += [rv_err_value]
							rv_mod['ErrVel(m/s)'] += [rv_err_value]
							rv_mod['errvel'] += [rv_err_value]
							rv_mod['sn26'] += [-99999]                        
							rv_mod['ccf_noise'] += [-99999]
							rv_mod['ccf_noise_err'] += [-99999]
							rv_mod['drift_used'] += [-99999]
							rv_mod['drift_used_err'] += [-99999]
							rv_mod['ccf_asym'] += [-99999]
							rv_mod['ccf_asym_err'] += [-99999]
							rv_mod['date_night'] += ['']  
							rv_mod['raw_file'] += ['none']
							rv_mod['prog_id'] += [' ']
							rv_mod['th_ar'] += [-99999]
							rv_mod['th_ar1'] += [-99999]
							rv_mod['th_ar2'] += [-99999]
							rv_mod['target'] += [target]                
							rv_mod['RA'] += [hdul[0].header['QRA']]
							rv_mod['DEC'] += [hdul[0].header['QDEC']]                                
						else:
							bjd_data = float(hdul[12].header['ccfjdsum'])
							rjd_data = bjd_data - 2400000.0
							t_data = bjd_data - 2450000.0                
							rv_mod['BJD'] = [bjd_data]
							rv_mod['rjd'] = [rjd_data]
							rv_mod['bjd'] = [bjd_data]                
							rv_mod['t'] = [t_data]                
							rv_mod['berv'] = [-99999]                       
							rv_mod['berv_err'] = [-99999]
							rv_mod['bispan'] = [-99999]                        
							rv_mod['bispan_err'] = [-99999]
							rv_mod['drift_noise'] = [-99999]
							rv_mod['drift_noise_err'] = [-99999]
							rv_mod['texp'] = [-99999]
							rv_mod['texp_err'] = [-99999]
							rv_mod['cal_therror'] = [-99999] 
							rv_mod['cal_therror_err'] = [-99999]
							rv_mod['protm08'] = [-99999] 
							rv_mod['protm08_err'] = [-99999] 
							rv_mod['protn84'] = [-99999] 
							rv_mod['protn84_err'] = [-99999]              
							rv_mod['fwhm'] = [-99999]                        
							rv_mod['fwhm_err'] = [-99999]
							rv_mod['spectroFluxSn20'] = [-99999]
							rv_mod['spectroFluxSn20_err'] = [-99999]
							rv_mod['cal_thfile'] = [-99999] 
							rv_mod['rhk_err'] = [-99999]
							rv_mod['rhk'] = [rhk]
							rv_mod['drs_qc'] = [True]
							rv_mod['ins_name'] = ['NEID']
							rv_mod['Telescope'] = ['NEID']                        
							rv_mod['ins_mode'] = ['']
							rv_mod['mask'] = ['']                
							rv_mod['contrast'] = [-99999]
							rv_mod['contrast_err'] = [-99999]
							rv_mod['spectroFluxSn50'] = [-99999]
							rv_mod['spectroFluxSn50_err'] = [-99999]
							rv_mod['naindex'] = [-99999]
							rv_mod['naindex_err'] = [-99999]
							rv_mod['snca2'] = [-99999]
							rv_mod['snca2_err'] = [-99999]
							rv_mod['public'] = [True]
							rv_mod['pub_reference'] = ['none']
							rv_mod['pub_bibcode'] = ['']
							rv_mod['sindex'] = [-99999]
							rv_mod['sindex_err'] = [-99999]
							rv_mod['haindex'] = [-99999]
							rv_mod['haindex_err'] = [-99999]
							rv_mod['drs_version'] = ['Pub']  
							rv_mod['caindex'] = [-99999]
							rv_mod['caindex_err'] = [-99999]
							rv_mod['rjd_err'] = [-99999]                
							rv_value = float(hdul[12].header['ccfrvmod'])*1000.0
							rv_mod['rv'] = [rv_value]
							rv_mod['Vel(m/s)'] = [rv_value]
							rv_mod['vel'] = [rv_value]                
							rv_err_value = float(hdul[12].header['dvrms'])*1000.0
							rv_mod['rv_err'] = [rv_err_value]
							rv_mod['ErrVel(m/s)'] = [rv_err_value]
							rv_mod['errvel'] = [rv_err_value]
							rv_mod['sn26'] = [-99999]                        
							rv_mod['ccf_noise'] = [-99999]
							rv_mod['ccf_noise_err'] = [-99999]
							rv_mod['drift_used'] = [-99999]
							rv_mod['drift_used_err'] = [-99999]
							rv_mod['ccf_asym'] = [-99999]
							rv_mod['ccf_asym_err'] = [-99999]
							rv_mod['date_night'] = ['']  
							rv_mod['raw_file'] = ['none']
							rv_mod['prog_id'] = [' ']
							rv_mod['th_ar'] = [-99999]
							rv_mod['th_ar1'] = [-99999]
							rv_mod['th_ar2'] = [-99999]
							rv_mod['target'] = [self.target]                
							rv_mod['RA'] = [hdul[0].header['QRA']]
							rv_mod['DEC'] = [hdul[0].header['QDEC']] 
							created=True
						rv_sum_neid = rv_sum_neid + rv_value
						rv_err_sum_neid = rv_err_sum_neid + rv_err_value
						rv_count_neid +=1                
					hdul.close
		if rv_count_neid > 0:
			average_rv_neid = rv_sum_neid / rv_count_neid
			average_rv_err_neid = rv_err_sum_neid / rv_count_neid
			average_SNR_neid = rv_sum_neid / rv_err_sum_neid
			format_SNR_neid = "{:.5f}".format(average_SNR_neid)
			print(f'Number of NEID observations is: {rv_count_neid}')
			print(f'Average NEID is: {format_SNR_neid}')
			print()
		rv_mod_df=pd.DataFrame.from_dict(rv_mod)                        
		return rv_mod_df
###################################################################################################################
def query_hires(target):
		target_alias,ra,dec = query_simbad(target)
		gamma = query_nasa(target)
		rhk = 3.0        
		rv_mod={}
		rv_sum_hires = 0.0
		rv_err_sum_hires = 0.0
		rv_count_hires = 0
		print("enter hires")
		print()
		created=False        
		alias_count=0
		alias_max=len(target_alias)
		aliasnotfound=True    
		simbad_alias=''
		while alias_count<alias_max and aliasnotfound:
			simbad_alias=target_alias[alias_count][0]               
			data = Download('prv.cookies', './')
			targetstring = simbad_alias
			targetlist = targetstring.partition(' ')
			if targetlist[0] == "HD":     
				filestub=targetlist[2]            
			else:
				if len(targetlist) > 0:
					filestub=targetlist[0]+targetlist[2]
				else:
					filestub=simbad_alias        
			rtn = data.rvcurve(filestub)
			try:
				num_lines = sum(1 for _ in open('vst'+filestub+'.csv'))
			except:
				num_lines = 0
			alias_count+=1            
			if num_lines > 1:
				aliasnotfound=False                
				with open('vst'+filestub+'.csv', 'r') as file:                
					header_lines=1
					count_lines=0                
					for line in file:
						if count_lines < header_lines:
							count_lines+=1
						else:                            
							values = line.split(',')    
							if created:                      
								bjd_data = float(values[1])+2450000.0
								rjd_data = bjd_data - 2400000.0
								t_data = bjd_data - 2450000.0                
								rv_mod['BJD'] += [bjd_data]
								rv_mod['rjd'] += [rjd_data]
								rv_mod['bjd'] += [bjd_data]                
								rv_mod['t'] += [t_data]                
								rv_mod['berv'] += [-99999]                       
								rv_mod['berv_err'] += [-99999]
								rv_mod['bispan'] += [-99999]                        
								rv_mod['bispan_err'] += [-99999]
								rv_mod['drift_noise'] += [-99999]
								rv_mod['drift_noise_err'] += [-99999]
								rv_mod['texp'] += [-99999]
								rv_mod['texp_err'] += [-99999]
								rv_mod['cal_therror'] += [-99999] 
								rv_mod['cal_therror_err'] += [-99999]
								rv_mod['protm08'] += [-99999] 
								rv_mod['protm08_err'] += [-99999] 
								rv_mod['protn84'] += [-99999] 
								rv_mod['protn84_err'] += [-99999]              
								rv_mod['fwhm'] += [-99999]                        
								rv_mod['fwhm_err'] += [-99999]
								rv_mod['spectroFluxSn20'] += [-99999]
								rv_mod['spectroFluxSn20_err'] += [-99999]
								rv_mod['cal_thfile'] += [-99999] 
								rv_mod['rhk_err'] += [-99999]
								rv_mod['rhk'] += [rhk]
								rv_mod['drs_qc'] += [True]
								rv_mod['ins_name'] += ['Keck HIRES']
								rv_mod['Telescope'] += ['Keck HIRES']                        
								rv_mod['ins_mode'] += ['']
								rv_mod['mask'] += ['']                
								rv_mod['contrast'] += [-99999]
								rv_mod['contrast_err'] += [-99999]
								rv_mod['spectroFluxSn50'] += [-99999]
								rv_mod['spectroFluxSn50_err'] += [-99999]
								rv_mod['naindex'] += [-99999]
								rv_mod['naindex_err'] += [-99999]
								rv_mod['snca2'] += [-99999]
								rv_mod['snca2_err'] += [-99999]
								rv_mod['public'] += [True]
								rv_mod['pub_reference'] += ['none']
								rv_mod['pub_bibcode'] += ['']
								rv_mod['sindex'] += [-99999]
								rv_mod['sindex_err'] += [-99999]
								rv_mod['haindex'] += [-99999]
								rv_mod['haindex_err'] += [-99999]
								rv_mod['drs_version'] += ['Pub']  
								rv_mod['caindex'] += [-99999]
								rv_mod['caindex_err'] += [-99999]
								rv_mod['rjd_err'] += [-99999]                
								rv_value = (float(values[2]))*1000.0
								rv_mod['rv'] += [rv_value]
								rv_mod['Vel(m/s)'] += [rv_value]
								rv_mod['vel'] += [rv_value]                
								rv_err_value = (float(values[3]))*1000.0
								rv_mod['rv_err'] += [rv_err_value]
								rv_mod['ErrVel(m/s)'] += [rv_err_value]
								rv_mod['errvel'] += [rv_err_value]
								rv_mod['sn26'] += [-99999]                        
								rv_mod['ccf_noise'] += [-99999]
								rv_mod['ccf_noise_err'] += [-99999]
								rv_mod['drift_used'] += [-99999]
								rv_mod['drift_used_err'] += [-99999]
								rv_mod['ccf_asym'] += [-99999]
								rv_mod['ccf_asym_err'] += [-99999]
								rv_mod['date_night'] += ['']  
								rv_mod['raw_file'] += ['none']
								rv_mod['prog_id'] += [' ']
								rv_mod['th_ar'] += [-99999]
								rv_mod['th_ar1'] += [-99999]
								rv_mod['th_ar2'] += [-99999]
								rv_mod['target'] += [target]                
								rv_mod['RA'] += [ra]
								rv_mod['DEC'] += [dec]                                
							else:
								bjd_data = float(values[1])+2450000.0
								rjd_data = bjd_data - 2400000.0
								t_data = bjd_data - 2450000.0                
								rv_mod['BJD'] = [bjd_data]
								rv_mod['rjd'] = [rjd_data]
								rv_mod['bjd'] = [bjd_data]                
								rv_mod['t'] = [t_data]                
								rv_mod['berv'] = [-99999]                       
								rv_mod['berv_err'] = [-99999]
								rv_mod['bispan'] = [-99999]                        
								rv_mod['bispan_err'] = [-99999]
								rv_mod['drift_noise'] = [-99999]
								rv_mod['drift_noise_err'] = [-99999]
								rv_mod['texp'] = [-99999]
								rv_mod['texp_err'] = [-99999]
								rv_mod['cal_therror'] = [-99999] 
								rv_mod['cal_therror_err'] = [-99999]
								rv_mod['protm08'] = [-99999] 
								rv_mod['protm08_err'] = [-99999] 
								rv_mod['protn84'] = [-99999] 
								rv_mod['protn84_err'] = [-99999]              
								rv_mod['fwhm'] = [-99999]                        
								rv_mod['fwhm_err'] = [-99999]
								rv_mod['spectroFluxSn20'] = [-99999]
								rv_mod['spectroFluxSn20_err'] = [-99999]
								rv_mod['cal_thfile'] = [-99999] 
								rv_mod['rhk_err'] = [-99999]
								rv_mod['rhk'] = [rhk]
								rv_mod['drs_qc'] = [True]
								rv_mod['ins_name'] = ['Keck HIRES']
								rv_mod['Telescope'] = ['Keck HIRES']                        
								rv_mod['ins_mode'] = ['']
								rv_mod['mask'] = ['']                
								rv_mod['contrast'] = [-99999]
								rv_mod['contrast_err'] = [-99999]
								rv_mod['spectroFluxSn50'] = [-99999]
								rv_mod['spectroFluxSn50_err'] = [-99999]
								rv_mod['naindex'] = [-99999]
								rv_mod['naindex_err'] = [-99999]
								rv_mod['snca2'] = [-99999]
								rv_mod['snca2_err'] = [-99999]
								rv_mod['public'] = [True]
								rv_mod['pub_reference'] = ['none']
								rv_mod['pub_bibcode'] = ['']
								rv_mod['sindex'] = [-99999]
								rv_mod['sindex_err'] = [-99999]
								rv_mod['haindex'] = [-99999]
								rv_mod['haindex_err'] = [-99999]
								rv_mod['drs_version'] = ['Pub']  
								rv_mod['caindex'] = [-99999]
								rv_mod['caindex_err'] = [-99999]
								rv_mod['rjd_err'] = [-99999]                
								rv_value = (float(values[2]))*1000.0
								rv_mod['rv'] = [rv_value]
								rv_mod['Vel(m/s)'] = [rv_value]
								rv_mod['vel'] = [rv_value]                
								rv_err_value = (float(values[3]))*1000.0
								rv_mod['rv_err'] = [rv_err_value]
								rv_mod['ErrVel(m/s)'] = [rv_err_value]
								rv_mod['errvel'] = [rv_err_value]
								rv_mod['sn26'] = [-99999]                        
								rv_mod['ccf_noise'] = [-99999]
								rv_mod['ccf_noise_err'] = [-99999]
								rv_mod['drift_used'] = [-99999]
								rv_mod['drift_used_err'] = [-99999]
								rv_mod['ccf_asym'] = [-99999]
								rv_mod['ccf_asym_err'] = [-99999]
								rv_mod['date_night'] = ['']  
								rv_mod['raw_file'] = ['none']
								rv_mod['prog_id'] = [' ']
								rv_mod['th_ar'] = [-99999]
								rv_mod['th_ar1'] = [-99999]
								rv_mod['th_ar2'] = [-99999]
								rv_mod['target'] = [target]                
								rv_mod['RA'] = [ra]
								rv_mod['DEC'] = [dec]  
								self.created=True
							rv_sum_hires = rv_sum_hires + rv_value
							rv_err_sum_hires = rv_err_sum_hires + rv_err_value
							rv_count_hires +=1
		if rv_count_hires > 0:
			average_rv_hires = rv_sum_hires / rv_count_hires
			average_rv_err_hires = rv_err_sum_hires / rv_count_hires
			average_SNR_hires = rv_sum_hires / rv_err_sum_hires
			format_SNR_hires = "{:.5f}".format(average_SNR_hires)
			print(f'Number of HIRES observations is: {rv_count_hires}')
			print(f'Average HIRES is: {format_SNR_hires}')
			print()
		rv_mod_df=pd.DataFrame.from_dict(rv_mod)                        
		return rv_mod_df
###################################################################################################################
def query_sophie(target):
		target_alias,ra,dec = query_simbad(target)
		gamma = query_nasa(target) 
		rhk = 3.0
		ravariation = 2/3600        
		rv_mod={}
		rv_sum_sophie = 0.0
		rv_err_sum_sophie = 0.0
		rv_count_sophie = 0
		print("enter sophie")
		print()
		created=False        
#SOPHIE uses SIMBAD to confirm targets so no need for the SIMBAD resolution for this archive        
		targetlist = target.partition(' ')
		if len(targetlist) > 0:
			filestub=targetlist[0]+targetlist[2]
		else:
			filestub=self.target        
		webUrl=urllib.request.urlopen('http://atlas.obs-hp.fr/sophie/sophie.cgi?n=sophiescc&c=o&of=1,leda,simbad&nra=l,simbad,d&d=seq%2C%24link%5Bobjname%5D%2Cdate%2C%24link%5Bslen%5D%2Cmask%2Cccf_offline%2C%24field%5Brv%5D%2Cfwhm%2Cspan%2Cmaxcpp%2Ccontrast%2Clines%2C%24link%5B%27view_head%27%5D%2C%24link%5B%27get_ccf%27%5D%2C%24link%5B%27view_ccf%27%5D%2Cra%2Cdec%2Csseq%2Cdvrms%2Cobstype%2Cexpno%2Cnexp%2Csn26%2Ctargetname%2C%22target_sptype%22%2C%22target_mv%22%2Cmjd%2Cbjd&o=WASP-77Ab&a=csv%5B%2C%5D')
		csvdata=webUrl.read()
		sophiedata=csvdata.decode(encoding='utf-8', errors='strict')
		sophielines=sophiedata.split("\n")        
		header_lines = 42
		count_lines = 0        
		for line in sophielines:
			if count_lines < header_lines:
				count_lines+=1
				continue
			sophiefields = line.split(',')            
			colcount=len(sophiefields)
			if colcount == 29:
				if sophiefields[28] != '\n' and sophiefields[1] != "SUN" and sophiefields[1] != "MOON" and sophiefields[16] != '' and sophiefields[17] != '':
					object_name = sophiefields[24]                                   
			try:
				num_lines = len(sophielines)-3 #trailer lines are in the returned data
			except:
				num_lines = 0            
			if num_lines > header_lines:                                               
				if count_lines < num_lines:
					if count_lines < header_lines:
						count_lines+=1
					else:
						rasophie = float(sophiefields[16])
						decsophie = float(sophiefields[17])
						if abs(rasophie-ra) < ravariation and abs(decsophie-dec) < ravariation:                        
							if created:                      
								bjd_data = float(sophiefields[28])
								rjd_data = bjd_data - 2400000.0
								t_data = bjd_data - 2450000.0                
								rv_mod['BJD'] += [bjd_data]
								rv_mod['rjd'] += [rjd_data]
								rv_mod['bjd'] += [bjd_data]                
								rv_mod['t'] += [t_data]                
								rv_mod['berv'] += [-99999]                       
								rv_mod['berv_err'] += [-99999]
								try:
									rv_mod['bispan'] += [float(sophiefields[9])]
								except:
									rv_mod['bispan'] += [-99999]
								rv_mod['bispan_err'] += [-99999]
								rv_mod['drift_noise'] += [-99999]
								rv_mod['drift_noise_err'] += [-99999]
								rv_mod['texp'] += [-99999]
								rv_mod['texp_err'] += [-99999]
								rv_mod['cal_therror'] += [-99999] 
								rv_mod['cal_therror_err'] += [-99999]
								rv_mod['protm08'] += [-99999] 
								rv_mod['protm08_err'] += [-99999] 
								rv_mod['protn84'] += [-99999] 
								rv_mod['protn84_err'] += [-99999]              
								try:
									rv_mod['fwhm'] += [float(sophiefields[8])]
								except:
									rv_mod['fwhm'] += [-99999]
								rv_mod['fwhm_err'] += [-99999]
								rv_mod['spectroFluxSn20'] += [-99999]
								rv_mod['spectroFluxSn20_err'] += [-99999]
								rv_mod['cal_thfile'] += [-99999] 
								rv_mod['rhk_err'] += [-99999]
								rv_mod['rhk'] += [rhk]
								rv_mod['drs_qc'] += [True]
								rv_mod['ins_name'] += ['SOPHIE']
								rv_mod['Telescope'] += ['SOPHIE']                        
								rv_mod['ins_mode'] += ['']
								rv_mod['mask'] += ['']                
								rv_mod['contrast'] += [-99999]
								rv_mod['contrast_err'] += [-99999]
								rv_mod['spectroFluxSn50'] += [-99999]
								rv_mod['spectroFluxSn50_err'] += [-99999]
								rv_mod['naindex'] += [-99999]
								rv_mod['naindex_err'] += [-99999]
								rv_mod['snca2'] += [-99999]
								rv_mod['snca2_err'] += [-99999]
								rv_mod['public'] += [True]
								rv_mod['pub_reference'] += ['none']
								rv_mod['pub_bibcode'] += ['']
								rv_mod['sindex'] += [-99999]
								rv_mod['sindex_err'] += [-99999]
								rv_mod['haindex'] += [-99999]
								rv_mod['haindex_err'] += [-99999]
								rv_mod['drs_version'] += ['Pub']  
								rv_mod['caindex'] += [-99999]
								rv_mod['caindex_err'] += [-99999]
								rv_mod['rjd_err'] += [-99999]                
								rv_value = (float(sophiefields[6])-gamma)*1000.0
								rv_mod['rv'] += [rv_value]
								rv_mod['Vel(m/s)'] += [rv_value]*1000.0
								rv_mod['vel'] += [rv_value]                
								rv_err_value = float(sophiefields[19])
								rv_mod['rv_err'] += [rv_err_value]
								rv_mod['ErrVel(m/s)'] += [rv_err_value]
								rv_mod['errvel'] += [rv_err_value]
								rv_mod['sn26'] += [float(sophiefields[23])]                        
								rv_mod['ccf_noise'] += [-99999]
								rv_mod['ccf_noise_err'] += [-99999]
								rv_mod['drift_used'] += [-99999]
								rv_mod['drift_used_err'] += [-99999]
								rv_mod['ccf_asym'] += [-99999]
								rv_mod['ccf_asym_err'] += [-99999]
								rv_mod['date_night'] += ['']  
								rv_mod['raw_file'] += ['none']
								rv_mod['prog_id'] += [' ']
								rv_mod['th_ar'] += [-99999]
								rv_mod['th_ar1'] += [-99999]
								rv_mod['th_ar2'] += [-99999]
								rv_mod['target'] += [target]                
								rv_mod['RA'] += [sophiefields[16]]
								rv_mod['DEC'] += [sophiefields[17]]                                
							else:
								bjd_data = float(sophiefields[28])
								rjd_data = bjd_data - 2400000.0
								t_data = bjd_data - 2450000.0                
								rv_mod['BJD'] = [bjd_data]
								rv_mod['rjd'] = [rjd_data]
								rv_mod['bjd'] = [bjd_data]                
								rv_mod['t'] = [t_data]                
								rv_mod['berv'] = [-99999]                       
								rv_mod['berv_err'] = [-99999]
								try:
									rv_mod['bispan'] = [float(sophiefields[9])]
								except:
									rv_mod['bispan'] = [-99999]                        
								rv_mod['bispan_err'] = [-99999]
								rv_mod['drift_noise'] = [-99999]
								rv_mod['drift_noise_err'] = [-99999]
								rv_mod['texp'] = [-99999]
								rv_mod['texp_err'] = [-99999]
								rv_mod['cal_therror'] = [-99999] 
								rv_mod['cal_therror_err'] = [-99999]
								rv_mod['protm08'] = [-99999] 
								rv_mod['protm08_err'] = [-99999] 
								rv_mod['protn84'] = [-99999] 
								rv_mod['protn84_err'] = [-99999]              
								try:
									rv_mod['fwhm'] = [float(sophiefields[8])]
								except:
									rv_mod['fwhm'] = [-99999]                        
								rv_mod['fwhm_err'] = [-99999]
								rv_mod['spectroFluxSn20'] = [-99999]
								rv_mod['spectroFluxSn20_err'] = [-99999]
								rv_mod['cal_thfile'] = [-99999] 
								rv_mod['rhk_err'] = [-99999]
								rv_mod['rhk'] = [rhk]
								rv_mod['drs_qc'] = [True]
								rv_mod['ins_name'] = ['SOPHIE']
								rv_mod['Telescope'] = ['SOPHIE']                        
								rv_mod['ins_mode'] = ['']
								rv_mod['mask'] = ['']                
								rv_mod['contrast'] = [-99999]
								rv_mod['contrast_err'] = [-99999]
								rv_mod['spectroFluxSn50'] = [-99999]
								rv_mod['spectroFluxSn50_err'] = [-99999]
								rv_mod['naindex'] = [-99999]
								rv_mod['naindex_err'] = [-99999]
								rv_mod['snca2'] = [-99999]
								rv_mod['snca2_err'] = [-99999]
								rv_mod['public'] = [True]
								rv_mod['pub_reference'] = ['none']
								rv_mod['pub_bibcode'] = ['']
								rv_mod['sindex'] = [-99999]
								rv_mod['sindex_err'] = [-99999]
								rv_mod['haindex'] = [-99999]
								rv_mod['haindex_err'] = [-99999]
								rv_mod['drs_version'] = ['Pub']  
								rv_mod['caindex'] = [-99999]
								rv_mod['caindex_err'] = [-99999]
								rv_mod['rjd_err'] = [-99999]                
								rv_value = (float(sophiefields[6])-gamma)*1000.0
								rv_mod['rv'] = [rv_value]
								rv_mod['Vel(m/s)'] = [rv_value]
								rv_mod['vel'] = [rv_value]                
								rv_err_value = float(sophiefields[19])*1000.0
								rv_mod['rv_err'] = [rv_err_value]
								rv_mod['ErrVel(m/s)'] = [rv_err_value]
								rv_mod['errvel'] = [rv_err_value]
								rv_mod['sn26'] = [float(sophiefields[23])]                        
								rv_mod['ccf_noise'] = [-99999]
								rv_mod['ccf_noise_err'] = [-99999]
								rv_mod['drift_used'] = [-99999]
								rv_mod['drift_used_err'] = [-99999]
								rv_mod['ccf_asym'] = [-99999]
								rv_mod['ccf_asym_err'] = [-99999]
								rv_mod['date_night'] = ['']  
								rv_mod['raw_file'] = ['none']
								rv_mod['prog_id'] = [' ']
								rv_mod['th_ar'] = [-99999]
								rv_mod['th_ar1'] = [-99999]
								rv_mod['th_ar2'] = [-99999]
								rv_mod['target'] = [target]                
								rv_mod['RA'] = [sophiefields[16]]
								rv_mod['DEC'] = [sophiefields[17]]   
								created=True
							rv_sum_sophie = rv_sum_sophie + rv_value
							rv_err_sum_sophie = rv_err_sum_sophie + rv_err_value
							rv_count_sophie +=1
						count_lines+=1
		if rv_count_sophie > 0:
			average_rv_sophie = rv_sum_sophie / rv_count_sophie
			average_rv_err_sophie = rv_err_sum_sophie / rv_count_sophie
			average_SNR_sophie = rv_sum_sophie / rv_err_sum_sophie
			format_SNR_sophie = "{:.5f}".format(average_SNR_sophie)
			print(f'Number of SOPHIE observations is: {rv_count_sophie}')
			print(f'Average SOPHIE is: {format_SNR_sophie}')
			print()
		rv_mod_df=pd.DataFrame.from_dict(rv_mod)                        
		return rv_mod_df                                                               

class RadialVelocityDatabase:

	# create list of archives to query
	query_fn = [query_dace]
	def __init__(self,target):       
		self.target=target       
	def get_data(self,target):
		self.target=target
## Instruments
# set an empty array for the instruments from which observations will be successfully collected
		self.instruments=[]
## Indicators
# stellar activity indicators are used to detect effects on radial velocity that could cause false positives
# for planetary perturbations
		self.indicators = ['rhk']            # rhk is the default        
		first_time=True        
		self.rv_data = {}
		self.rhk_value = 3.0
#		self.query_fn = [setup_rv_mod,query_dace,query_harps,query_espresso,query_neid,query_hires,query_sophie] 
		self.query_fn = [query_dace,query_harps,query_espresso,query_neid,query_hires,query_sophie]         
		for query_fn  in self.query_fn:
			query_df = query_fn(self.target)
			if first_time and len(query_df)>0:
				self.rv_data_df = query_df
				first_time = False                
			else:                
				self.rv_data_df = pd.concat([self.rv_data_df, query_df])                
		self.rv_data = self.rv_data_df.to_dict('list')      
		try:
			date_count = len(self.rv_data['rjd'])
		except:
			date_count = 0            
#sort combined sources by rjd    
		df_rv_data = pd.DataFrame.from_dict(self.rv_data)
		if date_count > 0:    
			df_rv_data.sort_values(by=['rjd'],inplace=True)
			df_rv_data = df_rv_data.reset_index(drop=True)
			sorted_rv_data = df_rv_data.to_dict(orient='list')
			nt = date_count
		else:
			print("nothing to sort for this target: ",self.target)
			nt = date_count
			sorted_rv_data = df_rv_data
		try:
			nt = len(self.rv_data['rjd'])
		except:
			nt = 0    
		if nt == 0:
			print(f'nt is: {nt}')
			print(f'No data was found for this target - {self.target}. Do not proceed.')
		else:
			print(f"average of rv_err = {mean(sorted_rv_data['rv_err'])}")
		for key in sorted_rv_data:
			try:
				if key == 'drs_qc':
					dsorted_rv_data[key] = np.array(sorted_rv_data[key], dtype=bool)
				else:
					sorted_rv_data[key] = np.array(sorted_rv_data[key], dtype=float)
			except:
				sorted_rv_data[key] = np.array(sorted_rv_data[key])
				pass
		print('Number of points before cleaning:', nt)
		if nt > 0:
# DRS quality check + nans
			for key in sorted_rv_data:
				if sorted_rv_data[key].dtype == float:
					sorted_rv_data[key][sorted_rv_data[key]==-99999] = np.nan
			keep_crit = sorted_rv_data['drs_qc']
			for key in ['rv', 'rv_err'] + self.indicators:
				keep_crit = keep_crit & (sorted_rv_data[key] == sorted_rv_data[key])
# Instruments
			if self.instruments:
				in_inst_list = False
				for inst in self.instruments:
					in_inst_list = in_inst_list | (sorted_rv_data['ins_name'] == inst)
				keep_crit = keep_crit & in_inst_list
# Apply filtering
			for key in sorted_rv_data:
				if len(sorted_rv_data[key]) == nt:
					sorted_rv_data[key] = sorted_rv_data[key][keep_crit]
# Remove empty instruments
			self.instruments = np.unique(sorted_rv_data['ins_name'])
			ninst = len(self.instruments)
# Sort by increasing time
			ksort = np.argsort(sorted_rv_data['rjd'])
			for key in sorted_rv_data:
				sorted_rv_data[key] = sorted_rv_data[key][ksort]
			nt = len(sorted_rv_data['rjd'])
			print('Number of points after cleaning:', nt)
			print('\nList of instruments kept after cleaning:')
			print(self.instruments)
		print()
		print("Exit RadialVelocityDatabase")
		print()        
		return df_rv_data    
	