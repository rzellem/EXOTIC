## This Python code is meant to be used as a replacement for radial velocity (RV) data input coding in DACE, RadVel and 
## NASA/JPL's Orbit Solver. In its current form, repositories of RV data are extracted from public sources: literature, 
## ESO-HARPS, ESO-ESPRESSO, Butler et al's KECK data, Butler et al's updated KECK data from HIRES Keck Precision Radial
## Velocity Exoplanet Survey, Tifonov et al's HARPS RV Bank v1, and the SOPHIE archive at Observatoire de Haute-Provence,
## Universit ́e d’Aix-Marseille & CNRS, F-04870 Saint Michel lObservatoire, France.
##
## Other sources of RV data and other radial velocity Python reduction environments will be considered.
##
## The RV data extracted from the sources listed above is meant to be stored on a local computer drive to allow for repeated
## usage and study without overloading the archive repositories listed and avoiding significant time loss due to repeated
## internet searches and retrievals.
##
## For simplicity and flexibiity, each source is developed as a stand alone method within the overall class. Since some
## RV comsuming code uses the full BJD format (245nnnn.nnn), some use the RJD format (245nnnn.nnn - 2400000.0) and others 
## use a four digit date (245nnnn.nnn - 2450000.0), all three versions of the time are included in the normalized data format. 
## Additional flexibility is allowed with upper and lower case instances of BJD/bjd. This means that the consuming program can 
## use BJD, bjd, rjd, or t to refer to the same Julian date without the consuming code needing to understand the distinction. 
## Similar flexibility is provided in alias versions of the radial velocity, radial velocity error and instrument. Since these 
## alias versions are "one-time" adaptations, it was considered that this cost borne by the developers of this code would be 
## saved many times over by the users of this code with their various RV consuming products.
##
## The invocation instruction contains a complex set of parameters that are meant to be set and used freely to explore
## the use of multiple RV sources simultaneously.
##
## In addition, there are quality controls built into this module: 
##      1. stellar indicators ensures that the use of Rhk, Sindex, Haindex, etc. for software that requires these indicators
##         can eliminate observations where these indicators are not present. 
##      2. instrument bucketing can assure that observations can be separated based upon spectrograph
##      3. targets can be indentified in multiple target environments which allows for larger on-computer repositories 
##      4. targets can be aliased using Simbad and fuzzy logic to accommodate those sources that use identifications that are
##         not standardized
##      5. observations can be trimmed from the input data based on observation date:
##         - remove all upper and/or lower values based upon specified date ranges
##         - remove all values within an upper and lower range (there can be several independent ranges)
##         - remove one or more unique date values as specified
##
## The calling sequence to invoke this class is:
##
##    new_rv_data = RV_Normalizer(target,instruments,indicators,rhk_value=rhk_value,
##    discovery_harps=True,discovery_coralie=True,
##    butler=False,butler_2017=False,trifonov=False,eso_harps=True,
##    eso_espresso=False,sophie=True,trim=True,trim_lower=0.0,
##    trim_upper=99999.0,remove_date_ranges=remove_date_ranges,remove_date_singles=remove_date_singles)
##
## The variables "target, instruments, indicators, rhk_value, remove_date_ranges, and
## remove_date_singles" must be initialzed in the calling program before invoking the class. Examples settings are listed
## below:
##
##    target = 'WASP-77A'
##    instruments = []
##    indicators = ['rhk']            
##    rhk_value = -4.46
##    remove_date_ranges = [[58428.620,58428.710],[57317.460,57317.600],[57302.000,57302.999]]
##    remove_date_singles = [55916.6680999999, 55832.87051999988, 55832.88853000011, 55832.89600999979, 57303.85785999987]
##
## After establishing the class, the following code will invoke the data input function:
##
##    df = new_rv_data.get_data()  # This instruction will create a dataframe to store the data
##
## Alternatively, this code can be added after the above instruction: 
##
##    rv_data = df.to_dict(orient='list') # This instruction will convert the dataframe to a  Python dictionary
## 

# class for data input and trimming

import os
import fnmatch
from astropy.io import fits
import pandas as pd
import numpy as np
from math import sin,cos,tan,asin,acos,atan,radians,degrees,pi,log,isnan
import math
from statistics import mean
from dace_query.spectroscopy import Spectroscopy
from astroquery.simbad import Simbad
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from fuzzywuzzy import process, fuzz
from datetime import datetime
import numbers

class RV_Mapper:

	def __init__(self,target,instruments,indicators,rhk_value, 
		dace,trifonov,
		butler_2017,butler_2020,eso_harps,
		eso_espresso,sophie,uves_2019,
		trim,trim_lower,trim_upper,remove_date_ranges,remove_date_singles,clip_snr,gamma):
		self.target=target
		print("self.target in normalizer ", self.target," target in normalizer ",target)
		print()
		self.instruments=instruments
		self.indicators=indicators
		self.rhk_value=rhk_value
		self.dace=dace        
#		self.discovery_TOI481=discovery_TOI481
		self.trifonov=trifonov        
		self.butler_2017=butler_2017
		self.butler_2020=butler_2020        
		self.eso_harps=eso_harps
		self.eso_espresso=eso_espresso
		self.sophie=sophie
		self.uves_2019=uves_2019        
		self.trim=trim
		self.trim_lower=trim_lower
		self.trim_upper=trim_upper
		self.remove_date_ranges=remove_date_ranges
		self.remove_date_singles=remove_date_singles
		self.clip_snr=clip_snr        
		self.created=False
		self.rv_data = {}
		self.dace_data = {}        
		self.gamma = gamma
#		self.gamma = 10.366151179475 #K2-132
#		self.today = datetime.today().strftime('%Y%m%d%H%M')
#		self.logname = 'normalizerlog'+self.today+'.txt'
#		self.startstring = 'begin log for '+self.today
#		with open(self.logname, 'w') as norm_log:
#			norm_log.write(self.startstring)
#		self.obsname = 'observationslog'+self.today+'.txt'
#		self.obsstring = 'begin log for '+self.today
#		with open(self.obsname, 'w') as observations_log:
#			observations_log.write(self.obsstring)            
	# normalizerlogYYMMDDhhmm.txt documents the fuzzy results for records to confirm that they are for the target specified 
	# observationslogYYMMDDhhmm.txt documents the accepted obserations for each source in a form for printing with Overleaf
	def get_data(self):
		self.created=False
		self.rv_sum_dace = 0.0
		self.rv_err_sum_dace = 0.0
		self.rv_count_dace = 0
#		self.rv_sum_discovery_TOI481 = 0.0
#		self.rv_err_sum_discovery_TOI481 = 0.0
#		self.rv_count_discovery_TOI481 = 0        
		self.rv_sum_discovery_hires = 0.0
		self.rv_err_sum_discovery_hires = 0.0
		self.rv_count_discovery_hires = 0        
		self.rv_sum_trifonov = 0.0
		self.rv_err_sum_trifonov = 0.0
		self.rv_count_trifonov = 0
		self.rv_sum_butler_2017 = 0.0
		self.rv_err_sum_butler_2017 = 0.0
		self.rv_count_butler_2017 = 0
		self.rv_sum_butler_2020 = 0.0        
		self.rv_err_sum_butler_2020 = 0.0
		self.rv_count_butler_2020 = 0        
		self.rv_sum_eso_harps = 0.0
		self.rv_err_sum_eso_harps = 0.0
		self.rv_count_eso_harps = 0
		self.rv_sum_espresso = 0.0
		self.rv_err_sum_espresso = 0.0
		self.rv_count_espresso = 0
		self.rv_sum_sophie = 0.0
		self.rv_err_sum_sophie = 0.0
		self.rv_count_sophie = 0
		self.rv_sum_uves_2019 = 0.0        
		self.rv_err_sum_uves_2019 = 0.0
		self.rv_count_uves_2019 = 0        
# Target set-up
		i=0
		self.target_alias = Simbad.query_objectids(self.target)
		self.target_alias_nospace = self.target_alias
		try:
			alias_count = len(self.target_alias)
		except:
			alias_count = 0        
		print(alias_count)
#		line = "Aliases found in Simbad "+str(alias_count)
#		with open(self.logname, 'a') as norm_log:
#			norm_log.write('\n')
#			norm_log.write('\n')            
#			norm_log.write(line)        
		while i<alias_count:
			print("self target alias", self.target_alias[i][0]," type  ",type(self.target_alias[i][0]))
#			line = self.target_alias[i][0]
			nospace = self.target_alias_nospace[i][0]           
			self.target_alias_nospace[i][0]=nospace.replace(" ","")            
#			with open(self.logname, 'a') as norm_log:
#				norm_log.write('\n')
#				norm_log.write(line)                                              
			i+=1
#		with open(self.logname, 'a') as norm_log:
#			norm_log.write('\n')
		result_table = Simbad.query_object(self.target)
		try:
			result_count = len(result_table)
		except:
			result_count = 0
		if result_count > 0:    
			print()        
			print(result_table[0][0])
			print()        
			print(result_table[0][1])
			print()        
			print(result_table[0][2])              
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
			self.radegrees = (rahh*15) + (ramm*.25) + (rass*0.00417)            
			self.decdegrees = (decd) + (decm/60.0) + (decs/3600.0)
			print("Converted to decimal degrees. RA = ",self.radegrees, " and DEC = ",self.decdegrees)        
		else:
			self.radegrees = 0
			self.decdegrees = 0
		if self.dace:
			self._get_dace()
#		if self.discovery_TOI481:
#			self._get_discovery_TOI481()            
		if self.trifonov:
			self._get_trifonov()
		if self.butler_2017:
			self._get_butler_2017()
		if self.butler_2020:
			self._get_butler_2020()            
		if self.eso_harps:
			self._get_eso_harps()
		if self.eso_espresso:
			self._get_eso_espresso()
		if self.sophie:
			self._get_sophie()
		if self.uves_2019:
			self._get_uves_2019()            
#		m=0                                                                               #
		try:
			date_count = len(self.rv_data['rjd'])
		except:
			date_count = 0
#		print()                                                                           #
#		while m< date_count:                                                             #
#			if self.rv_data['rjd'][m] >=  55069 and self.rv_data['rjd'][m] <= 55919:       #
#				print('prior CORALIE time before trim ',self.rv_data['rjd'][m])               #
#			m+=1                                                                          #
		if self.trim:
			print('******* self.trim is on**********')            
			self._trim_data()
#		m=0                                                                               #
#		date_count = len(self.rv_data['rjd'])                                             #
#		try:
#			date_count = len(self.rv_data['rjd'])
#		except:
#			date_count = 0
#		print()                                                                           #
#		while m< date_count:                                                             #
#			if self.rv_data['rjd'][m] >=  55069 and self.rv_data['rjd'][m] <= 55919:       #
#				print('prior CORALIE time after trim ',self.rv_data['rjd'][m])                #
#			m+=1                                                                          #
		if self.rv_count_dace > 0:
			self.average_rv_dace = self.rv_sum_dace / self.rv_count_dace
			self.average_rv_err_dace = self.rv_err_sum_dace / self.rv_count_dace
			self.average_SNR_dace = self.rv_sum_dace / self.rv_err_sum_dace
			format_SNR_dace = "{:.5f}".format(self.average_SNR_dace)
			print(f'Number of DACE observations is: {self.rv_count_dace}')
			print(f'Average SNR DACE is: {format_SNR_dace}')
			print()
            
		if self.rv_count_discovery_hires > 0:
			self.average_rv_discovery_hires = self.rv_sum_discovery_hires / self.rv_count_discovery_hires
			self.average_rv_err_discovery_hires = self.rv_err_sum_discovery_hires / self.rv_count_discovery_hires
			self.average_SNR_discovery_hires = self.rv_sum_discovery_hires / self.rv_err_sum_discovery_hires
			format_SNR_discovery_hires = "{:.5f}".format(self.average_SNR_discovery_hires)
			print(f'Number of Discovery HIRES observations is: {self.rv_count_discovery_hires}')
			print(f'Average SNR Discovery HIRES is: {format_SNR_discovery_hires}')
			print()            

		if self.rv_count_trifonov > 0:
			self.average_rv_trifonov = self.rv_sum_trifonov / self.rv_count_trifonov
			self.average_rv_err_trifonov = self.rv_err_sum_trifonov / self.rv_count_trifonov
			self.average_SNR_trifonov = self.rv_sum_trifonov / self.rv_err_sum_trifonov
			format_SNR_trifonov = "{:.5f}".format(self.average_SNR_trifonov)
			print(f'Number of Trifonov observations is: {self.rv_count_trifonov}')
			print(f'Average SNR Trifonov is: {format_SNR_trifonov}')
			print()

		if self.rv_count_butler_2017 > 0:
			self.average_rv_butler_2017 = self.rv_sum_butler_2017 / self.rv_count_butler_2017
			self.average_rv_err_butler_2017 = self.rv_err_sum_butler_2017 / self.rv_count_butler_2017
			self.average_SNR_butler_2017 = self.rv_sum_butler_2017 / self.rv_err_sum_butler_2017
			format_SNR_butler_2017 = "{:.5f}".format(self.average_SNR_butler_2017)
			print(f'Number of HIRES (Butler post) observations is: {self.rv_count_butler_2017}')
			print(f'Average SNR HIRES (Butler 2017) is: {format_SNR_butler_2017}')
			print()
            
		if self.rv_count_butler_2020 > 0:
			self.average_rv_butler_2020 = self.rv_sum_butler_2020 / self.rv_count_butler_2020
			self.average_rv_err_butler_2020 = self.rv_err_sum_butler_2020 / self.rv_count_butler_2020
			try:
				self.average_SNR_butler_2020 = self.rv_sum_butler_2020 / self.rv_err_sum_butler_2020
			except:
				self.average_SNR_butler_2020 = 99999.99
			format_SNR_butler_2020 = "{:.5f}".format(self.average_SNR_butler_2020)
			print(f'Number of HIRES (Butler post) observations is: {self.rv_count_butler_2020}')
			print(f'Average SNR HIRES (Butler 2020) is: {format_SNR_butler_2020}')
			print()            

		if self.rv_count_eso_harps > 0:
			self.average_rv_eso_harps = self.rv_sum_eso_harps / self.rv_count_eso_harps
			self.average_rv_err_eso_harps = self.rv_err_sum_eso_harps / self.rv_count_eso_harps
			self.average_SNR_eso_harps = self.rv_sum_eso_harps / self.rv_err_sum_eso_harps
			format_SNR_eso_harps = "{:.5f}".format(self.average_SNR_eso_harps)
			print(f'Number of ESO HARPS observations is: {self.rv_count_eso_harps}')
			print(f'Average SNR ESO HARPS is: {format_SNR_eso_harps}')
			print()

		if self.rv_count_espresso > 0:
			self.average_rv_espresso = self.rv_sum_espresso / self.rv_count_espresso
			self.average_rv_err_espresso = self.rv_err_sum_espresso / self.rv_count_espresso
			self.average_SNR_espresso = self.rv_sum_espresso / self.rv_err_sum_espresso
			format_SNR_espresso = "{:.5f}".format(self.average_SNR_espresso)
			print(f'Number of ESO ESPRESSO observations is: {self.rv_count_espresso}')
			print(f'Average SNR ESO ESPRESSO is: {format_SNR_espresso}')
			print()

		if self.rv_count_sophie > 0:
			self.average_rv_sophie = self.rv_sum_sophie / self.rv_count_sophie
			self.average_rv_err_sophie = self.rv_err_sum_sophie / self.rv_count_sophie
			self.average_SNR_sophie = self.rv_sum_sophie / self.rv_err_sum_sophie
			format_SNR_sophie = "{:.5f}".format(self.average_SNR_sophie)
			print(f'Number of SOPHIE observations is: {self.rv_count_sophie}')
			print(f'Average SNR SOPHIE is: {format_SNR_sophie}')
			print()
            
		if self.rv_count_uves_2019 > 0:
			self.average_rv_uves_2019 = self.rv_sum_uves_2019 / self.rv_count_uves_2019
			self.average_rv_err_uves_2019 = self.rv_err_sum_uves_2019 / self.rv_count_uves_2019
			self.average_SNR_uves_2019 = self.rv_sum_uves_2019 / self.rv_err_sum_uves_2019
			format_SNR_uves_2019 = "{:.5f}".format(self.average_SNR_uves_2019)
			print(f'Number of UVES 2019 observations is: {self.rv_count_uves_2019}')
			print(f'Average SNR UVES 2019 is: {format_SNR_uves_2019}')
			print()            
            
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
#			pass
#		df_rv_data = df_rv_data.reset_index(drop=True)
#		sorted_rv_data = df_rv_data.to_dict(orient='list')
# The following code was inserted by Suber Corley
#		nt = len(sorted_rv_data['rjd'])
		try:
			nt = len(self.rv_data['rjd'])
		except:
			nt = 0    
		print()
		print(f'target: {self.target}')
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
			print('\nList of available fields:')
			print(sorted_rv_data.keys())

		return df_rv_data

	def _get_dace(self):
		print("enter DACE")
		print()
		self.dace_data = Spectroscopy.get_timeseries(self.target, sorted_by_instrument=False, output_format="pandas")         
		i=0
		nt = len(self.dace_data['rjd'])
		print("nt in normalizer get dace ",nt)
#		pd.set_option('display.max_columns', None)
#		print(self.dace_data)
		while i < nt:             
			bjd_data = self.dace_data['rjd'][i] + 2400000.0
			if self.created:
				self.rv_data['rjd'] += [self.dace_data['rjd'][i]]
				self.rv_data['bjd'] += [bjd_data]                    
				self.rv_data['BJD'] += [bjd_data]                   
				self.rv_data['t'] += [self.dace_data['rjd'][i]]                    
				self.rv_data['berv'] += [self.dace_data['berv'][i]]
				self.rv_data['berv_err'] += [-99999]
				self.rv_data['bispan'] += [self.dace_data['bispan'][i]]
				self.rv_data['bispan_err'] += [self.dace_data['bispan_err'][i]] 
				self.rv_data['drift_noise'] += [self.dace_data['drift_noise'][i]]
				self.rv_data['drift_noise_err'] += [-99999]
				self.rv_data['texp'] += [self.dace_data['texp'][i]]
				self.rv_data['texp_err'] += [-99999]
				self.rv_data['cal_therror'] += [self.dace_data['cal_therror'][i]]  
				self.rv_data['cal_therror_err'] += [-99999]
				self.rv_data['protm08'] += [self.dace_data['protm08'][i]] 
				self.rv_data['protm08_err'] += [self.dace_data['protm08_err'][i]] 
				self.rv_data['protn84'] += [self.dace_data['protn84'][i]] 
				self.rv_data['protn84_err'] += [self.dace_data['protn84_err'][i]] 
				self.rv_data['fwhm'] += [self.dace_data['fwhm'][i]]  
				self.rv_data['fwhm_err'] += [self.dace_data['fwhm_err'][i]]
				self.rv_data['spectroFluxSn20'] += [self.dace_data['spectroFluxSn20'][i]]
				self.rv_data['spectroFluxSn20_err'] += [-99999]
				self.rv_data['cal_thfile'] += [self.dace_data['cal_thfile'][i]] 
				self.rv_data['rhk_err'] += [self.dace_data['rhk_err'][i]]
#				print("self.dace_data[rhk][i] ", self.dace_data['rhk'][i]," ",i)
				if (math.isnan(self.dace_data['rhk'][i])) or (self.dace_data['rhk'][i] == -99999.0):     					self.rv_data['rhk'] += [self.rhk_value]
				else:
					self.rv_data['rhk'] += [self.dace_data['rhk'][i]]
#				print("self.rv_data[rhk] ", self.rv_data['rhk'])
				self.rv_data['drs_qc'] += [self.dace_data['drs_qc'][i]]
				self.rv_data['ins_name'] += [self.dace_data['ins_name'][i]]
				self.rv_data['Telescope'] += [self.dace_data['ins_name'][i]]                    
				self.rv_data['ins_mode'] += [self.dace_data['ins_mode'][i]]
				self.rv_data['mask'] += [self.dace_data['mask'][i]]
				self.rv_data['contrast'] += [self.dace_data['contrast'][i]]
				self.rv_data['contrast_err'] += [self.dace_data['contrast_err'][i]]
				self.rv_data['spectroFluxSn50'] += [self.dace_data['spectroFluxSn50'][i]]
				self.rv_data['spectroFluxSn50_err'] += [-99999]
				self.rv_data['naindex'] += [self.dace_data['naindex'][i]]
				self.rv_data['naindex_err'] += [self.dace_data['naindex_err'][i]]
				self.rv_data['snca2'] += [self.dace_data['snca2'][i]]
				self.rv_data['snca2_err'] += [-99999]
				self.rv_data['public'] += [self.dace_data['public'][i]]
				self.rv_data['pub_reference'] += [self.dace_data['pub_reference'][i]]
				self.rv_data['pub_bibcode'] += [self.dace_data['pub_bibcode'][i]]
				self.rv_data['sindex'] += [self.dace_data['sindex'][i]]
				self.rv_data['sindex_err'] += [self.dace_data['sindex_err'][i]]
				self.rv_data['haindex'] += [self.dace_data['haindex'][i]]
				self.rv_data['haindex_err'] += [self.dace_data['haindex_err'][i]]
				self.rv_data['drs_version'] += [self.dace_data['drs_version'][i]]  
				self.rv_data['caindex'] += [self.dace_data['caindex'][i]]
				self.rv_data['caindex_err'] += [self.dace_data['caindex_err'][i]]
				self.rv_data['rjd_err'] += [-99999]
				rv_value = self.dace_data['rv'][i]                    
				self.rv_data['rv'] += [rv_value]
				self.rv_data['Vel(m/s)'] += [rv_value]
				self.rv_data['vel'] += [rv_value]                   
				self.rv_data['rv_err'] += [self.dace_data['rv_err'][i]]
				self.rv_data['ErrVel(m/s)'] += [self.dace_data['rv_err'][i]]
				self.rv_data['errvel'] += [self.dace_data['rv_err'][i]]
				self.rv_data['sn26'] += [-99999]                    
				self.rv_data['ccf_noise'] += [self.dace_data['ccf_noise'][i]]
				self.rv_data['ccf_noise_err'] += [-99999]
				self.rv_data['drift_used'] += [self.dace_data['drift_used'][i]]
				self.rv_data['drift_used_err'] += [-99999]
				self.rv_data['ccf_asym'] += [self.dace_data['ccf_asym'][i]]
				self.rv_data['ccf_asym_err'] += [self.dace_data['ccf_asym_err'][i]]
				self.rv_data['date_night'] += [self.dace_data['date_night'][i]]  
				self.rv_data['raw_file'] += [self.dace_data['raw_file'][i]]
				self.rv_data['prog_id'] += [self.dace_data['prog_id'][i]]
				self.rv_data['th_ar'] += [self.dace_data['th_ar'][i]]
				self.rv_data['th_ar1'] += [self.dace_data['th_ar1'][i]]
				self.rv_data['th_ar2'] += [self.dace_data['th_ar2'][i]]                    
			else:
				self.rv_data['rjd'] = [self.dace_data['rjd'][i]]
				self.rv_data['bjd'] = [bjd_data]                    
				self.rv_data['BJD'] = [bjd_data]                   
				self.rv_data['t'] = [self.dace_data['rjd'][i]]                    
				self.rv_data['berv'] = [self.dace_data['berv'][i]]
				self.rv_data['berv_err'] = [-99999]
				self.rv_data['bispan'] = [self.dace_data['bispan'][i]]
				self.rv_data['bispan_err'] = [self.dace_data['bispan_err'][i]] 
				self.rv_data['drift_noise'] = [self.dace_data['drift_noise'][i]]
				self.rv_data['drift_noise_err'] = [-99999]
				self.rv_data['texp'] = [self.dace_data['texp'][i]]
				self.rv_data['texp_err'] = [-99999]
				self.rv_data['cal_therror'] = [self.dace_data['cal_therror'][i]]  
				self.rv_data['cal_therror_err'] = [-99999]
				self.rv_data['protm08'] = [self.dace_data['protm08'][i]] 
				self.rv_data['protm08_err'] = [self.dace_data['protm08_err'][i]] 
				self.rv_data['protn84'] = [self.dace_data['protn84'][i]] 
				self.rv_data['protn84_err'] = [self.dace_data['protn84_err'][i]] 
				self.rv_data['fwhm'] = [self.dace_data['fwhm'][i]]  
				self.rv_data['fwhm_err'] = [self.dace_data['fwhm_err'][i]]
				self.rv_data['spectroFluxSn20'] = [self.dace_data['spectroFluxSn20'][i]]
				self.rv_data['spectroFluxSn20_err'] = [-99999]
				self.rv_data['cal_thfile'] = [self.dace_data['cal_thfile'][i]] 
				self.rv_data['rhk_err'] = [self.dace_data['rhk_err'][i]]
				if (math.isnan(self.dace_data['rhk'][i])) or (self.dace_data['rhk'][i] == -99999.0):
					self.rv_data['rhk'] = [self.rhk_value]
				else:
					self.rv_data['rhk'] = [self.dace_data['rhk'][i]]                        
				self.rv_data['drs_qc'] = [self.dace_data['drs_qc'][i]]
				self.rv_data['ins_name'] = [self.dace_data['ins_name'][i]]
				self.rv_data['Telescope'] = [self.dace_data['ins_name'][i]]                    
				self.rv_data['ins_mode'] = [self.dace_data['ins_mode'][i]]
				self.rv_data['mask'] = [self.dace_data['mask'][i]]
				self.rv_data['contrast'] = [self.dace_data['contrast'][i]]
				self.rv_data['contrast_err'] = [self.dace_data['contrast_err'][i]]
				self.rv_data['spectroFluxSn50'] = [self.dace_data['spectroFluxSn50'][i]]
				self.rv_data['spectroFluxSn50_err'] = [-99999]
				self.rv_data['naindex'] = [self.dace_data['naindex'][i]]
				self.rv_data['naindex_err'] = [self.dace_data['naindex_err'][i]]
				self.rv_data['snca2'] = [self.dace_data['snca2'][i]]
				self.rv_data['snca2_err'] = [-99999]
				self.rv_data['public'] = [self.dace_data['public'][i]]
				self.rv_data['pub_reference'] = [self.dace_data['pub_reference'][i]]
				self.rv_data['pub_bibcode'] = [self.dace_data['pub_bibcode'][i]]
				self.rv_data['sindex'] = [self.dace_data['sindex'][i]]
				self.rv_data['sindex_err'] = [self.dace_data['sindex_err'][i]]
				self.rv_data['haindex'] = [self.dace_data['haindex'][i]]
				self.rv_data['haindex_err'] = [self.dace_data['haindex_err'][i]]
				self.rv_data['drs_version'] = [self.dace_data['drs_version'][i]]  
				self.rv_data['caindex'] = [self.dace_data['caindex'][i]]
				self.rv_data['caindex_err'] = [self.dace_data['caindex_err'][i]]
				self.rv_data['rjd_err'] = [-99999]
				rv_value = self.dace_data['rv'][i]                    
				self.rv_data['rv'] = [rv_value]
				self.rv_data['Vel(m/s)'] = [rv_value]
				self.rv_data['vel'] = [rv_value]                   
				self.rv_data['rv_err'] = [self.dace_data['rv_err'][i]]
				self.rv_data['ErrVel(m/s)'] = [self.dace_data['rv_err'][i]]
				self.rv_data['errvel'] = [self.dace_data['rv_err'][i]]
				self.rv_data['sn26'] = [-99999]                    
				self.rv_data['ccf_noise'] = [self.dace_data['ccf_noise'][i]]
				self.rv_data['ccf_noise_err'] = [-99999]
				self.rv_data['drift_used'] = [self.dace_data['drift_used'][i]]
				self.rv_data['drift_used_err'] = [-99999]
				self.rv_data['ccf_asym'] = [self.dace_data['ccf_asym'][i]]
				self.rv_data['ccf_asym_err'] = [self.dace_data['ccf_asym_err'][i]]
				self.rv_data['date_night'] = [self.dace_data['date_night'][i]]  
				self.rv_data['raw_file'] = [self.dace_data['raw_file'][i]]
				self.rv_data['prog_id'] = [self.dace_data['prog_id'][i]]
				self.rv_data['th_ar'] = [self.dace_data['th_ar'][i]]
				self.rv_data['th_ar1'] = [self.dace_data['th_ar1'][i]]
				self.rv_data['th_ar2'] = [self.dace_data['th_ar2'][i]]
				self.created=True
			self.rv_sum_dace = self.rv_sum_dace + self.rv_data['rv'][i]
			self.rv_err_sum_dace = self.rv_err_sum_dace + self.dace_data['rv_err'][i]
			self.rv_count_dace += 1
			i+=1
                       
                     
	def _get_discovery_hires(self):
		print("enter HIRES for WASP77A")
		print()        
		myfile21 = open ('WASP_77_A_Discovery_HIRES.csv', 'r')
		for line in myfile21:
			values = line.split(',')
			object_name = values[7]
			i=0
			imax=0
			max_fuzzy=0
			while i < len(self.target_alias):
				String_Matched = fuzz.ratio(object_name, self.target_alias[i][0])
				if String_Matched > max_fuzzy:
					max_fuzzy=String_Matched
					imax=i
				i+=1
			if max_fuzzy > 80:
#				line = "String Matched % "+str(max_fuzzy)+"target object"+object_name+ "target_max "+ self.target_alias[imax][0]
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)                
				t_data = float(values[3]) - 2450000.0
				if self.created:
					self.rv_data['rjd'] += [float(values[6])]
					self.rv_data['bjd'] += [float(values[3])]                    
					self.rv_data['BJD'] += [float(values[3])]                   
					self.rv_data['t'] += [t_data]                    
					self.rv_data['berv'] += [-99999]
					self.rv_data['berv_err'] += [-99999]
					self.rv_data['bispan'] += [-99999]
					self.rv_data['bispan_err'] += [-99999] 
					self.rv_data['drift_noise'] += [-99999]
					self.rv_data['drift_noise_err'] += [-99999]
					self.rv_data['texp'] += [-99999]
					self.rv_data['texp_err'] += [-99999]
					self.rv_data['cal_therror'] += [-99999]  
					self.rv_data['cal_therror_err'] += [-99999]
					self.rv_data['protm08'] += [-99999] 
					self.rv_data['protm08_err'] += [-99999] 
					self.rv_data['protn84'] += [-99999] 
					self.rv_data['protn84_err'] += [-99999] 
					self.rv_data['fwhm'] += [-99999]  
					self.rv_data['fwhm_err'] += [-99999]
					self.rv_data['spectroFluxSn20'] += [-99999]
					self.rv_data['spectroFluxSn20_err'] += [-99999]
					self.rv_data['cal_thfile'] += [-99999] 
					self.rv_data['rhk_err'] += [-99999]
					self.rv_data['rhk'] += [self.rhk_value]
					self.rv_data['drs_qc'] += [True]
					self.rv_data['ins_name'] += [values[0]]
					self.rv_data['Telescope'] += [values[0]]                    
					self.rv_data['ins_mode'] += ['']
					self.rv_data['mask'] += ['']
					self.rv_data['contrast'] += [-99999]
					self.rv_data['contrast_err'] += [-99999]
					self.rv_data['spectroFluxSn50'] += [-99999]
					self.rv_data['spectroFluxSn50_err'] += [-99999]
					self.rv_data['naindex'] += [-99999]
					self.rv_data['naindex_err'] += [-99999]
					self.rv_data['snca2'] += [-99999]
					self.rv_data['snca2_err'] += [-99999]
					self.rv_data['public'] += [True]
					self.rv_data['pub_reference'] += ['Maxted']
					self.rv_data['pub_bibcode'] += ['']
					self.rv_data['sindex'] += [-99999]
					self.rv_data['sindex_err'] += [-99999]
					self.rv_data['haindex'] += [-99999]
					self.rv_data['haindex_err'] += [-99999]
					self.rv_data['drs_version'] += ['Pub']  
					self.rv_data['caindex'] += [-99999]
					self.rv_data['caindex_err'] += [-99999]
					self.rv_data['rjd_err'] += [-99999]
					rv_value = (float(values[1]))*1000.0                    
#					rv_value = float(values[4])*1.0
					self.rv_data['rv'] += [rv_value]
					self.rv_data['Vel(m/s)'] += [rv_value]
					self.rv_data['vel'] += [rv_value]                   
					rv_err_value = float(values[5])*1000.0
					self.rv_data['rv_err'] += [rv_err_value]
					self.rv_data['ErrVel(m/s)'] += [rv_err_value]
					self.rv_data['errvel'] += [rv_err_value]
					self.rv_data['sn26'] += [-99999]                    
					self.rv_data['ccf_noise'] += [-99999]
					self.rv_data['ccf_noise_err'] += [-99999]
					self.rv_data['drift_used'] += [-99999]
					self.rv_data['drift_used_err'] += [-99999]
					self.rv_data['ccf_asym'] += [-99999]
					self.rv_data['ccf_asym_err'] += [-99999]
					self.rv_data['date_night'] += ['']  
					self.rv_data['raw_file'] += ['Maxted']
					self.rv_data['prog_id'] += [' ']
					self.rv_data['th_ar'] += [-99999]
					self.rv_data['th_ar1'] += [-99999]
					self.rv_data['th_ar2'] += [-99999]                    
				else:
					self.rv_data['rjd'] = [float(values[6])]                    
					self.rv_data['bjd'] = [float(values[3])]                    
					self.rv_data['BJD'] = [float(values[3])]                    
					self.rv_data['t'] = [t_data]                   
					self.rv_data['berv'] = [-99999]
					self.rv_data['berv_err'] = [-99999]
					self.rv_data['bispan'] = [-99999]
					self.rv_data['bispan_err'] = [-99999] 
					self.rv_data['drift_noise'] = [-99999]
					self.rv_data['drift_noise_err'] = [-99999]
					self.rv_data['texp'] = [-99999]
					self.rv_data['texp_err'] = [-99999]
					self.rv_data['cal_therror'] = [-99999]  
					self.rv_data['cal_therror_err'] = [-99999]
					self.rv_data['protm08'] = [-99999] 
					self.rv_data['protm08_err'] = [-99999] 
					self.rv_data['protn84'] = [-99999] 
					self.rv_data['protn84_err'] = [-99999] 
					self.rv_data['fwhm'] = [-99999]  
					self.rv_data['fwhm_err'] = [-99999]
					self.rv_data['spectroFluxSn20'] = [-99999]
					self.rv_data['spectroFluxSn20_err'] = [-99999]
					self.rv_data['cal_thfile'] = [-99999] 
					self.rv_data['rhk_err'] = [-99999]
					self.rv_data['rhk'] = [self.rhk_value]
					self.rv_data['drs_qc'] = [True]
					self.rv_data['ins_name'] = [values[0]]
					self.rv_data['Telescope'] = [values[0]]                    
					self.rv_data['ins_mode'] = ['']
					self.rv_data['mask'] = ['']
					self.rv_data['contrast'] = [-99999]
					self.rv_data['contrast_err'] = [-99999]
					self.rv_data['spectroFluxSn50'] = [-99999]
					self.rv_data['spectroFluxSn50_err'] = [-99999]
					self.rv_data['naindex'] = [-99999]
					self.rv_data['naindex_err'] = [-99999]
					self.rv_data['snca2'] = [-99999]
					self.rv_data['snca2_err'] = [-99999]
					self.rv_data['public'] = [True]
					self.rv_data['pub_reference'] = ['Maxted']
					self.rv_data['pub_bibcode'] = ['']
					self.rv_data['sindex'] = [-99999]
					self.rv_data['sindex_err'] = [-99999]
					self.rv_data['haindex'] = [-99999]
					self.rv_data['haindex_err'] = [-99999]
					self.rv_data['drs_version'] = ['Pub']  
					self.rv_data['caindex'] = [-99999]
					self.rv_data['caindex_err'] = [-99999]
					self.rv_data['rjd_err'] = [-99999]
					rv_value = (float(values[1]))*1000.0                    
#					rv_value = float(values[4])*1.0
					self.rv_data['rv'] = [rv_value]
					self.rv_data['Vel(m/s)'] = [rv_value]
					self.rv_data['vel'] = [rv_value]                    
					rv_err_value = float(values[5])*1.0
					self.rv_data['rv_err'] = [rv_err_value]
					self.rv_data['ErrVel(m/s)'] = [rv_err_value]
					self.rv_data['errvel'] = [rv_err_value]
					self.rv_data['sn26'] = [-99999]                    
					self.rv_data['ccf_noise'] = [-99999]
					self.rv_data['ccf_noise_err'] = [-99999]
					self.rv_data['drift_used'] = [-99999]
					self.rv_data['drift_used_err'] = [-99999]
					self.rv_data['ccf_asym'] = [-99999]
					self.rv_data['ccf_asym_err'] = [-99999]
					self.rv_data['date_night'] = ['']  
					self.rv_data['raw_file'] = ['Maxted']
					self.rv_data['prog_id'] = [' ']
					self.rv_data['th_ar'] = [-99999]
					self.rv_data['th_ar1'] = [-99999]
					self.rv_data['th_ar2'] = [-99999]                    
					self.created=True
				self.rv_sum_discovery_hires = self.rv_sum_discovery_hires + float(values[4])
				self.rv_err_sum_discovery_hires = self.rv_err_sum_discovery_hires + float(values[5])
				self.rv_count_discovery_hires += 1
#				obsstring = "HARPS" + " & " + values[3] + " & " + values[4] + " & " + values[5] + " \\\\"               
#				with open(self.obsname, 'a') as observations_log:
#					observations_log.write('\n')                    
#					observations_log.write(obsstring)                 
		myfile21.close()    

 # Note that time is in rjd format not bjd
 # The file must be in .csv format and will have the Trifonov file format
 # file is loaded from the same directory as the notebook
 
	def _get_trifonov(self):
		print("enter HARPS03")
		print()        
		myfile1 = open ('HARPS_RVBank_v1.csv', 'r')
		for line in myfile1:
			values = line.split(',')
			object_name = values[0]
			i=0
			imax=0
			max_fuzzy=0
			while i < len(self.target_alias):
				String_Matched = fuzz.ratio(object_name, self.target_alias[i][0])
				if String_Matched > max_fuzzy:
					max_fuzzy=String_Matched
					imax=i
				i+=1
			if max_fuzzy > 92:
				line = "String Matched % "+str(max_fuzzy)+"target object"+object_name+ "target_max "+ self.target_alias[imax][0]
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)
				print(line)
				rjd_data = float(values[28]) - 2400000.0
				t_data = float(values[28]) - 2450000.0                
				if self.created:
					self.rv_data['rjd'] += [rjd_data]
					self.rv_data['bjd'] += [float(values[28])]
					self.rv_data['BJD'] += [float(values[28])]
					self.rv_data['t'] += [t_data]                     
					self.rv_data['berv'] += [float(values[29])]
					self.rv_data['berv_err'] += [-99999]
					self.rv_data['bispan'] += [float(values[25])]
					self.rv_data['bispan_err'] += [-99999] 
					self.rv_data['drift_noise'] += [float(values[31])]
					self.rv_data['drift_noise_err'] += [float(values[32])]
					self.rv_data['texp'] += [float(values[37])]
					self.rv_data['texp_err'] += [-99999]
					self.rv_data['cal_therror'] += [-99999]  
					self.rv_data['cal_therror_err'] += [-99999]
					self.rv_data['protm08'] += [-99999] 
					self.rv_data['protm08_err'] += [-99999] 
					self.rv_data['protn84'] += [-99999] 
					self.rv_data['protn84_err'] += [-99999] 
					self.rv_data['fwhm'] += [float(values[23])]  
					self.rv_data['fwhm_err'] += [-99999]
					self.rv_data['spectroFluxSn20'] += [-99999]
					self.rv_data['spectroFluxSn20_err'] += [-99999]
					self.rv_data['cal_thfile'] += [-99999] 
					self.rv_data['rhk_err'] += [-99999]
					self.rv_data['rhk'] += [self.rhk_value]
					self.rv_data['drs_qc'] += [True]
					self.rv_data['ins_name'] += ['HARPS03']
					self.rv_data['Telescope'] += ['HARPS03']                    
					self.rv_data['ins_mode'] += ['']
					self.rv_data['mask'] += [values[42]]
					self.rv_data['contrast'] += [float(values[24])]
					self.rv_data['contrast_err'] += [-99999]
					self.rv_data['spectroFluxSn50'] += [-99999]
					self.rv_data['spectroFluxSn50_err'] += [-99999]
					self.rv_data['naindex'] += [float(values[18])]
					self.rv_data['naindex_err'] += [float(values[19])]
					self.rv_data['snca2'] += [-99999]
					self.rv_data['snca2_err'] += [-99999]
					self.rv_data['public'] += [True]
					self.rv_data['pub_reference'] += ['Trifonov et al. 2020']
					self.rv_data['pub_bibcode'] += ['']
					self.rv_data['sindex'] += [-99999]
					self.rv_data['sindex_err'] += [-99999]
					self.rv_data['haindex'] += [float(values[16])]
					self.rv_data['haindex_err'] += [float(values[17])]
					self.rv_data['drs_version'] += ['Pub']  
					self.rv_data['caindex'] += [-99999]
					self.rv_data['caindex_err'] += [-99999]
					self.rv_data['rjd_err'] += [-99999]
					self.rv_data['rv'] += [float(values[4])]
					self.rv_data['Vel(m/s)'] += [float(values[4])]
					self.rv_data['vel'] += [float(values[4])]                    
					self.rv_data['rv_err'] += [float(values[5])]
					self.rv_data['ErrVel(m/s)'] += [float(values[5])]
					self.rv_data['errvel'] += [float(values[5])]
					self.rv_data['sn26'] += [-99999]                    
					self.rv_data['ccf_noise'] += [-99999]
					self.rv_data['ccf_noise_err'] += [-99999]
					self.rv_data['drift_used'] += [-99999]
					self.rv_data['drift_used_err'] += [-99999]
					self.rv_data['ccf_asym'] += [-99999]
					self.rv_data['ccf_asym_err'] += [-99999]
					self.rv_data['date_night'] += ['']  
					self.rv_data['raw_file'] += ['Trifonov et al. 2020']
					self.rv_data['prog_id'] += [' ']
					self.rv_data['th_ar'] += [-99999]
					self.rv_data['th_ar1'] += [-99999]
					self.rv_data['th_ar2'] += [-99999]                    
				else:
					self.rv_data['rjd'] = [rjd_data]
					self.rv_data['bjd'] = [float(values[28])]
					self.rv_data['BJD'] = [float(values[28])]
					self.rv_data['t'] = [t_data]                    
					self.rv_data['berv'] = [float(values[29])]
					self.rv_data['berv_err'] = [-99999]
					self.rv_data['bispan'] = [float(values[25])]
					self.rv_data['bispan_err'] = [-99999] 
					self.rv_data['drift_noise'] = [float(values[31])]
					self.rv_data['drift_noise_err'] = [float(values[32])]
					self.rv_data['texp'] = [float(values[37])]
					self.rv_data['texp_err'] = [-99999]
					self.rv_data['cal_therror'] = [-99999]  
					self.rv_data['cal_therror_err'] = [-99999]
					self.rv_data['protm08'] = [-99999] 
					self.rv_data['protm08_err'] = [-99999] 
					self.rv_data['protn84'] = [-99999] 
					self.rv_data['protn84_err'] = [-99999] 
					self.rv_data['fwhm'] = [float(values[23])]  
					self.rv_data['fwhm_err'] = [-99999]
					self.rv_data['spectroFluxSn20'] = [-99999]
					self.rv_data['spectroFluxSn20_err'] = [-99999]
					self.rv_data['cal_thfile'] = [-99999] 
					self.rv_data['rhk_err'] = [-99999]
					self.rv_data['rhk'] = [self.rhk_value]
					self.rv_data['drs_qc'] = [True]
					self.rv_data['ins_name'] = ['HARPS03']
					self.rv_data['Telescope'] = ['HARPS03']                    
					self.rv_data['ins_mode'] = ['']
					self.rv_data['mask'] = [values[42]]
					self.rv_data['contrast'] = [float(values[24])]
					self.rv_data['contrast_err'] = [-99999]
					self.rv_data['spectroFluxSn50'] = [-99999]
					self.rv_data['spectroFluxSn50_err'] = [-99999]
					self.rv_data['naindex'] = [float(values[18])]
					self.rv_data['naindex_err'] = [float(values[19])]
					self.rv_data['snca2'] = [-99999]
					self.rv_data['snca2_err'] = [-99999]
					self.rv_data['public'] = [True]
					self.rv_data['pub_reference'] = ['Trifonov et al. 2020']
					self.rv_data['pub_bibcode'] = ['']
					self.rv_data['sindex'] = [-99999]
					self.rv_data['sindex_err'] = [-99999]
					self.rv_data['haindex'] = [float(values[16])]
					self.rv_data['haindex_err'] = [float(values[17])]
					self.rv_data['drs_version'] = ['Pub']  
					self.rv_data['caindex'] = [-99999]
					self.rv_data['caindex_err'] = [-99999]
					self.rv_data['rjd_err'] = [-99999]
					self.rv_data['rv'] = [float(values[4])]
					self.rv_data['Vel(m/s)'] = [float(values[4])]
					self.rv_data['vel'] = [float(values[4])]                    
					self.rv_data['rv_err'] = [float(values[5])]
					self.rv_data['ErrVel(m/s)'] = [float(values[5])]
					self.rv_data['errvel'] = [float(values[5])]
					self.rv_data['sn26'] = [-99999]                    
					self.rv_data['ccf_noise'] = [-99999]
					self.rv_data['ccf_noise_err'] = [-99999]
					self.rv_data['drift_used'] = [-99999]
					self.rv_data['drift_used_err'] = [-99999]
					self.rv_data['ccf_asym'] = [-99999]
					self.rv_data['ccf_asym_err'] = [-99999]
					self.rv_data['date_night'] = ['']  
					self.rv_data['raw_file'] = ['Trifonov et al. 2020']
					self.rv_data['prog_id'] = [' ']
					self.rv_data['th_ar'] = [-99999]
					self.rv_data['th_ar1'] = [-99999]
					self.rv_data['th_ar2'] = [-99999]                    
					self.created=True
				print('trifonov date ',rjd_data,'rv  ',float(values[4]))
				print()
				self.rv_sum_trifonov = self.rv_sum_trifonov + float(values[4])
				self.rv_err_sum_trifonov = self.rv_err_sum_trifonov + float(values[5])
				self.rv_count_trifonov +=1
#				obsstring = "HARPS03" + " & " + values[28] + " & " + values[4] + " & " + values[5] + " \\\\"               
#				with open(self.obsname, 'a') as observations_log:
#					observations_log.write('\n')                    
#					observations_log.write(obsstring)                
		myfile1.close()

# hardener: add RHK and DRS_QC since they were not filled   

	def _get_butler_2017(self):
		print("enter HIRES2017")
		print()        
		file_name = 'HIRES_Keck_Precision_Radial_Velocity_Exoplanet_Survey.csv'
		myfile3 = open (file_name, 'r')
		for line in myfile3:
			values = line.split(',')
			object_name = values[0]
			i=0
			imax=0
			max_fuzzy=0
			try:
				alias_count = len(self.target_alias)
			except:
				alias_count = 0
			while i < alias_count:
				String_Matched = fuzz.ratio(object_name, self.target_alias[i][0])
				if String_Matched > max_fuzzy:
					max_fuzzy=String_Matched
					imax=i
				i+=1
			if max_fuzzy > 99:
				line = "String Matched % "+str(max_fuzzy)+"target object"+object_name+ "target_max "+ self.target_alias[imax][0]
				print(line)
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)
				rjd_data = float(values[1]) - 2400000.0
				t_data = float(values[1]) - 2450000.0                
				if self.created:
					self.rv_data['rjd'] += [rjd_data]
					self.rv_data['bjd'] += [float(values[1])]
					self.rv_data['BJD'] += [float(values[1])]
					self.rv_data['t'] += [t_data]                    
					self.rv_data['berv'] += [-99999]
					self.rv_data['berv_err'] += [-99999]
					self.rv_data['bispan'] += [-99999]
					self.rv_data['bispan_err'] += [-99999] 
					self.rv_data['drift_noise'] += [-99999]
					self.rv_data['drift_noise_err'] += [-99999]
					self.rv_data['texp'] += [-99999]
					self.rv_data['texp_err'] += [-99999]
					self.rv_data['cal_therror'] += [-99999]  
					self.rv_data['cal_therror_err'] += [-99999]
					self.rv_data['protm08'] += [-99999] 
					self.rv_data['protm08_err'] += [-99999] 
					self.rv_data['protn84'] += [-99999] 
					self.rv_data['protn84_err'] += [-99999] 
					self.rv_data['fwhm'] += [-99999]  
					self.rv_data['fwhm_err'] += [-99999]
					self.rv_data['spectroFluxSn20'] += [-99999]
					self.rv_data['spectroFluxSn20_err'] += [-99999]
					self.rv_data['cal_thfile'] += [-99999]
					self.rv_data['rhk_err'] += [-99999]
					self.rv_data['rhk'] += [self.rhk_value]
					self.rv_data['drs_qc'] += [True]
					self.rv_data['ins_name'] += ['HIRES']
					self.rv_data['Telescope'] += ['HIRES']                    
					self.rv_data['ins_mode'] += ['']
					self.rv_data['mask'] += ['']
					self.rv_data['contrast'] += [-99999]
					self.rv_data['contrast_err'] += [-99999]
					self.rv_data['spectroFluxSn50'] += [-99999]
					self.rv_data['spectroFluxSn50_err'] += [-99999]
					self.rv_data['naindex'] += [-99999]
					self.rv_data['naindex_err'] += [-99999]
					self.rv_data['snca2'] += [-99999]
					self.rv_data['snca2_err'] += [-99999]
					self.rv_data['public'] += [True]
					self.rv_data['pub_reference'] += ['Butler et al. 2017']
					self.rv_data['pub_bibcode'] += ['']
					self.rv_data['sindex'] += [float(values[4])]
					self.rv_data['sindex_err'] += [-99999]
					self.rv_data['haindex'] += [float(values[5])]
					self.rv_data['haindex_err'] += [-99999]
					self.rv_data['drs_version'] += ['Pub']  
					self.rv_data['caindex'] += [-99999]
					self.rv_data['caindex_err'] += [-99999]
					self.rv_data['rjd_err'] += [-99999]
					rv_value = (float(values[2])-self.gamma)*1000.0                        
					self.rv_data['rv'] += [rv_value]                    
					self.rv_data['Vel(m/s)'] += [rv_value]
					self.rv_data['vel'] += [rv_value]                    
					self.rv_data['rv_err'] += [float(values[3])]
					self.rv_data['ErrVel(m/s)'] += [float(values[3])]
					self.rv_data['errvel'] += [float(values[3])]
					self.rv_data['sn26'] += [-99999]                    
					self.rv_data['ccf_noise'] += [-99999]
					self.rv_data['ccf_noise_err'] += [-99999]
					self.rv_data['drift_used'] += [-99999]
					self.rv_data['drift_used_err'] += [-99999]
					self.rv_data['ccf_asym'] += [-99999]
					self.rv_data['ccf_asym_err'] += [-99999]
					self.rv_data['date_night'] += ['']  
					self.rv_data['raw_file'] += ['Butler et al. 2017']
					self.rv_data['prog_id'] += [' ']
					self.rv_data['th_ar'] += [-99999]
					self.rv_data['th_ar1'] += [-99999]
					self.rv_data['th_ar2'] += [-99999]                    
				else:
					self.rv_data['rjd'] = [rjd_data]
					self.rv_data['bjd'] = [float(values[1])]
					self.rv_data['BJD'] = [float(values[1])]
					self.rv_data['t'] = [t_data]                     
					self.rv_data['berv'] = [-99999]
					self.rv_data['berv_err'] = [-99999]
					self.rv_data['bispan'] = [-99999]
					self.rv_data['bispan_err'] = [-99999] 
					self.rv_data['drift_noise'] = [-99999]
					self.rv_data['drift_noise_err'] = [-99999]
					self.rv_data['texp'] = [-99999]
					self.rv_data['texp_err'] = [-99999]
					self.rv_data['cal_therror'] = [-99999]  
					self.rv_data['cal_therror_err'] = [-99999]
					self.rv_data['protm08'] = [-99999] 
					self.rv_data['protm08_err'] = [-99999] 
					self.rv_data['protn84'] = [-99999] 
					self.rv_data['protn84_err'] = [-99999] 
					self.rv_data['fwhm'] = [-99999]  
					self.rv_data['fwhm_err'] = [-99999]
					self.rv_data['spectroFluxSn20'] = [-99999]
					self.rv_data['spectroFluxSn20_err'] = [-99999]
					self.rv_data['cal_thfile'] = [-99999] 
					self.rv_data['rhk_err'] = [-99999]
					self.rv_data['rhk'] = [self.rhk_value]
					self.rv_data['drs_qc'] = [True]
					self.rv_data['ins_name'] = ['HIRES']
					self.rv_data['Telescope'] = ['HIRES']                    
					self.rv_data['ins_mode'] = ['']
					self.rv_data['mask'] = ['']
					self.rv_data['contrast'] = [-99999]
					self.rv_data['contrast_err'] = [-99999]
					self.rv_data['spectroFluxSn50'] = [-99999]
					self.rv_data['spectroFluxSn50_err'] = [-99999]
					self.rv_data['naindex'] = [-99999]
					self.rv_data['naindex_err'] = [-99999]
					self.rv_data['snca2'] = [-99999]
					self.rv_data['snca2_err'] = [-99999]
					self.rv_data['public'] = [True]
					self.rv_data['pub_reference'] = ['Butler et al. 2017']
					self.rv_data['pub_bibcode'] = ['']
					self.rv_data['sindex'] = [float(values[4])]
					self.rv_data['sindex_err'] = [-99999]
					self.rv_data['haindex'] = [float(values[5])]
					self.rv_data['haindex_err'] = [-99999]
					self.rv_data['drs_version'] = ['Pub']  
					self.rv_data['caindex'] = [-99999]
					self.rv_data['caindex_err'] = [-99999]
					self.rv_data['rjd_err'] = [-99999]
					self.rv_data['rv'] = [float(values[2])]
					rv_value = (float(values[2])-self.gamma)*1000.0                        
					self.rv_data['rv'] = [rv_value]                    
					self.rv_data['Vel(m/s)'] = [rv_value]
					self.rv_data['vel'] = [rv_value]                     
					self.rv_data['rv_err'] = [float(values[3])]
					self.rv_data['ErrVel(m/s)'] = [float(values[3])]
					self.rv_data['errvel'] = [float(values[3])]
					self.rv_data['sn26'] = [-99999]                    
					self.rv_data['ccf_noise'] = [-99999]
					self.rv_data['ccf_noise_err'] = [-99999]
					self.rv_data['drift_used'] = [-99999]
					self.rv_data['drift_used_err'] = [-99999]
					self.rv_data['ccf_asym'] = [-99999]
					self.rv_data['ccf_asym_err'] = [-99999]
					self.rv_data['date_night'] = ['']  
					self.rv_data['raw_file'] = ['Butler et al. 2017']
					self.rv_data['prog_id'] = [' ']
					self.rv_data['th_ar'] = [-99999]
					self.rv_data['th_ar1'] = [-99999]
					self.rv_data['th_ar2'] = [-99999]                    
					self.created=True
				self.rv_sum_butler_2017 = self.rv_sum_butler_2017 + float(values[2])
				self.rv_err_sum_butler_2017 = self.rv_err_sum_butler_2017 + float(values[3])
				self.rv_count_butler_2017 +=1
#				obsstring = "HIRES2017" + " & " + values[1] + " & " + values[2] + " & " + values[3] + " \\\\"               
#				with open(self.obsname, 'a') as observations_log:
#					observations_log.write('\n')                    
#					observations_log.write(obsstring)
		myfile3.close()

	def _get_butler_2020(self):
		print("enter HIRES2020 Fred's Counts")
		print()        
		file_name = 'KeckVels_Coordinates_Aliases_with_RV_counts.csv'
		myfilefred = open (file_name, 'r')
		header_lines = 1
		count_lines = 0
		for line in myfilefred:
			if count_lines < header_lines:
				count_lines+=1
				continue        
			values = line.split(',')
			object_name = values[0]                                
			ravariation = 2/3600
			try:
				ravalue = float(values[4])
			except:
				ravalue = 0.0                
			try:
				decvalue = float(values[5])
			except:
				decvalue = 0.0
			if abs(ravalue-self.radegrees) < ravariation and abs(decvalue-self.decdegrees) < ravariation:
				print('target ', object_name)
				print("RA ",ravalue, " DEC ",decvalue)
				print("RA variation",abs(ravalue-self.radegrees)," DEC variation ",abs(decvalue-self.decdegrees), " variation allowed ",ravariation)
				self.rv_sum_butler_2020 = self.rv_sum_butler_2020 + 0
				self.rv_err_sum_butler_2020 = self.rv_err_sum_butler_2020 + 0
				try:
					rvcount = float(values[1])
				except:
					rvcount = 0.0
				self.rv_count_butler_2020 += rvcount
		myfilefred.close()                
        
	def _get_butler_2020a(self):
		print("enter get_butler_2020")
		print()        
		for filename in os.listdir('keck_vels'):
			if fnmatch.fnmatch(filename, '*.vels'):
				target = filename.split('_',1)
				object_name = target[0]
#				print("object name ", object_name," filename ",filename)                
				i=0
				try:
					nospace_count = len(self.target_alias_nospace)
				except:
					nospace_count = 0
				while i < nospace_count:
#					if object_name == self.target or object_name == self.target_alias[i][0] or object_name == self.target_alias_nospace[i][0]:            
					if object_name == self.target_alias[i][0] or object_name == self.target_alias_nospace[i][0]:               
						fileopen = 'keck_vels/'+filename
						print("file name ",filename)                        
						myfile5 = open (fileopen, 'r')
						for line in myfile5:
							values = line.split()                    
							rjd_data = float(values[0]) - 2400000.0
							t_data = float(values[0]) - 2450000.0                    
							if self.created:
								self.rv_data['rjd'] += [rjd_data]
								self.rv_data['bjd'] += [float(values[0])]
								self.rv_data['BJD'] += [float(values[0])]    
								self.rv_data['t'] += [t_data]                                                
								self.rv_data['berv'] += [-99999]                       
								self.rv_data['berv_err'] += [-99999]
								self.rv_data['bispan'] += [-99999]                        
								self.rv_data['bispan_err'] += [-99999]                        
								self.rv_data['drift_noise'] += [-99999]
								self.rv_data['drift_noise_err'] += [-99999]
								self.rv_data['texp'] += [-99999]
								self.rv_data['texp_err'] += [-99999]
								self.rv_data['cal_therror'] += [-99999] 
								self.rv_data['cal_therror_err'] += [-99999]
								self.rv_data['protm08'] += [-99999] 
								self.rv_data['protm08_err'] += [-99999] 
								self.rv_data['protn84'] += [-99999] 
								self.rv_data['protn84_err'] += [-99999] 
								self.rv_data['fwhm'] += [-99999]                        
								self.rv_data['fwhm_err'] += [-99999]                        
								self.rv_data['spectroFluxSn20'] += [-99999]
								self.rv_data['spectroFluxSn20_err'] += [-99999]
								self.rv_data['cal_thfile'] += [-99999] 
								self.rv_data['rhk_err'] += [-99999]
								self.rv_data['rhk'] += [self.rhk_value]
								self.rv_data['drs_qc'] += [True]
								self.rv_data['ins_name'] += ['HIRES2020']
								self.rv_data['Telescope'] += ['HIRES2020']                        
								self.rv_data['ins_mode'] += ['']
								self.rv_data['mask'] += ['']
								self.rv_data['contrast'] += [-99999]
								self.rv_data['contrast_err'] += [-99999]                        
								self.rv_data['spectroFluxSn50'] += [-99999]
								self.rv_data['spectroFluxSn50_err'] += [-99999]
								self.rv_data['naindex'] += [-99999]
								self.rv_data['naindex_err'] += [-99999]
								self.rv_data['snca2'] += [-99999]
								self.rv_data['snca2_err'] += [-99999]
								self.rv_data['public'] += [True]
								self.rv_data['pub_reference'] += ['none']
								self.rv_data['pub_bibcode'] += ['']
								self.rv_data['sindex'] += [-99999]
								self.rv_data['sindex_err'] += [-99999]
								self.rv_data['haindex'] += [-99999]
								self.rv_data['haindex_err'] += [-99999]
								self.rv_data['drs_version'] += ['Pub']  
								self.rv_data['caindex'] += [-99999]
								self.rv_data['caindex_err'] += [-99999]
								self.rv_data['rjd_err'] += [-99999]
								rv_value = float(values[1])*1000.0
								self.rv_data['rv'] += [rv_value]
								self.rv_data['Vel(m/s)'] += [rv_value]
								self.rv_data['vel'] += [rv_value]                         
								rv_err_value = float(values[2])*1000.0
								self.rv_data['rv_err'] += [rv_err_value]
								self.rv_data['ErrVel(m/s)'] += [rv_err_value]
								self.rv_data['errvel'] += [rv_err_value]
								self.rv_data['sn26'] += [-99999]                        
								self.rv_data['ccf_noise'] += [-99999]
								self.rv_data['ccf_noise_err'] += [-99999]
								self.rv_data['drift_used'] += [-99999]
								self.rv_data['drift_used_err'] += [-99999]
								self.rv_data['ccf_asym'] += [-99999]
								self.rv_data['ccf_asym_err'] += [-99999]
								self.rv_data['date_night'] += ['']  
								self.rv_data['raw_file'] += ['none']
								self.rv_data['prog_id'] += [' ']
								self.rv_data['th_ar'] += [-99999]
								self.rv_data['th_ar1'] += [-99999]
								self.rv_data['th_ar2'] += [-99999]                                
							else:
								self.rv_data['rjd'] = [rjd_data]
								self.rv_data['bjd'] = [float(values[0])]
								self.rv_data['BJD'] = [float(values[0])]
								self.rv_data['t'] = [t_data]                                                
								self.rv_data['berv'] = [-99999]
								self.rv_data['berv_err'] = [-99999]
								self.rv_data['berv_err'] = [-99999]
								self.rv_data['bispan'] = [-99999]                        
								self.rv_data['bispan_err'] = [-99999] 
								self.rv_data['drift_noise'] = [-99999]
								self.rv_data['drift_noise_err'] = [-99999]
								self.rv_data['texp'] = [-99999]
								self.rv_data['texp_err'] = [-99999]
								self.rv_data['cal_therror'] = [-99999]  
								self.rv_data['cal_therror_err'] = [-99999]
								self.rv_data['protm08'] = [-99999] 
								self.rv_data['protm08_err'] = [-99999] 
								self.rv_data['protn84'] = [-99999] 
								self.rv_data['protn84_err'] = [-99999] 
								self.rv_data['fwhm'] = [-99999]                        
								self.rv_data['fwhm_err'] = [-99999]
								self.rv_data['spectroFluxSn20'] = [-99999]
								self.rv_data['spectroFluxSn20_err'] = [-99999]
								self.rv_data['cal_thfile'] = [-99999] 
								self.rv_data['rhk_err'] = [-99999]
								self.rv_data['rhk'] = [self.rhk_value]
								self.rv_data['drs_qc'] = [True]
								self.rv_data['ins_name'] = ['HIRES2020']
								self.rv_data['Telescope'] = ['HIRES2020']                        
								self.rv_data['ins_mode'] = ['']
								self.rv_data['mask'] = ['']
								self.rv_data['contrast'] = [-99999]
								self.rv_data['contrast_err'] = [-99999]
								self.rv_data['spectroFluxSn50'] = [-99999]
								self.rv_data['spectroFluxSn50_err'] = [-99999]
								self.rv_data['naindex'] = [-99999]
								self.rv_data['naindex_err'] = [-99999]
								self.rv_data['snca2'] = [-99999]
								self.rv_data['snca2_err'] = [-99999]
								self.rv_data['public'] = [True]
								self.rv_data['pub_reference'] = ['none']
								self.rv_data['pub_bibcode'] = ['']
								self.rv_data['sindex'] = [-99999]
								self.rv_data['sindex_err'] = [-99999]
								self.rv_data['haindex'] = [-99999]
								self.rv_data['haindex_err'] = [-99999]
								self.rv_data['drs_version'] = ['Pub']  
								self.rv_data['caindex'] = [-99999]
								self.rv_data['caindex_err'] = [-99999]
								self.rv_data['rjd_err'] = [-99999]
								rv_value = float(values[1])*1000.0
								self.rv_data['rv'] = [rv_value]
								self.rv_data['Vel(m/s)'] = [rv_value]
								self.rv_data['vel'] = [rv_value]                        
								rv_err_value = float(values[2])*1000.0
								self.rv_data['rv_err'] = [rv_err_value]
								self.rv_data['ErrVel(m/s)'] = [rv_err_value]
								self.rv_data['errvel'] = [rv_err_value]
								self.rv_data['sn26'] = [-99999]                        
								self.rv_data['ccf_noise'] = [-99999]
								self.rv_data['ccf_noise_err'] = [-99999]
								self.rv_data['drift_used'] = [-99999]
								self.rv_data['drift_used_err'] = [-99999]
								self.rv_data['ccf_asym'] = [-99999]
								self.rv_data['ccf_asym_err'] = [-99999]
								self.rv_data['date_night'] = ['']  
								self.rv_data['raw_file'] = ['none']
								self.rv_data['prog_id'] = [' ']
								self.rv_data['th_ar'] = [-99999]
								self.rv_data['th_ar1'] = [-99999]
								self.rv_data['th_ar2'] = [-99999]                                 
								self.created=True
							self.rv_sum_butler_2020 = self.rv_sum_butler_2020 + rv_value
							self.rv_err_sum_butler_2020 = self.rv_err_sum_butler_2020 + rv_err_value
							self.rv_count_butler_2020 +=1
#							obsstring = "BUTLER 2020" + " & " + values[0] + " & " + values[1] + " & " + values[2] + " \\\\"   
#							with open(self.obsname, 'a') as observations_log:
#								observations_log.write('\n')                        
#								observations_log.write(obsstring)
						myfile5.close() 
					i+=1        

	def _get_eso_harps(self):
		print("enter HARPS15")
		print()        
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
				if abs(HARPS_Table[k][35]-self.radegrees) < ravariation and abs(HARPS_Table[k][36]-self.decdegrees) < ravariation and bjd_data >= 2457387.00000:
					print('target ', object_name)
					print("RA ",HARPS_Table[k][35], " DEC ",HARPS_Table[k][36]," Program ID ",HARPS_Table[k][2])               
					rjd_data = bjd_data - 2400000.0
					t_data = bjd_data - 2450000.0                    
					if self.created:
						self.rv_data['rjd'] += [rjd_data]
						self.rv_data['bjd'] += [float(HARPS_Table[k][13])]
						self.rv_data['BJD'] += [float(HARPS_Table[k][13])]
						self.rv_data['t'] += [t_data]                       
						berv_value = float(HARPS_Table[k][16])*1000.0
						self.rv_data['berv'] += [berv_value]
						self.rv_data['berv_err'] += [-99999]
						self.rv_data['bispan'] += [-99999]
						self.rv_data['bispan_err'] += [-99999] 
						self.rv_data['drift_noise'] += [-99999]
						self.rv_data['drift_noise_err'] += [-99999]
						self.rv_data['texp'] += [-99999]
						self.rv_data['texp_err'] += [-99999]
						self.rv_data['cal_therror'] += [-99999]  
						self.rv_data['cal_therror_err'] += [-99999]
						self.rv_data['protm08'] += [-99999] 
						self.rv_data['protm08_err'] += [-99999] 
						self.rv_data['protn84'] += [-99999] 
						self.rv_data['protn84_err'] += [-99999] 
						fwhm_value = float(HARPS_Table[k][23])*1000.0 
						self.rv_data['fwhm'] += [fwhm_value]
						self.rv_data['fwhm_err'] += [-99999]
						self.rv_data['spectroFluxSn20'] += [-99999]
						self.rv_data['spectroFluxSn20_err'] += [-99999]
						self.rv_data['cal_thfile'] += [-99999] 
						self.rv_data['rhk_err'] += [-99999]
						self.rv_data['rhk'] += [self.rhk_value]
						self.rv_data['drs_qc'] += [True]
						self.rv_data['ins_name'] += ['HARPS15']
						self.rv_data['Telescope'] += ['HARPS15']                         
						self.rv_data['ins_mode'] += ['']
						self.rv_data['mask'] += ['']
						contrast_value = float(HARPS_Table[k][22])                     
						self.rv_data['contrast'] += [contrast_value]
						self.rv_data['contrast_err'] += [-99999]
						self.rv_data['spectroFluxSn50'] += [-99999]
						self.rv_data['spectroFluxSn50_err'] += [-99999]
						self.rv_data['naindex'] += [-99999]
						self.rv_data['naindex_err'] += [-99999]
						self.rv_data['snca2'] += [-99999]
						self.rv_data['snca2_err'] += [-99999]
						self.rv_data['public'] += [True]
						self.rv_data['pub_reference'] += ['none']
						self.rv_data['pub_bibcode'] += ['']
						self.rv_data['sindex'] += [-99999]
						self.rv_data['sindex_err'] += [-99999]
						self.rv_data['haindex'] += [-99999]
						self.rv_data['haindex_err'] += [-99999]
						self.rv_data['drs_version'] += ['Pub']  
						self.rv_data['caindex'] += [-99999]
						self.rv_data['caindex_err'] += [-99999]
						self.rv_data['rjd_err'] += [-99999]
#						rv_value = float(HARPS_Table[k][19])*1000.0
						rv_value = (float(HARPS_Table[k][19])-self.gamma)*1000.0                        
						self.rv_data['rv'] += [rv_value]
						self.rv_data['Vel(m/s)'] += [rv_value]
						self.rv_data['vel'] += [rv_value]                         
						rv_err_value = float(HARPS_Table[k][24])*1000.0
#						rv_err_value = (float(HARPS_Table[k][24])-self.gamma)*1000.0                    
						self.rv_data['rv_err'] += [rv_err_value]
						self.rv_data['ErrVel(m/s)'] += [rv_err_value]
						self.rv_data['errvel'] += [rv_err_value]
						self.rv_data['sn26'] += [-99999]                        
						self.rv_data['ccf_noise'] += [-99999]
						self.rv_data['ccf_noise_err'] += [-99999]
						self.rv_data['drift_used'] += [-99999]
						self.rv_data['drift_used_err'] += [-99999]
						self.rv_data['ccf_asym'] += [-99999]
						self.rv_data['ccf_asym_err'] += [-99999]
						self.rv_data['date_night'] += ['']  
						self.rv_data['raw_file'] += [HARPS_Table[k][1]]
						self.rv_data['prog_id'] += [HARPS_Table[k][2]]
						self.rv_data['th_ar'] += [-99999]
						self.rv_data['th_ar1'] += [-99999]
						self.rv_data['th_ar2'] += [-99999]                        
					else:
						self.rv_data['rjd'] = [rjd_data]
						self.rv_data['bjd'] = [float(HARPS_Table[k][13])]
						self.rv_data['BJD'] = [float(HARPS_Table[k][13])]
						self.rv_data['t'] = [t_data]                        
						berv_value = float(HARPS_Table[k][16])*1000.0
						self.rv_data['berv'] = [berv_value]
						self.rv_data['berv_err'] = [-99999]
						self.rv_data['bispan'] = [-99999]
						self.rv_data['bispan_err'] = [-99999] 
						self.rv_data['drift_noise'] = [-99999]
						self.rv_data['drift_noise_err'] = [-99999]
						self.rv_data['texp'] = [-99999]
						self.rv_data['texp_err'] = [-99999]
						self.rv_data['cal_therror'] = [-99999] 
						self.rv_data['cal_therror_err'] = [-99999]
						self.rv_data['protm08'] = [-99999] 
						self.rv_data['protm08_err'] = [-99999] 
						self.rv_data['protn84'] = [-99999] 
						self.rv_data['protn84_err'] = [-99999] 
						fwhm_value = float(HARPS_Table[k][23])*1000.0 
						self.rv_data['fwhm'] = [fwhm_value]  
						self.rv_data['fwhm_err'] = [-99999]
						self.rv_data['spectroFluxSn20'] = [-99999]
						self.rv_data['spectroFluxSn20_err'] = [-99999]
						self.rv_data['cal_thfile'] = [-99999] 
						self.rv_data['rhk_err'] = [-99999]
						self.rv_data['rhk'] = [self.rhk_value]
						self.rv_data['drs_qc'] = [True]
						self.rv_data['ins_name'] = ['HARPS15']
						self.rv_data['Telescope'] = ['HARPS15']                        
						self.rv_data['ins_mode'] = ['']
						self.rv_data['mask'] = ['']
						contrast_value = float(HARPS_Table[k][22])*1000.0 
						self.rv_data['contrast'] = [contrast_value]
						self.rv_data['contrast_err'] = [-99999]
						self.rv_data['spectroFluxSn50'] = [-99999]
						self.rv_data['spectroFluxSn50_err'] = [-99999]
						self.rv_data['naindex'] = [-99999]
						self.rv_data['naindex_err'] = [-99999]
						self.rv_data['snca2'] = [-99999]
						self.rv_data['snca2_err'] = [-99999]
						self.rv_data['public'] = [True]
						self.rv_data['pub_reference'] = ['none']
						self.rv_data['pub_bibcode'] = ['']
						self.rv_data['sindex'] = [-99999]
						self.rv_data['sindex_err'] = [-99999]
						self.rv_data['haindex'] = [-99999]
						self.rv_data['haindex_err'] = [-99999]
						self.rv_data['drs_version'] = ['Pub']  
						self.rv_data['caindex'] = [-99999]
						self.rv_data['caindex_err'] = [-99999]
						self.rv_data['rjd_err'] = [-99999]
#						rv_value = float(HARPS_Table[k][19])*1000.0
						rv_value = (float(HARPS_Table[k][19])-self.gamma)*1000.0                        
						self.rv_data['rv'] = [rv_value]
						self.rv_data['Vel(m/s)'] = [rv_value]
						self.rv_data['vel'] = [rv_value]                        
						rv_err_value = float(HARPS_Table[k][24])*1000.0
#						rv_err_value = (float(HARPS_Table[k][24])-self.gamma)*1000.0                    
						self.rv_data['rv_err'] = [rv_err_value]
						self.rv_data['ErrVel(m/s)'] = [rv_err_value]
						self.rv_data['errvel'] = [rv_err_value]
						self.rv_data['sn26'] = [-99999]                        
						self.rv_data['ccf_noise'] = [-99999]
						self.rv_data['ccf_noise_err'] = [-99999]
						self.rv_data['drift_used'] = [-99999]
						self.rv_data['drift_used_err'] = [-99999]
						self.rv_data['ccf_asym'] = [-99999]
						self.rv_data['ccf_asym_err'] = [-99999]
						self.rv_data['date_night'] = ['']  
						self.rv_data['raw_file'] = [HARPS_Table[k][1]]
						self.rv_data['prog_id'] = [HARPS_Table[k][2]]
						self.rv_data['th_ar'] = [-99999]
						self.rv_data['th_ar1'] = [-99999]
						self.rv_data['th_ar2'] = [-99999]
						self.created=True                        
					self.rv_sum_eso_harps = self.rv_sum_eso_harps + rv_value
					self.rv_err_sum_eso_harps = self.rv_err_sum_eso_harps + rv_err_value
					self.rv_count_eso_harps +=1
					print("HARPS15"," & ",bjd_data," & ",rv_value," & ",rv_err_value," \\\\")
					print("HARPS15 col 6 object ",HARPS_Table[k][6])                   
					print("HARPS15 col 16 berv ", HARPS_Table[k][16]," col 18 Bary RV corr ",HARPS_Table[k][18],"col 19 Bary RV no corr ",HARPS_Table[k][19]," col 20 BIS rv ",HARPS_Table[k][19])
					print("HARPS15 col 31 simbad ID ", HARPS_Table[k][31])                    
#					obsstring = "HARPS15" + " & " + str(bjd_data) + " & " + str(rv_value) + " & " + str(rv_err_value) + " \\\\"             
#					with open(self.obsname, 'a') as observations_log:
#						observations_log.write('\n')                        
#						observations_log.write(obsstring)                    
				k+=1

	def _get_eso_espresso(self):
		print("enter eso espresso")
		print()
		for filename in os.listdir('ESPRESSO_Full_Data_Set/archive'):
			if fnmatch.fnmatch(filename, '*.fits'):
				hdul = fits.open('ESPRESSO_Full_Data_Set/archive/'+filename)
				object_name = hdul[0].header['HIERARCH ESO OBS TARG NAME']  
				ravariation = 2/3600            
				if abs(hdul[0].header['RA']-self.radegrees) < ravariation and abs(hdul[0].header['DEC']-self.decdegrees) < ravariation:
					print('target ', object_name)
					print("RA ",hdul[0].header['RA'], " DEC ",hdul[0].header['DEC'])
					bjd_data = float(hdul[0].header['HIERARCH ESO QC BJD'])
					rjd_data = bjd_data - 2400000.0                    
					t_data = bjd_data - 2450000.0                    
					if self.created:
						self.rv_data['rjd'] += [rjd_data]
						self.rv_data['bjd'] += [bjd_data]
						self.rv_data['BJD'] += [bjd_data]    
						self.rv_data['t'] += [t_data]                         
						berv_value = float(hdul[0].header['HIERARCH ESO QC BERV'])*1000.0                       
						self.rv_data['berv'] += [berv_value]                       
						self.rv_data['berv_err'] += [-99999]
						bis_value = float(hdul[0].header['HIERARCH ESO QC CCF BIS SPAN'])*1000.0
						self.rv_data['bispan'] += [bis_value]                        
						bis_err_value = float(hdul[0].header['HIERARCH ESO QC CCF BIS SPAN ERROR'])*1000.0
						self.rv_data['bispan_err'] += [bis_err_value]                        
						self.rv_data['drift_noise'] += [-99999]
						self.rv_data['drift_noise_err'] += [-99999]
						self.rv_data['texp'] += [-99999]
						self.rv_data['texp_err'] += [-99999]
						self.rv_data['cal_therror'] += [-99999] 
						self.rv_data['cal_therror_err'] += [-99999]
						self.rv_data['protm08'] += [-99999] 
						self.rv_data['protm08_err'] += [-99999] 
						self.rv_data['protn84'] += [-99999] 
						self.rv_data['protn84_err'] += [-99999] 
						fwhm_value = float(hdul[0].header['HIERARCH ESO QC CCF FWHM'])*1000.0
						self.rv_data['fwhm'] += [fwhm_value]                        
						fwhm_err_value = float(hdul[0].header['HIERARCH ESO QC CCF FWHM ERROR'])*1000.0
						self.rv_data['fwhm_err'] += [fwhm_err_value]                        
						self.rv_data['spectroFluxSn20'] += [-99999]
						self.rv_data['spectroFluxSn20_err'] += [-99999]
						self.rv_data['cal_thfile'] += [-99999] 
						self.rv_data['rhk_err'] += [-99999]
						self.rv_data['rhk'] += [self.rhk_value]
						self.rv_data['drs_qc'] += [True]
						self.rv_data['ins_name'] += ['ESPRESSO']
						self.rv_data['Telescope'] += ['ESPRESSO']                        
						self.rv_data['ins_mode'] += ['']
						self.rv_data['mask'] += ['']
						contrast_value = float(hdul[0].header['HIERARCH ESO QC CCF CONTRAST'])*1000.0
						self.rv_data['contrast'] += [contrast_value]
						contrast_err_value = float(hdul[0].header['HIERARCH ESO QC CCF CONTRAST ERROR'])*1000.0
						self.rv_data['contrast_err'] += [contrast_err_value]                        
						self.rv_data['spectroFluxSn50'] += [-99999]
						self.rv_data['spectroFluxSn50_err'] += [-99999]
						self.rv_data['naindex'] += [-99999]
						self.rv_data['naindex_err'] += [-99999]
						self.rv_data['snca2'] += [-99999]
						self.rv_data['snca2_err'] += [-99999]
						self.rv_data['public'] += [True]
						self.rv_data['pub_reference'] += ['none']
						self.rv_data['pub_bibcode'] += ['']
						self.rv_data['sindex'] += [-99999]
						self.rv_data['sindex_err'] += [-99999]
						self.rv_data['haindex'] += [-99999]
						self.rv_data['haindex_err'] += [-99999]
						self.rv_data['drs_version'] += ['Pub']  
						self.rv_data['caindex'] += [-99999]
						self.rv_data['caindex_err'] += [-99999]
						self.rv_data['rjd_err'] += [-99999]
						rv_value = float(hdul[0].header['HIERARCH ESO QC CCF RV'])*1000.0
						self.rv_data['rv'] += [rv_value]
						self.rv_data['Vel(m/s)'] += [rv_value]
						self.rv_data['vel'] += [rv_value]                         
						rv_err_value = float(hdul[0].header['HIERARCH ESO QC CCF RV ERROR'])*1000.0
						self.rv_data['rv_err'] += [rv_err_value]
						self.rv_data['ErrVel(m/s)'] += [rv_err_value]
						self.rv_data['errvel'] += [rv_err_value]
						self.rv_data['sn26'] += [-99999]                        
						self.rv_data['ccf_noise'] += [-99999]
						self.rv_data['ccf_noise_err'] += [-99999]
						self.rv_data['drift_used'] += [-99999]
						self.rv_data['drift_used_err'] += [-99999]
						self.rv_data['ccf_asym'] += [-99999]
						self.rv_data['ccf_asym_err'] += [-99999]
						self.rv_data['date_night'] += ['']  
						self.rv_data['raw_file'] += ['none']
						self.rv_data['prog_id'] += [' ']
						self.rv_data['th_ar'] += [-99999]
						self.rv_data['th_ar1'] += [-99999]
						self.rv_data['th_ar2'] += [-99999]                         
					else:
						self.rv_data['rjd'] = [rjd_data]
						self.rv_data['bjd'] = [bjd_data]
						self.rv_data['BJD'] = [bjd_data]
						self.rv_data['t'] = [t_data]                        
						berv_value = float(hdul[0].header['HIERARCH ESO QC BERV'])*1000.0                        
						self.rv_data['berv'] = [berv_value]
						self.rv_data['berv_err'] = [-99999]
						bis_value = float(hdul[0].header['HIERARCH ESO QC CCF BIS SPAN'])*1000.0
						self.rv_data['bispan'] = [bis_value]                        
						bis_err_value = float(hdul[0].header['HIERARCH ESO QC CCF BIS SPAN ERROR'])*1000.0
						self.rv_data['bispan_err'] = [bis_err_value] 
						self.rv_data['drift_noise'] = [-99999]
						self.rv_data['drift_noise_err'] = [-99999]
						self.rv_data['texp'] = [-99999]
						self.rv_data['texp_err'] = [-99999]
						self.rv_data['cal_therror'] = [-99999]  
						self.rv_data['cal_therror_err'] = [-99999]
						self.rv_data['protm08'] = [-99999] 
						self.rv_data['protm08_err'] = [-99999] 
						self.rv_data['protn84'] = [-99999] 
						self.rv_data['protn84_err'] = [-99999] 
						fwhm_value = float(hdul[0].header['HIERARCH ESO QC CCF FWHM'])*1000.0
						self.rv_data['fwhm'] = [fwhm_value]                        
						fwhm_err_value = float(hdul[0].header['HIERARCH ESO QC CCF FWHM ERROR'])*1000.0
						self.rv_data['fwhm_err'] = [fwhm_err_value]
						self.rv_data['spectroFluxSn20'] = [-99999]
						self.rv_data['spectroFluxSn20_err'] = [-99999]
						self.rv_data['cal_thfile'] = [-99999] 
						self.rv_data['rhk_err'] = [-99999]
						self.rv_data['rhk'] = [self.rhk_value]
						self.rv_data['drs_qc'] = [True]
						self.rv_data['ins_name'] = ['ESPRESSO']
						self.rv_data['Telescope'] = ['ESPRESSO']                        
						self.rv_data['ins_mode'] = ['']
						self.rv_data['mask'] = ['']
						contrast_value = float(hdul[0].header['HIERARCH ESO QC CCF CONTRAST'])*1000.0
						self.rv_data['contrast'] = [contrast_value]
						contrast_err_value = float(hdul[0].header['HIERARCH ESO QC CCF CONTRAST ERROR'])*1000.0
						self.rv_data['contrast_err'] = [contrast_err_value]
						self.rv_data['spectroFluxSn50'] = [-99999]
						self.rv_data['spectroFluxSn50_err'] = [-99999]
						self.rv_data['naindex'] = [-99999]
						self.rv_data['naindex_err'] = [-99999]
						self.rv_data['snca2'] = [-99999]
						self.rv_data['snca2_err'] = [-99999]
						self.rv_data['public'] = [True]
						self.rv_data['pub_reference'] = ['none']
						self.rv_data['pub_bibcode'] = ['']
						self.rv_data['sindex'] = [-99999]
						self.rv_data['sindex_err'] = [-99999]
						self.rv_data['haindex'] = [-99999]
						self.rv_data['haindex_err'] = [-99999]
						self.rv_data['drs_version'] = ['Pub']  
						self.rv_data['caindex'] = [-99999]
						self.rv_data['caindex_err'] = [-99999]
						self.rv_data['rjd_err'] = [-99999]
						rv_value = float(hdul[0].header['ESO QC CCF RV'])*1000.0
						self.rv_data['rv'] = [rv_value]
						self.rv_data['Vel(m/s)'] = [rv_value]
						self.rv_data['vel'] = [rv_value]                        
						rv_err_value = float(hdul[0].header['ESO QC CCF RV ERROR'])*1000.0
						self.rv_data['rv_err'] = [rv_err_value]
						self.rv_data['ErrVel(m/s)'] = [rv_err_value]
						self.rv_data['errvel'] = [rv_err_value]
						self.rv_data['sn26'] = [-99999]                        
						self.rv_data['ccf_noise'] = [-99999]
						self.rv_data['ccf_noise_err'] = [-99999]
						self.rv_data['drift_used'] = [-99999]
						self.rv_data['drift_used_err'] = [-99999]
						self.rv_data['ccf_asym'] = [-99999]
						self.rv_data['ccf_asym_err'] = [-99999]
						self.rv_data['date_night'] = ['']  
						self.rv_data['raw_file'] = ['none']
						self.rv_data['prog_id'] = [' ']
						self.rv_data['th_ar'] = [-99999]
						self.rv_data['th_ar1'] = [-99999]
						self.rv_data['th_ar2'] = [-99999]                        
						self.created=True
					self.rv_sum_espresso = self.rv_sum_espresso + rv_value
					self.rv_err_sum_espresso = self.rv_err_sum_espresso + rv_err_value
					self.rv_count_espresso +=1
#					obsstring = "ESPRESSO" + " & " + str(bjd_data) + " & " + str(rv_value) + " & " + str(rv_err_value) + " \\\\"             
#					with open(self.obsname, 'a') as observations_log:
#						observations_log.write('\n')                        
#						observations_log.write(obsstring)
				hdul.close 

#stubs for future development 
	def _get_weiss(self):
		myfile9 = ''

	def _get_carmenes(self):
		myfile5 = ''

	def _get_sophie(self):
		print("enter SOPHIE")
		print()        
		myfile4 = open ('sophiecc_1686978863.txt', 'r')
		rjd_data = 0.0
		header_lines = 38
		count_lines = 0
		for line in myfile4:
			if count_lines < header_lines:
				count_lines+=1
				continue
			values = line.split(',')
			colcount=len(values)
			if colcount == 29:
				if values[28] != '\n' and values[1] != "SUN" and values[1] != "MOON" and values[16] != '' and values[17] != '':
					object_name = values[24]                                
					ravariation = 2250/3600
					if abs(float(values[16])-self.radegrees) < ravariation and abs(float(values[17])-self.decdegrees) < ravariation:
						print('target ', object_name)
						print("RA ",values[16], " DEC ",values[17])
						print("RA variation",abs(float(values[16])-self.radegrees)," DEC variation ",abs(float(values[17])-self.decdegrees), " variation allowed ",ravariation)                        
						new_rjd_data = float(values[28]) - 2400000.0
						t_data = float(values[28]) - 2450000.0 
#						print('line values for sophie ',values)
						if rjd_data != new_rjd_data and values[7] !='':                          # de-dup
							rjd_data = new_rjd_data
							if self.created:
								self.rv_data['rjd'] += [rjd_data]
								self.rv_data['bjd'] += [float(values[28])]
								self.rv_data['BJD'] += [float(values[28])]
								self.rv_data['t'] += [t_data]                            
								self.rv_data['berv'] += [-99999]
								self.rv_data['berv_err'] += [-99999]
								self.rv_data['bispan'] += [-99999]
								self.rv_data['bispan_err'] += [-99999] 
								self.rv_data['drift_noise'] += [-99999]
								self.rv_data['drift_noise_err'] += [-99999]
								self.rv_data['texp'] += [-99999]
								self.rv_data['texp_err'] += [-99999]
								self.rv_data['cal_therror'] += [-99999]  
								self.rv_data['cal_therror_err'] += [-99999]
								self.rv_data['protm08'] += [-99999] 
								self.rv_data['protm08_err'] += [-99999] 
								self.rv_data['protn84'] += [-99999] 
								self.rv_data['protn84_err'] += [-99999] 
								self.rv_data['fwhm'] += [float(values[8])]  
								self.rv_data['fwhm_err'] += [-99999]
								self.rv_data['spectroFluxSn20'] += [-99999]
								self.rv_data['spectroFluxSn20_err'] += [-99999]
								self.rv_data['cal_thfile'] += [-99999] 
								self.rv_data['rhk_err'] += [-99999]
								self.rv_data['rhk'] += [self.rhk_value]
								self.rv_data['drs_qc'] += [True]
								self.rv_data['ins_name'] += ['SOPHIE']
								self.rv_data['Telescope'] += ['SOPHIE']                            
								self.rv_data['ins_mode'] += ['']
								self.rv_data['mask'] += [values[4]]
								self.rv_data['contrast'] += [-99999]
								self.rv_data['contrast_err'] += [-99999]
								self.rv_data['spectroFluxSn50'] += [-99999]
								self.rv_data['spectroFluxSn50_err'] += [-99999]
								self.rv_data['naindex'] += [-99999]
								self.rv_data['naindex_err'] += [-99999]
								self.rv_data['snca2'] += [-99999]
								self.rv_data['snca2_err'] += [-99999]
								self.rv_data['public'] += [True]
								self.rv_data['pub_reference'] += ['SOPHIE']
								self.rv_data['pub_bibcode'] += ['']
								self.rv_data['sindex'] += [-99999]
								self.rv_data['sindex_err'] += [-99999]
								self.rv_data['haindex'] += [-99999]
								self.rv_data['haindex_err'] += [-99999]
								self.rv_data['drs_version'] += ['Pub']  
								self.rv_data['caindex'] += [-99999]
								self.rv_data['caindex_err'] += [-99999]
								self.rv_data['rjd_err'] += [-99999]
								rv_value = (float(values[6])-self.gamma)*1000.0
								self.rv_data['rv'] += [rv_value]
								self.rv_data['Vel(m/s)'] += [rv_value]
								self.rv_data['vel'] += [rv_value]
#								print('sophie rv ',self.rv_data['rv'])
								if values[7] =='':
									self.rv_data['rv_err'] += [-99999]
									self.rv_data['ErrVel(m/s)'] += [-99999]
									self.rv_data['errvel'] += [-99999]                                
								else:
									rv_err_value = float(values[7])*1000.0
									self.rv_data['rv_err'] += [rv_err_value]
									self.rv_data['ErrVel(m/s)'] += [rv_err_value]
									self.rv_data['errvel'] += [rv_err_value]
#								print('sophie rv err ',self.rv_data['rv_err'])    
								self.rv_data['sn26'] += [float(values[23])]
								self.rv_data['ccf_noise'] += [-99999]
								self.rv_data['ccf_noise_err'] += [-99999]
								self.rv_data['drift_used'] += [-99999]
								self.rv_data['drift_used_err'] += [-99999]
								self.rv_data['ccf_asym'] += [-99999]
								self.rv_data['ccf_asym_err'] += [-99999]
								self.rv_data['date_night'] += ['']  
								self.rv_data['raw_file'] += ['SOPHIE']
								self.rv_data['prog_id'] += [' ']
								self.rv_data['th_ar'] += [-99999]
								self.rv_data['th_ar1'] += [-99999]
								self.rv_data['th_ar2'] += [-99999]                                
							else:
								self.rv_data['rjd'] = [rjd_data]
								self.rv_data['bjd'] = [float(values[28])]
								self.rv_data['BJD'] = [float(values[28])]
								self.rv_data['t'] = [t_data]                            
								self.rv_data['berv'] = [-99999]
								self.rv_data['berv_err'] = [-99999]
								self.rv_data['bispan'] = [-99999]
								self.rv_data['bispan_err'] = [-99999] 
								self.rv_data['drift_noise'] = [-99999]
								self.rv_data['drift_noise_err'] = [-99999]
								self.rv_data['texp'] = [-99999]
								self.rv_data['texp_err'] = [-99999]
								self.rv_data['cal_therror'] = [-99999]  
								self.rv_data['cal_therror_err'] = [-99999]
								self.rv_data['protm08'] = [-99999] 
								self.rv_data['protm08_err'] = [-99999] 
								self.rv_data['protn84'] = [-99999] 
								self.rv_data['protn84_err'] = [-99999] 
								self.rv_data['fwhm'] = [float(values[8])]  
								self.rv_data['fwhm_err'] = [-99999]
								self.rv_data['spectroFluxSn20'] = [-99999]
								self.rv_data['spectroFluxSn20_err'] = [-99999]
								self.rv_data['cal_thfile'] = [-99999] 
								self.rv_data['rhk_err'] = [-99999]
								self.rv_data['rhk'] = [self.rhk_value]
								self.rv_data['drs_qc'] = [True]
								self.rv_data['ins_name'] = ['SOPHIE']
								self.rv_data['Telescope'] = ['SOPHIE']                            
								self.rv_data['ins_mode'] = ['']
								self.rv_data['mask'] = [values[4]]
								self.rv_data['contrast'] = [-99999]
								self.rv_data['contrast_err'] = [-99999]
								self.rv_data['spectroFluxSn50'] = [-99999]
								self.rv_data['spectroFluxSn50_err'] = [-99999]
								self.rv_data['naindex'] = [-99999]
								self.rv_data['naindex_err'] = [-99999]
								self.rv_data['snca2'] = [-99999]
								self.rv_data['snca2_err'] = [-99999]
								self.rv_data['public'] = [True]
								self.rv_data['pub_reference'] = ['SOPHIE']
								self.rv_data['pub_bibcode'] = ['']
								self.rv_data['sindex'] = [-99999]
								self.rv_data['sindex_err'] = [-99999]
								self.rv_data['haindex'] = [-99999]
								self.rv_data['haindex_err'] = [-99999]
								self.rv_data['drs_version'] = ['Pub']  
								self.rv_data['caindex'] = [-99999]
								self.rv_data['caindex_err'] = [-99999]
								self.rv_data['rjd_err'] = [-99999]
								rv_value = (float(values[6])-self.gamma)*1000.0
								self.rv_data['rv'] = [rv_value]
								self.rv_data['Vel(m/s)'] = [rv_value]
								self.rv_data['vel'] = [rv_value]                            
								if values[7] =='':
									self.rv_data['rv_err'] = [-99999]
									self.rv_data['ErrVel(m/s)'] = [-99999]
									self.rv_data['errvel'] = [-99999]                                
								else:
									rv_err_value = float(values[7])
									self.rv_data['rv_err'] = [rv_err_value]*1000.0
									self.rv_data['ErrVel(m/s)'] = [rv_err_value]
									self.rv_data['errvel'] = [rv_err_value]
								self.rv_data['sn26'] = [float(values[23])]                            
								self.rv_data['ccf_noise'] = [-99999]
								self.rv_data['ccf_noise_err'] = [-99999]
								self.rv_data['drift_used'] = [-99999]
								self.rv_data['drift_used_err'] = [-99999]
								self.rv_data['ccf_asym'] = [-99999]
								self.rv_data['ccf_asym_err'] = [-99999]
								self.rv_data['date_night'] = ['']  
								self.rv_data['raw_file'] = ['SOPHIE']
								self.rv_data['prog_id'] = [' ']
								self.rv_data['th_ar'] = [-99999]
								self.rv_data['th_ar1'] = [-99999]
								self.rv_data['th_ar2'] = [-99999]                                
								self.created=True
							self.rv_sum_sophie = self.rv_sum_sophie + rv_value
							self.rv_err_sum_sophie = self.rv_err_sum_sophie + rv_err_value
							self.rv_count_sophie +=1
							string = values[28].strip()
#							obsstring = "SOPHIE" + " & " + string + " & " + str(rv_value) + " & " + str(rv_err_value) + " \\\\"   
#							with open(self.obsname, 'a') as observations_log:
#								observations_log.write('\n')                
#								observations_log.write(obsstring)    
		myfile4.close()
                    
	def _get_uves_2019(self):
		print("enter get_uves_2019")
		print()        
		for filename in os.listdir('uves_vels_2019'):
			if fnmatch.fnmatch(filename, '*.vels'):
				target = filename.split('_',1)
				object_name = target[0]
#				print("object name ", object_name," filename ",filename)                
				i=0
				try:
					nospace_count = len(self.target_alias_nospace)
				except:
					nospace_count = 0
				while i < nospace_count:
					if object_name == self.target_alias_nospace[i][0]:
						fileopen = 'uves_vels_2019/'+filename
						myfile6 = open (fileopen, 'r')
						for line in myfile6:
							values = line.split()                    
							rjd_data = float(values[0]) - 2400000.0
							t_data = float(values[0]) - 2450000.0                    
							if self.created:
								self.rv_data['rjd'] += [rjd_data]
								self.rv_data['bjd'] += [float(values[0])]
								self.rv_data['BJD'] += [float(values[0])]    
								self.rv_data['t'] += [t_data]                                                
								self.rv_data['berv'] += [-99999]                       
								self.rv_data['berv_err'] += [-99999]
								self.rv_data['bispan'] += [-99999]                        
								self.rv_data['bispan_err'] += [-99999]                        
								self.rv_data['drift_noise'] += [-99999]
								self.rv_data['drift_noise_err'] += [-99999]
								self.rv_data['texp'] += [-99999]
								self.rv_data['texp_err'] += [-99999]
								self.rv_data['cal_therror'] += [-99999] 
								self.rv_data['cal_therror_err'] += [-99999]
								self.rv_data['protm08'] += [-99999] 
								self.rv_data['protm08_err'] += [-99999] 
								self.rv_data['protn84'] += [-99999] 
								self.rv_data['protn84_err'] += [-99999] 
								self.rv_data['fwhm'] += [-99999]                        
								self.rv_data['fwhm_err'] += [-99999]                        
								self.rv_data['spectroFluxSn20'] += [-99999]
								self.rv_data['spectroFluxSn20_err'] += [-99999]
								self.rv_data['cal_thfile'] += [-99999] 
								self.rv_data['rhk_err'] += [-99999]
								self.rv_data['rhk'] += [self.rhk_value]
								self.rv_data['drs_qc'] += [True]
								self.rv_data['ins_name'] += ['UVES2019']
								self.rv_data['Telescope'] += ['UVES2019']                        
								self.rv_data['ins_mode'] += ['']
								self.rv_data['mask'] += ['']
								self.rv_data['contrast'] += [-99999]
								self.rv_data['contrast_err'] += [-99999]                        
								self.rv_data['spectroFluxSn50'] += [-99999]
								self.rv_data['spectroFluxSn50_err'] += [-99999]
								self.rv_data['naindex'] += [-99999]
								self.rv_data['naindex_err'] += [-99999]
								self.rv_data['snca2'] += [-99999]
								self.rv_data['snca2_err'] += [-99999]
								self.rv_data['public'] += [True]
								self.rv_data['pub_reference'] += ['none']
								self.rv_data['pub_bibcode'] += ['']
								self.rv_data['sindex'] += [-99999]
								self.rv_data['sindex_err'] += [-99999]
								self.rv_data['haindex'] += [-99999]
								self.rv_data['haindex_err'] += [-99999]
								self.rv_data['drs_version'] += ['Pub']  
								self.rv_data['caindex'] += [-99999]
								self.rv_data['caindex_err'] += [-99999]
								self.rv_data['rjd_err'] += [-99999]
								rv_value = float(values[1])*1000.0
								self.rv_data['rv'] += [rv_value]
								self.rv_data['Vel(m/s)'] += [rv_value]
								self.rv_data['vel'] += [rv_value]                         
								rv_err_value = float(values[2])*1000.0
								self.rv_data['rv_err'] += [rv_err_value]
								self.rv_data['ErrVel(m/s)'] += [rv_err_value]
								self.rv_data['errvel'] += [rv_err_value]
								self.rv_data['sn26'] += [-99999]                        
								self.rv_data['ccf_noise'] += [-99999]
								self.rv_data['ccf_noise_err'] += [-99999]
								self.rv_data['drift_used'] += [-99999]
								self.rv_data['drift_used_err'] += [-99999]
								self.rv_data['ccf_asym'] += [-99999]
								self.rv_data['ccf_asym_err'] += [-99999]
								self.rv_data['date_night'] += ['']  
								self.rv_data['raw_file'] += ['none']
								self.rv_data['prog_id'] += [' ']
								self.rv_data['th_ar'] += [-99999]
								self.rv_data['th_ar1'] += [-99999]
								self.rv_data['th_ar2'] += [-99999]                                
							else:
								self.rv_data['rjd'] = [rjd_data]
								self.rv_data['bjd'] = [float(values[0])]
								self.rv_data['BJD'] = [float(values[0])]
								self.rv_data['t'] = [t_data]                                                
								self.rv_data['berv'] = [-99999]
								self.rv_data['berv_err'] = [-99999]
								self.rv_data['berv_err'] = [-99999]
								self.rv_data['bispan'] = [-99999]                        
								self.rv_data['bispan_err'] = [-99999] 
								self.rv_data['drift_noise'] = [-99999]
								self.rv_data['drift_noise_err'] = [-99999]
								self.rv_data['texp'] = [-99999]
								self.rv_data['texp_err'] = [-99999]
								self.rv_data['cal_therror'] = [-99999]  
								self.rv_data['cal_therror_err'] = [-99999]
								self.rv_data['protm08'] = [-99999] 
								self.rv_data['protm08_err'] = [-99999] 
								self.rv_data['protn84'] = [-99999] 
								self.rv_data['protn84_err'] = [-99999] 
								self.rv_data['fwhm'] = [-99999]                        
								self.rv_data['fwhm_err'] = [-99999]
								self.rv_data['spectroFluxSn20'] = [-99999]
								self.rv_data['spectroFluxSn20_err'] = [-99999]
								self.rv_data['cal_thfile'] = [-99999] 
								self.rv_data['rhk_err'] = [-99999]
								self.rv_data['rhk'] = [self.rhk_value]
								self.rv_data['drs_qc'] = [True]
								self.rv_data['ins_name'] = ['UVES2019']
								self.rv_data['Telescope'] = ['UVES2019']                        
								self.rv_data['ins_mode'] = ['']
								self.rv_data['mask'] = ['']
								self.rv_data['contrast'] = [-99999]
								self.rv_data['contrast_err'] = [-99999]
								self.rv_data['spectroFluxSn50'] = [-99999]
								self.rv_data['spectroFluxSn50_err'] = [-99999]
								self.rv_data['naindex'] = [-99999]
								self.rv_data['naindex_err'] = [-99999]
								self.rv_data['snca2'] = [-99999]
								self.rv_data['snca2_err'] = [-99999]
								self.rv_data['public'] = [True]
								self.rv_data['pub_reference'] = ['none']
								self.rv_data['pub_bibcode'] = ['']
								self.rv_data['sindex'] = [-99999]
								self.rv_data['sindex_err'] = [-99999]
								self.rv_data['haindex'] = [-99999]
								self.rv_data['haindex_err'] = [-99999]
								self.rv_data['drs_version'] = ['Pub']  
								self.rv_data['caindex'] = [-99999]
								self.rv_data['caindex_err'] = [-99999]
								self.rv_data['rjd_err'] = [-99999]
								rv_value = float(values[1])*1000.0
								self.rv_data['rv'] = [rv_value]
								self.rv_data['Vel(m/s)'] = [rv_value]
								self.rv_data['vel'] = [rv_value]                        
								rv_err_value = float(values[2])*1000.0
								self.rv_data['rv_err'] = [rv_err_value]
								self.rv_data['ErrVel(m/s)'] = [rv_err_value]
								self.rv_data['errvel'] = [rv_err_value]
								self.rv_data['sn26'] = [-99999]                        
								self.rv_data['ccf_noise'] = [-99999]
								self.rv_data['ccf_noise_err'] = [-99999]
								self.rv_data['drift_used'] = [-99999]
								self.rv_data['drift_used_err'] = [-99999]
								self.rv_data['ccf_asym'] = [-99999]
								self.rv_data['ccf_asym_err'] = [-99999]
								self.rv_data['date_night'] = ['']  
								self.rv_data['raw_file'] = ['none']
								self.rv_data['prog_id'] = [' ']
								self.rv_data['th_ar'] = [-99999]
								self.rv_data['th_ar1'] = [-99999]
								self.rv_data['th_ar2'] = [-99999]                                
								self.created=True
							self.rv_sum_uves_2019 = self.rv_sum_uves_2019 + rv_value
							self.rv_err_sum_uves_2019 = self.rv_err_sum_uves_2019 + rv_err_value
							self.rv_count_uves_2019 +=1
#							obsstring = "UVES 2019" + " & " + values[0] + " & " + values[1] + " & " + values[2] + " \\\\"   
#							with open(self.obsname, 'a') as observations_log:
#								observations_log.write('\n')                        
#								observations_log.write(obsstring)
						myfile6.close() 
					i+=1                     
        

# hardener: add RHK and DRS-QC since they were not included
# hardener: set up "try - except block" for the file read to cover stars not represented in HIRES
	def _shrink_data(self,i):
		a = np.array(self.rv_data['berv'])
		a = np.delete(a, i)
		self.rv_data['berv']=a
		a = np.array(self.rv_data['berv_err'])
		a = np.delete(a, i)
		self.rv_data['berv_err']=a
		a = np.array(self.rv_data['bispan'])
		a = np.delete(a, i)
		self.rv_data['bispan']=a
		a = np.array(self.rv_data['bispan_err'])
		a = np.delete(a, i)
		self.rv_data['bispan_err']=a
		a = np.array(self.rv_data['drift_noise'])
		a = np.delete(a, i)
		self.rv_data['drift_noise']=a
		a = np.array(self.rv_data['drift_noise_err'])
		a = np.delete(a, i)
		self.rv_data['drift_noise_err']=a
		a =np.array(self.rv_data['texp'])
		a = np.delete(a, i)
		self.rv_data['texp']=a
		a = np.array(self.rv_data['texp_err'])
		a = np.delete(a, i)
		self.rv_data['texp_err']=a
		a = np.array(self.rv_data['cal_therror'])
		a = np.delete(a, i)
		self.rv_data['cal_therror']=a
		a = np.array(self.rv_data['cal_therror_err'])
		a = np.delete(a, i)
		self.rv_data['cal_therror_err']=a
		a = np.array(self.rv_data['fwhm'])
		a = np.delete(a, i)
		self.rv_data['fwhm']=a
		a = np.array(self.rv_data['fwhm_err'])
		a = np.delete(a, i)
		self.rv_data['fwhm_err']=a
		a = np.array(self.rv_data['protm08'])
		a = np.delete(a, i)
		self.rv_data['protm08']=a
		a = np.array(self.rv_data['protm08_err'])
		a = np.delete(a, i)
		self.rv_data['protm08_err']=a
		a = np.array(self.rv_data['protn84'])
		a = np.delete(a, i)
		self.rv_data['protn84']=a
		a = np.array(self.rv_data['protn84_err'])
		a = np.delete(a, i)
		self.rv_data['protn84_err']=a
		a = np.array(self.rv_data['public'])
		a = np.delete(a, i)
		self.rv_data['public']=a
		a = np.array(self.rv_data['spectroFluxSn20'])
		a = np.delete(a, i)
		self.rv_data['spectroFluxSn20']=a
		a = np.array(self.rv_data['spectroFluxSn20_err'])
		a = np.delete(a, i)
		self.rv_data['spectroFluxSn20_err']=a
		a = np.array(self.rv_data['ins_mode'])
		a = np.delete(a, i)
		self.rv_data['ins_mode']=a
		a = np.array(self.rv_data['rhk'])
		a = np.delete(a, i)
		self.rv_data['rhk']=a
		a = np.array(self.rv_data['rhk_err'])
		a = np.delete(a, i)
		self.rv_data['rhk_err']=a
		a = np.array(self.rv_data['ins_name'])
		a = np.delete(a, i)
		self.rv_data['ins_name']=a
		a = np.array(self.rv_data['Telescope'])
		a = np.delete(a, i)
		self.rv_data['Telescope']=a        
		a = np.array(self.rv_data['mask'])
		a = np.delete(a, i)
		self.rv_data['mask']=a
		a = np.array(self.rv_data['contrast'])
		a = np.delete(a, i)
		self.rv_data['contrast']=a
		a = np.array(self.rv_data['contrast_err'])
		a = np.delete(a, i)
		self.rv_data['contrast_err']=a
		a = np.array(self.rv_data['cal_thfile'])
		a = np.delete(a, i)
		self.rv_data['cal_thfile']=a
		a = np.array(self.rv_data['spectroFluxSn50'])
		a = np.delete(a, i)
		self.rv_data['spectroFluxSn50']=a
		a = np.array(self.rv_data['spectroFluxSn50_err'])
		a = np.delete(a, i)
		self.rv_data['spectroFluxSn50_err']=a
		a = np.array(self.rv_data['naindex'])
		a = np.delete(a, i)
		self.rv_data['naindex']=a
		a = np.array(self.rv_data['naindex_err'])
		a = np.delete(a, i)
		self.rv_data['naindex_err']=a
		a = np.array(self.rv_data['snca2'])
		a = np.delete(a, i)
		self.rv_data['snca2']=a
		a = np.array(self.rv_data['snca2_err'])
		a= np.delete(a, i)
		self.rv_data['snca2_err']=a
		a = np.array(self.rv_data['pub_reference'])
		a = np.delete(a, i)
		self.rv_data['pub_reference']=a
		a = np.array(self.rv_data['pub_bibcode'])
		a = np.delete(a, i)
		self.rv_data['pub_bibcode']=a
		a = np.array(self.rv_data['sindex'])
		a = np.delete(a, i)
		self.rv_data['sindex']=a
		a = np.array(self.rv_data['sindex_err'])
		a = np.delete(a, i)
		self.rv_data['sindex_err']=a
		a = np.array(self.rv_data['drs_qc'])
		a = np.delete(a, i)
		self.rv_data['drs_qc']=a
		a = np.array(self.rv_data['haindex'])
		a = np.delete(a, i)
		self.rv_data['haindex']=a
		a = np.array(self.rv_data['haindex_err'])
		a = np.delete(a, i)
		self.rv_data['haindex_err']=a
		a = np.array(self.rv_data['drs_version'])
		a = np.delete(a, i)
		self.rv_data['drs_version']=a
		a = np.array(self.rv_data['caindex'])
		a = np.delete(a, i)
		self.rv_data['caindex']=a
		a = np.array(self.rv_data['caindex_err'])
		a = np.delete(a, i)
		self.rv_data['caindex_err']=a
		a = np.array(self.rv_data['rjd'])
		a = np.delete(a, i)
		self.rv_data['rjd']=a
		a = np.array(self.rv_data['bjd'])
		a = np.delete(a, i)
		self.rv_data['bjd']=a
		a = np.array(self.rv_data['BJD'])
		a = np.delete(a, i)
		self.rv_data['BJD']=a
		a = np.array(self.rv_data['t'])
		a = np.delete(a, i)
		self.rv_data['t']=a
		a = np.array(self.rv_data['rjd_err'])
		a = np.delete(a, i)
		self.rv_data['rjd_err']=a
		a = np.array(self.rv_data['rv'])
		a = np.delete(a, i)
		self.rv_data['rv']=a 
		a = np.array(self.rv_data['Vel(m/s)'])
		a = np.delete(a, i)
		self.rv_data['Vel(m/s)']=a
		a = np.array(self.rv_data['vel'])
		a = np.delete(a, i)
		self.rv_data['vel']=a        
		a = np.array(self.rv_data['rv_err'])
		a = np.delete(a, i)
		self.rv_data['rv_err']=a
		a = np.array(self.rv_data['ErrVel(m/s)'])
		a = np.delete(a, i)
		self.rv_data['ErrVel(m/s)']=a
		a = np.array(self.rv_data['errvel'])
		a = np.delete(a, i)
		self.rv_data['errvel']=a 
		a = np.array(self.rv_data['sn26'])
		a = np.delete(a, i)
		self.rv_data['sn26']=a        
		a = np.array(self.rv_data['ccf_noise'])
		a = np.delete(a, i)
		self.rv_data['ccf_noise']=a
		a = np.array(self.rv_data['ccf_noise_err'])
		a= np.delete(a, i)
		self.rv_data['ccf_noise_err']=a
		a = np.array(self.rv_data['drift_used'])
		a = np.delete(a, i)
		self.rv_data['drift_used']=a
		a = np.array(self.rv_data['drift_used_err'])
		a = np.delete(a, i)
		self.rv_data['drift_used_err']=a
		a = np.array(self.rv_data['ccf_asym'])
		a = np.delete(a, i)
		self.rv_data['ccf_asym']=a
		a = np.array(self.rv_data['ccf_asym_err'])
		a = np.delete(a, i)
		self.rv_data['ccf_asym_err']=a
		a = np.array(self.rv_data['date_night'])
		a = np.delete(a, i)
		self.rv_data['date_night']=a
		a = np.array(self.rv_data['raw_file'])
		a = np.delete(a, i)
		self.rv_data['raw_file']=a
		a = np.array(self.rv_data['prog_id'])
		a = np.delete(a, i)
		self.rv_data['prog_id']=a
		a = np.array(self.rv_data['th_ar'])
		a = np.delete(a, i)
		self.rv_data['th_ar']=a
		a = np.array(self.rv_data['th_ar1'])
		a = np.delete(a, i)
		self.rv_data['th_ar1']=a
		a = np.array(self.rv_data['th_ar2'])
		a = np.delete(a, i)
		self.rv_data['th_ar2']=a        
	def _trim_data(self):     
		k=len(self.rv_data['rjd'])       
		i=0
#		print()
#		with open(self.logname, 'a') as norm_log:
#			norm_log.write('\n')        
		while i < k:
#			print('********* rjd ',self.rv_data['rjd'][i],'**************')             
			if (self.rv_data['rjd'][i] < self.trim_lower or self.rv_data['rjd'][i] > self.trim_upper):
#				print('between upper and lower trim ',self.rv_data['rjd'][i])
#				line = 'between upper and lower trim ' + self.rv_data['rjd'][i]
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)            
				self._shrink_data(i)
				i=0
				k=len(self.rv_data['rjd'])
			elif (self.rv_data['rjd'][i] > self.remove_date_ranges[0][0] and self.rv_data['rjd'][i] < self.remove_date_ranges[0][1]):
#				print('between upper and lower range1 ',self.rv_data['rjd'][i])
#				line = 'between upper and lower range1 ' + str(self.rv_data['rjd'][i])
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)                
				self._shrink_data(i)
				i=0
				k=len(self.rv_data['rjd'])                
			elif (self.rv_data['rjd'][i]> self.remove_date_ranges[1][0] and self.rv_data['rjd'][i] < self.remove_date_ranges[1][1]):
#				print('between upper and lower range2 ',self.rv_data['rjd'][i])
#				line = 'between upper and lower range2 ' + str(self.rv_data['rjd'][i])
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)                
				self._shrink_data(i)
				i=0
				k=len(self.rv_data['rjd'])
			elif (self.rv_data['rjd'][i]> self.remove_date_ranges[2][0] and self.rv_data['rjd'][i] < self.remove_date_ranges[2][1]):
#				print('between upper and lower range3 ',self.rv_data['rjd'][i])
#				line = 'between upper and lower range3 ' + str(self.rv_data['rjd'][i])
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)                
				self._shrink_data(i)
				i=0
				k=len(self.rv_data['rjd'])  
			elif (self.rv_data['rjd'][i]> self.remove_date_ranges[3][0] and self.rv_data['rjd'][i] < self.remove_date_ranges[3][1]):
#				print('between upper and lower range4 ',self.rv_data['rjd'][i])
#				line = 'between upper and lower range4 ' + str(self.rv_data['rjd'][i])
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)                
				self._shrink_data(i)
				i=0
				k=len(self.rv_data['rjd']) 
			elif (self.rv_data['rjd'][i]> self.remove_date_ranges[4][0] and self.rv_data['rjd'][i] < self.remove_date_ranges[4][1]):
#				print('between upper and lower range5 ',self.rv_data['rjd'][i])
#				line = 'between upper and lower range5 ' + str(self.rv_data['rjd'][i])
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)                
				self._shrink_data(i)
				i=0
				k=len(self.rv_data['rjd'])                
			elif (self.rv_data['rjd'][i] in self.remove_date_singles):
#				print('in remove dates ',self.rv_data['rjd'][i])
#				line = 'in remove dates ' + str(self.rv_data['rjd'][i])
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)                
				self._shrink_data(i)                
				i=0
				k=len(self.rv_data['rjd'])
			elif (self.rv_data['sn26'][i] < self.clip_snr) and self.rv_data['sn26'][i] > -99999:
#				print(' sn26 in remove dates ',self.rv_data['rjd'][i], 'sn26 ', self.rv_data['sn26'][i], 'snr_clip ', self.clip_snr)
#				line = ' sn26 in remove dates ' + str(self.rv_data['rjd'][i]) +  'sn26 ' +  str(self.rv_data['sn26'][i]) +  'snr_clip ' +  str(self.clip_snr)
#				with open(self.logname, 'a') as norm_log:
#					norm_log.write('\n')
#					norm_log.write(line)                
				self._shrink_data(i)                
				i=0
				k=len(self.rv_data['rjd'])                
			else:
				i+=1
#		print()
#		with open(self.logname, 'a') as norm_log:
#			norm_log.write('\n')
#hardener: can use this method to delete outlier points (remove_dates)
 