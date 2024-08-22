import csv
import pandas as pd
class rv_file_extract:
	def __init__(self,target,exclude_instruments,exclude_times,harps_fiber_cutoff_time):       
		self.target=target
		self.exclude_instruments = exclude_instruments
		self.exclude_times = exclude_times 
		self.harps_fiber_cutoff_time = harps_fiber_cutoff_time        
	def get_data(self,target):
		self.target=target.replace(" ","")
		self.rv_data={}
		self.created=False 
		header_lines = 1
		count_lines = 0        
		myfile = open('radial_velocity/data/'+self.target+'.csv', 'r')
		for line in myfile:
			if count_lines < header_lines:
				count_lines+=1
				continue
			values = line.split(',')
			rvtime = values[1]            
			vel = values[2]
			err = values[3]
			instrument = values[4]          
			radegrees = values[5]
			decdegrees = values[6]
			rjdtime = float(values[7])            
			if rjdtime not in self.exclude_times and instrument not in self.exclude_instruments:                 
				if instrument== "HARPScat":
					if rjdtime > self.harps_fiber_cutoff_time:                    
						if self.created:
							self.rv_data['target']+=[values[0]]
							self.rv_data['BJD']+=[float(rvtime)]
							self.rv_data['Vel(m/s)']+=[float(vel)]
							self.rv_data['ErrVel(m/s)']+=[float(err)]
							self.rv_data['Telescope']+=[values[4]]                
							self.rv_data['RA']+=[float(radegrees)]
							self.rv_data['DEC']+=[float(decdegrees)]
							self.rv_data['rjd']+=[float(rjdtime)]
							self.rv_data['rv']+=[float(vel)]
							self.rv_data['ins_name']+=[values[4]]                
						else:
							self.rv_data['target']=[values[0]]
							self.rv_data['BJD']=[float(rvtime)]
							self.rv_data['Vel(m/s)']=[float(vel)]
							self.rv_data['ErrVel(m/s)']=[float(err)]
							self.rv_data['Telescope']=[values[4]]                
							self.rv_data['RA']=[float(radegrees)]
							self.rv_data['DEC']=[float(decdegrees)]                
							self.rv_data['rjd']=[float(rjdtime)]
							self.rv_data['rv']=[float(vel)]
							self.rv_data['ins_name']=[values[4]]                
							self.created=True
				else:
					if self.created:
						self.rv_data['target']+=[values[0]]
						self.rv_data['BJD']+=[float(rvtime)]
						self.rv_data['Vel(m/s)']+=[float(vel)]
						self.rv_data['ErrVel(m/s)']+=[float(err)]
						self.rv_data['Telescope']+=[values[4]]                
						self.rv_data['RA']+=[float(radegrees)]
						self.rv_data['DEC']+=[float(decdegrees)]
						self.rv_data['rjd']+=[float(rjdtime)]
						self.rv_data['rv']+=[float(vel)]
						self.rv_data['ins_name']+=[values[4]]                
					else:
						self.rv_data['target']=[values[0]]
						self.rv_data['BJD']=[float(rvtime)]
						self.rv_data['Vel(m/s)']=[float(vel)]
						self.rv_data['ErrVel(m/s)']=[float(err)]
						self.rv_data['Telescope']=[values[4]]                
						self.rv_data['RA']=[float(radegrees)]
						self.rv_data['DEC']=[float(decdegrees)]                
						self.rv_data['rjd']=[float(rjdtime)]
						self.rv_data['rv']=[float(vel)]
						self.rv_data['ins_name']=[values[4]]                
						self.created=True                    
		myfile.close()
#		print(self.rv_data)        
		rv_data_df=pd.DataFrame.from_dict(self.rv_data)
		return rv_data_df
   