import csv
import pandas as pd
class rv_file_extract:
	def __init__(self,target):       
		self.target=target              
	def get_data(self,target):
		self.target=target
		self.rv_data={}
		self.created=False 
		header_lines = 1
		count_lines = 0        
		myfile = open('/Users/subercorley/EXOTIC/examples/radial_velocity/data/'+self.target+'.csv', 'r')
		for line in myfile:
			if count_lines < header_lines:
				count_lines+=1
				continue
			values = line.split(',')
			rvtime = values[1]            
			vel = values[2]
			err = values[3]
			radegrees = values[5]
			decdegrees = values[6]
			rjdtime = values[7]            
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
   