def do_FRP_AC(filda_id, lut_id, geos_fp=False, geos_fp_id=None):

	'''
	Python function to conduct the first order atmospheric correction for FRP retrieval
	----------------
	Paramters:
	filda_id: handler of a VIIRS FILDA-2 file
	lut_id: handler of the look-up table file (downloadable from xxxxx)
	geos_fp: whether to correct for water vapor absorption using GEOS-FP
	geos_fp_id: handler of a GEOS-FP file that temporally matches the FILDA-2 file
	----------------
	Returns:
	FP_Power_AC: Atmosheric corrected FRP (unit: MW)
	'''
	
	import numpy as np
	import pandas as pd
	from netCDF4 import Dataset
	
	# read FILDA-2 file
	print('Reading FILDA-2 file: ', filda_id)
	f = Dataset(filda_id, 'r')
	FP_Latitude = f.variables['FP_Latitude'][:] # fire pixel latitude (degree)
	FP_Longitude = f.variables['FP_Longitude'][:] # fire pixel longitude (degree)
	FP_Power = f.variables['FP_Power'][:] # uncorrected FRP (MW)
	Sensor_Zenith = f.variables['Sensor_Zenith'][:] # sensor zenith angle (degree)
	f.close()
	
	# read look-up table used for atmospheric correction
	print('Reading look-up table file: ', lut_id)
	LUT = pd.read_csv(lut_id) # rows are sensor zenith angle and columns are precipitable water
	PW_array = np.array(LUT.keys()) # get a list of precipitable array (mm)
	
	# read GEOS-FP file 
	if geos_fp:
		if not geos_fp_id:
			raise ValueError('geos_fp_id is mandatory when geos_fp is True')
		else: 
			print('Reading GEOS-FP file: ', geos_fp_id)
			f = Dataset(geos_fp_id, 'r')
			GEOS_PW = f.variables['TQV'][:] # extract modeled precipitable water
			GEOS_Latitude = f.variables['lat'][:] # read GEOS-FP lat and lon to later match with FRP spatially
			GEOS_Longitude = f.variables['lon'][:]
			lat_interval = GEOS_Latitude[1] - GEOS_Latitude[0] # 0.25 deg interval in lat direction
			lon_interval = GEOS_Longitude[1] - GEOS_Longitude[0] # 0.3125 deg interval in the lon direction		
			lat_min = np.min(GEOS_Latitude)
			lon_min = np.min(GEOS_Longitude)
			f.close()
			
	# perform atmospheric correction
	FP_Power_AC = np.full(FP_Power.shape, np.nan)		
	for i in range(len(FP_Power)):
		if geos_fp: # if using precipitable water from GEOS-FP
			lat_idx = int((FP_Latitude[i]-(lat_min))//lat_interval)
			lon_idx = int((FP_Longitude[i]-(lon_min))//lon_interval)	
			PW = GEOS_PW[0,lat_idx,lon_idx] # kg/m2 (or mm)	
		else: # if not using GEOS-FP
			PW = 30 # assume a typical value of precipitable water (mm)
		
		PW_idx = np.argmin(np.abs(PW - PW_array.astype(float))) # find the closest PW entry in the look-up table
		VZA_idx = Sensor_Zenith[i].astype(int) # find the closest VZA entry in the Look-up table
		tau = LUT[PW_array[PW_idx]][VZA_idx] # read transmittance based on PW and VZA		
		FP_Power_AC[i] = FP_Power[i]/tau # correct FRP

# 		# checking:
# 		print('fire pixel index is ', i)
# 		print('FRP uncorrected is ', FP_Power[i]) 
# 		print('view zenith angle is ', Sensor_Zenith[i]) 		
# 		print('fire pixel latitude is ', FP_Latitude[i])
# 		print('fire pixel longitude is ', FP_Longitude[i])
# 		print('lat index for GEOS-FP is ', lat_idx)
# 		print('lon index for GEOS-FP is ', lon_idx)
# 		print('precipitable water is ', PW)
# 		print('PW index in LUT is ', PW_idx)
# 		print('VZA index in LUT is ', VZA_idx)
# 		print('transmitance is ', tau)
# 		print('FRP corrected is', FP_Power_AC[i])
# 		print('')

	return FP_Power_AC
	
filda_id = './VNP47IMG.A2014365.2354.002.2024304000833.nc'
lut_id = './LUT_VNP.csv'
geos_fp = True
geos_fp_id = './GEOS-FP/M2I1NXASM_Single-Level_Diagnostics/GEOS.fp.asm.inst3_2d_asm_Nx.20141231_2100.V01.nc4'
FP_Power_AC = do_FRP_AC(filda_id, lut_id, geos_fp, geos_fp_id)
