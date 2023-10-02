import numpy as np
import netCDF4 as nc4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os

from Functions import *
from Get_GLM import *
from Get_ABI import *
from ABI_grid import *

YYYY = '2023'
Sat = '16'
HOURS = np.arange(0,24,1)
i = 238

LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

lat_min = CONFIG.LAT_MIN - CONFIG.LAT_RES/2.
lat_max = CONFIG.LAT_MAX + CONFIG.LAT_RES/2.
lon_min = CONFIG.LON_MIN - CONFIG.LON_RES/2.
lon_max = CONFIG.LON_MIN + CONFIG.LON_RES/2.


#Convert Day of Year to Month Day
MMDD = (doy_to_mmdd(YYYY,i))
print('Starting for %s'%MMDD)
MM = MMDD[:2]
DD = MMDD[2:]

# Add the satellite zenith angle
lat_sat, lon_sat, area_sat = get_abi_grid()
satzen,satazm = read_satzen()

# Path of GLM
GLM_LNG_PATH = f'/Users/nmesa/Desktop/NASA_GLM_Project/gitrepo/nasaglm/Case_Studies/AL102023/{i}/' 

for HI, HOUR in enumerate(HOURS):
    #Get all the files over this hour (not always the same)
    LIGHTNING_LIST = GLM_LNG_PATH+str(HOUR).zfill(2)+'/'
    print('Processing %i GLM files for %iZ'%(len(LIGHTNING_LIST),HOUR))
    #REMEMBER QC SET TO FALSE
    lats_l, lons_l, ener_l, area_l, flag_l, flid_l, lats_lf, lons_lf, ener_lf, area_lf, flag_lf = read_glm_filelist(LIGHTNING_LIST, apply_qc=False,output_qc=False)
    quick_display(ener_lf)
                

print('test')

#2023/0826/OR_GLM-L2-LCFA_G16_s2023
