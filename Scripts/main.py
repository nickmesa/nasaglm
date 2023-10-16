import numpy as np
import netCDF4 as nc4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os

from Functions import *
from Get_GLM import *
from Get_ABI import *
from ABI_grid import *
from mesa_functions import *

YYYY = '2021'
Sat = '16'
HOURS = np.arange(0,24,1)
i = 238
#ATCF_storm_number = 'AL092021'
time_interval_min = 5

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

# Path of GLM
GLM_LNG_PATH = '/Users/nasarosesnas/glmtc/GLM/L2/'  + YYYY + '/' + MMDD +'/'
#GLM_LNG_PATH = f'/Users/nmesa/Desktop/NASA_GLM_Project/gitrepo/nasaglm/Case_Studies/AL102023/{i}/' 

# Add the satellite zenith angle
lat_sat, lon_sat, area_sat = get_abi_grid()
satzen,satazm = read_satzen()

for HI, HOUR in enumerate(HOURS):
  #Get all the files over this hour (not always the same)
    #LIGHTNING_LIST = GLM_LNG_PATH+str(HOUR).zfill(2) + '/'
    print('Processing %i GLM files for %iZ'%(len(GLM_LNG_PATH),HOUR))
    min_interval_filelist = time_interval_filelist(time_interval_min, GLM_LNG_PATH) #LIGHTNING_LIST) #)
    for min in range(0,len(min_interval_filelist),1):
        print('Processing for %i to %i minute interval in hour %iZ'%(min,(min*time_interval_min) + time_interval_min,HOUR))
        #REMEMBER QC SET TO FALSE
        lats_l, lons_l, ener_l, area_l, flag_l, flid_l, lats_lf, lons_lf, ener_lf, area_lf, flag_lf = read_glm_filelist(min_interval_filelist[min], GLM_LNG_PATH, apply_qc=False,output_qc=False) #

        #HURDAT2 best track center position
        center_lat, center_lon, vmax, mslp = get_hurdat('ida',2021,8,26,12,00)

        #SHIPS Shear Data
        

        ax = plt.axes(projection = ccrs.PlateCarree())
        #ax.add_feature(cfeature.COASTLINE)
        ax.set_extent([-60, -110, 0, 40])
        ax.scatter(lons_lf, lats_lf)
        plt.plot(-86.10,20.80, color = 'red', marker = 'x')
        plt.savefig(f'/Users/nmesa/Desktop/NASA_GLM_Project/gitrepo/nasaglm/Output_Plots/FlashLoc{YYYY}{MMDD}{HOUR}{min*time_interval_min}to{(min*time_interval_min) + time_interval_min}')
    break







