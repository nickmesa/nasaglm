#Code to regrid ABI brightness temperature data onto 0.02-degree regular grid

#Import necessary modules
import numpy as np
from scipy.signal import convolve
from netCDF4 import Dataset
import os
import glob
import shutil
import datetime
from scipy.interpolate import griddata as GridData
import warnings
warnings.filterwarnings("ignore")

from Functions import *
from ABI_grid import *


def read_abi(filename):
    ds=Dataset(filename)

    #Code from Alex from Kyle to convert Radiance to brightness temperature
    rad=ds.variables['Rad'][...]
    c1=ds.variables['planck_fk1'][0]
    c2=ds.variables['planck_fk2'][0]
    bc1=ds.variables['planck_bc1'][0]
    bc2=ds.variables['planck_bc2'][0]
    bt=c2/np.log((c1/rad)+1.)
    amap=(bt-bc1)/bc2

    DQF=ds.variables['DQF'][...]
    ds.close()
    return amap, DQF


def QC_abi(amap, DQF, area):
    #Input amap= brightness temperatures
    #Input DQF = Data Quality Flags
    #Input area = ABI area for each grid cell

    #First check the percentage of pixels to decide if they are good/bad
    pix_cold=np.where((amap.data>0) & (amap.data<183))
    pix_warm=np.where((amap.data<600) & (amap>313))
    pix_tot=np.where((amap.data>0) & (amap.data<600))

    if pix_tot==0: #no good pixels
        amap[:,:] = np.nan
        DQF[:,:]  = 9 #Bad fraction exceeds 10%
        return amap, DQF


    #Calculate the fraction of bad data to total data
    bad_frac=(len(amap.data[pix_cold])+len(amap.data[pix_warm]))/len(amap.data[pix_tot])

    if bad_frac>=.1: #if more than 10% of the data is bad, then we shouldn't use this file
        amap[:,:] = np.nan
        DQF[:,:]  = 9 #Bad fraction exceeds 10%
        return amap, DQF

    else: #If we have lots of good data 
        #Nan values where amap.mask==True
        abi_value=np.where(amap.mask==False,amap,np.nan)

        #Nan values where DQF flag value is equal to 1, 2, 3, or 4 (These indicate suspect data)
        abi_value=np.where(DQF>=1,np.nan,abi_value)

        #Nan values where pixel area is greater than 40 km
        abi_value=np.where(area > 40,np.nan,abi_value)
        DQF      =np.where(area > 40, 7 , DQF)

        #Nan values where the brightness temperature is less than 150 K or greater than 350 K
        abi_value=np.where(abi_value<150,np.nan,abi_value)
        abi_value=np.where(abi_value>350,np.nan,abi_value)
        DQF      =np.where(abi_value<150,5,DQF)
        DQF      =np.where(abi_value>350,5,DQF)


        return abi_value, DQF




def window_abi(var,lat,lon,latmin,latmax,lonmin,lonmax,return_latlon=True):
    #Uses the lat and lon to cutoff data in var that is outside the defined box
    
    if return_latlon!=True:
        return var[(lat>=latmin)&(lat<=latmax)&(lon>=lonmin)&(lon<=lonmax)]

    else:
        new_lat = lat[(lat>=latmin)&(lat<=latmax)&(lon>=lonmin)&(lon<=lonmax)]
        new_lon = lon[(lat>=latmin)&(lat<=latmax)&(lon>=lonmin)&(lon<=lonmax)]
        new_var = var[(lat>=latmin)&(lat<=latmax)&(lon>=lonmin)&(lon<=lonmax)]

        return new_lat, new_lon, new_var


def GET_ABI(YYYY,DOY,HH,sat_):
    PATHTOGLM = CONFIG.ABI_PATH+'G'+sat_+'/'

    #Convert Day of year to MMDD
    MMDD = doy_to_mmdd(YYYY,DOY)

    #Load the abi_lat, abi_lon, and pixel area data
    abi_lat, abi_lon, abi_area = get_abi_grid()

    #Grabe the new grid to be used
    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    #This is the meshed grid
    LAT_M_ARRAY, LON_M_ARRAY = np.meshgrid(LAT_ARRAY,LON_ARRAY)

    #Add a buffer to file boundaries to include in windo, helps the nearest neighbor calculation
    lat_min = CONFIG.LAT_MIN - CONFIG.LAT_RES*2.
    lat_max = CONFIG.LAT_MAX + CONFIG.LAT_RES*2.
    lon_min = CONFIG.LON_MIN - CONFIG.LON_RES*2.
    lon_max = CONFIG.LON_MAX + CONFIG.LON_RES*2.

    #File list to grab
    FILELIST = sorted(glob.glob(PATHTOGLM+YYYY+'/'+MMDD+'/OR_ABI-L1b-RadF-M*C13_G'+sat_+'_s'+YYYY+DOY+HH+'*'))

    #Define the arrays of data to save the ABI data to (another array is for the data quality flag
    WHOURLY_ABI = np.zeros((len(FILELIST),LON_ARRAY.size,LAT_ARRAY.size))
    WHOURLY_ABI[:,:,:] = np.nan
    WHOURLY_DQF = np.zeros((len(FILELIST),LON_ARRAY.size,LAT_ARRAY.size))
    WHOURLY_DQF[:,:,:] = np.nan

    print('Processing %i ABI files'%len(FILELIST))

    for i,FILE_ in enumerate(FILELIST): #Loop through the files
        ABI_VAL, ABI_DQF = read_abi(FILE_) #This code grabs the ABI data and converts to a brightness temperature

        try:
            ABI_VAL, ABI_DQF = QC_abi(ABI_VAL, ABI_DQF, abi_area) #This function QCs the data (if more than 10% bad then everything is nan'd
        except:
            print('Unable to QC ABI data')
            continue
        #Now we don't need to regrid all the data, so window the data to get only our domain
        ABI_LAT, ABI_LON, ABI_VAL = window_abi(ABI_VAL,abi_lat,abi_lon,lat_min,lat_max,lon_min,lon_max,return_latlon=True)
        ABI_DQF = window_abi(ABI_DQF,abi_lat,abi_lon,lat_min,lat_max,lon_min,lon_max,return_latlon=False)

        #Now use the subsetted data from the window to regrid the data using a nearest neighbor technique
        ABI_VAR_G = GridData((ABI_LON,ABI_LAT),ABI_VAL,(LON_M_ARRAY,LAT_M_ARRAY),method='nearest')
        ABI_DQF_G = GridData((ABI_LON,ABI_LAT),ABI_DQF,(LON_M_ARRAY,LAT_M_ARRAY),method='nearest')

        #Now save the regridded data to defined arrays
        WHOURLY_ABI[i,:,:] = ABI_VAR_G
        WHOURLY_DQF[i,:,:] = ABI_DQF_G

    try: 
        #Calculate the minimum, mean, std, max TB with the max DQF to send to output
        MINM_BT = (np.nanmin(WHOURLY_ABI,axis=0))
        MEAN_BT = (np.nanmean(WHOURLY_ABI,axis=0))
        STDV_BT = (np.nanstd(WHOURLY_ABI,axis=0))
        MAXM_BT = (np.nanmax(WHOURLY_ABI,axis=0))
        MAX_DQF = (np.nanmax(WHOURLY_DQF,axis=0))
    
    except:
        #If there is no ABI data, we can just output the field of nans
        MINM_BT = np.full((LON_ARRAY.size,LAT_ARRAY.size),np.nan)
        MEAN_BT = np.full((LON_ARRAY.size,LAT_ARRAY.size),np.nan)
        STDV_BT = np.full((LON_ARRAY.size,LAT_ARRAY.size),np.nan)
        MAXM_BT = np.full((LON_ARRAY.size,LAT_ARRAY.size),np.nan)
        MAX_DQF = np.full((LON_ARRAY.size,LAT_ARRAY.size),np.nan)


    return MINM_BT, MEAN_BT, STDV_BT, MAXM_BT, MAX_DQF







if __name__ == '__main__':
    #Start the timer
    start = datetime.datetime.utcnow()

    #Define year(yyyy) and day of year (ddd) to process
    yyyy='2018'
    ddd=list(np.linspace(200,200,1).astype(int))

    #Select the hours of the day for which to process data
    hours=[str(val).zfill(2) for val in np.arange(0,24,1)]


    #Load the abi_lat, abi_lon, and pixel area data
    abi_lat, abi_lon, abi_area = get_abi_grid()

    #Grabe the new grid to be used
    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    #This is the meshed grid
    LAT_M_ARRAY, LON_M_ARRAY = np.meshgrid(LAT_ARRAY,LON_ARRAY)

    lat_min = CONFIG.LAT_MIN - CONFIG.LAT_RES*2.
    lat_max = CONFIG.LAT_MAX + CONFIG.LAT_RES*2.
    lon_min = CONFIG.LON_MIN - CONFIG.LON_RES*2.
    lon_max = CONFIG.LON_MAX + CONFIG.LON_RES*2.


    #Loop through all days that are being processed
    for dddi in ddd:
        dddi2 = str(dddi).zfill(3)

        #Loop through the hours to process.
        for hh in hours:
            MINM_BT, MEAN_BT, STDV_BT, MAXM_BT, MAX_DQF = GET_ABI(yyyy,dddi2,hh)


    print("Finished in %s mins"%((datetime.datetime.utcnow()-start).total_seconds()/60.))

