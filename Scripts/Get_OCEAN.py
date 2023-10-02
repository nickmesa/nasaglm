import numpy as np
from netCDF4 import Dataset
import gzip
import glob
import datetime
from scipy.interpolate import griddata as GridData
from pathlib import Path
#-----------------------------------
# Import Functions from Other Scripts
from Regrid_DTL import *
from Get_GLM import *
from Functions import *
'''

This code has a few functions that reads in SST and OHC. These variables are left in a daily format and regridded to the domain of interest.

Functions are also provided for a distance to land file and a landmask file. Only the distance to land data is output along with the OHC and SST 

The main function that should be called is Regrid_OCEAN('year','day of year') to output DTL, SST, and OHC

'''

def get_landmask():
    mask_file="../Data/lsmask.oisst.v2.nc"
    lsmask=Dataset(mask_file)
    #This  must be changed based on the SST domain.
    m_lat = lsmask.variables['lat'][:]
    m_lon = lsmask.variables['lon'][:]
    mask = lsmask.variables['lsmask'][0][:][:].data
    lsmask.close()

    return m_lon, m_lat, mask



def get_sst(yyyy,doy):
    mmdd = doy_to_mmdd(yyyy,doy)
    #Load in the sst file
    SST_data = np.loadtxt(CONFIG.SST_PATH+yyyy+'/dsst_glb_'+yyyy+mmdd+'.dat',skiprows=1)

    #Acquire the longitude, latitude, and SST data from the file.
    lon_sst =np.unique(SST_data[:,0])
    lat_sst =np.unique(SST_data[:,1])

    #Reshape the SST data to 2D array.
    SST=np.reshape(SST_data[:,2],(lat_sst.size,lon_sst.size))
    SST=np.array(SST,dtype=np.float32)

    return lon_sst, lat_sst, SST


def get_ohc(yyyy,doy):
    mmdd = doy_to_mmdd(yyyy,doy)

    OHC_name=CONFIG.OHC_PATH + yyyy + '/eohc_glb_' + yyyy + mmdd + '.dat.gz'
    
    if Path(OHC_name).is_file()==False:
        print('File Does not exist: trying to use unzipped file')    
        OHC_name=CONFIG.OHC_PATH + yyyy + '/eohc_glb_' + yyyy + mmdd + '.dat'

    lon = []
    lat = []
    ohc_    = []

    if OHC_name[-3:]=='.gz':
        #The OHC data is gzipped, so we need to open those files and read line by line
        with gzip.open(OHC_name,'rt') as f:
            for i,line in enumerate(f.readlines()):
                if i==0: #This is the header
                    continue
                #Parse out the line
                data_line = np.array([(x.strip(' ').rstrip()) for x in line.split()])
                #save the data to lists
                lon.append(float(data_line[0]))
                lat.append(float(data_line[1]))
                ohc_.append(float(data_line[17]))

    elif OHC_name[-3:]=='dat':
        #The OHC data is gzipped, so we need to open those files and read line by line
        with open(OHC_name,'rt') as f:
            for i,line in enumerate(f.readlines()):
                if i==0: #This is the header
                    continue
                #Parse out the line
                data_line = np.array([(x.strip(' ').rstrip()) for x in line.split()])
                #save the data to lists
                lon.append(float(data_line[0]))
                lat.append(float(data_line[1]))
                ohc_.append(float(data_line[17]))
    else:
        print('Unknown File Extension for OHC')


    ohc_lon = np.unique(lon)
    ohc_lat = np.unique(lat)

    ohc_ = np.array(ohc_).reshape(ohc_lat.size,ohc_lon.size)

    #return as numpy arrays
    return np.array(ohc_lon),np.array(ohc_lat), ohc_.T


def Regrid_OCEAN(YYYY,DOY):
    
    #Read in the grid, get the meshed grid, and specify the boundaries for faster processing
    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    LAT_M_ARRAY, LON_M_ARRAY = np.meshgrid(LAT_ARRAY,LON_ARRAY)

    lat_min = CONFIG.LAT_MIN - CONFIG.LAT_RES/2.
    lat_max = CONFIG.LAT_MAX + CONFIG.LAT_RES/2.
    lon_min = CONFIG.LON_MIN - CONFIG.LON_RES/2.
    lon_max = CONFIG.LON_MAX + CONFIG.LON_RES/2.

    DOY= str(DOY).zfill(3)
    MMDD = doy_to_mmdd(YYYY,DOY)

#    NOTE-- REMOVED LANDMASK FROM CALCULATION (Redundant with DTL)
#    #----------------------------------
#    #  Grab and regrid the Landmask
#    #
#    #----------------------------------
#
#    #First lets grab the landmask data
#    m_lon, m_lat, mask = get_landmask()
#    #Now we will isolate only the data within the domain that we want
#    m_lon, m_lat, mask = window_var(mask,m_lat,m_lon,lat_min,lat_max,lon_min,lon_max,ret_latlon=True)
#
#    lonm_mesh,latm_mesh=np.meshgrid(m_lon,m_lat)
#    MASK =  GridData((lonm_mesh.flatten(),latm_mesh.flatten()),mask.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')

    #----------------------------------
    #  Grab and regrid the SST
    #
    #----------------------------------
    #Now lest get the SST data
    try:
        sst_lon, sst_lat, sst_ = get_sst(YYYY,DOY)
        #Now we will isolate only the data within the domain that we want
        sst_lon, sst_lat, sst_ = window_var(sst_,sst_lat,sst_lon,lat_min,lat_max,lon_min,lon_max,ret_latlon=True)
        #REgrid the data

        lon_mesh,lat_mesh=np.meshgrid(sst_lon,sst_lat)
        SST =  GridData((lon_mesh.flatten(),lat_mesh.flatten()),sst_.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')

    except: #Something wrong with the SST data
        SST = np.zeros(LAT_M_ARRAY.shape)
        SST[:] = np.nan

    #----------------------------------
    #  Grab and regrid the OHC
    #
    #----------------------------------
    try:
        ohc_lon, ohc_lat, ohc_ = get_ohc(YYYY,DOY)
        #Now we will isolate only the data within the domain that we want
        ohc_lon, ohc_lat, ohc_ = window_var(ohc_,ohc_lat,ohc_lon,lat_min,lat_max,lon_min,lon_max,ret_latlon=True)

        lonc_mesh,latc_mesh=np.meshgrid(ohc_lon,ohc_lat)
        OHC =  GridData((lonc_mesh.flatten(),latc_mesh.flatten()),ohc_.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')

    except: #Something wrong with the SST data
        OHC = np.zeros(LAT_M_ARRAY.shape)
        OHC[:] = np.nan

    #----------------------------------
    #  Grab and regrid the DTL
    #
    #----------------------------------
    dtl_lon,dtl_lat,dtl_ = get_01res_dtl()

    #Cut off the data outside domain of interest
    dtl_lon, dtl_lat, dtl_ = window_var(dtl_,dtl_lat,dtl_lon,lat_min,lat_max,lon_min,lon_max,ret_latlon=True)

    #Turn the list of unique lat/lons into a grid
    lond_mesh,latd_mesh=np.meshgrid(dtl_lon,dtl_lat)
    DTL =  GridData((lond_mesh.flatten(),latd_mesh.flatten()),dtl_.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')


    return SST, OHC, DTL








if __name__ == '__main__':

    #Define year(yyyy) and day of year (ddd) to process
    yyyy='2021'
    ddd=np.arange(100,110,1)


    for dddi in ddd:
#        lon,lat,ohc = get_ohc(yyyy,dddi)
#        print(np.nanmin(lon),np.nanmax(lon),np.nanmin(lat),np.nanmax(lat),np.nanmax(ohc))
        SST, OHC, DTL = Regrid_OCEAN(yyyy,dddi)
        print(np.nanmax(OHC))   


    
            



