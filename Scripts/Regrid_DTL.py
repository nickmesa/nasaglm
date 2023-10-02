import numpy as np
from netCDF4 import Dataset
import os
import glob
import netCDF4

#-------------------------------------------------------------------------
#
# First We will Read in the Original Data
#
#-------------------------------------------------------------------------

def get_001res_dtl():
    #Open the 0.01-degree file
    Data='/Users/nmesa/Desktop/NASA_GLM_Project/climatology_glm-main/Data/dist2coast_1deg.nc'
    Data_grid=Dataset(Data,mode='r')

    #Lat, Lon, and Dist variables
    lat=Data_grid.variables['lat'][:]
    lon=Data_grid.variables['lon'][:]
    dist=Data_grid.variables['dist'][:,:]

    #Adjust lon to be the same as the 0.02 degree grid (0 to 360, not -180 to 180)
#    lon2=np.where(lon<0,lon+360,lon)

    Data_grid.close() #close the file

    #Transopose the array so that it is lon,lat
    return lon[::10],lat[::10],dist.T[::10,::10]


def get_01res_dtl():
    #This is a function to read in the newly created netcdf file
    
    #Open the 0.1-degree file
    Data='/Users/nmesa/Desktop/NASA_GLM_Project/climatology_glm-main/Data/dist2coast_01deg.nc'
    Data_grid=Dataset(Data,mode='r')

    #Lat, Lon, and Dist variables
    lat=Data_grid.variables['lat'][:]
    lon=Data_grid.variables['lon'][:]
    dist=Data_grid.variables['Distance'][:,:]

    Data_grid.close() #close the file

    lon2=np.where(lon<0,lon+360,lon)

    return lon2,lat,dist


#-------------------------------------------------------------------------
#
# First We will Read in the Original Data
#
#-------------------------------------------------------------------------

if __name__ == '__main__':

    lon, lat, dist = get_001res_dtl()


    ncfile =  netCDF4.Dataset('/Users/nmesa/Desktop/NASA_GLM_Project/climatology_glm-main/Data/dist2coast_01deg.nc',mode='w') #open the file to write to

    #Create the dimensions of lat lons
    ncfile.createDimension('lons', len(lon))
    ncfile.createDimension('lats', len(lat))

    lat_ = ncfile.createVariable('lat', 'f', ('lats'),fill_value=-999.0)
    lat_[:] = lat

    lon_ = ncfile.createVariable('lon', 'f', ('lons'),fill_value=-999.0)
    lon_[:] = lon

    dtl_ = ncfile.createVariable('Distance', 'f', ('lons','lats'),fill_value=-999.0)
    dtl_.setncatts({'long_name': u"Distance to Land",\
                    'units': u"km", \
                    'var_desc': u"Regridded from The Pacific Islands Ocean Observing System (PacIOOS)"})
    dtl_[:] = dist

    ncfile.close() #Make sure that we close the netcdf file




