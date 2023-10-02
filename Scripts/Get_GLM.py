import datetime
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from netCDF4 import Dataset
import matplotlib
import glob
import os

from Functions import *


def window_var(var,lat,lon,latmin,latmax,lonmin,lonmax,ret_latlon=False):
    #Uses the lat and lon to cutoff data in var that is outside the defined box
    
    #conditionals for lat and lon to apply to the var
    lon_con = np.squeeze([(lon>=lonmin)&(lon<=lonmax)])
    lat_con = np.squeeze([(lat>=latmin)&(lat<=latmax)])

     
    if ret_latlon!=False: #If we want to return the lat/lon with conditions
        try: #If you get an error here, then you may need to switch the lat/lon placement
            return lon[lon_con], lat[lat_con], var.T[lat_con][:,lon_con]
        except:
            return lon[lon_con], lat[lat_con], var[lat_con][:,lon_con]

    else: #Just return the conditioned variable
        try: #If you get an error here, then you may need to switch the lat/lon placement
            return var.T[lat_con][:,lon_con]
        except:
            return var[lat_con][:,lon_con]



def read_glm_filelist(filelist, apply_qc=True,output_qc=False):
    '''
    This function will read in the filelist of glm netcdf files and then output arrays of the lightning groups





    '''

    #For Groups
    lats_l = []
    lons_l = []
    ener_l = []
    area_l = []
    flag_l = []    
    flid_l = []

    #For Flashes
    lats_lf = []
    lons_lf = []
    ener_lf = []
    area_lf = []
    flag_lf = []   
    flnum_lf = []

    if len(filelist)<1:
        return np.array(lats_l), np.array(lons_l), np.array(ener_l), np.array(area_l), np.array(flag_l), np.array(flid_l), np.array(lats_lf), np.array(lons_lf), np.array(ener_lf), np.array(area_lf), np.array(flag_lf)


    for filename in filelist:
        GLM_file = Dataset(filename)

        #Area units were changed in 2019 from km2 to m2 so convert to m2
        area_unit = GLM_file.variables['group_area'].units
        
        if area_unit=='m2':
            area_conv = 1.
        elif area_unit=='km2':
            area_conv=1000.**2

        temp_lats_g = GLM_file.variables['group_lat'][:]
        temp_lons_g = GLM_file.variables['group_lon'][:] + 360. #add 360 to convert to similar format as other grids
        temp_flag_g = GLM_file.variables['group_quality_flag'][:]

        temp_lats_f = GLM_file.variables['flash_lat'][:]
        temp_lons_f = GLM_file.variables['flash_lon'][:] + 360. #add 360 to convert to similar format as other grids
        temp_flag_f = GLM_file.variables['flash_quality_flag'][:]


        cond_g = np.squeeze([((temp_flag_g==0)&(temp_lats_g<=CONFIG.LAT_MAX)&(temp_lats_g>=CONFIG.LAT_MIN)&(temp_lons_g<=CONFIG.LON_MAX)&(temp_lons_g>=CONFIG.LON_MIN))])
        cond_f = np.squeeze([((temp_flag_f==0)&(temp_lats_f<=CONFIG.LAT_MAX)&(temp_lats_f>=CONFIG.LAT_MIN)&(temp_lons_f<=CONFIG.LON_MAX)&(temp_lons_f>=CONFIG.LON_MIN))])

        #For the group quantities
        lats_l.append(list(temp_lats_g[cond_g]))
        lons_l.append(list(temp_lons_g[cond_g]))
        ener_l.append(list(GLM_file.variables['group_energy'][:][cond_g]))
        area_l.append(list(area_conv*GLM_file.variables['group_area'][:][cond_g]))
        flag_l.append(list(temp_flag_g[cond_g]))
        flid_l.append(list(GLM_file.variables['group_parent_flash_id'][:][cond_g]))

        #For flash information
        lats_lf.append(list(temp_lats_f[cond_f]))
        lons_lf.append(list(temp_lons_f[cond_f]))   
        ener_lf.append(list(GLM_file.variables['flash_energy'][:][cond_f]))
        area_lf.append(list(area_conv*GLM_file.variables['flash_area'][:][cond_f]))
        flag_lf.append(list(temp_flag_f[cond_f]))
        flnum_lf.append(list(GLM_file.variables['flash_id'][:][cond_f]))


        GLM_file.close()
    
    #Reformat the group variables as numpy arrays
    lats_l = np.array(flatten(lats_l))
    lons_l = np.array(flatten(lons_l)) 
    ener_l = np.array(flatten(ener_l))
    area_l = np.array(flatten(area_l))
    flag_l = np.array(flatten(flag_l))
    flid_l = np.array(flatten(flid_l))

    #Reformat the flash variables as numpy arrays
    lats_lf = np.array(flatten(lats_lf))
    lons_lf = np.array(flatten(lons_lf)) 
    ener_lf = np.array(flatten(ener_lf))
    area_lf = np.array(flatten(area_lf))
    flag_lf = np.array(flatten(flag_lf))
    flnum_lf = np.array(flatten(flnum_lf))

    #if we need to qc the data
    if apply_qc==True:
        #Defined Exponential Function to eliminate questionable data
        AA = 0.9266284080291266
        BB = 2.1140353739008016e+20
        
        #Using the exponential, given the energy, what should the threshhold be for area
        Area_thresh = 10.**(AA*np.log10(ener_l)+np.log10(BB))
        Area_threshF = 10.**(AA*np.log10(ener_lf)+np.log10(BB))
        
        #Energy minimum to use as cutoff
        Energy_min = 1.9*1e-15
        
        #These are the QC conditions, above the area thresh and energy min, flag must be zero, and exclude non physical lat/lons
        condition_not = np.squeeze([((area_l<Area_thresh)|(ener_l<Energy_min))])
#        condition_ = np.squeeze([((area_l>=Area_thresh)&(ener_l>=Energy_min)&(flag_l==0)&(lats_l<=CONFIG.LAT_MAX)&(lats_l>=CONFIG.LAT_MIN)&(lons_l<=CONFIG.LON_MAX)&(lons_l>=CONFIG.LON_MIN))])

         #Get the unique flash id for all bad flashes
        bad_flash = flid_l[condition_not]

        # get the condition for only non bad identified groups using flash id
        condition_ = np.squeeze([(flid_l!=bad_flash)])

        condition_f = np.array([i for i,val in enumerate(flnum_lf) if val not in np.unique(bad_flash)]) 
        # now condition the flashes using bad flash ids
#        condition_f = np.squeeze([(flnum_lf==valid_flash)])

#        condition_f = np.squeeze([((area_lf>=Area_threshF)&(ener_lf>=Energy_min)&(flag_lf==0)&(lats_lf<=CONFIG.LAT_MAX)&(lats_lf>=CONFIG.LAT_MIN)&(lons_lf<=CONFIG.LON_MAX)&(lons_lf>=CONFIG.LON_MIN))])

        #Uncomment the line below to not include the QC equation
#        condition_ = np.squeeze([(flag_l==0)&(lats_l<=CONFIG.LAT_MAX)&(lats_l>=CONFIG.LAT_MIN)&(lons_l<=CONFIG.LON_MAX)&(lons_l>=CONFIG.LON_MIN)])

        lats_ = np.array(flatten(lats_l[condition_]))
        lons_ = np.array(flatten(lons_l[condition_]))
        ener_ = np.array(flatten(ener_l[condition_]))
        area_ = np.array(flatten(area_l[condition_]))
        flag_ = np.array(flatten(flag_l[condition_]))    
        flid_ = np.array(flatten(flid_l[condition_]))    

        lats_f = np.array((lats_lf[condition_f]))
        lons_f = np.array((lons_lf[condition_f])) 
        ener_f = np.array((ener_lf[condition_f]))
        area_f = np.array((area_lf[condition_f]))
        flag_f = np.array((flag_lf[condition_f]))

     
        return lats_, lons_, ener_, area_, flag_, flid_, lats_f, lons_f, ener_f, area_f, flag_f

    else: #If no QC is applied to the lightning data

        return lats_l, lons_l, ener_l, area_l, flag_l, flid_l, lats_lf, lons_lf, ener_lf, area_lf, flag_lf


def get_d2l(filepath,lat_,lon_,lat_min,lat_max,lon_min,lon_max):
    #uses the path for the d2l file, then converts nans to zeros and shrinks the domain
    
    Dist_land0=np.load(filepath) #load the data
    Dist_land=Dist_land0['dist_land_002'] #get the data from the variable name
    Dist_land[np.isnan(Dist_land)]=0   #convert nans to zeros
    
    #Now output the shrunken data
    return window_var(Dist_land,lat_,lon_,lat_min,lat_max,lon_min,lon_max)



if __name__ == '__main__':
    #-------------------
    
    YYYY = '2019'
    Sat = '16'
    #-------------------

    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    lat_min = CONFIG.LAT_MIN - CONFIG.LAT_RES/2.
    lat_max = CONFIG.LAT_MAX + CONFIG.LAT_RES/2.
    lon_min = CONFIG.LON_MIN - CONFIG.LON_RES/2.
    lon_max = CONFIG.LON_MAX + CONFIG.LON_RES/2.



    #Get the distance to land in the same window frame as the gridded lightning data

    #Loop through the days of the year
    for i in np.arange(100,101):
        MMDD = (doy_to_mmdd(YYYY,i))
        print('Starting for %s'%MMDD)

        #---------------------------------------
        # PATHS TO SET 
        #---------------------------------------        
        GLM_LNG_PATH = CONFIG.GLM_PATH+YYYY+'/'+MMDD+'/OR_GLM-L2-LCFA_G'+Sat+'_s'+YYYY+str(i).zfill(3)
        #---------------------------------------

        
        #Get the daily gridded file lists for the day

        start = 0
        for HOUR in range(23,24):
           
            LIGHTNING_LIST = glob.glob(GLM_LNG_PATH+str(HOUR).zfill(2)+'*')
            print('Processing %i individual files'%len(LIGHTNING_LIST))

            #read in teh lightning group data for single file
            #lats_l, lons_l, ener_l, area_l, flag_l, lats_lF, lons_lF, ener_lF, area_lF, flag_lF = read_glm_filelist(LIGHTNING_LIST, apply_qc=True,output_qc=True)
            lats_l, lons_l, ener_l, area_l, flag_l, flid_l, lats_lf, lons_lf, ener_lf, area_lf, flag_lf = read_glm_filelist(LIGHTNING_LIST, apply_qc=True,output_qc=False)

            LCOUNT      = bin_statistic(lats_l,lons_l,ener_l,LAT_BINS,LON_BINS,'count')
            LFCOUNT     = bin_statistic(lats_lf,lons_lf,ener_lf,LAT_BINS,LON_BINS,'count')

            MEAN_ENERGY = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'mean')
            STDV_ENERGY = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'std')
            MINM_ENERGY = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'min')
            MAXM_ENERGY = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'max')
            MEAN_AREA   = bin_statistic(lats_l[~np.isnan(area_l)],lons_l[~np.isnan(area_l)],area_l[~np.isnan(area_l)],LAT_BINS,LON_BINS,'mean')
            STDV_AREA   = bin_statistic(lats_l[~np.isnan(area_l)],lons_l[~np.isnan(area_l)],area_l[~np.isnan(area_l)],LAT_BINS,LON_BINS,'std')
            MINM_AREA   = bin_statistic(lats_l[~np.isnan(area_l)],lons_l[~np.isnan(area_l)],area_l[~np.isnan(area_l)],LAT_BINS,LON_BINS,'min')
            MAXM_AREA   = bin_statistic(lats_l[~np.isnan(area_l)],lons_l[~np.isnan(area_l)],area_l[~np.isnan(area_l)],LAT_BINS,LON_BINS,'max')

    
    






