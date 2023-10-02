import numpy as np
from netCDF4 import Dataset
import datetime
import glob
from matplotlib import pyplot as plt
from scipy.interpolate import griddata as GridData
#import custom functions
from Functions import *
from Get_GLM import bin_statistic
#-------------------------------
#
# This code has the functions to read in and process LIS data
# Ben Trabing, CIRA/CSU/NHC
# Ben.Trabing@noaa.gov
#
#-------------------------------

def read_LIS(filename_,domain_condition=True,show_image=False):
    # If domain_condition =True then we remove lightning outside the domain set in the Config file
    if type(filename_) is not list:
        filename_ = [filename_]

    print('Processing %i LIS Files'%len(filename_))    

    for filename in filename_:
        ds = Dataset(filename)

        g_lat=ds.variables['lightning_group_lat'][:].data
        g_lon=ds.variables['lightning_group_lon'][:].data
        g_rad=ds.variables['lightning_group_radiance'][:].data   #units of "uJ/sr/m2/um" 
        g_foot=ds.variables['lightning_group_footprint'][:].data  #units of km2
        g_time=ds.variables['lightning_group_TAI93_time'][:].data
        g_flag=ds.variables['lightning_group_alert_flag'][:].data

        v_lat=ds.variables['viewtime_lat'][:].data
        v_lon=ds.variables['viewtime_lon'][:].data
        v_start=ds.variables['viewtime_TAI93_start'][:].data
        v_end=ds.variables['viewtime_TAI93_end'][:].data
        
        if show_image!=False:
            #If show_image is not false, then read in and show the raster_image
            img=ds.variables['raster_image'][:].data
            fig = plt.figure(figsize=(12,8))
            plt.imshow(img)
            plt.show()

        #Close the netcdf file
        ds.close()
        
        return g_lat, g_lon, g_time, g_rad, g_foot,g_flag, v_lat, v_lon, v_start, v_end

def convert_nd_datetime(arr):
    #if the time arrays have any nans, convert them to numbers
    arr[np.isnan(arr)]=0

    start = datetime.datetime(1993,1,1,0,0,0) #start of LIS data
    
    arrx = np.copy(arr) #make a copy of the times
    arrx = arrx.astype('datetime64[s]') #convert to datetime units
    
    #Now we need to loop through the data and add the start datetime with the seconds since variable
    # datetimes don't play nice with numpy arrays, which is why we need a loop here
    if len(arr.shape)==1:
        for i in range(arr.shape[0]):
            arrx[i] = start + datetime.timedelta(seconds=int(arr[i]))

    #if we have a 2d array, we will loop through the dates in each dimension
    if len(arr.shape)==2:
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                arrx[i,j] = start + datetime.timedelta(seconds=int(arr[i,j]))

    return arrx #return the timing array in a datetime format


def LIS_path(filenames,YYYY,MM,DD,HH):
    #
    # Note that this data is being converted to hourly
    #
    
    #Read in the grid data in namelist
    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    LAT_M_ARRAY, LON_M_ARRAY = np.meshgrid(LAT_ARRAY,LON_ARRAY)
    
    #define the time values for window using the HH
    window_start = datetime.datetime(YYYY,MM,DD,HH,0,0)
    window_end = datetime.datetime(YYYY,MM,DD,HH,0,0)+datetime.timedelta(hours=1)


    ####
    # Define a .5deg x .5 deg grid that is the resolution of viewtime variables
    LIS_LATS=np.arange(-65,65.5,.5)
    LIS_LONS=np.arange(0,360,.5)
    LIS_MX,LIS_MY = np.meshgrid(LIS_LONS[:-1],LIS_LATS[:-1])
    
    ####
    
    Total_MASK = np.zeros(LON_M_ARRAY.shape)
    Total_MASK[:]= -1 #set mask to negative ones (i.e. no overpass)
    Total_LGRP = np.zeros(LON_M_ARRAY.shape)
    
    for filename in filenames:
        #Open the file
        ##################
        g_lat, g_lon, g_time, g_rad, g_foot,g_flag, v_lat, v_lon, v_start, v_end= read_LIS(filename,show_image=False)

        #fix the longitude data which is from -180->+180
        v_lon[v_lon<0] = v_lon[v_lon<0]+360.
        g_lon[g_lon<0] = g_lon[g_lon<0]+360.

        #convert times to datetimes
        g_time = (convert_nd_datetime(g_time))
        
        # create the timing condition for the lightning data to only get data inside the time window we want
        condition_L = np.squeeze([(g_time<window_end)&(g_time>=window_start)])
        
        try:
            #Calculate the group count within the config grid
            Swath_groups = bin_statistic(g_lat[condition_L],g_lon[condition_L],g_time[condition_L],LAT_BINS,LON_BINS,'count')
#            Total_LGRP[Swath_groups>0]+=Swath_groups[Swath_groups>0]
        except:
            #no lightning data so nothing to compute here and add to the array
            pass
        
        

        #Bin the viewtime locations in space using the .5 x.5 degree grid
        dat = bin_statistic(v_lat,v_lon,v_lon,LIS_LATS,LIS_LONS,'count')
        dat[dat<=0]=0 #Values of zero indicate no overpass
        dat[dat>0]=1 #values greater than zero means grid viewed


        #Convert .5x.5 degree LIS grid to .1x.1 degree using nearest neighbor
        VIEW_MASK =  GridData((LIS_MX.flatten(),LIS_MY.flatten()),dat.T.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='nearest')

        #linearly interpolate the ending time
        timdataS =  GridData((v_lon.flatten(),v_lat.flatten()),v_start.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')
        timdataE =  GridData((v_lon.flatten(),v_lat.flatten()),v_end.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')

        timdataS[VIEW_MASK<1] = np.nan #We dont want times where there was no overpass
        timdataE[VIEW_MASK<1] = np.nan #We dont want times where there was no overpass

        TIMEDATAS = convert_nd_datetime(timdataS) #convert the seconds since 1993 to datetimes
        TIMEDATAE = convert_nd_datetime(timdataE) #convert the seconds since 1993 to datetimes

        #Get the total duriation of the overpass over each point
        DURATION = (timdataE-timdataS)
        
        #Scale Data to be hourly using Duration- 
        OFFSET = 3600./DURATION
        Swath_groups = Swath_groups*OFFSET

        #Add the group counts with the offset to hourly
        Total_LGRP[Swath_groups>0]+=Swath_groups[Swath_groups>0]


        # create the condition for the data to only get data inside the time window we want
        condition_T = np.squeeze([(TIMEDATAS<window_end)&(TIMEDATAE>=window_start)])

        #convert data outside the time window to nan
        VIEW_MASK[~condition_T] = np.nan
        
        #Now take the masked data and if there was an overpass then add in the path data
        Total_MASK[VIEW_MASK==0]=0
        Total_MASK[VIEW_MASK>0]=1

    #Now that we have a final mask use it to remove zero group counts where no observations were made
    Total_LGRP[Total_MASK<0]= -1

    #return mask and lightning group count within the hour
    return Total_MASK,Total_LGRP




if __name__ == '__main__':
    #-------------------

    YYYY   = '2019'
    DOYS   = range(260,266)
    HOURS_ = range(0,3)
    #-------------------

    #Grab the grid we will use defined in the Namelist
    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    #Loop through the days of the year
    for i in DOYS:
        MMDD = (doy_to_mmdd(YYYY,i))
        print('Starting for %s'%MMDD)
        MM = MMDD[:2]
        DD = MMDD[2:]

        try: #This will not work for the first or last day of the year 
            MMDD_y = (doy_to_mmdd(YYYY,i-1))
            #Get the list of files from yesterday and grab the latest 2
            Filelist_yesterday = sorted(glob.glob(CONFIG.LIS_PATH+YYYY+'/ISS_LIS_SC_V1.0_'+YYYY+MMDD_y+'*'))[-2:]

            MMDD_t = (doy_to_mmdd(YYYY,i+1))
            #Get the list of files from tomorrow and grab the first 2
            Filelist_tomorrow = sorted(glob.glob(CONFIG.LIS_PATH+YYYY+'/ISS_LIS_SC_V1.0_'+YYYY+MMDD_t+'*'))[:2]

            #Combine the list of files
            Filelist = Filelist_yesterday+sorted(glob.glob(CONFIG.LIS_PATH+YYYY+'/ISS_LIS_SC_V1.0_'+YYYY+MMDD+'*')) +Filelist_tomorrow

        except:
            Filelist = sorted(glob.glob(CONFIG.LIS_PATH+YYYY+'/ISS_LIS_SC_V1.0_'+YYYY+MMDD+'*'))



  
        for HH in HOURS_:


            VMASK,LIS_GRP = LIS_path(Filelist,int(YYYY),int(MM),int(DD),int(HH))

            print(VMASK.shape,LIS_GRP.shape)




