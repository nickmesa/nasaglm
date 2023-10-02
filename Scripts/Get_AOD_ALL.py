import numpy as np
from netCDF4 import Dataset
import pygrib as pyg
import glob
import datetime
from scipy.interpolate import griddata as GridData
#-----------------------------------
from Get_GLM import *
from Functions import *

#At 1200 UTC August 16, 2019, the AOD resolution goes from 1 degree to 0.5 degrees. 


def AOD_Filepath(YYYY,DOY,HH,AODType):
#    AODTypes = ['','_du','_sa','_sm','_su']
    MMDD = doy_to_mmdd(YYYY,DOY)
    
    AOD_PATH = CONFIG.AOD_PATH

    #Data prior to this day is in a different resolution (different filenames)
    File_String = 'US058GMET-GR1mdl.*_*'

    filename = glob.glob(AOD_PATH+YYYY+'/'+YYYY+MMDD+'/'+File_String+YYYY+MMDD+HH+'_0001_000000-000000aero_optdep'+AODType)
    
    if len(filename)==1:
        return filename[0]




def pull_aod_data(filename):
    AOD_data = pyg.open(filename)
    grb = AOD_data.select()[0]
    aod = grb.values
    lat_, lon_ = grb.latlons()
    return lat_, lon_, aod


def get_utctime_sandwhich(current_dt):
    '''
    Returns the most recent utc time at the 6 hour mark
    
    provide a datetime object
    #Requires
    import datetime
    import numpy as np
    '''
    #Get current time rounded down to nearest hour
    current = current_dt.replace(microsecond=0, second=0, minute=0)

    #Declare the 4 potential utc 6 hour intervals to choose from 
    ptimes = np.array([current.replace(hour=0),current.replace(hour=6),current.replace(hour=12),current.replace(hour=18)])
    
    previous_time = ptimes[ptimes<=current].max()
    next_time = previous_time + datetime.timedelta(hours=6)
    return previous_time, next_time

def Get_Interp_AOD(YYYY,DOY,HOUR,AODType='_du'):
#    AODTypes = ['','_du','_sa','_sm','_su']
    
    #Read in the data
    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    LAT_M_ARRAY, LON_M_ARRAY = np.meshgrid(LAT_ARRAY,LON_ARRAY)

    HH = str(HOUR).zfill(2)
    DOY= str(DOY).zfill(3)

    print('Getting AOD Data and Interpolating')
    #First check the hours, if it is synoptic , then no interpolation is needed
    if (int(HOUR)%6==0)|(int(HOUR)%6==6):
#        Cfilename=aod_filepath(YYYY,DOY,HH) #Get the filename (weird structure in alternating strings for different times
        Cfilename=AOD_Filepath(YYYY,DOY,HH,AODType)
        try:
            C_LAT, C_LON, C_AOD = pull_aod_data(Cfilename) #Pull the data
        except: #if there is no data 12 hours into the past then say nan
            AOD_FILL = np.zeros(LAT_M_ARRAY.shape)
            AOD_FILL[:] = np.nan
            return AOD_FILL


        #Now regrid the data to defined grid in Nameslist
        C_AOD_G = GridData((C_LON.flatten(),C_LAT.flatten()),C_AOD.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')
        return C_AOD_G

    else: #Here we will have to interpolate the data
        # First we have to get the current time as a datetime to accurately get previous times
        Ctime = datetime.datetime.strptime(str(YYYY)+DOY+HH,'%Y%j%H')
        P_UTC, N_UTC = get_utctime_sandwhich(Ctime) #Gets last synoptic and next time

        try: #Tries to get the most recent file
            #From the data pull the year, day of year and hour
            Pyyyy = P_UTC.strftime('%Y')
            Pdddi = P_UTC.strftime('%j')
            Phh   = P_UTC.strftime('%H')

#            Pfilename=aod_filepath(Pyyyy,Pdddi,Phh) #Get the filename
            Pfilename=AOD_Filepath(Pyyyy,Pdddi,Phh,AODType) #Get the filename

            P_LAT, P_LON, P_AOD = pull_aod_data(Pfilename) #Get the AOD data
            #REgrid the data
            P_AOD_G =  GridData((P_LON.flatten(),P_LAT.flatten()),P_AOD.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')

        except: #If there is an error getting the 6 hour old file, then try and grab the 12 hour prior file
            P_UTC = P_UTC - datetime.timedelta(hours=6)

            #From the data pull the year, day of year and hour
            Pyyyy = P_UTC.strftime('%Y')
            Pdddi = P_UTC.strftime('%j')
            Phh   = P_UTC.strftime('%H')
                
#            Pfilename=aod_filepath(Pyyyy,Pdddi,Phh) #Get the filename
            Pfilename=AOD_Filepath(Pyyyy,Pdddi,Phh,AODType) #Get the filename

            try:
                P_LAT, P_LON, P_AOD = pull_aod_data(Pfilename) #Get the AOD data
            except: #if there is no data 12 hours into the past then say nan
                AOD_FILL = np.zeros(LAT_M_ARRAY.shape)
                AOD_FILL[:] = np.nan
                return AOD_FILL

            
            #REgrid the data
            P_AOD_G =  GridData((P_LON.flatten(),P_LAT.flatten()),P_AOD.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')

        try: #Now lets get the AOD for the next synoptic time

            #From the data pull the year, day of year and hour
            Nyyyy = N_UTC.strftime('%Y')
            Ndddi = N_UTC.strftime('%j')
            Nhh   = N_UTC.strftime('%H')
  
            #Since the AOD files have a weird naming convention
#            Pfilename=aod_filepath(Pyyyy,Pdddi,Phh)
            Nfilename=AOD_Filepath(Nyyyy,Ndddi,Nhh) #Get the filename

            # Now we have the data
            N_LAT, N_LON, N_AOD = pull_aod_data(Nfilename)

            #Lets regrid the data first, otherwise we will get errors on the day where the resolution changed
            N_AOD_G = GridData((N_LON.flatten(),N_LAT.flatten()),N_AOD.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')

        except: #Now lets get the AOD for the next synoptic time
            N_UTC = N_UTC + datetime.timedelta(hours=6)
            #From the data pull the year, day of year and hour
            Nyyyy = N_UTC.strftime('%Y')
            Ndddi = N_UTC.strftime('%j')
            Nhh   = N_UTC.strftime('%H')

            #Since the AOD files have a weird naming convention
#            Pfilename=aod_filepath(Pyyyy,Pdddi,Phh)
            Nfilename=AOD_Filepath(Nyyyy,Ndddi,Nhh,AODType) #Get the filename

            # Now we have the data
            try:
                N_LAT, N_LON, N_AOD = pull_aod_data(Nfilename)
            except: #if there is no data 12 hours in the future or current then say nan
                AOD_FILL = np.zeros(LAT_M_ARRAY.shape)
                AOD_FILL[:] = np.nan
                return AOD_FILL


            #Lets regrid the data first, otherwise we will get errors on the day where the resolution changed
            N_AOD_G = GridData((N_LON.flatten(),N_LAT.flatten()),N_AOD.flatten(),(LON_M_ARRAY,LAT_M_ARRAY),method='linear')


        #Now lets calculate the interpolated AOD data
        try:
            DeltaT    = np.abs((P_UTC - N_UTC).total_seconds()/60./60.)   #time between 2 periods in hours (should be 6 if no data is missing
            Timweight = int(HOUR)%int(CONFIG.HR_RES) #The modulo of hte hour in 6 hourly intervals
            Coef      = (DeltaT-Timweight)/DeltaT #Get the linear weighting coefficient
            
            return  Coef*P_AOD_G + (1.-Coef)* N_AOD_G

        except: #If there was an error calculating the mean, it was probably because other data was missing
            print('AOD Data is missing')
            AOD_FILL = np.zeros(LAT_M_ARRAY.shape)
            AOD_FILL[:] = np.nan
            return AOD_FILL










if __name__ == '__main__':

    #Define year(yyyy) and day of year (ddd) to process
    yyyy='2018'
    ddd=np.arange(100,110,1)

    #We add a buffer around the edge to make sure the data at the edges doesn't exclude data
    lat_min = CONFIG.LAT_MIN - CONFIG.LAT_RES/2.
    lat_max = CONFIG.LAT_MAX + CONFIG.LAT_RES/2.
    lon_min = CONFIG.LON_MIN - CONFIG.LON_RES/2.
    lon_max = CONFIG.LON_MAX + CONFIG.LON_RES/2.


    AODTypes = ['','_du','_sa','_sm','_su']
    AODNames = ['Total','Dust','Salt','Smoke','Sulfate']


    for dddi in ddd:

        for h in range(0,13):
            data = Get_Interp_AOD(yyyy,dddi,h)
            print(np.nanmin(data),np.nanmean(data),np.nanmax(data))
            print(yyyy,dddi,h,data.shape)


    
            



