#Here is my code to process WWLLN data, as well as acquiring SST and AOD information associated with each lightning strike
import numpy as np
import pandas as pd

from Functions import *
from Get_GLM import bin_statistic

def read_WWLLN_Aloc(filename,domain_condition=True):
    # This function reads in WWLLN daily lightning file and outputs the Hour, Lat and Lon of data
    # If domain_condition =True then we remove lightning outside the domain set in the Config file

    #read the csv file
    WWLLN=pd.read_csv(filename,delimiter=',')

    WWLLN = WWLLN.to_numpy() #convert to a numpy array

    #Date=WWLLN[:,0] #The dates are all the same since the files are daily
    Time=WWLLN[:,1]
    Hours = np.array([int(TT[:2]) for TT in Time]) #get only the hour component of time
    Lat=WWLLN[:,2].astype(float)
    Lon=WWLLN[:,3].astype(float)
    Lon=np.where(Lon<0,360+Lon,Lon) #Convert LonW from -180 to 180 to 0 to 360
    #ErrW0=WWLLN[:,4] #time residual less than 30 microseconds
    #NumW0=WWLLN[:,5] #number of stations observing groups

    if domain_condition==True: #Isolates the data within the domain set by the CONFIG file
        #Define the boundaries, extend the grid slightly to allow for centered counts within the grid boxes
        lat_min = CONFIG.LAT_MIN - CONFIG.LAT_RES/2.
        lat_max = CONFIG.LAT_MAX + CONFIG.LAT_RES/2.
        lon_min = CONFIG.LON_MIN - CONFIG.LON_RES/2.
        lon_max = CONFIG.LON_MAX + CONFIG.LON_RES/2.

        #Set the condition to apply to the returned data which is soley based on the domain
        condition_ = np.squeeze([((Lat<=lat_max)&(Lat>=lat_min)&(Lon<=lon_max)&(Lon>=lon_min))])

    else:
        condition_ = np.squeeze([Lat>-999]) #Should be all values (i.e. should be global data)

    return Hours[condition_],Lat[condition_],Lon[condition_]





if __name__ == '__main__':
    #-------------------

    YYYY   = '2020'
    DOYS   = range(360,366)
    HOURS_ = range(10,12)
    #-------------------

    #Grab the grid we will use defined in the Namelist
    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    #Loop through the days of the year
    for i in DOYS:
        MMDD = (doy_to_mmdd(YYYY,i))
        print('Starting for %s'%MMDD)


        if int(YYYY)<=2022:
            WHOURS,WLAT,WLON = read_WWLLN_Aloc(CONFIG.WLN_PATH+YYYY+'/A'+YYYY+MMDD+'.loc')
        #!! Files moved to use the path above
#        elif int(YYYY)==2020:
#            WHOURS,WLAT,WLON =read_WWLLN_Aloc(CONFIG.WLN_PATH+'/A'+YYYY+MMDD+'.loc')


        for HOUR in HOURS_:
            

            l_lat=WLAT[WHOURS==int(HOUR)]
            l_lon=WLON[WHOURS==int(HOUR)]

            #This will grid the data and sum the total in the bins
            COUNT = bin_statistic(l_lat,l_lon,l_lon,LAT_BINS,LON_BINS,'count')
            print('Total number of groups %f'%np.nansum(COUNT))




