#!/usr/bin/env python
'''

These Functions are intended to help with the processing and evaluation of tropical cyclones (TCs) within
the GLM climatology.

Written by Ben Trabing, March 23, 2022
Contact me: Ben.Trabing@noaa.gov



'''
import numpy as np
import pandas as pd
import subprocess
import re
import datetime
from Functions import *



def read_hur_storm(storm_name,BASIN='AL'):
    #Input the storm name and it will return the best track data as an array for that specific storm
    
    if BASIN=='AL':
        BEST = '../Data/hurdat2-1851-2021-041922.txt'
    elif BASIN=='EP':
        BEST = '../Data/hurdat2-nepac-1949-2021-042522.txt'


    #Names elaborated
    #A      Year Month and Day
    #B      Hours in UTC
    #C      Record Identifier  
    #D      Status of System
    #E      Latitude
    #F      Longitude
    #G      Maximum sustained winds in knots
    #H      Minimum Pressure in millibars
    #I      34 kt wind radii maximum in the northeastern quad in nautical miles
    #J      34 kt wind radii maximum in the southeastern quad
    #K      34 kt wind radii maximum in the southwestern quad
    #L      34 kt wind radii maximum in the northwestern quad
    #M      50 kt wind radii maximum in the northeastern quad
    #N      50 kt wind radii maximum in the southeastern quad
    #O      50 kt wind radii maximum in the southwestern quad
    #P      50 kt wind radii maximum in the northwestern quad
    #Q      64 kt wind radii maximum in the northeastern quad
    #R      64 kt wind radii maximum in the southeastern quad
    #S      64 kt wind radii maximum in the southwestern quad
    #T      64 kt wind radii maximum in the northwestern quad

    #If you are using a best track after 2020, you need to have a 'U' to account for new 'RMW' value
    hur = pd.read_csv(BEST,na_filter=False,names=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U'])


    #use the listed storm name to find the specified storm
    output = hur[hur['A'].str.contains(storm_name)]

    #Find the index number of the desired storm for the larger dataframe
    i = int(output.index.values)+1

    #using the listed number of files, find the last indexed value for the desired storm
    e = i + int(output['C'])

    #create a dataframe of the desired storm
    storm = hur[i:e]

    #Add the date and time columns together to easily put into datetime format
    storm.loc[:,'A'] = storm['A'] + storm['B']


    #simplify the longitude and latitude to floats
    storm.loc[:,'E'] = storm['E'].str[:5].astype(float)
    storm.loc[:,'F'] = storm['F'].str[:6].astype(float)*-1

    return storm.values

def get_stormnames_hurdat(BASIN='AL'):
    #This code will get all the unique storm names within the hurdat2
    if BASIN=='AL':
#        BEST = '../Data/hurdat2-1851-2020-052921.txt'
        BEST = '../Data/hurdat2-1851-2021-041922.txt'
    elif BASIN=='EP':
#        BEST = '../Data/hurdat2-nepac-1949-2020-043021a.txt'
        BEST = '../Data/hurdat2-nepac-1949-2021-042522.txt'


    #If you are using a best track after 2020, you need to have a 'U' to account for new 'RMW' value
    hur = pd.read_csv(BEST,na_filter=False,names=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U'])
 
    namemixdt = hur['A'].values #get the first column where the names are
    
    names = []
    for sname in namemixdt:
        if (sname[0:2].strip()=='AL') & (BASIN=='AL'): #Only if the line has a Storm Identifier
            names.append(sname.strip()) #Add the name to a list

        elif (sname[0:2].strip()=='EP')& ((BASIN=='EP') | (BASIN=='CP')) :
            names.append(sname.strip()) #Add the name to a list

    
    return np.unique(names) #Output the unique names of 



def output_best_track_dataset(BASIN='AL'):
    
    list_of_storms = get_stormnames_hurdat(BASIN) #Get the Atlantc storm names
    
    
    for i,name_ in enumerate(list_of_storms): #loop through every storm name
        
        sdata = read_hur_storm(name_,BASIN) #get the data for each storm name
        
        #Since the lines don't include storm names, we need to add it
        sname_array = np.full((sdata.shape[0]),'TEMP9999') #Create an emtpy array with temp name
        sname_array[:] = name_ #add the stormname to the array
        
        sdata_ = np.concatenate((sname_array.reshape(-1,1),sdata),axis=1) #add the name array to the data array
        
        if i==0:
            DATA=sdata_
        else:
            DATA = np.concatenate((DATA,sdata_),axis=0)
    
        
    return DATA #Output the data which should have everything in the best tracks
    
    

#NHC_BEST_AL = output_best_track_dataset('AL')
#NHC_BEST_EP = output_best_track_dataset('EP')

def Total_Best():
    # Combines AL and EP best track datasets
    # Call as:
    # NHC_BEST = Total_Best()

    NHC_BEST_AL = output_best_track_dataset('AL')
    NHC_BEST_EP = output_best_track_dataset('EP')

    NHC_BEST = np.concatenate((NHC_BEST_AL,NHC_BEST_EP),axis=0)
    
    return NHC_BEST


def get_storms_bytime(DATASET,YYYY,MM,DD,HH,Daily=False):
    '''
    Variable to Set
       Precision_hours : number of hours prior to and after the time of interest to output
    
    
    INPUT
       DATASET : the Best-Track Dataset- use: DATASET = output_best_track_dataset('AL')
       YYYY    : the year of interest
       MM      : the month of interest
       DD      : the day of interest
       HH      : the hour of interest
       Daily   : conditional- if set to True, then the hourly input is not considered
                    and the function will output all positions in the specified day.
                    Note that to get exactly 1 day, the Precision_hours must be set to zero.
                    Otherwise +/- the precision hours from the start and end of the specified 
                    day will be used.
    
    Returns
       Array   : the lines of the Best-Track Dataset within a specified time window
    
    
    '''    
    
    #Number of hours before and after to get best-track positions
    Precision_hours = 7.
    
    if Daily==True:
        #This is 12 hours plus the day of interest specified in the function call
        toi = datetime.datetime(int(YYYY),int(MM),int(DD),12,0)
    else: 
        #This is the time of interest specified in the function call
        toi = datetime.datetime(int(YYYY),int(MM),int(DD),int(HH),0)
    
    #Calculate from the DATASET file the time deviation in hours from the time of interest
    time_array_  = np.abs([(datetime.datetime.strptime(val,'%Y%m%d %H%M')-toi).total_seconds()/60./60. for val in DATASET[:,1]])

    if Daily==True:
        # return the best-track lines that are within +/- 12 hours plus the precision hour window
        # this will allow for the full day to be analyzed
        return DATASET[time_array_<=(12.+Precision_hours),:]

    else:
        #return the best-track lines that are within +/- the precision hour window
        return DATASET[time_array_<=Precision_hours,:]
    
    


def interp_track(best):
    '''
    This code takes the output of best track positions and interpolates to hourly track positions
    
    The best output from get_storms_bytime() should be used since it has multiple times to allow
    for the interpolation between points
    
    returns array of with shape of [#storms,(Times,Lats,Lons,Ints),#hours available up to 36]


    
    '''
    
    to_output = [] #the output list
    name_out  = []
    # Loop through the number of unique storms (do not want to interpolate points with multiple storms)
    for sname in np.unique(best[:,0]):
        
        #Make sure we only get the positions for one TC
        lines = best[best[:,0]==sname,:]
        if lines.shape[0]<1:
            print('No positions to interpolate around')
            continue
        elif lines.shape[0]==1:
            print('Only 1 Position, No interpolation can be done')
#            print(lines,lines[1])
#            to_output.append([datetime.datetime.strptime(lines[1],'%Y%m%d %H%M'),lines[5],lines[6],lines[7]])
            continue
        
        # New hours to interpolate to
        fhr_interp = np.arange(0,37,1)

        # First get the true times in the best track
        time_array_ = np.array([datetime.datetime.strptime(val,'%Y%m%d %H%M') for val in lines[:,1]])
        # Now get the times in hours since the first available
        hour_array_ = np.array([(datetime.datetime.strptime(val,'%Y%m%d %H%M')-time_array_[0]).total_seconds()/60./60. for val in lines[:,1]])

        # Compute new valid times 
        newtime_array_ = np.array([(time_array_[0]+datetime.timedelta(hours=float(hh))) for hh in fhr_interp])

        #The input lat/lon/intensities
        lons_in = lines[:,6]
        lats_in = lines[:,5]
        ints_in = lines[:,7]

        #The output- interpolated to hourly lat/lons
        lats_out = np.interp(fhr_interp.astype(float), hour_array_.astype(float), lats_in.astype(float),left=np.nan,right=np.nan)
        lons_out = np.interp(fhr_interp.astype(float), hour_array_.astype(float), lons_in.astype(float),left=np.nan,right=np.nan)
        ints_out = np.interp(fhr_interp.astype(float), hour_array_.astype(float), ints_in.astype(float),left=np.nan,right=np.nan)
        name_out.append(sname)
        #Output only non-nan values to the list to save
        to_output.append([newtime_array_[~np.isnan(lats_out)],lats_out[~np.isnan(lats_out)],lons_out[~np.isnan(lats_out)],ints_out[~np.isnan(lats_out)]])
    
    return np.array(to_output), name_out




def cosd(v):
    #function that does cosine in degrees
    return np.cos(np.deg2rad(v))

def sind(v):
    #function that does sine in degrees
    return np.sin(np.deg2rad(v))

def haversine(lat1,lon1,lat2,lon2):
    # haversine function from http://rosettacode.org/wiki/Haversine_formula#Julia to calculate distance between points
    # Written using numpy to allow for the calculation on grids
    # outputs units in km

    return 2.* 6372.8 * np.arcsin(np.sqrt(sind((lat2-lat1)/2.)**2. + cosd(lat1)* cosd(lat2)* sind((lon2 - lon1)/2.)**2.))


def get_distance_fromcenter(latc,lonc,LAT_M=False,LON_M=False):
    # Uses the haversine function to calculate the distance from storm center
    # This function's sole purpose is to redefine the config lat/lon grid if not already available
    
    #If provide, calculate the distance field
    if LAT_M is not None:
        return haversine(LAT_M,LON_M,latc,lonc)
    
    #If not provided, get the grid from the GLM grid 
    else:
        LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

        LAT_M, LON_M = np.meshgrid(LAT_ARRAY,LON_ARRAY)


        return haversine(LAT_M,LON_M,latc,lonc)


def get_storm_mask(DATASET,YYYY,MM,DD,HH,Daily=False):
    '''
    Wrapper script to produce one grid with the Great-Circle Distance from all TCs that are present
    at a specified time.
    
    INPUT
       DATASET : the Best-Track Dataset- use: DATASET = output_best_track_dataset('AL')
       YYYY    : the year of interest
       MM      : the month of interest
       DD      : the day of interest
       HH      : the hour of interest
       Daily   : conditional- if set to True, then the hourly input is not considered
                    and the function will output the interpolated positions at the 12 hour point
    
    Returns
       2D Array   : The Great Circle Distance to the nearest TC if present
       1D Array   : TC Vmax
       1D Array   : TC Intensity Change
       1D Array   : Lat
       1D Array   : Lon
       1D Array   : Stormid
       1D Array   : Datetime

       
    '''
    #Get the domain here (probably needs to be supplied via function)
    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    LAT_M, LON_M = np.meshgrid(LAT_ARRAY,LON_ARRAY)

    
    if Daily==True:
        #This is 12 hours plus the day of interest specified in the function call
        toi = datetime.datetime(int(YYYY),int(MM),int(DD),12,0).strftime('%Y%m%d %HZ')
    else: 
        #This is the time of interest specified in the function call
        toi = datetime.datetime(int(YYYY),int(MM),int(DD),int(HH),0).strftime('%Y%m%d %HZ')
        toim24 = (datetime.datetime(int(YYYY),int(MM),int(DD),int(HH),0)+datetime.timedelta(hours=12)).strftime('%Y%m%d %HZ')
        YYYY2 = toim24[:4]
        MM2 = toim24[4:6]
        DD2 = toim24[6:8]
        HH2 = toim24[9:11]


    #Get the storm positions that match the approximate input time
    BEST_DATA = get_storms_bytime(DATASET,YYYY,MM,DD,HH,Daily=False)
    BEST_DATA2 = get_storms_bytime(DATASET,YYYY2,MM2,DD2,HH2,Daily=False)

    print('TCs Found on %s:'%toi)
#    print(np.unique(BEST_DATA[:,0]))

    #Use the approximate best-track times to interpolate to hourly data
    TC_Tracks,NAMES = interp_track(BEST_DATA)
    TC_Tracks2,NAMES2 = interp_track(BEST_DATA2)

    #TC Tracks will have shape of [#storms,(Times,Lats,Lons,Ints),#hours available up to 36]
    
    # Create a list of distances to save each storm distance array to
    # Distances are calculated for each TC and then the minimum is applied since there
    # can be multiple TCs active at any one point
    DISTANCE_ARRAYS = []
    VMAX_ARRAYS     = []
    VCHANGE_ARRAYS  = []
    TCLAT_ARRAYS  = []
    TCLON_ARRAYS  = []
    STMID_ARRAYS  = []
    DTIME_ARRAYS  = []
 
    #Only do this if we have storms at valid time
    if len(TC_Tracks)>0:
        
        for i in range(0,np.array(TC_Tracks).shape[0]):
            name_ = (NAMES[i])
            Track = TC_Tracks[i]
            Times = np.array([val.strftime('%Y%m%d %HZ') for val in Track[0]])
            Lat  = Track[1].astype(float)[Times==toi]
            Lon  = Track[2].astype(float)[Times==toi]+360.
            Int  = Track[3].astype(float)[Times==toi]
            
            Intm24 = np.nan
            for i2 in range(0,np.array(TC_Tracks2).shape[0]):
                name_2 = NAMES2[i2]
                Track2 = TC_Tracks2[i2]
                Times2 = np.array([val.strftime('%Y%m%d %HZ') for val in Track2[0]])
                Int2  = Track2[3].astype(float)[Times2==toim24]
                if (len(Int2)>0)&(name_==name_2):
                    Intm24 = Int2
                
            
            if len(Lat)>0:   
                SS_D = get_distance_fromcenter(Lat,Lon,LAT_M=LAT_M, LON_M=LON_M) 
                # Calculate the distance from the center and save to the main list
                DISTANCE_ARRAYS.append(np.array(SS_D))
                VMAX_ARRAYS.append(Int)
                TCLAT_ARRAYS.append(Lat)
                TCLON_ARRAYS.append(Lon)
                STMID_ARRAYS.append(name_)
                DTIME_ARRAYS.append(toi)

                if np.isfinite(Intm24):                    
                    VCHANGE_ARRAYS.append(-(Int - Intm24))
                else:
                    VCHANGE_ARRAYS.append((np.nan))

            else:
                No_TCs = np.copy(LAT_M)
                No_TCs[:,:] = 99999
                DISTANCE_ARRAYS.append(No_TCs)
                VMAX_ARRAYS.append(np.nan)
                VCHANGE_ARRAYS.append((np.nan))
                TCLAT_ARRAYS.append(np.nan)
                TCLON_ARRAYS.append(np.nan)
                STMID_ARRAYS.append(np.nan)
                DTIME_ARRAYS.append(np.nan)

        #return the minimum distance array
        return np.array(DISTANCE_ARRAYS),np.array(VMAX_ARRAYS),np.array(VCHANGE_ARRAYS),np.array(TCLAT_ARRAYS),np.array(TCLON_ARRAYS),np.array(STMID_ARRAYS),np.array(DTIME_ARRAYS)
    
    else: #If we have no storms at the valid time
        No_TCs = np.copy(LAT_M)
        No_TCs[:,:] = 99999
        return No_TCs.reshape(1,LAT_M.shape[0],LAT_M.shape[1]),np.array(VMAX_ARRAYS),np.array(VCHANGE_ARRAYS),np.array(TCLAT_ARRAYS),np.array(TCLON_ARRAYS),np.array(STMID_ARRAYS),np.array(DTIME_ARRAYS)

    
    
        
def get_storm_min_mask(DATASET,YYYY,MM,DD,HH,Daily=False):
    '''
    Wrapper script to produce one grid with the Great-Circle Distance from all TCs that are present
    at a specified time.

    INPUT
       DATASET : the Best-Track Dataset- use: DATASET = output_best_track_dataset('AL')
       YYYY    : the year of interest
       MM      : the month of interest
       DD      : the day of interest
       HH      : the hour of interest
       Daily   : conditional- if set to True, then the hourly input is not considered
                    and the function will output the interpolated positions at the 12 hour point

    Returns
       2D Array   : The minimum Great Circle Distance to the nearest TC if present

    '''
    #Get the domain here (probably needs to be supplied via function)
    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    LAT_M, LON_M = np.meshgrid(LAT_ARRAY,LON_ARRAY)
 
    if Daily==True:
        #This is 12 hours plus the day of interest specified in the function call
        toi = datetime.datetime(int(YYYY),int(MM),int(DD),12,0).strftime('%Y%m%d %HZ')
    else:
        #This is the time of interest specified in the function call
        toi = datetime.datetime(int(YYYY),int(MM),int(DD),int(HH),0).strftime('%Y%m%d %HZ')

    #Get the storm positions that match the approximate input time
    BEST_DATA = get_storms_bytime(DATASET,YYYY,MM,DD,HH,Daily=False)

    print('TCs Found on %s:'%toi)
#    print(np.unique(BEST_DATA[:,0]))

    #Use the approximate best-track times to interpolate to hourly data
    TC_Tracks,NAMES = interp_track(BEST_DATA)

    #TC Tracks will have shape of [#storms,(Times,Lats,Lons,Ints),#hours available up to 36]

    # Create a list of distances to save each storm distance array to
    # Distances are calculated for each TC and then the minimum is applied since there
    # can be multiple TCs active at any one point
    DISTANCE_ARRAYS = []

    #Only do this if we have storms at valid time
    if len(TC_Tracks)>0:
        for i in range(0,np.array(TC_Tracks).shape[0]):
            name_ = (NAMES[i])
            Track = TC_Tracks[i]
            Times = np.array([val.strftime('%Y%m%d %HZ') for val in Track[0]])
            Lat  = Track[1].astype(float)[Times==toi]
            Lon  = Track[2].astype(float)[Times==toi]+360.
            Int  = Track[3].astype(float)[Times==toi]


            if len(Lat)>0:
                SS_D = get_distance_fromcenter(Lat,Lon,LAT_M=LAT_M, LON_M=LON_M)
                # Calculate the distance from the center and save to the main list
                DISTANCE_ARRAYS.append(np.array(SS_D))

            else:
                No_TCs = np.copy(LAT_M)
                No_TCs[:,:] = 99999
                DISTANCE_ARRAYS.append(No_TCs)


        #return the minimum distance array
        return np.nanmin(np.array(DISTANCE_ARRAYS),axis=0)

    else: #If we have no storms at the valid time
        No_TCs = np.copy(LAT_M)
        No_TCs[:,:] = 99999
        return No_TCs

    
    
 
def bin_byic(alldata,intchange):
    bw=5
    lowerb = np.arange(-52.5,50,bw)
    upperb = lowerb+bw
    
    freqs = []
    
    for lb,ub in zip(lowerb,upperb):
        data = alldata[(intchange>=lb)&(intchange<ub)]
        if len(data)==0:
            freqs.append(np.nan)
        else:
            freqs.append(data)
        
        
    return lowerb+bw/2.,freqs


def bin_byvmax(alldata_,vmax_):
    bw=10
    lowerb = np.arange(25,140,bw)
    upperb = lowerb+bw
    
    freqs = []
    
    for lb,ub in zip(lowerb,upperb):
        tdata = alldata_[(vmax_>=lb)&(vmax_<ub)]
        if len(tdata)==0:
            freqs.append(np.nan)
        else:
            freqs.append(tdata)
        
        
    return lowerb+bw/2.,freqs




















