import numpy as np
import datetime
from scipy import stats


class CONFIG:
    #Read in the namelist file and create a CONFIG class to be called by individual functions
    config_file = np.genfromtxt('../Config/Namelist', dtype=str)

    for line in config_file: #Go through the file and execute the commands
        exec(''.join(line)) #What this does is essentially creat a global variable that is found when importing this code


def GET_GRID():
    #Get the grid defined in the Namelist using the CONFIG class

    #Create the lat,lon arrays 
    LAT_ARRAY =np.arange(CONFIG.LAT_MIN,CONFIG.LAT_MAX+CONFIG.LAT_RES,CONFIG.LAT_RES)
    LON_ARRAY =np.arange(CONFIG.LON_MIN,CONFIG.LON_MAX+CONFIG.LON_RES,CONFIG.LON_RES)

    #Create bins which extend beyond the lat min/max
    LAT_BINS = np.append(LAT_ARRAY,(LAT_ARRAY[-1]+CONFIG.LAT_RES))
    LAT_BINS = LAT_BINS-CONFIG.LAT_RES/2.

    #Create bins which extend beyond the lon min/max
    LON_BINS = np.append(LON_ARRAY,(LON_ARRAY[-1]+CONFIG.LON_RES))
    LON_BINS = LON_BINS-CONFIG.LON_RES/2.


    return LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS


def doy_to_mmdd(YYYY,JJJ):
    # Day of year to Month-Day converter
    #YYYY is the year as a string or int
    #JJJ is the day of the year as string or int

    #Check if the
    if (len(str(JJJ))<3):
        JJJ=str(JJJ).zfill(3)

    #Return string in format of 'MMDD'
    return datetime.datetime.strptime(str(YYYY)+str(JJJ),'%Y%j').strftime('%m%d')

def flatten(t):
    #removes nested lists
    return [item for sublist in t for item in sublist]


def bin_statistic(lat,lon,val,lat_b,lon_b,STAT_TYPE):
    #Need
    #from scipy import stats
    #This function will bin ungridded data to the desired grid and calculate a given statistic within the bins
    #STAT_TYPE = 'count' 'mean' 'std' 'min' 'max'

    statistic, xedges, yedges, binnumber = stats.binned_statistic_2d(
        lon, lat, values=val, statistic=STAT_TYPE, bins=[lon_b, lat_b])

    return statistic


def quick_display(var,**kwargs):
    #Define a figure
    fig = plt.figure(figsize=(18,10))
    ax = plt.subplot(111) #create the subplot

    #make a colormap with white set under
    cmap_ = plt.cm.rainbow
    cmap_.set_under('white')
    
    #plot some data
    im = ax.pcolor(var,cmap=cmap_,**kwargs)

    plt.colorbar(im) #add colorbar

    plt.show() #show it




def get_nearest_min(xlats,xlons,gridlat,gridlon,gridvar):
    '''
    Input:
    xlats:    list or float of center latitude point
    ylats:    list or float of center longitude point
    gridlat: 1D array of latitudes corresponding to second gridvar dimension
    gradlon: 1D array of longitudes corresponding to first gridvar dimension 
    gridvar: 2D array of [lon,lat] that will be resampled and returned
    
    Output:
    List of maximum values of gridvar at within the number of buffer points surrounding (xlat,xlon)
    
    '''    
    # How many points to resample in the x/y direction
    # The buffer will be split, so 5 will give the nearest point and 2 points above and below the value
    buffer_ = 2
    
    if not ((type(xlats)==list) | (type(xlats)==np.ndarray)):
        xlats = [xlats]
        xlons = [xlons]
        
    output = []

    for xlat,xlon in zip(xlats,xlons):

        xis = np.argpartition((gridlat-xlat)**2,buffer_)
        yis = np.argpartition((gridlon-xlon)**2,buffer_)

        xi_mesh,yi_mesh = np.meshgrid(sorted(xis[:buffer_]),sorted(yis[:buffer_]))

        output.append(np.nanmin(gridvar[yi_mesh,xi_mesh]))

    return output



def get_nearest_max(xlats,xlons,gridlat,gridlon,gridvar):
    '''
    Input:
    xlats:    list or float of center latitude point
    ylats:    list or float of center longitude point
    gridlat: 1D array of latitudes corresponding to second gridvar dimension
    gradlon: 1D array of longitudes corresponding to first gridvar dimension 
    gridvar: 2D array of [lon,lat] that will be resampled and returned
    
    Output:
    List of maximum values of gridvar at within the number of buffer points surrounding (xlat,xlon)
    
    '''    
    # How many points to resample in the x/y direction
    # The buffer will be split, so 5 will give the nearest point and 2 points above and below the value
    buffer_ = 2
    
    if not ((type(xlats)==list) | (type(xlats)==np.ndarray)):
        xlats = [xlats]
        xlons = [xlons]
        
    output = []

    for xlat,xlon in zip(xlats,xlons):

        xis = np.argpartition((gridlat-xlat)**2,buffer_)
        yis = np.argpartition((gridlon-xlon)**2,buffer_)

        if len(xis.shape)>1:
            xi_mesh = xis
            yi_mesh = yis
        else:
            xi_mesh,yi_mesh = np.meshgrid(sorted(xis[:buffer_]),sorted(yis[:buffer_]))

        output.append(np.nanmax(gridvar[yi_mesh,xi_mesh]))

    return output

def get_nearest_mean(xlats,xlons,gridlat,gridlon,gridvar):
    '''
    Input:
    xlats:    list or float of center latitude point
    ylats:    list or float of center longitude point
    gridlat: 1D array of latitudes corresponding to second gridvar dimension
    gradlon: 1D array of longitudes corresponding to first gridvar dimension 
    gridvar: 2D array of [lon,lat] that will be resampled and returned
    
    Output:
    List of mean values of gridvar at within the number of buffer points surrounding (xlat,xlon)
    
    '''    
    # How many points to resample in the x/y direction
    # The buffer will be split, so 5 will give the nearest point and 2 points above and below the value
    buffer_ = 2
    
    if not ((type(xlats)==list) | (type(xlats)==np.ndarray)):
        xlats = [xlats]
        xlons = [xlons]
        
    output = []

    for xlat,xlon in zip(xlats,xlons):
        if len(gridlat.shape)>1:
            dist =np.sqrt(((gridlon-xlon)**2)**2+((gridlat-xlat)**2)**2)
            pts = np.unravel_index(np.nanargmin(dist),dist.shape)
            
            xis_arr = [pts[0]-1,pts[0],pts[0]+1]
            yis_arr = [pts[1]-1,pts[1],pts[1]+1]
            
            xi_mesh,yi_mesh = np.meshgrid(xis_arr,yis_arr)

        
        else:
            xis = np.argpartition((gridlat-xlat)**2,buffer_)
            yis = np.argpartition((gridlon-xlon)**2,buffer_)

            xi_mesh,yi_mesh = np.meshgrid(sorted(xis[:buffer_]),sorted(yis[:buffer_]))


        output.append(np.nanmean(gridvar[yi_mesh,xi_mesh]))

    return output
