import numpy as np

def get_abi_grid():
    '''
    This function will read in a binary file to get the abi grid

    Outputs the gridded lat, lon and area of each grid

    Note that the transition of G16 from the checkout location will be misrepresented

    '''

    filename="/Users/nmesa/Desktop/NASA_GLM_Project/climatology_glm-main/Data/G16_east_grid.bin"

    f = open(filename,'rb') #Open the file

    #get the dimensions in the x and y directions
    ydim = np.fromfile(f,dtype=np.int32,count=1)[0]
    xdim = np.fromfile(f,dtype=np.int32,count=1)[0]

    count = ydim*xdim #Total number of grid points
    shape = (ydim,xdim) #shape of the arrays (it is a square so order doesn't matter

    #Get the lat lon and area from the file
    lat  = np.fromfile(f,dtype=np.float32,count=count).reshape(shape)
    lon  = np.fromfile(f,dtype=np.float32,count=count).reshape(shape)
    area = np.fromfile(f,dtype=np.float32,count=count).reshape(shape)

    f.close() #close the file

    lon = lon+360.

    return lat, lon, area


    
def read_satzen():

    '''
    This function will read in a binary file to get the abi satellite zenith angle

    Outputs the gridded zenith angle and azimuth angle

    Note that the transition of G16 from the checkout location will be misrepresented

    Use in coordination with get_abi_grid
    '''

    filename="/Users/nmesa/Desktop/NASA_GLM_Project/climatology_glm-main/Data/G16_east_SATZEN.bin"

    f = open(filename,'rb') #Open the file

    #get the dimensions in the x and y directions
    ydim = np.fromfile(f,dtype=np.int32,count=1)[0]
    xdim = np.fromfile(f,dtype=np.int32,count=1)[0]
    
    #read the solar zenith angle and azimuth angle
    satzen = np.fromfile(f,dtype=np.float32,count=ydim*xdim).reshape((ydim,xdim))
    satazm = np.fromfile(f,dtype=np.float32,count=ydim*xdim).reshape((ydim,xdim))
    f.close()
    return satzen,satazm


if __name__ == '__main__':

    lat, lon, area = get_abi_grid()

    satzen,satazm = read_satzen()

    print(lon.shape,np.nanmin(lon),np.nanmax(lon))
    print(lat.shape,np.nanmin(lat),np.nanmax(lat))
    print(area.shape,np.nanmin(area),np.nanmax(area))

