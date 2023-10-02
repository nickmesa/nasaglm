import matplotlib.pyplot as plt
import numpy as np
import glob
from netCDF4 import Dataset
import os
from Functions import *


def Count_Groups(filename):
    '''
    Requires:
    from netCDF4 import Dataset

    Input: filename of GLM Climatology netCDF file
    Output: Sum of all groups in the file

    '''
    GLM_file = Dataset(filename)
    
    COUNT= GLM_file.variables['GROUP_COUNT'][:]
    
    GLM_file.close()

    Total_Flashed = np.nansum(COUNT)
    return Total_Flashed
    
def Proportion_Groups(filename):
    '''
    This code finds how many hours have lightning somewhere in the domain
    Use to help remove partial files

    Requires:
    from netCDF4 import Dataset

    Input: filename of GLM Climatology netCDF file

    Output: Number of hours with lightning

    '''
    GLM_file = Dataset(filename)

    COUNT= GLM_file.variables['GROUP_COUNT'][:]

    GLM_file.close()
    temp = []
    for i in range(0,24):
        temp.append(np.nansum(COUNT[i]))

#    Total_Flashed = np.nansum(np.nansum(COUNT,dtype=np.float32,axis=1),axis=1)
#    print(Total_Flashed)    
    temp = np.array(temp)
    return len(temp[temp>0])

if __name__ == '__main__':
    #-------------------

    file_list = sorted(glob.glob('../Output/2020/*'))

    p = 0
    #now loop thorugh all the files
    for i,file_ in enumerate(file_list):
        ngroup = Count_Groups(file_)
        nhours = Proportion_Groups(file_)

        filename = os.path.split(file_)[-1]
        #If there is no lightning in the file, then delete the file
        if ngroup==0:
            print('There is no lightning for %s. Deleting File'%filename)
            os.remove(file_)
        if nhours <23:
            print('There is only %s hours with lightning for %s. Deleting File'%(nhours,filename))
            os.remove(file_)    




