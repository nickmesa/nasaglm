#!/usr/bin/env python
import datetime
from scipy import stats
from netCDF4 import Dataset
import netCDF4
import glob
import os
#os.environ["OMP_NUM_THREADS"] = "2" # export OMP_NUM_THREADS=4
#os.environ["OPENBLAS_NUM_THREADS"] = "2" # export OPENBLAS_NUM_THREADS=4 import numpy as np

import multiprocessing 
#------------------------
# Import Custom Functions
#------------------------
from Functions import *
from Get_AOD_ALL import *
from Get_GLM import *
from Get_ABI import *
from Get_OCEAN import *
from Get_WWLLN import *
from Get_LIS import *

import warnings
warnings.filterwarnings("ignore")


def main(i):
    # The i is the day of the year for this fucntion to be run on
    #-------------------
    
    YYYY = '2023'
    Sat = '16'
    HOURS = np.arange(0,24,1)
    #-------------------

    LAT_ARRAY,LON_ARRAY,LAT_BINS,LON_BINS = GET_GRID()

    lat_min = CONFIG.LAT_MIN - CONFIG.LAT_RES/2.
    lat_max = CONFIG.LAT_MAX + CONFIG.LAT_RES/2.
    lon_min = CONFIG.LON_MIN - CONFIG.LON_RES/2.
    lon_max = CONFIG.LON_MIN + CONFIG.LON_RES/2.


    #Convert Day of Year to Month Day
    MMDD = (doy_to_mmdd(YYYY,i))
    print('Starting for %s'%MMDD)
    MM = MMDD[:2]
    DD = MMDD[2:]

    #---------------------------------------
    # OPEN NETCDF For Outputing the Data
    #--------------------------------------- 

    #---------------------------------------
    # Path of GLM
    GLM_LNG_PATH = CONFIG.GLM_PATH + YYYY + '/' + MMDD + '/OR_GLM-L2-LCFA_G' + Sat + '_s' + YYYY + str(i).zfill(3)

    # Output 
    OUTPUT_FILE = '/Users/nmesa/Desktop/NASA_GLM_Project/climatology_glm-main/Output/'+YYYY+'/GLM_CLIM_G'+Sat+'Gridded_'+YYYY+'_'+MMDD+'.nc'
    #---------------------------------------
    
    #Check if the file already exists (if it exists then we don't do anything
    if os.path.isfile(OUTPUT_FILE):
        return
   
    #If the file does not exist, then we will create a new file and write to it
    NC_FILE = netCDF4.Dataset(OUTPUT_FILE,mode='w') #open the file to write to 

    NC_FILE.createDimension('LON', len(LON_ARRAY))
    NC_FILE.createDimension('LAT', len(LAT_ARRAY))
    NC_FILE.createDimension('TIME', len(HOURS))

    LON_     = NC_FILE.createVariable('LON', 'f', ('LON'),fill_value=np.nan)
    LON_[:]  = LON_ARRAY
    LAT_     = NC_FILE.createVariable('LAT', 'f', ('LAT'),fill_value=np.nan)
    LAT_[:]  = LAT_ARRAY
    TIME_    = NC_FILE.createVariable('TIME', 'f', ('TIME'),fill_value=np.nan)
    TIME_[:] = HOURS
    TIME_.setncatts({'long_name': u"Hours in UTC from filename timestamp",\
                'units': u"UTC"})

    #Define the Variables Associated with the Ocean
    SST_ = NC_FILE.createVariable('SST', 'f', ('LON','LAT'),fill_value=np.nan)
    OHC_ = NC_FILE.createVariable('OHC', 'f', ('LON','LAT'),fill_value=np.nan)
    DTL_ = NC_FILE.createVariable('DTL', 'f', ('LON','LAT'),fill_value=np.nan)
    SST_.setncatts({'long_name': u"Reynolds Surface Temperature",\
                'units': u"C"})
    OHC_.setncatts({'long_name': u"Ocean Heat Content",\
                'units': u"J/C"})
    DTL_.setncatts({'long_name': u"Distance to Land",\
                'units': u"km"})

    SZA_ = NC_FILE.createVariable('SZA', 'f', ('LON','LAT'),fill_value=np.nan)
    SZA_.setncatts({'long_name': u"Satellite Zenith Angle",\
                'units': u"degrees"})


    #Define the Variables Associated with GLM
    LCOUNT_     = NC_FILE.createVariable('GROUP_COUNT', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    LCOUNT_.setncatts({'long_name': u"Count of GLM  Lightning Groups",\
                'units': u"Number"})
    LFCOUNT_     = NC_FILE.createVariable('FLASH_COUNT', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    LFCOUNT_.setncatts({'long_name': u"Count of GLM Lightning Flashes",\
                'units': u"Number"})
#    LCOUNT_FLAG_= NC_FILE.createVariable('COUNT_FLAG', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
#    LCOUNT_FLAG_.setncatts({'long_name': u"Count of Removed Lightning Groups from QC",\
#                    'units': u"Number"})
    MEAN_ENERGY_= NC_FILE.createVariable('MEAN_ENERGY', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    MEAN_ENERGY_.setncatts({'long_name': u"Mean Energy of Lightning Groups",\
                    'units': u"J"})
    STDV_ENERGY_= NC_FILE.createVariable('STD_ENERGY', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    STDV_ENERGY_.setncatts({'long_name': u"Standard Deviation of Energy in Grid Cell",\
                    'units': u"J"})
    MINM_ENERGY_= NC_FILE.createVariable('MIN_ENERGY', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    MINM_ENERGY_.setncatts({'long_name': u"Minimum Energy in Grid Cell",\
                'units': u"J"})
    MAXM_ENERGY_= NC_FILE.createVariable('MAX_ENERGY', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    MAXM_ENERGY_.setncatts({'long_name': u"Max Energy in Grid Cell",\
                'units': u"J"})
    TOTL_ENERGY_= NC_FILE.createVariable('TOTAL_ENERGY', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    TOTL_ENERGY_.setncatts({'long_name': u"Total Energy in Grid Cell",\
                'units': u"J"})
    MEAN_AREA_  = NC_FILE.createVariable('MEAN_AREA', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    MEAN_AREA_.setncatts({'long_name': u"Mean Group Area in Grid Cell",\
                'units': u"m2"})
    STDV_AREA_  = NC_FILE.createVariable('STD_AREA', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    STDV_AREA_.setncatts({'long_name': u"Area Standard Deviation in Grid Cell",\
                 'units': u"m2"})
    MINM_AREA_  = NC_FILE.createVariable('MIN_AREA', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    MINM_AREA_.setncatts({'long_name': u"Min Group Area in Grid Cell",\
                'units': u"m2"})
    MAXM_AREA_  = NC_FILE.createVariable('MAX_AREA', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    MAXM_AREA_.setncatts({'long_name': u"Max Group Area in Grid Cell",\
                'units': u"m2"})    

    #Define QC Variables Associated with GLM
    QC_LCOUNT_     = NC_FILE.createVariable('QC_GROUP_COUNT', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    QC_LCOUNT_.setncatts({'long_name': u"Count of GLM  Lightning Groups Removed by QC",\
                'units': u"Number"})
    QC_LFCOUNT_     = NC_FILE.createVariable('QC_FLASH_COUNT', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    QC_LFCOUNT_.setncatts({'long_name': u"Count of GLM Lightning Flashes Removed by QC",\
                'units': u"Number"})
    QC_MEAN_ENERGY_= NC_FILE.createVariable('QC_MEAN_ENERGY', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    QC_MEAN_ENERGY_.setncatts({'long_name': u"Mean Energy of Lightning Groups Removed by QC",\
                    'units': u"J"})
    QC_MEAN_AREA_  = NC_FILE.createVariable('QC_MEAN_AREA', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    QC_MEAN_AREA_.setncatts({'long_name': u"Mean Group Area in Grid Cell Removed by QC",\
                'units': u"m2"})

    #Now Create the Variables for the ABI Variable
    MINM_BT_  = NC_FILE.createVariable('MIN_BT', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    MINM_BT_.setncatts({'long_name': u"Min ABI IR Brightness Temperature in Grid Cell",\
                'units': u"K"})
    STDV_BT_  = NC_FILE.createVariable('STD_BT', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    STDV_BT_.setncatts({'long_name': u"Standard Deviation of ABI IR Brightness Temperature in Grid Cell",\
                'units': u"K"})
    MEAN_BT_  = NC_FILE.createVariable('MEAN_BT', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    MEAN_BT_.setncatts({'long_name': u"Mean ABI IR Brightness Temperature in Grid Cell",\
                'units': u"K"})
#    MAXM_BT_  = NC_FILE.createVariable('MAX_BT', 'f', ('TIME', 'LON','LAT'),fill_value=-np.nan)
#    MAXM_BT_.setncatts({'long_name': u"Max ABI IR Brightness Temperature in Grid Cell",\
#                'units': u"K"})
    MAX_DQF_  = NC_FILE.createVariable('MAX_DQF', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    MAX_DQF_.setncatts({'long_name': u"Max QC Flag in Grid Cell"})


    #Define the Variables Associated with WWLLN
    WLCOUNT_     = NC_FILE.createVariable('WWLLN_COUNT', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    WLCOUNT_.setncatts({'long_name': u"Count of WWLLN Lightning Groups",\
                'units': u"Number"})


    #Define the Variables Associated with LIS
    LIS_GRC_     = NC_FILE.createVariable('LIS_COUNT', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    LIS_GRC_.setncatts({'long_name': u"Count of LIS Lightning Groups",\
                'units': u"Number"})
    LIS_VW_     = NC_FILE.createVariable('LIS_SWATH', 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)
    LIS_VW_.setncatts({'long_name': u"Mask for LIS Swath",\
                'units': u"Number"})


    #Now create the variables for AOD
    AODTypes = ['','_du','_sa','_sm','_su']
    AODNames = ['Total','Dust','Salt','Smoke','Sulfate']

    # We will loop through the types and names to create the variables
    for AODT,AODN in zip(AODTypes,AODNames):
        name_ = 'AOD'+AODT
        lname_ = AODN+" Aerosal Optical Depth"

        #Define the variable in the Netcdf file
        exec("AOD_%s  = NC_FILE.createVariable(name_, 'f', ('TIME', 'LON','LAT'),fill_value=np.nan)"%AODT)
        exec("AOD_%s.setncatts({'long_name': lname_})"%AODT)





    #---------------------------------------
    # End Variable Creation in netCDF File
    #
    # Begin Output Data to NETCDF
    #--------------------------------------- 

    # Add the satellite zenith angle
    lat_sat, lon_sat, area_sat = get_abi_grid()
    satzen,satazm = read_satzen()
    SZA = bin_statistic(lat_sat[~np.isnan(satzen)],lon_sat[~np.isnan(satzen)],satzen[~np.isnan(satzen)],LAT_BINS,LON_BINS,'mean').astype(np.float32)
    SZA_[:] = SZA  #send to netcdf




    #----------------------------------------------------------------
    #
    # Processing for Daily Variables (SST/OHC) and Distance to Land
    #
    #----------------------------------------------------------------

    #Get the SST Data
    SST, OHC, DTL  = Regrid_OCEAN(YYYY,i) 
    SST_[:] = SST.astype(np.float32) #send to netcdf
    OHC_[:] = OHC.astype(np.float32) #send to netcdf
    DTL_[:] = DTL.astype(np.float32) #send to netcdf

    #---------------------------------------
    # Check Files for Missing Data
    #--------------------------------------- 

    # For WWWLN 
    # The files are daily so they are read here and if processing files then file is updated with np.nan
    try:
        #Data structure is different if the files are recent, this year will need to be updated/checked
        if int(YYYY)<2022: 
            WHOURS,WLAT,WLON = read_WWLLN_Aloc(CONFIG.WLN_PATH+YYYY+'/A'+YYYY+MMDD+'.loc')
         #Files moved
#        elif int(YYYY)==2020:
#            WHOURS,WLAT,WLON = read_WWLLN_Aloc(CONFIG.WLN_PATH+'A'+YYYY+MMDD+'.loc')
        missing_wwlln = False
    except:
        print('WWLLN: Missing MMDD in YYYY, Filling with NaNs')
        for HI,HOUR in enumerate(HOURS):
            temp = np.zeros((LON_ARRAY.size,LAT_ARRAY.size))
            temp[:,:] = np.nan
            missing_wwlln = True
            WLCOUNT_[HI,:,:]       = temp


    # Create the file list for LIS data files
    try: #This will not work for the first or last day of the year 
        MMDD_y = (doy_to_mmdd(YYYY,i-1))
        #Get the list of files from yesterday and grab the latest 2
        LIS_Flist_yesterday = sorted(glob.glob(CONFIG.LIS_PATH+YYYY+'/ISS_LIS_SC_V1.0_'+YYYY+MMDD_y+'*'))[-2:]

        MMDD_t = (doy_to_mmdd(YYYY,i+1))
        #Get the list of files from tomorrow and grab the first 2
        LIS_Flist_tomorrow = sorted(glob.glob(CONFIG.LIS_PATH+YYYY+'/ISS_LIS_SC_V1.0_'+YYYY+MMDD_t+'*'))[:2]

        #Combine the list of files
        LIS_Flist = LIS_Flist_yesterday+sorted(glob.glob(CONFIG.LIS_PATH+YYYY+'/ISS_LIS_SC_V1.0_'+YYYY+MMDD+'*')) +LIS_Flist_tomorrow

    except:
        LIS_Flist = sorted(glob.glob(CONFIG.LIS_PATH+YYYY+'/ISS_LIS_SC_V1.0_'+YYYY+MMDD+'*'))



         
    #Now lets loop through the hours of the day
    start = 0
    for HI,HOUR in enumerate(HOURS):
        #----------------------------------------------------------------
        #
        # Processing for ABI Data
        #
        #----------------------------------------------------------------


        #Grab the ABI data
        MINM_BT, MEAN_BT, STDV_BT, MAXM_BT, MAX_DQF = GET_ABI(YYYY,str(i).zfill(3),str(HOUR).zfill(2),Sat)
        MINM_BT_[HI,:,:] = MINM_BT.astype(np.float32)
        MEAN_BT_[HI,:,:] = MEAN_BT.astype(np.float32)
        STDV_BT_[HI,:,:] = STDV_BT.astype(np.float32)
        MAX_DQF_[HI,:,:] = MAX_DQF.astype(np.float32)

        #----------------------------------------------------------------
        #
        # Processing for GLM Data
        #
        #----------------------------------------------------------------

        #Get all the files over this hour (not always the same)
        LIGHTNING_LIST = glob.glob(GLM_LNG_PATH+str(HOUR).zfill(2)+'*')
        print('Processing %i GLM files for %iZ'%(len(LIGHTNING_LIST),HOUR))

        #read in the lightning group data for single file
        #lats_l, lons_l, ener_l, area_l, flag_l, lats_lF, lons_lF, ener_lF, area_lF, flag_lF = read_glm_filelist(LIGHTNING_LIST, apply_qc=True,output_qc=False)
        lats_l, lons_l, ener_l, area_l, flag_l, flid_l, lats_lf, lons_lf, ener_lf, area_lf, flag_lf = read_glm_filelist(LIGHTNING_LIST, apply_qc=True,output_qc=False)
        print(ener_l.size)
        #IF there is no data (or QC removed it, then convert to zeros (this is the center of the bins
        if int(ener_l.size)<1:
            LCOUNT_[HI,:,:]       = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            LFCOUNT_[HI,:,:]      = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            MEAN_ENERGY_[HI,:,:]  = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            STDV_ENERGY_[HI,:,:]  = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            MINM_ENERGY_[HI,:,:]  = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            MAXM_ENERGY_[HI,:,:]  = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            MEAN_AREA_[HI,:,:]    = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            STDV_AREA_[HI,:,:]    = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            MINM_AREA_[HI,:,:]    = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            MAXM_AREA_[HI,:,:]    = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            QC_LCOUNT_[HI,:,:]       = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            QC_LFCOUNT_[HI,:,:]      = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            QC_MEAN_ENERGY_[HI,:,:]  = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            QC_MEAN_AREA_[HI,:,:]    = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)



        else: #Here we have data so lets bin everything
            LCOUNT__       = bin_statistic(lats_l,lons_l,ener_l,LAT_BINS,LON_BINS,'count')
            LFCOUNT__      = bin_statistic(lats_lf,lons_lf,ener_lf,LAT_BINS,LON_BINS,'count')
            MEAN_ENERGY__  = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'mean').astype(np.float32)
            STDV_ENERGY__  = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'std').astype(np.float32)
            MINM_ENERGY__  = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'min').astype(np.float32)
            MAXM_ENERGY__  = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'max').astype(np.float32)
            TOTL_ENERGY__  = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'sum').astype(np.float32)
            MEAN_AREA__    = bin_statistic(lats_l[~np.isnan(area_l)],lons_l[~np.isnan(area_l)],area_l[~np.isnan(area_l)],LAT_BINS,LON_BINS,'mean').astype(np.float32)
            STDV_AREA__    = bin_statistic(lats_l[~np.isnan(area_l)],lons_l[~np.isnan(area_l)],area_l[~np.isnan(area_l)],LAT_BINS,LON_BINS,'std').astype(np.float32)
            MINM_AREA__    = bin_statistic(lats_l[~np.isnan(area_l)],lons_l[~np.isnan(area_l)],area_l[~np.isnan(area_l)],LAT_BINS,LON_BINS,'min').astype(np.float32)
            MAXM_AREA__    = bin_statistic(lats_l[~np.isnan(area_l)],lons_l[~np.isnan(area_l)],area_l[~np.isnan(area_l)],LAT_BINS,LON_BINS,'max').astype(np.float32)

#            QC_LCOUNT_[HI,:,:]       = bin_statistic(lats_l,lons_l,ener_l,LAT_BINS,LON_BINS,'count')
#            QC_LFCOUNT_[HI,:,:]      = bin_statistic(lats_lf,lons_lf,ener_lf,LAT_BINS,LON_BINS,'count')
#            QC_MEAN_ENERGY_[HI,:,:]  = bin_statistic(lats_l[~np.isnan(ener_l)],lons_l[~np.isnan(ener_l)],ener_l[~np.isnan(ener_l)],LAT_BINS,LON_BINS,'mean').astype(np.float32)
#            QC_MEAN_AREA_[HI,:,:]    = bin_statistic(lats_l[~np.isnan(area_l)],lons_l[~np.isnan(area_l)],area_l[~np.isnan(area_l)],LAT_BINS,LON_BINS,'mean').astype(np.float32)
            QC_LCOUNT__       = np.copy(LCOUNT__)
            QC_LFCOUNT__      = np.copy(LFCOUNT__)
            QC_MEAN_ENERGY__  = np.copy(MEAN_ENERGY__)
            QC_MEAN_AREA__    = np.copy(MEAN_AREA__)

            print('Applying QC')

            # Here we will do the Additional QC of the variables on the gridded products
            # Threshholds included in this QC
            Log_E = np.log10(MEAN_ENERGY__)
            Log_A = np.log10(MEAN_AREA__)
            Var_A = np.log10(STDV_AREA__**2)
            Var_E = np.log10(STDV_ENERGY__**2)

            THRESH_AREA   = Log_E * .54 + 15.75
            THRESH_ENERGY = Log_A * 2.3 - 32.0
            
            #Condition for Bad Groups
            CONDITION_BAD = ((((SZA>=50.)&(Log_A<THRESH_AREA))|(Log_E>THRESH_ENERGY)|(Var_E<-30.)|(Var_A<=13.5))|((MINM_BT>273.)&(Log_A>0.)))
            CONDITION_GOOD = ~((((SZA>=50.)&(Log_A<THRESH_AREA))|(Log_E>THRESH_ENERGY)|(Var_E<-30.)|(Var_A<=13.5))|((MINM_BT>273.)&(Log_A>0.)))
          

            #For QC variable, we need to set good values to nan
            QC_LCOUNT__[CONDITION_GOOD]       = np.nan 
            QC_LFCOUNT__[CONDITION_GOOD]      = np.nan
            QC_MEAN_ENERGY__[CONDITION_GOOD]  = np.nan
            QC_MEAN_AREA__[CONDITION_GOOD]    = np.nan

            # For regular variables set bad values to nan
            LCOUNT__[CONDITION_BAD]       = np.nan
            LFCOUNT__[CONDITION_BAD]      = np.nan
            MEAN_ENERGY__[CONDITION_BAD]  = np.nan
            STDV_ENERGY__[CONDITION_BAD]  = np.nan
            MINM_ENERGY__[CONDITION_BAD]  = np.nan
            MAXM_ENERGY__[CONDITION_BAD]  = np.nan
            TOTL_ENERGY__[CONDITION_BAD]  = np.nan
            MEAN_AREA__[CONDITION_BAD]    = np.nan
            STDV_AREA__[CONDITION_BAD]    = np.nan
            MINM_AREA__[CONDITION_BAD]    = np.nan
            MAXM_AREA__[CONDITION_BAD]    = np.nan

            # Now put into netcdf
            QC_LCOUNT_[HI,:,:]       = QC_LCOUNT__.astype(np.float32)
            QC_LFCOUNT_[HI,:,:]      = QC_LFCOUNT__.astype(np.float32)
            QC_MEAN_ENERGY_[HI,:,:]  = QC_MEAN_ENERGY__.astype(np.float32)
            QC_MEAN_AREA_[HI,:,:]    = QC_MEAN_AREA__.astype(np.float32)

            # For regular variables set bad values to nan
            LCOUNT_[HI,:,:]        = LCOUNT__.astype(np.float32)
            LFCOUNT_[HI,:,:]       = LFCOUNT__.astype(np.float32)
            MEAN_ENERGY_[HI,:,:]   = MEAN_ENERGY__.astype(np.float32)
            STDV_ENERGY_[HI,:,:]   = STDV_ENERGY__.astype(np.float32)
            MINM_ENERGY_[HI,:,:]   = MINM_ENERGY__.astype(np.float32)
            MAXM_ENERGY_[HI,:,:]   = MAXM_ENERGY__.astype(np.float32)
            TOTL_ENERGY_[HI,:,:]   = TOTL_ENERGY__.astype(np.float32)
            MEAN_AREA_[HI,:,:]     = MEAN_AREA__.astype(np.float32)
            STDV_AREA_[HI,:,:]     = STDV_AREA__.astype(np.float32)
            MINM_AREA_[HI,:,:]     = MINM_AREA__.astype(np.float32)
            MAXM_AREA_[HI,:,:]     = MAXM_AREA__.astype(np.float32)



# No longer add the count of removed data from QC
#        if lats_lF.size>0: #If there are groups removed by the QC, count them
#            LCOUNT_FLAG_[HI,:,:]  = bin_statistic(lats_lF,lons_lF,ener_lF,LAT_BINS,LON_BINS,'count').astype(np.float32)
#        else:
#            LCOUNT_FLAG_[HI,:,:]  = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)


        #----------------------------------------------------------------
        #
        # Processing for AOD Data
        #
        #----------------------------------------------------------------

        #Grab the AOD Data, 
        for AODT,AODN in zip(AODTypes,AODNames):
                exec("AOD_%s[HI,:,:]  = Get_Interp_AOD(YYYY,i,HOUR,AODT).astype(np.float32)"%AODT)

        #----------------------------------------------------------------
        #
        # Processing for WWLLN Data
        #
        #----------------------------------------------------------------
        if missing_wwlln==False: #If the data is missing then skip
            #Get the hourly lightning totals
            l_lat=WLAT[WHOURS==int(HOUR)]
            l_lon=WLON[WHOURS==int(HOUR)]

            #IF there is no data (or QC removed it, then convert to zeros (this is the center of the bins
            if l_lat.size<1:
                WLCOUNT_[HI,:,:]       = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            else:
                WLCOUNT_[HI,:,:]       = bin_statistic(l_lat,l_lon,l_lon,LAT_BINS,LON_BINS,'count').astype(np.float32)




        #----------------------------------------------------------------
        #
        # Processing for LIS Data
        #
        #----------------------------------------------------------------
        try:
            VMASK,LIS_GRP    = LIS_path(LIS_Flist,int(YYYY),int(MM),int(DD),int(HOUR))
            LIS_GRC_[HI,:,:] = LIS_GRP.astype(np.float32)
            LIS_VW_[HI,:,:]  = VMASK.astype(np.float32)
        except:
            missing = np.zeros((LON_ARRAY.size,LAT_ARRAY.size)).astype(np.float32)
            missing[:,:]     = np.nan
            LIS_GRC_[HI,:,:] = missing
            LIS_VW_[HI,:,:]  = missing

    print('Completed  %s'%MMDD)
    NC_FILE.close()




if __name__ == '__main__':
    DOYs = np.arange(1,366,1)

    p = multiprocessing.Pool(5)
    p.map(main, DOYs)
    p.close()


#DOYs = np.arange(1,366,1)

#main(DOYs)



