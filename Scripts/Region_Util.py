import numpy as np
import matplotlib.patches as patches



class CFG_REG:
    #Read in the namelist file and create a CONFIG class to be called by individual functions
    config_file = np.genfromtxt('../Config/Regions', dtype=str)

    for line in config_file: #Go through the file and execute the commands
        exec(''.join(line)) #What this does is essentially creat a global variable that is found when importing this code


def get_RX(var,lat,lon,REGNAME,ret_latlon=True):
    lon_min = eval('CFG_REG.%s_lon_min'%REGNAME)
    lat_min = eval('CFG_REG.%s_lat_min'%REGNAME)
    lon_max = eval('CFG_REG.%s_lon_max'%REGNAME)
    lat_max = eval('CFG_REG.%s_lat_max'%REGNAME)

    lon_con = np.squeeze([(lon>=lon_min)&(lon<lon_max)])
    lat_con = np.squeeze([(lat>=lat_min)&(lat<lat_max)])

    if ret_latlon==True:
        return lat[lat_con], lon[lon_con],var.T[lat_con][:,lon_con]
    else:
        return var.T[lat_con][:,lon_con]

def add_reg_boxes(ax,REGNAMES,PROJ):

    for REGNAME in REGNAMES:
        R_name  = eval('CFG_REG.%s_NAME'%REGNAME)
        lon_min = eval('CFG_REG.%s_lon_min'%REGNAME)
        lat_min = eval('CFG_REG.%s_lat_min'%REGNAME)
        lon_max = eval('CFG_REG.%s_lon_max'%REGNAME)
        lat_max = eval('CFG_REG.%s_lat_max'%REGNAME)

        #Add region patches
        rect1 = patches.Rectangle((lon_min, lat_min), np.abs(lon_max-lon_min), np.abs(lat_max-lat_min),transform=PROJ, linewidth=1, edgecolor='k', facecolor='none')
        ax.add_patch(rect1)


        #Add region labels
        ax.text(lon_min+np.abs(lon_max-lon_min)/2, lat_min+np.abs(lat_max-lat_min)/2, R_name,
            horizontalalignment='center',fontsize=8,
            verticalalignment='center',transform=PROJ)





