import netCDF4 as nc4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

file = '/Users/nasarosesnas/glmtc/GLM/L2/2017/0828/OR_GLM-L2-LCFA_G16_s20172400002000_e20172400002200_c20172400002226.nc'
data = Dataset(file, mode = 'r')

flash_lat = data.variables['flash_lat'][:]
flash_lon = data.variables['flash_lon'][:]
flash_energy = data.variables['flash_energy'][:]

#for i in flash_energy:

fig = plt.figure()
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
ax.set_extent([-70, -80, 20, 30])
ax.scatter(flash_lon, flash_lat)
#plt.plot(flash_lon, flash_lat)
#map = ax.pcolormesh([flash_lon[:], flash_lat[:]],)# flash_energy[:],)
plt.show()



print('stop')
