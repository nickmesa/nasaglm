import netCDF4 as nc4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

file = '/Users/nmesa/Desktop/NASA_GLM_Project/GLM-L2-LCFA-2019-244-16-OR_GLM-L2-LCFA_G16_s20192441600000_e20192441600200_c20192441600227.nc'
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
