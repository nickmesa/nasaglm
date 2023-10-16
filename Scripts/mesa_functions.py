import pandas as pd
import json
import datetime as dt
import matplotlib.pyplot as plt
from tropycal import tracks
import numpy as np

def get_ships(storm_name, year, month, day, hour):
    basin = tracks.TrackDataset(basin='north_atlantic',include_btk=False)
    storm = basin.get_storm((storm_name,year))
    ships = storm.get_ships(dt.datetime(year, month, day, hour))
    shearmag = ships.shear_kt[0]
    sheardir = ships.shear_dir[0]
    return shearmag, sheardir


def get_hurdat(storm_name, year, month, day, hour, minute):
    basin = tracks.TrackDataset(basin='north_atlantic',include_btk=False)
    storm = basin.get_storm((storm_name,year))
    time = storm.time
    toi = np.where(time == dt.datetime(year,month,day,hour,minute))
    lat = storm.lat[toi[0]][0]
    lon = storm.lon[toi[0]][0]
    vmax = storm.vmax[toi[0]][0]
    mslp = storm.mslp[toi[0]][0]
    return lat, lon, vmax, mslp


print(get_hurdat('ida',2021,8,26,12,0))
print(get_ships('ida',2021,8,27,18))