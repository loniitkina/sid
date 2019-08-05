from datetime import datetime
import numpy as np
import pyresample as pr
from sid_func import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

inpath = '/Data/sim/polona/sid/cs/'
inpath_b = '/Data/sim/polona/sid/buoys/'
outpath = '/Data/sim/polona/sid/deform/plots/'
metfile = '../data/10minute_nounits.csv'
radius = 50000


#Take corners of the SAR region area analyzed by sid_deform.py and check if based on OSI-SAF any of the corners appear to be FYI
#Take coordinates on 18. January 2015
# Create the figure and basemap object
maptime = datetime(2015,1,18,18,7,0)

#Lance postion (from Lance's met system)
mettime = getColumn(metfile,0)
dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
mi = np.argmin(abs(np.asarray(dtb)-maptime))
Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]

print(Lance_lat,Lance_lon)

#Lance centered projection
from pyproj import Proj, transform
inProj = Proj(init='epsg:4326')
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))
xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)

#corners
crnx = np.empty((4),dtype=float)
crny = np.empty((4),dtype=float)
crnx[0],crny[0] = (xlp-radius,ylp-radius)
crnx[1],crny[1] = (xlp+radius,ylp-radius)
crnx[2],crny[2] = (xlp+radius,ylp+radius)
crnx[3],crny[3] = (xlp-radius,ylp+radius)
print(crnx,crny)

for i in range(0,crnx.size):
    #transform back to latlon
    lon,lat = transform(outProj,inProj,crnx[i], crny[i])
    print(lon,lat)
    exit()



#If this works do same aslo for the CSM image corners
