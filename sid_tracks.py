from pyproj import Proj, transform
import numpy as np
from datetime import datetime

import csv
def getColumn(filename, column, delimiter=',', skipinitialspace=False, skipheader=1):
    results = csv.reader(open(filename),delimiter=delimiter,skipinitialspace=skipinitialspace)
    while skipheader>0:
        next(results, None)
        skipheader=skipheader-1
    return [result[column] for result in results]


#take one track and create several neighboring tracks at dfined distance in cardinal/corner directions from the original track

shipfile = '../../downloads/position_leg3_nh-track.csv'	#leg3 (and transition to leg 4 until 6 June)
shipfile = 'src/coord_trans/dshipextracts/transect_legs/position_leg3_nh-track.csv'

#ship postion
mettime = getColumn(shipfile,0)
dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
lon_dtb = np.asarray(getColumn(shipfile,1),dtype=float)
lat_dtb = np.asarray(getColumn(shipfile,2),dtype=float)

#define projections
inProj = Proj(init='epsg:4326')
#Use same projection as in pattern matching part of the sea ice drift algorithm
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +ellps=WGS84 +units=m' % (90, 10))
xlp,ylp = transform(inProj,outProj,lon_dtb, lat_dtb)

#get those points
radius=100000	#100km
xlp_nw=xlp-radius; ylp_nw=ylp+radius
xlp_sw=xlp-radius; ylp_sw=ylp-radius
xlp_se=xlp+radius; ylp_se=ylp-radius
xlp_ne=xlp+radius; ylp_ne=ylp+radius

#convert back to lon, lat
lon_nw,lat_nw = transform(outProj,inProj,xlp_nw, ylp_nw)
lon_sw,lat_sw = transform(outProj,inProj,xlp_sw, ylp_sw)
lon_se,lat_se = transform(outProj,inProj,xlp_se, ylp_se)
lon_ne,lat_ne = transform(outProj,inProj,xlp_ne, ylp_ne)

print(lon_ne,lat_ne)

#write out files
tt = [mettime,lon_nw,lat_nw]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_nw.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")

tt = [mettime,lon_sw,lat_sw]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_sw.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")
        
tt = [mettime,lon_se,lat_se]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_se.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")
        
tt = [mettime,lon_ne,lat_ne]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_ne.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")


