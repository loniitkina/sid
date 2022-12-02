import os
import sys
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


#take one track and create several neighboring tracks at defined distance in cardinal/corner directions from the original track

shipfile = '../../downloads/position_leg3_nh-track.csv'	#leg3 (and transition to leg 4 until 6 June)
#shipfile = 'src/coord_trans/dshipextracts/transect_legs/position_leg3_nh-track.csv'

shipfile = '../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0.csv'
#shipfile = '../../downloads/data_master-solution_mosaic-leg2-20191214-20200224-floenavi-refstat-v1p0.csv'

#copy this file to region 'c' - central
f = open(shipfile,"r")
copy = open(shipfile.split('.csv')[0]+'_c.csv',"wt")
line = f.read()
copy.write(str(line))
f.close()
copy.close()

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
#EW: extra-wide swath mode is 400km wide
#if we want to collate several tiles/sample from different image pairs - use large radius here > 200km
radius=100000	#100km
radius=200000
xlp_nw=xlp-radius; ylp_nw=ylp+radius
xlp_w=xlp-radius; ylp_w=ylp
xlp_sw=xlp-radius; ylp_sw=ylp-radius
xlp_s=xlp; ylp_s=ylp-radius
xlp_se=xlp+radius; ylp_se=ylp-radius
xlp_e=xlp+radius; ylp_e=ylp
xlp_ne=xlp+radius; ylp_ne=ylp+radius
xlp_n=xlp; ylp_n=ylp+radius

#convert back to lon, lat
lon_nw,lat_nw = transform(outProj,inProj,xlp_nw, ylp_nw)
lon_w,lat_w = transform(outProj,inProj,xlp_w, ylp_w)
lon_sw,lat_sw = transform(outProj,inProj,xlp_sw, ylp_sw)
lon_s,lat_s = transform(outProj,inProj,xlp_s, ylp_s)
lon_se,lat_se = transform(outProj,inProj,xlp_se, ylp_se)
lon_e,lat_e = transform(outProj,inProj,xlp_e, ylp_e)
lon_ne,lat_ne = transform(outProj,inProj,xlp_ne, ylp_ne)
lon_n,lat_n = transform(outProj,inProj,xlp_n, ylp_n)

#print(lon_ne,lat_ne)

#write out files
tt = [mettime,lon_dtb,lat_dtb]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_c_'+str(int(radius/1000))+'km.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")

tt = [mettime,lon_nw,lat_nw]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_nw_'+str(int(radius/1000))+'km.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")
        
tt = [mettime,lon_w,lat_w]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_w_'+str(int(radius/1000))+'km.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")

tt = [mettime,lon_sw,lat_sw]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_sw_'+str(int(radius/1000))+'km.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")
        
tt = [mettime,lon_s,lat_s]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_s_'+str(int(radius/1000))+'km.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")
        
tt = [mettime,lon_se,lat_se]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_se_'+str(int(radius/1000))+'km.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")
        
tt = [mettime,lon_e,lat_e]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_e_'+str(int(radius/1000))+'km.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")
               
tt = [mettime,lon_ne,lat_ne]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_ne_'+str(int(radius/1000))+'km.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")

tt = [mettime,lon_n,lat_n]
table = list(zip(*tt))
outname = shipfile.split('.csv')[0]+'_n_'+str(int(radius/1000))+'km.csv'
print(outname)
with open(outname, 'wb') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")
