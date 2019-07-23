from datetime import datetime
from glob import glob
import numpy as np
from sid_func import *
import matplotlib.pyplot as plt


radius = 5000
file_name_end = '_50km_more.csv'

#-------------------------------------------------------------------
inpath = '/Data/sim/polona/sid/drift_full_stp1/'
outpath_def = '/Data/sim/polona/sid/deform/'
outpath = outpath_def+'plots/'
metfile = '../data/10minute_nounits.csv'

#get Lance location
start = datetime(2015,2,22)
#Lance postion (from Lance's met system)
mettime = getColumn(metfile,0)
dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
mi = np.argmin(abs(np.asarray(dtb)-start))
Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]

#Lance location in geographical coordinates
from pyproj import Proj, transform
inProj = Proj(init='epsg:4326')
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))
xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)
        
print(xlp,ylp)

#create a grid of spacing of 5km in each direction until 25km away from the ship
xbp = np.arange(xlp-radius, xlp+radius, 1000)
ybp = np.arange(ylp-radius, ylp+radius, 1000)

x_buoy,y_buoy = np.meshgrid(xbp,ybp)

print(x_buoy.shape)



#read in the drift product and calculate displacements from starting point until the first break in the SAR data (1 week)
fl = sorted(glob(inpath+'*.npz'))
for i in fl:
    #read in all the data
    print(i)
    container = np.load(i)
    u = container['upm']
    v = container['vpm'] 
    lat = container['lat1']
    lon = container['lon1']
    x,y = transform(inProj,outProj,lon,lat)  #probably the arrays first need to be flattened before the transform can be used!
    
    print(lat)
    print(lon)
    
    print(x)
    print(y)
    
    outlon,outlat = transform(outProj,inProj,x,y)
    print(outlon)
    print(outlat)
    exit()
    
    
    print(np.max(x))
    print(np.min(x))
    print(np.max(y))
    print(np.min(y))
    #print(np.max(lat))
    #print(np.min(lat))
    
    #get time difference
    date1 = i.split('_')[-2]
    date2 = i.split('_')[-1].split('.')[0]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")

    diff = (dt2-dt1).seconds + (dt2-dt1).days*24*60*60
    #print(diff)

    
    #find closest location to the individual virtual buoy
    for m in range(0,x_buoy.shape[0]):
        for n in range(0,y_buoy.shape[1]):
            print(x_buoy[m,n])
            print(y_buoy[m,n])
            mask = (x>x_buoy[m,n]-10000) & (x<x_buoy[m,n]+10000) #& (y<y_buoy[m,n]-10000) & (y<y_buoy[m,n]+10000)
            sample = u[mask]
            print(sample)
            print(sample.shape)
            exit()
        
            #calculate displacement and new location for this buoy


#save the output locations (in latlon) to be plotted on top of divergence and shear maps by sid_defrom.py
