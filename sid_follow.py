from datetime import datetime
from glob import glob
import numpy as np
import pyresample as pr
from sid_func import *
import matplotlib.pyplot as plt


radius = 30000
file_name_end = '_50km_more.csv'

#-------------------------------------------------------------------
inpath = '/Data/sim/polona/sid/drift_full_stp1/'
outpath_data = '/Data/sim/polona/sid/deform/'
outpath = outpath_data+'plots/'
metfile = '../data/10minute_nounits.csv'
reg = 'leg1'
proj = reg

#get Lance location
start = datetime(2015,1,21,6,54,4)                 #this coincides with the 1st SAR image date/time
#Lance postion (from Lance's met system)
mettime = getColumn(metfile,0)
dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
mi = np.argmin(abs(np.asarray(dtb)-start))
Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]

#Lance start location in geographical coordinates
from pyproj import Proj, transform
inProj = Proj(init='epsg:4326')
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=km' % (90, 10))
xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)

#create a grid of spacing of 5km in each direction until half of the radius from the ship
#Lance will be in the center
xbf = np.arange(xlp-(radius/2), xlp+(radius/2)+5000, 5000)
ybf = np.arange(ylp-(radius/2), ylp+(radius/2)+5000, 5000)
x_buoy,y_buoy = np.meshgrid(xbf,ybf)

#read in the drift product and calculate displacements from starting point until the first break in the SAR data (1 week)
fl = sorted(glob(inpath+'*.npz'))

#empty arrays to store bouy drifts
x_path = np.zeros((7,x_buoy.shape[0],x_buoy.shape[1]))
y_path = np.zeros((7,x_buoy.shape[0],x_buoy.shape[1]))

#start with initial coordinates
x_path[0,:,:] = x_buoy
y_path[0,:,:] = y_buoy

date = [start]

for i in range(0,6):        #just first 6 day-pairs, then 1 week gap)
    #read in all the data
    print(fl[i])
    container = np.load(fl[i])
    lat1 = container['lat1']     #these arrays are as big as both scenes together (nan-where no overlap), ~5000x5000 of 40m-pixels = 200x200km
    lon1 = container['lon1']
    lat2 = container['lat2']
    lon2 = container['lon2']
    
    x1,y1 = transform(inProj,outProj,lon1,lat1)
    x2,y2 = transform(inProj,outProj,lon2,lat2)
    
    date2 = fl[i].split('_')[-1].split('.')[0]
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    date.append(dt2)
    
    #find closest location to the individual virtual buoy
    for m in range(0,x_buoy.shape[0]):
        for n in range(0,y_buoy.shape[1]):
             
            mask = (x1>x_buoy[m,n]-60) & (x1<x_buoy[m,n]+60) & (y1>y_buoy[m,n]-60) & (y1<y_buoy[m,n]+60)    #40m resolution data, give some margin
            #read the displacement from the end location
            #update VB location for the next step in time
            x_buoy[m,n] = np.mean(np.ma.masked_invalid(x2[mask]))
            y_buoy[m,n] = np.mean(np.ma.masked_invalid(y2[mask]))
            #step 0 is initial condition, step 1 is reached after one day etc.
            x_path[i+1,m,n] = x_buoy[m,n]
            y_path[i+1,m,n] = y_buoy[m,n]
            
print(x_path[:,0,0])
#save the output locations (in latlon) to be plotted on top of divergence and shear maps by sid_defrom.py
lon_path,lat_path = transform(outProj,inProj,x_path,y_path)

#dump data into numpy file
out_file = outpath_data+'VB.npz'
np.savez(out_file,x_path=x_path,y_path=y_path,lon_path=lon_path, lat_path=lat_path, date=date)

#make some plots
out_file = outpath_data+'VB.npz'
container = np.load(out_file)
print(container.files)
lon_path = container['lon_path']
lat_path = container['lat_path']
date = container['date']

radius_proj=55000
from pyresample.geometry import AreaDefinition
#Using a projection dictionary
area_id = 'around Lance'
description = 'North Pole LAEA Europe'
proj_id = 'lance'
proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'a':6378137, 'b':6356752.3142, 'units':'m'}
width = radius_proj*2/100 #100 m spacing
height = radius_proj*2/100 #100 m spacing
area_extent = (xlp-radius_proj,ylp-radius_proj,xlp+radius_proj,ylp+radius_proj)
area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)

print(area_def)
m = pr.plot.area_def2basemap(area_def)

#scale
#m.drawmapscale(Lance_lon, Lance_lat-.3, Lance_lon+8, Lance_lat-.2, 50, units='km', barstyle='fancy',fontsize=14)
m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])


color=iter(plt.cm.rainbow(np.linspace(0,1,lon_path.shape[0])))
   
for i in range(0,lon_path.shape[0]):
    fig1    = plt.figure(figsize=(10,10))
    ax      = fig1.add_subplot(111)

    x,y = m(lon_path[i,:,:],lat_path[i,:,:])

    #nods
    cl = next(color)
    ax.plot(x,y,'o',linewidth=2,color=cl)
    
    #triangulate betwen the points
    x = x.flatten()
    y = y.flatten()
    pts = np.zeros((len(x),2))
    pts[:,0]=x; pts[:,1]=y
    if i==0:
        from scipy.spatial import Delaunay
        tri = Delaunay(pts)
    ax.triplot(pts[:,0], pts[:,1], tri.simplices.copy())
    #print(tri.simplices.copy())

    #Lance initial position
    xl, yl = m(Lance_lon, Lance_lat)
    ax.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    
    #Lance moving with the VBs
    mi = np.argmin(abs(np.asarray(dtb)-date[i]))
    Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
    Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]
    xl, yl = m(Lance_lon, Lance_lat)
    ax.plot(xl,yl,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')

    outname='virtual_buoys_'+str(i)
    fig1.savefig(outpath+outname,bbox_inches='tight')
    plt.close()


#convert -delay 50 /Data/sim/polona/sid/deform/plots/virtual_buoys_* /Data/sim/polona/sid/deform/plots/vb_anim.gif
