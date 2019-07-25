from datetime import datetime
from glob import glob
import numpy as np
import pyresample as pr
from sid_func import *
import matplotlib.pyplot as plt


radius = 25000
file_name_end = '_50km_more.csv'

#-------------------------------------------------------------------
inpath = '/Data/sim/polona/sid/drift_full_stp1/'
outpath_data = '/Data/sim/polona/sid/deform/'
outpath = outpath_data+'plots/'
metfile = '../data/10minute_nounits.csv'
reg = 'leg1'
proj = reg

#get Lance location
start = datetime(2015,1,21)                 #this does not coincide with the 1st SAR image date/time
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

print(Lance_lon, Lance_lat)       
print(xlp,ylp)

        
#from pyresample.geometry import AreaDefinition
##Using a projection dictionary
#area_id = 'around Lance'
#description = 'North Pole LAEA Europe'
#proj_id = 'lance'
#proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'a':6378137, 'b':6356752.3142, 'units':'m'}
#width = radius*2/100 #100 m spacing
#height = radius*2/100 #100 m spacing
#area_extent = (xlp-radius,ylp-radius,xlp+radius,ylp+radius)
#area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)

#print(area_def)
#m = pr.plot.area_def2basemap(area_def)
            
#xlf, ylf = m(Lance_lon, Lance_lat)
#create a grid of spacing of 2km in each direction until half of the radius from the ship
xbf = np.arange(xlp-(radius/2), xlp+(radius/2), 5000)
ybf = np.arange(ylp-(radius/2), ylp+(radius/2), 5000)




x_buoy,y_buoy = np.meshgrid(xbf,ybf)

print(x_buoy.shape)
print(xbf[0])
print(x_buoy[0,0])
#exit()



#read in the drift product and calculate displacements from starting point until the first break in the SAR data (1 week)
fl = sorted(glob(inpath+'*.npz'))
#print(fl)
#exit()

#empty arrays to store bouy drifts
x_path = np.zeros((7,x_buoy.shape[0],x_buoy.shape[1]))
y_path = np.zeros((7,x_buoy.shape[0],x_buoy.shape[1]))

x_path[0,:,:] = x_buoy
y_path[0,:,:] = y_buoy

for i in range(0,5):        #just first 5 days (max 7, then 1 week gap)
    #read in all the data
    print(fl[i])
    container = np.load(fl[i])
    u = container['upm']            #these arrays are as big as both scenes together (nan-where no overlap), ~5000x5000 of 40m-pixels = 200x200km
    v = container['vpm'] 
    lat = container['lat1']
    lon = container['lon1']
    
    #x, y = m(lon, lat)
    
    x,y = transform(inProj,outProj,lon,lat)
    
    #print(lat)
    #print(lon)
    
    
    ##plt.imshow(u)
    ##plt.show()
    
    #print(x)
    #print(x.shape)
    #print(np.argmin(np.ma.masked_invalid(lat),axis=0))
    #print(np.argmin(np.ma.masked_invalid(lat),axis=1))
    #mask = np.ma.masked_invalid(lat) > 0.
    #print(lat[mask])
    #print(x[mask])
    #print(lon[100,3128],lat[100,3128])
    #print(x[100,3128],y[100,3128])
    
    #outlon,outlat = transform(outProj,inProj,x[100,3128],y[100,3128])
    #print(outlon)
    #print(outlat)
    #exit()
    
    
    ##print(np.max(x))
    ##print(np.min(x))
    ##print(np.max(y))
    ##print(np.min(y))
    ##print(np.max(lat))
    ##print(np.min(lat))
    
    #get time difference
    date1 = fl[i].split('_')[-2]
    date2 = fl[i].split('_')[-1].split('.')[0]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")

    diff = (dt2-dt1).seconds + (dt2-dt1).days*24*60*60
    #print(diff)

    
    #find closest location to the individual virtual buoy
    for m in range(0,x_buoy.shape[0]):
        for n in range(0,y_buoy.shape[1]):
            print(x_buoy[m,n])
            print(y_buoy[m,n])
            
            mask = (x>x_buoy[m,n]-60) & (x<x_buoy[m,n]+60) & (y>y_buoy[m,n]-60) & (y<y_buoy[m,n]+60)    #40m resolution data, give some margin
            print(x[mask])
            print(y[mask])
            #exit()
            
            usample = np.ma.masked_invalid(u[mask])
            vsample = np.ma.masked_invalid(v[mask])
            print(usample.shape)
            #exit()
            
            #calculate displacement
            dx = diff*np.mean(usample)
            dy = diff*np.mean(vsample)
            print(dx,dy)
            
            #new location
            x_buoy[m,n] = x_buoy[m,n]+dx
            y_buoy[m,n] = y_buoy[m,n]+dy
            
            x_path[i,m,n] = x_buoy[m,n]
            y_path[i,m,n] = y_buoy[m,n]
            
print(x_path[:,0,0])
#save the output locations (in latlon) to be plotted on top of divergence and shear maps by sid_defrom.py


lon_path,lat_path = transform(outProj,inProj,x_path,y_path)

#dump data into numpy file
out_file = outpath_data+'VB.npz'
np.savez(out_file,x_path=x_path,y_path=y_path,lon_path=lon_path, lat_path=lat_path)

#make some plots
out_file = outpath_data+'VB.npz'
container = np.load(out_file)
print(container.files)
lon_path = container['lon_path'][:5,:,:]
lat_path = container['lat_path'][:5,:,:]

#print(lon_path[4,:,:])
#print(lat_path[4,:,:])


#area_def = pr.utils.load_area('area.cfg', proj)  


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

    #Lance
    xl, yl = m(Lance_lon, Lance_lat)
    ax.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

    outname='virtual_buoys_'+str(i)
    fig1.savefig(outpath+outname,bbox_inches='tight')
    plt.close()


#convert -delay 50 /Data/sim/polona/sid/deform/plots/virtual_buoys_* /Data/sim/polona/sid/deform/plots/vb_anim.gif
