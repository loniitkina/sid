from datetime import datetime
from glob import glob
import numpy as np
import pyresample as pr
from sid_func import *
import matplotlib.pyplot as plt


radius = 7000


#-------------------------------------------------------------------
inpath_drift = '../sidrift/data/40m_combo/'
inpath = '../sidrift/data/80m_stp10_single_filter/'
outpath_data = inpath
outpath = '../sidrift/plots/'
outpath = inpath
metfile = '../sidrift/data/10minute_nounits.csv'

#get Lance location and seed the parcels
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
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))
xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)

#create a grid of spacing of 400m in each direction in the radius from the ship
#Lance will be in the center
xbf = np.arange(xlp-radius, xlp+radius, 100)
ybf = np.arange(ylp-radius, ylp+radius, 100)
x_buoy,y_buoy = np.meshgrid(xbf,ybf)

#read in the drift product and calculate displacements from starting point until the first break in the SAR data (1 week)
fl = sorted(glob(inpath_drift+'SeaIceDrift*.npz'))

fl_dmg = sorted(glob(inpath+'Damage*.npz'))

#empty arrays to store bouy drifts and damage
x_path = np.zeros((7,x_buoy.shape[0],x_buoy.shape[1]))
y_path = np.zeros((7,x_buoy.shape[0],x_buoy.shape[1]))
damage = np.zeros((7,x_buoy.shape[0],x_buoy.shape[1]))

#start with initial coordinates
x_path[0,:,:] = x_buoy
y_path[0,:,:] = y_buoy


date = [start]

for i in range(0,6):        #just first 6 day-pairs, then 1 week gap)
    #read in all the data
    break
    print(fl[i])
    container = np.load(fl[i])
    lat1 = container['lat1']     #these arrays are as big as both scenes together (nan-where no overlap), ~5000x5000 of 80m-pixels = 400x400km
    lon1 = container['lon1']
    u = container['upm']
    v = container['vpm']
    
    
    #also get damage data!!!
    container = np.load(fl_dmg[i])
    lat_dmg = container['lat']
    lon_dmg = container['lon']
    dmg = container['d']
    

    
    x1,y1 = transform(inProj,outProj,lon1,lat1)
    xd,yd = transform(inProj,outProj,lon_dmg,lat_dmg)
    
    
    #print(lat_dmg)
    #lon_dmg1,lat_dmg1 = transform(outProj,inProj,xd,yd)
    #print(lat_dmg1)
    #exit()
    
    
    #get displacements
    #get time difference
    date1 = fl[i].split('_')[-2]
    date2 = fl[i].split('_')[-1].split('.')[0]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    
    diff = (dt2-dt1).seconds + (dt2-dt1).days*24*60*60
    dx = diff*u; dy = diff*v
    x2 = x1 + dx
    y2 = y1 + dy
    


    date.append(dt2)
    
    #find closest location to the individual virtual buoy
    for m in range(0,x_buoy.shape[0]):
        for n in range(0,y_buoy.shape[1]):
             
            mask = (x1>x_buoy[m,n]-500) & (x1<x_buoy[m,n]+500) & (y1>y_buoy[m,n]-500) & (y1<y_buoy[m,n]+500)    #40m resolution data, give some margin
                        
            mask_d = (xd>x_buoy[m,n]-900) & (xd<x_buoy[m,n]+900) & (yd>y_buoy[m,n]-900) & (yd<y_buoy[m,n]+900)

            #read the displacement from the end location
            #update VB location for the next step in time
            
            #aaa = np.ma.masked_invalid(x2[mask])
            #if aaa.shape[0] > 0:
                #print(aaa.shape)
                #print(np.mean(aaa))
            ##exit()
            
            x_buoy[m,n] = np.mean(np.ma.masked_invalid(x2[mask]))
            y_buoy[m,n] = np.mean(np.ma.masked_invalid(y2[mask]))
            #step 0 is initial condition, step 1 is reached after one day etc.
            x_path[i+1,m,n] = x_buoy[m,n]
            y_path[i+1,m,n] = y_buoy[m,n]
            
            #check if there was any deformation
            dd = np.mean(np.ma.masked_invalid(dmg[mask_d]))
            damage[i+1,m,n] = dd
            

    #generation of new parcels in empty spaces - all are damaged (lead ice)
    #we can create new grid every time step and interpolate the old one onto it? Then we loose the track and this not lagrangian anymore
    
    #triangulate between all original parcels
    #if there is no parcel inside the triangle, put there a new point with damage=1
    
    #other options is to track triangles (or centroids of triangles)
    #if area and minangle of tringles become small - merge and re-triangulate that part
    #if area becomes too large, add new centroids and re-triangulate

#final damage is the sum of all time steps and can be fraction


            

##save the output locations (in latlon) to be plotted on top of divergence and shear maps by sid_defrom.py
#lon_path,lat_path = transform(outProj,inProj,x_path,y_path)
#print(lat_path[:,0,0])

##dump data into numpy file
#out_file = outpath_data+'VB.npz'
#np.savez(out_file,x_path=x_path,y_path=y_path,lon_path=lon_path, lat_path=lat_path, damage = damage)

##save dates separatelly
#out_file = outpath_data+'VB_dates.npz'
#np.savez(out_file,date=date)

#make some plots
out_file = outpath_data+'VB.npz'
container = np.load(out_file)
print(container.files)
lon_path = container['lon_path']
lat_path = container['lat_path']
damage = container['damage']

out_file = outpath_data+'VB_dates.npz'
container = np.load(out_file, allow_pickle=True)
date = container['date']

radius_proj=55000
radius_proj=radius
from pyresample.geometry import AreaDefinition
#Using a projection dictionary
area_id = 'around Lance'
description = 'North Pole LAEA Europe'
proj_id = 'lance'
proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'a':6378137, 'b':6356752.3142, 'units':'m'}
width = radius_proj*2/400 #100 m spacing
height = radius_proj*2/400 #100 m spacing
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
    print(i)
    fig1    = plt.figure(figsize=(10,10))
    ax      = fig1.add_subplot(111)

    #x,y = m(lon_path[i,:,:],lat_path[i,:,:])
    #ax.pcolormesh(x,y,damage[i])
    
    
    
    dd = np.ma.masked_where(~(np.isfinite(lat_path[i,:,:])),damage[i]); dd = np.ma.compressed(dd)
    lop = np.ma.masked_where(~(np.isfinite(lat_path[i,:,:])),lon_path[i,:,:]); lop = np.ma.compressed(lop)
    print(lop)
    lap = np.ma.masked_invalid(lat_path[i,:,:]); lap = np.ma.compressed(lap)
    print(lap)
    print(lop.size,lap.size)
    x,y = m(lop,lap)
    print(x)
    print(lap)
    print(dd)
    ax.scatter(x,y,c=dd)



    cl = next(color)
    ##Lance initial position
    #xl, yl = m(Lance_lon, Lance_lat)
    #ax.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    
    #Lance moving with the VBs
    mi = np.argmin(abs(np.asarray(dtb)-date[i]))
    Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
    Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]
    xl, yl = m(Lance_lon, Lance_lat)
    ax.plot(xl,yl,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')

    outname='virtual_buoys_'+str(i)
    fig1.savefig(outpath+outname,bbox_inches='tight')
    plt.close()


#total damage plot
date = date[1]
damage = np.ma.masked_invalid(damage).filled(fill_value=0)
dtotal = np.sum(damage[:2,:,:],axis=0)

fig2    = plt.figure(figsize=(10,10))
ax      = fig2.add_subplot(111)


m = pr.plot.area_def2basemap(area_def)

m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

#x,y = m(lon_path[-1,:,:],lat_path[-1,:,:])
#ax.pcolormesh(x,y,dtotal)


dd = np.ma.masked_where(~(np.isfinite(lat_path[1,:,:])),dtotal); dd = np.ma.compressed(dd)
lop = np.ma.masked_where(~(np.isfinite(lat_path[1,:,:])),lon_path[1,:,:]); lop = np.ma.compressed(lop)
lap = np.ma.masked_invalid(lat_path[1,:,:]); lap = np.ma.compressed(lap)
x,y = m(lop,lap)
ax.scatter(x,y,c=dd)



mi = np.argmin(abs(np.asarray(dtb)-date))
Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]
xl, yl = m(Lance_lon, Lance_lat)
ax.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

##scale
#m.drawmapscale(Lance_lon, Lance_lat-.3, Lance_lon+8, Lance_lat-.2, 50, units='km', barstyle='fancy',fontsize=14)


outname='virtual_buoys_damage_21-26Jan2015'
fig2.savefig(outpath+outname,bbox_inches='tight')
plt.close()


#save this a geotiff
geotiff_file = outpath+outname+'.tiff'
save_geotiff(dtotal, area_def, geotiff_file)


















#convert -delay 50 /Data/sim/polona/sid/deform/plots/virtual_buoys_* /Data/sim/polona/sid/deform/plots/vb_anim.gif
