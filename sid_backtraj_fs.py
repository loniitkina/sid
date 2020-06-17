from datetime import datetime, timedelta
import numpy as np
#from glob import glob
#from netCDF4 import Dataset
#from sid_func import *
import pyresample as pr
import matplotlib.pyplot as plt

import xarray as xr
from pyproj import Proj, transform

outpath = '../sidrift/plots/'

#search radius for OSI-SAF drift (in km)
sr = 34  #~half resolution of OSI-SAF sea ice drift product
#search radius for OSI-SAF ice concentration
sr_ic = 6    #5km resolution data

##mooring position
inProj = Proj(init='epsg:4326')
#OSI-SAF proj: proj4_string = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
outProj = Proj('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45 +units=km')

#Fram Strait moorings
#[ -8, 78.83, 0, 79 ]   #approx bounding box
moor_lat = [80,80,80,80,80]; moor_lon = [-8,-6,-4,-2,0]
xmoor_start,ymoor_start = transform(inProj,outProj,moor_lon, moor_lat)

#N-ICE region
moor_lat = [82,82,82,82]; moor_lon = [5,10,15,20]
xmoor_start,ymoor_start = transform(inProj,outProj,moor_lon, moor_lat)


print(xmoor_start,ymoor_start)

#OSI_SAF drift
ds_drift =xr.open_dataset('https://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/drift_lr_nh_agg')

#OSI_SAF ice concentration
ds_ic =xr.open_dataset('https://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/conc_nh_pol_agg')

#Set up the plot
fig1    = plt.figure(figsize=(20,20))
ax      = fig1.add_subplot(111)
area_def = pr.utils.load_area('area.cfg', 'arctic')  
m = pr.plot.area_def2basemap(area_def)
m.drawmapboundary(fill_color='#9999FF')
m.drawcoastlines()
m.fillcontinents(color='#ddaa66',lake_color='#9999FF')
##Draw parallels and meridians
m.drawparallels(np.arange(60.,86.,5),labels=[1,0,0,0], fontsize=16,latmax=85.)
m.drawmeridians(np.arange(-180.,180.,20.),latmax=85.,labels=[0,0,0,1], fontsize=16)

years = np.arange(2010,2021)
#years = [2020]
cls = iter(plt.cm.rainbow_r(np.linspace(0,1,len(years)+1)))

for yr in years:
    start = datetime(yr,4,1,0,0,0)
    xmoor = xmoor_start.copy(); ymoor = ymoor_start.copy()      #use .copy() or xmoor_start will become reference of xmoor
    cl = next(cls)
    
    #do backtrajectories for each mooring point
    for i in range(0,len(xmoor)):
        print('######################################################################### Mooring: '+str(i))
        date = start
        bt_lon = [moor_lon[i]]
        bt_lat = [moor_lat[i]]
        print(date,moor_lon[i],moor_lat[i])
        
        ice=True
        
        while ice==True:
            ic = ds_ic.ice_conc.sel(time=date,xc=xmoor[i],yc=ymoor[i],method='nearest')
            #ic = ds_ic.ice_conc.sel(time=date).sel(xc=xmoor[i],yc=ymoor[i],method='nearest',tolerance=sr_ic)
            #print(ic.values)
            icm = np.mean(np.ma.masked_invalid(ic.values))
            
            if icm < 70:
                ice=False
                print(date)
                print('Ice concentration bellow treshold: '+str(icm));continue

            #find the corresponding OSI-SAF drift file/date
            #dx = ds_drift.dX.sel(time=date,method='nearest').sel(xc=xmoor[i],yc=ymoor[i],method='pad')
            dx = ds_drift.dX.sel(time=date,xc=xmoor[i],yc=ymoor[i],method='bfill')
            
            #print(dx)
            #print(dx.values)
            #exit()
            
            
            dxm = np.mean(np.ma.masked_invalid(dx.values))/2         #this is displacement for 2 days
            
            dy = ds_drift.dY.sel(time=date,method='bfill',xc=xmoor[i],yc=ymoor[i])
            #dy = ds_drift.dX.sel(time=slice(date - timedelta(days=3),date - timedelta(days=3)))\
                #.sel(xc=xmoor[i],yc=ymoor[i],method='pad')
            dym = np.mean(np.ma.masked_invalid(dy.values))/2         #this is displacement for 2 days
            
            #the early version had different sign, data was never reprocessed
            if yr < 2016:
                dym = dym*-1
            
            #print(dx.values,dy.values)
            
            if abs(dxm)*abs(dym) > 0:
                #calculate position back in time
                xmoor[i] = xmoor[i] - dxm
                ymoor[i] = ymoor[i] - dym
                
                #transform back to latlon and to the OSI-SAF coordinates
                lon,lat = transform(outProj,inProj,xmoor[i], ymoor[i])
                #store back-trajectory
                bt_lon.append(lon)
                bt_lat.append(lat)
                
                #walk back in time by 1 day
                date = date - timedelta(days=1)
                
                #just stop when reaching freeze-up:
                if date < datetime(yr-1,11,1,0,0,0): print('freeze up'); ice=False
                
                
            else:
                ice=False
                print('No more ice!')
                print(date)
            
        #check the end coordinates
        print(ice); print(xmoor)
        print(bt_lon[-1])
        print(bt_lat[-1])
        
        #plot the trajectory
        x,y = m(bt_lon,bt_lat)
        if i == 0:
            #ax.plot(x,y,'o',label='ULS %s' %(i),c=cl)
            ax.plot(x,y,'o',label=str(yr),c=cl,alpha=0.2)
        else:
            ax.plot(x,y,'o',c=cl,alpha=0.2)

ax.legend(loc='upper left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
outname='backtraj_nice_2010-2020'
fig1.savefig(outpath+outname,bbox_inches='tight')
plt.close()
