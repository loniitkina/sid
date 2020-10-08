from datetime import datetime, timedelta
import numpy as np
#from glob import glob
#from netCDF4 import Dataset
#from sid_func import *
import pyresample as pr
import matplotlib.pyplot as plt
import os
from glob import glob

import xarray as xr
from pyproj import Proj, transform

outpath = '../sidrift/plots/'
outpath_bt = '../sidrift/data/backtraj_fs/'


#search radius for OSI-SAF drift (in km)
sr = 34  #~half resolution of OSI-SAF sea ice drift product
sr = 64  #full resolution will get more data
sr = 100 #get even more data
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

##N-ICE region
#moor_lat = [82,82,82,82]; moor_lon = [5,10,15,20]
#xmoor_start,ymoor_start = transform(inProj,outProj,moor_lon, moor_lat)


print(xmoor_start,ymoor_start)

#OSI_SAF drift
ds_drift =xr.open_dataset('https://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/drift_lr_nh_agg')

#OSI_SAF ice concentration
ds_ic =xr.open_dataset('https://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/conc_nh_pol_agg')

years = np.arange(2016,2021)
#years = [2010]
months = [1,2,3,4,5,12,11]
cls = iter(plt.cm.rainbow_r(np.linspace(0,1,len(months)+1)))

for yr in years:
    print('year',yr)
    
    #Set up a new plot for every year
    fig1    = plt.figure(figsize=(20,20))
    ax      = fig1.add_subplot(111)
    area_def = pr.utils.load_area('area.cfg', 'fram_strait')  
    m = pr.plot.area_def2basemap(area_def)
    m.drawmapboundary(fill_color='#9999FF')
    m.drawcoastlines()
    m.fillcontinents(color='#ddaa66',lake_color='#9999FF')
    ##Draw parallels and meridians
    m.drawparallels(np.arange(60.,86.,5),labels=[1,0,0,0], fontsize=16,latmax=85.)
    m.drawmeridians(np.arange(-180.,180.,20.),latmax=85.,labels=[0,0,0,1], fontsize=16)
    
    for mon in months:
        print('month',mon)
        
        if mon==2:
            days=range(1,29)
        elif mon==4 or mon==11:
            days=range(1,31)
        else:
            days=range(1,32)
        
        #get new color for each month
        cl = next(cls)        
        
        for day in days:
            print('day',day)
        
            start = datetime(yr,mon,day,0,0,0)
            xmoor = xmoor_start.copy(); ymoor = ymoor_start.copy()      #use .copy() or xmoor_start will become reference of xmoor
            
            #do backtrajectories for each mooring point
            for i in range(0,len(xmoor)):
                print('######################################################################### Mooring: '+str(i))
                date = start
                bt_lon = [moor_lon[i]]
                bt_lat = [moor_lat[i]]
                bt_date = [start]
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
                    #this will get only one value (if that value is nan, there is nothing to mask for average calculation)
                    #dx = ds_drift.dX.sel(time=date,xc=xmoor[i],yc=ymoor[i],method='bfill')
                    #this will get several inside the spatial slice
                    #the interval at yc has to be turned around!!!
                    dx = ds_drift.dX.sel(time=date,method='nearest').sel(xc=slice(xmoor[i]-sr,xmoor[i]+sr),yc=slice(ymoor[i]+sr,ymoor[i]-sr))   #IS THIS INTERVAL FOR YC ONLY SUCH UNTIL 2016???
                    #print(dx.values)
                    dxm = np.mean(np.ma.masked_invalid(dx.values))/2         #this is displacement for 2 days
                    
                    #dy = ds_drift.dY.sel(time=date,method='bfill',xc=xmoor[i],yc=ymoor[i])
                    dy = ds_drift.dY.sel(time=date,method='nearest').sel(xc=slice(xmoor[i]-sr,xmoor[i]+sr),yc=slice(ymoor[i]+sr,ymoor[i]-sr))
                    dym = np.mean(np.ma.masked_invalid(dy.values))/2         #this is displacement for 2 days
                    
                    #the early version had different sign, data was never reprocessed
                    #sign changed already in November, December 2015???
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
                        bt_date.append(date)
                        
                        #walk back in time by 1 day
                        date = date - timedelta(days=1)
                        
                        #just stop when reaching freeze-up:
                        if date < datetime(yr-1,11,1,0,0,0): print('freeze up', date); ice=False
                        
                        #also stop if the ice dirfts north of 87N (SAR polar hole)
                        if lat>87: print('polar hole', date); ice=False
                        
                    else:
                        #if we start in December or November (same year) freeze-up will not be identified and the loop will end here
                        ice=False
                        print('No more ice!', date)
                        #print(dx.values) #just nans
                    
                #print whre we end up
                print('trajectory of days: ',len(bt_lat), 'and final location: ', bt_lat[-1], bt_lon[-1])
                
                #write out text files with dates and coordinates
                tt = [bt_date, bt_lon, bt_lat]
                table = zip(*tt)
                #adjusted to python3:
                table = list(zip(*tt))

                output = outpath_bt + 'bt_'+start.strftime('%Y%m%d')+'_'+str(i)+'.csv'
                #remove the output from previous script runs (or it will be appended) - if exists
                rlist = glob(output)
                for fn in rlist:
                    os.remove(output)
                with open(output, 'ab') as f:
                    np.savetxt(f, table, fmt="%s", delimiter=",")
                
                
                #plot all the trajectories
                #label just the first mooring color every month
                x,y = m(bt_lon,bt_lat)
                if day==1 and i == 0:
                    ax.plot(x,y,'o',label=str(mon),c=cl,alpha=0.2)
                else:
                    ax.plot(x,y,'o',c=cl,alpha=0.2)


    #year level
        #month level
            #day lavel
                #mooring level
    ax.legend(loc='upper left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
    outname='backtraj_fs_'+str(yr)
    fig1.savefig(outpath+outname,bbox_inches='tight')
    plt.close()
        
    #exit()
