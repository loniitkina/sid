from datetime import datetime, timedelta
import numpy as np
from glob import glob
from netCDF4 import Dataset
from sid_func import *
import pyresample as pr
import matplotlib.pyplot as plt

inpath = '/Data/sim/data/OSISAF_ice_drift/'
inpath_ic = '/Data/sim/data/OSISAF_ice_conc/polstere/'
inpath_summer = '/Data/sim/polona/sid/nsidc/'    #Downloaded as daily values for a whole year (2014) from https://nsidc.org/data/nsidc-0116
outpath = '/Data/sim/polona/sid/deform/plots/'
metfile = '../data/10minute_nounits.csv'
radius = 50000
#search radius for OSI-SAF drift
sr = 35000  #~half resolution of OSI-SAF sea ice drift product (NSIDC drift is on 25km grid)
#search radius for OSI-SAF ice concentration
sr_ic = 10000

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
#outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))
#OSI-SAF proj: proj4_string = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
outProj = Proj('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)

#corners
crnx = np.empty((4),dtype=float)
crny = np.empty((4),dtype=float)
crnx[0],crny[0] = (xlp-radius,ylp-radius)
crnx[1],crny[1] = (xlp+radius,ylp-radius)
crnx[2],crny[2] = (xlp+radius,ylp+radius)
crnx[3],crny[3] = (xlp-radius,ylp+radius)
#print(crnx,crny)

#OSI_SAF drift grid
fn = inpath+'2015/01/ice_drift_nh_polstere-625_multi-oi_201501161200-201501181200.nc'
f = Dataset(fn)
#lats = f.variables['lat'][:]
#lons = f.variables['lon'][:]
xc = f.variables['xc'][:]*1000  #convert from km to m
yc = f.variables['yc'][:]*1000
f.close()
#make grid
xcm,ycm = np.meshgrid(xc,yc)

#OSI_SAF ice concentration grid
#proj4_string = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
fn = inpath_ic+'2015_nh_polstere/ice_conc_nh_polstere-100_multi_201501011200.nc'
f = Dataset(fn)
xc = f.variables['xc'][:]*1000  #convert from km to m
yc = f.variables['yc'][:]*1000
f.close()
#make grid
xcm_ic,ycm_ic = np.meshgrid(xc,yc)

#NSIDC for summer
#Tschudi, M. A., Meier, W. N., and Stewart, J. S.: An enhancement to sea ice motion and age products, The Cryosphere Discuss., https://doi.org/10.5194/tc-2019-40, in review, 2019.
#proj4text = "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"
fn = inpath_summer+'icemotion_daily_nh_25km_20140101_20141231_v4.1.nc'
f = Dataset(fn)
#lats = f.variables['latitude'][:]
#lons = f.variables['longitude'][:]
x = f.variables['x'][:]
y = f.variables['y'][:]
u = f.variables['u'][:]
v = f.variables['v'][:]
f.close()
#keep
nsidcProj = Proj('+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs')
#make grid
xcm_s,ycm_s = np.meshgrid(x,y)

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

#do backtrajectories for each corner point
for i in range(0,crnx.size):
    print(i)
    ice=True
    date = maptime
    lon,lat = transform(outProj,inProj,crnx[i], crny[i])
    bt_lon = [lon]
    bt_lat = [lat]
    ex=0
    
    while ice==True:
        #check if this place has enough of ice at all!
        #find the corresponding OSI-SAF ice concentration file/date
        ymd = datetime.strftime(date,'%Y%m%d')
        try:            #some files are missing >> in such case contine with the old displacements
            fn_ic = glob(inpath_ic+str(date.year)+'_nh_polstere/ice_conc_nh_polstere-100_multi_*'+ymd+'1200.nc')[0]
            #print(fn)
            f = Dataset(fn_ic)
            ic = f.variables['ice_conc'][0,:,:]
        except:
            print(date)
            print('Use old ice concentration data: '+fn_ic)
        
        mask = (xcm_ic>crnx[i]-sr_ic) & (xcm_ic<crnx[i]+sr_ic) & (ycm_ic>crny[i]-sr_ic) & (ycm_ic<crny[i]+sr_ic)
        icm = np.mean(np.ma.masked_invalid(ic[mask]))
        print(icm)
        #exit()
        
        if icm < 50:
            print(icm)
            ice=False
            print('Ice concentration bellow treshold.');continue

        #find the corresponding OSI-SAF drift file/date
        try:            #some files are missing >> in such case contine with the old displacements
            ym = datetime.strftime(date,'%Y/%m')
            fn = glob(inpath+ym+'/ice_drift_nh_polstere*'+ymd+'1200.nc')[0]
            #print(fn)
            f = Dataset(fn)
            dx = f.variables['dX'][0,:,:]*1000          #convert from km to m
            dy = f.variables['dY'][0,:,:]*1000
            #lon2 = f.variables['lon1'][0,:,:]
            #lat2 = f.variables['lat1'][0,:,:]
        except:
            print(date)
            #print('Use old data: '+fn)
        

        #find closest displacement
        mask = (xcm>crnx[i]-sr) & (xcm<crnx[i]+sr) & (ycm>crny[i]-sr) & (ycm<crny[i]+sr)
        dxm = np.mean(np.ma.masked_invalid(dx[mask]))/2         #this is displacement for 2 days
        dym = np.mean(np.ma.masked_invalid(dy[mask]))/2
        dym = -dym           #account for the OSI-SAF sign convention
        #print(dx[mask].shape)
        
        ##alternative displacement estimate
        #xcmm = np.ma.masked_invalid(xcm[mask])
        #ycmm = np.ma.masked_invalid(ycm[mask])
        #lon2m = np.ma.masked_invalid(lon2[mask])
        #lat2m = np.ma.masked_invalid(lat2[mask])
        #xcm2m,ycm2m = transform(inProj,outProj,lon2m, lat2m)
        #dx2 = np.mean(xcm2m-xcmm)/2
        #dy2 = np.mean(ycm2m-ycmm)/2
        ##print(dxm,dym)
        ##print(dx2,dy2)      #The distance estimate is different by ~100cm --OK, dY has inverted sign.
        #dxm = dx2
        #dym = -dy2           #account for the OSI-SAF sign convention
                
        if (date.month<11) & (date.month>4):
            if date.year<2014:
                print('Too far back!')
                ice=False
            summer=True
            #print('Summer: no OSI-SAF data - use alternative!')
            #use NSIDC drift to bridge the summer
            #get the rigt julian date
            idx = date.timetuple().tm_yday -1   #account for python indexing
            #print(idx)
            uj = u[idx,:,:]
            vj = v[idx,:,:]
            
            #Transform coordinates
            crnx[i],crny[i] = transform(outProj,nsidcProj,crnx[i],crny[i])
            
            #find closest grid point
            mask = (xcm_s>crnx[i]-sr) & (xcm_s<crnx[i]+sr) & (ycm_s>crny[i]-sr) & (ycm_s<crny[i]+sr)
            invalid = uj[mask]==-9999.
            dxm = np.mean(np.ma.array(uj[mask],mask=invalid))*24*60*60/100   #convert from velocity (cm/s) to displacement (m)
            dym = np.mean(np.ma.array(vj[mask],mask=invalid))*24*60*60/100
            #print(uj[mask].shape)
            #print(uj[mask])
            #print(dxm,dym)
            
        else:
            summer=False
            
        #stop where no more drift data (future: concentration falls below 50%)
        if abs(dxm) > 0:
            #calculate position back in time
            crnx[i] = crnx[i] - dxm
            crny[i] = crny[i] - dym
            
            #transform back to latlon and to the OSI-SAF coordinates
            if summer:
                lon,lat = transform(nsidcProj,inProj,crnx[i], crny[i])
                print(lon,lat)
                crnx[i],crny[i] = transform(nsidcProj,outProj,crnx[i],crny[i])
                #exit()
            else:
                lon,lat = transform(outProj,inProj,crnx[i], crny[i])
            #print(lon,lat)
            bt_lon.append(lon)
            bt_lat.append(lat)
            
            #walk back in time by 1 day
            date = date - timedelta(days=1)
        else:
            ice=False
            print('No more ice!')
            print(date)
        
    #save the backtrajectory for this corner
    print(bt_lon[-1])
    print(bt_lat[-1])
    
    #plot the trajectory
    x,y = m(bt_lon,bt_lat)
    ax.plot(x,y,'o',label='corner %s' %(i))

ax.legend()
outname='backtraj'
fig1.savefig(outpath+outname,bbox_inches='tight')
plt.close()



#If this works do same aslo for the CSM image corners
