from datetime import datetime, timedelta
import numpy as np
from glob import glob
from netCDF4 import Dataset
from sid_func import *
import pyresample as pr
import matplotlib.pyplot as plt

inpath = '/Data/sim/data/OSISAF_ice_drift/'
outpath = '/Data/sim/polona/sid/deform/plots/'
metfile = '../data/10minute_nounits.csv'
radius = 50000
sr = 35000  #search radius for OSI-SAF drift

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
print(crnx,crny)

#OSI_SAF drift grid
fn = '/Data/sim/data/OSISAF_ice_drift/2015/01/ice_drift_nh_polstere-625_multi-oi_201501161200-201501181200.nc'
f = Dataset(fn)
#lats = f.variables['lat'][:]
#lons = f.variables['lon'][:]
xc = f.variables['xc'][:]*1000  #convert from km to m
yc = f.variables['yc'][:]*1000
f.close()

xcm,ycm = np.meshgrid(xc,yc)

#Set up the plot
fig1    = plt.figure(figsize=(20,20))
ax      = fig1.add_subplot(111)
area_def = pr.utils.load_area('area.cfg', 'leg1')  
m = pr.plot.area_def2basemap(area_def)
m.drawmapboundary(fill_color='#9999FF')
m.drawcoastlines()
m.fillcontinents(color='#ddaa66',lake_color='#9999FF')
##Draw parallels and meridians
m.drawparallels(np.arange(80.,86.,1),labels=[1,0,0,0], fontsize=16,latmax=85.)
m.drawmeridians(np.arange(-5.,90.,5.),latmax=85.,labels=[0,0,0,1], fontsize=16)



for i in range(0,crnx.size):
    print(i)
    ice=True
    date = maptime
    lon,lat = transform(outProj,inProj,crnx[i], crny[i])
    bt_lon = [lon]
    bt_lat = [lat]
    while ice==True:
        #find the corresponding OSI-SAF drift file/date
        try:            #some files are missing >> in such case contine with the old displacements
            ym = datetime.strftime(date,'%Y/%m')
            ymd = datetime.strftime(date,'%Y%m%d')
            fn = glob(inpath+ym+'/ice_drift_nh_polstere*'+ymd+'1200.nc')[0]
            print(fn)
            f = Dataset(fn)
            dx = f.variables['dX'][0,:,:]*1000          #convert from km to m
            dy = f.variables['dY'][0,:,:]*1000
        except:
            print('Use old data.')

        #find closest desplacement
        mask = (xcm>crnx[i]-sr) & (xcm<crnx[i]+sr) & (ycm>crny[i]-sr) & (ycm<crny[i]+sr)
        dxm = np.mean(np.ma.masked_invalid(dx[mask]))/2         #this is displacement for 2 days
        dym = np.mean(np.ma.masked_invalid(dy[mask]))/2
        
        print(dx[mask])
        print(dx[mask].shape)
        print(dxm)
        
        #stop where no more drift data (future: concentration falls below 50%)
        if abs(dxm) > 0:
            #calculate position back in time
            crnx[i] = crnx[i] - dxm
            crny[i] = crny[i] + dym
            
            #transform back to latlon
            lon,lat = transform(outProj,inProj,crnx[i], crny[i])
            print(lon,lat)
            bt_lon.append(lon)
            bt_lat.append(lat)
            
            #walk back in time by 1 day
            date = date - timedelta(days=1)
        
        else:
            ice=False
        
    #save the backtrajectory for this corner
    print(bt_lon)
    print(bt_lat)
    
    #plot the trajectory
    x,y = m(bt_lon,bt_lat)
    ax.plot(x,y,'o')

outname='backtraj'
fig1.savefig(outpath+outname,bbox_inches='tight')
plt.close()



#If this works do same aslo for the CSM image corners
