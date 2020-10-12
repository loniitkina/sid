from datetime import datetime, timedelta
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import transform as shapely_transform
import pyresample as pr
import matplotlib.pyplot as plt
import os
from glob import glob
from pyproj import Proj, transform
from sid_func import *


inpath_bt = '../sidrift/data/backtraj_fs/'
inpath_poly = '../sidrift/data/RADARSAT-2/polygons/'
outpath = '../sidrift/data/RADARSAT-2/sequences/'

#prepare for conversion to geographic projection
inProj = Proj(init='epsg:4326')
#OSI-SAF proj: proj4_string = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
outProj = Proj('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45 +units=km')


#get the list of all trajectory data and sort them by date
bt_fl = sorted(glob(inpath_bt+'*_2.csv'))

#run through the trajectory and search for a matching polygon at same day
for fn in bt_fl:
    print(fn)
    
    #define which mooring this is
    mn = fn.split('_')[-1].split('.')[0]
    print(mn)
    
    rs2_fit_list = []
    mindist = []
    
    tmp = getColumn(fn,0)
    date = [ datetime.strptime(tmp[i], "%Y-%m-%d %H:%M:%S") for i in range(len(tmp)) ]
    lon = np.asarray(getColumn(fn,1),dtype=float)
    lat = np.asarray(getColumn(fn,2),dtype=float)
    
    #print(date)
    #print(lat)
    
    #prepare a sample trajectory plot
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

    #allow up to 5-day-long gaps in time
    gap=0    
    
    for i in range(0,len(date)):
        print(date[i])
        print(gap)
        
        
            
        #get daily bt coordintes
        xbt,ybt = transform(inProj,outProj,lon[i],lat[i])
        
        #make a list of all RS-2 meta files on that day
        #construct metafile name
        sd = date[i].strftime('%Y%m%d')
        fn_rs2 = glob(inpath_poly+'RS2_'+sd+'*.geojson')
        #print(fn_rs2)
        
        
        #open those files and check if the centorid of that polygon is close enough to the lat,lon
        distance = []
        tmp_store_cent_x = []
        tmp_store_cent_y = []
        tmp_store_poly = []
        if len(fn_rs2) > 0:
            for j in fn_rs2:
                #print(j)
                import json
                with open(j) as jsonfile:
                    data = json.load(jsonfile)
                
                poly = data['features'][0]['geometry']['coordinates'][0][0]
                #print(poly)
                poly = Polygon(poly)
                tmp_store_poly.append(poly)
                
                #get centroid
                cent = poly.centroid
                tmp_store_cent_x.append(cent.x)
                tmp_store_cent_y.append(cent.y)
                #print(cent)                
                
                #transform coordinates
                xc,yc = transform(inProj,outProj,cent.x,cent.y)
                #print(xc,yc)
                #print(xbt,ybt)
                
                dx = xc-xbt
                dy = yc-ybt
                
                #units are km!
                d = np.sqrt(dx**2+dy**2)
                
                #print(d)
                
                #collect all centroids
                distance.append(d)
            
            #select shortest distance
            #RS-2 scanSARwide is 500x500km
            #print(distance)
            mi = np.argmin(distance)
            #print(mi)
            if distance[mi] < 200:
                rs2_fit = fn_rs2[mi].split('/')[-1]
                
                #plot some of this
                #plot BT point for that day
                xbt,ybt = m(lon[i],lat[i])
                ax.scatter(xbt,ybt) #point
                #plot RS-2 polygon
                poly = tmp_store_poly[mi]
                x,y = poly.exterior.xy
                xrs2,yrs2 = m(x,y)
                ax.plot(xrs2,yrs2)  #line 
                #plot RS-2 centroid too
                xrs2c,yrs2c = m(tmp_store_cent_x[mi],tmp_store_cent_y[mi])
                ax.scatter(xrs2c,yrs2c,marker='*') #point
                
                
                #reset gap
                gap=0
                
                
                
            else:
                rs2_fit = 'Day with no good scene overlap'
                gap=gap+1
            
            mindist.append(distance[mi])

                
        else:
            #There are no scenes on the weekends
            rs2_fit = 'Day with no RS-2 scenes' 
            mindist.append('-')
            gap=gap+1
        
        if gap < 6:
            print(rs2_fit)
            rs2_fit_list.append(rs2_fit)
        else:
            break
        
        
    #print(rs2_fit_list)
    #print(mindist)
    
    if len(rs2_fit_list) > 30:
        print('Long time series!',len(rs2_fit_list), 'days')
            
        #save figure
        outname='sequence_'+date[0].strftime('%Y%m%d')+'_'+mn
        plt.savefig(outpath+outname)
        
        #save data
        #write out files with time, trajectory coordinates, distance and RS-2 scene name
        #exit()
        
    else:
        #close figure without saving
        plt.close()
    

#as a rule scenes follow the ice edge (as they were ordered to make ice charts)
#sea ice drifts mainly in N-S direction and the track is lost (reaches width of one RS-2 scene=500km in one month)

#as the scene is 500km wide, it is very likely that all moorings will be covered by same scenes
#so, to follow one mooring location (e.g. 2) is enough!

#usually when a trajectory from certain is well covered, there are several other days around it, covered by same sequence of scenes
#beginning of 2012 there is 2 months of such data (slow drift?)
#February 2013
#February 2014
