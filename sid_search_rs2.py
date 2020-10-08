from datetime import datetime, timedelta
import numpy as np
import pyresample as pr
import matplotlib.pyplot as plt
import os
from glob import glob
from pyproj import Proj, transform
from sid_func import *


inpath_bt = '../sidrift/data/backtraj_fs/'
inpath_poly = '../sidrift/data/RADARSAT-2/polygons/'
outpath = '../sidrift/data/RADARSAT-2/sequences/'


#read in RS-2 polygon data





#prepare for conversion to geographic projection
inProj = Proj(init='epsg:4326')
#OSI-SAF proj: proj4_string = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
outProj = Proj('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45 +units=km')


#get the list of all trajectory data and sort them by date
bt_fl = sorted(glob(inpath_bt+'*'))

#run through the trajectory and search for a matching polygon at same day
for fn in bt_fl:
    print(fn)
    
    tmp = getColumn(fn,0)
    date = [ datetime.strptime(tmp[i], "%Y-%m-%d %H:%M:%S") for i in range(len(tmp)) ]
    lon = np.asarray(getColumn(fn,1),dtype=float)
    lat = np.asarray(getColumn(fn,2),dtype=float)
    
    print(date)
    print(lat)
    
    
    
    for i in range(0,len(date)):
        print(date[i])
        
        #get daily bt coordintes
        xbt,ybt = transform(inProj,outProj,lon[i],lat[i])
        
        #make a list of all RS-2 meta files on that day
        #construct metafile name
        sd = date[i].strftime('%Y%m%d')
        fn_rs2 = glob(inpath_poly+'RS2_'+sd+'*.geojson')
        print(fn_rs2)
        
        
        #open those files and check if the centorid of that polygon is close enough to the lat,lon
        distance = []
        if len(fn_rs2) > 0:
            for j in fn_rs2:
                print(j)
                import json
                with open(j) as jsonfile:
                    data = json.load(jsonfile)
                
                poly = data['features'][0]['geometry']['coordinates'][0][0]
                print(poly)
                from shapely.geometry import Polygon
                poly = Polygon(poly)
                
                #get centroid
                cent = poly.centroid
                print(cent)                
                
                #transform coordinates
                xc,yc = transform(inProj,outProj,cent.x,cent.y)
                print(xc,yc)
                print(xbt,ybt)
                
                dx = xc-xbt
                dy = yc-ybt
                
                #units are km!
                d = np.sqrt(dx**2+dy**2)
                
                print(d)
                
                #collect all centroids
                distance.append(d)
            
            #select shortest distance
            #max distance should be half the RS-2 scene: 200km?
            mi = np.argmin(distance)
            print(mi)
            if distance[mi] < 200:
                rs2_fit = fn_rs2[mi].split('/')[-1]
            else:
                rs2_fit = 'No scene overlap'
            print(rs2_fit)

                
        else:
            print('Day with no RS-2 scenes')
        exit()
    
    
#allow 5 day time gaps
#get polygon centroid
#measure distanc between centroid and backtrajectory point
#store names of RS-2 scenes

#write out files with time, trajectory coordinates, distance and RS-2 scene name
