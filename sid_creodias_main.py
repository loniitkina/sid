import numpy as np
from datetime import datetime
from sid_func import getColumn
from pyproj import Proj, transform
from sid_creodias_func import download_from_polygon
from glob import glob

#This script will make a polygon of 100x100km around Polarstern and download the S-1 scenes to out_folder
out_folder='/Data/pit000/ResearchData/IFT/EarthObservation/MOSAIC/SAR/Sentinel-1'

#read the PS positions
#ps_file='../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0.csv'
#ps_file='../../downloads/data_master-solution_mosaic-leg2-20191214-20200224-floenavi-refstat-v1p0.csv'
ps_file='../../downloads/position_leg3_nh-track.csv'

#cover all quadrants
ps_files=sorted(glob('../../downloads/position_leg3_nh-track_[e,w,n,s,se,sw,nw,ne].csv')+glob('../../downloads/position_leg3_nh-track_[se,sw,nw,ne]?.csv'))
print(ps_files)

for ps_file in ps_files:
    print(ps_file)
    dt = getColumn(ps_file,0)
    lon = getColumn(ps_file,1)
    lat = getColumn(ps_file,2)
    lon = np.array(lon,dtype=np.float)
    lat = np.array(lat,dtype=np.float)
    dt = [ datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in dt ]

    #get 12UTC positions
    full_hour= np.where(np.array([ x.minute for x in dt ])==0,1,0)
    noon=np.where(np.array([ x.hour for x in dt ])==12,1,0)*full_hour

    lat_noon=np.ma.array(lat,mask=noon==0).compressed()
    lon_noon=np.ma.array(lon,mask=noon==0).compressed()
    dt_noon=np.ma.array(dt,mask=noon==0).compressed()

    outProj = Proj(init='epsg:4326')
    FloeNaviProj = Proj('+proj=stere +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +ellps=WGS84'%(0,0))
    radius=50000    #50km

    fns=[]
    dts=[]

    #cycle through every day
    for i in range(0,len(dt_noon)):
        coverage=0.8
        
        #convert to geogaphical coordinates
        x,y = transform(outProj,FloeNaviProj,lon_noon[i],lat_noon[i])
        
        #get corners
        ul_corner_x=x-radius
        ul_corner_y=y+radius
        lr_corner_x=x+radius
        lr_corner_y=y-radius

        #convert back to lat,lon
        ul_lon,ul_lat = transform(FloeNaviProj,outProj,ul_corner_x,ul_corner_y)
        lr_lon,lr_lat = transform(FloeNaviProj,outProj,lr_corner_x,lr_corner_y)

        #extend the polygon south of 85N (search querry limit)
        if lr_lat > 85:
            lr_lat=84
            coverage=0.1

        #write polygon with UL and LR coordinates
        #format [lat1,lat2,lon1,lon2]
        polarstern_poly = [ul_lat,lr_lat,ul_lon,lr_lon]
        print(polarstern_poly)

        #get the date
        start_date=datetime.strftime(dt_noon[i], '%Y-%m-%d')
        end_date=start_date
        
        print(start_date)
        #continue

        output = download_from_polygon(polarstern_poly, start_date, end_date, out_folder, start_hour = '00', end_hour = '23', coverage=coverage)
        if output: 
            fname=output[1]
            print(fname)
            fns.append(fname)
            dts.append(dt_noon[i])

    tt = [dts,fns]
    table = list(zip(*tt))
    outname = ps_file.split('.csv')[0]+'-fnames.csv'
    print(outname)
    with open(outname, 'wb') as f:
            np.savetxt(f, table, fmt="%s", delimiter=",")

    
#alternatives: get mosaics from sentinelhub
