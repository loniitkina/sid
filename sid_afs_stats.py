import os
from glob import glob
from datetime import datetime
import numpy as np
from shapely.geometry import Point, MultiPoint
from shapely.geometry import Polygon as Shapely_Polygon
import pickle
from sid_func import * 
import itertools
import matplotlib.pyplot as plt


inpath = '../sidrift/data/stp10_asf/'
outpath = inpath
reg = 'Lance'

file_name_end = '_25km_afs.cvs'

#time series of afs satistics
fig1    = plt.figure(figsize=(12,10))
ax      = fig1.add_subplot(411)
ax.set_title('floe number',fontsize=14, loc='left')
#ax.set_xlabel(r"Length scale (km)",fontsize=25)
#ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
ax.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 16))

bx      = fig1.add_subplot(412)
bx.set_title('floe area',fontsize=14, loc='left')
bx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 16))

cx      = fig1.add_subplot(413)
cx.set_title('floe roundness',fontsize=14, loc='left')
cx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 16))

dx      = fig1.add_subplot(414)
dx.set_title('floe fragmentation',fontsize=14, loc='left')
dx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 16))

outname_asf = 'asf_'+reg+file_name_end
rlist = glob(outpath+outname_asf)
for fn in rlist:
    os.remove(fn)



fl = sorted(glob(inpath+'*poly*'))
print(fl)

for i in fl:
    print(i)
    
    #get date
    date = i.split('_')[-1]
    dt = datetime.strptime(date, "%Y%m%dT%H%M%S")

    with open(i, "rb") as poly_file:
        poly = pickle.load(poly_file)
        #poly = np.load(poly_file, allow_pickle=True)

    #get some stats
    asf_date = []
    asf_num = []
    asf_area=[]
    asf_rr=[]
    asf_fr=[]

    if poly.geom_type == 'MultiPolygon':
        for geom in poly.geoms:
            #polygon area
            
            a = geom.area
            #print('area',a)
            #shortest radius
            centroid = geom.centroid
            rs = centroid.distance(geom.boundary)
            #print(rs)
            #longest radius
            #The Hausdorff distance between two geometries is the furthest distance that a point on either geometry can be from the nearest point to it on the other geometry.
            rl = centroid.hausdorff_distance(geom.boundary)
            #print(rl)
            #radius ratio ('roundness')
            rr = rs/rl
            #print('roundness',rr)
            
            #fragmentation ratio (fractuaction of the floe shape will be high if this floe is actually a conglomerate that can not be properly separated)
            #perimeter/area ratio
            bl = geom.boundary.length
            fr = bl/a
            #print('fragmentation',fr)
            
            #does this floe include smaller polygons? (fragmets of boundary)
            
            #fraction of area not covered by floes (LKF fraction)
            asf_area.append(a)
            asf_rr.append(rr)
            asf_fr.append(fr)
            
            asf_date.append(dt)
            asf_num.append(len(poly))
        
        #print(asf_date,asf_num,asf_area,asf_rr,asf_fr)
        
        #write individual floe stats
        tt = [asf_date,asf_num,asf_area,asf_rr,asf_fr]
        table = zip(*tt)
        #adjusted to python3:
        table = list(zip(*tt))

        output = outpath + outname_asf
        with open(output, 'ab') as f:
            np.savetxt(f, table, fmt="%s", delimiter=",")

        #time series of afs satistics
        ax.scatter(asf_date,asf_num)
        bx.scatter(asf_date,asf_area)
        cx.scatter(asf_date,asf_rr)
        dx.scatter(asf_date,asf_fr)

    else:
        print('Just a single polygon')

#save afs time series
fig1.tight_layout()
fig1.savefig(outpath+'asf_ts_25km',bbox_inches='tight')

