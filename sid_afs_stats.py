import os
from glob import glob
from datetime import datetime
import numpy as np
from shapely.geometry import Point, MultiPoint
from shapely.geometry import Polygon as Shapely_Polygon
from shapely.ops import unary_union
import pickle
from sid_func import * 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

#buffer size in m (max distance between two nods to get connected)
bf = 2000
#alpha for triangulation in concave hull (0.0002 will give max 5km (1/alpha) triangles inside a concave hull)
alpha = 0.0005

inpath = '../sidrift/data/stp10_asf/'
#canberra
#inpath = 'data/stp10_afs/'
outpath = inpath
reg = 'Lance'
metfile = '../sidrift/data/10minute_nounits.csv'

file_name_end = '_60km'

#time series of afs satistics
fig1    = plt.figure(figsize=(14,10))
ax      = fig1.add_subplot(611)
ax.set_title('floe number',fontsize=14, loc='left')
#ax.set_xlabel(r"Length scale (km)",fontsize=25)
#ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
ax.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

bx      = fig1.add_subplot(612)
bx.set_title('floe area',fontsize=14, loc='left')
bx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

cx      = fig1.add_subplot(613)
cx.set_title('floe roundness',fontsize=14, loc='left')
cx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

dx      = fig1.add_subplot(614)
dx.set_title('floe fragmentation',fontsize=14, loc='left')
dx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

ex      = fig1.add_subplot(615)
ex.set_title('LKF area',fontsize=14, loc='left')
ex.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

fx      = fig1.add_subplot(616)
fx.set_title('distance between LKF',fontsize=14, loc='left')
fx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

#scatter plot for floe numbers
fig2    = plt.figure(figsize=(6,6))
aax      = fig2.add_subplot(111)
aax.set_xlabel('Sampling Area Diameter (km)',fontsize=25)
aax.set_ylabel('Number of Floes',fontsize=25)

outname_asf = 'asf_'+reg+file_name_end+'.csv'
rlist = glob(outpath+outname_asf)
for fn in rlist:
    os.remove(fn)

fl = sorted(glob(inpath+'Afs*'+file_name_end+'*.npz'))

for i in fl:
    print(i)
    
    #get date
    date = i.split('_')[2]
    dt = datetime.strptime(date, "%Y%m%dT%H%M%S")

    #load data
    container = np.load(i, allow_pickle=True)
    pindex = container['pindex']
    tripts = container['tripts']
    
    #get all nods of triangles in lkfs
    lkf_tri = [ tripts[p] for p in pindex ]
    lkf_nods = [val for sublist in lkf_tri for val in sublist]
    print('This pair has # LKF nods',len(lkf_nods)) 
    #lkf_nods_extra = lkf_nods.copy()
        
    #a region for all area where there is data (following image pair edges)
    all_nods = [val for sublist in tripts for val in sublist]
    region = MultiPoint(all_nods).convex_hull
    
    #buffer size in meters
    print('buffer size',bf)
    
    #alpha for triangulation in concave hull
    print('max triangulation distance',1/alpha)
    
    #make polygons of all triangles
    tmp = [ Shapely_Polygon(p) for p in lkf_tri ]
    poly_tri = unary_union(tmp)
    
    #buffer around them in hope they will merge
    poly_buff = poly_tri.buffer(bf)
    
    #check what is contained in individual polygons of this multipolygon
    #make concave hull of those nods
    all_lkf=[]
    holes=[]
    if poly_buff.geom_type == 'MultiPolygon':
        for geom in poly_buff.geoms:
            poly_nods=[]
            for j in lkf_nods:
                if geom.contains(Point(j)):
                    poly_nods.append(j)
    
            #concave hull
            if len(poly_nods) > 3:
                try:
                    poly1a, edge_points = alpha_shape(poly_nods, alpha)
                    all_lkf.append(poly1a)
                
                    #check if this polygon has a hole inside
                    #make negative buffer polygon
                    bs = poly1a.buffer(-1*bf)
                    #make difference between whole polygon and negative  buffer
                    if bs.area > 0:
                        hole = poly1a.intersection(bs)
                        if hole.geom_type == 'MultiPolygon':
                            for geom in hole.geoms: 
                                #also check if there is LKFs inside these holes
                                whole = geom.difference(poly_buff)
                                if whole.geom_type == 'MultiPolygon':
                                    for geom in whole.geoms:
                                        #if there is significant surface left, attach this difference polygon as new floe
                                        if geom.area > 1e6:
                                            print('found hole with area:',geom.area)
                                            holes.append(geom.buffer(bf))          #get this buffer area back
                                else:    
                                    if whole.area > 1e6:
                                            print('found hole with area:',whole.area)
                                            holes.append(whole.buffer(bf))
                        else:
                            whole = hole.difference(poly_buff)
                            if whole.area > 1e7:
                                print('found hole with area:',whole.area)
                                holes.append(whole.buffer(bf))

                except:
                        print('problematic hull, likely too small anyway')  
    
    else:
        continue
    
    poly_lkf = unary_union(all_lkf)
    
    #leads wider then 1/alpha are lost by hulling, add them back!
    poly_lkf = poly_lkf.union(poly_tri)
    
    #how to connect some more LKFs
    #calculate area of each polygon and select only the small ones 
    #calculate distance from centroid to all vertex, select top 10% vertex
    #make buffers around these selected vertexes, one by one
    #check if there is another polygon in that buffer (small or big)
    #if yes: connect the vertex to closest point in that buffer, draw a line, make buffer around that line, make that into a polygon, unify all 3 polygons

    
    #add a small frame (500m) along the edges of the region
    #this will close any floes that run accros the region edge
    frame = region.boundary.buffer(500)
    poly4 = poly_lkf.union(frame)

    #get difference of both = floes!
    floes = region.difference(poly4)
    
    #add holes (as floes)
    holes = unary_union(holes)
    floes = floes.union(holes)
    
    #save this polygons
    with open(outpath+'afs_poly_'+date+'_'+file_name_end, "wb") as poly_file:
        pickle.dump(floes, poly_file, pickle.HIGHEST_PROTOCOL)
        
    #simplify polygons
    #break to lines
    #calculate angles between them
    
    
    #Plotting
    fig6    = plt.figure(figsize=(10,10))
    px      = fig6.add_subplot(111)
    
    #map area definition
    area_def_file = glob(inpath+'area_def_'+date+'*'+file_name_end+'*.pkl')[0]
    with open(area_def_file, "rb") as pkl:
        area_def = pickle.load(pkl)
    m = pr.plot.area_def2basemap(area_def)
    
    m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
                    
    #plot the union polygons (typically multipolygons)
    if holes.geom_type == 'Polygon':
        xg, yg = holes.exterior.xy 
        px.fill(xg, yg, alpha=1, fc='none', ec='b')
    if holes.geom_type == 'MultiPolygon':
        for geom in holes.geoms:  
            xg, yg = geom.exterior.xy    
            px.fill(xg, yg, alpha=1, fc='none', ec='b')
    
    #plot the union polygons (typically multipolygons)
    if poly_lkf.geom_type == 'Polygon':
        xg, yg = poly_lkf.exterior.xy 
        px.fill(xg, yg, alpha=0.5, fc='none', ec='r')
    if poly_lkf.geom_type == 'MultiPolygon':
        for geom in poly_lkf.geoms:  
            xg, yg = geom.exterior.xy    
            px.fill(xg, yg, alpha=0.5, fc='none', ec='r')
    
    xg, yg = region.exterior.xy 
    px.fill(xg, yg, alpha=0.5, fc='none', ec='g')
        
    #plot difference
    if floes.geom_type == 'Polygon':
        xg, yg = floes.exterior.xy 
        px.fill(xg, yg, alpha=0.5)
    if floes.geom_type == 'MultiPolygon':
        for geom in floes.geoms:  
            xg, yg = geom.exterior.xy    
            px.fill(xg, yg, alpha=0.5)

    #plot LKF triangles over
    patches_p = []
    for k in pindex:
        patch = Polygon(tripts[k])
        patches_p.append(patch)
        
    p = PatchCollection(patches_p, ec= 'g', fc=None, alpha=1)
    px.add_collection(p)

    fig6.savefig(outpath+'test_afs_'+str(int(bf/1000))+'km'+date.split('T')[0]+'_'+file_name_end,bbox_inches='tight')
    
    print('Done with ASF analysis for', date.split('T')[0])
    #exit()
    
    #with open(i, "rb") as poly_file:
        #poly = pickle.load(poly_file)

    #get some stats
    asf_date = []
    asf_num = []
    asf_area = []
    asf_rr = []
    asf_fr = []
    asf_lkfa = []
    asf_lkfd_min = []
    asf_lkfd_max = []

    if floes.geom_type == 'MultiPolygon':
        
        #fraction of area not covered by floes (LKF fraction)
        a_lkf=0
        if poly_tri.geom_type == 'MultiPolygon':
            for geom in poly_tri.geoms:
                #polygon area
                a_lkf = a_lkf  + geom.area
        if poly_tri.geom_type == 'Polygon':
            a_lkf = geom.area

        #fore every floe
        for geom in floes.geoms:
            #polygon area
            a = geom.area

            #shortest radius
            centroid = geom.centroid
            boundary = geom.boundary
            rs = centroid.distance(boundary)
            #longest radius
            #The Hausdorff distance between two geometries is the furthest distance that a point on either geometry can be from the nearest point to it on the other geometry.
            rl = boundary.hausdorff_distance(centroid)
            #radius ratio ('roundness')
            rr = rs/rl
                        
            #fragmentation ratio (will be high if this floe is actually a conglomerate that can not be properly separated by LKF)
            #perimeter/radius ratio
            bl = geom.boundary.length
            fr = bl/rs
                        
            #save data in lists
            asf_area.append(a)
            asf_rr.append(rr)
            asf_fr.append(fr)
            asf_lkfa.append(a_lkf)
            
            asf_date.append(dt)
            asf_num.append(len(floes))
            
            #min and max distance between LKFs
            asf_lkfd_min.append(rs*2)
            asf_lkfd_max.append(rl*2)
        
        #print(asf_date,asf_num,asf_area,asf_rr,asf_fr)
        
        #write individual floe stats
        tt = [asf_date,asf_num,asf_area,asf_rr,asf_fr,asf_lkfa,asf_lkfd_min,asf_lkfd_max]
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
        ex.scatter(asf_date,asf_lkfa)
        fx.scatter(asf_date,asf_lkfd_min)
        fx.scatter(asf_date,asf_lkfd_max)

    else:
        print('Just a single polygon')
        
    #get radius content
    
    #radius area (area over which we observe deformaiton)
    ra = np.arange(5000,60000,2000) #5km increase from 7km to 60km 
    
    n=9 # number of samples
    ra=np.exp(np.linspace(np.log(5000),np.log(60000),n))
    ra = ra.astype(int)

    print(ra)
    #exit()
    
    for i in ra:
        print(i)
        
        #center point of the image/polygon area
        #xl, yl = m(Lance_lon, Lance_lat)
        #px.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
        #Lance postion (from Lance's met system)
        mettime = getColumn(metfile,0)
        dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
        mi = np.argmin(abs(np.asarray(dtb)-dt))
        Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
        Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]
        xl, yl = m(Lance_lon, Lance_lat)
        
        #make buffer (or square) around Lance
        lance=Point(xl,yl)
        r_region=lance.buffer(i)
        
        #count how many floes we have inside this buffer
        #it counts if they are partly inside
        if floes.geom_type == 'MultiPolygon':
        
            floe_counter=0
            
            #fore every floe
            for geom in floes.geoms:
                if r_region.intersects(geom):
                    floe_counter=floe_counter+1
        
        else:
            floe_counter=1
            
        print(floe_counter)
        
        #plot
        aax.scatter(i*2/1000,floe_counter,alpha=.5,color='royalblue')
        
        #store to make mean value for each radius
    
    #exit()

#save afs time series
fig1.tight_layout()
fig1.savefig(outpath+'asf_ts_'+file_name_end,bbox_inches='tight')

#save radius scatter plots
fig2.tight_layout()
fig2.savefig(outpath+'asf_ra_'+file_name_end,bbox_inches='tight')
