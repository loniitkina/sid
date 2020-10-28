import os
from glob import glob
from datetime import datetime
import numpy as np
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from shapely.geometry import Point, MultiPoint
from shapely.geometry import Polygon as Shapely_Polygon
from shapely.ops import unary_union
import pickle
from sid_func import * 
import itertools
import matplotlib.pyplot as plt

#buffer size in m (max distance between two nods to get connected)
#6 and 7 km join a critical connection in the main 'N', but also too much of spaces between LKFs are filled
#2km keeps most fo the LKF bits separate
bf = 3000
#alpha for triangulation in concave hull (0.0002 will give max 5km (1/alpha) triangles inside a concave hull)
alpha = 0.0003


inpath = '../sidrift/data/stp10_asf/'
#canberra
inpath = 'data/stp10_afs/'
outpath = inpath
reg = 'Lance'

file_name_end = '_20km'

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

outname_asf = 'asf_'+reg+file_name_end+'.csv'
rlist = glob(outpath+outname_asf)
for fn in rlist:
    os.remove(fn)



#fl = sorted(glob(inpath+'*poly*'))

fl = sorted(glob(inpath+'Afs*'+file_name_end+'*.npz'))
fl = sorted(glob(inpath+'Afs*.npz'))
print(fl)

for i in fl:
    print(i)
    
    #load data
    container = np.load(i)
    pindex = container['pindex']
    tripts = container['tripts']
    xlist = container['xlist']
    ylist = container['ylist']

    #get all nods of triangles in lkfs
    lkf_tri = [ tripts[p] for p in pindex ]
    lkf_nods = [val for sublist in lkf_tri for val in sublist]
    print('This pair has # LKF nods',len(lkf_nods)) 
    #lkf_nods_extra = lkf_nods.copy()
    
    region = Shapely_Polygon([[xlist[0],ylist[0]],[xlist[1],ylist[1]],[xlist[2],ylist[2]],[xlist[3],ylist[3]]])
    
    #buffer size in meters
    print('buffer size',bf)
    
    #alpha for triangulation in concave hull
    print('triangle radius up to',1/alpha)
    
    #for every nod draw a buffer around of size bf
    all_poly2=[]
    all_poly3=[]
    for i in lkf_nods:
        point = Point(i)
        circle = Shapely_Polygon(point.buffer(bf).exterior)
        
        #check if any other points are inside this buffer
        close_nods=[i]
        for j in lkf_nods:
            if circle.contains(Point(j)):
                close_nods.append(j)
        
        #convex hull
        #poly2 = MultiPoint(close_nods).convex_hull
        
        #concave_hull should be better as it wont fill in bends in lkf shape
        #it will not work if there are too few points
        if len(close_nods) < 3:
            continue    #just ignore these cases - not sure how this is possible... huge triangles?
        else:
            #it also wont work if points are few and co-planar
            try:
                poly3, edge_points = alpha_shape(close_nods, alpha)
            except:
                continue    #also ignore those cases
        
        #collect all such polygons
        #all_poly2.append(poly2)
        all_poly3.append(poly3)
    
    #make union
    #poly2 = unary_union(all_poly2)
    poly3 = unary_union(all_poly3)
    
    #can we thin such polygons to lines (erosion, simplify)
    #extend them by some km
    #check if we meet another line
    #marge polygons
    
    #add a small frame (500m) along the edges of the region
    #this will close any floes that run accros the region edge
    frame = region.boundary.buffer(500)
    poly4 = poly3.union(frame)

    
    #get difference of both
    poly5 = region.difference(poly4)
    
    #save this polygons
    with open(outpath_def+'afs_poly_'+date1+'_'+file_name_end, "wb") as poly_file:
        pickle.dump(poly5, poly_file, pickle.HIGHEST_PROTOCOL)
        #np.save(poly_file,poly5)


    #plot filtered/remaining triangles
    #check if they are dis/connected
    fig6    = plt.figure(figsize=(10,10))
    px      = fig6.add_subplot(111)
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
    
    #Lance
    xl, yl = m(Lance_lon, Lance_lat)
    px.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

    #patches = []
    #for k in range(tripts_afs.shape[0]):
        #patch = Polygon(tripts_afs[k,:,:])
        #patches.append(patch)

    ##plot filled triangles
    #p = PatchCollection(patches, cmap=plt.cm.bwr, alpha=1)
    #p.set_array(floe_id)
    ##interval = [0, 2]
    #p.set_clim(interval)
    #px.add_collection(p)
    
    #plot the union polygons (typically multipolygons)
    if poly3.geom_type == 'Polygon':
        xg, yg = poly3.exterior.xy 
        px.fill(xg, yg, alpha=0.5, fc='none', ec='r')
    if poly3.geom_type == 'MultiPolygon':
        for geom in poly3.geoms:  
            xg, yg = geom.exterior.xy    
            px.fill(xg, yg, alpha=0.5, fc='none', ec='r')
    
    xg, yg = region.exterior.xy 
    px.fill(xg, yg, alpha=0.5, fc='none', ec='b')
    
    xg, yg = frame.exterior.xy 
    px.fill(xg, yg, alpha=0.5, fc='none', ec='g')
    
    #plot difference
    if poly5.geom_type == 'Polygon':
        xg, yg = poly5.exterior.xy 
        px.fill(xg, yg, alpha=0.5)
    if poly5.geom_type == 'MultiPolygon':
        for geom in poly5.geoms:  
            xg, yg = geom.exterior.xy    
            px.fill(xg, yg, alpha=0.5)

    #plot LKF triangles over
    patches_p = []
    for k in pindex:
        patch = Polygon(tripts[k])
        patches_p.append(patch)
        
    p = PatchCollection(patches_p, ec= 'g', fc=None, alpha=1)
    px.add_collection(p)

    fig6.savefig(outpath+'test_afs_'+str(int(bf/1000))+'km'+date1.split('T')[0]+'_'+str(int(radius/1000))+'km',bbox_inches='tight')
    
    print('Done with ASF analysis for', date1.split('T')[0])
    
    #loop through the time series and exit when done
    continue    

    
    
    
    
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
fig1.savefig(outpath+'asf_ts_'+file_name_end,bbox_inches='tight')

