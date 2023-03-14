import os
from glob import glob
from datetime import datetime
import numpy as np
from shapely.geometry import Point, MultiPoint, MultiPolygon
from shapely.geometry import Polygon as Shapely_Polygon
from shapely.ops import unary_union
import pickle
from sid_func import * 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable


#coherent dynamics features/element
#floes - summer idea of separate floes might confuse people

#buffer size in m (max distance between two nods to get connected)
#there is a trade-off between connecting ends of two LKFs and across a breath of the two intersecting LKFs
bf = 8000
#alpha for triangulation in concave hull/filling in small floes (0.0002 will give max 5km (1/alpha) triangulation distance inside a concave hull)
#small distance prevents 'micro floes'
alpha = 0.0005      #max 2km
#frame buffer - half of max distance between two LKF polygons to be connected
fbf = 1000
#additional buffer that increases connectivity of LKF/decreases floe area
margin_bf = 1700
#minimal size of the hole to be recognized as floe - the larger the 'bf', smaller holes are worth marking
#min_hole_area = 1e9
min_hole_area = 1e7     #good for 6km buffer
min_hole_area = 1e6 

inpath = '/scratch/pit000/results/sid/afs/'
inpath = '/scratch/pit000/results/sid/parcels/'
inpath = '/scratch/pit000/results/sid/deform200km/'
outpath = inpath
plotpath = '/scratch/pit000/results/sid/plots200km/'
#leg1
#shipfile = '../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0.csv'
#leg2
#shipfile = '../../downloads/data_master-solution_mosaic-leg2-20191214-20200224-floenavi-refstat-v1p0.csv'
#leg3 (and transition to leg 4 until 6 June)
#shipfile = '../../downloads/position_leg3_nh-track.csv'

shipfile = '../../downloads/MOSAiC_all_legs.csv'


reg = 'ship'
#file_name_date = '2019'	#leg1 and leg2
#file_name_date = '2020'	#leg 2 and leg3
file_name_date = '2019-2020'    #all legs
file_name_end = '_120km'

#time series of afs satistics
fig1    = plt.figure(figsize=(14,12))
ax      = fig1.add_subplot(611)
ax.set_title('Floe count',fontsize=20, loc='left')
#ax.set_xlabel(r"Length scale (km)",fontsize=25)
#ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
#ax.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

bx      = fig1.add_subplot(612)
bx.set_title('Floe area',fontsize=20, loc='left')
bx.set_ylabel('$km^2$',fontsize=20)
bx.set_ylim(0, 10000)
#bx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

cx      = fig1.add_subplot(613)
cx.set_title('Floe roundness',fontsize=20, loc='left')
#cx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

dx      = fig1.add_subplot(614)
dx.set_title('Floe fragmentation',fontsize=20, loc='left')
dx.set_ylim(0, 500)
#dx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

ex      = fig1.add_subplot(615)
ex.set_title('LKF area',fontsize=20, loc='left')
ex.set_ylabel('$km^2$',fontsize=20)
#ex.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

fx      = fig1.add_subplot(616)
#min and max floe radius
fx.set_title('Distance between LKFs',fontsize=20, loc='left')
fx.set_ylabel('$km$',fontsize=20)
fx.set_ylim(0, 200)
#fx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

#scatter plot for floe numbers
fig2    = plt.figure(figsize=(6,6))
aax      = fig2.add_subplot(111)
aax.set_xlabel('Sampling Area Diameter (km)',fontsize=25)
aax.set_ylabel('Number of Floes',fontsize=25)

#previously stored stats
outname_afs = 'afs_'+reg+file_name_date+file_name_end+'.csv'
rlist = glob(outpath+outname_afs)
for fn in rlist:
    os.remove(fn)

#prepare lists to store the data
afs_date = []
afs_num = []
afs_lkfa = []
afs_area_m = []
afs_area_s = []
afs_rr_m = []
afs_rr_s = []
afs_fr_m = []
afs_fr_s = []
afs_lkfd_min_m = []
afs_lkfd_min_s = []
afs_lkfd_max_m = []
afs_lkfd_max_s = []

fl = sorted(glob(inpath+'Afs*'+file_name_end+'*_tiled.npz'))#[:-2]
#fl = ['/scratch/pit000/results/sid/deform200km/Afs_20200401T114043_20200402T122140_120km.csv_tiled.npz']
#fl = ['/scratch/pit000/results/sid/deform200km/Afs_20200316T121329_20200317T111605_120km.csv_tiled.npz']
#fl = ['/scratch/pit000/results/sid/deform200km/Afs_20200415T112421_20200416T120518_120km.csv_tiled.npz']
#fl = ['/scratch/pit000/results/sid/deform200km/Afs_20200414T122141_20200415T112421_120km.csv_tiled.npz']

print(fl)


for i in fl:
    print(i)
    
    #get date
    date = i.split('_')[1]
    print(date)
    dt = datetime.strptime(date, "%Y%m%dT%H%M%S")
    
    #load data
    container = np.load(i, allow_pickle=True)
    pindex = container['pindex']
    tripts = container['tripts']
    region = container['region']    #This is typically a MultiPolygon (or sometimes still matplotlib Polygon)
    corner_nods = container['corners']
    
    #get all nods of triangles in lkfs
    lkf_tri = [ tripts[p] for p in pindex ]
    lkf_nods = [val for sublist in lkf_tri for val in sublist]
    print('This pair has # LKF nods',len(lkf_nods)) 
    #lkf_nods_extra = lkf_nods.copy()
        
    ##a region for all area where there is data (following image pair edges)
    ##this does not work well for tiled regions - empty corners are part of the region - read from container instead
    #all_nods = [val for sublist in tripts for val in sublist]
    #region = MultiPoint(all_nods).convex_hull
    
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
                    poly1a = poly1a.buffer(bf+margin_bf)    #add some buffer
                    all_lkf.append(poly1a)
                
                    #check if this polygon has a hole inside
                    #make negative buffer polygon
                    bs = poly1a.buffer(-1*bf)
                    #make difference between whole polygon and negative  buffer
                    if bs.area > min_hole_area:
                        hole = poly1a.intersection(bs)
                        if hole.geom_type == 'MultiPolygon':
                            for geom in hole.geoms: 
                                #also check if there is LKFs inside these holes
                                whole = geom.difference(poly_buff)
                                if whole.geom_type == 'MultiPolygon':
                                    for geom in whole.geoms:
                                        #if there is significant surface left, attach this difference polygon as new floe
                                        if geom.area > min_hole_area:
                                            print('found hole with area:',geom.area)
                                            holes.append(geom.buffer(bf-margin_bf))          #get this buffer area back
                                else:    
                                    if whole.area > min_hole_area:
                                            print('found hole with area:',whole.area)
                                            holes.append(whole.buffer(bf-margin_bf))
                        else:
                            whole = hole.difference(poly_buff)
                            if whole.area > min_hole_area:
                                print('found hole with area:',whole.area)
                                holes.append(whole.buffer(bf-margin_bf))

                except:
                        print('problematic hull, likely too small anyway')  
    
    else:
        continue
    
    ##make additional buffer around LKFs to merge more of them
    #all_lkf = all_lkf.buffer(1000)
    
    poly_lkf = unary_union(all_lkf)
    
    #leads wider then 1/alpha are lost by hulling, add them back!
    #also add a small buffer
    poly_lkf = poly_lkf.union(poly_tri.buffer(fbf*3))
    
    #how to connect some more LKFs
    #calculate area of each polygon and select only the small ones 
    #calculate distance from centroid to all vertex, select top 10% vertex
    #make buffers around these selected vertexes, one by one
    #check if there is another polygon in that buffer (small or big)
    #if yes: connect the vertex to closest point in that buffer, draw a line, make buffer around that line, make that into a polygon, unify all 3 polygons

    #get the region polygon
    #for some reason storing in numpy will preserve Polyongs, but not MultiPolyong types
    #print(region)
    try:
        region=unary_union(region)
        #try to merge two distant tiles (max 10km? appart)
        #add some buffer and see if there is any intersection
        if region.geom_type == 'MultiPolygon':
            tmp = region.buffer(5000)
            region=unary_union(tmp)
            region = region.buffer(-5000)
    except:
        print('matplotlib Polygon')    
        #WARNING:this is a temporary hack - it produces no desired frame!!!
        all_nods = [val for sublist in tripts for val in sublist]
        region = MultiPoint(all_nods).convex_hull 
    
    #the totoal coverage is sometimes much larger than the map we make here/region analysed for AFS - limit it to the picture frame
    plot_corners = Shapely_Polygon(corner_nods)
    region=region.intersection(plot_corners)
    #print(region)
    
    #add a small frame (~500m) along the edges of the region
    #this will close any floes that run accros the region edge
    frame = region.boundary.buffer(fbf)
    poly4 = poly_lkf.union(frame)

    #get difference of both = floes!
    floes1 = region.difference(poly4)
        
    #the LKFs are now too wide (floes are too skinny)
    #we cant increase/fatten the floe inside the union as they will merge
    #so we increase them one by one instead
    floes=[]
    if floes1.geom_type == 'MultiPolygon':
        used_space=floes1
        for geom in floes1.geoms:  
            fat = geom.buffer(bf)   #bloated floe
            
            #erode this polygon by already used space
            used_space = used_space.difference(geom)    #all other floes in the union
            keep = fat.difference(used_space)           #bloated floe without the area already covered by others
            used_space = used_space.union(keep)
            floes.append(keep)
    else:
        floes.append(floes1)
    
    #add holes (as floes)
    holes = unary_union(holes)
    if holes.geom_type == 'Polygon':
        floes.append(holes)
    if holes.geom_type == 'MultiPolygon':
        for geom in holes.geoms:  
            slim = geom.difference(used_space)  #is there any area of this 'hole' already covered now by the fat floe?
            used_space = used_space.union(slim) #add this 'hole' to the used space
            floes.append(slim)
    
    #convert list to MultiPolygon
    floes = MultiPolygon(floes)
    
    #save these polygons
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
    
    if region.geom_type == 'Polygon':
        xg, yg = region.exterior.xy 
        px.fill(xg, yg, alpha=0.5, fc='none', ec='g')
    if region.geom_type == 'MultiPolygon':
        for geom in region.geoms:  
            xg, yg = geom.exterior.xy    
            px.fill(xg, yg, alpha=0.5, fc='none', ec='g')
    
    #xg, yg = region.exterior.xy 
    #px.fill(xg, yg, alpha=0.5, fc='none', ec='g')
        
    #plot floes
    if floes.geom_type == 'Polygon':
        xg, yg = floes.exterior.xy 
        px.fill(xg, yg, alpha=0.5)
    if floes.geom_type == 'MultiPolygon':
        
        print('Floe number: ',len(floes.geoms))
        colors=iter(plt.cm.tab20(np.linspace(0,1,len(floes.geoms))))

        for geom in floes.geoms:  
            #geom = geom.buffer(bf)          #PROBLEM: holes should not be increased...
            xg, yg = geom.exterior.xy    
            cl = next(colors)
            px.fill(xg, yg, alpha=0.5, c=cl)

    #plot LKF triangles over
    patches_p = []
    for k in pindex:
        patch = Polygon(tripts[k])
        patches_p.append(patch)
        
    p = PatchCollection(patches_p, ec= 'k')#, fc=None, alpha=1)
    px.add_collection(p)
    
    #plot the ship
    mettime = getColumn(shipfile,0)
    dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
    mi = np.argmin(abs(np.asarray(dtb)-dt))
    ship_lon = np.asarray(getColumn(shipfile,1),dtype=float)[mi]
    ship_lat = np.asarray(getColumn(shipfile,2),dtype=float)[mi]
    xl, yl = m(ship_lon, ship_lat)
    px.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

    fig6.savefig(plotpath+'test_afs_'+str(int(bf/1000))+'km'+date.split('T')[0]+'_'+file_name_end,bbox_inches='tight')
    
    print('Done with AFS analysis for', date.split('T')[0])
    #exit()
    
    #with open(i, "rb") as poly_file:
        #poly = pickle.load(poly_file)

    #get some stats
    afs_area = []
    afs_rr = []
    afs_fr = []
    afs_lkfd_min = []
    afs_lkfd_max = []

    if floes.geom_type == 'MultiPolygon':
        
        #fraction of area not covered by floes (LKF fraction)
        a_lkf=0
        if poly_tri.geom_type == 'MultiPolygon':
            for geom in poly_tri.geoms:
                #polygon area
                a_lkf = a_lkf  + geom.area
        if poly_tri.geom_type == 'Polygon':
            a_lkf = geom.area

        #for every floe
        for geom in floes.geoms:
            
            #polygon area
            a = geom.area
            a = a /1000000  #convert to km^2

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
            #perimeter/radius
            bl = geom.boundary.length
            fr = bl/rs
                        
            #lists for this image pair
            afs_area.append(a)
            afs_rr.append(rr)
            afs_fr.append(fr)
            
            #min and max distance between LKFs
            afs_lkfd_min.append(rs*2/1000)
            afs_lkfd_max.append(rl*2/1000)

        
        #time series lists
        afs_area_m.append(np.mean(afs_area))
        afs_area_s.append(np.std(afs_area))
        
        afs_rr_m.append(np.mean(afs_rr))
        afs_rr_s.append(np.std(afs_rr))
        
        afs_fr_m.append(np.mean(afs_fr))
        afs_fr_s.append(np.std(afs_fr))
        
        afs_lkfd_min_m.append(np.mean(afs_lkfd_min))
        afs_lkfd_min_s.append(np.std(afs_lkfd_min))
        
        afs_lkfd_max_m.append(np.mean(afs_lkfd_max))
        afs_lkfd_max_s.append(np.std(afs_lkfd_max))
         
        afs_date.append(dt)
        afs_num.append(len(floes))
        afs_lkfa.append(a_lkf/1000000)
        
    else:
        print('Just a single polygon')
        
    #get radius content
    
    #radius area (area over which we observe deformaiton)
    ra = np.arange(5000,60000,2000) #5km increase from 7km to 60km 
    
    n=9 # number of samples
    ra=np.exp(np.linspace(np.log(5000),np.log(60000),n))
    ra = ra.astype(int)

    #print(ra)
    #exit()
    
    for i in ra:
        #print(i)
        
        #center point of the image/polygon area
        mettime = getColumn(shipfile,0)
        dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
        mi = np.argmin(abs(np.asarray(dtb)-dt))
        ship_lon = np.asarray(getColumn(shipfile,1),dtype=float)[mi]
        ship_lat = np.asarray(getColumn(shipfile,2),dtype=float)[mi]
        xl, yl = m(ship_lon, ship_lat)
        
        #make buffer (or square) around ship
        ship=Point(xl,yl)
        r_region=ship.buffer(i)
        
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
            
        #print('Floe counter: ',floe_counter)
        
        #plot
        aax.scatter(i*2/1000,floe_counter,alpha=.5,color='royalblue')
        
        #store to make mean value for each radius
    
    #exit()


#time series of statistics
ax.scatter(afs_date,afs_num)

bx.errorbar(afs_date,afs_area_m,afs_area_s,linestyle='None',marker='o')

cx.errorbar(afs_date,afs_rr_m,afs_rr_s,linestyle='None',marker='o')

dx.errorbar(afs_date,afs_fr_m,afs_fr_s,linestyle='None',marker='o')

ex.scatter(afs_date,afs_lkfa)

fx.errorbar(afs_date,afs_lkfd_min_m,afs_lkfd_min_s,linestyle='None',marker='o')
fx.errorbar(afs_date,afs_lkfd_max_m,afs_lkfd_max_s,linestyle='None',marker='o')


#save afs timeseries
#use the final date in the name
fig1.tight_layout()
fig1.savefig(plotpath+'afs_ts_'+date+'_'+str(int(bf/1000))+'km'+file_name_end,bbox_inches='tight')

#save radius scatter plots
fig2.tight_layout()
fig2.savefig(plotpath+'afs_ra_'+date+'_'+str(int(bf/1000))+'km'+file_name_end,bbox_inches='tight')

#save output for replotting
tt = [afs_date,afs_num,afs_area_m,afs_area_s,afs_rr_m,afs_rr_s,afs_fr_m,afs_fr_s,afs_lkfa,afs_lkfd_min_m,afs_lkfd_min_s,afs_lkfd_max_m,afs_lkfd_max_s]
table = zip(*tt)
table = list(zip(*tt))

output = outpath + outname_afs
print(output)
with open(output, 'ab') as f:
    #header
    f.write(b'date,DE count,DE area_m,DE area_s, DE rr_m,DE rr_s, DE fr_m,DE fr_s,LKF area,LKF d_min_m,LKF d_min_s,LKF d_max_m,LKF d_max_s\n')
    np.savetxt(f, table, fmt="%s", delimiter=",")

