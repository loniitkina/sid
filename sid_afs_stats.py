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
from matplotlib.patches import Polygon, PathPatch
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import AnchoredText

#coherent dynamic features/element
#floes - summer idea of separate floes might confuse people

#buffer size in m (max distance between two nods to get connected)
#there is a trade-off between connecting ends of two LKFs and across a breath of the two intersecting LKFs
#large buffer means rounded edges CDCs and hughe LKFs 
bf = 8000   #MOSAiC
#bf = 12000
#alpha for triangulation in concave hull/filling in small floes (0.0002 will give max 5km (1/alpha) triangulation distance inside a concave hull)
#small distance prevents 'micro floes'
alpha = 0.0005      #max 2km (mosaic)
#alpha = 0.00025     #max 4km
#frame buffer
fbf = 1000
#additional buffer that increases connectivity of LKF/decreases floe area
#margin_bf = 1700    #MOSAiC
#margin_bf = 2500
margin_bf = 3400
margin_bf = 4000; bf=margin_bf
#minimal size of the hole to be recognized as floe - the larger the 'bf', smaller holes are worth marking
#min_floe_area = 1e9
min_floe_area = 1e7     #good for 6km buffer
#min_floe_area = 1e6     #now only used to exclude tiny floes from statistics
lkf_max_tri_area = 800*800*6 #12 normal triangles, to exclude large triangles in the MIZ from the LKF total area

inpath = '/scratch/pit000/results/sid/cde/'
inpath = '/scratch/pit000/results/sid/parcels/'
inpath = '/scratch/pit000/results/sid/deform200km/'
outpath = inpath
plotpath = '/scratch/pit000/results/sid/plots200km/'
#plotpath = '/scratch/pit000/results/sid/plots_mosaic/'
plotpath = '/scratch/pit000/results/sid/plots200km_revision/'

#MOSAiC
#leg1
#shipfile = '../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0.csv'
#leg2
#shipfile = '../../downloads/data_master-solution_mosaic-leg2-20191214-20200224-floenavi-refstat-v1p0.csv'
#leg3 (and transition to leg 4 until 6 June)
#shipfile = '../../downloads/position_leg3_nh-track.csv'
#shipfile = '../../downloads/MOSAiC_all_legs.csv'
#year='2020'

#N-ICE
#WARNING: This expedition was close to the ice edge - recommended radius for CDC is ~100km.
shipfile = '../../downloads/lance_leg1_c_200km.csv' #original file is time,lat,lon, we need time, lon, lat!!!
year='2015'

##CIRFA
#shipfile = '../../downloads/CIRFA_cruise_stationM.csv'
#year='2022'

reg = 'ship'
#file_name_date = '2019'	#leg1 and leg2
#file_name_date = '2020'	#leg 2 and leg3
#file_name_date = '2019-2020'    #all legs
file_name_date = '2015' #N-ICE
#file_name_date = '2022' #CIRFA
#file_name_end = '_120km'
#file_name_end = '_200km'

#naming and input
radius = 200000
rname = '_'+str(int(radius/1000))+'km'
tname='_thfilter'
lname='_lkffilter'
file_name_end = rname+tname+lname
print('Your input is: ',file_name_end)

#time series of cde satistics
fig1    = plt.figure(figsize=(14,12))
ax      = fig1.add_subplot(611)
ax.set_title('CDC density',fontsize=20, loc='left')
ax.set_ylabel("$km^{-1}$",fontsize=20)
#ax.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

bx      = fig1.add_subplot(612)
bx.set_title('CDC area',fontsize=20, loc='left')
bx.set_ylabel('$km^2$',fontsize=20)
#bx.set_ylim(0, 50000)
#bx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

cx      = fig1.add_subplot(613)
cx.set_title('CDC roundness',fontsize=20, loc='left')
#cx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

dx      = fig1.add_subplot(614)
dx.set_title('CDC fragmentation',fontsize=20, loc='left')
#dx.set_ylim(0, 500)
#dx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

ex      = fig1.add_subplot(615)
ex.set_title('LKF fraction',fontsize=20, loc='left')
ex.set_ylabel('$\%$',fontsize=20)
#ex.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

fx      = fig1.add_subplot(616)
#min and max floe radius
fx.set_title('Distance between LKFs',fontsize=20, loc='left')
fx.set_ylabel('$km$',fontsize=20)
#fx.set_ylim(0, 200)
#fx.set_xlim(datetime(2015, 1, 21), datetime(2015, 2, 9))

#scatter plot for floe numbers
fig2    = plt.figure(figsize=(6,6))
aax      = fig2.add_subplot(111)
aax.set_xlabel('Sampling Area Diameter (km)',fontsize=25)
aax.set_ylabel('Number of Floes',fontsize=25)

#previously stored stats
outname_cde = 'cde_'+str(int(bf/1000))+'km_'+reg+file_name_date+file_name_end+'.csv'
rlist = glob(outpath+outname_cde)
for fn in rlist:
    os.remove(fn)

#prepare lists to store the data (_m is mean, _s is standard deviaition)
cde_date = []
cde_num = []
cde_lkfa = []
cde_lkfa_max = []
cde_area_m = []
cde_area_s = []
cde_rr_m = []
cde_rr_s = []
cde_fr_m = []
cde_fr_s = []
cde_cc_m = []
cde_cc_s = []
cde_lkfd_min_m = []
cde_lkfd_min_s = []
cde_lkfd_max_m = []
cde_lkfd_max_s = []

fl = sorted(glob(inpath+'Afs_'+year+'*'+file_name_end+'_tiled.npz'))    #this option does to grab data with 'stp1' suffix
#fl = sorted(glob(inpath+'Afs_'+year+'*'+file_name_end+'*_tiled.npz'))[:4] 
#fl = sorted(glob(inpath+'Afs_'+year+'0125*'+file_name_end+'*_tiled.npz'))

print(fl)

color=iter(plt.cm.jet_r(np.linspace(0,1,len(fl)+1)))

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
    #print('max triangulation distance',1/alpha)
    
    #make polygons of all triangles
    tmp = [ Shapely_Polygon(p) for p in lkf_tri ]
    poly_tri = unary_union(tmp)
    
    ###buffer around them in hope they will merge
    #tmp = MultiPolygon(tmp)
    #poly_buff = tmp.buffer(bf)
    
    ##WARNING: this part is actually never used!!!
    ##check what is contained in individual polygons of this multipolygon
    ##make concave hull of those nods
    #all_lkf=[]
    #all_hulls=[]
    #holes=[]
    #if poly_buff.geom_type == 'MultiPolygon':
        #for geom in poly_buff.geoms:
            #poly_nods=[]
            #for j in lkf_nods:
                #if geom.contains(Point(j)):
                    #poly_nods.append(j)
    
            ##concave hull
            #if len(poly_nods) > 3:
                #try:
                    #poly1, edge_points = alpha_shape(poly_nods, alpha)
                    #all_hulls.append(poly1)
                    #print(poly1)
                    
                    #poly2 = poly1.buffer(bf+margin_bf)    #add some buffer
                    
                    #print(poly2)
                    #exit()
                    #all_lkf.append(poly2)
                
                    ###WARNING: this part seems to be redundant
                    ###check if this polygon has a hole inside
                    ###make negative buffer polygon
                    ##bs = poly2.buffer(-1*bf)
                    ###make difference between whole polygon and negative  buffer
                    ##if bs.area > min_floe_area:
                        ##hole = poly2.intersection(bs)
                        ##if hole.geom_type == 'MultiPolygon':
                            ##for geom in hole.geoms: 
                                ###also check if there is LKFs inside these holes
                                ##whole = geom.difference(poly_buff)
                                ##if whole.geom_type == 'MultiPolygon':
                                    ##for geom in whole.geoms:
                                        ###if there is significant surface left, attach this difference polygon as new floe
                                        ##if geom.area > min_floe_area:
                                            ##print('found hole with area:',geom.area)
                                            ##holes.append(geom.buffer(bf-margin_bf))          #get this buffer area back
                                ##else:    
                                    ##if whole.area > min_floe_area:
                                            ##print('found hole with area:',whole.area)
                                            ##holes.append(whole.buffer(bf-margin_bf))
                        ##else:
                            ##whole = hole.difference(poly_buff)
                            ##if whole.area > min_floe_area:
                                ##print('found hole with area:',whole.area)
                                ##holes.append(whole.buffer(bf-margin_bf))

                #except:
                        #continue
                        #print('problematic hull, likely too small anyway')  
    
    #else:
        #continue
    
    ##make additional buffer around LKFs to merge more of them
    #all_lkf = all_lkf.buffer(1000)
    
    #poly_hull = unary_union(all_hulls)
    #poly_lkf = unary_union(all_lkf)
    
    #leads wider then 1/alpha are lost by hulling, add them back, and include a buffer
    #poly_lkf = poly_lkf.union(poly_tri.buffer(fbf*3))
    #poly_lkf = poly_lkf.union(poly_tri.buffer(margin_bf))
    poly_buff = poly_tri.buffer(margin_bf)
    poly_lkf = unary_union(poly_buff)
    
    #how to connect some more LKFs
    #calculate area of each polygon and select only the small ones 
    #calculate distance from centroid to all vertex, select top 10% vertex
    #make buffers around these selected vertexes, one by one
    #check if there is another polygon in that buffer (small or big)
    #if yes: connect the vertex to closest point in that buffer, draw a line, make buffer around that line, make that into a polygon, unify all 3 polygons

    #get the region polygon
    region = Shapely_Polygon(region)

    #If N-ICE avoid the MIZ 
    if year=='2015':
        #replace the corner nodes with manual coordinates
        #limit the data in image coordinates
        ##square around Lance
        #minx = 100000; maxx = 300000
        #miny = 100000; maxy = 300000
        ##whole area but the MIZ (south of Lance)
        minx = 0; maxx = 400000
        miny = 100000; maxy = 400000
        corner_nods = np.array([[minx,miny],
                    [minx,maxy],
                    [maxx,maxy],
                    [maxx,miny],
                    [minx,miny]])

    else:
        minx = 0; maxx = minx+2*radius
        miny = 0; maxy = miny+2*radius

    ##for some reason storing in numpy will preserve Polygons, but not MultiPolyon types
    ##print(region)
    #try:
        #region=unary_union(region)
        ##try to merge two distant tiles (max 20km? appart)
        ##add some buffer and see if there is any intersection
        #if region.geom_type == 'MultiPolygon':
            #tmp = region.buffer(20000)
            #region=unary_union(tmp)
            #region = region.buffer(-20000)
    #except:
        #print('matplotlib Polygon')    
        ##WARNING:this is a temporary hack - it produces no desired frame!!!
        #all_nods = [val for sublist in tripts for val in sublist]
        #region = MultiPoint(all_nods).convex_hull 
    
    #the total coverage is sometimes much larger than the map we make here/region analysed for AFS - limit it to the picture frame
    plot_corners = Shapely_Polygon(corner_nods)
    region=region.intersection(plot_corners)
    
    #add a small frame (~500m) along the edges of the region
    #this will close any floes that run accros the region edge
    frame = region.boundary.buffer(fbf)
    poly_lkf_union = poly_lkf.union(frame)
    #prepare the frame to reduce the spill over the region boundary later
    frame2 = region.boundary.buffer(bf)
    frame2 = frame2.difference(region)

    #region minus LKFs = floes!
    floes_tmp = region.difference(poly_lkf_union)
    
    #WARNING: this part seems to be redundant
    #the LKFs are now too wide (floes are too skinny)
    #we cant increase/fatten the floe inside the union as they will merge
    #so we increase them one by one instead
    floes=[]
    if floes_tmp.geom_type == 'MultiPolygon':
        used_space=unary_union(floes_tmp)
        for geom in floes_tmp.geoms:  
            fat = geom.buffer(margin_bf)   #bloated floe
            
            #floes.append(fat)
            
            #but dont get this over the region boundary
            fat = fat.difference(frame2)
            
            #erode this polygon by already used space
            used_space = used_space.difference(geom)    #all other floes in the union
            keep = fat.difference(used_space)           #bloated floe without the area already covered by others
            used_space = used_space.union(keep)
            
            if keep.area < min_floe_area:
                print('tiny floe - do not keep')
            else:
                if keep.geom_type == 'MultiPolygon':
                    for geom in keep.geoms:
                        #floes.append(geom)  #used just for concept figure
                        floes.append(Shapely_Polygon(geom.exterior.coords)) #by this we only keep the exterior part of the polygon and remove any holes
                else:    
                    #floes.append(keep)  #used just for concept figure
                    floes.append(Shapely_Polygon(keep.exterior.coords)) #by this we only keep the exterior part of the polygon and remove any holes
                    
    else:
        floes.append(floes_tmp)
    
    ##WARNING: this part seems to be redundant
    ##add holes (as floes)
    #holes = unary_union(holes)
    #if holes.geom_type == 'Polygon':
        #floes.append(holes)
    #if holes.geom_type == 'MultiPolygon':
        #for geom in holes.geoms:  
            #slim = geom.difference(used_space)  #is there any area of this 'hole' already covered now by the fat floe?
            #used_space = used_space.union(slim) #add this 'hole' to the used space
            #if slim.geom_type == 'MultiPolygon':
                #for geom in slim.geoms:
                    #floes.append(geom)
            #else:    
                #floes.append(slim)
    
    #floes_tmp=[]
    #for floe in floes:
        #floe = floe.difference(frame)
        #floes_tmp.append(floe)
    
    #convert list to MultiPolygon
    #print(floes)
    floes = MultiPolygon(floes)

    
    #ARE THERE MULTIPLE POLYGONS OVER THE SAME SPACE?
    #lkfs = region.symmetric_difference_all(floes)
    ufloes = unary_union(floes)
    lkfs = region.difference(ufloes)
    
    #poly_buffer = floes.buffer(-margin_bf)
    
    #poly_buffer = region.difference(floes)
    
    ###poly_buffer = region.difference(used_space)
    #ufloes = unary_union(floes)
    #poly_buffer = region.difference(ufloes)
    #poly_buffer = poly_buffer.intersection(floes)
    
    #lkfs=[]
    #if floes.geom_type == 'MultiPolygon':
        
        #for geom in floes.geoms:  
            #tmp = region.difference(geom)
            
            #if tmp.geom_type == 'MultiPolygon':
                #for geom in tmp.geoms:
                    #lkfs.append(geom)
            #else:
                #lkfs.append(tmp)
            

    #else:
        #lkfs.append(region.difference(floes))
    
    #poly_buffer = MultiPolygon(lkfs)
    
    #ufloes = unary_union(floes)
    #poly_buffer = region.difference(ufloes)
    #poly_buffer = poly_buffer.intersection(floes)
    
    #save these polygons
    with open(outpath+'cde_poly_'+date+'_'+file_name_end, "wb") as poly_file:
        pickle.dump(floes, poly_file, pickle.HIGHEST_PROTOCOL)
        
    #simplify polygons
    #break to lines
    #calculate angles between them
    
    #######################################################################################################################################################
    
    #get some stats
    cde_area = []
    cde_rr = []
    cde_fr = []
    cde_cc = []
    cde_lkfd_min = []
    cde_lkfd_max = []

    if floes.geom_type == 'MultiPolygon':
        
        #total area
        ta = (maxx-minx)*(maxy-miny)
        
        #fraction of area not covered by floes (LKF fraction)
        a_lkf=0
        if poly_tri.geom_type == 'MultiPolygon':
            for geom in poly_tri.geoms:
                #polygon area
                #do not allow any huge polygons (MIZ and artifacts)
                if geom.area < lkf_max_tri_area:
                    a_lkf = a_lkf  + geom.area
        if poly_tri.geom_type == 'Polygon':
            if geom.area < lkf_max_tri_area:
                a_lkf = geom.area
        #LKF fraction in %
        a_lkf = (a_lkf/ta)*100
        
        #extended fraction - everything that is not included into CDCs
        #this includes also undefined areas, artefacts, MIZ and large triangles
        a_lkf_max=0
        if lkfs.geom_type == 'MultiPolygon':
            for geom in lkfs.geoms:
                #polygon area
                a_lkf_max = a_lkf_max  + geom.area
        if lkfs.geom_type == 'Polygon':
            a_lkf_max = geom.area
        #LKF fraction in %
        a_lkf_max = (a_lkf_max/ta)*100

        #for every floe
        for geom in floes.geoms:
            
            #polygon area
            am = geom.area
            if am < min_floe_area:
                print('tiny floe');continue
            a = am /1000000  #convert to km^2

            #shortest radius
            centroid = geom.centroid
            boundary = geom.boundary
            rs = centroid.distance(boundary)
            #longest radius
            #The Hausdorff distance between two geometries is the furthest distance that a point on either geometry can be from the nearest point to it on the other geometry.
            rl = boundary.hausdorff_distance(centroid)
            #radius ratio ('roundness')
            rr = rs/rl
                        
            #fragmentation/complexity ratio (will be high if this floe is actually a conglomerate that can not be properly separated by LKF)
            #perimeter/diameter
            #The lowest values should be close to 2 in case of a 'stick'
            bl = geom.boundary.length
            diameter = np.mean([rs,rl])*2
            fr = bl/diameter
            
            #circularity
            hull = geom.convex_hull
            bl_hull = hull.boundary.length
            cc = (4*np.pi*am)/bl_hull**2
            
            #lists for this image pair
            cde_area.append(a)
            cde_rr.append(rr)
            cde_fr.append(fr)
            cde_cc.append(cc)
            
            #min and max distance between LKFs
            cde_lkfd_min.append(rs*2/1000)
            cde_lkfd_max.append(rl*2/1000)

        
        #time series lists
        cde_area_m.append(np.mean(cde_area))
        cde_area_s.append(np.std(cde_area))
        
        cde_rr_m.append(np.mean(cde_rr))
        cde_rr_s.append(np.std(cde_rr))
        
        cde_fr_m.append(np.mean(cde_fr))
        cde_fr_s.append(np.std(cde_fr))

        cde_cc_m.append(np.mean(cde_cc))
        cde_cc_s.append(np.std(cde_cc))
        
        cde_lkfd_min_m.append(np.mean(cde_lkfd_min))
        cde_lkfd_min_s.append(np.std(cde_lkfd_min))
        
        cde_lkfd_max_m.append(np.mean(cde_lkfd_max))
        cde_lkfd_max_s.append(np.std(cde_lkfd_max))
         
        cde_date.append(dt)
        cde_num.append(len(floes)/ta*1e6)   #floe density (per km^2)
        cde_lkfa.append(a_lkf)
        cde_lkfa_max.append(a_lkf_max)
        
    else:
        print('Just a single polygon')
    
    ##################################################################################################################################################
    #Plotting
    #Magnification of the sub-region for the concept explanation
    fig3    = plt.figure(figsize=(20,20))
    sx      = fig3.add_subplot(111)
    
    area_def_file = glob(inpath+'area_def_'+date+'*'+file_name_end+'*.pkl')[0]
    with open(area_def_file, "rb") as pkl:
        area_def = pickle.load(pkl)
    m = pr.plot.area_def2basemap(area_def)
    ##limit map extent
    sx.set_xlim(minx-2000,maxx+2000)
    sx.set_ylim(miny-2000,maxy+2000)
    
    ##The main purpose is to map all the intermediate steps in the method
            
    #plot triangle polygons
    if poly_tri.geom_type == 'Polygon':
        xg, yg = poly_tri.exterior.xy 
        sx.fill(xg, yg, alpha=.8, fc='gold', ec='k')
    if poly_tri.geom_type == 'MultiPolygon':
        for geom in poly_tri.geoms:  
            xg, yg = geom.exterior.xy    
            sx.fill(xg, yg, alpha=.8, fc='gold', ec='k')
    
    ##buffer - has exteriors and interiors (holes)
    #if poly_buff.geom_type == 'Polygon':
        #xg, yg = poly_buff.exterior.xy 
        #sx.fill(xg, yg, alpha=0.8, fc='none', ec='purple', lw=3)

        
        #for inner in poly_buff.interiors:
            #xi, yi = zip(*inner.coords[:])
            #sx.fill(xi, yi, alpha=0.8, fc='none', ec='purple', lw=3)
        
        
    #if poly_buff.geom_type == 'MultiPolygon':
        #for geom in poly_buff.geoms:  
            #xg, yg = geom.exterior.xy    
            #sx.fill(xg, yg, alpha=0.8, fc='none', ec='purple', lw=3)
            
            #for inner in geom.interiors:
                #xi, yi = zip(*inner.coords[:])
                #sx.fill(xi, yi, alpha=0.8, fc='none', ec='purple', lw=3)
    #CDCs
    if floes.geom_type == 'Polygon':
        xg, yg = floes.exterior.xy 
        sx.fill(xg, yg, alpha=.5,fc='0.7', ec='k')
    if floes.geom_type == 'MultiPolygon':
        for geom in floes.geoms:  
            xg, yg = geom.exterior.xy    
            sx.fill(xg, yg, alpha=.5, fc='0.7', ec='k')
    
    #LKFs - has exteriors and interiors (holes)
    #interiors need to be without color
    if lkfs.geom_type == 'Polygon':
        patches_p = []
        
        path = Path.make_compound_path(
            Path(np.asarray(lkfs.exterior.coords)[:, :2]),
            *[Path(np.asarray(ring.coords)[:, :2]) for ring in lkfs.interiors])
        
        patch = PathPatch(path)
        
        patches_p.append(patch)
        
        p = PatchCollection(patches_p, ec= 'none', fc='red', alpha=.2)
        sx.add_collection(p)

    if lkfs.geom_type == 'MultiPolygon':
        for geom in lkfs.geoms:  
            patches_p = []
            
            path = Path.make_compound_path(
            Path(np.asarray(geom.exterior.coords)[:, :2]),
            *[Path(np.asarray(ring.coords)[:, :2]) for ring in geom.interiors])
        
            patch = PathPatch(path)
            
            patches_p.append(patch)
            
            p = PatchCollection(patches_p, ec= 'none', fc='red', alpha=.2)
            sx.add_collection(p)
            
    ##region
    #xg, yg = region.exterior.xy 
    #sx.fill(xg, yg, alpha=0.5, fc='none', ec='g', lw=10)
    
    #xg, yg = plot_corners.exterior.xy 
    #sx.fill(xg, yg, alpha=0.2, fc='none', ec='k')
    
    #plot LKF triangles
    patches_p = []
    for k in pindex:
        patch = Polygon(tripts[k])
        patches_p.append(patch)
        
    p = PatchCollection(patches_p, ec= 'k', fc='none', alpha=1)
    sx.add_collection(p)
            
    plotname = plotpath+'CDC_concept'+str(int(bf/1000))+'km'+date.split('T')[0]+'_'+file_name_end+'_4'
    fig3.savefig(plotname,bbox_inches='tight')
    plt.close(fig3)
    print('saved: ',plotname)
    
    ####################################################################
    #Full map with CDCs
    fig6    = plt.figure(figsize=(20,10))
    p1x     = fig6.add_subplot(131)
    
    #annotate the date
    title_date=datetime.strftime(dt, "%Y/%m/%d")
    at = AnchoredText(title_date, prop=dict(size=12), frameon=True, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    p1x.add_artist(at)
    
    #annotate the data
    title_name='CDC size'
    at = AnchoredText(title_name, prop=dict(size=12), frameon=True, loc='upper right')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    p1x.add_artist(at)
    
    #map area definition
    area_def_file = glob(inpath+'area_def_'+date+'*'+file_name_end+'*.pkl')[0]
    with open(area_def_file, "rb") as pkl:
        area_def = pickle.load(pkl)
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,10.),latmax=90.,labels=[0,0,1,0])
    m.drawparallels(np.arange(82.,90.,1),labels=[0,1,0,0])
    
    p3x      = fig6.add_subplot(132)
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,10.),latmax=90.,labels=[0,0,1,0])
    m.drawparallels(np.arange(82.,90.,1),labels=[0,1,0,0])
    
    #annotate the data
    title_name='CDC circularity'
    at = AnchoredText(title_name, prop=dict(size=12), frameon=True, loc='upper right')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    p3x.add_artist(at)
    
    p2x      = fig6.add_subplot(133)
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,10.),latmax=90.,labels=[0,0,1,0])
    m.drawparallels(np.arange(82.,90.,1),labels=[0,1,0,0])
    
    #annotate the data
    title_name='CDC complexity'
    at = AnchoredText(title_name, prop=dict(size=12), frameon=True, loc='upper right')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    p2x.add_artist(at)
    
    ##region should not be multipolygon anymore...
    #xg, yg = region.exterior.xy 
    #p1x.fill(xg, yg, alpha=0.5, fc='none', ec='g')
    #p2x.fill(xg, yg, alpha=0.5, fc='none', ec='g')
    #p3x.fill(xg, yg, alpha=0.5, fc='none', ec='g')
        
    #plot floes
    if floes.geom_type == 'Polygon':
        xg, yg = floes.exterior.xy 
        p1x.fill(xg, yg, alpha=0.5,c='0.7', ec='k')
        p2x.fill(xg, yg, alpha=0.5,c='0.7', ec='k')
        p3x.fill(xg, yg, alpha=0.5,c='0.7', ec='k')
        
    if floes.geom_type == 'MultiPolygon':
        
        print('CDC count: ',len(floes.geoms),' with areas: ',sorted(cde_area))
            
        #use the colors for something meaningful
        interval = [0,50000]
        patches = []
        aaa = []
        for geom in floes.geoms:
            patch = Polygon(np.array(geom.exterior.xy).T)   #convert to matplotlib Polygons
            patches.append(patch)
    
            #the order of polygons is random and it changes every time
            #this is a temporal fix for the plotting
            aaa.append(geom.area/1000000)
        
        pl = PatchCollection(patches, cmap=plt.cm.magma_r, alpha=.5, ec='k')
        pl.set_array(np.array(aaa))
        pl.set_clim(interval)
        mp1 = p1x.add_collection(pl)
        
        cb = plt.colorbar(mp1, ax=p1x, pad=.01, orientation="horizontal")
        
        #print(sorted(cde_fr))
        interval = [2,7]
        patches = []
        aaa = []
        for geom in floes.geoms:
            patch = Polygon(np.array(geom.exterior.xy).T)   #convert to matplotlib Polygons
            patches.append(patch)
            
            centroid = geom.centroid
            boundary = geom.boundary
            rs = centroid.distance(boundary)
            rl = boundary.hausdorff_distance(centroid)
            bl = geom.boundary.length
            diameter = np.mean([rs,rl])*2
            fr = bl/diameter
            aaa.append(fr)
        
        pl = PatchCollection(patches, cmap=plt.cm.magma_r, alpha=.5, ec='k')
        pl.set_array(np.array(aaa))
        pl.set_clim(interval)
        mp2 = p2x.add_collection(pl)
        
        cb = plt.colorbar(mp2, ax=p2x, pad=.01, orientation="horizontal")
        
        interval = [.4,1]
        patches = []
        aaa = []
        for geom in floes.geoms:
            patch = Polygon(np.array(geom.exterior.xy).T)   #convert to matplotlib Polygons
            patches.append(patch)
            
            am = geom.area
            hull = geom.convex_hull
            bl_hull = hull.boundary.length
            cc = (4*np.pi*am)/bl_hull**2
            aaa.append(cc)
        
        pl = PatchCollection(patches, cmap=plt.cm.magma_r, alpha=.5, ec='k')
        pl.set_array(np.array(aaa))
        pl.set_clim(interval)
        mp3 = p3x.add_collection(pl)
        
        cb = plt.colorbar(mp3, ax=p3x, pad=.01, orientation="horizontal")

    #plot LKF triangles over
    patches_p = []
    for k in pindex:
        patch = Polygon(tripts[k])
        patches_p.append(patch)
        
    p = PatchCollection(patches_p, ec= 'k')#, fc=None, alpha=1)
    p1x.add_collection(p)
    
    patches_p = []
    for k in pindex:
        patch = Polygon(tripts[k])
        patches_p.append(patch)
        
    p = PatchCollection(patches_p, ec= 'k')#, fc=None, alpha=1)
    p2x.add_collection(p)
    
    patches_p = []
    for k in pindex:
        patch = Polygon(tripts[k])
        patches_p.append(patch)
        
    p = PatchCollection(patches_p, ec= 'k')#, fc=None, alpha=1)
    p3x.add_collection(p)
    
    #limit map extent
    p1x.set_xlim(minx,maxx)
    p1x.set_ylim(miny,maxy)
    p2x.set_xlim(minx,maxx)
    p2x.set_ylim(miny,maxy)
    p3x.set_xlim(minx,maxx)
    p3x.set_ylim(miny,maxy)
    
    #plot the ship
    mettime = getColumn(shipfile,0)
    dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
    mi = np.argmin(abs(np.asarray(dtb)-dt))
    ship_lon = np.asarray(getColumn(shipfile,1),dtype=float)[mi]
    ship_lat = np.asarray(getColumn(shipfile,2),dtype=float)[mi]
    xl, yl = m(ship_lon, ship_lat)
    print(xl,yl)
    p1x.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    p2x.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    p3x.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

    plotname = plotpath+'CDC_'+str(int(bf/1000))+'km'+date.split('T')[0]+'_'+file_name_end
    fig6.savefig(plotname,bbox_inches='tight')
    plt.close(fig6)
    print('saved: ',plotname)
    #exit()
    #########################################################################################################################################
    
    print('Done with CDC analysis for', date.split('T')[0])
    #exit()
    
    #with open(i, "rb") as poly_file:
        #poly = pickle.load(poly_file)


    #########################################################################################################################################    
    #get radius content
    #radius area (area over which we observe deformaiton)
    #ra = np.arange(5000,radius,2000) #5km increase from 7km to 60km 
    
    fig7    = plt.figure(figsize=(10,10))
    rx      = fig7.add_subplot(111)
    
    #map area definition
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,10.),latmax=90.,labels=[0,0,1,0])
    m.drawparallels(np.arange(82.,90.,1),labels=[0,1,0,0])
    
    #plot floes
    if floes.geom_type == 'Polygon':
        xg, yg = floes.exterior.xy 
        rx.fill(xg, yg, alpha=0.5,c='0.7', ec='k')
        
    if floes.geom_type == 'MultiPolygon':
        for geom in floes.geoms:
            xg, yg = geom.exterior.xy 
            rx.fill(xg, yg, alpha=0.5,c='0.7', ec='k')
    
    n=9 # number of samples
    ra=np.exp(np.linspace(np.log(2000),np.log(radius),n))
    ra = ra.astype(int)
    
    #ship location
    #print(ship_lon,ship_lat)
    #xl, yl = m(ship_lon, ship_lat)
    ship=Point(xl,yl)
    print('Centerpoint: ',ship)
    rx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    
    cl=next(color)
    
    for i in ra:
        print('Radius: ',i)
        
        #make a buffer polygon around ship
        r_region=ship.buffer(i)
        #print(r_region)
        xg, yg = r_region.exterior.xy 
        rx.fill(xg, yg, alpha=0.5, ec='g',fc=None)
        
        #count how many floes we have inside this buffer
        if floes.geom_type == 'MultiPolygon':
        
            floe_counter=0
            #fore every floe
            for geom in floes.geoms:
                #print(geom)
                #it counts if they are partly inside
                if r_region.intersects(geom):
                    floe_counter=floe_counter+1
        else:
            floe_counter=1
            
        print('CDC counter: ',floe_counter)
        
        #plot
        aax.scatter(i*2/1000,floe_counter,alpha=.5,color=cl)
        
        #store to make mean value for each radius
        
    plotname = plotpath+'ra_'+str(int(bf/1000))+'km'+date.split('T')[0]+'_'+file_name_end
    fig7.savefig(plotname,bbox_inches='tight')
    plt.close(fig7)
    print('saved: ',plotname)

#time series of statistics
ax.scatter(cde_date,cde_num)

bx.errorbar(cde_date,cde_area_m,cde_area_s,linestyle='None',marker='o')

cx.errorbar(cde_date,cde_rr_m,cde_rr_s,linestyle='None',marker='o')

dx.errorbar(cde_date,cde_fr_m,cde_fr_s,linestyle='None',marker='o')

ex.scatter(cde_date,cde_lkfa)
ex.scatter(cde_date,cde_lkfa_max)

fx.errorbar(cde_date,cde_lkfd_min_m,cde_lkfd_min_s,linestyle='None',marker='o')
fx.errorbar(cde_date,cde_lkfd_max_m,cde_lkfd_max_s,linestyle='None',marker='o')

#plt.show()

#save cde timeseries
#use the final start/end date in the name
start=datetime.strftime(cde_date[0], "%Y%m%d")
end=datetime.strftime(cde_date[-1], "%Y%m%d")
print(start,end)
plotname=plotpath+'cde_ts_'+start+'_'+end+'_'+str(int(bf/1000))+'km'+file_name_end
print(plotname)
fig1.savefig(plotname,bbox_inches='tight')

#save radius scatter plots
plotname=plotpath+'cde_ra_'+start+'_'+end+'_'+str(int(bf/1000))+'km'+file_name_end
print(plotname)
fig2.savefig(plotname,bbox_inches='tight')

#save output for replotting
tt = [cde_date,cde_num,cde_area_m,cde_area_s,cde_rr_m,cde_rr_s,cde_fr_m,cde_fr_s,cde_cc_m,cde_cc_s,cde_lkfa,cde_lkfa_max,cde_lkfd_min_m,cde_lkfd_min_s,cde_lkfd_max_m,cde_lkfd_max_s]
table = zip(*tt)
table = list(zip(*tt))

output = outpath + outname_cde
print(output)
with open(output, 'ab') as f:
    #header
    f.write(b'date,DE count,DE area_m,DE area_s, DE rr_m,DE rr_s, DE fr_m,DE fr_s,DE cc_m, DE cc_s,LKF area, LKF max area,LKF d_min_m,LKF d_min_s,LKF d_max_m,LKF d_max_s\n')
    np.savetxt(f, table, fmt="%s", delimiter=",")

#convert -delay 100 /scratch/pit000/results/sid/plots_mosaic/CDC_12km202004*__200km_thfilter_lkffilter.png /scratch/pit000/results/sid/plots_mosaic/afs_april2020_anim.gif
