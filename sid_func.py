#! /usr/bin/python

import numpy as np
import csv
import osr, gdal
import matplotlib.pyplot as plt
import pyresample as pr
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.interpolate import interpn
from shapely.ops import cascaded_union, polygonize
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, Polygon
from scipy.spatial import Delaunay
from shapely import geometry
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

def minang_tri(vert):
  #find smallest angle to filter out too acute triangles
  
  #sides of triangle
  ta = np.sqrt((vert[1,0]-vert[0,0])**2+(vert[1,1]-vert[0,1])**2)
  tb = np.sqrt((vert[2,0]-vert[1,0])**2+(vert[2,1]-vert[1,1])**2)
  tc = np.sqrt((vert[2,0]-vert[0,0])**2+(vert[2,1]-vert[0,1])**2)
  sides = [ta,tb,tc]
  srt = min(sides)
  srtid = np.argmin(sides)
  lng = np.delete(sides,srtid)
  
  #law of cosine
  #minang(ii) = acosd((sides(1)^2+sides(2)^2-shortest^2)/2/sides(1)/sides(2))
  minang = np.degrees(np.arccos((lng[0]**2+lng[1]**2-srt**2)/2/lng[0]/lng[1]))

  return(minang)

def deformation(vert,uvert,vvert):

  #smallest angle of triangle
  minang = minang_tri(vert)
  
  ##area of triangle
  #Hutchings etal, 2012 (ISPOL)
  area = .5* (vert[0,0]*vert[1,1] - vert[0,1]*vert[1,0] + vert[1,0]*vert[2,1] - vert[1,1]*vert[2,0] + vert[2,0]*vert[0,1] - vert[2,1]*vert[0,0])
  if area < 0: print('not counter-clockwise oriented triangle!'); exit()
  
  #some filtering to throw out the erroneous GPS positioning. If the accuracy is 2m.
  if area < 1: print(area); return 0,0,0,0,0,0	#Hutchings etal, 2011 (SEDNA)
    
  #calculate line integrals
  dux = (.5 / area) * (
    (uvert[1]+uvert[0])*(vert[1,1]-vert[0,1])+
    (uvert[2]+uvert[1])*(vert[2,1]-vert[1,1])+
    (uvert[0]+uvert[2])*(vert[0,1]-vert[2,1]) )
  
  duy = -(.5 / area) * (
    (uvert[1]+uvert[0])*(vert[1,0]-vert[0,0])+
    (uvert[2]+uvert[1])*(vert[2,0]-vert[1,0])+
    (uvert[0]+uvert[2])*(vert[0,0]-vert[2,0]) )
  
  dvx = (.5 / area) * (
    (vvert[1]+vvert[0])*(vert[1,1]-vert[0,1])+
    (vvert[2]+vvert[1])*(vert[2,1]-vert[1,1])+
    (vvert[0]+vvert[2])*(vert[0,1]-vert[2,1]) )
    
  dvy = -(.5 / area) * (
    (vvert[1]+vvert[0])*(vert[1,0]-vert[0,0])+
    (vvert[2]+vvert[1])*(vert[2,0]-vert[1,0])+
    (vvert[0]+vvert[2])*(vert[0,0]-vert[2,0]) )
  
  return dux,duy,dvx,dvy,minang,area

def getColumn(filename, column, delimiter=',', header=True):
    results = csv.reader(open(filename),delimiter=delimiter)
    if header==True:
        next(results, None)
    return [result[column] for result in results]

def logfit(xdata,ydata):
    logx=np.log10(xdata)
    logy=np.log10(ydata)

    # fit a curve to the data using a least squares 1st order polynomial fit
    z = np.polyfit(logx,logy,1)
    p = np.poly1d(z)
    fit = p(logx)

    # get the coordinates for the fit curve
    c_y = [np.min(fit),np.max(fit)]
    c_x = [np.min(logx),np.max(logx)]

    # predict y values of origional data using the fit
    p_y = z[0] * logx + z[1]

    # calculate the y-error (residuals)
    y_err = logy -p_y

    # create series of new test x-values to predict for
    p_x = np.arange(np.min(logx),np.max(logx)+.01,.01)

    # now calculate confidence intervals for new test x-series
    mean_x = np.mean(logx)         # mean of x
    n = len(logx)              # number of samples in origional fit
    tval = stats.t.ppf(1-0.005, n)		# appropriate t value for 2-tailed distribution
    s_err = np.sum(np.power(y_err,2))   # sum of the squares of the residuals

    confs = tval * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((p_x-mean_x),2)/
                ((np.sum(np.power(logx,2)))-n*(np.power(mean_x,2))))))

    # now predict y based on test x-values
    p_y = z[0]*p_x+z[1]

    # get lower and upper confidence limits based on predicted y and confidence intervals
    lower = p_y - abs(confs)
    upper = p_y + abs(confs)

    #get them back on the exponential scale
    k=z[0]
    loga=z[1]
    a=10.0**loga
    ciy_low = 10.0**lower
    ciy_upp = 10.0**upper
    cix = 10.0**p_x

    return(a,k,cix,ciy_upp,ciy_low)

#consider this option too:
#https://stackoverflow.com/questions/26851533/fit-a-power-law-function-to-the-data-with-both-x-and-y-errors

def density(x,y):
    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    return(x,y,z)


#log-spaced bins
def density_lsb( x,y,n=20):
    """
    Scatter plot colored by 2d histogram
    """
    
    vmin = np.min(x); vmax = np.max(x)
    lsbx = np.logspace(np.log10(vmin),np.log10(vmax), n)
    vmin = np.min(y); vmax = np.max(y)
    lsby = np.logspace(np.log10(vmin),np.log10(vmax), n)
    data , x_e, y_e = np.histogram2d( x, y, bins=[lsbx,lsby])
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    return(x,y,z)


def convertXY(xy_source, inproj, outproj):
    # function to convert coordinates

    shape = xy_source[0,:,:].shape
    size = xy_source[0,:,:].size

    # the ct object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inproj, outproj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))

    xx = xy_target[:,0].reshape(shape)
    yy = xy_target[:,1].reshape(shape)

    return xx, yy


def plot_def(area_def,tripts,deform,outname,label,interval,cmap,Lance_lon,Lance_lat,radius):
    from matplotlib.patches import Polygon
    fig3    = plt.figure(figsize=(20,20))
    cx      = fig3.add_subplot(111)
    #cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
    #cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

    #area_def = pr.utils.load_area('area.cfg', proj)  
    m = pr.plot.area_def2basemap(area_def)
    
    #scale
    #should be at the lower edge of the plot and have relative coordinates (not fixed lat,lon)
    #size depends on radius (in km)
    lonll,latll = area_def.get_lonlat(-1,0)
    m.drawmapscale(lonll, latll, lonll, latll, radius/1000, units='km', barstyle='fancy',fontsize=14)
    m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])


    #nods
    #cx.plot(x,y,'o',linewidth=2,color='purple')
    #mesh
    #cx.triplot(pts[:,0], pts[:,1], tri.simplices.copy(), color=clm, alpha=alpha, label=str(j))
    
    ##speed
    #sc = cx.scatter(x,y,s=50,c=speed,cmap=cmap, vmin=.06, vmax=.085)         #add colorbar and remove extreme values
    #cbar = plt.colorbar(sc)
    #cbar.ax.set_ylabel('Drift speed (m/s)',size=22)

    
    #triangles
    patches = []
    for k in range(deform.shape[0]):
        patch = Polygon(tripts[k,:,:], edgecolor='orchid', alpha=.2, fill=False)
        #plt.gca().add_patch(patch)
        patches.append(patch)
    
    #plot filled triangles
    p = PatchCollection(patches, cmap=cmap, alpha=1)
    p.set_array(np.array(deform))
    p.set_clim(interval)
    cx.add_collection(p)
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(cx)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    cbar = plt.colorbar(p, cax=cax, orientation='horizontal')
    cbar.set_label(label,size=16)

    #Lance
    xl, yl = m(Lance_lon, Lance_lat)
    cx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    
    fig3.savefig(outname,bbox_inches='tight')
    plt.close('all')
    
    return
    
def coarse_grain(tripts,tripts_seed,div,shr,limit,minang):
    #make area-weighted mean of all triangles that are in these triangles
    #weighted-mean deformation - calculating for triangles one by one
    #https://stackoverflow.com/questions/14697442/faster-way-of-polygon-intersection-with-shapely/14804366
    from shapely.geometry import Polygon
    from shapely.strtree import STRtree            
    polys = [Polygon(tt) for tt in tripts[:]]
    #print(len(polys)); print(len(div))
    s = STRtree(polys)
    
    #indexes for each trangle in the tree
    index_by_id = dict((id(pt), i)for i, pt in enumerate(polys))
    
            
    div_seed = []
    shr_seed = []
    area_seed = []
    minang_seed = []
    id_seed = []
    for t in range(0,len(tripts_seed)): 
        qg = Polygon(tripts_seed[t])
        #get area and minangle for the seeded triangle
        ars = qg.area
        mas = minang_tri(tripts_seed[t])
        #get list of interection area
        aa = [o.intersection(qg).area for o in s.query(qg) if o.intersects(qg)]
        #get indexes of idividual polygons
        bb = [(index_by_id[id(o)], o.wkt) for o in s.query(qg) if o.intersects(qg)]
        idxs = [ i[0] for i in bb ]
        
        #check size
        too_big = np.array(aa)>limit
        #check their min angles
        mm = minang[idxs]
        too_sharp = mm<15
        #mask out bad triangle
        too_bad = too_big | too_sharp
        aa = np.ma.array(aa,mask=too_bad)
        aa = np.ma.compressed(aa)
        idxs = np.ma.array(idxs,mask=too_bad)
        idxs = np.ma.compressed(idxs).tolist()

        #get weights for the weighted means
        weights = aa/ars
        
        #extract deformation values for weighted means
        dd = div[idxs]
        ss = shr[idxs]
        
        ##collect all lkf IDs of triangles
        #idd = lkf_id[idxs]
        ##get rid of doubles
        #iddd = list(dict.fromkeys(idd))
        ##print(iddd)
        ##get 
        ##print(len(iddd))
        ##exit()
        
        #check that at least 50% of the seeded triangle is covered by the original small triangles (their intersections)
        #get total area covered by small triangles in this seeded traingle
        #none of these values are masked
        aas = np.sum(aa)
        coverage = aas/ars 
        if coverage > .5:
            
            #store area, minang and number of different LKFs contributing to this seed triangle
            area_seed.append(aas)
            minang_seed.append(mas)
            #id_seed.append(len(iddd))
            
            #get weighted means
            ds = np.sum(weights*dd)
            sss = np.sum(weights*ss)
                        
            div_seed.append(ds)
            shr_seed.append(sss)

    div_seed = np.array(div_seed)
    shr_seed = np.array(shr_seed)
    area_seed = np.array(area_seed)
    minang_seed = np.array(minang_seed)
    #id_seed = np.array(id_seed)

    return(div_seed,shr_seed,area_seed,minang_seed)

#find in triangle centroid
def centroid(vertexes):
    _x_list = [vertex [0] for vertex in vertexes]
    _y_list = [vertex [1] for vertex in vertexes]
    _len = len(vertexes)
    _x = sum(_x_list) / _len
    _y = sum(_y_list) / _len
    return(_x, _y)

def get_lkf_id(tri,threshold,pindex):
    n_list = [] #hold record of all triangles that were used in the region
    lkf_id = np.zeros_like(threshold,dtype=np.int)
    #print(lkf_id.shape)
    
    lkf_counter = 1
    print('counter',lkf_counter)
    
    for p in pindex:
        if p in n_list: #if we used this seed before, just skip to next
            continue
        
        else:
            n_list.append(p)
            #start new LKF group for this seed
            lkf_idx = []  #hold record of all triangles used in this LKF
            lkf_idx.append(p)

            #get neighbors of this seed
            n = tri.neighbors[p]
            
            #iterate until neighbors of neighbors are found
            same_lkf=True
            while same_lkf:
                #just use the ones above threshold and never used before
                used=[]
                for nn in n:
                    used.append(nn in n_list)
                n = np.array(n)
                nmask = (n == -1) | used | threshold[n]
                n = np.ma.array(n,mask=nmask); n = np.ma.compressed(n)
        
                if len(n)>0:
                    ni=[]
                    for i in n:
                        lkf_idx.append(i)
                        n_list.append(i)
                        #accummulate next neighbors of all current ones                       
                        ni.extend(tri.neighbors[i])
                    n = ni
                else:
                    same_lkf=False
                
        #assign same ID to all triangles in this LKF group
        #take only really long ones
        if len(lkf_idx)>10:
            lkf_id[lkf_idx] = lkf_counter
            print(lkf_idx, lkf_counter)
            #increase the counter by one
            lkf_counter = lkf_counter+1
        
    return(lkf_id)

def sort_and_simplify(zigzag,origin=Point(0,0),tolerance=1000):
    #get line coordinates
    xl,yl=zigzag.xy
    
    origin=origin
    #if it has many nods (zig-zag) >> sort in relaton to a point far away (origin)
    if len(xl)> 20:
        dist=[]
        for i in range(0,len(xl)):
            d = origin.distance(Point(xl[i],yl[i]))
            dist.append(d)
    
        #sort nods
        line=zip(xl,yl)
        line = [line for _,line in sorted(zip(dist,line))]
        
        #subsample nods (50% of points is enough) 
        line = line[::2]

        #simplify again
        line = LineString(line)
        simple=line.simplify(tolerance=tolerance, preserve_topology=False)
        simple=simple.simplify(tolerance=tolerance, preserve_topology=False)
        simple=simple.simplify(tolerance=tolerance*2, preserve_topology=False)

    else:
        simple=zigzag.simplify(tolerance=tolerance*2, preserve_topology=False)   #not a zigzag actually, use larger tolerance
        
    return(simple)

def lines_angle(l1, l2):
    import math
    #α0=atan2(cy−ay,cx−ax)
    #α1=atan2(dy−by,dx−bx)
    #angle=α1−α0

    #get line nods
    xl1,yl1=l1.xy
    xl2,yl2=l2.xy
    
    #calculate both angles to the x-asis and then between themselves
    alpha1 = math.atan2(yl1[0]-yl1[-1],xl1[0]-xl1[-1])
    alpha2 = math.atan2(yl2[0]-yl2[-1],xl2[0]-xl2[-1])
    alpha = math.degrees(alpha2-alpha1)
    
    return(alpha)

def get_lkf_angle(tri,tripts,threshold,pindex):
    from shapely.ops import unary_union
    from shapely.ops import split
    from shapely.affinity import scale, rotate
    
    #get a list to store LKF lines
    lkfs=[]
    
    #Get cetroids of all high deformation triangles
    ctrd = []
    for p in pindex:
        ctrd.append(Polygon(tripts[p]).centroid)
    
    #sort them so they are ordered from south-west
    origin=Point(0,0)
    dist=[]
    cx=[]
    cy=[]
    for i in range(0,len(ctrd)):
        d = origin.distance(ctrd[i])
        dist.append(d)
        tmp1,tmp2=ctrd[i].xy
        cx.append(np.array(tmp1)[0]);cy.append(np.array(tmp1)[0])
    
    ctrd = [ctrd for _,ctrd in sorted(zip(dist,ctrd))]
    
    #make buffers around these centroids and unify
    multipoint = MultiPoint(ctrd)
    ctrd_buff = unary_union(multipoint.buffer(400))
    
    #long and complicated lines will be lost, shorter are easier to work with
    #can we split this polygons to several features?
    x_min=np.min(cx);x_max=np.max(cx)
    y_min=np.min(cy);y_max=np.max(cy)
    #construct a split line
    split_line = LineString([Point((x_max-x_min)/4,y_min), Point((x_max-x_min)/4,y_max), 
                             Point((x_max-x_min)/2,y_max), Point((x_max-x_min)/2,y_min),
                             Point(3*(x_max-x_min)/4,y_min), Point(3*(x_max-x_min)/4,y_max), 
                             
         Point(x_max,y_max), Point(x_max,(y_max-y_min)/4),Point(x_min,(y_max-y_min)/4),
                             Point(x_min,(y_max-y_min)/2),Point(x_max,(y_max-y_min)/2),
                             Point(x_max,3*(y_max-y_min)/4),Point(x_min,3*(y_max-y_min)/4) ])
    
    split_line_b = split_line.buffer(50)
    
    #split the polygons
    ctrd_buff = ctrd_buff.difference(split_line_b)
    
    #for each polygon get all centroids inside an connect to a line
    lkfs=[]
    if ctrd_buff.geom_type == 'MultiPolygon':
        for geom in ctrd_buff.geoms:
            lkf_ctrd=[]
            for j in ctrd:
                if geom.contains(Point(j)):
                    lkf_ctrd.append(j)
            
            if len(lkf_ctrd)>10:
                line = LineString(lkf_ctrd)
                
                #simplify
                simple=line.simplify(tolerance=500, preserve_topology=False)
                
                #check if line has a lot of nods
                #rotate the origin and do several iterations
                #incrasing tolerance will reduce most of the lines to two-point features
                simple = sort_and_simplify(simple,origin=Point(0,200000),tolerance=500)
                simple = sort_and_simplify(simple,origin=Point(0,0),tolerance=1000)
                simple = sort_and_simplify(simple,origin=Point(0,200000),tolerance=2000)
                simple = sort_and_simplify(simple,origin=Point(0,0),tolerance=4000)             #finish off with this origin as most features are aligned this way
                
                #store
                lkfs.append(simple)
    
    #prepare LKF buffer to work with (larger than the centroid buffer)
    lkf_buff = unary_union(multipoint.buffer(3000))
    
    #divide the lines at all nods
    lkfs_long=[]
    for i in range(0,len(lkfs)):
        
        #spliting the line on the nods
        xl,yl=lkfs[i].xy
        mp = [ Point(p) for p in zip(xl,yl) ]
        mp = MultiPoint(mp)
        segments = split(lkfs[i], mp)
    
        #throw away all short lines
        i=0
        for segment in segments.geoms:
            #print(segment.length)
            if segment.length>7000:     #5km is a good min lenght for a 60km radius area
                
                #keep only lines that are inside the initial buffers
                if lkf_buff.contains(segment):
                    
                    ##uncomment if no line buffering is used
                    ##extend this segment
                    #segment = scale(segment, xfact=2.5,yfact=2.5)
                    ##store
                    #lkfs_long.append(segment)

                    #uncomment if line buffering is used
                    if i==1:
                        #check if we already have a very similar line
                        if not segment_buffer.contains(segment):     #or touches if just lines
                            #extend this segment
                            segment = scale(segment, xfact=2.5,yfact=2.5)
                            #store
                            lkfs_long.append(segment)
                            
                            #extend the buffer of known lines
                            segment_buffer = unary_union([segment_buffer,segment.buffer(2000)])
                            #segment_buffer = unary_union([segment_buffer,segment])   #for lines
                            
                    else:
                        #extend this segment
                        segment = scale(segment, xfact=3,yfact=3)
                        #store
                        lkfs_long.append(segment)

                        #create the buffer of known lines
                        segment_buffer = segment.buffer(2000)
                        #segment_buffer = segment   #for lines
                        
    ##import ipdb; ipdb.set_trace()
    
    #measure all intersection angles
    angles=[]
    for l1 in lkfs_long:
        for l2 in lkfs_long:
            if l1.intersects(l2) and not l1.__eq__(l2):
                #print(l1)
                #print(l2)
                angle=lines_angle(l1,l2)
                
                #get smallest postive angle
                angle = abs(angle)
                if angle > 180:
                    angle = angle-180
                
                #filter out very small angles (angles that are close to 0 or 180)
                if angle > 10 and angle < 170:
                    angles.append(angle)

    #convert to plottable multistring/multilines
    lkfs=MultiLineString(lkfs_long)
    #lkfs=MultiLineString(lkfs)
    
    return(lkfs, lkf_buff, split_line, angles)

def get_distance_order(tri,tripts,pindex,lines):
    from shapely.ops import unary_union
    
    #Get centroids of all high deformation triangles
    ctrd = []
    for p in pindex:
        ctrd.append(Polygon(tripts[p]).centroid)
    
    #make buffer
    multipoint = MultiPoint(ctrd)
    ctrd_buff = unary_union(multipoint.buffer(200))
    
    #Get centroids of all triangles
    ctrd_all = []
    for i in range(0,len(tripts)):
        ctrd_all.append(Polygon(tripts[i]).centroid)
    
    #collect all centroids inside the buffer
    lkfs=[]
    if ctrd_buff.geom_type == 'MultiPolygon':
        for geom in ctrd_buff.geoms:
            lkf_ctrd=[]
            for j in ctrd_all:
                if geom.contains(Point(j)):
                    lkf_ctrd.append(j)
    
    #get minimal distance from centroid buffer
    dist1=[]
    for j in ctrd_all:
        d = ctrd_buff.distance(Point(j))
        dist1.append(d)
    
    dist2=0
    #uncoment this if you really want dist2
    #currently the angle lines are not good enough for this to be useful!!!
    #print(lines)
    #poly=[]
    #for geom in lines:
        #print(geom)
        #poly.append(geom.buffer(200))
    #lkf_buff = unary_union(poly)
    
    #print(lkf_buff)
    
    
    ##get minimal distance from angle line buffer
    #dist2=[]
    #for j in ctrd_all:
        #d = lkf_buff.distance(Point(j))
        #dist2.append(d)
    
    #print(dist2)
    
    
    return(dist1,dist2)
    



def save_geotiff(raster_array, area_def, export_path):
    # set manually the number of bands to one, because we know you have only one layer
    bands_number = 1

    # Create gdal geotiff driver
    gtiff_driver = gdal.GetDriverByName('GTiff')

    # Pick up the numbers format, use gdal.GDT_Float32 for floats
    gtiff_format = gdal.GDT_Float64
    gtiff_options=["COMPRESS=LZW", "PREDICTOR=2", "TILED=YES"]
    gtiff_options = ["COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6", "INTERLEAVE=BAND"]
    gtiff_options = []
        
    # Create output file (empty at this point) and define projection and extent parameters
    gtiff_dataset = gtiff_driver.Create(export_path,
                                                int(area_def.x_size),
                                                int(area_def.y_size),
                                                bands_number,
                                                gtiff_format,
                                                gtiff_options)
          
    # Define area extent for the Geotiff dataset        
    geometry_list = (area_def.area_extent[0],
                            area_def.pixel_size_x,
                            0,
                            area_def.area_extent[3],
                            0,
                            area_def.pixel_size_y * -1)

    # Set projection parameters
    gtiff_dataset.SetGeoTransform(geometry_list)
    srs = osr.SpatialReference()
    #print(srs)
    #print(area_def.proj4_string.encode('ascii'))
    #srs.ImportFromProj4(area_def.proj4_string.encode('ascii'))
    
    
    #tmp='+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6378137 +b=6356752.3142 +units=m +no_defs +type=crs'
    #srs.ImportFromProj4(tmp)
    srs.ImportFromProj4(area_def.proj4_string)
    
    gtiff_dataset.SetProjection(srs.ExportToWkt())

    # Get the empty band from the dataset, so gdal knows where to write the data	
    gtiff_band = gtiff_dataset.GetRasterBand(1)

    # Write the layer (your raster array) data into the geotiff dataset
    #raster_array = np.flipud(raster_array)
    gtiff_band.WriteArray(raster_array)

    gtiff_dataset = None

#from https://gist.github.com/dwyerk/10561690

def alpha_shape(points, alpha):
    """
    Compute the alpha shape (concave hull) of a set
    of points.
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """
    if len(points) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(points).convex_hull

    #coords = np.array([point.coords[0] for point in points])
    coords = np.zeros((len(points),2))
    xs = [ c[0] for c in points ]
    ys = [ c[1] for c in points ]
    coords[:,0]=xs; coords[:,1]=ys
    #print(coords)
    
    tri = Delaunay(coords)
    
    triangles = coords[tri.vertices]
    a = ((triangles[:,0,0] - triangles[:,1,0]) ** 2 + (triangles[:,0,1] - triangles[:,1,1]) ** 2) ** 0.5
    b = ((triangles[:,1,0] - triangles[:,2,0]) ** 2 + (triangles[:,1,1] - triangles[:,2,1]) ** 2) ** 0.5
    c = ((triangles[:,2,0] - triangles[:,0,0]) ** 2 + (triangles[:,2,1] - triangles[:,0,1]) ** 2) ** 0.5
    s = ( a + b + c ) / 2.0
    areas = (s*(s-a)*(s-b)*(s-c)) ** 0.5
    circums = a * b * c / (4.0 * areas)
    
    #print(circums)
    #print(1.0 / alpha)
    
    filtered = triangles[circums < (1.0 / alpha)]
    edge1 = filtered[:,(0,1)]
    edge2 = filtered[:,(1,2)]
    edge3 = filtered[:,(2,0)]
    edge_points = np.unique(np.concatenate((edge1,edge2,edge3)), axis = 0).tolist()
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return cascaded_union(triangles), edge_points

#create random log normal distributon with prescribed means and standard deviation
def generate_lognormal_samples(mean, stdev, n=1):
    """
    from: https://pythonhealthcare.org/2019/02/07/120-generating-log-normal-samples-from-provided-arithmetic-mean-and-standard-deviation-of-original-population/
    Returns n samples taken from a lognormal distribution, based on mean and
    standard deviation calculated from the original non-logged population.
    
    Converts mean and standard deviation to underlying lognormal distribution
    mu and sigma based on calculations desribed at:
        https://blogs.sas.com/content/iml/2014/06/04/simulate-lognormal-data-
        with-specified-mean-and-variance.html
        
    Returns a numpy array of floats if n > 1, otherwise return a float
    """
    
    # Calculate mu and sigma of underlying lognormal distribution
    phi = (stdev ** 2 + mean ** 2) ** 0.5
    mu = np.log(mean ** 2 / phi)
    sigma = (np.log(phi ** 2 / mean ** 2)) ** 0.5
    
    # Generate lognormal population
    generated_pop = np.random.lognormal(mu, sigma , n)
    
    # Convert single sample (if n=1) to a float, otherwise leave as array
    generated_pop = \
        generated_pop[0] if len(generated_pop) == 1 else generated_pop
        
    return generated_pop


#some notes for code development/debugging
#use 'import ipdb; ipdb.set_trace()' above problematic spot
#use 'dir(object)' to get all options for that object
