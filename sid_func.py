#! /usr/bin/python

import numpy as np
import csv
import osr, gdal
import matplotlib.pyplot as plt
import pyresample as pr
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

def getColumn(filename, column, delimiter=','):
    results = csv.reader(open(filename),delimiter=delimiter)
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
    from scipy import stats
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
    from scipy.stats import gaussian_kde
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
    from scipy.interpolate import interpn
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


def plot_def(area_def,tripts,deform,outname,label,interval,cmap,Lance_lon,Lance_lat):
    from matplotlib.patches import Polygon
    fig3    = plt.figure(figsize=(20,20))
    cx      = fig3.add_subplot(111)
    #cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
    #cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

    #area_def = pr.utils.load_area('area.cfg', proj)  
    m = pr.plot.area_def2basemap(area_def)
    
    #scale
    m.drawmapscale(Lance_lon, Lance_lat-.3, Lance_lon+8, Lance_lat-.2, 50, units='km', barstyle='fancy',fontsize=14)
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
    
    ##virtual buoys
    #xb,yb = m(lon_path[i,:,:],lat_path[i,:,:])
    #cx.plot(xb,yb,'o',linewidth=2,color='purple')
    
    ##mark most frequent values (estimated by sid_pl3d.py)
    #interval = [.1,1]   #rhomboids on 21/1/2015 at L1
    #comm = PatchCollection(patches, cmap=plt.cm.Greens, alpha=0.4)
    #cv = td*1e6
    #mask = (cv<interval[0]) |  (cv > interval[1])
    #cv = np.ma.array(cv,mask=mask)
    #comm.set_array(cv)
    ##comm.set_clim(interval)
    #cx.add_collection(comm)
    #outname = 'map_shr_'+reg+'_L'+str(j)+'_'+date1+'common'

    fig3.savefig(outname,bbox_inches='tight')
    plt.close('all')
    
    return
    
def coarse_grain(tripts,tripts_seed,div,shr):
    #make area-weighted mean of all triangles that are in these triangles
    #weighted-mean deformation - calculating for triangles one by one
    #https://stackoverflow.com/questions/14697442/faster-way-of-polygon-intersection-with-shapely/14804366
    from shapely.geometry import Polygon
    from shapely.strtree import STRtree            
    polys = [Polygon(tt) for tt in tripts[:]]
    print(len(polys)); print(len(div))
    s = STRtree(polys)
    
    #indexes for each trangle in the tree
    index_by_id = dict((id(pt), i)for i, pt in enumerate(polys))
            
    div_seed = []
    shr_seed = []
    area_seed = []
    minang_seed = []
    for t in range(0,len(tripts_seed)): 
        qg = Polygon(tripts_seed[t])
        ars = qg.area
        area_seed.append(ars)
        mas = minang_tri(tripts_seed[t])
        minang_seed.append(mas)
        #get list of interection area
        aa = [o.intersection(qg).area for o in s.query(qg) if o.intersects(qg)] #this is a list
        aa = np.array(aa)
        #get weights for the weighted means
        weights = aa/ars
        #get indexes of idividual polygons
        bb = [(index_by_id[id(o)], o.wkt) for o in s.query(qg) if o.intersects(qg)]
        idxs = [ i[0] for i in bb ]
        #extract deformation values for weighted means
        dd = div[idxs]
        ss = shr[idxs]
        #check we have min 50% coverage
        coverage = np.sum(np.ma.array(weights,mask=dd.mask))
        if coverage < 0.5:
            #print('too small area: '+str(coverage))
            div_seed.append(np.nan); shr_seed.append(np.nan)
            continue
        #get weighted means
        ds = np.sum(weights*dd)
        sss = np.sum(weights*ss)
        
        div_seed.append(ds)
        shr_seed.append(sss)
        
    div_seed = np.array(div_seed)
    shr_seed = np.array(shr_seed)
    area_seed = np.array(area_seed)
    minang_seed = np.array(minang_seed)
    
    div_seed = np.ma.masked_invalid(div_seed)
    shr_seed = np.ma.masked_invalid(shr_seed)

    return(div_seed,shr_seed,area_seed,minang_seed)

#find in triangle centroid
def centroid(vertexes):
    _x_list = [vertex [0] for vertex in vertexes]
    _y_list = [vertex [1] for vertex in vertexes]
    _len = len(vertexes)
    _x = sum(_x_list) / _len
    _y = sum(_y_list) / _len
    return(_x, _y)
