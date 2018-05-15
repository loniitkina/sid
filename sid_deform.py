from glob import glob
from datetime import datetime
import numpy as np
import pyresample as pr
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from sid_func import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

inpath = '../data/Sentinel1_drift/'
outpath_def = '../data/Sentinel1_def/'
outpath = '../plots/'
metfile = '../data/10minute_nounits.csv'
reg = 'leg1'

fl = sorted(glob(inpath+'*_upm.npy'))
for i in fl:
    #read in all the data
    u = np.load(i)
    v = np.load(i.split('_upm')[0]+'_vpm.npy')          #in m/s???
    rpm = np.load(i.split('_upm')[0]+'_rpm.npy')       #cross-correlation matrix
    lat = np.load(i.split('_upm')[0]+'_lat2pm.npy')      #consider using starting coordinates instead
    lon = np.load(i.split('_upm')[0]+'_lon2pm.npy')
    
    #Lance postion (from Lance's met system)
    name2 = i.split('/')[-1]
    date = name2.split('_')[4]
    print date
    dt = datetime.strptime(date, "%Y%m%dT%H%M%S")
    mettime = getColumn(metfile,0)
    dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
    if dtb[0]>dt: continue
    if dtb[-1]<dt: continue
    mi = np.argmin(abs(np.asarray(dtb)-dt))
    Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
    Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]
    if np.isnan(Lance_lon): continue
    
    #check time differences - they need to be similar for the data to be comparable!
    
    #get rid of all the nans, the resulting arrays are flattened
    u = np.ma.masked_invalid(u)
    u = np.ma.compressed(u)
    v = np.ma.masked_invalid(v)
    v = np.ma.compressed(v)
    lat = np.ma.masked_invalid(lat)
    lat = np.ma.compressed(lat)
    lon = np.ma.masked_invalid(lon)
    lon = np.ma.compressed(lon)
    rpm = np.ma.masked_invalid(rpm)
    rpm = np.ma.compressed(rpm)
    #print lon.shape
    
    #mask out all poor quality data: rpm > 0.4
    mask = rpm > 0.3
    u = np.ma.array(u,mask=~mask)
    u = np.ma.compressed(u)
    v = np.ma.array(v,mask=~mask)
    v = np.ma.compressed(v)
    lat = np.ma.array(lat,mask=~mask)
    lat = np.ma.compressed(lat)
    lon = np.ma.array(lon,mask=~mask)
    lon = np.ma.compressed(lon)
    
    #project the coordinates (units of distance have to be meters)
    #use QGIS for the first guess about the coordinates
    area_def = pr.utils.load_area('area.cfg', reg)
    m = pr.plot.area_def2basemap(area_def)
    #reproject vertices
    x, y = m(lon, lat)

    #triangulate betwen the points
    pts = np.zeros((len(x),2))
    pts[:,0]=x; pts[:,1]=y
    tri = Delaunay(pts)
    ##simple plot check
    #plt.triplot(pts[:,0], pts[:,1], tri.simplices.copy())
    #plt.show()
    
    tripts = pts[tri.simplices]
    upts = u[tri.simplices]
    vpts = v[tri.simplices]
    
    #calculate deformation - calculating for triangles one by one?
    dux=[];duy=[];dvx=[];dvy=[]
    for t in range(0,len(tripts)):
        vert = np.asarray(tripts[t])
        uvert = upts[t]
        vvert = vpts[t]
    
        a,b,c,d,minang=deformation(vert,uvert,vvert)
        dux.append(a);duy.append(b);dvx.append(c);dvy.append(d)
        #break
    
    dux = np.array(dux)
    duy = np.array(duy)
    dvx = np.array(dvx)
    dvy = np.array(dvy)
    
    div = dux + dvy
    shr = .5*np.sqrt((dux-dvy)**2+(duy+dvx)**2)
        
    #storing data
    
    #plotting 
    deform = div*1e6
    outname4 = 'map_div_'+date
    label = r'Divergence (10$^6$s$^{-1}$)'
    interval = [-50, 50]
    cmap=plt.cm.bwr
    #title = 'c'
    
    deform = shr*1e6
    outname4 = 'map_shr_'+date
    label = r'Shear (10$^6$s$^{-1}$)'
    interval = [0, 50]
    cmap=plt.cm.Reds
    

    #deformation plots
    fig3    = plt.figure(figsize=(10,10))
    cx      = fig3.add_subplot(111)
    #cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
    #cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

    ##nods
    #cx.plot(x,y,'o',linewidth=2,color='purple')
    
    ##speed
    #speed = np.sqrt(u**2+v**2)
    #sc = cx.scatter(x,y,s=50,c=speed,cmap=plt.cm.Reds)
    
    #triangles
    patches = []
    for i in range(deform.shape[0]):
        patch = Polygon(tripts[i,:,:], edgecolor='orchid', alpha=1, fill=False)
        #plt.gca().add_patch(patch)
        patches.append(patch)
    
    #plot filled triangles
    p = PatchCollection(patches, cmap=cmap, alpha=0.4)
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
    xa, ya = m(Lance_lon, Lance_lat)
    cx.plot(xa,ya,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

    fig3.savefig(outpath+outname4,bbox_inches='tight')
    #exit()
    
    


