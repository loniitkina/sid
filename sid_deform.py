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

#select lenght scale
lscale = [10,25,50,100,200,500][1]
print(lscale)

#-------------------------------------------------------------------
inpath = '../output/drift_'+str(lscale)+'/'
outpath_def = '../output/def_'+str(lscale)+'/'
outpath = outpath_def+'plots/'
metfile = '../data/10minute_nounits.csv'
reg = 'leg1'

outname_td = 'td_'+reg+'_L'+str(lscale)+'_15km.csv'
td_list=[]
ls_list=[]
ang_list=[]
time_list=[]
date_list=[]

fl = sorted(glob(inpath+'*_upm.npy'))
for i in fl:
    #read in all the data
    print(i)
    u = np.load(i)
    v = np.load(i.split('_upm')[0]+'_vpm.npy')          #in m/s???
    rpm = np.load(i.split('_upm')[0]+'_rpm.npy')       #cross-correlation matrix
    lat = np.load(i.split('_upm')[0]+'_lat1pm.npy')
    lon = np.load(i.split('_upm')[0]+'_lon1pm.npy')
    
    lat2 = np.load(i.split('_upm')[0]+'_lat2pm.npy')
    lon2 = np.load(i.split('_upm')[0]+'_lon2pm.npy')
    
    #get time difference
    date1 = i.split('/')[-1].split('_')[1]
    date2 = i.split('/')[-1].split('_')[2]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    
    if dt1 > datetime(2015,2,12): continue

    diff = (dt2-dt1).seconds + (dt2-dt1).days*24*60*60
    
    #Lance postion (from Lance's met system)
    mettime = getColumn(metfile,0)
    dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
    if dtb[0]>dt1: continue
    if dtb[-1]<dt1: continue
    mi = np.argmin(abs(np.asarray(dtb)-dt1))
    Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
    Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]
    if np.isnan(Lance_lon): continue
    
    #select a region around Lance that corresponds to the buoy array area
    #radius of 100 km
    #otherwise to much sea ice deep in the Arctic is included (reduced scaling slope)
    
    #convert Lance position and sea ice drift array to projected coordinates
    #project the coordinates (units of distance have to be meters)
    #use QGIS for the first guess about the coordinates
    area_def = pr.utils.load_area('area.cfg', reg)
    m = pr.plot.area_def2basemap(area_def)
        
    #reproject vertices
    x, y = m(lon, lat)
    
    #recalculate u,v
    x2, y2 = m(lon2, lat2)
    u = (x2-x)/diff
    v = (y2-y)/diff
    
    #reproject Lance position
    xl, yl = m(Lance_lon, Lance_lat)


    #cut out region
    radius = 15000
    mask = (x<xl-radius) | (x>xl+radius) | (y<yl-radius) | (y>yl+radius)
    
    x = np.ma.array(x,mask=mask)
    x = np.ma.compressed(x)
    y = np.ma.array(y,mask=mask)
    y = np.ma.compressed(y)
    u = np.ma.array(u,mask=mask)
    u = np.ma.compressed(u)
    v = np.ma.array(v,mask=mask)
    v = np.ma.compressed(v)
    rpm = np.ma.array(rpm,mask=mask)
    rpm = np.ma.compressed(rpm)
   
    #get rid of all the nans, the resulting arrays are flattened
    x = np.ma.masked_where(~(np.isfinite(u)),x)
    x = np.ma.compressed(x)
    y = np.ma.masked_where(~(np.isfinite(u)),y)
    y = np.ma.compressed(y)

    
    u = np.ma.masked_invalid(u)
    u = np.ma.compressed(u)
    v = np.ma.masked_invalid(v)
    v = np.ma.compressed(v)
    rpm = np.ma.masked_invalid(rpm)
    rpm = np.ma.compressed(rpm)
    #print lon.shape
    
    #mask out all poor quality data: rpm > 0.4
    mask = rpm > 0.3
    try:
        u = np.ma.array(u,mask=~mask)       #sometimes problems with the mask - as if gridding of lat2pm is not always perfect!
    except:
        continue
    u = np.ma.compressed(u)
    v = np.ma.array(v,mask=~mask)
    v = np.ma.compressed(v)
    x = np.ma.array(x,mask=~mask)
    x = np.ma.compressed(x)
    y = np.ma.array(y,mask=~mask)
    y = np.ma.compressed(y)
    

    #triangulate betwen the points
    pts = np.zeros((len(x),2))
    pts[:,0]=x; pts[:,1]=y
    try:
        tri = Delaunay(pts)
    except:
        continue
    ##simple plot check
    #plt.triplot(pts[:,0], pts[:,1], tri.simplices.copy())
    #plt.show()
    #exit()
    
    tripts = pts[tri.simplices]
    upts = u[tri.simplices]
    vpts = v[tri.simplices]
    
    #calculate deformation - calculating for triangles one by one?
    dux=[];duy=[];dvx=[];dvy=[];minang=[];area=[]
    for t in range(0,len(tripts)):
        vert = np.asarray(tripts[t])
        uvert = upts[t]
        vvert = vpts[t]
    
        try:
            a,b,c,d,e,f=deformation(vert,uvert,vvert)
        except:
            continue
        dux.append(a);duy.append(b);dvx.append(c);dvy.append(d);minang.append(e);area.append(f)
        
    #print area
    dux = np.array(dux)
    duy = np.array(duy)
    dvx = np.array(dvx)
    dvy = np.array(dvy)
    minang = np.array(minang)
    area = np.array(area)
    
    
    div = dux + dvy
    shr = .5*np.sqrt((dux-dvy)**2+(duy+dvx)**2)
    td = np.sqrt(div**2 + shr**2)
    ls = np.sqrt(area)
        
    #storing data for the scatter plots
    td_list.extend(td.tolist())
    ls_list.extend(ls.tolist())
    ang_list.extend(minang.tolist())
    
    #print(td)
    #print(diff_tri)
    #print(dt_tri.astype(datetime))
    #exit()
    dt_tri = np.full_like(dux,np.datetime64(dt1))
    diff_tri = np.ones_like(dux)*diff
    
    date_list.extend(dt_tri.tolist())
    time_list.extend(diff_tri.tolist())
    
    #print(len(diff_tri.tolist()))
    #exit()
    
    ##plotting 
    #deform = div*1e6
    #outname = 'map_div_'+str(lscale)+'_'+date1
    #label = r'Divergence (10$^6$s$^{-1}$)'
    #interval = [-50, 50]
    #cmap=plt.cm.bwr
    ##title = 'c'
    
    ##deform = shr*1e6
    ##outname = 'map_shr_'+str(lscale)+'_'+date1
    ##label = r'Shear (10$^6$s$^{-1}$)'
    ##interval = [0, 50]
    ##cmap=plt.cm.Reds

    #print(outname)
    ##deformation plots
    #fig3    = plt.figure(figsize=(10,10))
    #cx      = fig3.add_subplot(111)
    ##cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
    ##cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

    #area_def = pr.utils.load_area('area.cfg', reg)
    #m = pr.plot.area_def2basemap(area_def)
    
    ##scale
    #m.drawmapscale(15, 82, 25, 82, 100, units='km', barstyle='fancy',fontsize=14)
    #m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    #m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])


    ##nods
    #cx.plot(x,y,'o',linewidth=2,color='purple')
    
    ###speed
    ##speed = np.sqrt(u**2+v**2)
    ##sc = cx.scatter(x,y,s=50,c=speed,cmap=plt.cm.Reds)
    
    ##triangles
    #patches = []
    #for i in range(deform.shape[0]):
        #patch = Polygon(tripts[i,:,:], edgecolor='orchid', alpha=1, fill=False)
        ##plt.gca().add_patch(patch)
        #patches.append(patch)
    
    ##plot filled triangles
    #p = PatchCollection(patches, cmap=cmap, alpha=0.4)
    #p.set_array(np.array(deform))
    #p.set_clim(interval)
    #cx.add_collection(p)
    
    ## create an axes on the right side of ax. The width of cax will be 5%
    ## of ax and the padding between cax and ax will be fixed at 0.05 inch.
    #divider = make_axes_locatable(cx)
    #cax = divider.append_axes("bottom", size="5%", pad=0.1)
    #cbar = plt.colorbar(p, cax=cax, orientation='horizontal')
    #cbar.set_label(label,size=16)

    ##Lance
    ##xl, yl = m(Lance_lon, Lance_lat)
    #cx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    

    #fig3.savefig(outpath+outname,bbox_inches='tight')
    ##exit()
    

#write out lists into csv file
tt = [date_list, ls_list, time_list, td_list, ang_list]
table = zip(*tt)
#adjusted to python3:
table = list(zip(*tt))

output = outpath_def + outname_td

with open(output, 'wb') as f:
  #header
  f.write(b'date, length scale, time difference, total deformation, min angle\n')
  np.savetxt(f, table, fmt="%s", delimiter=",")
