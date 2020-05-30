from glob import glob
from datetime import datetime
import numpy as np
import pyresample as pr
from pyproj import Proj, transform
from pyresample.geometry import AreaDefinition
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from sid_func import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

#How long do you want it to run?
first_week=True
#first_week=False    #This will make it run for all the data

#after_storm=True
after_storm=False

#Do you want all output in figures?
image=True
#image=False

#select lenght scale
radius = 20000
file_name_end = '_20km'

#create log-spaced vector and convert it to integers
n=8 # number of samples
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)

##for large radius (over 50km)
#radius = 120000
#file_name_end = '_120km'
#n=9
#stp=np.exp(np.linspace(np.log(1),np.log(800),n))
#stp = stp.astype(int)

##check just the largest triangles
#stp = stp[-2:]

if first_week==True:
    file_name_end = file_name_end+'FW'
elif after_storm==True:
    file_name_end = file_name_end+'AS'
    
file_name_end = file_name_end+'.csv'

#-------------------------------------------------------------------
#inpath = '../sidrift/data/'
inpath = '../sidrift/data/40m_combo/'
outpath_def = inpath
#outpath = '../sidrift/plots/'
outpath = inpath
metfile = '../sidrift/data/10minute_nounits.csv'
reg = 'leg1'
proj = reg
#reg = 'leg1_FYI'
#reg = 'leg1_SYI'
#reg = 'fixed'

##virtual buoys
#out_file = outpath_def+'VB.npz'
#container = np.load(out_file)
#print(container.files)
#lon_path = container['lon_path'][:6,:,:]   #first step is just vb grid initialization
#lat_path = container['lat_path'][:6,:,:]


#remove all the output text files from previous script runs (or they will be appended)
#outname_td = 'td_'+reg+'_L'+str(j)+file_name_end
#outname_ts = 'ts_'+reg+'_L'+str(j)+file_name_end
import os
rlist = glob(outpath_def+'td_*'+reg+'*'+file_name_end)
for fn in rlist:
    os.remove(fn)
rlist = glob(outpath_def+'ts_*'+reg+'*'+file_name_end)
for fn in rlist:
    os.remove(fn)


#map of frames and valid data
#make a map of the whole area for witch the drift data is processed
inProj = Proj(init='epsg:4326')
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))

#Using a projection dictionary
area_id = 'around Lance'
description = 'North Pole LAEA Europe'
proj_id = 'lance'
proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'a':6378137, 'b':6356752.3142, 'units':'m'}
regn = 84; regs = 81
regw = 10; rege = 30
xc1,yc1 = transform(inProj,outProj,regw, regs)
xc2,yc2 = transform(inProj,outProj,rege, regn)
area_extent = (xc1,yc1,xc2,yc2)
area_def2 = AreaDefinition(area_id, description, proj_id, proj_dict, 300, 1000, area_extent)

fig1    = plt.figure(figsize=(10,10))
ax      = fig1.add_subplot(111)
ma = pr.plot.area_def2basemap(area_def2)
ma.drawmapscale(regw+2, regs+.2, regw+10, regs+.3, 50, units='km', barstyle='fancy',fontsize=14)
ma.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
ma.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])





fl = sorted(glob(inpath+'*.npz'))
#fl = fl[:2]

#colors
color=iter(plt.cm.jet_r(np.linspace(0,1,len(fl)+1)))



for i in range(0,len(fl)):
    #read in all the data
    print(fl[i])
    container = np.load(fl[i])
    #print(container.files)
    u = container['upm']#[::j, ::j]  #this is a regular grid/matrix and such indexing produces a resonably well distributed mesh to construct non-acute triangles
    v = container['vpm']#[::j, ::j] 
    #rpm = container['rpm'][::j, ::j]   #cross-correlation matrix (measure of drift quality worse than hpm) - not used
    hpm = container['hpm']##[::j, ::j]    #hessian
    lat = container['lat1']#[::j, ::j] 
    lon = container['lon1']#[::j, ::j] 
    #lat2 = container['lat2'][::j, ::j] #can be used to recalculate drift from displacements - not used
    #lon2 = container['lon2'][::j, ::j] 
    
    ##at L close to the radius, very few trinagles is selected in this way, here an alternative method with random choice:
    #full = container['upm']
    #print(full.shape)
    ##full = np.ma.masked_invalid(full)
    ##u = np.ma.compressed(u)
    #sample = int(radius*2/j/40/10)
    #print(sample)
    #idx = np.random.choice(np.arange(len(full)), sample, replace=False)     #choose 10 in every row
    #print(idx)
    #u = full[idx]
    #print(u.shape)
    
    #v = container['vpm'][idx] 
    ##rpm = container['rpm'][idx]   #cross-correlation matrix (measure of drift quality worse than hpm) - not used
    #hpm = container['hpm'][idx]    #hessian
    #lat = container['lat1'][idx] 
    #lon = container['lon1'][idx] 

    #from L65 on the area dosnt get fully coverd anymore - corners are missing.
    #simple solution would be to insert those corners...
    #clip option of np.take?
    
    #if j > 0:
        ##insert last non NaN row and column
        #u = get_limits(container['upm'],j)
        #v = get_limits(container['vpm'],j)
        #hpm = get_limits(container['hpm'],j)
        #lat = get_limits(container['lat1'],j)
        #lon = get_limits(container['lon1'],j)
    ##else:   
        ##continue
    
    
    #exit()
    #print(u.shape)
    #print(u.size)
    
    #third option - flatten the array and then index
    #u = container['upm'].flatten()[::j]
    
    print('Size of input matrix:')
    print(u.size)
    #continue
    
    #get time difference
    date1 = fl[i].split('_')[-2]
    date2 = fl[i].split('_')[-1].split('.')[0]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    
    #if we want just the data until the week-long gap
    if (dt1 > datetime(2015,1,28)) & (first_week==True): print('First week only'); continue


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

    if reg == 'fixed':   #keep positions from previous step
        if i==0:
            llato = Lance_lat
            llono = Lance_lon        
            
        else:
            Lance_lat = llato
            Lance_lon = llono

    
    #if we want just data during/after storm
    if (dt1 < datetime(2015,2,6)) & (after_storm==True): continue

    
    #select a region around Lance that corresponds to the buoy array area
    #radius of 100 km
    #otherwise to much sea ice deep in the Arctic is included (reduced scaling slope)
    
    #convert Lance position and sea ice drift array to projected image coordinates
    #project the coordinates (units of distance have to be meters)
    #use QGIS for the first guess about the coordinates
    #area_def = pr.utils.load_area('area.cfg', proj)
    
    #Lance centered projection
    #inProj = Proj(init='epsg:4326')
    #outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))
    xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)
    
    ##Using a projection dictionary (same stuff is same as above)
    #area_id = 'around Lance'
    #description = 'North Pole LAEA Europe'
    #proj_id = 'lance'
    #proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'a':6378137, 'b':6356752.3142, 'units':'m'}
    width = radius*2/40 #x m spacing
    height = radius*2/40 #1x m spacing
    area_extent = (xlp-radius,ylp-radius,xlp+radius,ylp+radius)
    area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
    
    m = pr.plot.area_def2basemap(area_def)
        
    #reproject vertices
    x, y = m(lon, lat)
    
    #x= np.ma.masked_invalid(x)
    #print(np.min(x),np.max(x))
    #print(np.unravel_index(np.argmin(x),x.shape),np.unravel_index(np.argmax(x),x.shape))
    
    #y= np.ma.masked_invalid(y)
    #print(np.min(y),np.max(y))
    #print(np.unravel_index(np.argmin(y),x.shape),np.unravel_index(np.argmax(y),x.shape))
    
    
    #exit()

    
    ##recalculate u,v
    #x2, y2 = m(lon2, lat2)
    #u = (x2-x)/diff
    #v = (y2-y)/diff
            
    ##possible fix for old version velocities
    #u = u/diff
    #v = v/diff
    
    print(u.size)
    
            
    #reproject Lance position
    xl, yl = m(Lance_lon, Lance_lat)
    xla, yla = ma(Lance_lon, Lance_lat)
    
    
    #cut out
    if reg == 'leg1_FYI':   #shift region 25km southwards into the pure FYI zone
        #yl = yl-5000
        xl = xl+20000

    if reg == 'leg1_SYI':   #shift region 25km northwards into the pure SYI zone
        yl = yl+20000
        #xl = xl-20000
        
        
        
    
    #cut out region
    mask = (x<xl-radius) | (x>xl+radius) | (y<yl-radius) | (y>yl+radius)
    
    #find corners
    #add all 4 corner points and centers of the sides from j==1
    
    ##make radius little shorter first
    #rad12=radius-500
    #rad34=radius-500
    
    #idcr1 = np.unravel_index(np.argmin(np.abs(x-(xl-rad12))*np.abs(y-(yl-rad12))),x.shape)
    #idcr2 = np.unravel_index(np.argmin(np.abs(x-(xl-rad12))*np.abs(y-(yl+rad12))),x.shape)
    #idcr3 = np.unravel_index(np.argmin(np.abs(x-(xl+rad34))*np.abs(y-(yl+rad34))),x.shape)
    #idcr4 = np.unravel_index(np.argmin(np.abs(x-(xl+rad34))*np.abs(y-(yl-rad34))),x.shape)
    #print(idcr1,idcr2,idcr3,idcr4)
    ##idcr1 = list(idcr1)
    ##idcr1[0] = idcr1[0]+40
    ##idcr1[1] = idcr1[1]+40
    ##idcr1 = tuple(idcr1)
    #print(idcr1,idcr2,idcr3,idcr4)
    #print(x[idcr1],x[idcr2],x[idcr3],x[idcr4])
    #print(y[idcr1],y[idcr2],y[idcr3],y[idcr4])
    
    #print(u[idcr1],u[idcr2],u[idcr3],u[idcr4])
    #print(v[idcr1],v[idcr2],v[idcr3],v[idcr4])
    
    
    #exit()
    
    #this will make the code much faster
    x = np.ma.array(x,mask=mask,fill_value=np.nan)
    y = np.ma.array(y,mask=mask,fill_value=np.nan)
    u = np.ma.array(u,mask=mask,fill_value=np.nan)
    v = np.ma.array(v,mask=mask,fill_value=np.nan)
    hpm = np.ma.array(hpm,mask=mask,fill_value=np.nan)
    lon = np.ma.array(lon,mask=mask,fill_value=np.nan)
    lat = np.ma.array(lat,mask=mask,fill_value=np.nan)
    
    
    #make radius little shorter first
    rad12=radius-100
    rad34=radius-100
    
    idcr1 = np.unravel_index(np.argmin(np.sqrt((x-(xl-rad12))**2+(y-(yl-rad12))**2)),x.shape)
    idcr2 = np.unravel_index(np.argmin(np.sqrt((x-(xl-rad12))**2+(y-(yl+rad12))**2)),x.shape)
    idcr3 = np.unravel_index(np.argmin(np.sqrt((x-(xl+rad34))**2+(y-(yl+rad34))**2)),x.shape)
    idcr4 = np.unravel_index(np.argmin(np.sqrt((x-(xl+rad34))**2+(y-(yl-rad34))**2)),x.shape)
    print(idcr1,idcr2,idcr3,idcr4)
    print(x[idcr1],x[idcr2],x[idcr3],x[idcr4])
    print(y[idcr1],y[idcr2],y[idcr3],y[idcr4])
    
    print(u[idcr1],u[idcr2],u[idcr3],u[idcr4])
    print(v[idcr1],v[idcr2],v[idcr3],v[idcr4])
    #exit()
 

    
    ##################################################################3
    ##plot on overview map
    ##reproject vertices
    #xa, ya = ma(lon[::50,::50], lat[::50,::50])
    #ax.plot(xa,ya,'o',color='purple')
    
    
    ##################################################################3
    #mesh test plot
    fig4    = plt.figure(figsize=(10,10))
    gx      = fig4.add_subplot(111)

    colormesh=iter(plt.cm.jet_r(np.linspace(0,1,len(stp)+1)))
    ##################################################################3
    
    #now masking begin that will collapse the grid, so spatial array will disappers - do spatial sampling here
    for j in stp:
        print('Step: '+str(j))
        outname_td = 'td_'+reg+'_L'+str(j)+file_name_end
        outname_ts = 'ts_'+reg+'_L'+str(j)+file_name_end
        td_list=[]
        ls_list=[]
        ang_list=[]
        time_list=[]
        date_list=[]
        date_ts=[]
        mdiv=[]
        mpdiv=[]
        mndiv=[]
        mshr=[]
        
        print(u.shape)
        
        ##instead of samplig by indexing, select points by reducing resolution (witsh, height) in projection and select nearest neigbors with valid numbers
        #width = radius*2/40/j
        #height = radius*2/40/j
        #area_def3 = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
        
        #ms = pr.plot.area_def2basemap(area_def3)
            
        ##reproject vertices
        #xs, ys = m(lon, lat)
        
        ##now find the closest neighbors to these points in x and y
        #cn = np.unravel_index(np.argmin(np.abs(x-xs)*np.abs(y-ys)),x.shape)
        
        #print(cs.shape)
        #print(x.shape)
        
        ##instead: take min and max x,y - divide into equal steps
        #print(np.min(x),np.max(x))
        #print(np.argmin(x),np.argmax(x))
        #print(x[0,0])
        #exit()
        

        
        
        #get subsamples for each lenght scale step
        us = u.filled()[::j,::j]
        vs = v.filled()[::j,::j]
        xs = x.filled()[::j,::j]
        ys = y.filled()[::j,::j]
        hpms = hpm.filled()[::j,::j]

        lons = lon.filled()[::j,::j]
        lats = lat.filled()[::j,::j]
        
        ##get subsamples for each lenght scale step
        #us = u[::j,::j]
        #vs = v[::j,::j]
        #xs = x[::j,::j]
        #ys = y[::j,::j]
        #hpms = hpm[::j,::j]

        #lons = lon[::j,::j]
        #lats = lat[::j,::j]
        
        ##use corners to get subsamples
        ##radius is not perfectly aligned with this projection
        #us = u.filled()[icrn1::j,::j]
        #vs = v.filled()[::j,::j]
        #xs = x.filled()[::j,::j]
        #ys = y.filled()[::j,::j]
        #hpms = hpm.filled()[::j,::j]

        #lons = lon.filled()[::j,::j]
        #lats = lat.filled()[::j,::j]
        
        
        
        print(us.shape)
        
        
        print(idcr1,idcr2,idcr3,idcr4)
        print(x[idcr1],x[idcr2],x[idcr3],x[idcr4])
        print(y[idcr1],y[idcr2],y[idcr3],y[idcr4])
    
        print(u[idcr1],u[idcr2],u[idcr3],u[idcr4])
        print(v[idcr1],v[idcr2],v[idcr3],v[idcr4])
    
    
    
        #exit()

        
        clm = next(colormesh)
        if j > 11:
            gx.plot(xs,ys,'o',c=clm)
            gx.plot(x[idcr1],y[idcr1],'*')
            gx.plot(x[idcr2],y[idcr2],'*')
            gx.plot(x[idcr3],y[idcr3],'*')
            gx.plot(x[idcr4],y[idcr4],'*')
            #Lance
            gx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
        
        #if reg == 'leg1_FYI':   #shift region 25km southwards into the pure FYI zone
            ##yl = yl-5000
            #xl = xl+20000

        #if reg == 'leg1_SYI':   #shift region 25km northwards into the pure SYI zone
            #yl = yl+20000
            ##xl = xl-20000
            
            
            
        
        ##cut out region
        #mask = (xs<xl-radius) | (xs>xl+radius) | (ys<yl-radius) | (ys>yl+radius)
        
        ##find corners
        #if j == 1:
            ##add all 4 corner points and centers of the sides from j==1
            #idcr1 = np.unravel_index(np.argmin(np.abs(x-(xl-radius))*np.abs(y-(yl-radius))),x.shape)
            #idcr2 = np.unravel_index(np.argmin(np.abs(x-(xl-radius))*np.abs(y-(yl+radius))),x.shape)
            #idcr3 = np.unravel_index(np.argmin(np.abs(x-(xl+radius))*np.abs(y-(yl+radius))),x.shape)
            #idcr4 = np.unravel_index(np.argmin(np.abs(x-(xl+radius))*np.abs(y-(yl-radius))),x.shape)
            ##print(idcr1,idcr2,idcr3,idcr4)
            ##print(u[idcr1])

        
        #xs = xs[~mask]
        #ys = ys[~mask]
        #us = us[~mask]
        #vs = vs[~mask]
        #hpms = hpms[~mask]
        
        #print(us.size)
        #print(us.shape)
        
        ##try to get the masked area as matrix again
        #uw = np.zeros(u.shape) + np.nan
    
    
        #uw[mask] = us   
        
        #print(uw.shape)
        
        #exit()

        #lons = lons[~mask]
        #lats = lats[~mask]
       
            
        
        #mask out all poor quality data: rpm < 0.4
        gpi = hpms > 9    #this maskes out artefacts like images edges and wierd lines of scnes, but still leaves in in the 'rhomboids'
        us = us[gpi]
        vs = vs[gpi]
        xs = xs[gpi]
        ys = ys[gpi]


        lons = lons[gpi]
        lats = lats[gpi]


        print(us.size)

        #add corners to the flattend data
        if j > 5:
            xlist = [x[idcr1],x[idcr2],x[idcr3],x[idcr4]]
            ylist = [y[idcr1],y[idcr2],y[idcr3],y[idcr4]]
            ulist = [u[idcr1],u[idcr2],u[idcr3],u[idcr4]]
            vlist = [v[idcr1],v[idcr2],v[idcr3],v[idcr4]]
            xs = np.append(xs,xlist)
            ys = np.append(ys,ylist)
            us = np.append(us,ulist)
            vs = np.append(vs,vlist)

        
        
        
        #get rid of all the nans, the resulting arrays are flattened
        xs = np.ma.masked_where(~(np.isfinite(us)),xs)
        xs = np.ma.compressed(xs)
        ys = np.ma.masked_where(~(np.isfinite(us)),ys)
        ys = np.ma.compressed(ys)

        if j==5:
            lons = np.ma.masked_where(~(np.isfinite(us)),lons)
            lons = np.ma.compressed(lons)
            lats = np.ma.masked_where(~(np.isfinite(us)),lats)
            lats = np.ma.compressed(lats)

        
        #rpm = np.ma.masked_invalid(rpm)
        #rpm = np.ma.compressed(rpm)
        #hpm = np.ma.masked_where(~(np.isfinite(u)),hpm)
        #hpm = np.ma.compressed(hpm)
        #print lon.shape
        
        us = np.ma.masked_invalid(us)
        us = np.ma.compressed(us)
        vs = np.ma.masked_invalid(vs)
        vs = np.ma.compressed(vs)

            
            
        #check how many values are left in the region
        print(us.size)
    
        #print('hello1')
        
        #################################################################3
        if j==5:
            #plot on overview map
            #reproject vertices
            xa, ya = ma(lons, lats)
            cl = next(color)
            ax.plot(xa,ya,'.',color=cl, alpha=.3)
            #Lance
            ax.plot(xla,yla,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')
        
        
        #################################################################3


        
        #triangulate betwen the points
        pts = np.zeros((len(xs),2))
        pts[:,0]=xs; pts[:,1]=ys
        try:
            tri = Delaunay(pts)
        except:
            continue
        
        #mesh test plot
        #clm = next(colormesh)
        if j > 11:
            alpha=1
        else:
            alpha=0.5
        gx.triplot(pts[:,0], pts[:,1], tri.simplices.copy(), color=clm, alpha=alpha, label=str(j))


        tripts = pts[tri.simplices]
        upts = us[tri.simplices]
        vpts = vs[tri.simplices]
        
        #define thershold value
        #at factor=1 and step=10, we get 400m distance between grid points
        #distance between nods==400m or triangle area is 400*400/2 
        #there is one nod, where displacement difference is 40m more than at the other two
        #example velocities nod1=1000m/diff, nod2=1000m/diff, nod3=1040m/diff
        #diff is different for every scene pair and needs to be recalculated
        
        print(tripts[0])
        print(upts[0])
        print(vpts[0])
        
        dummy_vert = np.array([[0,0],[400,0],[0,400]])
        dummy_uvert = np.array([1000/diff,1000/diff,1040/diff])
        dummy_vvert = np.array([1000/diff,1000/diff,1000/diff])
        print(dummy_vert)
        print(dummy_uvert)
        dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_f=deformation(dummy_vert,dummy_uvert,dummy_vvert)
        
        print(dummy_a)
        print(dummy_b)
        print(dummy_e)
        print(dummy_f)
        dummy_div = dummy_a+dummy_b
        print(dummy_div)
        
        #exit()
        
        #calculate deformation - calculating for triangles one by one?
        dux=[];duy=[];dvx=[];dvy=[];minang=[];area=[]
        for t in range(0,len(tripts)):
            vert = np.asarray(tripts[t])
            uvert = upts[t]
            vvert = vpts[t]
        
            try:
                a,b,c,d,e,f=deformation(vert,uvert,vvert)
                
                ##try detecting and filtering out rhomboids
                ##at factor=1 and step=10, we get 400m distance between grid points
                ##if distance between nods==400m or triangle area is 400*400/2 & dux==
                #dummy_vert = #exact spacing as in [::j,::j] step
                #dummy_uvert = #three values, one of them has displacement 40m more
                #dummy_vvert =
                #dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,f=deformation(dummy_vert,dummy_uvert,dummy_vvert)
                
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
        
        #print(minang.shape)
        #print(minang)
        ##exit()
        
        
        div = dux + dvy
        shr = .5*np.sqrt((dux-dvy)**2+(duy+dvx)**2)
        
        ##use threshold
        #threshold = np.abs(div)<abs(dummy_div)
        #div = np.ma.array(div,mask=threshold)
        #shr = np.ma.array(shr,mask=threshold)
        
        td = np.sqrt(div**2 + shr**2)
        ls = np.sqrt(area)
        #ls = np.ma.array(np.sqrt(area),mask=threshold)
        #minang = np.ma.array(minang,mask=threshold)
        
            
        ##can we mask out somehow the values cause by the speed stepfunction?
        ##how do we get a relative displacement?
        ##if known we could filter out all that are exactly the smallest step...
        ##lets start by plotting a pdf of divergence...
        #ax = plt.subplot(131)
        #ddd = div*1e6
        #weights = np.ones_like(ddd) / (len(ddd))
        #dbins = np.arange(-3,3,.1)
        #n, bins, patches = ax.hist(ddd, dbins, weights=weights)
        
        #bx = plt.subplot(132)
        #sss = shr*1e6
        #n, bins, patches = bx.hist(sss, dbins, weights=weights)     
        
        #cx = plt.subplot(133)
        #ttt = td*1e6
        #n, bins, patches = cx.hist(ttt, dbins, weights=weights)
        
        #plt.show()
        ##exit()


            
        #storing data for the scatter plots
        td_list.extend(np.ma.compressed(td).tolist())
        ls_list.extend(np.ma.compressed(ls).tolist())
        ang_list.extend(np.ma.compressed(minang).tolist())
        
        #print(td)
        #print(diff_tri)
        #print(dt_tri.astype(datetime))
        #exit()
        dt_tri = np.full_like(dux,np.datetime64(dt1))
        diff_tri = np.ones_like(dux)*diff
        
        #dt_tri = np.ma.array(dt_tri,mask=threshold)
        #diff_tri = np.ma.array(diff_tri,mask=threshold)
        
        date_list.extend(np.ma.compressed(dt_tri).tolist())
        time_list.extend(np.ma.compressed(diff_tri).tolist())
        
        #print(len(diff_tri.tolist()))
        
        #save the data for time series
        date_ts.append(dt1)
        posd = np.ma.array(div,mask=div<0)
        negd = np.ma.array(div,mask=div>0)
        mpdiv.append(np.mean(posd))
        mndiv.append(np.mean(negd))
        mdiv.append(np.mean(div))
        mshr.append(np.mean(shr))
        
        ##save for plotting of Figure 1
        #if date1 == '20150126T070245':        #hack to give me just this one day (for plotting of Figure 1)
            #print(date1)
            #out_file = outpath_def+'Fig1_data_'+date1+'.npz'
            #np.savez(out_file,tripts = tripts,shr = shr, div = div)
            #exit()
        
        #if (i < 6) & (j==1):  #only for the first 6 image pairs and highest resolution
        if image == True:
            #continue
            ##Plotting 
            deform = div*1e6
            #print(deform)
            outname = 'map_div'+reg+'_L'+str(j)+'_'+date1
            label = r'Divergence (10$^6$s$^{-1}$)'
            interval = [-5, 5]
            cmap=plt.cm.bwr
            ##title = 'c'
            
            #deform = shr*1e6
            #outname = 'map_shr_'+reg+'_L'+str(j)+'_'+date1
            #label = r'Shear (10$^6$s$^{-1}$)'
            #interval = [0, 10]
            #cmap=plt.cm.Reds

            #speed = np.sqrt(u**2+v**2)
            #outname = 'map_speed_'+reg+'_L'+str(j)+'_'+date1
            #label = r'Speed (m/s)'
            #interval = [0, 10000]
            #cmap=plt.cm.nipy_spectral

                        
            print(outname)
            #deformation plots
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
            
            ##speed
            #sc = cx.scatter(x,y,s=50,c=speed,cmap=cmap, vmin=.06, vmax=.085)         #add colorbar and remove extreme values
            #cbar = plt.colorbar(sc)
            #cbar.ax.set_ylabel('Drift speed (m/s)',size=22)

            
            #triangles
            patches = []
            for k in range(deform.shape[0]):
                patch = Polygon(tripts[k,:,:], edgecolor='orchid', alpha=1, fill=False)
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
            #xl, yl = m(Lance_lon, Lance_lat)
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


            fig3.savefig(outpath+outname,bbox_inches='tight')
            plt.close('all')
            #exit()
        


        #continue
        #write out lists into csv file
        tt = [date_list, ls_list, time_list, td_list, ang_list]
        table = zip(*tt)
        #adjusted to python3:
        table = list(zip(*tt))

        output = outpath_def + outname_td
        with open(output, 'ab') as f:
            #header
            #f.write(b'date, length scale, time difference, total deformation, min angle\n')
            np.savetxt(f, table, fmt="%s", delimiter=",")

        #write data for time series 
        tt = [date_ts, mpdiv, mndiv, mdiv, mshr]
        table = zip(*tt)
        table = list(zip(*tt))

        output = outpath_def + outname_ts
        with open(output, 'ab') as f:
            #header
            #f.write(b'date, pos. divergence, neg. divergence, mean divergence, mean shear\n')
            np.savetxt(f, table, fmt="%s", delimiter=",")
            
    outname='overview_mesh'+date1
    gx.legend()
    fig4.savefig(outpath+outname)

            
fig1.savefig(outpath+'overview_map_'+reg+'_'+str(int(radius/1000.)),bbox_inches='tight')
    #exit()
