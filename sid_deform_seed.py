from glob import glob
from datetime import datetime
import numpy as np
import pyresample as pr
from pyproj import Proj, transform
from pyresample.geometry import AreaDefinition
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from sid_func import *
import itertools
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

#just save level 1 data and exit
parcel=True
parcel=False

#select lenght scale
radius = 7000
file_name_end = '_7km'

#create log-spaced vector and convert it to integers
n=8 # number of samples
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)
#stp=[1]


if first_week==True:
    file_name_end = file_name_end+'FW'
elif after_storm==True:
    file_name_end = file_name_end+'AS'
    
file_name_end = file_name_end+'.csv'

#-------------------------------------------------------------------
#also set step aand factor for each input data

inpath = '../sidrift/data/40m_combo/'
outpath_def = '../sidrift/data/80m_stp10/'
step = 10
factor = 40    #(factor in sid_drift 0.5)

#inpath = '/media/polona/Polona/s1data/data/stp1_afternoon/'
#outpath_def = '../sidrift/data/80m_stp1/'
#step = 1
#factor = 80     #(factor in sid_drift 0.5)



#outpath = '../sidrift/plots/'
outpath = outpath_def

metfile = '../sidrift/data/10minute_nounits.csv'
reg = 'Lance'
proj = reg
#reg = 'FYI'
#reg = 'SYI'
#reg = 'fixed'

##virtual buoys
#out_file = outpath_def+'VB.npz'
#container = np.load(out_file)
#print(container.files)
#lon_path = container['lon_path'][:6,:,:]   #first step is just vb grid initialization
#lat_path = container['lat_path'][:6,:,:]


#remove all the output text files from previous script runs (or they will be appended)
import os
rlist = glob(outpath_def+'td_*'+reg+'*'+file_name_end)
for fn in rlist:
    os.remove(fn)
rlist = glob(outpath_def+'ts_*'+reg+'*'+file_name_end)
for fn in rlist:
    os.remove(fn)

outname_dummy = 'dummy_'+reg+file_name_end
rlist = glob(outpath_def+outname_dummy)
for fn in rlist:
    os.remove(fn)



#supress warnings
import warnings
warnings.filterwarnings("ignore")

#map of frames and valid data = overview map
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

#get all drift pair files
fl = sorted(glob(inpath+'SeaIceDrift*.npz'))
#colors for overview map
color=iter(plt.cm.jet_r(np.linspace(0,1,len(fl)+1)))

for i in range(0,len(fl)):
    #read in all the data
    print(fl[i])
    container = np.load(fl[i])
    u = container['upm']#[::5,::5]
    v = container['vpm']#[::5,::5]
    hpm = container['hpm']#[::5,::5]    #hessian
    lat = container['lat1']#[::5,::5]
    lon = container['lon1']#[::5,::5]
    
    #print('Size of input matrix:')
    #print(u.size)
        
    #get time difference
    date1 = fl[i].split('_')[-2]
    date2 = fl[i].split('_')[-1].split('.')[0]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    
    #if we want just the data until the week-long gap
    if (dt1 > datetime(2015,1,28)) & (first_week==True): print('First week only'); break

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

    
    #select a region around Lance
    
    #convert Lance position and sea ice drift array to projected image coordinates
    #project the coordinates (units of distance have to be meters)
    #use QGIS for the first guess about the coordinates
    #area_def = pr.utils.load_area('area.cfg', proj)
    
    #Lance centered projection
    #inProj = Proj(init='epsg:4326')
    #outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))
    xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)
    
    if reg == 'FYI':   #shift region into the pure FYI zone
        ylp = ylp-5000
        xlp = xlp+25000

    if reg == 'SYI':   #shift region into the pure SYI zone
        ylp = ylp+5000
        xlp = xlp-25000

    
    ##Using a projection dictionary (same stuff is same as above)
    width = radius*2/40 # m spacing
    height = radius*2/40 # m spacing
    area_extent = (xlp-radius,ylp-radius,xlp+radius,ylp+radius)
    area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
    
    m = pr.plot.area_def2basemap(area_def)
        
    #reproject vertices
    x, y = m(lon, lat)
        
    ##recalculate u,v
    #x2, y2 = m(lon2, lat2)
    #u = (x2-x)/diff
    #v = (y2-y)/diff
            
    ##possible fix for old version velocities
    #u = u/diff
    #v = v/diff
                
    #reproject Lance position
    xl, yl = m(Lance_lon, Lance_lat)
    xla, yla = ma(Lance_lon, Lance_lat)
    
    
    #needs to be cut out from projected space again or coordinates run far beyond the image boundaries and slow down the calculation!
    if reg == 'FYI':   #shift region southwards into the pure FYI zone
        yl = yl-5000
        xl = xl+25000

    if reg == 'SYI':   #shift region northwards into the pure SYI zone
        yl = yl+5000
        xl = xl-25000

    #cut out region
    mask = (x<xl-radius) | (x>xl+radius) | (y<yl-radius) | (y>yl+radius)
    
    #find corners
    idcr1 = np.unravel_index(np.argmin(np.sqrt((x-(xl-radius))**2+(y-(yl-radius))**2)),x.shape)
    idcr2 = np.unravel_index(np.argmin(np.sqrt((x-(xl-radius))**2+(y-(yl+radius))**2)),x.shape)
    idcr3 = np.unravel_index(np.argmin(np.sqrt((x-(xl+radius))**2+(y-(yl+radius))**2)),x.shape)
    idcr4 = np.unravel_index(np.argmin(np.sqrt((x-(xl+radius))**2+(y-(yl-radius))**2)),x.shape)

    #mask the region and fill in nans
    x = np.ma.array(x,mask=mask,fill_value=np.nan)
    y = np.ma.array(y,mask=mask,fill_value=np.nan)
    u = np.ma.array(u,mask=mask,fill_value=np.nan)
    v = np.ma.array(v,mask=mask,fill_value=np.nan)
    hpm = np.ma.array(hpm,mask=mask,fill_value=np.nan)
    lon = np.ma.array(lon,mask=mask,fill_value=np.nan)
    lat = np.ma.array(lat,mask=mask,fill_value=np.nan)
    
    us = u.filled()
    vs = v.filled()
    xs = x.filled()
    ys = y.filled()
    hpms = hpm.filled()
    lons = lon.filled()
    lats = lat.filled()
    
    #mask out all poor quality data: rpm < 0.4
    gpi = hpms > 4    #this maskes out artefacts like images edges and wierd lines on scenes, but still leaves in in the 'rhomboids'
    us = us[gpi]      #it also removes nods at LKFs and makes mesh less dense there >> very few small triangles are left
    vs = vs[gpi]      #9 is very high value, 4 is quite low
    xs = xs[gpi]
    ys = ys[gpi]
    lons = lons[gpi]
    lats = lats[gpi]
    
    #get rid of all the nans, the resulting arrays are flattened
    xs = np.ma.masked_where(~(np.isfinite(us)),xs)
    xs = np.ma.compressed(xs)
    ys = np.ma.masked_where(~(np.isfinite(us)),ys)
    ys = np.ma.compressed(ys)

    
    us = np.ma.masked_invalid(us)
    us = np.ma.compressed(us)
    vs = np.ma.masked_invalid(vs)
    vs = np.ma.compressed(vs)

    #check how many values are left in the region
    #print('Values left for this step:')
    #print(us.size)
    if us.size < 500: continue
    
    #################################################################3
    #plot on overview map
    #reproject vertices
    xa, ya = ma(lons, lats)
    cl = next(color)
    ax.plot(xa,ya,'.',color=cl, alpha=.3)
    #Lance
    ax.plot(xla,yla,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')
    ##################################################################3

    #triangulate betwen the points
    pts = np.zeros((len(xs),2))
    pts[:,0]=xs; pts[:,1]=ys
    tri = Delaunay(pts)

    tripts = pts[tri.simplices]
    upts = us[tri.simplices]
    vpts = vs[tri.simplices]
    
    #################################################################3
    #mesh test plot
    fig4    = plt.figure(figsize=(10,10))
    gx      = fig4.add_subplot(111)
    
    #Lance
    gx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    #full mesh
    gx.triplot(pts[:,0], pts[:,1], tri.simplices.copy(), color='k', alpha=0.5, label='full')
    #colors for the nods
    colormesh=iter(plt.cm.jet_r(np.linspace(0,1,len(stp))))
    ##################################################################3

    
    #define threshold value
    #at factor=1 and step=10, we get 400m distance between grid points
    #distance between nods==400m or triangle area is 400*400/2 
    #there is one nod, where displacement difference is 40m more than at the other two
    #example velocities nod1=1000m/diff, nod2=1000m/diff, nod3=1040m/diff
    #exaggerate the factor as the possible displacements are discrete steps as 40, 80, 120 m...
    #diff is different for every scene pair and threshold needs to be recalculated
    
    #calculate dummy values for all scales, all days and store in a file
    #use on power law plots to show detection limits...
    #should be compared to error estimates by Hutchings et al, 2012: 
    #deformation vealues become reasonably noise free with L >> 340m (A >> 8*3**2*dx**2, dx=40m)
    
    #mean displacement
    dxmean = np.mean(us)*diff
    dymean = np.mean(vs)*diff

    dst = step*factor
    dummy_td_all = []
    dummy_ls_all = []
    for j in stp:
        distance = dst*j
        #print(distance)
        
        dummy_vert = np.array([[0,0],[distance,0],[0,distance]])
        dummy_uvert = np.array([dxmean/diff,dxmean/diff,(dxmean+factor)/diff])      #here use exact grid spacing
        dummy_vvert = np.array([dymean/diff,dymean/diff,dymean/diff])             #this will take out also cases with step with v direction and motion at several nods
        
        dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_f=deformation(dummy_vert,dummy_uvert,dummy_vvert)
        
        dummy_div = dummy_a+dummy_b
        dummy_shr = .5*np.sqrt((dummy_a-dummy_d)**2+(dummy_b+dummy_c)**2)
        dummy_td = np.sqrt(dummy_div**2 + dummy_shr**2)
        dummy_ls = np.sqrt(dummy_f)
    
        dummy_td_all.append(dummy_td)
        dummy_ls_all.append(dummy_ls)
    
    ##store for filtering/masking
    #print(dummy_td_all)
    #print(dummy_td_all[0])
    #exit()
    
    
    #store in a file
    tt = [dummy_ls_all, dummy_td_all]
    table = zip(*tt)
    #adjusted to python3:
    table = list(zip(*tt))

    output = outpath_def + outname_dummy
    with open(output, 'ab') as f:
        np.savetxt(f, table, fmt="%s", delimiter=",")
        
    #now get a slightly higher value for the LKF filtering
    distance = step*factor
    exag_fac = factor+10
        
    dummy_vert = np.array([[0,0],[distance,0],[0,distance]])
    dummy_uvert = np.array([dxmean/diff,dxmean/diff,(dxmean+exag_fac)/diff])    #use average velocities and add a min (exaggerated) step at only one node
    dummy_vvert = np.array([dymean/diff,dymean/diff,dymean/diff])
    
    dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_f=deformation(dummy_vert,dummy_uvert,dummy_vvert)
    dummy_div = dummy_a+dummy_b



    #***************************************************************************************
    #calculate deformation - calculating for triangles one by one
    dux=[];duy=[];dvx=[];dvy=[];minang=[];area=[]
    for t in range(0,len(tripts)):
        vert = np.asarray(tripts[t])
        uvert = upts[t]
        vvert = vpts[t]
    
        #try:
        a,b,c,d,e,f=deformation(vert,uvert,vvert)
        #except:
        #continue
        dux.append(a);duy.append(b);dvx.append(c);dvy.append(d);minang.append(e);area.append(f)
        
    dux = np.array(dux)
    duy = np.array(duy)
    dvx = np.array(dvx)
    dvy = np.array(dvy)
    minang = np.array(minang)
    area = np.array(area)
        
    div = dux + dvy
    shr = .5*np.sqrt((dux-dvy)**2+(duy+dvx)**2)
    
    #use threshold and triangle size criteria to detect noise due to step function in speed
    threshold = (np.abs(div)<abs(dummy_div)) #| (area <= (distance)**2/2)
    
    print(dummy_div)
    
    #threshold = 0
        
    #div = np.ma.array(div,mask=threshold) 
    #shr = np.ma.array(shr,mask=threshold)
    #area = np.ma.array(area,mask=threshold)
    #minang = np.ma.array(minang,mask=threshold)
    
    #store all this for parcel tracking!
    if parcel:
        #find in triangle centroid
        ctrdx = np.zeros_like(div);ctrdy = np.zeros_like(div) 
        for i in range(0,len(tripts)):
            ctrdx[i],ctrdy[i] = centroid(tripts[i])
            
        #convert to lat,lon
        ctrd_lon,ctrd_lat = m(ctrdx,ctrdy,inverse=True)
        
        #get 0,1 array with all damaged triangles as 1
        damage = np.where(~threshold,1,0)
        
        #dump damage data into numpy file
        out_file = outpath_def+'Damage_'+date1+'_'+date2+'.npz'
        np.savez(out_file,lon = ctrd_lon,lat = ctrd_lat, d = damage)

        print('Storing data: ',out_file)
        continue

    ##apply LKF filter
    #prepare place to store filtered data
    
    div_f = np.ma.array(div, mask=threshold)
    shr_f = np.ma.array(shr, mask=threshold)
    
    ##fill with some value or coarse graining will not work
    ##get mean values bellow threshold for mask and use that as fill_value
    #low_div = np.ma.array(div,mask=~threshold)
    #low_mean_div = np.mean(low_div)
    
    #low_shr = np.ma.array(shr,mask=~threshold)
    #low_mean_shr = np.mean(low_shr)
    
    ##or use zeros
    #low_mean_div = 0; low_mean_shr = 0
    
    #print('Fill value bellow threshold:')
    #print(low_mean_div)
    #div_f2 = np.ones_like(div)*low_mean_div
    #shr_f2 = np.ones_like(div)*low_mean_shr
    
    #or just dont touch the low values
    div_f2 = div
    shr_f2 = shr
    
    ##non-masked triangles
    pindex = np.arange(0,len(tri.vertices))
    pindex = np.ma.array(pindex, mask=threshold)
    pindex = np.ma.compressed(pindex)
    
    #plot LKF##################################################################################################3
    fig5    = plt.figure(figsize=(20,10))
    
    #plot original divergence field
    nx      = fig5.add_subplot(121)
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
    
    patches = []
    for k in range(div.shape[0]):
        patch = Polygon(tripts[k,:,:])
        patches.append(patch)

    #plot filled triangles
    p = PatchCollection(patches, cmap=plt.cm.bwr, alpha=1)
    p.set_array(div*1e6)
    interval = [-5, 5]
    p.set_clim(interval)
    nx.add_collection(p)
    
    #set up the plot for the flitered field
    lx      = fig5.add_subplot(122)
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
    #plot basic field
    p = PatchCollection(patches, cmap=plt.cm.bwr, alpha=1)
    p.set_array(div*1e6)
    interval = [-5, 5]
    p.set_clim(interval)
    lx.add_collection(p)
    ###########################################################################################################3
    
    n_list = []
    for p in pindex:
        #get a seed triangle
        #check if this triangle was already part of some other LKF
        if p in n_list: continue
        
        #start new LKF group
        lkf = []
        lkf_div = []
        lkf_shr = []
        lkf_idx = []
                
        lkf.append(tripts[p])
        lkf_div.append(div[p])
        lkf_shr.append(shr[p])  
        lkf_idx.append(p)
        n_list.append(p)                #use pop function of list instead!!!
        
        #get neighbors of this seed, max number is 3 - it should be much more if single vertice neighbours were taken into account
        #use vertex_neighbor_vertices
        #This is not necessary if resolution is higer. Then all triengles in LKF should have a connection along a side.
        n = tri.neighbors[p]
        
        #cycle through the neighbors until no more new are found
        same_lkf = True
        while same_lkf:
            #no neighbor (-1), bellow threshold, already used elements and acute triangles will be masked  
            used=[]
            for nn in n:
                used.append(nn in n_list)
            nmask = (n == -1) | (div_f.mask[n]) | used | (minang[n] < 15)
            n = np.ma.array(n,mask=nmask); n = np.ma.compressed(n)
            if len(n) > 0:
                for i in n:
                    lkf.append(tripts[i])
                    lkf_div.append(div[i])
                    lkf_div.append(shr[i])
                    lkf_idx.append(i)
                    n_list.append(i)
                
                #get next neighbors
                n = tri.neighbors[n]
                #flatten
                n = list(itertools.chain.from_iterable(n))
                #get rid of doubles
                n = list(dict.fromkeys(n))
            else:
               same_lkf = False
                
        #problems in the MIZ: whole belts of LKFs at the image edges...
        #MIZ is area where there is no small triangles... can we filter that out?
                
        #make LKF mean value
        if len(lkf) < 2: 
            #print('Single triangle LKF')
            #these are mainly frame edges (somehow minang does not filter them out)
            #print(minang[lkf_idx])
            #add these triangles to the mask
            threshold[lkf_idx] = True
        else:
            lkf_mdiv = np.mean(lkf_div)
            lkf_mshr = np.mean(lkf_shr)
            #print('Triangle # in LKF: '); print(len(lkf)); print(lkf_mdiv)
            
            #get data to the arrays
            div_f2[lkf_idx] = lkf_mdiv
            shr_f2[lkf_idx] = lkf_mshr
        
            #plot LKF###########################################################################################################3
            patches = []
            for k in range(0,len(lkf)):
                patch = Polygon(lkf[k])
                patches.append(patch)
            
            #plot filled triangles
            p = PatchCollection(patches, cmap=plt.cm.bwr, alpha=1, edgecolor='k')
            lkf_ma = np.ones((len(lkf_div)))*lkf_mdiv*1e6
            p.set_array(lkf_ma)
            interval = [-5, 5]
            p.set_clim(interval)
            lx.add_collection(p)

    fig5.savefig(outpath+'test',bbox_inches='tight')
    #print('Figure saved!')###########################################################################################################3
    #exit()   
    
    #div = div_f2
    #shr = shr_f2
    
    #course-grain these results:
    ddidx = 0
    for j in stp:
        print('Step: '+str(j))
        outname_td = 'td_seed_f_'+reg+'_L'+str(j)+file_name_end
        outname_ts = 'ts_seed_f_'+reg+'_L'+str(j)+file_name_end
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
        
        if j > 1:
            #get seeding points for each lenght scale step
            #flatten the array
            xs = x.filled()[::j,::j]
            ys = y.filled()[::j,::j]
            hpms = hpm.filled()[::j, ::j]
            
            #mask out all poor quality data
            #this maskes out artefacts like images edges and wierd lines on scenes, but still leaves in in the 'rhomboids'
            #it also removes nods at LKFs and makes mesh less dense there
            gpi = hpms > 5    
            xs = xs[gpi]
            ys = ys[gpi]
            
            #add corner points for higher steps
            if j > 5:
                xlist = [x[idcr1],x[idcr2],x[idcr3],x[idcr4]]
                ylist = [y[idcr1],y[idcr2],y[idcr3],y[idcr4]]
                xs = np.append(xs,xlist)
                ys = np.append(ys,ylist)

            #keep only valid data
            xs = np.ma.masked_invalid(xs)
            xs = np.ma.compressed(xs)
            ys = np.ma.masked_invalid(ys)
            ys = np.ma.compressed(ys)

            #check how many nods do we have left
            if len(xs) < 3:
                print('No nods for triangles left!')
                continue

            #triangulate between these seeding points
            pts_seed = np.zeros((len(xs),2))
            pts_seed[:,0]=xs; pts_seed[:,1]=ys
            tri_seed = Delaunay(pts_seed)
            tripts_seed = pts_seed[tri_seed.simplices]
            
            ###############################################################3
            #continue mesh plot
            if j > 11:
                alpha=1
            else:
                alpha=0.5
            clm = next(colormesh)
            gx.triplot(pts_seed[:,0], pts_seed[:,1], tri_seed.simplices.copy(), color=clm, alpha=alpha, label=str(j))
            ###############################################################3
            
            #do area-weighted coarse-graining
            div_seed,shr_seed,area_seed,minang_seed = coarse_grain(tripts,tripts_seed,div_f2,shr_f2)
            
            #print(np.min(np.abs(div_seed)))
            ##print(np.min(area_seed))
            ##print(np.min(minang_seed))
            
            
            ##this will mask out all the long scales...
            ##now mask the triangles that have values exactly like the threshold/all traingles incoming traingles have this threshold value
            #threshold_seed = (np.abs(div_seed)==abs(dummy_div)) #& (area <= (410)**2/2)
            #div_seed = np.ma.array(div_seed,mask=threshold_seed)
            #shr_seed = np.ma.array(shr_seed,mask=threshold_seed)
            
            #print(np.min(np.abs(div_seed)))
            ###exit()
            
            
            #mask with threshold that depends on the lenght scale
            dyn_dummy = dummy_td_all[ddidx]
            threshold_seed = (np.abs(div_seed)<abs(dyn_dummy))
            ddidx = ddidx + 1
            
            div_seed = np.ma.array(div_seed,mask=threshold_seed)
            shr_seed = np.ma.array(shr_seed,mask=threshold_seed)            
            area_seed = np.ma.array(area_seed, mask=threshold_seed)
            minang_seed = np.ma.array(minang_seed, mask=threshold_seed)            
            

        else:
            div_seed=np.ma.array(div_f2, mask=threshold)
            shr_seed=np.ma.array(shr_f2, mask=threshold)
            tripts_seed=tripts
            area_seed = np.ma.array(area, mask=threshold)
            minang_seed = np.ma.array(minang, mask=threshold)
            
            #div_seed=div_f2
            #shr_seed=shr_f2
            #tripts_seed=tripts
            #area_seed = area
            #minang_seed = minang
            
        
        ####################################################################3
        if image == True:
            ##Plotting deformation 
            deform = div_seed*1e6
            outname = outpath+'map_div_seed_f_'+reg+'_L'+str(j)+'_'+date1
            label = r'Divergence (10$^6$s$^{-1}$)'
            interval = [-5, 5]
            cmap=plt.cm.bwr
            
            #deform = shr_seed*1e6
            #outname = outpath+'map_shr_'+reg+'_L'+str(j)+'_'+date1
            #label = r'Shear (10$^6$s$^{-1}$)'
            #interval = [0, 10]
            #cmap=plt.cm.Reds
            
            #if j > 1:
                ##otherwise masked values are plotted
                #deform = deform.filled(fill_value=0)

            
            #print(outname)
            plot_def(area_def,tripts_seed,deform,outname,label,interval,cmap,Lance_lon,Lance_lat)
        #####################################################################3

        #do I still need to mask out everything that is bellow the threshold???
        
        #should threshold depend on the lenght scale???
        #otherwise some deformation values will be masked too 
        #this will not help, they merge anyway as seen on the scatter plots...
        
        #what about if seeding detects how much values bellow threshold were in the triangle
        #for example, if WHOLE triangle is masked, then keep masked
        
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

        #this will cause artificially too steep slope as we gradually incorporate more and more low values into the averages of large triangles...
        
        #this can not be fixed by up-sampling. The values of detectable and un-detectable deformation will still be separate clouds (artificially!)
        #just the threshold will be better and there will be no noise, just where we have really high deformation
        #and no real relativelly low deformtion will be masked out as noise (LKFs will be better connected)
        #this is important for parcel tracking!!!
        
        #so this is perfect for maps or location of LKFs, but not for the values and power law
        #except for the relative comparions (not accros products)
        
        #FYI and SYI could still be compared...
        #the scalling could be different if deformation at highest spatial scale is different or if generally drift speed is different (different threshold value)
        #threshold values should be recored!!!
        
        #what happens if we start at much coarser resolution? E.g. 2km???
        
        #The means of values bellow thershold are too low, because there is still motion/deformation inside those, but is completelly undetected=close to zero
        #the step function detects only resultant cumulative motion/deformation. If averaged with near zero values for the rest of the triangles, that is too low.
        #The only way is to analyze sepatatelly places resolutions where we can separate these two (bellow 2km)
        #and resolutions where we can analyze separatelly completelly merged values (over 10km)
        #but values over 10km are very few as we have to work with image pairs...



        #storing data for the scatter plots
        #mask = div_seed.mask
        
        
        
        
        td = np.sqrt(div_seed**2 + shr_seed**2)
        ls = np.sqrt(area_seed)
        minang = np.array(minang_seed)
        
        
        td_list.extend(np.ma.compressed(td).tolist())
        ls_list.extend(np.ma.compressed(ls).tolist())
        ang_list.extend(np.ma.compressed(minang).tolist())
        
        #time handing
        dt_tri = np.full_like(td,np.datetime64(dt1))
        diff_tri = np.ones_like(td)*diff
        
        
        
        date_list.extend(np.ma.compressed(dt_tri).tolist())
        time_list.extend(np.ma.compressed(diff_tri).tolist())
                
        ##save the data for time series
        #date_ts.append(dt1)
        #posd = np.ma.array(div_seed,mask=div_seed<0)
        #negd = np.ma.array(div_seed,mask=div_seed>0)
        #mpdiv.append(np.mean(posd))
        #mndiv.append(np.mean(negd))
        #mdiv.append(np.mean(div_seed))
        #mshr.append(np.mean(shr_seed))


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

        ##write data for time series 
        #tt = [date_ts, mpdiv, mndiv, mdiv, mshr]
        #table = zip(*tt)
        #table = list(zip(*tt))

        #output = outpath_def + outname_ts
        #with open(output, 'ab') as f:
            ##header
            ##f.write(b'date, pos. divergence, neg. divergence, mean divergence, mean shear\n')
            #np.savetxt(f, table, fmt="%s", delimiter=",")
            
    outname='overview_mesh_seed'+date1
    gx.legend()
    fig4.savefig(outpath+outname)

    
            
fig1.savefig(outpath+'overview_map_'+reg+'_'+str(int(radius/1000.)),bbox_inches='tight')
