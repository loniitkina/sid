from glob import glob
from datetime import datetime
import numpy as np
import pyresample as pr
from pyproj import Proj, transform
from pyresample.geometry import AreaDefinition
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
#from shapely.geometry import Point, MultiPoint
#from shapely.geometry import Polygon as Shapely_Polygon
#from shapely.ops import unary_union
import pickle
from sid_func import * 
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

#WARNING: shapely can be unstable, lots of segfaults
#try: pip uninstall shapely; pip install --no-binary :all: shapely
#works better, but not good

#Conda will have problems with some libraries installed by pip (if shared resources), better if all installed by conda
#try: pip uninstall shapely; conda install shapely
#no, that does not work use this instead: https://anaconda.org/conda-forge/shapely
#does not work :(


#How long do you want it to run?
first_week=True
#first_week=False    #This will make it run for all the data

#after_storm=True
after_storm=False

#Do you want all output in figures?
image=True
#image=False

#active floe size
#afs=True
afs=False

##for parcel tracking we need to have consequtive data: second scene in pair needs to be first scene in the next pair! (combo option is not possible here)
#just save level 1 data and exit
parcel=True
parcel=False

#time series
#time_series=True
time_series=False

#
#no_threshold=True
no_threshold=False

#
LKF_filter=True
#LKF_filter=False

#select lenght scale
radius = 60000
file_name_end = '_60km'

#create log-spaced vector and convert it to integers
n=9 # number of samples
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
#inpath = '../sidrift/data/40m_stp1_afternoon/' #files too big to open???
#canberra
#inpath = 'data/40m_combo/'

outpath_def = '../sidrift/data/80m_stp10_single_filter/'
outpath_def = '../sidrift/data/80m_stp10_nofilter/'

#parcels 
inpath = '../sidrift/data/stp10_parcels_f1/'
#inpath = '../sidrift/data/stp10_parcels_f1_rs2/'
#outpath_def = inpath

#afs 
#outpath_def = '../sidrift/data/stp10_asf/'
#canberra
#outpath_def = 'data/stp10_afs/'

step = 10
factor = 40    #(factor in sid_drift 1)
extra_margin = 20   #20 in test5 gace best results so far (also lets in the last artifacts in LKF filter)

##inpath = '/media/polona/Polona/s1data/data/stp1_afternoon/'
#outpath_def = '../sidrift/data/80m_stp1/'
#inpath = outpath_def
#step = 2
#factor = 40     #(factor in sid_drift 1)
#extra_margin = 10

#inpath = '../sidrift/data/drift_full_time/'
#outpath_def = '../sidrift/data/80m_stp10_time/'
#step = 10
#factor = 80    #(factor in sid_drift 0.5)
#extra_margin = 20                           #this should get the sub-daily data under control (otherwise very noisy)

#outpath = '../sidrift/plots/'
outpath = outpath_def

metfile = '../sidrift/data/10minute_nounits.csv'
#canberra
#metfile = 'data/10minute_nounits.csv'
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
#Use same projection as in pattern matching part of the sea ice drift algorithm
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +ellps=WGS84 +units=m' % (90, 10))

#Using a projection dictionary
area_id = 'around Lance'
description = 'North Pole LAEA Europe'
proj_id = 'lance'
#proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'a':6378137, 'b':6356752.3142, 'units':'m'}
#SRS used in sea ice drift:
#srs = '+proj=laea lat_0=%f lon_0=%f +datum=WGS84 +ellps=WGS84 +no_defs' % (90, 10)
#because displacements and velocities are in this projection we should use it here too!
proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'datum':'WGS84', 'ellps':'WGS84', 'units':'m'}

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
    u = container['upm']#[::10,::10]     
    v = container['vpm']#[::10,::10]
    hpm = container['hpm']#[::10,::10]    #hessian
    lat = container['lat1']#[::10,::10]
    lon = container['lon1']#[::10,::10]
    
    #print('Size of input matrix:')
    #print(u.size)
        
    #get time difference
    date1 = fl[i].split('_')[-2]
    date2 = fl[i].split('_')[-1].split('.')[0]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    
    #if we want just the data until the week-long gap
    if (dt1 > datetime(2015,1,27)) & (first_week==True): print('First week only') ;break

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
        #ylp = ylp-5000
        xlp = xlp+20000

    if reg == 'SYI':   #shift region into the pure SYI zone
        #ylp = ylp+5000
        xlp = xlp-20000

    
    ##Using a projection dictionary (same stuff is same as above)
    width = radius*2/factor # m spacing
    height = radius*2/factor # m spacing
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
                
    #reproject Lance position    ##mask all very small or big triangles
    ##if not masked the range of the ls is big and has several clouds (expected ls, twice the ls and all kinds of smaller ls)
    #center = np.mean(ls)
    ##center = stats.mode(ls)[0][0]                      #this takes too much time
    #print(center)
    #minlen = center-margin[i]; maxlen = center+margin[i]
    ##minlen = center-margin; maxlen = center+margin
    #mask = (ls<minlen) | (ls>maxlen)
    #ls = np.ma.array(ls,mask=mask)
    #td = np.ma.array(td,mask=mask)
    #ls = np.ma.compressed(ls)
    #td = np.ma.compressed(td) #*1e6   

    xl, yl = m(Lance_lon, Lance_lat)
    xla, yla = ma(Lance_lon, Lance_lat)
    
    
    #needs to be cut out from projected space again or coordinates run far beyond the image boundaries and slow down the calculation!
    if reg == 'FYI':   #shift region southwards into the pure FYI zone
        #yl = yl-5000
        xl = xl+20000

    if reg == 'SYI':   #shift region northwards into the pure SYI zone
        #yl = yl+5000
        xl = xl-20000

    #cut out region
    mask = (x<xl-radius) | (x>xl+radius) | (y<yl-radius) | (y>yl+radius)
    
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
    gpi = hpms > 9    #this maskes out artefacts like images edges and wierd lines on scenes, but still leaves in in the 'rhomboids'
    us = us[gpi]      #it also removes nods at LKFs and makes mesh less dense there >> very few small triangles are left
    vs = vs[gpi]      #9 is very high value, 4 is quite low, 4 was used for all initial runs
    xs = xs[gpi]
    ys = ys[gpi]
    lons = lons[gpi]
    lats = lats[gpi]
    
    #find corners
    idcr1 = np.unravel_index(np.argmin(np.sqrt((x-(xl-radius))**2+(y-(yl-radius))**2)),x.shape)
    idcr2 = np.unravel_index(np.argmin(np.sqrt((x-(xl-radius))**2+(y-(yl+radius))**2)),x.shape)
    idcr3 = np.unravel_index(np.argmin(np.sqrt((x-(xl+radius))**2+(y-(yl+radius))**2)),x.shape)
    idcr4 = np.unravel_index(np.argmin(np.sqrt((x-(xl+radius))**2+(y-(yl-radius))**2)),x.shape)
    
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
    if us.size < 100: print(us.size);continue
    
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
        #print(j)
        distance = dst*j
        exag_fac = factor+extra_margin
        
        dummy_vert = np.array([[0,0],[distance,0],[0,distance]])
        dummy_uvert = np.array([dxmean/diff,dxmean/diff,(dxmean+exag_fac)/diff])      #because this is only one nod, there will be divergence (increase in area) and shearing (change in shape, angles)
        dummy_vvert = np.array([dymean/diff,dymean/diff,dymean/diff])             
        
        dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_f=deformation(dummy_vert,dummy_uvert,dummy_vvert)
        
        ##alternative calculation of dummy_b (the rest are zero anyway)
        #duy = .5 * factor * distance / diff / area
        #print(duy)
        #print(dummy_a,dummy_b,dummy_c,dummy_d)
        
        dummy_div = dummy_a+dummy_b
        dummy_shr = .5*np.sqrt((dummy_a-dummy_d)**2+(dummy_b+dummy_c)**2)
        dummy_td = np.sqrt(dummy_div**2 + dummy_shr**2)
        dummy_ls = np.sqrt(dummy_f)
    
        dummy_td_all.append(dummy_td)
        dummy_ls_all.append(dummy_ls)
        
        #print(dummy_div)
        #print(dummy_shr)
        #print(dummy_td) #little higher than div, because shr is low
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
    exag_fac = factor+extra_margin
        
    dummy_vert = np.array([[0,0],[distance,0],[0,distance]])
    dummy_uvert = np.array([dxmean/diff,dxmean/diff,(dxmean+exag_fac)/diff])    #use average velocities and add a min (exaggerated) step at only one node
    dummy_vvert = np.array([dymean/diff,dymean/diff,dymean/diff])
    
    dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_f=deformation(dummy_vert,dummy_uvert,dummy_vvert)
    dummy_div = dummy_a+dummy_b

    print(dummy_div)
    dummy_shr = .5*np.sqrt((dummy_a-dummy_d)**2+(dummy_b+dummy_c)**2)
    dummy_td = np.sqrt(dummy_div**2 + dummy_shr**2)
    print(dummy_td)
    
    if no_threshold==True:
        dummy_td = 0    
    
    ##compare this value to the strain rate error (noise/signal ratio < 0.5)
    #sigma_x = factor
    #mean_area = distance**2/2
    #mean_area2 = (distance*5)**2/2
    #sigma_area = 2*np.sqrt(2)*3*np.sqrt(mean_area2)*sigma_x
    #limit_area = 8*3**2*sigma_x**2
    #mean_displacement = np.sqrt(dxmean**2+dymean**2)
    #print(mean_area)            #this area is way too small to give reliable results (400 m spacing)
    #print(mean_area2)           #but it becomes much better at the next step        (800 m spacing)
    #print(sigma_area)
    #print(limit_area)
    #print(mean_displacement)
        
    ##time error is negligable
    #noise_ratio = 2*np.sqrt((4*(sigma_x**2)/mean_area) + (2*(sigma_x**2)/(mean_displacement**2)) + (sigma_area**2/mean_area**2))
    #print(noise_ratio)
    
    ##only at step=5, noise to signal ration falls bellow 0.5 >> lenght scales over 1 km
    #noise_ratio = 2*np.sqrt((4*sigma_x**2/mean_area2) + (2*sigma_x**2/mean_displacement**2) + (sigma_area**2/mean_area2**2))
    #print(noise_ratio)
    
    #print(4*(sigma_x**2)/mean_area2)
    #print(2*(sigma_x**2)/(mean_displacement**2))
    #print(sigma_area**2/mean_area2**2)                   #this is the problematic term!
    #exit()


    #***************************************************************************************
    #alternative for quad  mesh generation: https://scicomp.stackexchange.com/questions/530/unstructured-quad-mesh-generation
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
    td = np.sqrt(div**2 + shr**2)
    
    #use threshold and triangle size criteria to detect noise due to step function in speed
    #hessian filter got rid of the image artifacts and bad data in the shattered zones (MIZ etc), all large triangle left are of good qualities
    #step function artefacts are all small traingles
    #WARNING: still not all large triangles are part of this!
    threshold = ~((td>dummy_td) | (area > (distance*1.1)**2/2))
    
    #for high resolution data only:
    #if increassing hessian mask to 8, we get larger triangles aroud the LKFs
    #threshold = ~((np.abs(div)>abs(dummy_div)) | (area > (distance*2.5)**2/2))
    
    
    ##threshold average
    #threshold_average_div = 
    #threshold_average_shr =
                
    ##apply LKF filter
    #prepare place to store filtered data
    div_f2 = div.copy()
    shr_f2 = shr.copy()
    
    if LKF_filter==True:
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
        
        #Lance
        xl, yl = m(Lance_lon, Lance_lat)
        nx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

        
        patches_all = []
        for k in range(div.shape[0]):
            patch = Polygon(tripts[k,:,:])
            patches_all.append(patch)

        #plot filled triangles
        p = PatchCollection(patches_all, cmap=plt.cm.bwr, alpha=1)
        p.set_array(div*1e6)
        interval = [-10, 10]
        p.set_clim(interval)
        nx.add_collection(p)
        
        ##check what we are filtering
        #patches_p = []
        #for k in pindex:
            #patch = Polygon(tripts[k])
            #patches_p.append(patch)

        ##plot filled triangles
        #p = PatchCollection(patches_p, ec= 'g', fc=None, alpha=1)
        #nx.add_collection(p)
        
        #set up the plot for the flitered field
        lx      = fig5.add_subplot(122)
        m = pr.plot.area_def2basemap(area_def)
        m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
        m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])        
        ###########################################################################################################3
        
        #n_list = []
        for p in pindex:
            #start new LKF group for each seed triangle
            n_list = []; n_list.append(p)
            lkf = []; lkf.append(tripts[p])
            lkf_div = []; lkf_div.append(div[p])
            lkf_shr = []; lkf_shr.append(shr[p])
            lkf_area = []; lkf_area.append(area[p])
            lkf_idx = []; lkf_idx.append(p)
            
            #get neighbors of this seed
            n = tri.neighbors[p]
            
            #get number of sides of seed triangle occupied by neighbors (max 3) - this will be used to limit the kernel size
            nmask = (n == -1) | threshold[n] | (minang[n] < 5)
            n = np.ma.array(n,mask=nmask); n = np.ma.compressed(n)
            
            #for each side of the seed triangle that has a neighbor
            for nn in n:
                same_side = True
                side = []; side.append(tripts[nn])
                side_div = []; side_div.append(div[nn])
                side_shr = []; side_shr.append(shr[nn])
                side_area = []; side_area.append(area[nn])
                side_idx = []; side_idx.append(nn)
                n_list.append(nn)
                
                while same_side:                    
                    #get next neighbors
                    #nn is growing in every step of the while loop
                    nn = tri.neighbors[nn]
                    
                    if (len(nn)==0):
                        same_side = False
                    else:
                        #flatten if list of lists
                        try: 
                            nn = list(itertools.chain.from_iterable(nn))
                        except:
                            #print('Not list of lists')
                            pass
                        #get rid of doubles
                        nn = list(dict.fromkeys(nn))
                         
                        #mask
                        #no neighbor (-1), bellow threshold, already used elements and acute triangles (near co-planar at the region edge) will be masked  
                        used=[]
                        for nnn in nn:
                            used.append(nnn in n_list)
                        nmask = (nn == -1) | threshold[nn] | used | (minang[nn] < 5)
                        nn = np.ma.array(nn,mask=nmask); nn = np.ma.compressed(nn)
                        
                        #store values
                        for j in nn:
                            #check if we still did not reach the max size of kernel for this side
                            #and if we still have some unused neighbors
                            #3 is recommended max value by Buollion
                            #use < 3 to get to at least 3 (at least one more will be added after this)
                            if (len(side)<3) & (len(nn)>0):
                                side.append(tripts[j])
                                side_div.append(div[j])
                                side_shr.append(shr[j])
                                side_area.append(area[j])
                                side_idx.append(j)
                                n_list.append(j)
                            else:
                                same_side = False
                            
                #collect all values from this side
                lkf.extend(side)
                lkf_div.extend(side_div)
                lkf_shr.extend(side_shr)
                lkf_area.extend(side_area)
                lkf_idx.extend(side_idx)

            ##check what we have:
            #print(div[p])
            #print(lkf_div)
            #print(lkf_idx)
            
                    
            #make LKF mean value
            if len(lkf) < 3: 
                #print('Too short LKF')
                #these are mainly frame edges (somehow minang does not filter them out) and crossings of two step-function artefacts
                #if not filtered out, they will stick out as very high values compared to the averaged ones
                #add these triangles to the mask
                threshold[p] = True
            else:
                #scale all values so that they correspond to the the size of the 'p' triangle
                weights = np.array(lkf_area)/lkf_area[0]    #this is a Q&D way of scaling, ideally we need to use the Oikkonen power law
                lkf_mdiv = np.mean(lkf_div*weights)
                lkf_mshr = np.mean(lkf_shr*weights)
                #print('Triangle # in LKF: '); print(len(lkf)); print(lkf_mdiv)
                #get data to the arrays
                div_f2[p] = lkf_mdiv
                shr_f2[p] = lkf_mshr
                
                #plot LKF###########################################################################################################3   
                
                patches = []
                #for k in lkf_idx:
                patch = Polygon(tripts[p])
                patches.append(patch)
                
                ##plot filled triangles
                pc = PatchCollection(patches, cmap=plt.cm.bwr, alpha=1, edgecolor='k')
                pc.set_array(np.ones((1))*lkf_mdiv*1e6) #has to be an array
                pc.set_clim(interval)
                lx.add_collection(pc)
                
        ##plot all (filtered and bellow threshold) triangles
        #pc = PatchCollection(patches_all, cmap=plt.cm.bwr, alpha=1)
        #pc.set_array(div_f2*1e6)
        #interval = [-5, 5]
        #pc.set_clim(interval)
        #lx.add_collection(pc)

        fig5.savefig(outpath+'test',bbox_inches='tight')
        plt.close(fig5)
        print('Figure saved!')###########################################################################################################3
        #exit()   
    
        #update pindex (triangle index) for afs and lkf_id
        pindex = np.arange(0,len(tri.vertices));pindex = np.ma.array(pindex, mask=threshold); pindex = np.ma.compressed(pindex)
        
        #get region polygon
        xlist = [x[idcr1],x[idcr2],x[idcr3],x[idcr4]]
        ylist = [y[idcr1],y[idcr2],y[idcr3],y[idcr4]]

        if afs==True:
            #dump active floe size input data into numpy file
            out_file = outpath_def+'Afs_'+date1+'_'+date2+file_name_end+'.npz'
            np.savez(out_file,pindex = pindex, tripts = tripts, xlist = xlist, ylist = ylist)

            #for maping
            out_area_def = outpath_def+'area_def_'+date1+'_'+date2+file_name_end+'.pkl'
            f = open(out_area_def,'wb')
            pickle.dump(area_def,f)
            f.close()

            print('Storing data: ',out_file)
            continue
        
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
            lead = np.where((div_f2>0)&(~threshold),1,0)
            ridge = np.where((div_f2<0)&(~threshold),1,0)
            
            #how can i track the area change over this triangle???
            #keep triangulation but apply to lon2,lat2
            #final coordinates can be stored by drift script or calculated from displacements
            
            #dump damage data into numpy file
            out_file = outpath_def+'Damage_'+date1+'_'+date2+'.npz'
            np.savez(out_file,lon = ctrd_lon,lat = ctrd_lat, d = damage, l=lead, r=ridge)

            print('Storing data: ',out_file)
            continue
        
        #get LKF IDs
        #result is an array same shape as div, with ID=0 for all sub-threshold triangles and unique ID for each individual LKF 
        #lkf_id = get_lkf_angle(tri,tripts,threshold,pindex)
        #print(lkf_id)
        
        
        ##plot them all##################################################################################################3
        #fig6    = plt.figure(figsize=(10,10))
        
        ##plot original divergence field
        #xx      = fig6.add_subplot(111)
        #m = pr.plot.area_def2basemap(area_def)
        #m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
        #m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
        
        #patches = []
        #for k in range(threshold.shape[0]):
            #patch = Polygon(tripts[k,:,:])
            #patches.append(patch)

        ##plot filled triangles
        #p = PatchCollection(patches, cmap=plt.cm.bwr, alpha=1)
        #p.set_array(lkf_id)
        #xx.add_collection(p)
        #fig6.savefig(outpath+'LKF_id'+date1,bbox_inches='tight')
        #print('LKF ID figure saved!')##################################################################################3
        ##exit()
        
        #this only works at times when LKF are straight lines
        lkfs, buff, sline = get_lkf_angle(tri,tripts,threshold,pindex)     
        
        #Plotting
        fig6    = plt.figure(figsize=(10,10))
        px      = fig6.add_subplot(111)
        
        #map area definition
        #area_def_file = glob(inpath+'area_def_'+date+'*'+file_name_end+'*.pkl')[0]
        #with open(area_def_file, "rb") as pkl:
            #area_def = pickle.load(pkl)
        m = pr.plot.area_def2basemap(area_def)
        
        m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
        m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
        
        ##Lance
        #xl, yl = m(Lance_lon, Lance_lat)
        #px.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
                    
        #plot the lkfs
        for geom in lkfs:  
            xg, yg = geom.xy  
            px.plot(xg, yg, c='b')
            
        #plot the split line 
        xg, yg = sline.xy  
        px.plot(xg, yg, c='g')
        
        #plot the buffers (typically multipolygons)
        if buff.geom_type == 'Polygon':
            xg, yg = buff.exterior.xy 
            px.fill(xg, yg, alpha=1, fc='none', ec='r')
        if buff.geom_type == 'MultiPolygon':
            for geom in buff.geoms:  
                xg, yg = geom.exterior.xy    
                px.fill(xg, yg, alpha=1, fc='none', ec='r')
                
        ##plot the buffers (typically multipolygons)
        #if sline.geom_type == 'Polygon':
            #xg, yg = sline.exterior.xy 
            #px.fill(xg, yg, alpha=1, fc='none', ec='g')
        #if sline.geom_type == 'MultiPolygon':
            #for geom in sline.geoms:  
                #xg, yg = geom.exterior.xy    
                #px.fill(xg, yg, alpha=1, fc='none', ec='g')
                

        #plot LKF triangles over
        patches_p = []
        for k in pindex:
            patch = Polygon(tripts[k])
            patches_p.append(patch)
            
        p = PatchCollection(patches_p, ec= 'g', fc=None, alpha=1)
        px.add_collection(p)

        fig6.savefig(outpath+'test_lkf_lines_'+date1.split('T')[0],bbox_inches='tight')
        plt.close(fig6)
        
        print('Done with LKF analysis for', date1.split('T')[0])
        #exit()
        continue
        
        
    
    #course-grain these results:
    ddidx = 0
    for j in stp:
        print('Step: '+str(j))
        outname_td = 'td_seed_f_'+reg+'_L'+str(j)+file_name_end
        outname_ts = 'ts_seed_f_'+reg+'_L'+str(j)+file_name_end
        td_list=[]
        ls_list=[]
        ang_list=[]
        id_list=[]
        time_list=[]
        date_list=[]
        date_ts=[]
        mdiv=[]
        mpdiv=[]
        mndiv=[]
        mshr=[]
        mdiv_sd=[]
        mpdiv_sd=[]
        mndiv_sd=[]
        mshr_sd=[]
        
        
        if j > 1:
            #time series dont need corse-graining
            if time_series==True:
                print('time series');continue
            
            
            #get seeding points for each lenght scale step
            #flatten the array
            xs = x.filled()[::j,::j]
            ys = y.filled()[::j,::j]
            hpms = hpm.filled()[::j, ::j]
            
            #mask out all poor quality data
            #5 is sufficient as these are just the seeding polygons
            gpi = hpms > 9    
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
            div_seed,shr_seed,area_seed,minang_seed,id_seed = coarse_grain(tripts,tripts_seed,div_f2,shr_f2,lkf_id)
            
            #mask with threshold that depends on the lenght scale (mask out large triangles)
            #this is done in the plotting, so not necessary here...
            td_seed = np.sqrt(div_seed**2 + shr_seed**2)
            dyn_dummy = dummy_td_all[ddidx]
            ls_dummy = dummy_ls_all[ddidx]
            if no_threshold==True:
                dyn_dummy = 0
            #threshold_seed = (td_seed<dyn_dummy) #| (area_seed > ls_dummy**2)
            threshold_seed = (td_seed<dyn_dummy)
            ddidx = ddidx + 1
            
            #masking
            div_seed = np.ma.array(div_seed,mask=threshold_seed)
            shr_seed = np.ma.array(shr_seed,mask=threshold_seed)            
            area_seed = np.ma.array(area_seed, mask=threshold_seed)
            minang_seed = np.ma.array(minang_seed, mask=threshold_seed)
            id_seed = np.ma.array(id_seed,mask=threshold_seed)
            
            #keep all the data
            
            

        else:
            if LKF_filter==False:
                lkf_id = np.zeros_like(area)
            
            if time_series==False:
                div_seed=np.ma.array(div_f2, mask=threshold)
                shr_seed=np.ma.array(shr_f2, mask=threshold)
                tripts_seed=tripts
                area_seed = np.ma.array(area, mask=threshold)
                minang_seed = np.ma.array(minang, mask=threshold)
                id_seed = np.ma.array(lkf_id,mask=threshold)
            else:
                #keep all the data
                div_seed=div_f2
                shr_seed=shr_f2
                tripts_seed=tripts
                area_seed = area
                minang_seed = minang
                id_seed = lkf_id
            
        
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
            plot_def(area_def,tripts_seed,deform,outname,label,interval,cmap,Lance_lon,Lance_lat,radius)
        #####################################################################3

        ##write out lists into csv file
        #tt = [date_list, ls_list, time_list, td_list, ang_list]
        #table = zip(*tt)
        ##adjusted to python3:
        #table = list(zip(*tt))

        #output = outpath_def + outname_td
        #with open(output, 'ab') as f:
            ##header
            ##f.write(b'date, length scale, time difference, total deformation, min angle\n')
            #np.savetxt(f, table, fmt="%s", delimiter=",")

        #this will cause artificially too steep slope as we gradually incorporate more and more low values into the averages of large triangles...
        
        #at high resolution deformation data we cant avoid of having two separate clouds (detectable and un-detectable deformation) - that happens when our method cant detect deformation anymore
        
        #this can not be fixed by filtering and the line in power law will break even if we use the Sylvain filter
        
        #still the filter removes the artifical opening-clossing mess in shear zones (and removes some high values in that zones), so it makes it useful for the mapping of larger features such as LKFs and leads.
        



        #storing data for the scatter plots
        
        td = np.sqrt(div_seed**2 + shr_seed**2)
        ls = np.sqrt(area_seed)
        minang = np.array(minang_seed)
        id_num = np.array(id_seed)
        
        
        td_list.extend(np.ma.compressed(td).tolist())
        ls_list.extend(np.ma.compressed(ls).tolist())
        ang_list.extend(np.ma.compressed(minang).tolist())
        id_list.extend(np.ma.compressed(id_num).tolist())
        
        #time handing
        dt_tri = np.full_like(td,np.datetime64(dt1))
        diff_tri = np.ones_like(td)*diff
             
        date_list.extend(np.ma.compressed(dt_tri).tolist())
        time_list.extend(np.ma.compressed(diff_tri).tolist())
                
        #save the data for time series
        #these are spatial averages and some value need to be used for all those masked areas - average???
        date_ts.append(dt1)
        posd = np.ma.array(div_seed,mask=div_seed<0)
        negd = np.ma.array(div_seed,mask=div_seed>0)
        mpdiv.append(np.ma.filled(np.mean(posd), fill_value=-999))                     #add fill values (-999)
        mpdiv_sd.append(np.ma.filled(np.std(posd), fill_value=-999))
        mndiv.append(np.ma.filled(np.mean(negd), fill_value=-999))
        mndiv_sd.append(np.ma.filled(np.std(negd), fill_value=-999))
        mdiv.append(np.ma.filled(np.mean(div_seed), fill_value=-999))
        mdiv_sd.append(np.ma.filled(np.std(div_seed), fill_value=-999))
        mshr.append(np.ma.filled(np.mean(shr_seed), fill_value=-999))
        mshr_sd.append(np.ma.filled(np.std(shr_seed), fill_value=-999))

        #continue
        #write out lists into csv file
        tt = [date_list, ls_list, time_list, td_list, ang_list, id_list]
        table = zip(*tt)
        #adjusted to python3:
        table = list(zip(*tt))

        output = outpath_def + outname_td
        with open(output, 'ab') as f:
            #header
            #f.write(b'date, length scale, time difference, total deformation, min angle\n')
            np.savetxt(f, table, fmt="%s", delimiter=",")

        #write data for time series 
        tt = [date_ts, mpdiv, mndiv, mdiv, mshr, mpdiv_sd, mndiv_sd, mdiv_sd, mshr_sd]
        table = zip(*tt)
        table = list(zip(*tt))

        output = outpath_def + outname_ts
        with open(output, 'ab') as f:
            #header
            #f.write(b'date, pos. divergence, neg. divergence, mean divergence, mean shear\n')
            np.savetxt(f, table, fmt="%s", delimiter=",")
            
    outname='overview_mesh_seed'+date1
    gx.legend()
    fig4.savefig(outpath+outname)
    plt.close(fig4)
    
            
fig1.savefig(outpath+'overview_map_'+reg+'_'+str(int(radius/1000.)),bbox_inches='tight')
plt.close(fig1)
