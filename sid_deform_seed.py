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
first_week=False    #This will make it run for all the data

after_storm=True
after_storm=False

#Do you want all output in figures?
image=True
#image=False

#active floe size
afs=True
#afs=False

##for parcel tracking we need to have consequtive data: second scene in pair needs to be first scene in the next pair! (combo option is not possible here)
#just save level 1 data and exit
parcel=True
parcel=False

#time series
time_series=True        #this will not scale the data - useful also for the angle PDF!
time_series=False

#
no_threshold=True
no_threshold=False

#
kernel=3    #3 is recommend by Sylvain (but for coarser resolution)
LKF_filter=True
#LKF_filter=False

#low mean adjustment
low_mean=True
low_mean=False

#min angle limit, best set very low or zero when wanting to detect all features (will include MIZ!)
#can be set to 15, when working with scalling
minang_limit=5

#select area size
radius = 120000
file_name_end = '_120km'

#expected divergence values (div in s-1 * 10e6)
interval = [-1, 1]

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
#also set step and factor for each input data

inpath = '../../results/sid/drift/'
outpath_def = '../../results/sid/afs/'



#parcels 
#inpath = '../sidrift/data/stp10_parcels_f1/'
#inpath = '../sidrift/data/stp10_parcels_f1_rs2/'
#outpath_def = inpath

#afs 
#outpath_def = '../sidrift/data/stp10_asf/'
#canberra
#outpath_def = 'data/stp10_afs/'

step = 10
factor = 80    #(factor in sid_drift, default factor in FT is 0.5, which translates to 80 m here)
extra_margin = 20   #20 in test5 gave best results so far (but also lets in the last artifacts in LKF filter)

#20 was a good margin for factor=80, seems like 25% was a good estimate
#no margin necessary for scaling, but use 10 to get clean LKFs for other analysis where 'less is more' (e.g. parcels)
#extra_margin = 0


##CHECK IF WE GET OLD DATA NOW
#step = 10
#factor = 80    #(factor in sid_drift 1)
#extra_margin = 20   #20 in test5 gave best results so far (but also lets in the last artifacts in LKF filter)


##inpath = '/media/polona/Polona/s1data/data/stp1_afternoon/'
#outpath_def = '../sidrift/data/80m_stp1/'
#inpath = outpath_def
#step = 2
#factor = 40     #(factor in sid_drift 1)
#extra_margin = 10

#inpath = '../sidrift/data/drift_full_time/'
#outpath_def = '../sidrift/data/80m_stp10_time/'
##canberra
#inpath = '../../s1data/data/drift_full_factor1/'
#outpath_def = 'data/40m_stp10_time_margin_m10_k4/'
#step = 10
#factor = 40    #(factor in sid_drift 1)
#extra_margin = -10                           #20m margin cuts out too much, especially at shorter than 12h

outpath = '../../results/sid/plots/'

shipfile = '../../downloads/position_leg3_nh-track.csv'
reg = 'ship'
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
area_id = 'around ship'
description = 'North Pole LAEA Europe'
proj_id = 'ship'
#proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'a':6378137, 'b':6356752.3142, 'units':'m'}
#SRS used in sea ice drift:
#srs = '+proj=laea lat_0=%f lon_0=%f +datum=WGS84 +ellps=WGS84 +no_defs' % (90, 10)
#because displacements and velocities are in this projection we should use it here too!
proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'datum':'WGS84', 'ellps':'WGS84', 'units':'m'}

#regn = 84; regs = 81
#regw = 10; rege = 30
#specify region
lon_diff = 15
ship_lon=17.147909; ship_lat=87.132429      #March/April event start
regn = ship_lat+.1; regs = ship_lat-4
regw = ship_lon-lon_diff; rege = ship_lon+lon_diff

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
    
    #print(np.ma.masked_invalid(u).compressed())
    
    print('Size of input matrix:')
    print(u.size)
    
        
    #get time difference
    date1 = fl[i].split('_')[-2]
    date2 = fl[i].split('_')[-1].split('.')[0]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    
    ##if we want just the data until the week-long gap
    #if (dt1 > datetime(2015,1,27)) & (first_week==True): print('First week only') ;break

    diff = (dt2-dt1).seconds + (dt2-dt1).days*24*60*60
    
    #ship postion (from ship's met system)
    mettime = getColumn(shipfile,0)
    dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
    if dtb[0]>dt1: continue
    if dtb[-1]<dt1: continue
    mi = np.argmin(abs(np.asarray(dtb)-dt1))
    ship_lon = np.asarray(getColumn(shipfile,1),dtype=float)[mi]
    ship_lat = np.asarray(getColumn(shipfile,2),dtype=float)[mi]
    if np.isnan(ship_lon): continue

    print('ship at: ',ship_lon,ship_lat,dtb[mi])

    #if reg == 'fixed':   #keep positions from previous step
        #if i==0:
            #llato = ship_lat
            #llono = ship_lon        
            
        #else:
            #ship_lat = llato
            #ship_lon = llono

    
    ##if we want just data during/after storm
    #if (dt1 < datetime(2015,2,3)) & (after_storm==True): continue

    #select a region around ship
    
    #convert ship position and sea ice drift array to projected image coordinates
    #project the coordinates (units of distance have to be meters)
    #use QGIS for the first guess about the coordinates
    #area_def = pr.utils.load_area('area.cfg', proj)
    
    #ship centered projection
    #inProj = Proj(init='epsg:4326')
    #outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))
    xlp,ylp = transform(inProj,outProj,ship_lon,ship_lat)
    
    #if reg == 'FYI':   #shift region into the pure FYI zone
        ##ylp = ylp-5000
        #xlp = xlp+20000

    #if reg == 'SYI':   #shift region into the pure SYI zone
        ##ylp = ylp+5000
        #xlp = xlp-20000

    ##Using a projection dictionary (same stuff as above)
    width = radius*2/factor # m spacing
    height = radius*2/factor # m spacing
    area_extent = (xlp-radius,ylp-radius,xlp+radius,ylp+radius)
    area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
    
    m = pr.plot.area_def2basemap(area_def)
        
    #reproject vertices
    x, y = m(lon, lat)
        
    ###recalculate u,v
    ##x2, y2 = m(lon2, lat2)
    ##u = (x2-x)/diff
    ##v = (y2-y)/diff
            
    ##possible fix for old version velocities - check if this is not done already in sid_drift.py
    #u = u/diff
    #v = v/diff
                
    #reproject ship position    ##mask all very small or big triangles
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

    xl, yl = m(ship_lon, ship_lat)
    print(xl,yl)
    
    
    
    
    ##needs to be cut out from projected space again or coordinates run far beyond the image boundaries and slow down the calculation!
    #if reg == 'FYI':   #shift region southwards into the pure FYI zone
        ##yl = yl-5000
        #xl = xl+20000

    #if reg == 'SYI':   #shift region northwards into the pure SYI zone
        ##yl = yl+5000
        #xl = xl-20000

    #cut out region
    mask = ~((x<xl-radius) | (x>xl+radius) | (y<yl-radius) | (y>yl+radius))
    
    ##mask the region and fill in nans
    #x = np.ma.array(x,mask=mask,fill_value=np.nan)
    #y = np.ma.array(y,mask=mask,fill_value=np.nan)
    #u = np.ma.array(u,mask=mask,fill_value=np.nan)
    #v = np.ma.array(v,mask=mask,fill_value=np.nan)
    #hpm = np.ma.array(hpm,mask=mask,fill_value=np.nan)
    #lon = np.ma.array(lon,mask=mask,fill_value=np.nan)
    #lat = np.ma.array(lat,mask=mask,fill_value=np.nan)
    
    #us = u.filled()
    #vs = v.filled()
    #xs = x.filled()
    #ys = y.filled()
    #hpms = hpm.filled()
    #lons = lon.filled()
    #lats = lat.filled()
    
    ##dummies
    #us = u
    #vs = v
    #xs = x
    #ys = y
    #hpms = hpm
    #lons = lon
    #lats = lat
    
    us = u[mask]
    vs = v[mask]
    xs = x[mask]
    ys = y[mask]
    hpms = hpm[mask]
    lons = lon[mask]
    lats = lat[mask]
    
    
    
    
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


    print(xs,ys)
    #exit()

    lons = (np.ma.masked_where(~(np.isfinite(us)),lons))
    lats = (np.ma.masked_where(~(np.isfinite(us)),lats))
    
    us = np.ma.masked_invalid(us)
    us = np.ma.compressed(us)
    vs = np.ma.masked_invalid(vs)
    vs = np.ma.compressed(vs)

    #check how many values are left in the region
    print('Values left for this step:')
    print(us.size)#; exit()
    if us.size < 100: continue
    
    #################################################################3
    #plot on overview map
    #reproject vertices
    xa, ya = ma(lons, lats)
    #reproject ship coordinates
    xla, yla = ma(ship_lon, ship_lat)
    cl = next(color)
    ax.plot(xa,ya,'.',color=cl, alpha=.3)
    #ship
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
    
    #ship
    gx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    #full mesh
    gx.triplot(pts[:,0], pts[:,1], tri.simplices.copy(), color='k', alpha=0.5, label='full')
    #colors for the nods
    colormesh=iter(plt.cm.jet_r(np.linspace(0,1,len(stp))))
    
    #otherwise it never gets saved...
    outname='overview_mesh_seed'+date1
    gx.legend()
    fig4.savefig(outpath+outname)
    plt.close(fig4)

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
    dummy_max_td_all = []
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
        #print(dummy_td)
        
        #print(dummy_div)
        #print(dummy_shr)
        #print(dummy_td) #little higher than div, because shr is low
        #exit()
        
        #now also make some estimates of what is the max deformation/displacement that we measure
        exag_fac = factor*100                                                        #40*100=4km total displacement!
        dummy_uvert = np.array([dxmean/diff,dxmean/diff,(dxmean+exag_fac)/diff])
                     
        dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_f=deformation(dummy_vert,dummy_uvert,dummy_vvert)
                
        dummy_div = dummy_a+dummy_b
        dummy_shr = .5*np.sqrt((dummy_a-dummy_d)**2+(dummy_b+dummy_c)**2)
        dummy_td = np.sqrt(dummy_div**2 + dummy_shr**2)
    
        dummy_max_td_all.append(dummy_td)
        #print(dummy_td)
        #exit()
        
    #store in a file
    tt = [dummy_ls_all, dummy_td_all, dummy_max_td_all]
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
    print(div)
    print(np.max(np.abs(div)),np.min(np.abs(div)))
    
    #use threshold and triangle size criteria to detect noise due to step function in speed
    #hessian filter got rid of the image artifacts and bad data in the shattered zones (MIZ etc), all large triangle left are of good qualities
    #step function artefacts are all small traingles
    #WARNING: still not all large triangles are part of this!
    threshold = ~( (td>dummy_td) | (area > (distance*1.1)**2/2) & (minang>minang_limit) )
    
    #for high resolution data only:
    #if increassing hessian mask to 8, we get larger triangles aroud the LKFs
    #threshold = ~((np.abs(div)>abs(dummy_div)) | (area > (distance*2.5)**2/2))
    
    #first check for MOSAiC
    threshold = td<dummy_td
    
        
    ##apply LKF filter
    #prepare place to store filtered data
    div_f2 = div.copy()
    shr_f2 = shr.copy()
    minang_f2 = minang.copy()   #this will be used later for coarsning
    
    if LKF_filter==True:
        print('starting LKF filter')
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
        
        #ship
        xl, yl = m(ship_lon, ship_lat)
        nx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

        
        patches_all = []
        for k in range(div.shape[0]):
            patch = Polygon(tripts[k,:,:])
            patches_all.append(patch)

        #plot filled triangles
        p = PatchCollection(patches_all, cmap=plt.cm.bwr, alpha=1)
        p.set_array(div*1e6)
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
            nmask = (n == -1) | threshold[n] | (minang[n] < minang_limit) 
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
                        nmask = (nn == -1) | threshold[nn] | used | (minang[nn] < minang_limit)
                        nn = np.ma.array(nn,mask=nmask); nn = np.ma.compressed(nn)
                        
                        #store values
                        for j in nn:
                            #check if we still did not reach the max size of kernel for this side
                            #and if we still have some unused neighbors
                            #3 is recommended max value by Buollion
                            #use < 3 to get to at least 3 (at least one more will be added after this)
                            if (len(side)<kernel) & (len(nn)>0):
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
        print('Test figure saved!')###########################################################################################################3
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
        
        #get LKF angles        
        #this only works at times when LKF are straight lines
        #it will not work during extensive leads/free drift, but again when the situation calms down again!!!
        if len(pindex) > 0:
            lkfs, buff, sline, angles = get_lkf_angle(tri,tripts,threshold,pindex)     
        
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
            
            #ship
            xl, yl = m(ship_lon, ship_lat)
            px.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
                        
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
                                    
            #plot LKF triangles over
            patches_p = []
            for k in pindex:
                patch = Polygon(tripts[k])
                patches_p.append(patch)
                
            p = PatchCollection(patches_p, ec= 'g', fc=None, alpha=1)
            px.add_collection(p)

            fig6.savefig(outpath+'test_lkf_lines_'+date1.split('T')[0]+file_name_end+'.png',bbox_inches='tight')
            plt.close(fig6)
            
            print('Done with LKF analysis for', date1.split('T')[0])
        else:
            angles=[]
            
        if low_mean==True:
            #under-threshold triangles means
            #the mean of calculated sub-t values is very small and messes up the scaling (steepens the power law)
            #instead we calculate the mean of the low (sub-threshold) values based on the comparison to the thresholded vs all ship radar data
            #this is based on emperical values from ship radar total deformaton for the whole period
            #we assume that the ratios between td, div and shr are same
            #WARNING: in some cases the estimated low mean can be above the threshold!!! >>> this has to be wrong!!!
            #WARNING: especially for small radius this low means can only be stimated once the deformation values are filtered by LKF filter


            a_ratio = 1.5310070655098624
            k_ratio =-0.14707155728256027

            #picking a short ls means high ratio=low td, and high standard deviation
            ls_low_mean = np.mean(np.sqrt(np.ma.array(area,mask=~threshold)))/1000  #power law is in km
            #ls_low_mean = np.mean(np.sqrt(area))/1000                               #mean over all triangles
            #ls_low_mean = np.mean(np.sqrt(np.ma.array(area,mask=threshold)))/1000  #mean over deforming triangles (those that count for #WARNING
            
            ratio = a_ratio*ls_low_mean**k_ratio
            print(ratio)    #since the non-deforming triangles are same big for all days, this ratio will not change actually!
            
            
            n1 = np.ma.compressed(np.ma.array(np.abs(div),mask=threshold)).size
            n2 = div.size - n1
            print(n1,n2)
            
            td = np.sqrt(div_f2**2 + shr_f2**2)
            td_high_mean = np.mean(np.ma.array(td,mask=threshold))
            td_mean = td_high_mean/ratio                                #ratio was estimated based on several days, here used on a single day
            td_low_mean = (td_mean*(n1+n2) - td_high_mean*n1 ) / n2     #maybe this would work if eare was large (and diverse) enough
            
            step_low = np.mean(np.ma.array(td,mask=td>dummy_td))
            
            
            a_std = 2.9243933850630908e-06
            k_std = -0.6661488984596151
            sample_std = a_std*ls_low_mean**k_std
            #sample = generate_lognormal_samples(td_mean, sample_std, n=n1+n2)
            sample = generate_lognormal_samples(td_mean, sample_std, n=n1*5+n2)    #WARNING: this is a hack to compensate to the large area covered by large tri 
            #sample = sample[:n1+n2]                                                #same hack, wee below
            
            #instead: compare the area in the deforming and undeforming triangles to estimate how many small ones would fit into deforming ones...
            
            #also size sample should be smaller for the days when deforming triangles are few - how to adjust for that???
            
            
            
            
            
            #distribution is a bit too spiky (several days)
            #this can be checked in sid_pl3.py
            
            #too spiky distribution could be due to too low sigma (shape parameter)
            #sigma could be too low due to too high mean - again the use of ratio is questionable!!! (see above)
            
            
            #sea ice deformation is multifractal, so all moments are equally important and it will only scale if all moments are there!
            #we could get the right moment values from the higher scales of SAR
            
            #Yes, Marsan states: if we randomly reshufle the velocity values in RGP, there is no more power law. this means that the spatial structure of the values is the one leading to power law. 
            #for us that means: only values that are correctly spatially aranged on the map will scale to power law
            #this also means: we can generate a random sample for a certain lenght scale, but it will never scale accross all different scales unless we distribute it properly in space!!!
            #this means that we can generate a map for one scale and compute averages
            #and yes, we can scale means and standard deviation for the time series!
            
            #this also means that we would need a construct a filter to spatially concentrate the deformation and mimic natural distribution... That is our plan for the snow topography...
            
            #sample_low=np.ma.array(sample,mask=sample>dummy_td)
            #sample_low=np.ma.compressed(sample_low)
            #print(sample_low.size)
            
            
            
            #sort the sample
            #print(sample)
            #print(sample.shape)
            sample = np.sort(sample)
            sample = sample[:n1+n2]                                                #part of same hack as above
            
            sample_low=sample[:n2]
            
            #this sample is ordered and so is the triangulation
            #therefore the maps are blobby (highest values where last triangles were drawn)
            #should we shuffle these values again???
            
            
            print('check the new low means - total deformation')        
            print('threshold:', dummy_td)
            print('high:', td_high_mean)
            print('est mean:', td_mean)
            print('sample mean:', np.mean(sample))
            print('step low:', np.mean(np.ma.array(td,mask=td>dummy_td)))
            print('est low:', np.mean(sample_low))  #these values should still be lower than threshold
            print('old mean:', np.mean(td))
            
            #td[threshold] = td_low_mean
            td[threshold] = sample_low
            
            ##sample_low has still many high values
            #print(sample_low[0])
            #print(sample_low[-1])
            #print(sample[-1])
            ##exit()

            
            
            #we should also spatially order the reamining high values in sample_low
            #triangles that are inside the buffer of the high deformation should get priority
            #or use triangles that are inside the buffer of the angles lines
            
            #get distances of all triangles (centroids) from the buffer/lines
            dist1,dist2 = get_distance_order(tri,tripts,pindex,lkfs)
            
            
            #sort sample from highest to lowest
            sample = np.array(sorted(sample,reverse=True))
            
            #td[threshold] and dist1,dist2 have same order
            #sample_low needs to be sorted in a same way
            #indexes that would sort the array from min to max: indexes to arrive to sorted state
            dist=np.array(dist1)
            idx = np.argsort(dist)
            #print(idx)
                        
            #we need indexes to arrive from sorted back to unsorted
            idxu = range(0,len(dist))
            idxuu = [s for _,s in sorted(zip(dist,idxu))]
                        
            sample_a=np.zeros_like(dist)
            sample_a[idxuu] = sample
            sample_low = sample_a[threshold]
            
            td[threshold] = sample_low
            
            ###angle lines
            #dist=np.array(dist2)
            #idx = np.argsort(dist)
            #print(idx)
                        
            ##we need indexes to arrive from sorted back to unsorted
            #idxu = range(0,len(dist))
            #idxuu = [s for _,s in sorted(zip(dist,idxu))]
                        
            #sample_a=np.zeros_like(dist)
            #sample_a[idxuu] = sample
            #sample_low = sample_a[threshold]
            
            #td[threshold] = sample_low

            
            
            
            
            
            
            print('new mean:', np.mean(td)) #this should be same as estimated mean!
            #exit()
            
            
            ##this is the problem!!!
            #ratio2 = np.mean(np.abs(div_f2))/np.mean(shr)
            ##print(ratio2)
            ##exit()
            #div_low_mean = td_low_mean /np.sqrt(2) * ratio2
            #if np.mean(div_f2) < 0:
                #div_low_mean = div_low_mean *-1
            #shr_low_mean = td_low_mean * (1-ratio2)
            
            #print('check the new low means - divergence')
            #print('old mean div:', np.mean(np.abs(div_f2)))
            #print('threshold:', dummy_div)
            #print('est low:', div_low_mean)
            
            
            
            #what if we simply replace the low values by the low ship radar data at this initial scale
            #we already know that the high values are co-inciding...
            #and generate a random sample from them
            
            #why do we actually need there maps for?
            #just to recreate a perfect power law?
            
            #to keep them for later to generate topography???
            
            #is our distribution really log-normal?
            #or is it logarithmic???
            
            #or what if we have to make log-normal distribution for real mean and then use the dummy to get only the low values???
            
            #div_f2[threshold] = td_low_mean/np.sqrt(2)
            
            
            ##can we have a randoly distributed values instead?
            ##log-normal distribution with known mean and standard deviation?
            #sample_mean = td_low_mean/np.sqrt(2)
            #sample_mean = np.mean(np.ma.array(div_f2,mask=~threshold))
            #sample_mean = np.mean(np.ma.array(np.abs(div_f2),mask=td>dummy_td))     #this mean should be much lower!!!
            #a_std = 1.6299983312588754e-07
            #k_std = -0.9780933554534802
            #sample_std = a_std*ls_low_mean**k_std
            #sample = generate_lognormal_samples(sample_mean, sample_std, n=n2)
            
            sample_div = sample_low/np.sqrt(2)
            
            #then also randomly multiply by -1???
            #for example: randrange(0,1), then change all 0 to -1
            sign = np.random.random_integers(0, 1, size=n2)
            sign = np.where(sign==0,-1,sign)
            print(sign)
            sample_div = sample_div*sign
            #exit()
            
            div_f2[threshold] = sample_div
            print('new mean div:', np.mean(np.abs(div_f2)))
                
            print('check the new low means - shear')        
            print('old mean shr:', np.mean(shr_f2))
            print('threshold:', dummy_shr)
            #print('est low:', shr_low_mean)
            
            #shr_f2[threshold] = td_low_mean/np.sqrt(2)
            sample_shr = sample_low/np.sqrt(2)
            shr_f2[threshold] = sample_shr
            print('new mean shr:', np.mean(shr_f2))  
            
            test = np.sqrt(np.mean(shr_f2)**2+np.mean(div_f2)**2)       #this ones comes quite low...
            print(test)
            
            #exit()

    
    #course-grain these results:
    ddidx = 0
    for j in stp:
        print('Step: '+str(j))
        outname_td = 'td_seed_f_'+reg+'_L'+str(j)+file_name_end
        outname_ts = 'ts_seed_f_'+reg+'_L'+str(j)+file_name_end
        outname_angles = 'angle_seed_f_'+reg+'_L'+str(j)+file_name_end
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
        angles_list=[]
        
        
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
            
            ##add corner points for higher steps
            #if j > 5:
                #xlist = [x[idcr1],x[idcr2],x[idcr3],x[idcr4]]
                #ylist = [y[idcr1],y[idcr2],y[idcr3],y[idcr4]]
                #xs = np.append(xs,xlist)
                #ys = np.append(ys,ylist)

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
            
            #WARNING: these are triangles of various sizes at LKFs, in MIZ and at the scene edges
            #allow max 2x larger triangles than optimal for smallest scale for seeding
            max_area = (distance*2)**2/2
            #do area-weighted coarse-graining
            div_seed,shr_seed,area_seed,minang_seed = coarse_grain(tripts,tripts_seed,div_f2,shr_f2,max_area,minang_f2)
            
            td_seed = np.sqrt(div_seed**2 + shr_seed**2)
            
            if low_mean==False:
                #mask with threshold that depends on the lenght scale (mask out large triangles)
                #this is done in the plotting, so not necessary here...
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
                #id_seed = np.ma.array(id_seed,mask=threshold_seed)
            
        else:
            #if LKF_filter==False:
                #lkf_id = np.zeros_like(area)
            
            if low_mean==False:
                div_seed=np.ma.array(div_f2, mask=threshold)
                shr_seed=np.ma.array(shr_f2, mask=threshold)
                tripts_seed=tripts
                area_seed = np.ma.array(area, mask=threshold)
                minang_seed = np.ma.array(minang, mask=threshold)
                #id_seed = np.ma.array(lkf_id,mask=threshold)
                                
            else:
                #keep all the data
                div_seed=div_f2
                shr_seed=shr_f2
                tripts_seed=tripts
                area_seed = area
                minang_seed = minang
                #id_seed = lkf_id
            
        
        ####################################################################3
        if image == True:
            ##Plotting deformation 
            deform = div_seed*1e6
            outname = outpath+'map_div_seed_f_'+reg+'_L'+str(j)+'_'+date1+'_'+date2
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
            plot_def(area_def,tripts_seed,deform,outname,label,interval,cmap,ship_lon,ship_lat,radius)
        #####################################################################3
            #exit()

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

        
        
        #at high resolution deformation data we cant avoid of having two separate clouds (detectable and un-detectable deformation) - that happens when our method cant detect deformation anymore
        
        #this can not be fixed by filtering and the line in power law will break even if we use the Sylvain filter
        
        #still the filter removes the artifical opening-clossing mess in shear zones (and removes some high values in that zones), so it makes it useful for the mapping of larger features such as LKFs and leads.
        



        #storing data for the scatter plots
        
        td = np.sqrt(div_seed**2 + shr_seed**2)
        ls = np.sqrt(area_seed)
        minang = np.array(minang_seed)
        #id_num = np.array(id_seed)
        
        
        td_list.extend(np.ma.compressed(td).tolist())
        ls_list.extend(np.ma.compressed(ls).tolist())
        ang_list.extend(np.ma.compressed(minang).tolist())
        #id_list.extend(np.ma.compressed(id_num).tolist())
        
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
        tt = [date_ts, mpdiv, mndiv, mdiv, mshr, mpdiv_sd, mndiv_sd, mdiv_sd, mshr_sd]
        table = zip(*tt)
        table = list(zip(*tt))

        output = outpath_def + outname_ts
        with open(output, 'ab') as f:
            #header
            #f.write(b'date, pos. divergence, neg. divergence, mean divergence, mean shear\n')
            np.savetxt(f, table, fmt="%s", delimiter=",")
        
        #angle data
        if LKF_filter:
            angles_list.extend(angles)
            #print(angles_list)
            
            #write data for angle PDFs
            output = outpath_def + outname_angles
            with open(output, 'ab') as f:
                np.savetxt(f, angles_list, fmt="%s", delimiter=",")
        
    outname='overview_mesh_seed'+date1
    gx.legend()
    fig4.savefig(outpath+outname)
    plt.close(fig4)
    
            
fig1.savefig(outpath+'overview_map_'+reg+'_'+str(int(radius/1000.)),bbox_inches='tight')
plt.close(fig1)
