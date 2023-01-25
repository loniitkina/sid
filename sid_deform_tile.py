import os
from glob import glob
from datetime import datetime
import numpy as np
import pyresample as pr
from pyproj import Proj, transform
from pyresample.geometry import AreaDefinition
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from shapely.geometry import Point, MultiPoint, MultiPolygon
from shapely.geometry import Polygon as Shapely_Polygon
from shapely.ops import unary_union
import pickle
from sid_func import * 
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

#supress warnings
import warnings
warnings.filterwarnings("ignore")

#WARNING: shapely can be unstable, lots of segfaults
#try: pip uninstall shapely; pip install --no-binary :all: shapely
#works better, but not good

#Conda will have problems with some libraries installed by pip (if shared resources), better if all installed by conda
#try: pip uninstall shapely; conda install shapely
#no, that does not work use this instead: https://anaconda.org/conda-forge/shapely
#does not work :(

#Do you want all output in figures?
image=True
#image=False

#active floe size
afs=True
#afs=False

##for parcel tracking we need to have consequtive data: second scene in pair needs to be first scene in the next pair! (combo option is not possible here)
#just save level 1 data and exit
parcel=True
#parcel=False

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
minang_limit=5  #this will leave out some triangles over wide LKFs (dense points at the edges of the feature make narrow triangles!)
minang_limit=2

#select area size
radius = 120000
file_name_end = '_120km'

##gaps between tiles at _120km
#radius = 150000
#file_name_end = '_150km'

#expected divergence values (div in s-1 * 10e6)
interval = [-1, 1]

#for scaling
#create log-spaced vector and convert it to integers
n=9 # number of samples
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)
#stp=[1]

file_name_end = file_name_end+'.csv'

#-------------------------------------------------------------------INPUT
#threshold detection values: step and factor differ for each for each input data (usually put in different folders/paths)
#resolution=factor X step
#at factor=1 and step=10, we get 400m distance between grid points
#distance between nods==400m or triangle area is 400*400/2 
#there is one nod, where displacement difference is 40m more than at the other two
#example velocities nod1=1000m/diff, nod2=1000m/diff, nod3=1040m/diff
#exaggerate the factor as the possible displacements are discrete steps as 40, 80, 120 m...
#diff is different for every scene pair and threshold needs to be recalculated

extra_margin = 20   #some allowance #20 in test5 gave best results so far (but also lets in the last artifacts in LKF filter)
#20 was a good margin for factor=80, seems like 25% was a good estimate
#no margin necessary for scaling, but use 10 to get clean LKFs for other analysis where 'less is more' (e.g. parcels)
extra_margin = 30   #this is a good option for 'less is more': removed most of the artifacts (and some data) should work well for parcels and AFS!

#800 m resolution data (pixel/factor=80m, step=10)
inpath = '/scratch/pit000/results/sid/drift/stp10_factor05/'
inpath = '/scratch/pit000/results/sid/drift/stp10_factor05_200km/'
step = 10      #grid resolution of drift 
factor = 80    #pixel size (factor in sid_drift, default factor in FT is 0.5, which translates to 80 m here)

#200 m resolution data
#inpath = '/scratch/pit000/results/sid/drift/stp5_factor1/'
#step = 5
#factor = 40

#------------------------------------------------------------------OUTPUT
#outpath_def = '/scratch/pit000/results/sid/afs/'
#outpath_def = '/scratch/pit000/results/sid/parcels/'

outpath_def = '/scratch/pit000/results/sid/deform200km/'


#parcels 
#inpath = '../sidrift/data/stp10_parcels_f1/'
#inpath = '../sidrift/data/stp10_parcels_f1_rs2/'
#outpath_def = inpath

#afs 
#outpath_def = '../sidrift/data/stp10_asf/'
#canberra
#outpath_def = 'data/stp10_afs/'

outpath = '/scratch/pit000/results/sid/plots/'
outpath = '/scratch/pit000/results/sid/plots200km/'

#central tile - for mapping/projection
#shipfile = '../../downloads/position_leg3_nh-track.csv'	#leg3 (and transition to leg 4 until 6 June)
shipfile = '../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0.csv'
#shipfile = '../../downloads/lance_leg1.csv'
reg = shipfile.split('_')[1]
proj = 'ship'

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
#WARNING: this is still manually adjusted!
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

#remove all the output text files from previous script runs (or they will be appended)
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


##should we instead
#-loop through all dates
#-loop throu all regions
#-make daily tiles

#regions have to be arranged CW  - or convex hull of frames wont work!
#start with c and continue with most likely good coverage (e.g. s)
regions=['c','sw','s','se','e','ne','n','nw','w']
regions=['c','s','w','e','n','sw','se','nw','ne']

#WARNING - this is also hard-coded....
#check how many days there are in the central tile/main track file (the one with the ship)
#main_trackfile='../../downloads/position_leg3_nh-track_c_200km-fnames.csv'
main_trackfile='../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0_c_200km-fnames.csv'
#main_trackfile='../../downloads/lance_leg1_c_200km-fnames.csv'

noons = getColumn(main_trackfile,0,header=False)
noons = [ datetime.strptime(noons[i], "%Y-%m-%d %H:%M:%S") for i in range(len(noons)) ]
days = [ datetime.strftime(noons[i], "%Y%m%d") for i in range(len(noons)) ]
print(days)

#colors for overview map
color=iter(plt.cm.jet_r(np.linspace(0,1,len(days)+1)))

#days=['20200401','20200402','20200403','20200511']
#days=['20200402']
#days=['20200323']
#days=['20200307']
#days=['20200414']

for day in days:
    print(day,'******************************************************************************************************************************')
    
    #color for the overview map
    cl = next(color)
    
    #check if any data was collected - important for saving at the end of this loop
    anydata=False
    
    #keep all the timestamps and time differences
    timestamp1=[]
    timestamp2=[]
    timediff=[]

    for region in regions:
        print(region)
        try:
            fname=sorted(glob(inpath+'SeaIceDrift_'+day+'*_'+region+'.npz'))[0]
            print(fname)
        except:
            fname='empty'
        #no point in going on if there is nothing in the central tile   ==== central file should always be there: reprocess!!!!!!!
        if region=='c' and fname=='empty': print('no central tile');break
        if fname=='empty': 
            print(region,' - no tile',day)
            timestamp1.append('notile');timestamp2.append('notile');timediff.append('notile')
            continue
            
        #read in all the data
        container = np.load(fname)
        u = container['upm']     
        v = container['vpm']
        hpm = container['hpm']    #hessian
        lat = container['lat1']
        lon = container['lon1']
        print('Size of input matrix:', u.size)
            
        #get time difference
        date1 = fname.split('_')[-3]
        date2 = fname.split('_')[-2].split('.')[0]
        dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
        dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
        diff = (dt2-dt1).seconds + (dt2-dt1).days*24*60*60
        
        timestamp1.append(date1)
        timestamp2.append(date2)
        timediff.append(diff)
        
        #tile center - from track file
        shipfile_tile = shipfile.split('.csv')[0]+'_'+region+'_200km.csv'
        print(shipfile_tile)
        mettime = getColumn(shipfile_tile,0)
        dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
        if dtb[0]>dt1: continue
        if dtb[-1]<dt1: continue
        mi = np.argmin(abs(np.asarray(dtb)-dt1))
        ship_lon = np.asarray(getColumn(shipfile_tile,1),dtype=float)[mi]
        ship_lat = np.asarray(getColumn(shipfile_tile,2),dtype=float)[mi]
        if np.isnan(ship_lon): continue
        print('center at: ',ship_lon,ship_lat,dtb[mi])

        #ship centered projection
        #work in geographical projection, where units are meters!
        #also keep the dates from the central tile >> parcels and AFS
        if region=='c':
            
            xlp,ylp = transform(inProj,outProj,ship_lon,ship_lat)
            #Using a projection dictionary (same stuff as for overview map)
            width = radius*2/factor # m spacing
            height = radius*2/factor # m spacing
            area_extent = (xlp-radius,ylp-radius,xlp+radius,ylp+radius)
            area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
            m = pr.plot.area_def2basemap(area_def)
            
            date1_c=date1
            date2_c=date2
            
        #reproject center position
        xl, yl = m(ship_lon, ship_lat)

        #reproject vertices
        x, y = m(lon, lat)
        
        #cut out region
        mask = ~((x<xl-radius) | (x>xl+radius) | (y<yl-radius) | (y>yl+radius))
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
        
        #find corners of valid data
        idcr1 = np.unravel_index(np.argmin(np.sqrt((x-(xl-radius))**2+(y-(yl-radius))**2)),x.shape)
        idcr2 = np.unravel_index(np.argmin(np.sqrt((x-(xl-radius))**2+(y-(yl+radius))**2)),x.shape)
        idcr3 = np.unravel_index(np.argmin(np.sqrt((x-(xl+radius))**2+(y-(yl+radius))**2)),x.shape)
        idcr4 = np.unravel_index(np.argmin(np.sqrt((x-(xl+radius))**2+(y-(yl-radius))**2)),x.shape)
        
        #get rid of all the nans, the resulting arrays are flattened
        xs = np.ma.masked_where(~(np.isfinite(us)),xs)
        xs = np.ma.compressed(xs)
        ys = np.ma.masked_where(~(np.isfinite(us)),ys)
        ys = np.ma.compressed(ys)
        lons = (np.ma.masked_where(~(np.isfinite(us)),lons))
        lats = (np.ma.masked_where(~(np.isfinite(us)),lats))
        us = np.ma.masked_invalid(us)
        us = np.ma.compressed(us)
        vs = np.ma.masked_invalid(vs)
        vs = np.ma.compressed(vs)

        #check how many values are left in the region
        print('Values left for this step:')
        print(us.size)
        if region=='c' and us.size < 100: #if no data in the central tile - stop here
            print('empty central tile');break
        elif us.size < 100:
            continue            
        
        #################################################################3
        #plot on overview map
        #reproject vertices
        xa, ya = ma(lons, lats)
        #reproject ship coordinates
        xla, yla = ma(ship_lon, ship_lat)
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
        if region=='c':
            #mesh test plot
            fig4    = plt.figure(figsize=(10,10))
            gx      = fig4.add_subplot(111)
            print('gx created')
    
        #center
        gx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
        #full mesh
        gx.triplot(pts[:,0], pts[:,1], tri.simplices.copy(), color='k', alpha=0.5, label='full')
        #colors for the nods
        colormesh=iter(plt.cm.jet_r(np.linspace(0,1,len(stp))))

        ##################################################################3

        #define threshold value
        #due to different time difference between pairs this is going to be different for every tile
        
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
            distance = dst*j
            exag_fac = factor+extra_margin
            
            dummy_vert = np.array([[0,0],[distance,0],[0,distance]])
            dummy_uvert = np.array([dxmean/diff,dxmean/diff,(dxmean+exag_fac)/diff])
            dummy_vvert = np.array([dymean/diff,dymean/diff,dymean/diff])             
            dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_f=deformation(dummy_vert,dummy_uvert,dummy_vvert)
            
            #because this is only one nod, there will be divergence (increase in area) and shearing (change in shape, angles)
            dummy_div = dummy_a+dummy_b
            dummy_shr = .5*np.sqrt((dummy_a-dummy_d)**2+(dummy_b+dummy_c)**2)
            dummy_td = np.sqrt(dummy_div**2 + dummy_shr**2)
            dummy_ls = np.sqrt(dummy_f)
        
            #now also make some estimates of what is the max deformation/displacement that we measure
            exag_fac = factor*100                                                        #40*100=4km total displacement!
            dummy_uvert = np.array([dxmean/diff,dxmean/diff,(dxmean+exag_fac)/diff])
            dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_f=deformation(dummy_vert,dummy_uvert,dummy_vvert)
            dummy_div = dummy_a+dummy_b
            dummy_shr = .5*np.sqrt((dummy_a-dummy_d)**2+(dummy_b+dummy_c)**2)
            dummy_td = np.sqrt(dummy_div**2 + dummy_shr**2)
            
            #store in a list
            dummy_td_all.append(dummy_td)
            dummy_ls_all.append(dummy_ls)
            dummy_max_td_all.append(dummy_td)
            
        #store in a file
        tt = [dummy_ls_all, dummy_td_all, dummy_max_td_all]
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
        #print(dummy_div)
        
        dummy_shr = .5*np.sqrt((dummy_a-dummy_d)**2+(dummy_b+dummy_c)**2)
        dummy_td = np.sqrt(dummy_div**2 + dummy_shr**2)
        #print(dummy_td)
        
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
        dux=[];duy=[];dvx=[];dvy=[];minang=[];area=[];utri=[];vtri=[]
        for t in range(0,len(tripts)):
            vert = np.asarray(tripts[t])
            uvert = upts[t]
            vvert = vpts[t]
        
            #deformation:
            a,b,c,d,e,f=deformation(vert,uvert,vvert)
            
            #mean drift (for parcels):
            uu = np.mean(uvert)
            vv = np.mean(vvert)
            
            #store
            dux.append(a);duy.append(b);dvx.append(c);dvy.append(d);minang.append(e);area.append(f);utri.append(uu);vtri.append(vv)
            
        dux = np.array(dux)
        duy = np.array(duy)
        dvx = np.array(dvx)
        dvy = np.array(dvy)
        minang = np.array(minang)
        area = np.array(area)
        
        #triangle averaged drift and time difference information for the parcels
        utri = np.array(utri)
        vtri = np.array(vtri)
        dt = np.ones_like(utri)*diff
            
        div = dux + dvy
        shr = .5*np.sqrt((dux-dvy)**2+(duy+dvx)**2)
        td = np.sqrt(div**2 + shr**2)
        #print(div)
        print('Dummy value for divergence: ',dummy_div)
        print('Max and Min divergence values: ',np.max(np.abs(div)),np.min(np.abs(div)))
        
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
        threshold = (td<dummy_td) | (minang<minang_limit) #| (area > (distance*5)**2/2)
        
            
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
            if region=='c':
                fig5    = plt.figure(figsize=(20,10))
                
                #plot original divergence field
                nx      = fig5.add_subplot(121)
                mn = pr.plot.area_def2basemap(area_def)
                mn.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
                mn.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
                
                #set up the plot for the flitered field
                lx      = fig5.add_subplot(122)
                ml = pr.plot.area_def2basemap(area_def)
                ml.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
                ml.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])        

                #ship
                xl, yl = mn(ship_lon, ship_lat)
                nx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
                lx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
                
            else:
                #tile center
                xl, yl = mn(ship_lon, ship_lat)
                nx.plot(xl,yl,'o',markeredgewidth=2,color='hotpink',markersize=10,markeredgecolor='k')
                
            

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
                    
                    ###plot LKF###########################################################################################################3   
                    
                    ##patches = []
                    ###for k in lkf_idx:
                    ##patch = Polygon(tripts[p])
                    ##patches.append(patch)
                    
                    ####plot filled triangles
                    ##pc = PatchCollection(patches, cmap=plt.cm.bwr, alpha=1, edgecolor='k')
                    ##pc.set_array(np.ones((1))*lkf_mdiv*1e6) #has to be an array
                    ##pc.set_clim(interval)
                    ##lx.add_collection(pc)
                    
            print('LKF filter completed')    
                    
            #update pindex
            pindex = np.arange(0,tripts.shape[0]);pindex = np.ma.array(pindex, mask=threshold); pindex = np.ma.compressed(pindex)

            print('Start the tiling')
            #a region for all area where there is data (following image pair edges)
            #get triangle centroids
            ctrdx = np.zeros_like(div);ctrdy = np.zeros_like(div)
            for p in range(0,tripts.shape[0]):
                ctrdx[p],ctrdy[p] = centroid(tripts[p])
                
            #convert to format that shapely can work with
            ctr_nods = [ np.array([i,j]) for i,j in zip(ctrdx,ctrdy) ]

            print('Get initial tile')
            tile = MultiPoint(ctr_nods).convex_hull
            
            #get rid of large triangles at the tile edge - artefacts - 20km
            #centroid has to lay inside this buffer polygon
            good_nods=[]
            eroded=tile.buffer(-20000)    #this erodes the polygon by buffer width
            for i in range(0,len(ctr_nods)):
                if eroded.disjoint(Point(ctr_nods[i])) and ( area[i] > (distance*2)**2/2 or minang[i] < 10 ):
                    good_nods.append(False)
                else:
                    good_nods.append(True)
            
            ctrdx = ctrdx[good_nods]
            ctrdy = ctrdy[good_nods]
            tripts = tripts[good_nods]
            div_f2 = div_f2[good_nods]
            shr_f2 = shr_f2[good_nods]
            minang_f2 = minang_f2[good_nods]
            threshold = threshold[good_nods]
            utri = utri[good_nods]
            vtri = vtri[good_nods]
            dt = dt[good_nods]
            
            print('Get refined tile')
            ctr_nods = [ np.array([i,j]) for i,j in zip(ctrdx,ctrdy) ]
            tile = MultiPoint(ctr_nods).convex_hull        
            
            #get rid of potential lines and points
            if tile.geom_type == 'GeometryCollection':
                print(tile.geom_type)
                tmp=[]
                for geom in tile.geoms:
                    if geom.geom_type == 'Polygon': #if produces by one image pair, there should always be just one polygon
                        tmp.append(geom)
                if len(tmp)>1:
                    tile = MultiPolygon(tmp)
                elif len(tmp)==0:
                    continue
                else:
                    tile=tmp[0]
            
            print(tile.geom_type)
                
            #in case this is just a single point or line (still not Polygon or MultiPolygon), there is no point in continuing here
            if tile.geom_type != 'Polygon':# or tile.geom_type != 'MultiPolygon':
                print('No useful tile polygons in ',region)
                if region=='c':
                    print('empty central tile');break
                else:    
                    continue
            
            #store and plot
            if region=='c':
                #store all nods
                tiled_cover=tile
                anydata=True
                
                tiled_tripts=tripts.copy()
                tiled_div_f2=div_f2.copy()
                tiled_shr_f2=shr_f2.copy()
                tiled_minang_f2=minang_f2.copy()
                tiled_threshold=threshold.copy()
                tiled_utri=utri.copy()
                tiled_vtri=vtri.copy()
                tiled_dt=dt.copy()
                
                #plot the area covered
                xg, yg = tiled_cover.exterior.xy 
                lx.fill(xg, yg, fc='none', ec='green')
                
                #get plot frame polygon
                xlist = [xl-radius,xl-radius,xl+radius,xl+radius,xl-radius]
                ylist = [yl-radius,yl+radius,yl+radius,yl-radius,yl-radius]
                corner_nods = [ np.array([i,j]) for i,j in zip(xlist,ylist) ]
                #plot_corners = Shapely_Polygon(corner_nods)
                #print(plot_corners)

            else:
                #add new nods
                #difference: a.difference(b) is what is covered by a, but not by b
                diff_tile = tile.difference(tiled_cover)
                new_nods=[]
                for i in range(0,len(ctr_nods)):
                    if diff_tile.contains(Point(ctr_nods[i])):
                        new_nods.append(True)
                    else:
                        new_nods.append(False)

                tiled_tripts=np.concatenate((tiled_tripts, tripts[new_nods]), axis=0)
                tiled_div_f2=np.concatenate((tiled_div_f2, div_f2[new_nods]), axis=0)
                tiled_shr_f2=np.concatenate((tiled_shr_f2, shr_f2[new_nods]), axis=0)
                tiled_minang=np.concatenate((tiled_minang_f2, minang_f2[new_nods]), axis=0)
                tiled_threshold=np.concatenate((tiled_threshold,threshold[new_nods]), axis=0)
                tiled_utri=np.concatenate((tiled_utri,utri[new_nods]), axis=0)
                tiled_vtri=np.concatenate((tiled_vtri,vtri[new_nods]), axis=0)
                tiled_dt=np.concatenate((tiled_dt,dt[new_nods]), axis=0)
                
                #plot the new area covered and extend coverage - dont do convex hull yet or intenral area are filled in!
                if diff_tile.geom_type == 'MultiPolygon':
                    for geom in diff_tile.geoms: 
                        
                        #plot
                        xg, yg = geom.exterior.xy 
                        lx.fill(xg, yg, fc='none', ec='green')
                        
                    #extend
                    tiled_cover=tiled_cover.union(diff_tile)
                
                #sometimes the difference may be also point and lines and not just polygons
                elif diff_tile.geom_type == 'GeometryCollection':
                    for geom in diff_tile.geoms:
                        if geom.geom_type == 'Polygon':
                            #plot
                            xg, yg = geom.exterior.xy 
                            lx.fill(xg, yg, fc='none', ec='green')
                            
                            #extend
                            tiled_cover=tiled_cover.union(geom)
                
                else:
                    #plot
                    xg, yg = diff_tile.exterior.xy 
                    lx.fill(xg, yg, fc='none', ec='green')
                    
                    #extend
                    tiled_cover=tiled_cover.union(diff_tile)
                    
                
    #store combined data and figures for all regions/tiles
    if anydata:
        print('We have data for date: ',day)
        #update pindex (triangle index) for afs and lkf_id
        pindex = np.arange(0,tiled_tripts.shape[0]);pindex = np.ma.array(pindex, mask=tiled_threshold); pindex = np.ma.compressed(pindex)
        
        ##tiles can go beyond the mapped region (plot extent)
        ##limit the cover by corner coordinates
        #tiled_cover_frame=tiled_cover.intersection(plot_corners)
        #print(tiled_cover_frame)
        
        #store all pindex/divergence/shear/tripts/area/minang data
        if afs:
            #dump active floe size input data into numpy file
            out_file = outpath_def+'Afs_'+date1_c+'_'+date2_c+file_name_end+'_tiled.npz'
            np.savez(out_file,pindex = pindex, tripts = tiled_tripts, region = tiled_cover, corners = corner_nods)

            #for maping
            out_area_def = outpath_def+'area_def_'+date1_c+'_'+date2_c+file_name_end+'.pkl'
            f = open(out_area_def,'wb')
            pickle.dump(area_def,f)
            f.close()

            print('Storing data: ',out_file)
            #continue
        
        #store all this for parcel tracking!
        if parcel:
            #find in triangle centroid
            ctrdx = np.zeros_like(tiled_div_f2);ctrdy = np.zeros_like(tiled_div_f2) 
            for i in range(0,len(tiled_tripts)):
                ctrdx[i],ctrdy[i] = centroid(tiled_tripts[i])
                
            #convert to lat,lon
            ctrd_lon,ctrd_lat = m(ctrdx,ctrdy,inverse=True)
            
            #get 0,1 array with all damaged triangles as 1
            damage = np.where(~tiled_threshold,1,0)
            lead = np.where((tiled_div_f2>0)&(~tiled_threshold),1,0)
            ridge = np.where((tiled_div_f2<0)&(~tiled_threshold),1,0)
            
            #how can i track the area change over this triangle???
            #keep triangulation but apply to lon2,lat2
            #final coordinates can be stored by drift script or calculated from displacements
            
            #dump damage & drift data into numpy file
            out_file = outpath_def+'Damage_'+date1_c+'_'+date2_c+'_tiled.npz'
            np.savez(out_file,lon = ctrd_lon,lat = ctrd_lat, d = damage, l=lead, r=ridge, u=tiled_utri, v=tiled_vtri, dt=tiled_dt)

            print('Storing data: ',out_file)
            #continue

        #store regional timestamps and differences
        tt = [regions,timestamp1,timestamp2,timediff]
        table = list(zip(*tt))
        outname=outpath_def+'timestamps_'+day+'_'+str(int(radius/1000))+'km.csv'
        print(outname)
        with open(outname, 'wb') as f:
                np.savetxt(f, table, fmt="%s", delimiter=",")

        
        #figures
        outname='overview_mesh_'+day
        gx.legend()
        fig4.savefig(outpath+outname)
        plt.close(fig4)

        #plot all tile coverage
        #here it is safe to do convex hull - no more new regions are coming
        tiled_cover=unary_union(tiled_cover)
        if tiled_cover.geom_type == 'MultiPolygon':
            for geom in tiled_cover.geoms:
                xg, yg = geom.exterior.xy 
                lx.fill(xg, yg, fc='none', ec='purple', ls='--')
        else:
            xg, yg = tiled_cover.exterior.xy 
            lx.fill(xg, yg, fc='none', ec='purple', ls='--')
        
        #plot tiled and filterd data   
        patches = []
        for p in pindex:
            patch = Polygon(tiled_tripts[p])
            patches.append(patch)
        
        ##plot filled triangles
        pc = PatchCollection(patches, cmap=plt.cm.bwr, alpha=1, edgecolor='k')
        pc.set_array(tiled_div_f2[~tiled_threshold])
        pc.set_clim(interval)
        lx.add_collection(pc)

        fig5.savefig(outpath+'filter_check'+date1_c.split('T')[0]+file_name_end+'.png',bbox_inches='tight')
        plt.close(fig5)
        print('Filter CHECK figure saved!: ',outpath+'filter_check'+date1.split('T')[0]+file_name_end+'.png')
        #exit()
    ###########################################################################################################################################3

    
    
    
    


#end-of-functional-code
        
        ##get LKF angles        
        ##this only works at times when LKF are straight lines
        ##it will not work during extensive leads/free drift, but again when the situation calms down again!!!
        #if len(pindex) > 0:
            #lkfs, buff, sline, angles = get_lkf_angle(tri,tripts,threshold,pindex)     
        
            ##Plotting
            #fig6    = plt.figure(figsize=(10,10))
            #px      = fig6.add_subplot(111)
            
            ##map area definition
            ##area_def_file = glob(inpath+'area_def_'+date+'*'+file_name_end+'*.pkl')[0]
            ##with open(area_def_file, "rb") as pkl:
                ##area_def = pickle.load(pkl)
            #m = pr.plot.area_def2basemap(area_def)
            
            #m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
            #m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
            
            ##ship
            #xl, yl = m(ship_lon, ship_lat)
            #px.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
                        
            ##plot the lkfs
            #for geom in lkfs:  
                #xg, yg = geom.xy  
                #px.plot(xg, yg, c='b')
                
            ##plot the split line 
            #xg, yg = sline.xy  
            #px.plot(xg, yg, c='g')
            
            ##plot the buffers (typically multipolygons)
            #if buff.geom_type == 'Polygon':
                #xg, yg = buff.exterior.xy 
                #px.fill(xg, yg, alpha=1, fc='none', ec='r')
            #if buff.geom_type == 'MultiPolygon':
                #for geom in buff.geoms:  
                    #xg, yg = geom.exterior.xy    
                    #px.fill(xg, yg, alpha=1, fc='none', ec='r')
                                    
            ##plot LKF triangles over
            #patches_p = []
            #for k in pindex:
                #patch = Polygon(tripts[k])
                #patches_p.append(patch)
                
            #p = PatchCollection(patches_p, ec= 'g', fc=None, alpha=1)
            #px.add_collection(p)

            #fig6.savefig(outpath+'test_lkf_lines_'+date1.split('T')[0]+file_name_end+'.png',bbox_inches='tight')
            #plt.close(fig6)
            
            #print('Done with LKF analysis for', date1.split('T')[0])
        #else:
            #angles=[]
            
        #if low_mean==True:
            ##under-threshold triangles means
            ##the mean of calculated sub-t values is very small and messes up the scaling (steepens the power law)
            ##instead we calculate the mean of the low (sub-threshold) values based on the comparison to the thresholded vs all ship radar data
            ##this is based on emperical values from ship radar total deformaton for the whole period
            ##we assume that the ratios between td, div and shr are same
            ##WARNING: in some cases the estimated low mean can be above the threshold!!! >>> this has to be wrong!!!
            ##WARNING: especially for small radius this low means can only be stimated once the deformation values are filtered by LKF filter


            #a_ratio = 1.5310070655098624
            #k_ratio =-0.14707155728256027

            ##picking a short ls means high ratio=low td, and high standard deviation
            #ls_low_mean = np.mean(np.sqrt(np.ma.array(area,mask=~threshold)))/1000  #power law is in km
            ##ls_low_mean = np.mean(np.sqrt(area))/1000                               #mean over all triangles
            ##ls_low_mean = np.mean(np.sqrt(np.ma.array(area,mask=threshold)))/1000  #mean over deforming triangles (those that count for #WARNING
            
            #ratio = a_ratio*ls_low_mean**k_ratio
            #print(ratio)    #since the non-deforming triangles are same big for all days, this ratio will not change actually!
            
            
            #n1 = np.ma.compressed(np.ma.array(np.abs(div),mask=threshold)).size
            #n2 = div.size - n1
            #print(n1,n2)
            
            #td = np.sqrt(div_f2**2 + shr_f2**2)
            #td_high_mean = np.mean(np.ma.array(td,mask=threshold))
            #td_mean = td_high_mean/ratio                                #ratio was estimated based on several days, here used on a single day
            #td_low_mean = (td_mean*(n1+n2) - td_high_mean*n1 ) / n2     #maybe this would work if eare was large (and diverse) enough
            
            #step_low = np.mean(np.ma.array(td,mask=td>dummy_td))
            
            
            #a_std = 2.9243933850630908e-06
            #k_std = -0.6661488984596151
            #sample_std = a_std*ls_low_mean**k_std
            ##sample = generate_lognormal_samples(td_mean, sample_std, n=n1+n2)
            #sample = generate_lognormal_samples(td_mean, sample_std, n=n1*5+n2)    #WARNING: this is a hack to compensate to the large area covered by large tri 
            ##sample = sample[:n1+n2]                                                #same hack, wee below
            
            ##instead: compare the area in the deforming and undeforming triangles to estimate how many small ones would fit into deforming ones...
            
            ##also size sample should be smaller for the days when deforming triangles are few - how to adjust for that???
            
            
            
            
            
            ##distribution is a bit too spiky (several days)
            ##this can be checked in sid_pl3.py
            
            ##too spiky distribution could be due to too low sigma (shape parameter)
            ##sigma could be too low due to too high mean - again the use of ratio is questionable!!! (see above)
            
            
            ##sea ice deformation is multifractal, so all moments are equally important and it will only scale if all moments are there!
            ##we could get the right moment values from the higher scales of SAR
            
            ##Yes, Marsan states: if we randomly reshufle the velocity values in RGP, there is no more power law. this means that the spatial structure of the values is the one leading to power law. 
            ##for us that means: only values that are correctly spatially aranged on the map will scale to power law
            ##this also means: we can generate a random sample for a certain lenght scale, but it will never scale accross all different scales unless we distribute it properly in space!!!
            ##this means that we can generate a map for one scale and compute averages
            ##and yes, we can scale means and standard deviation for the time series!
            
            ##this also means that we would need a construct a filter to spatially concentrate the deformation and mimic natural distribution... That is our plan for the snow topography...
            
            ##sample_low=np.ma.array(sample,mask=sample>dummy_td)
            ##sample_low=np.ma.compressed(sample_low)
            ##print(sample_low.size)
            
            
            
            ##sort the sample
            ##print(sample)
            ##print(sample.shape)
            #sample = np.sort(sample)
            #sample = sample[:n1+n2]                                                #part of same hack as above
            
            #sample_low=sample[:n2]
            
            ##this sample is ordered and so is the triangulation
            ##therefore the maps are blobby (highest values where last triangles were drawn)
            ##should we shuffle these values again???
            
            
            #print('check the new low means - total deformation')        
            #print('threshold:', dummy_td)
            #print('high:', td_high_mean)
            #print('est mean:', td_mean)
            #print('sample mean:', np.mean(sample))
            #print('step low:', np.mean(np.ma.array(td,mask=td>dummy_td)))
            #print('est low:', np.mean(sample_low))  #these values should still be lower than threshold
            #print('old mean:', np.mean(td))
            
            ##td[threshold] = td_low_mean
            #td[threshold] = sample_low
            
            ###sample_low has still many high values
            ##print(sample_low[0])
            ##print(sample_low[-1])
            ##print(sample[-1])
            ###exit()

            
            
            ##we should also spatially order the reamining high values in sample_low
            ##triangles that are inside the buffer of the high deformation should get priority
            ##or use triangles that are inside the buffer of the angles lines
            
            ##get distances of all triangles (centroids) from the buffer/lines
            #dist1,dist2 = get_distance_order(tri,tripts,pindex,lkfs)
            
            
            ##sort sample from highest to lowest
            #sample = np.array(sorted(sample,reverse=True))
            
            ##td[threshold] and dist1,dist2 have same order
            ##sample_low needs to be sorted in a same way
            ##indexes that would sort the array from min to max: indexes to arrive to sorted state
            #dist=np.array(dist1)
            #idx = np.argsort(dist)
            ##print(idx)
                        
            ##we need indexes to arrive from sorted back to unsorted
            #idxu = range(0,len(dist))
            #idxuu = [s for _,s in sorted(zip(dist,idxu))]
                        
            #sample_a=np.zeros_like(dist)
            #sample_a[idxuu] = sample
            #sample_low = sample_a[threshold]
            
            #td[threshold] = sample_low
            
            ####angle lines
            ##dist=np.array(dist2)
            ##idx = np.argsort(dist)
            ##print(idx)
                        
            ###we need indexes to arrive from sorted back to unsorted
            ##idxu = range(0,len(dist))
            ##idxuu = [s for _,s in sorted(zip(dist,idxu))]
                        
            ##sample_a=np.zeros_like(dist)
            ##sample_a[idxuu] = sample
            ##sample_low = sample_a[threshold]
            
            ##td[threshold] = sample_low

            
            
            
            
            
            
            #print('new mean:', np.mean(td)) #this should be same as estimated mean!
            ##exit()
            
            
            ###this is the problem!!!
            ##ratio2 = np.mean(np.abs(div_f2))/np.mean(shr)
            ###print(ratio2)
            ###exit()
            ##div_low_mean = td_low_mean /np.sqrt(2) * ratio2
            ##if np.mean(div_f2) < 0:
                ##div_low_mean = div_low_mean *-1
            ##shr_low_mean = td_low_mean * (1-ratio2)
            
            ##print('check the new low means - divergence')
            ##print('old mean div:', np.mean(np.abs(div_f2)))
            ##print('threshold:', dummy_div)
            ##print('est low:', div_low_mean)
            
            
            
            ##what if we simply replace the low values by the low ship radar data at this initial scale
            ##we already know that the high values are co-inciding...
            ##and generate a random sample from them
            
            ##why do we actually need there maps for?
            ##just to recreate a perfect power law?
            
            ##to keep them for later to generate topography???
            
            ##is our distribution really log-normal?
            ##or is it logarithmic???
            
            ##or what if we have to make log-normal distribution for real mean and then use the dummy to get only the low values???
            
            ##div_f2[threshold] = td_low_mean/np.sqrt(2)
            
            
            ###can we have a randoly distributed values instead?
            ###log-normal distribution with known mean and standard deviation?
            ##sample_mean = td_low_mean/np.sqrt(2)
            ##sample_mean = np.mean(np.ma.array(div_f2,mask=~threshold))
            ##sample_mean = np.mean(np.ma.array(np.abs(div_f2),mask=td>dummy_td))     #this mean should be much lower!!!
            ##a_std = 1.6299983312588754e-07
            ##k_std = -0.9780933554534802
            ##sample_std = a_std*ls_low_mean**k_std
            ##sample = generate_lognormal_samples(sample_mean, sample_std, n=n2)
            
            #sample_div = sample_low/np.sqrt(2)
            
            ##then also randomly multiply by -1???
            ##for example: randrange(0,1), then change all 0 to -1
            #sign = np.random.random_integers(0, 1, size=n2)
            #sign = np.where(sign==0,-1,sign)
            #print(sign)
            #sample_div = sample_div*sign
            ##exit()
            
            #div_f2[threshold] = sample_div
            #print('new mean div:', np.mean(np.abs(div_f2)))
                
            #print('check the new low means - shear')        
            #print('old mean shr:', np.mean(shr_f2))
            #print('threshold:', dummy_shr)
            ##print('est low:', shr_low_mean)
            
            ##shr_f2[threshold] = td_low_mean/np.sqrt(2)
            #sample_shr = sample_low/np.sqrt(2)
            #shr_f2[threshold] = sample_shr
            #print('new mean shr:', np.mean(shr_f2))  
            
            #test = np.sqrt(np.mean(shr_f2)**2+np.mean(div_f2)**2)       #this ones comes quite low...
            #print(test)
            
            ##exit()

    
    ##course-grain these results:
    #ddidx = 0
    #for j in stp:
        #print('Step: '+str(j))
        #outname_td = 'td_seed_f_'+reg+'_L'+str(j)+file_name_end
        #outname_ts = 'ts_seed_f_'+reg+'_L'+str(j)+file_name_end
        #outname_angles = 'angle_seed_f_'+reg+'_L'+str(j)+file_name_end
        #td_list=[]
        #ls_list=[]
        #ang_list=[]
        #id_list=[]
        #time_list=[]
        #date_list=[]
        #date_ts=[]
        #mdiv=[]
        #mpdiv=[]
        #mndiv=[]
        #mshr=[]
        #mdiv_sd=[]
        #mpdiv_sd=[]
        #mndiv_sd=[]
        #mshr_sd=[]
        #angles_list=[]
        
        
        #if j > 1:
            ##time series dont need corse-graining
            #if time_series==True:
                #print('time series');continue

            ##get seeding points for each lenght scale step
            ##flatten the array
            #xs = x.filled()[::j,::j]
            #ys = y.filled()[::j,::j]
            #hpms = hpm.filled()[::j, ::j]
            
            ##mask out all poor quality data
            ##5 is sufficient as these are just the seeding polygons
            #gpi = hpms > 9    
            #xs = xs[gpi]
            #ys = ys[gpi]
            
            ###add corner points for higher steps
            ##if j > 5:
                ##xlist = [x[idcr1],x[idcr2],x[idcr3],x[idcr4]]
                ##ylist = [y[idcr1],y[idcr2],y[idcr3],y[idcr4]]
                ##xs = np.append(xs,xlist)
                ##ys = np.append(ys,ylist)

            ##keep only valid data
            #xs = np.ma.masked_invalid(xs)
            #xs = np.ma.compressed(xs)
            #ys = np.ma.masked_invalid(ys)
            #ys = np.ma.compressed(ys)

            ##check how many nods do we have left
            #if len(xs) < 3:
                #print('No nods for triangles left!')
                #continue

            ##triangulate between these seeding points
            #pts_seed = np.zeros((len(xs),2))
            #pts_seed[:,0]=xs; pts_seed[:,1]=ys
            #tri_seed = Delaunay(pts_seed)
            #tripts_seed = pts_seed[tri_seed.simplices]
            
            ################################################################3
            ##continue mesh plot
            #if j > 11:
                #alpha=1
            #else:
                #alpha=0.5
            #clm = next(colormesh)
            #gx.triplot(pts_seed[:,0], pts_seed[:,1], tri_seed.simplices.copy(), color=clm, alpha=alpha, label=str(j))
            ################################################################3
            
            ##WARNING: these are triangles of various sizes at LKFs, in MIZ and at the scene edges
            ##allow max 2x larger triangles than optimal for smallest scale for seeding
            #max_area = (distance*2)**2/2
            ##do area-weighted coarse-graining
            #div_seed,shr_seed,area_seed,minang_seed = coarse_grain(tripts,tripts_seed,div_f2,shr_f2,max_area,minang_f2)
            
            #td_seed = np.sqrt(div_seed**2 + shr_seed**2)
            
            #if low_mean==False:
                ##mask with threshold that depends on the lenght scale (mask out large triangles)
                ##this is done in the plotting, so not necessary here...
                #dyn_dummy = dummy_td_all[ddidx]
                #ls_dummy = dummy_ls_all[ddidx]
                #if no_threshold==True:
                    #dyn_dummy = 0
                ##threshold_seed = (td_seed<dyn_dummy) #| (area_seed > ls_dummy**2)
                #threshold_seed = (td_seed<dyn_dummy)
                #ddidx = ddidx + 1
                
                ##masking
                #div_seed = np.ma.array(div_seed,mask=threshold_seed)
                #shr_seed = np.ma.array(shr_seed,mask=threshold_seed)            
                #area_seed = np.ma.array(area_seed, mask=threshold_seed)
                #minang_seed = np.ma.array(minang_seed, mask=threshold_seed)
                ##id_seed = np.ma.array(id_seed,mask=threshold_seed)
            
        #else:
            ##if LKF_filter==False:
                ##lkf_id = np.zeros_like(area)
            
            #if low_mean==False:
                #div_seed=np.ma.array(div_f2, mask=threshold)
                #shr_seed=np.ma.array(shr_f2, mask=threshold)
                #tripts_seed=tripts
                #area_seed = np.ma.array(area, mask=threshold)
                #minang_seed = np.ma.array(minang, mask=threshold)
                ##id_seed = np.ma.array(lkf_id,mask=threshold)
                                
            #else:
                ##keep all the data
                #div_seed=div_f2
                #shr_seed=shr_f2
                #tripts_seed=tripts
                #area_seed = area
                #minang_seed = minang
                ##id_seed = lkf_id
            
        
        #####################################################################3
        #if image == True:
            ###Plotting deformation 
            #deform = div_seed*1e6
            #outname = outpath+'map_div_seed_f_'+reg+'_L'+str(j)+'_'+date1+'_'+date2
            #label = r'Divergence (10$^6$s$^{-1}$)'
            #interval = [-5, 5]
            #cmap=plt.cm.bwr
            
            ##deform = shr_seed*1e6
            ##outname = outpath+'map_shr_'+reg+'_L'+str(j)+'_'+date1
            ##label = r'Shear (10$^6$s$^{-1}$)'
            ##interval = [0, 10]
            ##cmap=plt.cm.Reds
            
            ##if j > 1:
                ###otherwise masked values are plotted
                ##deform = deform.filled(fill_value=0)

            
            ##print(outname)
            #plot_def(area_def,tripts_seed,deform,outname,label,interval,cmap,ship_lon,ship_lat,radius)
        ######################################################################3
            ##exit()

        ###write out lists into csv file
        ##tt = [date_list, ls_list, time_list, td_list, ang_list]
        ##table = zip(*tt)
        ###adjusted to python3:
        ##table = list(zip(*tt))

        ##output = outpath_def + outname_td
        ##with open(output, 'ab') as f:
            ###header
            ###f.write(b'date, length scale, time difference, total deformation, min angle\n')
            ##np.savetxt(f, table, fmt="%s", delimiter=",")

        
        
        ##at high resolution deformation data we cant avoid of having two separate clouds (detectable and un-detectable deformation) - that happens when our method cant detect deformation anymore
        
        ##this can not be fixed by filtering and the line in power law will break even if we use the Sylvain filter
        
        ##still the filter removes the artifical opening-clossing mess in shear zones (and removes some high values in that zones), so it makes it useful for the mapping of larger features such as LKFs and leads.
        



        ##storing data for the scatter plots
        
        #td = np.sqrt(div_seed**2 + shr_seed**2)
        #ls = np.sqrt(area_seed)
        #minang = np.array(minang_seed)
        ##id_num = np.array(id_seed)
        
        
        #td_list.extend(np.ma.compressed(td).tolist())
        #ls_list.extend(np.ma.compressed(ls).tolist())
        #ang_list.extend(np.ma.compressed(minang).tolist())
        ##id_list.extend(np.ma.compressed(id_num).tolist())
        
        ##time handing
        #dt_tri = np.full_like(td,np.datetime64(dt1))
        #diff_tri = np.ones_like(td)*diff
             
        #date_list.extend(np.ma.compressed(dt_tri).tolist())
        #time_list.extend(np.ma.compressed(diff_tri).tolist())
                
        ##save the data for time series
        ##these are spatial averages and some value need to be used for all those masked areas - average???
        #date_ts.append(dt1)
        #posd = np.ma.array(div_seed,mask=div_seed<0)
        #negd = np.ma.array(div_seed,mask=div_seed>0)
        #mpdiv.append(np.ma.filled(np.mean(posd), fill_value=-999))                     #add fill values (-999)
        #mpdiv_sd.append(np.ma.filled(np.std(posd), fill_value=-999))
        #mndiv.append(np.ma.filled(np.mean(negd), fill_value=-999))
        #mndiv_sd.append(np.ma.filled(np.std(negd), fill_value=-999))
        #mdiv.append(np.ma.filled(np.mean(div_seed), fill_value=-999))
        #mdiv_sd.append(np.ma.filled(np.std(div_seed), fill_value=-999))
        #mshr.append(np.ma.filled(np.mean(shr_seed), fill_value=-999))
        #mshr_sd.append(np.ma.filled(np.std(shr_seed), fill_value=-999))

        ##continue
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

        ##write data for time series 
        #tt = [date_ts, mpdiv, mndiv, mdiv, mshr, mpdiv_sd, mndiv_sd, mdiv_sd, mshr_sd]
        #table = zip(*tt)
        #table = list(zip(*tt))

        #output = outpath_def + outname_ts
        #with open(output, 'ab') as f:
            ##header
            ##f.write(b'date, pos. divergence, neg. divergence, mean divergence, mean shear\n')
            #np.savetxt(f, table, fmt="%s", delimiter=",")
        
        ##angle data
        #if LKF_filter:
            #angles_list.extend(angles)
            ##print(angles_list)
            
            ##write data for angle PDFs
            #output = outpath_def + outname_angles
            #with open(output, 'ab') as f:
                #np.savetxt(f, angles_list, fmt="%s", delimiter=",")
        
    #outname='overview_mesh_seed'+date1
    #gx.legend()
    #fig4.savefig(outpath+outname)
    #plt.close(fig4)
    
            
#fig1.savefig(outpath+'overview_map_'+reg+'_'+str(int(radius/1000.)),bbox_inches='tight')
#plt.close(fig1)
