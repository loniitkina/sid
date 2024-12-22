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
from matplotlib.offsetbox import AnchoredText

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

#WARNING: tiles are just overlaid and drift between the timestamps is not taken into account. There will be discontinuities in positions between the regions. 
#TODO: implement uniform correction for the mean dirift difference between central region and the rest

#using the treshold to discard the data with poor signal/noise ratio
threshold_filter=True

#filtering for better values along the LKFs as recommended for SAR data by Buillion et al
lkf_filter=True
kernel=3    #3 is recommend by Sylvain (but for coarser resolution)
interval = [-1, 1]  #expected divergence values (div in s-1 * 10e6) for filter plot

#minang_limit is used in threshold_filter and lkf_filter
#best set very low or zero when wanting to detect all features (will include MIZ!)
#can be set to 15, when working with scalling
minang_limit=5  #this will leave out some triangles over wide LKFs (dense points at the edges of the feature make narrow triangles!)
minang_limit=2

#min amount of data to do the triangulation
min_nod_number=100

#select area size
#radius = 120000 #will leave gaps between the 200-km tiles
radius = 200000 #some hard-coded naming depends on this

#should we triangulate second or first SAR scene coordinates
#second works works OK for deformation and CDE products, but not for parcels yet
second_latlon=False

#naming and output
rname = '_'+str(int(radius/1000))+'km'
if threshold_filter:
    tname='_thfilter'
    #such output is only sensible if threshold filter is on
    filter_fig=True
    #active floe size
    afs=True
    ##for parcel tracking we need to have consequtive data: second scene in pair needs to be first scene in the next pair! (combo option is not possible here)
    #just save level 1 data and exit
    parcel=True
    #LKF angles
    lkf_angles=False
else:
    tname='_nothfilter'
    afs=False
    parcel=False
    lkf_angles=False
    
if lkf_filter:
    lname='_lkffilter'
else:
    lname='_nolkfilter'

file_name_end = rname+tname+lname
print('Your output will be: ',file_name_end)

file_name_end_csv = file_name_end+'.csv'

scaling=False
#create log-spaced vector and convert it to integers
n=9 # number of samples
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)

#-------------------------------------------------------------------INPUT
#threshold detection values: step and factor differ for each for each input data (usually put in different folders/paths)
#resolution=factor X step
#at factor=1 and step=10, we get 400m distance between grid points
#distance between nods==400m or triangle area is 400*400/2 
#there is one nod, where displacement difference is 40m more than at the other two
#example velocities nod1=1000m/diff, nod2=1000m/diff, nod3=1040m/diff
#exaggerate the factor as the possible displacements are discrete steps as 40, 80, 120 m...
#diff is different for every scene pair and threshold needs to be recalculated

#extra_margin = 20   #some allowance #20 in test5 gave best results so far (but also lets in the last artifacts in LKF filter)
#20 was a good margin for factor=80, seems like 25% was a good estimate
#no margin necessary for scaling, but use 10 to get clean LKFs for other analysis where 'less is more' (e.g. parcels)
extra_margin = 30   #this is a good option for 'less is more': removed most of the artifacts (and some data) should work well for parcels and AFS!

#800 m resolution data (pixel/factor=80m, step=10)
#inpath = '/scratch/pit000/results/sid/drift/stp10_factor05/'
inpath = '/scratch/pit000/results/sid/drift/stp10_factor05'+rname+'/'
step = 10      #grid resolution of drift 
factor = 80    #pixel size (factor in sid_drift, default factor in FT is 0.5, which translates to 80 m here)

##200 m resolution data
#inpath = '/scratch/pit000/results/sid/drift/stp5_factor1/'
#step = 5
#factor = 40

##MOSAiC
#reg = 'mosaic'   #used in filenames
##for mapping/projection
#shipfile = '../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0.csv'
##shipfile = '../../downloads/data_master-solution_mosaic-leg2-20191214-20200224-floenavi-refstat-v1p0.csv'
##shipfile = '../../downloads/position_leg3_nh-track.csv'	#leg3 (and transition to leg 4 until 6 June)
##Region limits for the overview map
#lon_diff = 15
#ship_lon=17.147909; ship_lat=87.132429      #March/April event start
#regn = ship_lat+.1; regs = ship_lat-4
#regw = ship_lon-lon_diff; rege = ship_lon+lon_diff

#N-ICE
reg = 'lance_leg1'   #used in filenames
shipfile = '../../downloads/lance_leg1.csv'
#Region limits for the overview map
regn = 85; regs = 79
regw = 0; rege = 60

##CIRFA cruise
#reg = 'cirfa'   #used in filenames
#shipfile = '../../downloads/CIRFA_cruise_stationM.csv'
##Region limits for the overview map
#regn = 82; regs = 77
#regw = -25; rege = 5

#for getting timestampts of drift files
main_trackfile_tail = '_c'+rname+'-fnames.csv'
main_trackfile=shipfile.split('.csv')[0]+main_trackfile_tail
print('Your study area is: ',main_trackfile)

#------------------------------------------------------------------OUTPUT
outpath_def = '/scratch/pit000/results/sid/deform200km/'
#outpath_def = '/scratch/pit000/results/sid/deform200km_stp5_res40/'

#outpath = '/scratch/pit000/results/sid/plots/'
outpath = '/scratch/pit000/results/sid/plots200km/'
#outpath = '/scratch/pit000/results/sid/plots200km_stp5_res40/'

#-------------------------------------------------------------------------------------------------------------------end of parameters
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

xc1,yc1 = transform(inProj,outProj,regw, regs)
xc2,yc2 = transform(inProj,outProj,rege, regn)
area_extent = (xc1,yc1,xc2,yc2)
area_def2 = AreaDefinition(area_id, description, proj_id, proj_dict, 300, 1000, area_extent)

fig1    = plt.figure(figsize=(10,10))
ax      = fig1.add_subplot(111)
ma = pr.plot.area_def2basemap(area_def2, resolution='h')
ma.drawmapscale(regw+5, regs+.6, 10, 82, 100, units='km', barstyle='fancy',fontsize=14)    #drawmapscale(lon, lat, lon0, lat0, length, **kwargs)
ma.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
ma.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
ma.drawcoastlines()
ma.fillcontinents()
#ma.etopo()
#ma.drawmapboundary()#color='k', linewidth=1.0, fill_color=None, zorder=None, ax=None)

#remove all the output text files from previous script runs (or they will be appended)
rlist = glob(outpath_def+'td_*'+reg+'*'+file_name_end_csv)
for fn in rlist:
    os.remove(fn)
rlist = glob(outpath_def+'ts_*'+reg+'*'+file_name_end_csv)
for fn in rlist:
    os.remove(fn)

outname_dummy = 'dummy_'+reg+rname+'.csv'
rlist = glob(outpath_def+outname_dummy)
for fn in rlist:
    os.remove(fn)

#regions have to be arranged CW  - or convex hull of frames wont work!
#start with c and continue with most likely good coverage (e.g. s)
regions=['c','msw','sw','s','mse','se','e','mne','ne','n','mnw','nw','w']
regions=['c']

noons = getColumn(main_trackfile,0,header=False)
noons = [ datetime.strptime(noons[i], "%Y-%m-%d %H:%M:%S") for i in range(len(noons)) ]
days = [ datetime.strftime(noons[i], "%Y%m%d") for i in range(len(noons)) ]
print(days)

#colors for overview map
color=iter(plt.cm.rainbow_r(np.linspace(0,1,len(days)+1)))

#days=['20200401','20200402','20200403','20200511']
#days=['20191114']
#days = ['20150115']
#days=['20150119','20150120','20150121','20150122','20150123','20150124','20150125','20150126','20150127','20150203']

#for plotting
ship_lon_daily=[]
ship_lat_daily=[]

for day in days:
    print(day,'******************************************************************************************************************************')

    #color for the overview map
    cl = next(color)
    
    #check if any data was collected - important for saving at the end of this loop
    anydata=False
    
    #for backward overlaying of sigma data
    sigma_data=[]
    xs1_data=[]
    ys1_data=[]

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
        
        #pick the right reference date
        if second_latlon:
            ref_date=dt2
        else:
            ref_date=dt1
        
        #do we need to tile the SAR intensities here?
        
        
        if region=='c':
            #tile center - from track file
            shipfile_tile = shipfile.split('.csv')[0]+'_'+region+'_200km.csv'
            print(shipfile_tile)
            mettime = getColumn(shipfile_tile,0)
            dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
            if dtb[0]>ref_date: continue
            if dtb[-1]<ref_date: continue
            mi = np.argmin(abs(np.asarray(dtb)-ref_date))
            ship_lon = np.asarray(getColumn(shipfile_tile,1),dtype=float)[mi]
            ship_lat = np.asarray(getColumn(shipfile_tile,2),dtype=float)[mi]
            if np.isnan(ship_lon): continue
            print('center at: ',ship_lon,ship_lat,dtb[mi])

            #ship centered projection
            #work in geographical projection, where units are meters!
            
        
            
            xlp,ylp = transform(inProj,outProj,ship_lon,ship_lat)
            ##Using a projection dictionary (same stuff as for overview map)
            #width = radius*2/factor # m spacing
            #height = radius*2/factor # m spacing
            #area_extent = (xlp-radius,ylp-radius,xlp+radius,ylp+radius)
            #area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
            #m = pr.plot.area_def2basemap(area_def)
            
            #also keep the dates from the central tile for filenaming
            date1_c=date1
            date2_c=date2
            
            ##reproject center position to image coordinates
            #xl, yl = m(ship_lon, ship_lat)

        #WARNING: here we start working with the second image coordinates
        if second_latlon:
            #reproject vertices to geo coordinates
            x1, y1 = transform(inProj,outProj,lon,lat)
            
            #get lat,lon for the second image in the SAR pair
            #get displacements
            dx = diff*u; dy = diff*v
            #estimate positions after drift
            x2 = x1 + dx
            y2 = y1 + dy
            lon,lat = transform(outProj,inProj,x2,y2) 
        
        #reproject vertices to image coordinates
        x, y = transform(inProj,outProj,lon,lat)
        
        #cut out region
        mask = ~((x<xlp-radius) | (x>xlp+radius) | (y<ylp-radius) | (y>ylp+radius))
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
        
        ##find corners of valid data
        #idcr1 = np.unravel_index(np.argmin(np.sqrt((x-(xl-radius))**2+(y-(yl-radius))**2)),x.shape)
        #idcr2 = np.unravel_index(np.argmin(np.sqrt((x-(xl-radius))**2+(y-(yl+radius))**2)),x.shape)
        #idcr3 = np.unravel_index(np.argmin(np.sqrt((x-(xl+radius))**2+(y-(yl+radius))**2)),x.shape)
        #idcr4 = np.unravel_index(np.argmin(np.sqrt((x-(xl+radius))**2+(y-(yl-radius))**2)),x.shape)
        
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

        #make timestemp and time difference data per point
        dt1s = np.full_like(xs,np.datetime64(dt1))
        dt2s = np.full_like(xs,np.datetime64(dt2))
        diffs = np.ones_like(xs)*diff

        #check how many values are left in the region
        print('Values left for this step:')
        print(us.size)
        if region=='c' and us.size < min_nod_number: #if no data in the central tile - stop here
            print('empty central tile');break
        elif us.size < min_nod_number:
            continue            

        print('Start tiling region: ',region)
        #a region of all area where there is data (following image pair edges)
        #convert to format that shapely can work with
        ctr_nods = [ np.array([i,j]) for i,j in zip(xs,ys) ]

        #get initial tile for this region
        tile = MultiPoint(ctr_nods).convex_hull
        
        #get rid of large triangles at the tile edge - artefacts - 20km
        #centroid has to lay inside this buffer polygon
        good_nods=[]
        eroded=tile.buffer(-20000)    #this erodes the polygon by buffer width
        for i in range(0,len(ctr_nods)):
            if eroded.disjoint(Point(ctr_nods[i])):# and ( area[i] > (distance*2)**2/2 or minang[i] < 10 ):
                good_nods.append(False)
            else:
                good_nods.append(True)
        
        xs = xs[good_nods]
        ys = ys[good_nods]
        lats = lats[good_nods]
        lons = lons[good_nods]
        us = us[good_nods]
        vs = vs[good_nods]
        dt1s = dt1s[good_nods]
        dt2s = dt2s[good_nods]
        diffs = diffs[good_nods]
        
        print('Get refined tile for this region')
        ctr_nods = [ np.array([i,j]) for i,j in zip(xs,ys) ]
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
        
        #print(tile.geom_type)
            
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
            
            tiled_xs=xs.copy()
            tiled_ys=ys.copy()
            tiled_lats=lats.copy()
            tiled_lons=lons.copy()
            tiled_us=us.copy()
            tiled_vs=vs.copy()
            tiled_dt1s=dt1s.copy()
            tiled_dt2s=dt2s.copy()
            tiled_diffs=diffs.copy()
            
            ##plot the area covered
            #xg, yg = tiled_cover.exterior.xy 
            #ax.fill(xg, yg, fc='none', ec=cl)
            
            #get plot frame polygon
            xlist = [xlp-radius,xlp-radius,xlp+radius,xlp+radius,xlp-radius]
            ylist = [ylp-radius,ylp+radius,ylp+radius,ylp-radius,ylp-radius]
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

            tiled_xs=np.concatenate((tiled_xs, xs[new_nods]), axis=0)
            tiled_ys=np.concatenate((tiled_ys, ys[new_nods]), axis=0)
            tiled_lons=np.concatenate((tiled_lons, lons[new_nods]), axis=0)
            tiled_lats=np.concatenate((tiled_lats, lats[new_nods]), axis=0)
            tiled_us=np.concatenate((tiled_us, us[new_nods]), axis=0)
            tiled_vs=np.concatenate((tiled_vs, vs[new_nods]), axis=0)
            tiled_dt1s=np.concatenate((tiled_dt1s, dt1s[new_nods]), axis=0)
            tiled_dt2s=np.concatenate((tiled_dt2s, dt2s[new_nods]), axis=0)
            tiled_diffs=np.concatenate((tiled_diffs, diffs[new_nods]), axis=0)
            
            
            #plot the new area covered and extend coverage - dont do convex hull yet or internal area are filled in!
            if diff_tile.geom_type == 'MultiPolygon':
                #for geom in diff_tile.geoms: 
                    
                    ##plot
                    #xg, yg = geom.exterior.xy 
                    #ax.fill(xg, yg, fc='none', ec=cl)
                    
                #extend
                tiled_cover=tiled_cover.union(diff_tile)
            
            #sometimes the difference may be also point and lines and not just polygons
            elif diff_tile.geom_type == 'GeometryCollection':
                for geom in diff_tile.geoms:
                    if geom.geom_type == 'Polygon':
                        #plot
                        xg, yg = geom.exterior.xy 
                        gx.fill(xg, yg, fc='none', ec=cl)
                        
                        #extend
                        tiled_cover=tiled_cover.union(geom)
            else:
                ##plot
                #xg, yg = diff_tile.exterior.xy 
                #ax.fill(xg, yg, fc='none', ec=cl)
                
                #extend
                tiled_cover=tiled_cover.union(diff_tile)
                    
    #store combined data and figures for all regions/tiles
    if anydata:
        print('We have data for date: ',day)
        
        #here it is safe to do convex hull - no more new regions are coming
        tiled_cover = tiled_cover.convex_hull
        print(tiled_cover)
        #store as vertices - This should not be MultiPolygon
        xg, yg = tiled_cover.exterior.xy
        tiled_cover = [ np.array([i,j]) for i,j in zip(xg,yg) ]
        
        #plot on the overview map
        ##points
        #xa, ya = ma(tiled_lons, tiled_lats)
        #ax.plot(xa,ya,'.',ms=.1,color=cl, alpha=.1)
        #ship
        xla, yla = ma(ship_lon, ship_lat)
        ax.plot(xla,yla,'*',markeredgewidth=1,color=cl,markersize=15,markeredgecolor='k')
        #tile
        xg_lon,yg_lat = transform(outProj,inProj,xg,yg)
        xg,yg=ma(xg_lon,yg_lat)
        frames = ax.fill(xg, yg, fc='none', ec=cl)
        
        #store the data
        out_file = inpath+'Tile_'+str(step)+'_'+str(factor)+'_'+date1_c+'_'+date2_c+file_name_end+'.npz'
        np.savez(out_file, xs=tiled_xs, ys=tiled_ys, lons=tiled_lons, lats=tiled_lats, us=tiled_us,
                    vs=tiled_vs, dt1s=tiled_dt1s, dt2s=tiled_dt2s, diffs=tiled_diffs, region=tiled_cover, corners=corner_nods)
        
        ship_lon_daily.append(ship_lon)
        ship_lat_daily.append(ship_lat)

#put the entire ship track over
xla, yla = ma(ship_lon_daily, ship_lat_daily)
ax.plot(xla,yla,c='k',lw=2)

cb = plt.colorbar(frames, ax=ax, pad=.01)  # draw colorbar
#cb.set_label(label='Intensity (dB)',fontsize=20)
cb.ax.set_yticklabels(days.strftime('%d %b %Y'))
cb.ax.tick_params(labelsize=20)

outname = outpath+'tiling_map_'+reg+file_name_end
print(outname)
plt.show()
fig1.savefig(outname,bbox_inches='tight')
plt.close(fig1)

