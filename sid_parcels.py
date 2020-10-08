from datetime import datetime
from glob import glob
import numpy as np
import pyresample as pr
from sid_func import *
import matplotlib.pyplot as plt

#for parcel tracking we need to have consequtive data: second scene in pair needs to be first scene in the next pair! (combo option is not possible here)



radius = 60000          #use a radius shorter than in the parcel producing script - or lots of boundary parcels will be lost fast!
spacing= 600            #at spacing shorter than 500m, there will be lots of nans...


#for output
fp = 5  #final image pair

#-------------------------------------------------------------------
inpath = '../sidrift/data/stp10_parcels_f1/'
inpath = '../sidrift/data/stp10_parcels_f1_rs2/'    #alternative endings
inpath_drift = inpath
outpath_data = inpath
outpath = inpath

#inpath = '../sidrift/data/80m_stp10_single_filter/'
#inpath_drift = '../sidrift/data/40m_combo/'
#outpath_data = inpath
#outpath = inpath

metfile = '../sidrift/data/10minute_nounits.csv'

#get Lance location and seed the parcels
start = datetime(2015,1,21,6,54,4)                 #this coincides with the 1st SAR image date/time
#Lance postion (from Lance's met system)
mettime = getColumn(metfile,0)
dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
mi = np.argmin(abs(np.asarray(dtb)-start))
llon_dtb = np.asarray(getColumn(metfile,2),dtype=float)
llat_dtb = np.asarray(getColumn(metfile,1),dtype=float)
Lance_lon = llon_dtb[mi]
Lance_lat = llat_dtb[mi]

#Lance start location in geographical coordinates
from pyproj import Proj, transform
inProj = Proj(init='epsg:4326')
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 0))
xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)

#create a grid of spacing of 100m in each direction in the radius from the ship's initial position
#Lance will be in the center
xbf = np.arange(xlp-radius, xlp+radius, spacing)    #initial conditions are same as the final grid!
ybf = np.arange(ylp-radius, ylp+radius, spacing)
x_buoy,y_buoy = np.meshgrid(xbf,ybf)
#flatten
x_buoy = x_buoy.flatten()
y_buoy = y_buoy.flatten()

#read in the drift product and calculate displacements from starting point until the first break in the SAR data (1 week)
fl = sorted(glob(inpath_drift+'SeaIceDrift*.npz'))

fl_dmg = sorted(glob(inpath+'Damage*.npz'))

#empty arrays to store parcel drifts and damage
x_path = np.zeros((len(fl_dmg)+1,x_buoy.shape[0]))
y_path = np.zeros((len(fl_dmg)+1,x_buoy.shape[0]))
damage = np.zeros((len(fl_dmg)+1,x_buoy.shape[0]))
leads = np.zeros((len(fl_dmg)+1,x_buoy.shape[0]))
ridges = np.zeros((len(fl_dmg)+1,x_buoy.shape[0]))
 
#start with initial coordinates
x_path[0,:] = x_buoy
y_path[0,:] = y_buoy
date = [start]

for i in range(0,len(fl_dmg)):
    #break
    #read in the drift data
    print(fl[i])
    container = np.load(fl[i])
    lat1 = container['lat1']     #these arrays are as big as both scenes together (nan-where no overlap), ~5000x5000 of 80m-pixels = 400x400km
    lon1 = container['lon1']
    u = container['upm']
    v = container['vpm']
    hpm = container['hpm'] 
    
    #quality filter 
    gpi = hpm > 9
    u = u[gpi]
    v = v[gpi]
    lon1 = lon1[gpi]
    lat1 = lat1[gpi]
    
    #get damage data
    print(fl_dmg[i])
    container = np.load(fl_dmg[i])
    lat_dmg = container['lat']      #triangle centroids
    lon_dmg = container['lon']
    dmg = container['d']
    rid = container['r']
    lea = container['l']
    
    #project data coordinates
    x1,y1 = transform(inProj,outProj,lon1,lat1)
    xd,yd = transform(inProj,outProj,lon_dmg,lat_dmg)
    
    #get time difference for drift
    date1 = fl[i].split('_')[-2]
    date2 = fl[i].split('_')[-1].split('.')[0]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    diff = (dt2-dt1).total_seconds()
    #print(diff)
    
    #get displacements
    dx = diff*u; dy = diff*v
    
    #print(u.shape)
    #print(x1.shape)
    #print(dx.shape)
    #print(np.min(np.ma.masked_invalid(dx)))
    #print(np.max(np.ma.masked_invalid(dx)))
    #print(np.mean(np.ma.masked_invalid(dx)))
    
    #check how much Lance drifted in this time and apply one number correction
    #this helps keeping the deformed ice in a line (relative location), but it deteriorates the absolute positions
    mi = np.argmin(abs(np.asarray(dtb)-dt1))
    Lance_lon1 = llon_dtb[mi]
    Lance_lat1 = llat_dtb[mi]

    mi = np.argmin(abs(np.asarray(dtb)-dt2))
    Lance_lon2 = llon_dtb[mi]
    Lance_lat2 = llat_dtb[mi]
    
    lx1,ly1 = transform(inProj,outProj,Lance_lon1, Lance_lat1)
    lx2,ly2 = transform(inProj,outProj,Lance_lon2, Lance_lat2)
    
    ldx = lx2-lx1
    ldy = ly2-ly1
        
    #find local displacement (closest coordinate) and correct uniformly for that displacement
    mask_l = (x1>lx1-spacing) & (x1<lx1+spacing) & (y1>ly1-spacing) & (y1<ly1+spacing)
    #print(x1[mask_l], y1[mask_l])
    dx_near = np.mean(np.ma.masked_invalid(dx[mask_l]))
    dy_near = np.mean(np.ma.masked_invalid(dy[mask_l]))
    #print(dx_near, dy_near)
    
    x_correction = dx_near-ldx
    y_correction = dy_near-ldy
    print(x_correction, y_correction)
    
    dx = dx-x_correction
    dy = dy-y_correction
    
    #estimate positions after drift
    x2 = x1 + dx
    y2 = y1 + dy
    
    #save the date 
    date.append(dt2)
    
    #find closest drift/damage coordinate for each parcel
    print(x_buoy.shape[0])
    for m in range(0,x_buoy.shape[0]):
        #print(m)
        #sea ice drift/displacement mask for data with 400m resolution, 600m spacing between buoys, 300 m drift precission (based on buoys)    
        mask = (x1>x_buoy[m]-spacing) & (x1<x_buoy[m]+spacing) & (y1>y_buoy[m]-spacing) & (y1<y_buoy[m]+spacing)
        #damage mask for data with ~400m resolution            
        mask_d = (xd>x_buoy[m]-spacing/2) & (xd<x_buoy[m]+spacing/2) & (yd>y_buoy[m]-spacing/2) & (yd<y_buoy[m]+spacing/2)
        
        #get the mean location for parcel (and update parcel location for the search in next step)
        #print(np.ma.compressed(np.ma.masked_invalid(x2[mask])))
        x_buoy[m] = np.mean(np.ma.masked_invalid(x2[mask]))
        y_buoy[m] = np.mean(np.ma.masked_invalid(y2[mask]))
        x_path[i+1,m] = x_buoy[m]   #step 0 is initial condition, step 1 is reached after one day etc.
        y_path[i+1,m] = y_buoy[m]
        
        #check if there was any deformation
        #print(np.ma.compressed(np.ma.masked_invalid(dmg[mask_d])))
        dd = np.mean(np.ma.masked_invalid(dmg[mask_d])) #this is just 0 and 1 data, if empty, we get a nan and loose the buoy...
        damage[i+1,m] = dd
        ll = np.mean(np.ma.masked_invalid(lea[mask_d]))
        leads[i+1,m] = ll
        rr = np.mean(np.ma.masked_invalid(rid[mask_d]))
        ridges[i+1,m] = rr
        
        ##how close are other parcels?
        ##check in all 4 directions if we have a neighbor at less than 2 spacing away. if not create new neighbor at 1 spacig away
        #if np.min(np.abs(x_buoy-x_buoy[m])) > 2*spacing:
            ##check that this is not out of region
            #print('empty space')
            ##make new parcel and attach to the end of the list
            
            ##how long is x_buoy.shape[0]?
            #print(x_buoy.shape[0])
        
            
    #plt.plot(damage[i+1,:],label='damage')
    #plt.plot(leads[i+1,:],label='leads')
    #plt.plot(ridges[i+1,:],label='ridges')
    #plt.legend()
    #plt.show()
    #exit()
    
    #controlling parcels: adding and removing 
    #mark the bounary VBs and dont search for new parcles next to them
    #flatten the array into a list, also results need to be flat. add ID to keep the track of individual parcels as their relative position may change!
    #check for every buoy how far the others are and add new one as soon as the distance in x or y direction < 2x grid step: leads
    #remove buoys if they have distance of less than 50m: ridges
    

            

#save the output locations (in latlon) to be plotted on top of divergence and shear maps by sid_defrom.py
lon_path,lat_path = transform(outProj,inProj,x_path,y_path)
print(x_path[:,0])
print(lat_path[:,0]) #does 90 come from nan???

#dump data into numpy file
out_file = outpath_data+'VB.npz'
np.savez(out_file,x_path=x_path,y_path=y_path,lon_path=lon_path, lat_path=lat_path, damage = damage, leads=leads, ridges=ridges)

#save dates separatelly
out_file = outpath_data+'VB_dates.npz'
np.savez(out_file,date=date)

#make some plots
out_file = outpath_data+'VB.npz'
container = np.load(out_file)
print(container.files)
lon_path = container['lon_path']
lat_path = container['lat_path']
damage = container['damage']
leads = container['leads']
ridges = container['ridges']


print(lat_path[:,10])        #no more data on image 12
#exit()


out_file = outpath_data+'VB_dates.npz'
container = np.load(out_file, allow_pickle=True)
date = container['date']

#radius_proj=55000           #big enough to contain all the moving ice
radius_proj=radius
from pyresample.geometry import AreaDefinition
#Using a projection dictionary
area_id = 'around Lance'
description = 'North Pole LAEA'
proj_id = 'lance'
proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':0, 'a':6378137, 'b':6356752.3142, 'units':'m'}
width = radius_proj*2/100 #100 m spacing    (should have same resolution as the final gridded map)
height = radius_proj*2/100 #100 m spacing
area_extent = (xlp-radius_proj,ylp-radius_proj,xlp+radius_proj,ylp+radius_proj)
area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
#print(area_def)
m = pr.plot.area_def2basemap(area_def)

color=iter(plt.cm.rainbow(np.linspace(0,1,lon_path.shape[0])))

for i in range(1,fp+1):#lon_path.shape[0]):
    print(i)
    fig1    = plt.figure(figsize=(30,10))
    ax      = fig1.add_subplot(131)
    bx      = fig1.add_subplot(132)
    cx      = fig1.add_subplot(133)
    
    ax.set_title('damage')
    bx.set_title('leads')
    cx.set_title('ridges')
    
    #x,y = m(lon_path[i,:,:],lat_path[i,:,:])
    #ax.pcolormesh(x,y,damage[i])
    
    #uncomment here if you want it with geographical coordinates
    #m = pr.plot.area_def2basemap(area_def)
    #m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    #m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

    lop = np.ma.masked_where(~(np.isfinite(lat_path[i,:])),lon_path[i,:]); lop = np.ma.compressed(lop)
    lap = np.ma.masked_invalid(lat_path[i,:]); lap = np.ma.compressed(lap)
    x,y = m(lop,lap)
    
    dd = np.ma.masked_where(~(np.isfinite(lat_path[i,:])),damage[i,:]); dd = np.ma.compressed(dd)
    rr = np.ma.masked_where(~(np.isfinite(lat_path[i,:])),ridges[i,:]); rr = np.ma.compressed(rr)
    ll = np.ma.masked_where(~(np.isfinite(lat_path[i,:])),leads[i,:]); ll = np.ma.compressed(ll)
    
    ax.scatter(x,y,c=dd,lw=5)         #marker size should be radius dependant
    bx.scatter(x,y,c=ll,lw=5)
    cx.scatter(x,y,c=rr,lw=5)

    cl = next(color)
    
    #Lance moving with the VBs
    mi = np.argmin(abs(np.asarray(dtb)-date[i]))
    Lance_lon = llon_dtb[mi]
    Lance_lat = llat_dtb[mi]
    xl, yl = m(Lance_lon, Lance_lat)
    ax.plot(xl,yl,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')
    bx.plot(xl,yl,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')
    cx.plot(xl,yl,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')

    outname='virtual_buoys_'+date[i].strftime('%Y%m%d')
    print(outname)
    fig1.savefig(outpath+outname,bbox_inches='tight')
    plt.close()


#total damage plot
#date = date[-1]
damage = np.ma.masked_invalid(damage).filled(fill_value=0)
dtotal = np.sum(damage[:fp+1,:],axis=0)

leads = np.ma.masked_invalid(leads).filled(fill_value=0)
ltotal = np.sum(leads[:fp+1,:],axis=0)

ridges = np.ma.masked_invalid(ridges).filled(fill_value=0)
rtotal = np.sum(ridges[:fp+1,:],axis=0)


classified = np.where((ltotal<.3)&(rtotal>.3),1,0)
classified = np.where((ltotal>.3)&(rtotal<.3),2,classified)
classified = np.where((ltotal>.3)&(rtotal>.3),3,classified)

radius_proj=radius+1000  #extend to account of VBs that drifted out
xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)
area_extent = (xlp-radius_proj,ylp-radius_proj,xlp+radius_proj,ylp+radius_proj) #just the last scene, use last Lance location from previous scene
area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)

print(date[fp])
print(Lance_lon, Lance_lat)
print(xlp,ylp)

#m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
#m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
fig2    = plt.figure(figsize=(40,10))
ax      = fig2.add_subplot(141)
bx      = fig2.add_subplot(142)
cx      = fig2.add_subplot(143)
dx      = fig2.add_subplot(144)

ax.set_title('total damage')
bx.set_title('leads')
cx.set_title('ridges')
dx.set_title('classified damage')


#m = pr.plot.area_def2basemap(area_def)

#m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
#m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

#x,y = m(lon_path[-1,:,:],lat_path[-1,:,:])
#ax.pcolormesh(x,y,dtotal)

lop = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),lon_path[fp,:]); lop = np.ma.compressed(lop)
lap = np.ma.masked_invalid(lat_path[fp,:]); lap = np.ma.compressed(lap)
x,y = m(lop,lap)

dd = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),dtotal); dd = np.ma.compressed(dd)
ax.scatter(x,y,c=dd,lw=5)

ll = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),ltotal); ll = np.ma.compressed(ll)
bx.scatter(x,y,c=ll,lw=5)

rr = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),rtotal); rr = np.ma.compressed(rr)
cx.scatter(x,y,c=rr,lw=5)

cc = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),classified); cc = np.ma.compressed(cc)
scatter = dx.scatter(x,y,c=cc,lw=5)

# produce a legend with the unique colors from the scatter
legend1 = dx.legend(*scatter.legend_elements(),
                    loc="lower right", title="Classes")
dx.add_artist(legend1)


xl, yl = m(Lance_lon, Lance_lat)
ax.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
bx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
cx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
dx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

##scale
#m.drawmapscale(Lance_lon, Lance_lat-.3, Lance_lon+8, Lance_lat-.2, 50, units='km', barstyle='fancy',fontsize=14)
start = date[0].strftime('%Y%m%d')
end = date[fp].strftime('%Y%m%d%H%M%S')
outname='virtual_buoys_classified_'+start+'_'+end+'_'+str(int(radius/1000))+'km'
fig2.savefig(outpath+outname,bbox_inches='tight')
plt.close()

#interpolate these parcels to a regular grid, so that they can be saved as geotiff
swath_def = pr.geometry.SwathDefinition(lons=lop, lats=lap)
classified_map = pr.kd_tree.resample_gauss(swath_def, cc,area_def, radius_of_influence=1000, sigmas=500)

#save this a geotiff
geotiff_file = outpath+outname+'.tiff'
save_geotiff(classified_map, area_def, geotiff_file)


#convert -delay 50 /Data/sim/polona/sid/deform/plots/virtual_buoys_* /Data/sim/polona/sid/deform/plots/vb_anim.gif

