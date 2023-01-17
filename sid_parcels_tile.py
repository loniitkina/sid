from datetime import datetime
from glob import glob
import numpy as np
import pyresample as pr
from sid_func import *
import matplotlib.pyplot as plt

#for parcel tracking we need to have consequtive data: second scene in pair needs to be first scene in the next pair! (combo option is not possible here)

#WARNING: how accurate are the input S-1 data? Geolocation in GRD is still approximate... Use snap ESA to imporve geolocaton for s-1
#This option not available for RS-2, ask KSAT

#TASK: feed your own tiffs to Sea_ice_drift (change geolocation, apply calibration, apply noise correction, use DB instead of linar scale for intensity)
#algorithm works with open CV for feature detection

#TASK: make flowchart for parcels

#WARNING: CHECK what is the time difference between the parts of the mosaic and use different drift products for different parts of the mosaic!!!

radius = 60000          #use a radius shorter than in the parcel producing script - or lots of boundary parcels will be lost fast!
radius = 100000
spacing= 300            #at spacing shorter than 500m, there will be lots of nans...
spacing= 800

#for output
#fp = 6  #final image pair
#fp = 24	#in MOSAiC leg 3 events there is 24 pairs

#-------------------------------------------------------------------
inpath = '/scratch/pit000/results/sid/parcels/'
inpath = '/scratch/pit000/results/sid/deform200km/'

inpath_drift = '/scratch/pit000/results/sid/drift/stp10_factor05/'
outpath_data = inpath
outpath = inpath
outpath = inpath

resolution=str(spacing)	#for naming of output, this is combination of drift/parcel input and grid spacing

#get ship location
shipfile = '../../downloads/position_leg3_nh-track.csv'	#leg3 (and transition to leg 4 until 6 June)
ship_time = getColumn(shipfile,0)
ship_time = [ datetime.strptime(ship_time[i], "%Y-%m-%d %H:%M:%S") for i in range(len(ship_time)) ]
ship_lon = np.asarray(getColumn(shipfile,1),dtype=float)
ship_lat = np.asarray(getColumn(shipfile,2),dtype=float)

#projection conversion parameters
from pyproj import Proj, transform
inProj = Proj(init='epsg:4326')
#Use same projection as in pattern matching part of the sea ice drift algorithm
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +ellps=WGS84 +units=m' % (90, 10))
ship_x,ship_y = transform(inProj,outProj,ship_lon, ship_lat)

#read in the drift product and calculate displacements from starting point until the first break in the SAR data (1 week)
#fl = sorted(glob(inpath_drift+'SeaIceDrift_2020*.npz'))[1:]

##some dates
#start = datetime(2015,1,21,6,54,4)              #this coincides with the 1st SAR image date/time
#start = datetime(2020,3,14,12,29,54)		#MOSAiC leg 3 events
#start = datetime(2020,3,15,11,32,28)		#MOSAiC leg 3 events - day with better initial coverage   
#start = datetime(2020,3,12,11,32,28)		#MOSAiC leg 3 events - day with reasonable initial coverage  

fl = sorted(glob(inpath+'Damage_2020*.npz'))[5:65]  #this should work for leg 3
#print(fl)
#exit()
start = fl[0].split('_')[-3]
start = datetime.strptime(start, "%Y%m%dT%H%M%S")
print(start)

#seed the parcels from ship position
def seed_parcels(target_date,base_time,base_x,base_y,spacing):
    mi = np.argmin(abs(np.asarray(base_time)-target_date))
    target_x = base_x[mi]
    target_y = base_y[mi]

    #create a grid of spacing of X m in each direction in the radius from the ship's initial position
    xbf = np.arange(target_x-radius, target_x+radius+spacing, spacing)    #add some space over the edge (will be nan immediately at first step: to keep space for new parcels)
    ybf = np.arange(target_y-radius, target_y+radius, spacing)    #to increase resolution use fractions of spacing: e.g. spacing/3!
    x_buoy,y_buoy = np.meshgrid(xbf,ybf)
    #flatten
    x_buoy = x_buoy.flatten()
    y_buoy = y_buoy.flatten()
    
    return(x_buoy,y_buoy)

x_buoy,y_buoy = seed_parcels(start,ship_time,ship_x,ship_y,spacing)



#empty arrays to store parcel drifts and damage
x_path = np.zeros((len(fl)+1,x_buoy.shape[0]))
y_path = np.zeros((len(fl)+1,x_buoy.shape[0]))
damage = np.zeros((len(fl)+1,x_buoy.shape[0]))
leads = np.zeros((len(fl)+1,x_buoy.shape[0]))
ridges = np.zeros((len(fl)+1,x_buoy.shape[0]))
age = np.zeros((len(fl)+1,x_buoy.shape[0]))
 
#start with initial coordinates
x_path[0,:] = x_buoy
y_path[0,:] = y_buoy
date = [start]

for i in range(0,len(fl)):
    #break
    ##read in the drift data
    #print(fl[i])
    #container = np.load(fl[i])
    #lat1 = container['lat1']     #these arrays are as big as both scenes together (nan-where no overlap), ~5000x5000 of 80m-pixels = 400x400km
    #lon1 = container['lon1']
    #u = container['upm']
    #v = container['vpm']
    #hpm = container['hpm'] 
        
    ##quality filter 
    ##strict quality filter (hpm > 9) will leave holes in deforming areas
    #gpi = hpm > 4
    #u = u[gpi]
    #v = v[gpi]
    #lon1 = lon1[gpi]
    #lat1 = lat1[gpi]
    
    #get damage data
    print(fl[i])
    container = np.load(fl[i])
    lat_dmg = container['lat']      #triangle centroids
    lon_dmg = container['lon']
    dmg = container['d']
    rid = container['r']
    lea = container['l']
    
    u = container['u']
    v = container['v']
    dt = container['dt']
    
    #project data coordinates
    #x1,y1 = transform(inProj,outProj,lon1,lat1)
    x1,y1 = transform(inProj,outProj,lon_dmg,lat_dmg)
    
    #get dates
    date1 = fl[i].split('_')[-3]
    date2 = fl[i].split('_')[-2].split('.')[0]
    dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    #diff = (dt2-dt1).total_seconds()
    #print(diff)

    #get displacements
    dx = dt*u; dy = dt*v
        
    #estimate positions after drift
    x2 = x1 + dx
    y2 = y1 + dy
    
    #save the date 
    date.append(dt2)
    
    #find closest drift/damage coordinate for each parcel
    print(x_buoy.shape[0])
    for m in range(0,x_buoy.shape[0]):
        #print(m)
        
        #skip this buoy if already nan
        if np.isnan(x_buoy[m]) and np.isnan(damage[i,m]):
            continue
        
        ##sea ice drift/displacement mask for data with 400m resolution, 600m spacing between buoys, 300 m drift precission (based on buoys)    
        #mask = (x1>x_buoy[m]-spacing) & (x1<x_buoy[m]+spacing) & (y1>y_buoy[m]-spacing) & (y1<y_buoy[m]+spacing)
        #damage mask for data with ~400m resolution            
        mask = (x1>x_buoy[m]-spacing/2) & (x1<x_buoy[m]+spacing/2) & (y1>y_buoy[m]-spacing/2) & (y1<y_buoy[m]+spacing/2)
        
        #get the mean location for parcel (and update parcel location for the search in next step)
        #print(np.ma.compressed(np.ma.masked_invalid(x2[mask])))
        x_buoy[m] = np.mean(np.ma.masked_invalid(x2[mask]))
        y_buoy[m] = np.mean(np.ma.masked_invalid(y2[mask]))
        x_path[i+1,m] = x_buoy[m]   #step 0 is initial condition, step 1 is reached after one day etc.
        y_path[i+1,m] = y_buoy[m]
        
        #check if there was any deformation
        #print(np.ma.compressed(np.ma.masked_invalid(dmg[mask_d])))
        dd = np.mean(np.ma.masked_invalid(dmg[mask])) #this is just 0 and 1 data, if empty, we get a nan and loose the buoy...
        damage[i+1,m] = dd
        ll = np.mean(np.ma.masked_invalid(lea[mask]))
        leads[i+1,m] = ll
        rr = np.mean(np.ma.masked_invalid(rid[mask]))
        ridges[i+1,m] = rr
        
        #parcel age
        #age[i+1,m] = age[i,m]+1
        
    ##fill empty area with undamaged parcels - leads
    ##make a new mesh and check for every new buoy if it has any buoy in the radius
    ##if not add this buoy
    ##track the parcel age
    #x_seed,y_seed = seed_parcels(dt1,ship_time,ship_x,ship_y,spacing)  
    #for m in range(0,x_seed.shape[0]):
        
        #mask = (x_buoy>x_seed[m]-spacing/2) & (x_buoy<x_seed[m]+spacing/2) & (y_buoy>y_seed[m]-spacing/2) & (y_buoy<y_seed[m]+spacing/2)
        
        #if len(mask)==0:
            #print('add new parcel')
            ##we add coordinates for this step. All previous time steps (i) are nan or masked
            
            
            #TASK: extend m dimension - this means that for next time step i, the m dimension (number of parcels will be longer)
            ##concatenate path, damage, leads, ridges, age
            
            #parcel age
            #age[i+1,m] = 0

        
        #how can i track the area change over this triangle???
        #for that i need to always use the same nods and keep record of the area
        
        #this is only possible if output points from one drift pair are kept as input for the other.
        #and some mew points are added
        
        #or simlpy, why not claculating deformaiton again the parcel script, there cells are tracked...
        
        #then track volume of lead ice and ridge ice in every parcel...
        #no, parcels are just points, how should i distribute the values from triangles to them???
        
        #use some simple themodynamical model for sea ice growth inside the parcel
        #dynamic components based on area change creates leads and ridges
        #these new, ridge ice needs to be distributed somehow inside the parcels
        #what can we learn from Albeldyl et al, TC???
        
            
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


#print(lat_path[:,10])        #no more data on image 12
#exit()


out_file = outpath_data+'VB_dates.npz'
container = np.load(out_file, allow_pickle=True)
date = container['date']

#radius_proj=55000           #big enough to contain all the moving ice
radius_proj=radius+2000  #extend to account of VBs that drifted out
from pyresample.geometry import AreaDefinition
#Using a projection dictionary
area_id = 'around ship'
description = 'North Pole LAEA'
proj_id = 'lance'
proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'datum':'WGS84', 'ellps':'WGS84', 'units':'m'}
width = radius_proj*2/600 #100 m spacing    (should have same resolution as the final gridded map)
height = radius_proj*2/600 #100 m spacing
#area_extent = (xlp-radius_proj,ylp-radius_proj,xlp+radius_proj,ylp+radius_proj)
#area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
#print(area_def)
#m = pr.plot.area_def2basemap(area_def)

color=iter(plt.cm.rainbow(np.linspace(0,1,lon_path.shape[0])))

for i in range(0,lon_path.shape[0]):
    print(i)

    #moving projection with the ship
    mi = np.argmin(abs(np.asarray(dtb)-date[i]))
    ship_lon = llon_dtb[mi]
    ship_lat = llat_dtb[mi]
    xl, yl = transform(inProj,outProj,ship_lon, ship_lat)


    #Using a projection dictionary
    area_extent = (xl-radius_proj,yl-radius_proj,xl+radius_proj,yl+radius_proj)
    area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
    #print(area_def)
    m = pr.plot.area_def2basemap(area_def)


    fig1    = plt.figure(figsize=(30,10))
    ax      = fig1.add_subplot(131)
    #uncomment here if you want it with geographical coordinates (comment to get figure coordinates - and all data)
    m = pr.plot.area_def2basemap(area_def)                                      #note: this is static/same for all days!!!
    m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

    bx      = fig1.add_subplot(132)
    #uncomment here if you want it with geographical coordinates
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

    cx      = fig1.add_subplot(133)
    #uncomment here if you want it with geographical coordinates
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

    
    ax.set_title('damage')
    bx.set_title('leads')
    cx.set_title('ridges')
    
    lop = np.ma.masked_where(~(np.isfinite(lat_path[i,:])),lon_path[i,:]); lop = np.ma.compressed(lop)
    lap = np.ma.masked_invalid(lat_path[i,:]); lap = np.ma.compressed(lap)
    x,y = m(lop,lap)
    
    dd = np.ma.masked_where(~(np.isfinite(lat_path[i,:])),damage[i,:]); dd = np.ma.compressed(dd)
    rr = np.ma.masked_where(~(np.isfinite(lat_path[i,:])),ridges[i,:]); rr = np.ma.compressed(rr)
    ll = np.ma.masked_where(~(np.isfinite(lat_path[i,:])),leads[i,:]); ll = np.ma.compressed(ll)
    
    ax.scatter(x,y,c=dd,lw=2)         #marker size should be radius dependant
    bx.scatter(x,y,c=ll,lw=2)
    cx.scatter(x,y,c=rr,lw=2)

    cl = next(color)
    
    #ship moving with the VBs
    ax.plot(xl,yl,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')
    bx.plot(xl,yl,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')
    cx.plot(xl,yl,'*',markeredgewidth=2,color=cl,markersize=20,markeredgecolor='k')

    outname='virtual_buoys_'+date[i].strftime('%Y%m%d')
    print(outname)
    fig1.savefig(outpath+outname,bbox_inches='tight')
    plt.close()


#total damage plot - best coverage is x days before the end
fp=-3	#March April event: 6 April
fp=-33	#whole spring: 30 April there is about half parcels left
print(date[fp])

damage = np.ma.masked_invalid(damage).filled(fill_value=0)
dtotal = np.sum(damage[:fp,:],axis=0)

leads = np.ma.masked_invalid(leads).filled(fill_value=0)
ltotal = np.sum(leads[:fp,:],axis=0)

ridges = np.ma.masked_invalid(ridges).filled(fill_value=0)
rtotal = np.sum(ridges[:fp,:],axis=0)


classified = np.where((ltotal<.3)&(rtotal>.3),1,0)
classified = np.where((ltotal>.3)&(rtotal<.3),2,classified)
classified = np.where((ltotal>.3)&(rtotal>.3),3,classified)

radius_proj=radius_proj+10000  #extend to account of VBs that drifted out
width = radius_proj*2/600 #100 m spacing    (should have same resolution as the final gridded map)
height = radius_proj*2/600 #100 m spacing

#xlp,ylp = transform(inProj,outProj,ship_lon, ship_lat)
#area_extent = (xl-radius_proj,ylp-radius_proj,xl+radius_proj,ylp+radius_proj) #just the last scene, use last ship location from previous scene
#area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)

#print(date[fp])
#print(ship_lon, ship_lat)
#print(xlp,ylp)

#m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
#m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
fig2    = plt.figure(figsize=(40,10))
ax      = fig2.add_subplot(141)

mi = np.argmin(abs(np.asarray(dtb)-date[fp]))
ship_lon = llon_dtb[mi]
ship_lat = llat_dtb[mi]
xl, yl = transform(inProj,outProj,ship_lon, ship_lat)
#Using a projection dictionary
area_extent = (xl-radius_proj,yl-radius_proj,xl+radius_proj,yl+radius_proj)
area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)


m = pr.plot.area_def2basemap(area_def)
m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

bx      = fig2.add_subplot(142)
m = pr.plot.area_def2basemap(area_def)
m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

cx      = fig2.add_subplot(143)
m = pr.plot.area_def2basemap(area_def)
m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

dx      = fig2.add_subplot(144)
m = pr.plot.area_def2basemap(area_def)
m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])


ax.set_title('total damage')
bx.set_title('leads')
cx.set_title('ridges/rafted lead ice')
dx.set_title('classified damage')



#x,y = m(lon_path[-1,:,:],lat_path[-1,:,:])
#ax.pcolormesh(x,y,dtotal)

lop = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),lon_path[fp,:]); lop = np.ma.compressed(lop)
lap = np.ma.masked_invalid(lat_path[fp,:]); lap = np.ma.compressed(lap)
x,y = m(lop,lap)

dd = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),dtotal); dd = np.ma.compressed(dd)
heat = ax.scatter(x,y,c=dd,lw=2)
cb = plt.colorbar(heat, ax=ax, pad=.01, orientation="horizontal")
cb.set_label(label='Number of days',fontsize=20)
cb.ax.tick_params(labelsize=20)

ll = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),ltotal); ll = np.ma.compressed(ll)
heat = bx.scatter(x,y,c=ll,lw=2)
cb = plt.colorbar(heat, ax=bx, pad=.01, orientation="horizontal")
cb.set_label(label='Number of days',fontsize=20)
cb.ax.tick_params(labelsize=20)

rr = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),rtotal); rr = np.ma.compressed(rr)
heat = cx.scatter(x,y,c=rr,lw=2)
cb = plt.colorbar(heat, ax=cx, pad=.01, orientation="horizontal")
cb.set_label(label='Number of days',fontsize=20)
cb.ax.tick_params(labelsize=20)

cc = np.ma.masked_where(~(np.isfinite(lat_path[fp,:])),classified); cc = np.ma.compressed(cc)
scatter = dx.scatter(x,y,c=cc,lw=2)

# produce a legend with the unique colors from the scatter
legend1 = dx.legend(*scatter.legend_elements(),
                    loc="lower right", title="Classes")
dx.add_artist(legend1)


xl, yl = m(ship_lon, ship_lat)
ax.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
bx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
cx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
dx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

##scale
#m.drawmapscale(ship_lon, ship_lat-.3, ship_lon+8, ship_lat-.2, 50, units='km', barstyle='fancy',fontsize=14)
start = date[0].strftime('%Y%m%d')
end = date[fp].strftime('%Y%m%d%H%M%S')
outname='virtual_buoys_classified_wholespring_'+resolution+'_m_'+start+'_'+end+'_'+str(int(radius/1000))+'km'
print(outpath+outname)
fig2.savefig(outpath+outname,bbox_inches='tight')
plt.close()

#just classified parcels for Jack's YRT proposal

#a bit of extra fiddling with the projection
#radius_proj=radius+2000
#area_extent = (xlp-radius-10000,ylp-radius-10000,xlp+radius+20000,ylp+radius+10000)
#proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'datum':'WGS84', 'ellps':'WGS84', 'units':'m'}
#area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)

fig2a = plt.figure(figsize=(10,10))
ax    = fig2a.add_subplot(111)
m = pr.plot.area_def2basemap(area_def)
m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

#try to use a color map where 'no deformation' is white

scatter = ax.scatter(x,y,c=cc,lw=1,s=20)
handles,labels=scatter.legend_elements()
legend1 = ax.legend(handles=handles,labels=['No deformation','Mainly convergence','Mainly divergence','Mixed deformation'],
                    loc="upper left",fontsize=12,title='Deformation 15 Mar-30 Apr 2020',title_fontsize=14)
ax.add_artist(legend1)
print(outpath+'final_map')
fig2a.savefig(outpath+'final_map'+resolution+'_m_wholespring',bbox_inches='tight')
plt.close()
exit()

#interpolate these parcels to a regular grid, so that they can be saved as geotiff
swath_def = pr.geometry.SwathDefinition(lons=lop, lats=lap)
#classified_map = pr.kd_tree.resample_gauss(swath_def, cc,area_def, radius_of_influence=600, sigmas=500)
##this will give interpolated values between the classes
##reclasify back to integer values
##classified_map = np.where(classified_map<0.1,0,classified_map)
#classified_map = np.where((classified_map>0) &(classified_map<1.1),1,classified_map)
#classified_map = np.where((classified_map>1.1) & (classified_map<2.1),2,classified_map)
#classified_map = np.where(classified_map>2.1,3,classified_map)

#rather use the nearest neighbor
classified_map = pr.kd_tree.resample_nearest(swath_def, cc,area_def, radius_of_influence=600)

#save as png
fig3    = plt.figure(figsize=(10,10))
ax      = fig3.add_subplot(111)
ax.imshow(classified_map)
fig3.savefig(outpath+outname+'_nn',bbox_inches='tight')

#save as geotiff
geotiff_file = outpath+outname+'.tiff'
save_geotiff(classified_map, area_def, geotiff_file)


#convert -delay 50 /Data/sim/polona/sid/deform/plots/virtual_buoys_* /Data/sim/polona/sid/deform/plots/vb_anim.gif

