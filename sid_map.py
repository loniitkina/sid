from datetime import datetime
import numpy as np
import pyresample as pr
from sid_func import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

inpath = '/Data/sim/polona/sid/buoys/'
outpath = '/Data/sim/polona/sid/deform/plots/'
metfile = '../data/10minute_nounits.csv'
radius = 51000

maptime = datetime(2015,1,27,12,0,0)

#Lance postion (from Lance's met system)
mettime = getColumn(metfile,0)
dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
mi = np.argmin(abs(np.asarray(dtb)-maptime))
Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]

print(Lance_lat,Lance_lon)

#Lance centered projection
from pyproj import Proj, transform
inProj = Proj(init='epsg:4326')
outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))
xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)

from pyresample.geometry import AreaDefinition
#Using a projection dictionary
area_id = 'around Lance'
description = 'North Pole LAEA Europe'
proj_id = 'lance'
proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'a':6378137, 'b':6356752.3142, 'units':'m'}
width = radius*2/100 #100 m spacing
height = radius*2/100 #100 m spacing
area_extent = (xlp-radius,ylp-radius,xlp+radius,ylp+radius)
area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)
#print(area_def)

#Plotting
fig1    = plt.figure(figsize=(20,20))
ax      = fig1.add_subplot(111)
#cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
#cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

#area_def = pr.utils.load_area('area.cfg', proj)  
m = pr.plot.area_def2basemap(area_def)

#scale
m.drawmapscale(Lance_lon-.5, Lance_lat-.08, Lance_lon+.5, Lance_lat-.1, 10, units='km', barstyle='fancy',fontsize=14)
m.drawmeridians(np.arange(0.,360.,2.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,.5),labels=[1,0,0,0])

#Lance
xl, yl = m(Lance_lon, Lance_lat)
ax.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k',label='RV Lance')

#box for SAR data
#get 100km box around Lance
xbox = [xl-50000, xl+50000, xl+50000, xl-50000]
ybox = [yl-50000, yl-50000, yl+50000, yl+50000]
xybox = zip(xbox,ybox)
xybox = list(xybox)     #work around of Python 3 code
poly = Polygon( xybox, edgecolor='k', alpha=1, fill=False, facecolor='.2', linewidth=10, ls='--', label='SAR')
plt.gca().add_patch(poly)

#get 7.5km box in the SYI
ys = yl+25000   #shift center 25km north
xbox = [xl-7500, xl+7500, xl+7500, xl-7500]
ybox = [ys-7500, ys-7500, ys+7500, ys+7500]
xybox = zip(xbox,ybox)
xybox = list(xybox)     #work around of Python 3 code
poly = Polygon( xybox, edgecolor='darkred', alpha=.5, fill=True, facecolor='darkred', linewidth=10, label='SAR SYI')
plt.gca().add_patch(poly)

#get 7.5km box in the FYI
yf = yl-25000   #shift center 25km south
xbox = [xl-7500, xl+7500, xl+7500, xl-7500]
ybox = [yf-7500, yf-7500, yf+7500, yf+7500]
xybox = zip(xbox,ybox)
xybox = list(xybox)     #work around of Python 3 code
poly = Polygon( xybox, edgecolor='gold', alpha=.5, fill=True, facecolor='gold', linewidth=10, label='SAR FYI')
plt.gca().add_patch(poly)

#box for ship radar data
#get 14km box around Lance
xbox = [xl-7500, xl+7500, xl+7500, xl-7500]
ybox = [yl-7500, yl-7500, yl+7500, yl+7500]
xybox = zip(xbox,ybox)
xybox = list(xybox)     #work around of Python 3 code
poly = Polygon( xybox, edgecolor='purple', alpha=.5, fill=True, facecolor='purple', linewidth=6, label='Ship Radar')
plt.gca().add_patch(poly)

#buoys
fn = inpath+'buoys.csv'
blon = np.asarray(getColumn(fn,2),dtype=float)
blat = np.asarray(getColumn(fn,3),dtype=float)
print(blon,blat)
xb,yb = m(blon,blat)
ax.plot(xb,yb,'o',markeredgewidth=2,color='royalblue',markersize=20,markeredgecolor='k')
ax.plot(xb[0],yb[0],'o',markeredgewidth=2,color='royalblue',markersize=20,markeredgecolor='k', label='Buoys')

#use deformation plot for background
            ###Plotting 
            ##deform = div*1e6
            ###print(deform)
            ##outname = 'map_div'+reg+'_L'+str(j)+'_'+date1
            ##label = r'Divergence (10$^6$s$^{-1}$)'
            ##interval = [-5, 5]
            ##cmap=plt.cm.bwr
            ###title = 'c'
            
            #deform = shr*1e6
            #outname = 'map_shr_'+reg+'_L'+str(j)+'_'+date1
            #label = r'Shear (10$^6$s$^{-1}$)'
            #interval = [0, 10]
            #cmap=plt.cm.Reds

            ##speed = np.sqrt(u**2+v**2)
            ##outname = 'map_speed_'+reg+'_L'+str(j)+'_'+date1
            ##label = r'Speed (m/s)'
            ##cmap=plt.cm.nipy_spectral
                        
            #print(outname)
            ##deformation plots
            #fig3    = plt.figure(figsize=(20,20))
            #cx      = fig3.add_subplot(111)
            ##cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
            ##cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

            ##area_def = pr.utils.load_area('area.cfg', proj)  
            #m = pr.plot.area_def2basemap(area_def)
            
            ##scale
            #m.drawmapscale(Lance_lon, Lance_lat-.3, Lance_lon+8, Lance_lat-.2, 50, units='km', barstyle='fancy',fontsize=14)
            #m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
            #m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])


            ##nods
            ##cx.plot(x,y,'o',linewidth=2,color='purple')
            
            ###speed
            ##sc = cx.scatter(x,y,s=50,c=speed,cmap=cmap, vmin=.06, vmax=.085)         #add colorbar and remove extreme values
            ##cbar = plt.colorbar(sc)
            ##cbar.ax.set_ylabel('Drift speed (m/s)',size=22)

            
            ##triangles
            #patches = []
            #for k in range(deform.shape[0]):
                #patch = Polygon(tripts[k,:,:], edgecolor='orchid', alpha=1, fill=False)
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
            
            ##virtual buoys
            #xb,yb = m(lon_path[i,:,:],lat_path[i,:,:])
            #cx.plot(xb,yb,'o',linewidth=2,color='purple')


#legend
ax.legend(loc=0,fancybox=True,fontsize='xx-large')

fig1.savefig(outpath+'sid_map',bbox_inches='tight')
  
