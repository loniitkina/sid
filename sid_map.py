from datetime import datetime
import numpy as np
import pyresample as pr
from sid_func import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

outpath_def = '/Data/sim/polona/sid/deform/'
outpath = outpath_def+'plots/'
metfile = '../data/10minute_nounits.csv'
radius = 75000


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




#make map of Lance, SAR image box, ship radar box and buoy positions
#also make polygons of FYI and SYI 

fig1    = plt.figure(figsize=(20,20))
ax      = fig1.add_subplot(111)
#cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
#cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

#area_def = pr.utils.load_area('area.cfg', proj)  
m = pr.plot.area_def2basemap(area_def)

#scale
m.drawmapscale(Lance_lon, Lance_lat-.3, Lance_lon+8, Lance_lat-.2, 50, units='km', barstyle='fancy',fontsize=14)
m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])

#Lance
xl, yl = m(Lance_lon, Lance_lat)
ax.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

#box for SAR data
#get 100km box around Lance
xbox = [xl-50000, xl+50000, xl+50000, xl-50000]
ybox = [yl-50000, yl-50000, yl+50000, yl+50000]
xybox = zip(xbox,ybox)
xybox = list(xybox)     #work around of Python 3 code
poly = Polygon( xybox, edgecolor='k', alpha=1, fill=False, linewidth=3)
plt.gca().add_patch(poly)

#box for ship radar data
#get 14km box around Lance
xbox = [xl-7000, xl+7000, xl+7000, xl-7000]
ybox = [yl-7000, yl-7000, yl+7000, yl+7000]
xybox = zip(xbox,ybox)
xybox = list(xybox)     #work around of Python 3 code
poly = Polygon( xybox, edgecolor='purple', alpha=1, fill=False, linewidth=3)
plt.gca().add_patch(poly)

#buoys



fig1.savefig(outpath+'sid_map',bbox_inches='tight')
  
