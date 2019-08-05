from datetime import datetime
import numpy as np
import pyresample as pr
from sid_func import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

inpath = '/Data/sim/polona/sid/cs/'
inpath_b = '/Data/sim/polona/sid/buoys/'
outpath = '/Data/sim/polona/sid/deform/plots/'
metfile = '../data/10minute_nounits.csv'
radius = 8000

##simple preview
#from osgeo import gdal
#ds = gdal.Open(inpath+'CSM_20150118.tif').ReadAsArray()
#im = plt.imshow(ds)
#plt.show()

#based on this example: https://stackoverflow.com/questions/20488765/plot-gdal-raster-using-matplotlib-basemap
import osr, gdal

# This is the cropped tiff. The original is 2-3 times the size.
# Read the data and metadata
ds = gdal.Open(inpath+'CSM_20150118.tif')

data = ds.ReadAsArray()
print(data.shape)
gt = ds.GetGeoTransform()
proj = ds.GetProjection()

xres = gt[1]
yres = gt[5]
print(xres,yres)

# get the edge coordinates and add half the resolution 
# to go to center coordinates
xmin = gt[0] + xres * 0.5
xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
ymax = gt[3] - yres * 0.5

print(xmin,ymin)
print(xmax,ymax)
ds = None

##alternative tif
#ds = gdal.Open(inpath+'cosmo7.tif')

#data = ds.ReadAsArray()
#print(data.shape)
#gt = ds.GetGeoTransform()
#proj = ds.GetProjection()

#xres = gt[1]
#yres = gt[5]
#print(xres,yres)

## get the edge coordinates and add half the resolution 
## to go to center coordinates
#xmin = gt[0] + xres * 0.5
#xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
#ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
#ymax = gt[3] - yres * 0.5

#print(xmin,ymin)
#print(xmax,ymax)
#ds = None

#exit()

# create a grid of xy coordinates in the original projection
xy_source = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]

print(proj)
print(xy_source)


# Create the figure and basemap object
maptime = datetime(2015,1,18,18,7,0)

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
m.drawmapscale(Lance_lon-.1, Lance_lat-.05, Lance_lon+.3, Lance_lat-.1, 10, units='km', barstyle='fancy',fontsize=14)
m.drawmeridians(np.arange(0.,360.,1.),latmax=90.,labels=[0,0,0,1,])
m.drawparallels(np.arange(79.,90.,.2),labels=[1,0,0,0])

#Lance
xl, yl = m(Lance_lon, Lance_lat)
ax.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

# Create the projection objects for the convertion
# original (Albers)
inproj = osr.SpatialReference()
inproj.ImportFromWkt(proj)

# Get the target projection from the basemap object
outproj = osr.SpatialReference()
outproj.ImportFromProj4(m.proj4string)

# Convert from source projection to basemap projection
xx, yy = convertXY(xy_source, inproj, outproj)

# plot the sat image
im1 = m.pcolormesh(xx, yy, data.T, cmap=plt.cm.Greys_r)

#box for SAR data
#get 100km box around Lance
xbox = [xl-7500, xl+7500, xl+7500, xl-7500]
ybox = [yl-7500, yl-7500, yl+7500, yl+7500]
xybox = zip(xbox,ybox)
xybox = list(xybox)     #work around of Python 3 code
poly = Polygon( xybox, edgecolor='k', alpha=1, fill=False, facecolor='.2', linewidth=10, label='SAR')
plt.gca().add_patch(poly)

#box for ship radar data
#get 14km box around Lance
xbox = [xl-7000, xl+7000, xl+7000, xl-7000]
ybox = [yl-7000, yl-7000, yl+7000, yl+7000]
xybox = zip(xbox,ybox)
xybox = list(xybox)     #work around of Python 3 code
poly = Polygon( xybox, edgecolor='purple', alpha=1, fill=False, facecolor='purple', linewidth=10, label='Ship radar')
plt.gca().add_patch(poly)

#buoys
fn = inpath_b+'buoys.csv'
blon = np.asarray(getColumn(fn,2),dtype=float)
blat = np.asarray(getColumn(fn,3),dtype=float)
print(blon,blat)
xb,yb = m(blon,blat)
ax.plot(xb,yb,'o',markeredgewidth=2,color='royalblue',markersize=20,markeredgecolor='k')
ax.plot(xb[0],yb[0],'o',markeredgewidth=2,color='royalblue',markersize=20,markeredgecolor='k', label='Buoys')

#legend
ax.legend(loc=0,fancybox=True,fontsize='xx-large')

fig1.savefig(outpath+'sid_map_cs',bbox_inches='tight')
  
