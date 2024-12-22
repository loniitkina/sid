import os
from glob import glob
from datetime import datetime
import numpy as np
from shapely.geometry import Point, MultiPoint, MultiPolygon
from shapely.geometry import Polygon as Shapely_Polygon
from shapely.ops import unary_union
import pickle
from sid_func import * 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import AnchoredText

#plotting TS created by sid_afs_stats.py

inpath = '/scratch/pit000/results/sid/deform200km/'
plotpath = '/scratch/pit000/results/sid/plots200km_revision/'

#MOSAiC

#N-ICE
#fname='/scratch/pit000/results/sid/deform200km/cde_ship2015_200km_thfilter_lkffilter.csv'
#fname='/scratch/pit000/results/sid/deform200km/cde_8km_ship2015_200km_thfilter_lkffilter.csv'
fname='/scratch/pit000/results/sid/deform200km/cde_4km_ship2015_200km_thfilter_lkffilter.csv'
start=datetime(2015, 1, 16)
end=datetime(2015, 2, 27)

#CIRFA

#load data
date = getColumn(fname,0)
cde_date = [ datetime.strptime(date[i], "%Y-%m-%d %H:%M:%S") for i in range(len(date)) ]

cde_num = np.asarray(getColumn(fname,1),dtype=float)*10000
cde_area_m = np.asarray(getColumn(fname,2),dtype=float)
cde_area_s = np.asarray(getColumn(fname,3),dtype=float)
cde_rr_m = np.asarray(getColumn(fname,4),dtype=float)
cde_rr_s = np.asarray(getColumn(fname,5),dtype=float)
cde_fr_m = np.asarray(getColumn(fname,6),dtype=float)
cde_fr_s = np.asarray(getColumn(fname,7),dtype=float)
cde_cc_m = np.asarray(getColumn(fname,8),dtype=float)
cde_cc_s = np.asarray(getColumn(fname,9),dtype=float)
cde_lkfa = np.asarray(getColumn(fname,10),dtype=float)
cde_lkfa_max = np.asarray(getColumn(fname,11),dtype=float)
cde_lkfd_min_m = np.asarray(getColumn(fname,12),dtype=float)
cde_lkfd_min_s = np.asarray(getColumn(fname,13),dtype=float)
cde_lkfd_max_m = np.asarray(getColumn(fname,14),dtype=float)
cde_lkfd_max_s = np.asarray(getColumn(fname,15),dtype=float)

#WARNING:invert this temporary (will be fixed in next version of sid_afs_stats.py)
cde_rr_m=(cde_rr_m)**(-1)
cde_rr_s=(cde_rr_s)**(-1)

#time series of cde satistics
fig1    = plt.figure(figsize=(14,12))
ax      = fig1.add_subplot(611)
ax.set_title('a) CDC density',fontsize=18, loc='left')
ax.set_ylabel('$1/{100km^2}$',fontsize=16)
ax.set_ylim(0, 20)
ax.set_xlim(start,end)

bx      = fig1.add_subplot(612)
bx.set_title('b) CDC area',fontsize=18, loc='left')
bx.set_ylabel('$km^2$',fontsize=16)
bx.set_ylim(0, 50000)
bx.set_xlim(start,end)

cx      = fig1.add_subplot(613)
cx.set_title('c) CDC circularity/roundness',fontsize=18, loc='left')
cx.set_ylim(0.45, 1)
cx.set_xlim(start,end)
cx1 = cx.twinx()
cx1.set_ylim(1, 5)

dx      = fig1.add_subplot(614)
dx.set_title('d) CDC complexity',fontsize=18, loc='left')
#dx.set_ylabel('CDC\nfragmentation',fontsize=16)
dx.set_ylim(2.3, 6.5)
dx.set_xlim(start,end)

ex      = fig1.add_subplot(615)
ex.set_title('e) LKF fraction',fontsize=18, loc='left')
ex.set_ylabel('$\%$',fontsize=16)
ex.set_ylim(0, 3)
ex.set_xlim(start,end)
#twin axis
ex1 = ex.twinx()
ex1.set_ylim(0, 50)

fx      = fig1.add_subplot(616)
#min and max floe radius
fx.set_title('f) LKF distance',fontsize=18, loc='left')
fx.set_ylabel('$km$',fontsize=16)
fx.set_ylim(0, 40)
fx.set_xlim(start,end)
#twin axis
fx1 = fx.twinx()
fx1.set_ylim(0, 200)

#highlight MAJOR storms (all based on Lana's storm table, except the first one which is based on temperature above -20)
plots = [ax,bx,cx,dx,ex,fx]
[pp.axvspan(datetime(2015,1,21,15,0), datetime(2015,1,22,15,0), facecolor='0.7', alpha=0.4, linewidth=0) for pp in plots]
[pp.axvspan(datetime(2015,2,3,11,0), datetime(2015,2,8,21,0), facecolor='0.7', alpha=0.4, linewidth=0) for pp in plots]
[pp.axvspan(datetime(2015,2,15,12,0), datetime(2015,2,21,4,0), facecolor='0.7', alpha=0.4, linewidth=0) for pp in plots]

#annotate the storms
ax.text(datetime(2015,1,21,18,0),15, 'A', fontsize=20)
ax.text(datetime(2015,2,3,13,0),15, 'B', fontsize=20)
ax.text(datetime(2015,2,15,15,0),15, 'C', fontsize=20)

#time series of statistics
ax.scatter(cde_date,cde_num, color='indigo')

bx.plot(cde_date,cde_area_m,linestyle='None',marker='o',color='purple')
bx.fill_between(cde_date, cde_area_m-cde_area_s, cde_area_m+cde_area_s,alpha=.2, color='purple')

cx.plot(cde_date,cde_cc_m,linestyle='None',marker='o',color='teal')
cx.fill_between(cde_date, cde_cc_m-cde_cc_s, cde_cc_m+cde_cc_s,alpha=.2, color='teal')

cx1.plot(cde_date,cde_rr_m,linestyle='None',marker='o',color='0.7')
cx1.fill_between(cde_date, cde_rr_m-cde_rr_s, cde_rr_m+cde_rr_s,alpha=.2, color='0.7')

dx.plot(cde_date,cde_fr_m,linestyle='None',marker='o',color='royalblue')
dx.fill_between(cde_date, cde_fr_m-cde_fr_s, cde_fr_m+cde_fr_s,alpha=.2, color='royalblue')

ex.scatter(cde_date,cde_lkfa, color='gold')
ex1.scatter(cde_date,cde_lkfa_max, color='hotpink')

fx.plot(cde_date,cde_lkfd_min_m,linestyle='None',marker='o',color='gold')
fx.fill_between(cde_date, cde_lkfd_min_m-cde_lkfd_min_s, cde_lkfd_min_m+cde_lkfd_min_s,alpha=.2, color='gold')
#twin axis
fx1.plot(cde_date,cde_lkfd_max_m,linestyle='None',marker='o',color='hotpink')
fx1.fill_between(cde_date, cde_lkfd_max_m-cde_lkfd_max_s, cde_lkfd_max_m+cde_lkfd_max_s,alpha=.2, color='hotpink')

fig1.autofmt_xdate()
fig1.tight_layout(h_pad=2)  #make space for titles

#ax.xaxis.set_major_formatter(ticker.NullFormatter())
#ax.xaxis.set_minor_formatter(dates.DateFormatter('%b'))

plt.show()


#save cde timeseries
ss = datetime.strftime(start, "%Y-%m-%d")
ee = datetime.strftime(end, "%Y-%m-%d")
plotname=plotpath+'cde_ts_'+ss+'_'+ee+'_replot1'
print(plotname)
fig1.savefig(plotname,bbox_inches='tight')




