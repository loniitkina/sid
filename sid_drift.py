# modified from simple.py (part of SeaIceDrift by Anton Korosov)
import os
import sys
import glob
import unittest
import inspect

import csv
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('Agg')

from nansat import Nansat, Domain, NSR

from sea_ice_drift import SeaIceDrift

# download Sentinel-1 data with sentinelsat
#make json point files with qgis (or write with a script)
#in QGIS: create layer - new temporary scratch layer - draw point - save as
#create user on scihub: https://scihub.copernicus.eu/
#download data to the data partition!
#sentinelsat --user polona --password jank0vna --geometry /mnt/data/sea_ice_drift/lance.geojson --sentinel 1 --producttype GRD --start 20150128 --end 20150202 --download --path /mnt/data/sea_ice_drift/Sentinel1

#retrieve drift
#plot drift and write text output

#retrieve sea ice deformation (next script)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
outpath_drift = '../output/drift_50/el/'
outpath = '../plots/drift_50/el/'
#selection mode (in time)
#all: all images are analyzed
#ee: only image from early morning and early afternoon (same date) - ee and el will often be mixed up (lots of missing images etc) - that is not a problem!
#el: only image from early morning and late afteroon (same date)
#24h: only early mornming images (2 different dates)
mode = ['all', 'ee', 'el', '24h'][2]
print(mode)


# ==== ICE DRIFT RETRIEVAL ====
#inpath = '../data/Sentinel1/'
inpath = '/Data/sim/data/Sentinel1/'

#specify region
regn = 84
regs = 82
regw = 10
rege = 25
#and number of grid points in each direction
#used in post-processing:
#lsc_list = [25,50,100,200,500]   #not ls but number of nominal grid points
#minlen = [4,2,1,.5,.2]
#maxlen = [6,3,1.5,.75,.3]
gridp_we = 50
gridp_sn = gridp_we*1.5         #the region is elongated in NS direction (it needs more points to get nicer triangles/squares)


#show Lance position
def getColumn(filename, column):
    results = csv.reader(open(filename))
    next(results, None)
    return [result[column] for result in results]
metfile = '../data/10minute_nounits.csv'

#file list
fl = sorted(glob.glob(inpath+'S1A_EW_GRDM_1SDH_*.zip'))
print(fl)

#date list
dl = []
for i in range(0,len(fl)):
    tmp = fl[i].split('/')[-1].split('.zip')[0].split('_')[4]
    date = datetime.strptime(tmp, "%Y%m%dT%H%M%S")
    dl.append(date)

print(dl)

for i in range(0,len(fl)):    

    dt1 = dl[i]
    if mode != 'all':
        #check if the file is early morning
        if dt1.hour > 8: continue

        #select dt2, depending on the mode
        dt_diff = np.array(dl) - dt1
        #ee: expected time difference: 5-6 h
        if mode == 'ee':
            match = np.argmin(abs(dt_diff - timedelta(hours=5)))
        #el: expected time difference: 7-9 h
        elif mode == 'el':
            match = np.argmin(abs(dt_diff - timedelta(hours=8)))
        #24h: expected time difference: 22-24 h 
        elif mode =='24h':
            match = np.argmin(abs(dt_diff - timedelta(hours=22)))
        
        print(match)
    
    else:
        match = i+1
    
    #find where Lance is
    mettime = getColumn(metfile,0)
    dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
    if dtb[0]>dt1: continue
    if dtb[-1]<dt1: continue
    mi = np.argmin(abs(np.asarray(dtb)-dt1))
    Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
    Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]
    if np.isnan(Lance_lon): continue

    # open files, read 'sigma0_HV' band and convert to UInt8 image
    f1 = fl[i]
    f2 = fl[match]
    print(f1)
    print(f2)
    if f1 == f2: print('same file'); continue
    
    print(datetime.strftime(dl[i], "%Y%m%dT%H%M%S"))
    print(datetime.strftime(dl[match], "%Y%m%dT%H%M%S"))

    #uncoment to check the pairs before processing
    #continue
    
    
    sid = SeaIceDrift(f1, f2)

    # apply Feature Tracking algorithm and retrieve ice drift speed
    # and starting/ending coordinates of matched keypoints
    uft, vft, lon1ft, lat1ft, lon2ft, lat2ft = sid.get_drift_FT()

    # user defined grid of points:
    lon1pm, lat1pm = np.meshgrid(np.linspace(regw, rege, gridp_we), np.linspace(regs, regn, gridp_sn))

    # apply Pattern Matching and find sea ice drift speed
    # for the given grid of points
    upm, vpm, apm, rpm, hpm, lon2pm, lat2pm = sid.get_drift_PM(lon1pm, lat1pm, lon1ft, lat1ft, lon2ft, lat2ft)

    #dump the PM data into numpy files
    name1 = fl[i].split('/')[-1].split('.zip')[0]
    np.save(outpath_drift+name1+'_upm',upm)
    np.save(outpath_drift+name1+'_vpm',vpm)
    np.save(outpath_drift+name1+'_apm',apm)
    np.save(outpath_drift+name1+'_rpm',rpm)
    np.save(outpath_drift+name1+'_hpm',hpm)
    np.save(outpath_drift+name1+'_lon1pm',lon1pm)
    np.save(outpath_drift+name1+'_lat1pm',lat1pm)
    np.save(outpath_drift+name1+'_lon2pm',lon2pm)
    np.save(outpath_drift+name1+'_lat2pm',lat2pm)

    # ==== PLOTTING ====
    date1 = datetime.strftime(dl[i], "%Y%m%dT%H%M%S")
    date2 = datetime.strftime(dl[match], "%Y%m%dT%H%M%S")
    # get coordinates of SAR scene borders
    lon1, lat1 = sid.n1.get_border()
    lon2, lat2 = sid.n2.get_border()

    # prepare projected images with sigma0_HV
    sid.n1.reproject(Domain(NSR().wkt, '-te -10 82 25 84 -ts 1000 1000'))
    s01 = sid.n1['sigma0_HV']
    sid.n2.reproject(Domain(NSR().wkt, '-te -10 82 25 84 -ts 1000 1000'))
    s02 = sid.n2['sigma0_HV']

    # plot the projected image from the first SAR scene
    plt.imshow(s01, extent=[regw, rege, regs, regn], cmap='gray', aspect=12)
    # plot vectors of sea ice drift from Feature Tracking
    plt.quiver(lon1ft, lat1ft, uft, vft, color='r',
            angles='xy', scale_units='xy', scale=0.5)
    # plot border of the second SAR scene
    plt.plot(lon2, lat2, '.-r')
    # set X/Y limits of figure
    plt.xlim([regw, rege])
    plt.ylim([regs, regn])
    plt.savefig(outpath+date1+'_'+date2+'_sea_ice_drift_FT_img1.png', dpi=500, bbox_inches='tight', pad_inches=0)
    plt.close('all')

    # plot the projected image from the second SAR scene
    plt.imshow(s02, extent=[regw, rege, regs, regn], cmap='gray', aspect=12)
    # filter only high quality pixels
    gpi = rpm > 0.4
    # plot vectors of sea ice drift from Feature Tracking, color by MCC
    plt.quiver(lon1pm[gpi], lat1pm[gpi], upm[gpi], vpm[gpi], rpm[gpi],
            angles='xy', scale_units='xy', scale=0.5)
    # plot border of the first SAR scene
    plt.plot(lon1, lat1, '.-r')
    # set X/Y limits of figure
    plt.xlim([regw, rege])
    plt.ylim([regs, regn])
    #plot Lance
    plt.plot(Lance_lon, Lance_lat, '*', color='purple', markersize=6)
    
    plt.savefig(outpath+date1+'_'+date2+'_sea_ice_drift_PM_img2.png', dpi=500, bbox_inches='tight', pad_inches=0)
    plt.close('all')

