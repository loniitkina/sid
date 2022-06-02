# modified from simple.py (part of SeaIceDrift by Anton Korosov)
import os
import sys
from glob import glob
import unittest
import inspect

import csv
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('Agg')

from nansat import Nansat, Domain, NSR

from sea_ice_drift import SeaIceDrift


def getColumn(filename, column, delimiter=',', header=True):
    results = csv.reader(open(filename),delimiter=delimiter)
    if header==True:
        next(results, None)
    return [result[column] for result in results]

# download Sentinel-1 data with sentinelsat
#make json point files with qgis (or write with a script)
#in QGIS: create layer - new temporary scratch layer - draw point - save as
#create user on scihub: https://scihub.copernicus.eu/
#download data to the data partition!
#create custom geojson by http://geojson.io/#map=2/45.0/4.2
#sentinelsat --user polona --password jank0vna --geometry /mnt/data/sea_ice_drift/lance.geojson --sentinel 1 --producttype GRD --start 20150128 --end 20150202 --download --path /mnt/data/sea_ice_drift/Sentinel1
#sentinelsat --user polona --password jank0vna --geometry nice.geojson --sentinel 1 --producttype GRD --start 20150128 --end 20150202 --download --path ../sidrift/data/Sentinel1_new

#retrieve drift
#plot drift and write text output
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#selection mode (in time)
#all: all images are analyzed
#ee: only image from early morning and early afternoon (same date) - ee and el will often be mixed up (lots of missing images etc) - that is not a problem!
#el: only image from early morning and late afteroon (same date)
#24h: only early mornming images (2 different dates)
#all drift files will be in same folder files with same time steps will be overwritten - finally there will be no double files!
mode = ['all', 'ee', 'el', '24h', '2d', '3d'][0]
print(mode)

# full resolution (pixel spacing) is 80m (or even 40m - something to check)
# template for PM is set to 35x35 pixels, so at 10 step resolution, overlap is substational - how high is it sensible to go?
#stp=1: only 24h mode needs such high spatial resolution.
#stp=10: sufficient for the time scalling law
#stp = 5; factor=1;        #200m step
stp = 10; factor=0.5    #default run with 800m step (80m averaged pixel, sampled every 100 points)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
outpath_drift = '../../results/sid/drift/stp10_factor05/'
#outpath_drift = '../../results/sid/drift/stp5_factor1/'
outpath = '../../results/sid/plots/'

# ==== ICE DRIFT RETRIEVAL ====
#inpath = '/Data/pit000/ResearchData/IFT/EarthObservation/MOSAIC/SAR/Sentinel-1/'
inpath = '../../data/'  #make ln -s of all files from remote server above: ln -s /Data/pit000/ResearchData/IFT/EarthObservation/MOSAIC/SAR/Sentinel-1/* ../../data/

#show ship/CO position
#shipfile = '../sidrift/data/10minute_nounits.csv'
shipfile = '../../downloads/position_leg3_nh-track.csv'
#shipfile = '../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0.csv'

#cover all quadrants
ps_files=sorted(glob('../../downloads/position_leg3_nh-track_[e,w,n,s,se,sw,nw,ne].csv')+glob('../../downloads/position_leg3_nh-track_[se,sw,nw,ne]?.csv'))
print(ps_files)

for shipfile in ps_files:
    print(shipfile)
    region=shipfile.split('.csv')[0].split('_')[-1]
    print(region)
    
    #read filelist for SAR images from track filename list
    trackfile=shipfile.split('.csv')[0]+'-fnames.csv'
    print(trackfile)
    fl = getColumn(trackfile,1)
    fl = [ inpath+i.split('.SAFE')[0]+'.zip' for i in fl ]
    print(fl)

    ##file list
    #fl = sorted(glob(inpath+'S1*.zip'))
    ##print(fl)
    ##fl = fl[:-80]   #leg 1
    ##fl = fl[56:-80] #leg 1 after processing error (file 47 was last)
    ##fl = fl[-80:-55]    #process all March/April event period
    ##fl = fl[-55:]	#process all the way
    #print(fl)
    ##exit()

    #date list
    dl = []
    for i in range(0,len(fl)):
        tmp = fl[i].split('/')[-1].split('.zip')[0].split('_')[4]
        date = datetime.strptime(tmp, "%Y%m%dT%H%M%S")
        dl.append(date)
    #print(dl)

    for i in range(0,len(fl)):    

        dt1 = dl[i]
        if mode != 'all':
            #check if the file is early morning
            if dt1.hour > 7: print('not early morning');continue

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
                if (dl[match].hour > 7): print('second file not early morning');continue
                if (dl[match].hour == 7) & (dl[match].minute > 58): print('second file not early morning');continue
            #2days
            elif mode =='2d':
                match = np.argmin(abs(dt_diff - timedelta(hours=46)))
                if (dl[match].hour > 7): print('second file not early morning');continue
                if (dl[match].hour == 7) & (dl[match].minute > 58): print('second file not early morning');continue
            #3days
            elif mode =='3d':
                match = np.argmin(abs(dt_diff - timedelta(hours=70)))
                if (dl[match].hour > 7): print('second file not early morning');continue
                if (dl[match].hour == 7) & (dl[match].minute > 58): print('second file not early morning');continue
        else:
            match = i+1
        
        #find where ship/CO is
        mettime = getColumn(shipfile,0)
        dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
        if dtb[0]>dt1: continue
        if dtb[-1]<dt1: continue
        mi = np.argmin(abs(np.asarray(dtb)-dt1))
        ship_lon = np.asarray(getColumn(shipfile,1),dtype=float)[mi]
        ship_lat = np.asarray(getColumn(shipfile,2),dtype=float)[mi]
        if np.isnan(ship_lon): continue

        print('ship at: ',ship_lon,ship_lat)
        #exit()

        # open files, read 'sigma0_HV' band and convert to UInt8 image
        f1 = fl[i]
        f2 = fl[match]
        print(f1)
        print(f2)
        if f1 == f2: print('same file'); continue
        
        date1 = datetime.strftime(dl[i], "%Y%m%dT%H%M%S")
        date2 = datetime.strftime(dl[match], "%Y%m%dT%H%M%S")
        print(date1,date2)
        
        #get time difference
        timediff = (dl[match] - dl[i]).total_seconds()
        print('Time difference: '+str(timediff)+' seconds')

        #get sea ice drift
        #up- or down-sample the images, default factor=0.5 downsamples from 40 to 80m pixels
        sid = SeaIceDrift(f1, f2)

        # apply Feature Tracking algorithm and retrieve ice drift speed
        # and starting/ending coordinates of matched keypoints
        print('*****************************Feature Tracking*****************************')
        uft, vft, lon1ft, lat1ft, lon2ft, lat2ft = sid.get_drift_FT()

        #this could be substituted from final positions in previous image pair (in time)
        #IF the location was precise!!!


        # user defined grid of points:
        srs = '+proj=laea lat_0=%f lon_0=%f +datum=WGS84 +ellps=WGS84 +no_defs' % (90, 10)
        #d = Domain(srs, '-te -100000 -100000 100000 100000 -tr 1000 1000')
        #lon1pm, lat1pm = d.get_geolocation_grids()
        lon1pm, lat1pm = sid.n1.get_geolocation_grids()

        # subsample lon,lat with even steps
        lon1pm, lat1pm = lon1pm[::stp, ::stp], lat1pm[::stp, ::stp] 
        
        # subsample around ship
        lonlat_shape = lon1pm.shape
        lon_diff = 15
        lat_diff = 2
        near_lance_pix = ((lon1pm > (ship_lon - lon_diff )) *
                        (lon1pm < (ship_lon + lon_diff )) * 
                        (lat1pm > (ship_lat - lat_diff )) * 
                        (lat1pm < (ship_lat + lat_diff )))
        # that will return lon,lat as vectors, not as grids   
        lon1pm, lat1pm  = lon1pm[near_lance_pix] , lat1pm[near_lance_pix] 
        
        # apply Pattern Matching and find sea ice drift speed
        # for the given grid of points
        #if no spatial overlap this will give an error: UnboundLocalError: local variable 'lon_pm2_grd' referenced before assignment
        print('*****************************Pattern Matching*****************************')
        try:
            #try other template sizes/img_size of 1st image. Values > 20 should be fine. Margin = search distance on 2nd image.
            upm, vpm, apm, rpm, hpm, lon2pm, lat2pm = sid.get_drift_PM(lon1pm, lat1pm, lon1ft, lat1ft, lon2ft, lat2ft, img_size=45, margin=0, srs=srs)
        except:
            continue

        #calculate velocity from displacements
        upm = upm/timediff
        vpm = vpm/timediff

        # create grids from vectors
        upm2 = np.zeros(lonlat_shape) + np.nan
        vpm2 = np.zeros(lonlat_shape) + np.nan
        rpm2 = np.zeros(lonlat_shape) + np.nan
        apm2 = np.zeros(lonlat_shape) + np.nan
        hpm2 = np.zeros(lonlat_shape) + np.nan
        lon2pm2 = np.zeros(lonlat_shape) + np.nan
        lat2pm2 = np.zeros(lonlat_shape) + np.nan
        lon1pm2 = np.zeros(lonlat_shape) + np.nan
        lat1pm2 = np.zeros(lonlat_shape) + np.nan
        
        
        upm2[near_lance_pix] = upm
        vpm2[near_lance_pix] = vpm
        apm2[near_lance_pix] = apm
        rpm2[near_lance_pix] = rpm
        hpm2[near_lance_pix] = hpm
        #lon2pm2[near_lance_pix] = lon2pm
        #lat2pm2[near_lance_pix] = lat2pm
        lon1pm2[near_lance_pix] = lon1pm
        lat1pm2[near_lance_pix] = lat1pm
        
        #dump the PM data into numpy file
        out_file = outpath_drift+'SeaIceDrift_'+date1+'_'+date2+'_'+region+'.npz'
        #np.savez(out_file,upm = upm2,vpm = vpm2, apm = apm2, rpm = rpm2, hpm = hpm2, lon1 = lon1pm2, lat1 = lat1pm2, lon2 = lon2pm2, lat2 = lat2pm2)
        np.savez(out_file,upm = upm2,vpm = vpm2, apm = apm2, rpm = rpm2, hpm = hpm2, lon1 = lon1pm2, lat1 = lat1pm2)
        
        container = np.load(out_file)
        print(container.files)
        print(container['upm'])


        ## ==== PLOTTING ====
        ## get coordinates of SAR scene borders
        lon1, lat1 = sid.n1.get_border()
        lon2, lat2 = sid.n2.get_border()
        
        #specify region
        regn = ship_lat+lat_diff; regs = ship_lat-lat_diff
        regw = ship_lon-lon_diff; rege = ship_lon+lon_diff

        ## prepare projected images with sigma0_HV
        #sid.n1.reproject(Domain(NSR().wkt, '-te -10 82 25 84 -ts 1000 1000'))
        #s01 = sid.n1['sigma0_HV']
        region_text = '-te %f %f %f %f -ts 1000 1000' % (regw,regs,rege,regn)
        print(region_text)
        sid.n2.reproject(Domain(NSR().wkt, region_text ))
        s02 = sid.n2['sigma0_HV']

        ## plot the projected image from the first SAR scene
        #plt.imshow(s01, extent=[regw, rege, regs, regn], cmap='gray', aspect=12)
        ## plot vectors of sea ice drift from Feature Tracking
        #plt.quiver(lon1ft, lat1ft, uft, vft, color='r',
                #angles='xy', scale_units='xy', scale=0.5)
        ## plot border of the second SAR scene
        #plt.plot(lon2, lat2, '.-r')
        ## set X/Y limits of figure
        #plt.xlim([regw, rege])
        #plt.ylim([regs, regn])
        #plt.savefig(outpath+'SeaIceDrift_FT_img1_'+date1+'_'+date2+'.png', dpi=500, bbox_inches='tight', pad_inches=0)
        #plt.close('all')

        # plot the projected image from the second SAR scene
        plt.imshow(s02, extent=[regw, rege, regs, regn], cmap='gray', aspect=12)
        # filter only high quality pixels
        gpi = rpm > 0.4
        # plot vectors of sea ice drift from Feature Tracking, color by MCC
        plt.quiver(lon1pm[gpi], lat1pm[gpi], upm[gpi], vpm[gpi], hpm[gpi],
                angles='xy', scale_units='xy', scale=1)
        # plot border of the first SAR scene
        plt.plot(lon1, lat1, '.-r')
        # set X/Y limits of figure
        plt.xlim([regw, rege])
        plt.ylim([regs, regn])
        #plot ship
        plt.plot(ship_lon, ship_lat, '*', color='purple', markersize=6)
        
        plt.savefig(outpath+'SeaIceDrift_PM_img2_'+date1+'_'+date2+'_'+region+'.png', dpi=500, bbox_inches='tight', pad_inches=0)
        plt.close('all')
