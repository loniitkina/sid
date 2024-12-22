# modified from simple.py (part of SeaIceDrift by Anton Korosov)
import os
from os.path import exists
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
stp = 10; factor=0.5    #default run with 800m step (80m averaged pixel, sampled every 10 points)

ext = '_%s_%s_' %(stp,factor)   #for plot naming

#subsampling area around the ship (in degrees)
lon_diff = 15
lat_diff = 2

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
outpath_drift = '../../results/sid/drift/stp10_factor05/'
outpath_drift = '/scratch/pit000/results/sid/drift/stp10_factor05_200km/'
#outpath_drift = '/scratch/pit000/results/sid/drift/stp5_factor1/'
outpath = '/scratch/pit000/results/sid/plots/'

# ==== ICE DRIFT RETRIEVAL ====
#inpath = '/Data/pit000/ResearchData/IFT/EarthObservation/MOSAIC/SAR/Sentinel-1/'
inpath = '../../data/Sentinel-1/'  #make ln -s of all files from remote server above: ln -s /Data/pit000/ResearchData/IFT/EarthObservation/MOSAIC/SAR/Sentinel-1/* ../../data/Sentinel-1/
    #and remove all those duplicates by:
#rm ../../data/Sentinel-1/*\(1\).zip*

#show ship/CO position
#shipfile = '../sidrift/data/10minute_nounits.csv'
#shipfile = '../../downloads/position_leg3_nh-track.csv'
#shipfile = '../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0.csv'
    
#cover all tiles

#leg1
#ps_files=sorted(glob('../../downloads/data_master-solution_mosaic-leg1*_200km.csv'))
ps_files=sorted(glob('../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0_m*_200km.csv'))
ps_file_base = '../../downloads/data_master-solution_mosaic-leg1-20191016-20191213-floenavi-refstat-v1p0_%s_200km.csv'

#leg2
#ps_files=sorted(glob('../../downloads/data_master-solution_mosaic-leg2*_200km.csv'))
#ps_files=sorted(glob('../../downloads/data_master-solution_mosaic-leg2-20191214-20200224-floenavi-refstat-v1p0_m*_200km.csv'))
ps_file_base = '../../downloads/data_master-solution_mosaic-leg2-20191214-20200224-floenavi-refstat-v1p0_%s_200km.csv'

#leg3
#ps_files=sorted(glob('../../downloads/position_leg3_nh-track_[c,e,w,n,s,se,sw,nw,ne].csv')+glob('../../downloads/position_leg3_nh-track_[se,sw,nw,ne]?.csv'))
#ps_files=sorted(glob('../../downloads/position_leg3_nh-track_[c,w,n,s,se,sw,nw,ne]_200km.csv')+glob('../../downloads/position_leg3_nh-track_[se,sw,nw,ne]?_200km.csv'))
#ps_files=sorted(glob('../../downloads/position_leg3_nh-track_m*_200km.csv'))
ps_file_base = '../../downloads/position_leg3_nh-track_%s_200km.csv'

#N-ICE 2015
#ps_file_base = '../../downloads/lance_leg1_%s_200km.csv'

##CIRFA
#ps_file_base='../../downloads/CIRFA_cruise_stationM_%s_200km.csv'

regions=['c','msw','sw','s','mse','se','e','mne','ne','n','mnw','nw','w']
regions=['se','e','mne','ne','n','mnw','nw','w']    #finishing leg 3 run after a crach in mse
#regions=['sw','s','se','ne','n','nw','w']

ps_files=[]
for rr in regions:
    rfile = glob(ps_file_base %(rr))[0]
    #print(rfile)
    ps_files.append(rfile)

print(ps_files)

for shipfile in ps_files:
    print('coordinates from file: ',shipfile)
    region=shipfile.split('.csv')[0].split('_')[-2]
    print('region: ',region)
    
    #read filelist for SAR images from track filename list
    trackfile=shipfile.split('.csv')[0]+'-fnames.csv'
    print(trackfile)
    fl = getColumn(trackfile,1,header=False)
    fl = [ inpath+i.split('.SAFE')[0]+'.zip' for i in fl ]
    print(fl)

    #date list
    dl = []
    for i in range(0,len(fl)):
        tmp = fl[i].split('/')[-1].split('.zip')[0]
        if tmp == 'missing file':
            dl.append(tmp)
        else:
            tmp=tmp.split('_')[4]
            date = datetime.strptime(tmp, "%Y%m%dT%H%M%S")
            dl.append(date)
    #print(dl)

    for i in range(0,len(fl)-1):    

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
        
        #check what files we have        
        f1 = fl[i]
        f2 = fl[match]
        print('We have files:')
        print(f1)
        print(f2)
        
        #sometimes there is no scene for that day
        #then f[i] or f[match] == 'missing file'
        j=i
        while f1==inpath+'missing file.zip':
            j=j+1
            f1 = fl[j]
            try:
                f2 = fl[j+1]
            except:
                #at some point there will be no more files
                print('no more f1 files in region: ', region)
                f1='dummy'
                continue
            print(f1)
            print(f2)

        k=j
        while f2==inpath+'missing file.zip':
            k=k+1
            try:
                f2 = fl[k+1]
            except:
                print('no more f2 files in region: ', region)
                f2='dummy'
                continue
            print(f2)
        
        #go to next region
        if (f1=='dummy') or (f2=='dummy'):
            print('moving to next region: ',f1,f2); continue

        #check that files really exist
        file_exists = exists(f1)
        if file_exists:
            print(f1)
        #sometimes there is a small difference in last part of the name
        else:
            #../../data/S1B_EW_GRDM_1SDH_20191019T053003_20191019T053107_018540_022EED_E941.zip
            parts=f1.split('_')
            first_part=parts[0]+'_'+parts[1]+'_'+parts[2]+'_'+parts[3]+'_'+parts[4]+'_'+parts[5]+'_'
            try:
                f1=glob(first_part+'*')[0]
            except:
                #or just take another scene from same date
                #it can be both S1B or S1A - but it needs to be EW
                alternative = inpath+'S1*_EW_GRDM_1SDH_'+parts[4].split('T')[0]+'*.zip'
                f1=glob(alternative)[0]
            print(f1)
        
        #same for the second file
        file_exists = exists(f2)
        if file_exists:
            print(f2)
        #sometimes there is a small difference in last part of the name
        else:
            #../../data/S1B_EW_GRDM_1SDH_20191019T053003_20191019T053107_018540_022EED_E941.zip
            parts=f2.split('_')
            first_part=parts[0]+'_'+parts[1]+'_'+parts[2]+'_'+parts[3]+'_'+parts[4]+'_'+parts[5]+'_'
            try:
                f2=glob(first_part+'*')[0]
            except:
                #or just take another scene from same date
                #it can be both S1B or S1A - but it needs to be EW
                alternative = inpath+'S1*_EW_GRDM_1SDH_'+parts[4].split('T')[0]+'*.zip'
                print(alternative)
                f2=glob(alternative)[0]
            print(f2)
        
        print('We are working with files: ', region)
        print(f1)
        print(f2)
        
        #WARNING: the file has to be EW, IW wont work with Nansat: ValueError: Cannot find band {'name': 'sigma0_HV'}! band_number is from 1 to 13
        
        if f1 == f2: print('same file'); continue
        
        date1 = f1.split('/')[-1].split('.zip')[0].split('_')[4]       
        date2 = f2.split('/')[-1].split('.zip')[0].split('_')[4]
        dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
        dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
        print(date1,date2)
            
        #find where ship/CO is
        mettime = getColumn(shipfile,0)
        dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
        if dtb[0]>dt1: continue
        if dtb[-1]<dt1: continue
        mi = np.argmin(abs(np.asarray(dtb)-dt1))
        ship_lon = np.asarray(getColumn(shipfile,1),dtype=float)[mi]
        ship_lat = np.asarray(getColumn(shipfile,2),dtype=float)[mi]
        if np.isnan(ship_lon): continue

        print('Ship at: ',ship_lon,ship_lat)

        #get time difference
        timediff = (dt2 - dt1).total_seconds()
        print('Time difference: '+str(timediff)+' seconds')

        #get sea ice drift
        #open files, read 'sigma0_HV' band and convert to UInt8 image
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
        
        #container = np.load(out_file)
        #print(container.files)
        #print(container['upm'])


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

        #import ipdb; ipdb.set_trace()
        #store the intensities and coordinates for a mosaic plot
        #15,19,21 Jan are good dates
        #dump the PM data into numpy file
        lon2_all, lat2_all = sid.n2.get_geolocation_grids()     #get all coordinates for the second image
        out_file = outpath_drift+'S1_'+date1+'_'+date2+'_'+region+ext+'.npz'
        np.savez(out_file,sigma = s02,lon2=lon2_all,lat2=lat2_all)  #17MB = reasonable size
        print(out_file, 'stored.')
        
        
