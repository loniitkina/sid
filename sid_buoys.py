#! /usr/bin/python
# -*- coding: utf-8 -*-

#triangulate the buoy cordinates and calculate the deformation

#use WGS84 ellipsoid for the distances
#https://pypi.python.org/pypi/LatLon

from datetime import date, timedelta, datetime
import csv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from glob import glob
import os
from pyproj import Proj
from scipy.spatial import Delaunay
import pandas as pd

from sid_func import *

##get all the buoy files
#################################################Pick only one block, comment out the others
dep='1'
path = '../../buoys/data/deployment1/'
path_out = '../data/buoys/'
##1st SAR period
start = datetime(2015, 1, 21, 0, 0)
end = datetime(2015, 2, 5, 0, 0)
outname_sh = 'nice1_in_SAR1_'
###2nd SAR period
#start = datetime(2015, 2, 5, 0, 0)
#end = datetime(2015, 2, 16, 0, 0)
#outname_sh = 'nice1_in_SAR2_'
##pick which buoys to use: 
##all buoys
flist= glob(path+'*rev.csv')
#inner ring
flist=glob(path+'CALIB*rev.csv')+glob(path+'SVP*rev.csv')+glob(path+'SIMBA_2015[c-f]*rev.csv')+glob(path+'SNOW_2015a_rev.csv')
##outer ring
#flist=glob(path+'IMB*rev.csv')+glob(path+'WAVE*rev.csv')+glob(path+'SNOW_2015a_rev.csv')+glob(path+'SIMBA[a-b]*rev.csv')
##selection of buoys
#flist= glob(path+'IMB*rev.csv')+glob(path+'WAVE*rev.csv')+glob(path+'SNOW*_rev.csv')+glob(path+'SIMBA_2015*rev.csv')
####################################################end

##1h -- dont use this with all buoys - use just the ones with valid hourly output
##leave out CALIBs and bad SVPs with 3h positions
#flist= glob(path+'SNOW*rev.csv')+glob(path+'SIMBA*rev.csv')+glob(path+'WAVE*rev.csv')+glob(path+'IMB*rev.csv')+glob(path+'SVP_2015[a-b]*rev.csv')
##inner ring
#flist=glob(path+'SVP*[a-b]*rev.csv')+glob(path+'SIMBA_2015[c-f]*rev.csv')+glob(path+'SNOW_2015a_rev.csv')
#idx1=2;idx2=1

##3h
#idx1=2;idx2=3

##6h
#idx1=2;idx2=6

##12h
#idx1=2;idx2=12

#1 day
#pick an entry every day, start with 2, time step 24h
idx1=2;idx2=24

##5 days
##pick an entry every 5 days, start with 43, time step 120h
#idx1=43;idx2=120

outname2 = outname_sh+str(idx2)+'h'
print outname2
####################################################end

#min angle of triangle that will be used for calculation
minang_val = 15
####################################################end

#empty lists to collect the values for power law scatter plots
dr_list_l=[]
ls_list_l=[]

#for idx1 in range(2,3):
#sliding picking of the corner points (to increase the sample number)
for idx1 in range(1,idx2+1):
    
    print '###############################################################################SLIDER: '+str(idx1)
    
    #get the dates
    time = getColumn(flist[0],0)
    dt = [ datetime.strptime(time[x], "%Y-%m-%d %H:%M:%S") for x in range(len(time)) ][idx1::idx2]
    #get just the predefined date/time interval
    si = np.argmin(abs(np.asarray(dt)-start))
    ei = np.argmin(abs(np.asarray(dt)-end))
    dates = dt[si:ei]

    #make an empty array for all the buoy data
    data = np.zeros((len(dates),4,len(flist)))

    #open all files and store all the buoy data in one array
    i=0
    for buoy in flist:
        #print buoy

        #remove duplicated positions
        #if dep=='1':
            #if buoy == '../../buoys/data/deployment1/SIMBA_2015i_rev.csv':print 'match';continue	#french buoy with 2h positions
            #if buoy == '../../buoys/data/deployment1/SIMBA_2015h_edge_dist_rev.csv':print 'match';continue	#french buoy with 2h positions
            #if buoy == '../../buoys/data/deployment1/SIMBA_2015d_rev.csv':print 'match';continue	#deployed on floe2 in March
            #if buoy == '../../buoys/data/deployment1/SIMBA_2015b_rev.csv':print 'match';continue	#deployed together with IMB_2015a, almost no GPS data
            #if buoy == '../../buoys/data/deployment1/SNOW_2015d_rev.csv':print 'match';continue		#deployed on floe2 in March
            #if buoy == '../../buoys/data/deployment1/WAVE_2015a_rev.csv':print 'match';continue		#deployed together with SIMBA_2015a
            #if buoy == '../../buoys/data/deployment1/SIMBA_2015f_rev.csv':print 'match';continue	#SIMBA deployed next to SNOW_2015a, almost no GPS data
            
            #outer ring
            #if buoy == '../../buoys/data/deployment1/WAVE_2015a_rev.csv':print 'match';continue
            #if buoy == '../../buoys/data/deployment1/WAVE_2015b_rev.csv':print 'match';continue
            #if buoy == '../../buoys/data/deployment1/WAVE_2015c_rev.csv':print 'match';continue
            #if buoy == '../../buoys/data/deployment1/WAVE_2015d_rev.csv':print 'match';continue
            
        x = np.asarray(getColumn(buoy,4),dtype=float)[idx1::idx2]
        y = np.asarray(getColumn(buoy,3),dtype=float)[idx1::idx2]
        
        #calculate velocity from positions
        x = np.ma.array(x,mask=x==0)
        xsh = x[1:]
        up = (xsh-x[:-1])/float(3600*idx2)
        
        y = np.ma.array(y,mask=y==0)
        ysh = y[1:]
        vp = (ysh-y[:-1])/float(3600*idx2)
        
        ##th_mask = ((up[si:ei]<0.02) & (up[si:ei]>-0.02)) | (vp[si:ei]<0.02) & (vp[si:ei]>-0.02)	#exclude all velocities under error treshold
        th_mask = (up[si:ei]==0)| (vp[si:ei]==0)
        up=np.ma.array(up[si:ei]/1000*3600*24,mask=th_mask)
        vp=np.ma.array(vp[si:ei]/1000*3600*24,mask=th_mask)
        x=np.ma.array(x[si:ei]/1000,mask=x[si:ei]==0)
        y=np.ma.array(y[si:ei]/1000,mask=y[si:ei]==0)
        
        #no interpolations
        data[:,0,i] =x
        data[:,1,i] =y
        data[:,2,i] =up
        data[:,3,i] =vp
        
        i=i+1

        del x,y,up,vp

    data = np.where(np.isnan(data),0,data)
    data = np.ma.array(data,mask=data==0)

    count = 0
    tottri = 0

    #for all time steps
    for k in range(1,data.shape[0]):
        #print '******************************************************************************************k:',k
        #triangulate all the positions
        #get the x and y coordinates, skip the masked values and fit them into pairs
        pts = data[k,:2,:].T.compressed().reshape(-1,2)
        #latlon = data[k,4:,:].T.compressed().reshape(-1,2)
        
        #skip the time steps with less than 3 points
        #print 'number of points:', pts.shape[0]
        if pts.shape[0] < 3: continue
        
        #velocities
        u = data[k,2,:].compressed()
        v = data[k,3,:].compressed()
        
        if u.shape!=data[k,0,:].compressed().shape:
            #print 'this can be a problem!'
            continue

        #get all possible combinations of 3 points - just keep unique combinations
        from itertools import combinations,chain
        
        tri = []
        tri.extend(combinations(pts, 3))
        trin = len(tri)
        
        utri = []
        utri.extend(combinations(u, 3))
        
        vtri = []
        vtri.extend(combinations(v, 3))
        
        dux = np.zeros(trin)
        duy = np.zeros(trin)
        dvx = np.zeros(trin)
        dvy = np.zeros(trin)
        tris = np.zeros(trin)
        trid = np.zeros(trin)
        
        lsa = np.zeros(trin)
        
        t=-1
        ##for each triangle
        for vert in combinations(pts, 3):
            t = t + 1
            #transform tuple of arrays into an array
            vert = np.asarray(vert)

            #sorting the vertices so that they are always counter-clockwise
            from scipy.spatial import ConvexHull
            try:
                hull = ConvexHull(vert)
            except:
                print 'bad triangle!!!'; continue
            vert = vert[hull.vertices]
        
            #skip all triangles with the too short or too long sides
            
            #% sides of triangle
            ta = np.sqrt((vert[1,0]-vert[0,0])**2+(vert[1,1]-vert[0,1])**2)
            tb = np.sqrt((vert[2,0]-vert[1,0])**2+(vert[2,1]-vert[1,1])**2)
            tc = np.sqrt((vert[2,0]-vert[0,0])**2+(vert[2,1]-vert[0,1])**2)
            #if ta < 20 or tb < 20 or tc < 20: continue		#shorter than 100m     
            #if ta > 200 or tb > 200 or tc > 200: continue	#longer than 500km
            
            #find smallest angle to filter out too acute triangles
            sides = [ta,tb,tc]
            srt = min(sides)
            srtid = np.argmin(sides)
            lng = np.delete(sides,srtid)
            #minang(ii) = acosd((sides(1)^2+sides(2)^2-shortest^2)/2/sides(1)/sides(2)); % law of cosines
            minang = np.degrees(np.arccos((lng[0]**2+lng[1]**2-srt**2)/2/lng[0]/lng[1]))
            #print minang; exit()
            if minang < minang_val: continue
            
            ##area of triangle
            ##A = |0.5 * (x1(y2-y3)+x2(y3-y1)+x3(y1-y2))|
            #the absolute value not necessary: vertices should be oriented counter-clockwise, the formula only gives neg. values if this is not so
            #Hutchings etal, 2012 (ISPOL)
            area = .5* (vert[0,0]*vert[1,1] - vert[0,1]*vert[1,0] + vert[1,0]*vert[2,1] - vert[1,1]*vert[2,0] + vert[2,0]*vert[0,1] - vert[2,1]*vert[0,0])
            if area < 0: print 'not counter-clockwise oriented triangle!'; exit()
            
            #some filtering to throw out the erroneous GPS positioning. If the accuracy is 2m.
            if area < 1: print area; continue	#Hutchings etal, 2011 (SEDNA)
        
            ##also filter out all the data where the area changes too fast (Jennifer Hutchings)
            ##we cant do this since we dont track the individual triangles, but we can filter out too high and too low velocities instead
            ##change 1km**2/10min == 12km/day
            #this unevenly filters out high deformation at large triangles
            
            #print t
            #print vert
            #print len(utri)
            #print utri
            uvert = np.array(utri[t])
            try: 
                vvert = np.array(vtri[t])
            except:
                print uvert
                print vvert
                continue
            uvert = np.array(uvert[hull.vertices])
            vvert = np.array(vvert[hull.vertices])
            
            maxts = np.max(np.sqrt(uvert**2+vvert**2))	#unit is km/day
            mints = np.min(np.sqrt(uvert**2+vvert**2))
            tsdiff = np.abs(maxts-mints)/np.sqrt(area)
            #if tsdiff > 6: print '%s%s' %('too high speed difference: ',tsdiff) ;continue
        
        

            #print '%s%s' %('calculating triangle',t)
            dux[t] = (.5 / area) * (
            (uvert[1]+uvert[0])*(vert[1,1]-vert[0,1])+
            (uvert[2]+uvert[1])*(vert[2,1]-vert[1,1])+
            (uvert[0]+uvert[2])*(vert[0,1]-vert[2,1]) )
            
            duy[t] = -(.5 / area) * (
            (uvert[1]+uvert[0])*(vert[1,0]-vert[0,0])+
            (uvert[2]+uvert[1])*(vert[2,0]-vert[1,0])+
            (uvert[0]+uvert[2])*(vert[0,0]-vert[2,0]) )
            
            dvx[t] = (.5 / area) * (
            (vvert[1]+vvert[0])*(vert[1,1]-vert[0,1])+
            (vvert[2]+vvert[1])*(vert[2,1]-vert[1,1])+
            (vvert[0]+vvert[2])*(vert[0,1]-vert[2,1]) )
            
            dvy[t] = -(.5 / area) * (
            (vvert[1]+vvert[0])*(vert[1,0]-vert[0,0])+
            (vvert[2]+vvert[1])*(vert[2,0]-vert[1,0])+
            (vvert[0]+vvert[2])*(vert[0,0]-vert[2,0]) )

            tris[t] = (np.sqrt(uvert[0]**2+vvert[0]**2)+
                    np.sqrt(uvert[1]**2+vvert[1]**2)+
                    np.sqrt(uvert[2]**2+vvert[2]**2) )
            
            dx = np.mean(uvert)
            dy = np.mean(vvert)
            
            trid[t] = np.degrees(np.arctan2(dy,dx))
            #arctan2 gives the angle from origin at [1,0] - that is 90 deg. clockwise of the N direction
            #the wind angle is given for the direction the wind is coming from: add 180 to be comparable
            trid[t] = 90 - trid[t] + 180
            
            #skip very high speeds (helicopter deployment/errors):
            #if tris[t] > 80: dux[t]=0 ; continue
            
            #plot deformation rate vs. length scale
            #divergence
            a = (dux[t] + dvy[t])
            #shear
            b = .5*np.sqrt((dux[t]-dvy[t])**2+(duy[t]+dvx[t])**2)
            #deformation rate
            dr = np.sqrt(a**2+b**2)
            ls = np.sqrt(area)
            lsa[t] = ls
            
            dr_list_l.append(dr)
            ls_list_l.append(ls)
            
            tottri = tottri+1

        del dux, duy, dvx, dvy, tris, pts, u, v, lsa
        count = count + 1

    print 'total number of all valid triangles t:',tottri, 'in time steps n:',count, 'over hours h:', k


#save all the data for the power law
ls_list_l = np.array(ls_list_l)
dr_list_l = np.array(dr_list_l)
ls_list_l.dump(path_out+'ls_'+outname2)
dr_list_l.dump(path_out+'dr_'+outname2)

print '%s%s%s' %('This calculation is based on positons of ',len(flist),' buoys.')
