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
end = datetime(2015, 1, 27, 0, 0)
outname2 = 'nice1_comb_SAR1_1h'
##2nd SAR period
start = datetime(2015, 2, 3, 0, 0)
end = datetime(2015, 2, 16, 0, 0)
outname2 = 'nice1_comb_SAR2_1h'
#pick which buouys to use: 
##all buoys
flist= glob(path+'*rev.csv')
##inner ring
#flist=glob(path+'CALIB*rev.csv')+glob(path+'SVP*rev.csv')+glob(path+'SIMBA_2015[c-f]*rev.csv')+glob(path+'SNOW_2015a_rev.csv')+glob(path+'IMB*rev.csv')+glob(path+'SIMBA_2015a_rev.csv')
#print len(flist)
####results for time series (leave out CALIBs and bad SVPs for the time series)
flist= glob(path+'SNOW*rev.csv')+glob(path+'SIMBA*rev.csv')+glob(path+'WAVE*rev.csv')+glob(path+'IMB*rev.csv')+glob(path+'SVP_2015[a-b]*rev.csv')
scales = [0,100]
####################################################end

##pick the scales (borders of the bins over which averages are made) for the power law scatter plot
##the power law is very scale sensitive! Results can be very different depending on what bin margins are chosen

#1h -- dont use this with all buoys - use juse the ones with valid hourly output
idx1=2;idx2=1

##3h
#idx1=2;idx2=3

##6h
#idx1=2;idx2=6

##12h
#idx1=2;idx2=12

##1 day
##pick an entry every day, start with 2, time step 24h
#idx1=2;idx2=24

##5 days
##pick an entry every 5 days, start with 43, time step 120h
#idx1=43;idx2=120

#min angle of triangle that will be used for calculation
minang_val = 15
####################################################end

print flist

#get the dates
time = getColumn(flist[0],0)
dt = [ datetime.strptime(time[x], "%Y-%m-%d %H:%M:%S") for x in range(len(time)) ][idx1::idx2]
#get just the predefined date/time interval
si = np.argmin(abs(np.asarray(dt)-start))
ei = np.argmin(abs(np.asarray(dt)-end))
dates = dt[si:ei]
print dates

#write out the array
dta = np.array(dates)
dta.dump('dates'+outname2)

#meteorological data (from Lance's met system)
metfile = '../data/10minute_nounits.csv'
mettime = getColumn(metfile,0)
metdates = [ datetime.strptime(mettime[x], "%Y-%m-%d %H:%M:%S") for x in range(len(mettime)) ]
metspeed = np.asarray(getColumn(metfile,10),dtype=float)
metspeed = np.ma.masked_invalid(metspeed)
metdir = np.asarray(getColumn(metfile,11),dtype=float)
metdir = np.ma.masked_invalid(metdir)
mettemp = np.asarray(getColumn(metfile,6),dtype=float)
mettemp = np.ma.masked_invalid(mettemp)

#average the wind speed for every hour
col = len(metspeed)/6
tmp = metspeed.reshape(col,6)
ttmp = np.mean(tmp,axis=1)

tmp = mettemp.reshape(col,6)
atmp = np.mean(tmp,axis=1)

dmp = metdir.reshape(col,6)

def circmean(alpha,axis=None):
  #To convert from radians to degrees, multiply by (180o/(PI))
  tod = 180/np.pi
  tor = np.pi/180

  sa = np.mean(np.sin(alpha*tor),axis)
  ca  = np.mean(np.cos(alpha*tor),axis)
  mean_angle = np.arctan2(sa,ca)*tod
  mean_angle = np.where(mean_angle<0,mean_angle+360,mean_angle)
  return mean_angle

tdmp = circmean(dmp,axis=1)

#get the hourly dates
dtmp = metdates[48*6::120*6]
wsi = np.argmin(abs(np.asarray(dtmp)-start))
wei = np.argmin(abs(np.asarray(dtmp)-end))

wspeed = ttmp[wsi:wei]
wdir = tdmp[wsi:wei]
wtemp = atmp[wsi:wei]

#plt.plot(metdates,metdir)
#plt.plot(dates,wdir)
#plt.xlim(start,end)
#plt.show()

#always dump the data on disc, so that the time series can be done separately with another script
wspeed.dump('wspeed'+outname2)
wdir.dump('wdir'+outname2)
wtemp.dump('wtemp'+outname2)

#make an empty array for all the buoy data
data = np.zeros((len(dates),4,len(flist)))

#open all files and store all the buoy data in one array
i=0
for buoy in flist:
  print buoy

  #remove duplicated positions
  if dep=='1':
    if buoy == '../data/deployment1/SIMBA_2015i_rev.csv':print 'match';continue	#french buoy with 2h positions
    if buoy == '../data/deployment1/SIMBA_2015h_rev.csv':print 'match';continue	#french buoy with 2h positions
    #if buoy == '../data/deployment1/SIMBA_2015d_rev.csv':print 'match';continue	#deployed on floe2 in March
    #if buoy == '../data/deployment1/SIMBA_2015b_rev.csv':print 'match';continue	#deployed together with IMB_2015a, almost no GPS data
    #if buoy == '../data/deployment1/SNOW_2015d_rev.csv':print 'match';continue		#deployed on floe2 in March
    #if buoy == '../data/deployment1/WAVE_2015a_rev.csv':print 'match';continue		#deployed together with SIMBA_2015a
    #if buoy == '../data/deployment1/SIMBA_2015f_rev.csv':print 'match';continue	#SIMBA deployed next to SNOW_2015a, almost no GPS data
    
    #outer ring
    #if buoy == '../data/deployment1/WAVE_2015a_rev.csv':print 'match';continue
    #if buoy == '../data/deployment1/WAVE_2015b_rev.csv':print 'match';continue
    if buoy == '../data/deployment1/WAVE_2015c_rev.csv':print 'match';continue
    if buoy == '../data/deployment1/WAVE_2015d_rev.csv':print 'match';continue
    

  if dep=='2':
    if buoy == '../data/deployment2/IMB-B_2015a_rev.csv':print 'match';continue	#all BAS IMBs were deployed next to AFARs
    if buoy == '../data/deployment2/IMB-B_2015b_rev.csv':print 'match';continue
    if buoy == '../data/deployment2/IMB-B_2015c_rev.csv':print 'match';continue
    if buoy == '../data/deployment2/SNOW_2015b_rev.csv':print 'match';continue	#deployed next to SIMBA_2015c
    #if buoy == '../data/deployment2/SIMBA_2015d_rev.csv':print 'match';continue	#deployed on floe2 in March
    #if buoy == '../data/deployment2/RIDGE_2015a_rev.csv':print 'match';continue


  x = np.asarray(getColumn(buoy,4),dtype=float)[idx1::idx2]
  y = np.asarray(getColumn(buoy,3),dtype=float)[idx1::idx2]
  #up = np.asarray(getColumn(buoy,6),dtype=float)[idx1::idx2]
  #vp = np.asarray(getColumn(buoy,5),dtype=float)[idx1::idx2]
  #lat = np.asarray(getColumn(buoy,1),dtype=float)
  #lon = np.asarray(getColumn(buoy,2),dtype=float)
  #speed = np.asarray(getColumn(buoy,7),dtype=float)
  
  #calculate velocity from positions
  x = np.ma.array(x,mask=x==0)
  xsh = x[1:]
  up = (xsh-x[:-1])/float(3600*idx2)
  
  y = np.ma.array(y,mask=y==0)
  ysh = y[1:]
  vp = (ysh-y[:-1])/float(3600*idx2)
  
  #print x
  #print y
  #print up
  #print vp
  #exit() 
  
  ##interpolation from 3h to 1h
  #mask=up==0
  #up = np.ma.array(up,mask=mask)
  #sc = up.compressed()
  #dt = np.ma.array(dt,mask=mask)
  #dc = dt.compressed()
  #import pandas as pd
  #ts = pd.Series(sc, index=dc)
  #tsmv = pd.Series(up, index=dt)
  #tc_h = ts.reindex(tsmv.index, method='bfill')
  
  #plt.plot(dt,up)
  #plt.show()
  #exit()
  ##up = tc_h.values
  #plt.plot(dt,up)
  
  #vp = np.ma.array(vp,mask=mask)
  #sc = vp.compressed()
  #dt = np.ma.array(dt,mask=mask)
  #dc = dt.compressed()
  #import pandas as pd
  #ts = pd.Series(sc, index=dc)
  #tsmv = pd.Series(vp, index=dt)
  #tc_h = ts.reindex(tsmv.index, method='bfill')
  #vp = tc_h.values
  
  #x = np.ma.array(x,mask=mask)
  #sc = x.compressed()
  #dt = np.ma.array(dt,mask=mask)
  #dc = dt.compressed()
  #import pandas as pd
  #ts = pd.Series(sc, index=dc)
  #tsmv = pd.Series(x, index=dt)
  #tc_h = ts.reindex(tsmv.index, method='bfill')
  #x = tc_h.values
  
  #y = np.ma.array(y,mask=mask)
  #sc = y.compressed()
  #dt = np.ma.array(dt,mask=mask)
  #dc = dt.compressed()
  #import pandas as pd
  #ts = pd.Series(sc, index=dc)
  #tsmv = pd.Series(y, index=dt)
  #tc_h = ts.reindex(tsmv.index, method='bfill')
  #y = tc_h.values
  
  
  #print up[500:550]
  #plt.show()
  #exit()
  
  
  ##calculate velocities from the mean distance covered in the previous and next hour
  #un = np.concatenate((up[1:],[0]))
  #u = (np.ma.array(up,mask=up==0,fill_value=0)+np.ma.array(un,mask=un==0,fill_value=0))/2/1000*3600*24
  #u = u.filled()
  
  #vn = np.concatenate((vp[1:],[0]))
  #v = (np.ma.array(vp, mask=vp==0,fill_value=0)+np.ma.array(vn,mask=np==0,fill_value=0))/2/1000*3600*24
  #v = v.filled()
  
  ##th_mask = ((up[si:ei]<0.02) & (up[si:ei]>-0.02)) | (vp[si:ei]<0.02) & (vp[si:ei]>-0.02)	#exclude all velocities under error treshold
  th_mask = (up[si:ei]==0)| (vp[si:ei]==0)
  up=np.ma.array(up[si:ei]/1000*3600*24,mask=th_mask)
  vp=np.ma.array(vp[si:ei]/1000*3600*24,mask=th_mask)
  x=np.ma.array(x[si:ei]/1000,mask=x[si:ei]==0)
  y=np.ma.array(y[si:ei]/1000,mask=y[si:ei]==0)
  ##make daily averages (or leave as hourly)
  #ts = pd.Series(x, index=dates)
  #ts_d = ts.resample(pp, how=method)
  #ts_d = ts_d.interpolate()
  #data[:,0,i] = ts_d.values#[1:]	#FOR FRAMZY
  #ts = pd.Series(y, index=dates)
  #ts_d = ts.resample(pp, how=method)
  #ts_d = ts_d.interpolate()
  #data[:,1,i] = ts_d.values#[1:]
  #ts = pd.Series(up, index=dates)
  #ts_d = ts.resample(pp, how=method) 
  #ts_d = ts_d.interpolate()
  #data[:,2,i] = ts_d.values#[1:]
  #ts = pd.Series(vp, index=dates)
  #ts_d = ts.resample(pp, how=method)
  #ts_d = ts_d.interpolate()
  #data[:,3,i] = ts_d.values#[1:]
    
  #plt.plot(ts_d)
  #print ts_d.values[500:505]
  #plt.show()
  #exit()
  
  #no interpolations
  data[:,0,i] =x
  data[:,1,i] =y
  data[:,2,i] =up
  data[:,3,i] =vp
  
  
  
  i=i+1

del x,y,up,vp

data = np.where(np.isnan(data),0,data)
data = np.ma.array(data,mask=data==0)

#empty arrays for the results
div = np.zeros((data.shape[0],len(scales)-1))
shr = np.zeros((data.shape[0],len(scales)-1))
vor = np.zeros((data.shape[0],len(scales)-1))
num = np.zeros((data.shape[0],len(scales)-1))
speed = np.zeros((data.shape[0],len(scales)-1))
angle = np.zeros((data.shape[0],len(scales)-1))

##create scalling plot
#fig1 = plt.figure(figsize=(9,9))
#aax = fig1.add_subplot(111)
##aax.set_title('Scalling effect'+'\n'+title,fontsize=22)
#aax.set_title(title,fontsize=40)
#aax.set_ylabel(r"Total deformation (d$^{-1}$)",fontsize=28)
#aax.set_xlabel(r"Length scale (km)",fontsize=28)
##set the limits for the axis
#aax.set_xlim(1,200)
#aax.set_ylim(1e-5,10)
#aax.grid('on')
#aax.set_axis_bgcolor('.9')
#aax.set_xscale('log')
#aax.set_yscale('log')
##aviod scientific format of the axis thicks
#from matplotlib.ticker import ScalarFormatter
#for axis in [aax.xaxis, aax.yaxis]:
    #axis.set_major_formatter(ScalarFormatter())
#aax.tick_params(axis='both',which='major',labelsize=22)
#fig1.tight_layout()

###speed vs. deformation rate
##fig3 = plt.figure(figsize=(9,9))
##bbx = fig3.add_subplot(111)
###bbx.set_title('Scalling effect'+'\n'+title,fontsize=22)
##bbx.set_title(title,fontsize=40)
##bbx.set_ylabel(r"Total deformation (d$^{-1}$)",fontsize=28)
##bbx.set_xlabel(r"Speed (km/day)",fontsize=28)
###set the limits for the axis
##bbx.set_xlim(1,100)
##bbx.set_ylim(1e-4,100)
##bbx.grid('on')
##bbx.set_axis_bgcolor('.9')
##bbx.set_xscale('log')
##bbx.set_yscale('log')
##from matplotlib.ticker import ScalarFormatter
##for axis in [bbx.xaxis, bbx.yaxis]:
    ##axis.set_major_formatter(ScalarFormatter())
##bbx.tick_params(axis='both',which='major',labelsize=22)
##fig3.tight_layout()

#empty lists to collect the values for binning (power law scatter plot)
dr_list_h=[]
ls_list_h=[]
dr_list_l=[]
ls_list_l=[]

count = 0
tottri = 0
valn = 0

#for all time steps
for k in range(1,data.shape[0]):
  print '******************************************************************************************k:',k
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
    print 'this can be a problem!'; continue

  #get all possible combinations of 3 points - just keep unique combinations
  from itertools import combinations,chain
  
  tri = []
  tri.extend(combinations(pts, 3))
  trin = len(tri)
  
  utri = []
  utri.extend(combinations(u, 3))
  
  vtri = []
  vtri.extend(combinations(v, 3))
  
  #latlontri = []
  #latlontri.extend(combinations(latlon, 3))
  
  ##plot the triangles
  #fig = plt.figure()
  #plt.gca().set_aspect('equal')
  #import itertools
  #plt.plot(
    #*zip(*itertools.chain.from_iterable(tri)),
    #color='royalblue', marker = 'o')
  ##plt.fill(*zip(*itertools.chain.from_iterable(tri)),
    ##fill='lightblue', hatch='\\')
  #plt.title('All Triangles')
  ##fig.show()
  #tname = '%s%s' % ('triangle_comb_nice1_',k)
  #fig.savefig(path_outt+tname)
  ##exit()

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
    #ta = sqrt((ttxy(2,1)-ttxy(1,1))^2+(ttxy(2,2)-ttxy(1,2))^2);
    #tb = sqrt((ttxy(3,1)-ttxy(2,1))^2+(ttxy(3,2)-ttxy(2,2))^2);
    #tc = sqrt((ttxy(3,1)-ttxy(1,1))^2+(ttxy(3,2)-ttxy(1,2))^2);
    ta = np.sqrt((vert[1,0]-vert[0,0])**2+(vert[1,1]-vert[0,1])**2)
    tb = np.sqrt((vert[2,0]-vert[1,0])**2+(vert[2,1]-vert[1,1])**2)
    tc = np.sqrt((vert[2,0]-vert[0,0])**2+(vert[2,1]-vert[0,1])**2)
    #if ta < 20 or tb < 20 or tc < 20: continue		#shorter than 100m
    #print ta
    #print tb
    #print tc
    
    ##calculate the side lenghts from the latlon
    #llvert = np.asarray(latlontri[t])
    #print llvert
    
    #from LatLon import LatLon
    #palmyra = LatLon(llvert[0,0],llvert[0,1]) # Location of Palmyra Atoll
    #honolulu = LatLon(llvert[1,0],llvert[1,1]) # Location of Honolulu, HI
    #ta = palmyra.distance(honolulu) # WGS84 distance in km  
    
    #palmyra = LatLon(llvert[2,0],llvert[2,1]) # Location of Palmyra Atoll
    #honolulu = LatLon(llvert[1,0],llvert[1,1]) # Location of Honolulu, HI
    #tb = palmyra.distance(honolulu) # WGS84 distance in km  
    
    #palmyra = LatLon(llvert[0,0],llvert[0,1]) # Location of Palmyra Atoll
    #honolulu = LatLon(llvert[2,0],llvert[2,1]) # Location of Honolulu, HI
    #tc = palmyra.distance(honolulu) # WGS84 distance in km  
     
    #if ta > 200 or tb > 200 or tc > 200: continue		#longer than 500km
    
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
  
  
    ##example of a valid triangle
    #fig = plt.figure()
    #plt.gca().set_aspect('equal')
    #import itertools
    #plt.plot(
      #*zip(*itertools.chain.from_iterable(itertools.combinations(pts, 3))),
      #color = 'royalblue', marker = 'o')
    #plt.plot(vert[0][0],vert[0][1],'o',color = 'purple')
    #plt.plot(vert[1][0],vert[1][1],'o',color = 'purple')
    #plt.plot(vert[2][0],vert[2][1],'o',color = 'purple')
    #plt.title('All Triangles')
    ##fig.show()
    #tname = '%s%s' % ('triangle_comb_nice1_',k)
    #fig.savefig(path_outt+tname)
    ##exit()

    
    ##du/dx (calculate line integrals)
    #from Gunnar's matlab script (RGPS_MITgcm_lagr_deform.m)
    #rdux(ii) = 0.5/triarea(ii) * ( ...
	#(ru(2)+ru(1))*(ttxy(2,2)-ttxy(1,2)) + ...
	#(ru(3)+ru(2))*(ttxy(3,2)-ttxy(2,2)) + ...
	#(ru(1)+ru(3))*(ttxy(1,2)-ttxy(3,2)) );
    #rduy(ii) = -0.5/triarea(ii) * ( ...
	#(ru(2)+ru(1))*(ttxy(2,1)-ttxy(1,1)) + ...
	#(ru(3)+ru(2))*(ttxy(3,1)-ttxy(2,1)) + ...
	#(ru(1)+ru(3))*(ttxy(1,1)-ttxy(3,1)) );

    #rdvx(ii) = 0.5/triarea(ii) * ( ...
	#(rv(2)+rv(1))*(ttxy(2,2)-ttxy(1,2)) + ...
	#(rv(3)+rv(2))*(ttxy(3,2)-ttxy(2,2)) + ...
	#(rv(1)+rv(3))*(ttxy(1,2)-ttxy(3,2)) );
    #rdvy(ii) = -0.5/triarea(ii) * ( ...
	#(rv(2)+rv(1))*(ttxy(2,1)-ttxy(1,1)) + ...
	#(rv(3)+rv(2))*(ttxy(3,1)-ttxy(2,1)) + ...
	#(rv(1)+rv(3))*(ttxy(1,1)-ttxy(3,1)) );

    #also Hutchings et al, 2012 (ISPOL paper)
    
    #from Jenny Hutching, SEDNA routines (IDL code), strain_nogaps_la_v2.pro
    #; Estimate strain rates
    #;
	#ux=dblarr(npts)*0.
	#uy=ux
	#vx=ux
	#vy=ux

	#for nn=0,nbuoys-2 do begin
	  #ux = ux + (u(nn+1,0,*)+u(nn,0,*)) * (xytrack(nn+1,1,*)-xytrack(nn,1,*)) 
	  #uy = uy + (u(nn+1,0,*)+u(nn,0,*)) * (xytrack(nn+1,0,*)-xytrack(nn,0,*))
	  #vx = vx + (u(nn+1,1,*)+u(nn,1,*)) * (xytrack(nn+1,1,*)-xytrack(nn,1,*))
	  #vy = vy + (u(nn+1,1,*)+u(nn,1,*)) * (xytrack(nn+1,0,*)-xytrack(nn,0,*))
	#endfor
	#ux = ux $ 
	#+ (u(0,0,*)+u(nbuoys-1,0,*)) * (xytrack(0,1,*)-xytrack(nbuoys-1,1,*))
	#uy = uy $
	#+ (u(0,0,*)+u(nbuoys-1,0,*)) * (xytrack(0,0,*)-xytrack(nbuoys-1,0,*))
	#vx = vx $
	#+ (u(0,1,*)+u(nbuoys-1,1,*)) * (xytrack(0,1,*)-xytrack(nbuoys-1,1,*))
	#vy = vy $
	#+ (u(0,1,*)+u(nbuoys-1,1,*)) * (xytrack(0,0,*)-xytrack(nbuoys-1,0,*))
	#ux=ux*0.5/area 
	#uy=uy*0.5/area 
	#vx=vx*0.5/area 
	#vy=vy*0.5/area 

    #;
    #; Principle components of strain rate
    #;
	#e1=ux+vy
	#e2=0.5*sqrt((ux-vy)*(ux-vy)+(uy+vx)*(uy+vx))
      
    #;
    #; Strain rate components
    #;
	#vorticity = vx-uy
	#shear = uy+vx
	#normal = ux-vy

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
    
    #if ls > 80: print '%s%s%s%s%s%s' %('side lenghts: ',ta,' ',tb,' ',tc)
    #separate the high and low speeds and keep all values for binning and averaging
    #if dr > 1000: continue		#skip wierd values (SEDNA)

    #sort by speed 
    #if tris[t] < .02*wspeed[k]*3600*24/1000:
    #if tris[t] > 8.64: 	#0.1m/s
    #aax.scatter(ls,dr,color='royalblue',alpha=.4)
    dr_list_l.append(dr)
    ls_list_l.append(ls)
    #print dr
    #else:
      #aax.scatter(ls,dr,color='purple',alpha=.4)
      #dr_list_h.append(dr)
      #ls_list_h.append(ls)
    
    ##sort by lenght scale
    ##complete lists for speed vs. deformation
    #if ls < scales[1]:
      #bbx.scatter(tris[t],dr,color='b',alpha=.2)
    #elif ls < scales[2]:
      #bbx.scatter(tris[t],dr,color='g',alpha=.2)
    #elif ls < scales[3]:
      #bbx.scatter(tris[t],dr,color='r',alpha=.2)
    #else:
      #bbx.scatter(tris[t],dr,color='k',alpha=.2)
 
    count = count + 1

  #mask = (dux==0)
  #get timeseries for lenght scale
  for i in range(0,len(scales)-1):
    mask = ((scales[i] < lsa) & (lsa > scales[i+1]) | (dux==0))
        
    #divergence
    tridiv = np.ma.array(dux,mask=mask) + np.ma.array(dvy,mask=mask)    
    div[k,i] = np.mean(tridiv)	# /day
        
    #number of valid triangles
    num[k,i] = np.ma.count(tridiv)
 
    #shear (pure and normal)
    pshear = np.mean(np.ma.array(dux,mask=mask)-np.ma.array(dvy,mask=mask))
    nshear = np.mean(np.ma.array(duy,mask=mask)+np.ma.array(dvx,mask=mask))
    shr[k,i] = .5 * np.sqrt(pshear**2+nshear**2)	# /day
  
    #vorticity
    trivor = np.ma.array(dvx,mask=mask) - np.ma.array(duy,mask=mask)
    vor[k,i] = np.mean(trivor)	# /day
    
    #speed
    speed[k,i] = np.mean(np.ma.array(tris,mask=mask))	#km/day
    
    #direction
    sin = np.sin(np.radians(np.ma.array(trid,mask=mask)))
    sm = np.mean(sin)
    cos = np.cos(np.radians(np.ma.array(trid,mask=mask)))
    cm = np.mean(cos)

    if np.invert(np.ma.is_masked(sm)):
      angle[k,i] = np.degrees(np.arctan2(sm,cm))
      angle[k,i] = np.where(angle[k,i]<0,angle[k,i]+360,angle[k,i])
      
    tottri = tottri + num[k]
    valn = valn+1
  
  del tridiv, pshear, nshear, trivor, dux, duy, dvx, dvy, tris, pts, u, v, lsa

print 'total number of all valid triangles t:',tottri, 'in time steps n:',valn, 'over hours h:', k

#binning and averaging
#ls_list_h = np.array(ls_list_h)
#dr_list_h = np.array(dr_list_h)

#ls_list_l = np.array(ls_list_l)
#dr_list_l = np.array(dr_list_l)

#from scipy import stats
#binval = 5
#binval = scales

##sortls_l = sorted(ls_list_l)
##sortdr_l = [x for (y,x) in sorted(zip(dr_list_l,ls_list_l))]

##sortls_h = sorted(ls_list_h)
##sortdr_h = [x for (y,x) in sorted(zip(dr_list_h,ls_list_h))]


#print len(ls_list_l)
#print len(ls_list_h)

##binning the data
#bin_means_l, bin_edges, binnumber = stats.binned_statistic(ls_list_l,
                #dr_list_l, statistic='mean', bins=binval)
###aax.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=5,
           ###label='high bins')
##bin_means_h, bin_edges, binnumber = stats.binned_statistic(ls_list_h,
                ##dr_list_h, statistic='mean', bins=binval)

##finding the power law
#x = (bin_edges[1:] + bin_edges[:-1])/2
#print x

#y_l = bin_means_l
##y_h = bin_means_h
#print dr_list_l
#print y_l
##print y_h

#logx=np.log(x)
#logyl=np.log(y_l)
##logyh=np.log(y_h)
#p=np.polyfit(logx,logyl,1)
#aax.plot(x,y_l,'*',color='royalblue',ms=10,mew=1,label='compact ice')
#k=p[0]
#loga=p[1]
#a=np.exp(loga)
#aax.plot(x,a*x**k,'--k',linewidth=2,label='y=%.2f*x^%.2f' %(a,k))
#aax.plot(x,a*x**-.2,'--k',linewidth=1,label='y=%.2f*x^%.2f' %(a,-.2))	#Marsan et al, 2004

##p=np.polyfit(logx,logyh,1)
##aax.plot(x,y_h,'*',color='purple',ms=10,mew=1,label='free drift')
##k=p[0]
##loga=p[1]
##a=np.exp(loga)
##aax.plot(x,a*x**k,'-.k',linewidth=3,label='y=%.2f*x^%.2f' %(a,k))

#aax.legend(loc='lower left',prop={'size':15}, fancybox=True, framealpha=0.5,numpoints=1)
###plt.show()

#fig1.savefig(path_out+outname2)

#plot speed vs. deformation rate
#fig3.savefig(path_out+outname3)

#save all the data for the time series
div.dump(path_out+'div_'+outname2)
shr.dump(path_out+'shr_'+outname2)
vor.dump(path_out+'vor_'+outname2)
speed.dump(path_out+'speed_'+outname2)
angle.dump(path_out+'angle_'+outname2)
num.dump(path_out+'num_'+outname2)

#and for the power law
ls_list_l = np.array(ls_list_l)
dr_list_l = np.array(dr_list_l)
ls_list_l.dump(path_out+'ls_'+outname2)
dr_list_l.dump(path_out+'dr_'+outname2)

print '%s%s%s' %('This calculation is based on positons of ',len(flist),' buoys.')
