#! /usr/bin/python

import numpy as np
import csv

def deformation(vert,uvert,vvert):

  #skip all triangles with the too short or too long sides
  #sides of triangle
  ta = np.sqrt((vert[1,0]-vert[0,0])**2+(vert[1,1]-vert[0,1])**2)
  tb = np.sqrt((vert[2,0]-vert[1,0])**2+(vert[2,1]-vert[1,1])**2)
  tc = np.sqrt((vert[2,0]-vert[0,0])**2+(vert[2,1]-vert[0,1])**2)
  #if ta < 20 or tb < 20 or tc < 20: continue		#shorter than 100m
  
  #find smallest angle to filter out too acute triangles
  sides = [ta,tb,tc]
  srt = min(sides)
  srtid = np.argmin(sides)
  lng = np.delete(sides,srtid)
  #minang(ii) = acosd((sides(1)^2+sides(2)^2-shortest^2)/2/sides(1)/sides(2)); % law of cosines
  minang = np.degrees(np.arccos((lng[0]**2+lng[1]**2-srt**2)/2/lng[0]/lng[1]))
  #if minang < 2: return 0,0,0,0,0
  
  ##area of triangle
  #Hutchings etal, 2012 (ISPOL)
  area = .5* (vert[0,0]*vert[1,1] - vert[0,1]*vert[1,0] + vert[1,0]*vert[2,1] - vert[1,1]*vert[2,0] + vert[2,0]*vert[0,1] - vert[2,1]*vert[0,0])
  if area < 0: print 'not counter-clockwise oriented triangle!'; exit()
  
  #some filtering to throw out the erroneous GPS positioning. If the accuracy is 2m.
  if area < 1: print area; return 0,0,0,0,0	#Hutchings etal, 2011 (SEDNA)
    
  #calculate line integrals
  dux = (.5 / area) * (
    (uvert[1]+uvert[0])*(vert[1,1]-vert[0,1])+
    (uvert[2]+uvert[1])*(vert[2,1]-vert[1,1])+
    (uvert[0]+uvert[2])*(vert[0,1]-vert[2,1]) )
  
  duy = -(.5 / area) * (
    (uvert[1]+uvert[0])*(vert[1,0]-vert[0,0])+
    (uvert[2]+uvert[1])*(vert[2,0]-vert[1,0])+
    (uvert[0]+uvert[2])*(vert[0,0]-vert[2,0]) )
  
  dvx = (.5 / area) * (
    (vvert[1]+vvert[0])*(vert[1,1]-vert[0,1])+
    (vvert[2]+vvert[1])*(vert[2,1]-vert[1,1])+
    (vvert[0]+vvert[2])*(vert[0,1]-vert[2,1]) )
    
  dvy = -(.5 / area) * (
    (vvert[1]+vvert[0])*(vert[1,0]-vert[0,0])+
    (vvert[2]+vvert[1])*(vert[2,0]-vert[1,0])+
    (vvert[0]+vvert[2])*(vert[0,0]-vert[2,0]) )

  return dux,duy,dvx,dvy,minang

def getColumn(filename, column, delimiter=','):
    results = csv.reader(open(filename),delimiter=delimiter)
    next(results, None)
    return [result[column] for result in results]
