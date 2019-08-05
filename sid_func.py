#! /usr/bin/python

import numpy as np
import csv
import osr, gdal

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
  if area < 0: print('not counter-clockwise oriented triangle!'); exit()
  
  #some filtering to throw out the erroneous GPS positioning. If the accuracy is 2m.
  if area < 1: print(area); return 0,0,0,0,0	#Hutchings etal, 2011 (SEDNA)
    
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
  
  #B&R, 2015 filtering
  

  return dux,duy,dvx,dvy,minang,area

def getColumn(filename, column, delimiter=','):
    results = csv.reader(open(filename),delimiter=delimiter)
    next(results, None)
    return [result[column] for result in results]

def logfit(xdata,ydata):
    logx=np.log10(xdata)
    logy=np.log10(ydata)

    # fit a curve to the data using a least squares 1st order polynomial fit
    z = np.polyfit(logx,logy,1)
    p = np.poly1d(z)
    fit = p(logx)

    # get the coordinates for the fit curve
    c_y = [np.min(fit),np.max(fit)]
    c_x = [np.min(logx),np.max(logx)]

    # predict y values of origional data using the fit
    p_y = z[0] * logx + z[1]

    # calculate the y-error (residuals)
    y_err = logy -p_y

    # create series of new test x-values to predict for
    p_x = np.arange(np.min(logx),np.max(logx)+.01,.01)

    # now calculate confidence intervals for new test x-series
    mean_x = np.mean(logx)         # mean of x
    n = len(logx)              # number of samples in origional fit
    from scipy import stats
    tval = stats.t.ppf(1-0.005, n)		# appropriate t value for 2-tailed distribution
    s_err = np.sum(np.power(y_err,2))   # sum of the squares of the residuals

    confs = tval * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((p_x-mean_x),2)/
                ((np.sum(np.power(logx,2)))-n*(np.power(mean_x,2))))))

    # now predict y based on test x-values
    p_y = z[0]*p_x+z[1]

    # get lower and upper confidence limits based on predicted y and confidence intervals
    lower = p_y - abs(confs)
    upper = p_y + abs(confs)

    #get them back on the exponential scale
    k=z[0]
    loga=z[1]
    a=10.0**loga
    ciy_low = 10.0**lower
    ciy_upp = 10.0**upper
    cix = 10.0**p_x

    return(a,k,cix,ciy_upp,ciy_low)

#consider this option too:
#https://stackoverflow.com/questions/26851533/fit-a-power-law-function-to-the-data-with-both-x-and-y-errors

def density(x,y):
    from scipy.stats import gaussian_kde
    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    return(x,y,z)


#log-spaced bins
def density_lsb( x,y,n=20):
    """
    Scatter plot colored by 2d histogram
    """
    from scipy.interpolate import interpn
    vmin = np.min(x); vmax = np.max(x)
    lsbx = np.logspace(np.log10(vmin),np.log10(vmax), n)
    vmin = np.min(y); vmax = np.max(y)
    lsby = np.logspace(np.log10(vmin),np.log10(vmax), n)
    data , x_e, y_e = np.histogram2d( x, y, bins=[lsbx,lsby])
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    return(x,y,z)


def convertXY(xy_source, inproj, outproj):
    # function to convert coordinates

    shape = xy_source[0,:,:].shape
    size = xy_source[0,:,:].size

    # the ct object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inproj, outproj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))

    xx = xy_target[:,0].reshape(shape)
    yy = xy_target[:,1].reshape(shape)

    return xx, yy
