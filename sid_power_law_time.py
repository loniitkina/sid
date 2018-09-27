from sid_func import *
import matplotlib.pyplot as plt

#plotting
fig1 = plt.figure(figsize=(9,9))
fig1 = plt.figure(figsize=(7.5,7))
ax = fig1.add_subplot(111, axisbg='.9')
#title = 'ship_radar+buoys+SAR'
title = 'ship_radar+buoys'
ax.set_title(title,fontsize=29, loc='left')
ax.set_xlabel(r"Time scale (h)",fontsize=25)
ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
ax.set_xscale('log')
ax.set_yscale('log')

meanls_list=[]
meantd_list=[]

meanls_list_sar=[]
meantd_list_sar=[]

#ship radar data
inpath = '../data/ship_radar/L6/'
outpath = '../plots/'
fname_start = 'F1_Deformation_L6_'
lsc_list = [1,3,6,12,24]

for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print scale
    fname = inpath+fname_start+str(scale)+'h.txt'
    
    ls = scale
    td = getColumn(fname,1, delimiter=' ')
    #ls = np.array(ls,dtype=np.float)/1000  #convert from m to km
    td = np.array(td,dtype=np.float)/3600      #convert from h to s
    
    #calculate and store averages
    meanls=ls
    meantd=np.mean(td)
    meanls_list.append(meanls)
    meantd_list.append(meantd)
    
    #plot all the data
    #ax.scatter(ls, td, lw=0, alpha=.2)  # Data
    ax.plot(meanls,meantd,'s',markersize=7,markeredgecolor='k')

#buoy data - do we have enough of 1 day data (should be enough for the entire leg 1)
#scales 2-100km
inpath = '../data/buoys/'
fname_start = 'nice1_in_SAR'
lsc_list = [1,3,6,12,24]

for i in lsc_list:
    scale = i
    print scale
    fname_td1 = inpath+'dr_'+fname_start+'1_'+str(scale)+'h'
    fname_ls1 = inpath+'ls_'+fname_start+'1_'+str(scale)+'h'
    fname_td2 = inpath+'dr_'+fname_start+'2_'+str(scale)+'h'
    fname_ls2 = inpath+'ls_'+fname_start+'2_'+str(scale)+'h'
    
    td = np.append(np.load(fname_td1),np.load(fname_td2))/24/60/60      #convert from days to s
    ls = np.append(np.load(fname_ls1),np.load(fname_ls2))
    
    print len(ls)
    
    #throw away very high deformation rates (unrealistic values)
    mask = td>10e-3
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)   

    #limit the lengh scale to 4-5km
    minlen = 3
    maxlen = 5

    mask = (ls<minlen) | (ls>maxlen)
    #ls_class = np.ma.array(ls,mask=mask)
    td_class = np.ma.array(td,mask=mask)
    #ls_class = np.ma.compressed(ls_class)
    td_class = np.ma.compressed(td_class)  
    
    #print td_class
    #exit()
    
    #calculate and store averages
    meanls=i
    meantd=np.mean(td_class)
    meanls_list.append(meanls)
    meantd_list.append(meantd)
    
    #print meanls, meantd
    
    #plot all the data
    #ax.scatter(ls_class, td_class, lw=0, alpha=.2)  # Data
    ax.plot(meanls,meantd,'o',markersize=7,markeredgecolor='k')

##SAR data
#inpath = '../data/Sentinel1_def_24h_'
#outpath = '../plots/'
#fname_start = 'td_leg1_L'
#lsc_list = [7,10,25,50,100,200,500]   #not ls but number of nominal grid points
#minlen = [20,10,4,2,1,.5,.2]
#maxlen = [50,18,6,3,1.5,.75,.3]


#for i in range(0,len(lsc_list)):
    #scale = lsc_list[i]
    #print scale
    #fname = inpath+str(scale)+'/'+fname_start+str(scale)+'.csv'
    
    #ls = getColumn(fname,0, delimiter=',')
    #td = getColumn(fname,1, delimiter=',')
    #ang = getColumn(fname,2)
    #ls = np.array(ls,dtype=np.float)
    #td = np.array(td,dtype=np.float)
    #ang = np.array(ang,dtype=np.float)
    
    ##mask out all the acute triangles
    #mask=ang<15
    #ls = np.ma.array(ls,mask=mask)
    #td = np.ma.array(td,mask=mask)
    #ls = np.ma.compressed(ls)/1000  #convert from m to km
    #td = np.ma.compressed(td)
    
    ##mask all very small or big triangles
    ##if not masked the renge of the ls is big and has several clouds (expected ls, twide the ls and all kinds of smaller ls)
    #mask = (ls<minlen[i]) | (ls>maxlen[i])
    #ls = np.ma.array(ls,mask=mask)
    #td = np.ma.array(td,mask=mask)
    #ls = np.ma.compressed(ls)
    #td = np.ma.compressed(td)   
    
    ##calculate and store averages
    #meanls=np.mean(ls)
    #meantd=np.mean(td)
    #meanls_list_sar.append(meanls)
    #meantd_list_sar.append(meantd)
    
    ##plot all the data
    #ax.scatter(ls, td, lw=0, alpha=.2)  # Data
    #ax.plot(meanls,meantd,'*',markersize=10,markeredgecolor='k')


##fit the line
#xdata = meanls_list
#ydata = meantd_list
#logx=np.log10(xdata)
#logy=np.log10(ydata)

## fit a curve to the data using a least squares 1st order polynomial fit
#z = np.polyfit(logx,logy,1)
#p = np.poly1d(z)
#fit = p(logx)

## get the coordinates for the fit curve
#c_y = [np.min(fit),np.max(fit)]
#c_x = [np.min(logx),np.max(logx)]

## predict y values of origional data using the fit
#p_y = z[0] * logx + z[1]

## calculate the y-error (residuals)
#y_err = logy -p_y

## create series of new test x-values to predict for
#p_x = np.arange(np.min(logx),np.max(logx)+.01,.01)

## now calculate confidence intervals for new test x-series
#mean_x = np.mean(logx)         # mean of x
#n = len(logx)              # number of samples in origional fit
#from scipy import stats
#tval = stats.t.ppf(1-0.005, n)		# appropriate t value for 2-tailed distribution
#s_err = np.sum(np.power(y_err,2))   # sum of the squares of the residuals

#confs = tval * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((p_x-mean_x),2)/
            #((np.sum(np.power(logx,2)))-n*(np.power(mean_x,2))))))

## now predict y based on test x-values
#p_y = z[0]*p_x+z[1]

## get lower and upper confidence limits based on predicted y and confidence intervals
#lower = p_y - abs(confs)
#upper = p_y + abs(confs)

##get them back on the exponential scale
#k=z[0]
#loga=z[1]
#a=10.0**loga
#ciy_low = 10.0**lower
#ciy_upp = 10.0**upper
#cix = 10.0**p_x

##dummy x data for plotting
#x = np.arange(min(xdata), max(xdata), 1)

#ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k))
#ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
#ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

##and separate for SAR
##fit the line
#xdata = meanls_list_sar
#ydata = meantd_list_sar
#logx=np.log10(xdata)
#logy=np.log10(ydata)

## fit a curve to the data using a least squares 1st order polynomial fit
#z = np.polyfit(logx,logy,1)
#p = np.poly1d(z)
#fit = p(logx)

## get the coordinates for the fit curve
#c_y = [np.min(fit),np.max(fit)]
#c_x = [np.min(logx),np.max(logx)]

## predict y values of origional data using the fit
#p_y = z[0] * logx + z[1]

## calculate the y-error (residuals)
#y_err = logy -p_y

## create series of new test x-values to predict for
#p_x = np.arange(np.min(logx),np.max(logx)+.01,.01)

## now calculate confidence intervals for new test x-series
#mean_x = np.mean(logx)         # mean of x
#n = len(logx)              # number of samples in origional fit
#from scipy import stats
#tval = stats.t.ppf(1-0.005, n)		# appropriate t value for 2-tailed distribution
#s_err = np.sum(np.power(y_err,2))   # sum of the squares of the residuals

#confs = tval * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((p_x-mean_x),2)/
            #((np.sum(np.power(logx,2)))-n*(np.power(mean_x,2))))))

## now predict y based on test x-values
#p_y = z[0]*p_x+z[1]

## get lower and upper confidence limits based on predicted y and confidence intervals
#lower = p_y - abs(confs)
#upper = p_y + abs(confs)

##get them back on the exponential scale
#k=z[0]
#loga=z[1]
#a=10.0**loga
#ciy_low = 10.0**lower
#ciy_upp = 10.0**upper
#cix = 10.0**p_x

##dummy x data for plotting
#x = np.arange(min(xdata), max(xdata), 1)

#ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k))
#ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
#ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)



#ax.set_xlim(1,400)
#ax.set_ylim(1e-4,100)
ax.grid('on')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())    
ax.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig1.tight_layout()

fig1.savefig(outpath+'power_law_4-5km_'+title)
