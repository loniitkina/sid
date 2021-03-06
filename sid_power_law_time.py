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

meantd_list_sr=[]
meantd_list_b=[]

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
    meantd_list_sr.append(meantd)
    
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
    #meanls_list.append(meanls)
    meantd_list_b.append(meantd)
    
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


#ship radar
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,meantd_list_sr)

#dummy x data for plotting
x = np.arange(min(meanls_list), max(meanls_list), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='m')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#buoys
#fit the line to bins
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list[1:],meantd_list_b[1:])

ax.loglog(x[1:],a*x[1:]**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='g')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)


ax.grid('on')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())    
ax.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig1.tight_layout()

fig1.savefig(outpath+'power_law_3-5km_'+title)
