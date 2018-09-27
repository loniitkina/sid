from sid_func import *
import matplotlib.pyplot as plt

#plotting
fig1 = plt.figure(figsize=(9,9))
fig1 = plt.figure(figsize=(7.5,7))
ax = fig1.add_subplot(111, axisbg='.9')
#title = 'ship_radar+buoys+SAR_select'
#title = 'ship_radar+buoys'
title = 'ship_radar+buoys+SAR_in'
ax.set_title(title,fontsize=29, loc='left')
ax.set_xlabel(r"Length scale (km)",fontsize=25)
ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
ax.set_xscale('log')
ax.set_yscale('log')

meanls_list_sr=[]
meantd_list_sr=[]

meanls_list_sar=[]
meantd_list_sar=[]

meanls_list_b=[]
meantd_list_b=[]


#ship radar data
inpath = '../data/ship_radar/24h/'
outpath = '../plots/'
fname_start = 'F1_Deformation_L'
lsc_list = range(1,7)

for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print scale
    fname = inpath+fname_start+str(scale)+'_24h.txt'
    
    ls = getColumn(fname,0, delimiter=' ')
    td = getColumn(fname,1, delimiter=' ')
    ls = np.array(ls,dtype=np.float)/1000  #convert from m to km
    td = np.array(td,dtype=np.float)/3600      #convert from h to s
    
    #calculate and store averages
    meanls=np.mean(ls)
    meantd=np.mean(td)
    meanls_list_sr.append(meanls)
    meantd_list_sr.append(meantd)
    
    #plot all the data
    ax.scatter(ls, td, lw=0, alpha=.2)  # Data
    ax.plot(meanls,meantd,'s',markersize=7,markeredgecolor='k')

#buoy data - do we have enough of 1 day data (should be enough for the entire leg 1)
#scales 2-100km
inpath = '../data/buoys/'
outpath = '../plots/'
#fname_start = 'nice1_comb_SAR'
fname_start = 'nice1_in_SAR'
#fname_start = 'nice1_outer_SAR'
#fname_start = 'nice1_late_SAR'
#fname_start = 'nice1_early_SAR'
#fname_start = 'nice1_select_SAR'
fname_td1 = inpath+'dr_'+fname_start+'1_'+'24h'
fname_ls1 = inpath+'ls_'+fname_start+'1_'+'24h'
fname_td2 = inpath+'dr_'+fname_start+'2_'+'24h'
fname_ls2 = inpath+'ls_'+fname_start+'2_'+'24h'
 
td = np.append(np.load(fname_td1),np.load(fname_td2))/24/60/60      #convert from days to s
ls = np.append(np.load(fname_ls1),np.load(fname_ls2))

#throw away very high deformation rates (unrealistic values)
mask = td>10e-3
ls = np.ma.array(ls,mask=mask)
td = np.ma.array(td,mask=mask)
ls = np.ma.compressed(ls)
td = np.ma.compressed(td)   

##divide in lengh scale classes
#lsc_list = [1]   #order number
#minlen = [1]
#maxlen = [200]

#inner ring
lsc_list = [1,2,3]   #order number
minlen = [2,4,7]
maxlen = [4,7,10]


for i in range(0,len(lsc_list)):
    print i
    mask = (ls<minlen[i]) | (ls>maxlen[i])
    ls_class = np.ma.array(ls,mask=mask)
    td_class = np.ma.array(td,mask=mask)
    ls_class = np.ma.compressed(ls_class)
    td_class = np.ma.compressed(td_class)   
      
    #calculate and store averages
    meanls=np.mean(ls_class)
    meantd=np.mean(td_class)
    meanls_list_b.append(meanls)
    meantd_list_b.append(meantd)
    
    #plot all the data
    ax.scatter(ls_class, td_class, lw=0, alpha=.2)  # Data
    ax.plot(meanls,meantd,'o',markersize=7,markeredgecolor='k')

#SAR data
inpath = '../data/Sentinel1_def_24h_'
outpath = '../plots/'
fname_start = 'td_leg1_L'
lsc_list = [25,50,100,200,500]   #not ls but number of nominal grid points
minlen = [4,2,1,.5,.2]
maxlen = [6,3,1.5,.75,.3]


for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print scale
    fname = inpath+str(scale)+'/'+fname_start+str(scale)+'_15km.csv'
    
    ls = getColumn(fname,0, delimiter=',')
    td = getColumn(fname,1, delimiter=',')
    ang = getColumn(fname,2)
    ls = np.array(ls,dtype=np.float)
    td = np.array(td,dtype=np.float)
    ang = np.array(ang,dtype=np.float)
    
    #mask out all the acute triangles
    mask=ang<15
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)/1000  #convert from m to km
    td = np.ma.compressed(td)
    
    #mask all very small or big triangles
    #if not masked the renge of the ls is big and has several clouds (expected ls, twide the ls and all kinds of smaller ls)
    mask = (ls<minlen[i]) | (ls>maxlen[i])
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)   
    
    #calculate and store averages
    meanls=np.mean(ls)
    meantd=np.mean(td); print meantd
    meanls_list_sar.append(meanls)
    meantd_list_sar.append(meantd)
    
    #plot all the data
    ax.scatter(ls, td, lw=0, alpha=.2)  # Data
    ax.plot(meanls,meantd,'*',markersize=10,markeredgecolor='k')

#ship radar
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sr,meantd_list_sr)

#dummy x data for plotting
x = np.arange(min(meanls_list_sr), max(meanls_list_sr), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='m')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#buoys
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_b,meantd_list_b)

#dummy x data for plotting
x = np.arange(min(meanls_list_b), max(meanls_list_b), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='g')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#and separate for SAR
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar,meantd_list_sar)

#dummy x data for plotting
x = np.arange(min(meanls_list_sar), max(meanls_list_sar), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='purple')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

ax.grid('on')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())    
ax.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig1.tight_layout()

fig1.savefig(outpath+'power_law_24h_'+title)
