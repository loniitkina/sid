#https://www.thoughtco.com/what-are-moments-in-statistics-3126234
#Statistical moments are basic calculations that can be used to find a probability distribution's mean, variance, and skewness.
#Suppose that we have a set of data with a total of n discrete points. One important calculation, which is actually several numbers, is called the sth moment. 
#The sth moment of the data set with values x1, x2, x3, . . . , xn is given by the formula:
#(x1s + x2s + x3s + . . . + xns)/n
#s=1, 1st moment=mean
#s=2, 2nd moment=variance
#s=3, 3rd moment=skewness - if normalized (symetricity)
#s=4, 4th moment=kurtosis - if normalized (length of tail)

from sid_func import *
import matplotlib.pyplot as plt

#plotting
fig1 = plt.figure(figsize=(10,5))
ax = fig1.add_subplot(121)
title = 'moments_new'
ax.set_title(title,fontsize=29, loc='left')
ax.set_xlabel(r"Length scale (km)",fontsize=25)
ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
ax.set_xscale('log')
ax.set_yscale('log')

#ship radar data
inpath = '../data/ship_radar/24h/'
outpath = '../plots/'
fname_start = 'F1_Deformation_L'
lsc_list = range(1,7)

meanls_list=[]
meantd_list=[]
mom2_list=[]
mom3_list=[]
mom4_list=[]

for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print(scale)
    fname = inpath+fname_start+str(scale)+'_24h.txt'
    
    ls = getColumn(fname,0, delimiter=' ')
    td = getColumn(fname,1, delimiter=' ')
    ls = np.array(ls,dtype=np.float)/1000  #convert from m to km
    td = np.array(td,dtype=np.float)/3600      #convert from h to s
    
    #calculate and store moments
    meanls=np.mean(ls)
    meantd=np.mean(td)
    print(meantd)
    from scipy.stats import moment
    moments = moment(td,moment=[1,2,3,4])  #Calculate the nth moment about the mean for a sample. Attention: this is already a central moment, 1st moment==0
    #if ordinary moments are wanted, they must be calculated by own function: np.mean(x^n)

    meanls_list.append(meanls)
    meantd_list.append(meantd)
    
    
    #plot all the data
    #ax.scatter(ls, td, lw=0, alpha=.2)  # Data
    ax.plot(meanls,meantd,'s',markersize=7,markeredgecolor='k',color='r')
    
    ax.plot(meanls,moments[1],'s',markersize=7,markeredgecolor='k',color='r')
    ax.plot(meanls,moments[2],'s',markersize=7,markeredgecolor='k',color='r')
    ax.plot(meanls,moments[3],'s',markersize=7,markeredgecolor='k',color='r')
    
    mom2_list.append(moments[1])
    mom3_list.append(moments[2])
    mom4_list.append(moments[3])

slopes_rad=[]

a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,meantd_list)
#dummy x data for plotting
x = np.arange(min(meanls_list), max(meanls_list), 1)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='orange',ls='-')
slopes_rad.append(k*-1)

a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,mom2_list)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='salmon',ls='-')
slopes_rad.append(k*-1)
    
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,mom3_list)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='royalblue',ls='-')
slopes_rad.append(k*-1)
    
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,mom4_list)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='darkred',ls='-')
slopes_rad.append(k*-1)

    
#buoy data - do we have enough of 1 day data (should be enough for the entire leg 1)
#scales 2-100km
inpath = '../data/buoys/'
outpath = '../plots/'
fname_start = 'nice1_in_SAR'
fname_td1 = inpath+'dr_'+fname_start+'1_'+'24h'
fname_ls1 = inpath+'ls_'+fname_start+'1_'+'24h'
fname_td2 = inpath+'dr_'+fname_start+'2_'+'24h'
fname_ls2 = inpath+'ls_'+fname_start+'2_'+'24h'
 
td = np.append(np.load(fname_td1,encoding='latin1'),np.load(fname_td2,encoding='latin1'))/24/60/60      #convert from days to s
ls = np.append(np.load(fname_ls1,encoding='latin1'),np.load(fname_ls2,encoding='latin1'))

#throw away very high deformation rates (unrealistic values)
mask = td>10e-3
ls = np.ma.array(ls,mask=mask)
td = np.ma.array(td,mask=mask)
ls = np.ma.compressed(ls)
td = np.ma.compressed(td)   

#inner ring
lsc_list = [1,2,3]   #order number
minlen = [2,4,7]
maxlen = [4,7,10]

meanls_list=[]
meantd_list=[]
mom2_list=[]
mom3_list=[]
mom4_list=[]

for i in range(0,len(lsc_list)):
    print(i)
    mask = (ls<minlen[i]) | (ls>maxlen[i])
    ls_class = np.ma.array(ls,mask=mask)
    td_class = np.ma.array(td,mask=mask)
    ls_class = np.ma.compressed(ls_class)
    td_class = np.ma.compressed(td_class)   
      
    #calculate and store averages
    meanls=np.mean(ls_class)
    meantd=np.mean(td_class)
    meanls_list.append(meanls)
    meantd_list.append(meantd)
    
    #moments
    moments = moment(td,moment=[1,2,3,4])  #Calculate the nth moment about the mean for a sample. Attention: this is already a central moment, 1st moment==0
    #if ordinary moments are wanted, they must be calculated by own function: np.mean(x^n)

    
    #plot all the data
    #ax.scatter(ls_class, td_class, lw=0, alpha=.2)  # Data
    ax.plot(meanls,meantd,'o',markersize=7,markeredgecolor='k',color='g')
    
    ax.plot(meanls,moments[1],'o',markersize=7,markeredgecolor='k',color='g')
    ax.plot(meanls,moments[2],'o',markersize=7,markeredgecolor='k',color='g')
    ax.plot(meanls,moments[3],'o',markersize=7,markeredgecolor='k',color='g')
    
    mom2_list.append(moments[1])
    mom3_list.append(moments[2])
    mom4_list.append(moments[3])

a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,meantd_list)
#dummy x data for plotting
x = np.arange(min(meanls_list), max(meanls_list), 1)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='orange',ls=':')

a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,mom2_list)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='salmon',ls=':')
    
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,mom3_list)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='royalblue',ls=':')
    
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,mom4_list)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='darkred',ls=':')
    

#SAR data
inpath = '../output/def_full/'
outpath = '../plots/'
fname_start = 'td_leg1_L'
lsc_list = [1,2,3,5,7,10,15,25,40]
#lsc_list = [1,2,3,5,7,10,15,25,40,60]  #large steps dont have enough triangles to be representative
minlen = [0,.2,.3,.5,.7,1,1.5,2.5,4]
maxlen = [.15,.3,.4,.6,.9,1.5,2,3.5,5]

meanls_list=[]
meantd_list=[]
mom2_list=[]
mom3_list=[]
mom4_list=[]

for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print(scale)
    fname = inpath+fname_start+str(scale)+'_15km.csv'
    print(fname)
    
    ls = getColumn(fname,1, delimiter=',')
    tls = getColumn(fname,2, delimiter=',')
    td = getColumn(fname,3, delimiter=',')
    ang = getColumn(fname,4)
    ls = np.array(ls,dtype=np.float)/1000  #convert from m to km
    tls = np.array(tls,dtype=np.float)/60/60    #convert from s to h
    td = np.array(td,dtype=np.float)
    ang = np.array(ang,dtype=np.float)
    
    #get only 24h data
    mask = (tls<22) | (tls>25)
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ang = np.ma.array(ang,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)
    ang = np.ma.compressed(ang)
    
    #mask out all the acute triangles
    mask = ang<15
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)
        
    #throw away very low deformation rates (pixel/template edge noise)
    mask = td<.5e-7
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
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
    meantd=np.mean(td); print(meantd)
    meanls_list.append(meanls)
    meantd_list.append(meantd)
    
    #moments
    moments = moment(td,moment=[1,2,3,4])  #Calculate the nth moment about the mean for a sample. Attention: this is already a central moment, 1st moment==0
    #if ordinary moments are wanted, they must be calculated by own function: np.mean(x^n)

    
    #plot all the data
    #ax.scatter(ls_class, td_class, lw=0, alpha=.2)  # Data
    ax.plot(meanls,meantd,'*',markersize=9,markeredgecolor='k',color='purple')
    
    ax.plot(meanls,moments[1],'*',markersize=9,markeredgecolor='k',color='purple')
    ax.plot(meanls,moments[2],'*',markersize=9,markeredgecolor='k',color='purple')
    ax.plot(meanls,moments[3],'*',markersize=9,markeredgecolor='k',color='purple')
    
    mom2_list.append(moments[1])
    mom3_list.append(moments[2])
    mom4_list.append(moments[3])
    

slopes_sar=[]
    
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,meantd_list)
#dummy x data for plotting
x = np.arange(min(meanls_list), max(meanls_list), 1)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='orange',ls='--')
slopes_sar.append(k*-1)

a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,mom2_list)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='salmon',ls='--')
slopes_sar.append(k*-1)
    
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,mom3_list)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='royalblue',ls='--')
slopes_sar.append(k*-1)
    
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,mom4_list)
ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='darkred',ls='--')
slopes_sar.append(k*-1)

ax.grid('on')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())    
#ax.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)

#quadratic fit of the slopes
bx = fig1.add_subplot(122)

moms = np.arange(1,5)
bx.plot(moms,slopes_sar,'o', label='SAR')
bx.plot(moms,slopes_rad,'o', label='ship radar')

z1 = np.polyfit(moms, slopes_sar, 2)
z2 = np.polyfit(moms, slopes_rad, 2)
print(z1)

x=np.linspace(0,5,100)
y1=z1[0]*x**2+z1[1]*x+z1[2]
y2=z2[0]*x**2+z2[1]*x+z2[2]

bx.plot(x,y1,'-')
bx.plot(x,y2,'-')

bx.legend()


fig1.tight_layout()

fig1.savefig(outpath+'power_law_24h_'+title)
