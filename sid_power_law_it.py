from sid_func import *
import matplotlib.pyplot as plt

#plotting
fig1 = plt.figure(figsize=(9,9))
fig1 = plt.figure(figsize=(7.5,7))
ax = fig1.add_subplot(111)
title = 'SAR'
ax.set_title(title,fontsize=29, loc='left')
ax.set_xlabel(r"Length scale (km)",fontsize=25)
ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
ax.set_xscale('log')
ax.set_yscale('log')

meanls_list_sar=[]
meantd_list_sar=[]
meanls_list_fyi=[]
meantd_list_fyi=[]

ls_list_sar=[]
td_list_sar=[]
ls_list_fyi=[]
td_list_fyi=[]

#SAR at Lance
inpath = '../output/def_'
outpath = '../plots/'
fname_start = 'td_leg1_L'
lsc_list = [25,50,100,200,500,1000]   #not ls but number of nominal grid points
minlen = [4,2,1,.5,.2,0]
maxlen = [6,3,1.5,.75,.3,.15]

#lsc_list = [50,100,200,500,1000]   #not ls but number of nominal grid points
#minlen = [2,1,.5,.2,0]
#maxlen = [3,1.5,.75,.3,.15]


for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print(scale)
    fname = inpath+str(scale)+'/'+fname_start+str(scale)+'_15km.csv'
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
    
    #throw away very high deformation rates (unrealistic values)
    mask = td>10e-3
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)   
    
    
    #throw away very low deformation rates (noise)
    mask = td<10e-8
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
    meanls_list_sar.append(meanls)
    meantd_list_sar.append(meantd)
    
    ls_list_sar.extend(ls)
    td_list_sar.extend(td)
        
    #plot all the data
    ax.scatter(ls, td, lw=0, alpha=.2, color='k')  # Data
    ax.plot(meanls,meantd,'*',markersize=10,markeredgecolor='w',color='k')

#SAR in FYI
inpath = '../output/def_'
outpath = '../plots/'
fname_start = 'td_leg1_FYI_L'
lsc_list = [25,50,100,200,500,1000]   #not ls but number of nominal grid points
minlen = [4,2,1,.5,.2,0]
maxlen = [6,3,1.5,.75,.3,.15]

#lsc_list = [50,100,200,500,1000]   #not ls but number of nominal grid points
#minlen = [2,1,.5,.2,0]
#maxlen = [3,1.5,.75,.3,.15]



for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print(scale)
    fname = inpath+str(scale)+'/'+fname_start+str(scale)+'_15km.csv'
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
    
    #throw away very high deformation rates (unrealistic values)
    mask = td>10e-3
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)   
    
    
    #throw away very low deformation rates (noise)
    mask = td<10e-8
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
    meanls_list_fyi.append(meanls)
    meantd_list_fyi.append(meantd)
    
    ls_list_fyi.extend(ls)
    td_list_fyi.extend(td)
        
    #plot all the data
    ax.scatter(ls, td, lw=0, alpha=.2, color='r')  # Data
    ax.plot(meanls,meantd,'*',markersize=10,markeredgecolor='w',color='r')



#Lance
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar,meantd_list_sar)
#no binning
#a,k,cix,ciy_upp,ciy_low = logfit(ls_list_sar,td_list_sar)

#dummy x data for plotting
x = np.arange(min(meanls_list_sar), max(meanls_list_sar), 1)
x = np.arange(min(ls_list_sar), max(ls_list_sar), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='purple')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#FYI
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_fyi,meantd_list_fyi)
#no binning
#a,k,cix,ciy_upp,ciy_low = logfit(ls_list_sar,td_list_sar)

#dummy x data for plotting
x = np.arange(min(meanls_list_fyi), max(meanls_list_fyi), 1)
x = np.arange(min(ls_list_fyi), max(ls_list_fyi), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='royalblue')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)






ax.grid('on')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())    
ax.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig1.tight_layout()

fig1.savefig(outpath+'power_law_24h_it_'+title)
