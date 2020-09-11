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

#SYI
inpath = '../sidrift/data/40m_combo/'
inpath = '../sidrift/data/80m_stp10/'
inpath = '../sidrift/data/80m_stp10_single_filter/'

outpath = inpath
fname_start = 'td_leg1_SYI_L'
fname_start = 'td_seed_f_SYI_L'

#create log-spaced vector and convert it to integers
n=8 # number of samples
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)
#size envelope also needs to increase (from 10m to 3km)
margin = np.exp(np.linspace(np.log(.15),np.log(3),n))


for i in range(0,len(stp)-4):                           #the last two steps are off the curve, try removing them
    scale = stp[i]
    fname = inpath+fname_start+str(scale)+'_7km_n9.csv'
    #fname = inpath+fname_start+str(scale)+'_25kmFW.csv'
    #fname = inpath+fname_start+str(scale)+'_20kmFW.csv'
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
        
    ##throw away very low deformation rates (pixel edge noise)
    #mask = td<.5e-7
    #ls = np.ma.array(ls,mask=mask)
    #td = np.ma.array(td,mask=mask)
    #ls = np.ma.compressed(ls)
    #td = np.ma.compressed(td)   

    ###throw away very high deformation rates (image edges and other artifacts)
    #mask = td>1e-4
    #ls = np.ma.array(ls,mask=mask)
    #td = np.ma.array(td,mask=mask)
    #ls = np.ma.compressed(ls)
    #td = np.ma.compressed(td)   


    #mask all very small or big triangles
    #if not masked the range of the ls is big and has several clouds (expected ls, twice the ls and all kinds of smaller ls)
    center = np.mean(ls)
    minlen = center-margin[i]; maxlen = center+margin[i]
    mask = (ls<minlen) | (ls>maxlen)
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
    ax.scatter(ls, td, lw=0, alpha=.5, color='darkred')  # Data
    ax.plot(meanls,meantd,'*',markersize=10,markeredgecolor='w',color='brown')

#FYI
fname_start = 'td_leg1_FYI_L'
fname_start = 'td_seed_f_FYI_L'


for i in range(0,len(stp)-4):                           #the last two steps are off the curve, try removing them
    scale = stp[i]
    fname = inpath+fname_start+str(scale)+'_7km_n9.csv'
    #fname = inpath+fname_start+str(scale)+'_25kmFW.csv'
    #fname = inpath+fname_start+str(scale)+'_20kmFW.csv'
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
        
    ##throw away very low deformation rates (template noise)
    #mask = td<.5e-7
    #ls = np.ma.array(ls,mask=mask)
    #td = np.ma.array(td,mask=mask)
    #ls = np.ma.compressed(ls)
    #td = np.ma.compressed(td)   
    
    #throw away very high deformation rates (image edges and other artifacts)
    mask = td>1e-4
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)   


    #mask all very small or big triangles
    #if not masked the range of the ls is big and has several clouds (expected ls, twice the ls and all kinds of smaller ls)
    center = np.mean(ls)
    minlen = center-margin[i]; maxlen = center+margin[i]
    mask = (ls<minlen) | (ls>maxlen)
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
    ax.scatter(ls, td, lw=0, alpha=.5, color='gold')  # Data
    ax.plot(meanls,meantd,'*',markersize=10,markeredgecolor='w',color='darkorange')



#SYI
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar,meantd_list_sar)

#dummy x data for plotting
x = np.arange(min(meanls_list_sar), max(meanls_list_sar), 1)
x = np.arange(min(ls_list_sar), max(ls_list_sar), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$ (SYI)' %(a*10e6,k),c='darkred')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#FYI
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_fyi,meantd_list_fyi)

#dummy x data for plotting
x = np.arange(min(meanls_list_fyi), max(meanls_list_fyi), 1)
x = np.arange(min(ls_list_fyi), max(ls_list_fyi), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$ (FYI)' %(a*10e6,k),c='darkorange')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

ax.grid(True)
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())    
ax.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig1.tight_layout()

fig1.savefig(outpath+'power_law_24h_it_7km_seed_f_masked_n9'+title)
#fig1.savefig(outpath+'power_law_24h_it_7km_'+title)
