from sid_func import *
import matplotlib.pyplot as plt

#plotting
fig1 = plt.figure(figsize=(9,9))
fig1 = plt.figure(figsize=(7.5,7))
ax = fig1.add_subplot(111)
#title = 'ship_radar+buoys+SAR_faithful'
#title = 'ship_radar+buoys+SAR_bold2'
title = 'ship_radar+buoys+SAR_filter'
ax.set_title(title,fontsize=29, loc='left')
ax.set_xlabel(r"Time scale (h)",fontsize=25)
ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
ax.set_xscale('log')
ax.set_yscale('log')

meanls_list=[]

meantd_list_sr=[]
meantd_list_b=[]

#ship radar data
inpath = '../sidrift/data/ship_radar/time_sliced_data/'
name1 = 'Annu Oikkonen - Period'
name2 = '_Deformation_L6_'
lsc_list = [1,3,6,12,24]

for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print(scale)

    #period 1
    fname = inpath+name1+'1'+name2+str(scale)+'h.txt'
    td = getColumn(fname,1, delimiter=' ')

    ls = scale
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
inpath = '../sidrift/data/buoys/'
fname_start = 'nice1_in_SAR'
lsc_list = [1,3,6,12,24]

for i in lsc_list:
    scale = i
    print(scale)
    
    fname_td1 = inpath+'dr_'+fname_start+'1_'+str(scale)+'h'
    fname_ls1 = inpath+'ls_'+fname_start+'1_'+str(scale)+'h'
    
    td = np.load(fname_td1,encoding='latin1',allow_pickle=True)/24/60/60      #convert from days to s
    ls = np.load(fname_ls1,encoding='latin1',allow_pickle=True)
    
    print(len(ls))
    
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
    ls_class = np.ma.array(ls,mask=mask)
    td_class = np.ma.array(td,mask=mask)
    ls_class = np.ma.compressed(ls_class)
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

#SAR data
reg = 'leg1'
lscale = 26
inpath = '../sidrift/data/80m_stp10_time/'
outpath = inpath
fname_start = 'td_seed_f_Lance_L'

minlen=0
maxlen=10

tc_min = [1,2,15,20,25,50]
tc_max = [2,9,20,25,50,80]


#tc_min = [1,4,6,15,20,25,50]
#tc_max = [2,6,7,20,25,50,80]

#tc_min = [20,25,50]
#tc_max = [25,50,80]


meanls_list_sar=[]
meantd_list_sar=[]



f = inpath+fname_start+str(lscale)+'_7kmFW.csv'
#f = inpath+'td_'+reg+'_L'+str(lscale)+'_25km.csv'
#print(f)

#extract time lengh scale from the name of the file
sls = getColumn(f,1, delimiter=',')
ls = getColumn(f,2, delimiter=',')
td = getColumn(f,3, delimiter=',')
ang = getColumn(f,4)
sls = np.array(sls,dtype=np.float)
ls = np.array(ls,dtype=np.float)
td = np.array(td,dtype=np.float)
ang = np.array(ang,dtype=np.float)



#mask out all the acute triangles
mask=ang>15
ls = ls[mask]/60/60  #convert from s to h
td = td[mask]
sls = sls[mask]/1000  #convert from m to km

#np.set_printoptions(threshold=np.nan)
#print(ls)
#exit()

##mask outliers
#mask = td<1
#ls = ls[mask]
#td = td[mask]
#sls = sls[mask]

##mask all very small or big triangles
##if not masked the renge of the ls is big and has several clouds (expected ls, twice the ls and all kinds of smaller ls)
#mask = (sls>minlen) & (sls<maxlen)
#ls = ls[mask]
#td = td[mask]
#sls = sls[mask]

#print(ls)
#print(sls)
#exit()


#plot all the data
ax.scatter(ls, td, lw=0, alpha=.2)  # Data

#cluster triangles by time difference
#some clusters will be empty!
for i in range(0,len(tc_min)):
    mask = (ls>tc_min[i]) & (ls<tc_max[i])
    lsm = ls[mask]
    tdm = td[mask]
    
    #calculate and store averages
    meanls=np.mean(lsm)
    meantd=np.mean(tdm)
    meanls_list_sar.append(meanls)
    meantd_list_sar.append(meantd)
    
    #plot means
    ax.plot(meanls,meantd,'*',markersize=10,markeredgecolor='k', color='r')
    
    print(meanls)
    print(meantd)


#dummy x data for plotting
tc_min = [.4,1,4,6,7,15,20,25,50]
tc_max = [1,2,6,7,9,20,25,50,80]

x = (np.array(tc_min)+np.array(tc_max))/2




#ship radar
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list,meantd_list_sr)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='m')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#buoys
#fit the line to bins
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list[1:],meantd_list_b[1:])

ax.loglog(x[1:],a*x[1:]**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='g')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

##sar 
##fit the line
#a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar,meantd_list_sar)



ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f L^{%.2f}$' %(a,k),c='royalblue')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)



ax.grid('on')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())    
ax.legend(loc='upper right',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig1.tight_layout()

fig1.savefig(outpath+'power_law_3-5km_'+title)
