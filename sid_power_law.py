from sid_func import *
from scipy import stats
import matplotlib.pyplot as plt

#plotting
fig1 = plt.figure(figsize=(9,9))
fig1 = plt.figure(figsize=(7.5,7))
ax = fig1.add_subplot(111)
title = 'ship_radar+buoys+SAR_UiT_seed_f_test'
#title = 'ship_radar+buoys+SAR_UiT_seed_nofilter_test'
#radius = '_7km_hes9FW.csv'
#radius = '_7km_hes9.csv'
#radius = '_7km_ttdFW.csv'
#radius = '_7km_n9.csv'
#radius = '_20kmFW.csv'
#radius = '_20km_hes9.csv'
#radius = '_20km_hes9_1km.csv'
#radius = '_20km_testFW.csv'
#radius = '_50kmFW.csv'
#radius = '_20km.csv'
#radius = '_100kmFW.csv'
#radius = '_100km_hes9_1km.csv'
radius = '_100km_ttd.csv'
#radius = '_100km_n9.csv'
#radius = '_120km.csv'
rr = radius.split('.')[0]
name = 'ship_radar+buoys+SAR'
ax.set_title(name,fontsize=29, loc='left')
ax.set_xlabel(r"Length scale (km)",fontsize=25)
ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
#ax.set_xscale('log')
#ax.set_yscale('log')

#ax.set_xlim(.1,10)

meanls_list_sr=[]
meantd_list_sr=[]
meanls_list_sar=[]
meantd_list_sar=[]
meanls_list_b=[]
meantd_list_b=[]

ls_list_sr=[]
td_list_sr=[]
ls_list_b=[]
td_list_b=[]
ls_list_sar=[]
td_list_sar=[]


#SAR data
#inpath = '../sidrift/data/whole_series_10stp_factor_def/data/'
#outpath = '../sidrift/data/whole_series_10stp_factor_def/plots/'

inpath = '../sidrift/data/40m_combo/'
inpath = '../sidrift/data/80m_stp10/'
inpath = '../sidrift/data/80m_stp10_canberra/'
inpath = '../sidrift/data/80m_stp10_nofilter/'
inpath = '../sidrift/data/80m_stp10_single_filter/'
outpath = inpath

fname_start = 'td_leg1_L'
fname_start = 'td_seed_leg1_L'
fname_start = 'td_seed_f_Lance_L'
#lsc_list = [1,2,3,5,7,10,15,25,40]
#lsc_list = [1,2,3,5,7,10,15,25,40,60]  #large steps dont have enough triangles to be representative
#minlen = [0,  .2,.3,.5,.7,1,  1.5,2.5,0]
#maxlen = [.15,.3,.4,.6,.9,1.5,2,  3.5,20]

#create log-spaced vector and convert it to integers
n=8 # number of samples
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)
print(stp)

##for 50km radius
#n=9
#stp=np.exp(np.linspace(np.log(1),np.log(800),n))
#stp = stp.astype(int)
#margin = np.exp(np.linspace(np.log(.2),np.log(3),n))


##recalculate these lenghts to sqrt of area (lenght scale)
##assume an average triangle is unilateral
##change unit from step to 40m (.04 km) pixel size
#ls_stp = np.sqrt(np.sqrt(3)/4*stp**2*.04)
#print(ls_stp)

#size envelope also needs to increase (from 10m to 3km)
margin = np.exp(np.linspace(np.log(.2),np.log(10),n))

#colors
color=iter(plt.cm.Blues_r(np.linspace(0,1,len(stp)+1)))

#put dummy values (detection limits) on the scatter plots
fname = inpath+'dummy_Lance'+radius
print(fname)

dum1 = getColumn(fname,0, delimiter=',', header=False)
dum2 = getColumn(fname,1, delimiter=',', header=False)
dum1 = np.array(dum1,dtype=np.float)/1000  #convert from m to km
dum2 = np.array(dum2,dtype=np.float)

ax.plot(dum1, dum2, 'x', color='k', alpha=.2)

##for the LKF plot
#fig2 = plt.figure(figsize=(9,9))
#bx = fig2.add_subplot(111)

for i in range(0,len(stp)-0):                           
#for i in range(0,len(stp)):    
    scale = stp[i]
    print(scale)
    fname = inpath+fname_start+str(scale)+radius
    print(fname)
    
    ls = getColumn(fname,1, delimiter=',')
    tls = getColumn(fname,2, delimiter=',')
    td = getColumn(fname,3, delimiter=',')
    ang = getColumn(fname,4)
    ls = np.array(ls,dtype=np.float)/1000  #convert from m to km
    tls = np.array(tls,dtype=np.float)/60/60    #convert from s to h
    td = np.array(td,dtype=np.float)
    ang = np.array(ang,dtype=np.float)
    
    #lkf_ids = getColumn(fname,5)
    #lkf_ids = np.array(lkf_ids,dtype=np.float)
    #print(lkf_ids)
        
    #get only 24h data
    mask = (tls<23) | (tls>25)
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
    #only for smallest scales, where we have underdetected deformation
    mask = (td<dum2[i]) #& (ls < 1)
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)  
    
    ###throw away very high deformation rates (image edges and other artifacts)
    #mask = td>1.5e-4
    #ls = np.ma.array(ls,mask=mask)
    #td = np.ma.array(td,mask=mask)
    #ls = np.ma.compressed(ls)
    #td = np.ma.compressed(td)   


    #mask all very small or big triangles
    #if not masked the range of the ls is big and has several clouds (expected ls, twice the ls and all kinds of smaller ls)
    center = np.mean(ls)
    #center = stats.mode(ls)[0][0]                      #this takes too much time
    print('mean lenght')
    print(center)
    print(margin[i])
    minlen = center-margin[i]; maxlen = center+margin[i]
    #minlen = center-margin; maxlen = center+margin
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
    cl = next(color)
    ax.scatter(ls, td, lw=0, alpha=.1, color=cl)  # Data
    ax.plot(meanls,meantd,'*',markersize=10,markeredgecolor='w',color=cl)


    ##plot example threshold
    #ax.axhline(y=2.1040826060411214e-06)
    ##why are there values bellow this line at L=1 ??????????????????????
    ##low values due to shear...
    #ax.axhline(y=4.0344082184489904e-07)
    #ax.axhline(y=7.827424274136162e-09)



    ##plot a pdf of lkf id lenght (len1=only no LKF triangles, this should be masked actually at all ls)
    ##this will not really work as some LKFs are sometimes connected and sometimes not...
    #if i > 0:   #first step is special data
        #weights = np.ones_like(lkf_ids) / (len(lkf_ids))
        #n, bins, patches = bx.hist(lkf_ids, range(1,7),label=str(stp[i]), density=True, alpha=1, histtype='step',lw=5, weights=weights)

#bx.legend()
#fig2.savefig(outpath+'lkf_ids_len'+title+rr)

#ship radar data
#inpath = '../data/ship_radar/24h/'
inpath = '../sidrift/data/ship_radar/time_sliced_data/'
name1 = 'Annu Oikkonen - Period'
name2 = '_Deformation_L'

lsc_list = range(1,7)

#colors
color=iter(plt.cm.Purples_r(np.linspace(0,1,len(lsc_list)+1)))

#dummy values for ship radar
dum2r = [dum2[0],dum2[1],dum2[1],dum2[2],dum2[3],dum2[3]]


for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print(scale)
    
    #period 1
    fname = inpath+name1+'1'+name2+str(scale)+'_24h.txt'
    ls1 = getColumn(fname,0, delimiter=' ')
    td1 = getColumn(fname,1, delimiter=' ')
    
    #period 2
    fname = inpath+name1+'2'+name2+str(scale)+'_24h.txt'
    ls2 = getColumn(fname,0, delimiter=' ')
    td2 = getColumn(fname,1, delimiter=' ')
    
    #combine
    ls = np.array(np.append(ls1,ls2),dtype=np.float)/1000  #convert from m to km
    td = np.array(np.append(td1,td2),dtype=np.float)/3600  #convert from h to s   
    
    ##use only first part
    #ls = np.array(ls1,dtype=np.float)/1000  #convert from m to km
    #td = np.array(td1,dtype=np.float)/3600  #convert from h to s 
    
    #mask with the dummy values from SAR
    mask = (td<dum2r[i])
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)
    
    #calculate and store averages
    meanls=np.mean(ls)
    meantd=np.mean(td)
    meanls_list_sr.append(meanls)
    meantd_list_sr.append(meantd)
    ls_list_sr.extend(ls)
    td_list_sr.extend(td)

    
    #plot all the data
    cl = next(color)
    ax.scatter(ls, td, marker='s', lw=0, alpha=.2, color=cl)  # Data
    ax.plot(meanls,meantd,'s',markersize=7,markeredgecolor='orange', color=cl)

#buoy data - do we have enough of 1 day data (should be enough for the entire leg 1)
#scales 2-100km
inpath = '../sidrift/data/buoys/'
fname_start = 'nice1_comb_SAR'
fname_start = 'nice1_in_SAR'        #only inner ring of the buoys
#fname_start = 'nice1_short_SAR'    
fname_td1 = inpath+'dr_'+fname_start+'1_'+'24h'
fname_ls1 = inpath+'ls_'+fname_start+'1_'+'24h'
fname_td2 = inpath+'dr_'+fname_start+'2_'+'24h'
fname_ls2 = inpath+'ls_'+fname_start+'2_'+'24h'
 
td = np.append(np.load(fname_td1,encoding='latin1',allow_pickle=True),np.load(fname_td2,encoding='latin1',allow_pickle=True))/24/60/60      #convert from days to s
ls = np.append(np.load(fname_ls1,encoding='latin1',allow_pickle=True),np.load(fname_ls2,encoding='latin1',allow_pickle=True))

##use only first part
#td = np.load(fname_td1,encoding='latin1',allow_pickle=True)/24/60/60      #convert from days to s
#ls = np.load(fname_ls1,encoding='latin1',allow_pickle=True)



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

#colors
color=iter(plt.cm.Greens_r(np.linspace(0,1,len(lsc_list)+1)))

#dummy values for buoys
dum2b = [dum2[3],dum2[4],dum2[4]]


for i in range(0,len(lsc_list)):
    print(i)
    mask = (ls<minlen[i]) | (ls>maxlen[i])
    ls_class = np.ma.array(ls,mask=mask)
    td_class = np.ma.array(td,mask=mask)
    ls_class = np.ma.compressed(ls_class)
    td_class = np.ma.compressed(td_class) 
    
    #mask with the dummy values from SAR
    mask = (td_class<dum2b[i])
    ls_class = np.ma.array(ls_class,mask=mask)
    td_class = np.ma.array(td_class,mask=mask)
    ls_class = np.ma.compressed(ls_class)
    td_class = np.ma.compressed(td_class)

      
    #calculate and store averages
    meanls=np.mean(ls_class)
    meantd=np.mean(td_class)
    meanls_list_b.append(meanls)
    meantd_list_b.append(meantd)
    ls_list_b.extend(ls)
    td_list_b.extend(td)

    
    #plot all the data
    cl = next(color)
    ax.scatter(ls_class, td_class, marker='o', lw=0, alpha=.2, color=cl)  # Data
    ax.plot(meanls,meantd,'o',markersize=7,markeredgecolor='yellow', color=cl)


#ship radar
#fit the line
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sr,meantd_list_sr)
#a,k,cix,ciy_upp,ciy_low = logfit(ls_list_sr,td_list_sr)

#dummy x data for plotting
x = np.arange(min(meanls_list_sr), max(meanls_list_sr), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='orange')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#buoys
#fit the line to bins
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_b,meantd_list_b)
#a,k,cix,ciy_upp,ciy_low = logfit(ls_list_b,td_list_b)

#dummy x data for plottingoutpath+'power_law_24h_'+title+rr
x = np.arange(min(meanls_list_b), max(meanls_list_b), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='y')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#and separate for SAR
#fit the line
#a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar,meantd_list_sar)
#no binning
#a,k,cix,ciy_upp,ciy_low = logfit(ls_list_sar,td_list_sar)
#skip the last 2 lenghts
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar[:-4],meantd_list_sar[:-4])

#dummy x data for plotting
x = np.arange(min(meanls_list_sar), max(meanls_list_sar), 1)
x = np.arange(min(ls_list_sar), max(ls_list_sar), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='darkred')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

ax.grid('on')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())    
ax.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig1.tight_layout()

fig1.savefig(outpath+'power_law_24h_'+title+rr)
print(outpath+'power_law_24h_'+title+rr)


#some extra LKF statistics
