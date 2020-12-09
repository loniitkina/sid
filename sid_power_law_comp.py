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
#radius = '_7km_test5.csv'
#radius = '_20km.csv'
#radius = '_20km_hes9.csv'
#radius = '_20km_hes9_1km.csv'
#radius = '_20km_test.csv'
#radius = '_20km_test6.csv'
#radius = '_50kmFW.csv'
#radius = '_20km.csv'
#radius = '_100kmFW.csv'
#radius = '_100km_hes9_1km.csv'
#radius = '_100km_ttd.csv'
radius = '_100km_n9.csv'
#radius = '_120km.csv'
rr = radius.split('.')[0]
name = 'ship_radar+buoys+SAR'
ax.set_title(name,fontsize=29, loc='left')
ax.set_xlabel(r"Length scale (km)",fontsize=25)
ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
#ax.set_xscale('log')
#ax.set_yscale('log')

#ax.set_xlim(.1,10)


#scale total divergence to average displacement (in one direction)
#dux = n * sigma_x**2 / (2* A * dt).
fig3 = plt.figure(figsize=(7.5,7))
cx = fig3.add_subplot(111)

cx.set_xlabel(r"Length scale (km)",fontsize=25)
cx.set_ylabel(r"Total displacement (m)",fontsize=25)
cx.set_xscale('log')
cx.set_yscale('log')
cx.set_xlim(.001,100)



meanls_list_sr=[]
meantd_list_sr=[]
meanls_list_sr_all=[]
meantd_list_sr_all=[]


meanls_list_sar=[]
meantd_list_sar=[]
stdtd_list_sar=[]

meanls_list_b=[]
meantd_list_b=[]

ls_list_sr=[]
td_list_sr=[]
ls_list_b=[]
td_list_b=[]
ls_list_sar=[]
td_list_sar=[]

mean_dist_sar=[]
std_dist_sar=[]
mean_dist_sr=[]
std_dist_sr=[]

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
n=9 # number of samples
n=8
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
margin = np.exp(np.linspace(np.log(.25),np.log(10),n))

#colors
color=iter(plt.cm.Blues_r(np.linspace(0,1,len(stp)+1)))

#put dummy values (detection limits) on the scatter plots
#radius = '_7km_test1.csv'               #use dummy values from different run
fname = inpath+'dummy_Lance'+radius
print(fname)

dum1 = getColumn(fname,0, delimiter=',', header=False)
dum2 = getColumn(fname,1, delimiter=',', header=False)
dum1 = np.array(dum1,dtype=np.float)/1000  #convert from m to km
dum2 = np.array(dum2,dtype=np.float)

ax.plot(dum1, dum2, 'x', color='k', alpha=.2)

#get a power law fit on these dummy values
a_dum,k_dum,cix,ciy_upp,ciy_low = logfit(dum1[:6],dum2[:6])
x = np.arange(0.1, 100, 1)
ax.loglog(x,a_dum*x**k_dum,linewidth=1,c='k')



##for the LKF plot
#fig2 = plt.figure(figsize=(9,9))
#bx = fig2.add_subplot(111)

print('***********************sar***************************')

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
    tls = np.ma.array(tls,mask=mask); tls = np.ma.compressed(tls)
    
    #mask out all the acute triangles
    mask = ang<15
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)
    tls = np.ma.array(tls,mask=mask); tls = np.ma.compressed(tls)
        
    ##throw away very low deformation rates (pixel/template edge noise)
    no_all = len(td)
    #mask = (td<dum2[i])                 #this is only from the first image pair!
    mask = td < a_dum*ls**k_dum
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)
    tls = np.ma.array(tls,mask=mask); tls = np.ma.compressed(tls)
    no_part = len(td)
    fraction = no_part/no_all
    print('fraction: '+str(fraction)+'from: '+str(no_all))
    
    ##replace low deformation rates by averages
    #mask = (td<dum2[i])
    #fv = np.mean(np.ma.array(td,mask=~mask))
    #print('replacing bellow', str(dum2[i]), 'by', str(fv) )
    #td[mask] = fv
    
    
    
    
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
    tls = np.ma.array(tls,mask=mask); tls = np.ma.compressed(tls)
    
    #calculate and store averages
    meanls=np.mean(ls)
    meantd=np.mean(td); print('mean td :',meantd, 'from', td.shape[0], 'values')
    meanls_list_sar.append(meanls)
    meantd_list_sar.append(meantd)
    stdtd_list_sar.append(np.std(td))
    
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


    #convert td to displacement
    #dux = n * sigma_x**2 / (2* A * dt)
    #sigma_x = np.sqrt(dux * (2* A * dt) /n )
    tls = tls *3600 #convert back to seconds
    ls = ls *1000 #convert back to meters
    sigma_x = np.sqrt(td * (2* ls**2 * tls) /3 )
    sigma_x_mean = np.mean(sigma_x)
    mean_dist_sar.append(sigma_x_mean)
    std_dist_sar.append(np.std(sigma_x))

    cx.scatter(ls/1000, sigma_x, lw=0, alpha=.1, color=cl)
    cx.plot(meanls,sigma_x_mean,'*',markersize=10,markeredgecolor='w',color=cl)

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

print('***********************ship_radar***************************')
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
    
    #calculate and store averages
    meanls=np.mean(ls)
    meantd=np.mean(td)
    meanls_list_sr_all.append(meanls)
    meantd_list_sr_all.append(meantd)
    #ls_list_sr.extend(ls)
    #td_list_sr.extend(td)
    
    #plot averages for all the data
    ax.plot(meanls,meantd,'s',markersize=7,markeredgecolor='orange', color=cl)
    #all data in grey
    ax.scatter(ls, td, marker='s', lw=0, alpha=.2, color='0.8')
    
    #convert td to displacement
    tls = 24 *3600 #convert back to seconds
    lsc = ls *1000 #convert back to meters
    sigma_x = np.sqrt(td * (2* lsc**2 * tls) /3 )
    sigma_x_mean = np.mean(sigma_x)
    cx.scatter(ls, sigma_x, lw=0, alpha=.1, color=cl)
    cx.plot(meanls,sigma_x_mean,'*',markersize=10,markeredgecolor='w',color=cl)
    
    mean_dist_sr.append(sigma_x_mean)
    std_dist_sr.append(np.std(sigma_x))
    
    #mask with the dummy values from SAR
    no_all = len(td)
    #mask = (td<dum2r[i])
    mask = td < a_dum*ls**k_dum
    ls = np.ma.array(ls,mask=mask)
    td = np.ma.array(td,mask=mask)
    ls = np.ma.compressed(ls)
    td = np.ma.compressed(td)
    no_part = len(td)
    fraction = no_part/no_all
    print('fraction: '+str(fraction)+'from: '+str(no_all))
    
    #calculate and store averages
    meanls=np.mean(ls)
    meantd=np.mean(td)
    meanls_list_sr.append(meanls)
    meantd_list_sr.append(meantd)
    #ls_list_sr.extend(ls)
    #td_list_sr.extend(td)

    
    #plot all the data
    cl = next(color)
    ax.scatter(ls, td, marker='s', lw=0, alpha=.2, color=cl)  # Data
    #ax.plot(meanls,meantd,'s',markersize=7,markeredgecolor='orange', color=cl)
    
    


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

print('***********************buoys***************************')
for i in range(0,len(lsc_list)):
    print(i)
    mask = (ls<minlen[i]) | (ls>maxlen[i])
    ls_class = np.ma.array(ls,mask=mask)
    td_class = np.ma.array(td,mask=mask)
    ls_class = np.ma.compressed(ls_class)
    td_class = np.ma.compressed(td_class) 
    
    #mask with the dummy values from SAR
    #mask = (td_class<dum2b[i])
    mask = td_class < a_dum*ls_class**k_dum
    ls_class = np.ma.array(ls_class,mask=mask)
    td_class = np.ma.array(td_class,mask=mask)
    ls_class = np.ma.compressed(ls_class)
    td_class = np.ma.compressed(td_class)

      
    #calculate and store averages
    meanls=np.mean(ls_class)
    meantd=np.mean(td_class)
    meanls_list_b.append(meanls)
    meantd_list_b.append(meantd)
    #ls_list_b.extend(ls)
    #td_list_b.extend(td)

    
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

#ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='orange')
#ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
#ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#all the ship radar data
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sr_all[:],meantd_list_sr_all[:])
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
##fit the line
#a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar,meantd_list_sar)
#no binning
#a,k,cix,ciy_upp,ciy_low = logfit(ls_list_sar,td_list_sar)
#skip the last 2 lenghts
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar[:-3],meantd_list_sar[:-3])

#dummy x data for plotting
x = np.arange(min(meanls_list_sar), max(meanls_list_sar), 1)
x = np.arange(min(ls_list_sar), max(ls_list_sar), 1)

ax.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='darkred')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

##and just SAR long scales
#a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar[-5:],meantd_list_sar[-5:])
#ax.loglog(x,a*x**k,linewidth=2, ls='--',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='darkred')

#and combined SR for short scales and SAR for long scales
comb_ls=[]
comb_ls.extend(meanls_list_sr_all[:-2])
comb_ls.extend(meanls_list_sar[-4:])

comb_td=[]
comb_td.extend(meantd_list_sr_all[:-2])
comb_td.extend(meantd_list_sar[-4:])

a,k,cix,ciy_upp,ciy_low = logfit(comb_ls,comb_td)
ax.loglog(x,a*x**k,linewidth=2, ls=':',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='darkred')
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)




#and Ship radar distances
x = np.arange(.001, max(ls_list_sar), 1)

ad,kd,cix,ciy_upp,ciy_low = logfit(meanls_list_sr,mean_dist_sr)
cx.loglog(x,ad*x**kd,linewidth=2, ls='--',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(ad*10e6,kd),c='purple')
#ast,kst,cix,ciy_upp,ciy_low = logfit(meanls_list_sr,std_dist_sr)
#cx.loglog(x,ast*x**kst,linewidth=2, ls=':',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(ast*10e6,kst),c='purple')

#and SAR distances
ad,kd,cix,ciy_upp,ciy_low = logfit(meanls_list_sar[-5:],mean_dist_sar[-5:])
cx.loglog(x,ad*x**kd,linewidth=2, ls='--',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(ad*10e6,kd),c='darkred')
##SAR standard deviation
#ast,kst,cix,ciy_upp,ciy_low = logfit(meanls_list_sar[0:],std_dist_sar[0:])
#cx.loglog(x,ast*x**kst,linewidth=2, ls=':',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(ast*10e6,kst),c='darkred')


#combined
comb_dist=[]
comb_dist.extend(mean_dist_sr[:-2])
comb_dist.extend(mean_dist_sar[-4:])
comb_std=[]
comb_std.extend(std_dist_sr[:-2])
comb_std.extend(std_dist_sar[-4:])

ad,kd,cix,ciy_upp,ciy_low = logfit(comb_ls,comb_dist)
cx.loglog(x,ad*x**kd,linewidth=2, ls=':',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(ad*10e6,kd),c='darkred')
ast,kst,cix,ciy_upp,ciy_low = logfit(comb_ls,comb_std)

#extend the x
dist=[]
dist.extend([.001,.03,.2])
dist.extend(meanls_list_sar)

#create random log normal distributon with such means and standar deviation
def generate_lognormal_samples(mean, stdev, n=1):
    """
    from: https://pythonhealthcare.org/2019/02/07/120-generating-log-normal-samples-from-provided-arithmetic-mean-and-standard-deviation-of-original-population/
    Returns n samples taken from a lognormal distribution, based on mean and
    standard deviation calaculated from the original non-logged population.
    
    Converts mean and standard deviation to underlying lognormal distribution
    mu and sigma based on calculations desribed at:
        https://blogs.sas.com/content/iml/2014/06/04/simulate-lognormal-data-
        with-specified-mean-and-variance.html
        
    Returns a numpy array of floats if n > 1, otherwise return a float
    """
    
    # Calculate mu and sigma of underlying lognormal distribution
    phi = (stdev ** 2 + mean ** 2) ** 0.5
    mu = np.log(mean ** 2 / phi)
    sigma = (np.log(phi ** 2 / mean ** 2)) ** 0.5
    
    # Generate lognormal population
    generated_pop = np.random.lognormal(mu, sigma , n)
    
    # Convert single sample (if n=1) to a float, otherwise leave as array
    generated_pop = \
        generated_pop[0] if len(generated_pop) == 1 else generated_pop
        
    return generated_pop


for i in range(0,len(dist)):
    #distances
    mean = ad*dist[i]**kd
    #stdev = std_dist_sar[i]                                       #here we need to use something else... scale this from the last 4 too!!!
    stdev = ast*dist[i]**kst
    generated_pop = generate_lognormal_samples(mean, stdev, 1000)
    print ('Mean:', generated_pop.mean())
    print ('Standard deviation:', generated_pop.std())

    sample_x = np.ones_like(generated_pop)*dist[i]

    cx.scatter(sample_x, generated_pop, marker='o', lw=0, alpha=.2, color='r')
    cx.plot(dist[i],generated_pop.mean(),'*',markersize=10,markeredgecolor='w',color='r')
    
    ##total deformation
    #mean = a*meanls_list_sar[i]**k
    #stdev = stdtd_list_sar[i]                                       #here we need to use something else...
    
    #generated_pop = generate_lognormal_samples(mean, stdev, 1000)
    #print ('Mean:', generated_pop.mean())
    #print ('Standard deviation:', generated_pop.std())

    #sample_x = np.ones_like(generated_pop)*meanls_list_sar[i]

    #ax.scatter(sample_x, generated_pop, marker='o', lw=0, alpha=.2, color='r')
    #ax.plot(meanls_list_sar[i],generated_pop.mean(),'*',markersize=10,markeredgecolor='w',color='r')


ax.grid('on')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())    
ax.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig1.tight_layout()

fig1.savefig(outpath+'power_law_24h_comp_'+title+rr)
print(outpath+'power_law_24h_comp_'+title+rr)


#storing displacements
cx.grid('on')
cx.xaxis.set_major_formatter(ScalarFormatter())   
#cx.minorticks_off()
#cx.set_yscale('log', subsy=[.01, .1, 1, 10, 100, 1000, 10000]) #basey, 'subsy', nonposy'
#cx.yaxis.set_major_formatter(ScalarFormatter())
cx.set_yticks([.001,.01, .1, 1, 10, 100, 1000, 10000, 100000],[0.001,0.01, 0.1, 1, 10, 100, 1000, 10000, 100000])


#cx.yaxis.ticker.Formatter.IndexFormatter(labels=[0.001,0.01, 0.1, 1, 10, 100, 1000, 10000, 100000])
cx.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig3.tight_layout()

fig3.savefig(outpath+'power_law_24h_disp_'+title+rr)




#difference plots for short time steps

x = np.array(meanls_list_sar)
sar_est = a*x**k


diff_sar = np.array(meantd_list_sar) - sar_est


x = np.array(meanls_list_sr)
sh_est = a*x**k

diff_sr = np.array(meantd_list_sr) - sh_est

diff_sr_all = np.array(meantd_list_sr_all) - sh_est

#print(meanls_list_sr)
#print(meanls_list_sr_all)

print(diff_sar)
print(diff_sr)









fig2 = plt.figure(figsize=(7.5,7))
bx = fig2.add_subplot(111)

bx.set_xlabel(r"Length scale (km)",fontsize=25)
bx.set_ylabel(r"Total deformation difference (s$^{-1}$)",fontsize=25)

bx.loglog(meanls_list_sar,diff_sar,linewidth=2, marker='o')
bx.loglog(meanls_list_sr,diff_sr,linewidth=2, marker='o')
bx.loglog(meanls_list_sr,diff_sr_all,linewidth=2, marker='o')

#dummy line
a_dum,k_dum,cix,ciy_upp,ciy_low = logfit(dum1[:6],dum2[:6])
bx.loglog(x,a_dum*x**k_dum,linewidth=1,c='k',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k))

#SAR
x = np.arange(min(meanls_list_sar[:-4]), max(meanls_list_sar[:-4]), 1)
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar[:-4],diff_sar[:-4])
bx.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='y')
#bx.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
#bx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#ship radar
x = np.arange(min(meanls_list_sr), max(meanls_list_sr), 1)
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sr[:-2],diff_sr[:-2])
bx.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='r')
#bx.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
#bx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#ship radar all
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sr[:-3],diff_sr_all[:-3])
bx.loglog(x,a*x**k,linewidth=2,label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k),c='b')


bx.grid('on')
bx.xaxis.set_major_formatter(ScalarFormatter())    
bx.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig2.tight_layout()

fig2.savefig(outpath+'power_law_24h_comp_diff_'+title+rr)






