import numpy as np
from sid_func import getColumn, logfit
from scipy import stats
from glob import glob
import pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from datetime import datetime

reg = 'lance_leg1'

filtering=False
filtering=True
#This is threshold filtering.
#All data here are stored by sid_deform_tile.py and are tiled and lkf-filtered/kernel-filtered

if filtering:
    naming='filter'
else:
    naming='nofilter'

#file_name_end and output
radius = 200000
rname = '_'+str(int(radius/1000))+'km'

if filtering:
    #threshold_filter:
    tname='_thfilter'
    #WARNING: choose the right option here!
    #lkf_filter:
    lname='_lkffilter'
    #lname='_nolkfilter'
else:
    tname='_nothfilter'
    #WARNING: choose the right option here!
    #lkf_filter:
    lname='_lkffilter'
    lname='_nolkfilter'

file_name_end = rname+tname+lname

print('Doing: ',file_name_end, reg)

powerlaw = 'powerlaw_'+naming+'_'+reg+'_fit.png'
outpath = '/scratch/pit000/results/sid/plots200km/'


#Plotting
fig1 = plt.figure(figsize=(21,7))
ax = fig1.add_subplot(131)
#name = 'ship_radar+buoys+SAR'
#ax.set_title(name,fontsize=29, loc='left')
ax.set_xlabel(r"Length scale (km)",fontsize=25)
ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlim(.1,10)
ax.set_ylim(1e-9,1e-3)

bx = fig1.add_subplot(132)
bx.set_xlabel(r"Length scale (km)",fontsize=25)
bx.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
bx.set_xscale('log')
bx.set_yscale('log')
bx.set_ylim(1e-9,1e-3)


#scale total divergence to average displacement (in one direction)
#dux = n * sigma_x**2 / (2* A * dt).
cx = fig1.add_subplot(133)
cx.set_xlabel(r"Length scale (m)",fontsize=25)
cx.set_ylabel(r"Total displacement (m)",fontsize=25)
cx.set_xscale('log')
cx.set_yscale('log')
#cx.set_xlim(.001,100)
#cx.set_ylim(1e-9,1e-3)


#threshold
#dummy_ls_all, dummy_td_all, dummy_max_td_all
inpath = '/scratch/pit000/results/sid/deform200km/'
date='2015'
#threshold_file = glob(inpath+'dummy_*'+date+'*c'+file_name_end+'.csv')[0]   #use the values from the center tile
threshold_files = glob(inpath+'dummy_*'+date+'*'+rname+'.csv')   #use the values from the center tile

dummy_td_mean=[]
dummy_max=[]

#LS are the same for all files 
dummy_ls = getColumn(threshold_files[0],0,header=False)
dummy_ls = np.array(dummy_ls,dtype=np.float)[:-3]/1000      #dont take the last 3, no data there anyway
#print(dummy_ls)

#make empty lists to store the values
dls_name = [str(int(i)) for i in dummy_ls]
dtd_lists = {ls:[] for ls in dls_name}
mtd_lists = {ls:[] for ls in dls_name}

#These values depend on the time difference and regions: very different from day to day, region to region...
for tf in threshold_files:
    #print(tf)
    dummy_td = getColumn(tf,1,header=False)
    max_td = getColumn(tf,2,header=False)
    dummy_td = np.array(dummy_td,dtype=np.float)[:-3]
    max_td = np.array(max_td,dtype=np.float)[:-3]/3      #found this factor empirically!!!
    
    #collect all values per LS
    for i in range(0,len(dls_name)):
        #print(dls_name[i])
        dtd_lists[dls_name[i]].append(dummy_td[i])
        mtd_lists[dls_name[i]].append(max_td[i])
    
#get one mean/max for every LS
for ii in dls_name:   
    #print(np.mean(dtd_lists[ii]))
    dummy_td_mean.append(np.min(dtd_lists[ii]))
    dummy_max.append(np.max(mtd_lists[ii]))

##plot
#bx.scatter(dummy_ls,dummy_td_mean,marker='x')
#bx.scatter(dummy_ls,dummy_max,marker='.')


    
#get a power law fit on these dummy values
x_all = np.arange(0.1,np.max(dummy_ls))


a_dum1,k_dum1,cix,ciy_upp,ciy_low = logfit(dummy_ls,dummy_td_mean)
ax.loglog(x_all,a_dum1*x_all**k_dum1,linewidth=1,c='k',ls=':')
bx.loglog(x_all,a_dum1*x_all**k_dum1,linewidth=1,c='k',ls=':')

a_dum2,k_dum2,cix,ciy_upp,ciy_low = logfit(dummy_ls,dummy_max)
ax.loglog(x_all,a_dum2*x_all**k_dum2,linewidth=1,c='k',ls=':')
bx.loglog(x_all,a_dum2*x_all**k_dum2,linewidth=1,c='k',ls=':')

#convert to displacements
#get time difference
#/scratch/pit000/results/sid/deform200km/dummy_lance_leg1_20150206T062144_20150207T133434_se_200km.csv
date1=tf.split('_')[3]
date2=tf.split('_')[4]
dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
tls = (dt2-dt1).seconds + (dt2-dt1).days*24*60*60

#convert detection limit to displacements
dummy_ls = dummy_ls*1000    #convert to m
sigma_x = np.sqrt(dummy_td_mean * (2* (dummy_ls)**2 * tls) /3 )
#cx.scatter(dummy_ls,sigma_x,marker='x')

#convert dummy max to displacements
sigma_x_max = np.sqrt(dummy_max * (2* (dummy_ls)**2 * tls) /3 )
#cx.scatter(dummy_ls,sigma_x_max,marker='.')

#power law fits
#make a dummy that starts at 1m
x_all1m = np.arange(1,np.max(dummy_ls))

#lower detection limit
a_dum,k_dum,cix,ciy_upp,ciy_low = logfit(dummy_ls,sigma_x)
cx.loglog(x_all1m,a_dum*x_all1m**k_dum,linewidth=1,c='k',ls=':')

#max detection limit
a_dum,k_dum,cix,ciy_upp,ciy_low = logfit(dummy_ls,sigma_x_max)
cx.loglog(x_all1m,a_dum*x_all1m**k_dum,linewidth=1,c='k',ls=':')
 
#max deformation rates also follows a powerlaw - with apperantly same slope as the threshold!!! with a constant offset between 1 and 2-orders of magnitude!!!
#is this resolution and product independent? Is there any limit in feature tracking?
#exag_fac_max = factor*100                   #80*100=8km total displacement!
#do the SAR, ship radar and buoy max values coincide? ship radar and buoys do, also SAR with LS_eff does...The offest here is about 1 order of magnitude and slope exactly as threshold
#this indicates that max displacement is limited by this method
#check if the max displacements for ship radar on 1-hr and 24-hr are the same - Oikkonen et al, indicated that they are not - so, this could be a real signal

#NOTE: I m not sure if our scaling is any good. In the buoy paper we got differnet spatial scales by calculating deformation from drift at at every scale and not to averaging it like here... Here it is same as in Sylvain's paper.

#step 1 is not resampled to a fresh mesh and has different point cloud - very narrow LS and very high values. Try re-meshing those results too...


print('***********************SAR***************************')

#empty lists to collect the data
meanls_list_sar=[]
meantd_list_sar=[]
stdtd_list_sar=[]

ls_list_sar=[]
td_list_sar=[]

mean_dist_sar=[]
std_dist_sar=[]

inpath = '/scratch/pit000/results/sid/deform200km/'
#infile = inpath+'td_'+reg+file_name_end+'_1000.npz'     #allowing all time differences up to 33 hours
infile = inpath+'td_'+reg+file_name_end+'_strict.npz'   #only allowing 22 to 26h time difference - not working as this is c-tile only...
infile = inpath+'td_'+reg+file_name_end+'_strict30.npz' #now also with minang=30 limitation to see with max_def limit changes...
print(infile)
container = np.load(infile, allow_pickle=True)
ls_lists = container['ls_lists'].tolist()
ls_lists_eff = container['ls_lists_eff'].tolist()
td_lists = container['td_lists'].tolist()
#tls_list = container['tls_list'].tolist()
#print(ls_lists)

n=9
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)
stp_name = [str(i) for i in stp]

#colors
color=iter(plt.cm.Blues_r(np.linspace(0,1,len(stp)+1)))

for jj in stp_name:
    
    #TLS is different for every day and tile
    
    ##flatten the list of arrays
    print(len(td_lists[jj]))
    if len(td_lists[jj])>0:
        td = np.concatenate(td_lists[jj]).ravel()
        ls = np.concatenate(ls_lists[jj]).ravel()/1000 #convert to km
        ls_eff = np.concatenate(ls_lists_eff[jj]).ravel()/1000 #convert to km
        
        #calculate and store averages
        meanls=np.mean(ls); print('mean ls :',meanls)
        meanls_eff=np.mean(ls_eff); print('mean ls EFF :',meanls_eff)
        meantd=np.mean(td); print('mean td :',meantd, 'from', td.shape[0], 'values')
        #meanls_list_sar.append(meanls)
        meanls_list_sar.append(meanls_eff)
        meantd_list_sar.append(meantd)
        stdtd_list_sar.append(np.std(td))
        
        ls_list_sar.extend(ls)
        td_list_sar.extend(td)
            
        #plot all the data
        cl = next(color)
        #ax.scatter(ls, td, lw=0, alpha=.1, color=cl)  # Data
        bx.scatter(ls_eff, td, lw=0, alpha=.1, color=cl)  # Data
        #bx.plot(meanls,meantd,'*',markersize=10,markeredgecolor='w',color=cl)
        bx.plot(meanls_eff,meantd,'*',markersize=10,markeredgecolor='w',color=cl)
        
        
        #convert td to displacement
        #dux = n * sigma_x**2 / (2* A * dt)
        #sigma_x = np.sqrt(dux * (2* A * dt) /n )
        tls = 24 *3600  #approx 1 day in seconds (this needs to be stored too!)
        ls_eff = ls_eff *1000 #convert back to meters
        sigma_x = np.sqrt(td * (2* ls_eff**2 * tls) /3 )
        sigma_x_mean = np.mean(sigma_x)
        mean_dist_sar.append(sigma_x_mean)
        #std_dist_sar.append(np.std(sigma_x))

        cx.scatter(ls_eff, sigma_x, lw=0, alpha=.1, color=cl)
        cx.plot(meanls_eff*1000,sigma_x_mean,'*',markersize=10,markeredgecolor='w',color=cl)

        
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar[:5],meantd_list_sar[:5])   #skip the last LS

#dummy x data for plotting
#x = np.arange(min(meanls_list_sar), max(meanls_list_sar), 1)
x = np.arange(meanls_list_sar[0], meanls_list_sar[5], 1)

bx.loglog(x,a*x**k,linewidth=3,c='b')
bx.loglog(x_all,a*x_all**k,linewidth=3,c='b',ls=':',label=r'$D=%.2f*10^{-6} L^{%.2f}$ (DL & LKFF)' %(a*10e6,k))
bx.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
bx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)        

#introduce other SAR powerlaws in other colors of blues - produced in a separate run of this script!
#color = plt.cm.Blues_r(np.linspace(0,1,len(stp)+1))

a = 1.4631134586019913e-06
k = -0.5949872436341559
bx.loglog(x,a*x**k,linewidth=3,c='hotpink')
bx.loglog(x_all,a*x_all**k,linewidth=3,c='hotpink',ls=':',label=r'$D=%.2f*10^{-6} L^{%.2f}$ (DL)' %(a*10e6,k))

a =4.794292292427179e-07
k =-0.12618239257819433
bx.loglog(x,a*x**k,linewidth=3,c='purple')
bx.loglog(x_all,a*x_all**k,linewidth=3,c='purple',ls=':',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k))

#displacements
meanls_list_sar = np.array(meanls_list_sar)*1000#convert to m
x = x*1000
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sar[:5],mean_dist_sar[:5])
cx.loglog(x,a*x**k,linewidth=3,c='b')
x1 = np.arange(1,30000)
cx.loglog(x1,a*x1**k,linewidth=3,c='b',ls=':',label=r'$D=%.1f L^{%.2f}$ (DL & LKFF)' %(a,k))
cx.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
cx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)
        
#===why are there low values in thresholded data? 
#The threshold values are tile per tile, while dummy is just for the central tile
#also time difference (diff or tls) can be very different inside the same tiled information. This is not implemented yet.
#for now TLS is approximated to exactly 24h (the large gaps over 32h are taken out). This makes the total displacement calculations a bit wrong and causes the scatter bellow/over the dummy lines

#TODO: make experiement, where data is limited to 24h (with only 2 hours deviation) = this will be nice, but SAR-SH-Buouy data wont be over same time = magnitude difference???
#also experiment with the 15 degree angle rule to see if any changes are visible at the max deformation rates...


#===LS differences: 9 classes, why do we have different LS for filtered data??? in non-filtered they match the dummy LS, in filtered they dont!
#problem was that the effective area of seeded triangles was stored and not the nominal (only the area covered by the values)
#effective area for LS gives a cloud of values and nominal gives a line at exact value.
#consequnetly the powerlaw now breaks at 10km -this demonstrates clearly how this is the measurement problem.
#Weiss shows that the power law in temporal scaling does not break. This means we make good measurements of intermittency, but lack goos spatial measurements...



#===max deformation rates are liklely connected to how far we can track the features (e.g. template size at SAR drift and ship radar radius at ship radar)
#feature trackingn is done on 800m regular grid
#img_size=45 (template) in pattern matching: 45x80=3.6km
#detection limit-200m, max deformation about 2km
#ship radar useful radius at N-ICE: ??? likley up to 10km, so allowing large max deformation
#in deformation only input size of triangles is limited, not output.
#This could be max deformation inside winter ice pack (obviously at coasts e.g. North of Greenland in 20xx, summer ice pack, polynyas this can be larger!)






#COPY TO THE PAPER
#results and discussion
#no thersholding and no filtering gives nice power law slope, but much too gentle and wrong values compared to ship radar and buoys
#no thersholding with LKF filtering gives higher values, but still very wrong slope
#thersholding improves the slope
#LKF filtering makes the values higher
#buoy data gets a better magnitude by thresholding, but the slope decreases. However, the confidence interval is wide as there are only 3 classes...
#this kind of SAR thresholding and filtering means lots of missing data at large spatial scales. This means we can coarsen our data beyond 10km. But if we dont do that we get very wrong data anyway - spatial scales break.
#Hypothesis: This happens also because ~10km is typical LKF distance (LKFs are resolved by ship radar, SAR and buoys, but not smaller features). After 10km we start to have 2 LKFs in triangles and this makes our values artifically high. This statement needs further checking...
#About displacements: ship radar all data arrived to 2m at scale 1m. This is extremly close to the mean sea ice thickness (1.5m) measured at N-ICE and to the theoretical limits described by Weiss et al, 2017. SAR data aligns totaly to the 'high masked' ship radar data. Bouy data has some unexpected behavious, but has also wide confidence intervals. Max displacements at 1m are about 50m. This is based on all 3 datasets. This is a max displacement in winter pack ice away from coast (and polynyas) expected within one day at a scale of 1m. BUT what meaning does this have if we know at at MOSAiC met tower moved by 600m in a few hours??? Same for the eastern lead on N-ICE that opened on 7 February. Those large motions were not localized within one feature then? Right, all those triangles were above the spatial scales and filtered out??? NO, it is only the input triangle size that matter and here the smallest were 200m. Does it matter that this 50m is same as one pixes of S-1? Does this mean we could detect small cracks/ridges if they are of such size? No, this could be the limit of the 15 degree angles... And the limit of this method... This is why is scales so nicely. If we want to measure anything else, we need different shapes...
#Buoy data: inner ring is taken as this is the only complete ring - larger rings are asymmetric towards the MIZ and magnitudes are higher there.

#bottomline: 
#we have two detection levels: lower determied by the spatial resolution and upper determined by the method (minimal angle)
#filtering works and values inside the detection limits look good. Different input data give same/similar results.
#SR and B data is not LKF filtered... LKF filtering makes the SAR values more simlar to buoy values: LKF filterign has double effect.
#Besides its original purpose (smoothing and noise removal along LKF) it also removes remaining artifacts (rhomboid corners) and engances the lower detection limit filtering. So: it removes, the artifacts, but the mean values inside LKFs do not change.
#B data filtering is not possible - not enough data
#SH LKF filtering is recommended if any mapps are to be used, but for statistical comparison like we do here, this should not be essential. (Sylvain paper is not agreeing with on this!!!)
#limiting SAR to 24-h is difficult (has to be done per tile in sid_deform_tile.py). If only c-tile is filtered, we get even steeper PL.
#limiting to 30 degrees: no visible effect, no steeping of the powerlaw. Likley because the triangles are already very regular, typicall right angle triangles. #This shows that 15 degrees is a good measure (buoys and ship radar have it) and we arrive to same results as with SAR and 45 degrees??? 


print('***********************ship_radar***************************')
meanls_list_sr=[]
meantd_list_sr=[]
meanls_list_sr_all=[]
meantd_list_sr_all=[]
meanls_list_sr_low=[]
meantd_list_sr_low=[]
stdtd_list_sr_low=[]
stdtd_list_sr_all=[]

ls_list_sr=[]
td_list_sr=[]

r_mean_list=[]
r_n_list=[]

mean_dist_sr=[]
mean_dist_sr_high=[]

#ship radar data
inpath = '/home/pit000/data/ship_radar/'
name1 = 'Annu Oikkonen - Period'
name2 = '_Deformation_L'

lsc_list = range(1,7)

#colors
color=iter(plt.cm.Purples_r(np.linspace(0,1,len(lsc_list)+1)))

#if threshold:
    ##dummy values for ship radar
    #dum2r = [dum2[0],dum2[1],dum2[1],dum2[2],dum2[3],dum2[3]]

for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print(scale)
    
    cl = next(color)
    
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
    stdtd_list_sr_all.append(np.std(td))
    #ls_list_sr.extend(ls)
    #td_list_sr.extend(td)
    
    #plot averages for all the data
    ax.plot(meanls,meantd,'s',markersize=7,markeredgecolor='orange', color=cl)
    #all data
    ax.scatter(ls, td, marker='s', lw=0, alpha=.2, color=cl)
    
    #mask with the dummy values from SAR
    no_all = len(td)
    mask = td < a_dum1*ls**k_dum1
    
    #store low values mean
    ls_low = np.ma.array(ls,mask=~mask)
    td_low = np.ma.array(td,mask=~mask)
    ls_low = np.ma.compressed(ls_low)
    td_low = np.ma.compressed(td_low)
    meanls_list_sr_low.append(np.mean(ls_low))
    meantd_list_sr_low.append(np.mean(td_low))
    stdtd_list_sr_low.append(np.std(td_low))
    #ax.plot(np.mean(ls_low),np.mean(td_low),'s',markersize=7,markeredgecolor='orange', color=cl)
    n1=td_low.size
    
    #high values
    ls_high = np.ma.array(ls,mask=mask)
    td_high = np.ma.array(td,mask=mask)
    ls_high = np.ma.compressed(ls_high)
    td_high = np.ma.compressed(td_high)
    no_high = len(td_high)
    fraction = no_high/no_all
    print('fraction: '+str(fraction)+'from: '+str(no_all))
    
    #calculate and store averages
    meanls=np.mean(ls_high)
    meantd=np.mean(td_high)
    meanls_list_sr.append(meanls)
    meantd_list_sr.append(meantd)
    #ls_list_sr.extend(ls)
    #td_list_sr.extend(td)
    
    #plot the data
    bx.plot(meanls,meantd,'s',markersize=7,markeredgecolor='orange', color=cl)
    #ax.scatter(ls, td, marker='s', lw=0, alpha=.2, color=cl)  # Data

    #convert td to displacement
    tls = 24 *3600 #convert back to seconds
    lsc = ls *1000 #convert back to meters
    sigma_x = np.sqrt(td * (2* lsc**2 * tls) /3 )
    sigma_x_mean = np.mean(sigma_x)
    #cx.scatter(ls, sigma_x, lw=0, alpha=.1, color=cl)
    #cx.plot(meanls,sigma_x_mean,'*',markersize=10,markeredgecolor='w',color=cl)
    mean_dist_sr.append(sigma_x_mean)
    
    lsc = ls_high *1000 #convert back to meters
    sigma_x = np.sqrt(td_high * (2* lsc**2 * tls) /3 )
    sigma_x_mean = np.mean(sigma_x)
    mean_dist_sr_high.append(sigma_x_mean)

#x data for plotting
x = np.arange(min(meanls_list_sr_all), max(meanls_list_sr_all), 1)

#all the ship radar data
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sr_all[:],meantd_list_sr_all[:])
ax.loglog(x,a*x**k,linewidth=2,c='orange')
ax.loglog(x_all,a*x_all**k,linewidth=2,c='orange',ls='--',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k))
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#high ship radar data
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sr[:],meantd_list_sr[:])
bx.loglog(x,a*x**k,linewidth=2,c='orange')
bx.loglog(x_all,a*x_all**k,linewidth=2,c='orange',ls=':',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k))
bx.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
bx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#displacements
x=x*1000
meanls_list_sr = np.array(meanls_list_sr)*1000
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sr,mean_dist_sr)
cx.loglog(x,a*x**k,linewidth=3,c='orange')
cx.loglog(x1,a*x1**k,linewidth=3,c='orange',ls='--',label=r'$D=%.1f L^{%.2f}$' %(a,k))
cx.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
cx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_sr,mean_dist_sr_high)
cx.loglog(x,a*x**k,linewidth=3,c='orange')
cx.loglog(x1,a*x1**k,linewidth=3,c='orange',ls=':',label=r'$D=%.1f L^{%.2f}$' %(a,k))
cx.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
cx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

print('***********************buoys***************************')

ls_list_b=[]
td_list_b=[]

meanls_list_b_all=[]
meantd_list_b_all=[]

meanls_list_b=[]
meantd_list_b=[]

mean_dist_b=[]
mean_dist_b_high=[]

#buoy data - do we have enough of 1 day data? (should be enough for the entire leg 1)
#scales 2-100km
inpath = '/home/pit000/data/buoys/'
#fname_start = 'nice1_comb_SAR'
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
#mask = td>10e-6
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

#if threshold:
    ##dummy values for buoys
    #dum2b = [dum2[3],dum2[4],dum2[4]]

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
    meanls_list_b_all.append(meanls)
    meantd_list_b_all.append(meantd)
    #ls_list_b.extend(ls)
    #td_list_b.extend(td)

    #plot all the data
    cl = next(color)
    ax.scatter(ls_class, td_class, marker='o', lw=0, alpha=.2, color=cl)  # Data
    ax.plot(meanls,meantd,'o',markersize=7,markeredgecolor='yellow', color=cl)

    #mask with the dummy values from SAR
    no_all = len(td_class)
    mask = td_class < a_dum1*ls_class**k_dum1
    
    #high values
    ls_high = np.ma.array(ls_class,mask=mask)
    td_high = np.ma.array(td_class,mask=mask)
    ls_high = np.ma.compressed(ls_high)
    td_high = np.ma.compressed(td_high)
    no_high = len(td_high)
    fraction = no_high/no_all
    print('fraction: '+str(fraction)+'from: '+str(no_all))
    
    #calculate and store averages
    meanls=np.mean(ls_high)
    meantd=np.mean(td_high)
    meanls_list_b.append(meanls)
    meantd_list_b.append(meantd)
    
    #plot the data
    bx.plot(meanls,meantd,'o',markersize=7,markeredgecolor='g', color=cl)
    
    #convert td to displacement
    tls = 24 *3600 #convert back to seconds
    lsc = ls_class *1000 #convert back to meters
    sigma_x = np.sqrt(td_class * (2* lsc**2 * tls) /3 )
    sigma_x_mean = np.mean(sigma_x)
    #cx.scatter(ls, sigma_x, lw=0, alpha=.1, color=cl)
    #cx.plot(meanls,sigma_x_mean,'*',markersize=10,markeredgecolor='w',color=cl)
    mean_dist_b.append(sigma_x_mean)
    
    lsc = ls_high *1000 #convert back to meters
    sigma_x = np.sqrt(td_high * (2* lsc**2 * tls) /3 )
    sigma_x_mean = np.mean(sigma_x)
    mean_dist_b_high.append(sigma_x_mean)

#x data for plotting
x = np.arange(min(meanls_list_b_all), max(meanls_list_b_all), 1)

#all the buoy data
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_b_all[:],meantd_list_b_all[:])
ax.loglog(x,a*x**k,linewidth=2,c='g')
ax.loglog(x_all,a*x_all**k,linewidth=2,c='g',ls='--',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k))
ax.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
ax.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#high buoy data
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_b[:],meantd_list_b[:])
bx.loglog(x,a*x**k,linewidth=2,c='g')
bx.loglog(x_all,a*x_all**k,linewidth=2,c='g',ls=':',label=r'$D=%.2f*10^{-6} L^{%.2f}$' %(a*10e6,k))
bx.plot(cix,ciy_low,'--', c= 'r',linewidth=1,label=r'$99\%\,confidence\,band$')
bx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

#displacements
x=x*1000
meanls_list_b = np.array(meanls_list_b)*1000
a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_b,mean_dist_b)
cx.loglog(x,a*x**k,linewidth=3,c='g')
cx.loglog(x1,a*x1**k,linewidth=3,c='g',ls='--',label=r'$D=%.1f L^{%.2f}$' %(a,k))
cx.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
cx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)

a,k,cix,ciy_upp,ciy_low = logfit(meanls_list_b,mean_dist_b_high)
cx.loglog(x,a*x**k,linewidth=3,c='g')
cx.loglog(x1,a*x1**k,linewidth=3,c='g',ls=':',label=r'$D=%.1f L^{%.2f}$' %(a,k))
cx.plot(cix,ciy_low,'--', c= 'r',linewidth=1)
cx.plot(cix,ciy_upp,'--', c= 'r',linewidth=1)


###################################################Finish the plots
ax.grid('on')
ax.xaxis.set_major_formatter(ScalarFormatter()) 
bx.grid('on')
bx.xaxis.set_major_formatter(ScalarFormatter())
cx.grid('on')
cx.xaxis.set_major_formatter(ScalarFormatter()) 
cx.yaxis.set_major_formatter(ScalarFormatter()) 

ax.legend(loc='upper right',prop={'size':13}, fancybox=True, framealpha=0.5,numpoints=1)
bx.legend(loc='lower left',prop={'size':13}, fancybox=True, framealpha=0.5,numpoints=1)
cx.legend(loc='lower right',prop={'size':13}, fancybox=True, framealpha=0.5,numpoints=1)

fig1.tight_layout()

plt.show()
fig1.savefig(outpath+'power_law_tile_'+naming)
print(outpath+'power_law_tile_'+naming)
