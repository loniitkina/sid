from sid_func import *
from scipy import stats
import matplotlib.pyplot as plt

title = 'density_lkf'
radius = '_20km.csv'
radius = '_20kmFW.csv'
#radius = '_7kmFW.csv'

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


#Plotting
fig1 = plt.figure(figsize=(20,6.5))
ax = fig1.add_subplot(131)
name = 'Ship radar'
ax.set_title(name,fontsize=29, loc='left')
ax.set_xlabel(r"Length scale (km)",fontsize=25)
ax.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
ax.set_xscale('log')
ax.set_yscale('log')

#PDF plot
fig2 = plt.figure(figsize=(20,6.5))
zx = fig2.add_subplot(131)
zx.set_xscale('log')
zx.set_title('Ship radar')

xx = fig2.add_subplot(132)
xx.set_xscale('log')
xx.set_title('Buoys')

yx = fig2.add_subplot(133)
yx.set_xscale('log')
yx.set_title('SAR')

logbins = np.logspace(np.log10(1e-10),np.log10(1e-3), 100)

#ymax=.1
#ax.set_ylim(0,ymax)
#ax.set_xlabel('Snow depth (m)')
#srbins = np.arange(0,.8,.01)



#ship radar data
#inpath = '../data/ship_radar/24h/'
inpath = '../sidrift/data/ship_radar/time_sliced_data/'
name1 = 'Annu Oikkonen - Period'
name2 = '_Deformation_L'

lsc_list = range(1,7)

#colors
color=iter(plt.cm.Purples_r(np.linspace(0,1,len(lsc_list)+1)))



for i in range(0,len(lsc_list)):
    scale = lsc_list[i]
    print(scale)
    
    #period 1
    fname = inpath+name1+'1'+name2+str(scale)+'_24h.txt'
    ls1 = getColumn(fname,0, delimiter=' ')
    td1 = getColumn(fname,1, delimiter=' ')
    
    ##period 2
    #fname = inpath+name1+'2'+name2+str(scale)+'_24h.txt'
    #ls2 = getColumn(fname,0, delimiter=' ')
    #td2 = getColumn(fname,1, delimiter=' ')
    
    ##combine
    #ls = np.array(np.append(ls1,ls2),dtype=np.float)/1000  #convert from m to km
    #td = np.array(np.append(td1,td2),dtype=np.float)/3600  #convert from h to s   
    
    ls = np.array(ls1,dtype=np.float)/1000  #convert from m to km
    td = np.array(td1,dtype=np.float)/3600  #convert from h to s
    
    
    #calculate and store averages
    meanls=np.mean(ls)
    meantd=np.mean(td)
    meanls_list_sr.append(meanls)
    meantd_list_sr.append(meantd)
    ls_list_sr.extend(ls)
    td_list_sr.extend(td)

    #density plots
    #x, y, z = density(ls,td)           #density classes will not be log-spaced (low td values will be all in a single wide color class)
    x, y, z = density_lsb(ls,td,n=100)  #log-spaced ls and td bins
    ax.scatter(x, y, c=z, s=50, edgecolor='',cmap=plt.cm.jet, alpha=.2)
    #plot the means
    cl = next(color)
    ax.plot(meanls,meantd,'s',markersize=7,markeredgecolor='orange', color=cl)
    
    #plot PDFs for each ls
    weights = np.ones_like(td) / (len(td))
    n, bins, patches = zx.hist(td, logbins, lw=3, alpha=0.5, weights=weights, label=str(scale), histtype='step')





bx = fig1.add_subplot(132)
name = 'Buoys'
bx.set_title(name,fontsize=29, loc='left')
bx.set_xlabel(r"Length scale (km)",fontsize=25)
#bx.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
bx.set_xscale('log')
bx.set_yscale('log')


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
 
#td = np.append(np.load(fname_td1,encoding='latin1',allow_pickle=True),np.load(fname_td2,encoding='latin1',allow_pickle=True))/24/60/60      #convert from days to s
#ls = np.append(np.load(fname_ls1,encoding='latin1',allow_pickle=True),np.load(fname_ls2,encoding='latin1',allow_pickle=True))

#use only first part
td = np.load(fname_td1,encoding='latin1',allow_pickle=True)/24/60/60      #convert from days to s
ls = np.load(fname_ls1,encoding='latin1',allow_pickle=True)


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
color=iter(plt.cm.Blues_r(np.linspace(0,1,len(lsc_list)+1)))


for i in range(0,len(lsc_list)):
    print(i)
    mask = (ls<minlen[i]) | (ls>maxlen[i])
    ls_class = np.ma.array(ls,mask=mask)
    td_class = np.ma.array(td,mask=mask)
    ls_class = np.ma.compressed(ls_class)
    td_class = np.ma.compressed(td_class)   
      
    #calculate and store averages
    meanls=np.mean(ls_class); print(meanls)
    meantd=np.mean(td_class); print(meantd)
    meanls_list_b.append(meanls)
    meantd_list_b.append(meantd)
    ls_list_b.extend(ls)
    td_list_b.extend(td)

    #density plots
    x, y, z = density_lsb(ls_class,td_class,n=100)  #log-spaced ls and td bins
    bx.scatter(x, y, c=z, s=50, edgecolor='',cmap=plt.cm.jet, alpha=.2)
    cl = next(color)
    bx.plot(meanls,meantd,'o',markersize=7,markeredgecolor='yellow', color=cl)
    
    #plot PDFs for each ls
    weights = np.ones_like(td_class) / (len(td_class))
    n, bins, patches = xx.hist(td_class, logbins, lw=3, alpha=0.5, weights=weights, label=str(i), histtype='step')


cx = fig1.add_subplot(133)
name = 'SAR'
cx.set_title(name,fontsize=29, loc='left')
cx.set_xlabel(r"Length scale (km)",fontsize=25)
#cx.set_ylabel(r"Total deformation (s$^{-1}$)",fontsize=25)
cx.set_xscale('log')
cx.set_yscale('log')


#SAR data
inpath = '../sidrift/data/40m_combo/'
outpath = inpath
fname_start = 'td_leg1_L'
fname_start = 'td_seed_f_Lance_L'

n=8 # number of samples
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)

##for 50km radius
#n=9
##size envelope also needs to increase (from 10m to 3km)
#stp=np.exp(np.linspace(np.log(1),np.log(800),n))
#stp = stp.astype(int)

margin = np.exp(np.linspace(np.log(.1),np.log(3),n))

for i in range(0,len(stp)-5):                           #the last two steps are off the curve, try removing them
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
        
    ##mask all very small or big triangles
    ##if not masked the range of the ls is big and has several clouds (expected ls, twice the ls and all kinds of smaller ls)
    #center = np.mean(ls)
    ##center = stats.mode(ls)[0][0]                      #this takes too much time
    #print(center)
    #minlen = center-margin[i]; maxlen = center+margin[i]
    ##minlen = center-margin; maxlen = center+margin
    #mask = (ls<minlen) | (ls>maxlen)
    #ls = np.ma.array(ls,mask=mask)
    #td = np.ma.array(td,mask=mask)
    #ls = np.ma.compressed(ls)
    #td = np.ma.compressed(td) #*1e6   
    
    ##throw away very low deformation rates (template noise)
    #mask = td<.6e-6
    #ls = np.ma.array(ls,mask=mask)
    #td = np.ma.array(td,mask=mask)
    #ls = np.ma.compressed(ls)
    #td = np.ma.compressed(td)   

    
    #calculate and store averages
    meanls=np.mean(ls)
    meantd=np.mean(td); print(meantd)
    meanls_list_sar.append(meanls)
    meantd_list_sar.append(meantd)
    
    ls_list_sar.extend(ls)
    td_list_sar.extend(td)
        
    #density plots
    x, y, z = density_lsb(ls,td,n=100)  #log-spaced ls and td bins
    #x = np.ma.masked_where(~(np.isfinite(z)),x)
    #x = np.ma.compressed(x)
    #y = np.ma.masked_where(~(np.isfinite(z)),y)
    #y = np.ma.compressed(y)
    #z = np.ma.masked_where(~(np.isfinite(z)),z)
    #z = np.ma.compressed(z)
    #print(z.shape)
    #print(np.max(z))
    #print(np.min(z))    #how can z have negative value??? - this is just a random interval of number in which histogram is scaled? not number of hits?
    #exit()
    cx.scatter(x, y, c=z, s=50, edgecolor='',cmap=plt.cm.jet, alpha=.2)
    #means
    cx.plot(meanls,meantd,'*',markersize=10,markeredgecolor='w',color='k')

    #plot PDFs for each ls
    weights = np.ones_like(td) / (len(td))
    n, bins, patches = yx.hist(td, logbins, lw=3, alpha=0.5, weights=weights, label=str(scale), histtype='step')



ax.grid(True)
bx.grid(True)
cx.grid(True)

ax.set_xlim(.01,25.)
ax.set_ylim(1e-10,1e-3)
bx.set_xlim(.01,25.)
bx.set_ylim(1e-10,1e-3)
cx.set_xlim(.01,25.)
cx.set_ylim(1e-10,1e-3)
#cx.set_ylim(1e-10*1e6,1e-3*1e6)


from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())  
bx.xaxis.set_major_formatter(ScalarFormatter())
cx.xaxis.set_major_formatter(ScalarFormatter())
#ax.legend(loc='lower left',prop={'size':16}, fancybox=True, framealpha=0.5,numpoints=1)
fig1.tight_layout()

fig1.savefig(outpath+'power_law_24h_'+title)

zx.legend()
xx.legend()
yx.legend()
fig2.savefig(outpath+'power_law_pdf_lkf')



##density plots
#from scipy.stats import gaussian_kde

## Generate fake data
#x = ls
#y = td

## Calculate the point density
#xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)

## Sort the points by density, so that the densest points are plotted last
#idx = z.argsort()
#x, y, z = x[idx], y[idx], z[idx]

##fig, ax = plt.subplots()
#ax.scatter(x, y, c=z, s=50, edgecolor='')
#plt.show()
