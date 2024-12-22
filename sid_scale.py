from glob import glob
from datetime import datetime
import numpy as np
import pyresample as pr
from scipy.spatial import Delaunay
import pickle
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from sid_func import coarse_grain, getColumn

#supress warnings
import warnings
warnings.filterwarnings("ignore")

#WARNING: filtered data has large no-data areas. Adjustment of cov_limit may be required to get some coase scale data
#Potentially implementation of cov_limit=0 should be done here...

filtering=True
#filtering=False
#This is threshold filtering.
#All data here are stored by sid_deform_tile.py and are already lkf-filtered/kernel-filtered

reg = 'lance_leg1'

inpath = '/scratch/pit000/results/sid/deform200km/'
inpath = '/scratch/pit000/results/sid/deform200km_stp5_factor1/'

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
    #lname='_nolkfilter'

file_name_end = rname+tname+lname

print('Doing: ',file_name_end, reg)

#outname_td = inpath+'td_'+reg+file_name_end+'_1000.npz'
#powerlaw = 'powerlaw_test_'+reg+file_name_end+'_1000.png'

minang_limit=30

#strict time difference filtering
outname_td = inpath+'td_'+reg+file_name_end+'_strict30.npz'
powerlaw = 'powerlaw_test_'+reg+file_name_end+'_strict30.png'

outpath = '/scratch/pit000/results/sid/plots200km/'
outpath = '/scratch/pit000/results/sid/plots200km_stp5_factor1/'
interval = [-1, 1]  #expected divergence values (div in s-1 * 10e6) for filter plot

plt.rcParams['agg.path.chunksize'] = 10000
fig1 = plt.figure(figsize=(7.5,7))
ax = fig1.add_subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')

#The coarse-grained data for scaling have a limited lenght scale - defined by these steps
#For filtered data this is shorter
#This is because the coverage needs to be at least 50%
#create log-spaced vector and convert it to integers
n=9 # number of samples
minscale=1  #shortest scale = resolution x minscale - 1x800m=800m
maxscale=1000#longest scale = resolution x maxscale - 300x800m=24km, 1000x80=80km
#1=800m, 300=24km
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)

#make empty lists to store the values for the powerlaw plots
stp_name = [str(i) for i in stp]
td_lists = {sn:[] for sn in stp_name}
ls_lists = {sn:[] for sn in stp_name}
ls_lists_eff = {sn:[] for sn in stp_name}
tls_list = []

#input produced by sid_deform_tile.py
fnames = sorted(glob(inpath+'Scaling_*.npz'))
fnames = sorted(glob(inpath+'Scaling_800_2015*'+file_name_end+'_tiled.npz'))
fnames = sorted(glob(inpath+'Scaling_200_2015*'+file_name_end+'_stp1_tiled.npz'))

for fname in fnames:
    print(fname)
    #get date
    #date = fname.split('_')[2]
    date = fname.split('_')[4]
    dt = datetime.strptime(date, "%Y%m%dT%H%M%S")
    
    #check that we have max 26h time difference (discard the data spanning over missing days) - name come just from c-tile, so this is not ideal...
    #date2 = fname.split('_')[3]
    date2 = fname.split('_')[5]
    dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
    diff = (dt2-dt).seconds + (dt2-dt).days*24*60*60
    #if diff> 33*60*60:  #morning afternoon difference can be 32 hours - very liberal setting to cover as much of time in January/February (similar to TS data)
    if (diff< 22*60*60)|(diff> 26*60*60):  #strict setting to get only daily data
        print('Too short/long time difference: ',date,date2)
        continue
        
    #get nominal resolution (distance between nods)
    #dst = int(fname.split('_')[1])
    dst = int(fname.split('_')[3])
    #dst=800
    
    #load data
    container = np.load(fname, allow_pickle=True)
    tripts = container['tripts']
    div = container['div']
    shr = container['shr']
    area = container['area']
    minang = container['minang']
    max_area = container['max_area']        #this seems too large, and it should increase by steps!
    corner_nods = container['corners']
    threshold = container['threshold']
    #diff = container['dt']
    
    print(np.sqrt(max_area))
    max_area = (dst*1.3)**2/2  #this will also reduce wide leads and shear zones, but otherwise we get various separate clusters
    print(np.sqrt(max_area))
    
    print(div.shape)
    
    if filtering:
        tripts = tripts[~threshold]
        div = div[~threshold]
        shr = shr[~threshold]
        area = area[~threshold]
        minang = minang[~threshold]
        #diff = diff[~threshold]
    
    print(div.shape)
    
    #reduce pointless analysis
    #check size and shape of input triangles - if filtering this was already done as part of thresholding
    #too_big = area>max_area
    too_big = area>max_area*10
    #check their min angles
    too_sharp = minang<minang_limit
    #mask out bad triangle
    too_bad = too_big | too_sharp
    tripts = tripts[~too_bad]
    div = div[~too_bad]
    shr = shr[~too_bad]
    area = area[~too_bad]
    minang = minang[~too_bad]
    #diff = diff[~too_bad]
    
    print(div.shape)    #ok, we get less ad less values, this is OK.
    
    #threshold
    #dummy_ls_all, dummy_td_all, dummy_max_td_all
    #threshold_file = glob(inpath+'dummy_*'+date+'*c'+file_name_end+'.csv')[0]   #use the values from the center tile
    threshold_file = glob(inpath+'dummy_'+reg+'_'+date+'*c'+rname+'.csv')[0]    #there should only be one such file!
    print(threshold_file)
    dummy_ls = getColumn(threshold_file,0,header=False)
    dummy_td = getColumn(threshold_file,1,header=False)
    max_td = getColumn(threshold_file,2,header=False)
    dummy_ls = np.array(dummy_ls,dtype=np.float)
    dummy_td = np.array(dummy_td,dtype=np.float)  
    max_td = np.array(max_td,dtype=np.float)
    
    #Plot the maps - we need also area_def
    #map area definition
    area_def_file = glob(inpath+'area_def_'+date+'*'+rname+'*.pkl')[0]
    with open(area_def_file, "rb") as pkl:
        area_def = pickle.load(pkl)
    
    #fig2, bx = plt.subplots(3, 3,figsize=(20,20))
    fig2    = plt.figure(figsize=(20,20))
    #bx = fig2.add_subplot(3,3,1)
    
    #m = pr.plot.area_def2basemap(area_def)
    #m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    #m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
    
    ##bx[0].plot(corner_nods[0][0],corner_nods[0][1],'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    
    #patches_all = []
    #for k in range(div.shape[0]):
        #patch = Polygon(tripts[k,:,:])
        #patches_all.append(patch)

    ##plot filled triangles
    #p = PatchCollection(patches_all, cmap=plt.cm.bwr, alpha=1)
    #p.set_array(div*1e6)
    #p.set_clim(interval)
    #bx.add_collection(p)
    
    #td = np.sqrt(div**2+shr**2)
    #ls = np.sqrt(area)
    
    #print('means for this L: ', np.mean(ls), np.mean(td))
    
    #ls_m = np.int(np.mean(ls))
    #bx.set_title(('L'+str(ls_m)))#title
    
    #ax.scatter(ls,td,alpha=.5)
    #ax.scatter(np.mean(ls),np.mean(td),marker='*')
    
    #ax.scatter(dummy_ls,dummy_td,marker='x')
    #ax.scatter(dummy_ls,max_td,marker='.')
    
    #tls_list.append(diff)   #only appended for one scale (the rest are the same) - no, they need to be scaled too!
    
    #or comment this out and do remeshing for all scales/steps
    #td_lists[stp_name[0]].append(td)
    #ls_lists[stp_name[0]].append(ls)
    #ls_lists_eff[stp_name[0]].append(ls)

    #idx=2
    #for j in stp[1:]:
    idx=1   #subplot index
    for j in stp:   
        print('step: ',j)
        spacing=dst*j

        #use corner nods to create a regular grid of seed nods
        #print(corner_nods)
        xs, ys = np.mgrid[corner_nods[0][0]:corner_nods[2][0]:spacing, corner_nods[0][1]:corner_nods[2][1]:spacing]
        pts_seed = np.zeros((len(xs.flatten()),2))
        pts_seed[:,0]=xs.flatten(); pts_seed[:,1]=ys.flatten()
        
        #triangulate between these seeding points
        tri_seed = Delaunay(pts_seed)
        tripts_seed = pts_seed[tri_seed.simplices]
        
        #coarsen
        div_coarse,shr_coarse,area_coarse,area_coarse_eff,minang_coarse,tripts_coarse = coarse_grain(tripts,tripts_seed,div,shr,minang,size_limit=max_area,ang_limit=minang_limit, cov_limit=.1)
        
        #get stats
        print('points left: ',len(div_coarse))
        
        if len(div_coarse)>0:
        
            td = np.sqrt(div_coarse**2+shr_coarse**2)
            ls = np.sqrt(area_coarse)
            ls_eff = np.sqrt(area_coarse_eff)
            
            #plot these on the powerlaw plot too
            ax.scatter(ls,td,alpha=.5)
            ax.scatter(np.mean(ls),np.mean(td),marker='*')
            
            print('means for this L: ', np.mean(ls), np.mean(ls_eff), np.mean(td))
            
            td_lists[str(j)].append(td)
            ls_lists[str(j)].append(ls)
            ls_lists_eff[str(j)].append(ls_eff)
            
            ls_m = np.round(np.mean(ls)/1000,1)
            print(ls_m)
            
            bx = fig2.add_subplot(3,3,idx)
            bx.set_title(('$\lambda$'+str(ls_m)+' km'))#title
    
            m = pr.plot.area_def2basemap(area_def)
            m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
            m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
            
            #plot grey background, so that the missing data is better visible
            bx.set_facecolor('0.9')
                
            #draw coarse triangles with colors
            patches_all = []
            for k in range(div_coarse.shape[0]):
                patch = Polygon(tripts_coarse[k,:,:])
                patches_all.append(patch)

            #plot filled triangles
            p = PatchCollection(patches_all, cmap=plt.cm.bwr, alpha=1)
            p.set_array(div_coarse*1e6)
            p.set_clim(interval)
            bx.add_collection(p)
            
            if ls_m > 1:
                #draw original shortest triangles over this as empty polygons with black lines
                patches = []
                for k in range(div.shape[0]):
                    patch = Polygon(tripts[k,:,:])
                    patches.append(patch)
                
                p = PatchCollection(patches, ec= '0.5', fc=None, alpha=0.5, lw=.5)
                bx.add_collection(p)
            
        idx=idx+1
    
    #plt.show()
    multiplot = 'multiplot_'+file_name_end+'_'+reg+'_'+date+'_'+file_name_end+'_1000.png'
    
    fig2.savefig(outpath+multiplot,bbox_inches='tight')
    print('Multiplot figure saved to: '+outpath+multiplot)
    
    #Plot over the dummies
    ax.scatter(dummy_ls,dummy_td,marker='x')    #interesting: there is a huge jump at 10km even in dummies, this must be some artifact...
    ax.scatter(dummy_ls,max_td,marker='.')

fig1.savefig(outpath+powerlaw)

##flatten this list of lists
#td_list = np.concatenate(td_list).ravel()
#ls_list = np.concatenate(ls_list).ravel()

##write data for scatter plots
#print('Storing data for the scatter plots')
#tt = [ls_list,td_list]
#table = list(zip(*tt))

#print(outname_td)
#with open(outname_td, 'wb') as f:
    ##header
    #f.write(b'lenght scale, total deformation\n')
    #np.savetxt(f, table, fmt="%s", delimiter=",")

#write data for scatter plots
np.savez(outname_td, ls_lists=ls_lists, ls_lists_eff=ls_lists_eff, td_lists=td_lists, tls_list=tls_list)
print('Stored data for the scatter plots to: ', outname_td)
