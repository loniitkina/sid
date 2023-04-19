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

filtering=False

reg = 'lance_leg1'

inpath = '/scratch/pit000/results/sid/deform200km/'
file_name_end = '_200km'
outname_td = inpath+'td_nofilter_'+reg+file_name_end+'_1000.csv'
powerlaw = 'powerlaw_test_nofilter_'+reg+'_1000.png'


outpath = '/scratch/pit000/results/sid/plots200km/'
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

fnames = sorted(glob(inpath+'Scaling_*.npz'))
fnames = sorted(glob(inpath+'Scaling_800_2015*.npz'))

td_list=[]
ls_list=[]

for fname in fnames:
    print(fname)

    #get date
    date = fname.split('_')[2]
    #print(date)
    dt = datetime.strptime(date, "%Y%m%dT%H%M%S")
    
    #get nominal resolution (distance between nods)
    dst = int(fname.split('_')[1])
    #dst=800
    
    #load data
    container = np.load(fname, allow_pickle=True)
    tripts = container['tripts']
    div = container['div']
    shr = container['shr']
    area = container['area']
    minang = container['minang']
    max_area = container['max_area']
    corner_nods = container['corners']
    threshold = container['threshold']
    
    if filtering:
        tripts = tripts[~threshold]
        div = div[~threshold]
        shr = shr[~threshold]
        area = area[~threshold]
        minang = minang[~threshold]
        
    #threshold
    #dummy_ls_all, dummy_td_all, dummy_max_td_all
    threshold_file = glob(inpath+'dummy_*'+date+'*c'+file_name_end+'.csv')[0]   #use the values from the center tile
    print(threshold_file)
    dummy_ls = getColumn(threshold_file,0,header=False)
    dummy_td = getColumn(threshold_file,1,header=False)
    max_td = getColumn(threshold_file,2,header=False)
    dummy_ls = np.array(dummy_ls,dtype=np.float)
    dummy_td = np.array(dummy_td,dtype=np.float)  
    max_td = np.array(max_td,dtype=np.float)
    
    #Plot the maps - we need also area_def
    #map area definition
    area_def_file = glob(inpath+'area_def_'+date+'*'+file_name_end+'*.pkl')[0]
    with open(area_def_file, "rb") as pkl:
        area_def = pickle.load(pkl)
    
    #fig2, bx = plt.subplots(3, 3,figsize=(20,20))
    fig2    = plt.figure(figsize=(20,10))
    bx = fig2.add_subplot(3,3,1)
    
    m = pr.plot.area_def2basemap(area_def)
    m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
    m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
    
    #bx[0].plot(corner_nods[0][0],corner_nods[0][1],'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
    
    patches_all = []
    for k in range(div.shape[0]):
        patch = Polygon(tripts[k,:,:])
        patches_all.append(patch)

    #plot filled triangles
    p = PatchCollection(patches_all, cmap=plt.cm.bwr, alpha=1)
    p.set_array(div*1e6)
    p.set_clim(interval)
    bx.add_collection(p)
    
    td = np.sqrt(div**2+shr**2)
    ls = np.sqrt(area)
    
    ax.scatter(ls,td,alpha=.5)
    ax.scatter(dummy_ls,dummy_td,marker='x')
    ax.scatter(dummy_ls,max_td,marker='.')
    
    td_list.append(td)
    ls_list.append(ls)

    idx=2
    for j in stp[1:]:
        print(j)
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
        div_coarse,shr_coarse,area_coarse,minang_coarse,tripts_coarse = coarse_grain(tripts,tripts_seed,div,shr,max_area,minang)
        
        #get stats
        print(len(div_coarse))
        
        if len(div_coarse)>0:
        
            td = np.sqrt(div_coarse**2+shr_coarse**2)
            ls = np.sqrt(area_coarse)
            
            td_list.append(td)
            ls_list.append(ls)
            
            
            bx = fig2.add_subplot(3,3,idx)
    
            m = pr.plot.area_def2basemap(area_def)
            m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
            m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])
            
            patches_all = []
            for k in range(div_coarse.shape[0]):
                patch = Polygon(tripts_coarse[k,:,:])
                patches_all.append(patch)

            #plot filled triangles
            p = PatchCollection(patches_all, cmap=plt.cm.bwr, alpha=1)
            p.set_array(div_coarse*1e6)
            p.set_clim(interval)
            bx.add_collection(p)
            
        idx=idx+1
    
    multiplot = 'multiplot_nofilter'+reg+'_'+date+'_'+file_name_end+'_1000.png'
    fig2.savefig(outpath+multiplot)

fig1.savefig(outpath+powerlaw)

#flatten this list of lists
td_list = np.concatenate(td_list).ravel()
ls_list = np.concatenate(ls_list).ravel()

#write data for scatter plots
print('Storing data for the scatter plots')
tt = [ls_list,td_list]
table = list(zip(*tt))

print(outname_td)
with open(outname_td, 'wb') as f:
    #header
    f.write(b'lenght scale, total deformation\n')
    np.savetxt(f, table, fmt="%s", delimiter=",")

