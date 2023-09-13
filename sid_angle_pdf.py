import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from sid_func import getColumn

#PDFs for winter heat budget
inpath = '/scratch/pit000/results/sid/deform200km/'
outpath = '/scratch/pit000/results/sid/plots200km/'

outname = 'pdf_angles.png'

fnames = glob(inpath+'angle_mosaic2019*_200km.csv')
#fnames = glob(inpath+'angle_*2015*_200km.csv')

#cols = plt.cm.rainbow(np.linspace(0, 1, len(dates)))


srbins = np.arange(0,90,1)

#PDFs
fig1 = plt.figure(figsize=(10,10))
#fig1.suptitle(title, fontsize=30)
ax = fig1.add_subplot(111)
#ax.text(-.13, .355, "a", ha="center", va="center", size=45)  #make simple figure annotation
ax.set_xlabel('Intersection angle', fontsize=25)
ax.set_ylabel('Probability', fontsize=25)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
#ax.set_xlim(0,1)
#ax.set_ylim(0,0.35)

angles=[]
ss_mean=[]
ii_mode=[]
for fn in fnames:
        
    aa = getColumn(fn,0)
    aa = np.array(aa,dtype=np.float)
    
    angles.extend(aa)

#print(angles)
print(len(angles))    

angles = np.array(angles)
angles = np.where(angles>90,angles-90,angles)
        
weights = np.ones_like(angles) / (len(angles))
n, bins, patches = ax.hist(angles, srbins, histtype='step', color='royalblue', linewidth=4, alpha=.5, weights=weights)
#ax.plot([mn,mn],[0,ymax],c='b',ls='--', label='mean = '+str(round(mn,2))+' m')
    
#for i in range(0,len(locs)):
    #print(locs[i])
    #mn=ss_mean[i]
    #ax.plot([mn,mn],[0,ymax],c='r',ls=':')
    
    #mo=ii_mode[i]
    #bx.plot([mo,mo],[0,ymax],c='r',ls=':')


#ax.legend(fontsize=20)

ax.text(60, .04, '$N$= '+str(np.round(len(angles),2)), ha="left", va="center", size=20)

print(outpath+outname)
fig1.savefig(outpath+outname,bbox_inches='tight')
