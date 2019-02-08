from glob import glob
from datetime import datetime
import numpy as np
import pyresample as pr
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from sid_func import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

#select lenght scale
lscale = 'full'
print(lscale)

#-------------------------------------------------------------------
inpath = '../output/drift_'+str(lscale)+'/'
inpath = '/Data/sim/polona/sid/test/'
outpath_def = '../output/def_'+str(lscale)+'/'
outpath = outpath_def+'plots/'
metfile = '../data/10minute_nounits.csv'
reg = 'leg1'
proj = reg

fl = sorted(glob(inpath+'*size35.npz'))
i = fl[0]


#read in all the data
print(i)
container = np.load(i)
print(container.files)
u = container['upm'] 
v = container['vpm'] 
rpm = container['rpm'] #cross-correlation matrix
hpm = container['hpm'] #hessian

plt.imshow(rpm)
plt.savefig(inpath+'rpm_35')

plt.imshow(hpm)
plt.savefig(inpath+'hpm_35')

plt.imshow(hpm*rpm)
plt.savefig(inpath+'rpm_hpm_35')

speed = np.sqrt(u**2+v**2)
plt.imshow(speed)
plt.savefig(inpath+'speed_35')
    
