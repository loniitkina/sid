import numpy as np
from glob import glob
from datetime import datetime
from sid_func import getColumn 
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.dates import DateFormatter, date2num
import locale

locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')


radius = 200000 #some hard-coded naming depends on this
rname = '_'+str(int(radius/1000))+'km'
shipfile = '../../downloads/lance_leg1.csv'
inpath = '/scratch/pit000/results/sid/drift/stp10_factor05'+rname+'/'

#for getting timestampts of drift files
main_trackfile_tail = '_c'+rname+'-fnames.csv'
main_trackfile=shipfile.split('.csv')[0]+main_trackfile_tail
print('Your study area is: ',main_trackfile)

noons = getColumn(main_trackfile,0,header=False)
days = [ datetime.strptime(noons[i], "%Y-%m-%d %H:%M:%S") for i in range(len(noons)) ]
print(days)

N = len(days)
df = pd.DataFrame({
           'x': np.random.randn(N),
           'y': np.random.randn(N),
           'z': pd.date_range(start=days[0], periods=N, freq='D')},index=days)

fig, ax = plt.subplots(figsize=(10,4))

smap = ax.scatter(df.x, df.y, s=50,
                   c=[date2num(i.date()) for i in df.z], cmap=plt.cm.jet_r) # <==

fig.gca().set_visible(False)
#cb = fig.colorbar(orientation="horizontal", cax=cax)

cb = fig.colorbar(smap, orientation='horizontal',format=DateFormatter('%d %b')) # <==
cb.ax.tick_params(labelsize=20)
#cb = fig.colorbar(smap, orientation='vertical')
#cb.ax.set_yticklabels(df.index.strftime('%d %b'))

#plt.show()
fig.savefig('tiling_map_colorbar',bbox_inches='tight')
