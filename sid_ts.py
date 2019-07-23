from sid_func import *
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt

#wind direction data: (derivative) -d does the directin change
outname2 = 'nice1_comb_ts'

wspeed = np.load('wspeed'+outname2)
wdir = np.load('wdir'+outname2)
atemp = np.load('wtemp'+outname2)
dates = np.load('dates'+outname2)


#SAR data
inpath = '/Data/sim/polona/sid/deform/'
outpath = '../plots/'
fname_start = 'ts_leg1_L'

#create log-spaced vector and convert it to integers
n=8 # number of samples
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)
print(stp)
#get a grip on the lengh scales
print(stp*.04)

#for 50km radius
n=9
stp=np.exp(np.linspace(np.log(1),np.log(800),n))
stp = stp.astype(int)
print(stp)

#distance in m
lenght = stp*40
print(lenght)
#exit()

radius = 50000
file_name_end = '_50km_more.csv'

for i in range(0,len(stp)-2): #skip the last two scales (no data anyway)  
    scale = stp[i]
    print(scale)
    fname = inpath+fname_start+str(scale)+file_name_end
    print(fname)
    
    date = getColumn(fname,0, delimiter=',')
    pdiv = getColumn(fname,1, delimiter=',')
    ndiv = getColumn(fname,2, delimiter=',')
    div = getColumn(fname,3, delimiter=',')
    shr = getColumn(fname,4, delimiter=',')

    pdiv = np.array(pdiv,dtype=np.float)
    ndiv = np.array(ndiv,dtype=np.float)
    div = np.array(div,dtype=np.float)
    shr = np.array(shr,dtype=np.float)

    time = [datetime.strptime(date[x], '%Y-%m-%d %H:%M:%S') for x in range(0,len(date))]

    #make pandas time series
    df = pd.DataFrame({ 'Divergence scale:%i m' % (lenght[i]) : div,
                        'Pos. div. scale:%i m' % (lenght[i]) : pdiv,
                        'Neg.div. scale:%i m' % (lenght[i]) : ndiv,
                        'Shear scale:%i m' % (lenght[i]): shr,
                        }, index=time)
    
    #merge with shorter lenghts
    if i == 0:
        big = df
    else:
        big = big.join(df, how='outer')


print(big)
#plot
fig, axes = plt.subplots(nrows=2, ncols=1,figsize=(10,8))

#colors
color=plt.cm.rainbow_r(np.linspace(0,1,len(stp)))

#divergence
ax = big.iloc[:,0::4].plot(ax=axes[0],title='Sea ice deformation',color=color[:-2])
ax.set_ylabel(r'Divergence (s$^-1$)')

ax = big.iloc[:,1::4].plot(ax=axes[0],color=color[:-2], linestyle='--', label=False)
ax = big.iloc[:,2::4].plot(ax=axes[0],color=color[:-2], linestyle='--', label=False)

#shear
bx = big.iloc[:,3::4].plot(ax=axes[1],color=color[:-2])
bx.set_ylabel(r'Shear (s$^-1$)')

fig.savefig(outpath+'ts.png',bbox_inches='tight')


