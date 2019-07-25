from sid_func import *
from datetime import datetime, timedelta
import pandas as pd
import matplotlib.pyplot as plt

#wind direction data: (derivative) -detect the direction change
metfile = '../data/10minute_nounits.csv'

mettime = getColumn(metfile,0)[::6]
dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
dtb = [ dtb[i]-timedelta(seconds=5*60) for i in range(len(mettime)) ]  
metspeed = np.asarray(getColumn(metfile,10),dtype=float)
metspeed = np.ma.masked_invalid(metspeed)
metdir = np.asarray(getColumn(metfile,11),dtype=float)
metdir = np.ma.masked_invalid(metdir)
mettemp = np.asarray(getColumn(metfile,6),dtype=float)
mettemp = np.ma.masked_invalid(mettemp)

#average the wind speed for every hour
col = int(len(metspeed)/6)
print(col)
tmp = metspeed.reshape(col,6)
ws = np.mean(tmp,axis=1)

tmp = mettemp.reshape(col,6)
at = np.mean(tmp,axis=1)

tmp = metdir.reshape(col,6)

def circmean(alpha,axis=None):
  #To convert from radians to degrees, multiply by (180o/(PI))
  tod = 180/np.pi
  tor = np.pi/180

  sa = np.mean(np.sin(alpha*tor),axis)
  ca  = np.mean(np.cos(alpha*tor),axis)
  mean_angle = np.arctan2(sa,ca)*tod
  mean_angle = np.where(mean_angle<0,mean_angle+360,mean_angle)
  return mean_angle

wd = circmean(tmp,axis=1)

dir_dt = (wd[:-1]-wd[1:])/60/60 #rate of wind direction change per second (60 min steps)
speed_dt = (ws[:-1]-ws[1:])/60/60

#make pandas time series
dfw = pd.DataFrame({ 'Wind speed' : ws[1:],
                    'Wind direction' : wd[1:],
                    'Air temperature' : at[1:],
                    'Wind speed change' : speed_dt,
                    }, index=dtb[1:])

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

    time = [ datetime.strptime(date[x], '%Y-%m-%d %H:%M:%S') for x in range(0,len(date)) ]
    print(time)
    #time = [ time[x].date()+timedelta(seconds=7*60*60) for x in range(0,len(date)) ]

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
        
    #smallest scale in enough
    break

#big = big.join(dfw, how='outer')
#print(big['2015-01-22'])

#make daily time series of deformation
big = big.resample('D').mean()
print(big)

start = datetime(2015,1,21)
end = datetime(2015,2,16)

#plot
fig, axes = plt.subplots(nrows=3, ncols=1,figsize=(10,8))

#colors
#color=plt.cm.rainbow_r(np.linspace(0,1,len(stp)))

#divergence
ax = big.iloc[:,0].plot(ax=axes[0],title='Sea ice deformation',xlim=(start,end),label=True)
ax.set_ylabel(r'Divergence (s$^-1$)')

ax = big.iloc[:,1].plot(ax=axes[0], linestyle='--',xlim=(start,end))
ax = big.iloc[:,2].plot(ax=axes[0], linestyle='--',xlim=(start,end))

#shear
bx = big.iloc[:,3].plot(ax=axes[1],xlim=(start,end))
bx.set_ylabel(r'Shear (s$^-1$)')

#wind direction
cx = dfw.plot(ax=axes[2],xlim=(start,end))
cx.set_ylabel(r'Wind direction change')

fig.savefig(outpath+'ts.png',bbox_inches='tight')


