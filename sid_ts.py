from sid_func import *
from datetime import datetime, timedelta
import pandas as pd
import matplotlib.pyplot as plt

#wind direction data: (derivative) -detect the direction change
metfile = '../sidrift/data/10minute_nounits.csv'

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
#inpath = '../sidrift/data/80m_stp10_single_filter/'
inpath = '../sidrift/data/80m_stp10_nofilter/'
outpath = inpath
fname_start = 'ts_seed_f_Lance_L1'

radius = 7000
file_name_end = '_7km.csv'

fname = inpath+fname_start+file_name_end


date = getColumn(fname,0, delimiter=',')
pdiv = getColumn(fname,1, delimiter=',')
ndiv = getColumn(fname,2, delimiter=',')
div = getColumn(fname,3, delimiter=',')
shr = getColumn(fname,4, delimiter=',')

#print(date)

#convert from s-1 to hours-1
pdiv = np.array(pdiv,dtype=np.float);pdiv = np.ma.array(pdiv,mask=pdiv==-999) 
ndiv = np.array(ndiv,dtype=np.float);ndiv = np.ma.array(ndiv,mask=ndiv==-999) 
div = np.array(div,dtype=np.float);div = np.ma.array(div,mask=div==-999)
shr = np.array(shr,dtype=np.float);shr = np.ma.array(shr,mask=shr==-999)

time = [ datetime.strptime(date[x], '%Y-%m-%d %H:%M:%S') for x in range(0,len(date)) ]
#print(time)
#time = [ time[x].date()+timedelta(seconds=7*60*60) for x in range(0,len(date)) ]

#make pandas time series
df = pd.DataFrame({ 'Divergence' : div,
                    'Pos. div.'  : pdiv,
                    'Neg.div.'   : ndiv,
                    'Shear'      : shr,
                    }, index=time)

#make daily time series of deformation
df_d = df.resample('D').mean()


#get data from Oikkonen et al, 2017
inpath_radar = '../sidrift/plots/'
fname = inpath_radar+'Oikkonen_ts_shear.csv'
tmp = getColumn(fname,0, delimiter=',')
time_radar = [ datetime.strptime(tmp[x], '%Y/%m/%d/%H') for x in range(0,len(tmp)) ]
shear_radar = getColumn(fname,1, delimiter=',')
shear_radar = np.array(shear_radar,dtype=np.float)/3600 #chnage units to s-1

#some dates are not sorted
time_radar = [x for x,_ in sorted(zip(time_radar,shear_radar))]

shear_radar = [y for _,y in sorted(zip(time_radar,shear_radar))]

#use temporal scaling law exponent to scale from hourly to daily data
#there will be just some 5 points to compare!!! >> is this worth doing???


#print(time_radar)
#exit()

dfr = pd.DataFrame({ 'Shear_radar' : shear_radar
                    }, index=time_radar)


print(dfr)

#make daily time series of deformation
dfr_d = dfr.resample('D').mean()

print(dfr_d)
#exit()

big = df_d.join(dfr_d, how='outer')
print(big)



start = datetime(2015,1,21)
end = datetime(2015,2,9)

#plot
fig, axes = plt.subplots(nrows=3, ncols=1,figsize=(10,8))

#colors
#color=plt.cm.rainbow_r(np.linspace(0,1,len(stp)))

#divergence
ax = big.iloc[:,0].plot(ax=axes[0],title='Sea ice deformation',xlim=(start,end),label=True)
ax.set_ylabel(r'Divergence (s$^-1$)')

ax = big.iloc[:,1].plot(ax=axes[0], linestyle='--',xlim=(start,end),ylim=(-1e-4,1e-4))
ax = big.iloc[:,2].plot(ax=axes[0], linestyle='--',xlim=(start,end))

#shear
#bx = big.iloc[:,3].plot(ax=axes[1])#,xlim=(start,end),ylim=(0,1e-4))
#bx = big.iloc[:,4].plot(ax=axes[1])
bx = df.iloc[:,3].plot(ax=axes[1])
bx = dfr.plot(ax=axes[1])

bx.set_ylabel(r'Shear (s$^-1$)')

#wind direction
cx = dfw.plot(ax=axes[2],xlim=(start,end))
cx.set_ylabel(r'Wind direction change')

fig.savefig(outpath+'ts.png',bbox_inches='tight')


