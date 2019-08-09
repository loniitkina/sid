import json
from glob import glob
from datetime import datetime, timedelta
import numpy as np

#set date/time we want to extract
maptime = datetime(2015,1,27,12,0,0)

#get all files
path = '/Data/sim/polona/sid/buoys/'
fl = glob(path+'*.json')

print(fl)

names = []
dates = []
lons = []
lats = []

#open all files
for i in fl:
    print(i)
    with open(i, "r") as read_file:
        data = json.load(read_file)
        name = data['properties']['buoy']
        print(name)
        dt = data['properties']['measured']
        
        #search for the right date
        date = [ datetime.strptime(dt[k], "%Y-%m-%dT%H:%M:%SZ") for k in range(len(dt)) ]
        mi = np.argmin(abs(np.asarray(date)-maptime))
        lonlat = data['geometry']['coordinates'][mi]
        print(date[mi])
        
        #if date is much later than the maptime (spring deployment), skip this buoy
        #else store the data
        if (date[mi] < maptime-timedelta(hours=5)) | (date[mi] > maptime+timedelta(hours=5)):
            print('Bad buoy or spring deployment.'); continue
        else:
            names.append(name); dates.append(date[mi]); lons.append(lonlat[0]); lats.append(lonlat[1])

print(names)
print(dates)
print(lons,lats)
    
#write csv file with buoy name, date, time, lat, lon
tt = [names, dates, lons, lats]
table = list(zip(*tt))

output = path + 'buoys_fw.csv'
with open(output, 'wb') as f:
    #header
    f.write(b'names, dates, lons, lats\n')
    np.savetxt(f, table, fmt="%s", delimiter=",")
