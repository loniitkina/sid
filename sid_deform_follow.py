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

#How long do you want it to run?
first_week=True
#first_week=False    #This will make it run for all the data

#Do you want all output in figures?
image=True

#select lenght scale
lscale = 'full'
print(lscale)
##steps between the image pixels (1 pixel = 40m)
##for 25km radius
#radius = 25000
#file_name_end = '_25km.csv'
#for 7.5km radius
radius = 7500
file_name_end = '_7km.csv'

#create log-spaced vector and convert it to integers
#start with ~500m step
n=8 # number of samples
stp=np.exp(np.linspace(np.log(1),np.log(300),n))
stp = stp.astype(int)

##for 50km radius
#radius = 50000
#file_name_end = '_50km.csv'
#n=9
#stp=np.exp(np.linspace(np.log(1),np.log(800),n))
#stp = stp.astype(int)

#distance in m
print(stp*40)
print(stp*40*10)

##check just the largest triangles
#stp = stp[-2:]


#-------------------------------------------------------------------
#inpath = '../sidrift/data/'
#inpath = '../sidrift/data/test1/'
inpath = '../sidrift/data/whole_series_10stp_factor_def/data/'
outpath_def = inpath
#outpath = '../sidrift/plots/'
outpath = '../sidrift/data/whole_series_10stp_factor_def/plots/'
metfile = '../sidrift/data/10minute_nounits.csv'
reg = 'leg1'
proj = reg
#reg = 'leg1_FYI'
#reg = 'leg1_SYI'

##virtual buoys
#out_file = outpath_def+'VB.npz'
#container = np.load(out_file)
#print(container.files)
#lon_path = container['lon_path'][:6,:,:]   #first step is just vb grid initialization
#lat_path = container['lat_path'][:6,:,:]


for j in stp:
    print('Step: '+str(j))
    outname_td = 'td_'+reg+'_L'+str(j)+file_name_end
    td_list=[]
    ls_list=[]
    ang_list=[]
    time_list=[]
    date_list=[]
    date_ts=[]
    mdiv=[]
    mpdiv=[]
    mndiv=[]
    mshr=[]

    fl = sorted(glob(inpath+'*.npz'))
    fl = fl[1:]
    for i in range(0,len(fl)):
        #read in all the data
        print(fl[i])
        container = np.load(fl[i])
        print(container.files)
        u = container['upm'][::j, ::j]  #this is a regular grid/matrix and such indexing produces a resonably well distributed mesh to construct non-acute triangles
        print(u.shape)
        #continue
        v = container['vpm'][::j, ::j] 
        rpm = container['rpm'][::j, ::j] #cross-correlation matrix
        hpm = container['hpm'][::j, ::j] #hessian
        lat = container['lat1'][::j, ::j] 
        lon = container['lon1'][::j, ::j] 
        
        lat2 = container['lat2'][::j, ::j] 
        lon2 = container['lon2'][::j, ::j] 
        
        #get time difference
        date1 = fl[i].split('_')[-2]
        date2 = fl[i].split('_')[-1].split('.')[0]
        dt1 = datetime.strptime(date1, "%Y%m%dT%H%M%S")
        dt2 = datetime.strptime(date2, "%Y%m%dT%H%M%S")
        
        #if we want just the data until the week-long gap
        if (dt1 > datetime(2015,1,28)) & (first_week==True): continue

        #if date1 != '20150126T070245': continue        #hack to give me just this one day (for plotting of Figure 1)

        diff = (dt2-dt1).seconds + (dt2-dt1).days*24*60*60
        
        #Lance postion (from Lance's met system)
        mettime = getColumn(metfile,0)
        dtb = [ datetime.strptime(mettime[i], "%Y-%m-%d %H:%M:%S") for i in range(len(mettime)) ]
        if dtb[0]>dt1: continue
        if dtb[-1]<dt1: continue
        mi = np.argmin(abs(np.asarray(dtb)-dt1))
        Lance_lon = np.asarray(getColumn(metfile,2),dtype=float)[mi]
        Lance_lat = np.asarray(getColumn(metfile,1),dtype=float)[mi]
        if np.isnan(Lance_lon): continue
        
        #select a region around Lance that corresponds to the buoy array area
        #radius of 100 km
        #otherwise to much sea ice deep in the Arctic is included (reduced scaling slope)
        
        #convert Lance position and sea ice drift array to projected image coordinates
        #project the coordinates (units of distance have to be meters)
        #use QGIS for the first guess about the coordinates
        #area_def = pr.utils.load_area('area.cfg', proj)
        

        
        #Lance centered projection
        from pyproj import Proj, transform
        inProj = Proj(init='epsg:4326')
        outProj = Proj('+proj=laea +lat_0=%f +lon_0=%f +a=6378137 +b=6356752.3142 +units=m' % (90, 10))
        xlp,ylp = transform(inProj,outProj,Lance_lon, Lance_lat)
        
        from pyresample.geometry import AreaDefinition
        #Using a projection dictionary
        area_id = 'around Lance'
        description = 'North Pole LAEA Europe'
        proj_id = 'lance'
        proj_dict = {'proj':'laea', 'lat_0':90, 'lon_0':10, 'a':6378137, 'b':6356752.3142, 'units':'m'}
        width = radius*2/40 #x m spacing
        height = radius*2/40 #1x m spacing
        area_extent = (xlp-radius,ylp-radius,xlp+radius,ylp+radius)
        area_def = AreaDefinition(area_id, description, proj_id, proj_dict, width, height, area_extent)

        print(area_def)
        #exit()

        
        m = pr.plot.area_def2basemap(area_def)
            
        #reproject vertices
        #lon = lon.reshape
        x, y = m(lon, lat)
        
        ##recalculate u,v
        #x2, y2 = m(lon2, lat2)
        #u = (x2-x)/diff
        #v = (y2-y)/diff
        
        u = u/diff
        v = v/diff        
        
        print(u)
        #exit()
               
        #reproject Lance position
        xl, yl = m(Lance_lon, Lance_lat)
        
        if reg == 'leg1_FYI':   #shift region 25km southwards into the pure FYI zone
            yl = yl-25000

        if reg == 'leg1_SYI':   #shift region 25km northwards into the pure SYI zone
            yl = yl+25000
            
        #cut out region
        mask = (x<xl-radius) | (x>xl+radius) | (y<yl-radius) | (y>yl+radius)
        
        x = x[~mask]
        y = y[~mask]
        u = u[~mask]
        v = v[~mask]
        rpm = rpm[~mask]
        hpm = hpm[~mask]
       
        #get rid of all the nans, the resulting arrays are flattened
        #x = np.ma.masked_where(~(np.isfinite(u)),x)
        x = np.ma.compressed(x)
        #y = np.ma.masked_where(~(np.isfinite(u)),y)
        y = np.ma.compressed(y)

        
        #u = np.ma.masked_invalid(u)
        u = np.ma.compressed(u)
        #v = np.ma.masked_invalid(v)
        v = np.ma.compressed(v)
        #rpm = np.ma.masked_invalid(rpm)
        rpm = np.ma.compressed(rpm)
        #hpm = np.ma.masked_invalid(hpm)
        hpm = np.ma.compressed(hpm)
        #print lon.shape
        
        #check how many values are left in the region
        print(u)
        print(u.shape)
        #continue
        
        ##mask out all poor quality data: rpm < 0.4
        #gpi = hpm > 5
        #u = u[gpi]
        #v = v[gpi]
        #x = x[gpi]
        #y = y[gpi]
        
        print(u.shape)
        
        if i == 0:
            print('First time step')
            #triangulate betwen the points
            pts = np.zeros((len(x),2))
            pts[:,0]=x; pts[:,1]=y
            try:
                tri = Delaunay(pts)
            except:
                continue
            #simple plot check
            fig1    = plt.figure(figsize=(10,10))
            ax      = fig1.add_subplot(111)
            ax.triplot(pts[:,0], pts[:,1], tri.simplices.copy())
            grid_filename = outpath+'grid'+reg+'_L'+str(j)+'_'+date1
            fig1.savefig(grid_filename)
            ##exit()
            
            tripts = pts[tri.simplices]
            
            #setup displacement that acummulates along the time
            pts2 = np.zeros((len(x),2))
            pts2[:,0]=x; pts2[:,1]=y
            
        else:
            #coordinates of the nods need to be updated from the last step
            #add dispacement from previous step
            pts2[:,0]=pts2[:,0]+dx; pts2[:,1]=pts2[:,1]+dy
            
            print(pts2)
            print(pts2.shape)
            #print(pts)
            #print(pts.shape)
            #print(u.shape)

            #check plot
            fig1    = plt.figure(figsize=(10,10))
            ax      = fig1.add_subplot(111)
            ax.triplot(pts2[:,0], pts2[:,1], tri.simplices.copy())
            grid_filename = outpath+'grid'+reg+'_L'+str(j)+'_'+date1
            fig1.savefig(grid_filename)
            #exit()

            tripts = pts2[tri.simplices]
        
        #calculate displacement for the next step
        dx = u*diff; dy = v*diff
        
        #but some u,v data can be really bad
        #if dx or dy is greater than the step - then what???
        #some interpolation?
        #new triangulation for those spots???
        
        
        
        #print(dx)
        #print(u)
        #exit()
        
        #velocities
        upts = u[tri.simplices]
        vpts = v[tri.simplices]
        
        #calculate deformation - calculating for triangles one by one?
        dux=[];duy=[];dvx=[];dvy=[];minang=[];area=[]
        for t in range(0,len(tripts)):
            vert = np.asarray(tripts[t])
            uvert = upts[t]
            vvert = vpts[t]
        
            try:
                a,b,c,d,e,f=deformation(vert,uvert,vvert)
            except:
                continue
            dux.append(a);duy.append(b);dvx.append(c);dvy.append(d);minang.append(e);area.append(f)
            
        #print area
        dux = np.array(dux)
        duy = np.array(duy)
        dvx = np.array(dvx)
        dvy = np.array(dvy)
        minang = np.array(minang)
        area = np.array(area)
        
        print(minang.shape)
        print(minang)
        #exit()
        
        
        div = dux + dvy
        shr = .5*np.sqrt((dux-dvy)**2+(duy+dvx)**2)
        td = np.sqrt(div**2 + shr**2)
        ls = np.sqrt(area)
            
        ##can we mask out somehow the values cause by the speed stepfunction?
        ##how do we get a relative displacement?
        ##if known we could filter out all that are exactly the smallest step...
        ##lets start by plotting a pdf of divergence...
        #ax = plt.subplot(131)
        #ddd = div*1e6
        #weights = np.ones_like(ddd) / (len(ddd))
        #dbins = np.arange(-3,3,.1)
        #n, bins, patches = ax.hist(ddd, dbins, weights=weights)
        
        #bx = plt.subplot(132)
        #sss = shr*1e6
        #n, bins, patches = bx.hist(sss, dbins, weights=weights)     
        
        #cx = plt.subplot(133)
        #ttt = td*1e6
        #n, bins, patches = cx.hist(ttt, dbins, weights=weights)
        
        #plt.show()
        ##exit()

            
            
        #storing data for the scatter plots
        td_list.extend(td.tolist())
        ls_list.extend(ls.tolist())
        ang_list.extend(minang.tolist())
        
        #print(td)
        #print(diff_tri)
        #print(dt_tri.astype(datetime))
        #exit()
        dt_tri = np.full_like(dux,np.datetime64(dt1))
        diff_tri = np.ones_like(dux)*diff
        
        date_list.extend(dt_tri.tolist())
        time_list.extend(diff_tri.tolist())
        
        #print(len(diff_tri.tolist()))
        
        #save the data for time series
        date_ts.append(dt1)
        posd = np.ma.array(div,mask=div<0)
        negd = np.ma.array(div,mask=div>0)
        mpdiv.append(np.mean(posd))
        mndiv.append(np.mean(negd))
        mdiv.append(np.mean(div))
        mshr.append(np.mean(shr))
        
        ##save for plotting of Figure 1
        #if date1 == '20150126T070245':        #hack to give me just this one day (for plotting of Figure 1)
            #print(date1)
            #out_file = outpath_def+'Fig1_data_'+date1+'.npz'
            #np.savez(out_file,tripts = tripts,shr = shr, div = div)
            #exit()
        
        #if (i < 6) & (j==1):  #only for the first 6 image pairs and highest resolution
        if image == True:
            #continue
            ##Plotting 
            #deform = div*1e6
            ##print(deform)
            #outname = 'map_div'+reg+'_L'+str(j)+'_'+date1
            #label = r'Divergence (10$^6$s$^{-1}$)'
            #interval = [-5, 5]
            #cmap=plt.cm.bwr
            ##title = 'c'
            
            deform = shr*1e6
            outname = 'map_shr_'+reg+'_L'+str(j)+'_'+date1
            label = r'Shear (10$^6$s$^{-1}$)'
            interval = [0, 10]
            cmap=plt.cm.Reds

            #speed = np.sqrt(u**2+v**2)
            #outname = 'map_speed_'+reg+'_L'+str(j)+'_'+date1
            #label = r'Speed (m/s)'
            #interval = [0, 10000]
            #cmap=plt.cm.nipy_spectral

                        
            print(outname)
            #deformation plots
            fig3    = plt.figure(figsize=(20,20))
            cx      = fig3.add_subplot(111)
            #cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
            #cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

            #area_def = pr.utils.load_area('area.cfg', proj)  
            m = pr.plot.area_def2basemap(area_def)
            
            #scale
            m.drawmapscale(Lance_lon, Lance_lat-.3, Lance_lon+8, Lance_lat-.2, 50, units='km', barstyle='fancy',fontsize=14)
            m.drawmeridians(np.arange(0.,360.,5.),latmax=90.,labels=[0,0,0,1,])
            m.drawparallels(np.arange(79.,90.,1),labels=[1,0,0,0])


            #nods
            #cx.plot(x,y,'o',linewidth=2,color='purple')
            
            ##speed
            #sc = cx.scatter(x,y,s=50,c=speed,cmap=cmap, vmin=.06, vmax=.085)         #add colorbar and remove extreme values
            #cbar = plt.colorbar(sc)
            #cbar.ax.set_ylabel('Drift speed (m/s)',size=22)

            
            #triangles
            patches = []
            for k in range(deform.shape[0]):
                patch = Polygon(tripts[k,:,:], edgecolor='orchid', alpha=1, fill=False)
                #plt.gca().add_patch(patch)
                patches.append(patch)
            
            #plot filled triangles
            p = PatchCollection(patches, cmap=cmap, alpha=0.4)
            p.set_array(np.array(deform))
            p.set_clim(interval)
            cx.add_collection(p)
            
            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            divider = make_axes_locatable(cx)
            cax = divider.append_axes("bottom", size="5%", pad=0.1)
            cbar = plt.colorbar(p, cax=cax, orientation='horizontal')
            cbar.set_label(label,size=16)

            #Lance
            #xl, yl = m(Lance_lon, Lance_lat)
            cx.plot(xl,yl,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')
            
            ##virtual buoys
            #xb,yb = m(lon_path[i,:,:],lat_path[i,:,:])
            #cx.plot(xb,yb,'o',linewidth=2,color='purple')

            print(outpath+outname)
            fig3.savefig(outpath+outname,bbox_inches='tight')
            #exit()
        else:
            exit() #just make these plots and exit
                

    #continue
    #write out lists into csv file
    tt = [date_list, ls_list, time_list, td_list, ang_list]
    table = zip(*tt)
    #adjusted to python3:
    table = list(zip(*tt))

    output = outpath_def + outname_td
    with open(output, 'wb') as f:
        #header
        f.write(b'date, length scale, time difference, total deformation, min angle\n')
        np.savetxt(f, table, fmt="%s", delimiter=",")

    #write data for time series 
    tt = [date_ts, mpdiv, mndiv, mdiv, mshr]
    table = zip(*tt)
    table = list(zip(*tt))

    outname_ts = 'ts_'+reg+'_L'+str(j)+file_name_end
    output = outpath_def + outname_ts
    with open(output, 'wb') as f:
        #header
        f.write(b'date, pos. divergence, neg. divergence, mean divergence, mean shear\n')
        np.savetxt(f, table, fmt="%s", delimiter=",")
