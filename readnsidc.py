import numpy as np # for array handling
import scipy as sp # scientific python
import matplotlib.pyplot as plt # for basic plotting
from subprocess import call # for doing system calls - not really needed
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap # for maps
import datetime as datetime

plt.close("all")

basepath = '/home/rkm/work/DATA/NSIDC/passivemicro/north/timeseries/'
fbase = 'seaice_conc_monthly_nh_'
timp = '197811-201212'

timeidx=346 # index 346 s/b 9/2007

#cmapc = plt.cm.gist_yarg
cmapc = plt.cm.YlGnBu

baseparams = dict(projection='npstere',boundinglat=45,lon_0=
                  310,resolution='c')
climparams = dict(shading='flat',latlon=True,cmap=cmapc)
diffparams = dict(shading='flat',latlon=True,cmap=plt.cm.RdBu,vmin=-.05,vmax=.05)


# ======= first plot GODDARD data

dtyp = 'godd'  # or 'nsidc'
fname = basepath + fbase + timp + '_' + dtyp +  '_timeseries.nc'

ncfile = Dataset(fname,'r')

lat = ncfile.variables['latitude'][:]
lon = ncfile.variables['longitude'][:]
gsi = ncfile.variables['goddard_merged_seaice_conc_monthly'][timeidx,:,:]
dattime = ncfile.variables['time'][timeidx]
  #double time(time) ;
  #  time:standard_name = "time" ;
  #  time:long_name = "ANSI date" ;
  #  time:units = "days since 1601-01-01 00:00:00" ;
  #  time:calendar = "standard" ;

datdate = datetime.date(1601, 1, 1) + datetime.timedelta(int(dattime))

# guessing at lon_0 for NSIDC projection (=50W?)
plt.figure()
bm = Basemap(**baseparams)
pc = bm.pcolormesh(lon,lat,gsi,**climparams)
bm.drawcoastlines()
#bm.drawmapboundary(fill_color='#99ffff')
bm.drawlsmask()
bm.drawparallels(np.arange(-90.,120.,30.))
bm.drawmeridians(np.arange(0.,360.,30.))
plt.title(datdate)

# add colorbar.
cbar = bm.colorbar(pc,location='bottom',pad="5%")
cbar.set_label('fraction')

# ===== now plot NSIDC data
dtyp = 'nsidc'
fname = basepath + fbase + timp + '_' + dtyp +  '_timeseries.nc'

ncfile = Dataset(fname,'r')
nsi = ncfile.variables['seaice_conc_monthly_cdr'][timeidx,:,:]

# guessing at lon_0 for NSIDC projection (=50W?)
plt.figure()
bm = Basemap(projection='npstere', boundinglat=45,lon_0=310,resolution='c')
pc = bm.pcolormesh(lon,lat,nsi,**climparams)
bm.drawcoastlines()
#bm.drawmapboundary(fill_color='#99ffff')
bm.drawlsmask()
bm.drawparallels(np.arange(-90.,120.,30.))
bm.drawmeridians(np.arange(0.,360.,30.))
plt.title(datdate)

# add colorbar.
cbar = bm.colorbar(pc,location='bottom',pad="5%")
cbar.set_label('fraction')


# ==== plot NSIDC - GODD

sidiff = nsi-gsi
# Note that plt.cm.RdBu_r is the reverse of RdBu colormap

plt.figure()
bm = Basemap(projection='npstere', boundinglat=45,lon_0=310,resolution='c')
pc = bm.pcolormesh(lon,lat,sidiff,**diffparams)
bm.drawcoastlines()
#bm.drawmapboundary(fill_color='#99ffff')
bm.drawlsmask()
bm.drawparallels(np.arange(-90.,120.,30.))
bm.drawmeridians(np.arange(0.,360.,30.))
plt.title(datdate)

# add colorbar.
cbar = bm.colorbar(pc,location='bottom',pad="5%")
cbar.set_label('fraction')
#cbar.set_clim(vmin=-.1,vmax=.1) # This leaves the colorbar lims as the original and looks weird

# TO SEE ALL PYPLOT COLORMAPS
#execfile('pycmaps.py')



