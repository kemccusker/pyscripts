# plotnsidc_allmonths.py
#    1/31/2014
#    plot the difference in two time periods for each month of the year
#    Start with 2000-2010 minus 1979-1989 to mimic the CanESM2 output
#
#    

import numpy as np # for array handling
import scipy as sp # scientific python
import matplotlib.pyplot as plt # for basic plotting
from subprocess import call # for doing system calls - not really needed
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap # for maps
import datetime as datetime

plt.close("all")

# @@ put in a constants module later?
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
printtofile=0

basepath = '/home/rkm/work/DATA/NSIDC/passivemicro/north/timeseries/'
fbase = 'seaice_conc_monthly_nh_'
timpb = '1979-1989' # base time period
timp  = '2000-2010' # 


# seaice_conc_monthly_nh_2000-2010climo_godd.nc

dtyp = 'godd'  # or 'nsidc'
field = 'goddard_merged_seaice_conc_monthly'
fnameb = basepath + fbase + timpb + 'climo_' + dtyp +  '.nc'
fname = basepath + fbase + timp + 'climo_' + dtyp +  '.nc'



#cmapc = plt.cm.gist_yarg
cmapc = plt.cm.YlGnBu # climo colormap

baseparams = dict(projection='npstere',boundinglat=45,lon_0=
                  310,resolution='c') # for basemap
climparams = dict(shading='flat',latlon=True,cmap=cmapc) # for pcolormesh
diffparams = dict(shading='flat',latlon=True,cmap=plt.cm.RdBu,vmin=-.1,vmax=.1) # for pcolormesh

# example of subplots:
#   # Four axes, returned as a 2-d array
#     f, axarr = plt.subplots(2, 2)
#     axarr[0, 0].plot(x, y)
#     axarr[0, 0].set_title('Axis [0,0]')
#     axarr[0, 1].scatter(x, y)
#     axarr[0, 1].set_title('Axis [0,1]')
#     axarr[1, 0].plot(x, y ** 2)
#     axarr[1, 0].set_title('Axis [1,0]')
#     axarr[1, 1].scatter(x, y ** 2)
#     axarr[1, 1].set_title('Axis [1,1]')
#
#     #Fine-tune figure; hide x ticks for top plots and y ticks for right plots
#     plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
#     plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
#
#   Four polar axes
#   plt.subplots(2, 2, subplot_kw=dict(polar=True))
#   plt.show()


# ======= first plot GODDARD data

ncfileb = Dataset(fnameb,'r') # base time period
ncfile = Dataset(fname,'r')   # later time period

lat = ncfile.variables['latitude'][:]
lon = ncfile.variables['longitude'][:]
fld = ncfile.variables[field][...]
fldb = ncfileb.variables[field][...]

flddiff = fld - fldb # does this work? 

#fig = plt.figure()
##fig.subplots_adjust(wspace=0,hspace=0)

  #fig,axarr = plt.subplots(2,6)
fig, ax = plt.subplots(2,6)
fig.set_size_inches(12,6)
fig.subplots_adjust(hspace=0.1)

for midx in range(len(months)):
#for midx in range(1,3):
    plotfld = flddiff[midx,:,:]

    plt.subplot(2,6,midx)
    
    bm = Basemap(**baseparams)
    pc = bm.pcolormesh(lon,lat,plotfld,**diffparams)
    bm.drawcoastlines()
    #bm.drawmapboundary(fill_color='#99ffff')
    bm.drawlsmask()
    bm.drawparallels(np.arange(-90.,120.,30.))
    bm.drawmeridians(np.arange(0.,360.,30.))

#cbar = fig.colorbar(pc,location='bottom',pad="5%")
#cbar.set_label('fraction')

#Add an axes at position rect [left, bottom, width, height]
#where all quantities are in fractions of figure width and height.
cbar_ax = fig.add_axes([.95,.15,.02,.7])
fig.colorbar(pc,cax=cbar_ax)
fig.suptitle('sea ice concentration')

plt.show()

if printtofile:
    plt.savefig('monthlyseaice.pdf',dpi=100)



# Positioning one colorbar with a set of subplots:
#
# from: http://stackoverflow.com/questions/13784201/matplotlib-2-subplots-1-colorbar
#
##import numpy as np
##import matplotlib.pyplot as plt

##fig, axes = plt.subplots(nrows=2, ncols=2)
##for ax in axes.flat:
##    im = ax.imshow(np.random.random((10,10)), vmin=0, vmax=1)

##fig.subplots_adjust(right=0.8)
##cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
##fig.colorbar(im, cax=cbar_ax)

##plt.show()

### More code example related to subplots:
#     from same webpage
#
##fig = plt.figure()
##fig.subplots_adjust(wspace=0,hspace=0)

### subplot number 1
##ax1 = fig.add_subplot(1,2,1,aspect='equal',xlim=[-2,2],ylim=[-2,2])
##plt.title(r"$\gamma_{1}$",fontsize="18")
##plt.xlabel(r"x ($\theta_{E}$)",fontsize="15")
##plt.ylabel(r"y ($\theta_{E}$)",rotation='horizontal',fontsize="15")
##plt.xticks([-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5])
##plt.xticks([-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5])
##plt.imshow(g1out,extent=(-2,2,-2,2))
##plt.axhline(y=0,linewidth=2,color='k',linestyle="--")
##plt.axvline(x=0,linewidth=2,color='k',linestyle="--")
##e1 = patches.Ellipse((0,0),2,2,color='white')
##ax1.add_patch(e1)

### subplot number 2
##ax2 = fig.add_subplot(1,2,2,sharey=ax1,xlim=[-2,2],ylim=[-2,2])
##plt.title(r"$\gamma_{2}$",fontsize="18")
##plt.xlabel(r"x ($\theta_{E}$)",fontsize="15")
##ax2.yaxis.set_major_formatter( NullFormatter() )
##plt.axhline(y=0,linewidth=2,color='k',linestyle="--")
##plt.axvline(x=0,linewidth=2,color='k',linestyle="--")
##plt.imshow(g2out,extent=(-2,2,-2,2))
##e2 = patches.Ellipse((0,0),2,2,color='white')
##ax2.add_patch(e2)
