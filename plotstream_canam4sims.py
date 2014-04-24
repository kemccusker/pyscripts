# plotstream_canam4sims.py
#    2/24/2014
#    Streamfunction plots for CanAM4 sea ice simulations
#    (taken from plotvert_canam4sims.py)
#    

import numpy as np # for array handling
import numpy.ma as ma # masked array
import scipy as sp # scientific python
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt # for basic plotting
import matplotlib.cm as cm
from subprocess import call # for doing system calls - not really needed
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap # for maps
import datetime as datetime
import matplotlib.colors as col
import platform as platform
import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import matplotlib.font_manager as fm


# while I'm still creating these modules, have to reload to get changes
cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)

plt.close("all")
plt.ion()

printtofile=1

model = 'CanAM4'
ftype = 'ts'    # timeseries
#ftype = 'climo'  # 12-month climatology
climo=0
timeavg = 'SON'

if timeavg in ('DJF','MAM','NDJ','JJA','SON'):
    diffcontsn=[-21,-18,-15,-12,-10,-8,-6,-4,-2,-1] # negative
    diffcontsp=[0,1,2,4,6,8,10,12,15,18,21]      # positive
    diffcontsnz=diffcontsn    # negative zoom
    diffcontspz=diffcontsp    # positive zoom
    contsn=[-8,-6,-4,-2]
    contsp=[0,2,4,6,8,10,12,14,16,18,20]
    contsnz=[-5,-4,-3,-2,-1.5,-1,-.5]
    contspz=[0,.5,1,1.5,2,3]
 
# # # ########### set Simulations #############
# Control run
casename = 'kemctl1'
#timstr = '001-061'
timstr = '001-111'
styr = 2             # skip year 1
enyr = 61 

# Pert run
casenamep1 = 'kem1pert1'  # 2002-2012 sic and sit
casenamep2 = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
casenamep3 = 'kem1pert3'  # 2002-2012 sic, adjusted sst. control sit
#timstrp = '001-061'
timstrp = '001-111'
styrp = 2             # skip year 1
enyrp = 61

cmap = 'RdBu_r'
cmapclimo = 'Spectral_r'
conv = 1   # conversion factor to convert units, etc

# # # ######## set Field info ###############
field = 'zmpsi'
ncfield = 'ZMPSI'
units = '10^10 kg/s' # @@
conv=1/1e10
## cminc = 200
## cmaxc = 310
## cmin = -.5
## cmax = .5
## cminm = -.8 # for monthly
## cmaxm = .8  # for monthly
## cminsc = -2.5 # as Screen et al paper
## cmaxsc = 2.5  # as Screen et al paper



# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'

fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_' + ftype + '.nc'
fnamep1 = basepath + casenamep1 + subdir + casenamep1 + '_' + field + '_' + timstrp + '_' + ftype + '.nc'
fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp + '_' + ftype + '.nc'
fnamep3 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp + '_' + ftype + '.nc'

# Get the data
ncfilec = Dataset(fnamec,'r') # control
ncfilep1 = Dataset(fnamep1,'r') # pert1
ncfilep2 = Dataset(fnamep2,'r') # pert2
ncfilep3 = Dataset(fnamep3,'r') # pert3

lat = ncfilec.variables['lat'][:]
#lon = ncfilec.variables['lon'][:]
lev = ncfilec.variables['plev'][:]

# time, lev, lat

if ftype == 'ts':
    fldc = ncfilec.variables[ncfield][(styr-1)*12:(enyr*12+1),:,:]*conv # time start year to end
    fldp1 = ncfilep1.variables[ncfield][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
    fldp2 = ncfilep2.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
    fldp3 = ncfilep3.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
else:
    fldc = ncfilec.variables[ncfield][...]*conv 
    fldp1 = ncfilep1.variables[ncfield][...]*conv
    fldp2 = ncfilep2.variables[ncfield][...]*conv # time start year to end
    fldp3 = ncfilep3.variables[ncfield][...]*conv # time start year to end

    
nt,nlev,nlat = fldc.shape # if nt == 12 then it's a climo

#  ############# set which perturbation run! ######
casenamep = casenamep2
fldp = fldp2
#  ###############################################

seasfldc = np.squeeze(cutl.seasonalize_monthlyts(fldc,timeavg,climo=climo))
seasfldp = np.squeeze(cutl.seasonalize_monthlyts(fldp,timeavg,climo=climo))

print seasfldc.shape


lats,levs = np.meshgrid(lat,lev)

tstat,pval = sp.stats.ttest_ind(seasfldp,seasfldc,axis=0)
print tstat.shape
print timeavg
pval = ma.masked_invalid(pval) # pcolormesh needs masked_array

#   SEASONAL MEAN CLIMO (CONTROL)
plotfld = np.average(seasfldc,0)


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
CS1 = plt.contour(lats,levs/100,plotfld,contsn,\
            colors='k',linestyles='dashed')
plt.clabel(CS1,fmt = '%2.1f',inline=0,fontsize=10)
CS1 = plt.contour(lats,levs/100,plotfld,contsp,\
            colors='k',linestyles='solid')
plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
ax1.set_xlim(-90,90)
ax1.set_ylim(10,1000)
ax1.invert_yaxis()
ax1.set_yscale('log')
ax1.set_yticks([1000,800, 500, 300, 100, 10])
ax1.set_yticklabels((1000,800,500,300,100,10))
ax1.set_ylabel('Pressure (hPa)')
ax1.set_xlabel('Latitude')
ax1.set_xticks([-45, 0, 45])
ax1.set_xticklabels((-45, 0, 45))
ax1.set_title(ncfield + ' (' + units + ')')


#      Zoom on NH: shallow top (200hPa)
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
CS3 = plt.contour(lats,levs/100,plotfld,\
                 contsn,\
                 colors='k',linestyles='dashed')
plt.clabel(CS3,fmt = '%2.1f',inline=0,fontsize=10)
CS3 = plt.contour(lats,levs/100,plotfld,\
                 contsp,\
                 colors='k',linestyles='solid')
plt.clabel(CS3,fmt = '%2.1f',inline=1,fontsize=10)
ax3.set_xlim(0,90)
ax3.set_ylim(200,1000)
ax3.invert_yaxis()
ax3.set_yscale('log')
ax3.set_ylabel('Pressure (hPa)')
ax3.set_xlabel('Latitude')
ax3.set_yticks([1000,800, 500, 300, 200])
ax3.set_yticklabels((1000,800,500,300,200))
ax3.set_xticks([0, 30, 60, 80])
ax3.set_xticklabels((0, 30, 60, 80))
ax3.set_title(ncfield + ' (' + units + ')')
if printtofile:
    plt.savefig('MMC_' + timeavg + \
                '_' + casename + '_NH.pdf')
    plt.savefig('MMC_' + timeavg +\
                '_' + casename + '_NH.png')


#      Zoom on NH HIGH LAT (>30N): shallow top (200hPa)
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed' # @@?? isn't working! s/b default
fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
CS5 = plt.contour(lats,levs/100,plotfld,\
                 levels=contsnz,\
                 colors='k',linestyles='dashed')
plt.clabel(CS5,fmt = '%2.1f',inline=0,fontsize=10)
CS5 = plt.contour(lats,levs/100,plotfld,\
                 levels=contspz,\
                 colors='k',linestyles='solid')
plt.clabel(CS5,fmt = '%2.2f',inline=1,fontsize=10)
ax5.set_xlim(30,90)
ax5.set_ylim(200,1000)
ax5.invert_yaxis()
ax5.set_yscale('log')
ax5.set_ylabel('Pressure (hPa)')
ax5.set_xlabel('Latitude')
ax5.set_yticks([1000,800, 500, 300, 200])
ax5.set_yticklabels((1000,800,500,300,200))
ax5.set_xticks([30, 60, 80])
ax5.set_xticklabels((30, 60, 80))
ax5.set_title(ncfield + ' (' + units + ')')
if printtofile:
    plt.savefig('MMC_' + timeavg + \
                '_' + casename + '_NHzoom.pdf')
    plt.savefig('MMC_' + timeavg +\
                '_' + casename + '_NHzoom.png')

#   SEASONAL MEAN DIFF
#      full height of atmosphere
plotfld = np.average(seasfldp - seasfldc,0)
plotfld = plotfld*100
units = "10^8 kg/s"

fig = plt.figure()
ax = fig.add_subplot(111)
CS = plt.contour(lats,levs/100,plotfld,diffcontsn,\
                 colors='k',linestyles='dashed')
plt.clabel(CS,fmt = '%2.1f',inline=0,fontsize=10)
CS = plt.contour(lats,levs/100,plotfld,diffcontsp,\
                 colors='k',linestyles='solid')
plt.clabel(CS,fmt = '%2.1f',inline=1,fontsize=10)
ax.invert_yaxis()
ax.set_xlim(-90,90)
ax.set_yscale('log')
ax.set_ylabel('Pressure (hPa)')
ax.set_xlabel('Latitude')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_xticks([-45, 0, 45])
ax.set_xticklabels((-45, 0, 45))
ax.set_title(ncfield + ' (' + units + ')')
pc = cplt.addtsig(ax,pval,lat,lev/100,type='color')
cb = plt.colorbar(pc)# doesn't work? ,boundaries=(0,.05))


#      shallow top (100hPa)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
CS2 = plt.contour(lats,levs/100,plotfld,diffcontsn,\
                 colors='k',linestyles='dashed')
plt.clabel(CS2,fmt = '%2.1f',inline=0,fontsize=10)
CS2 = plt.contour(lats,levs/100,plotfld,diffcontsp,\
                 colors='k',linestyles='solid')
plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)
ax2.set_xlim(-90,90)
ax2.set_ylim(100,1000)
ax2.invert_yaxis()
ax2.set_yscale('log')
ax2.set_ylabel('Pressure (hPa)')
ax2.set_xlabel('Latitude')
ax2.set_yticks([1000,800, 500, 300, 100])
ax2.set_yticklabels((1000,800,500,300,100))
ax2.set_xticks([-45, 0, 45])
ax2.set_xticklabels((-45, 0, 45))
ax2.set_title(ncfield + ' (' + units + ')')
cplt.addtsig(ax2,pval,lat,lev/100,type='color')
if printtofile:
    plt.savefig('MMCsig_' + timeavg + '_' + casenamep +\
                '_v_' + casename + '_shall.pdf')
    plt.savefig('MMCsig_' + timeavg + '_' + casenamep +\
                '_v_' + casename + '_shall.png')


#      Zoom on NH: shallow top (200hPa)
fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
CS4 = plt.contour(lats,levs/100,plotfld,\
                 diffcontsn,\
                 colors='k',linestyles='dashed')
plt.clabel(CS4,fmt = '%2.1f',inline=0,fontsize=10)
CS4 = plt.contour(lats,levs/100,plotfld,\
                 diffcontsp,\
                 colors='k',linestyles='solid')
plt.clabel(CS4,fmt = '%2.1f',inline=1,fontsize=10)
ax4.set_xlim(0,90)
ax4.set_ylim(200,1000)
ax4.invert_yaxis()
ax4.set_yscale('log')
ax4.set_ylabel('Pressure (hPa)')
ax4.set_xlabel('Latitude')
ax4.set_yticks([1000,800, 500, 300, 200])
ax4.set_yticklabels((1000,800,500,300,200))
ax4.set_xticks([0, 30, 60, 80])
ax4.set_xticklabels((0, 30, 60, 80))
ax4.set_title(ncfield + ' (' + units + ')')
cplt.addtsig(ax4,pval,lat,lev/100,type='color')
if printtofile:
    plt.savefig('MMCsig_' + timeavg + '_' + casenamep +\
                '_v_' + casename + '_NH.pdf')
    plt.savefig('MMCsig_' + timeavg + '_' + casenamep +\
                '_v_' + casename + '_NH.png')


#      Zoom on NH HIGH LAT (>30N): shallow top (200hPa)
fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
CS6 = plt.contour(lats,levs/100,plotfld,\
                 diffcontsnz,\
                 colors='k',linestyles='dashed')
plt.clabel(CS6,fmt = '%2.1f',inline=0,fontsize=10)
CS6 = plt.contour(lats,levs/100,plotfld,\
                 diffcontspz,\
                 colors='k',linestyles='solid')
plt.clabel(CS6,fmt = '%2.1f',inline=1,fontsize=10)
ax6.set_xlim(30,90)
ax6.set_ylim(200,1000)
ax6.invert_yaxis()
ax6.set_yscale('log')
ax6.set_ylabel('Pressure (hPa)')
ax6.set_xlabel('Latitude')
ax6.set_yticks([1000,800, 500, 300, 200])
ax6.set_yticklabels((1000,800,500,300,200))
ax6.set_xticks([30, 60, 80])
ax6.set_xticklabels((30, 60, 80))
ax6.set_title(ncfield + ' (' + units + ')')
cplt.addtsig(ax6,pval,lat,lev/100,type='color')
if printtofile:
    plt.savefig('MMCsig_' + timeavg + '_' + casenamep +\
                '_v_' + casename + '_NHzoom.pdf')
    plt.savefig('MMCsig_' + timeavg + '_' + casenamep +\
                '_v_' + casename + '_NHzoom.png')




# TODO: zonal mean vprime*tprime
#       This eddy heat flux drives the indirect MMC (ferrell)
#       zonal mean vprime*uprime
#       This eddy momentum flux also drive ferrell
#       @@ stat sig on these
#       @@ save seasonal figures w/ stat sig. Add TEMP/U?

########################## END #######################@@
#     ALL MONTHS
months = con.get_mon()

midx=0
fig4,ax4 = plt.subplots(2,6) 
fig4.set_size_inches(12,6)
fig4.subplots_adjust(hspace=.15,wspace=.05)
for ax in ax4.flat:
    plotfld = fldpzm[midx,:,:] - fldczm[midx,:,:]
    pc = ax.pcolormesh(lats,levs/100,plotfld,\
                       cmap= plt.cm.get_cmap(cmap),shading='gouraud',\
                       vmin=cminm,vmax=cmaxm)
    ax.set_title(months[midx])
    ax.set_xlim(-90,90)
    ax.set_ylim(10,1000)
    ax.invert_yaxis()
    ax.set_yscale('log')
    if midx == 0 or midx == 6:
        ax.set_ylabel('Pressure (hPa)')
        ax.set_yticks([1000,800, 500, 300, 100, 10])
        ax.set_yticklabels((1000,800,500,300,100,10))
    else:
        ax.set_yticklabels('')

    if midx in range(6,12):
        ax.set_xlabel('Latitude')
        ax.set_xticks([-45, 0, 45])
        ax.set_xticklabels((-45, 0, 45))
    else:
        ax.set_xticklabels('')
    midx = midx+1

#cbar_ax = fig4.add_axes([.2, .02, .7, .03])
#fig4.colorbar(pc,cax=cbar_ax,orientation='horizontal')
cbar_ax = fig4.add_axes([.91,.15, .02,.7])
fig4.colorbar(pc,cax=cbar_ax)
plt.suptitle(ncfield + ' (' + units + ')')
if printtofile:
    plt.savefig(field + 'VERTzm_' + casenamep +\
                '_v_' + casename + '_shall_allmos.pdf')
    plt.savefig(field + 'VERTzm_' + casenamep +\
                '_v_' + casename + '_shall_allmos.png')


#     ALL MONTHS As Screen et al 2013, ClimDyn

midx=0
fig5,ax5 = plt.subplots(2,6) 
fig5.set_size_inches(12,4.5)
fig5.subplots_adjust(hspace=.15,wspace=.05)
for ax in ax5.flat:
    plotfld = fldpzm[midx,:,:] - fldczm[midx,:,:]
    pc = ax.pcolormesh(lats,levs/100,plotfld,\
                       cmap= plt.cm.get_cmap(cmap),shading='gouraud',\
                       vmin=cminsc,vmax=cmaxsc)
    ax.set_title(months[midx])
    ax.set_xlim(20,90)
    ax.set_ylim(300,1000)
    ax.invert_yaxis()
    #ax.set_yscale('log')
    if midx == 0 or midx == 6:
        ax.set_ylabel('Pressure (hPa)')
        ax.set_yticks([900,700, 500, 300])
        ax.set_yticklabels((900,700,500,300))
    else:
        ax.set_yticklabels('')

    if midx in range(6,12):
        ax.set_xlabel('Latitude')
        ax.set_xticks([40, 60, 80])
        ax.set_xticklabels((40,60,80))
    else:
        ax.set_xticklabels('')
    midx = midx+1

#cbar_ax = fig4.add_axes([.2, .02, .7, .03])
#fig4.colorbar(pc,cax=cbar_ax,orientation='horizontal')
cbar_ax = fig5.add_axes([.91,.15, .02,.7])
fig5.colorbar(pc,cax=cbar_ax)
plt.suptitle(ncfield + ' (' + units + ')')
if printtofile:
    plt.savefig(field + 'VERTzm_'+ casenamep +\
                '_v_' + casename + '_screen_allmos.pdf')
    plt.savefig(field + 'VERTzm_' + casenamep +\
                '_v_' + casename + '_screen_allmos.png')
