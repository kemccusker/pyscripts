# plotvert_canam4sims2.py
#    2/19/2014
#    Vertical plots for CanAM4 sea ice simulations
#    3/3/2014: taken from first version. Try and use new
#              cccmaNC module to read and average netcdf data
#    

#import numpy as np # for array handling
#import scipy as sp # scientific python
import scipy.stats
#import matplotlib.pyplot as plt # for basic plotting
import matplotlib.cm as cm
from subprocess import call # for doing system calls - not really needed
#from netCDF4 import Dataset
#from mpl_toolkits.basemap import Basemap # for maps
import datetime as datetime
import matplotlib.colors as col
import platform as platform
#import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import matplotlib.font_manager as fm
#import cccmaNC as cnc
import cccmacmaps

#import cdo as cdo; cdo = cdo.Cdo()
#import os

# https://code.zmaw.de/projects/cdo/wiki/Cdo%7Brbpy%7D

# while I'm still creating these modules, have to reload to get changes
## cplt = reload(cplt)
## con = reload(con)
## cutl = reload(cutl)
## cnc = reload(cnc)

#os.system('rm -rf /tmp/cdoPy*')
plt.close("all")
plt.ion()

printtofile=1
allmos=1 # make monthly figures
bimos=0  # bi-monthly figures
seasonal=0 # seasonal figures
singleplots=0  # seasonal climo and mean diff
addcontlines=1 # want contour lines in addition to colors on figs
obssims=1    # this will overrule the simulation settings and set to kemhad*
obssimscomp=0 # compare hadpert to pert2 @@ not ready yet

# # # ######## set Field info ###############
field = 'gz'  # t, u, gz


seasons = 'DJF','MAM','JJA','SON'

model = 'CanAM4'
ftype = 'ts'    # timeseries
## @@ really shouldn't use anymore now that code is all worked out #ftype = 'climo'  # 12-month climatology
season = 'ANN'
 
# # # ########### set Simulations #############
# Control run
casename = 'kemctl1'
timstr = '001-061'
timstr2 = '062-111'

# Pert run
casenamep1 = 'kem1pert1'  # 2002-2012 sic and sit
casenamep2 = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
casenamep3 = 'kem1pert3'  # 2002-2012 sic, adjusted sst. control sit
timstrp = '001-061'
timstrp2 = '062-111'

######## set pert run ############
casenamep = casenamep2

######## set new control run ###########
casename = casenamep3


########### OR SET OTHER RUNS ##############
if obssims:
    casename = 'kemhadctl'
    casenamep = 'kemhadpert'
    timstr2 = '062-121'
    timstrp2 = '062-121'
elif obssimscomp: # can't do this yet...need had runs to finish
    casename = 'kemhadpert'
    casenamep = 'kem1pert2'

#cmap = 'RdBu_r'
cmap = 'blue2red_20'
cmapclimo = 'Spectral_r'
conv = 1   # conversion factor to convert units, etc


if field == 't':
    ncfield = 'TEMP'
    units = 'K' # @@
    cminc = 200; cmaxc = 310
    cmin = -.5; cmax = .5
    cminm = -.8; cmaxm = .8  # for monthly
    cminsc = -2.5; cmaxsc = 2.5  # as Screen et al paper
    cminp = -.5; cmaxp=.5 # for when a pert is 'control'
elif field == 'u':
    ncfield = 'U'
    units = 'm/s' #@@
    cminc = -40; cmaxc = 40
    cmin = -.5; cmax = .5
    cminm = -1; cmaxm = 1
    cminsc=cminm; cmaxsc=cmaxm
    cminp = -.5; cmaxp = .5
elif field == 'gz':
    ncfield = 'PHI'
    units = 'm' # @@
    conv = 1/con.get_g()
    cminc = 0; cmaxc = 48000  # climo
    cmin = -10; cmax = 10  # annual mean
    cminm = -15; cmaxm = 15  # for monthly
    cminsc = -25; cmaxsc = 25  # as Screen et al paper
    cminp = -10; cmaxp=10
else:
    print 'no such field: ' + field


# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'

# The 3D data is split into 2 timeseries
fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_' + ftype + '.nc'
fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_' + ftype + '.nc'
## fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp + '_' + ftype + '.nc'
## fnamep3 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp + '_' + ftype + '.nc'

fnamec2 = basepath + casename + subdir + casename + '_' + field + '_' + timstr2 + '_' + ftype + '.nc'
fnamep2 = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp2 + '_' + ftype + '.nc'
## fnamep22 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp2 + '_' + ftype + '.nc'
## fnamep32 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp2 + '_' + ftype + '.nc'



## ############################
## ########## set pert! #######
## #fnamep = fnamep2; fnamep2 = fnamep22; casenamep = casenamep2
## #fnamep = fnamep1; fnamep2 = fnamep12; casenamep = casenamep1
## fnamep = fnamep3; fnamep2 = fnamep32; casenamep = casenamep3

print 'CONTROL: ' + casename
print 'PERT: ' + casenamep

if casename != 'kemctl1' and casename != 'kemhadctl':
    cmin=cminp; cmax=cmaxp
    cminm=cminp; cmaxm=cmaxp
    cminsc=cminp; cmaxsc=cmaxp

lat = cnc.getNCvar(fnamec,'lat')
lev = cnc.getNCvar(fnamec,'plev')

# get the data already zonal-meaned
#fldczm = cnc.getNCvar(fnamec,ncfield,calc='zm',timechunk=((styr-1)*12,enyr*12+1) )*conv
#fldpzm = cnc.getNCvar(fnamep2,ncfield,calc='zm',timechunk=((styr-1)*12,enyr*12+1) )*conv

seasfldczm = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',seas=season,calc='zm')*conv,
               cnc.getNCvar(fnamec2,ncfield,seas=season,calc='zm')*conv,
               axis=0)
seasfldpzm = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',seas=season,calc='zm')*conv,
               cnc.getNCvar(fnamep2,ncfield,seas=season,calc='zm')*conv,
               axis=0)



#nt,nlev,nlat,nlon = fldc.shape # if nt == 12 then it's a climo
nt,nlev,nlat = seasfldczm.shape
print nt


# seasonalize
#seasfldczm = cutl.seasonalize_monthlyts(fldczm,timeavg)
#seasfldpzm = cutl.seasonalize_monthlyts(fldpzm,timeavg)

if ftype=='ts':
    tstat,pval = sp.stats.ttest_ind(seasfldpzm,seasfldczm,axis=0)

# time-mean
seasfldczmtm=np.mean(seasfldczm,0)
seasfldpzmtm=np.mean(seasfldpzm,0)


lats,levs = np.meshgrid(lat,lev)


############################## SEASONAL MEAN ##########################
if singleplots:
    #   SEASONAL MEAN CLIMO (CONTROL)
    cmlen=15
    incr = (cmaxc-cminc) / (cmlen)
    conts = np.arange(cminc,cmaxc+incr,incr)

    plotfld = seasfldczmtm

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ## plt.pcolormesh(lats,levs/100,plotfld,\
    ##                cmap= plt.cm.get_cmap(cmapclimo),shading='gouraud',\
    ##                vmin=cminc,vmax=cmaxc)
    cf = plt.contourf(lats,levs/100,plotfld,
                 cmap= plt.cm.get_cmap(cmapclimo),levels=conts,
                 vmin=cminc,vmax=cmaxc,extend='both')
    if addcontlines:
        plt.contour(lats,levs/100,plotfld,levels=conts,colors='.3')
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
    cbar = fig1.colorbar(cf)

    if printtofile:
        plt.savefig(field + 'VERTzm_' + season + \
                    '_' + casename + '_shall2.pdf')
        plt.savefig(field + 'VERTzm_' + season +\
                    '_' + casename + '_shall2.png')



    #   SEASONAL MEAN DIFF
    #      full height of atmosphere
    cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    plotfld = seasfldpzmtm - seasfldczmtm

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ## pc = plt.pcolormesh(lats,levs/100,plotfld,\
    ##                     cmap= plt.cm.get_cmap(cmap),shading='gouraud',\
    ##                     vmin=cmin,vmax=cmax)
    pc = plt.contourf(lats,levs/100,plotfld,
                      cmap=plt.cm.get_cmap(cmap),levels=conts,
                      vmin=cmin,vmax=cmax,extend='both')
    if addcontlines:
        plt.contour(lats,levs/100,plotfld,levels=conts,colors='.3')
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
    cbar = fig.colorbar(pc)
    if ftype == 'ts':
        # add sig
        cplt.addtsig(ax,pval,lat,lev/100,type='hatch')

    if printtofile:
        if ftype=='ts':
            plt.savefig(field + 'VERTzmdiff_' + season + '_' + casenamep + '_v_' + casename + '2.pdf')
            plt.savefig(field + 'VERTzmdiffsig_' + season + '_' + casenamep + '_v_' + casename + '2.png')
        else:
            plt.savefig(field + 'VERTzmdiff_' + season + '_' + casenamep + '_v_' + casename + '2.pdf')
            plt.savefig(field + 'VERTzmdiff_' + season + '_' + casenamep + '_v_' + casename + '2.png')



    #      shallow top (100hPa)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ## plt.pcolormesh(lats,levs/100,plotfld,\
    ##                cmap= plt.cm.get_cmap(cmap),shading='gouraud',\
    ##                vmin=cmin,vmax=cmax)
    cf = plt.contourf(lats,levs/100,plotfld,
                 cmap=plt.cm.get_cmap(cmap),levels=conts,
                 vmin=cmin,vmax=cmax,extend='both')
    if addcontlines:
        plt.contour(lats,levs/100,plotfld,levels=conts,colors='.3')
    ax2.set_xlim(-90,90)
    ax2.set_ylim(10,1000)
    ax2.invert_yaxis()
    ax2.set_yscale('log')
    ax2.set_ylabel('Pressure (hPa)')
    ax2.set_xlabel('Latitude')
    ax2.set_yticks([1000,800, 500, 300, 100, 10])
    ax2.set_yticklabels((1000,800,500,300,100,10))
    ax2.set_xticks([-45, 0, 45])
    ax2.set_xticklabels((-45, 0, 45))
    ax2.set_title(ncfield + ' (' + units + ')')
    cbar = fig2.colorbar(cf)
    if ftype == 'ts':
        # add sig
        cplt.addtsig(ax2,pval,lat,lev/100,type='hatch')

    if printtofile:
        if ftype=='ts':
            plt.savefig(field + 'VERTzmdiff_' + season + '_' + casenamep +\
                        '_v_' + casename + '_shall2.pdf') # sig doesn't plot to pdf
            plt.savefig(field + 'VERTzmdiffsig_' + season + '_' + casenamep +\
                        '_v_' + casename + '_shall2.png')
        else:
            plt.savefig(field + 'VERTzmdiff_' + season + '_' + casenamep +\
                        '_v_' + casename + '_shall2.pdf')
            plt.savefig(field + 'VERTzmdiff_' + season + '_' + casenamep +\
                        '_v_' + casename + '_shall2.png')

######################## ALL MONTHS #######################
if allmos:
    #     ALL MONTHS
    cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
    incr = (cmaxm-cminm) / (cmlen)
    conts = np.arange(cminm,cmaxm+incr,incr)

    months = con.get_mon()
    # these are suitable for plotting, but not tstat calcs. they come back as climos
    ## fldczm = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',seas='climo',calc='zm')*conv,
    ##                cnc.getNCvar(fnamec2,ncfield,seas='climo',calc='zm')*conv,
    ##                axis=0)
    ## fldpzm = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',seas='climo',calc='zm')*conv,
    ##                cnc.getNCvar(fnamep2,ncfield,seas='climo',calc='zm')*conv,
    ##                axis=0)

    tstat = np.zeros((12,nlev,nlat))
    pval = np.zeros((12,nlev,nlat))
    fldczmmos = np.zeros((12,nlev,nlat))
    fldpzmmos = np.zeros((12,nlev,nlat))

    midx=0
    fig4,ax4 = plt.subplots(2,6) 
    fig4.set_size_inches(12,6)
    fig4.subplots_adjust(hspace=.15,wspace=.05)
    for ax in ax4.flat:

        fldczmmo = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',seas=midx+1,calc='zm')*conv,
                         cnc.getNCvar(fnamec2,ncfield,seas=midx+1,calc='zm')*conv,
                         axis=0)
        fldpzmmo = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',seas=midx+1,calc='zm')*conv,
                         cnc.getNCvar(fnamep2,ncfield,seas=midx+1,calc='zm')*conv,
                         axis=0)
        tstat[midx,:,:],pval[midx,:,:] = sp.stats.ttest_ind(fldpzmmo,fldczmmo,axis=0)

        fldczmmos[midx,:,:] = np.mean(fldczmmo,axis=0)
        fldpzmmos[midx,:,:] = np.mean(fldpzmmo,axis=0)
        
        # plotfld = fldpzm[midx,:,:] - fldczm[midx,:,:]
        plotfld = fldpzmmos[midx,:,:] - fldczmmos[midx,:,:]

        ## pc = ax.pcolormesh(lats,levs/100,plotfld,\
        ##                    cmap= plt.cm.get_cmap(cmap),shading='gouraud',\
        ##                    vmin=cminm,vmax=cmaxm)
        pc = ax.contourf(lats,levs/100,plotfld,
                         cmap=plt.cm.get_cmap(cmap),levels=conts,
                         vmin=cminm,vmax=cmaxm,extend='both')
        if addcontlines:
            ax.contour(lats,levs/100,plotfld,levels=conts,colors='.3')
        if ftype == 'ts':
            # add sig
            cplt.addtsig(ax,pval[midx,:,:],lat,lev/100,type='hatch')

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
        if ftype=='ts':
            plt.savefig(field + 'VERTzmdiff_' + casenamep +\
                        '_v_' + casename + '_shall_allmos2.pdf') # sig doesn't plot to pdf
            plt.savefig(field + 'VERTzmdiffsig_' + casenamep +\
                        '_v_' + casename + '_shall_allmos2.png')
        else:
            plt.savefig(field + 'VERTzmdiff_' + casenamep +\
                        '_v_' + casename + '_shall_allmos2.pdf')
            plt.savefig(field + 'VERTzmdiff_' + casenamep +\
                        '_v_' + casename + '_shall_allmos2.png')


    #     ALL MONTHS As Screen et al 2013, ClimDyn
    incr = (cmaxsc-cminsc) / (cmlen)
    conts = np.arange(cminsc,cmaxsc+incr,incr)

    midx=0
    fig5,ax5 = plt.subplots(2,6) 
    fig5.set_size_inches(12,4.5)
    fig5.subplots_adjust(hspace=.15,wspace=.05)
    for ax in ax5.flat:
        plotfld = fldpzmmos[midx,:,:] - fldczmmos[midx,:,:]
        
        ## pc = ax.pcolormesh(lats,levs/100,plotfld,\
        ##                    cmap= plt.cm.get_cmap(cmap),shading='gouraud',\
        ##                    vmin=cminsc,vmax=cmaxsc)
        pc = ax.contourf(lats,levs/100,plotfld,
                         cmap=plt.cm.get_cmap(cmap),levels=conts,
                         vmin=cminsc,vmax=cmaxsc, extend='both')
        if addcontlines:
            ax.contour(lats,levs/100,plotfld,levels=conts,colors='.3')
        if ftype == 'ts':
            # add sig
            cplt.addtsig(ax,pval[midx,:,:],lat,lev/100,type='hatch')

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
        if ftype == 'ts':
            plt.savefig(field + 'VERTzmdiff_'+ casenamep +\
                        '_v_' + casename + '_screen_allmos2.pdf') # sig doesn't plot to pdf
            plt.savefig(field + 'VERTzmdiffsig_' + casenamep +\
                        '_v_' + casename + '_screen_allmos2.png')
        else:
            plt.savefig(field + 'VERTzmdiff_'+ casenamep +\
                        '_v_' + casename + '_screen_allmos2.pdf')
            plt.savefig(field + 'VERTzmdiff_' + casenamep +\
                        '_v_' + casename + '_screen_allmos2.png')



if bimos:
    print 'bimos'


if seasonal:

    tstat = np.zeros((len(seasons),nlev,nlat))
    pval = np.zeros((len(seasons),nlev,nlat))
    fldczmallseas = np.zeros((len(seasons),nlev,nlat))
    fldpzmallseas = np.zeros((len(seasons),nlev,nlat))
    
    incr = (cmaxsc-cminsc) / (cmlen)
    conts = np.arange(cminsc,cmaxsc+incr,incr)

    midx=0
    fig6,ax6 = plt.subplots(1,4) 
    fig6.set_size_inches(12,2.5)
    fig6.subplots_adjust(hspace=.15,wspace=.05)
    for ax in ax6.flat:


        fldczmsea = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',
                                       seas=seasons[midx],calc='zm')*conv,
                          cnc.getNCvar(fnamec2,ncfield,seas=seasons[midx],calc='zm')*conv,
                          axis=0)
        fldpzmsea = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',
                                       seas=seasons[midx],calc='zm')*conv,
                          cnc.getNCvar(fnamep2,ncfield,seas=seasons[midx],calc='zm')*conv,
                          axis=0)
        
        tstat[midx,:,:],pval[midx,:,:] = sp.stats.ttest_ind(fldpzmsea,fldczmsea,axis=0)
        fldczmallseas[midx,:,:] = np.mean(fldczmsea,axis=0)
        fldpzmallseas[midx,:,:] = np.mean(fldpzmsea,axis=0)
        
        plotfld = fldpzmallseas[midx,:,:] - fldczmallseas[midx,:,:]

        pc = ax.contourf(lats,levs/100,plotfld,
                         cmap=plt.cm.get_cmap(cmap),levels=conts,
                         vmin=cminsc,vmax=cmaxsc, extend='both')
        if addcontlines:
            ax.contour(lats,levs/100,plotfld,levels=conts,colors='.3')
        if ftype == 'ts':
            # add sig
            cplt.addtsig(ax,pval[midx,:,:],lat,lev/100,type='hatch')

        ax.set_title(seasons[midx])
        ax.set_xlim(20,90)
        ax.set_ylim(300,1000)
        ax.invert_yaxis()
        #ax.set_yscale('log')
        if midx == 0:
            ax.set_ylabel('Pressure (hPa)')
            ax.set_yticks([900,700, 500, 300])
            ax.set_yticklabels((900,700,500,300))
        else:
            ax.set_yticklabels('')

        ax.set_xlabel('Latitude')
        ax.set_xticks([40, 60, 80])
        ax.set_xticklabels((40,60,80))

        midx = midx+1

    #cbar_ax = fig4.add_axes([.2, .02, .7, .03])
    #fig4.colorbar(pc,cax=cbar_ax,orientation='horizontal')
    cbar_ax = fig6.add_axes([.91,.15, .02,.7])
    fig6.colorbar(pc,cax=cbar_ax)
    plt.suptitle(ncfield + ' (' + units + ')')
    if printtofile:
        if ftype == 'ts':
            plt.savefig(field + 'VERTzmdiff_'+ casenamep +\
                        '_v_' + casename + '_screen_seas2.pdf') # sig doesn't plot to pdf
            plt.savefig(field + 'VERTzmdiffsig_' + casenamep +\
                        '_v_' + casename + '_screen_seas2.png')
        else:
            plt.savefig(field + 'VERTzmdiff_'+ casenamep +\
                        '_v_' + casename + '_screen_seas2.pdf')
            plt.savefig(field + 'VERTzmdiff_' + casenamep +\
                        '_v_' + casename + '_screen_seas2.png')


    incr = (cmaxm-cminm) / (cmlen)
    conts = np.arange(cminm,cmaxm+incr,incr)

    midx=0
    fig6,ax6 = plt.subplots(1,4) 
    fig6.set_size_inches(12,3.5)
    fig6.subplots_adjust(hspace=.15,wspace=.05)
    for ax in ax6.flat:


        ## fldczmsea[midx,:,:] = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',
        ##                                              seas=seasons[midx],calc='zm')*conv,
        ##                      cnc.getNCvar(fnamec2,ncfield,seas=seasons[midx],calc='zm')*conv,
        ##                      axis=0)
        ## fldpzmsea[midx,:,:] = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',
        ##                                              seas=seasons[midx],calc='zm')*conv,
        ##                      cnc.getNCvar(fnamep2,ncfield,seas=seasons[midx],calc='zm')*conv,
        ##                      axis=0)
        ## tstat[midx,:,:],pval[midx,:,:] = sp.stats.ttest_ind(fldpzmsea,fldczmsea,axis=0)

        plotfld = fldpzmallseas[midx,:,:] - fldczmallseas[midx,:,:]

        pc = ax.contourf(lats,levs/100,plotfld,
                         cmap=plt.cm.get_cmap(cmap),levels=conts,
                         vmin=cminm,vmax=cmaxm, extend='both')
        if addcontlines:
            ax.contour(lats,levs/100,plotfld,levels=conts,colors='.3')
        if ftype == 'ts':
            # add sig
            cplt.addtsig(ax,pval[midx,:,:],lat,lev/100,type='hatch')

        ax.set_title(seasons[midx])
        ax.set_xlim(-90,90)
        ax.set_ylim(10,1000)
        ax.invert_yaxis()
        #ax.set_yscale('log')
        if midx == 0:
            ax.set_ylabel('Pressure (hPa)')
            ax.set_yticks([1000,800, 500, 300, 100, 10])
            ax.set_yticklabels((1000,800,500,300,100,10))
        else:
            ax.set_yticklabels('')

        ax.set_xlabel('Latitude')
        ax.set_xticks([-45, 0, 45])
        ax.set_xticklabels((-45, 0, 45))

        midx = midx+1

    #cbar_ax = fig4.add_axes([.2, .02, .7, .03])
    #fig4.colorbar(pc,cax=cbar_ax,orientation='horizontal')
    cbar_ax = fig6.add_axes([.91,.15, .02,.7])
    fig6.colorbar(pc,cax=cbar_ax)
    plt.suptitle(ncfield + ' (' + units + ')')
    if printtofile:
        if ftype == 'ts':
            plt.savefig(field + 'VERTzmdiff_'+ casenamep +\
                        '_v_' + casename + '_shall_seas2.pdf') # sig doesn't plot to pdf
            plt.savefig(field + 'VERTzmdiffsig_' + casenamep +\
                        '_v_' + casename + '_shall_seas2.png')
        else:
            plt.savefig(field + 'VERTzmdiff_'+ casenamep +\
                        '_v_' + casename + '_shall_seas2.pdf')
            plt.savefig(field + 'VERTzmdiff_' + casenamep +\
                        '_v_' + casename + '_shall_seas2.png')



    
