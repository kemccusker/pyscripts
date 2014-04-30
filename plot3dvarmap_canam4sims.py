""" plot3dvarmap_canam4sims.py
#    4/8/2014
#    Taken from plotvert_canam4sims2.py
#    This script will plot a map of a level of a 3d atmos var
#
#    2/19/2014
#    Vertical plots for CanAM4 sea ice simulations
#    3/3/2014: taken from first version. Try and use new
#              cccmaNC module to read and average netcdf data
"""   

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
#import cccmaplots as cplt # in startup now: ~/.config/ipython/profile_default/startup/00-startup.py
import constants as con
import cccmautils as cutl
import matplotlib.font_manager as fm
#import cccmaNC as cnc # in startup
#import cccmacmaps # in startup

#import cdo as cdo; cdo = cdo.Cdo()
#import os

# https://code.zmaw.de/projects/cdo/wiki/Cdo%7Brbpy%7D

# while I'm still creating these modules, have to reload to get changes
# don't need to do anymore b/c added to configuration (autoreload)
#cplt = reload(cplt)
#con = reload(con)
#cutl = reload(cutl)
#cnc = reload(cnc)

#os.system('rm -rf /tmp/cdoPy*')
plt.close("all")
plt.ion()

printtofile=1
allmos=1 # make monthly figures
bimos=0  # bi-monthly figures. not implemented yet 4/29/14
seasonal=1 # seasonal figures
singleplots=1  # seasonal climo and mean diff
obssims=1    # this will overrule the simulation settings and set to kemhad*
thickness=0; level2=70000 # calc thickness between level and level2 for GZ only, typically 1000-700hPa
#level = 100000 # for thickness calc

#level = 30000
#level = 50000 # 500hPa
level = 70000

sigtype = 'cont' # significance: 'cont' or 'hatch' which is default

"""  plev = 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, 12500, 
    15000, 17500, 20000, 22500, 25000, 30000, 35000, 40000, 45000, 50000, 
    55000, 60000, 65000, 70000, 75000, 77500, 80000, 82500, 85000, 87500, 
    90000, 92500, 95000, 97500, 100000 ;
"""

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
casenamep = casenamep1

########### OR SET OTHER RUNS ##############
if obssims:
    casename = 'kemhadctl'
    casenamep = 'kemhadpert'
    timstr2 = '062-121'
    timstrp2 = '062-121'


cmap = 'blue2red_w20'
cmapclimo = 'Spectral_r'
conv = 1   # conversion factor to convert units, etc

# # # ######## set Field info ###############
field = 'gz'  # t, u, gz


if sigtype=='cont':
    suff='pdf'
else:
    suff='png'
    
if field == 't':
    ncfield = 'TEMP'
    units = 'K' # @@
    if level == 30000:
        cminc = 215; cmaxc = 245        
    elif level == 70000:
        cminc = 245; cmaxc = 285
        
    if level == 30000:
        cmin = -.3; cmax = .3
        cminm = -.8; cmaxm = .8  # for monthly
        cminsea = -.5; cmaxsea = .5
    elif level == 70000:
        cmin = -.3; cmax = .3
        cminm = -1; cmaxm = 1  # for monthly
        cminsea = -.5; cmaxsea = .5
    
elif field == 'u':
    ncfield = 'U'
    units = 'm/s' #@@
    if level==50000:
        cminc=-25; cmaxc=25
    elif level==70000:
        cminc=-15; cmaxc=15
    elif level == 30000:
        cminc=-40; cmaxc=40

    if level == 30000:
        cmin = -2; cmax = 2
        cminm = -5; cmaxm = 5
        cminsea = -3; cmaxsea = 3
    else:
        cmin = -1; cmax = 1
        cminm = -3; cmaxm = 3
        cminsea = -1; cmaxsea = 1

    cmapclimo='blue2red_20'
elif field == 'gz':
    ncfield = 'PHI'
    units = 'm' # @@
    conv = 1/con.get_g()
    if level==50000:
        cminc = 5200; cmaxc = 5900  # climo 500hPa
    elif level==70000:
        cminc=2800; cmaxc = 3200
    elif level==30000:
        cminc=8600; cmaxc = 9800
    elif level==100000: # use this for thickness calc 1000-700
        cminc=2650; cmaxc = 3050
        
    cmin = -8 # annual mean
    cmax = 8  # annual mean
    
    if level==30000:
        cmin = -15; cmax = 15
        cminsea=-20; cmaxsea = 20
        cminm = -30; cmaxm = 30  # for monthly
    else:
        cminsea = -15; cmaxsea = 15
        cminm = -20; cmaxm = 20  # for monthly

        
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

print field + ' level ' + str(level/100)
if field=='gz' and thickness==1: print ' to level ' + str(level2/100)
print 'CONTROL: ' + casename
print 'PERT: ' + casenamep

lat = cnc.getNCvar(fnamec,'lat')
lon = cnc.getNCvar(fnamec,'lon')


seasfldc = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',levsel=level,seas=season)*conv,
               cnc.getNCvar(fnamec2,ncfield,levsel=level,seas=season)*conv,
               axis=0)
seasfldp = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',levsel=level,seas=season)*conv,
               cnc.getNCvar(fnamep2,ncfield,levsel=level,seas=season)*conv,
               axis=0)
if field=='gz' and thickness==1:
    seasfldc2 = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',levsel=level2,seas=season)*conv,
                   cnc.getNCvar(fnamec2,ncfield,levsel=level2,seas=season)*conv,
                   axis=0)
    seasfldp2 = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',levsel=level2,seas=season)*conv,
                   cnc.getNCvar(fnamep2,ncfield,levsel=level2,seas=season)*conv,
                   axis=0)
    seasfldc=seasfldc2 - seasfldc
    seasfldp=seasfldp2 - seasfldp



#nt,nlev,nlat,nlon = fldc.shape # if nt == 12 then it's a climo
nt,nlat,nlon = seasfldc.shape
print nt


# seasonalize
#seasfldczm = cutl.seasonalize_monthlyts(fldczm,timeavg)
#seasfldpzm = cutl.seasonalize_monthlyts(fldpzm,timeavg)

if ftype=='ts':
    tstat,pval = sp.stats.ttest_ind(seasfldp,seasfldc,axis=0)

# time-mean
seasfldctm=np.mean(seasfldc,0)
seasfldptm=np.mean(seasfldp,0)


############################## SEASONAL (ANN) MEAN ##########################
if singleplots:
    #   SEASONAL MEAN CLIMO (CONTROL)
    if field=='gz' and thickness==1:
        title = season + ' ' + field + ' ' + str(level2/100) + '-' + str(level/100) + ': ' + casename
    else:
        title = season + ' ' + field + ' ' + str(level/100) + ': ' + casename

    plotfld = seasfldctm

    fig1 = plt.figure()
    bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cminc,cmax=cmaxc,cmap=cmapclimo,type='nh',\
                        title=title,units=units)
    if printtofile:
        if field=='gz' and thickness==1:
            fig1.savefig(field + 'THK' + str(level2/100) + '-' + str(level/100) + '_' + season + \
                         '_' + casename + '_nh.' + suff)
        else:
            fig1.savefig(field + 'LEV' + str(level/100) + '_' + season + \
                         '_' + casename + '_nh.' + suff)


    #   SEASONAL MEAN DIFF
    #      full height of atmosphere
    plotfld = seasfldptm - seasfldctm

    fig = plt.figure()
    bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,type='nh',\
                        title=title,units=units)
    cplt.addtsigm(bm,pval,lat,lon,type=sigtype)    

    if printtofile:
        if field=='gz' and thickness==1:
            fig.savefig(field + 'DIFFTHK' + str(level2/100) + '-' + str(level/100) + sigtype + '_' + season + \
                        '_' + casenamep + '_v_' + casename + '_nh.' + suff)
        else:
            fig.savefig(field + 'DIFFLEV' + str(level/100) + sigtype + '_' + season + \
                        '_' + casenamep + '_v_' + casename + '_nh.' + suff)



######################## ALL MONTHS #######################
if allmos:

    fldc = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',levsel=level)*conv,
               cnc.getNCvar(fnamec2,ncfield,levsel=level)*conv,
               axis=0)
    fldp = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',levsel=level)*conv,
               cnc.getNCvar(fnamep2,ncfield,levsel=level)*conv,
               axis=0)
    if field == 'gz' and thickness==1:
        fldc2 = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',levsel=level2)*conv,
                   cnc.getNCvar(fnamec2,ncfield,levsel=level2)*conv,
                   axis=0)
        fldp2 = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',levsel=level2)*conv,
                   cnc.getNCvar(fnamep2,ncfield,levsel=level2)*conv,
                   axis=0)
        fldc = fldc2 - fldc
        fldp = fldp2 - fldp
    
    #     ALL MONTHS
    months = con.get_mon()

    
#    sigs = np.zeros((12,nlat,nlon))


    midx=0
    fig4,ax4 = plt.subplots(2,6) 
    fig4.set_size_inches(12,6)
    fig4.subplots_adjust(hspace=.15,wspace=.05)
    for ax in ax4.flat:
        
        monfldc = fldc[midx::12,:,:]
        monfldp = fldp[midx::12,:,:]

        tstat,pval = sp.stats.ttest_ind(monfldp,monfldc,axis=0)
        #sigs[midx,:,:] = ma.masked_where(pval>siglevel,pval)

        plotfld = np.mean(monfldp-monfldc,0)

        bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',\
                     title=months[midx],axis=ax,suppcb=1)
        ax.set_title(months[midx])
        cplt.addtsigm(bm,pval,lat,lon,type=sigtype)

        midx = midx+1

    #cbar_ax = fig4.add_axes([.2, .02, .7, .03])
    #fig4.colorbar(pc,cax=cbar_ax,orientation='horizontal')
    cbar_ax = fig4.add_axes([.91,.25, .02,.5])
    fig4.colorbar(pc,cax=cbar_ax)
    if field=='gz' and thickness==1:
        plt.suptitle(ncfield + ' ' + str(level2/100) + '-' + str(level/100)+ ' (' + units + ')')
    else:
        plt.suptitle(ncfield + ' ' + str(level/100) + ' (' + units + ')')
    if printtofile:
        if field=='gz' and thickness==1:
            fig4.savefig(field + 'DIFFTHK' + str(level2/100) + '-' + str(level/100) + sigtype + '_' + casenamep +\
                    '_v_' + casename + '_allmos_nh.' + suff)
        else:
            fig4.savefig(field + 'DIFFLEV' + str(level/100) + sigtype + '_' + casenamep +\
                    '_v_' + casename + '_allmos_nh.' + suff)


if bimos:
    print 'bimos'

if seasonal:

    tstat = np.zeros((len(seasons),nlat,nlon))
    pval = np.zeros((len(seasons),nlat,nlon))
    fldcallseas = np.zeros((len(seasons),nlat,nlon))
    fldpallseas = np.zeros((len(seasons),nlat,nlon))
    
    midx=0
    fig6,ax6 = plt.subplots(1,4) 
    fig6.set_size_inches(12,2.5)
    fig6.subplots_adjust(hspace=.15,wspace=.05)
    for ax in ax6.flat:

        fldcsea = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',
                                       levsel=level,seas=seasons[midx])*conv,
                          cnc.getNCvar(fnamec2,ncfield,levsel=level,seas=seasons[midx])*conv,
                          axis=0)
        fldpsea = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',
                                       levsel=level,seas=seasons[midx])*conv,
                          cnc.getNCvar(fnamep2,ncfield,levsel=level,seas=seasons[midx])*conv,
                          axis=0)

        if field == 'gz' and thickness==1:
            fldcsea2 = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0061-12-31',levsel=level2,seas=seasons[midx])*conv,
                       cnc.getNCvar(fnamec2,ncfield,levsel=level2,seas=seasons[midx])*conv,
                       axis=0)
            fldpsea2 = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0061-12-31',levsel=level2,seas=seasons[midx])*conv,
                       cnc.getNCvar(fnamep2,ncfield,levsel=level2,seas=seasons[midx])*conv,
                       axis=0)
            fldcsea = fldcsea2 - fldcsea
            fldpsea = fldpsea2 - fldpsea

        tstat[midx,:,:],pval[midx,:,:] = sp.stats.ttest_ind(fldpsea,fldcsea,axis=0)
        fldcallseas[midx,:,:] = np.mean(fldcsea,axis=0)
        fldpallseas[midx,:,:] = np.mean(fldpsea,axis=0)
        
        plotfld = fldpallseas[midx,:,:] - fldcallseas[midx,:,:]

        bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cminsea,cmax=cmaxsea,cmap=cmap,type='nh',\
                        axis=ax,suppcb=1)
        ax.set_title(seasons[midx])
        cplt.addtsigm(bm,pval[midx,:,:],lat,lon,type=sigtype)

        midx = midx+1

    #cbar_ax = fig4.add_axes([.2, .02, .7, .03])
    #fig4.colorbar(pc,cax=cbar_ax,orientation='horizontal')
    cbar_ax = fig6.add_axes([.91,.25, .02,.5])
    fig6.colorbar(pc,cax=cbar_ax)
    #plt.suptitle(ncfield + ' (' + units + ')')
    if printtofile:
        if field=='gz' and thickness==1:
            fig6.savefig(field + 'DIFFTHK' + str(level2/100) + '-' + str(level/100) + sigtype + '_' + casenamep +\
                         '_v_' + casename + '_seas_nh.' + suff)
        else:
            fig6.savefig(field + 'DIFFLEV' + str(level/100) + sigtype + '_' + casenamep +\
                         '_v_' + casename + '_seas_nh.' + suff)
 
