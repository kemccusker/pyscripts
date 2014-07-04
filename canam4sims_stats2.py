"""
 canam4sims_stats.py
    2/20/2014: taken from plot_canam4sims_hists.py: 
               calculate & plot statistical properties of the runs

"""

#import numpy as np # for array handling
import numpy.ma as ma
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
#import cccmaplots as cplt    # my module
import constants as con      # my module
import cccmautils as cutl    # my module
import matplotlib.font_manager as fm
import copy
#import cccmacmaps as ccm
#import cccmaNC as cnc

# while I'm still creating these modules, have to reload to get changes
## cplt = reload(cplt)
## con = reload(con)
## cutl = reload(cutl)
## ccm = reload(ccm)
## cnc = reload(cnc)

plt.close("all")
plt.ion()

printtofile=1
plotann=0    # annual average
plotallmos=0 # each month separately
bimos=0 # averages every 2 mos (JF, MA, MJ, JA, SO, ND) @@ add
seasonal=1 # averages seasons (DJF, MAM, JJA, SON)
obssims=0  # override settings to do observed runs (kemhad*)

sigtype = 'cont' # significance: 'cont' or 'hatch' which is default

seasons = 'SON','DJF','MAM','JJA'

model = 'CanAM4'

# # # ########### set Simulations #############
# Control run
casename = 'kemctl1'
timstr = '001-121'
timesel = '0002-01-01,0121-12-31'

# Pert run
#casenamep1 = 'kem1pert1'  # 2002-2012 sic and sit
casenamep1 = 'kem1pert1b'  # 2002-2012 sic and sit

casenamep2 = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
casenamep3 = 'kem1pert3'  # 2002-2012 sic, adjusted sst. control sit
timstrp = '001-121'
casenamepra = 'kem1rcp85a' # 2022-2032 sic, adjusted sst, sit from RCP8.5



####### SET PERT RUN ############
casenamep = casenamep3
####### SET NEW CTL RUN #########
#casename = casenamep3

########## SET OBSERVED RUNS #######
if obssims==1:
    casename = 'kemhadctl'
    casenamep = 'kemhadpert'
    timstr = '001-121'
    timstrp = timstr
    timesel = '0002-01-01,0121-12-31'

    
print 'CONTROL IS ' + casename
print 'PERT IS ' + casenamep


# # # ######## set Field info ###################
# st, sicn, sic, gt, pmsl, pcp, hfl, hfs, turb, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)
field = 'st'

cmap = 'blue2red_w20' # default cmap
cmapclimo = 'Spectral_r'
pct = 0 # if 1, do calculation as a percent
skipnmin=1 # if 1, do not calc and plot Nmin




# # # ###########################################
#   Shouldn't have to mod below....

if field == 'st':
    units = 'K'
    conv = 1  # no conversion
    cmin = -2; cmax = 2  # for anomaly plots
    cminp=-.5; cmaxp=.5 # for when pert is 'ctl'
    cminm = -3; cmaxm = 3   # monthly
    ## print 'small clim!'
    ## cmin = -1; cmax = 1  # for anomaly plots
    ## cminm = -1.5; cmaxm = 1.5   # monthly
    
    cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
    cmap = 'blue2red_w20'

elif field == 'sicn':
    units = 'frac'
    conv=1
    cmin=-.15; cmax=.15
    cminp=-.15; cmaxp=.15
    cminm=-.15; cmaxm=.15
    cminmp=-.15; cmaxmp=.15
    cmap = 'red2blue_w20'
elif field == 'sic':
    units='m'
    conv=1/913.
    cmin=-.5
    cmax=.5
    cminm=-.5
    cmaxm=.5
    cmap = 'red2blue_w20'
elif field == 'gt':
    units='K'
    conv=1
    cmin=-2; cmax=2
    cminm=-3; cmaxm=3
    cminp=-.5; cmaxp=.5 # for when pert is 'ctl'
    cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
    cmap = 'blue2red_w20'
elif field == 'pmsl':
    units = 'hPa' # pretty sure hpa @@double check
    conv = 1
    cmin = -1; cmax = 1  # for anomaly plots
    cminm=-2; cmaxm=2  # for monthly maps
    cminp=cmin; cmaxp=cmax # for when pert is 'ctl'
    cminmp=cminm; cmaxmp=cmaxm
    cmap = 'blue2red_20'
elif field == 'pcp':
#    pct=1; units = '%'
    units = 'mm/day' # original: kg m-2 s-1
    conv = 86400  # convert from kg m-2 s-1 to mm/day
    cmin = -.2; cmax = .2  # for anomaly plots
    cminp=-.15; cmaxp=.15
    cminm = -.4; cmaxm = .4
    #cmap = 'PuOr'
    cmap = 'brown2blue_16w'
    cminpct=-12; cmaxpct=12
    cminmpct=-20; cmaxmpct=20
    cminmp =-.25; cmaxmp=.25
    cminpctp=-8; cmaxpctp=8
    cminpctmp=-12; cmaxpctmp=12
elif field == 'hfl': # sfc upward LH flux
    units = 'W/m2'
    conv = 1
    cmin = -5
    cmax = 5
    cminm = -8
    cmaxm = 8

elif field == 'hfs': # sfc upward SH flux
    units = 'W/m2'
    conv = 1
    cmin = -5
    cmax = 5
    cminm = -8
    cmaxm = 8

elif field == 'turb': # combine hfl and hfs
    units = 'W/m2'
    conv=1
    cmin=-10
    cmax=10
    cminm = -20
    cmaxm = 20
    cmap='blue2red_20'
    
elif field == 'net': # net of all sfc fluxes
    print " 'net ' not yet implemented! @@"
    
elif field == 'flg': # net downward LW at the sfc. Positive down?
    units = 'W/m2'
    conv = 1
    cmin = -5
    cmax = 5
    cminm = -8
    cmaxm = 8

elif field == 'fsg': # net (absorbed) solar downard at sfc
    units = 'W/m2'
    conv = 1
    cmin = -5
    cmax = 5
    cminm = -8
    cmaxm = 8

elif field == 'fn': # snow fraction
    units = '%'
    conv=100
    cmin = -5
    cmax = 5
    cminm = -5
    cmaxm = 5
    cmap = 'red2blue_w20'

elif field == 'pcpn': # snowfall rate (water equivalent, kg/m2/s)
    #pct = 1; units='%'
    units = 'mm/day'
    conv = 86400  # convert from kg m-2 s-1 to mm/day (I think it's same as pcp) @@
    cmap = 'brown2blue_16w'
    cmin = -.1 # for anomaly plots
    cmax = .1  # for anomaly plots
    cminm = -.15
    cmaxm = .15
    cminpct=-12
    cmaxpct=12
    cminmpct=-25
    cmaxmpct=25

elif field == 'zn': # snow depth (m)
#    pct=1; units='%'
    units = 'cm'
    conv = 100; # convert to cm
    cmap = 'brown2blue_16w'
    cmin = -2
    cmax = 2
    cminm = -3
    cmaxm = 3
    cminpct=-10
    cmaxpct=10
    cminmpct=-10
    cmaxmpct=10

elif field == 'su':
    units = 'm/s'
    conv = 1;
    cmap = 'blue2red_20'
    cmin = -1; cmax = 1
    cminm = -1; cmaxm = 1
    cminp = -.5; cmaxp=.5
    cminmp = -.5; cmaxmp=.5
elif field == 'sv':
    units = 'm/s'
    conv = 1;
    cmap = 'blue2red_20'
    cmin = -.5
    cmax = .5
    cminm = -.5
    cmaxm = .5
    
else:
    print 'No settings for ' + field


# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'

if field=='turb':
    field='hfl'; fieldb='hfs'
    fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
    fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_ts.nc'
    fnamecb = basepath + casename + subdir + casename + '_' + fieldb + '_' + timstr + '_ts.nc'
    fnamepb = basepath + casenamep + subdir + casenamep + '_' + fieldb + '_' + timstrp + '_ts.nc'

    fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv + \
           cnc.getNCvar(fnamecb,fieldb.upper(),timesel=timesel)*conv
    fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv+ \
           cnc.getNCvar(fnamepb,fieldb.upper(),timesel=timesel)*conv
    ## fldp2 = cnc.getNCvar(fnamep2,field.upper(),timesel='0002-01-01,0121-12-31')*conv+cnc.getNCvar(fnamep2b,fieldb.upper(),timesel='0002-01-01,0121-12-31')*conv
    ## fldp3 = cnc.getNCvar(fnamep3,field.upper(),timesel='0002-01-01,0121-12-31')*conv+cnc.getNCvar(fnamep3b,fieldb.upper(),timesel='0002-01-01,0121-12-31')*conv
    field='turb'
else:
    fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
    fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_ts.nc'
    ## fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp + '_ts.nc'
    ## fnamep3 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp + '_ts.nc'

    fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv
    fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv
    ## fldp2 = cnc.getNCvar(fnamep2,field.upper(),timesel='0002-01-01,0121-12-31')*conv
    ## fldp3 = cnc.getNCvar(fnamep3,field.upper(),timesel='0002-01-01,0121-12-31')*conv

# Get the data
## ncfilec = Dataset(fnamec,'r') # control
## ncfilep1 = Dataset(fnamep1,'r') # pert1
## ncfilep2 = Dataset(fnamep2,'r') # pert2
## ncfilep3 = Dataset(fnamep3,'r') # pert3
#   time period 2
## ncfilec2 = Dataset(fnamec2,'r') # control
## ncfilep12 = Dataset(fnamep12,'r') # pert1
## ncfilep22 = Dataset(fnamep22,'r') # pert2
## ncfilep32 = Dataset(fnamep32,'r') # pert3

#lat = ncfilec.variables['lat'][:]
#lon = ncfilec.variables['lon'][:]

lat = cnc.getNCvar(fnamec,'lat')
lon = cnc.getNCvar(fnamec,'lon')



## fldc = ncfilec.variables[field.upper()][(styr-1)*12:,:,:]*conv # time start year to end
## #fldc = np.append(fldc, ncfilec2.variables[field.upper()][...]*conv, axis=0)

## fldp1 = ncfilep1.variables[field.upper()][(styrp-1)*12:,:,:]*conv # time start year to end
## #fldp1 = np.append(fldp1,ncfilep12.variables[field.upper()][...]*conv, axis=0)

## fldp2 = ncfilep2.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
## #fldp2 = np.append(fldp2,ncfilep22.variables[field.upper()][...]*conv, axis=0)

## fldp3 = ncfilep3.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
## #fldp3 = np.append(fldp3,ncfilep32.variables[field.upper()][...]*conv, axis=0)


# annual time-series (3d)
anntsc = cutl.annualize_monthlyts(fldc)
anntsp = cutl.annualize_monthlyts(fldp)
#anntsp2 = cutl.annualize_monthlyts(fldp2)
#anntsp3 = cutl.annualize_monthlyts(fldp3)

nt,nlev,nlat = anntsc.shape # @@the var names are "wrong" but work fine in the script as written

##### Set which pert run! ######
## if casenamep == casenamep1:
##     fnamep = fnamep1
##     fnamepb = fnamep1b
##     anntsp = anntsp1
##     fldp = fldp1
## elif casenamep == casenamep2:
##     fnamep = fnamep2
##     fnamepb = fnamep2b
##     anntsp = anntsp2
##     fldp = fldp2
## elif casenamep == casenamep3:
##     fnamep = fnamep3
##     fnamepb = fnamep3b
##     anntsp = anntsp3
##     fldp = fldp3

if casename != 'kemctl1' and casename != 'kemhadctl':
    cmin=cminp; cmax=cmaxp
    cminm=cminmp; cmaxm=cmaxmp
    cminpct=cminpctp; cmaxpct=cmaxpctp
    cminpctm=cminpctmp; cmaxpctm=cmaxpctmp


timeavg = 'ANN'

if sigtype=='cont':
    suff='pdf'
else:
    suff='png'
    
if skipnmin==0:
    nminANN = cutl.calc_Nmin(anntsp,anntsc)

tstat,pval = sp.stats.ttest_ind(anntsp,anntsc,axis=0)
#tstatb,pvalb = sp.stats.ttest_ind(anntsp,anntsc,axis=0,equal_var=False) # basically the same as above
# Note that NaN is returned for zero variance (I think..from googling..)
# If that is the case, pcolormesh() needs a masked_array rather than ndarray (??)
#  : http://stackoverflow.com/questions/7778343/pcolormesh-with-missing-values

if plotann:
    title = timeavg + " " + field + ": " + casenamep + "-" + casename
    if pct:
        anntsctm = np.mean(anntsc,0)
        anntsctm = ma.masked_where(anntsctm<=0.01,anntsctm)

        plotfld = np.mean(anntsp-anntsc,0) / anntsctm * 100
        cmin=cminpct
        cmax=cmaxpct
    else:
        plotfld = np.mean(anntsp,0)-np.mean(anntsc,0)

    fig1 = plt.figure()
    bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,type='nh',\
                        title=title,units=units)
    cplt.addtsigm(bm,pval,lat,lon,type=sigtype) # add significance info (hatching. for contour, type='contour')
    if printtofile:            
        if pct:
            fig1.savefig(field + 'pctdiffsig' + sigtype + '_' + casenamep +\
                         '_v_' + casename + '_' + timeavg + '_nh.' + suff )
        else:
            fig1.savefig(field + 'diffsig' + sigtype + '_' + casenamep +\
                         '_v_' + casename + '_' + timeavg + '_nh.' + suff)



#   http://easycalculation.com/statistics/critical-t-test.php
# T critical (2-tailed 0.05): 110 dof: 1.9818, 60 dof: 2.0003, (109 dof: 1.982, 59 dof: 2.001 for DJF, NDJ)
# Or is dof: N+M-2: 218 (upper lim 200) dof: 1.9719, 118 dof: 1.9803 (216 dof:  116 dof: 1.9806 for DJF, NDJ)
#
#  Nmin = 2tc^2*(sp/(xbar-ybar))^2
#  sp = sqrt(  sum1_to_n[ (xi-xbar)^2 ] + sum1_to_m[ (yi-ybar)^2 ] / (n+m-2)  )
## tcritbig = 1.9719
## tcritsmall = 1.9803
## tcritsmallw = 1.9806
## sp = cutl.pooledstd(fldc,fldp)

nmins = np.zeros((12,fldc.shape[1],fldc.shape[2]))
sigs = np.ones((12,fldc.shape[1],fldc.shape[2]))


months=con.get_mon()

if plotallmos:
    title = field + ": " + casenamep + "-" + casename
    midx=0
    fig, spax = plt.subplots(2,6)
    #fig.set_size_inches(12,6)
    fig.set_size_inches(12,4.5)
    fig.subplots_adjust(hspace=0,wspace=0)

    for ax in spax.flat:

        monfldc = fldc[midx::12,:,:]
        monfldp = fldp[midx::12,:,:]

        tstat,pval = sp.stats.ttest_ind(monfldp,monfldc,axis=0)
        sigs[midx,:,:] = ma.masked_where(pval>0.05,pval) 
        if skipnmin==0:
            nmins[midx,:,:]  = cutl.calc_Nmin(monfldp,monfldc)

        if pct:
            monfldctm = np.mean(monfldc,0)
            monfldctm = ma.masked_where(monfldctm<=0.01,monfldctm)
            plotfld = np.mean(monfldp-monfldc,0) / monfldctm *100
            cminm=cminmpct
            cmaxm=cmaxmpct
        else:
            plotfld = np.mean(monfldp,0)-np.mean(monfldc,0)

        bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',\
                     title=months[midx],axis=ax,suppcb=1)
        ax.set_title(months[midx])
        cplt.addtsigm(bm,pval,lat,lon,type=sigtype)

        midx = midx+1

    cbar_ax = fig.add_axes([.91,.25, .02,.5])
    fig.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    plt.suptitle(title)
    if printtofile:
        if pct:
            fig.savefig(field + 'pctdiffsig' + sigtype + '_' + casenamep +\
                        '_v_' + casename + '_allmos_nh.' + suff)
        else:
            fig.savefig(field + 'diffsig' + sigtype + '_' + casenamep +\
                    '_v_' + casename + '_allmos_nh.' + suff)


    if skipnmin != 1:

        midx=0
        fig, spax = plt.subplots(2,6)
        fig.set_size_inches(12,6)
        fig.subplots_adjust(hspace=0,wspace=0)

        for ax in spax.flat:

            monfldc = fldc[midx::12,:,:]
            monfldp = fldp[midx::12,:,:]

            #tstat,pval = sp.stats.ttest_ind(monfldp,monfldc,axis=0)
            #nmins[midx,:,:]  = cutl.calc_Nmin(monfldp,monfldc)

            plotfld = nmins[midx,:,:]
            #plotfld = plotfld[sigs[midx,:,:]<=.05 # @@@ can't figure out how to mask out insignificant areas! ?

            bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=-100,cmax=100,cmap='gist_heat',type='nh',\
                         title=months[midx],axis=ax,suppcb=1)
            ax.set_title(months[midx])
            cplt.addtsigm(bm,pval,lat,lon,type=sigtype)

            midx = midx+1

        cbar_ax = fig.add_axes([.91,.25, .02,.5])
        fig.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
        plt.suptitle('NMIN: ' + title)
        if printtofile:
            fig.savefig(field + 'diff_NMIN_' + casenamep +\
                    '_v_' + casename + '_allmos_nh.' + suff)

# done with if plotallmos

if bimos:
    print 'plotbimos not implemented!'

if seasonal:


    cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division

    
    tstat = np.zeros((len(seasons),nlev,nlat))
    pval = np.zeros((len(seasons),nlev,nlat))
    fldcallseas = np.zeros((len(seasons),nlev,nlat))
    fldpallseas = np.zeros((len(seasons),nlev,nlat))
    
    incr = (cmaxm-cminm) / (cmlen)
    conts = np.arange(cminm,cmaxm+incr,incr)

    midx=0
    fig6,ax6 = plt.subplots(1,4) 
    fig6.set_size_inches(12,3)
    fig6.subplots_adjust(hspace=.15,wspace=.05)

    for ax in ax6.flat:

        
        if field=='turb':
            field='hfl'; fieldb='hfs'
            fldcsea = cnc.getNCvar(fnamec,field.upper(),timesel='0002-01-01,0121-12-31',
                                           seas=seasons[midx])*conv + cnc.getNCvar(fnamecb,fieldb.upper(),
                                           timesel='0002-01-01,0121-12-31',seas=seasons[midx])*conv
            fldpsea = cnc.getNCvar(fnamep,field.upper(),timesel='0002-01-01,0121-12-31',
                                           seas=seasons[midx])*conv + cnc.getNCvar(fnamepb,fieldb.upper(),
                                           timesel='0002-01-01,0121-12-31',seas=seasons[midx])*conv 
            field='turb'
        else:
            fldcsea = cnc.getNCvar(fnamec,field.upper(),timesel='0002-01-01,0121-12-31',
                                           seas=seasons[midx])*conv # @@ returns only 109 indices?, oh prob for winter...
            fldpsea = cnc.getNCvar(fnamep,field.upper(),timesel='0002-01-01,0121-12-31',
                                           seas=seasons[midx])*conv

        tstat[midx,:,:],pval[midx,:,:] = sp.stats.ttest_ind(fldpsea,fldcsea,axis=0)
        fldcallseas[midx,:,:] = np.mean(fldcsea,axis=0)
        fldpallseas[midx,:,:] = np.mean(fldpsea,axis=0)

        if pct:
            plotfld = (fldpallseas[midx,:,:]-fldcallseas[midx,:,:]) / fldcallseas[midx,:,:] *100
            cminm=cminmpct
            cmaxm=cmaxmpct
        else:
            plotfld = fldpallseas[midx,:,:] - fldcallseas[midx,:,:]

        bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',\
                        axis=ax,suppcb=1)
        ax.set_title(seasons[midx])
        cplt.addtsigm(bm,pval[midx,:,:],lat,lon,type=sigtype)

        midx = midx+1

    cbar_ax = fig6.add_axes([.91,.25, .02,.5])
    fig6.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    #plt.suptitle(title)
    if printtofile:
        if pct: # version 2 has seasons in new order
            fig6.savefig(field + 'pctdiffsig' + sigtype + '_' + casenamep +\
                        '_v_' + casename + '_seas_nh2.' + suff)
        else:
            fig6.savefig(field + 'diffsig' + sigtype + '_' + casenamep +\
                    '_v_' + casename + '_seas_nh2.' + suff)

