# new figs for WAIS paper, with stat sig
#  5/9/2015

import cccmautils as cutl
import numpy.ma as ma

plt.close('all')

printtofile=False
atmos=True
ocean=False

basepath = '/Users/kelly/School/DATA/'
basepath = '/Volumes/MyPassport2TB/DATA/ccsm4/'
casenamec = 'b40.20th.track1.1deg.006'; timeperc = '1970-1999'
casenamep = 'geo2035ensavg'; timeperp = '2045-2054'
casenamep2 = 'rcp8_5GHGrem1850'; timeperp = '2045-2054' # or 'rcp8_5GHGrem1850'
casenamep3 = 'b40.rcp8_5.1deg.006'; timeperp = '2045-2054' # or 'rcp8_5GHGrem1850'

subdir=''


if atmos:

    siglev=0.1
    var='T'; cmin=-3; cmax=3; units='$^\circ$C'
    fs=14 # fonstize

    # ############# ATMOS #####################
    # VERTICAL PLOT ===========================
    # like: U.b40.20th.track1.1deg.006.cam2.h0.1970-1999_timeseries.nc
    # and:  U.geo2035ensavg.cam2.h0.2045-2054_timeseries.nc

    fnamec =  basepath + casenamec + '/atm_proc/' + var + '/' + var + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc'
    fnamep =  basepath + casenamep + '/atm_proc/' + var + '/' + var + '.' + casenamep + '.cam2.h0.' + timeperp + '_timeseries.nc'

    lev = cnc.getNCvar(fnamec,'lev')
    lat = cnc.getNCvar(fnamep,'lat')

    fldc = cnc.getNCvar(fnamec,var).mean(axis=3) # for ttest
    fldp = cnc.getNCvar(fnamep,var).mean(axis=3) # for ttest
    plotdiff = fldp.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals=cutl.ttest_ind(fldp,fldc)

    var='U'; units='hPa' # ===========

    fnamec =  basepath + casenamec + '/atm_proc/' + var + '/' + var + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc'
    fnamep =  basepath + casenamep + '/atm_proc/' + var + '/' + var + '.' + casenamep + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldc2 = cnc.getNCvar(fnamec,var).mean(axis=3) # for ttest
    fldp2 = cnc.getNCvar(fnamep,var).mean(axis=3) # for ttest
    plotdiff2 = fldp2.mean(axis=0)-fldc2.mean(axis=0) # for the plot

    tstats,pvals2=cutl.ttest_ind(fldp2,fldc2)

    fig,axs=plt.subplots(2,1)
    fig.set_size_inches(6,10.5)
    ax = axs[0]
    cplt.vert_plot(plotdiff,lev,lat,cmin=cmin,cmax=cmax, title='$\Delta$ T ($^\circ$C)', 
                   units=units,cmap='blue2red_20',levlim=10,hPa=True,axis=ax,suppcb=True)
    cplt.addtsig(ax, pvals,lat,lev,type='hatch',siglevel=siglev)
    ax.annotate('a)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 

    ax = axs[1]
    cf = cplt.vert_plot(plotdiff2,lev,lat,cmin=cmin,cmax=cmax, title='$\Delta$ U (m/s)', 
                   units=units,cmap='blue2red_20',levlim=10,hPa=True,axis=ax,suppcb=True)
    cplt.addtsig(ax, pvals2,lat,lev,type='hatch',siglevel=siglev)
    ax.annotate('b)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 

    #cbar_ax = fig.add_axes([.91,.25, .02,.5])
    cbar_ax = fig.add_axes([.15,.02, .7,.02])   
    cb=fig.colorbar(cf,cax=cbar_ax,orientation='horizontal') # or do bm.colorbar....
    cb.set_ticks(np.arange(-3,4,1))
    cb.set_ticklabels(np.arange(-3,4,1))

    if printtofile:
        fig.savefig('WAISpaper_TUzm_sig' + str((1-siglev)*100) + '.png',dpi=400)


    # SH MAPS ==============================
    printtofile=False
    var='TREFHT'; units='$^\circ$C'; cmin=-4; cmax=4 # ===========

    pparams= {'cmin':cmin,'cmax':cmax,'type':'sh','suppcb':True}

    fnamec =  basepath + casenamec + '/atm_proc/2D/' + var + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc'
    fnamep3 =  basepath + casenamep3 + '/atm_proc/2D/' + var + '.' + casenamep3 + '.cam2.h0.' + timeperp + '_timeseries.nc'

    lon = cnc.getNCvar(fnamec,'lon')
    lons, lats = np.meshgrid(lon,lat)

    fldc = cnc.getNCvar(fnamec,var) # for ttest
    fldp3 = cnc.getNCvar(fnamep3,var) # for ttest
    plotdiff3 = fldp3.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals3=cutl.ttest_ind(fldp3,fldc)

    fig,axs=plt.subplots(2,3)
    fig.set_size_inches(9,5)
    ax=axs[0][0]
    bm,pc = cplt.kemmap(plotdiff3,lat,lon,title='RCP8.5',axis=ax,**pparams)
    cplt.addtsigm(bm,pvals3,lat,lon,type='hatch',siglevel=siglev)
    ax.annotate('a)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 

    cbar_ax = fig.add_axes([.35,.55, .02, .35])   
    cb=fig.colorbar(pc,cax=cbar_ax) 
    cb.set_ticks(np.arange(-4,5,2))
    cb.set_ticklabels(np.arange(-4,5,2))

    cmin=-2; cmax=2 
    fnamep =  basepath + casenamep + '/atm_proc/2D/' + var + '.' + casenamep + '.cam2.h0.' + timeperp + '_timeseries.nc'

    pparams= {'cmin':cmin,'cmax':cmax,'type':'sh','suppcb':True}

    fldp = cnc.getNCvar(fnamep,var) # for ttest
    plotdiff = fldp.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals=cutl.ttest_ind(fldp,fldc)

    ax=axs[0][1]
    bm,pc = cplt.kemmap(plotdiff,lat,lon,axis=ax,title='Sulf',**pparams)
    cplt.addtsigm(bm,pvals,lat,lon,type='hatch',siglevel=siglev)
    ax.annotate('b)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 

    fnamep2 =  basepath + casenamep2 + '/atm_proc/2D/' + var + '.' + casenamep2 + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldp2 = cnc.getNCvar(fnamep2,var) # for ttest
    plotdiff2 = fldp2.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals2=cutl.ttest_ind(fldp2,fldc)

    ax=axs[0][2]
    bm,pc = cplt.kemmap(plotdiff2,lat,lon,axis=ax,title='GHGrem',**pparams)
    cplt.addtsigm(bm,pvals2,lat,lon,type='hatch',siglevel=siglev)
    ax.annotate('c)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 

    cbar_ax = fig.add_axes([.9,.55, .02, .35])   
    cb = fig.colorbar(pc,cax=cbar_ax) #,orientation='horizontal')
    cb.set_ticks(np.arange(-2,3,1))

    var='TAUX'; units='N/m$^2$'; cmin=-0.03; cmax=0.03 # ===========
    varb='ICEFRAC'; unitsb='frac'; lw=1; ls='-'; ithresh=0.01; pcol='0.5' # linewidths, linestyles

    pparams= {'cmin':cmin,'cmax':cmax,'type':'sh','suppcb':True}

    fnamec =  basepath + casenamec + '/atm_proc/2D/' + var + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc'
    fnamep3 =  basepath + casenamep3 + '/atm_proc/2D/' + var + '.' + casenamep3 + '.cam2.h0.' + timeperp + '_timeseries.nc'
    fnamecb =  basepath + casenamec + '/atm_proc/2D/' + varb + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc'
    fnamep3b =  basepath + casenamep3 + '/atm_proc/2D/' + varb + '.' + casenamep3 + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldc = cnc.getNCvar(fnamec,var) # for ttest
    fldp3 = cnc.getNCvar(fnamep3,var) # for ttest
    plotdiff3 = fldp3.mean(axis=0)-fldc.mean(axis=0) # for the plot
    tstats,pvals3=cutl.ttest_ind(fldp3,fldc)

    # ice frac
    fldcb = cnc.getNCvar(fnamecb,varb).mean(axis=0)
    fldcb = ma.masked_where(fldcb<ithresh,fldcb)
    fldp3b = cnc.getNCvar(fnamep3b,varb).mean(axis=0)
    fldp3b = ma.masked_where(fldp3b<ithresh,fldp3b)

    ax=axs[1][0]
    bm,pc = cplt.kemmap(plotdiff3*-1,lat,lon,axis=ax,**pparams)
    cplt.addtsigm(bm,pvals3,lat,lon,type='hatch',siglevel=siglev)
    ax.annotate('d)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 
    bm.contour(lons,lats,fldcb,levels=[0.15, 0.15],colors='k',linewidths=lw,latlon=True)
    bm.contour(lons,lats,fldp3b,levels=[0.15, 0.15],colors=pcol,linewidths=lw,latlon=True,linestyles=ls)

    fnamep =  basepath + casenamep + '/atm_proc/2D/' + var + '.' + casenamep + '.cam2.h0.' + timeperp + '_timeseries.nc'
    fnamepb =  basepath + casenamep + '/atm_proc/2D/' + varb + '.' + casenamep + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldp = cnc.getNCvar(fnamep,var) # for ttest
    plotdiff = fldp.mean(axis=0)-fldc.mean(axis=0) # for the plot
    tstats,pvals=cutl.ttest_ind(fldp,fldc)

    # ice frac
    fldpb = cnc.getNCvar(fnamepb,varb).mean(axis=0)
    fldpb = ma.masked_where(fldpb<ithresh,fldpb)

    ax=axs[1][1]
    bm,pc = cplt.kemmap(plotdiff*-1,lat,lon,axis=ax,**pparams)
    cplt.addtsigm(bm,pvals,lat,lon,type='hatch',siglevel=siglev)
    ax.annotate('e)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 
    bm.contour(lons,lats,fldcb,levels=[0.15, 0.15],colors='k',linewidths=lw,latlon=True)
    bm.contour(lons,lats,fldpb,levels=[0.15, 0.15],colors=pcol,linewidths=lw,latlon=True,linestyles=ls)

    fnamep2 =  basepath + casenamep2 + '/atm_proc/2D/' + var + '.' + casenamep2 + '.cam2.h0.' + timeperp + '_timeseries.nc'
    fnamep2b =  basepath + casenamep2 + '/atm_proc/2D/' + varb + '.' + casenamep2 + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldp2 = cnc.getNCvar(fnamep2,var) # for ttest
    plotdiff2 = fldp2.mean(axis=0)-fldc.mean(axis=0) # for the plot
    tstats,pvals2=cutl.ttest_ind(fldp2,fldc)

    # ice frac
    fldp2b = cnc.getNCvar(fnamep2b,varb).mean(axis=0)
    fldp2b = ma.masked_where(fldp2b<ithresh,fldp2b)

    ax=axs[1][2]
    bm,pc = cplt.kemmap(plotdiff2*-1,lat,lon,axis=ax,**pparams)
    cplt.addtsigm(bm,pvals2,lat,lon,type='hatch',siglevel=siglev)
    ax.annotate('f)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 
    bm.contour(lons,lats,fldcb,levels=[0.15, 0.15],colors='k',linewidths=lw,latlon=True)
    bm.contour(lons,lats,fldpb,levels=[0.15, 0.15],colors=pcol,linewidths=lw,latlon=True,linestyles=ls)

    cbar_ax = fig.add_axes([.9,.1, .02, .35])   
    cb = fig.colorbar(pc,cax=cbar_ax) #,orientation='horizontal')
    cb.set_ticks(np.arange(-0.03, 0.039, .01))

    if printtofile:
        fig.savefig('WAISpaper_SHmaps_sig' + str((1-siglev)*100) + '.png',dpi=400)


# annual average with nc tools:
# # Annual average (use the feature of 'Duration')
#      ncra -O --mro -d time,"1956-01-01 00:00:0.0","2005-12-31 23:59:9.9",12,12
# my version failed b/c nco version too old: no mro option
# ncra -O --mro -d time "1970-01-01 00:00:0.0","1979-12-31 23:59:9.9",12,12 b40.20th.track1.1deg.006.pop.h.WISOP.197001-197912.nc b40.20th.track1.1deg.006.pop.h.WISOP.1970-1979_nctimeseries.nc
# not using cdo because it gave me 11 outputs! there should only be 10....


# ############# OCEAN #####################
if ocean:
    filenamec = basepath + casenamec + '/' + casenamec + '.pop.ANN.1970-1999.nc'
    filenamep = basepath + casenamep + '/' + casenamep + '.pop.ANN.2045-2054.nc'

    print filenamec
    print filenamep

