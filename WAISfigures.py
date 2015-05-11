# new figs for WAIS paper, with stat sig
#  5/9/2015

import cccmautils as cutl

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

    var='T'; cmin=-3; cmax=3; units='$^\circ$C'

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
    cplt.addtsig(ax, pvals,lat,lev,type='hatch')

    ax = axs[1]
    cf = cplt.vert_plot(plotdiff2,lev,lat,cmin=cmin,cmax=cmax, title='$\Delta$ U (m/s)', 
                   units=units,cmap='blue2red_20',levlim=10,hPa=True,axis=ax,suppcb=True)
    cplt.addtsig(ax, pvals2,lat,lev,type='hatch')

    #cbar_ax = fig.add_axes([.91,.25, .02,.5])
    cbar_ax = fig.add_axes([.15,.02, .7,.02])
    
    fig.colorbar(cf,cax=cbar_ax,orientation='horizontal') # or do bm.colorbar....
    cbar_ax.set_xticks(np.arange(-3,4,1))
    cbar_ax.set_xticklabels(np.arange(-3,4,1))
    print 'ADD Panel Letters @@@ !!'
    if printtofile:
        fig.savefig('WAISpaper_TUzm.png',dpi=400)


    # SH MAPS ==============================
    printtofile=True
    var='TREFHT'; units='$^\circ$C'; cmin=-4; cmax=4 # ===========

    fnamec =  basepath + casenamec + '/atm_proc/2D/' + var + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc'
    fnamep3 =  basepath + casenamep3 + '/atm_proc/2D/' + var + '.' + casenamep3 + '.cam2.h0.' + timeperp + '_timeseries.nc'

    lon = cnc.getNCvar(fnamec,'lon')

    fldc = cnc.getNCvar(fnamec,var) # for ttest
    fldp3 = cnc.getNCvar(fnamep3,var) # for ttest
    plotdiff3 = fldp3.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals3=cutl.ttest_ind(fldp3,fldc)

    fig,axs=plt.subplots(2,3)
    fig.set_size_inches(10,5)
    ax=axs[0][0]
    bm,pc = cplt.kemmap(plotdiff3,lat,lon,cmin=cmin,cmax=cmax,axis=ax,type='sh')
    cplt.addtsigm(bm,pvals3,lat,lon,type='hatch')

    cmin=-2; cmax=2 
    fnamep =  basepath + casenamep + '/atm_proc/2D/' + var + '.' + casenamep + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldp = cnc.getNCvar(fnamep,var) # for ttest
    plotdiff = fldp.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals=cutl.ttest_ind(fldp,fldc)

    ax=axs[0][1]
    bm,pc = cplt.kemmap(plotdiff,lat,lon,cmin=cmin,cmax=cmax,axis=ax,type='sh')
    cplt.addtsigm(bm,pvals,lat,lon,type='hatch')

    fnamep2 =  basepath + casenamep2 + '/atm_proc/2D/' + var + '.' + casenamep2 + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldp2 = cnc.getNCvar(fnamep2,var) # for ttest
    plotdiff2 = fldp2.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals2=cutl.ttest_ind(fldp2,fldc)

    ax=axs[0][2]
    bm,pc = cplt.kemmap(plotdiff2,lat,lon,cmin=cmin,cmax=cmax,axis=ax,type='sh')
    cplt.addtsigm(bm,pvals2,lat,lon,type='hatch')



    var='TAUX'; units='N/m$^2$'; cmin=-0.03; cmax=0.03 # ===========

    fnamec =  basepath + casenamec + '/atm_proc/2D/' + var + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc'
    fnamep3 =  basepath + casenamep3 + '/atm_proc/2D/' + var + '.' + casenamep3 + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldc = cnc.getNCvar(fnamec,var) # for ttest
    fldp3 = cnc.getNCvar(fnamep3,var) # for ttest
    plotdiff3 = fldp3.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals3=cutl.ttest_ind(fldp3,fldc)

    ax=axs[1][0]
    bm,pc = cplt.kemmap(plotdiff3*-1,lat,lon,cmin=cmin,cmax=cmax,axis=ax,type='sh')
    cplt.addtsigm(bm,pvals3,lat,lon,type='hatch')

    fnamep =  basepath + casenamep + '/atm_proc/2D/' + var + '.' + casenamep + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldp = cnc.getNCvar(fnamep,var) # for ttest
    plotdiff = fldp.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals=cutl.ttest_ind(fldp,fldc)

    ax=axs[1][1]
    bm,pc = cplt.kemmap(plotdiff*-1,lat,lon,cmin=cmin,cmax=cmax,axis=ax,type='sh')
    cplt.addtsigm(bm,pvals,lat,lon,type='hatch')

    fnamep2 =  basepath + casenamep2 + '/atm_proc/2D/' + var + '.' + casenamep2 + '.cam2.h0.' + timeperp + '_timeseries.nc'

    fldp2 = cnc.getNCvar(fnamep2,var) # for ttest
    plotdiff2 = fldp2.mean(axis=0)-fldc.mean(axis=0) # for the plot

    tstats,pvals2=cutl.ttest_ind(fldp2,fldc)

    ax=axs[1][2]
    bm,pc = cplt.kemmap(plotdiff2*-1,lat,lon,cmin=cmin,cmax=cmax,axis=ax,type='sh')
    cplt.addtsigm(bm,pvals2,lat,lon,type='hatch')

    print 'ADD SIC!! and panel letters @@'

    if printtofile:
        fig.savefig('WAISpaper_SHmaps.png',dpi=400)


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

