# new figs for WAIS paper, with stat sig
#  5/9/2015

import cccmautils as cutl
import numpy.ma as ma
import cccmacmaps as ccm

plt.close('all')

printtofile=False
atmos=False
atmosts=True
ocean=False
oceanwvel=False
oceants=False

basepath2 = '/Users/kelly/School/DATA/'
basepath = '/Volumes/MyPassport2TB/DATA/ccsm4/'
casenamec = 'b40.20th.track1.1deg.006'; timeperc = '1970-1999'
casenamep = 'geo2035ensavg'; timeperp = '2045-2054'
casenamep2 = 'rcp8_5GHGrem1850'; timeperp = '2045-2054' 
casenamep3 = 'b40.rcp8_5.1deg.006'; timeperp = '2045-2054' 

subdir=''

siglev=0.1


if atmos:

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

    # ================= ICE FRACTION MAPS ======================

    printtofile=False
    cmin=-10; cmax=10
    pparams= {'cmin':cmin,'cmax':cmax,'type':'sh','suppcb':False, 'cmap': 'red2blue_w20'}

    fldcb = cnc.getNCvar(fnamecb,varb)*100 # for ttest
    fldcb=ma.masked_where(fldcb<=0,fldcb)
    fldpb = cnc.getNCvar(fnamepb,varb)*100 # for ttest
    fldpb=ma.masked_where(fldpb<=0,fldpb)
    plotdiffb = fldpb.mean(axis=0)-fldcb.mean(axis=0) # for the plot
    tstats,pvalsb=cutl.ttest_ind(fldpb,fldcb)

    # ice frac
    #fldcb = ma.masked_where(fldcb<ithresh,fldcb)
    #fldp3b = ma.masked_where(fldp3b<ithresh,fldp3b)

    fig,axs=plt.subplots(1,2)
    fig.set_size_inches(8,4)
    ax=axs[0]
    bm,pc = cplt.kemmap(plotdiffb,lat,lon,title='Sulf',axis=ax,**pparams)
    cplt.addtsigm(bm,pvalsb,lat,lon,type='hatch',siglevel=siglev)
    ax.annotate('a)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 

    fldp2b = cnc.getNCvar(fnamep2b,varb)*100 # for ttest
    fldp2b=ma.masked_where(fldp2b<=0,fldp2b)
    plotdiff2b = fldp2b.mean(axis=0)-fldcb.mean(axis=0) # for the plot
    tstats,pvals2b=cutl.ttest_ind(fldp2b,fldcb)

    ax=axs[1]
    bm,pc = cplt.kemmap(plotdiff2b,lat,lon,title='GHGrem',axis=ax,**pparams)
    cplt.addtsigm(bm,pvals2b,lat,lon,type='hatch',siglevel=siglev)
    ax.annotate('b)',xy=(0,1.02),xycoords='axes fraction',fontsize=fs) 

    if printtofile:
        fig.savefig('WAISpaper_ICEFRACmaps_sig' + str((1-siglev)*100) + '.png',dpi=400)

    printtofile=False

    # ================= ZONAL WIND STRESS ======================

    var='TAUX' # wind stress

    fnamec =  basepath + casenamec + '/atm_proc/2D/' + var + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc'
    fnamep =  basepath + casenamep + '/atm_proc/2D/' + var + '.' + casenamep + '.cam2.h0.' + timeperp + '_timeseries.nc'
    fnamep2 =  basepath + casenamep2 + '/atm_proc/2D/' + var + '.' + casenamep2 + '.cam2.h0.' + timeperp + '_timeseries.nc'
    fnamep3 =  basepath + casenamep3 + '/atm_proc/2D/' + var + '.' + casenamep3 + '.cam2.h0.' + timeperp + '_timeseries.nc'


    # Now plot zonal mean wind and wind stress curl
    # On challenger: /home/disk/eos11/ccsmruns/compute_windstress (or something close)
    basepath='/Users/kelly/School/DATA/wscurl/'

    atmfile=basepath + '../b40.20th.track1.1deg.006_ANN_climo.nc' # for landfrac
    lfrac=np.squeeze(cnc.getNCvar(atmfile,'LANDFRAC'))
    ofrac=1-lfrac

    #fname=basepath + 'geo2035ensavg_v_b40.20th.track1.1deg.006_ANN_climo_wscurl.nc' # ??
    # use regular climos instead: geo2035ensavg_ANN_climo_wscurl.nc
    # field = 'wscurl'
    # @@ need to get landfrac

    fname = basepath + 'wscurl.b40.20th.track1.1deg.006.cam2.1970-1999.ANN.nc'
    fnameclim = basepath + 'b40.20th.track1.1deg.006_ANN_climo_wscurl.nc'

    # climo test
    fldcclim = np.squeeze(cnc.getNCvar(fnameclim,'wscurl'))*-1 
    fldcclimwt = fldcclim*ofrac
    fldcclimwtzm = np.mean(fldcclimwt,axis=1)

    fldc = np.squeeze(cnc.getNCvar(fname,'wscurl'))*-1 
    lat = cnc.getNCvar(fname,'lat')
    lon = cnc.getNCvar(fname,'lon')
    print fldc.shape

    (nt,nlt,nln)=fldc.shape
    ofract=np.tile(ofrac,(nt,1,1))
    fldcwt = fldc*ofract
    fldcwtzm = np.mean(fldcwt,axis=2)


    #--- wind stress files
    fnamews = basepath + '../b40.20th.track1.1deg.006_ANN_climo.nc' # 1970-1999
    wsfldc = np.squeeze(cnc.getNCvar(fnamec,'TAUX'))*-1 # (also mult by -1 ?) Yes, atmos field

    # tile ofroc with time
    (nt,nlt,nln)=wsfldc.shape
    ofract=np.tile(ofrac,(nt,1,1))
    wsfldcwt = wsfldc*ofract
    wsfldcwtzm = np.mean(wsfldcwt,axis=2)

    fnames = {'Sulf': fnamep, # 2045-2054
              'GHGrem': fnamep2, # 2045-2054
              'RCP8.5': fnamep3} # 2045-2054

    # wind stress curl files
    wscfnames = {'Sulf': basepath + 'wscurl.geo2035ensavg.cam2.2045-2054.ANN.nc',
              'GHGrem': basepath + 'wscurl.rcp8_5GHGrem1850.cam2.2045-2054.ANN.nc',
              'RCP8.5': basepath + 'wscurl.b40.rcp8_5.1deg.006.cam2.2045-2054.ANN.nc' }
    wscclimfnames = {'Sulf': basepath + 'geo2035ensavg_ANN_climo_wscurl.nc',
                     'GHGrem': basepath + 'rcp8_5GHGrem1850_ANN_climo_wscurl.nc',
                     'RCP8.5': basepath + 'b40.rcp8_5.1deg.006_ANN_climo_wscurl.nc'}


    coldt = {'RCP8.5': ccm.get_linecolor('firebrick1'),
             'Sulf': ccm.get_linecolor('mediumblue'),
             'GHGrem': ccm.get_linecolor('darkolivegreen3') }

    printtofile=True

    fig,axs = plt.subplots(1,2,sharex=True)
    fig.set_size_inches(14,3.2) # match zonal mean ocean TEMP
    ax=axs[0] # wind stress

    savepvalswsdt={}
    wssaveplotdt={}
    for fname in ('RCP8.5','Sulf','GHGrem'):
        fld = np.squeeze(cnc.getNCvar(fnames[fname],'TAUX'))*-1         

        #fldd = fld.mean(axis=0)-wsfldc.mean(axis=0)
        #flddwt = fldd*ofrac
        #flddwtzm = np.mean(flddwt,axis=1)
        
        # do stats now (Nothing is sig!)
        (nt,nlt,nln)=fld.shape
        ofract=np.tile(ofrac,(nt,1,1))
        fldwt = fld*ofract
        fldwtzm = np.mean(fldwt,axis=2)

        tstat,pvalsws = cutl.ttest_ind(fldwtzm,wsfldcwtzm)
        savepvalswsdt[fname] = pvalsws
        plotfld=fldwtzm.mean(axis=0)-wsfldcwtzm.mean(axis=0)
        wssaveplotdt[fname] = plotfld

        ax.plot(lat,plotfld,color=coldt[fname],linewidth=2)

    for fname in ('RCP8.5','Sulf','GHGrem'):
        pvalsws=savepvalswsdt[fname]
        plotfld = wssaveplotdt[fname]
        plotfldm = ma.masked_where(pvalsws>siglev, plotfld)
        ax.plot(lat,plotfldm, color=coldt[fname],linewidth=6,alpha=0.7)

    print '@@@ fix legend...'
    ax.set_xlim(-75,-40)
    ax.axhline(y=0,color='k')
    ax.set_xticks(np.arange(-75,-35,5))
    ax.set_xticklabels(['','70$^\circ$S','','60$^\circ$S','',\
                        '50$^\circ$S','','40$^\circ$S'],fontsize=18)
    ax.set_yticks(np.arange(-0.01,0.025,0.005))
    ax.set_yticklabels([-0.01,'',0,'', 0.01, '', 0.02],fontsize=18)
    ax.set_ylim(-0.01,0.02)
    ax.set_title(r'$\Delta$ TAUx (N/m$^2$)',fontsize=18)
    ax.legend(('RCP8.5','Sulf','GHGrem'),loc='upper right',fancybox=True,framealpha=0.5,fontsize=15)

    cmintest=-1e-7; cmaxtest=1e-7

    ax=axs[1] # wind stress curl
    savedt={}
    savepvaldt={}
    for fname in ('RCP8.5','Sulf','GHGrem'): #'RCP8.5','Sulf','GHGrem'
        print wscfnames[fname]

        if fname=='RCP8.5': conv=1 # already multiplied (CMIP5 output)
        else: conv=-1

        fld = np.squeeze(cnc.getNCvar(wscfnames[fname],'wscurl'))*conv

        # climo test: climo looks good.

        """fldclim = np.squeeze(cnc.getNCvar(wscclimfnames[fname],'wscurl'))*-1 # get climo for testing
        fldclimwt = fldclim*ofrac
        fldclimwtzm = np.mean(fldclimwt,axis=1)

        plotfldclimmap = fldclim - fldcclim
        plotfldclim = fldclimwtzm - fldcclimwtzm

        fig,axtest = plt.subplots(1,1)
        cplt.kemmap(plotfldclimmap,lat,lon,type='sh',axis=axtest,title=fname,cmin=cmintest,cmax=cmaxtest) # @@@
        
        plt.figure(); plt.plot(lat,plotfldclim); plt.title(fname); plt.ylim(-2e-8,2e-8);
        """

        #fldd = fld-fldc
        #flddwt = fldd*ofrac
        #flddwtzm = np.mean(flddwt,axis=1)

        print '@@ RCP85 data is wrong -- much too large, pattern seems shifted...'

        # do stats now (Nothing is sig!)
        (nt,nlt,nln)=fld.shape
        ofract=np.tile(ofrac,(nt,1,1))
        fldwt = fld*ofract
        fldwtzm = np.mean(fldwt,axis=2)

        tstat,pvalswsc = cutl.ttest_ind(fldwtzm,fldcwtzm)
        plotfldwsc=fldwtzm.mean(axis=0)-fldcwtzm.mean(axis=0)

        

        ax.plot(lat,plotfldwsc,color=coldt[fname],linewidth=2)
        plotfldwscm = ma.masked_where(pvalswsc>siglev, plotfldwsc)
        ax.plot(lat,plotfldwscm, color=coldt[fname],linewidth=6,alpha=0.7)

        savedt[fname] = plotfldwscm # @@@
        savepvaldt[fname] = pvalswsc

    ax.set_xlim(-75,-40)
    ax.axhline(y=0,color='k')
    ax.set_xticklabels(['','70$^\circ$S','','60$^\circ$S','',\
                        '50$^\circ$S','','40$^\circ$S'],fontsize=18)
    ax.set_yticks(np.arange(-3e-8,2.5e-8,.5e-8))
    ax.set_yticklabels([-3,'',-2,'',-1,'',0,'',1,'',2],fontsize=18)
    ax.set_title(r'$\Delta$ ( $\nabla$ x TAU ) (10$^{-8}$ N/m$^3$)',fontsize=18)

    if printtofile:
        fig.savefig('TAUX_curlTAUX_allruns_zonmean_ANN.pdf')

    printtofile=False




# annual average with nc tools:
# # Annual average (use the feature of 'Duration')
#      ncra -O --mro -d time,"1956-01-01 00:00:0.0","2005-12-31 23:59:9.9",12,12
# my version failed b/c nco version too old: no mro option
# ncra -O --mro -d time "1970-01-01 00:00:0.0","1979-12-31 23:59:9.9",12,12 b40.20th.track1.1deg.006.pop.h.WISOP.197001-197912.nc b40.20th.track1.1deg.006.pop.h.WISOP.1970-1979_nctimeseries.nc
# not using cdo because it gave me 11 outputs! there should only be 10....

# #### ======= Atmos SAT with time ======
if atmosts:

    var='TREFHT'
    timeperp='2035-2094'
    timeperc='1970-1999'
    fnamec =  basepath + casenamec + '/atm_proc/2D/' + var + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc' # ann average
    fnamep =  basepath + casenamep + '/atm_proc/2D/' + var + '.' + casenamep + '.cam2.' + timeperp + '.nc'

    lon = cnc.getNCvar(fnamec,'lon')
    lat = cnc.getNCvar(fnamep,'lat')

    fldc = cnc.getNCvar(fnamec,var).mean(axis=2) # for ttest
    fldp = cutl.annualize_monthlyts(cnc.getNCvar(fnamep,var))
    fldp = fldp.mean(axis=2) # for ttest
    plotdiff = fldp-fldc.mean(axis=0) # for the plot


    tt=5
    pvals=np.zeros((fldp.shape[0]-tt,len(lat)))
    for ii in np.arange(tt,fldp.shape[0]):
        # use the surrounding 10 years to estimate that grid point's significance
        test = fldpstats[ii-tt:ii+tt+1,:] #make the year the midpoint of 11-year period

        tstats,pvals[ii-tt,:]=cutl.ttest_ind(fldp[ii-tt:ii+tt+1,:],fldc)


    xx=np.arange(2035,2095)
    xxpvals=np.arange(2035+tt,2095)
    
    lats,times = np.meshgrid(lat,xx)

    cmax=2; cmin=-2
    cmap='blue2red_20'; cmlen=20.
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)


    fig,ax = plt.subplots(1,1)
    fig.set_size_inches(8,5)

    CF = ax.contourf(times,lats,plotdiff,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
    cplt.addtsig(ax,pvals.T,xxpvals,lat,siglevel=siglev,type='hatch')

    ax.set_xlabel('Years',fontsize=14)
    ax.set_ylabel('Latitude',fontsize=14)
    ax.set_yticks(np.arange(-90,100,30))

    cbar_ax = fig.add_axes([.91,.15, .02,.7])
    fig.colorbar(CF,cax=cbar_ax)

    if printtofile:
        fig.savefig('WAISpap_SuppFig_TREFHTtimeserieswithlat.png',dpi=300)






# ############# OCEAN #####################


def ocnzonalmean_ccsm4(fld,tarea,kmt):
    """ calculates an area-weighted zonal average with depth.
                fld must have the last 3 dims be lev x lat x lon (can handle time x lev x lat x lon)
                tarea has dims of lat x lon
                kmt has dims of just lev
    """

    ndims=fld.ndim
    if ndims==4:
        (ntime,nlev,nlat,nlon) = fld.shape

        # put depth/level first dim for calcs, so now nlev,ntime,nlat,nlon
        fldsave=fld
        fld=np.transpose(fld,(1,0,2,3))
        tareat = np.tile(tarea,(nlev,ntime,1,1))
        kmtt = np.tile(kmt,(ntime,1,1))
    elif ndims==3:
        (nlev,nlat,nlon)=fld.shape
        ntime=None
        tareat = np.tile(tarea,(nlev,1,1))
        kmtt=kmt
    else:
        print 'fld.shape is <3 or >4!'
        return -1


    for lii in np.arange(0,nlev):
        # mask out levels below sea floor
        #print fld.shape
        #print kmtt.shape
        fld[lii,...] = ma.masked_where(kmtt <= lii,fld[lii,...])
        #tareat[lii,...] = ma.masked_where(kmt <= lii,tareat[lii,...]) # @@ not working??

    tareat = ma.masked_where(fld.mask,tareat)

    totzonalarea= ma.sum(tareat,axis=ndims-1) # sum over last dimension (lon)

    if ndims==4:
        tiledims=(tareat.shape[ndims-1],1,1,1)
        trans=(1,2,3,0)
    else:
        tiledims=(tareat.shape[ndims-1],1,1)
        trans=(1,2,0)

    totzonalareat=np.tile(totzonalarea,tiledims) # tile total area over lon dim to make weights

    totzonalareat=np.transpose(totzonalareat,trans)
    fullzonalwgts= tareat/totzonalareat
    fullzonalwgts=ma.masked_where(tareat.mask,fullzonalwgts) # remember have to use fld.mask !

    fld=np.squeeze(ma.mean(fld,axis=ndims-1))
    # now put time back in first dim if necessary:
    if ndims==4:
        fld=np.transpose(fld,(1,0,2))

    return fld


#@@def ocnregmean_ccsm4(fld,tarea,kmt):



if ocean:


    # /Volumes/MyPassport2TB/DATA/ccsm4/b40.20th.track1.1deg.006/ocn_proc/TEMP/TEMP.b40.20th.track1.1deg.006.pop.h.1970-1999.nc
    var='TEMP'
    reg='SH' # SH or PIG
    lonlims = [230,280]; strlims='80W-120W' # for PIG

    basepath='/Volumes/MyPassport2TB/DATA/ccsm4/'
    basepath2='/Users/kelly/School/DATA/'

    casenamec='b40.20th.track1.1deg.006'
    casenamep='geo2035ensavg'
    casenamep2='rcp8_5GHGrem1850'
    casenamep3='b40.rcp8_5.1deg.006'

    timeperc='1970-1999'
    timeperp='2045-2054'
    regstr=''

    if reg=='PIG':
        reg='SH'
        pig=True
    else:
        pig=False

    filenameclim = basepath2 + casenamec + '/' + casenamec + '.pop.ANN.' + timeperc + '.nc' 

    filenamec = basepath + casenamec + '/ocn_proc/' + var + '/' + var + '.' + casenamec + '.pop.' + timeperc + '.' + reg + '.nc' 
    filenamep = basepath + casenamep + '/ocn_proc/' + var + '/' + var + '.' + casenamep  + '.pop.' + timeperp + '.' + reg + '.nc' 
    filenamep2 = basepath + casenamep2 + '/ocn_proc/' + var + '/' + var + '.' + casenamep2  + '.pop.' + timeperp + '.' + reg + '.nc' 
    filenamep3 = basepath + casenamep3 + '/ocn_proc/' + var + '/' + var + '.' + casenamep3  + '.pop.' + timeperp + '.' + reg + '.nc' 


    print filenamec
    print filenamep

    kmt = cnc.getNCvar(filenameclim,'KMT') # get from climo
    tarea=cnc.getNCvar(filenameclim,'TAREA') # get from climo
    tlatc=cnc.getNCvar(filenameclim,'TLAT')
    tlonc=cnc.getNCvar(filenameclim,'TLONG')

    zt = cnc.getNCvar(filenamec, 'z_t')
    tlat = cnc.getNCvar(filenamec,'TLAT') # these are the shape of the var
    tlon = cnc.getNCvar(filenamec,'TLONG') 

    fldc = cnc.getNCvar(filenamec,var)
    fldp = cnc.getNCvar(filenamep,var)
    fldp2 = cnc.getNCvar(filenamep2,var)

    (ntime,nlev,nlat,nlon) = fldc.shape

    if reg=='SH' and not pig:
        #tareareg = tarea[tlat<0].reshape(fldc[0,0,...].shape)
        tareareg = tarea[tlatc<0].reshape(fldc[0,0,...].shape)
        kmtreg = kmt[tlatc<0].reshape(fldc[0,0,...].shape)
        #kmtreg = kmt[tlat<0].reshape(fldc[0,0,...].shape)
    elif pig:
        # not sure cdo selection for PIG was working correctly. do it here instead:
        # first do the data
        tmp = fldc[:,:,np.logical_and(tlon>=lonlims[0],tlon<=lonlims[1])]
        (nt,nl,space)=tmp.shape # unknown number of lons
        nlon=space/nlat
        fldc=tmp.reshape((ntime,nlev,nlat,nlon))
        tmp = fldp[:,:,np.logical_and(tlon>=lonlims[0],tlon<=lonlims[1])]
        (nt,nl,space)=tmp.shape # unknown number of lons
        fldp=tmp.reshape((nt,nlev,nlat,nlon))
        tmp = fldp2[:,:,np.logical_and(tlon>=lonlims[0],tlon<=lonlims[1])]
        (nt,nl,space)=tmp.shape # unknown number of lons
        fldp2=tmp.reshape((nt,nlev,nlat,nlon))

        # now fix the coords
        tlat=tlat[np.logical_and(tlon>=lonlims[0],tlon<=lonlims[1])].reshape((nlat,nlon))
        tlon=tlon[np.logical_and(tlon>=lonlims[0],tlon<=lonlims[1])].reshape((nlat,nlon))

        # here we get lons between 230 and 280 at the same time as lats<0
        tareareg = tarea[np.logical_and(np.logical_and(tlonc>=lonlims[0],tlonc<=lonlims[1]), tlatc<0)].reshape(fldc[0,0,...].shape)
        kmtreg = kmt[np.logical_and(np.logical_and(tlonc>=lonlims[0],tlonc<=lonlims[1]), tlatc<0)].reshape(fldc[0,0,...].shape)


    if pig:
        reg='PIG'
        regstr='PIG '

    fldczm=ocnzonalmean_ccsm4(fldc,tareareg,kmtreg)
    fldpzm=ocnzonalmean_ccsm4(fldp,tareareg,kmtreg)
    fldp2zm=ocnzonalmean_ccsm4(fldp2,tareareg,kmtreg)

    tstats,pvals=cutl.ttest_ind(fldpzm,fldczm)
    tstats,pvals2=cutl.ttest_ind(fldp2zm,fldczm)
    
    # =========== paper figure ====================

    tlats,zlevs = np.meshgrid(np.squeeze(tlat[:,1]),zt/100.)

    plotfld=fldpzm.mean(axis=0)-fldczm.mean(axis=0)
    plotfld2=fldp2zm.mean(axis=0)-fldczm.mean(axis=0)

    
    ylim=800
    xlims=(-77,-50)
    ylims=(0,ylim)

    cmlen=float(20)

    cmap='blue2red_w20'


    fig,axs = plt.subplots(1,2,sharey=True)
    ax=axs[0]
    #fig.set_size_inches(10,3)
    fig.set_size_inches(14,3) # to match PIG region

    cmax=.5; cmin=-.5
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    CF = ax.contourf(tlats,zlevs,plotfld,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
    cplt.addtsig(ax,pvals,np.squeeze(tlat[:,1]),zt/100.,siglevel=siglev,type='cont')

    ax.set_ylim(ylims)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(0,900,100))
    ax.set_yticklabels([0,'',200,'',400,'',600,'',800],fontsize=18)
    ax.set_xlim(xlims)
    ax.set_xticks(np.arange(-75,-45,5))
    ax.set_xticklabels(['75$^\circ$S', '70$^\circ$S', '65$^\circ$S', \
                        '60$^\circ$S', '55$^\circ$S', '50$^\circ$S'],fontsize=18)
    ax.set_title(regstr + 'Sulf',fontsize=18)
    ax.set_ylabel('Depth (m)',fontsize=18)
    if pig:
        # add vert line
        ax.axvline(x=-65,linestyle='--',color='k') # @@@ the vert line is to show the area averaged in the VHT plots
    
    cbar_ax = fig.add_axes([.485,.15, .02,.7])
    fig.colorbar(CF,cax=cbar_ax)

    ax2=axs[1]
    cmax=.5; cmin=-.5
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    CF2 = ax2.contourf(tlats,zlevs,plotfld2,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
    cplt.addtsig(ax2,pvals2,np.squeeze(tlat[:,1]),zt/100.,siglevel=siglev,type='cont')

    ax2.set_ylim(ylims)
    ax2.invert_yaxis()
    ax2.set_yticks(np.arange(0,900,100))
    ax2.set_yticklabels([0,'',200,'',400,'',600,'',800],fontsize=18)
    ax2.set_xlim(xlims)
    ax2.set_xticks(np.arange(-75,-45,5))
    ax2.set_xticklabels(['75$^\circ$S', '70$^\circ$S', '65$^\circ$S', \
                        '60$^\circ$S', '55$^\circ$S', '50$^\circ$S'],fontsize=18)
    if pig:
        ax2.axvline(x=-65,linestyle='--',color='k') # @@@ the vert line is to show the area averaged in the VHT plots

    #ax2.set_title(casenamep2)
    ax2.set_title(regstr + 'GHGrem',fontsize=18)
    #cbar_ax = fig.add_axes([.91,.15, .02,.7])
    #fig.colorbar(CF,cax=cbar_ax)

    if printtofile:
        #fig.savefig('TEMPanom_subplot' + reg + '_ylim' + str(ylims[1]) + 'xlim' + str(xlims[1]) + '_c_sig' + str(1-siglev) + '.png',dpi=400)
        fig.savefig('TEMPanom_subplot' + reg + '_ylim' + str(ylims[1]) + 'xlim' + str(xlims[1]) + '_c_sig' + str(1-siglev) + 'cont.pdf')



################# WVEL averages ##################
if oceanwvel:

    climo=True # test with climo files first.


    var='WVEL'
    varb='WISOP'
    vart='TEMP'
    reg='SH'; pig=True # because the cdo files w/ PIG already selected seemed incorrect


    if climo:
        filenamec = filenameclim = filenamect = filenamecb = basepath2 + casenamec + '/' + casenamec + '.pop.ANN.' + timeperc + '.nc' 
        filenamep = filenamept = filenamepb = basepath2 + casenamep + '/' + casenamep + '.pop.ANN.' + timeperp + '.nc' 
        filenamep2 = filenamept2 = filenamepb2 = basepath2 + casenamep2 + '/' + casenamep2 + '.pop.ANN.' + timeperp + '.nc' 
        filenamep3 = filenamept3 = filenamepb3 = basepath2 + casenamep3 + '/' + casenamep3 + '.pop.ANN.' + timeperp + '.nc' 
    else:
        filenameclim = basepath2 + casenamec + '/' + casenamec + '.pop.ANN.' + timeperc + '.nc' 

        # TEMP
        filenamect = basepath + casenamec + '/ocn_proc/' + vart + '/' + vart + '.' + casenamec + '.pop.' + timeperc + '.' + reg + '.nc' 
        filenamept = basepath + casenamep + '/ocn_proc/' + vart + '/' + vart + '.' + casenamep  + '.pop.' + timeperp + '.' + reg + '.nc' 
        filenamept2 = basepath + casenamep2 + '/ocn_proc/' + vart + '/' + vart + '.' + casenamep2  + '.pop.' + timeperp + '.' + reg + '.nc' 

        # WVEL
        filenamec = basepath + casenamec + '/ocn_proc/' + var + '/' + var + '.' + casenamec + '.pop.' + timeperc + '.' + reg + '.nc' 
        filenamep = basepath + casenamep + '/ocn_proc/' + var + '/' + var + '.' + casenamep  + '.pop.' + timeperp + '.' + reg + '.nc' 
        filenamep2 = basepath + casenamep2 + '/ocn_proc/' + var + '/' + var + '.' + casenamep2  + '.pop.' + timeperp + '.' + reg + '.nc' 

        # WISOP
        filenamecb = basepath + casenamec + '/ocn_proc/' + varb + '/' + varb + '.' + casenamec + '.pop.' + timeperc + '.' + reg + '.nc' 
        filenamepb = basepath + casenamep + '/ocn_proc/' + varb + '/' + varb + '.' + casenamep  + '.pop.' + timeperp + '.' + reg + '.nc' 
        filenamepb2 = basepath + casenamep2 + '/ocn_proc/' + varb + '/' + varb + '.' + casenamep2  + '.pop.' + timeperp + '.' + reg + '.nc' 



    # get from climo: standard coords (standard shapes)
    kmt = cnc.getNCvar(filenameclim,'KMT') 
    tarea=cnc.getNCvar(filenameclim,'TAREA') 
    tlatc=cnc.getNCvar(filenameclim,'TLAT')
    tlonc=cnc.getNCvar(filenameclim,'TLONG')

    # get from var file: these are in the shape of the var
    zt = cnc.getNCvar(filenamec, 'z_t')
    zw = cnc.getNCvar(filenamec,'z_w')
    tlat = cnc.getNCvar(filenamec,'TLAT') 
    tlon = cnc.getNCvar(filenamec,'TLONG') 

    wvelc = np.squeeze(cnc.getNCvar(filenamec,var))/100. # convert to m/s
    wisopc = np.squeeze(cnc.getNCvar(filenamec,varb))/100.
    tempc = np.squeeze(cnc.getNCvar(filenamec,vart))

    wvelp = np.squeeze(cnc.getNCvar(filenamep,var))/100.
    wisopp = np.squeeze(cnc.getNCvar(filenamep,varb))/100.
    tempp = np.squeeze(cnc.getNCvar(filenamep,vart))

    if climo:
        (nlev,nlat,nlon) = tempc.shape # on T grid
        (nlevw,nlatw,nlonw) = wvelc.shape # on U grid
    else:
        (ntime,nlev,nlat,nlon) = tempc.shape # on T grid
        (ntimew,nlevw,nlatw,nlonw) = wvelc.shape # on U grid

    print tempc.shape # @@

    # Create the mask for the region: PIG
    lonlims = [230,280]; region = 'PIG'; strlims='80W-120W'
    rmaskout = np.logical_and(tlon>lonlims[0], tlon<lonlims[1]) # region mask! this masks OUT the region itself
    rmask = np.logical_or(tlon<=lonlims[0],tlon>=lonlims[1]) # use this one for averaging. keep only the region

    # tile the mask
    rmask = np.tile(rmask,(len(zt),1,1))
    print rmask.shape # @@

    # tile area with depth
    tareat = np.tile(tarea,(len(zt),1,1))




    # TAKE INTO ACCOUNT the time dimension @@@@@

    # for each level, calc wprime*dTvar/dz
    # Sulf
    wtotprime = (wvelp+wisopp)-(wvelc+wisopc)
    wvelprime = wvelp-wvelc
    wisopprime = wisopp-wisopc

    # GHGrem
    #wprime2= (wvp2+wip2)-(wvc+wic)
    #wvelprime2 = wvp2-wvc
    #wisopprime2 = wip2-wic


    tempcreg = ma.masked_where(rmask,tempc)
    temppreg = ma.masked_where(rmask,tempp)
    tareatreg = ma.masked_where(rmask,tareat)

    # PUT INTO A FUNCTION @@@
    # @@ if swap dims so lev is first, don't need to know # of dims
    # now also mask out cells below ocean floor:
    for lii,zz in enumerate(zt):
        # mask out levels below sea floor
        tempcreg[lii,...] = ma.masked_where(kmt <= lii,tempcreg[lii,...]) # @@@ do the masks combine? hope so.
        temppreg[lii,...] = ma.masked_where(kmt <= lii,temppreg[lii,...])
        tareatreg[lii,...] = ma.masked_where(kmt <= lii, tareatreg[lii,...])
        wtotprime[lii,...] = ma.masked_where(kmt <= lii,wtotprime[lii,...])
        wvelprime[lii,...] = ma.masked_where(kmt <= lii,wvelprime[lii,...])
        wisopprime[lii,...] = ma.masked_where(kmt <= lii,wisopprime[lii,...])

    wtotprimereg = ma.masked_where(rmask,wtotprime)
    wvelprimereg = ma.masked_where(rmask,wvelprime)
    wisopprimereg = ma.masked_where(rmask,wisopprime)

    #zonal mean: may not need the weights
    regzonalarea= ma.sum(tareatreg,axis=2) # only want to sum where there isn't land
    regzonalareat=np.tile(regzonalarea,(tareatreg.shape[2],1,1))
    regzonalareat=np.transpose(regzonalareat,(1,2,0))
    regzonalwgts= tareatreg/regzonalareat
    regzonalwgts=ma.masked_where(tareatreg.mask,regzonalwgts)

    tempcreg=np.squeeze(ma.average(tempcreg,axis=2,weights=regzonalwgts))
    temppreg=np.squeeze(ma.average(temppreg,axis=2,weights=regzonalwgts))

    # @@@ NOT weighted
    wtotprimereg = np.squeeze(ma.mean(wtotprimereg,axis=2))
    wvelprimereg = np.squeeze(ma.mean(wvelprimereg,axis=2))
    wisopprimereg = np.squeeze(ma.mean(wisopprimereg,axis=2))

    tbarreg = tempcreg# already zonal meaned MEAN T
    dtbarreg = ma.diff(tbarreg,axis=0) # delta of MEAN T with height

    # thickness of each layer
    dzt = np.diff(zt/100.) # convert to m
    tottransreg = ma.zeros((len(dzt),wtotprimereg.shape[1])) # initialize total heat trans (used to be wtransreg)
    wveltransreg = ma.zeros((len(dzt),wtotprimereg.shape[1])) # heat trans due to wvel (used to be wtranswvreg)
    wisoptransreg = ma.zeros((len(dzt),wtotprimereg.shape[1])) # (used to be wtranswireg)

    # -------------- Calc heat trans----
    # calc heat transport (heating rate) for each level
    for lii,dz in enumerate(dzt):
        #print 'ind: ' + str(lii) + ', dz: ' + str(dz)

        #  W prime * (dTbar / dz)
        tottransreg[lii,...] = wtotprimereg[lii,...]*(dtbarreg[lii,...]/dz)
        wveltransreg[lii,...] = wvelprimereg[lii,...]*(dtbarreg[lii,...]/dz) # WVEL only
        wisoptransreg[lii,...] = wisopprimereg[lii,...]*(dtbarreg[lii,...]/dz) # WISOP only



    onelat = tlat[:,1] # note this is actually one lon. values from -79 to +72

    rho_sw=cnc.getNCvar(filenamec,'rho_sw')
    cp_sw = cnc.getNCvar(filenamec,'cp_sw')
    rhocp = 1e-1*cp_sw*rho_sw # [J/K/m^3]

    # tile the layer thickness on T grid
    dzttile = np.tile(dzt,(tottransreg.shape[1],1))
    dzttile = np.transpose(dzttile)

    # multiply by layer thickness to get W/m2
    totwreg = tottransreg*rhocp*dzttile
    totwvreg = wveltransreg*rhocp*dzttile
    totwireg = wisoptransreg*rhocp*dzttile

    print totwreg.shape


    # ---------- average the PIG region with depth ------
    
     # trying to avoid cells too far south 
    Nlim=-65   # 65 to 74 is good for pig
    Slim=-74

    tareay = tarea[:,0] # because we are dealing w/ SH only, doesn't matter what lon we choose
    # now tile tareay for each depth
    tareayt = np.tile(tareay,(len(zt)-1,1)) # for transport (or when a dt or dz is involved)
    # create weights:
    totareay = ma.sum(tareayt[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
    totareayt = np.tile(totareay,(len(tareayt[0,np.logical_and(onelat<= Nlim,onelat>Slim)]),1))
    totareayt = np.transpose(totareayt,(1,0))
    wgts = tareayt[:,np.logical_and(onelat<= Nlim,onelat>Slim)] / totareayt

    #tareaytw = np.tile(tareay,(len(zt),1)) # for all 60 levels
    tareaytw = tareat[...,0] # already tile to 60 levels @@@@
    totareayw = ma.sum(tareaytw[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1) # @@@@ change to ma.sum
    totareaytw = np.tile(totareayw,(len(tareaytw[0,np.logical_and(onelat<= Nlim,onelat>Slim)]),1))
    totareaytw = np.transpose(totareaytw,(1,0))
    wgtsw = tareaytw[:,np.logical_and(onelat<= Nlim,onelat>Slim)] / totareaytw


    # meridional average with depth: vertical heat transport
    totwreg = ma.average(totwreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwvreg = ma.average(totwvreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwireg = ma.average(totwireg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)


    # meridional average w PRIME with depth
    wavgreg=ma.average(wtotprimereg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    wvavgreg=ma.average(wvelprimereg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    wiavgreg=ma.average(wisopprimereg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    
    # climo dTbar, averaged meridionally
    dTbaravgreg=ma.average(dtbarreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)




    # PLOT HEAT TRANS FOR PAPER ====================
    s2day = 60*60*24
    sec2yr = s2day*365

    mediumblue = ccm.get_linecolor('mediumblue') # Sulf
    dodgerblue = ccm.get_linecolor('dodgerblue') # GHGrem
    darkolivegreen3 = ccm.get_linecolor('darkolivegreen3') # GHGrem
    firebrick = ccm.get_linecolor('firebrick') # RCP8.5

    ylim=500

    # VERSION 2 has dTbar instead of vert heat through a layer ======= PAPER
    fig2,axs = plt.subplots(1,3,sharey=True)
    fig2.set_size_inches(14,4)
    ax = axs[2] #fig2.add_subplot(131,sharey=True)

    ax.plot(-1*totwreg,zt[1:]/100.,color=mediumblue,linewidth=4)
    #ax.plot(-1*totw2reg,zt[1:]/100.,color=darkolivegreen3,linewidth=4) # GHGrem
    ax.plot(-1*totwvreg,zt[1:]/100.,color=mediumblue,linewidth=2,linestyle='--')
    #ax.plot(-1*totwv2reg,zt[1:]/100.,color=darkolivegreen3,linewidth=2,linestyle='--')
    ax.plot(-1*totwireg,zt[1:]/100.,color=mediumblue,linewidth=2)#,linestyle=':')
    #ax.plot(-1*totwi2reg,zt[1:]/100.,color=darkolivegreen3,linewidth=2)#3,linestyle=':')

    yticks=np.arange(0,ylim,100)
    ax.plot([0,0],[0,1000],'k')
    ax.legend(('Sulf','GHGrem'),loc='best',fancybox=True,framealpha=0.5)
    ax.set_ylim((0,ylim))
    ax.set_xlim(-.15,.15)
    ax.set_xticks([-0.15,-0.10,-0.05, 0, 0.05, 0.1, 0.15])
    ax.set_xticklabels([-0.15,'',-0.05, 0, .05,'', 0.15], fontsize=18)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks,fontsize=18)
    #ax.set_title('VHT (W/m$^2$)',fontsize=18)
    ax.set_title('$\Delta$w*d$\overline{T}$/dz (W/m$^2$)',fontsize=18)
    ax.set_ylabel('Depth (m)',fontsize=18)
    #plt.title(region + ' Avg vert heat trans (W/m2) ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
    ax.invert_yaxis()

    ax2=axs[0] #fig2.add_subplot(132,sharey=True)
    ax2.plot(wavgreg*sec2yr,zt/100.,color=mediumblue,linewidth=4)
   # ax2.plot(wavg2reg*sec2yr,zt/100.,color=darkolivegreen3,linewidth=4)
    ax2.plot(wvavgreg*sec2yr,zt/100.,color=mediumblue,linewidth=2,linestyle='--')
    #ax2.plot(wvavg2reg*sec2yr,zt/100.,color=darkolivegreen3,linewidth=2,linestyle='--')
    ax2.plot(wiavgreg*sec2yr,zt/100.,color=mediumblue,linewidth=2)#,linestyle='-.')
    #ax2.plot(wiavg2reg*sec2yr,zt/100.,color=darkolivegreen3,linewidth=2)#,linestyle='-.')

    ax2.plot([0,0],[0,1000],'k')
    ax2.set_ylim((0,ylim))
    ax2.set_xlim(-10,8)
    ax2.set_xticks([-10,-8,-6,-4,-2, 0, 2, 4, 6, 8, 10])
    ax2.set_xticklabels([-10,'',-6,'',-2,0,2, '',6,'',10], fontsize=18)
    ax2.set_title('$\Delta$w (m/yr)',fontsize=18)
    ax2.invert_yaxis()

    ax3=axs[1] #fig2.add_subplot(133)
#    # POSITIVE means getting warmer with depth
#    plt.plot(dTbaravgreg,zt[1:]/100.,color='k',linewidth=3)
    # NEGATIVE means getting warmer with depth (if mult by -1)
    ax3.plot(-1*dTbaravgreg,zt[1:]/100.,color='k',linewidth=3) 

    ax3.plot([0,0],[0,1000],'k')
    ax3.set_ylim((0,ylim))
    #plt.xlim(-10,8)
    ax3.set_xticks([-0.3,-0.2,-0.1,0,0.1]) # for when dT/dz is neg for warm with depth
    ax3.set_xticklabels([-0.3,-0.2,-0.1,0,0.1], fontsize=18)
    ax3.set_title('d$\overline{T}$/dz ($^\circ$C/m)',fontsize=18)
    ax3.invert_yaxis()

    if printtofile:
        fig2.savefig('WAISpap_vertheattrans_wvels_dTbarnegwrmwithdepth_' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 
                     'S_ylim' + str(ylim) + '_' + region + '.pdf')


if oceants:


    # /Volumes/MyPassport2TB/DATA/ccsm4/b40.20th.track1.1deg.006/ocn_proc/TEMP/TEMP.b40.20th.track1.1deg.006.pop.h.1970-1999.nc
    var='TEMP'
    reg='SH' # SH or PIG
    lonlims = [230,280]; strlims='80W-120W' # for PIG

    latlim=-50 # average south of 50s

    basepath='/Volumes/MyPassport2TB/DATA/ccsm4/'
    basepath2='/Users/kelly/School/DATA/'

    casenamec='b40.20th.track1.1deg.006'
    casenamep='geo2035ensavg'
    casenamep2='rcp8_5GHGrem1850'
    casenamep3='b40.rcp8_5.1deg.006'

    timeperc='1970-2005'
    timeperp='2035-2094'
    timeperp3='2005-2034'
    regstr=''

    if reg=='PIG':
        reg='SH'
        pig=True
    else:
        pig=False

    filenameclim = basepath2 + casenamec + '/' + casenamec + '.pop.ANN.1970-1999.nc' 
    
    filenamec = basepath + casenamec + '/ocn_proc/' + var + '/' + var + '.' + casenamec + '.pop.' + timeperc + '.' + reg + '.nc' # ctl
    filenamep = basepath + casenamep + '/ocn_proc/' + var + '/' + var + '.' + casenamep  + '.pop.' + timeperp + '.' + reg + '.nc' # sulf
    filenamep3 = basepath + casenamep3 + '/ocn_proc/' + var + '/' + var + '.' + casenamep3  + '.pop.' + timeperp3 + '.' + reg + '.nc' # rcp


    print filenamec
    print filenamep

    kmt = cnc.getNCvar(filenameclim,'KMT') # get from climo
    tarea=cnc.getNCvar(filenameclim,'TAREA') # get from climo
    tlatc=cnc.getNCvar(filenameclim,'TLAT')
    tlonc=cnc.getNCvar(filenameclim,'TLONG')

    zt = cnc.getNCvar(filenamec, 'z_t')
    tlat = cnc.getNCvar(filenamec,'TLAT') # these are the shape of the var
    tlon = cnc.getNCvar(filenamec,'TLONG') 
    onelat=tlat[:,0]
    
    fldclimo = cnc.getNCvar(filenameclim,var)
    fldc = cnc.getNCvar(filenamec,var)
    fldp = cnc.getNCvar(filenamep,var)
    fldp3 = cnc.getNCvar(filenamep3,var)

    fldclimoreg = fldclimo[:,onelat<latlim,:]
    fldcreg = fldc[:,:,onelat<latlim,:]
    fldpreg = fldp[:,:,onelat<latlim,:]
    fldp3reg = fldp3[:,:,onelat<latlim,:]
    tareareg = tarea[onelat<latlim,:]
    kmtreg = kmt[onelat<latlim]
    
    (ntime,nlev,nlat,nlon) = fldcreg.shape
    (ntimep,nlevp,nlatp,nlonp) = fldpreg.shape
    (ntimep3,nlevp3,nlatp3,nlonp3) = fldp3reg.shape

    onelatarea = tareareg[onelat<latlim,0]
    totarea = np.sum(onelatarea)
    wgts = onelatarea / totarea
    # tile wgts with depth
    wgtst = np.tile(wgts,(nlev,1))
    
    fldclimozm=ocnzonalmean_ccsm4(fldclimoreg,tareareg,kmtreg)
    fldczm=ocnzonalmean_ccsm4(fldcreg,tareareg,kmtreg)
    fldpzm=ocnzonalmean_ccsm4(fldpreg,tareareg,kmtreg)
    fldp3zm=ocnzonalmean_ccsm4(fldp3reg,tareareg,kmtreg)    

    # now average the zonal means with latitude. tile wgts with time for each
    fldclimoavg = np.average(fldclimozm,axis=1,weights=wgtst)
    fldcavg = np.average(fldczm,axis=2,weights=np.tile(wgtst,(ntime,1,1)))
    fldpavg = np.average(fldpzm,axis=2,weights=np.tile(wgtst,(ntimep,1,1)))
    fldp3avg = np.average(fldp3zm,axis=2,weights=np.tile(wgtst,(ntimep3,1,1)))

    print fldcavg.shape

    # CALCULATE SIGNIFICANCE
    fldcstats = fldcavg[:-6,...]
    #fldpstats = np.concatenate((fldp3avg,fldpavg),axis=0)
    fldpstats = np.concatenate((fldcavg,fldp3avg[1:,:],fldpavg),axis=0)

    tt=5
    pvalsall=np.zeros((fldpstats.shape[0]-tt,60))
    for ii in np.arange(tt,fldpstats.shape[0]):
        # use the surrounding 10 years to estimate that grid point's significance
        test = fldpstats[ii-tt:ii+tt+1,:] #make the year the midpoint of 11-year period
        print ii
        print test.shape

        tstats,pvalsall[ii-tt,:]=cutl.ttest_ind(fldpstats[ii-tt:ii+tt+1,:],fldcstats)
    
   # pvals=np.zeros((60-tt,60))
   # pvals3=np.zeros((30,60))

   # for ii in np.arange(tt,60):
   #     # use the surrounding 10 years to estimate that grid point's significance
   #     tstats,pvals[ii-tt,:]=cutl.ttest_ind(fldpavg[ii-tt:ii+tt,:],fldcstats)
   # for ii in np.arange(0,30):
   #     tstats,pvals3[ii,:]=cutl.ttest_ind(fldp3avg[ii,:],fldcstats)

    # concat the timeseries
    plotfld = np.concatenate((fldcavg,fldp3avg[1:,:],fldpavg),axis=0) - fldclimoavg
    #plotpvals = np.concatenate((pvals3,pvals),axis=0)
    plotpvals=pvalsall

    # ========= PLOT Southern OCN TEMP with time ========
    xx=np.arange(1970,2095)
    #xxpvals=np.arange(2005,2095)
    #xxpvals=np.arange(2005+tt,2095)
    xxpvals=np.arange(1970+tt,2095)

    zlevs,times = np.meshgrid(zt/100.,xx)
    zlevspv,timespv = np.meshgrid(zt/100.,xxpvals)

    cmax=.5; cmin=-.5
    cmap='blue2red_20'; cmlen=20.
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)


    fig,axs = plt.subplots(2,1)
    fig.set_size_inches(8,5)
    ax=axs[0]

    CF = ax.contourf(times,zlevs,plotfld,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
    cplt.addtsig(ax,plotpvals.T,xxpvals,zt/100.,siglevel=siglev,type='hatch')

    ax.set_ylim((0,1000))
    ax.invert_yaxis()
    ax.set_yticks(np.arange(0,1500,500))
    ax.set_yticklabels([0,0.5,1],fontsize=18)
    #ax.set_xlim(xlims)
    #ax.set_xticks(np.arange(-75,-45,5))
    ax.set_xticklabels("")
    #ax.set_title(regstr + 'Sulf',fontsize=18)
    ax.set_ylabel('Depth (m)',fontsize=18)
    
    ax=axs[1]

    CF = ax.contourf(times,zlevs,plotfld,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
    cplt.addtsig(ax,plotpvals.T,xxpvals,zt/100.,siglevel=siglev,type='hatch')

    ax.set_ylim((1000,5000))
    ax.invert_yaxis()
    ax.set_yticks(np.arange(2000,6000,1000))
    ax.set_yticklabels([2,3,4,5],fontsize=18)
    ax.set_xlabel('Years',fontsize=18)

    cbar_ax = fig.add_axes([.91,.15, .02,.7])
    fig.colorbar(CF,cax=cbar_ax)

    if printtofile:
        fig.savefig('WAISpap_SuppFig_shtocnTEMPwithtime.png',dpi=300)
