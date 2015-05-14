# new figs for WAIS paper, with stat sig
#  5/9/2015

import cccmautils as cutl
import numpy.ma as ma

plt.close('all')

printtofile=True
atmos=False
ocean=True

basepath2 = '/Users/kelly/School/DATA/'
basepath = '/Volumes/MyPassport2TB/DATA/ccsm4/'
casenamec = 'b40.20th.track1.1deg.006'; timeperc = '1970-1999'
casenamep = 'geo2035ensavg'; timeperp = '2045-2054'
casenamep2 = 'rcp8_5GHGrem1850'; timeperp = '2045-2054' # or 'rcp8_5GHGrem1850'
casenamep3 = 'b40.rcp8_5.1deg.006'; timeperp = '2045-2054' # or 'rcp8_5GHGrem1850'

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



    # ================= ZONAL WIND STRESS ======================

    var='TAUX' # wind stress

    fnamec =  basepath + casenamec + '/atm_proc/' + var + '/' + var + '.' + casenamec + '.cam2.h0.' + timeperc + '_timeseries.nc'
    fnamep =  basepath + casenamep + '/atm_proc/' + var + '/' + var + '.' + casenamep + '.cam2.h0.' + timeperp + '_timeseries.nc'
    fnamep2 =  basepath + casenamep2 + '/atm_proc/' + var + '/' + var + '.' + casenamep2 + '.cam2.h0.' + timeperp + '_timeseries.nc'
    fnamep3 =  basepath + casenamep3 + '/atm_proc/' + var + '/' + var + '.' + casenamep3 + '.cam2.h0.' + timeperp + '_timeseries.nc'


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

    fname = basepath + 'wscurl.b40.20th.track1.1deg.006.cam2.1970-1999.nc'

    fldc = np.squeeze(cnc.getNCvar(fname,'wscurl'))*-1 
    lat = cnc.getNCvar(fname,'lat')
    lon = cnc.getNCvar(fname,'lon')
    print fldc.shape

    # wind stress files
    fnamews = basepath + '../b40.20th.track1.1deg.006_ANN_climo.nc' # 1970-1999
    wsfldc = np.squeeze(cnc.getNCvar(fnamews,'TAUX'))*-1 # (also mult by -1 ?) Yes, atmos field

    fnames = {'Sulf': fnamep, # 2045-2054
              'GHGrem': fnamep2, # 2045-2054
              'RCP8.5': fnamep3} # 2045-2054

    # wind stress curl files
    wscfnames = {'Sulf': basepath + 'wscurl.geo2035ensavg.cam2.2045-2054.nc',
              'GHGrem': basepath + 'wscurl.rcp8_5GHGrem1850.cam2.2045-2054.nc',
              'RCP8.5': basepath + 'wscurl.b40.rcp8_5.1deg.006.cam2.2045-2054.nc' }


    coldt = {'RCP8.5': ccm.get_linecolor('firebrick1'),
             'Sulf': ccm.get_linecolor('mediumblue'),
             'GHGrem': ccm.get_linecolor('darkolivegreen3') }


    # @@@@@@@@@ ended here. Need to add compute significance and add to figs.

    fig,axs = plt.subplots(1,2,sharex=True)
    fig.set_size_inches(14,3.2) # match zonal mean ocean TEMP
    ax=axs[0] # wind stress
    for fname in ('RCP8.5','Sulf','GHGrem'):
        fld = np.squeeze(cnc.getNCvar(fnames[fname],'TAUX'))*-1 
        fldd = fld-wsfldc
        flddwt = fldd*ofrac
        flddwtzm = np.mean(flddwt,axis=1)

        ax.plot(lat,flddwtzm,color=coldt[fname],linewidth=3)

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

    ax=axs[1] # wind stress curl
    for fname in ('RCP8.5','Sulf','GHGrem'):
        fld = np.squeeze(cnc.getNCvar(wscfnames[fname],'wscurl'))*-1 
        fldd = fld-fldc
        flddwt = fldd*ofrac
        flddwtzm = np.mean(flddwt,axis=1)

        ax.plot(lat,flddwtzm,color=coldt[fname],linewidth=3)

    ax.set_xlim(-75,-40)
    ax.axhline(y=0,color='k')
    ax.set_xticklabels(['','70$^\circ$S','','60$^\circ$S','',\
                        '50$^\circ$S','','40$^\circ$S'],fontsize=18)
    ax.set_yticks(np.arange(-3e-8,2.5e-8,.5e-8))
    ax.set_yticklabels([-3,'',-2,'',-1,'',0,'',1,'',2],fontsize=18)
    ax.set_title(r'$\Delta$ ( $\nabla$ x TAU ) (10$^{-8}$ N/m$^3$)',fontsize=18)

    if printtofile:
        fig.savefig('TAUX_curlTAUX_allruns_zonmean_ANN.pdf')





# annual average with nc tools:
# # Annual average (use the feature of 'Duration')
#      ncra -O --mro -d time,"1956-01-01 00:00:0.0","2005-12-31 23:59:9.9",12,12
# my version failed b/c nco version too old: no mro option
# ncra -O --mro -d time "1970-01-01 00:00:0.0","1979-12-31 23:59:9.9",12,12 b40.20th.track1.1deg.006.pop.h.WISOP.197001-197912.nc b40.20th.track1.1deg.006.pop.h.WISOP.1970-1979_nctimeseries.nc
# not using cdo because it gave me 11 outputs! there should only be 10....


# ############# OCEAN #####################
if ocean:


    def ocnzonalmean_ccsm4(fld,tarea,kmt):


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
            # first mask out levels below sea floor
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





    # /Volumes/MyPassport2TB/DATA/ccsm4/b40.20th.track1.1deg.006/ocn_proc/TEMP/TEMP.b40.20th.track1.1deg.006.pop.h.1970-1999.nc
    var='TEMP'
    reg='PIG' # SH or PIG
    lonlims = [230,280]; strlims='80W-120W' # for PIG

    basepath='/Volumes/MyPassport2TB/DATA/ccsm4/'
    basepath2='/Users/kelly/School/DATA/'
    casenamec='b40.20th.track1.1deg.006'
    casenamep='geo2035ensavg'

    timeperc='1970-1999'
    timeperp='2045-2054'

    if reg=='PIG':
        reg='SH'
        pig=True
    else:
        pig=False

    filenamec = basepath + casenamec + '/ocn_proc/' + var + '/' + var + '.' + casenamec + '.pop.' + timeperc + '.' + reg + '.nc' 
    filenameclim = basepath2 + casenamec + '/' + casenamec + '.pop.ANN.' + timeperc + '.nc' 

    filenamep = basepath + casenamep + '/ocn_proc/' + var + '/' + var + '.' + casenamep  + '.pop.' + timeperp + '.' + reg + '.nc' 
    filenamep2 = basepath + casenamep2 + '/ocn_proc/' + var + '/' + var + '.' + casenamep2  + '.pop.' + timeperp + '.' + reg + '.nc' 

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
    cplt.addtsig(ax,pvals,np.squeeze(tlat[:,1]),zt/100.,siglevel=siglev)

    ax.set_ylim(ylims)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(0,900,100))
    ax.set_yticklabels([0,'',200,'',400,'',600,'',800],fontsize=18)
    ax.set_xlim(xlims)
    ax.set_xticks(np.arange(-75,-45,5))
    ax.set_xticklabels(['75$^\circ$S', '70$^\circ$S', '65$^\circ$S', \
                        '60$^\circ$S', '55$^\circ$S', '50$^\circ$S'],fontsize=18)
    ax.set_title('Sulf',fontsize=18)
    ax.set_ylabel('Depth (m)',fontsize=18)

    ax2=axs[1]
    cmax=.5; cmin=-.5
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    CF2 = ax2.contourf(tlats,zlevs,plotfld2,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
    cplt.addtsig(ax2,pvals2,np.squeeze(tlat[:,1]),zt/100.,siglevel=siglev)#,type='cont')

    ax2.set_ylim(ylims)
    ax2.invert_yaxis()
    ax2.set_yticks(np.arange(0,900,100))
    ax2.set_yticklabels([0,'',200,'',400,'',600,'',800],fontsize=18)
    ax2.set_xlim(xlims)
    ax2.set_xticks(np.arange(-75,-45,5))
    ax2.set_xticklabels(['75$^\circ$S', '70$^\circ$S', '65$^\circ$S', \
                        '60$^\circ$S', '55$^\circ$S', '50$^\circ$S'],fontsize=18)
    #ax2.set_title(casenamep2)
    ax2.set_title('GHGrem',fontsize=18)
    cbar_ax = fig.add_axes([.91,.15, .02,.7])
    fig.colorbar(CF,cax=cbar_ax)

    if printtofile:
        fig.savefig('TEMPanom_subplot' + reg + '_ylim' + str(ylims[1]) + 'xlim' + str(xlims[1]) + '_c_sig' + str(1-siglev) + '.png',dpi=400)




