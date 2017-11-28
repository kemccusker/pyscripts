"""
    rlx_utils.py
    10/24/2016

    This class is for "nudging" or relaxation runs.

"""

import cccmautils as cutl
import numpy as np
import loadCanESM2rlxdata as lrlx
import cccmaplots as cplt


def rlx_diffs(ncdt,dotimmean=False):
    """ Take a dictionary of all the simulation data and return
        a dictionary of the differences as analyzed in paper (e.g. ICEcold, ICEwarm, etc)
    """
    keys=('ICEcold','CO2hi','Full','ICEwarm','CO2lo')
    diffcases={'ICEcold':('prei2xco2iceb','preipreiice'),
              'CO2hi': ('2xco2preiice','preipreiice'),
              'Full': ('2xco22xco2ice','preipreiice'),
              'ICEwarm':('2xco22xco2ice','2xco2preiice'),
              'CO2lo': ('2xco22xco2ice','prei2xco2iceb')}
    retdt={}
    for dkey in keys:
        sims = diffcases[dkey]
        if dotimmean:
            retdt[dkey] = ncdt[sims[0]].mean(axis=0) - ncdt[sims[1]].mean(axis=0)
        else:
            retdt[dkey]=ncdt[sims[0]]-ncdt[sims[1]]
        
    return retdt


def calc_linpatterncorr(ncdt,lin='one',northof=0,conv=1,vert=False):
    
    """ Computes the pattern correlation b/w the 'SUM' and 'FULL'
           for given 'lin'. Input data is seasonal avg timeseries of spatial data.
           
           computes the correlation for area north of 'northof' value.
           
           returns rval,pval
           
    """
    
    pico = ncdt['prei2xco2iceb']
    pipi = ncdt['preipreiice']
    copi = ncdt['2xco2preiice']
    coco = ncdt['2xco22xco2ice']
    lat = ncdt['lat']

    if lin=='one':
        ice = (np.mean(pico,axis=0) - np.mean(pipi,axis=0))*conv
        co2 = (np.mean(copi,axis=0) - np.mean(pipi,axis=0))*conv
        ctl = np.mean(pipi,axis=0)*conv        
        suff=''
        suff1='cold'; suff2='hi'
    elif lin=='two':
        ice = (np.mean(coco,axis=0) - np.mean(copi,axis=0))*conv
        co2 = (np.mean(coco,axis=0) - np.mean(pico,axis=0))*conv
        ctl = np.mean(coco,axis=0)*conv
        suff='2'
        suff1='warm'; suff2='lo'

    combo = ice+co2
    full = (np.mean(coco,axis=0) - np.mean(pipi,axis=0))*conv
    
    if vert: # lat in 2nd dim
        rval,pval = cutl.pattcorr_pearson(combo[:,lat>northof].flatten(),
                                          full[:,lat>northof].flatten())
    else: # lat in 1st dim
        rval,pval = cutl.pattcorr_pearson(combo[lat>northof,:].flatten(),
                                          full[lat>northof,:].flatten())
    #print cutl.pattcorr(combo[lat>northof,:].flatten(),full[lat>northof,:].flatten())
    
    return rval,pval

def calc_monthly_linpatterncorr(ncdt,lin='one',northof=0,conv=1,vert=False):
    
    """ Computes the pattern correlation b/w the 'SUM' and 'FULL'
           for given 'lin'. Input data is climatology of spatial data.
           
           computes the correlation for area north of 'northof' value for e/ month.
           
           returns rval,pval
           
    """
    
    pico = ncdt['prei2xco2iceb']
    pipi = ncdt['preipreiice']
    copi = ncdt['2xco2preiice']
    coco = ncdt['2xco22xco2ice']
    lat = ncdt['lat']
    
    if lin=='one':
        ice = (pico - pipi)*conv
        co2 = (copi - pipi)*conv
        ctl = pipi*conv        
        suff=''
        suff1='cold'; suff2='hi'
    elif lin=='two':
        ice = (coco - copi)*conv
        co2 = (coco - pico)*conv
        ctl = coco*conv
        suff='2'
        suff1='warm'; suff2='lo'

    combo = ice+co2
    full = (coco - pipi)*conv
    
    #print lat>northof
    rvals=np.zeros((12,))
    pvals=np.zeros_like(rvals)
    
    for moii in range(0,12):
        if vert: # lat in 2nd dim
            rvals[moii],pvals[moii] = cutl.pattcorr_pearson(combo[moii,:,lat>northof].flatten(),
                                                            full[moii,:,lat>northof].flatten())
        else: # lat in 1st dim
            rvals[moii],pvals[moii] = cutl.pattcorr_pearson(combo[moii,lat>northof,:].flatten(),
                                                            full[moii,lat>northof,:].flatten())
    
    return rvals,pvals
    
def calc_monthly_rmse(ncdt,lin='one',norm=True,northof=0,conv=1,vert=False):
    
    """ Computes the average RMSE b/w the 'SUM' and 'FULL'
           for given 'lin'. Input data is climatology of spatial data.
           
           computes the RMSE for area north of 'northof' value for e/ month.
           
           if norm=True: normalize by spatial variance
           
           returns rmse
           
    """
    
    pico = ncdt['prei2xco2iceb']
    pipi = ncdt['preipreiice']
    copi = ncdt['2xco2preiice']
    coco = ncdt['2xco22xco2ice']

    if lin=='one':
        ice = (pico - pipi)*conv
        co2 = (copi - pipi)*conv
        ctl = pipi*conv        
        suff=''
        suff1='cold'; suff2='hi'
    elif lin=='two':
        ice = (coco - copi)*conv
        co2 = (coco - pico)*conv
        ctl = coco*conv
        suff='2'
        suff1='warm'; suff2='lo'

    combo = ice+co2
    full = (coco - pipi)*conv
    
    rmse=np.zeros((12,))
    spvar = np.zeros_like(rmse)
    for moii in range(0,12):
        if vert: # lat in 2nd dim
            print 'how to weight vertical coord? dp?'
            #rmse[moii] = cutl.calc_areawgted_rmse(combo[moii,:,lat>northof].flatten(),
            #                                                full[moii,:,lat>northof].flatten())
        else: # lat in 1st dim
            #print 'testing area weighted var func @@@'
            spvar[moii] = cutl.calc_areawgted_var(full[moii,lat>northof,:-1],
                                            lat[lat>northof],lon[:-1],model='CanESM2')
            
            rmse[moii] = cutl.calc_areawgted_rmse(full[moii,lat>northof,:-1],
                                                  combo[moii,lat>northof,:-1],
                                                  lat[lat>northof],lon[:-1],model='CanESM2')
            if norm:
                rmse[moii] = rmse[moii]/spvar[moii]
            
    print 'spvar ' +str(spvar) # @@
            
    return rmse

def climatologize_ncdt(ncdt):
    
    climdt={}
    for ckey in ncdt.keys():
        dat = ncdt[ckey]
        cdat,_ = cutl.climatologize(dat)
        climdt[ckey] = cdat
        
    return climdt


##### plotting funcs

def plot_nc_linearity_maps(ncdt, axs, lin='one',cmin='',cmax='',cmind='',cmaxd='',cmin2='',cmax2='',
                           conv=1, ptype='nh',suppcb=False,cmap='blue2red_w20',subtime=None,
                           screen=False,vert=False,levlim=None,suppttl=False,addsig=False,sigtype='cont',
                           latlim=None,vertptype=None,addclimcont=False,ctlconts=None,
                           nosum=False,nofull=False,nolin=False,scaledlin=False,ncicedt=None,sigreverse=False,
                           northof=0,fsz=None):
    """ 
        cmin/cmax:  for ice
        cmin2/cmax2: for co2 & sum & full, if given. otherwise same as ice
        cmind/cmaxd: for lin combo subtraction from full
        
        lin = 'one' is ice and co2 (ICEcold, CO2hi)
            = 'two' is ice2 and co22 (ICEwarm, CO2lo)
            = 'sub' subtract lin 'one' from lin 'two' (shows how warm clim/low ice different from cold/high)
            
        nosum: include ice+co2 panel or not?
        nofull: include Full panel or not?
        nolin: include non-linearity panel or not?
        scaledlin: is the non-linearity panel to be scaled or not?
        ncicedt: must be set if scaledlin=True!
        
        sigreverse: if True, plot hatching/contours where NOT significant on
                    every map but non-lin. (only if addsig=True)
        northof: compute pattern correlation b/w SUM and FULL north of this latitude.
        fsz: fontsize for title
                    
        returns: plot handles (ph,ph2[,phd]) (ice,co2|sum|full[,difference/nonlin])
    """
    
    
    # @@ subtime not implemented yet. assume average over full time 
    #tmplen = ncdt['2xco2preiice'].shape[0]-1
    fmt='%2.1f' # clabel format
    #if ctlconts!=None:
    #    if (ctlconts < 1.).any() and (ctlconts > -1.).any():
    #    #if np.logical_and((ctlconts<1).any(),(ctlconts>-1).any()):
    #        print ctlconts
    #        fmt='%2.1f'
    
    pico = ncdt['prei2xco2iceb']
    pipi = ncdt['preipreiice']
    copi = ncdt['2xco2preiice']
    coco = ncdt['2xco22xco2ice']
    lat = ncdt['lat']
    
    if lin=='one':
        _,icepv = cutl.ttest_ind(pico, pipi,axis=0,effdof=False)        
        ice = (np.mean(pico,axis=0) - np.mean(pipi,axis=0))*conv
        #icep = ncfldzmdt['pi2xco2ipulse'] - ncfldzmdt['preipreiice']    
        _,co2pv = cutl.ttest_ind(copi, pipi,axis=0,effdof=False)   
        co2 = (np.mean(copi,axis=0) - np.mean(pipi,axis=0))*conv
        _,combopv = cutl.ttest_ind(pico-pipi+(copi-pipi),
                                   coco-pipi)
        ctl = np.mean(pipi,axis=0)*conv
        
        suff=''
        suff1='cold'; suff2='hi'
    elif lin=='two':
        _,icepv = cutl.ttest_ind(coco, copi,axis=0,effdof=False)
        ice = (np.mean(coco,axis=0) - np.mean(copi,axis=0))*conv
        _,co2pv = cutl.ttest_ind(coco, pico,axis=0,effdof=False)
        co2 = (np.mean(coco,axis=0) - np.mean(pico,axis=0))*conv
        _,combopv = cutl.ttest_ind(coco-copi+(coco-pico),
                                   coco-pipi)
        ctl = np.mean(coco,axis=0)*conv
        suff='2'
        suff1='warm'; suff2='lo'
    elif lin=='sub':
        #two - one: (coco-copi) - (pico-pipi)
        _,icepv = cutl.ttest_ind(coco-copi, pico-pipi,axis=0,effdof=False)
        ice = (np.mean(coco-copi,axis=0) - np.mean(pico-pipi,axis=0))*conv
        #two - one: (coco-pico) - (copi-pipi)
        _,co2pv = cutl.ttest_ind(coco-pico, copi-pipi,axis=0,effdof=False)
        co2 = (np.mean(coco-pico,axis=0) - np.mean(copi-pipi,axis=0))*conv
        suff='2-1'
        suff1='warm-cold'; suff2='lo-hi'
        
        
    _,fullpv = cutl.ttest_ind(coco, pipi,axis=0,effdof=False)
    full = (np.mean(coco,axis=0) - np.mean(pipi,axis=0))*conv
    
    # compute pattern correlation b/w SUM and FULL
    prval,ppval = calc_linpatterncorr(ncdt,lin=lin,northof=northof,vert=vert)
    pvarexp = (prval**2)*100
    
    # pparams for ice!
    # pparams2 for co2![sum|full]
    if vert:
        lev = ncdt['lev']
        
        pparams = {'cmin': cmin, 'cmax': cmax,
                   'suppcb': suppcb, 'screen': screen, 'cmap':cmap,'levlim':levlim,
                   'latlim':latlim,'ptype':vertptype}
        if cmin2 != '':
            pparams2 = {'cmin': cmin2, 'cmax': cmax2,
                       'suppcb': suppcb, 'screen': screen, 'cmap':cmap,'levlim':levlim,
                       'latlim':latlim,'ptype':vertptype}
        else:
            pparams2 = {'cmin': cmin, 'cmax': cmax,
                       'suppcb': suppcb, 'screen': screen, 'cmap':cmap,'levlim':levlim,
                       'latlim':latlim,'ptype':vertptype}
    else:
        lon = ncdt['lon']
        
        pparams = {'ptype': ptype, 'cmin': cmin, 'cmax': cmax,
                   'suppcb': suppcb, 'cmap':cmap}
        if cmin2 != '':
            pparams2 = {'ptype': ptype, 'cmin': cmin2, 'cmax': cmax2,
                       'suppcb': suppcb, 'cmap':cmap}
        else:
            pparams2 = {'ptype': ptype, 'cmin': cmin, 'cmax': cmax,
                       'suppcb': suppcb, 'cmap':cmap}
    
    aii=0
    ax=axs[aii]
    if suppttl: ttl=''
    else: ttl='ICE'+suff1
    if vert:
        ph = cplt.vert_plot(ice,lev,lat,axis=ax,title=ttl,**pparams)
        if addsig:
            cplt.addtsig(ax,icepv,lat,lev/100.,sigtype=sigtype,reverse=sigreverse)
        if addclimcont:
            cplt.add_contoursvert(ax,ctl,lat,lev,verb=True,clab=True,conts=ctlconts,fmt=fmt)
    else:
        bm,ph = cplt.kemmap(ice,lat,lon,axis=ax,title=ttl,**pparams)
        if addsig:
            cplt.addtsigm(bm,icepv,lat,lon,sigtype=sigtype,reverse=sigreverse)
    if fsz!=None:
        ax.set_title(ttl,fontsize=fsz)
    aii+=1    
    
    ax=axs[aii]
    if suppttl: ttl=''
    else: ttl='CO2'+suff2
    if vert:
        ph2 = cplt.vert_plot(co2,lev,lat,axis=ax,title=ttl,suppylab=True,**pparams2)
        if addsig:
            cplt.addtsig(ax,co2pv,lat,lev/100.,sigtype=sigtype,reverse=sigreverse)
        if addclimcont:
            cplt.add_contoursvert(ax,ctl,lat,lev,verb=True,clab=True,conts=ctlconts,fmt=fmt)
    else:
        bm,ph2 = cplt.kemmap(co2,lat,lon,axis=ax,title=ttl,**pparams2)
        if addsig:
            cplt.addtsigm(bm,co2pv,lat,lon,sigtype=sigtype,reverse=sigreverse)
    if fsz!=None:
        ax.set_title(ttl,fontsize=fsz)
    aii+=1
    
    if not nosum:
        ax=axs[aii]
        if suppttl: ttl=''
        else: ttl='Sum' #ice'+suff1+'+co2'+suff2
        if vert:
            ph2 = cplt.vert_plot(ice+co2,lev,lat,axis=ax,title=ttl,suppylab=True,**pparams2)
            if addclimcont:
                cplt.add_contoursvert(ax,ctl,lat,lev,verb=True,clab=True,conts=ctlconts,fmt=fmt)
        else:
            bm,ph2 = cplt.kemmap(ice+co2,lat,lon,axis=ax,title=ttl,**pparams2)
        if fsz!=None:
            ax.set_title(ttl,fontsize=fsz)            
        aii+=1

    if not nofull:
        ax=axs[aii]
        if suppttl: ttl=''
        else: ttl='Full'
        if vert:
            ph2 = cplt.vert_plot(full,lev,lat,axis=ax,title=ttl,suppylab=True,**pparams2)
            if addsig:
                cplt.addtsig(ax,fullpv,lat,lev/100.,sigtype=sigtype,reverse=sigreverse)
            if addclimcont:
                cplt.add_contoursvert(ax,ctl,lat,lev,verb=True,clab=True,conts=ctlconts,fmt=fmt)
        else:
            bm,ph2 = cplt.kemmap(full, lat,lon,axis=ax,title=ttl,**pparams2)
            if addsig:
                cplt.addtsigm(bm,fullpv,lat,lon,sigtype=sigtype,reverse=sigreverse)
            ax.annotate('$%.0f$'%(pvarexp) + '%', xy=(0.78,0.98),xycoords='axes fraction',
                       fontsize=16)
        if fsz!=None:
            ax.set_title(ttl,fontsize=fsz)
        aii+=1    
    
    if not nolin:
        ax=axs[aii] # the linearity test
        # @@@ add scaledlin here @@@
        # e.g.  c1*ICEcold + c2*CO2hi - Full
        #       c1*(pico-pipi) + c2*(copi-pipi) - (coco-pipi)
        #       solve c1, c2 with sea ice s.t. above eqn would give 0
        # @@@ Question: how to scale? how to get c1,c2? Use difference from target sea ice?
        if scaledlin:
            # compute target offset
            offs = compute_targoffset(ncicedt,dosia=True)
            # want to be a fraction so div by 100.
            # and for ice comparisons, we don't reach target so we need to bump it up (add 1)
            c1 = offs['ice'+suff1]/100. + 1
            # for co2 comparison, not sure how to scale, b/c want to remove ice effect, which isn't a scaling?
            c2 = offs['co2'+suff2]/100.
    
        if suppttl: ttl=''
        else: ttl='Full-Sum' #ttl='(ice'+suff +'+co2'+suff+')-full'
        if cmind=='': # just divide the other clims by 2
            if vert:
                phd=cplt.vert_plot(full-(ice+co2),lev,lat,axis=ax,title=ttl,
                            cmin=cmin/2,cmax=cmax/2,suppcb=suppcb,screen=screen,cmap=cmap,
                            suppylab=True,levlim=levlim,latlim=latlim,ptype=vertptype)
                if addsig:
                    cplt.addtsig(ax,combopv,lat,lev,sigtype=sigtype)
                if addclimcont: # want in the linearity fig?? @@@@
                    cplt.add_contoursvert(ax,ctl,lat,lev,verb=True,clab=True,conts=ctlconts,fmt=fmt)
            else:
                bm,phd=cplt.kemmap(full-(ice+co2),lat,lon,axis=ax,title=ttl,ptype=ptype,
                            cmin=cmin/2,cmax=cmax/2,suppcb=suppcb,cmap=cmap)
                if addsig:
                    cplt.addtsigm(bm,combopv,lat,lon,sigtype=sigtype)
        else:
            if vert:
                phd=cplt.vert_plot(full-(ice+co2),lev,lat,axis=ax,title=ttl,
                            cmin=cmind,cmax=cmaxd,suppcb=suppcb,screen=screen,cmap=cmap,
                            suppylab=True,levlim=levlim,latlim=latlim,ptype=vertptype)
                if addsig:
                    cplt.addtsig(ax,combopv,lat,lev,sigtype=sigtype)
                if addclimcont: # want in the linearity fig?? @@@@
                    cplt.add_contoursvert(ax,ctl,lat,lev,verb=True,clab=True,conts=ctlconts,fmt=fmt)
            else:
                bm,phd=cplt.kemmap(full-(ice+co2),lat,lon,axis=ax,title=ttl,ptype=ptype,
                            cmin=cmind,cmax=cmaxd,suppcb=suppcb,cmap=cmap)
                if addsig:
                    cplt.addtsigm(bm,combopv,lat,lon,sigtype=sigtype)
        if fsz!=None:
            ax.set_title(ttl,fontsize=fsz)
    
        return ph,ph2,phd
    else:
        return ph,ph2
    
    
def plot_sims_quad(nhsies,gms,toffsets,coconc=False,fullann=True,printtofile=False):
    dosia=True
    fsz=20
    
    # Now make the plot
    fig,ax = plt.subplots(1,1,figsize=(10,8))

    #print nhsies # @@@
    #print gms # @@@
    
    
    for ii,case in enumerate(casenames):
        if case=='pi2xco2ipulse':
            pass
        else:
            ax.plot(nhsies[case]/1e12,gms[case], marker='s',color=colors[case])
            if fullann:
                ax.annotate(shortnames[case],(nhsies[case]/1e12,gms[case]),
                            xytext=(4,10),textcoords='offset points',fontweight='bold')

    # add ICEwarm+CO2lo point to fig
    sumvecx = (nhsies['prei2xco2iceb']-nhsies['2xco22xco2ice']) +\
              (nhsies['2xco2preiice']-nhsies['2xco22xco2ice']) # CO2lo + ICEwarm
    #targoffset1 = sumvecx / nhsies['2xco22xco2ice'] * 100
    sumvecy = (gms['prei2xco2iceb']-gms['2xco22xco2ice']) +\
              (gms['2xco2preiice']-gms['2xco22xco2ice']) # CO2lo + ICEwarm
    
    #print 'CO2lo + ICEwarm (ice,sat) ' + str(sumvecx/1e12), str(sumvecy)
    
    ax.plot((nhsies['2xco22xco2ice']+sumvecx)/1e12,
            gms['2xco22xco2ice']+sumvecy,
            marker='s',color='r',alpha=0.5)
    if fullann:
        ax.annotate('ICEwarm+CO2lo',((nhsies['2xco22xco2ice']+sumvecx)/1e12,gms['2xco22xco2ice']+sumvecy),
                    xytext=(-30,-15),textcoords='offset points')
    
    # add ICEcold+CO2hi point to fig
    sumvecx = (nhsies['prei2xco2iceb']-nhsies['preipreiice']) +\
              (nhsies['2xco2preiice']-nhsies['preipreiice']) # ICEcold + CO2hi
    sumvecy = (gms['prei2xco2iceb']-gms['preipreiice']) +\
              (gms['2xco2preiice']-gms['preipreiice'])
    #targoffset2 = sumvecx / nhsies['preipreiice'] * 100
    
    #print 'CO2hi + ICEcold (ice,sat) ' + str(sumvecx/1e12), str(sumvecy)
    
    ax.plot((nhsies['preipreiice']+sumvecx)/1e12,
            gms['preipreiice']+sumvecy,
            marker='s',color='b',alpha=0.5)
    if fullann:
        ax.annotate('ICEcold+CO2hi',((nhsies['preipreiice']+sumvecx)/1e12,gms['preipreiice']+sumvecy),
                    xytext=(-30,-15),textcoords='offset points')
            
    if dosia:
        print '---Doing SEA ICE AREA---'
        
    print sea
    print 'gms ' + str(gms)
    print 'nhsies' + str(np.array(nhsies.values())/1e12)
    
    # print GM SAT differences b/w ctl and perturbed:
    gmtsdiffs={}
    gmtsdiffs['icecold'] = gms['prei2xco2iceb'] - gms['preipreiice']
    gmtsdiffs['icewarm'] = gms['2xco22xco2ice'] - gms['2xco2preiice']
    gmtsdiffs['co2lo'] = gms['2xco22xco2ice'] - gms['prei2xco2iceb']
    gmtsdiffs['co2hi'] = gms['2xco2preiice'] - gms['preipreiice']
    gmtsdiffs['full'] = gms['2xco22xco2ice'] - gms['preipreiice']
    
    siediffs={}
    siediffs['icecold'] = nhsies['prei2xco2iceb'] - nhsies['preipreiice']
    siediffs['icewarm'] = nhsies['2xco22xco2ice'] - nhsies['2xco2preiice']
    siediffs['co2lo'] = nhsies['2xco22xco2ice'] - nhsies['prei2xco2iceb']
    siediffs['co2hi'] = nhsies['2xco2preiice'] - nhsies['preipreiice']
    siediffs['full'] = nhsies['2xco22xco2ice'] - nhsies['preipreiice']
    
    print 'prei2xco2iceb - preipreiice GM SAT : ' + str(gms['prei2xco2iceb'] - gms['preipreiice'])
    print '2xco22xco2ice - 2xco2preiice GM SAT : ' + str(gms['2xco22xco2ice'] - gms['2xco2preiice'])
    
    # Add Ice annotation
    ax.annotate('',xy=(nhsies['preipreiice']/1e12,gms['preipreiice']), 
                xytext=(nhsies['prei2xco2iceb']/1e12,gms['prei2xco2iceb']),
                arrowprops=dict(width=1,headlength=12,headwidth=10,
                                facecolor='b',edgecolor='b',shrink=.01,alpha=0.5))    #dict(arrowstyle='->',connectionstyle='arc3',
                                #facecolor='k',edgecolor='k')
    ax.annotate('ICEcold',xy=((nhsies['preipreiice']/1e12 + nhsies['prei2xco2iceb']/1e12)/2.,
                          (gms['preipreiice']+gms['prei2xco2iceb'])/2.), xytext=(0,-30),
                textcoords='offset points', fontsize=fsz) # former Ice
    # Add CO2 annotation
    ax.annotate('',xy=(nhsies['preipreiice']/1e12,gms['preipreiice']), 
                xytext=(nhsies['2xco2preiice']/1e12,gms['2xco2preiice']),
                arrowprops=dict(width=1,headlength=12,headwidth=10,
                                facecolor='b',edgecolor='b',shrink=.01,alpha=0.5))    
    ax.annotate('CO2hi',xy=((nhsies['preipreiice']/1e12 + nhsies['2xco2preiice']/1e12)/2.,
                          (gms['preipreiice']+gms['2xco2preiice'])/2.), xytext=(5,0),
                textcoords='offset points', fontsize=fsz) # former CO2
    # Add Ice2 annotation
    ax.annotate('',xy=(nhsies['2xco2preiice']/1e12,gms['2xco2preiice']), 
                xytext=(nhsies['2xco22xco2ice']/1e12,gms['2xco22xco2ice']),
                arrowprops=dict(width=1,headlength=12,headwidth=10,
                                facecolor='r',edgecolor='r',shrink=.01,alpha=0.5))    
    ax.annotate('ICEwarm',xy=((nhsies['2xco2preiice']/1e12 + nhsies['2xco22xco2ice']/1e12)/2.,
                          (gms['2xco2preiice']+gms['2xco22xco2ice'])/2.), xytext=(0,5),
                textcoords='offset points', fontsize=fsz) # former Ice2
    # Add CO22 annotation
    ax.annotate('',xy=(nhsies['prei2xco2iceb']/1e12,gms['prei2xco2iceb']), 
                xytext=(nhsies['2xco22xco2ice']/1e12,gms['2xco22xco2ice']),
                arrowprops=dict(width=1,headlength=12,headwidth=10,
                                facecolor='r',edgecolor='r',shrink=.01,alpha=0.5))    
    ax.annotate('CO2lo',xy=((nhsies['prei2xco2iceb']/1e12 + nhsies['2xco22xco2ice']/1e12)/2.,
                          (gms['prei2xco2iceb']+gms['2xco22xco2ice'])/2.), xytext=(-65,0),
                textcoords='offset points', fontsize=fsz) # former CO22
    # Add FULL annotation
    ax.annotate('',xy=(nhsies['preipreiice']/1e12,gms['preipreiice']), 
                xytext=(nhsies['2xco22xco2ice']/1e12,gms['2xco22xco2ice']),
                arrowprops=dict(width=1,headlength=12,headwidth=10,
                                facecolor='purple',edgecolor='purple',shrink=.01,alpha=0.5))    
    ax.annotate('Full',xy=((nhsies['preipreiice']/1e12 + nhsies['2xco22xco2ice']/1e12)/2.,
                          (gms['preipreiice']+gms['2xco22xco2ice'])/2.), xytext=(5,5),
                textcoords='offset points', fontsize=fsz)
    
    # Add x marker for values interpolated to 'target ice'
    icem = (gms['preipreiice']-gms['prei2xco2iceb']) / (nhsies['preipreiice']/1e12-nhsies['prei2xco2iceb']/1e12)
    iceb = gms['prei2xco2iceb'] - (icem*nhsies['prei2xco2iceb']/1e12)
    targx = nhsies['2xco22xco2ice']/1e12
    targy = icem*targx + iceb
    #print icem, iceb
    #print targx,targy

    ax.plot(targx,targy,color='0.5',marker='x')
    ax.axvline(x=targx, color='0.5',linewidth=.5)
    
    ice2m = (gms['2xco22xco2ice']-gms['2xco2preiice']) / (nhsies['2xco22xco2ice']/1e12-nhsies['2xco2preiice']/1e12)
    ice2b = gms['2xco2preiice'] - (ice2m*nhsies['2xco2preiice']/1e12)
    targ2x = nhsies['preipreiice']/1e12
    targ2y = ice2m*targ2x + ice2b
    #print ice2m, ice2b
    #print targ2x,targ2y
    
    ax.plot(targ2x,targ2y,color='0.5',marker='x')
    ax.axvline(x=targ2x, color='0.5',linewidth=.5)
    
    # Add % offset information
    #   original just two values at the target ice markers
    mm1 = toffsets['icewarm']; pctmm1 = pcttoffsets['icewarm']
    mm2 = toffsets['icecold']; pctmm2 = pcttoffsets['icecold']
    if fullann:
        ax.annotate('%2.2f'%mm2+' (%2.1f'%pctmm2+'%)', xy=(targx,targy),
                    xytext=(-55,-15),textcoords='offset points')
        ax.annotate('%2.2f'%mm1+' (%2.1f'%pctmm1+'%)', xy=(targ2x,targ2y),
                    xytext=(10,-5),textcoords='offset points')
    
    #   try 4 values for each of the 'differences' e.g. ICEcold, etc
    # @@@@

    if dosia:
        ax.set_xlabel('Arctic Sea Ice Area (10$^6$ km$^2$)',fontsize=fsz)
        plab='sia'
    else:
        ax.set_xlabel('NH SIE (10$^6$ km$^2$)',fontsize=fsz)
        plab='sie'
    
    ax.set_ylabel('Global Mean Surface Air Temperature ($^\circ$C)',fontsize=fsz)
    ax.set_title(sea,fontsize=fsz)

    if sea=='ANN':
        ax.set_ylim((13,18))
        if dosia:
            ax.set_xlim((4,11))
        else:
            ax.set_xlim((10,15))
    elif sea=='JJA':
        ax.set_ylim((14.5,19.5))
        if dosia:
            ax.set_xlim((1.5,7))
        else:
            ax.set_xlim((6,11))
    elif sea=='DJF':
        ax.set_ylim((11,16))
        if dosia:
            ax.set_xlim((7,13.5))
        else:
            ax.set_xlim((10,16))
    elif sea=='SON':
        ax.set_ylim((13,18))
        if dosia:
            ax.set_xlim((-1,8))
        else:
            ax.set_xlim((-1,12))
    elif sea=='MAM':
        ax.set_ylim((13,17.5))
        if dosia:
            ax.set_xlim((9.5,13.5))
        else:
            pass

    if coconc:
        ax.set_ylabel('x CO$_2$ concentration',fontsize=16)
        cogmstr='coconc'
        ax.set_ylim((0.8,2.2))
    else:
        cogmstr='gmsat'
        
    ax.set_xticklabels(ax.get_xticks(), fontsize=fsz-2)
    ax.set_yticklabels(ax.get_yticks(), fontsize=fsz-2)
    
    if printtofile:
        if fullann: suffsuff=''
        else: suffsuff='noannot'
            
        fig.savefig('relaxruns_schematic_'+plab+'_' +cogmstr +'_' + sea + '_' + suff + suffsuff + '.pdf')
        fig.savefig('relaxruns_schematic_'+plab+'_' +cogmstr +'_' + sea + '_' + suff + suffsuff + '.jpg',dpi=300)

    return gmtsdiffs,siediffs


def compute_targoffset(ncicedt,dosia=True,dospatial=False,verb=True):
    """ dospatial supercedes dosia

    """
    sies={}
    for case in ncicedt.keys():
        if verb:
            print case
            
        icetmp = np.mean(ncicedt[case],axis=0)#/100.

            
        if dospatial: # this is for spatial scaling
            verb=False
            sies[case] = icetmp/100.
        else:
            if dosia:
                print 'icetmp.shape ' + str(icetmp.shape) # @@@@@@@@@@@@@@@22
                if last=='last200':
                    sies[case] = icetmp #@@ test
                else:
                    sies[case],_ = cutl.calc_totseaicearea(icetmp/100.,lat,lon,model=None,isarea=False)
            else:
                if last=='last200':
                    print 'SIE for last200 not implemented! @@@'
                else:
                    sies[case],_ = cutl.calc_seaiceextent(icetmp/100.,lat,lon,model=None)
            
    fullchange = sies['2xco22xco2ice'] - sies['preipreiice']
    
    pctoffsets={}; offsets={}
    offsets['icecold'] = (sies['2xco22xco2ice'] - sies['prei2xco2iceb'])/1e12
    pctoffsets['icecold'] = (sies['2xco22xco2ice'] - sies['prei2xco2iceb']) / fullchange * 100
    # for co2lo: how much deviation is the sea ice from climo (here an increase 
    #            --> except still means co2 isolation contaminated by sea ice loss)
    offsets['co2lo'] = (sies['2xco22xco2ice']-sies['prei2xco2iceb'])/1e12
    pctoffsets['co2lo'] = (sies['2xco22xco2ice']-sies['prei2xco2iceb'])/sies['2xco22xco2ice']*100
    #pctoffsets['co2lo'] = pctoffsets['icecold']
    offsets['icewarm'] = (sies['2xco2preiice'] - sies['preipreiice'])/1e12
    pctoffsets['icewarm'] = (sies['2xco2preiice'] - sies['preipreiice']) / fullchange * 100
    # for co2hi: how much deviation is the sea ice from climo (here a decrease)
    offsets['co2hi'] = (sies['2xco2preiice'] - sies['preipreiice'])/1e12
    pctoffsets['co2hi'] = (sies['2xco2preiice'] - sies['preipreiice'])/sies['preipreiice']*100
    #pctoffsets['co2hi'] = pctoffsets['icewarm']
    
    if verb:
        print '\n   CO2hi: CO2 isolation has a effect from ice melt (ice not grown enough): ' +\
            str(pctoffsets['co2hi']) + '% of climo'
        print '   CO2lo: CO2 isolation has a effect from ice melt (ice not melted enough): ' +\
            str(pctoffsets['co2lo']) + '% of climo'
        print '   ICEwarm: sea ice melt effect is underestimated (ice not grown enough): ' +\
            str(pctoffsets['icewarm']) + '% of full change'
        print '   ICEcold: sea ice melt effect is underestimated (ice not melted enough): ' +\
            str(pctoffsets['icecold']) + '% of full change\n'
         
    return offsets,pctoffsets



def calc_nc_agreement_map(ncdt, axs, lin='one',conv=1, subtime=None, magtype='abs'):
    """         
        lin = 'one' is ice and co2 (ICEcold, CO2hi)
            = 'two' is ice2 and co22 (ICEwarm, CO2lo)
            = 'sub' subtract lin 'one' from lin 'two' (shows how warm clim/low ice different from cold/high)
                                       
        magtype: if 'abs' then magnitude is np.abs(ice) + np.abs(co2)
                 else it is ice+co2
                 
        returns: mag=|ice|+|co2|, and agreement (mag with sign)
    """
       
    pico = ncdt['prei2xco2iceb']
    pipi = ncdt['preipreiice']
    copi = ncdt['2xco2preiice']
    coco = ncdt['2xco22xco2ice']
    
    if lin=='one':
        _,icepv = cutl.ttest_ind(pico, pipi,axis=0,effdof=False)        
        ice = (np.mean(pico,axis=0) - np.mean(pipi,axis=0))*conv
        #icep = ncfldzmdt['pi2xco2ipulse'] - ncfldzmdt['preipreiice']    
        _,co2pv = cutl.ttest_ind(copi, pipi,axis=0,effdof=False)   
        co2 = (np.mean(copi,axis=0) - np.mean(pipi,axis=0))*conv
        _,combopv = cutl.ttest_ind(pico-pipi+(copi-pipi),
                                   coco-pipi)
        ctl = np.mean(pipi,axis=0)*conv
        
        suff=''
        suff1='cold'; suff2='hi'
    elif lin=='two':
        _,icepv = cutl.ttest_ind(coco, copi,axis=0,effdof=False)
        ice = (np.mean(coco,axis=0) - np.mean(copi,axis=0))*conv
        _,co2pv = cutl.ttest_ind(coco, pico,axis=0,effdof=False)
        co2 = (np.mean(coco,axis=0) - np.mean(pico,axis=0))*conv
        _,combopv = cutl.ttest_ind(coco-copi+(coco-pico),
                                   coco-pipi)
        ctl = np.mean(coco,axis=0)*conv
        suff='2'
        suff1='warm'; suff2='lo'
    elif lin=='sub':
        #two - one: (coco-copi) - (pico-pipi)
        _,icepv = cutl.ttest_ind(coco-copi, pico-pipi,axis=0,effdof=False)
        ice = (np.mean(coco-copi,axis=0) - np.mean(pico-pipi,axis=0))*conv
        #two - one: (coco-pico) - (copi-pipi)
        _,co2pv = cutl.ttest_ind(coco-pico, copi-pipi,axis=0,effdof=False)
        co2 = (np.mean(coco-pico,axis=0) - np.mean(copi-pipi,axis=0))*conv
        suff='2-1'
        suff1='warm-cold'; suff2='lo-hi'
        
        
    _,fullpv = cutl.ttest_ind(coco, pipi,axis=0,effdof=False)
    full = (np.mean(coco,axis=0) - np.mean(pipi,axis=0))*conv

    
    if magtype == 'abs':
        mag = np.abs(ice) + np.abs(co2)
    else:
        mag = ice + co2 # could also try: ice + np.abs(co2)
    
    import copy as copy
    
    magsign=copy.copy(mag)
    
    magsign[np.sign(ice) != np.sign(co2)] = -1*mag[np.sign(ice) != np.sign(co2)] 
    
    return mag,magsign


def plot_nc_agreement_map(ncdt, ax, lin='one',magtype='abs',cmin='',cmax='',cmind='',cmaxd='',cmin2='',cmax2='',
                           conv=1, ptype='nh',suppcb=False,cmap='blue2red_w20',subtime=None,
                           screen=False,vert=False,levlim=None,suppttl=False,addsig=False,sigtype='cont',
                           latlim=None,vertptype=None, fsz=None,plab=None):
    """ 
        cmin/cmax:  for ice
        
        lin = 'one' is ice and co2 (ICEcold, CO2hi)
            = 'two' is ice2 and co22 (ICEwarm, CO2lo)
            = 'sub' subtract lin 'one' from lin 'two' (shows how warm clim/low ice different from cold/high)
                    
        fsz: fontsize for title
                    
        returns: plot handle (ph) 
    """
    
    
    mag,magsign = calc_nc_agreement_map(ncdt, axs, lin=lin,magtype=magtype,
                           conv=1, subtime=None)
    
    
    # @@ subtime not implemented yet. assume average over full time 
    #tmplen = ncdt['2xco2preiice'].shape[0]-1
    fmt='%2.1f' # clabel format
    #if ctlconts!=None:
    #    if (ctlconts < 1.).any() and (ctlconts > -1.).any():
    #    #if np.logical_and((ctlconts<1).any(),(ctlconts>-1).any()):
    #        print ctlconts
    #        fmt='%2.1f'
    
    if lin=='one':
        suff=''
        suff1='cold'; suff2='hi'
    elif lin=='two':
        suff='2'
        suff1='warm'; suff2='lo'
    
        # pparams for ice!
    # pparams2 for co2![sum|full]
    if vert:
        pparams = {'cmin': cmin, 'cmax': cmax,
                   'suppcb': suppcb, 'screen': screen, 'cmap':cmap,'levlim':levlim,
                   'latlim':latlim,'ptype':vertptype}
    else:
        pparams = {'ptype': ptype, 'cmin': cmin, 'cmax': cmax,
                   'suppcb': suppcb, 'cmap':cmap}
    

    if suppttl: ttl=''
    else: ttl='ICE'+suff1 +' & CO2'+suff2
    if vert:
        ph = cplt.vert_plot(magsign,lev,lat,axis=ax,title=ttl,**pparams)
        #if addsig:
        #    cplt.addtsig(ax,icepv,lat,lev/100.,sigtype=sigtype,reverse=sigreverse)
        #if addclimcont:
        #    cplt.add_contoursvert(ax,ctl,lat,lev,verb=True,clab=True,conts=ctlconts,fmt=fmt)
    else:
        bm,ph = cplt.kemmap(magsign,lat,lon,axis=ax,title=ttl,**pparams)
        #if addsig:
        #    cplt.addtsigm(bm,icepv,lat,lon,sigtype=sigtype,reverse=sigreverse)
    if fsz!=None:
        ax.set_title(ttl,fontsize=fsz)
        if plab!=None:
            ax.annotate(plab, xy=(0.02,1.02),xycoords='axes fraction',
                       fontsize=fsz-1,fontweight='bold')

    
    
    return ph
