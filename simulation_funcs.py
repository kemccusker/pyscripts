""" simulation_funcs.py
        9/9/2014: Goal is to use this to house functions to
                     prepare simulation data, and then plot data (at a higher
                     level than cccmaplots.py)
                     e.g. simulations and figures more specific to
                          sea ice paper.
                  In future, this need not be limited to canam4, so
                  I have left the filename general, and will update
                  the innards to be able to handle other model data.
                  (right now model='CanAM4' etc is hard-coded. 9/9/14)
                  
"""

import numpy as np
import numpy.ma as ma
import scipy as sp
import scipy.stats
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import cccmacmaps as ccm
import cccmaNC as cnc

cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)
ccm = reload(ccm)
cnc = reload(cnc)


def plot_seasonal_maps(fielddict,coords,sims,pparams,vert=False,loctimesel=None,info=None,printtofile=False,seas=None):
    """ plot_seasonal_maps(fielddict,coords,sims,pparams,vert=False,loctimesel=None,info=None,printtofile=False):
    
              info should be a dict of catch-all keywords/info
    """
    
    field=fielddict['field']
    ncfield=fielddict['ncfield']
    fieldstr=fielddict['fieldstr']
    conv=fielddict['conv']

    if seas != None:
        seasons=seas
    else: # standard seasons
        seasons=('SON','DJF','MAM','JJA')
        
    model=info['model'] #'CanAM4' # @@ move out of function
    sigtype=info['sigtype'] # specify sig marking type
    sigoff=info['sigoff']  # turn off significance marking?
    pct=info['pct']   # percentage change?
    nonstandardlev=fielddict['nonstandardlev']


    bp=con.get_basepath()
    basepath=bp['basepath'] + model + '/'; subdir=bp['subdir'] # @@ move out of function?


    
    lat=coords['lat']
    nlat=len(lat)
    if vert==True:
        lev=coords['lev']
        nlev=len(lev)
        theshape=(len(seasons),nlev,nlat)      
    else:
        lon=coords['lon']
        nlon=len(lon)
        theshape=(len(seasons),nlat,nlon)


    tstat = np.zeros(theshape)
    pval = np.zeros(theshape)
    fldcallseas = np.zeros(theshape)
    fldpallseas = np.zeros(theshape)

    cmap=pparams['cmap']
    cmin=pparams['cmin']; cmax=pparams['cmax']
    cmlen=float( plt.cm.get_cmap(cmap).N)
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    fig6,ax6 = plt.subplots(len(seasons),len(sims)) # 1 row for e/ of 5 ens members, plus mean, plus meanBC
    fig6.set_size_inches(12,8)  
    fig6.subplots_adjust(hspace=.15,wspace=.05)

    lastcol=len(sims)-1
    
    for colidx,sim in enumerate(sims): # traverse cols

        rowidx=0
        simpair = con.get_simpair(sim)
        simctl = simpair['ctl']['fullname']
        timstr = simpair['ctl']['timestr']
        simpt = simpair['pert']['fullname']
        timstrp = simpair['pert']['timestr']
        frootc = basepath + simctl + subdir + simctl + '_' + field + '_'
        frootp = basepath + simpt + subdir + simpt + '_' + field + '_'
        rowl=sim

        if loctimesel !=None: # e.g. for when want to look at half of a run...
            timesel = loctimesel
        else:
            timesel = simpair['ctl']['timesel']
        
        fnamec = frootc + timstr + '_ts.nc'
        fnamep = frootp + timstrp + '_ts.nc'
        # @@@ I *think* I don't need this anymore since processing the 3D files more
        # @@  although getting a nonstandard lev is not supported
        if nonstandardlev:
            print 'not yet supported for seasonal'
        ## if threed==0:
        ##     fnamec = frootc + timstr + '_ts.nc'
        ##     fnamep = frootp + timstrp + '_ts.nc'
        ## else:
        ##     print '@@ fix to use the level NC files'
        ##     fnamec = frootc + '001-061_ts.nc'
        ##     fnamec2 = frootc + '062-121_ts.nc'
        ##     fnamep = frootp + '001-061_ts.nc'
        ##     fnamep2 = frootp + '062-121_ts.nc'
                
        for sea in seasons: # traverse rows (or use map_allseas() ?? @@)
            #ax = ax6[ridx][cidx]
            ax = ax6[rowidx][colidx] # swapped row and col index positions in subplot
            
            if field=='turb':
                field='hfl'; fieldb='hfs'
                fldcsea = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                               seas=sea)*conv + cnc.getNCvar(fnamecb,fieldb.upper(),
                                               timesel=timesel,seas=sea)*conv
                fldpsea = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                               seas=seas)*conv + cnc.getNCvar(fnamepb,fieldb.upper(),
                                               timesel=timesel,seas=sea)*conv 
                field='turb'
            elif field=='net':
                print '@@ not implemented for seasonal maps'
            else:
                    
                fldcsea = cnc.getNCvar(fnamec,ncfield,timesel=timesel,
                                               seas=sea)*conv
                fldpsea = cnc.getNCvar(fnamep,ncfield,timesel=timesel,
                                               seas=sea)*conv
                ## if threed==0:
                ##     fldcsea = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                ##                                    seas=sea)*conv
                ##     fldpsea = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                ##                                    seas=sea)*conv
                ## else:
                ##     fldcsea = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,061-12-31',levsel=level,
                ##                                    seas=sea)*conv,
                ##                         cnc.getNCvar(fnamec2,ncfield,levsel=level,seas=sea)*conv,axis=0)
                ##     fldpsea = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,061-12-31',levsel=level,
                ##                                    seas=sea)*conv,
                ##                         cnc.getNCvar(fnamep2,ncfield,levsel=level,seas=sea)*conv,axis=0)
                    
            tstat[rowidx,:,:],pval[rowidx,:,:] = sp.stats.ttest_ind(fldpsea,fldcsea,axis=0)
            fldcallseas[rowidx,:,:] = np.mean(fldcsea,axis=0)
            fldpallseas[rowidx,:,:] = np.mean(fldpsea,axis=0)

            if pct:
                plotfld = (fldpallseas[rowidx,:,:]-fldcallseas[rowidx,:,:]) / fldcallseas[rowidx,:,:] *100
            else:
                plotfld = fldpallseas[rowidx,:,:] - fldcallseas[rowidx,:,:]

            pparams['axis']=ax
            
            if vert: # zonal mean with height
                pparams['suppcb'] = True
                pparams['levlim'] = levlim
                pparams['screen'] = screen
                pparams['addcontlines'] = True
                if colidx!=0: # if not the first column, suppress y labels
                    pparams['suppylab'] = True
                
                    
                pc = cplt.vert_plot(plotfld,lev,lat,**pparams)
                if sigoff==0:
                     cplt.addtsig(ax,pval[rowidx,...],lat,lev/100.,type=sigtype) # @@ dims?

                if colidx==lastcol:
                    # put season label on right side.
                    ax.yaxis.set_label_position("right")
                    ax.set_ylabel(sea)
                    
            else: # maps
                pparams['suppcb'] = 1
                bm,pc = cplt.kemmap(plotfld,lat,lon,**pparams)#@@

                if sigoff==0:
                    cplt.addtsigm(bm,pval[rowidx,:,:],lat,lon,type=sigtype)

                if colidx==0: # when col index is 0, set season
                    ax.set_ylabel(sea)

            if rowidx==0: # when row index is 0, set simulation
                ax.set_title(rowl)     

            rowidx = rowidx+1

    cbar_ax = fig6.add_axes([.91,.25, .02,.5])
    fig6.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    plt.suptitle(fieldstr)

    ## if printtofile:
    ##     if sigoff==0:
    ##         sigstr='sig' + sigtype
    ##     else:
    ##         sigstr=''

    ##     if latlim!= None:
    ##         latstr=str(latlim)
    ##     else:
    ##         latstr=''

    ##     if seasvert:
    ##         if screen:
    ##             style = 'screen'
    ##         else:
    ##             style = str(latlim) + 'N' + str(levlim) + 'hPa'
                
    ##         if pct:
    ##             fig6.savefig(fieldstr + 'pctdiff' + sigstr + '_enssubplot' + obsstr + ctstr +
    ##                          '_seas_' + style + '2.' + suff)
    ##         else:
    ##             fig6.savefig(fieldstr + 'diff' + sigstr + '_enssubplot' + obsstr + ctstr +
    ##                          '_seas_' + style + '2.' + suff)
    ##     else: # maps
    ##         if pct: # version 2 has new season order, new filename/key org
    ##             fig6.savefig(fieldstr + 'pctdiff' + sigstr + '_enssubplot' + obsstr + ctstr +
    ##                          '_seas_nh' + latstr + '2.' + suff)
    ##         else:
    ##             fig6.savefig(fieldstr + 'diff' + sigstr + '_enssubplot' + obsstr + ctstr + '_seas_nh'
    ##                          + latstr + '2.' + suff)


def calc_seasonal_cycle(fielddict,coords,sims,withlat=False,loctimesel=None,info=None):
    """ calc_seasonal_cycle(fielddict,coords,sims,loctimesel=None,info=None)

              returns a data blob:
                   blob['ctl'] = fldcdict
                   blob['pert'] = fldpdict
                   blob['tstat'] = tstatdict
                   blob['pval'] = pvaldict
                   blob['ctlstd'] = fldcstddict
                   blob['pertstd'] = fldpstddict
                   blob['diff'] = flddiffdict
                   blob['mask'] = flddmaskdict
    """
    
    # @@ switch to handling any region

    field=fielddict['field']
    ncfield=fielddict['ncfield']
    fieldstr=fielddict['fieldstr']
    conv=fielddict['conv']

    # ======  from canam4sims_analens......

    """    # get data for either zonal mean or sea cycle figures
    if plotzonmean==1 or pattcorrwithtime==1 or plotregmean==1:
        # seasons is defined above
        corrlim = 45
    elif plotseacyc==1:
        
        if field in (fluxes,'fsg','turb','net'):
            seacyclatlim=40
        # else leave seacyclatlim as set at top of script"""
    seacyc = con.get_mon()
    seacyclatlim=info['seacyclatlim'] # should be 40N for (fluxes,'fsg','turb','net')

    if field=='sia':
        sia=True
        field = 'sicn' # while getting the data...
        
    tstatdict = dict.fromkeys(sims,{}); pvaldict = dict.fromkeys(sims,{})
    fldcdict = dict.fromkeys(sims,{}); fldpdict = dict.fromkeys(sims,{})
    fldcstddict = dict.fromkeys(sims,{}); fldpstddict = dict.fromkeys(sims,{})
    flddiffdict = dict.fromkeys(sims,{}); flddmaskdict = dict.fromkeys(sims,{})
    fldpcorrdict = dict.fromkeys(sims,{});
    cidict = dict.fromkeys(sims,{})

    for ridx,sim in enumerate(sims):    
        seatstatdict=dict.fromkeys(seacyc); seapvaldict=dict.fromkeys(seacyc)
        seafldcdict=dict.fromkeys(seacyc); seafldpdict=dict.fromkeys(seacyc)
        seafldcstddict=dict.fromkeys(seacyc); seafldpstddict=dict.fromkeys(seacyc)
        seadiffdict=dict.fromkeys(seacyc); seadmaskdict=dict.fromkeys(seacyc)
        seapcorrdict=dict.fromkeys(seacyc)
        seacidict=dict.fromkeys(seacyc)
        # switch these to just an array? 0:12. @@


        simpair = con.get_simpair(sim)
        simctl = simpair['ctl']['fullname']
        timstr = simpair['ctl']['timestr']
        simpt = simpair['pert']['fullname']
        timstrp = simpair['pert']['timestr']
        frootc = basepath + simctl + subdir + simctl + '_' 
        frootp = basepath + simpt + subdir + simpt + '_' 

        if loctimesel !=None: # e.g. for when want to look at half of a run...
            timesel = loctimesel
        else:
            timesel = simpair['ctl']['timesel']

        """if sim=='kemhad' or sim=='kemnsidc': 
            frootc = basepath + sim + 'ctl' + subdir + sim + 'ctl' + '_' 
            frootp = basepath + sim + 'pert' + subdir + sim + 'pert' + '_' 
        elif sim in ('kem1pert1b','kem1pert3','kem1rcp85a'):
            frootc = basepath + 'kemctl1' + subdir + 'kemctl1' + '_' 
            frootp = basepath + sim + subdir + sim + '_'
        else:
            frootc =  basepath + bcasename + sim + subdir + bcasename + sim + '_'
            frootp = basepath + bcasenamep + sim + subdir + bcasenamep + sim + '_' """

        if threed and nonstandardlev:
            fnamec = frootc + field + '_' + '001-061_ts.nc'
            fnamec2 = frootc + field + '_' + '062-121_ts.nc'
            fnamep = frootp + field + '_' + '001-061_ts.nc'
            fnamep2 = frootp + field + '_' + '062-121_ts.nc'
        else:
            fnamec = frootc + field + '_' + timstr + '_ts.nc'
            fnamep = frootp + field + '_' + timstrp + '_ts.nc'
        
        for sii,sea in enumerate(seacyc):

            #if plotzonmean==1 or pattcorrwithtime==1 or plotregmean==1:
            #    ncparams = {'seas': sea}
            #elif plotseacyc==1:
            ncparams = {'monsel': sii+1}

            # Now get the data
            if field in ('turb','net'):
                #print 'not implemented @@'
                #print 'field is ' + field + '. getting hfl, hfs'
                fielda='hfl'; fieldb='hfs'
                fnamec = frootc + fielda + '_' + timstr + '_ts.nc'
                fnamep = frootp + fielda + '_' + timstrp + '_ts.nc'
                fnamecb = frootc + fieldb + '_' + timstr + '_ts.nc'
                fnamepb = frootp + fieldb + '_' + timstrp + '_ts.nc'

                fldc = cnc.getNCvar(fnamec,fielda.upper(),timesel=timesel,
                                             **ncparams)*conv + cnc.getNCvar(fnamecb,fieldb.upper(),
                                             timesel=timesel,**ncparams)*conv
                fldp = cnc.getNCvar(fnamep,fielda.upper(),timesel=timesel,
                                             **ncparams)*conv + cnc.getNCvar(fnamepb,fieldb.upper(),
                                             timesel=timesel,**ncparams)*conv
                if field=='net':
                    #print 'getting flg for net'
                    fieldb='flg'
                    conv=-1
                    fnamecb = frootc + fieldb + '_' + timstr + '_ts.nc'
                    fnamepb = frootp + fieldb + '_' + timstrp + '_ts.nc'
                    fldc = fldc + cnc.getNCvar(fnamecb,fieldb.upper(),
                                                 timesel=timesel,**ncparams)*conv
                    fldp = fldp + cnc.getNCvar(fnamepb,fieldb.upper(),
                                                 timesel=timesel,**ncparams)*conv
                    field='net' # return field name to what it should be
                    conv=1
                else:
                    field='turb' # return field name to what it should be
            # end if turb or net  
            else:
                """if threed==0:
                    fldczm = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                                   **ncparams)*conv
                    fldpzm = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                                   **ncparams)*conv

                else:
                    #print '@@ fix to use level NC files'"""
                    
                if threed and nonstandardlev:
                    ncparams['levsel'] = level
                    fldc = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                      cnc.getNCvar(fnamec2,ncfield,**ncparams)*conv,
                                      axis=0)
                    fldp = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                      cnc.getNCvar(fnamep2,ncfield,**ncparams)*conv,
                                      axis=0)
                else:
                    fldc = cnc.getNCvar(fnamec,ncfield,timesel=timesel,
                                        **ncparams)*conv
                    fldp = cnc.getNCvar(fnamep,ncfield,timesel=timesel,
                                        **ncparams)*conv
                        
                if sia==True:
                    fldc = cutl.calc_seaicearea(fldc,lat,lon)
                    fldp = cutl.calc_seaicearea(fldp,lat,lon)
            # end getting data initially
        
            # now do processing (regional means, tot SIA, etc)
            if field in (fluxes,'fsg','turb','net'):
                # mask out regions that are not ice in the control (as P&M 2014 JClim)
                sicnc = cnc.getNCvar(frootc + 'sicn_' + timstr + '_ts.nc','SICN',timesel=timesel,**ncparams)
                
                fldc = ma.masked_where(sicnc<.10,fldc)
                fldp = ma.masked_where(sicnc<.10,fldp)
                
            """if plotzonmean==1:
                fldczm = np.mean(fldczm[...,:-1],axis=2) # actually take zonal mean now
                fldpzm = np.mean(fldpzm[...,:-1],axis=2)
                
            elif pattcorrwithtime==1:
                # loop through each year
                # calc pattern corr either yearly or integrated
                years=np.arange(0,fldczm.shape[0])
                fldctm = np.mean(fldczm[:,lat>corrlim,...],axis=0)
                pcorr = np.zeros(len(years))
                for yr in years:
                    areas = cutl.calc_cellareas(lat,lon)
                    areas = areas[lat>corrlim,:]
                    weights = areas / np.sum(np.sum(areas,axis=1),axis=0)
                    if pattcorryr:
                        # yearly anomaly pattern corr w/ the time mean pattern
                        tmp = fldpzm[yr,lat>corrlim,...]-fldctm
                    else:
                        tmp = np.mean(fldpzm[:yr,lat>corrlim,...],axis=0)-fldctm # integrated anomaly pattern
                    tmpmean = np.mean(fldpzm[:,lat>corrlim,...],axis=0) - fldctm # end pattern to compare against
                    pcorr[yr] = cutl.pattcorr(tmp.flatten()*weights.flatten(),tmpmean.flatten()*weights.flatten())

                seapcorrdict[sea] = pcorr"""

            #elif plotseacyc==1:
            if withlat:
                # leave the latitude dimension intact
                # dims are time x lat x lon to start
                # take zonal mean
                fldc = np.mean(fldc,axis=2) # take zonal mean
                fldp = np.mean(fldp,axis=2)
            else:
                if sia:
                    #calc total area instead of average
                    fldc,sh = cutl.calc_totseaicearea(fldc,lat,lon) #np.sum(np.sum(fldczm[:,lat>0,:],axis=2),axis=1)
                    fldp,sh = cutl.calc_totseaicearea(fldp,lat,lon) #np.sum(np.sum(fldpzm[:,lat>0,:],axis=2),axis=1)
                else:
                    # consider masking out land for sfc fluxes...?
                    # @@@ switch to regional mean? which includes polar means....
                    fldc = cutl.polar_mean_areawgted3d(fldc,lat,lon,latlim=seacyclatlim)
                    fldp = cutl.polar_mean_areawgted3d(fldp,lat,lon,latlim=seacyclatlim)

            """elif plotregmean==1:
                
                #limsdict = con.get_regionlims(region)
                fldczm = cutl.calc_regmean(fldczm,lat,lon,region)#limsdict)
                fldpzm = cutl.calc_regmean(fldpzm,lat,lon,region)#limsdict) """

            # now save the processed fields
            seafldcstddict[sea] = np.std(fldc,axis=0)
            seafldpstddict[sea] = np.std(fldp,axis=0)
            ttmp,pvtmp = sp.stats.ttest_ind(fldp,fldc,axis=0)
            """if plotregmean==1:
                # calculate confidence interval on the regional mean
                # double-check the scale setting
                ci = sp.stats.t.interval(1-siglevel,len(fldpzm)-1,loc=np.mean(fldpzm,axis=0)-np.mean(fldczm,axis=0),
                                         scale=np.std(fldpzm,axis=0)/np.sqrt(len(fldpzm)))
                seacidict[sea] = ci
                #print ci # @@@ """
                
            seatstatdict[sea] = ttmp
            seapvaldict[sea] = pvtmp
            seafldcdict[sea] =  np.mean(fldc,axis=0) # time mean
            seafldpdict[sea] =  np.mean(fldp,axis=0)
            seadiffdict[sea] = np.mean(fldp,axis=0)- np.mean(fldc,axis=0)
            seadmaskdict[sea] = ma.masked_where(pvtmp>siglevel,seadiffdict[sea])

            # end loop through seacyc

        fldcstddict[sim] = seafldcstddict
        fldpstddict[sim] = seafldpstddict
        #if plotregmean==1:
        #    cidict[sim] = seacidict
        tstatdict[sim] = seatstatdict
        pvaldict[sim] = seapvaldict
        fldcdict[sim] = seafldcdict
        fldpdict[sim] = seafldpdict
        flddiffdict[sim] = seadiffdict
        flddmaskdict[sim] = seadmaskdict
        #if pattcorrwithtime==1:
        #    fldpcorrdict[sim] = seapcorrdict

        # end loop through simulations
    if sia:
        field = 'sia' # put back after getting the data

    blob = {}
    blob['ctl'] = fldcdict
    blob['pert'] = fldpdict
    blob['tstat'] = tstatdict
    blob['pval'] = pvaldict
    blob['ctlstd'] = fldcstddict
    blob['pertstd'] = fldpstddict
    blob['diff'] = flddiffdict
    blob['mask'] = flddmaskdict

    return blob # datablob

    # ======
