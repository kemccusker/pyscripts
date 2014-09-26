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
import matplotlib.font_manager as fm
import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import cccmacmaps as ccm
import cccmaNC as cnc
import pandas as pd

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
    savestr=info['savestr']

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

    latlim=pparams['latlim']
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
                #pparams['levlim'] = levlim
                #pparams['screen'] = screen
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

    if printtofile:
    ##     if sigoff==0:
    ##         sigstr='sig' + sigtype
    ##     else:
    ##         sigstr=''
        sigstr='sigcont' #@@ hard code
        suff='pdf'

        if latlim!= None:
            latstr=str(latlim)
        else:
            latstr=''

        if vert:
            if screen:
                style = 'screen'
            else:
                style = str(latlim) + 'N' + str(levlim) + 'hPa'
                
            if pct:
                fig6.savefig(fieldstr + 'pctdiff' + sigstr + '_enssubplot' + savestr +
                             '_seas_' + style + '2.' + suff)
            else:
                fig6.savefig(fieldstr + 'diff' + sigstr + '_enssubplot' + savestr +
                             '_seas_' + style + '2.' + suff)
        else: # maps
            if pct: # version 2 has new season order, new filename/key org
                fig6.savefig(fieldstr + 'pctdiff' + sigstr + '_enssubplot' + savestr +
                             '_seas_nh' + latstr + '2.' + suff)
            else:
                fig6.savefig(fieldstr + 'diff' + sigstr + '_enssubplot' + savestr + '_seas_nh'
                             + latstr + '2.' + suff)

def calc_seasons(fielddict,coords,sims,loctimesel=None,info=None,siglevel=0.05,
                 calctype=None):
    """ calc_seasons(fielddict,coords,sims,withlat=False,loctimesel=None,info=None,siglevel=0.05)

              info: a dict which specifies lots of things, but in particular, which region to average
              calctype: 'zonmean', 'regmean', 'pattcorrwithtime', 'pattcorrwithtimeyr', None
              
              returns a data blob:
                   blob['ctl'] = fldcdict
                   blob['pert'] = fldpdict
                   blob['tstat'] = tstatdict
                   blob['pval'] = pvaldict
                   blob['ci'] = cidict
                   blob['ctlstd'] = fldcstddict
                   blob['pertstd'] = fldpstddict
                   blob['diff'] = flddiffdict
                   blob['mask'] = flddmaskdict
                   blob['pcorr'] = fldpcorrdict # pattern corr with time only
    """

    seasons='SON','DJF','MAM','JJA'
    sia = False

    field=fielddict['field']
    ncfield=fielddict['ncfield']
    fieldstr=fielddict['fieldstr']
    conv=fielddict['conv']
    threed=fielddict['threed']
    isflux=fielddict['isflux']
    nonstandardlev=fielddict['nonstandardlev']

    model=info['model'] 
    pct=info['pct']   # percentage change?
    corrlim=info['corrlim'] # southern limit for pattern correlation
    
    bp=con.get_basepath()
    basepath=bp['basepath'] + model + '/'; subdir=bp['subdir'] # @@ move out of function?

    lat=coords['lat']
    nlat=len(lat)
    lon=coords['lon']
    nlon=len(lon)

    plotzonmean=False; plotregmean=False; pattcorrwithtime=False; pattcorryr=False

    if calctype!=None and calctype=='zonmean':
        plotzonmean=True
    elif calctype!=None and calctype=='regmean':
        plotregmean=True
        region=info['region']
    elif calctype!=None and calctype=='pattcorrwithtime':
         # note that this is pattern correlating a run's
         # pattern in time (cumulative avg) w/ its final pattern
        pattcorrwithtime=True
    elif  calctype!=None and calctype=='pattcorrwithtimeyr':
         # note that this is pattern correlating a run's
         # pattern each year in time w/ its final pattern, then sorted
        pattcorrwithtime=True
        pattcorryr=True
        
    # note that pattern corr with time will be a 2D processed field -->
    #    for each season, fld.shape = ntime
    # zonal mean will be a 2D processed field -->
    #    for each season, fld.shape = nlat
    # regional mean will be a 1D processed field -->
    #    for each season, fld.shape = 1 (regional mean)
    

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
        seatstatdict=dict.fromkeys(seasons); seapvaldict=dict.fromkeys(seasons)
        seafldcdict=dict.fromkeys(seasons); seafldpdict=dict.fromkeys(seasons)
        seafldcstddict=dict.fromkeys(seasons); seafldpstddict=dict.fromkeys(seasons)
        seadiffdict=dict.fromkeys(seasons); seadmaskdict=dict.fromkeys(seasons)
        seapcorrdict=dict.fromkeys(seasons)
        seacidict=dict.fromkeys(seasons)

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

        if threed and nonstandardlev:
            fnamec = frootc + field + '_' + '001-061_ts.nc'
            fnamec2 = frootc + field + '_' + '062-121_ts.nc'
            fnamep = frootp + field + '_' + '001-061_ts.nc'
            fnamep2 = frootp + field + '_' + '062-121_ts.nc'
        else:
            fnamec = frootc + field + '_' + timstr + '_ts.nc'
            fnamep = frootp + field + '_' + timstrp + '_ts.nc'
        
        for sii,sea in enumerate(seasons):

            ncparams = {'seas': sea}  

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
                    field='net'
                    conv=1
                else:
                    field='turb'               
            else:

                if threed and nonstandardlev:
                    ncparams['levsel'] = level
                    fldc = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                        cnc.getNCvar(fnamec2,ncfield,**ncparams)*conv,
                                        axis=0)
                    fldp = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                        cnc.getNCvar(fnamep2,ncfield,**ncparams)*conv,
                                        axis=0)
                else:
                    print fnamec
                    fldc = cnc.getNCvar(fnamec,ncfield,timesel=timesel,
                                          **ncparams)*conv
                    fldp = cnc.getNCvar(fnamep,ncfield,timesel=timesel,
                                          **ncparams)*conv
                        
                if sia:
                    fldc = cutl.calc_seaicearea(fldc,lat,lon)
                    fldp = cutl.calc_seaicearea(fldp,lat,lon)
                    

            if isflux: #field in (fluxes,'fsg','turb','net'):
                # mask out regions that are not ice in the control (as P&M 2014 JClim)
                sicnc = cnc.getNCvar(frootc + 'sicn_' + timstr + '_ts.nc','SICN',timesel=timesel,**ncparams)
                
                fldc = ma.masked_where(sicnc<.10,fldc)
                fldp = ma.masked_where(sicnc<.10,fldp)
                
            if plotzonmean:
                fldc = np.mean(fldc[...,:-1],axis=2)# take zonal mean, removing extra lon
                fldp = np.mean(fldp[...,:-1],axis=2)
                
            elif pattcorrwithtime:
                # loop through each year
                # calc pattern corr either yearly or integrated
                years=np.arange(0,fldc.shape[0])
                fldctm = np.mean(fldc[:,lat>corrlim,...],axis=0)
                pcorr = np.zeros(len(years))
                for yr in years:
                    areas = cutl.calc_cellareas(lat,lon)
                    areas = areas[lat>corrlim,:]
                    weights = areas / np.sum(np.sum(areas,axis=1),axis=0)
                    if pattcorryr:
                        # yearly anomaly pattern corr w/ the time mean pattern
                        tmp = fldp[yr,lat>corrlim,...]-fldctm
                    else:
                        tmp = np.mean(fldp[:yr,lat>corrlim,...],axis=0)-fldctm # integrated anomaly pattern
                        
                    tmpmean = np.mean(fldp[:,lat>corrlim,...],axis=0) - fldctm # end pattern to compare against
                    pcorr[yr] = cutl.pattcorr(tmp.flatten()*weights.flatten(),tmpmean.flatten()*weights.flatten())

                seapcorrdict[sea] = pcorr

            elif plotregmean:
                
                #limsdict = con.get_regionlims(region)
                fldc = cutl.calc_regmean(fldc,lat,lon,region)#limsdict)
                fldp = cutl.calc_regmean(fldp,lat,lon,region)#limsdict)
                
            else: # just calculate a polar mean
                if sia:
                    fldc,sh = cutl.calc_totseaicearea(fldc,lat,lon)
                    fldp,sh = cutl.calc_totseaicearea(fldp,lat,lon)
                else:
                    fldc = cutl.polar_mean_areawgted3d(fldc,lat,lon,latlim=seacyclatlim)
                    fldp = cutl.polar_mean_areawgted3d(fldp,lat,lon,latlim=seacyclatlim)

            seafldcstddict[sea] = np.std(fldc,axis=0)
            seafldpstddict[sea] = np.std(fldp,axis=0)
            ttmp,pvtmp = sp.stats.ttest_ind(fldp,fldc,axis=0)
 
            # calculate confidence interval
            # double-check the scale setting
            ci = sp.stats.t.interval(1-siglevel,len(fldp)-1,loc=np.mean(fldp,axis=0)-np.mean(fldc,axis=0),
                                     scale=np.std(fldp,axis=0)/np.sqrt(len(fldp)))
            seacidict[sea] = ci
                
            seatstatdict[sea] = ttmp
            seapvaldict[sea] = pvtmp
            seafldcdict[sea] =  np.mean(fldc,axis=0) # time mean
            seafldpdict[sea] =  np.mean(fldp,axis=0)
            seadiffdict[sea] = np.mean(fldp,axis=0)- np.mean(fldc,axis=0)
            seadmaskdict[sea] = ma.masked_where(pvtmp>siglevel,seadiffdict[sea])

            # end loop through seasons

        fldcstddict[sim] = seafldcstddict
        fldpstddict[sim] = seafldpstddict
        cidict[sim] = seacidict
        tstatdict[sim] = seatstatdict
        pvaldict[sim] = seapvaldict
        fldcdict[sim] = seafldcdict
        fldpdict[sim] = seafldpdict
        flddiffdict[sim] = seadiffdict
        flddmaskdict[sim] = seadmaskdict
        if pattcorrwithtime==1:
            fldpcorrdict[sim] = seapcorrdict

        # end loop through simulations
    if sia:
        field = 'sia' # put back after getting the data. prob not important in the function context

    blob = {}
    blob['ctl'] = fldcdict
    blob['pert'] = fldpdict
    blob['tstat'] = tstatdict
    blob['pval'] = pvaldict
    blob['ci'] = cidict
    blob['ctlstd'] = fldcstddict
    blob['pertstd'] = fldpstddict
    blob['diff'] = flddiffdict
    blob['mask'] = flddmaskdict
    if pattcorrwithtime:
        blob['pcorr'] = fldpcorrdict

    return blob # datablob


def calc_seasonal_cycle(fielddict,coords,sims,withlat=False,loctimesel=None,info=None,siglevel=0.05):
    """ calc_seasonal_cycle(fielddict,coords,sims,withlat=False,loctimesel=None,info=None,siglevel=0.05)

              info: a dict which specifies lots of things, but in particular, which region to average
              
              returns a data blob:
                   blob['ctl'] = fldcdict
                   blob['pert'] = fldpdict
                   blob['tstat'] = tstatdict
                   blob['pval'] = pvaldict
                   blob['ci'] = cidict
                   blob['ctlstd'] = fldcstddict
                   blob['pertstd'] = fldpstddict
                   blob['diff'] = flddiffdict
                   blob['mask'] = flddmaskdict
    """

    print 'calc_seasonal_cycle()'

    sia = False
    
    # make standard function that will be called by "seasonal_cycle", "regional mean", etc...
    #  this can be adapted for lots of processing rather than have multiple copies of code. @@
    
    field=fielddict['field']
    ncfield=fielddict['ncfield']
    fieldstr=fielddict['fieldstr']
    conv=fielddict['conv']
    threed=fielddict['threed']
    isflux=fielddict['isflux']

    model=info['model'] 
    pct=info['pct']   # percentage change?
    nonstandardlev=fielddict['nonstandardlev']
    
    bp=con.get_basepath()
    basepath=bp['basepath'] + model + '/'; subdir=bp['subdir'] # @@ move out of function?

    lat=coords['lat']
    nlat=len(lat)
    lon=coords['lon']
    nlon=len(lon)
        
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
    region = info['region'] # if region is set, it supercedes seacyclatlim! @@

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
        seatstat = np.zeros((12)); seapval = np.zeros((12))
        seafldc = np.zeros((12)); seafldp = np.zeros((12))
        seafldcstd = np.zeros((12)); seafldpstd = np.zeros((12))
        seadiff= np.zeros((12)); seadmask=np.zeros((12))
        seaci = np.zeros((12,2));#@@ might need new shape for ci.
        
        ## seatstatdict=dict.fromkeys(seacyc); seapvaldict=dict.fromkeys(seacyc)
        ## seafldcdict=dict.fromkeys(seacyc); seafldpdict=dict.fromkeys(seacyc)
        ## seafldcstddict=dict.fromkeys(seacyc); seafldpstddict=dict.fromkeys(seacyc)
        ## seadiffdict=dict.fromkeys(seacyc); seadmaskdict=dict.fromkeys(seacyc)
        ## seapcorrdict=dict.fromkeys(seacyc)
        ## seacidict=dict.fromkeys(seacyc)
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
        
        #for sii,sea in enumerate(seacyc):
        if 1: # @@ try getting rid of loop through months.

            print '@@ will have to remove ncparams'
            #@@ncparams = {'monsel': sii+1}

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
                    # @@ really have only tested this one...
                    fldc = cnc.getNCvar(fnamec,ncfield,timesel=timesel)*conv#,
                                        #**ncparams)*conv
                    fldp = cnc.getNCvar(fnamep,ncfield,timesel=timesel)*conv#,
                                        #**ncparams)*conv
                        
                if sia==True:
                    fldc = cutl.calc_seaicearea(fldc,lat,lon)
                    fldp = cutl.calc_seaicearea(fldp,lat,lon)
            # end getting data initially
        
            # now do processing (regional means, tot SIA, etc)
            if isflux:  #field in (fluxes,'fsg','turb','net'):
                # mask out regions that are not ice in the control (as P&M 2014 JClim)
                sicnc = cnc.getNCvar(frootc + 'sicn_' + timstr + '_ts.nc','SICN',timesel=timesel,**ncparams)
                
                fldc = ma.masked_where(sicnc<.10,fldc)
                fldp = ma.masked_where(sicnc<.10,fldp)

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
                    if region!=None:
                        fldc = cutl.calc_regmean(fldc,lat,lon,region)
                        fldp = cutl.calc_regmean(fldp,lat,lon,region)
                    else:
                        fldc = cutl.polar_mean_areawgted3d(fldc,lat,lon,latlim=seacyclatlim)
                        fldp = cutl.polar_mean_areawgted3d(fldp,lat,lon,latlim=seacyclatlim)

            # now save the processed fields
            # @@ got rid of loop through months so have to deal with that here:
            seafldcstd = cutl.calc_monthlystd(fldc)
            seafldpstd = cutl.calc_monthlystd(fldp)
            seatstat,seapval = cutl.calc_monthlytstat(fldp,fldc)
            
            #seafldcstddict[sea] = np.std(fldc,axis=0)
            #seafldpstddict[sea] = np.std(fldp,axis=0)
            #ttmp,pvtmp = sp.stats.ttest_ind(fldp,fldc,axis=0)
                
            #seatstat = ttmp
            #seapval = pvtmp
            seafldc,junk = cutl.climatologize(fldc)
            seafldp,junk = cutl.climatologize(fldp)
            seadiff = seafldp - seafldc

            
            #seafldcdict[sea] =  np.mean(fldc,axis=0) # time mean
            #seafldpdict[sea] =  np.mean(fldp,axis=0)
            #seadiffdict[sea] = np.mean(fldp,axis=0)- np.mean(fldc,axis=0)
            seadmask = ma.masked_where(seapval>siglevel,seadiff)

            # calculate confidence interval on the regional mean
            # double-check the scale setting
            # @@ can do conf interval for any mean, not just regional mean?
            
            for moidx in np.arange(0,12):
                seaci[moidx,:] = sp.stats.t.interval(1-siglevel,len(fldp[moidx::12])-1,
                                                  loc=np.mean(fldp[moidx::12],axis=0)-np.mean(fldc[moidx::12],axis=0),
                                                  scale=np.std(fldp[moidx::12],axis=0)/np.sqrt(len(fldp[moidx::12])))
            #seaci = ci
            """if plotregmean==1:
                # calculate confidence interval on the regional mean
                # double-check the scale setting
                ci = sp.stats.t.interval(1-siglevel,len(fldpzm)-1,loc=np.mean(fldpzm,axis=0)-np.mean(fldczm,axis=0),
                                         scale=np.std(fldpzm,axis=0)/np.sqrt(len(fldpzm)))
                seacidict[sea] = ci
                #print ci # @@@ """

            # end loop through seacyc

        fldcstddict[sim] = seafldcstd
        fldpstddict[sim] = seafldpstd
        cidict[sim] = seaci
        tstatdict[sim] = seatstat
        pvaldict[sim] = seapval
        fldcdict[sim] = seafldc
        fldpdict[sim] = seafldp
        flddiffdict[sim] = seadiff
        flddmaskdict[sim] = seadmask
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
    blob['ci'] = cidict
    blob['ctlstd'] = fldcstddict
    blob['pertstd'] = fldpstddict
    blob['diff'] = flddiffdict
    blob['mask'] = flddmaskdict

    return blob # datablob

    # ======

# ==================================================
def plot_seasonal_cycle(datablob,fielddict,sims,pparams=None,ptypes=('anom',),withlat=False,info=None,printtofile=False,figsize=(6,2.5)):
    """ plot_seasonal_cycle(datablob,pparams,ptypes=('anom',),info=None,printtofile=False,figsize=(6,2.5)):

                  ptypes: 'anom': anomaly sea cycle (pert - ctl)
                          'climo': pert and ctl climatological sea cycle
                          'stddev': pert and ctl monthly standard dev
                          'stdan': monthly std dev anomaly (pert - ctl)
                          
                  figsize = (6,2.5) This is the 'squatter' size

    """
        
    colordict=ccm.get_colordict()
    leglocs='best','best','best','best'
    if info !=None:
        if info['leglocs'] != None:
            leglocs=info['leglocs']

    field=fielddict['field']
    fieldstr=fielddict['fieldstr']

    shadeens = info['shadeens'] # list of ensembles
    savestr = info['savestr']
    seacyclatlim = info['seacyclatlim']
    region = info['region']

    if region!=None:
        regstr=region
    else:
        regstr= 'pol' + str(seacyclatlim) + 'N'

    simspdt = con.get_simpairsdict()
    
    months = con.get_mon()

    fontP = fm.FontProperties()
    fontP.set_size('small')

    # ===== from canam4sims_analens
    # Now plot seasonal cycle
    moidxs=np.arange(1,13)

    if withlat:
        print '@@ implement withlat'
        
        # plot sea cycle of std deviation over ens runs with lat and month
        # may want to eventually add a suplot of each ens run sea cycle of the
        # var (currently done separately in plottimebylat.py) @@
        cmlen=float( plt.cm.get_cmap(cmap).N)
        incr = (cmaxm-cminm) / (cmlen)
        conts = np.arange(cminm,cmaxm+incr,incr)

        fldcdf = pd.DataFrame(fldcdict)
        fldpdf = pd.DataFrame(fldpdict)
        
        tmpcdf = fldcdf.loc[mol]
        tmppdf = fldpdf.loc[mol]

        cestd = np.zeros((len(mol),len(lat)))
        pestd = np.zeros((len(mol),len(lat)))
        # have to do these loops b/c np.array(tmpcdf) comes out as an object of
        # size (mons x sims) with inner arrays of length lat, and can't take std of that
        for mii,mon in enumerate(months): 
            tmpc = np.zeros((5,len(lat)))
            tmpp = np.zeros((5,len(lat)))
            for simii,sim in enumerate(sims[0:5]):# accumulate all ens sims together @@hard-coded sims!
                tmpc[simii,:] = tmpcdf[sim][mon].values
                tmpp[simii,:] = tmppdf[sim][mon].values
            cestd[mii,:] = np.std(tmpc,axis=0)
            pestd[mii,:] = np.std(tmpp,axis=0)

        lats,mos = np.meshgrid(lat,np.arange(0,12))
        fig,axs = plt.subplots()
        fig.set_size_inches(6,5)

        cf = axs.contourf(mos,lats,cestd,cmap=plt.cm.get_cmap(cmap),
                          levels=conts,vmin=cminm,vmax=cmaxm,extend='both')

        axs.set_xlim(0,11)
        axs.set_xticks(range(0,12))
        axs.set_xticklabels(months)
        axs.set_ylim(0,90)
        axs.set_ylabel('Latitude')
        axs.set_xlabel('Month')

        cbar_ax = fig.add_axes([.91,.15, .02,.7])
        fig.colorbar(cf,cax=cbar_ax)
        if printtofile:
            fig.savefig(fieldstr + 'stddev_overens_monxlat_nh.' + suff)
        
    else: # regular seasonal cycle
        print '@@ implement printtofile!'
        
        if 'climo' in ptypes:
 
            fldcdf = pd.DataFrame(datablob['ctl'])
            fldpdf = pd.DataFrame(datablob['pert'])
            
            # climo
            fig,axs = plt.subplots()
            fig.set_size_inches(figsize)

            for skey in sims:
                # try removing [mol] b/c I think it's just an array of months now
                axs.plot(moidxs,fldcdf[skey],color=colordict[skey],linewidth=2)

            for skey in sims:
                axs.plot(moidxs,fldpdf[skey],color=colordict[skey],linewidth=2,linestyle='--')

            plt.legend(sims,leglocs[0], fancybox=True,prop=fontP,framealpha=0.5,ncol=2)
            plt.xlim((1,12))
            plt.gca().set_xticks(range(1,13))
            plt.gca().set_xticklabels(months)
            #plt.xlabel('Month')
            plt.ylabel(fieldstr)
            plt.title('Climos')

            if printtofile: 
                fig.savefig(fieldstr + '_ens_meanBC' + savestr + '_seacyc_' + regstr + '4.pdf')
        # end if 'climo'

        # set up ensemble dict to do max/min/stddev calcs later
        if ('anom' in ptypes) or ('stddev' in ptypes) or ('stdan' in ptypes):           
            simsplt=() #  initialize tuple
            callensdt={}; pallensdt={}
            
            if shadeens != None:
                for ensname in shadeens:
                    censdt={}; pensdt={};
                    print ensname
                    for skey in sims: # for each simulation check if it's in an ensemble to shade

                        if simspdt[skey]['pert']['ensname']==ensname:
                            # create an ensemble dict
                            censdt[skey] = datablob['ctl'][skey]
                            pensdt[skey] = datablob['pert'][skey]

                    callensdt[ensname] = censdt # dict of ens -> dict of sims in ens --> data
                    pallensdt[ensname] = pensdt
                # end loop through ens to shade
                    
                for skey in sims: # only want to do this once: the leftover non-ens sims
                    if simspdt[skey]['pert']['ensname'] not in shadeens:
                        simsplt = simsplt + (skey,)
                                            
                print simsplt
                

        if 'anom' in ptypes:
            print 'anom'
            
            flddiffdf = pd.DataFrame(datablob['diff'])
            fldmaskdf = pd.DataFrame(datablob['mask'])
            seacycylim = info['seacycylim']
            
            # differences
            fig,axs = plt.subplots()
            fig.set_size_inches(figsize)

            for skey in sims:
                axs.plot(moidxs,flddiffdf[skey],color=colordict[skey],linewidth=2)
            for skey in sims:
                axs.plot(moidxs,fldmaskdf[skey],linestyle='none',color=colordict[skey],marker='s')

            plt.legend(sims,leglocs[1], fancybox=True, prop=fontP, framealpha=0.5,ncol=2)
            plt.xlim((1,12))
            if seacycylim!=None:
                plt.ylim(seacycylim)
            plt.gca().set_xticks(range(1,13))
            plt.gca().set_xticklabels(months)
            #plt.xlabel('Month')
            plt.ylabel(fieldstr)
            plt.title('Anomalies')
            axylims = axs.get_ylim()
            if axylims[0]<=0 and axylims[1]>=0:
                axs.axhline(y=0,color='k',linewidth=.5) # @@ figure out how to make it first layer of plot...

            if printtofile: # version 2 loops through sims in order of melt
                fig.savefig(fieldstr + 'diff_ens_meanBC' + savestr + '_seacyc_' + regstr + '4.pdf')


            # sims are 'diffname' in simsdict (con.get_simsdict())
            # sims are also keys to simspairsdict (con.get_simpairsdict())
            #    see if part of ensemble, e.g.: simpairsdict['R1']['pert']['ensname']
            #    If no ensemble, then it's None. Otherwise, it's ensname.
            if shadeens != None:
           
                # differences SHADED
                fig,axs = plt.subplots()
                fig.set_size_inches(figsize)
                fc='0.3','y'
                for eii,ensname in enumerate(shadeens):
                    censdf=pd.DataFrame(callensdt[ensname])
                    pensdf=pd.DataFrame(pallensdt[ensname])
                    demin = np.min(pensdf-censdf,axis=1)
                    demax = np.max(pensdf-censdf,axis=1)
                    axs.fill_between(moidxs,demin,demax,facecolor=fc[eii],alpha=0.2)
                    #fc=fc-.3
                for skey in simsplt:
                    axs.plot(moidxs,flddiffdf[skey],color=colordict[skey],linewidth=2)
                for skey in simsplt:
                    axs.plot(moidxs,fldmaskdf[skey],linestyle='none',color=colordict[skey],marker='s')

                plt.legend(simsplt,leglocs[1], fancybox=True,prop=fontP,framealpha=0.5,ncol=2)
                plt.xlim((1,12))
                if seacycylim!=None:
                    plt.ylim(seacycylim)
                plt.gca().set_xticks(range(1,13))
                plt.gca().set_xticklabels(months)
                #plt.xlabel('Month')
                plt.ylabel(fieldstr)
                plt.title('Anomalies')
                axylims = axs.get_ylim()
                if axylims[0]<=0 and axylims[1]>=0:
                    axs.axhline(y=0,color='k',linewidth=.5) # @@ figure out how to make it first layer of plot...

                if printtofile: # version 2 loops through sims in order of melt
                    fig.savefig(fieldstr + 'diff_ens_meanBC' + savestr + '_seacyc_' + regstr + '4shade.pdf')
        # end if 'anom'

        if 'stddev' in ptypes:
            fldcstddf = pd.DataFrame(datablob['ctlstd'])
            fldpstddf = pd.DataFrame(datablob['pertstd'])
          
            # Standard deviation climos
            fig,axs = plt.subplots()
            fig.set_size_inches(figsize)

            for skey in sims:
                axs.plot(moidxs,fldcstddf[skey],color=colordict[skey],linewidth=2)

            for skey in sims:
                axs.plot(moidxs,fldpstddf[skey],color=colordict[skey],linestyle='--',linewidth=2)

            if shadeens != None:
                col=.3
                for ensname in shadeens:
                    censdf=pd.DataFrame(callensdt[ensname])
                    pensdf=pd.DataFrame(pallensdt[ensname])
                    cestd = np.std(censdf,axis=1)
                    pestd = np.std(pensdf,axis=1)                    

                    axs.plot(range(1,13),cestd,color=str(col),linewidth=3)
                    axs.plot(range(1,13),pestd,color=str(col),linewidth=3,linestyle='--')
                    col=col+.3

            plt.legend(sims,leglocs[2], fancybox=True,prop=fontP, framealpha=0.5, ncol=2)
            plt.xlim((1,12))
            plt.gca().set_xticks(range(1,13))
            plt.gca().set_xticklabels(months)
            #plt.xlabel('Month')
            plt.ylabel(fieldstr)
            plt.title('Sigma')
            axylims = axs.get_ylim()
            if axylims[0]<=0 and axylims[1]>=0:
                axs.axhline(y=0,color='k',linewidth=.5) # @@ figure out how to make it first layer of plot...

            if printtofile:
                fig.savefig(fieldstr + 'STD_ens_meanBC' + savestr + '_seacyc_' + regstr + '4.pdf')
        # end if 'stddev'

        if 'stdan' in ptypes:

            fldcstddf = pd.DataFrame(datablob['ctlstd'])
            fldpstddf = pd.DataFrame(datablob['pertstd'])
            
            # Difference in standard deviation
            fig,axs = plt.subplots()
            fig.set_size_inches(figsize)

            for skey in sims:
                axs.plot(moidxs,fldpstddf[skey]-fldcstddf[skey],color=colordict[skey],linewidth=2)
            col=0.3
            for ensname in shadeens:
                censdf=pd.DataFrame(callensdt[ensname])
                pensdf=pd.DataFrame(pallensdt[ensname])
                cestd = np.std(censdf,axis=1)
                pestd = np.std(pensdf,axis=1)                    
                axs.plot(moidxs,pestd-cestd,color=str(col),linewidth=3)
                col=col+.3

            plt.legend(sims,leglocs[3], fancybox=True, prop=fontP, framealpha=0.5,ncol=2)
            plt.xlim((1,12))
            plt.gca().set_xticks(range(1,13))
            plt.gca().set_xticklabels(months)
            #plt.xlabel('Month')
            plt.ylabel(fieldstr)
            plt.title('Sigma anomalies')
            axylims = axs.get_ylim()
            if axylims[0]<=0 and axylims[1]>=0:
                axs.axhline(y=0,color='k',linewidth=.5) # @@ figure out how to make it first layer of plot...

            if printtofile:
                fig.savefig(fieldstr + 'STDdiff_ens_meanBC' + savestr + '_seacyc_' + regstr + '4.pdf')
        # end if 'stdanom'

def plot_zonmean_byseas(datablob,fielddict,coords,sims,pparams=None,ptypes=('anom',),withlat=False,info=None,printtofile=False,figsize=(6,2.5)):


    seasons = 'SON','DJF','MAM','JJA'
    lat=coords['lat']
    nlat=len(lat)
    
    colordict=ccm.get_colordict()
    leglocs='best','best','best','best'
    if info !=None:
        if info['leglocs'] != None:
            leglocs=info['leglocs']

    field=fielddict['field']
    fieldstr=fielddict['fieldstr']

    shadeens = info['shadeens'] # list of ensembles
    savestr = info['savestr']

    simspdt = con.get_simpairsdict()
    
    fontP = fm.FontProperties()
    fontP.set_size('small')


    # =========================================
    ## seal = list(seasons)
    
    flddiffdf = pd.DataFrame(datablob['diff'])
    fldcdf = pd.DataFrame(datablob['ctl'])
    fldpdf = pd.DataFrame(datablob['pert'])
    fldmaskdf = pd.DataFrame(datablob['mask'])
    fldcstddf = pd.DataFrame(datablob['ctlstd'])
    fldpstddf = pd.DataFrame(datablob['pertstd'])   
    
    ## tmpcdf = fldcdf.loc[seal]
    ## tmppdf = fldpdf.loc[seal]

    ## cestd = np.zeros((4,len(lat)))
    ## pestd = np.zeros((4,len(lat)))
    ## cemin = np.zeros((4,len(lat))) # ctl ens min/max
    ## cemax = np.zeros((4,len(lat)))
    ## pemin = np.zeros((4,len(lat))) # pert ens min/max
    ## pemax = np.zeros((4,len(lat)))
    ## demin = np.zeros((4,len(lat))) # difference min/max
    ## demax = np.zeros((4,len(lat)))
    ## # have to do these loops b/c np.array(tmpcdf) comes out as an object of
    ## # size (seasons x sims) with inner arrays of length lat, and can't take std of that
    ## # also construct a min and max line of ens, for fill_between spread
    ## for sii,sea in enumerate(seasons):
    ##     tmpc = np.zeros((5,len(lat)))
    ##     tmpp = np.zeros((5,len(lat)))
    ##     print '@@ fix this because ens members are hardcoded into sims'
    ##     for simii,sim in enumerate(sims[0:5]):
    ##         # here we just accumulate all ens sims into one matrix
    ##         tmpc[simii,:] = tmpcdf[sim][sea].values
    ##         tmpp[simii,:] = tmppdf[sim][sea].values
    ##     # now do calcs: standard dev, min, max
    ##     cestd[sii,:] = np.std(tmpc,axis=0)
    ##     pestd[sii,:] = np.std(tmpp,axis=0)
    ##     cemax[sii,:] = np.max(tmpc,axis=0)
    ##     cemin[sii,:] = np.min(tmpc,axis=0)
    ##     pemax[sii,:] = np.max(tmpp,axis=0)
    ##     pemin[sii,:] = np.min(tmpp,axis=0)
    ##     demax[sii,:] = np.max(tmpp-tmpc,axis=0)
    ##     demin[sii,:] = np.min(tmpp-tmpc,axis=0)

    # ========================
    # set up ensemble dict to do max/min/stddev calcs later
    if ('anom' in ptypes) or ('stddev' in ptypes) or ('stdan' in ptypes):           
        simsplt=() #  initialize tuple for legends
        callensdt={}; pallensdt={} # contains all ensembles

        if shadeens != None:
            for ensname in shadeens:
                censdt={}; pensdt={}; # all runs w/in one ensemble
                print ensname
                for skey in sims: # for each simulation check if it's in an ensemble to shade

                    if simspdt[skey]['pert']['ensname']==ensname:
                        # create an ensemble dict
                        censdt[skey] = datablob['ctl'][skey]
                        pensdt[skey] = datablob['pert'][skey]

                callensdt[ensname] = censdt # dict of ens -> dict of sims in ens --> data
                pallensdt[ensname] = pensdt
            # end loop through ens to shade

            for skey in sims: # only want to do this once: the leftover non-ens sims
                if simspdt[skey]['pert']['ensname'] not in shadeens:
                    simsplt = simsplt + (skey,)

            print simsplt
    #====================
    
    # climo
    if 'climo' in ptypes:
        fig,axs = plt.subplots(4,1)
        fig.set_size_inches(6,12)
        fig.subplots_adjust(hspace=.2,wspace=.2)

        for sii,sea in enumerate(seasons):
            ax=axs[sii]

            for skey in sims:
                ax.plot(lat,fldcdf[skey][sea],color=colordict[skey],linewidth=2)

            for skey in sims:
                ax.plot(lat,fldpdf[skey][sea],color=colordict[skey],linewidth=2,linestyle='--')

            ax.set_xlim(-90,90)
            ax.set_title(sea)

        ax.set_xlabel('lat')
        if printtofile: # version 2 has new colors and seasons and sims are reordered
            fig.savefig(fieldstr + '_ens_meanBC' + savestr + '_allseassp_zonmean2.pdf')

    # standard deviation
    if 'stddev' in ptypes:

        fig,axs = plt.subplots(4,1)
        fig.set_size_inches(6,12)
        fig.subplots_adjust(hspace=.2,wspace=.2)
        for sii,sea in enumerate(seasons):
            ax=axs[sii]
            cols='.3','y'
            for eii,ensname in enumerate(shadeens):
                censdf = pd.DataFrame(callensdt[ensname]).loc[sea] # returns series of length sim. each element has nlat data
                pensdf = pd.DataFrame(pallensdt[ensname]).loc[sea]
                # hack to get data into a matrix
                tmpc=np.zeros((len(censdf),nlat))
                tmpp=np.zeros((len(pensdf),nlat))
                for ii,ekey in enumerate(censdf.keys()):
                    tmpc[ii,:] = censdf[ekey]
                    tmpp[ii,:] = pensdf[ekey]
                cestd = np.std(tmpc,axis=0)
                pestd = np.std(tmpp,axis=0)
                ax.plot(lat,cestd,color=cols[eii],linewidth=2)
                ax.plot(lat,pestd,color=cols[eii],linewidth=2,linestyle='--')
                
            for skey in sims:
                ax.plot(lat,fldcstddf[skey][sea],color=colordict[skey],linewidth=2)
            #ax.plot(lat,cestd,'k',linewidth=2)

            for skey in sims:
                ax.plot(lat,fldpstddf[skey][sea],color=colordict[skey],linewidth=2,linestyle='--')
            #ax.plot(lat,pestd,'k',linewidth=2,linestyle='--')

            ax.set_xlim(-90,90)
            ax.set_title(sea)

        ax.set_xlabel('lat')
        if printtofile:
            fig.savefig(fieldstr + 'STD_ens_meanBC' + savestr + '_allseassp_zonmean2.pdf')

    # differences
    if 'anom' in ptypes:
        fig,axs = plt.subplots(4,1)
        fig.set_size_inches(6,12)
        fig.subplots_adjust(hspace=.2,wspace=.2)
        for sii,sea in enumerate(seasons):
            ax=axs[sii]

            for skey in sims:
                ax.plot(lat,flddiffdf[skey][sea],color=colordict[skey],linewidth=2)

            for skey in sims:
                ax.plot(lat,fldmaskdf[skey][sea],color=colordict[skey],linestyle='none',marker='s')

            ax.set_xlim(0,90)
            ax.set_title(field + ': ' + sea)

        ax.set_xlabel('lat')
        ax.legend(sims,'upper left', prop=fontP,ncol=2,fancybox=True,framealpha=0.5)
        if printtofile:
            fig.savefig(fieldstr + 'diff_ens_meanBC' + savestr + '_allseassp_zonmean_nh2.pdf')

        # sims are 'diffname' in simsdict (con.get_simsdict())
        # sims are also keys to simspairsdict (con.get_simpairsdict())
        #    see if part of ensemble, e.g.: simpairsdict['R1']['pert']['ensname']
        #    If no ensemble, then it's None. Otherwise, it's ensname.
        if shadeens != None:
           
            # differences SHADED
            fig,axs = plt.subplots(4,1)
            fig.set_size_inches(6,12)
            fig.subplots_adjust(hspace=.2,wspace=.2)
            for sii,sea in enumerate(seasons):
                ax=axs[sii]
                fc='0.3','y'
                for eii,ensname in enumerate(shadeens):
                    censdf = pd.DataFrame(callensdt[ensname]).loc[sea] # returns series of length sim. each element has nlat data
                    pensdf = pd.DataFrame(pallensdt[ensname]).loc[sea]
                    # hack to get data into a matrix
                    tmpc=np.zeros((len(censdf),nlat))
                    tmpp=np.zeros((len(pensdf),nlat))
                    for ii,ekey in enumerate(censdf.keys()):
                        tmpc[ii,:] = censdf[ekey]
                        tmpp[ii,:] = pensdf[ekey]
                    demin = np.min(tmpp-tmpc,axis=0)
                    demax = np.max(tmpp-tmpc,axis=0)
                    ax.fill_between(lat,demin,demax,facecolor=fc[eii],alpha=0.2)

                #ax.fill_between(lat,demin[sii,...],demax[sii,...],facecolor='0.7',alpha=0.2)
                #print '@@ sims other than ens appear to be hard-coded. danger'
                for skey in simsplt: #sims[5:]: # @@ ack, hard-coded! dangerous
                    ax.plot(lat,flddiffdf[skey][sea],color=colordict[skey],linewidth=2)

                for skey in simsplt: #sims[5:]: # @@ danger hard-code again
                    ax.plot(lat,fldmaskdf[skey][sea],color=colordict[skey],linestyle='none',marker='s')

                ax.set_xlim(0,90)
                ax.set_title(field + ': ' + sea)

            ax.set_xlabel('lat')
            ax.legend(simsplt,'upper left', prop=fontP,ncol=2,fancybox=True,framealpha=0.5) # @@ hard-coded
            if printtofile:
                fig.savefig(fieldstr + 'diff_ens_meanBC' + savestr + '_allseassp_zonmean_nh2shade.pdf')

    # standard deviation differences
    if 'stdan' in ptypes:
                
        fig,axs = plt.subplots(4,1)
        fig.set_size_inches(6,12)
        fig.subplots_adjust(hspace=.2,wspace=.2)
        for sii,sea in enumerate(seasons):
            ax=axs[sii]
            cols='.3','y'
            for eii,ensname in enumerate(shadeens):
                censdf = pd.DataFrame(callensdt[ensname]).loc[sea] # returns series of length sim. each element has nlat data
                pensdf = pd.DataFrame(pallensdt[ensname]).loc[sea]
                # hack to get data into a matrix
                tmpc=np.zeros((len(censdf),nlat))
                tmpp=np.zeros((len(pensdf),nlat))
                for ii,ekey in enumerate(censdf.keys()):
                    tmpc[ii,:] = censdf[ekey]
                    tmpp[ii,:] = pensdf[ekey]
                cestd = np.std(tmpc,axis=0)
                pestd = np.std(tmpp,axis=0)
                ax.plot(lat,pestd-cestd,color=cols[eii],linewidth=2) # ensemble mean std dev diff
                
            for skey in sims:
                ax.plot(lat,fldpstddf[skey][sea]-fldcstddf[skey][sea],color=colordict[skey],linewidth=2)
            #ax.plot(lat,pestd-cestd,'k',linewidth=2) # ensemble mean std dev diff

            ax.set_xlim(0,90)
            ax.set_title(field + ': ' + sea)

        ax.set_xlabel('lat')

        if printtofile:
            fig.savefig(fieldstr + 'STDdiff_ens_meanBC' + savestr + '_allseassp_zonmean_nh2.pdf')

    # =========================================

def plot_regmean_byseas(datablob,fielddict,sims,pparams=None,info=None,printtofile=False):
    """
    def plot_regmean_byseas(datablob,fielddict,sims,pparams=None,info=None,printtofile=False):
        
        plot the regional mean by season, like a box and whisker plot, with confidence intervals
        
    """

    seasons='SON','DJF','MAM','JJA'
    
    colordict=ccm.get_colordict()

    field=fielddict['field']
    fieldstr=fielddict['fieldstr']

    shadeens = info['shadeens'] # list of ensembles
    savestr = info['savestr']
    region = info['region']

    simspdt = con.get_simpairsdict()
    
    months = con.get_mon()

    fontP = fm.FontProperties()
    fontP.set_size('small')

    flddiffdf = pd.DataFrame(datablob['diff']) 
    fldcdf = pd.DataFrame(datablob['ctl'])
    fldpdf = pd.DataFrame(datablob['pert']) 
    fldmaskdf = pd.DataFrame(datablob['mask'])
    fldcstddf = pd.DataFrame(datablob['ctlstd'])
    fldpstddf = pd.DataFrame(datablob['pertstd'])
    cidf = pd.DataFrame(datablob['ci'])

    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(10,8)
    for sii,sea in enumerate(seasons):

        ax=axs[sii]
        if sii==0:
            ax.set_title(fieldstr + ' ' + region)
            
        for skeyii,skey in enumerate(sims):
             
            val=flddiffdf[skey][sea]
            ci=cidf[skey][sea]
            
            ax.plot(skeyii,val,color=colordict[skey],marker='s',markersize=8)
            ax.plot((skeyii,skeyii),ci,color=colordict[skey],linewidth=2,marker='_',markersize=6)
            axylims = ax.get_ylim()
            if axylims[0]<=0 and axylims[1]>=0:
                ax.axhline(y=0,color='k',linewidth=.5) # @@ figure out how to make it first layer of plot...
            print sea + ' ' + skey + ' ' + str(val)

        ax.set_xticks(np.arange(0,len(sims)))
        ax.set_xticklabels(sims)
        ax.set_ylabel(sea)
        ax.set_xlim(-.5,len(sims)+.5)
        ax.grid()

    if printtofile:
        fig.savefig(fieldstr + 'diffCI_ens_meanBC' + savestr + '_allseassp_' + region + '.pdf')


def plot_pattcorrwithtime_byseas(datablob,fielddict,sims,pparams=None,info=None,calctype=None,printtofile=False):


    colordict=ccm.get_colordict()

    field=fielddict['field']
    fieldstr=fielddict['fieldstr']

    shadeens = info['shadeens'] # list of ensembles
    savestr = info['savestr']

    simspdt = con.get_simpairsdict()
    
    fontP = fm.FontProperties()
    fontP.set_size('small')

    if  calctype!=None and calctype=='pattcorrwithtimeyr':
         # note that this is pattern correlating a run's
         # pattern each year in time w/ its final pattern, then sorted
        pattcorryr=True
    else:
        pattcorryr=False
        
    # ==============================
    seasons='SON','DJF','MAM','JJA'
    seal=list(seasons)
    
    fldpcorrdf = pd.DataFrame(datablob['pcorr']) #fldpcorrdict)
    tmppcdf = fldpcorrdf.loc[seal] # put seasons in order


    # ========================
    # set up ensemble dict to do max/min/stddev calcs later
    #if ('anom' in ptypes) or ('stddev' in ptypes) or ('stdan' in ptypes):           
    simsplt=() #  initialize tuple for legends
    pcallensdt={} # contains all ensembles

    if shadeens != None:
        for ensname in shadeens:
            pcensdt={}; # all runs w/in one ensemble
            print ensname
            for skey in sims: # for each simulation check if it's in an ensemble to shade

                if simspdt[skey]['pert']['ensname']==ensname:
                    # create an ensemble dict
                    #censdt[skey] = datablob['ctl'][skey]
                    pcensdt[skey] = datablob['pcorr'][skey]

            #callensdt[ensname] = censdt # dict of ens -> dict of sims in ens --> data
            pcallensdt[ensname] = pcensdt
        # end loop through ens to shade

        for skey in sims: # only want to do this once: the leftover non-ens sims
            if simspdt[skey]['pert']['ensname'] not in shadeens:
                simsplt = simsplt + (skey,)

        print simsplt
    #====================

    ## pcemax = np.zeros((4,len(tmppcdf['r1']['DJF'])))# get any element to get time
    ## pcemin = np.zeros((4,len(tmppcdf['r1']['DJF'])))
    ## # have to do these loops b/c np.array(tmppcdf) comes out as an object of
    ## # size (seasons x sims) with inner arrays of length time, and can't take max/min of that
    ## for sii,sea in enumerate(seasons):

        
    ##     tmppc = np.zeros( (5,len(tmppcdf['r1']['DJF'])) )
    ##     for simii,sim in enumerate(simsplt):#[0:5]): # @@ hard-coded ens members!
    ##         # here we just accumulate all ens sims into one matrix
    ##         tmpsorted = tmppcdf[sim][sea].values
    ##         if sea in ('JJA','MAM','SON'):
    ##             if threed==1:
    ##                 tmpsorted=tmpsorted[:-2]
    ##             else:
    ##                 tmpsorted=tmpsorted[:-1] # all timeseries have to be same length
                
    ##         if pattcorryr:
    ##             # sort
    ##             tmpsorted.sort()
                
    ##         tmppc[simii,:] = tmpsorted

    ##     # now do calcs: min, max            
    ##     pcemax[sii,:] = np.max(tmppc,axis=0)
    ##     pcemin[sii,:] = np.min(tmppc,axis=0)
        
    
    ylims = 0,1
    if pattcorryr:
        ylims = -1,1

    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for ii,ax in enumerate(axs.flat):
        sea = seasons[ii]

        for simii,sim in enumerate(sims):
            
            if pattcorryr:
                sortfld = fldpcorrdf[sim][sea]
                sortfld.sort() # sort from smallest to largest patt corr
                ax.plot(sortfld,color=colordict[sim],linewidth=2)
            else:
                ax.plot(fldpcorrdf[sim][sea],color=colordict[sim],linewidth=2)
        ax.set_ylim(ylims)    
        ax.set_title(fieldstr + ': ' + sea)
        
    ax.set_xlabel('lat')
    ax.legend(sims,'lower right', prop=fontP,ncol=2)
    
    if printtofile:
        if pattcorryr:
            fig.savefig(fieldstr + 'diffpattcorryrly_ens_meanBC' + savestr + '_allseassp_nh.pdf')
        else:
            fig.savefig(fieldstr + 'diffpattcorr_ens_meanBC' + savestr + '_allseassp_nh.pdf')

    # with SHADING
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for ii,ax in enumerate(axs.flat):
        sea = seasons[ii]

        # ===
        fc='0.3','y'
        for eii,ensname in enumerate(shadeens):
            pcensdf = pd.DataFrame(pcallensdt[ensname]).loc[sea] # returns series of length sim. each element has ntime data
            fudge = pcensdf.keys()[0]
            ntime = len(pcensdf[fudge])
            
            #pensdf = pd.DataFrame(pallensdt[ensname]).loc[sea]
            # hack to get data into a matrix
            tmppc=np.zeros((len(pcensdf),ntime))
            #tmpp=np.zeros((len(pensdf),ntime))
            for ii,ekey in enumerate(pcensdf.keys()):
                if pattcorryr:
                    tmppc[ii,:] = pcensdf[ekey]
                    tmppc[ii,:].sort()
                else:
                    tmppc[ii,:] = pcensdf[ekey]
                    #tmpp[ii,:] = pensdf[ekey]
            pcemin = np.min(tmppc,axis=0)
            pcemax = np.max(tmppc,axis=0)
            ax.fill_between(np.arange(0,ntime),pcemin,pcemax,facecolor=fc[eii],alpha=0.2)

        # ===

        #ax.fill_between(np.arange(0,pcemin.shape[1]),pcemin[ii,...],pcemax[ii,...],facecolor='0.7',alpha=0.2)
        for simii,sim in enumerate(simsplt):
            
            if pattcorryr:
                sortfld = fldpcorrdf[sim][sea]
                sortfld.sort() # sort from smallest to largest patt corr
                ax.plot(sortfld,color=colordict[sim],linewidth=2)
            else:
                ax.plot(fldpcorrdf[sim][sea],color=colordict[sim],linewidth=2)
        ax.set_ylim(ylims)    
        ax.set_title(fieldstr + ': ' + sea)
        
    ax.set_xlabel('lat')
    ax.legend(simsplt,'lower right', prop=fontP,ncol=2)
    
    if printtofile:
        if pattcorryr:
            fig.savefig(fieldstr + 'diffpattcorryrly_ens_meanBC' + savestr + '_allseassp_nhshade.pdf')
        else:
            fig.savefig(fieldstr + 'diffpattcorr_ens_meanBC' + savestr + '_allseassp_nhshade.pdf')

    
