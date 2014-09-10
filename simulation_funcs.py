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
        
    model='CanAM4' # @@ move out of function
    bp=con.get_basepath()
    basepath=bp['basepath'] + model + '/'; subdir=bp['subdir'] # @@ move out of function
    nonstandardlev=False # @@ move out of function
    pct=False # @@ move out of function
    sigoff=False # @@ move out of function
    sigtype='cont' # @@ move out of function

    
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
                cminm=cminmpct
                cmaxm=cmaxmpct
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
