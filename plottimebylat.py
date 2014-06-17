""" plottimebylat.py
    4/16/2014: plot time (or N) by latitude, showing when significance emergences
    Used canam4sims_stats2.py as a template
    
"""
#import sys as sys
#sys.path.insert(0,'./classes/')

import copy
#import numpy as np
import numpy.ma as ma
#import scipy as sp # scientific python
import scipy.stats
import matplotlib.font_manager as fm
#import matplotlib.pyplot as plt
#import platform as platform
#import cccmaplots as cplt
import constants as con
import cccmautils as cutl
#import cccmaNC as cnc
import cccmacmaps as ccm

#cplt = reload(cplt)
#con = reload(con)
cutl = reload(cutl)
#ccm = reload(ccm)
#cnc = reload(cnc)
# commented out all the things I supposedly fixed in config
# and startup scripts: 4/28/2014

plt.close("all")
plt.ion()

printtofile=0
seasonal = 1
seacycle = 0
pattcorr=1 # do pattern correlation with time instead
pattcorryr=0 # yearly pattern corr with mean instead of time-integrated patt corr
dostd=0 # plot standard dev with time (just the pert run, not anomaly for now@@)

# version 2 uses control climo as baseline (rather than individual times),
#   and full timeseries for ttest
v2=1 # only for seasonal=1

sigtype = 'cont'

seasons = 'DJF','MAM','JJA','SON'
model = 'CanAM4'
threed = 0

############ set simulations #############
casename = 'kemctl1'
casenameh = 'kemhadctl'
timstr = '001-111'
timstr1 = '001-061' # for 3d vars
timstr2 = '062-111' # "
timesel = '0002-01-01,0111-12-31'

# Pert run
casenamep1 = 'kem1pert1'  # 2002-2012 sic and sit
casenamep2 = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
casenamep3 = 'kem1pert3'  # 2002-2012 sic, adjusted sst. control sit
casenameph = 'kemhadpert'
casenamep21 = 'kem1pert2r1'; 
casenamep22 = 'kem1pert2r2'
casenamep23 = 'kem1pert2r3'
casenamep24 = 'kem1pert2r4'
casenamep25 = 'kem1pert2r5'
ensruns = casenamep21, casenamep22, casenamep23, casenamep24, casenamep25
casenamep2e = 'kem1pert2ens'

timstrp = '001-111'
timstrp1 = '001-061' # for 3d vars
timstrp2 = '062-111' # "


casenamep = casenamep2e
rnum=1 # @@ set if casenamep is one of the ens members




# # # ######## set Field info ###################
# st, sicn, sic, gt, pmsl, pcp, hfl, hfs, turb, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)
# OR threed: 'gz','t','u'
field = 'pmsl'
level=30000
#level = 100000  # only for threed vars
thickness=0 # do thickness instead: just for gz
level2=70000 # for thickness calc: typically 1000-700hPa thickness

cmap = 'blue2red_w20' # default cmap
cmapclimo = 'Spectral_r'

if casenamep == casenameph:
    casename = casenameh
    timstr = '001-121'
    timstr2 = '062-121'
    timstrp = '001-121'
    timstrp2 = '062-121'
    timesel = '0002-01-01,0121-12-31'
elif casenamep in ensruns:
    casename = casename + 'r' + str(rnum)
    timstr = '001-121'
    timstr2 = '062-121'
    timstrp = '001-121'
    timstrp2 = '062-121'
    timesel = '0002-01-01,0121-12-31'
elif casenamep == casenamep2e:
    casename = 'kemctl1ens'
    timstr = '001-121'
    timstr2 = '062-121'
    timstrp = '001-121'
    timstrp2 = '062-121'
    timesel = '0002-01-01,0121-12-31'

print 'CONTROL IS ' + casename
print 'PERT IS ' + casenamep


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
    cminm = -2; cmaxm = 2   # monthly
    
    cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
    if dostd==1 and casenamep in casenamep2e:
        cminm=-1.5; cmaxm=1.5

    cmap = 'blue2red_w20'
    #cmap = 'blue2red_20'
elif field == 'gt':
    units = 'K'
    conv = 1  # no conversion
    cmin = -2; cmax = 2  # for anomaly plots
    cminp=-.5; cmaxp=.5 # for when pert is 'ctl'
    cminm = -3; cmaxm = 3   # monthly
    ## print 'small clim!'
    ## cmin = -1; cmax = 1  # for anomaly plots
    cminm = -2; cmaxm = 2   # monthly
    
    cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
    cmap = 'blue2red_w20'
elif field == 'sic':
    units='m'
    conv=1/913.
    cmin=-.5
    cmax=.5
    cminm=-.2
    cmaxm=.2
    cmap = 'red2blue_w20'
elif field == 'sicn':
    units='frac'
    conv = 1
    cmin=-.15
    cmax=.15
    cminm=-.1
    cmaxm=.1
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
    if dostd==1 and casenamep!=casenamep2e:
        cminm=-5; cmaxm=5
elif field == 'pcp':
#    pct=1; units = '%'
    units = 'mm/day' # original: kg m-2 s-1
    conv = 86400  # convert from kg m-2 s-1 to mm/day
    cmin = -.2; cmax = .2  # for anomaly plots
    cminp=-.15; cmaxp=.15
    cminm = -.2; cmaxm = .2
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
    cminm = -.1
    cmaxm = .1
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
    cminm = -1
    cmaxm = 1
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

elif field == 'gz':
    threed = 1

    ncfield = 'PHI'
    units = 'm' # @@
    conv = 1/con.get_g()
    if level==50000:
        cminc = 5200; cmaxc = 5900  # climo 500hPa
    elif level==70000:
        cminc=2800; cmaxc = 3200
    elif level==30000:
        cminc=8600; cmaxc = 9800
        
    cmin = -8 # annual mean
    cmax = 8  # annual mean
    
    if level==30000:
        cmin = -15; cmax = 15
        #cminsea=-20; cmaxsea = 20
        cminm = -15; cmaxm = 15  # for monthly
        if dostd==1 and casenamep!=casenamep2e:
            cminm=-60; cmaxm=60
        elif dostd==1 and casenamep==casenamep2e:
            cminm=-40; cmaxm=40
    else:
        #cminsea = -15; cmaxsea = 15
        cminm = -10; cmaxm = 10  # for monthly

    if thickness:
        cminm = -15; cmaxm = 15

elif field == 't':
    threed=1
    ncfield = 'TEMP'
    units = 'K' # @@
    if level == 30000:
        cminc = 215; cmaxc = 245        
    elif level == 70000:
        cminc = 245; cmaxc = 285
        
    if level == 30000:
        cmin = -.3; cmax = .3
        cminm = -.1; cmaxm = .1  # for monthly
        #cminsea = -.5; cmaxsea = .5
    elif level == 70000:
        cmin = -.3; cmax = .3
        cminm = -.1; cmaxm = .1  # for monthly
        #cminsea = -.5; cmaxsea = .5
    
elif field == 'u':
    threed = 1
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
        cminm = -.25; cmaxm = .25
        #cminsea = -3; cmaxsea = 3
    elif level == 85000:
        cminm = -.15; cmaxm = .15
    else:
        cmin = -1; cmax = 1
        cminm = -.25; cmaxm = .25
        #cminsea = -1; cmaxsea = 1

    cmapclimo='blue2red_20'
      
else:
    print 'No settings for ' + field

## fontP = fm.FontProperties()
## fontP.set_size('small')

# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'

if threed: # @@ try to put 3d vars in the same scripts
    if field == 'turb':
        print '3d turb not ready yet!'
    else:
        fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr1 + '_ts.nc'
        fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp1 + '_ts.nc'
        fnamec2 = basepath + casename + subdir + casename + '_' + field + '_' + timstr2 + '_ts.nc'
        fnamep2 = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp2 + '_ts.nc'
else:

    if field=='turb':
        field='hfl'; fieldb='hfs'
        fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
        fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_ts.nc'
        fnamecb = basepath + casename + subdir + casename + '_' + fieldb + '_' + timstr + '_ts.nc'
        fnamepb = basepath + casenamep + subdir + casenamep + '_' + fieldb + '_' + timstrp + '_ts.nc'

        ## fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv + \
        ##        cnc.getNCvar(fnamecb,fieldb.upper(),timesel=timesel)*conv
        ## fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv+ \
        ##        cnc.getNCvar(fnamepb,fieldb.upper(),timesel=timesel)*conv
        field='turb'
    else:
        fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
        fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_ts.nc'

        ## fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv
        ## fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv

lat = cnc.getNCvar(fnamec,'lat')
lon = cnc.getNCvar(fnamec,'lon')

## nt,nlat,nlon = fldc.shape
## nyr = nt/12.

if sigtype == 'hatch':
    suff = 'png'
else:
    suff = 'pdf'
    
if seasonal: # plot all 4 seasons in a subplot

    cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
    
    incr = (cmaxm-cminm) / (cmlen)
    conts = np.arange(cminm,cmaxm+incr,incr)

    midx=0
    fig,axs = plt.subplots(4,1) 
    fig.set_size_inches(6,8)
    fig.subplots_adjust(hspace=.15,wspace=.05)

    for ax in axs.flat:

        if threed:
            
                
            fldcsea = np.append(cnc.getNCvar(fnamec,ncfield,timesel=timesel,levsel=level,seas=seasons[midx])*conv,
                   cnc.getNCvar(fnamec2,ncfield,levsel=level,seas=seasons[midx])*conv,
                   axis=0)
            fldpsea = np.append(cnc.getNCvar(fnamep,ncfield,timesel=timesel,levsel=level,seas=seasons[midx])*conv,
                   cnc.getNCvar(fnamep2,ncfield,levsel=level,seas=seasons[midx])*conv,
                   axis=0)

            if field == 'gz' and thickness == 1:
                # get level2, then thickness is level2-level
                fldcsea2 = np.append(cnc.getNCvar(fnamec,ncfield,timesel=timesel,levsel=level2,seas=seasons[midx])*conv,
                       cnc.getNCvar(fnamec2,ncfield,levsel=level2,seas=seasons[midx])*conv,
                       axis=0)
                fldpsea2 = np.append(cnc.getNCvar(fnamep,ncfield,timesel=timesel,levsel=level2,seas=seasons[midx])*conv,
                       cnc.getNCvar(fnamep2,ncfield,levsel=level2,seas=seasons[midx])*conv,
                       axis=0)
                fldcsea = fldcsea2 - fldcsea
                fldpsea = fldpsea2 - fldpsea

                
        else:
            if field=='turb':
                field='hfl'; fieldb='hfs'
                fldcsea = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                               seas=seasons[midx])*conv + cnc.getNCvar(fnamecb,fieldb.upper(),
                                               timesel=timesel,seas=seasons[midx])*conv
                fldpsea = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                               seas=seasons[midx])*conv + cnc.getNCvar(fnamepb,fieldb.upper(),
                                               timesel=timesel,seas=seasons[midx])*conv 
                field='turb'
            else:
                fldcsea = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                               seas=seasons[midx])*conv # @@note winter returns 2 fewer elements
                fldpsea = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                               seas=seasons[midx])*conv

        nyr,nlat,nlon = fldcsea.shape
        years = np.arange(0,nyr)
        
        
        if pattcorr==1:

            fldcseazm = fldcsea # do not take zonal mean
            fldpseazm = fldpsea
            fldcseazmtm = np.mean(fldcseazm,axis=0) # time mean now
            fldpseazmtm = np.mean(fldpseazm,axis=0)

            plotd = np.zeros(nyr)
            tstat = np.zeros(nyr)
            pval = np.zeros(nyr)

            testd = np.zeros(nyr) # @@
        else:
            # take a zonal mean for each time
            fldcseazm = np.mean(fldcsea[...,:-1],axis=2) # time x lat
            fldpseazm = np.mean(fldpsea[...,:-1],axis=2)
            fldcseazmtm = np.mean(fldcseazm,axis=0) # time mean now
            fldpseazmtm = np.mean(fldpseazm,axis=0)

            plotd = np.zeros((nyr,nlat))
            tstat = np.zeros((nyr,nlat))
            pval = np.zeros((nyr,nlat))
        
        for yr in years:
            
            if pattcorr==1:
                # do a pattern correlation with time
                
                corrlim=45 # @@

                areas = cutl.calc_cellareas(lat,lon)
                areas = areas[lat>corrlim,:]
                #areas = ma.masked_where(lmask[lat>corrlim,:]==-1,areas)
                weights = areas / np.sum(np.sum(areas,axis=1),axis=0)
                
                # Could also do each year's pattern corr rather than cumulative mean DONE
                # OR, use the ens mean pattern as the pattern to compare against -- not useful I think,
                #        for what we want to know: is the pattern of response for each ensemble member random
                #        or dependent on the boundary condition.
                # @@ Also, would be nice to have multiple variables in one plot and/or multiple simulations

                if pattcorryr:
                    # yearly anomaly pattern corr w/ the time mean pattern
                    tmp = fldpseazm[yr,lat>corrlim,...] - fldcseazmtm[lat>corrlim,...]
                else:
                    # time-integrated anomaly pattern corr w/ the end anomaly pattern
                    tmp = np.mean(fldpseazm[0:yr,lat>corrlim,...],axis=0)-fldcseazmtm[lat>corrlim,...] 

                tmpmean = fldpseazmtm[lat>corrlim,...]-fldcseazmtm[lat>corrlim,...] # the end pattern

                tmpcorr = ma.corrcoef(tmp.flatten()*weights.flatten(),
                                      tmpmean.flatten()*weights.flatten())
                plotd[yr] = tmpcorr[0,1]
                testd[yr] = cutl.pattcorr(tmp.flatten()*weights.flatten(),
                                          tmpmean.flatten()*weights.flatten()) # @@ same result as built-in method
 

                """ from canam4sims_analens.py. modify for here
                ensmem = fldpdict[sim][moidx,lat>corrlim,...] - fldcdict[sim][moidx,lat>corrlim,...]
                obsbc = fldp2[moidx,lat>corrlim,...] - fldc2[moidx,lat>corrlim,...]

                # weight the fields by area
                areas = cutl.calc_cellareas(lat,lon)
                areas = areas[lat>corrlim,:]
                areas = ma.masked_where(lmask[lat>corrlim,:]==-1,areas)
                weights = areas / np.sum(np.sum(areas,axis=1),axis=0)

                ensmem = ma.masked_where(lmask[lat>corrlim,:]==-1,ensmem) # mask out land
                obsbc = ma.masked_where(lmask[lat>corrlim,:]==-1,obsbc) # mask out land

                # @@@ note have to use masked corrcoef!
                # @@ This is probably doing a centered pearson where I really
                # want uncentered (no central mean removed). Could calc myself from
                # http://www.stanford.edu/~maureenh/quals/html/ml/node53.html 
                tmpcorr = ma.corrcoef(ensmem.flatten()*weights.flatten(), obsbc.flatten()*weights.flatten())
                pcorr[moidx] = tmpcorr[0,1]
                @@ see above: doesn't seem to remove mean, gives same result as calc'ing the cov, etc.
                """

            elif dostd==1:
                                
                if yr>=5:
                    plotd[yr,:] = np.std(fldpseazm[0:yr,:],axis=0)#@@-np.std(fldcseazmtm,axis=0)
                    #tstat[yr,:],pval[yr,:] = sp.stats.ttest_ind(fldpseazm[0:yr,:],fldcseazm,axis=0)
                else:
                    #pval[yr,:] = np.ones((1,nlat))
                    plotd[yr,:] = np.NaN
            elif v2==1:
                # here, use the control climo as the baseline (rather than 1 year), and give
                #   the full control timeseries to the ttest every time
                plotd[yr,:] = np.mean(fldpseazm[0:yr,:]-fldcseazmtm,axis=0)
                if yr>5:
                    tstat[yr,:],pval[yr,:] = sp.stats.ttest_ind(fldpseazm[0:yr,:],fldcseazm,axis=0)
                else:
                    pval[yr,:] = np.ones((1,nlat))
            else:
                plotd[yr,:] = np.mean(fldpseazm[0:yr,:]-fldcseazm[0:yr,:],axis=0)
                if yr>5: # start doing stats after 5 years
                    tstat[yr,:],pval[yr,:] = sp.stats.ttest_ind(fldpseazm[0:yr,:],fldcseazm[0:yr,:],axis=0)
                else:
                    pval[yr,:] = np.ones((1,nlat))
            
        if pattcorr: # pattern corr by time
            
            if pattcorryr:
                plotd.sort()
                ax.plot(years,plotd,marker='*')
                ax.set_ylim(-1,1)
            else:
                ax.plot(years,plotd)
                #ax.plot(years,testd,color='r',marker='.') @@ this gives same result!
                ax.set_ylim(0,1)
        else: # lat by time
            lats,times = np.meshgrid(lat,years)
        
            plotfld = plotd #@@these years could be in any order...how take this into account??
        
            cf = ax.contourf(times,lats,plotfld,cmap=plt.cm.get_cmap(cmap),
                              levels=conts,vmin=cminm,vmax=cmaxm,extend='both')
            #cplt.addtsig(ax,pval,lat,years,type=sigtype) # @@ I don't think the stats are right.
            if dostd!=1:
                if sigtype == 'cont':
                    ax.contour(times,lats,pval,levels=[0.05,0.05],colors='k')
                elif sigtype == 'hatch':
                    ax.contourf(times,lats,pval,levels=[0,0.05],colors='none',hatches='.')

            
            ax.set_ylim(0,90)
            
        ax.set_ylabel(seasons[midx])
        if midx == 3:
            ax.set_xlabel('Time')
        ax.set_xlim(0,nyr)

        midx=midx+1
        # ending the loop through seasons

    if pattcorr!=1:
        cbar_ax = fig.add_axes([.91,.15, .02,.7])
        fig.colorbar(cf,cax=cbar_ax)

    if threed:
        if field == 'gz' and thickness==1:
            plt.suptitle(field + ' ' + str(level2/100) + '-' + str(level/100) + ' (' + units + ')')
        else:
            plt.suptitle(field + ' ' + str(level/100) + ' (' + units + ')')

    else:
        plt.suptitle(field + ' (' + units + ')')

        #plt.tight_layout() # makes it worse

    if printtofile:
        prfield = field # print field
        if threed:
            if field == 'gz' and thickness==1:
                prfield = field + str(level2/100) + '-' + str(level/100)
            else:
                prfield = field + str(level/100)
        
        if pattcorr==1:
            if pattcorryr==1: # yearly
                fig.savefig(prfield + 'diffpattcorryrly' + '_' + casenamep +\
                            '_v_' + casename + '_timexlat_seas_nh' + str(corrlim) + 'N.' + suff)
            else: # time-integrated
                fig.savefig(prfield + 'diffpattcorr' + '_' + casenamep +\
                            '_v_' + casename + '_timexlat_seas_nh' + str(corrlim) + 'N.' + suff)
        elif dostd==1:
             fig.savefig(prfield + 'stddev' + '_' + casenamep +\
                        '_timexlat_seas_nh.' + suff)
        elif v2: 
            fig.savefig(prfield + 'diffsig' + sigtype + '_' + casenamep +\
                        '_v_' + casename + '_timexlat_seas_nh_v2.' + suff)
        else:
            fig.savefig(prfield + 'diffsig' + sigtype + '_' + casenamep +\
                        '_v_' + casename + '_timexlat_seas_nh.' + suff)



if seacycle: # want month x lat (or height)

    cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
    
    incr = (cmaxm-cminm) / (cmlen)
    conts = np.arange(cminm,cmaxm+incr,incr)
    months = con.get_mon()
    
    if threed:
        fldc = np.append(cnc.getNCvar(fnamec,ncfield,timesel=timesel,levsel=level)*conv,
               cnc.getNCvar(fnamec2,ncfield,levsel=level)*conv,
               axis=0)
        fldp = np.append(cnc.getNCvar(fnamep,ncfield,timesel=timesel,levsel=level)*conv,
               cnc.getNCvar(fnamep2,ncfield,levsel=level)*conv,
               axis=0)
    else:
        if field=='turb':
            field='hfl'; fieldb='hfs'
            fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv + \
                   cnc.getNCvar(fnamecb,fieldb.upper(),timesel=timesel)*conv
            fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv + \
                   cnc.getNCvar(fnamepb,fieldb.upper(),timesel=timesel)*conv 
            field='turb'
        else:
            fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv
            fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv

    tstat = np.zeros((12,fldc.shape[1]))
    pval = np.zeros((12,fldc.shape[1]))
    monfldczmclimo = np.zeros((12,fldc.shape[1]))
    monfldpzmclimo = np.zeros((12,fldc.shape[1]))
        
    # loop through months, calcing mean pert and stat sig
    for midx in np.arange(0,12):

        if pattcorr:
            print 'do pattern corr @@'
        elif dostd:
            # taken zonal mean
            monfldczm = np.mean(fldc[midx::12,:,:-1],axis=2)
            monfldpzm = np.mean(fldp[midx::12,:,:-1],axis=2)
            # std dev
            monfldczmclimo[midx,:] = np.std(monfldczm,axis=0)
            monfldpzmclimo[midx,:] = np.std(monfldpzm,axis=0)
        else:
            # taken zonal mean
            monfldczm = np.mean(fldc[midx::12,:,:-1],axis=2)
            monfldpzm = np.mean(fldp[midx::12,:,:-1],axis=2)
            # time mean
            monfldczmclimo[midx,:] = np.mean(monfldczm,axis=0)
            monfldpzmclimo[midx,:] = np.mean(monfldpzm,axis=0)
        
            tstat[midx,:],pval[midx,:] = sp.stats.ttest_ind(monfldpzm,monfldczm,axis=0)

    lats,mos = np.meshgrid(lat,np.arange(0,12))

 
    fig2,ax = plt.subplots(1,1)
    fig2.set_size_inches(6, 5)

    if dostd==1:
        plotfld = monfldpzmclimo
    else:
        plotfld = monfldpzmclimo - monfldczmclimo

    cf = ax.contourf(mos,lats,plotfld,cmap=plt.cm.get_cmap(cmap),
                      levels=conts,vmin=cminm,vmax=cmaxm,extend='both')

    if dostd!=1:       
        if sigtype == 'cont':
            ax.contour(mos,lats,pval,levels=[0.05,0.05],colors='k')
        elif sigtype == 'hatch':
            ax.contourf(mos,lats,pval,levels=[0,0.05],colors='none',hatches='.')

    ax.set_xlim(0,11)
    ax.set_xticks(range(0,12))
    ax.set_xticklabels(months)
    ax.set_ylim(0,90)
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Month')

    cbar_ax = fig2.add_axes([.91,.15, .02,.7])
    fig2.colorbar(cf,cax=cbar_ax)

    if threed:
        if field == 'gz' and thickness==1:
            plt.suptitle(field + ' ' + str(level2/100) + '-' + str(level/100) + ' (' + units + ')')
        else:
            plt.suptitle(field + ' ' + str(level/100) + ' (' + units + ')')
    else:
        plt.suptitle(field + ' (' + units + ')')
    #plt.tight_layout() makes it worse

    if printtofile:
        prfield = field

        if threed:
            if field == 'gz' and thickness==1:
                prfield = field + str(level2/100) + '-' + str(level/100)
            else:
                prfield = field + str(level/100)

        if pattcorr==1:
            fig2.savefig(prfield + 'diffpattcorr' + sigtype + '_' + casenamep +\
                    '_v_' + casename + '_monxlat_nh.' + suff)
        elif dostd==1:
            fig2.savefig(prfield + 'stddev' + '_' + casenamep +\
                    '_monxlat_nh.' + suff)
        else:
            fig2.savefig(prfield + 'diffsig' + sigtype + '_' + casenamep +\
                         '_v_' + casename + '_monxlat_nh.' + suff)
