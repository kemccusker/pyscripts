"""
 canam4sims_analens.py
    6/2/2014: taken from canam4sims_stats2.py
              This script is specifically for the CanAM4 ensemble runs
              (kemctl1r? and kem1pert2r?)
              
    2/20/2014: taken from plot_canam4sims_hists.py: 
               calculate & plot statistical properties of the runs

"""

import numpy.ma as ma
import scipy.stats
import matplotlib.cm as cm
from subprocess import call # for doing system calls - not really needed
import datetime as datetime
import matplotlib.colors as col
#import platform as platform
import constants as con      # my module
import cccmautils as cutl    # my module
import matplotlib.font_manager as fm
import copy
import cccmacmaps as ccm
import pandas as pd

# while I'm still creating these modules, have to reload to get changes
cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)
ccm = reload(ccm)
cnc = reload(cnc)

plt.close("all")
plt.ion()

printtofile=True

plotann=0    # seasonal avg map, comparing ens runs and meanBC
plotallmos=0 # monthly maps (@@ not implemented)
seasonal=0 # seasonal maps (SON, DJF, MAM, JJA)
seasvert=0 # seasonal must =1. seasonal vertical zonal means (SON,DJF,MAM,JJA) instead of maps
screen=True # whether to have screen-style vertical zonal means

plotzonmean=0 # plotzonmean,plotseacyc,pattcorrwithtime are mutually exclusive
plotseacyc=0 # plotzonmean,plotseacyc,pattcorrwithtime are mutually exclusive
withlat=0 # plot the seasonal cycle with latitude dimension too (only for plotseacyc=1)@@for now just std over ens
squatseacyc=0 # plot seacycle figs as shorter than wide
squatterseacyc=1 # even shorter, for paper
pattcorrwithtime=0 # plot pattern correlation with time for each ens member
pattcorryr=0 # if 1, do a yearly anomaly pattern rather than time-integrated

plotregmean=1
#latlims=[70,89]; lonlims=[0,359]; region='polcap70' # Polar cap north of 70N
#latlims=[65,89]; lonlims=[0,359]; region='polcap65' # Polar cap north of 65N for NAM proxy
#latlims=[35,60]; lonlims=[40,120]; region='eurasia' # Eurasia 35-60N, 40E-120E
#latlims=[35,60]; lonlims=[240,280]; region='ntham' # North America 35-60N, 120W-80W
latlims=[35,60]; lonlims=[300,360]; region='nthatl' # North Atlantic 35-60N, 60W-0

testhadisst=0 # check which ens member most similar to hadisst
normbystd=0
halftime=False # get only the first 60yrs. make sure to set the other flag the opp
halftime2=False # get only the last 60yrs. make sure to set the other flag the opp

sensruns=False # sensruns only: addr4ct=1 and addsens=1. no meanBC, r mean, or obs
addobs=1 # add mean of kemhad* runs to line plots, seasonal maps. add nsidc if SIA/SIT (@@for now)
addr4ct=0 # add kem1pert2r4ct (constant thickness version of ens4)
addsens=0 # add sensitivity runs (kem1pert1b, kem1pert3)
simsforpaper=False # meanBC, HAD, NSIDC only. best for maps and zonal mean figs (not line plots)
    
latlim = None # None #45 # lat limit for NH plots. Set to None otherwise.
levlim= 100 # level limit for vertical ZM plots (in hPa). ignored if screen=True

sigtype = 'cont' # significance: 'cont' or 'hatch' which is default
sigoff=0 # if 1, don't add significance
siglevel=0.05

# # # ######## set Field info ###################
# gz, t, u, v, q (3D !)
# st, sic, sicn (sia), gt, pmsl, pcp, hfl, hfs, turb, net, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)
field = 'gz'

print field
timeavg = 'DJF'

# only for threed vars
#level = 30000
level = 50000 # 500hPa
#level = 70000
nonstandardlev=False # standards are 700,500,300



seasons = 'SON','DJF','MAM','JJA'

model = 'CanAM4'
threed=0
sia=0 # is the requested field sea ice area
ctstr=''

# # # ########### set Simulations #############
# Control run
bcasename = 'kemctl1'
casename = bcasename + 'ens'
timstr = '001-121'
timesel = '0002-01-01,0121-12-31'
bcasenamep = 'kem1pert2'
casenamep = bcasenamep + 'ens'
timstrp = '001-121'

# ######### set second set of sims (mean BC) ##########
casename2 = 'kemctl1'
casenamep2 = 'kem1pert2'
timstr2='001-121'

if halftime:
    timesel= '0002-01-01,0061-12-31'
elif halftime2:
    timesel='0062-01-01,0121-12-31'
        

if sensruns:
    addobs=0
    addr4ct=1
    addsens=1
    # but don't plot the meanBC or mean of Ens members


# HERE THE SECOND SET OF RUNS IS HadISST BC's
if testhadisst:
    casename2 = 'kemhadctl'
    casenamep2 = 'kemhadpert'
    timstr2=timstr
    field='sicn'


cmap = 'blue2red_w20' # default cmap
cmapclimo = 'Spectral_r'
pct = 0 # if 1, do calculation as a percent




# # # ###########################################
#   Shouldn't have to mod below....

"""  plev = 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, 12500, 
    15000, 17500, 20000, 22500, 25000, 30000, 35000, 40000, 45000, 50000, 
    55000, 60000, 65000, 70000, 75000, 77500, 80000, 82500, 85000, 87500, 
    90000, 92500, 95000, 97500, 100000 ;
"""

if field == 'st':
    units = 'K'
    conv = 1  # no conversion
    cmin = -2; cmax = 2  # for anomaly plots
    cminp=-.5; cmaxp=.5 # for when pert is 'ctl'
    cminm = -3; cmaxm = 3   # monthly
    ## print 'small clim!'
    ## cmin = -1; cmax = 1  for anomaly plots
    ## cminm = -1.5; cmaxm = 1.5   monthly
    
    cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
    cminn = -5; cmaxn = 5 # for norm by std
    cmap = 'blue2red_w20'

    if plotseacyc==1  and withlat==1:
        cminm=-.5; cmaxm=.5 # @@ will have to update this when add subplots
        
    leglocs = 'upper left', 'upper left', 'upper right', 'upper left'
    seacycylim=(-.2,1.1)
elif field == 'sic':
    units='m'
    conv=1/913.
    cmin=-.5
    cmax=.5
    cminm=-.5
    cmaxm=.5
    cmap = 'red2blue_w20'
    leglocs = 'lower left', 'lower left', 'upper left', 'upper left'
elif field == 'sicn' or field == 'sia':
    units = 'frac'
    conv=1
    cmin=-.15; cmax=.15
    cminm=-.15; cmaxm=.15
    cmap = 'red2blue_w20'
    leglocs = 'lower left', 'lower left', 'upper left', 'upper left'
    seacycylim=(-2e12,0)
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
    ## print 'new cmap and small clim! @@ '
    ## cmap = 'blue2red_w20'
    ## cminm=-1; cmaxm=1  for monthly maps, small clim
    
    cminn = -1; cmaxn = 1 # for norm by std

    if plotseacyc==1  and withlat==1:
        cminm=-1; cmaxm=1 # @@ will have to update this when add subplots
    leglocs = 'lower left', 'lower left', 'upper center', 'lower left'
    seacycylim=(-1,0.4)
elif field == 'pcp':
    units = 'mm/day' # original: kg m-2 s-1
    
    pct=1; units = '%'
    
    conv = 86400  # convert from kg m-2 s-1 to mm/day
    cmin = -.2; cmax = .2  # for anomaly plots
    cminp=-.15; cmaxp=.15
    cminm = -.2; cmaxm = .2

    ## print '@@ big clim!'
    ## cmin = -.75; cmax = .75 
    ## cminm = -.75; cmaxm = .75

    ## print '@@ medium clim!'
    ## cmin = -.5; cmax = .5 
    ## cminm = -.5; cmaxm = .5

    
    #cmap = 'PuOr'
    cmap = 'brown2blue_16w'
    cminpct=-12; cmaxpct=12
    cminmpct=-20; cmaxmpct=20
    cminmp =-.25; cmaxmp=.25
    cminpctp=-8; cmaxpctp=8
    cminpctmp=-12; cmaxpctmp=12
    leglocs = 'upper left', 'upper left', 'upper left', 'upper left'
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
    #print " 'net ' not yet implemented! @@"
    units = 'W/m2'
    conv=1
    cmin=-10
    cmax=10
    cminm = -20
    cmaxm = 20
    cmap='blue2red_20'
    leglocs = 'upper left', 'upper left', 'upper left', 'upper left'
    seacycylim=(-5,25)
elif field == 'flg': # net downward LW at the sfc.
    units = 'W/m2'
    conv = -1 # so positive heats atmos
    cmin = -5
    cmax = 5
    cminm = -8
    cmaxm = 8
    leglocs = 'upper left', 'lower left', 'upper left', 'upper left'
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
    leglocs = 'upper left', 'upper left', 'upper center', 'upper left'

elif field == 'zn': # snow depth (m)
#    pct=1; units='%'
    units = 'cm'
    conv = 100; # convert to cm
    cmap = 'brown2blue_16w'
    cmin = -2
    cmax = 2
    cminm = -2.5; cmaxm = 2.5
    cminpct=-10
    cmaxpct=10
    cminmpct=-10
    cmaxmpct=10
    leglocs = 'upper right', 'lower left', 'upper right', 'lower right'

elif field == 'su':
    units = 'm/s'
    conv = 1;
    cmap = 'blue2red_20'
    cmin = -1; cmax = 1
    cminm = -1; cmaxm = 1
    cminp = -.5; cmaxp=.5
    cminmp = -.5; cmaxmp=.5
    leglocs = 'upper left', 'lower left', 'upper left', 'upper left'
elif field == 'sv':
    units = 'm/s'
    conv = 1;
    cmap = 'blue2red_20'
    cmin = -.5
    cmax = .5
    cminm = -.5
    cmaxm = .5

elif field == 't':
    threed = 1

    conv=1
    ncfield = 'TEMP'
    units = 'K' # @@
    ## if level == 30000:
    ##     cminc = 215; cmaxc = 245        
    ## elif level == 70000:
    ##     cminc = 245; cmaxc = 285
        
    if level == 30000:
        cmin = -.3; cmax = .3
        cminm = -.5; cmaxm = .5  # for monthly
        #cminsea = -.5; cmaxsea = .5
    elif level == 70000:
        cmin = -.3; cmax = .3
        cminm = -.5; cmaxm = .5  # for monthly
        #cminsea = -.5; cmaxsea = .5

    if seasvert:
        if screen:
            cmin=-2.5; cmax=2.5
            cminm=-2.5; cmaxm=2.5
        else:
            cmin = -.5; cmax = .5
            cminm = -.8; cmaxm = .8 
        
elif field == 'u':
    threed = 1

    conv=1
    ncfield = 'U'
    units = 'm/s' #@@
    ## if level==50000:
    ##     cminc=-25; cmaxc=25
    ## elif level==70000:
    ##     cminc=-15; cmaxc=15
    ## elif level == 30000:
    ##     cminc=-40; cmaxc=40

    if level == 30000:
        cmin = -2; cmax = 2
        cminm = -3; cmaxm = 3
        #cminsea = -3; cmaxsea = 3
    else:
        cmin = -1; cmax = 1
        cminm = -1; cmaxm = 1
        #cminsea = -1; cmaxsea = 1

    if seasvert:
        cmin=-.5; cmax=.5
        cminm=-1; cmaxm=1
    #cmapclimo='blue2red_20'

elif field == 'gz':
    threed=1

    ncfield = 'PHI'
    units = 'm' # @@
    conv = 1/con.get_g()
    ## if level==50000:
    ##     cminc = 5200; cmaxc = 5900  # climo 500hPa
    ## elif level==70000:
    ##     cminc=2800; cmaxc = 3200
    ## elif level==30000:
    ##     cminc=8600; cmaxc = 9800
    ## elif level==100000: # use this for thickness calc 1000-700
    ##     cminc=2650; cmaxc = 3050
        
    cmin = -8 # annual mean
    cmax = 8  # annual mean
    
    if level==30000:
        cmin = -15; cmax = 15
        #cminsea=-20; cmaxsea = 20
        cminm = -20; cmaxm = 20  # for monthly
    else:
        #cminsea = -15; cmaxsea = 15
        cminm = -15; cmaxm = 15  # for monthly

    if seasvert:
        if screen:
            cmin=-10; cmax=10
            cminm=-25; cmaxm=25
        else:
            cmin = -10; cmax = 10
            cminm = -15; cmaxm = 15 

else:
    print 'No settings for ' + field

fluxes = 'hfl','hfs','flg' # LH, SH, LWdown

# # # ########## Read NC data ###############

bp=con.get_basepath()
basepath=bp['basepath'] + model + '/'; subdir=bp['subdir']


# set filename for getting lat,lon
if threed==1: # either map of a level, or a zonal mean
    if seasonal==1 and seasvert==1:
        fieldstr=field+'ZM'
        field=field+'ZM'
        #fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc' # for lat,lon,plev
    else:
        fieldstr=field+str(level/100) # figure names. want integer for string
        field=field+str(level) # read in files
        #fnamec = basepath + casename + subdir + casename + '_' + field + '_001-061_ts.nc' # for lat,lon
    fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc' # for lat,lon

else:
    fieldstr=field
    if field == 'sia':
        sia=1
        field='sicn' # only really sia for zonal mean and seasonal cycle
        
    fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
    if sia==1:
        field='sia'
    if field in ('turb','net'): # just for lat/lon
        fnamec = basepath + casename + subdir + casename + '_flg_' + timstr + '_ts.nc'

lat = cnc.getNCvar(fnamec,'lat')
lon = cnc.getNCvar(fnamec,'lon')
if threed:
    lev=cnc.getNCvar(fnamec,'plev')
    nlev=len(lev)
nlat = len(lat)
nlon = len(lon)

if sigtype=='cont' or sigoff==1:
    suff='pdf'
else:
    suff='png'

if simsforpaper:
    sims = 'kemhad','kemnsidc',''
    ctstr = '_forpap'
    addobs=0
    addr4ct=0
    addsens=0
    seasons = ('SON','DJF')
else:

    # order ens simulations in order of most ice loss in melt season to least. Then ens mean, PERT2, observations if requested
    sims = 'r1','r4','r3','r5','r2','ens','' # suffixes to bcasename and bcasenamep
    if addobs:
        #if field in ('sia','sicn'): # for now@@
        sims = sims + ('kemhad','kemnsidc')
        #print 'adding yrs 1-61 for nsidc for now@@'
        #else:
        #    sims = sims + ('kemhad',)
        #    print '@@ not adding nsidc yet for field: ' + field

    if sensruns: # add sensitivity runs. don't plot meanBC, mean of ens
        sims = sims[0:5] + ('r4ct','kem1pert1b','kem1pert3')
        ctstr = '_sensruns'
    else:
        if addr4ct:
            sims = sims + ('r4ct',)
            ctstr = '_r4ct' # for figure filenames
        if addsens:
            sims = sims + ('kem1pert1b','kem1pert3') # control is kemctl1 (or '' key)
            ctstr = ctstr + 'sens'
if halftime:
    ctstr = ctstr + '_60yrs' # @@
elif halftime2:
    ctstr = ctstr + '_60yrs2' # @@
    
print sims
#ensmems=np.arange(0,5)

if plotann:

    if field=='turb':
        print 'not fully implemented! @@, no second set of sims'

        field='hfl'; fieldb='hfs'
        fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
        fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_ts.nc'
        fnamecb = basepath + casename + subdir + casename + '_' + fieldb + '_' + timstr + '_ts.nc'
        fnamepb = basepath + casenamep + subdir + casenamep + '_' + fieldb + '_' + timstrp + '_ts.nc'

        fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv + \
               cnc.getNCvar(fnamecb,fieldb.upper(),timesel=timesel)*conv
        fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv+ \
               cnc.getNCvar(fnamepb,fieldb.upper(),timesel=timesel)*conv

        field='turb'
    else:
        if threed==0:
            if field == 'sia':
                sia=1
                field='sicn' # only really sia for zonal mean and seasonal cycle

            # get ensemble mean
            fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
            fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_ts.nc'

            fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv
            fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv

            # get original runs (mean BC)
            fnamec2 = basepath + casename2 + subdir + casename2 + '_' + field + '_' + timstr2 + '_ts.nc'
            fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstr2 + '_ts.nc'

            fldc2 = cnc.getNCvar(fnamec2,field.upper(),timesel=timesel)*conv
            fldp2 = cnc.getNCvar(fnamep2,field.upper(),timesel=timesel)*conv

            if sia==1:
                field='sia'
        else: # is threed
            # get ensemble mean

            if nonstandardlev: # leave this in in case want a different level, not tested w/ 3D zonal mean
                # read 3D variables at specified level, 2 timeseries files
                frootc = basepath + casename + subdir + casename + '_' + field + '_'
                frootp =  basepath + casenamep + subdir + casenamep + '_' + field + '_'
                fnamec = frootc + '001-061_ts.nc'

                fldc = np.append(cnc.getNCvar(fnamec,ncfield,
                                              timesel='0002-01-01,061-12-31',levsel=level)*conv,
                                 cnc.getNCvar(frootc+'062-121_ts.nc',ncfield,levsel=level)*conv,
                                 axis=0)
                fldp = np.append(cnc.getNCvar(frootp+'001-061_ts.nc',ncfield,
                                              timesel='0002-01-01,061-12-31',levsel=level)*conv,
                                 cnc.getNCvar(frootp+'062-121_ts.nc',ncfield,levsel=level)*conv,
                                 axis=0)

                # get original runs (mean BC)
                frootc2 = basepath + casename2 + subdir + casename2 + '_' + field + '_'
                frootp2 =  basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_'

                fldc2 = np.append(cnc.getNCvar(frootc2+'001-061_ts.nc',ncfield,
                                              timesel='0002-01-01,061-12-31',levsel=level)*conv,
                                 cnc.getNCvar(frootc2+'062-121_ts.nc',ncfield,levsel=level)*conv,
                                 axis=0)
                fldp2 = np.append(cnc.getNCvar(frootp2+'001-061_ts.nc',ncfield,
                                              timesel='0002-01-01,061-12-31',levsel=level)*conv,
                                 cnc.getNCvar(frootp2+'062-121_ts.nc',ncfield,levsel=level)*conv,
                                 axis=0)
            else: # read in a standard level from processed file (default) @@reorg since same as above almost
                fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
                fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_ts.nc'

                fldc = cnc.getNCvar(fnamec,ncfield,timesel=timesel)*conv
                fldp = cnc.getNCvar(fnamep,ncfield,timesel=timesel)*conv

                # get original runs (mean BC)
                fnamec2 = basepath + casename2 + subdir + casename2 + '_' + field + '_' + timstr2 + '_ts.nc'
                fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstr2 + '_ts.nc'

                fldc2 = cnc.getNCvar(fnamec2,ncfield,timesel=timesel)*conv
                fldp2 = cnc.getNCvar(fnamep2,ncfield,timesel=timesel)*conv                

    # annual time-series (3d)
    seastsc = cutl.seasonalize_monthlyts(fldc,timeavg)
    seastsp = cutl.seasonalize_monthlyts(fldp,timeavg)

    nt,nlev,nlat = seastsc.shape # @@the var names are "wrong" but work fine in the script as written

    tstat,pval = sp.stats.ttest_ind(seastsp,seastsc,axis=0)
    seastdc = np.std(seastsc,axis=0)
    seastdp = np.std(seastsp,axis=0)

    seastsc2 = cutl.seasonalize_monthlyts(fldc2,timeavg)
    seastsp2 = cutl.seasonalize_monthlyts(fldp2,timeavg)

    nt2,nlev,nlat = seastsc.shape # @@the var names are "wrong" but work fine in the script as written

    tstat2,pval2 = sp.stats.ttest_ind(seastsp2,seastsc2,axis=0)
    seastdc2 = np.std(seastsc2,axis=0)
    seastdp2 = np.std(seastsp2,axis=0)

    

    #tstatb,pvalb = sp.stats.ttest_ind(anntsp,anntsc,axis=0,equal_var=False) # basically the same as above
    # Note that NaN is returned for zero variance (I think..from googling..)
    # If that is the case, pcolormesh() needs a masked_array rather than ndarray (??)
    #  : http://stackoverflow.com/questions/7778343/pcolormesh-with-missing-values

#if plotann: @@ moving this might require defining some things again...test...6/17
    print timeavg
    if field=='sia':
        print 'Plotting maps of SICN instead of sia'
        field='sicn'

    if pct:
        seastsctm = np.mean(seastsc,0)
        seastsctm = ma.masked_where(seastsctm<=0.01,seastsctm)

        plotfld = np.mean(seastsp-seastsc,0) / seastsctm * 100

        seastsctm2 = np.mean(seastsc2,0)
        seastsctm2 = ma.masked_where(seastsctm2<=0.01,seastsctm2)

        plotfld2 = np.mean(seastsp2-seastsc2,0) / seastsctm2 * 100
        cmin=cminpct
        cmax=cmaxpct
    else:
        if normbystd:
            plotfld = (np.mean(seastsp,0)-np.mean(seastsc,0))/seastdc
            plotfld2 = (np.mean(seastsp2,0)-np.mean(seastsc2,0))/seastdc2
            cmin=cminn; cmax=cmaxn 
            sigoff=1
        else:
            plotfld = np.mean(seastsp,0)-np.mean(seastsc,0)
            plotfld2 = np.mean(seastsp2,0)-np.mean(seastsc2,0)

    # Plot mean of ens runs, the meanBC run, and their difference
    fig1,ax1 = plt.subplots(1,3) 
    ax = ax1[0]
    bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,type='nh',\
                        title='ens',units=units,axis=ax,suppcb=1)
    if sigoff==0:
        cplt.addtsigm(bm,pval,lat,lon,type=sigtype) # add significance info (hatching. for contour, type='contour')

    ax = ax1[1]
    bm,pc = cplt.kemmap(plotfld2,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,type='nh',\
                        title=casenamep2 + '-' + casename2,units=units,axis=ax,suppcb=1)
    if sigoff==0:
        cplt.addtsigm(bm,pval2,lat,lon,type=sigtype)

    ax = ax1[2]
    bm,pc = cplt.kemmap(plotfld-plotfld2,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,type='nh',\
                        title='ens-meanBC',units=units,axis=ax,suppcb=1)
    if sigoff==0:
        tstattmp,pvaltmp = sp.stats.ttest_ind(seastsp-seastsc,seastsp2-seastsc2,axis=0)
        cplt.addtsigm(bm,pvaltmp,lat,lon,type=sigtype)

    cbar_ax = fig1.add_axes([.91,.25, .02,.5])
    fig1.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    plt.suptitle(timeavg + ': ' + fieldstr)
    
    if printtofile:
        if sigoff==0:
            sigstr='sig' + sigtype
        else:
            if normbystd:
                sigstr='norm'
            else:
                sigstr=''
            
        if pct:
            fig1.savefig(fieldstr + 'pctdiff' + sigstr + '_ens_v_meanBC_' + timeavg + '_nh.' + suff )
        else:
            fig1.savefig(fieldstr + 'diff' + sigstr + '_ens_v_meanBC_' + timeavg + '_nh.' + suff)

    # ==================================================
    # do a subplot with each ens member, plus the mean
    ridx=0
    fig,spax = plt.subplots(2,3)
    for ax in spax.flat:

        if ridx==len(spax.flat)-1: # last spot is for the ens mean
            if threed==0:
                fnamec = basepath + bcasename + 'ens' + subdir + bcasename +\
                         'ens_' + field + '_' + timstr + '_ts.nc'
                fnamep = basepath + bcasenamep + 'ens' + subdir + bcasenamep +\
                         'ens_' + field + '_' + timstrp + '_ts.nc'
                seasfldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,seas=timeavg)*conv
                seasfldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,seas=timeavg)*conv
            else:
                # read 3D variables at specified level, 2 timeseries files
                if nonstandardlev:
                    frootc = basepath + bcasename + 'ens' + subdir + bcasename + 'ens_' + field + '_'
                    frootp =  basepath + bcasenamep + 'ens' + subdir + bcasenamep + 'ens_' + field + '_'
                    fnamec = frootc + '001-061_ts.nc'

                    seasfldc = np.append(cnc.getNCvar(fnamec,ncfield,
                                                  timesel='0002-01-01,061-12-31',levsel=level,seas=timeavg)*conv,
                                     cnc.getNCvar(frootc+'062-121_ts.nc',ncfield,levsel=level,seas=timeavg)*conv,
                                     axis=0)
                    seasfldp = np.append(cnc.getNCvar(frootp+'001-061_ts.nc',ncfield,
                                                  timesel='0002-01-01,061-12-31',levsel=level,seas=timeavg)*conv,
                                     cnc.getNCvar(frootp+'062-121_ts.nc',ncfield,levsel=level,seas=timeavg)*conv,
                                     axis=0)
                else: # read from already processed file (default). @@ should reorg b/c this is almost same as above
                    fnamec = basepath + bcasename + 'ens' + subdir + bcasename +\
                         'ens_' + field + '_' + timstr + '_ts.nc'
                    fnamep = basepath + bcasenamep + 'ens' + subdir + bcasenamep +\
                             'ens_' + field + '_' + timstrp + '_ts.nc'
                    seasfldc = cnc.getNCvar(fnamec,ncfield,timesel=timesel,seas=timeavg)*conv
                    seasfldp = cnc.getNCvar(fnamep,ncfield,timesel=timesel,seas=timeavg)*conv

            ttl = 'ens'
            
        else:
            if threed==0:
                fnamec = basepath + bcasename + 'r' + str(ridx+1) + subdir + bcasename +\
                         'r' + str(ridx+1) + '_' + field + '_' + timstr + '_ts.nc'
                fnamep = basepath + bcasenamep + 'r' + str(ridx+1) + subdir + bcasenamep +\
                         'r' + str(ridx+1) + '_' + field + '_' + timstrp + '_ts.nc'
                seasfldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,seas=timeavg)*conv
                seasfldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,seas=timeavg)*conv
            else:
                if nonstandardlev:
                    frootc = basepath + bcasename + 'r' + str(ridx+1) + subdir + bcasename +\
                             'r' + str(ridx+1) + '_' + field + '_'
                    frootp =  basepath + bcasenamep + 'r' + str(ridx+1) + subdir + bcasenamep +\
                             'r' + str(ridx+1) + '_' + field + '_'
                    fnamec = frootc + '001-061_ts.nc'

                    seasfldc = np.append(cnc.getNCvar(fnamec,ncfield,
                                                  timesel='0002-01-01,061-12-31',levsel=level,seas=timeavg)*conv,
                                     cnc.getNCvar(frootc+'062-121_ts.nc',ncfield,levsel=level,seas=timeavg)*conv,
                                     axis=0)
                    seasfldp = np.append(cnc.getNCvar(frootp+'001-061_ts.nc',ncfield,
                                                  timesel='0002-01-01,061-12-31',levsel=level,seas=timeavg)*conv,
                                     cnc.getNCvar(frootp+'062-121_ts.nc',ncfield,levsel=level,seas=timeavg)*conv,
                                     axis=0)
                else:
                    fnamec = basepath + bcasename + 'r' + str(ridx+1) + subdir + bcasename +\
                         'r' + str(ridx+1) + '_' + field + '_' + timstr + '_ts.nc'
                    fnamep = basepath + bcasenamep + 'r' + str(ridx+1) + subdir + bcasenamep +\
                             'r' + str(ridx+1) + '_' + field + '_' + timstrp + '_ts.nc'
                    seasfldc = cnc.getNCvar(fnamec,ncfield,timesel=timesel,seas=timeavg)*conv
                    seasfldp = cnc.getNCvar(fnamep,ncfield,timesel=timesel,seas=timeavg)*conv
            
            ttl = 'r' + str(ridx+1)

        seastdc = np.std(seasfldc,axis=0)
        seastdp = np.std(seasfldp,axis=0)

        tstat,pval = sp.stats.ttest_ind(seasfldp,seasfldc,axis=0)
        plotfld = np.mean(seasfldp,axis=0)-np.mean(seasfldc,axis=0)
        if normbystd:
            plotfld = plotfld/seastdc
            cmin=cminn; cmax=cmaxn
            sigoff=1

        bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,type='nh',\
                            title=ttl,units=units,axis=ax,suppcb=1)
        if sigoff==0:
            cplt.addtsigm(bm,pval,lat,lon,type=sigtype)
        
        ridx = ridx+1

    cbar_ax = fig.add_axes([.91,.25, .02,.5])
    fig.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    plt.suptitle(timeavg + ': ' + fieldstr)
    
    if printtofile:

        if sigoff==0:
            sigstr='sig' + sigtype
        else:
            if normbystd:
                sigstr='norm'
            else:
                sigstr=''           
        if pct:
            fig.savefig(fieldstr + 'pctdiff' + sigstr + '_enssubplot_' + timeavg + '_nh.' + suff )
        else:
            fig.savefig(fieldstr + 'diff' + sigstr + '_enssubplot_' + timeavg + '_nh.' + suff)


# end plotann

#sigs = np.ones((12,fldc.shape[1],fldc.shape[2])) unused

if addobs:
    obsstr = 'obs' #'had' # for figure file name
else:
    obsstr = ''

months=con.get_mon()

if plotallmos:
    print 'plotallmos not implemented'
    
    ## title = field + ": " + casenamep + "-" + casename
    ## midx=0
    ## fig, spax = plt.subplots(2,6)
    ## #fig.set_size_inches(12,6)
    ## fig.set_size_inches(12,4.5)
    ## fig.subplots_adjust(hspace=0,wspace=0)

    ## for ax in spax.flat:

    ##     monfldc = fldc[midx::12,:,:]
    ##     monfldp = fldp[midx::12,:,:]

    ##     tstat,pval = sp.stats.ttest_ind(monfldp,monfldc,axis=0)
    ##     sigs[midx,:,:] = ma.masked_where(pval>0.05,pval) 

    ##     if pct:
    ##         monfldctm = np.mean(monfldc,0)
    ##         monfldctm = ma.masked_where(monfldctm<=0.01,monfldctm)
    ##         plotfld = np.mean(monfldp-monfldc,0) / monfldctm *100
    ##         cminm=cminmpct
    ##         cmaxm=cmaxmpct
    ##     else:
    ##         plotfld = np.mean(monfldp,0)-np.mean(monfldc,0)

    ##     bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',\
    ##                  title=months[midx],axis=ax,suppcb=1)
    ##     ax.set_title(months[midx])
    ##     cplt.addtsigm(bm,pval,lat,lon,type=sigtype)

    ##     midx = midx+1

    ## cbar_ax = fig.add_axes([.91,.25, .02,.5])
    ## fig.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    ## plt.suptitle(title)
    ## if printtofile:
    ##     if pct:
    ##         fig.savefig(field + 'pctdiffsig' + sigtype + '_' + casenamep +\
    ##                     '_v_' + casename + '_allmos_nh.' + suff)
    ##     else:
    ##         fig.savefig(field + 'diffsig' + sigtype + '_' + casenamep +\
    ##                 '_v_' + casename + '_allmos_nh.' + suff)

# done with if plotallmos


if seasonal:

    if field=='sia':
        print 'Plotting maps of SICN instead of sia'
        field='sicn'

    cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division


    if seasvert:
        tstat = np.zeros((len(seasons),nlev,nlat))
        pval = np.zeros((len(seasons),nlev,nlat))
        fldcallseas = np.zeros((len(seasons),nlev,nlat))
        fldpallseas = np.zeros((len(seasons),nlev,nlat))
    else:
        tstat = np.zeros((len(seasons),nlat,nlon))
        pval = np.zeros((len(seasons),nlat,nlon))
        fldcallseas = np.zeros((len(seasons),nlat,nlon))
        fldpallseas = np.zeros((len(seasons),nlat,nlon))
    
    incr = (cmaxm-cminm) / (cmlen)
    conts = np.arange(cminm,cmaxm+incr,incr)

    fig6,ax6 = plt.subplots(len(seasons),len(sims)) # 1 row for e/ of 5 ens members, plus mean, plus meanBC
    fig6.set_size_inches(12,8)  
    fig6.subplots_adjust(hspace=.15,wspace=.05)

    lastcol=len(sims)-1
    
    for ridx,sim in enumerate(sims): # traverse cols

        cidx=0
        if sim in ('kemhad','kemnsidc'):
            frootc = basepath + sim + 'ctl' + subdir + sim + 'ctl' + '_' + field + '_'
            frootp = basepath + sim + 'pert' + subdir + sim + 'pert' + '_' + field + '_'
            if sim=='kemhad':
                rowl='had'
            else:
                rowl='nsidc'
        elif sim in ('kem1pert1b','kem1pert3'):
            frootc = basepath + 'kemctl1' + subdir + 'kemctl1' + '_' + field + '_'
            frootp = basepath + sim + subdir + sim + '_' + field + '_'
            if sim=='kem1pert1b':
                rowl='nosst'
            else:
                rowl='nosit'
        else:
            frootc =  basepath + bcasename + sim + subdir + bcasename + sim + '_' + field + '_'
            frootp = basepath + bcasenamep + sim + subdir + bcasenamep + sim + '_' + field + '_'
            rowl=sim


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
            ax = ax6[cidx][ridx] # swapped row and col index positions in subplot
            
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
                if threed:
                    getfld=ncfield
                else:
                    getfld=field.upper()
                    
                fldcsea = cnc.getNCvar(fnamec,getfld,timesel=timesel,
                                               seas=sea)*conv
                fldpsea = cnc.getNCvar(fnamep,getfld,timesel=timesel,
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
                    
            tstat[cidx,:,:],pval[cidx,:,:] = sp.stats.ttest_ind(fldpsea,fldcsea,axis=0)
            fldcallseas[cidx,:,:] = np.mean(fldcsea,axis=0)
            fldpallseas[cidx,:,:] = np.mean(fldpsea,axis=0)

            if pct:
                plotfld = (fldpallseas[cidx,:,:]-fldcallseas[cidx,:,:]) / fldcallseas[cidx,:,:] *100
                cminm=cminmpct
                cmaxm=cmaxmpct
            else:
                plotfld = fldpallseas[cidx,:,:] - fldcallseas[cidx,:,:]

            pparams = dict(cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',axis=ax,latlim=latlim)
            
            if seasvert: # zonal mean with height
                pparams['suppcb'] = True
                pparams['levlim'] = levlim
                pparams['screen'] = screen
                pparams['addcontlines'] = True
                if ridx!=0: # if not the first column, suppress y labels
                    pparams['suppylab'] = True
                
                    
                pc = cplt.vert_plot(plotfld,lev,lat,**pparams)
                if sigoff==0:
                     cplt.addtsig(ax,pval[cidx,...],lat,lev/100.,type=sigtype) # @@ dims?

                if ridx==lastcol:
                    # put season label on right side.
                    ax.yaxis.set_label_position("right")
                    ax.set_ylabel(sea)
                    
            else: # maps
                pparams['suppcb'] = 1
                bm,pc = cplt.kemmap(plotfld,lat,lon,**pparams)#@@

                if sigoff==0:
                    cplt.addtsigm(bm,pval[cidx,:,:],lat,lon,type=sigtype)

                if ridx==0: # when col index is 0, set season
                    ax.set_ylabel(sea)

            if cidx==0: # when row index is 0, set simulation
                ax.set_title(rowl)     

            cidx = cidx+1

    cbar_ax = fig6.add_axes([.91,.25, .02,.5])
    fig6.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    plt.suptitle(fieldstr)
    if printtofile:
        if sigoff==0:
            sigstr='sig' + sigtype
        else:
            sigstr=''

        if latlim!= None:
            latstr=str(latlim)
        else:
            latstr=''

        if seasvert:
            if screen:
                style = 'screen'
            else:
                style = str(latlim) + 'N' + str(levlim) + 'hPa'
                
            if pct:
                fig6.savefig(fieldstr + 'pctdiff' + sigstr + '_enssubplot' + obsstr + ctstr +
                             '_seas_' + style + '2.' + suff)
            else:
                fig6.savefig(fieldstr + 'diff' + sigstr + '_enssubplot' + obsstr + ctstr +
                             '_seas_' + style + '2.' + suff)
        else: # maps
            if pct: # version 2 has new season order, new filename/key org
                fig6.savefig(fieldstr + 'pctdiff' + sigstr + '_enssubplot' + obsstr + ctstr +
                             '_seas_nh' + latstr + '2.' + suff)
            else:
                fig6.savefig(fieldstr + 'diff' + sigstr + '_enssubplot' + obsstr + ctstr + '_seas_nh'
                             + latstr + '2.' + suff)



if plotzonmean==1 or plotseacyc==1 or pattcorrwithtime==1 or plotregmean==1:
    # get data for either zonal mean or sea cycle figures
    if plotzonmean==1 or pattcorrwithtime==1 or plotregmean==1:
        # seasons is defined above
        corrlim = 45
    elif plotseacyc==1:
        
        if field in (fluxes,'fsg','turb','net'):
            latlim=40
        else:
            latlim = 70 # for area averaging
        seasons = con.get_mon()

    if sia==1:
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

        if sim=='kemhad' or sim=='kemnsidc': 
            frootc = basepath + sim + 'ctl' + subdir + sim + 'ctl' + '_' 
            frootp = basepath + sim + 'pert' + subdir + sim + 'pert' + '_' 
        elif sim in ('kem1pert1b','kem1pert3'):
            frootc = basepath + 'kemctl1' + subdir + 'kemctl1' + '_' 
            frootp = basepath + sim + subdir + sim + '_'
        else:
            frootc =  basepath + bcasename + sim + subdir + bcasename + sim + '_'
            frootp = basepath + bcasenamep + sim + subdir + bcasenamep + sim + '_'

        if threed==0:
            #if sim=='kemnsidc':
            #    print '@@ temporary fix while nsidc runs go'
            #    fnamec = con.getBCfilenames('sicn',sim+'ctl')
            #    fnamep = con.getBCfilenames('sicn',sim+'pert')
            #else:
            fnamec = frootc + field + '_' + timstr + '_ts.nc'
            fnamep = frootp + field + '_' + timstrp + '_ts.nc'
        else:
            if nonstandardlev:               
                fnamec = frootc + field + '_' + '001-061_ts.nc'
                fnamec2 = frootc + field + '_' + '062-121_ts.nc'
                fnamep = frootp + field + '_' + '001-061_ts.nc'
                fnamep2 = frootp + field + '_' + '062-121_ts.nc'
            else:
                fnamec = frootc + field + '_' + timstr + '_ts.nc'
                fnamep = frootp + field + '_' + timstrp + '_ts.nc'
        
        for sii,sea in enumerate(seasons):

            if plotzonmean==1 or pattcorrwithtime==1 or plotregmean==1:
                ncparams = {'seas': sea}
            elif plotseacyc==1:
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

                fldczm = cnc.getNCvar(fnamec,fielda.upper(),timesel=timesel,
                                               **ncparams)*conv + cnc.getNCvar(fnamecb,fieldb.upper(),
                                               timesel=timesel,**ncparams)*conv
                fldpzm = cnc.getNCvar(fnamep,fielda.upper(),timesel=timesel,
                                               **ncparams)*conv + cnc.getNCvar(fnamepb,fieldb.upper(),
                                               timesel=timesel,**ncparams)*conv
                if field=='net':
                    #print 'getting flg for net'
                    fieldb='flg'
                    conv=-1
                    fnamecb = frootc + fieldb + '_' + timstr + '_ts.nc'
                    fnamepb = frootp + fieldb + '_' + timstrp + '_ts.nc'
                    fldczm = fldczm + cnc.getNCvar(fnamecb,fieldb.upper(),
                                                   timesel=timesel,**ncparams)*conv
                    fldpzm = fldpzm + cnc.getNCvar(fnamepb,fieldb.upper(),
                                                   timesel=timesel,**ncparams)*conv
                    #print fldpzm-fldczm
                    field='net'
                    conv=1
                else:
                    field='turb'

                
            else:
                if threed==0:
                    fldczm = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                                   **ncparams)*conv
                    fldpzm = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                                   **ncparams)*conv

                else:
                    #print '@@ fix to use level NC files'
                    
                    if nonstandardlev:
                        ncparams['levsel'] = level
                        fldczm = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                            cnc.getNCvar(fnamec2,ncfield,**ncparams)*conv,
                                            axis=0)
                        fldpzm = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                            cnc.getNCvar(fnamep2,ncfield,**ncparams)*conv,
                                            axis=0)
                    else:
                        fldczm = cnc.getNCvar(fnamec,ncfield,timesel=timesel,
                                              **ncparams)*conv
                        fldpzm = cnc.getNCvar(fnamep,ncfield,timesel=timesel,
                                              **ncparams)*conv
                        
                if sia==1:
                    #areas = cutl.calc_cellareas(lat,lon,repeat=fldczm.shape)
                    #print areas.shape
                    #fldczm = fldczm*areas
                    #fldpzm = fldpzm*areas
                    fldczm = cutl.calc_seaicearea(fldczm,lat,lon)
                    fldpzm = cutl.calc_seaicearea(fldpzm,lat,lon)
                    

            if field in (fluxes,'fsg','turb','net'):
                # mask out regions that are not ice in the control (as P&M 2014 JClim)
                sicnc = cnc.getNCvar(frootc + 'sicn_' + timstr + '_ts.nc','SICN',timesel=timesel,**ncparams)
                
                fldczm = ma.masked_where(sicnc<.10,fldczm)
                fldpzm = ma.masked_where(sicnc<.10,fldpzm)
                
            if plotzonmean==1:
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

                seapcorrdict[sea] = pcorr

            elif plotseacyc==1:
                if withlat==1:
                    # leave the latitude dimension intact
                    # dims are time x lat x lon to start
                    # take zonal mean
                    fldczm = np.mean(fldczm,axis=2) # take zonal mean
                    fldpzm = np.mean(fldpzm,axis=2)
                else:
                    #print field
                    if sia==1:
                        #calc total area instead of average
                        fldczm = np.sum(np.sum(fldczm[:,lat>0,:],axis=2),axis=1)
                        fldpzm = np.sum(np.sum(fldpzm[:,lat>0,:],axis=2),axis=1)
                    else:
                        # consider masking out land for sfc fluxes...?
                        fldczm = cutl.polar_mean_areawgted3d(fldczm,lat,lon,latlim=latlim)
                        fldpzm = cutl.polar_mean_areawgted3d(fldpzm,lat,lon,latlim=latlim)

            elif plotregmean==1:

                lons,lats = np.meshgrid(lon,lat)
                
                ntime = fldczm.shape[0]
                
                reglatsbool = np.logical_and(lat>latlims[0],lat<latlims[1])
                reglonsbool = np.logical_and(lon>lonlims[0],lon<lonlims[1])
                regmask = np.logical_or(
                    np.logical_or(lats<latlims[0],lats>latlims[1]), 
                    np.logical_or(lons<lonlims[0],lons>lonlims[1]))
                regmaskt = np.tile(regmask,(ntime,1,1)) # tiled regional mask

                areas = cutl.calc_cellareas(lat,lon)
                areasm = ma.masked_where(regmask,areas)
                weightsm = areasm / np.sum(np.sum(areasm,axis=1),axis=0) # weights masked
                weightsmt = np.tile(weightsm,(ntime,1,1)) # weights masked tiled
                
                # regional subset: control
                ## tmp = fldczm[:,reglatsbool,:] # this swaps the lat and time dims for some reason
                ## tmp = np.transpose(tmp,(1,0,2))
                ## tmp = tmp[:,:,reglonsbool]
                ## tmp = np.transpose(tmp,(1,0,2))
                ## tmpreg = np.sum(np.sum(tmp*weightsmt,axis=2),axis=1)
                ## fldczm = tmpreg # should be timeseries of regional mean

                tmp = ma.masked_where(regmaskt,fldczm)
                tmpreg = np.sum(np.sum(tmp*weightsmt,axis=2),axis=1)
                fldczm = tmpreg # should be timeseries of regional mean

                # regional subset: pert
                ## tmp = fldpzm[:,reglatsbool,:] # this swaps the lat and time dims for some reason
                ## tmp = np.transpose(tmp,(1,0,2))
                ## tmp = tmp[:,:,reglonsbool]
                ## tmp = np.transpose(tmp,(1,0,2))
                ## tmpreg = np.sum(np.sum(tmp*weightsmt,axis=2),axis=1)
                ## fldpzm = tmpreg # should be timeseries of regional mean

                tmp = ma.masked_where(regmaskt,fldpzm)
                tmpreg = np.sum(np.sum(tmp*weightsmt,axis=2),axis=1)
                fldpzm = tmpreg # should be timeseries of regional mean
 

            seafldcstddict[sea] = np.std(fldczm,axis=0)
            seafldpstddict[sea] = np.std(fldpzm,axis=0)
            ttmp,pvtmp = sp.stats.ttest_ind(fldpzm,fldczm,axis=0)
            if plotregmean==1:
                # double-check the scale setting
                ci = sp.stats.t.interval(1-siglevel,len(fldpzm)-1,loc=np.mean(fldpzm,axis=0)-np.mean(fldczm,axis=0),
                                         scale=np.std(fldpzm,axis=0)/np.sqrt(len(fldpzm)))
                seacidict[sea] = ci
                #print ci # @@@
                
            seatstatdict[sea] = ttmp
            seapvaldict[sea] = pvtmp
            seafldcdict[sea] =  np.mean(fldczm,axis=0) # time mean
            seafldpdict[sea] =  np.mean(fldpzm,axis=0)
            seadiffdict[sea] = np.mean(fldpzm,axis=0)- np.mean(fldczm,axis=0)
            seadmaskdict[sea] = ma.masked_where(pvtmp>siglevel,seadiffdict[sea])

            # end loop through seasons

        fldcstddict[sim] = seafldcstddict
        fldpstddict[sim] = seafldpstddict
        if plotregmean==1:
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
    if sia==1:
        field = 'sia' # put back after getting the data

fontP = fm.FontProperties()
fontP.set_size('small')

darkolivegreen1 = np.array([202, 255, 112])/255. # terrible
darkolivegreen3 = np.array([162, 205, 90])/255.
darkseagreen = np.array([143, 188, 143])/255.
darkseagreen4 = np.array([105, 139, 105])/255.
dodgerblue = np.array([30, 144, 255])/255. 
orangered4 = np.array([139, 37, 0])/255.
lightsteelblue3 = np.array([162, 181, 205])/255. # more grey looking
lightsteelblue4 = np.array([110, 123, 139])/255. # more grey looking
steelblue3 = np.array([79, 148, 205])/255.  # more blue looking
steelblue4 = np.array([54, 100, 139])/255.  # more blue looking

colors = darkseagreen,darkseagreen4,lightsteelblue3,lightsteelblue4,steelblue4,dodgerblue,orangered4,darkolivegreen3
## colordict = {'kemctl1r1': darkseagreen, 'kemctl1r2': darkseagreen4, 'kemctl1r3': lightsteelblue3,
##              'kemctl1r4': lightsteelblue4, 'kemctl1r5': steelblue4, 'kemctl1ens': dodgerblue,
##              'kemctl1': orangered4, 'kemhadctl': darkolivegreen3 }

# try to match the warmness of color to the amount of sea ice loss in the simulation
#  Warmer means more loss in melt season. Ensemble mean, PERT2, and obs (hadisst) have special colors
"""colordict = {'kemctl1r1': ccm.get_linecolor('chocolate4'), #'firebrick'),
             'kemctl1r4': ccm.get_linecolor('firebrick'), #'firebrick1'),
             'kemctl1r3': ccm.get_linecolor('firebrick1'), #'yelloworange'),#'chocolate1'),
             'kemctl1r5': ccm.get_linecolor('yelloworange'), #'darkyellow') #'skyblue'), #yelloworange'),
             'kemctl1r2': ccm.get_linecolor('darkyellow'), #'steelblue3'), #'darkgoldenrod1'),
             'kemctl1ens': ccm.get_linecolor('magenta'),
             'kemctl1': ccm.get_linecolor('dodgerblue'), #'mediumpurple1'), #darkyellow'),
             'kemhadctl': ccm.get_linecolor('darkolivegreen3')} """

"""colordict = {'kemctl1r1': ccm.get_linecolor('red1'), #'firebrick'),
             'kemctl1r4': ccm.get_linecolor('red2'), #'firebrick1'),
             'kemctl1r3': ccm.get_linecolor('red3'), #'yelloworange'),#'chocolate1'),
             'kemctl1r5': ccm.get_linecolor('red4'), #'darkyellow') #'skyblue'), #yelloworange'),
             'kemctl1r2': ccm.get_linecolor('red5'), #'steelblue3'), #'darkgoldenrod1'),
             'kemctl1ens': ccm.get_linecolor('magenta'),
             'kemctl1': ccm.get_linecolor('dodgerblue'), #'mediumpurple1'), #darkyellow'),
             'kemhadctl': ccm.get_linecolor('darkolivegreen3')}"""

"""colordict = {'r1': ccm.get_linecolor('warm1'), #'firebrick'),
             'r4': ccm.get_linecolor('warm2'), #'firebrick1'),
             'r3': ccm.get_linecolor('warm3'), #'yelloworange'),#'chocolate1'),
             'r5': ccm.get_linecolor('warm4'), #'darkyellow') #'skyblue'), #yelloworange'),
             'r2': ccm.get_linecolor('warm5'), #'steelblue3'), #'darkgoldenrod1'),
             'ens': ccm.get_linecolor('magenta'),
             '': ccm.get_linecolor('mediumblue'), #'mediumpurple1'), #darkyellow'),
             'kemhad': ccm.get_linecolor('deepskyblue'),
             'kemnsidc': ccm.get_linecolor('steelblue4'),
             'r4ct': ccm.get_linecolor('darkseagreen4'),
             'kem1pert1b': ccm.get_linecolor('darkolivegreen3'),
             'kem1pert3': ccm.get_linecolor('darkolivegreen1')}"""

colordict = ccm.get_colordict()


# now make the figures ==============================================
if plotzonmean:

    seal = list(seasons)
    
    flddiffdf = pd.DataFrame(flddiffdict)
    fldcdf = pd.DataFrame(fldcdict)
    fldpdf = pd.DataFrame(fldpdict)
    fldmaskdf = pd.DataFrame(flddmaskdict)
    fldcstddf = pd.DataFrame(fldcstddict)
    fldpstddf = pd.DataFrame(fldpstddict)   
    
    tmpcdf = fldcdf.loc[seal]
    tmppdf = fldpdf.loc[seal]

    cestd = np.zeros((4,len(lat)))
    pestd = np.zeros((4,len(lat)))
    cemin = np.zeros((4,len(lat))) # ctl ens min/max
    cemax = np.zeros((4,len(lat)))
    pemin = np.zeros((4,len(lat))) # pert ens min/max
    pemax = np.zeros((4,len(lat)))
    demin = np.zeros((4,len(lat))) # difference min/max
    demax = np.zeros((4,len(lat)))
    # have to do these loops b/c np.array(tmpcdf) comes out as an object of
    # size (seasons x sims) with inner arrays of length lat, and can't take std of that
    # also construct a min and max line of ens, for fill_between spread
    for sii,sea in enumerate(seasons):
        tmpc = np.zeros((5,len(lat)))
        tmpp = np.zeros((5,len(lat)))
        for simii,sim in enumerate(sims[0:5]):
            # here we just accumulate all ens sims into one matrix
            tmpc[simii,:] = tmpcdf[sim][sea].values
            tmpp[simii,:] = tmppdf[sim][sea].values
        # now do calcs: standard dev, min, max
        cestd[sii,:] = np.std(tmpc,axis=0)
        pestd[sii,:] = np.std(tmpp,axis=0)
        cemax[sii,:] = np.max(tmpc,axis=0)
        cemin[sii,:] = np.min(tmpc,axis=0)
        pemax[sii,:] = np.max(tmpp,axis=0)
        pemin[sii,:] = np.min(tmpp,axis=0)
        demax[sii,:] = np.max(tmpp-tmpc,axis=0)
        demin[sii,:] = np.min(tmpp-tmpc,axis=0)
        
    # climo
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    
    for sii,sea in enumerate(seasons):
        ax=axs[sii]
        
        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): #== '' or skey == 'ens' or skey=='kemhad':
                ax.plot(lat,fldcdf[skey][sea],color=colordict[skey],linewidth=3)
            else:
                ax.plot(lat,fldcdf[skey][sea],color=colordict[skey],linewidth=2)

        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                ax.plot(lat,fldpdf[skey][sea],color=colordict[skey],linewidth=3,linestyle='--')
            else:
                ax.plot(lat,fldpdf[skey][sea],color=colordict[skey],linewidth=2,linestyle='--')
            
        ax.set_xlim(-90,90)
        ax.set_title(sea)
        
    ax.set_xlabel('lat')
    if printtofile: # version 2 has new colors and seasons and sims are reordered
        fig.savefig(fieldstr + '_ens_meanBC' + obsstr + ctstr + '_allseassp_zonmean2.pdf')

    # standard deviation
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for sii,sea in enumerate(seasons):
        ax=axs[sii]
        
        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                ax.plot(lat,fldcstddf[skey][sea],color=colordict[skey],linewidth=3)
            else:
                ax.plot(lat,fldcstddf[skey][sea],color=colordict[skey],linewidth=2)
        ax.plot(lat,cestd[sii,...],'k',linewidth=2)

        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                ax.plot(lat,fldpstddf[skey][sea],color=colordict[skey],linewidth=3,linestyle='--')
            else:
                ax.plot(lat,fldpstddf[skey][sea],color=colordict[skey],linewidth=2,linestyle='--')
        ax.plot(lat,pestd[sii,...],'k',linewidth=2,linestyle='--')

        ax.set_xlim(-90,90)
        ax.set_title(sea)
        
    ax.set_xlabel('lat')
    if printtofile:
        fig.savefig(fieldstr + 'STD_ens_meanBC' + obsstr + ctstr + '_allseassp_zonmean2.pdf')

    # differences
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for sii,sea in enumerate(seasons):
        ax=axs[sii]
        
        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                ax.plot(lat,flddiffdf[skey][sea],color=colordict[skey],linewidth=3)
            else:
                ax.plot(lat,flddiffdf[skey][sea],color=colordict[skey],linewidth=2)

        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                ax.plot(lat,fldmaskdf[skey][sea],color=colordict[skey],linestyle='none',marker='s')
            else:
                ax.plot(lat,fldmaskdf[skey][sea],color=colordict[skey],linestyle='none',marker='s')

        ax.set_xlim(0,90)
        ax.set_title(field + ': ' + sea)
        
    ax.set_xlabel('lat')
    ax.legend(sims,'upper left', prop=fontP,ncol=2)
    if printtofile:
        fig.savefig(fieldstr + 'diff_ens_meanBC' + obsstr + ctstr + '_allseassp_zonmean_nh2.pdf')

    # differences SHADED
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for sii,sea in enumerate(seasons):
        ax=axs[sii]

        ax.fill_between(lat,demin[sii,...],demax[sii,...],facecolor='0.7',alpha=0.2)
        for skey in sims[5:]:# actually don't skip # skip the r4ct sim here
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                ax.plot(lat,flddiffdf[skey][sea],color=colordict[skey],linewidth=3)
            else:
                ax.plot(lat,flddiffdf[skey][sea],color=colordict[skey],linewidth=2)

        for skey in sims[5:]: #-1]:
            ax.plot(lat,fldmaskdf[skey][sea],color=colordict[skey],linestyle='none',marker='s')

        ax.set_xlim(0,90)
        ax.set_title(field + ': ' + sea)
        
    ax.set_xlabel('lat')
    ax.legend(sims[5:],'upper left', prop=fontP,ncol=2)
    if printtofile:
        fig.savefig(fieldstr + 'diff_ens_meanBC' + obsstr + ctstr + '_allseassp_zonmean_nh2shade.pdf')

    # standard deviation differences
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for sii,sea in enumerate(seasons):
        ax=axs[sii]
        
        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                ax.plot(lat,fldpstddf[skey][sea]-fldcstddf[skey][sea],color=colordict[skey],linewidth=3)
            else:
                ax.plot(lat,fldpstddf[skey][sea]-fldcstddf[skey][sea],color=colordict[skey],linewidth=2)
        ax.plot(lat,pestd[sii,...]-cestd[sii,...],'k',linewidth=2)

        ax.set_xlim(0,90)
        ax.set_title(field + ': ' + sea)
        
    ax.set_xlabel('lat')
    
    if printtofile:
        fig.savefig(fieldstr + 'STDdiff_ens_meanBC' + obsstr + ctstr + '_allseassp_zonmean_nh2.pdf')

if plotseacyc:

    # http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing
    
    ## Operation 	Syntax 	Result
    ## Select column 	df[col] 	Series
    ## Select row by label 	df.loc[label] 	Series
    ## Select row by integer location 	df.iloc[loc] 	Series
    ## Slice rows 	df[5:10] 	DataFrame
    ## Select rows by boolean vector 	df[bool_vec] 	DataFrame

    # To plot from DataFrame: @@
    #      import pandas as pd
    #   create DataFrame from dictionary:
    #      df = pd.DataFrame(flddiffdict)
    #   select all months for one simulation for plotting:
    #      #this gets the 12 months for the given simulation, but months are out of order
    #      seacycle = df[simkey]
    #   change months tuple to list for indexing df and getting in correct order
    #      mol = list(months)
    #      plt.plot(seacycle[mol]) # plots in correct order

    # @@ this should help with getting a range of ensemble runs (max,min) for shading

    mol = list(months) # use this list of strings for indexing the dataframe

    fontP = fm.FontProperties()
    fontP.set_size('small')

    # Now plot seasonal cycle
    moidxs=np.arange(1,13)

    if withlat:

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
            for simii,sim in enumerate(sims[0:5]):# accumulate all ens sims together
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
        squatfs=(6,3) # figsize
        squatterfs=(6,2.5) #figsize

        flddiffdf = pd.DataFrame(flddiffdict)
        fldcdf = pd.DataFrame(fldcdict)
        fldpdf = pd.DataFrame(fldpdict)
        fldmaskdf = pd.DataFrame(flddmaskdict)
        fldcstddf = pd.DataFrame(fldcstddict)
        fldpstddf = pd.DataFrame(fldpstddict)
        
        # climo
        fig,axs = plt.subplots()
        if squatseacyc:
            fsuff='short'
            fig.set_size_inches(squatfs)
        elif squatterseacyc:
            fsuff='shorter'
            fig.set_size_inches(squatterfs)
                
        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                axs.plot(moidxs,fldcdf[skey][mol],color=colordict[skey],linewidth=3)
            else:
                axs.plot(moidxs,fldcdf[skey][mol],color=colordict[skey],linewidth=2)

        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                axs.plot(moidxs,fldpdf[skey][mol],color=colordict[skey],linewidth=3,linestyle='--')
            else:
                axs.plot(moidxs,fldpdf[skey][mol],color=colordict[skey],linewidth=2,linestyle='--')

        plt.legend(sims,leglocs[0], prop=fontP,ncol=2)
        plt.xlim((1,12))
        plt.gca().set_xticks(range(1,13))
        plt.gca().set_xticklabels(months)
        #plt.xlabel('Month')
        plt.ylabel(fieldstr)
        plt.title('Climos')

        if printtofile: 
            fig.savefig(fieldstr + '_ens_meanBC' + obsstr + ctstr + '_seacyc_pol' + str(latlim) + 'N3' + fsuff + '.pdf')


        # differences
        fig,axs = plt.subplots()
        if squatseacyc:
            fsuff='short'
            fig.set_size_inches(squatfs)
        elif squatterseacyc:
            fsuff='shorter'
            fig.set_size_inches(squatterfs)
            
        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                axs.plot(moidxs,flddiffdf[skey][mol],color=colordict[skey],linewidth=3)
            else:
                axs.plot(moidxs,flddiffdf[skey][mol],color=colordict[skey],linewidth=2)
        for skey in sims:
            axs.plot(moidxs,fldmaskdf[skey][mol],linestyle='none',color=colordict[skey],marker='s')

        plt.legend(sims,leglocs[1], prop=fontP,ncol=2)
        plt.xlim((1,12))
        plt.ylim(seacycylim)
        plt.gca().set_xticks(range(1,13))
        plt.gca().set_xticklabels(months)
        #plt.xlabel('Month')
        plt.ylabel(fieldstr)
        plt.title('Anomalies')

        if printtofile: # version 2 loops through sims in order of melt
            fig.savefig(fieldstr + 'diff_ens_meanBC' + obsstr + ctstr + '_seacyc_pol' + str(latlim) + 'N3' + fsuff + '.pdf')

        # calc stddev over ensemble, and min/max for shading
        tmpcdf = fldcdf.loc[mol]
        tmppdf = fldpdf.loc[mol]
        ce = np.array(tmpcdf.loc[:,sims[0:5]]) # gives array of month x simulation in correct order
        pe = np.array(tmppdf.loc[:,sims[0:5]])
        #ce = np.array(tmpcdf.loc[:,sims[ensmems]]) # index as is (np array) doesn't work
        #pe = np.array(tmppdf.loc[:,sims[ensmems]])
        cestd = np.std(ce,axis=1) # std dev over ensemble, 1st 5 simulations
        pestd = np.std(pe,axis=1)
        cemax = np.max(ce,axis=1) # ens ctl max/min
        cemin = np.min(ce,axis=1)
        pemax = np.max(pe,axis=1) # ens pert max/min
        pemin = np.min(pe,axis=1)
        demax = np.max(pe-ce,axis=1) # ens diff max/min
        demin = np.min(pe-ce,axis=1)

        # differences SHADED
        fig,axs = plt.subplots()
        if squatseacyc:
            fsuff='short'
            fig.set_size_inches(squatfs)
        elif squatterseacyc:
            fsuff='shorter'
            fig.set_size_inches(squatterfs)
            
        for skey in sims[5:]:  #-1]:
            axs.fill_between(moidxs,demin,demax,facecolor='0.7',alpha=0.2)
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                axs.plot(moidxs,flddiffdf[skey][mol],color=colordict[skey],linewidth=3)
            else:
                axs.plot(moidxs,flddiffdf[skey][mol],color=colordict[skey],linewidth=2)
        for skey in sims[5:]: #-1]:
            axs.plot(moidxs,fldmaskdf[skey][mol],linestyle='none',color=colordict[skey],marker='s')

        plt.legend(sims[5:],leglocs[1], prop=fontP,ncol=2)
        plt.xlim((1,12))
        plt.ylim(seacycylim)
        plt.gca().set_xticks(range(1,13))
        plt.gca().set_xticklabels(months)
        #plt.xlabel('Month')
        plt.ylabel(fieldstr)
        plt.title('Anomalies')

        if printtofile: # version 2 loops through sims in order of melt
            fig.savefig(fieldstr + 'diff_ens_meanBC' + obsstr + ctstr + '_seacyc_pol' + str(latlim) + 'N3shade' + fsuff + '.pdf')


        # Standard deviation climos
        fig,axs = plt.subplots()
        if squatseacyc:
            fsuff='short'
            fig.set_size_inches(squatfs)
        elif squatterseacyc:
            fsuff='shorter'
            fig.set_size_inches(squatterfs)
            
        for skey in sims:
            
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                axs.plot(moidxs,fldcstddf[skey][mol],color=colordict[skey],linewidth=3)
            else:
                axs.plot(moidxs,fldcstddf[skey][mol],color=colordict[skey],linewidth=2)

        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): # == '' or skey == 'ens' or skey=='kemhad':
                axs.plot(moidxs,fldpstddf[skey][mol],color=colordict[skey],linestyle='--',linewidth=3)
            else:
                axs.plot(moidxs,fldpstddf[skey][mol],color=colordict[skey],linestyle='--',linewidth=2)

        axs.plot(range(1,13),cestd,color='k',linewidth=3)
        axs.plot(range(1,13),pestd,color='k',linewidth=3,linestyle='--')

        plt.legend(sims,leglocs[2], prop=fontP,ncol=2)
        plt.xlim((1,12))
        plt.gca().set_xticks(range(1,13))
        plt.gca().set_xticklabels(months)
        #plt.xlabel('Month')
        plt.ylabel(fieldstr)
        plt.title('Sigma')

        if printtofile:
            fig.savefig(fieldstr + 'STD_ens_meanBC' + obsstr + ctstr + '_seacyc_pol' + str(latlim) + 'N3' + fsuff + '.pdf')


        # Difference in standard deviation
        fig,axs = plt.subplots()
        if squatseacyc:
            fsuff='short'
            fig.set_size_inches(squatfs)
        elif squatterseacyc:
            fsuff='shorter'
            fig.set_size_inches(squatterfs)
            
        for skey in sims:
            if skey in ('','ens','kemhad','kemnsidc'): #  == '' or skey == 'ens' or skey=='kemhad':
                axs.plot(moidxs,fldpstddf[skey][mol]-fldcstddf[skey][mol],color=colordict[skey],linewidth=3)
            else:
                axs.plot(moidxs,fldpstddf[skey][mol]-fldcstddf[skey][mol],color=colordict[skey],linewidth=2)
        axs.plot(moidxs,pestd-cestd,color='k',linewidth=3)

        plt.legend(sims,leglocs[3], prop=fontP,ncol=2)
        plt.xlim((1,12))
        plt.gca().set_xticks(range(1,13))
        plt.gca().set_xticklabels(months)
        #plt.xlabel('Month')
        plt.ylabel(fieldstr)
        plt.title('Sigma anomalies')

        if printtofile:
            fig.savefig(fieldstr + 'STDdiff_ens_meanBC' + obsstr + ctstr + '_seacyc_pol' + str(latlim) + 'N3' + fsuff + '.pdf')


if plotregmean==1:

    flddiffdf = pd.DataFrame(flddiffdict)
    fldcdf = pd.DataFrame(fldcdict)
    fldpdf = pd.DataFrame(fldpdict)
    fldmaskdf = pd.DataFrame(flddmaskdict)
    fldcstddf = pd.DataFrame(fldcstddict)
    fldpstddf = pd.DataFrame(fldpstddict)
    cidf = pd.DataFrame(cidict)

    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(10,8)
    print '@@ add confidence intervals!'
    for sii,sea in enumerate(seasons):

        ax=axs[sii]
        if sii==0:
            ax.set_title(fieldstr + ' ' + region)
            
        for skeyii,skey in enumerate(sims):
             
            val=flddiffdf[skey][sea]
            ci=cidf[skey][sea]
            #print ci
            
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
        fig.savefig(fieldstr + 'diffCI_ens_meanBC' + obsstr + ctstr + '_allseassp_' + region + '.pdf')

if pattcorrwithtime==1:


    fldpcorrdf = pd.DataFrame(fldpcorrdict)
    tmppcdf = fldpcorrdf.loc[seal] # put seasons in order


    pcemax = np.zeros((4,len(tmppcdf['r1']['DJF'])))# get any element to get time
    pcemin = np.zeros((4,len(tmppcdf['r1']['DJF'])))
    # have to do these loops b/c np.array(tmppcdf) comes out as an object of
    # size (seasons x sims) with inner arrays of length time, and can't take max/min of that
    for sii,sea in enumerate(seasons):
        tmppc = np.zeros( (5,len(tmppcdf['r1']['DJF'])) )
        for simii,sim in enumerate(sims[0:5]):
            # here we just accumulate all ens sims into one matrix
            tmpsorted = tmppcdf[sim][sea].values
            if sea in ('JJA','MAM','SON'):
                if threed==1:
                    tmpsorted=tmpsorted[:-2]
                else:
                    tmpsorted=tmpsorted[:-1] # all timeseries have to be same length
                
            if pattcorryr:
                # sort
                tmpsorted.sort()
                
            tmppc[simii,:] = tmpsorted

        # now do calcs: min, max            
        pcemax[sii,:] = np.max(tmppc,axis=0)
        pcemin[sii,:] = np.min(tmppc,axis=0)
        
    
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
                sortfld = fldpcorrdict[sim][sea]
                sortfld.sort() # sort from smallest to largest patt corr
                ax.plot(sortfld,color=colordict[sim],linewidth=2)
            else:
                ax.plot(fldpcorrdict[sim][sea],color=colordict[sim],linewidth=2)
        ax.set_ylim(ylims)    
        ax.set_title(fieldstr + ': ' + sea)
        
    ax.set_xlabel('lat')
    ax.legend(sims,'lower right', prop=fontP,ncol=2)
    
    if printtofile:
        if pattcorryr:
            fig.savefig(fieldstr + 'diffpattcorryrly_ens_meanBC' + obsstr + ctstr + '_allseassp_nh.pdf')
        else:
            fig.savefig(fieldstr + 'diffpattcorr_ens_meanBC' + obsstr + ctstr + '_allseassp_nh.pdf')

    # with SHADING
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for ii,ax in enumerate(axs.flat):
        sea = seasons[ii]
        ax.fill_between(np.arange(0,pcemin.shape[1]),pcemin[ii,...],pcemax[ii,...],facecolor='0.7',alpha=0.2)
        for simii,sim in enumerate(sims[5:]):
            
            if pattcorryr:
                sortfld = fldpcorrdict[sim][sea]
                sortfld.sort() # sort from smallest to largest patt corr
                ax.plot(sortfld,color=colordict[sim],linewidth=2)
            else:
                ax.plot(fldpcorrdict[sim][sea],color=colordict[sim],linewidth=2)
        ax.set_ylim(ylims)    
        ax.set_title(fieldstr + ': ' + sea)
        
    ax.set_xlabel('lat')
    ax.legend(sims[5:],'lower right', prop=fontP,ncol=2)
    
    if printtofile:
        if pattcorryr:
            fig.savefig(fieldstr + 'diffpattcorryrly_ens_meanBC' + obsstr + ctstr + '_allseassp_nhshade.pdf')
        else:
            fig.savefig(fieldstr + 'diffpattcorr_ens_meanBC' + obsstr + ctstr + '_allseassp_nhshade.pdf')

    
if testhadisst:
    # here want to check which ensemble run is most similar to hadisst in terms of SICN

    # get HadISST runs
    fnamec2 = basepath + casename2 + subdir + casename2 + '_' + field + '_' + timstr2 + '_ts.nc'
    fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstr2 + '_ts.nc'

    fldc2 = cnc.getNCvar(fnamec2,field.upper(),timesel=timesel)*conv
    fldp2 = cnc.getNCvar(fnamep2,field.upper(),timesel=timesel)*conv
            
    corrlim=40 # the latitude north of which to consider the pattern correlation
    lmask = con.get_t63landmask()
    
    fldcdict = {}; fldpdict = {}
    diffdict = {}
    pcorrdict = {}

    # for weighting the pattern corr by area
    areas = cutl.calc_cellareas(lat,lon)
    areas = areas[lat>corrlim,:]
    areas = ma.masked_where(lmask[lat>corrlim,:]==-1,areas)
    weights = areas / np.sum(np.sum(areas,axis=1),axis=0)
    
    for ridx,sim in enumerate(sims[0:7]):
        
        frootc =  basepath + bcasename + sim + subdir + bcasename + sim + '_' + field + '_'
        frootp = basepath + bcasenamep + sim + subdir + bcasenamep + sim + '_' + field + '_'

        print sim
        fnamec = frootc + timstr + '_ts.nc'
        fnamep = frootp + timstrp + '_ts.nc'

        fldcdict[sim] = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv
        fldpdict[sim] = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv

        #cplt.map_allmonths((fldpdict[sim]-fldcdict[sim]) - (fldp2-fldc2),
        #                   lat,lon,cmin=-.1,cmax=.1,cmap='red2blue_w20',title=sim,type='nh')

        # root mean square error(?) @@not quite
        # want to take the difference b/w the ensrun and hadisst anomalies at each grid point,
        # squaring and square rooting will get the distance away ensrun is from hadisst
        # then take the polar mean of the distances
        #yrfld = np.zeros(12)
        tmpfld = np.sqrt(cutl.polar_mean_areawgted3d(np.square((fldpdict[sim][0:12,...]-
                                                                fldcdict[sim][0:12,...]) -
                                                               (fldp2[0:12,...]-fldc2[0:12,...])),
                                                     lat,lon,latlim=40))
                    
        diffdict[sim] = tmpfld

        # do correlation here:
        pcorr = np.zeros((12))
        for moidx in range(0,12):
            ensmem = fldpdict[sim][moidx,lat>corrlim,...] - fldcdict[sim][moidx,lat>corrlim,...]
            obsbc = fldp2[moidx,lat>corrlim,...] - fldc2[moidx,lat>corrlim,...] # @@where do these flds get set?

            ensmem = ma.masked_where(lmask[lat>corrlim,:]==-1,ensmem) # mask out land
            obsbc = ma.masked_where(lmask[lat>corrlim,:]==-1,obsbc) # mask out land

            # @@@ note have to use masked corrcoef!
            # @@ This is probably doing a centered pearson where I really
            # want uncentered (no central mean removed). Could calc myself from
            # http://www.stanford.edu/~maureenh/quals/html/ml/node53.html 
            tmpcorr = ma.corrcoef(ensmem.flatten()*weights.flatten(), obsbc.flatten()*weights.flatten())
            pcorr[moidx] = tmpcorr[0,1]
            
        pcorrdict[sim] = pcorr

    moidxs = range(1,13)
    fontP = fm.FontProperties()
    fontP.set_size('small')
    
    plt.figure(); # @@ could also just loop sim keys here
    plt.plot(moidxs,diffdict['r1'],color=colordict['r1'],linewidth=2)
    plt.plot(moidxs,diffdict['r2'],color=colordict['r2'],linewidth=2)
    plt.plot(moidxs,diffdict['r3'],color=colordict['r3'],linewidth=2)
    plt.plot(moidxs,diffdict['r4'],color=colordict['r4'],linewidth=2)
    plt.plot(moidxs,diffdict['r5'],color=colordict['r5'],linewidth=2)
    plt.plot(moidxs,diffdict['ens'],color=colordict['ens'],linewidth=3)
    plt.plot(moidxs,diffdict[''],color=colordict[''],linewidth=3)
    plt.legend(('r1','r2','r3','r4','r5','ens','meanBC'),prop=fontP)
    plt.xlim((1,12))
    plt.title('RMSE compared to HadISST')
    if printtofile:
        plt.savefig('sicnRMSE_ens_v_hadisst_seacycle.pdf')

    plt.figure();
    plt.plot(moidxs,pcorrdict['r1'],color=colordict['r1'],linewidth=2)
    plt.plot(moidxs,pcorrdict['r2'],color=colordict['r2'],linewidth=2)
    plt.plot(moidxs,pcorrdict['r3'],color=colordict['r3'],linewidth=2)
    plt.plot(moidxs,pcorrdict['r4'],color=colordict['r4'],linewidth=2)
    plt.plot(moidxs,pcorrdict['r5'],color=colordict['r5'],linewidth=2)
    plt.plot(moidxs,pcorrdict['ens'],color=colordict['ens'],linewidth=3)
    plt.plot(moidxs,pcorrdict[''],color=colordict[''],linewidth=3)
    plt.legend(('r1','r2','r3','r4','r5','ens','meanBC'),'upper left',prop=fontP)
    plt.xlim((1,12))
    plt.title('pattern correlation with HadISST')
    if printtofile:
        plt.savefig('sicnPatternCorr_ens_v_hadisst_seacycle.pdf')

    wgts = con.get_monweights()
    print 'r1: ' + str(np.average(diffdict['r1'],weights=wgts))
    print 'r2: ' + str(np.average(diffdict['r2'],weights=wgts))
    print 'r3: ' + str(np.average(diffdict['r3'],weights=wgts))
    print 'r4: ' + str(np.average(diffdict['r4'],weights=wgts))
    print 'r5: ' + str(np.average(diffdict['r5'],weights=wgts))
    print 'ens: ' + str(np.average(diffdict['ens'],weights=wgts))
    print 'meanBC: ' + str(np.average(diffdict[''],weights=wgts))

    print '===== PATTERN CORR ====='
    print 'r1: ' + str(np.average(pcorrdict['r1'],weights=wgts))
    print 'r2: ' + str(np.average(pcorrdict['r2'],weights=wgts))
    print 'r3: ' + str(np.average(pcorrdict['r3'],weights=wgts))
    print 'r4: ' + str(np.average(pcorrdict['r4'],weights=wgts))
    print 'r5: ' + str(np.average(pcorrdict['r5'],weights=wgts))
    print 'ens: ' + str(np.average(pcorrdict['ens'],weights=wgts))
    print 'meanBC: ' + str(np.average(pcorrdict[''],weights=wgts))

    # @@ ens/meanBC wins for both of these measures when considering all cells
    #        north of 0 (corr=0.419) or 40N (corr=0.404) (no masking, unweighted).
    #        R4 wins otherwise (corr=0.39,0.38) (no masking, unweighted)
    #  masking out land and weighting by area makes ens corr=0.422, R4 corr=0.407


"""
Pattern correlation calc?:
 From: http://www.ncl.ucar.edu/Support/talk_archives/2011/2603.html
> r_x = pattern_cor(X,O,1.0,1) ; r_x = 0.9998
> r_y = pattern_cor(Y,O,1.0,1) ; r_y= 0.99996
>
> # centered pattern correlation
> Avg_X = avg(X)
> Avg_Y = avg(Y)
> Avg_O = avg(O)
>
> Sum_X = 0.0
> Sum_y = 0.0
> Sum_dev_X = 0.0
> Sum_dev_Y = 0.0
> Sum_dev_O = 0.0
>
> do i=0, 67
> do j=0, 67
>
> Sum_x = Sum_x + (X(i,j) - Avg_X)*(O(i,j) - Avg_O)
> Sum_y = Sum_y + (Y(i,j) - Avg_Y)*(O(i,j) - Avg_O)
> Sum_dev_X = Sum_dev_X + (X(i,j) - Avg_X)*(X(i,j) - Avg_X)
> Sum_dev_Y = Sum_dev_Y +(Y(i,j) - Avg_Y)*(Y(i,j) - Avg_Y)
> Sum_dev_O = Sum_dev_O +(O(i,j) - Avg_O)*(O(i,j) - Avg_O)
>
> end do
> end do
>
> R_x = Sum_x/(sqrt(Sum_dev_X)*sqrt(Sum_dev_O)) ; R_x = 0.70
>
> R_y = Sum_y/(sqrt(Sum_dev_Y)*sqrt(Sum_dev_O)) ; R_y = 0.98
>
> 

Also try: Pearson product-moment correlation coefficient

"""
