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
import platform as platform
import constants as con      # my module
import cccmautils as cutl    # my module
import matplotlib.font_manager as fm
import copy

# while I'm still creating these modules, have to reload to get changes
cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)
#ccm = reload(ccm)
cnc = reload(cnc)

plt.close("all")
plt.ion()

printtofile=1

plotann=0    # seasonal avg map, comparing ens runs and meanBC
plotallmos=0 # monthly maps (@@ not implemented)
seasonal=0 # seasonal maps (DJF, MAM, JJA, SON)
plotzonmean=0 # plotzonmean and plotseacyc are mutually exclusive
plotseacyc=0 # plotzonmean and plotseacyc are mutually exclusive
testhadisst=1 # check which ens member most similar to hadisst
normbystd=0
addobs=1 # add mean of kemhad* runs to line plots, seasonal maps (so far@@)
latlim = None #45 # lat limit for NH plots. Set to None otherwise.

sigtype = 'cont' # significance: 'cont' or 'hatch' which is default
sigoff=0 # if 1, don't add significance
siglevel=0.05
seasons = 'DJF','MAM','JJA','SON'

model = 'CanAM4'
threed=0
sia=0 # is the requested field sea ice area


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
timstr2='001-111'


# # # ######## set Field info ###################
# gz, t, u, v, q (3D !)
# st, sic, sicn (sia), gt, pmsl, pcp, hfl, hfs, turb, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)
field = 't'
print field
timeavg = 'DJF'

# only for threed vars
#level = 30000
#level = 50000 # 500hPa
level = 70000

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
    ## cmin = -1; cmax = 1  # for anomaly plots
    ## cminm = -1.5; cmaxm = 1.5   # monthly
    
    cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
    cminn = -5; cmaxn = 5 # for norm by std
    cmap = 'blue2red_w20'

    leglocs = 'upper left', 'upper left', 'upper right', 'upper left'
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
    cminn = -1; cmaxn = 1 # for norm by std
    leglocs = 'lower left', 'upper left', 'upper center', 'upper left'
elif field == 'pcp':
    units = 'mm/day' # original: kg m-2 s-1
    
    pct=1; units = '%'
    
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
    cminm = -3
    cmaxm = 3
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

else:
    print 'No settings for ' + field


# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Volumes/MyPassport1TB/DATA/CanSISE/' + model + '/'
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'

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
            
        fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
        fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_ts.nc'

        fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)*conv
        fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)*conv

        fnamec2 = basepath + casename2 + subdir + casename2 + '_' + field + '_' + timstr2 + '_ts.nc'
        fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstr2 + '_ts.nc'

        fldc2 = cnc.getNCvar(fnamec2,field.upper(),timesel=timesel)*conv
        fldp2 = cnc.getNCvar(fnamep2,field.upper(),timesel=timesel)*conv

        if sia==1:
            field='sia'
    else:
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
        
        frootc2 = basepath + casename2 + subdir + casename2 + '_' + field + '_'
        frootp2 =  basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_'

        fldc2 = np.append(cnc.getNCvar(frootc2+'001-061_ts.nc',ncfield,
                                      timesel='0002-01-01,061-12-31',levsel=level)*conv,
                         cnc.getNCvar(frootc2+'062-111_ts.nc',ncfield,levsel=level)*conv,
                         axis=0)
        fldp2 = np.append(cnc.getNCvar(frootp2+'001-061_ts.nc',ncfield,
                                      timesel='0002-01-01,061-12-31',levsel=level)*conv,
                         cnc.getNCvar(frootp2+'062-111_ts.nc',ncfield,levsel=level)*conv,
                         axis=0)
        

lat = cnc.getNCvar(fnamec,'lat')
lon = cnc.getNCvar(fnamec,'lon')


if sigtype=='cont' or sigoff==1:
    suff='pdf'
else:
    suff='png'

# annual time-series (3d)
seastsc = cutl.seasonalize_monthlyts(fldc,timeavg)
seastsp = cutl.seasonalize_monthlyts(fldp,timeavg)
#anntsc = cutl.annualize_monthlyts(fldc)
#anntsp = cutl.annualize_monthlyts(fldp)

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

if threed==1:
    fieldstr=field+str(level/100)
else:
    fieldstr=field
    
#tstatb,pvalb = sp.stats.ttest_ind(anntsp,anntsc,axis=0,equal_var=False) # basically the same as above
# Note that NaN is returned for zero variance (I think..from googling..)
# If that is the case, pcolormesh() needs a masked_array rather than ndarray (??)
#  : http://stackoverflow.com/questions/7778343/pcolormesh-with-missing-values

if plotann:
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
            cmin=cminn; cmax=cmaxn # @@@ guessing
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
            
            ttl = 'r' + str(ridx+1)

        seastdc = np.std(seasfldc,axis=0)
        seastdp = np.std(seasfldp,axis=0)

        tstat,pval = sp.stats.ttest_ind(seasfldp,seasfldc,axis=0)
        plotfld = np.mean(seasfldp,axis=0)-np.mean(seasfldc,axis=0)
        if normbystd:
            plotfld = plotfld/seastdc
            cmin=cminn; cmax=cmaxn # @@@ guessing
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




sigs = np.ones((12,fldc.shape[1],fldc.shape[2]))


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

    
    tstat = np.zeros((len(seasons),nlev,nlat))
    pval = np.zeros((len(seasons),nlev,nlat))
    fldcallseas = np.zeros((len(seasons),nlev,nlat)) # @@ dont' need to do this anymore, not saving all
    fldpallseas = np.zeros((len(seasons),nlev,nlat))
    
    incr = (cmaxm-cminm) / (cmlen)
    conts = np.arange(cminm,cmaxm+incr,incr)

    if addobs:
        sims = bcasename + 'r1', bcasename+'r2',bcasename+'r3',bcasename+'r4', bcasename+'r5',bcasename+'ens',bcasename,'kemhadctl'
    else:
        sims = bcasename + 'r1', bcasename+'r2',bcasename+'r3',bcasename+'r4', bcasename+'r5',bcasename+'ens',bcasename
#    cidx=0 # traverse cols (seasons)
#    ridx=0 # traverse rows
#    fig6,ax6 = plt.subplots(7,4) # 1 row for e/ of 5 ens members, plus mean, plus meanBC
#    fig6.set_size_inches(8,12)
    fig6,ax6 = plt.subplots(4,len(sims)) # 1 row for e/ of 5 ens members, plus mean, plus meanBC
    fig6.set_size_inches(12,8)  
    fig6.subplots_adjust(hspace=.15,wspace=.05)
        
    for ridx in range(0,len(sims)): # traverse rows

        cidx=0

        if ridx==5: # 2nd to last row is for the ens mean
            frootc = basepath + bcasename + 'ens' + subdir + bcasename +\
                     'ens' + '_' + field + '_'
            frootp = basepath + bcasenamep + 'ens' + subdir + bcasenamep +\
                     'ens' + '_' + field + '_'
            if threed==0:
                fnamec =  frootc + timstr + '_ts.nc'
                fnamep =  frootp + timstrp + '_ts.nc'
            else:
                fnamec = frootc + '001-061_ts.nc'
                fnamec2 = frootc + '062-121_ts.nc'
                fnamep = frootp + '001-061_ts.nc'
                fnamep2 = frootp + '062-121_ts.nc'
            
            rowl = 'r mean' # row label
            
        elif ridx==6: # last row is for meanBC
            frootc = basepath + bcasename + subdir + bcasename +\
                     '_' + field + '_'
            frootp = basepath + bcasenamep + subdir + bcasenamep +\
                     '_' + field + '_'
            if threed==0:
                fnamec = frootc + timstr2 + '_ts.nc'
                fnamep = frootp + timstr2 + '_ts.nc'
            else:
                fnamec = frootc + '001-061_ts.nc'
                fnamec2 = frootc + '062-111_ts.nc'
                fnamep = frootp + '001-061_ts.nc'
                fnamep2 = frootp + '062-111_ts.nc'
            
            rowl = 'meanBC'

        elif ridx==7: # observation simulation
            frootc = basepath + 'kemhadctl' + subdir + 'kemhadctl' + '_' + field + '_'
            frootp = basepath + 'kemhadpert' + subdir + 'kemhadpert' + '_' + field + '_'
            if threed==0:
                fnamec =  frootc + timstr + '_ts.nc'
                fnamep =  frootp + timstrp + '_ts.nc'
            else:
                fnamec = frootc + '001-061_ts.nc'
                fnamec2 = frootc + '062-121_ts.nc'
                fnamep = frootp + '001-061_ts.nc'
                fnamep2 = frootp + '062-121_ts.nc'

            rowl = 'had'
            
        else:
            frootc =  basepath + bcasename + 'r' + str(ridx+1) + subdir + bcasename +\
                     'r' + str(ridx+1) + '_' + field + '_'
            frootp = basepath + bcasenamep + 'r' + str(ridx+1) + subdir + bcasenamep +\
                     'r' + str(ridx+1) + '_' + field + '_'
            if threed==0:
                fnamec = frootc + timstr + '_ts.nc'
                fnamep = frootp + timstrp + '_ts.nc'
            else:
                fnamec = frootc + '001-061_ts.nc'
                fnamec2 = frootc + '062-121_ts.nc'
                fnamep = frootp + '001-061_ts.nc'
                fnamep2 = frootp + '062-121_ts.nc'
            
            rowl = 'r' + str(ridx+1)
        
        for sea in seasons: # traverse cols
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
            else:
                if threed==0:
                    fldcsea = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                                   seas=sea)*conv
                    fldpsea = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                                   seas=sea)*conv
                else:
                    fldcsea = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,061-12-31',levsel=level,
                                                   seas=sea)*conv,
                                        cnc.getNCvar(fnamec2,ncfield,levsel=level,seas=sea)*conv,axis=0)
                    fldpsea = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,061-12-31',levsel=level,
                                                   seas=sea)*conv,
                                        cnc.getNCvar(fnamep2,ncfield,levsel=level,seas=sea)*conv,axis=0)
                    
            tstat[cidx,:,:],pval[cidx,:,:] = sp.stats.ttest_ind(fldpsea,fldcsea,axis=0)
            fldcallseas[cidx,:,:] = np.mean(fldcsea,axis=0)
            fldpallseas[cidx,:,:] = np.mean(fldpsea,axis=0)

            if pct:
                plotfld = (fldpallseas[cidx,:,:]-fldcallseas[cidx,:,:]) / fldcallseas[cidx,:,:] *100
                cminm=cminmpct
                cmaxm=cmaxmpct
            else:
                plotfld = fldpallseas[cidx,:,:] - fldcallseas[cidx,:,:]

            bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',\
                            axis=ax,suppcb=1,latlim=latlim)#@@
            if cidx==0: # when row index is 0, set simulation
                ax.set_title(rowl)

            if sigoff==0:
                cplt.addtsigm(bm,pval[cidx,:,:],lat,lon,type=sigtype)
                
            if ridx==0: # when col index is 0, set season
                ax.set_ylabel(sea)

            cidx = cidx+1

    cbar_ax = fig6.add_axes([.91,.25, .02,.5])
    fig6.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    plt.suptitle(fieldstr)
    if printtofile:
        if sigoff==0:
            sigstr='sig' + sigtype
        else:
            sigstr=''
        if addobs:
            obsstr = 'had'
        else:
            obsstr = ''
        if latlim!= None:
            latstr=str(latlim)
        else:
            latstr=''
            
        if pct:
            fig6.savefig(fieldstr + 'pctdiff' + sigstr + '_enssubplot' + obsstr +
                         '_seas_nh' + latstr + '.' + suff)
        else:
            fig6.savefig(fieldstr + 'diff' + sigstr + '_enssubplot' + obsstr + '_seas_nh'
                         + latstr + '.' + suff)



if plotzonmean==1 or plotseacyc==1:
    # get data for either zonal mean or sea cycle figures
    if plotzonmean:
        pass # seasons is defined above
    elif plotseacyc:
        latlim = 70 # for area averaging
        seasons = con.get_mon()

    if addobs:
        sims = bcasename + 'r1', bcasename+'r2',bcasename+'r3',bcasename+'r4', bcasename+'r5',bcasename+'ens',bcasename,'kemhadctl'
    else:
        sims = bcasename + 'r1', bcasename+'r2',bcasename+'r3',bcasename+'r4', bcasename+'r5',bcasename+'ens',bcasename

    if sia==1:
        field = 'sicn' # while getting the data...
        
    tstatdict = dict.fromkeys(sims,{}); pvaldict = dict.fromkeys(sims,{})
    fldcdict = dict.fromkeys(sims,{}); fldpdict = dict.fromkeys(sims,{})
    fldcstddict = dict.fromkeys(sims,{}); fldpstddict = dict.fromkeys(sims,{})
    flddiffdict = dict.fromkeys(sims,{}); flddmaskdict = dict.fromkeys(sims,{})

    for ridx in range(0,len(sims)): # # of simulations
        seatstatdict=dict.fromkeys(seasons); seapvaldict=dict.fromkeys(seasons)
        seafldcdict=dict.fromkeys(seasons); seafldpdict=dict.fromkeys(seasons)
        seafldcstddict=dict.fromkeys(seasons); seafldpstddict=dict.fromkeys(seasons)
        seadiffdict=dict.fromkeys(seasons); seadmaskdict=dict.fromkeys(seasons)
        
        if ridx==5: # 2nd to last sim is the ens mean
            frootc = basepath + bcasename + 'ens' + subdir + bcasename +\
                     'ens' + '_' + field + '_'
            frootp = basepath + bcasenamep + 'ens' + subdir + bcasenamep +\
                     'ens' + '_' + field + '_'
            if threed==0:
                fnamec =  frootc + timstr + '_ts.nc'
                fnamep =  frootp + timstrp + '_ts.nc'
            else:
                fnamec = frootc + '001-061_ts.nc'
                fnamec2 = frootc + '062-121_ts.nc'
                fnamep = frootp + '001-061_ts.nc'
                fnamep2 = frootp + '062-121_ts.nc'

            sim = bcasename + 'ens'

        elif ridx==6: # last sim is meanBC
            frootc = basepath + bcasename + subdir + bcasename +\
                     '_' + field + '_'
            frootp = basepath + bcasenamep + subdir + bcasenamep +\
                     '_' + field + '_'
            if threed==0:
                fnamec = frootc + timstr2 + '_ts.nc'
                fnamep = frootp + timstr2 + '_ts.nc'
            else:
                fnamec = frootc + '001-061_ts.nc'
                fnamec2 = frootc + '062-111_ts.nc'
                fnamep = frootp + '001-061_ts.nc'
                fnamep2 = frootp + '062-111_ts.nc'

            sim = bcasename
        elif ridx==7: # observation simulation
            frootc = basepath + 'kemhadctl' + subdir + 'kemhadctl' + '_' + field + '_'
            frootp = basepath + 'kemhadpert' + subdir + 'kemhadpert' + '_' + field + '_'
            if threed==0:
                fnamec =  frootc + timstr + '_ts.nc'
                fnamep =  frootp + timstrp + '_ts.nc'
            else:
                fnamec = frootc + '001-061_ts.nc'
                fnamec2 = frootc + '062-121_ts.nc'
                fnamep = frootp + '001-061_ts.nc'
                fnamep2 = frootp + '062-121_ts.nc'

            sim = 'kemhadctl'

        else:
            frootc =  basepath + bcasename + 'r' + str(ridx+1) + subdir + bcasename +\
                     'r' + str(ridx+1) + '_' + field + '_'
            frootp = basepath + bcasenamep + 'r' + str(ridx+1) + subdir + bcasenamep +\
                     'r' + str(ridx+1) + '_' + field + '_'
            if threed==0:
                fnamec = frootc + timstr + '_ts.nc'
                fnamep = frootp + timstrp + '_ts.nc'
            else:
                fnamec = frootc + '001-061_ts.nc'
                fnamec2 = frootc + '062-121_ts.nc'
                fnamep = frootp + '001-061_ts.nc'
                fnamep2 = frootp + '062-121_ts.nc'

            sim = bcasename + 'r' + str(ridx+1)

        for sii,sea in enumerate(seasons):

            if plotzonmean:
                ncparams = {'seas': sea}
            elif plotseacyc:
                ncparams = {'monsel': sii+1}

            # Now get the data
            if field=='turb':
                print 'not implemented @@'
                field='hfl'; fieldb='hfs'
                fldczm = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                               **ncparams)*conv + cnc.getNCvar(fnamecb,fieldb.upper(),
                                               timesel=timesel,**ncparams)*conv
                fldpzm = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                               **ncparams)*conv + cnc.getNCvar(fnamepb,fieldb.upper(),
                                               timesel=timesel,**ncparams)*conv 
                field='turb'
            else:
                if threed==0:
                    fldczm = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                                   **ncparams)*conv
                    fldpzm = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                                   **ncparams)*conv

                else:                
                    ncparams['levsel'] = level
                    fldczm = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                        cnc.getNCvar(fnamec2,ncfield,**ncparams)*conv,
                                        axis=0)
                    fldpzm = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                        cnc.getNCvar(fnamep2,ncfield,**ncparams)*conv,
                                        axis=0)
                if sia==1:                    
                    fldczm = cutl.calc_seaicearea(fldczm,lat,lon)
                    fldpzm = cutl.calc_seaicearea(fldpzm,lat,lon)
                    
                    

            if plotzonmean:
                fldczm = np.mean(fldczm[...,:-1],axis=2)
                fldpzm = np.mean(fldpzm[...,:-1],axis=2)

            elif plotseacyc:
                if field == 'sia':
                    #calc total area instead of average
                    fldczm = np.sum(np.sum(fldczm[lat>0,:],axis=2),axis=1)
                    fldpzm = np.sum(np.sum(fldpzm[lat>0,:],axis=2),axis=1)
                else:
                    fldczm = cutl.polar_mean_areawgted3d(fldczm,lat,lon,latlim=latlim)
                    fldpzm = cutl.polar_mean_areawgted3d(fldpzm,lat,lon,latlim=latlim)

            seafldcstddict[sea] = np.std(fldczm,axis=0)
            seafldpstddict[sea] = np.std(fldpzm,axis=0)
            ttmp,pvtmp = sp.stats.ttest_ind(fldpzm,fldczm,axis=0)
            seatstatdict[sea] = ttmp
            seapvaldict[sea] = pvtmp
            seafldcdict[sea] =  np.mean(fldczm,axis=0)
            seafldpdict[sea] =  np.mean(fldpzm,axis=0)
            seadiffdict[sea] = np.mean(fldpzm,axis=0)- np.mean(fldczm,axis=0)
            seadmaskdict[sea] = ma.masked_where(pvtmp>siglevel,seadiffdict[sea])

            # end loop through seasons

        fldcstddict[sim] = seafldcstddict
        fldpstddict[sim] = seafldpstddict
        tstatdict[sim] = seatstatdict
        pvaldict[sim] = seapvaldict
        fldcdict[sim] = seafldcdict
        fldpdict[sim] = seafldpdict
        flddiffdict[sim] = seadiffdict
        flddmaskdict[sim] = seadmaskdict

        # end loop through simulations
    if sia==1:
        field = 'sia' # put back after getting the data


darkolivegreen1 = np.array([202, 255, 112])/255 # terrible
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
colordict = {'kemctl1r1': darkseagreen, 'kemctl1r2': darkseagreen4, 'kemctl1r3': lightsteelblue3,
             'kemctl1r4': lightsteelblue4, 'kemctl1r5': steelblue4, 'kemctl1ens': dodgerblue,
             'kemctl1': orangered4, 'kemhadctl': darkolivegreen3 }
if addobs:
    skip=-3 # use this to only take std dev over ens members r*
else:
    skip=-2


# now make the figures ==============================================
if plotzonmean:

    tmpval = fldcdict.values()[0]

    fldcestd = np.zeros((len(seasons),tmpval.values()[0].shape[0]))
    fldpestd = np.zeros((len(seasons),tmpval.values()[0].shape[0]))

    # take standard deviation over ensemble members:
    for sii,enssea in enumerate(seasons):

        tmpc = np.zeros((5,tmpval.values()[0].shape[0])) # 5 sims
        tmpp = np.zeros((5,tmpval.values()[0].shape[0]))
        for eii,enssim in enumerate(sims[0:skip]): # skip last 2 sims
            tmpc[eii,...] = fldcdict[enssim][enssea] # gather all ens members for calc
            tmpp[eii,...] = fldpdict[enssim][enssea]

        fldcestd[sii,...] = np.std(tmpc,axis=0) # for each season, what is sigma
        fldpestd[sii,...] = np.std(tmpp,axis=0)

    print fldcestd.shape



    # climo
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for ii,ax in enumerate(axs.flat):
        sea = seasons[ii]
        ax.plot(lat,fldcdict[bcasename+'r1'][sea],color=colors[0])
        ax.plot(lat,fldcdict[bcasename+'r2'][sea],color=colors[1])
        ax.plot(lat,fldcdict[bcasename+'r3'][sea],color=colors[2])
        ax.plot(lat,fldcdict[bcasename+'r4'][sea],color=colors[3])
        ax.plot(lat,fldcdict[bcasename+'r5'][sea],color=colors[4])
        ax.plot(lat,fldcdict[bcasename+'ens'][sea],color=colors[5],linewidth=2)
        ax.plot(lat,fldcdict[bcasename][sea],color=colors[6],linewidth=2)
        if addobs:
            ax.plot(lat,fldcdict['kemhadctl'][sea],color=colordict['kemhadctl'],linewidth=2)

        ax.plot(lat,fldpdict[bcasename+'r1'][sea],color=colors[0],linestyle='--')
        ax.plot(lat,fldpdict[bcasename+'r2'][sea],color=colors[1],linestyle='--')
        ax.plot(lat,fldpdict[bcasename+'r3'][sea],color=colors[2],linestyle='--')
        ax.plot(lat,fldpdict[bcasename+'r4'][sea],color=colors[3],linestyle='--')
        ax.plot(lat,fldpdict[bcasename+'r5'][sea],color=colors[4],linestyle='--')
        ax.plot(lat,fldpdict[bcasename+'ens'][sea],color=colors[5],linestyle='--',linewidth=2)
        ax.plot(lat,fldpdict[bcasename][sea],color=colors[6],linestyle='--',linewidth=2)
        if addobs:
            ax.plot(lat,fldpdict['kemhadctl'][sea],color=colordict['kemhadctl'],linestyle='--',linewidth=2)

        ax.set_xlim(-90,90)
        ax.set_title(sea)
        
    ax.set_xlabel('lat')
    if printtofile:
        fig.savefig(fieldstr + '_ens_meanBC_allseassp_zonmean.pdf')

    # standard deviation
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for ii,ax in enumerate(axs.flat):
        sea = seasons[ii]
        ax.plot(lat,fldcstddict[bcasename+'r1'][sea],color=colors[0])
        ax.plot(lat,fldcstddict[bcasename+'r2'][sea],color=colors[1])
        ax.plot(lat,fldcstddict[bcasename+'r3'][sea],color=colors[2])
        ax.plot(lat,fldcstddict[bcasename+'r4'][sea],color=colors[3])
        ax.plot(lat,fldcstddict[bcasename+'r5'][sea],color=colors[4])
        ax.plot(lat,fldcstddict[bcasename+'ens'][sea],color=colors[5],linewidth=2)
        ax.plot(lat,fldcstddict[bcasename][sea],color=colors[6],linewidth=2)
        if addobs:
            ax.plot(lat,fldcstddict['kemhadctl'][sea],color=colordict['kemhadctl'],linewidth=2)
        ax.plot(lat,fldcestd[ii,...],color='k',linewidth=2)

        ax.plot(lat,fldpstddict[bcasename+'r1'][sea],color=colors[0],linestyle='--')
        ax.plot(lat,fldpstddict[bcasename+'r2'][sea],color=colors[1],linestyle='--')
        ax.plot(lat,fldpstddict[bcasename+'r3'][sea],color=colors[2],linestyle='--')
        ax.plot(lat,fldpstddict[bcasename+'r4'][sea],color=colors[3],linestyle='--')
        ax.plot(lat,fldpstddict[bcasename+'r5'][sea],color=colors[4],linestyle='--')
        ax.plot(lat,fldpstddict[bcasename+'ens'][sea],color=colors[5],linestyle='--',linewidth=2)
        ax.plot(lat,fldpstddict[bcasename][sea],color=colors[6],linestyle='--',linewidth=2)
        if addobs:
            ax.plot(lat,fldpstddict['kemhadctl'][sea],color=colordict['kemhadctl'],linestyle='--',linewidth=2)
        ax.plot(lat,fldpestd[ii,...],color='k',linestyle='--',linewidth=2)
        
        ax.set_xlim(-90,90)
        ax.set_title(sea)
        
    ax.set_xlabel('lat')
    if printtofile:
        fig.savefig(fieldstr + 'STD_ens_meanBC_allseassp_zonmean.pdf')

    # differences
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for ii,ax in enumerate(axs.flat):
        sea = seasons[ii]

        ax.plot(lat,flddiffdict[bcasename+'r1'][sea],color=colors[0])
        ax.plot(lat,flddiffdict[bcasename+'r2'][sea],color=colors[1])
        ax.plot(lat,flddiffdict[bcasename+'r3'][sea],color=colors[2])
        ax.plot(lat,flddiffdict[bcasename+'r4'][sea],color=colors[3])
        ax.plot(lat,flddiffdict[bcasename+'r5'][sea],color=colors[4])
        ax.plot(lat,flddiffdict[bcasename+'ens'][sea],color=colors[5],linewidth=2)
        ax.plot(lat,flddiffdict[bcasename][sea],color=colors[6],linewidth=2)
        if addobs:
            ax.plot(lat,flddiffdict['kemhadctl'][sea],color=colordict['kemhadctl'],linewidth=2)

        ax.plot(lat,flddmaskdict[bcasename+'r1'][sea],linestyle='none',marker='s',color=colors[0])
        ax.plot(lat,flddmaskdict[bcasename+'r2'][sea],linestyle='none',marker='s',color=colors[1])
        ax.plot(lat,flddmaskdict[bcasename+'r3'][sea],linestyle='none',marker='s',color=colors[2])
        ax.plot(lat,flddmaskdict[bcasename+'r4'][sea],linestyle='none',marker='s',color=colors[3])
        ax.plot(lat,flddmaskdict[bcasename+'r5'][sea],linestyle='none',marker='s',color=colors[4])
        ax.plot(lat,flddmaskdict[bcasename+'ens'][sea],linestyle='none',marker='s',color=colors[5])
        ax.plot(lat,flddmaskdict[bcasename][sea],linestyle='none',marker='s',color=colors[6])
        if addobs:
            ax.plot(lat,flddmaskdict['kemhadctl'][sea],color=colordict['kemhadctl'],linestyle='none',marker='s')

        ax.set_xlim(0,90)
        ax.set_title(field + ': ' + sea)
        
    ax.set_xlabel('lat')
    ax.legend(sims,'upper left', prop=fontP,ncol=2)
    if printtofile:
        fig.savefig(fieldstr + 'diff_ens_meanBC_allseassp_zonmean_nh.pdf')

    # standard deviation differences
    fig,axs = plt.subplots(4,1)
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    for ii,ax in enumerate(axs.flat):
        sea = seasons[ii]

        ax.plot(lat,fldpstddict[bcasename+'r1'][sea]-fldcstddict[bcasename+'r1'][sea],color=colors[0])
        ax.plot(lat,fldpstddict[bcasename+'r2'][sea]-fldcstddict[bcasename+'r2'][sea],color=colors[1])
        ax.plot(lat,fldpstddict[bcasename+'r3'][sea]-fldcstddict[bcasename+'r3'][sea],color=colors[2])
        ax.plot(lat,fldpstddict[bcasename+'r4'][sea]-fldcstddict[bcasename+'r4'][sea],color=colors[3])
        ax.plot(lat,fldpstddict[bcasename+'r5'][sea]-fldcstddict[bcasename+'r5'][sea],color=colors[4])
        ax.plot(lat,fldpstddict[bcasename+'ens'][sea]-fldcstddict[bcasename+'ens'][sea],color=colors[5],linewidth=2)
        ax.plot(lat,fldpstddict[bcasename][sea]-fldcstddict[bcasename][sea],color=colors[6],linewidth=2)
        if addobs:
            ax.plot(lat,fldpstddict['kemhadctl'][sea]-fldcstddict['kemhadctl'][sea],color=colordict['kemhadctl'],linewidth=2)

        ax.plot(lat,fldpestd[ii,...]-fldcestd[ii,...],color='k',linewidth=2)

        ax.set_xlim(0,90)
        ax.set_title(field + ': ' + sea)
        
    ax.set_xlabel('lat')
    
    if printtofile:
        fig.savefig(fieldstr + 'STDdiff_ens_meanBC_allseassp_zonmean_nh.pdf')

if plotseacyc:

    fldcestd = np.zeros((len(seasons))) # ensemble std dev
    fldpestd = np.zeros((len(seasons)))

    fldcestddict = dict.fromkeys(seasons)
    fldpestddict = dict.fromkeys(seasons)
    
    # take standard deviation over ensemble members:
    for sii,enssea in enumerate(seasons):

        tmpc = np.zeros(5) # 5 sims
        tmpp = np.zeros(5)
        for eii,enssim in enumerate(sims[0:skip]): # skip last 2 sims
            tmpc[eii] = fldcdict[enssim][enssea] # gather all ens members for calc
            tmpp[eii] = fldpdict[enssim][enssea]

        fldcestddict[enssea] = np.std(tmpc,axis=0) # for each season, what is sigma
        fldpestddict[enssea] = np.std(tmpp,axis=0)

    #print fldcestd.shape


    fontP = fm.FontProperties()
    fontP.set_size('small')
    
    seacycddict = dict.fromkeys(sims)
    seacyccdict = dict.fromkeys(sims)
    seacycpdict = dict.fromkeys(sims)
    seacycmaskdict = dict.fromkeys(sims)
    seacycdstddict = dict.fromkeys(sims)
    seacyccstddict = dict.fromkeys(sims)
    seacycpstddict = dict.fromkeys(sims)
    
    # first, reorganize the data: put all months into an array
    for skey,simval in flddiffdict.iteritems():
        cfld = fldcdict[skey]
        pfld = fldpdict[skey]
        
        seacycd = np.zeros(12)
        seacycc = np.zeros(12)
        seacycp = np.zeros(12)
        #print skey
        #sii=0
        #for seakey,seaval in simval.iteritems():
        for sii,seakey in enumerate(seasons): # have to do it this way to get in right order
            
            seacycd[sii] = simval[seakey]
            seacycc[sii] = cfld[seakey]
            seacycp[sii] = pfld[seakey]
            #sii=sii+1
            
        seacycddict[skey] = seacycd
        seacyccdict[skey] = seacycc
        seacycpdict[skey] = seacycp
        
    for skey,simval in flddmaskdict.iteritems():        
        seacyc = np.zeros(12)
        #print skey
        #sii=0
        #for seakey,seaval in simval.iteritems():
        for sii,seakey in enumerate(seasons): # have to do it this way to get in right order
            seacyc[sii] = simval[seakey]
            #sii=sii+1
            
        seacycmaskdict[skey] = seacyc

    # standard deviation
    for skey,simval in fldcstddict.iteritems():
        pstd = fldpstddict[skey]
        
        seacycd = np.zeros(12)
        seacycc = np.zeros(12)
        seacycp = np.zeros(12)
        #print skey
        #sii=0
        #for seakey,seaval in simval.iteritems():
        for sii,seakey in enumerate(seasons): # have to do it this way to get in right order
            pstdval = pstd[seakey]
            
            seacycd[sii] = pstdval - simval[seakey]
            seacycc[sii] = simval[seakey]
            seacycp[sii] = pstdval
            #sii=sii+1
            
        seacycdstddict[skey] = seacycd
        seacyccstddict[skey] = seacycc
        seacycpstddict[skey] = seacycp

    for sii,seakey in enumerate(seasons):
        fldcestd[sii] = fldcestddict[seakey]
        fldpestd[sii] = fldpestddict[seakey]

    # Now can plot the seasonal cycle
    mons=np.arange(1,13)

    # climo
    fig,axs = plt.subplots()
    for skey,val in seacyccdict.iteritems():

        if skey == 'kemctl1' or skey == 'kemctl1ens' or skey=='kemhadctl':
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=3)
        else:
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=2)

    for skey,val in seacycpdict.iteritems():
        if skey == 'kemctl1' or skey == 'kemctl1ens' or skey=='kemhadctl':
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=3,linestyle='--')
        else:
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=2,linestyle='--')
         
    plt.legend(seacyccdict.keys(),leglocs[0], prop=fontP,ncol=2)
    plt.xlim((1,12))
    plt.xlabel('Month')
    plt.ylabel(fieldstr)
    plt.title('Climos')

    if printtofile:
        fig.savefig(fieldstr + '_ens_meanBC_seacyc_pol' + str(latlim) + 'N.pdf')


    # differences
    fig,axs = plt.subplots()
    for skey,val in seacycddict.iteritems():
        if skey == 'kemctl1' or skey == 'kemctl1ens' or skey=='kemhadctl':
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=3)
        else:
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=2)
    # significance dots
    for skey,val in seacycmaskdict.iteritems():
        axs.plot(mons,val,linestyle='none',color=colordict[skey],marker='s')
         
    plt.legend(seacycddict.keys(),leglocs[1], prop=fontP,ncol=2)
    plt.xlim((1,12))
    plt.xlabel('Month')
    plt.ylabel(fieldstr)
    plt.title('Anomalies')

    if printtofile:
        fig.savefig(fieldstr + 'diff_ens_meanBC_seacyc_pol' + str(latlim) + 'N.pdf')


    # Standard deviation
    fig,axs = plt.subplots()
    for skey,val in seacyccstddict.iteritems():

        if skey == 'kemctl1' or skey == 'kemctl1ens' or skey=='kemhadctl':
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=3)
        else:
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=2)

    for skey,val in seacycpstddict.iteritems():
        if skey == 'kemctl1' or skey == 'kemctl1ens' or skey=='kemhadctl':
            axs.plot(range(1,13),val,color=colordict[skey],linestyle='--',linewidth=3)
        else:
            axs.plot(range(1,13),val,color=colordict[skey],linestyle='--',linewidth=2)

    axs.plot(range(1,13),fldcestd,color='k',linewidth=3)
    axs.plot(range(1,13),fldpestd,color='k',linewidth=3,linestyle='--')
             
    plt.legend(seacyccstddict.keys(),leglocs[2], prop=fontP,ncol=2)
    plt.xlim((1,12))
    plt.xlabel('Month')
    plt.ylabel(fieldstr)
    plt.title('Sigma')

    if printtofile:
        fig.savefig(fieldstr + 'STD_ens_meanBC_seacyc_pol' + str(latlim) + 'N.pdf')


    # Difference in standard deviation
    fig,axs = plt.subplots()
    for skey,val in seacycdstddict.iteritems():
        if skey == 'kemctl1' or skey == 'kemctl1ens' or skey=='kemhadctl':
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=3)
        else:
            axs.plot(range(1,13),val,color=colordict[skey],linewidth=2)
    axs.plot(range(1,13),fldpestd-fldcestd,color='k',linewidth=3)
             
    plt.legend(seacycdstddict.keys(),leglocs[3], prop=fontP,ncol=2)
    plt.xlim((1,12))
    plt.xlabel('Month')
    plt.ylabel(fieldstr)
    plt.title('Sigma anomalies')

    if printtofile:
        fig.savefig(fieldstr + 'STDdiff_ens_meanBC_seacyc_pol' + str(latlim) + 'N.pdf')



    
if testhadisst:
    # here want to check which ensemble run is most similar to hadisst in terms of SICN
    # can use fldc2 and fldp2 for hadisst data
    # need to get all the r* runs

    corrlim=40 # the latitude north of which to consider the pattern correlation
    lmask = con.get_t63landmask()
    
    fldcdict = {}; fldpdict = {}
    diffdict = {}
    pcorrdict = {}
    
    for ridx in range(0,7): # # of simulations

        if ridx==5: # 2nd to last sim is the ens mean
            frootc = basepath + bcasename + 'ens' + subdir + bcasename +\
                     'ens' + '_' + field + '_'
            frootp = basepath + bcasenamep + 'ens' + subdir + bcasenamep +\
                     'ens' + '_' + field + '_'

            fnamec =  frootc + timstr + '_ts.nc'
            fnamep =  frootp + timstrp + '_ts.nc'
            sim = bcasename + 'ens'

        elif ridx==6: # last sim is meanBC
            frootc = basepath + bcasename + subdir + bcasename +\
                     '_' + field + '_'
            frootp = basepath + bcasenamep + subdir + bcasenamep +\
                     '_' + field + '_'
            fnamec = frootc + '001-111_ts.nc'
            fnamep = frootp + '001-111_ts.nc'
            sim = bcasename

        else:
            frootc =  basepath + bcasename + 'r' + str(ridx+1) + subdir + bcasename +\
                     'r' + str(ridx+1) + '_' + field + '_'
            frootp = basepath + bcasenamep + 'r' + str(ridx+1) + subdir + bcasenamep +\
                     'r' + str(ridx+1) + '_' + field + '_'

            fnamec = frootc + timstr + '_ts.nc'
            fnamep = frootp + timstrp + '_ts.nc'
            sim = bcasename + 'r' + str(ridx+1)

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
            
        pcorrdict[sim] = pcorr

    mons = range(1,13)
    fontP = fm.FontProperties()
    fontP.set_size('small')
    
    plt.figure();
    plt.plot(mons,diffdict['kemctl1r1'],color=colordict['kemctl1r1'],linewidth=2)
    plt.plot(mons,diffdict['kemctl1r2'],color=colordict['kemctl1r2'],linewidth=2)
    plt.plot(mons,diffdict['kemctl1r3'],color=colordict['kemctl1r3'],linewidth=2)
    plt.plot(mons,diffdict['kemctl1r4'],color=colordict['kemctl1r4'],linewidth=2)
    plt.plot(mons,diffdict['kemctl1r5'],color=colordict['kemctl1r5'],linewidth=2)
    plt.plot(mons,diffdict['kemctl1ens'],color=colordict['kemctl1ens'],linewidth=3)
    plt.plot(mons,diffdict['kemctl1'],color=colordict['kemctl1'],linewidth=3)
    plt.legend(('r1','r2','r3','r4','r5','ens','meanBC'),prop=fontP)
    plt.xlim((1,12))
    plt.title('RMSE compared to HadISST')
    if printtofile:
        plt.savefig('sicnRMSE_ens_v_hadisst_seacycle.pdf')

    plt.figure();
    plt.plot(mons,pcorrdict['kemctl1r1'],color=colordict['kemctl1r1'],linewidth=2)
    plt.plot(mons,pcorrdict['kemctl1r2'],color=colordict['kemctl1r2'],linewidth=2)
    plt.plot(mons,pcorrdict['kemctl1r3'],color=colordict['kemctl1r3'],linewidth=2)
    plt.plot(mons,pcorrdict['kemctl1r4'],color=colordict['kemctl1r4'],linewidth=2)
    plt.plot(mons,pcorrdict['kemctl1r5'],color=colordict['kemctl1r5'],linewidth=2)
    plt.plot(mons,pcorrdict['kemctl1ens'],color=colordict['kemctl1ens'],linewidth=3)
    plt.plot(mons,pcorrdict['kemctl1'],color=colordict['kemctl1'],linewidth=3)
    plt.legend(('r1','r2','r3','r4','r5','ens','meanBC'),'upper left',prop=fontP)
    plt.xlim((1,12))
    plt.title('pattern correlation with HadISST')
    if printtofile:
        plt.savefig('sicnPatternCorr_ens_v_hadisst_seacycle.pdf')

    wgts = con.get_monweights()
    print 'r1: ' + str(np.average(diffdict['kemctl1r1'],weights=wgts))
    print 'r2: ' + str(np.average(diffdict['kemctl1r2'],weights=wgts))
    print 'r3: ' + str(np.average(diffdict['kemctl1r3'],weights=wgts))
    print 'r4: ' + str(np.average(diffdict['kemctl1r4'],weights=wgts))
    print 'r5: ' + str(np.average(diffdict['kemctl1r5'],weights=wgts))
    print 'ens: ' + str(np.average(diffdict['kemctl1ens'],weights=wgts))
    print 'meanBC: ' + str(np.average(diffdict['kemctl1'],weights=wgts))

    print '===== PATTERN CORR ====='
    print 'r1: ' + str(np.average(pcorrdict['kemctl1r1'],weights=wgts))
    print 'r2: ' + str(np.average(pcorrdict['kemctl1r2'],weights=wgts))
    print 'r3: ' + str(np.average(pcorrdict['kemctl1r3'],weights=wgts))
    print 'r4: ' + str(np.average(pcorrdict['kemctl1r4'],weights=wgts))
    print 'r5: ' + str(np.average(pcorrdict['kemctl1r5'],weights=wgts))
    print 'ens: ' + str(np.average(pcorrdict['kemctl1ens'],weights=wgts))
    print 'meanBC: ' + str(np.average(pcorrdict['kemctl1'],weights=wgts))

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
