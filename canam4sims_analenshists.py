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
import cccmacmaps as ccm
import matplotlib.font_manager as fm
import copy

# while I'm still creating these modules, have to reload to get changes
cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)
ccm = reload(ccm)
cnc = reload(cnc)

plt.close("all")
plt.ion()

printtofile=1


seasonal=0 # seasonal maps (DJF, MAM, JJA, SON)
addobs=1 # add mean of kemhad* runs to line plots, seasonal maps (so far@@)
latlim = 60 # lat limit for NH plots. Set to None otherwise.

maskland=False # To get ocean-lake-only hist: mutually exclusive
maskocean=False # To get land-only hist: mutually exclusive

seasonal=1 # mutually exclusive
monthly=0  # mutually exclusive

#seasons = 'DJF','MAM','JJA','SON'
seasons = 'SON','DJF','MAM','JJA'

model = 'CanAM4'
threed=0
sia=0 # is the requested field sea ice area
mtype=''

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


# # # ######## set Field info ###################
# gz, t, u, v, q (3D !)
# st, sic, sicn (sia), gt, pmsl, pcp, hfl, hfs, turb, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)
field = 'pmsl'
print field

# only for threed vars
#level = 30000
#level = 50000 # 500hPa
level = 70000

## # HERE THE SECOND SET OF RUNS IS HadISST BC's
## if testhadisst:
##     casename2 = 'kemhadctl'
##     casenamep2 = 'kemhadpert'
##     timstr2=timstr
##     field='sicn'


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
    xlims = -14,14; ylims = 0,0.6; ylims2=0,11 
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

    xlims = -25,25; ylims = 0,0.2; ylims2=0,6 
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
    xlims = -4,4; ylims=0,3; ylims2 = 0,20
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
bp=con.get_basepath()
basepath=bp['basepath'] + model + '/'; subdir=bp['subdir']

## plat = platform.system()

## if plat == 'Darwin':  # means I'm on my mac
##     basepath = '/Volumes/MyPassport1TB/DATA/CanSISE/' + model + '/'
##     subdir = '/'
## else:  # on linux workstation in Vic
##     basepath = '/home/rkm/work/DATA/' + model + '/'
##     subdir = '/ts/'

# prob don't need this next chunk if looping down below @@@
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
                         cnc.getNCvar(frootc2+'062-121_ts.nc',ncfield,levsel=level)*conv,
                         axis=0)
        fldp2 = np.append(cnc.getNCvar(frootp2+'001-061_ts.nc',ncfield,
                                      timesel='0002-01-01,061-12-31',levsel=level)*conv,
                         cnc.getNCvar(frootp2+'062-121_ts.nc',ncfield,levsel=level)*conv,
                         axis=0)
        

lat = cnc.getNCvar(fnamec,'lat')
lon = cnc.getNCvar(fnamec,'lon')


if threed==1:
    fieldstr=field+str(level/100)
else:
    fieldstr=field
    





months=con.get_mon()

## darkolivegreen1 = np.array([202, 255, 112])/255 # terrible
## darkolivegreen3 = np.array([162, 205, 90])/255.
## darkseagreen = np.array([143, 188, 143])/255.
## darkseagreen4 = np.array([105, 139, 105])/255.
## dodgerblue = np.array([30, 144, 255])/255. 
## orangered4 = np.array([139, 37, 0])/255.
## lightsteelblue3 = np.array([162, 181, 205])/255. # more grey looking
## lightsteelblue4 = np.array([110, 123, 139])/255. # more grey looking
## steelblue3 = np.array([79, 148, 205])/255.  # more blue looking
## steelblue4 = np.array([54, 100, 139])/255.  # more blue looking


## colordict = {'kemctl1r1': darkseagreen, 'kemctl1r2': darkseagreen4, 'kemctl1r3': lightsteelblue3,
##              'kemctl1r4': lightsteelblue4, 'kemctl1r5': steelblue4, 'kemctl1ens': dodgerblue,
##              'kemctl1': orangered4, 'kemhadctl': darkolivegreen3 }

colordict = {'r1': ccm.get_linecolor('warm1'), #'firebrick'),
             'r4': ccm.get_linecolor('warm2'), #'firebrick1'),
             'r3': ccm.get_linecolor('warm3'), #'yelloworange'),#'chocolate1'),
             'r5': ccm.get_linecolor('warm4'), #'darkyellow') #'skyblue'), #yelloworange'),
             'r2': ccm.get_linecolor('warm5'), #'steelblue3'), #'darkgoldenrod1'),
             'ens': ccm.get_linecolor('magenta'),
             '': ccm.get_linecolor('mediumblue'), #'mediumpurple1'), #darkyellow'), # @@ empty key?
             'kemhad': ccm.get_linecolor('deepskyblue')}


# Load data into dicts

# order ens simulations in order of most ice loss in melt season to least. Then ens mean, PERT2, observations if requested
## sims = bcasename+'r1', bcasename+'r4', bcasename+'r3', bcasename+'r5', bcasename+'r2', bcasename+'ens',bcasename
## if addobs:
##     sims = sims + ('kemhadctl',)

sims = 'r1', 'r4','r3','r5','r2','ens','' # suffixes to bcasename and bcasenamep
if addobs:
    sims = sims + ('kemhad',)



if sia==1:
    field = 'sicn' # while getting the data...

flddiffhistdict = dict.fromkeys(sims,{})
binedgedict = dict.fromkeys(sims,{})

#for ridx in range(0,len(sims)): # # of simulations
for ridx,sim in enumerate(sims):
    print sim
    seadiffhistdict=dict.fromkeys(seasons)
    seabinedgedict=dict.fromkeys(seasons)

    if sim == 'kemhad':
        frootc = basepath + sim + 'ctl' + subdir + sim + 'ctl' + '_' + field + '_'
        frootp = basepath + sim + 'pert' + subdir + sim + 'pert' + '_' + field + '_'
    else:
        frootc =  basepath + bcasename + sim + subdir + bcasename + sim + '_' + field + '_'
        frootp = basepath + bcasenamep + sim + subdir + bcasenamep + sim + '_' + field + '_'
    
    if threed==0:
        fnamec =  frootc + timstr + '_ts.nc'
        fnamep =  frootp + timstrp + '_ts.nc'
    else:
        fnamec = frootc + '001-061_ts.nc'
        fnamec2 = frootc + '062-121_ts.nc'
        fnamep = frootp + '001-061_ts.nc'
        fnamep2 = frootp + '062-121_ts.nc'

    print fnamec
    print fnamep

    ## I think can get rid of all this now, with all runs having
        # same timeseries length and organizing simulation names this way
    ## #if ridx==5: # 2nd to last sim is the ens mean
    ## if sim == bcasename+'ens':
    ##     frootc = basepath + bcasename + 'ens' + subdir + bcasename +\
    ##              'ens' + '_' + field + '_'
    ##     frootp = basepath + bcasenamep + 'ens' + subdir + bcasenamep +\
    ##              'ens' + '_' + field + '_'
    ##     if threed==0:
    ##         fnamec =  frootc + timstr + '_ts.nc'
    ##         fnamep =  frootp + timstrp + '_ts.nc'
    ##     else:
    ##         fnamec = frootc + '001-061_ts.nc'
    ##         fnamec2 = frootc + '062-121_ts.nc'
    ##         fnamep = frootp + '001-061_ts.nc'
    ##         fnamep2 = frootp + '062-121_ts.nc'

    ##     #sim = bcasename + 'ens'

    ## #elif ridx==6: # last sim is meanBC
    ## elif sim = bcasename:
    ##     frootc = basepath + bcasename + subdir + bcasename +\
    ##              '_' + field + '_'
    ##     frootp = basepath + bcasenamep + subdir + bcasenamep +\
    ##              '_' + field + '_'
    ##     if threed==0:
    ##         fnamec = frootc + timstr2 + '_ts.nc'
    ##         fnamep = frootp + timstr2 + '_ts.nc'
    ##     else:
    ##         fnamec = frootc + '001-061_ts.nc'
    ##         fnamec2 = frootc + '062-121_ts.nc'
    ##         fnamep = frootp + '001-061_ts.nc'
    ##         fnamep2 = frootp + '062-121_ts.nc'

    ##     #sim = bcasename
    ## #elif ridx==7: # observation simulation
    ## elif sim = 'kemhadctl':
    ##     frootc = basepath + 'kemhadctl' + subdir + 'kemhadctl' + '_' + field + '_'
    ##     frootp = basepath + 'kemhadpert' + subdir + 'kemhadpert' + '_' + field + '_'
    ##     if threed==0:
    ##         fnamec =  frootc + timstr + '_ts.nc'
    ##         fnamep =  frootp + timstrp + '_ts.nc'
    ##     else:
    ##         fnamec = frootc + '001-061_ts.nc'
    ##         fnamec2 = frootc + '062-121_ts.nc'
    ##         fnamep = frootp + '001-061_ts.nc'
    ##         fnamep2 = frootp + '062-121_ts.nc'

    ##     #sim = 'kemhadctl'

    ## else:
    ##     frootc =  basepath + sim + subdir + sim + '_' + field + '_'
    ##     # have to figure out per casename@@
    ##     frootp = basepath + bcasenamep + 'r' + str(ridx+1) + subdir + bcasenamep +\
    ##              'r' + str(ridx+1) + '_' + field + '_'
        
    ##     if threed==0:
    ##         fnamec = frootc + timstr + '_ts.nc'
    ##         fnamep = frootp + timstrp + '_ts.nc'
    ##     else:
    ##         fnamec = frootc + '001-061_ts.nc'
    ##         fnamec2 = frootc + '062-121_ts.nc'
    ##         fnamep = frootp + '001-061_ts.nc'
    ##         fnamep2 = frootp + '062-121_ts.nc'

    ##     #sim = bcasename + 'r' + str(ridx+1)

    for sii,sea in enumerate(seasons):
        print sea
        
        if seasonal:# @@ fix
            ncparams = {'seas': sea}
        elif monthly: # @@fix
            ncparams = {'monsel': sii+1}

        # Now get the data
        if field=='turb':
            print 'not implemented @@'
            field='hfl'; fieldb='hfs'
            fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                           **ncparams)*conv + cnc.getNCvar(fnamecb,fieldb.upper(),
                                           timesel=timesel,**ncparams)*conv
            fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                           **ncparams)*conv + cnc.getNCvar(fnamepb,fieldb.upper(),
                                           timesel=timesel,**ncparams)*conv 
            field='turb'
        else:
            if threed==0:
                fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,
                                               **ncparams)*conv
                fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,
                                               **ncparams)*conv

            else:                
                ncparams['levsel'] = level
                fldc = np.append(cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                    cnc.getNCvar(fnamec2,ncfield,**ncparams)*conv,
                                    axis=0)
                fldp = np.append(cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,061-12-31',**ncparams)*conv,
                                    cnc.getNCvar(fnamep2,ncfield,**ncparams)*conv,
                                    axis=0)
            if sia==1:                    
                fldc = cutl.calc_seaicearea(fldc,lat,lon)
                fldp = cutl.calc_seaicearea(fldp,lat,lon)


        # Prepare the data here: @@
        #  @@ from hists script:
        areas = cutl.calc_cellareas(lat,lon,repeat=fldc.shape)
        areas = areas[:,lat>latlim,:]
        fldcorig=copy.copy(fldc)
        fldporig=copy.copy(fldp)
        fldc = fldc[:,lat>latlim,:]
        fldp = fldp[:,lat>latlim,:]
        if maskland:
            lmask = con.get_t63landmask(repeat=fldc.shape)
            areas = ma.masked_where(lmask[:,lat>latlim,:]==-1,areas)
            fldc = ma.masked_where(lmask[:,lat>latlim,:]==-1,fldc)
            fldp = ma.masked_where(lmask[:,lat>latlim,:]==-1,fldp)
            mtype='ocn'
        elif maskocean:
            lmask = con.get_t63landmask(repeat=fldc.shape)
            areas = ma.masked_where(lmask[:,lat>latlim,:]!=-1,areas)
            fldc = ma.masked_where(lmask[:,lat>latlim,:]!=-1,fldc)
            fldp = ma.masked_where(lmask[:,lat>latlim,:]!=-1,fldp)
            mtype='lnd'

        totarea = np.sum(np.sum(areas,axis=2),axis=1)
        totarea = np.tile(totarea,(fldc.shape[1],fldc.shape[2],1))
        totarea = np.transpose(totarea,(2,0,1))
        weights = areas/totarea
        #weights = np.tile(weights,(len(lat),len(lon),fldc.shape[0]))
        #weights = np.transpose(weights,(2,0,1))
                          

        histfld = fldp-fldc
        print 'max val: ' + str(np.max(np.max(histfld)))
        print 'min val: ' + str(np.min(np.min(histfld)))
        if maskocean or maskland: # pass only the non-masked data to histogram
            (seadiffhistdict[sea],seabinedgedict[sea]) = np.histogram(histfld.compressed(),
                                                                      bins=100,weights=weights.compressed(),
                                                                      density=True)
        else:
            (seadiffhistdict[sea],seabinedgedict[sea]) = np.histogram(histfld.flatten(),
                                                                      bins=100,weights=weights.flatten(),
                                                                      density=True)
        #widths = np.diff(binedges)
        
        ## seafldcstddict[sea] = np.std(fldc,axis=0)
        ## seafldpstddict[sea] = np.std(fldpaxis=0)
        ## ttmp,pvtmp = sp.stats.ttest_ind(fldp,fldc,axis=0)
        ## seatstatdict[sea] = ttmp
        ## seapvaldict[sea] = pvtmp
        ## seafldcdict[sea] =  np.mean(fldc,axis=0)
        ## seafldpdict[sea] =  np.mean(fldp,axis=0)
        ## seadiffdict[sea] = np.mean(fldp,axis=0)- np.mean(fldc,axis=0)
        ## seadmaskdict[sea] = ma.masked_where(pvtmp>siglevel,seadiffdict[sea])

        # end loop through seasons

    flddiffhistdict[sim] = seadiffhistdict
    binedgedict[sim] = seabinedgedict
    
    ## fldcstddict[sim] = seafldcstddict
    ## fldpstddict[sim] = seafldpstddict
    ## tstatdict[sim] = seatstatdict
    ## pvaldict[sim] = seapvaldict
    ## fldcdict[sim] = seafldcdict
    ## fldpdict[sim] = seafldpdict
    ## flddiffdict[sim] = seadiffdict
    ## flddmaskdict[sim] = seadmaskdict

    # end loop through simulations
if sia==1:
    field = 'sia' # put back after getting the data

# @@ plot HISTOGRAMS HERE

widths=0.15 # np.diff(binedges)
sea='DJF'

fig = plt.figure()
for sim in sims:
    binedges = binedgedict[sim][sea]
    plt.bar(binedges[0:-1],flddiffhistdict[sim][sea],
            width=widths,color=colordict[sim],alpha=0.5)
#plt.xlim(xlims)
plt.grid()


fig = plt.figure()
for sim in sims:
    binedges = binedgedict[sim][sea]
    plotfld = flddiffhistdict[sim][sea]
    plotfld = plotfld / np.sum(plotfld)*100
    
    plt.bar(binedges[0:-1],plotfld,
            width=widths,color=colordict[sim],alpha=0.5)
#plt.xlim(xlims)
plt.grid()
plt.ylabel('% of total')
plt.title(sea + ' ' + fieldstr + '>' + str(latlim) + 'N')
if printtofile:
    fig.savefig(fieldstr + 'diffpcthist_' + str(latlim) + 'N' + mtype + '_' + sea + '.pdf')


#xlims = -15,15; ylims = 0,0.6; ylims2=0,11 # ST DJF @@
fig2,axs = plt.subplots(2,4)
for sii,ax in enumerate(axs.flat):
    sim = sims[sii]
    plotfld=flddiffhistdict[sim][sea]
    binedges = binedgedict[sim][sea]
    ax.bar(binedges[0:-1],plotfld,
           width=widths,color=colordict[sim],alpha=0.7,
           linewidth=0,edgecolor='none')
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_title(sim)
    ax.grid()
if printtofile:
    fig2.savefig(fieldstr + 'diffhist_' + str(latlim) + 'N' + mtype + '_' + sea + '_subplt.pdf')


fig2,axs = plt.subplots(2,4)
for sii,ax in enumerate(axs.flat):
    sim = sims[sii]
    plotfld=flddiffhistdict[sim][sea]
    plotfld=plotfld / np.sum(plotfld)*100
    binedges = binedgedict[sim][sea]
    ax.bar(binedges[0:-1],plotfld,
           width=widths,color=colordict[sim],alpha=0.7,
           linewidth=0,edgecolor='none')
    ax.set_xlim(xlims)
    ax.set_ylim(ylims2)
    ax.set_title(sim)
    if sii==0 or sii==4:
        ax.set_ylabel('% of total')
    ax.grid()
if printtofile:
    fig2.savefig(fieldstr + 'diffpcthist_' + str(latlim) + 'N' + mtype + '_' + sea + '_subplt.pdf')


# # Cumulative PDF
# @@ consider using 'ax.fill_between' for the r set of runs
# @@ also see example of cum here:
# http://matplotlib.org/examples/pylab_examples/histogram_demo_extended.html
fig = plt.figure()
for sim in sims:
    binedges = binedgedict[sim][sea]
    plotfld = flddiffhistdict[sim][sea]

    dcum = np.zeros(plotfld.shape)
    for el in np.arange(1,len(dcum)+1):
        dcum[el-1] = np.sum(plotfld[0:el])
        
    plotfld = dcum / np.sum(plotfld)
    
    plt.plot(binedges[0:-1],plotfld,
            color=colordict[sim],linewidth=2)
plt.xlim(xlims)
plt.grid()
plt.ylabel('Fraction of total')
plt.title(sea + ' ' + fieldstr + '>' + str(latlim) + 'N')
if printtofile:
    fig.savefig(fieldstr + 'difffracCDD_' + str(latlim) + 'N' + mtype + '_' + sea + '.pdf')

    





"""
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
                fnamec2 = frootc + '062-121_ts.nc'
                fnamep = frootp + '001-061_ts.nc'
                fnamep2 = frootp + '062-121_ts.nc'
            
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
"""
