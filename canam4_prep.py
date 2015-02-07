""" canam4_prep.py:
      9/9/2014: Goal is to prep the data for 1. retrieval/calcs and 2. plotting
                To be used in conjunction with simulation_funcs.py
                Splitting canam4sims_analens.py
                This is particular to canam4, but the functions being called
                should not be (in the end...) so can have prep scripts for CanESM2, etc
                
"""
import constants as con
import simulation_funcs as sfnc
import statshelper as sh

sh=reload(sh)
sfnc = reload(sfnc)
con = reload(con)
#ccm = reload(ccm)

plt.close("all")
plt.ion()

# Will use this script to set up all the stuff.
# for example, plot_seasonal_maps() requires:
       # simpair keys
       # fielddict with 'field','ncfield','fieldstr'
       # pparams dict with cmin, cmax, cmap, latlim, type, suppcb ....
       #    if vert=True, also screen, levlim, addcontlines
       # coords = {'lat': con.get_t63lat(), 'lon': con.get_t63lon()}


printtofile=True

field = 'st'
smclim=True
level=50000 # for threed

addcont=False # overlay map with contours
sigoff=False # if True, don't add significance
effdof=False # use effective deg of freedom or no. Set to False.
field2='gz'
#field2='pmsl'
level2=50000

# seasonalmap, seasonalvert, plotzonmean, plotseacyc, pattcorrwithtime, plotregmean,calcregmeanwithtime, calcregunccascade,timetosig, timetosigsuper, plotscatter
plottype='calcregunccascade' 
projtype='eastere' #'eastere' # 'nh','sh','sq','eastere','nastere'

# None, nh, polcap60, polcap65, polcap70, eurasia, eurasiamori, eurasiasth,eurasiathin,eurasiathinw,eurasiathine,ntham, nthatl, bks, bksmori, soo
region='polcap60' #'eurasia' #'eurasiamori'
screen=True # just for seasonalvert
seacyclatlim=60
withlat=False
pattcorryr=False # need anymore?
latlim = None # None #45 # lat limit for NH plots. Set to None otherwise. use 45 for BC-type maps
round=False # if latlim is not None, this specifies whether the polar stereographic figure should be square or round
levlim= 100 # level limit for vertical ZM plots (in hPa). ignored if screen=True
fallwin=False # just SON and DJF
bimon=True # do bi-montly seasons instead
figtrans=False # for maps/vert: make seasons the cols and sims the rows if True. auto True for simsforpaper


# Choose how to handle the data ==============
normbystd=False
pct = False # if True, do calculation as a percent (@@not fully implemented)
halftime=False # get only the first 60yrs. make sure to set the other flag the opp
halftime2=False # get only the last 60yrs. make sure to set the other flag the opp

# Choose what simulations to add =============
#  default is R1-5, ENS
canens=False # just the CAN ensemble (E1-E5) plus mean, plus mean of R ensemble. option to addobs only.
allens=False # this is ONLY the ensemble means, plus superensemble
sensruns=False # sensruns only: addr4ct=1,addsens=1. others=0 no meanBC, r mean, or obs
ivar=False # this will show ENS (TOT) and ENSE (ANTH) and their difference = internal var
simsforpaper=False # ANT, TOT, HAD, NSIDC only. best for maps and zonal mean figs (not line plots)
simsforpaperwace=False # Just R1,R5 (warm/cold Eurasia) or whatever it's set to below.
antcat=False # this is the concatenation of ens members within each ensemble (really only useful for ANT)
bothcat=False # can do concatenation of both ensembles if want to. These are useful for timetosig
onlyens=False # just do ensemble means ANT and TOT

addobs=True # add mean of kemhad* & kemnsidc* runs to line plots, seasonal maps. 
addr4ct=False # add kem1pert2r4ct (constant thickness version of ens4)
addsens=False # add sensitivity runs (kem1pert1b, kem1pert3)
addrcp=False # add kem1rcp85a simulation (and others if we do more)
addcanens=True # add "initial condition" ensemble of kemctl1/kem1pert2
addsuper=False # add superensemble mean


# not sure these flags are in use?
sigtype = 'cont' # significance: 'cont' or 'hatch' which is default
siglevel=0.05



ptparams={}
ptparams['smclim'] = smclim # for maps typically
ptparams['latlim'] = latlim
ptparams['round'] = round
ptparams['levlim'] = levlim
ptparams['region'] = region # plotregmean
ptparams['screen'] = screen # seasonalvert
ptparams['seacyclatlim'] = seacyclatlim # plotseacyc
ptparams['withlat'] = withlat # plotseacyc
ptparams['pattcorryr'] = pattcorryr # pattcorrwithtime

if plottype=='seasonalvert':
    seasonalvert=True
else:
    seasonalvert=False

## # Choose type of plot =========================
## seasonalmap=False # seasonal maps (SON, DJF, MAM, JJA)
## seasonalvert=False # seasonal vertical zonal means instead of maps
## screen=False # whether to have screen-style vertical zonal means

## plotzonmean=False # plotzonmean,plotseacyc,pattcorrwithtime are mutually exclusive

## plotseacyc=False # plotzonmean,plotseacyc,pattcorrwithtime are mutually exclusive
## seacyclatlim=60 # southern limit for plotting polar mean seasonal cycles (line plot)
## withlat=False # plot the seasonal cycle with latitude dimension too (only for plotseacyc=1)@@for now just std over ens
## #squatseacyc=False # plot seacycle figs as shorter than wide
## #squatterseacyc=True # even shorter, for paper

## pattcorrwithtime=False # plot pattern correlation with time for each ens member
## pattcorryr=False # if True, do a yearly anomaly pattern rather than time-integrated

## plotregmean=True
## region = 'ntham' # None, polcap60, polcap65, polcap70, eurasia, eurasiamori, ntham, nthatl, bks, bksmori, soo

testhadisst=0 # check which ens member most similar to hadisst

model='CanAM4'


# #################################################################
#   Probably don't need to modify below if everything goes well.
# #################################################################

# initialize things that don't get set otherwise
savestr='' # string for plot filenames
## threed=False # is the field three dimensional
## sia=False # is the requested field sea ice area
## conv=1
timesel=None
## isflux=False
shadeens=('histBC',)
corrlim=45 # southern lat limit for pattern correlation with time

# set up simulations and figure filename strings
sims = 'R1','R4','R3','R5','R2','ENS'#,'ENSE'#,'CAN' # R's in order of sea ice loss
defaultsims=sims
seasons = ('SON','DJF','MAM','JJA')
biseas = ('SO','ND','JF') # @@@ so far only these implemented. expecting to add all 11/25/14
#biseas = ('ND',) # @@@@


if simsforpaper: # best for maps only
    sims = ('HAD','NSIDC','ENSE','ENS')
    savestr = '_forpap4' # add ENS. # 4 means fig is transposed
    
    if bimon:
        seasons=biseas
        savestr = savestr + 'bimon'
    else:
        seasons=('SON','DJF')
    figtrans=True
    
elif simsforpaperwace:
    
    print '@@@ simsforpaperwace is WACE paper -- ND cold and warm extremes, R1, R5 only'
    savestr = '_forpapwace'; sims = ('R1','R5'); seasons=('ND',); figtrans=True # 95% stat sig
    #print '@@@ simsforpaperwace is WACE paper -- ND ENSEMBLE MEANS ONLY'
    #savestr = '_forpapwace_ensmean'; sims = ('ENS','ENSE'); seasons=('ND',); figtrans=True # 95% stat sig
    
    #print '@@@ simsforpaper is WACE B paper right now -- ND cold and warm extremes, E4, E1 only'
    #savestr = '_forpapwaceb'; sims = ('E4','E1'); seasons=('ND',); figtrans=False # 90% stat sig
    
elif sensruns: # add sensitivity runs. with Shaded ENS. don't plot meanBC, mean of ens
    sims = sims[0:5] + ('R4ct','CANnosst','CANnothk') # @@ change to E1nosst, etc?
    savestr = '_sensruns'
    if bimon:
        seasons=biseas
        savestr = savestr + 'bimon'
    if fallwin:
        seasons=('SON','DJF')
        savestr = savestr + '_SONDJF'
    
elif canens: # do canens instead of r ens. Useful for maps.
    sims = ('E1','E2','E3','E4','E5','ENSE','ENS')
    savestr = '_canensonly'
    shadeens=('histIC',)
    if addobs:
        sims = sims + ('HAD','NSIDC')
        savestr = savestr + 'obs'
    if addsuper:
        sims = sims + ('ESPR',)
        savestr = savestr + 'spr'
    if bimon:
        seasons=biseas
        savestr = savestr + 'bimon'   
    if fallwin:
        seasons=('SON','DJF')
        savestr = savestr + '_SONDJF' 
elif allens: # just do the ens means and superensemble mean
    sims = ('ENS','ENSE','ESPR')
    savestr = '_allens'
    if bimon:
        seasons=biseas
        savestr = savestr + 'bimon'  
    if fallwin:
        seasons=('SON','DJF')
        savestr = savestr + '_SONDJF'
    # addobs? @@@
elif ivar:
    sims = ('ENS','ENSE','IVAR')
    savestr = '_ivar'
    smclim=True # @@
    if bimon:
        seasons=biseas
        savestr = savestr + 'bimon'
    if fallwin:
        seasons=('SON','DJF')
        savestr = savestr + '_SONDJF'
elif antcat: # eventually add totcat @@@
    sims = ('ENSECAT',)
    savestr = '_antcat'
    figtrans=True
    if bimon:
        seasons=biseas
        savestr = savestr + 'bimon'
    if fallwin:
        seasons=('SON','DJF')
        savestr = savestr + '_SONDJF'
elif bothcat: # both antcat and totcat
    sims = ('ENSECAT','ENSCAT')
    savestr = '_anttotcat'
    figtrans=True
    if bimon:
        seasons=biseas
        savestr = savestr + 'bimon'
    if fallwin:
        seasons=('SON','DJF')
        savestr = savestr + '_SONDJF'

elif onlyens:
    sims = ('ENSE','ENS')
    savestr = '_onlyens'
    if bimon:
        seasons=biseas
        savestr = savestr + 'bimon'   
    if fallwin:
        seasons=('SON','DJF')
        savestr = savestr + '_SONDJF' 
else:
    savestr=savestr+'_ens' # default sims
    
    if addcanens:
        sims = sims + ('E1','E2','E3','E4','E5','ENSE') # E1=CAN
        savestr = savestr + 'canens'
        shadeens=shadeens+('histIC',)
    else:
        sims = sims + ('ENSE',)
    if addsuper:
        sims = sims + ('ESPR',)
        savestr = savestr + 'spr'
    if addobs:
        sims = sims + ('HAD','NSIDC')
        savestr = savestr + 'obs'
    if addr4ct:
        sims = sims + ('R4ct',)
        savestr = 'r4ct' # for figure filenames
    if addsens:
        sims = sims + ('CANnosst','CANnothk') # control is kemctl1 (or '' key)
        savestr = savestr + 'sens'
    if addrcp:
        sims = sims + ('RCPa',) # control is kemctl1
        savestr = savestr + 'rcpa'
    if bimon:
        seasons=biseas
        savestr = savestr + 'bimon'  
    if fallwin:
        seasons=('SON','DJF')
        savestr = savestr + '_SONDJF'
    
import load_fldmeta as ld
ld=reload(ld)
## ############################ start copy to load_fldmeta.py
## # # # ######## set Field info ###################
## # gz, t, u, v, q (3D !)
## # st, sic, sicn (sia), gt, pmsl, pcp, hfl, hfs, turb, net, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)


## print field

if halftime:
    timesel= '0002-01-01,0061-12-31'
elif halftime2:
    timesel='0062-01-01,0121-12-31'



## # # # ###########################################
## #   Shouldn't have to mod below....


## """  plev = 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, 12500, 
##     15000, 17500, 20000, 22500, 25000, 30000, 35000, 40000, 45000, 50000, 
##     55000, 60000, 65000, 70000, 75000, 77500, 80000, 82500, 85000, 87500, 
##     90000, 92500, 95000, 97500, 100000 ;
## """

## fdict = {'field': field, 'ncfield': None, 'fieldstr': None,
##          'units': None, 'conv': conv,
##          'nonstandardlev': nonstandardlev,
##          'threed': threed} # fielddict

## # reserved for expansion into the plotfunction call
## pparams = {'cmin': None, 'cmax': None, 'cmap': 'blue2red_20',              
##            'type':'nh', 'latlim': latlim} # plotparams
seacycylim=None
infodict ={'cmapclimo': 'Spectral_r','leglocs': None,
           'seacycylim': None, 'savestr': savestr,
           'model': model, 'sigtype': sigtype, 'sigoff': sigoff,
           'pct': pct, 'seacyclatlim': seacyclatlim, 'region': region,
           'shadeens': shadeens, 'corrlim': corrlim, 'figtrans':figtrans,
           'type': projtype} # random other info. projtype for maps only

fdict,pparams=ld.loadfldmeta(field,infodict,plottype,ptparams,level=level)

if addcont: # overlay map with another field in contours
    # start with just anomaly contours. @@later add option for climo
    fdict2,pparams2=ld.loadfldmeta(field2,infodict,plottype,ptparams,level=level2)

    
## if field == 'st':
##     fdict['units'] = 'K'
##     fdict['ncfield'] = field.upper()
##     fdict['fieldstr'] = field

##     pparams['cmin'] = -3; pparams['cmax'] = 3 # seasonal/monthly
##     pparams['cmap'] = 'blue2red_w20'
##     if smclim:
##         pparams['cmin'] = -1.5; pparams['cmax'] = 1.5 # seasonal/monthly
##         pparams['cmap'] = 'blue2red_20'
##         savestr= savestr + '_smclim'
        
##     #cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
##     #cminn = -5; cmaxn = 5 # for norm by std

##     if plotseacyc  and withlat:
##         pparams['cmin']=-.5;
##         pparams['cmax']=.5 # @@ will have to update this when add subplots
        
##     leglocs = 'upper left', 'upper left', 'upper right', 'upper left'
##     if region not in ('eurasia','ntham','nthatl'):
##         seacycylim=(-.5,4) # >70N
    
## elif field == 'sic':

##     fdict['units'] = 'm'
##     fdict['ncfield'] = field.upper()
##     fdict['fieldstr'] = field

##     pparams['cmin'] = -.5; pparams['cmax'] = .5 # seasonal/monthly

##     fdict['conv']=1/913.

##     pparams['cmap'] = 'red2blue_w20'
##     #pparams['leglocs'] = 'lower left', 'lower left', 'upper left', 'upper left'
    
## elif field == 'sicn' or field == 'sia':
##     fdict['units'] = 'frac'
##     fdict['ncfield'] = field.upper()
##     fdict['fieldstr'] = field

##     pparams['cmin'] = -.15; pparams['cmax'] = .15 # seasonal/monthly

##     pparams['cmap'] = 'red2blue_w20'
##     leglocs = 'lower left', 'lower left', 'upper left', 'upper left'

##     if field=='sia':
##         fdict['ncfield'] = 'SICN'
##         seacycylim=(-2e12,0) # for sia
##         sia=True
    
## elif field == 'gt':
##     units='K'
##     conv=1
##     cmin=-2; cmax=2
##     cminm=-3; cmaxm=3
##     cminp=-.5; cmaxp=.5 # for when pert is 'ctl'
##     cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
##     cmap = 'blue2red_w20'
## elif field == 'pmsl':
##     fdict['units'] = 'hPa'
##     fdict['ncfield'] = field.upper()
##     fdict['fieldstr'] = field

##     pparams['cmin'] = -2; pparams['cmax'] = 2 # seasonal/monthly
##     if smclim:
##         pparams['cmin'] = -1; pparams['cmax'] = 1 # seasonal/monthly
##         savestr= savestr + '_smclim'

##     pparams['cmap'] = 'blue2red_20'

##     if plotseacyc  and withlat:
##         pparams['cmin']=-1;
##         pparams['cmax']=1 # @@ will have to update this when add subplots

##     leglocs = 'lower left', 'lower left', 'upper center', 'lower left'
##     if region not in ('eurasia','ntham','nthatl'):
##         seacycylim=(-2,1.5) # >70N

    
##     ## units = 'hPa' # pretty sure hpa @@double check
##     ## conv = 1
##     ## cmin = -1; cmax = 1  # for anomaly plots
##     ## cminm=-2; cmaxm=2  # for monthly maps
##     ## cminp=cmin; cmaxp=cmax # for when pert is 'ctl'
##     ## cminmp=cminm; cmaxmp=cmaxm
##     ## cmap = 'blue2red_20'
##     ## ## print 'new cmap and small clim! @@ '
##     ## ## cmap = 'blue2red_w20'
##     ## ## cminm=-1; cmaxm=1  for monthly maps, small clim
    
##     ## cminn = -1; cmaxn = 1 # for norm by std

    
## elif field == 'pcp':
##     fdict['fieldstr'] = field
##     fdict['ncfield'] = field.upper()
    
##     fdict['units'] = 'mm/day' # original: kg m-2 s-1
##     fdict['conv'] = 86400 # convert from kg m-2 s-1 to mm/day

##     pparams['cmap'] = 'brown2blue_16w'

##     if pct:
##         pparams['cmin'] = -20; pparams['cmax'] = 20
##         fdict['units'] = '%'
##     else:
##         pparams['cmin'] = -0.2; pparams['cmax'] = 0.2
   
##     conv = 86400  # convert from kg m-2 s-1 to mm/day
##     cmin = -.2; cmax = .2  # for anomaly plots
##     cminp=-.15; cmaxp=.15
##     cminm = -.2; cmaxm = .2

##     ## print '@@ big clim!'
##     ## cmin = -.75; cmax = .75 
##     ## cminm = -.75; cmaxm = .75

##     ## print '@@ medium clim!'
##     ## cmin = -.5; cmax = .5 
##     ## cminm = -.5; cmaxm = .5

##     ## #cmap = 'PuOr'
##     ## cmap = 'brown2blue_16w'
##     ## cminpct=-12; cmaxpct=12
##     ## cminmpct=-20; cmaxmpct=20
##     ## cminmp =-.25; cmaxmp=.25
##     ## cminpctp=-8; cmaxpctp=8
##     ## cminpctmp=-12; cmaxpctmp=12
##     leglocs = 'upper left', 'upper left', 'upper left', 'upper left'
    
## elif field == 'hfl': # sfc upward LH flux
##     units = 'W/m2'
##     conv = 1
##     cmin = -5
##     cmax = 5
##     cminm = -8
##     cmaxm = 8

##     isflux=True
    
## elif field == 'hfs': # sfc upward SH flux
##     units = 'W/m2'
##     conv = 1
##     cmin = -5
##     cmax = 5
##     cminm = -8
##     cmaxm = 8

##     isflux=True
## elif field == 'turb': # combine hfl and hfs
##     units = 'W/m2'
##     conv=1
##     cmin=-10
##     cmax=10
##     cminm = -20
##     cmaxm = 20
##     cmap='blue2red_20'
    
## elif field == 'net': # net of all sfc fluxes
##     #print " 'net ' not yet implemented! @@"
##     units = 'W/m2'
##     conv=1
##     cmin=-10
##     cmax=10
##     cminm = -20
##     cmaxm = 20
##     cmap='blue2red_20'
##     leglocs = 'upper left', 'upper left', 'upper left', 'upper left'
##     seacycylim=(-5,25) # always >40N (where there is ice in CTL)

##     isflux=True
## elif field == 'flg': # net downward LW at the sfc.
##     units = 'W/m2'
##     conv = -1 # so positive heats atmos
##     cmin = -5
##     cmax = 5
##     cminm = -8
##     cmaxm = 8
##     leglocs = 'upper left', 'lower left', 'upper left', 'upper left'

##     isflux=True
## elif field == 'fsg': # net (absorbed) solar downard at sfc
##     units = 'W/m2'
##     conv = 1
##     cmin = -5
##     cmax = 5
##     cminm = -8
##     cmaxm = 8

##     isflux=True

## elif field == 'fn': # snow fraction
##     units = '%'
##     conv=100
##     cmin = -5
##     cmax = 5
##     cminm = -5
##     cmaxm = 5
##     cmap = 'red2blue_w20'

## elif field == 'pcpn': # snowfall rate (water equivalent, kg/m2/s)
##     #pct = 1; units='%'
##     units = 'mm/day'
##     conv = 86400  # convert from kg m-2 s-1 to mm/day (I think it's same as pcp) @@
##     cmap = 'brown2blue_16w'
##     cmin = -.1 # for anomaly plots
##     cmax = .1  # for anomaly plots
##     cminm = -.15
##     cmaxm = .15
##     cminpct=-12
##     cmaxpct=12
##     cminmpct=-25
##     cmaxmpct=25
##     leglocs = 'upper left', 'upper left', 'upper center', 'upper left'

## elif field == 'zn': # snow depth (m)
## #    pct=1; units='%'
##     units = 'cm'
##     conv = 100; # convert to cm
##     cmap = 'brown2blue_16w'
##     cmin = -2
##     cmax = 2
##     cminm = -2.5; cmaxm = 2.5
##     cminpct=-10
##     cmaxpct=10
##     cminmpct=-10
##     cmaxmpct=10
##     leglocs = 'upper right', 'lower left', 'upper right', 'lower right'

## elif field == 'su':
##     fdict['units'] = 'm/s'
##     fdict['ncfield'] = 'SU'    
##     fdict['conv'] = 1;
    
##     pparams['cmap'] = 'blue2red_20'
##     pparams['cmin'] = -1; pparams['cmax'] = 1
##     cminm = -1; cmaxm = 1
##     cminp = -.5; cmaxp=.5
##     cminmp = -.5; cmaxmp=.5
##     leglocs = 'upper left', 'lower left', 'upper left', 'upper left'
## elif field == 'sv':
##     units = 'm/s'
##     conv = 1;
##     cmap = 'blue2red_20'
##     cmin = -.5
##     cmax = .5
##     cminm = -.5
##     cmaxm = .5

## elif field == 't':

##     fdict['units'] = 'K'
##     fdict['ncfield'] = 'TEMP'

##     threed = True
##     fdict['conv'] = 1

##     pparams['cmap'] = 'blue2red_w20'
        
##     if level == 30000:
##         pparams['cmin'] = -.5; pparams['cmax'] = .5
##         #cminm = -.5; cmaxm = .5  # for monthly
##         #cminsea = -.5; cmaxsea = .5
##     elif level == 70000:
##         #cmin = -.3; cmax = .3
##         pparams['cmin'] = -.5; pparams['cmax'] = .5  # for monthly
##         #cminsea = -.5; cmaxsea = .5

##     if seasonalvert:
##         if screen:
##             pparams['levlim']=300
##             pparams['cmin']=-2.5; pparams['cmax']=2.5
##             #cminm=-2.5; cmaxm=2.5
##         else:
##             cmin = -.5; cmax = .5
##             pparams['cmin'] = -.8; pparams['cmax'] = .8 
## elif field == 'u':
##     threed = True
##     fdict['units'] = 'm/s'
##     fdict['ncfield'] = 'U'
##     fdict['conv'] = 1

##     if level == 30000:
##         #cmin = -2; cmax = 2
##         pparams['cmin'] = -3; pparams['cmax'] = 3
##         #cminsea = -3; cmaxsea = 3
##     else:
##         #cmin = -1; cmax = 1
##         pparams['cmin'] = -1; pparams['cmax'] = 1
##         #cminsea = -1; cmaxsea = 1

##     if seasonalvert:
##         #cmin=-.5; cmax=.5
##         pparams['cmin'] = -1; pparams['cmax'] = 1
##         if screen:
##             pparams['levlim']=300
        
##     pparams['cmap']='blue2red_20'

## elif field == 'gz':
    
##     fdict['units'] = 'm'
##     fdict['ncfield'] = 'PHI'

##     threed = True
##     fdict['conv'] = 1/con.get_g()
        
##     if level==30000:
##         pparams['cmin'] = -20; pparams['cmax'] = 20 # seasonal/monthly
##     else:
##         pparams['cmin'] = -15; pparams['cmax'] = 15 # seasonal/monthly
##         if region not in ('eurasia','ntham','nthatl'):
##             seacycylim=(-16,20) # >70N, 500hPa

##     if seasonalvert:
##         if screen:
##             #cmin=-10; cmax=10
##             #cminm=-25; cmaxm=25
##             pparams['cmin'] = -25; pparams['cmax'] = 25 # seasonal/monthly
##             pparams['levlim'] = 300
##         else:
##             #cmin = -10; cmax = 10
##             #cminm = -15; cmaxm = 15
##             pparams['cmin'] = -15; pparams['cmax'] = 15 # seasonal/monthly

##     pparams['cmap'] = 'blue2red_w20'
                
##     leglocs = 'upper left', 'upper left', 'upper right', 'upper left'

## else:
##     print 'No settings for ' + field

## if threed: # either map of a level, or a zonal mean
##     if seasonalvert:
##         fieldstr=field+'ZM'
##         field=field+'ZM'
##     else:
##         fieldstr=field+str(level/100) # figure names. want integer for string
##         field=field+str(level) # read in files
##     fdict['field']=field
##     fdict['fieldstr']=fieldstr
## else:
##     fdict['fieldstr']=field

## fluxes = 'hfl','hfs','flg' # LH, SH, LWdown

coords = {'lev': con.get_t63lev(), 'lat': con.get_t63lat(), 'lon': con.get_t63lon()}
## infodict['savestr'] = savestr
## infodict['leglocs'] = leglocs
## infodict['seacycylim'] = seacycylim
## fdict['isflux'] = isflux
## fdict['threed'] = threed
## ###################### end copy to load_fldmeta.py

# do an if elif elif ....
if plottype=='seasonalvert':
    infodict['screen']=screen
else:
    infodict['screen'] = None
        
if plottype in ('seasonalmap','seasonalvert'):

    print fdict
    print sims
    print pparams

    #pparams2['cmin']=-.6; pparams2['cmax']=.6 #@@@
        
    # this one does data processing and plotting together
    # some stuff in the function need to be removed or set differently.@@
    # marked in the function @@
    if addcont:
        addflds=(fdict2,)
        addpparams=(pparams2,)
    else:
        addflds=None
        addpparams=None

    thefig = sfnc.calc_plot_seasonal_maps(fdict,coords,sims,pparams,vert=seasonalvert,
                            loctimesel=timesel,info=infodict,seas=seasons,
                            printtofile=printtofile,addflds=addflds,addpparams=addpparams)

    if simsforpaperwace:
        # here I will adjust the figure for the WACE paper
        if figtrans:
            if projtype=='eastere':
                thefig.set_size_inches((5,8))
                thefig.subplots_adjust(hspace=.02,wspace=.02)
            else:
                thefig.set_size_inches((4,10))
            
        else:
            thefig.set_size_inches((10,5))
        theaxs = thefig.get_axes()
        ax1=theaxs[0]
        if 'ENS' in sims: # assume it is just the two ensemble means
            prstr='ANTTOT'
            ax1.set_title('a. Individual SIC forcings')
            ax1.set_ylabel('')
        else: # else assume it's the cooling and warming cases for paper
            prstr=''
            ax1.set_title('a. Significant cooling case')
            ax1.set_ylabel('Latitude',fontsize=12)
            if not figtrans:
                ax1.set_xlabel('Longitude',fontsize=12)

        ax2=theaxs[1]
        if 'ENS' in sims:  # assume it is just the two ensemble means
            ax2.set_title('b. Average SIC forcing')
            ax2.set_ylabel('')
        else:  # else assume it's the cooling and warming cases for paper
            ax2.set_title('b. Significant warming case')
            ax2.set_xlabel('Longitude',fontsize=12)
            if figtrans:
                ax2.set_ylabel('Latitude',fontsize=12)
            else:
                ax2.set_ylabel('',fontsize=12)

        if not figtrans: # one row
            axcb=theaxs[2] # colorbar
            axcb.set_position([.91,.25, .02,.5])
        else: # one col
            axcb=theaxs[2] # colorbar
            axcb.set_position([.15,0.07, 0.75, .02])
            if 'ENS' not in sims:
                axcb.set_xticks(np.arange(-1,1.2,0.2))
                axcb.set_xticklabels(('-1.0','','','','','0','','','','','1.0'))

        thefig.suptitle('')
        if addcont:
            thefig.savefig('wacefigure3_' + field + '_' + field2 + 'cont_trans_' + prstr + '.pdf')
            thefig.savefig('wacefigure3_' + field + '_' + field2 + 'cont_trans_' + prstr + '.eps')
        else:
            thefig.savefig('wacefigure3_' + field + '_trans_' + prstr + '.pdf')
            thefig.savefig('wacefigure3_' + field + '_trans_' + prstr + '.eps')
        
        
    

if plottype=='plotseacyc':

    dblob = sfnc.calc_seasonal_cycle(fdict,coords,sims,withlat=withlat,loctimesel=timesel,info=infodict)

    sfnc.plot_seasonal_cycle(dblob,fdict,sims,ptypes=('anom','stddev'),info=infodict,printtofile=printtofile)
    

if plottype=='plotzonmean':
    
    dblob = sfnc.calc_seasons(fdict,coords,sims,loctimesel=timesel,info=infodict,calctype='zonmean',effdof=effdof)
    sfnc.plot_zonmean_byseas(dblob,fdict,coords,sims,ptypes=('climo','anom','stddev','stdan'),info=infodict,printtofile=printtofile)

if plottype=='pattcorrwithtime':

    if pattcorryr==True:
        calctype='pattcorrwithtimeyr'
    else:
        calctype='pattcorrwithtime'
        
    dblob = sfnc.calc_seasons(fdict,coords,sims,loctimesel=timesel,info=infodict,calctype=calctype,effdof=effdof)
    sfnc.plot_pattcorrwithtime_byseas(dblob,fdict,sims,info=infodict,calctype=calctype,printtofile=printtofile)

if plottype=='plotregmean':

    dblob = sfnc.calc_seasons(fdict,coords,sims,seas=seasons,loctimesel=timesel,info=infodict,calctype='regmean',effdof=effdof)
    sfnc.plot_regmean_byseas(dblob,fdict,sims,info=infodict,seas=seasons,printtofile=printtofile)


    # DEBUG:
    ci=dblob['ci']
    pv=dblob['pval']


    import pandas as pd
    
    

    # calc whether ens means are sig different @@@@
    if shadeens != None: # right now, will never be none. Should change that. @@@

        # calc anom range vs mean anom for the given field within given ensemble
        allensdt,allensmdt = con.build_ensembles(shadeens,dblob,calctype='diff')
        for ens in shadeens:
            ensdf=pd.DataFrame(allensdt[ens])
            #ensmdf=pd.DataFrame(allensmdt[ens])
            ensrng = ensdf.max(axis=1)-ensdf.min(axis=1)
            ensrat = ensrng / ensdf.mean(axis=1)
            print ens + ' RANGE: '
            print str(ensrng)
            print ens + ' MEAN: '
            print str(ensdf.mean(axis=1))
            print ens + ' RATIO*100: '
            print str(ensrat*100)
            
        if len(shadeens)>1:
            sh.calc_ensemblestats(dblob,shadeens,seas=seasons)
    
if plottype=='timetosig' or plottype=='timetosigsuper':

    calctype='timetosig'
    pparams['cmin'] = 0; pparams['cmax'] = 120 # same clims for all runs/vars in this project
    pparams['cmap'] = 'Spectral_r' #'YlGnBu'
    
    if plottype=='timetosigsuper':
        calctype='timetosigsuper'
        pparams['cmax'] = 600
        
    dblob = sfnc.calc_seasons(fdict,coords,sims,loctimesel=timesel,info=infodict,
                              calctype=calctype,seas=seasons,effdof=effdof)
  
    sfnc.plot_seasonal_maps(dblob,fdict,coords,sims,pparams,plottype='timetosig',
                            vert=False,seas=seasons,info=infodict,printtofile=printtofile)

    
if plottype=='calcregmeanwithtime' or plottype=='calcregunccascade' or plottype=='calccomposites':
    
    dblob = sfnc.calc_seasons(fdict,coords,sims,seas=seasons,
                              loctimesel=timesel,info=infodict,calctype='regmeanwithtime',
                              effdof=effdof)

    infodict['region'] = 'eurasia'
    dblob2 = sfnc.calc_seasons(fdict,coords,sims,seas=seasons,
                               loctimesel=timesel,info=infodict,calctype='regmeanwithtime',
                               effdof=effdof)
    
    import pandas as pd

    # dblob should have regional means *with* time dimension
    # want to test sig different mean and variance b/w ctl and pert
    allstats=sh.calc_runstats(dblob,sims, seas=seasons,siglevel=siglevel)

    fpvaldf=pd.DataFrame(allstats['fpval']) # pval of f statistic for variance significance
    fpvaldft=fpvaldf.transpose()
    tpvaldf=pd.DataFrame(allstats['tpval']) # pval of t statistic for mean significance
    tpvaldft=tpvaldf.transpose()

    import cccmacmaps as ccm

    if plottype=='calcregunccascade':
        effdof=False

        col=('0.3',ccm.get_linecolor('firebrick'))
        if field=='st' and region in ('eurasia','eurasiamori'):
            #xlab = '$\Delta$ Eurasia SAT ($^\circ$C)'
            xlab = '($^\circ$C)'
        elif field=='st' and region=='polcap60':
            #xlab = '$\Delta$ >60$^\circ$N SAT ($^\circ$C)'
            xlab = '($^\circ$C)'
        else:
            xlab=None

        for sea in seasons:
            seas=(sea,)
            fig,axs = plt.subplots(2,1)
            fig.set_size_inches(4,8)
            fig.set_frameon(False)
            
            ax=axs[0] # first polar cap SAT
            if 'ND' in seas:
                annlab='a. Nov-Dec Polar SAT changes'
            else:
                annlab='a. Polar SAT changes'
            annloc=(0.02,0.94)
            sfnc.plot_uncertainty_cascade(dblob,fdict,coords,sims,pparams,
                                          info=infodict,seas=seas,ax=ax,xlab=xlab,annlab=annlab,annloc=annloc,
                                          color=col,effdof=effdof)
            fs=9
            ax.annotate('Ensemble mean',xy=(1.9,2.93),xycoords='data',fontsize=fs)
            ax.annotate('120-yr means', xy=(1.9,1.93),xycoords='data',fontsize=fs)
            ax.annotate('60-yr means', xy=(1.9,0.93),xycoords='data',fontsize=fs)
            ax.annotate('Individual SIC forcings',xy=(1.31,6.25),xycoords='data',fontsize=fs)
            ax.annotate('Average SIC forcing',xy=(1.31,3.25),xycoords='data',fontsize=fs)

            if 'ND' in seas:
                annlab='b. Nov-Dec Eurasian SAT changes'
            else:
                annlab='b. Eurasian SAT changes'
            annloc=(0.02,0.94)

            ax=axs[1] # second eurasian SAT
            sfnc.plot_uncertainty_cascade(dblob2,fdict,coords,sims,pparams,
                                          info=infodict,seas=seas,ax=ax,xlab=xlab,annlab=annlab,annloc=annloc,
                                          color=col,effdof=effdof)
            if printtofile:
                rg2 = infodict['region']
                if effdof:
                    fig.savefig(field + '_' + region + '_' + seas[0] + '_unccascade3_effdofobs.pdf')
                else:
                    fig.savefig(field + '_polcap60' + rg2 + '_' + seas[0] + '_unccascade7ann.pdf')
                    fig.savefig(field + '_polcap60' + rg2 + '_' + seas[0] + '_unccascade7ann.eps')

    elif plottype=='calccomposites':

        print 'here calculate high and low composites using a given threshold, say 1 sigma' # @@@@
        thresh=1 # the multiplier to sigma, ie 1sigma or 2sigma etc.
        # first do eurasia SAT anomalies.
        # assume I have a timeseries of anomalies, with a sigma associated with it.
        # select all times with val < sigma and with val > sigma
        #diffdf = pd.DataFrame(dblob['diff'])
        #std= diffdf.std()
        diffdt=dblob['diff']
        
        for sea in seasons:
            for skey in sims:
                diff = diffdt[skey][sea]
                std = diff.std()
                maxthresh = diff.mean()+thresh*std
                minthresh = diff.mean()-thresh*std

                high = ma.masked_where(diff<maxthresh,diff)
                low = ma.masked_where(diff>minthresh,diff)
                xx=np.arange(0,len(diff))
                fig,axs=plt.subplots(1,2)
                ax=axs[0] # timeseries
                ax.plot(xx,diff,color='k')
                ax.axhline(y=0,color='k',linewidth=.5)
                ax.axhline(y=maxthresh,color='k',linewidth=.5,linestyle='--')
                ax.axhline(y=minthresh,color='k',linewidth=.5,linestyle='--')
                ax.plot(xx,high,marker='o',linestyle='none',color='r')
                ax.plot(xx,low,marker='o',linestyle='none',color='b')
                #plt.fill_between(xx,y1=maxthresh,y2=diff,where=diff>=maxthresh,color='r')
                #plt.fill_between(xx,y1=minthresh,y2=diff,where=diff<=minthresh,color='b')
                ax.set_title(field + ' ' + region + ' ' + skey + ' ' + sea)

                ax=axs[1]
                nnhi,binshi,patcheshi = ax.hist(high.compressed(),color='r',alpha=0.5,normed=True,histtype='stepfilled')
                nnlo,binslo,patcheslo = ax.hist(low.compressed(),color='b',alpha=0.5,normed=True,histtype='stepfilled')
                nntot,binstot,patchestot = ax.hist(diff,color='0.5',alpha=0.5,normed=True,histtype='stepfilled')

    else:

        # Plot pvals for mean anomaly and standard dev
        # @@ for some reason if all sims are in dataframe, only some are plotted..??
        fig,axs=plt.subplots(1,2)
        fig.set_size_inches(12,4)
        ax=axs[0]
        tpvaldft.plot(linestyle='None',color=ccm.get_colordict().values(),marker='s',markersize=6,ax=ax)
        ax.axhline(y=siglevel,color='k')
        ax.set_ylabel('PVAL of T statistic for mean')
        ax.set_title(fdict['fieldstr'] + ' ' + region)
        ax.set_ylim((0,1))

        ax=axs[1]
        fpvaldft.plot(linestyle='None',colors=ccm.get_colordict().values(),marker='s',markersize=6,ax=ax)#make legend better@@
        ax.axhline(y=0.05,color='k')
        ax.set_ylabel('PVAL of F statistic for variance')
        ax.set_title(fdict['fieldstr'] + ' ' + region)
        ax.set_ylim((0,1))

        if printtofile:
            fig.savefig(fdict['fieldstr']+'_ftpvals' + savestr + '_' + region + '.pdf')

if plottype=='plotscatter':

    print 'Implement plotscatter! @@@'
    ## # get datablobs of two fields
    ## dblob = sfnc.calc_seasons(fdict,coords,sims,seas=seasons,
    ##                           loctimesel=timesel,info=infodict,calctype='regmeanwithtime',
    ##                           effdof=effdof)

    ## infodict['region'] = region2    
    ## dblob2 = sfnc.calc_seasons(fdict,coords,sims,seas=seasons,
    ##                            loctimesel=timesel,info=infodict,calctype='regmeanwithtime',
    ##                            effdof=effdof)

    ## sfnc.plot_scatter_regions()
        
if testhadisst:
    print '@@testhadisst not implemented'
