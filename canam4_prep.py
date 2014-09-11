""" canam4_prep.py:
      9/9/2014: Goal is to prep the data for 1. retrieval/calcs and 2. plotting
                To be used in conjunction with simulation_funcs.py
                Splitting canam4sims_analens.py
                This is particular to canam4, but the functions being called
                should not be (in the end...) so can have prep scripts for CanESM2, etc
                
"""
import constants as con
import simulation_funcs as sfnc

sfnc = reload(sfnc)
con = reload(con)

plt.close("all")
plt.ion()

# Will use this script to set up all the stuff.
# for example, plot_seasonal_maps() requires:
       # simpair keys
       # fielddict with 'field','ncfield','fieldstr'
       # pparams dict with cmin, cmax, cmap, latlim, type, suppcb ....
       #    if vert=True, also screen, levlim, addcontlines
       # coords = {'lat': con.get_t63lat(), 'lon': con.get_t63lon()}


printtofile=False

field = 'st'
smclim=False
level=50000 # for threed
nonstandardlev=False # standards are 700,500,300

# Choose type of plot =========================
seasonalmap=False # seasonal maps (SON, DJF, MAM, JJA)
seasonalvert=False # seasonal vertical zonal means instead of maps
screen=True # whether to have screen-style vertical zonal means

plotzonmean=False # plotzonmean,plotseacyc,pattcorrwithtime are mutually exclusive
plotseacyc=True # plotzonmean,plotseacyc,pattcorrwithtime are mutually exclusive
seacyclatlim=60 # southern limit for plotting polar mean seasonal cycles (line plot)
withlat=False # plot the seasonal cycle with latitude dimension too (only for plotseacyc=1)@@for now just std over ens
squatseacyc=False # plot seacycle figs as shorter than wide
squatterseacyc=True # even shorter, for paper

pattcorrwithtime=False # plot pattern correlation with time for each ens member
pattcorryr=False # if 1, do a yearly anomaly pattern rather than time-integrated

plotregmean=False
region = 'polcap60' # polcap60, polcap65, polcap70, eurasia, ntham, nthatl

testhadisst=0 # check which ens member most similar to hadisst

# Choose how to handle the data ==============
normbystd=False
pct = False # if 1, do calculation as a percent
halftime=False # get only the first 60yrs. make sure to set the other flag the opp
halftime2=False # get only the last 60yrs. make sure to set the other flag the opp

# Choose what simulations to add =============
sensruns=False # sensruns only: addr4ct=1,addsens=1. others=0 no meanBC, r mean, or obs
addobs=True # add mean of kemhad* & kemnsidc* runs to line plots, seasonal maps. 
addr4ct=False # add kem1pert2r4ct (constant thickness version of ens4)
addsens=False # add sensitivity runs (kem1pert1b, kem1pert3)
addrcp=True # add kem1rcp85a simulation (and others if we do more)
simsforpaper=False # meanBC, HAD, NSIDC only. best for maps and zonal mean figs (not line plots)

latlim = None # None #45 # lat limit for NH plots. Set to None otherwise. use 45 for BC-type maps
levlim= 100 # level limit for vertical ZM plots (in hPa). ignored if screen=True

# not sure these flags are in use?
sigtype = 'cont' # significance: 'cont' or 'hatch' which is default
sigoff=False # if True, don't add significance
siglevel=0.05
model='CanAM4'
# ??


# #################################################################
#   Probably don't need to modify below if everything goes well.
# #################################################################

# initialize things that don't get set otherwise
savestr='' # string for plot filenames
threed=False # is the field three dimensional
sia=False # is the requested field sea ice area
conv=1
timesel=None
isflux=False

# set up simulations and figure filename strings
sims = 'R1','R4','R3','R5','R2','ENS','CAN' # R's in order of sea ice loss

if simsforpaper: # best for maps only
    sims = ('HAD','NSIDC','CAN')
    savestr = '_forpap'
    seasons=('SON','DJF')
elif sensruns: # add sensitivity runs. with Shaded ENS. don't plot meanBC, mean of ens
    sims = sims[0:5] + ('R4ct','CANnosst','CANnothk')
    savestr = '_sensruns'
else:
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
    



# # # ######## set Field info ###################
# gz, t, u, v, q (3D !)
# st, sic, sicn (sia), gt, pmsl, pcp, hfl, hfs, turb, net, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)


print field

if halftime:
    timesel= '0002-01-01,0061-12-31'
elif halftime2:
    timesel='0062-01-01,0121-12-31'



# # # ###########################################
#   Shouldn't have to mod below....


"""  plev = 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, 12500, 
    15000, 17500, 20000, 22500, 25000, 30000, 35000, 40000, 45000, 50000, 
    55000, 60000, 65000, 70000, 75000, 77500, 80000, 82500, 85000, 87500, 
    90000, 92500, 95000, 97500, 100000 ;
"""

fdict = {'field': field, 'ncfield': None, 'fieldstr': None,
         'units': None, 'conv': conv,
         'nonstandardlev': nonstandardlev,
         'threed': threed} # fielddict

# reserved for expansion into the plotfunction call
pparams = {'cmin': None, 'cmax': None, 'cmap': 'blue2red_20',              
           'type':'nh', 'latlim': latlim} # plotparams

infodict ={'cmapclimo': 'Spectral_r','leglocs': None,
           'seacycylim': None, 'savestr': None,
           'model': model, 'sigtype': sigtype, 'sigoff': sigoff,
           'pct': pct, 'seacyclatlim': seacyclatlim} # random other info

if field == 'st':
    fdict['units'] = 'K'
    fdict['ncfield'] = field.upper()
    fdict['fieldstr'] = field

    pparams['cmin'] = -3; pparams['cmax'] = 3 # seasonal/monthly
    if smclim:
        pparams['cmin'] = -1.5; pparams['cmax'] = 1.5 # seasonal/monthly
        savestr= savestr + '_smclim'

    pparams['cmap'] = 'blue2red_w20'
        
    #cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
    #cminn = -5; cmaxn = 5 # for norm by std

    if plotseacyc  and withlat:
        pparams['cmin']=-.5;
        pparams['cmax']=.5 # @@ will have to update this when add subplots
        
    leglocs = 'upper left', 'upper left', 'upper right', 'upper left'
    seacycylim=(-.5,4) # >70N
    
elif field == 'sic':

    fdict['units'] = 'm'
    fdict['ncfield'] = field.upper()
    fdict['fieldstr'] = field

    pparams['cmin'] = -.5; pparams['cmax'] = .5 # seasonal/monthly

    fdict['conv']=1/913.

    pparams['cmap'] = 'red2blue_w20'
    #pparams['leglocs'] = 'lower left', 'lower left', 'upper left', 'upper left'
    
elif field == 'sicn' or field == 'sia':
    units = 'frac'
    conv=1
    cmin=-.15; cmax=.15
    cminm=-.15; cmaxm=.15
    cmap = 'red2blue_w20'
    leglocs = 'lower left', 'lower left', 'upper left', 'upper left'
    seacycylim=(-2e12,0) # for sia
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

    if plotseacyc  and withlat:
        cminm=-1; cmaxm=1 # @@ will have to update this when add subplots
    leglocs = 'lower left', 'lower left', 'upper center', 'lower left'
    seacycylim=(-2,1.5) # >70N
elif field == 'pcp':
    units = 'mm/day' # original: kg m-2 s-1
    
    #pct=1; units = '%'; print 'PCT'
    
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

    isflux=True
    
elif field == 'hfs': # sfc upward SH flux
    units = 'W/m2'
    conv = 1
    cmin = -5
    cmax = 5
    cminm = -8
    cmaxm = 8

    isflux=True
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
    seacycylim=(-5,25) # always >40N (where there is ice in CTL)

    isflux=True
elif field == 'flg': # net downward LW at the sfc.
    units = 'W/m2'
    conv = -1 # so positive heats atmos
    cmin = -5
    cmax = 5
    cminm = -8
    cmaxm = 8
    leglocs = 'upper left', 'lower left', 'upper left', 'upper left'

    isflux=True
elif field == 'fsg': # net (absorbed) solar downard at sfc
    units = 'W/m2'
    conv = 1
    cmin = -5
    cmax = 5
    cminm = -8
    cmaxm = 8

    isflux=True

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
    threed = True

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
    threed = True

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
    threed=True

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
        seacycylim=(-16,20) # >70N, 500hPa

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

coords = {'lev': con.get_t63lev(), 'lat': con.get_t63lat(), 'lon': con.get_t63lon()}
infodict['savestr'] = savestr
infodict['leglocs'] = leglocs
infodict['seacycylim'] = seacycylim
fdict['isflux'] = isflux
fdict['threed'] = threed


# do an if elif elif ....

if seasonalmap or seasonalvert:

    print fdict
    print sims
    print pparams

    if seasonalvert:
        pparams['screen']=screen
        
    # this one does data processing and plotting together
    # some stuff in the function need to be removed or set differently.@@
    # marked in the function @@
    sfnc.plot_seasonal_maps(fdict,coords,sims,pparams,vert=seasonalvert,loctimesel=timesel,info=infodict,seas=seasons)
    

elif plotseacyc:

    print 'test me @@'
    dblob = sfnc.calc_seasonal_cycle(fdict,coords,sims,withlat=withlat,loctimesel=timesel,info=infodict)

    #sfnc.plot_seasonal_cycle(dblob,fdict,sims,ptypes=('climo','anom'),info=infodict)
    

elif plotzonmean or pattcorrwithtime or plotregmean:

    # call a func that isn't there yet
    print '@@ not implemented'

elif testhadisst:
    print '@@testhadisst not implemented'
