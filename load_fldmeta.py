# the point here is to set up the info dict etc that
#   has all the field metadata in it for plotting and such
# basically strip that stuff from canam4_prep.py
#   11/20/2014
import constants as con

def loadfldmeta(field,infodict,plottype,ptparams,level=None):
    """ NOTE not all Fields are implemented into fdict and infodict, etc!
        11/20/2014
    """

    # check for nonstandardlev
    # set plottype
    # ptparams has specific settings for plottype

    print ptparams
    
    smclim=ptparams['smclim']
    latlim=ptparams['latlim']
    levlim=ptparams['levlim']
    region = ptparams['region'] 
    savestr=infodict['savestr'] # string for plot filenames
    seacycylim=infodict['seacycylim']
    
    threed=False # is the field three dimensional
    sia=False # is the requested field sea ice area
    conv=1
    isflux=False

    if level!=None:
        nonstandardlev = con.is_standardlev(level)
    else:
        nonstandardlev=False


    # initialize flags
    seasonalmap = seasonalvert = plotzonmean = plotseacyc = pattcorrwithtime = plotregmean = False

    if plottype=='seasonalmap':
        seasonalmap=True
    elif plottype=='seasonalvert':
        seasonalvert=True
        screen=ptparams['screen']
    elif plottype=='plotzonmean':
        plotzonmean=True
    elif plottype=='plotseacyc':
        plotseacyc=True
        seacyclatlim=ptparams['seacyclatlim']
        withlat = ptparams['withlat']

    elif plottype=='pattcorrwithtime':
        pattcorrwithtime=True
        pattcorryr = ptparams['pattcorryr']
    elif plottype=='plotregmean':
        plotregmean=True
    else:
        print 'Plottype not recognized!'
        return


    ############################ start copy to load_fldmeta.py
    # # # ######## set Field info ###################
    # gz, t, u, v, q (3D !)
    # st, sic, sicn (sia), gt, pmsl, pcp, hfl, hfs, turb, net, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)


    print field

    ## if halftime:
    ##     timesel= '0002-01-01,0061-12-31'
    ## elif halftime2:
    ##     timesel='0062-01-01,0121-12-31'



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
    ## seacycylim=None
    ## infodict ={'cmapclimo': 'Spectral_r','leglocs': None,
    ##            'seacycylim': None, 'savestr': None,
    ##            'model': model, 'sigtype': sigtype, 'sigoff': sigoff,
    ##            'pct': pct, 'seacyclatlim': seacyclatlim, 'region': region,
    ##            'shadeens': shadeens, 'corrlim': corrlim} # random other info

    if field == 'st':
        fdict['units'] = 'K'
        fdict['ncfield'] = field.upper()
        fdict['fieldstr'] = field

        pparams['cmin'] = -3; pparams['cmax'] = 3 # seasonal/monthly
        pparams['cmap'] = 'blue2red_w20'
        if smclim:
            pparams['cmin'] = -1.5; pparams['cmax'] = 1.5 # seasonal/monthly
            pparams['cmap'] = 'blue2red_20'
            savestr= savestr + '_smclim'

        #cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
        #cminn = -5; cmaxn = 5 # for norm by std

        if plotseacyc  and withlat:
            pparams['cmin']=-.5;
            pparams['cmax']=.5 # @@ will have to update this when add subplots

        leglocs = 'upper left', 'upper left', 'upper right', 'upper left'
        if region not in ('eurasia','ntham','nthatl'):
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
        fdict['units'] = 'frac'
        fdict['ncfield'] = field.upper()
        fdict['fieldstr'] = field

        pparams['cmin'] = -.15; pparams['cmax'] = .15 # seasonal/monthly

        pparams['cmap'] = 'red2blue_w20'
        leglocs = 'lower left', 'lower left', 'upper left', 'upper left'

        if field=='sia':
            fdict['ncfield'] = 'SICN'
            seacycylim=(-2e12,0) # for sia
            sia=True

    elif field == 'gt':
        units='K'
        conv=1
        cmin=-2; cmax=2
        cminm=-3; cmaxm=3
        cminp=-.5; cmaxp=.5 # for when pert is 'ctl'
        cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
        cmap = 'blue2red_w20'
    elif field == 'pmsl':
        fdict['units'] = 'hPa'
        fdict['ncfield'] = field.upper()
        fdict['fieldstr'] = field

        pparams['cmin'] = -2; pparams['cmax'] = 2 # seasonal/monthly
        if smclim:
            pparams['cmin'] = -1; pparams['cmax'] = 1 # seasonal/monthly
            savestr= savestr + '_smclim'

        pparams['cmap'] = 'blue2red_20'

        if plotseacyc  and withlat:
            pparams['cmin']=-1;
            pparams['cmax']=1 # @@ will have to update this when add subplots

        leglocs = 'lower left', 'lower left', 'upper center', 'lower left'
        if region not in ('eurasia','ntham','nthatl'):
            seacycylim=(-2,1.5) # >70N


        ## units = 'hPa' # pretty sure hpa @@double check
        ## conv = 1
        ## cmin = -1; cmax = 1  # for anomaly plots
        ## cminm=-2; cmaxm=2  # for monthly maps
        ## cminp=cmin; cmaxp=cmax # for when pert is 'ctl'
        ## cminmp=cminm; cmaxmp=cmaxm
        ## cmap = 'blue2red_20'
        ## ## print 'new cmap and small clim! @@ '
        ## ## cmap = 'blue2red_w20'
        ## ## cminm=-1; cmaxm=1  for monthly maps, small clim

        ## cminn = -1; cmaxn = 1 # for norm by std


    elif field == 'pcp':
        fdict['fieldstr'] = field
        fdict['ncfield'] = field.upper()

        fdict['units'] = 'mm/day' # original: kg m-2 s-1
        fdict['conv'] = 86400 # convert from kg m-2 s-1 to mm/day

        pparams['cmap'] = 'brown2blue_16w'

        if pct:
            pparams['cmin'] = -20; pparams['cmax'] = 20
            fdict['units'] = '%'
        else:
            pparams['cmin'] = -0.2; pparams['cmax'] = 0.2

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

        ## #cmap = 'PuOr'
        ## cmap = 'brown2blue_16w'
        ## cminpct=-12; cmaxpct=12
        ## cminmpct=-20; cmaxmpct=20
        ## cminmp =-.25; cmaxmp=.25
        ## cminpctp=-8; cmaxpctp=8
        ## cminpctmp=-12; cmaxpctmp=12
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
        fdict['units'] = 'm/s'
        fdict['ncfield'] = 'SU'    
        fdict['conv'] = 1;

        pparams['cmap'] = 'blue2red_20'
        pparams['cmin'] = -1; pparams['cmax'] = 1
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

        fdict['units'] = 'K'
        fdict['ncfield'] = 'TEMP'

        threed = True
        fdict['conv'] = 1

        pparams['cmap'] = 'blue2red_w20'

        if level == 30000:
            pparams['cmin'] = -.5; pparams['cmax'] = .5
            #cminm = -.5; cmaxm = .5  # for monthly
            #cminsea = -.5; cmaxsea = .5
        elif level == 70000:
            #cmin = -.3; cmax = .3
            pparams['cmin'] = -.5; pparams['cmax'] = .5  # for monthly
            #cminsea = -.5; cmaxsea = .5

        if seasonalvert:
            if screen:
                pparams['levlim']=300
                pparams['cmin']=-2.5; pparams['cmax']=2.5
                #cminm=-2.5; cmaxm=2.5
            else:
                cmin = -.5; cmax = .5
                pparams['cmin'] = -.8; pparams['cmax'] = .8 
    elif field == 'u':
        threed = True
        fdict['units'] = 'm/s'
        fdict['ncfield'] = 'U'
        fdict['conv'] = 1

        if level == 30000:
            #cmin = -2; cmax = 2
            pparams['cmin'] = -3; pparams['cmax'] = 3
            #cminsea = -3; cmaxsea = 3
        else:
            #cmin = -1; cmax = 1
            pparams['cmin'] = -1; pparams['cmax'] = 1
            #cminsea = -1; cmaxsea = 1

        if seasonalvert:
            #cmin=-.5; cmax=.5
            pparams['cmin'] = -1; pparams['cmax'] = 1
            if screen:
                pparams['levlim']=300

        pparams['cmap']='blue2red_20'

    elif field == 'gz':

        fdict['units'] = 'm'
        fdict['ncfield'] = 'PHI'

        threed = True
        fdict['conv'] = 1/con.get_g()

        if level==30000:
            pparams['cmin'] = -20; pparams['cmax'] = 20 # seasonal/monthly
        else:
            pparams['cmin'] = -15; pparams['cmax'] = 15 # seasonal/monthly
            if region not in ('eurasia','ntham','nthatl'):
                seacycylim=(-16,20) # >70N, 500hPa

        if seasonalvert:
            if screen:
                #cmin=-10; cmax=10
                #cminm=-25; cmaxm=25
                pparams['cmin'] = -25; pparams['cmax'] = 25 # seasonal/monthly
                pparams['levlim'] = 300
            else:
                #cmin = -10; cmax = 10
                #cminm = -15; cmaxm = 15
                pparams['cmin'] = -15; pparams['cmax'] = 15 # seasonal/monthly

        pparams['cmap'] = 'blue2red_w20'

        leglocs = 'upper left', 'upper left', 'upper right', 'upper left'

    else:
        print 'No settings for ' + field

    if threed: # either map of a level, or a zonal mean
        if seasonalvert:
            fieldstr=field+'ZM'
            field=field+'ZM'
        else:
            fieldstr=field+str(level/100) # figure names. want integer for string
            field=field+str(level) # read in files
        fdict['field']=field
        fdict['fieldstr']=fieldstr
    else:
        fdict['fieldstr']=field

    fluxes = 'hfl','hfs','flg' # LH, SH, LWdown

    infodict['savestr'] = savestr
    infodict['leglocs'] = leglocs
    infodict['seacycylim'] = seacycylim
    fdict['isflux'] = isflux
    fdict['threed'] = threed
    ###################### end copy to load_fldmeta.py

    return fdict,pparams
