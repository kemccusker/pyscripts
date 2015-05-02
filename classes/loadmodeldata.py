""" loadmodeldata.py
        6/24/2014: I should standardize loading the simulation and variable
                   data since I do it over and over again.
"""
import numpy as np
import constants as con
import cccmaNC as cnc
import cccmautils as cutl

def loadmondata(fields,simulations,model='CanAM4',levsel=None,meantype=None):
    """ loadmondata(fields, simulations, model='CanAM4',levsel=None, meantype=None)
        Calls loaddate()
            fields: tuple of field names
            simulations: tuple of simulation names
            model: for now 'CanAM4' is implemented
            levsel: select level in Pa (e.g. 50000 for 500hPa)
            meantype: 'time','zonal' @@for now. default None, but recommended to choose
                       one if loading multiple variables and multiple simulations at once.
                       It is assumed that time is the first dimension.
                       If 'time' is the meantype, datadict contains tuples of (mean,std)


            returns nested dictionary@@
                    FIELDS->SIMULATIONS->MONTHS
                       
            Load requested MONTHLY fields from requested CanAM4 simulations
                 into dictionaries (dataframes?)
    """
    return loaddata(fields,simulations,model,timefreq='monthly',levsel=levsel,meantype=meantype)

def loadseasdata(fields,simulations,model='CanAM4',levsel=None,meantype=None):
    """ loadseasdata(fields, simulations, model='CanAM4',levsel=None, meantype=None)
        Calls loaddata()
            fields: tuple of field names
            simulations: tuple of simulation names
            model: for now 'CanAM4' is implemented
            levsel: select level in Pa (e.g. 50000 for 500hPa)
            meantype: 'time','zonal' @@for now. default None, but recommended to choose
                       one if loading multiple variables and multiple simulations at once.
                       It is assumed that time is the first dimension.
                       If 'time' is the meantype, datadict contains tuples of (mean,std)


            returns nested dictionary@@
                    FIELDS->SIMULATIONS->MONTHS
                       
            Load requested SEASONAL fields from requested CanAM4 simulations
                 into dictionaries (dataframes?)
    """
    return loaddata(fields,simulations,model,timefreq='seasonal',levsel=levsel,meantype=meantype)


def loaddata(fields, simulations, model='CanAM4',timeper='001-121',timefreq=None, 
             levsel=None,meantype=None,filetype='diff',region=None):
    """ loaddata(fields, simulations,model='CanAM4',timeper='001-121',timefreq=None, 
                 levsel=None, meantype=None,filetype='diff',region=None)
    
            fields: tuple of field names
            simulations: tuple of simulation names (diff names='E1'...'ENS','NSIDC' etc)
            model: for now only 'CanAM4' is implemented
            timeper: time period of the data (used for filename). default '001-121'
            timefreq: time frequency TO RETURN. default all data
                      'monthly'|'seasonal'|'climo'|'ANN'|'DJF'|'JJA'|'NDJ'|'MAM'|'SON'|
                      1,2,3,4,5,6,7,8,9,10,11,12 
            levsel: select level in Pa (e.g. 50000 for 500hPa). default all levels
            meantype: 'time','zonal' @@for now. default None, but recommended to choose
                       one if loading multiple variables and multiple simulations at once.
                       It is assumed that time is the first dimension.
                       If 'time' is the meantype, datadict contains tuples of (mean,std)
            filetype: 'diff','ctl','pert'
                       Default is diff where both ctl and pert are read in and differenced.
            region:   any of the regions in constants -> region dict. default None.

            returns: nested dictionary@@
                    FIELDS->SIMULATIONS->TIMEFREQ
                       
            Load requested fields from requested CanAM4 simulations
                 into dictionaries (dataframes?). Function automatically
                 skips first year of simulation and gets a timeseries (or climo
                 if requested) of the rest (assumed through year 121).
                 3D data and 'turb' not yet implemented! @@
    """
    if model!='CanAM4':
        print 'model not supported!'
        return -1

    print '@@ probably should invert the order such that it is field, season, sim?'
    #bp=con.get_basepath()
    #basepath=bp['basepath'] + model + '/'; subdir=bp['subdir']
    timesel='0002-01-01,0121-12-31'

    seabool=False # set to True if requested time freq is one season or climo
    monbool=False # set to True if requested time freq is one month

    if timefreq=='monthly':
        # get all months
        tf = con.get_mon()
    elif timefreq=='seasonal':
        # get all 4 seasons
        tf = 'DJF','MAM','JJA','SON'
    elif timefreq in range(1,13):
        # choose an individual month
        tf=timefreq
        monbool=True
    elif timefreq in ('climo','ANN','DJF','JJA','NDJ','MAM','SON'):
        tf=timefreq
        seabool=True

    print tf # @@
    
    datadict = dict.fromkeys(fields,{})
    for field in fields:

        print field # @@
        # assume simulations entered are of the form: E1, R3, ENSE etc. Then
        # filetype input arg tells which simulation to get (or both) @@@@ NOT DONE.
        simdict = dict.fromkeys(simulations,{})

        for sim in simulations:
            print sim # @@
            timdict = dict.fromkeys(tf)

            # construct filename here @@
            fnamec,fnamep = con.build_filepathpair(sim,field)
            #fname = basepath+sim+subdir+sim+'_'+field+'_'+timeper+'_ts.nc'
            print fnamec 
            
            tfkey = tf
            #for tfkey in tf:
            #    print tfkey # @@
            #print 'too many levels in the dict...get rid of seasonal keys and just do one@@@ 5/1/2015'
            # @@ get the data with cnc.getNCvar() here
            ncparams = {}
            if monbool:
                ncparams = {'monsel': tfkey}
            elif seabool:
                ncparams = {'seas': tfkey}

            if levsel!=None:
                ncparams['levsel'] = levsel

            if meantype=='zonal':
                ncparams['calc'] = 'zm'

            if filetype=='diff':
                fld = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,**ncparams) -\
                      cnc.getNCvar(fnamec,field.upper(),timesel=timesel,**ncparams)
            elif filetype=='ctl':
                fld = cnc.getNCvar(fnamec,field.upper(),timesel=timesel,**ncparams)
            elif filetype=='pert':
                fld = cnc.getNCvar(fnamep,field.upper(),timesel=timesel,**ncparams)
            else:
                print "filetype not supported! ['diff'|'ctl'|'pert']"
                return -1

            if region != None:
                lat=cnc.getNCvar(fnamec,'lat'); lon=cnc.getNCvar(fnamec,'lon')
                fld = cutl.calc_regmean(fld,lat,lon,region)

            if meantype=='time':
                #fldstd = np.std(fld,axis=0)
                fld = np.mean(fld,axis=0)
                #timdict[tfkey] = fld #,fldstd
                timdict=fld
            else:
                #timdict[tfkey] = fld
                timdict=fld

            simdict[sim] = timdict

        datadict[field]=simdict
        # can I set attributes to a dictionary? @@ like for nfields, nsims, ntimes?



    return datadict
