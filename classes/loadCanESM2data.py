import numpy as np
import constants as con
import cccmaNC as cnc
import cccmautils as cutl
import scipy as sp
import scipy.signal 
import matplotlib.pyplot as plt

# from scipy import signal

""" This class is actually for the CanESM2 LE. Not general runs. CMIP style.

"""
bp=con.get_basepath(runtype='CMIP')

def get_timepers(casename):

    if casename=='piControl':
        timepers = ('201501-231012','231101-241012','241101-251012','251101-261012',
            '261101-271012','271101-281012','281101-291012','291101-301012')
    elif casename=='gregory_2xco2':
        timepers = ('245101-255012',)

    return timepers

def build_filenames(fielddict, casename, ftype=None,timesel=None,verb=True,local=False):
    """ ftype and timesel are not implemented here yet.

        ftype: type of filename to build. Right now just 'fullts' for full timeseries
              or 'fullclimo' for '1950-2020_climo' or, given timesel:
              for styr-enyr_climo. or 'ensmean'
        local: is the data on /raid/ra40 or ~/work/DATA (local=True)
              
        returns a list of all filenames 
        
    """
    print 'local ' + str(local)
    
    # for now, don't need these flags
    """if ftype=='fullts':
        suff='195001-202012'
    elif ftype=='ensmean':
        # assume fullts
        suff='195001-202012'
        ensmean=True
    elif ftype=='fullclimo':
        suff=ftype
    else: # climo with specified years
        if timesel==None:
            print 'timesel cannot be None if ftype=climo'
            return -1
        (timselst,timselen)=timesel.split(',')       
        (styear,stmon,stday) = timselst.split('-')
        (enyear,enmon,enday) = timselen.split('-')
        suff=str(styear)+ '-' + str(enyear) + '_climo'
    """

    if local: # as of 10/7/2016, no consistent local directory (ie on laptop or imac itself. just external drives)
        basepath = '/HOME/rkm/work/DATA/CanESM2/'
    else:
        print 'getting basepath'
        basepath=bp['basepath']

        
    field=fielddict['field']
    comp=fielddict['comp']
    
    basedir=basepath + casename

    flist=[]

    timepers = get_timepers(casename)

    for timeper in timepers:

        fname=basedir + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
           casename + '_r1i1p1_' + timeper + '.nc'
        if verb:
            print fname
        flist.append(fname)

    return flist


def load_data(fielddict, casename, seas=None, timesel=None,infodict=None,ftype='fullts',
              calctype=None, calcdict=None,conv=1, region=None,local=False,
              orig=None,verb=True, detrend=True):
    """ def load_data(fielddict, seas=('DJF','MAM','JJA','SON'), timesel=None, infodict=None, calctype=None, calcdict=None)

            fielddict: should have {'field': field, 'ncfield': ncfield, 'comp': comp} at minimum
            ftype: type of filename to build. Right now just 'fullts' for full timeseries
              or '1950-2020_climo' or 'ensmean'
            seas has to be a tuple
            rettype: 'dict' or 'ndarray' as return type. 
                      default is dict of DataFrames (to make a Panel). 
                      else, 'ndarray' is a 4D array
            local: is data in ~/work/DATA (local=True) or /raid/ra40
            detrend: removes any spurious trend (drift) found in data, for e/ grid point 
                     

            @@ add regional avg functionality?
            @@ region, calctype, calcdict? not implemented

            returns an object of type rettype
    """

    field = fielddict['field']
    ncfield=fielddict['ncfield']
    comp = fielddict['comp']

    flist = build_filenames(fielddict, casename,ftype=ftype,timesel=timesel,local=local,verb=verb)


    for fii,fname in enumerate(flist):

        if seas==None:
            # don't average, return all months
            fld = cnc.getNCvar(fname,ncfield,timesel=timesel)*conv
            #if region!=None:
            #    lat=cnc.getNCvar(fname,'lat')
            #    lon=cnc.getNCvar(fname,'lon')
            #    fld = cutl.calc_regmean(fld,lat,lon,region=region)
            
            # do calcs:
            if calctype!=None:
                print 'calcs not implemented yet!! @@@'

            fldseas=fld

            
        else:
            sea = seas[0]
            fld = cnc.getNCvar(fname,ncfield,timesel=timesel,seas=sea)*conv

            # do calcs
            if calctype!=None:
                print 'calcs not implemented yet!! @@@'
                
            fldseas=fld
            
        # hstack for 1D
        if fii==0:
            fldret=fldseas
        else:
            if fldseas.ndim==1:
                fldret = np.hstack((fldret,fldseas))
            else:
                fldret = np.vstack((fldret,fldseas))
        
    if detrend and seas!=None: # @@@ eventually implement for monthly too
        flddtr = sp.signal.detrend(fldret,axis=0,type='linear') + fldret.mean(axis=0)

        fldret = flddtr

    return fldret
