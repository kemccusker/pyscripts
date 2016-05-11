import cccmautils as cutl
import pandas as pd
import cccmaNC as cnc
import constants as con
import numpy as np

""" as loadLE.py but specifically for CESM1 (instead of CanESM2
    3/24/2016

"""
con=reload(con)

bp=con.get_LEbasepath(model='CESM1')

basepath=bp['basepath']
ensnum=30
prefix = 'b.e11'
res = 'f09_g16'

#  Make this a class! CESMLE

def get_hh(comp):

    if comp=='cam':
        return 'h0'
    else:
        raise Exception('loadCESMLE.get_hh(): Comp other than cam not implemented!')
    

def load_LEdata(fielddict, ens, seas=None, timesel=None,infodict=None,ftype='fullts',
                calctype=None, calcdict=None, rettype='dict',conv=1, region=None,local=False,
                subens=False,verb=True,zonal=False):
    """

            ftype: type of filename to build. Default 'fullts' for full timeseries (assumes 192001-202012)
                   Also can override default time period of file by setting ftype to timeperiod.
                   e.g. ftype = '192001-200512'

            seas has to be a tuple (but will only get one season)
            rettype: 'dict' or 'ndarray' as return type. 
                      default is dict of DataFrames (to make a Panel). 
                      else, 'ndarray' is a 4D array
            local: is data in ~/work/DATA (local=True) or /raid/rc40
          
            subens: overrides getting all ensemble members. Must be tuple (e.g. (5,10) get ens mems 5 thru 9)

            @@ add regional avg functionality?

            returns an object of type rettype
    """

    field = fielddict['field']
    ncfield=fielddict['ncfield']
    comp = fielddict['comp']

    flist = build_filenames(fielddict, ens,ftype=ftype,timesel=timesel,local=local,
                            subens=subens,verb=verb)

    fname1 = flist[0]
    if verb:
        print ' @@ fname1 ' + fname1
        
        
    if rettype=='dict':
        fldret={}
    elif rettype=='ndarray':
        if zonal:
            tmp = cnc.getNCvar(fname1,ncfield,timesel=timesel,calc='zm')*conv
        else:
            tmp = cnc.getNCvar(fname1,ncfield,timesel=timesel)*conv

        if len(tmp.shape)>3:
            print '4D variable not supported! @@@@'
            return -1
        elif len(tmp.shape)==3:
            # first dim is time, reshape the 2nd two into one spatial dim
            n0=tmp.shape[0]
            n1=tmp.shape[1]
            n2=tmp.shape[2]
            initshape=(len(flist), n0, n1*n2) # 3D, ens x time x space
            rshape=initshape[1:] # what to reshape individual data to
        elif len(tmp.shape)==2:
            initshape=list(tmp.shape)
            initshape.insert(0,len(flist))
            initshape=tuple(initshape) # 3D, ens x time x space
            rshape=None
        elif len(tmp.shape)==1:
            initshape = (len(flist),tmp.shape[0])
            rshape=None

        print initshape
        fldret = np.zeros(initshape)
        print fldret.shape # @@@
    else:
        raise Exception('return type not supported!! @@@')
        return -1

    numens = len(flist)
    #fldall = np.zeros(numens)
    #dlist=[]
    

    if rettype=='dict':
        timedim = cnc.getNCvar(fname1,'time',timesel=timesel)
        if verb: print timedim


    for eii,fname in enumerate(flist):
        
        if seas==None:
            # don't average, return all months
            if zonal:
                fld = cnc.getNCvar(fname,ncfield,timesel=timesel,calc='zm')*conv
            else:
                fld = cnc.getNCvar(fname,ncfield,timesel=timesel)*conv
            #if region!=None:
            #    lat=cnc.getNCvar(fname,'lat')
            #    lon=cnc.getNCvar(fname,'lon')
            #    fld = cutl.calc_regmean(fld,lat,lon,region=region)
            
            # do calcs:
            if calctype!=None:
                print 'calcs not implemented yet!! @@@'

            if rettype=='ndarray':
                fldseas=fld.reshape(rshape)
            else:
                if len(fld.shape)!=3:
                    raise Exception('load CESM LE with rettype dict will break if dim of var is not 3 @@@')
                    return -1

                n0,n1,n2=fld.shape
                fldseas=fld.reshape((n0,n1*n2)) # flattening
            
        else:
            fldseas={}
            for sea in seas:
                if zonal:
                    # @@@ not expecting this to work.
                    fld = cnc.getNCvar(fname,ncfield,timesel=timesel,seas=sea,calc='zm')*conv
                else:
                    fld = cnc.getNCvar(fname,ncfield,timesel=timesel,seas=sea)*conv

                # do calcs
                if calctype!=None:
                    print 'calcs not implemented yet!! @@@'
                
                fldseas[sea]=fld
            
        if rettype=='dict':
            fldret[fname]= pd.DataFrame(fldseas,index=timedim)
        elif rettype=='ndarray':
            if seas!=None:
                print 'ndarray return type with a seasonal avg does not work yet @@@@'
            fldret[eii,...] = fldseas

        #dlist.append(fldseas)

        
    return fldret



def build_filenames(fielddict, ens, ftype='fullts',timesel=None,verb=True,local=False,subens=False):
    """ 

        ens: name of runs: B20TRC5CNBDRD-BRCP85C5CNBDRD, B20TRC5CNBDRD, or BRCP85C5CNBDRD

        ftype: type of filename to build. Right now just 'fullts' for full timeseries
              or 'fullclimo' for '1950-2020_climo' or, given timesel:
              for styr-enyr_climo. or 'ensmean'. 
              Also can override default time period by setting ftype to timeperiod.
                   e.g. ftype = '195001-201012'
        local: is the data on /raid/rc40 or ~/work/DATA (local=True)
              
        subens: overrides getting all ensemble members. Must be tuple of ens nums (e.g. (5,10) to get ens 5 thru 9.)

        returns a list of all filenames.
        
    """

    ensmean=False
    if ftype=='fullts':
        if ens=='B20TRC5CNBDRD-BRCP85C5CNBDRD':
            suff='192001-202012'
        elif ens=='B20TRC5CNBDRD':
            suff='192001-200512'
        elif ens=='BRCP85C5CNBDRD':
            suff='200601-208012'
        else: raise Exception('Ensemble name not recognized')

    #elif ftype=='ensmean':
    #    # assume fullts
    #    suff='195001-202012'
    #    ensmean=True
    #elif ftype=='fullclimo':
    #    suff=ftype
    #elif 'climo' in ftype: # @@ not the right test. want to test for substr@@@ # climo with specified years
    #    if timesel==None:
    #        print 'timesel cannot be None if ftype=climo'
    #        return -1
    #    (timselst,timselen)=timesel.split(',')       
    #    (styear,stmon,stday) = timselst.split('-')
    #    (enyear,enmon,enday) = timselen.split('-')
    #    suff=str(styear)+ '-' + str(enyear) + '_climo'
    else:
        # assume ftype is a time period -- add Exception here @@@@@
        suff = ftype # a way to override default timeperiod


    if local:
        basepath = '/HOME/rkm/work/DATA/cesm1/'
    else:
        basepath=bp['basepath']

        
    field=fielddict['field']
    comp=fielddict['comp']
    
    flist=[]

    #if ensmean: # @@ note, should add the casename to ensmean filename!
    #    fname = basepath + ens + '-ens' + '/' + field + '/' + field + '_' + comp + '_CanESM2_ensmean_' +\
    #            suff + '.nc'
    #    if verb:
    #        print fname
    #    flist.append(fname)

    #else:
    if 1:
        if not subens:
            ensrange = range(1,ensnum+1)
        else:
            ensrange = range(subens[0],subens[1])

        hh = get_hh(comp)

        for eii in ensrange:
            if suff=='192001-200512' and eii==1: # first ens member of 20thC non-processed files is diff time per
                fname=basepath + '/' + field + '/' + prefix + '.' + ens + '.' + res + '.00' + str(eii) + \
                       '.' + comp + '.' + hh + '.' + field + '.185001-200512.nc'
            else:
                if eii<10:
                    # b.e11.B20TRC5CNBDRD-BRCP85C5CNBDRD.f09_g16.002.cam.h0.TREFHTeurasiamori.192001-202012.nc
                    fname=basepath + '/' + field + '/' + prefix + '.' + ens + '.' + res + '.00' + str(eii) + \
                           '.' + comp + '.' + hh + '.' + field + '.' + suff + '.nc'
                else:
                    fname=basepath + '/' + field + '/' + prefix + '.' + ens + '.' + res + '.0' + str(eii) + \
                           '.' + comp + '.' + hh + '.' + field + '.' + suff + '.nc'


            if verb:
                print fname
            flist.append(fname)

    return flist
        
def get_lat(local=False):
    
    fnames=build_filenames({'field':'Z3500','comp':'cam'},'B20TRC5CNBDRD-BRCP85C5CNBDRD',verb=False,local=local)
    return cnc.getNCvar(fnames[1],'lat')

def get_lon(local=False):
    
    fnames=build_filenames({'field':'Z3500','comp':'cam'},'B20TRC5CNBDRD-BRCP85C5CNBDRD',verb=False,local=local)
    return cnc.getNCvar(fnames[1],'lon')



    
