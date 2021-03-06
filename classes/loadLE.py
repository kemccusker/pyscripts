import cccmautils as cutl
import pandas as pd
import cccmaNC as cnc
import constants as con
import numpy as np

con=reload(con)

bp=con.get_LEbasepath()

basepath=bp['basepath']
ensnum=10

ensDT = {'historical':'historical',
         'historicalNat': 'historicalNat',
         'historicalAA': 'historicalMisc',
         'historicalSO': 'historicalMisc',
         'historicalLU': 'historicalMisc',
         'historicalSI': 'historicalMisc',
         'historicalGHG': 'historicalGHG'}
pnumDT = {'historical':'1', 'historicalNat': '1',
          'historicalAA': '4',
          'historicalSO':'6', 'historicalLU': '2',
          'historicalSI': '3','historicalGHG':'1'}
timeperDT={'le':'195001-202012',
           'cmip':'185001-201212',
           'historicalSO': '195001-200012'} # except for some reason SO is '195001-200012', '200101-205012'

#  Make this a class! CanESM2LE

def load_LEdataxr(fielddict, ens, seas=None, timesel=None,infodict=None,ftype='fullts',
                calctype=None, calcdict=None, rettype='dict',conv=1, region=None,local=True,
                orig=None,subens=False,verb=True,zonal=False):
    """ load LE data using xarray
    """

    print 'load_LEdataxr() not yet implemented' #@@@
    

    
def load_LEdata(fielddict, ens, seas=None, timesel=None,infodict=None,ftype='fullts',
                calctype=None, calcdict=None, rettype='dict',conv=1, region=None,local=True,
                orig=None,subens=False,verb=True,zonal=False):
    """ def loadLEdata(fielddict, seas=('DJF','MAM','JJA','SON'), timesel=None, infodict=None, calctype=None, calcdict=None)

            ftype: type of filename to build. Default 'fullts' for full timeseries.
              'fullclimo' for '1950-2020_climo' or, given timesel:
              for styr-enyr_climo. or 'ensmean'. 
              Also can override default time period of file by setting ftype to timeperiod.
                   e.g. ftype = '195001-201012'
            seas has to be a tuple
            rettype: 'dict' or 'ndarray' as return type. 
                      default is dict of DataFrames (to make a Panel). 
                      else, 'ndarray' is a 4D array
            local: is data in /Volumes/KellyDataDisk/home/work/DATA (local=True) or /raid/ra40
            orig='just' means just load the original 5 (historical+rcp85) runs (consistent w/ LE). 
            orig='just45' means just load original 5 (historical+rcp45) to the full LE.

            orig='add' means add the original 5 (historical+rcp85) to the full LE (consistent w/ LE)
            orig='add45' means add the original 5 (historical+rcp45) to the full LE.
          
            subens: overrides getting all sim groups in ensemble. must be tuple. e.g. ('historical-r1',)

            @@ add regional avg functionality?

            returns an object of type rettype
    """

    field = fielddict['field']
    ncfield=fielddict['ncfield']
    comp = fielddict['comp']

    flist = build_filenames(fielddict, ens,ftype=ftype,timesel=timesel,local=local,
                            orig=orig,subens=subens,verb=verb)

    fname1 = flist[0]
    if verb:
        print ' @@ fname1 ' + fname1
    timedim = cnc.getNCvar(fname1,'time',timesel=timesel)
    #print timedim

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
        print 'return type not supported!! @@@'
        return -1

    numens = len(flist)
    #fldall = np.zeros(numens)
    #dlist=[]
    
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
                    print 'load LE with rettype dict will break if dim of var is not 3 @@@'
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

def load_originalfive(fielddict, ens,seas=None, timesel=None,infodict=None,ftype='fullts',
                      calctype=None, calcdict=None, rettype='dict',conv=1, region=None,local=False,
                      orig='just',verb=True):
    """ orig='just' means just load the original 5 (historical+rcp85) runs (consistent w/ LE). 
        orig='just45' means just load original 5 (historical+rcp45) to the full LE.

        orig='add' means add the original 5 (historical+rcp85) to the full LE (consistent w/ LE)
        orig='add45' means add the original 5 (historical+rcp45) to the full LE.
    """

    return load_LEdata(fielddict, ens, seas=seas, timesel=timesel, infodict=infodict,
                       ftype=ftype, calctype=calctype, calcdict=calcdict, rettype=rettype,
                       conv=conv, region=region, local=local,orig=orig,verb=verb)
    


def build_filenames(fielddict, ens, ftype='fullts',timesel=None,verb=True,local=False,orig=None,
                        subens=False,etype='le'):
    """ here we know that each 'sim' has 10 sub-ensemble members

        ftype: type of filename to build. Right now just 'fullts' for full timeseries
              or 'fullclimo' for '1950-2020_climo' or, given timesel:
              for styr-enyr_climo. or 'ensmean'. 
              Also can override default time period by setting ftype to timeperiod.
                   e.g. ftype = '195001-201012'
        local: is the data on /raid/ra40 or /Volumes/KellyDataDisk/home/work/DATA (local=True)
              
        orig='just' means just load the original 5 (historical+rcp85) runs (consistent w/ LE). 
        orig='just45' means just load original 5 (historical+rcp45) to the full LE.

        orig='add' means add the original 5 (historical+rcp85) to the full LE (consistent w/ LE)
        orig='add45' means add the original 5 (historical+rcp45) to the full LE.

        subens: overrides getting all sim groups in ensemble. must be tuple. e.g. ('historical-r1',)

        etype: 'le', 'cmip': le are the 50-member large ensembles. cmip are the 5-member cmip-submitted ensembles.
               'cmip' w/ ens='historical' is the same as orig='just'
        returns a list of all filenames (50).
        
    """

    ensmean=False
    if ftype=='fullts':
        if ens=='historicalSO':
            suff='195001-200012'
        else:
            suff='195001-202012'
    elif ftype=='ensmean':
        # assume fullts
        suff='195001-202012'
        ensmean=True
    elif ftype=='fullclimo':
        suff=ftype
    elif 'climo' in ftype: # @@ not the right test. want to test for substr@@@ # climo with specified years
        if timesel==None:
            print 'timesel cannot be None if ftype=climo'
            return -1
        (timselst,timselen)=timesel.split(',')       
        (styear,stmon,stday) = timselst.split('-')
        (enyear,enmon,enday) = timselen.split('-')
        suff=str(styear)+ '-' + str(enyear) + '_climo'
    else:
        # assume ftype is a time period -- add Exception here @@@@@
        suff = ftype # a way to override default timeperiod

    #if ens in ('historicalMisc', 'historicalAA'):
    #    pnum='4'
    #elif ens in ('historicalSO',):
    #    pnum='6'
    #else: pnum='1'
    pnum = pnumDT[ens]

    if local:
        basepath = '/Volumes/KellyDataDisk/home/work/DATA/CanESM2/'
    else:
        basepath=bp['basepath']

    if etype=='le': basepath=basepath+'LE/'
        
    field=fielddict['field']
    comp=fielddict['comp']
    
    flist=[]
    if etype=='cmip' or (orig in ('just','just45')):

        #try:
        if ens=='historical':
            if orig=='just45':
                casename='historicalrcp45'
                timeext='185001-201212'
            elif orig=='just':
                casename='historicalrcp85'
                timeext='185001-202012'
            else:
                # etype=cmip
                casename='historical'
                timeext='185001-201212'
        else:
            casename=ensDT[ens]
            timeext=timeperDT[etype]
                #raise
        #except:
        #    raise Exception

        orignum=5
        basedir=basepath + casename

        for eii in np.arange(1,orignum+1):
            fname=basedir + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
                   casename + '_r' + str(eii) + 'i1p' + pnum + '_' + timeext + '.nc'
            if verb:
                print fname
            flist.append(fname)
    else:
        if ensmean: # @@ note, should add the casename to ensmean filename!
            fname = basepath + ens + '-ens' + '/' + field + '/' + field + '_' + comp + '_CanESM2_ensmean_' +\
                    suff + '.nc'
            if verb:
                print fname
            flist.append(fname)

        else:
            if subens:
                sims=subens
            else:
                sims = get_sims(ens)

            for sim in sims:
                for eii in range(1,ensnum+1):

                    fname=basepath + sim + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
                       sim + '_r' + str(eii) + 'i1p' + pnum + '_' + suff + '.nc'
                    if verb:
                        print fname
                    flist.append(fname)

            if orig in ('add','add45'):
                #try:
                print 'this is not going to be functional anymore (basepath wrong) 8/18/17'#@@@
                if ens=='historical':
                    if orig=='add45':
                        casename='historicalrcp45'
                        timeext='185001-201212'
                    else:
                        casename='historicalrcp85'
                        timeext='185001-202012'
                else:
                    print 'ens != historical. What to do?'
                        #raise
                #except:
                #    raise Exception

                orignum=5
                basedir='/HOME/rkm/work/DATA/CanESM2/' + casename

                for eii in np.arange(1,orignum+1):
                    fname=basedir + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
                           casename + '_r' + str(eii) + 'i1p' + pnum + '_' + timeext + '.nc'

                    if verb:
                        print fname
                    flist.append(fname)

    return flist
        
def get_lat(local=False):
    
    fnames=build_filenames({'field':'tas','comp':'Amon'},'historical',verb=False,local=local)
    return cnc.getNCvar(fnames[0],'lat')

def get_lon(local=False):
    
    fnames=build_filenames({'field':'tas','comp':'Amon'},'historical',verb=False,local=local)
    return cnc.getNCvar(fnames[0],'lon')


def get_sims(ens):

    if ens=='historical':
        return ('historical-r1','historical-r2','historical-r3','historical-r4','historical-r5')
    elif ens=='historicalNat':
        return ('historicalNat-r1','historicalNat-r2','historicalNat-r3','historicalNat-r4','historicalNat-r5')
    elif ens in ('historicalMisc', 'historicalAA','historicalSO'):
        return ('historicalMisc-r1','historicalMisc-r2','historicalMisc-r3','historicalMisc-r4','historicalMisc-r5')
    else:
        print 'ens not defined!! @@' # should throw an exception
        return -1


    
