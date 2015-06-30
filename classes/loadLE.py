import cccmautils as cutl
import pandas as pd
import cccmaNC as cnc
import constants as con
import numpy as np

con=reload(con)

bp=con.get_LEbasepath()

basepath=bp['basepath']
ensnum=10

#  Make this a class! CanESM2LE

def load_LEdata(fielddict, ens, seas=None, timesel=None,infodict=None,ftype='fullts',
                calctype=None, calcdict=None, rettype='dict',conv=1, region=None,local=False,
                orig=None,verb=True):
    """ def loadLEdata(fielddict, seas=('DJF','MAM','JJA','SON'), timesel=None, infodict=None, calctype=None, calcdict=None)

            ftype: type of filename to build. Right now just 'fullts' for full timeseries
              or '1950-2020_climo' or 'ensmean'
            seas has to be a tuple
            rettype: 'dict' or 'ndarray' as return type. 
                      default is dict of DataFrames (to make a Panel). 
                      else, 'ndarray' is a 4D array
            local: is data in ~/work/DATA (local=True) or /raid/ra40
                      
            @@ add regional avg functionality?

            returns an object of type rettype
    """

    field = fielddict['field']
    ncfield=fielddict['ncfield']
    comp = fielddict['comp']

    flist = build_filenames(fielddict, ens,ftype=ftype,timesel=timesel,local=local,orig=orig,verb=verb)

    fname1 = flist[0]
    if verb:
        print ' @@ fname1 ' + fname1
    timedim = cnc.getNCvar(fname1,'time',timesel=timesel)
    #print timedim

    if rettype=='dict':
        fldret={}
    elif rettype=='ndarray':
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
    """ orig='just' means just load the original 5 runs. 
        orig='add' means add the original 5 to the full LE.
    """

    return load_LEdata(fielddict, ens, seas=seas, timesel=timesel, infodict=infodict,
                       ftype=ftype, calctype=calctype, calcdict=calcdict, rettype=rettype,
                       conv=conv, region=region, local=local,orig=orig,verb=verb)
    


def build_filenames(fielddict, ens, ftype='fullts',timesel=None,verb=True,local=False,orig=None):
    """ here we know that each 'sim' has 10 sub-ensemble members

        ftype: type of filename to build. Right now just 'fullts' for full timeseries
              or 'fullclimo' for '1950-2020_climo' or, given timesel:
              for styr-enyr_climo. or 'ensmean'
        local: is the data on /raid/ra40 or ~/work/DATA (local=True)
              
        returns a list of all filenames (50). @@add original 5?!
        
    """
    if ens=='historicalMisc':
        pnum='4'
    else: pnum='1'

    ensmean=False
    if ftype=='fullts':
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

    if local:
        basepath = '/HOME/rkm/work/DATA/CanESM2/LE/'
    else:
        basepath=bp['basepath']

        
    field=fielddict['field']
    comp=fielddict['comp']
    
    flist=[]
    if orig=='just':
        if ens=='historical':
            casename='historicalrcp45'

        orignum=5
        basedir='/HOME/rkm/work/DATA/CanESM2/' + casename

        for eii in np.arange(1,orignum+1):
            fname=basedir + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
                   casename + '_r' + str(eii) + 'i1p' + pnum + '_185001-201212.nc'

            flist.append(fname)
    else:
        if ensmean: # @@ note, should add the casename to ensmean filename!
            fname = basepath + ens + '-ens' + '/' + field + '/' + field + '_' + comp + '_CanESM2_ensmean_' +\
                    suff + '.nc'
            if verb:
                print fname
            flist.append(fname)

        else:
            sims = get_sims(ens)

            for sim in sims:
                for eii in range(1,ensnum+1):

                    fname=basepath + sim + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
                       sim + '_r' + str(eii) + 'i1p' + pnum + '_' + suff + '.nc'
                    if verb:
                        print fname
                    flist.append(fname)
        if orig=='add':
            if ens=='historical':
                casename='historicalrcp45'
                
            orignum=5
            basedir='/HOME/rkm/work/DATA/CanESM2/' + casename

            for eii in np.arange(1,orignum+1):
                fname=basedir + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
                       casename + '_r' + str(eii) + 'i1p' + pnum + '_185001-201212.nc'
                
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
    elif ens=='historicalMisc':
        return ('historicalMisc-r1','historicalMisc-r2','historicalMisc-r3','historicalMisc-r4','historicalMisc-r5')
    else:
        print 'ens not defined!! @@' # should throw an exception
        return -1


    
