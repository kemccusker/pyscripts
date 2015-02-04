import cccmautils as cutl
import pandas as pd
import cccmaNC as cnc

basepath='/ra40/data/kem/CanSISE/CanESM2/LE/'
ensnum=10

#  Make this a class! CanESM2LE

def load_LEdata(fielddict, ens,seas=None, timesel=None,infodict=None,ftype='fullts',calctype=None, calcdict=None):
    """ def loadLEdata(fielddict, seas=('DJF','MAM','JJA','SON'), timesel=None, infodict=None, calctype=None, calcdict=None)

            ftype: type of filename to build. Right now just 'fullts' for full timeseries
              or '1950-2020_climo'
            seas has to be a tuple

            returns a dict of DataFrames (to make a Panel)
    """

    field = fielddict['field']
    ncfield=fielddict['ncfield']
    comp = fielddict['comp']

    flist = build_filenames(fielddict, ens,ftype=ftype)
    flddt = {}

    fname1 = flist[0]
    timedim = cnc.getNCvar(fname1,'time',timesel=timesel)
    print timedim

    numens = len(flist)
    #fldall = np.zeros(numens)
    #dlist=[]
    
    for fname in flist:
        
        if seas==None:
            # don't average, return all months
            fld = cnc.getNCvar(fname,ncfield,timesel=timesel)
            n0,n1,n2=fld.shape
            
            # do calcs:
            if calctype!=None:
                print 'calcs not implemented yet!! @@@'

            fldseas=fld.reshape((n0,n1*n2))
            
        else:
            fldseas={}
            for sea in seas:
                fld = cnc.getNCvar(fname,ncfield,timesel=timesel,seas=sea)

                # do calcs
                if calctype!=None:
                    print 'calcs not implemented yet!! @@@'
                
                fldseas[sea]=fld
            
        flddt[fname]= pd.DataFrame(fldseas,index=timedim)
        #dlist.append(fldseas)

        
    return flddt # dict of dataframes


def build_filenames(fielddict, ens, ftype='fullts'):
    """ here we know that each 'sim' has 10 sub-ensemble members

        ftype: type of filename to build. Right now just 'fullts' for full timeseries
              or '1950-2020_climo'
              
        returns a list of all filenames (50). @@add original 5?!
        
    """

    if ftype=='fullts':
        suff='195001-202012'
    else:
        suff=ftype
        
    field=fielddict['field']
    comp=fielddict['comp']
    
    sims = get_sims(ens)

    flist=[]
    for sim in sims:
        for eii in range(1,ensnum+1):

            fname=basepath + sim + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
               sim + '_r' + str(eii) + 'i1p1_' + suff + '.nc'
            print fname
            flist.append(fname)

    return flist
        
    


def get_sims(ens):

    if ens=='historical':
        return ('historical-r1','historical-r2','historical-r3','historical-r4','historical-r5')
    elif ens=='historicalNat':
        return ('historicalNat-r1','historicalNat-r2','historicalNat-r3','historicalNat-r4','historicalNat-r5')
    else:
        print 'ens not defined!! @@' # should throw an exception
        return -1

    
