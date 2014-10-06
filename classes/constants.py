"""
constants.py
Use this module to house access to constant data, variables
    
    e.g. get_t63landmask()

    2/11/2014

"""

import numpy as np
from netCDF4 import Dataset
import platform as platform
import cccmaNC as cnc


def get_basepath():
    """ get_basepath():
           returns path to data, depending on what machine we're running on
           output: bp['basepath'|'subdir'] (dictionary with basepath and subdir keys)
                     typically the pathtodata would be basepath+model+subdir
    """

    bp={}
    
    plat = platform.system()   
    if plat == 'Darwin':  # means I'm on my mac
        bp['basepath'] = '/Volumes/MyPassport2TB/DATA/CanSISE/'
        bp['subdir'] = '/'

    else:  # on linux workstation in Vic
        bp['basepath'] = '/home/rkm/work/DATA/'
        bp['subdir'] = '/ts/'

    return bp

def get_BCbasepath():
    """ get_basepath():
           returns path to data, depending on what machine we're running on
           output: bp['basepath'|'subdir'] (dictionary with basepath and subdir keys)
                     typically the pathtodata would be basepath+model+subdir
    """

    bp={}
    
    plat = platform.system()   
    if plat == 'Darwin':  # means I'm on my mac
        bp['basepath'] = '/Volumes/MyPassport2TB/DATA/CanSISE/CanAM4/BCs/'
        bp['subdir'] = '/'

    else:  # on linux workstation in Vic
        bp['basepath'] = '/home/rkm/work/BCs/'
        bp['subdir'] = '/'

    return bp
    
def get_t63landmask(repeat=None):
    """ Return ground cover.
        64x129
        -1=land, 0=open water, 1=sea ice (not present), 2=inland lake

        repeat: shape you want landmask as. assume last 2 dims are lat,lon
        
        """
    plat = platform.system()

    if plat == "Linux":
        basepath = '/home/rkm/work/DATA/CanAM4/constants/'
    else:
        basepath = '/Users/kelly/CCCma/CanSISE/DATA/constants/' #@@

    ncfile = Dataset(basepath + 't63_landmask.nc','r')
    #float GC(lat, lon) ;
    #           GC:long_name = "Ground cover (-1=land, 0=open water, +1=sea ice)" ;
    #           GC:min_avg_max = -1.f, -0.2944274f, 2.f ;
    # Also 2 is inland lake
    lmask = ncfile.variables['GC'][...]

    if repeat != None:
        nrep = repeat[0:-2] # leave off last 2 dims (lat, lon)
        nrep = nrep + (1,1)
        lmask = np.tile(lmask,nrep)
        
    return lmask
        
def get_t63lat():
    plat = platform.system()

    if plat == "Linux":
        basepath = '/home/rkm/work/DATA/CanAM4/constants/'
    else:
        basepath = '/Users/kelly/CCCma/CanSISE/DATA/constants/' #@@

    fname = basepath + 't63_landmask.nc'

    return cnc.getNCvar(fname,'lat')

def get_t63lon():

    plat = platform.system()

    if plat == "Linux":
        basepath = '/home/rkm/work/DATA/CanAM4/constants/'
    else:
        basepath = '/Users/kelly/CCCma/CanSISE/DATA/constants/' #@@

    fname = basepath + 't63_landmask.nc'
    
    return cnc.getNCvar(fname,'lon')

def get_t63lev(): # prob isn't tied to t63..

    plat = platform.system()

    if plat == "Linux":
        basepath = '/home/rkm/work/DATA/CanAM4/constants/'
    else:
        basepath = '/Users/kelly/CCCma/CanSISE/DATA/constants/' #@@

    fname = basepath + 'kem1rcp85a_v_001-061_climo.nc'
    
    return cnc.getNCvar(fname,'plev')


def get_monweights():

    monwgts = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    monwgts = monwgts / 365.
    return monwgts

def get_mon():
    months = 'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'
    return months

def get_earthrad():
    earthrad = 6.37122e6  # m
    return earthrad

def get_g():
    # this is the value used to convert PHI in model output to meters (divide PHI by g)
    g = 9.810616 
    return g

## def casenamec():
##     return 'kemctl1'

## def casenamep1():
##     return 'kem1pert1'

## def casenamep2():
##     return 'kem1pert2'

## def casenamep3():
##     return 'kem1pert3'

## def casenamech():
##     return 'kemhadctl'

## def casenameph():
##     return 'kemhadpert'

def get_simsdict():
    """ get_simsdict(): return a dictionary of all simulation dictionaries
              of metadata.
              
               # simsdict subdict of metadata:
               #   fullname, altname (possibly shorter/better for paper),
               #   timestr, timesel, stdctl, diffname (when referring to pert-ctl)

              
    """

    # simsdict subdict of metadata:
    #   fullname, altname (possibly shorter/better for paper),
    #   timestr, timesel, stdctl, diffname (when referring to pert-ctl)

    # these will be keys (shortnames)
    sims = ('ctl1', 'pert1', 'pert2', 'pert3', 'rcp85a',
            'ctl1r1', 'ctl1r2', 'ctl1r3', 'ctl1r4', 'ctl1r4', 'ctl1ens',
            'pert2r1','pert2r2', 'pert2r3','pert2r4','pert2r5','pert2ens',
            'pert2r4ct',
            'hadctl', 'hadpert',
            'nsidcctl', 'nsidcpert',
            'ctl1e2', 'ctl1e3', 'ctl1e4', 'ctl1e5', 'ctl1ense',
            'pert2e2', 'pert2e3','pert2e4','pert2e5', 'pert2ense',
            'ctl1espr',
            'pert2espr') # super ensemble mean

    simsdict = dict.fromkeys(sims,{})
    model='CanAM4' # @@@ all of these are CanAM4
    
    ## # kemctl1
    ## meta = {'fullname': 'kemctl1', 'altname': 'CANctl', 'timestr': '001-121',
    ##         'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
    ##         'ensname': 'histIC'}
    ## simsdict['ctl1'] = meta

    # kemctl1
    meta = {'fullname': 'kemctl1', 'altname': 'E1ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histIC','model':model} # aka CAN
    simsdict['ctl1'] = meta
    
    # kem1pert1b
    meta = {'fullname': 'kem1pert1b', 'altname': 'nosst', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1','diffname': 'CANnosst',
            'ensname': None,'model':model}
    simsdict['pert1'] = meta

    ## # kem1pert2
    ## meta = {'fullname': 'kem1pert2', 'altname': 'CANpt', 'timestr': '001-121',
    ##         'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1','diffname': 'CAN',
    ##         'ensname': 'histIC'}
    ## simsdict['pert2'] = meta

    # kem1pert2
    meta = {'fullname': 'kem1pert2', 'altname': 'E1pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1','diffname': 'E1',
            'ensname': 'histIC','model':model} # aka CAN
    simsdict['pert2'] = meta
    
    # kem1pert3
    meta = {'fullname': 'kem1pert3', 'altname': 'nothk', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1','diffname': 'CANnothk',
            'ensname': None,'model':model}
    simsdict['pert3'] = meta

    # kem1rcp85a
    meta = {'fullname': 'kem1rcp85a', 'altname': 'RCPpt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1','diffname': 'RCPa',
            'ensname': None,'model':model}
    simsdict['rcp85a'] = meta

    # kemctl1r1
    meta = {'fullname': 'kemctl1r1', 'altname': 'R1ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histBC','model':model}
    simsdict['ctl1r1'] = meta

    # kemctl1r2
    meta = {'fullname': 'kemctl1r2', 'altname': 'R2ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histBC','model':model}
    simsdict['ctl1r2'] = meta

    # kemctl1r3
    meta = {'fullname': 'kemctl1r3', 'altname': 'R3ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histBC','model':model}
    simsdict['ctl1r3'] = meta

    # kemctl1r4
    meta = {'fullname': 'kemctl1r4', 'altname': 'R4ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histBC','model':model}
    simsdict['ctl1r4'] = meta

    # kemctl1r5
    meta = {'fullname': 'kemctl1r5', 'altname': 'R5ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histBC','model':model}
    simsdict['ctl1r5'] = meta

    # kemctl1ens
    meta = {'fullname': 'kemctl1ens', 'altname': 'ENSctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None, 'diffname': None,
            'ensname': 'histBCmean','model':model} # @@
    simsdict['ctl1ens'] = meta

    # kem1pert2r1
    meta = {'fullname': 'kem1pert2r1', 'altname': 'R1pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1r1', 'diffname': 'R1',
            'ensname': 'histBC','model':model}
    simsdict['pert2r1'] = meta

    # kem1pert2r2
    meta = {'fullname': 'kem1pert2r2', 'altname': 'R2pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1r2', 'diffname': 'R2',
            'ensname': 'histBC','model':model}
    simsdict['pert2r2'] = meta

    # kem1pert2r3
    meta = {'fullname': 'kem1pert2r3', 'altname': 'R3pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1r3', 'diffname': 'R3',
            'ensname': 'histBC','model':model}
    simsdict['pert2r3'] = meta

    # kem1pert2r4
    meta = {'fullname': 'kem1pert2r4', 'altname': 'R4pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1r4', 'diffname': 'R4',
            'ensname': 'histBC','model':model}
    simsdict['pert2r4'] = meta

    # kem1pert2r5
    meta = {'fullname': 'kem1pert2r5', 'altname': 'R5pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1r5', 'diffname': 'R5',
            'ensname': 'histBC','model':model}
    simsdict['pert2r5'] = meta

    # kem1pert2ens
    meta = {'fullname': 'kem1pert2ens', 'altname': 'ENSpt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1ens', 'diffname': 'ENS',
            'ensname': 'histBCmean','model':model}# @@
    simsdict['pert2ens'] = meta

    # kem1pert2r4ct: constant thickness 'ct'
    meta = {'fullname': 'kem1pert2r4ct', 'altname': 'R4ctpt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1r4', 'diffname': 'R4ct',
            'ensname': None,'model':model}
    simsdict['pert2r4ct'] = meta

    # kemhadctl
    meta = {'fullname': 'kemhadctl', 'altname': 'HADctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': None,'model':model}
    simsdict['hadctl'] = meta

    # kemhadpert
    meta = {'fullname': 'kemhadpert', 'altname': 'HADpt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'hadctl','diffname': 'HAD',
            'ensname': None,'model':model}
    simsdict['hadpert'] = meta

    # kemnsidcctl
    meta = {'fullname': 'kemnsidcctl', 'altname': 'NSIDCctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': None,'model':model}
    simsdict['nsidcctl'] = meta

    # kemnsidcpert
    meta = {'fullname': 'kemnsidcpert', 'altname': 'NSIDCpt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'nsidcctl','diffname': 'NSIDC',
            'ensname': None,'model':model}
    simsdict['nsidcpert'] = meta

    # kemctl1e2
    meta = {'fullname': 'kemctl1e2', 'altname': 'E2ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histIC','model':model}
    simsdict['ctl1e2'] = meta

    # kemctl1e3
    meta = {'fullname': 'kemctl1e3', 'altname': 'E3ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histIC','model':model}
    simsdict['ctl1e3'] = meta

    # kemctl1e4
    meta = {'fullname': 'kemctl1e4', 'altname': 'E4ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histIC','model':model}
    simsdict['ctl1e4'] = meta

    # kemctl1e5
    meta = {'fullname': 'kemctl1e5', 'altname': 'E5ctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None,'diffname': None,
            'ensname': 'histIC','model':model}
    simsdict['ctl1e5'] = meta

    # kemctl1ense
    meta = {'fullname': 'kemctl1ense', 'altname': 'ENSEctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None, 'diffname': None,
            'ensname': 'histICmean','model':model} # @@
    simsdict['ctl1ense'] = meta

    # kemctl1ensspr
    meta = {'fullname': 'kemctl1ensspr', 'altname': 'ESPRctl', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': None, 'diffname': None,
            'ensname': 'histBCICmean','model':model} # @@
    simsdict['ctl1espr'] = meta

    # kem1pert2e2
    meta = {'fullname': 'kem1pert2e2', 'altname': 'E2pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1e2', 'diffname': 'E2',
            'ensname': 'histIC','model':model}
    simsdict['pert2e2'] = meta

    # kem1pert2e3
    meta = {'fullname': 'kem1pert2e3', 'altname': 'E3pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1e3', 'diffname': 'E3',
            'ensname': 'histIC','model':model}
    simsdict['pert2e3'] = meta

     # kem1pert2e4
    meta = {'fullname': 'kem1pert2e4', 'altname': 'E4pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1e4', 'diffname': 'E4',
            'ensname': 'histIC','model':model}
    simsdict['pert2e4'] = meta

     # kem1pert2e5
    meta = {'fullname': 'kem1pert2e5', 'altname': 'E5pt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1e5', 'diffname': 'E5',
            'ensname': 'histIC','model':model}
    simsdict['pert2e5'] = meta

    # kem1pert2ense
    meta = {'fullname': 'kem1pert2ense', 'altname': 'ENSEpt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1ense', 'diffname': 'ENSE',
            'ensname': 'histICmean','model':model}# @@
    simsdict['pert2ense'] = meta

    # kem1pert2ensspr
    meta = {'fullname': 'kem1pert2ensspr', 'altname': 'ESPRpt', 'timestr': '001-121',
            'timesel': '0002-01-01,0121-12-31','stdctl': 'ctl1espr', 'diffname': 'ESPR',
            'ensname': 'histBCICmean','model':model}# @@
    simsdict['pert2espr'] = meta

    return simsdict


def get_simmeta(simkey):
    """ get_simmeta(simkey): return the dictionary of metadata for given key
    """

    simsdict = get_simsdict()
    return simsdict[simkey]
    
def get_simpairsdict():
    """ get_simpairsdict(): returns a dictionary of dictionaries of sim metadata
              Layout is: pairdict --> ctldict and pertdict --> each of those
                            is a metadict (w/ sim metadata)
    """
    simsdict=get_simsdict()
    
    # build the pair dictionary. This matches a perturbation
    #  simulation with its standard control, and leaves out anything
    #  that isn't yet defined in the simsdict.
    pairdict = {}
    for skey in simsdict.iterkeys():

        #print skey
        simdict=simsdict[skey]
        

        if simdict == {}:
            pass
        elif simdict['diffname'] != None:
            #print simdict['diffname']
            # add the diffname to pairs dict as new key
            # set up the simpair
            simpt = simsdict[skey]
            diffkey = simpt['diffname']
            stdctl = simpt['stdctl']
            
            ctldict = simsdict[stdctl]
            ctldict['shortname'] = stdctl
            pertdict = simsdict[skey]
            pertdict['shortname'] = skey
            thepair = {'ctl': ctldict, 'pert': pertdict}
            
            pairdict[diffkey] = thepair

    return pairdict

def get_simpair(diffkey):
    """ get_simpair(diffkey): returns a dictionary of
                 ctl and pert simulation for given
                 diffkey

    """
    
    pairdict = get_simpairsdict()
    return pairdict[diffkey]

def build_filepath(sim,field,timeper=None):
    """ build_filename(sim,field,timeper=None)
            sim: simulation shortname (ie. 'ctl1')
            simtype: 'ctl' or 'pert'
            field: what variable
            timeper: filename time period. Default is to use
                     the time period specified in the sim metadata
                     
            returns filename
            
    """
    meta = get_simmeta(sim)
    simfull = meta['fullname']
    timstr = meta['timestr']
    model = meta['model']
    
    bp=get_basepath()
    basepath=bp['basepath'] + model + '/'; subdir=bp['subdir']   

    froot = basepath + simfull + subdir + simfull + '_' 

    fname = froot + field + '_' + timstr + '_ts.nc'

    return fname

def build_filepathpair(simpairkey,field,timeper=None):
    """ build_filenamepair(simpairkey,field,timeper=None)
            sim: simulation diffname/pairname (ie. 'R1' or 'ENSE')
            field: what variable
            timeper: filename time period. Default is to use
                     the time period specified in the sim metadata
                     
            returns tuple (fnamecontrol, fnamepert)
            
    """
    
    simpair = get_simpair(simpairkey)
    simctl = simpair['ctl']['shortname']
    simpert = simpair['pert']['shortname']
    
    fnamec = build_filepath(simctl,field,timeper)
    fnamep = build_filepath(simpert,field,timeper)

    return (fnamec,fnamep)


def getBCfilenames(field,sim=None):

    bp = get_BCbasepath()
    basepath=bp['basepath']
    
    fnames = {'kemctl1': basepath + 'CanESM2/canesm2_kemctl1_128x64_0001_0125_' + field + '.nc',
              'kem1pert1b': basepath + 'CanESM2/canesm2_kem1pert1b_128x64_0001_0125_' + field + '.nc',
              'kem1pert2': basepath + 'CanESM2/canesm2_kem1pert2_128x64_0001_0125_' + field + '.nc',
              'kem1pert3': basepath + 'CanESM2/canesm2_kem1pert3_128x64_0001_0125_' + field + '.nc',
              'kem1rcp85a': basepath + 'CanESM2/canesm2_kem1rcp85a_128x64_0001_0125_' + field + '.nc',
              'kemctl1r1':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kemctl1r1_128x64_0001_0125_' + field + '.nc',
              'kemctl1r2':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kemctl1r2_128x64_0001_0125_' + field + '.nc',
              'kemctl1r3':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kemctl1r3_128x64_0001_0125_' + field + '.nc',
              'kemctl1r4':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kemctl1r4_128x64_0001_0125_' + field + '.nc',
              'kemctl1r5':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kemctl1r5_128x64_0001_0125_' + field + '.nc',
              'kem1pert2r1':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kem1pert2r1_128x64_0001_0125_' + field + '.nc',
              'kem1pert2r2':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kem1pert2r2_128x64_0001_0125_' + field + '.nc',
              'kem1pert2r3':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kem1pert2r3_128x64_0001_0125_' + field + '.nc',
              'kem1pert2r4':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kem1pert2r4_128x64_0001_0125_' + field + '.nc',
              'kem1pert2r5':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kem1pert2r5_128x64_0001_0125_' + field + '.nc',
              'kem1pert2r4ct':  basepath + 'CanESM2/ENSMEMBERS/canesm2_kem1pert2r4ct_128x64_0001_0125_' + field + '.nc',
              'kemhadctl': basepath + 'HadISST/hadisst_kemhadctl_128x64_0001_0125_' + field + '.nc',
              'kemhadpert': basepath + 'HadISST/hadisst_kemhadpert_128x64_0001_0125_' + field + '.nc',
              'kemnsidcctl': basepath + 'NSIDC/nsidcbt_kemnsidcctl_128x64_0001_0125_' + field + '.nc',
              'kemnsidcpert': basepath + 'NSIDC/nsidcbt_kemnsidcpert_128x64_0001_0125_' + field + '.nc' }

    if sim==None:
        return fnames
    else:
        return fnames[sim]
    

def get_regiondict():
    """ get_regiondict(): return a dictionary of all defined regions

                 regiondict:
                          polcap70: latlims=[70,89]; lonlims=[0,359]
                          polcap65: latlims=[65,89]; lonlims=[0,359]
                          polcap60: latlims=[60,89]; lonlims=[0,359]
                          eurasia: latlims=[35,60]; lonlims=[40,120]
                          ntham: latlims=[35,60]; lonlims=[240,280]
                          nthatl: latlims=[35,60]; lonlims=[300,360]
    """
    regions = ('polcap70', 'polcap65', 'polcap60', 'eurasia',
               'ntham', 'nthatl')
    regdict = dict.fromkeys(regions, {})

    regdict['polcap70'] = {'latlims': [70,89], 'lonlims': [0,359]}
    regdict['polcap65'] = {'latlims': [65,89], 'lonlims': [0,359]}
    regdict['polcap60'] = {'latlims': [60,89], 'lonlims': [0,359]}
    regdict['eurasia'] = {'latlims': [35,60], 'lonlims': [40,120]}
    regdict['ntham'] = {'latlims': [35,60], 'lonlims': [240,280]}
    regdict['nthatl'] = {'latlims': [35,60], 'lonlims': [300,360]}
    
    return regdict
    
def get_regionlims(regname):
    """ get_regionlim(regname): Given a region name, return a dict of latlims and lonlims

                      regname options:
                          polcap70: latlims=[70,89]; lonlims=[0,359] # Polar cap north of 70N
                          polcap65: latlims=[65,89]; lonlims=[0,359] # Polar cap north of 65N for NAM proxy
                          polcap60: latlims=[60,89]; lonlims=[0,359] # Polar cap north of 60N to match pattern corrs
                          eurasia: latlims=[35,60]; lonlims=[40,120] # Eurasia 35-60N, 40E-120E
                          ntham: latlims=[35,60]; lonlims=[240,280]  # North America 35-60N, 120W-80W
                          nthatl: latlims=[35,60]; lonlims=[300,360] # North Atlantic 35-60N, 60W-0

    """


    regdict = get_regiondict()

    return regdict[regname]

    
def get_t63regionmask(regname,limsdict=None):
    """ get_t63regionmask(regname,limsdict=None):
                          Return a lat x lon mask for given region name.
                          (Masks everything BUT the region of interest)
                          If regname='other', use region defined in limsdict
    """

    lat=get_t63lat()
    lon=get_t63lon()
    lons,lats = np.meshgrid(lon,lat)    
        
    # create mask
    if regname=='other':
        pass # use given limsdict
    else:
        limsdict = get_regionlims(regname)
        
    latlims = limsdict['latlims']
    lonlims = limsdict['lonlims']

    reglatsbool = np.logical_and(lat>latlims[0],lat<latlims[1])
    reglonsbool = np.logical_and(lon>lonlims[0],lon<lonlims[1])

    # create mask of everything but the region of interest
    regmask = np.logical_or( 
                            np.logical_or(lats<latlims[0],lats>latlims[1]), 
                            np.logical_or(lons<lonlims[0],lons>lonlims[1]))

    return regmask


if __name__ == "__main__":
    plot_allregions()
    

