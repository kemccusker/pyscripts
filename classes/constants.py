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

def casenamec():
    return 'kemctl1'

def casenamep1():
    return 'kem1pert1'

def casenamep2():
    return 'kem1pert2'

def casenamep3():
    return 'kem1pert3'

def casenamech():
    return 'kemhadctl'

def casenameph():
    return 'kemhadpert'


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

                 @@ should define region dict upon import, not in a func...
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
    
