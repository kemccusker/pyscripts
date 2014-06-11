"""
constants.py
Use this module to house access to constant data, variables
    
    e.g. get_t63landmask()

    2/11/2014

"""

import numpy as np
from netCDF4 import Dataset
import platform as platform

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


