"""
cmip5meta.py
    Use this module to house metadata about CMIP5 models, simulations, years, etc
    
    12/4/2014

"""

import numpy as np
from netCDF4 import Dataset
import platform as platform
import cccmaNC as cnc


cmip5base='/rd40/data/CMIP5/' # @@@ check for what machine...

models = ('ACCESS1-0','ACCESS1-3','BNU-ESM','CCSM4','CESM1-BGC','CESM1-CAM5','CMCC-CMS',
          'CMCC-CM','CNRM-CM5','CSIRO-Mk3-6-0','CanESM2','EC-EARTH','FGOALS-g2','GFDL-CM3',
          'GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR',
          'IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC-ESM-CHEM','MIROC-ESM','MIROC5','MPI-ESM-LR',
          'MPI-ESM-MR','MRI-CGCM3','NorESM1-ME')

allmodeldt = dict.fromkeys(models)

# standard defaults (for historical sim)
histyears='185001-200512'
histstyr=1850
stmon=1
enmon=12
histlen=156

model='ACCESS1-0'
meta = {'model': model,
        'fyears': histyears, # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'numens': 1 # number of ensemble members
        }
allmodeldt[model] = meta


model='ACCESS1-3'
meta = {'model': model,
        'fyears': histyears, # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'numens': 1, # number of ensemble members
        }
allmodeldt[model] = meta
    

model='BNU-ESM'
meta = {'model': model,
        'fyears': histyears, # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'numens': 1 # number of ensemble members
        }
allmodeldt[model] = meta

model='CCSM4'
meta = {'model': model,
        'fyears': histyears, # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'numens': 6 # number of ensemble members
        }
allmodeldt[model] = meta

model='CESM1-BGC'
meta = {'model': model,
        'fyears': histyears, # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'numens': 1 # number of ensemble members
        }
allmodeldt[model] = meta

model='CESM1-CAM5'
meta = {'model': model,
        'fyears': histyears, # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'styr': histstyr,
        'stmon': stmon,
        'enmon': enmon,
        'flen': histlen, # length of each file in years
        'fileyrincr': histlen, # length of each file
        'spclfirstfincr': None, # special incremement of first file
        'spcllastfincr':None, # special increment of last file 
        'numens': 3 # number of ensemble members
        }
allmodeldt[model] = meta


model='CMCC-CM'
meta = {'model': model,
        'fyears': ('185001-185912','186001-186912','187001-187912','188001-188912',
                   '189001-189912','190001-190912','191001-191912','192001-192912',
                   '193001-193912','194001-194912','195001-195912','196001-196912',
                   '197001-197912','198001-198912','199001-199912','200001-200512',), # filename years
        'yrspan': (histstyr,2005), # span of years in all files
        'styr': histstyr,
        'stmon': stmon,
        'enmon': enmon,
        'flen': histlen, # length of all files in years
        'fileyrincr': 10, # length of each file in years
        'spclfirstfincr': None, # special incremement of first file
        'spcllastfincr':6, # special increment of last file 
        'numens': 1 # number of ensemble members
        }
allmodeldt[model] = meta

model='CMCC-CMS'
meta = {'model': model,
        'fyears': ('185001-185912','186001-186912','187001-187912','188001-188912',
                   '189001-189912','190001-190912','191001-191912','192001-192912',
                   '193001-193912','194001-194912','195001-195912','196001-196912',
                   '197001-197912','198001-198912','199001-199912','200001-200512',), # filename years
        'yrspan': (histstyr,2005), # span of years in all files
        'styr': histstyr,
        'stmon': stmon,
        'enmon': enmon,
        'flen': histlen, # length of all files in years
        'fileyrincr': 10, # length of each file in years
        'spclfirstfincr': None, # special incremement of first file
        'spcllastfincr':6, # special increment of last file 
        'numens': 1 # number of ensemble members
        }
allmodeldt[model] = meta

model='CMCC-CESM'
meta = {'model': model,
        'fyears': ('185001-185512','185601-186112','186201-186712','186801-187312',
                   '187401-187912','188001-188512','188601-189112','189201-189712',
                   '189801-190312','190401-190912','191001-191512','191601-192112',
                   '192201-192712','192801-193312','193401-193912','194001-194512',
                   '194601-195112','195201-195712','195801-196312','196401-196912',
                   '197001-197512','197601-198112','198201-198712','198801-199312',
                   '199401-199912','200001-200512'), # filename years
        'yrspan': (histstyr,2005), # span of years in all files
        'styr': histstyr,
        'stmon': stmon,
        'enmon': enmon,
        'flen': histlen, # length of all files in years
        'fileyrincr': 6, # length of each file in years
        'spclfirstfincr': None, # special incremement of first file
        'spcllastfincr':None, # special increment of last file 
        'numens': 1 # number of ensemble members
        }
allmodeldt[model] = meta

model='CNRM-CM5'
meta = {'model': model,
        'fyears': histyears, # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'styr': histstyr,
        'stmon': stmon,
        'enmon': enmon,
        'flen': histlen, # length of each file in years
        'fileyrincr': histlen, # length of each file
        'spclfirstfincr': None, # special incremement of first file
        'spcllastfincr':None, # special increment of last file 
        'numens': 1 # number of ensemble members
        }
allmodeldt[model] = meta


model='CSIRO-Mk3-6-0'
meta = {'model': model,
        'fyears': histyears, # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'styr': histstyr,
        'stmon': stmon,
        'enmon': enmon,
        'flen': histlen, # length of each file in years
        'fileyrincr': histlen, # length of each file
        'spclfirstfincr': None, # special incremement of first file
        'spcllastfincr':None, # special increment of last file 
        'numens': 10 # number of ensemble members
        }
allmodeldt[model] = meta

model='CanESM2'
meta = {'model': model,
        'fyears': histyears, # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'styr': histstyr,
        'stmon': stmon,
        'enmon': enmon,
        'flen': histlen, # length of each file in years
        'fileyrincr': histlen, # length of each file
        'spclfirstfincr': None, # special incremement of first file
        'spcllastfincr':None, # special increment of last file 
        'numens': 5 # number of ensemble members
        }
allmodeldt[model] = meta

model='EC-EARTH'
meta = {'model': model,
        'fyears':histyears, #('185001-185912','186001-186912','187001-187912','188001-188912',
                   #'189001-189912','190001-190912','191001-191912','192001-192912',
                   #'193001-193912','194001-194912','195001-195912','196001-196912',
                   #'197001-197912','198001-198912','199001-199912','200001-200512'), # filename years
        'yrspan': (histstyr,2005), # span of years in the file
        'styr': histstyr,
        'stmon': stmon,
        'enmon': enmon,
        'flen': histlen, # length of each file in years
        'fileyrincr': histlen, # length of each file
        'spclfirstfincr': None, # special incremement of first file
        'spcllastfincr':None, # special increment of last file 
        'numens': 5, # number of ensemble members
        'skipens': [1,4,7,8] # skip these ens members. weird years or whatever.@@
        }
allmodeldt[model] = meta


model='FGOALS-g2'
meta = {'model': model,
        'fyears': ('185001-185912','186001-186912','187001-187912','188001-188912',
                   '189001-189912','190001-190912','191001-191912','192001-192912',
                   '193001-193912','194001-194912','195001-195912','196001-196912',
                   '197001-197912','198001-198912','199001-199912','200001-200512'), # filename years
        'yrspan': (histstyr,2005), # span of years in all files
        'styr': histstyr,
        'stmon': stmon,
        'enmon': enmon,
        'flen': histlen, # length of all files in years
        'fileyrincr': 10, # length of each file in years
        'spclfirstfincr': None, # special incremement of first file
        'spcllastfincr':6, # special increment of last file 
        'numens': 5 # number of ensemble members
        }
allmodeldt[model] = meta



def get_modellist():
    
    return models

def get_cmip5basepath():
    return cmip5base
