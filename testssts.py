import numpy as np
import numpy.ma as ma
import scipy as sp # scientific python
import scipy.stats
import matplotlib.pyplot as plt
import platform as platform
import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import cccmaNC as cnc

cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)
cnc = reload(cnc)


# New Pert1 dataset that checks for open water colder than 271.2
pert1newf = 'GTfreezechkfor_historicalrcp452002-2012_BC_CanESM2_historical_1979-1989_1870010100-2010120100.nc'
lat = cnc.getNCvar(pert1newf,'lat')
lon = cnc.getNCvar(pert1newf,'lon')

pert1n = cnc.getNCvar(pert1newf,'GT')
pert1n = pert1n[0:12,...]




# original BC files for PERT1. I expected that the ssts would be below freezing in some places
#    where SICN was gone in PERT1 but not control. Masking out landmask==-1 and SICN>.15, I can't find any tho. ??
pert1f = '/home/rkm/work/BCs/CanESM2/GT_BC_CanESM2_historical1979-1989_1870010100-2011020100.nc'
pert1=cnc.getNCvar(pert1f,'GT')
pert1 = pert1[0:12,...]

pert1sicf = '/home/rkm/work/BCs/CanESM2/SICN_BC_CanESM2_historical2002-2012_1870010100-2011020100.nc'
pert1sic = cnc.getNCvar(pert1sicf,'SICN')
pert1sic = pert1sic[0:12,...]

landmask = con.get_t63landmask()
landmask = np.tile(landmask,(12,1,1))


pert1n = ma.masked_where(landmask==-1,pert1n)
pert1n = ma.masked_where(pert1sic>.15,pert1n)

pert1 = ma.masked_where(landmask==-1,pert1)
pert1 = ma.masked_where(pert1sic>.15,pert1)


cplt.map_allmonths(pert1n,lat,lon,cmin=271.2-35,cmax=271.2+35,cmap='blue2red_20',type='nh',conts=[271.2],climo=1)

cplt.map_allmonths(pert1,lat,lon,cmin=271.2-35,cmax=271.2+35,cmap='blue2red_20',type='nh',conts=[271.2],climo=1)

cplt.map_allmonths(pert1n-pert1,lat,lon,cmin=-.1,cmax=.1,cmap='blue2red_20',type='nh',climo=1)
#cplt.map_allmonths(landmask,lat,lon,cmap='blue2red_20',type='nh',climo=1)
