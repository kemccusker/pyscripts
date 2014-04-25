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

""" Note, other testing of this is done in compareOBS.py when plotseacycmag=1
    created 4/23/2014
"""
checkpert1s=0 # check new pert1 BC against original
checknewbcs=1 # check BCs for new runs: rcp85 2022-2032, etc

timeperp = '2022-2032'
#timeperp = '2042-2052'


plat = platform.system()
if plat == 'Darwin':  # means I'm on my mac
    basepath2 = '/Volumes/MyPassport1TB/DATA/CanSISE/'
    basepath = '/Users/kelly/CCCma/CanSISE/BoundaryConditionFiles/' # @@ not sure
    
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/BCs/'



if checkpert1s:

    # New Pert1 dataset that checks for open water colder than 271.2
    pert1newf = 'GTfrzchkhistoricalrcp452002-2012_BC_CanESM2_historical_1979-1989_1870010100-2010120100.nc'

    lat = cnc.getNCvar(pert1newf,'lat')
    lon = cnc.getNCvar(pert1newf,'lon')

    pert1n = cnc.getNCvar(pert1newf,'GT')
    pert1n = pert1n[0:12,...]




    # original BC files for PERT1. I expected that the ssts would be below freezing in some places
    #    where SICN was gone in PERT1 but not control. 
    pert1f = '/home/rkm/work/BCs/CanESM2/GT_BC_CanESM2_historical1979-1989_1870010100-2011020100.nc'
    pert1=cnc.getNCvar(pert1f,'GT')
    pert1 = pert1[0:12,...]

    pert1sicf = '/home/rkm/work/BCs/CanESM2/SICN_BC_CanESM2_historical2002-2012_1870010100-2011020100.nc'
    pert1sic = cnc.getNCvar(pert1sicf,'SICN')
    pert1sic = pert1sic[0:12,...]

    landmask = con.get_t63landmask()
    landmask = np.tile(landmask,(12,1,1))

    #pert1n[pert1n<271.2] = -1000 # @@ just for testing
    pert1n = ma.masked_where(landmask==-1,pert1n)
    pert1n = ma.masked_where(pert1sic>.15,pert1n)

    #pert1[pert1<271.2] = -1000 # @@ just for testing
    pert1 = ma.masked_where(landmask==-1,pert1)
    pert1 = ma.masked_where(pert1sic>.15,pert1)


    #cplt.map_allmonths(pert1n,lat,lon,cmin=271.2-35,cmax=271.2+35,cmap='blue2red_20',type='nh',conts=[271.2],climo=1,lmask=1)

    #cplt.map_allmonths(pert1,lat,lon,cmin=271.2-35,cmax=271.2+35,cmap='blue2red_20',type='nh',conts=[271.2],climo=1,lmask=1)

    cplt.map_allmonths(pert1n-pert1,lat,lon,cmin=-.1,cmax=.1,cmap='blue2red_20',type='nh',climo=1)
    #cplt.map_allmonths(landmask,lat,lon,cmap='blue2red_20',type='nh',climo=1)

if checknewbcs:


    deni = 913 # density of ice

    if timeperp=='2022-2032':
        cming = -5; cmaxg = 5 # GT
        cminc = -.4; cmaxc = .4 # SI concentration
        cmint = -.7*deni; cmaxt = .7*deni # SI thickness
    else:
        cming = -6; cmaxg = 6
        cminc = -.6; cmaxc = .6 # SI concentration
        cmint = -.8*deni; cmaxt = .8*deni # SI thickness
        
    ctlfile = basepath + 'CanESM2/GT_BC_CanESM2_historical_1979-1989_1870010100-2010120100.nc'
    pertfile = basepath + 'CanESM2/GTadjusted_BC_CanESM2_rcp85_' + timeperp + '_1870010100-2010120100_abs10thresh.nc'

    lat = cnc.getNCvar(ctlfile,'lat')
    lon = cnc.getNCvar(ctlfile,'lon')
    fldc = cnc.getNCvar(ctlfile,'GT')
    fldp = cnc.getNCvar(pertfile,'GT')

    #cplt.map_allmonths(fldp-fldc,lat,lon,cmin=cming,cmax=cmaxg,cmap='blue2red_20',type='nh')


    ctlfile = basepath + 'CanESM2/SICN_BC_CanESM2_historical_1979-1989_1870010100-2010120100.nc'
    pertfile = basepath + 'CanESM2/SICN_BC_CanESM2_rcp85_' + timeperp + '_1870010100-2010120100.nc'

    fldc = cnc.getNCvar(ctlfile,'SICN')
    fldp = cnc.getNCvar(pertfile,'SICN')

    #cplt.map_allmonths(fldp-fldc,lat,lon,cmin=cminc,cmax=cmaxc,cmap='red2blue_w20',type='nh')
    

    ctlfile = basepath + 'CanESM2/SIC_BC_CanESM2_historical_1979-1989_1870010100-2010120100.nc'
    pertfile = basepath + 'CanESM2/SIC_BC_CanESM2_rcp85_' + timeperp + '_1870010100-2010120100.nc'

    fldc = cnc.getNCvar(ctlfile,'SIC')
    fldp = cnc.getNCvar(pertfile,'SIC')

    cplt.map_allmonths(fldp-fldc,lat,lon,cmin=cmint,cmax=cmaxt,cmap='red2blue_w20',type='nh')
   
