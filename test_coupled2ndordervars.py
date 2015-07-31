""" bits taken from notebook of same name
    July 30 2015
"""

import cccmautils as cutl
import constants as con
import pandas as pd


zonal=False 
calctype='total'

domonth=False # otherwise do season
sea='ANN'
mo=1

if domonth:
    seasonalizedt={'mo':mo}
else:
    seasonalizedt={'season':sea}

basedir='/HOME/rkm/work/DATA/CanESM2/'
casename='prei2xco2icenh' 
basepath=basedir+casename+'/ts/'
timeper = '2021-3041'
ncfield='MLT'

if zonal:
    timesel='2031-01-01,3041-12-31' # skip first decade
    fldsuff='ZM'
else:
    timesel='3002-01-01,3041-12-31'
    fldsuff=''


Cp=1004 # specific heat at const pressure (J/K/kg)
Lv=2.5e6 # specific heat of condensation (latent heat of vapo) at 0C(?@@) (J/kg)
erad = con.get_earthrad() # m
grav = con.get_g() # m/s2

# 1 cal = 4.186J
# convert to cal / day
# J/s -> cal/day
W2calperday = (60*60*24)/4.186

lat = con.get_t63lat()
lon = con.get_t63lon()
#lev = con.get_t63lev() # 37 levs
fname=basepath + casename + '_u_' + timeper + '_ts.nc' # 22 levs!!
lev = cnc.getNCvar(fname,'plev')

dp = np.diff(lev)
nlon = len(lon)


def calc_PHT(fld, onelat, nlon, dp, calc='total'):
    """
         calculates zonal integration and vertical integration/average

            calc='total' calculates zonal and vertical integral
            calc='average' calculates zonal integral and vertical average

    """
    latrad = np.deg2rad(onelat)
    totlam=2*np.pi 
    dlam=totlam/np.float(nlon) # div by number of lons

    tmpfld = np.sum(fldtm,axis=1)*erad*np.cos(latrad)*dlam # zonal integration

    fldzint= (tmpfld[0:-1] + tmpfld[1:]) / 2 # interp in b/w levels

    if calc=='total':
        fldpint = np.sum(fldzint*dp) # vertically integrate
    elif calc=='average':
        fldpint = np.average(fldzint, weights=dp)
    else:
        print 'calc type not supported' # @@@
        return -1

    fldpint = fldpint/grav

    return fldpint


fdict={'MPEF': basepath + casename + '_meridional_potential_energy_flux' + fldsuff +'_' + timeper + '_ts.nc',
       'MSHF': basepath + casename + '_meridional_sens_heat_flux' + fldsuff + '_' + timeper + '_ts.nc',
       'MLHF': basepath + casename + '_meridional_latent_heat_flux' + fldsuff + '_' + timeper + '_ts.nc',
       'MFZM': basepath + casename + '_meridional_flux_zonal_mom' + fldsuff + '_' + timeper + '_ts.nc'}
flddt={'MPEF': 'MLT',
       'MSHF': 'MLT',
       'MLHF': 'MLT',
       'MFZM': 'MLT'}
condt={'MPEF': 1,
       'MSHF': Cp,
       'MLHF': Lv,
       'MFZM': 0.5} # conversion factors


fldintdt={}
for fkey in fdict:
    
    ncfield=flddt[fkey]
    
    fname=fdict[fkey]
    conv=condt[fkey] 

    fld=cnc.getNCvar(fname,ncfield,timesel=timesel)*conv

    fldint=np.zeros(len(lat))
    latidx=0
    for ll in lat: # do all lats
        latrad = np.deg2rad(ll)

        fldsub = np.squeeze(fld[:,:,latidx,:-1]) # get just one lat and remove extra lon
        fldtm = np.mean(cutl.seasonalize_monthlyts(fldsub,**seasonalizedt),axis=0) # time mean
    
        fldint[latidx] = calc_PHT(fldtm, latidx, nlon, dp, calc=calctype)
        
        latidx += 1
    
    fldintdt[fkey] = fldint


flddf = pd.DataFrame(fldintdt)
fldtot = flddf.sum(axis=1)


flddf.plot(x=lat)
plt.plot(lat,fldtot,color='k')
plt.axhline(y=0,color='0.5')
plt.xlim((-80,80))


flddf.plot(x=lat)
plt.plot(lat,fldtot,color='k')
plt.axhline(y=0,color='0.5')
plt.ylim((-1e16,1e16))
plt.xlim((-15,75))


(flddf*W2calperday/1e19).plot(x=lat)
plt.plot(lat,fldtot*W2calperday/1e19,color='k')
plt.axhline(y=0,color='0.5')
plt.ylim((-25,25))
plt.xlim((-15,75))
plt.title(str(seasonalizedt) + ' 1e19 cal/day')
