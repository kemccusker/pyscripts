""" bits taken from notebook of same name
    July 30 2015
"""

import cccmautils as cutl
import constants as con
import pandas as pd


printtofile=True

#zonal=False 
#calctype='vertint' # 'vertint' or 'vertavg', for calc'ing PHT in script
calctype='vertavg'

# this next setting is for calcs already done in nc file
#fldpre = '_vert_int'
fldpre = '_vert_avg'; 


domonth=False # otherwise do season
sea='ANN'
mo=1

if domonth:
    seasonalizedt={'mo':mo}
else:
    seasonalizedt={'season':sea}


if calctype=='vertint':
    ylimsW=(-10e15,10e15)
    ylimscal=(-25,25)
    ylimszaW2=(-2e14,1.5e14)
    ylimszaW=(-2e15,1.5e15)
elif calctype=='vertavg':
    ylimsW=(-1e11,1e11)
    ylimscal=(-.0004,.0004)
    ylimszaW2=(-2e10,2e10)
    ylimszaW=(-2e10,2e10)

basedir='/HOME/rkm/work/DATA/CanESM2/'
fldsuff=''

"""
casename='prei2xco2icenh' 
basepath=basedir+casename+'/ts/'
timeper = '2021-3041'
ncfield='MLT'
if zonal:
    # can't really do zonal mean if want zonal integration in PHT
    timesel='2031-01-01,3041-12-31' # skip first decade
    fldsuff='ZM'
else:
    timesel='3002-01-01,3041-12-31'
    fldsuff=''
"""

casename='preitest' 
basepath=basedir+casename+'/ts/'
timeper = '2922-2931'
ncfield='MLT'
timesel='2922-01-01,2931-12-31'



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
nlon = len(lon)-1 # removing extra lon


def calc_PHT(fld, onelat, nlon, dp, calc='vertint'):
    """
         calculates zonal integration and vertical integration/average

            calc='vertint' calculates zonal and vertical integral
            calc='vertavg' calculates zonal integral and vertical average

    """
    #print 'calc_PHT()'

    latrad = np.deg2rad(onelat)
    totlam=2*np.pi 
    dlam=totlam/np.float(nlon) # div by number of lons

    tmpfld = np.sum(fld,axis=1)*erad*np.cos(latrad)*dlam # zonal integration

    fldzint= (tmpfld[0:-1] + tmpfld[1:]) / 2 # interp in b/w levels

    if calc=='vertint':
        fldpint = np.sum(fldzint*dp) # vertically integrate
    elif calc=='vertavg':
        fldpint = np.average(fldzint, weights=dp)
    else:
        print 'calc type not supported' # @@@
        return -1

    fldpint = fldpint/grav

    return fldpint


def calc_PHT2(fld, onelat, nlon, dp=None, calc='vertint'):
    """
         calculates zonal integration and vertical integration/average

            calc='vertint' calculates zonal and vertical integral
            calc='vertavg' calculates zonal integral and vertical average
            calc='none' just zonal integral, b/c already vertically integrated(?)
    """
    #print 'calc_PHT2()'

    latrad = np.deg2rad(onelat)
    totlam=2*np.pi 
    dlam=totlam/np.float(nlon) # div by number of lons

    

    if calc != 'none':
        tmpfld = np.sum(fld,axis=1)*erad*np.cos(latrad)*dlam # zonal integration
        fldzint= (tmpfld[0:-1] + tmpfld[1:]) / 2 # interp in b/w levels
    else:
        tmpfld = np.sum(fld,axis=0)*erad*np.cos(latrad)*dlam # zonal integration
        fldzint = tmpfld

    if calc=='vertint':
        fldpint = np.sum(fldzint*dp) # vertically integrate
    elif calc=='vertavg':
        fldpint = np.average(fldzint, weights=dp)
    elif calc=='none':
        fldpint = fldzint
    else:
        print 'calc type not supported' # @@@
        return -1

    fldpint = fldpint/grav

    return fldpint


def calc_PHT3(fld, onelat, nlon, dp=None, calc='vertint'):
    """
         calculates zonal AVERAGE and vertical integration/average

            calc='vertint' calculates zonal and vertical integral
            calc='vertavg' calculates zonal integral and vertical average
            calc='none' just zonal integral, b/c already vertically integrated(?)
    """
    #print 'calc_PHT2()'

    latrad = np.deg2rad(onelat)
    totlam=2*np.pi 
    dlam=totlam/np.float(nlon) # div by number of lons

    

    if calc != 'none':
        tmpfld = np.sum(fld,axis=1)*erad*np.cos(latrad)*dlam # zonal integration
        tmpfld = tmpfld / totlam # zonal AVERAGE
        fldzint= (tmpfld[0:-1] + tmpfld[1:]) / 2 # interp in b/w levels
    else:
        tmpfld = np.sum(fld,axis=0)*erad*np.cos(latrad)*dlam # zonal integration
        tmpfld = tmpfld / totlam # zonal AVERAGE
        fldzint = tmpfld

    if calc=='vertint':
        fldpint = np.sum(fldzint*dp) # vertically integrate
    elif calc=='vertavg':
        fldpint = np.average(fldzint, weights=dp)
    elif calc=='none':
        fldpint = fldzint
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


fldintdt={}; fldintzadt={}
for fkey in fdict.keys():

    ncfield=flddt[fkey]
    
    fname=fdict[fkey]
    conv=condt[fkey] 

    fld=cnc.getNCvar(fname,ncfield,timesel=timesel)*conv

    fldint=np.zeros(len(lat))
    fldintza=np.zeros(len(lat))
    latidx=0
    for ll in lat: # do all lats
        latrad = np.deg2rad(ll)

        fldsub = np.squeeze(fld[:,:,latidx,:-1]) # get just one lat and remove extra lon
        fldtm = np.mean(cutl.seasonalize_monthlyts(fldsub,**seasonalizedt),axis=0) # time mean
    
        #fldint[latidx] = calc_PHT(fldtm, latidx, nlon, dp, calc=calctype)
        fldint[latidx] = calc_PHT2(fldtm, latidx, nlon, dp, calc=calctype) # zonal integration
        fldintza[latidx] = calc_PHT3(fldtm, latidx, nlon, dp, calc=calctype) # zonal average

        latidx += 1
    
    fldintdt[fkey] = fldint
    fldintzadt[fkey] = fldintza

flddf = pd.DataFrame(fldintdt)
fldtot = flddf.sum(axis=1)
fldzadf = pd.DataFrame(fldintzadt)
fldzatot = fldzadf.sum(axis=1)


flddf.plot(x=lat)
plt.plot(lat,fldtot,color='k')
plt.axhline(y=0,color='0.5')
plt.xlim((-80,80))
plt.xlabel('latitude')
plt.ylabel('poleward transport (W)')
plt.title(str(seasonalizedt) + ' ' + calctype + ' done in script')
if printtofile:
    plt.savefig('PHTglob_W_' + str(seasonalizedt.values()) + '_' + calctype + '_inscript2.pdf')

flddf.plot(x=lat)
plt.plot(lat,fldtot,color='k')
plt.axhline(y=0,color='0.5')
plt.ylim(ylimsW)
plt.xlim((-15,75))
plt.xlabel('latitude')
plt.ylabel('poleward transport (W)')
plt.title(str(seasonalizedt) + ' ' + calctype + ' done in script')
if printtofile:
    plt.savefig('PHTnh_W_' + str(seasonalizedt.values()) + '_' + calctype + '_inscript2.pdf')

(flddf*W2calperday/1e19).plot(x=lat)
plt.plot(lat,fldtot*W2calperday/1e19,color='k')
plt.axhline(y=0,color='0.5')
plt.ylim(ylimscal)
plt.ylabel('poleward transport 1e19 cal/day()')
plt.xlim((-15,75))
plt.title(str(seasonalizedt) + ' ' + calctype + ' done in script')
if printtofile:
    plt.savefig('PHTnh_calpday_' + str(seasonalizedt.values()) + '_' + calctype + '_inscript2.pdf')

# ------ zonal AVERAGE
fldzadf.plot(x=lat)
plt.plot(lat,fldzatot,color='k')
plt.axhline(y=0,color='0.5')
plt.ylim(ylimszaW)
plt.xlim((-15,75))
plt.xlabel('latitude')
plt.ylabel('(zonal average) poleward transport (W)')
plt.title(fldpre + ' ' + str(seasonalizedt) + ' calc in script')
if printtofile:
    plt.savefig('PHTnh_zonavg_W_' + str(seasonalizedt.values()) + calctype + '_inscript2.pdf')
# --------------


# =================
# Now test vertical int / avg outputs


if fldpre == '_vert_int':
    flddt2 = flddt
else:
    flddt2={'MPEF': 'DIV',
            'MSHF': 'DIV',
            'MLHF': 'DIV',
            'MFZM': 'DIV'}

fdict2={'MPEF': basepath + casename + fldpre + '_meridional_potential_energy_flux' + fldsuff +'_' + timeper + '_ts.nc',
       'MSHF': basepath + casename + fldpre + '_meridional_sens_heat_flux' + fldsuff + '_' + timeper + '_ts.nc',
       'MLHF': basepath + casename + fldpre + '_meridional_latent_heat_flux' + fldsuff + '_' + timeper + '_ts.nc',
       'MFZM': basepath + casename + fldpre + '_meridional_flux_zonal_mom' + fldsuff + '_' + timeper + '_ts.nc'}

fldintdt2={}; fldintzadt2={}
for fkey in fdict2:
    
    ncfield=flddt2[fkey]
    
    fname=fdict2[fkey]
    conv=condt[fkey] 

    fld=cnc.getNCvar(fname,ncfield,timesel=timesel)*conv

    fldint=np.zeros(len(lat))
    fldintza=np.zeros(len(lat))
    latidx=0
    for ll in lat: # do all lats
        latrad = np.deg2rad(ll)

        fldsub = np.squeeze(fld[:,latidx,:-1]) # get just one lat and remove extra lon
        fldtm = np.mean(cutl.seasonalize_monthlyts(fldsub,**seasonalizedt),axis=0) # time mean
    
        fldint[latidx] = calc_PHT2(fldtm, latidx, nlon, calc='none') # zonal integration
        fldintza[latidx] = calc_PHT3(fldtm, latidx, nlon, calc='none') # zonal average

        latidx += 1
    
    fldintdt2[fkey] = fldint
    fldintzadt2[fkey] = fldintza


flddf2 = pd.DataFrame(fldintdt2)
fldtot2 = flddf2.sum(axis=1)
fldzadf2 = pd.DataFrame(fldintzadt2)
fldzatot2 = fldzadf2.sum(axis=1)


flddf2.plot(x=lat)
plt.plot(lat,fldtot2,color='k')
plt.axhline(y=0,color='0.5')
plt.ylabel('poleward transport (W)')
plt.xlim((-80,80))
plt.title(fldpre+ ' ' + str(seasonalizedt) + ' calc in file')
if printtofile:
    plt.savefig('PHTglob_W_' + str(seasonalizedt.values()) + fldpre + '_infile.pdf')

flddf2.plot(x=lat)
plt.plot(lat,fldtot2,color='k')
plt.axhline(y=0,color='0.5')
#plt.ylim((-1e15,1e15))
plt.xlim((-15,75))
plt.xlabel('latitude')
plt.ylabel('poleward transport (W)')
plt.title(fldpre + ' ' + str(seasonalizedt) + ' calc in file')
if printtofile:
    plt.savefig('PHTnh_W_' + str(seasonalizedt.values()) + fldpre + '_infile.pdf')

(flddf2*W2calperday/1e19).plot(x=lat)
plt.plot(lat,fldtot2*W2calperday/1e19,color='k')
plt.axhline(y=0,color='0.5')
plt.ylim(ylimscal)
plt.xlim((-15,75))
plt.ylabel('(1e19 cal/day)')
plt.title(fldpre + ' ' + str(seasonalizedt) + ' calc in file')
if printtofile:
    plt.savefig('PHTglob_calpday_' + str(seasonalizedt.values()) + fldpre + '_infile.pdf')


# ------ zonal AVERAGE
fldzadf2.plot(x=lat)
plt.plot(lat,fldzatot2,color='k')
plt.axhline(y=0,color='0.5')
plt.ylim(ylimszaW2)
plt.xlim((-15,75))
plt.xlabel('latitude')
plt.ylabel('(zonal average) poleward transport (W)')
plt.title(fldpre + ' ' + str(seasonalizedt) + ' calc in file')
if printtofile:
    plt.savefig('PHTnh_zonavg_W_' + str(seasonalizedt.values()) + fldpre + '_infile.pdf')
# --------------


# vertint:====
# vert_int is smoother, looks much better. But magnitude is possible
#   x10 off? Yes, ANN total northward flux around 45N is maybe 3PW
#   in Oort. My figure shows maybe ~0.4PW.
# ACTUALLY, need it to be ZONAL AVERAGE: when comparing zonal average values:
#    the calc done in nc file is x100 smaller. This may be in units of pressure??
#    And here the calc in this script is off by x10. (compare to Fasullo & Trenberth for example)

# vert_avg looks like it still might have some elevation effects: bit
#   jaggety in similar spots to previous problem areas.
# ZONAL AVG and vert avg are very similar for calc in nc file and in script. 
# same as zonal integral.
