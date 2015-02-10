import cccmaplots as cplt
import cccmaNC as cnc
import cccmautils as cutl


cnc=reload(cnc)
cutl=reload(cutl)

region='eurasia'

basepath = '/Volumes/MyPassport2TB/DATA/OBSERVATIONS/'
file = 'gistemp1200_ERSST.nc'

fname = basepath + file

fld = cnc.getNCvar(fname,'tempanomaly', seas='ND')
lat = cnc.getNCvar(fname,'lat')
lon = cnc.getNCvar(fname,'lon')

dates = cnc.get_NCdates(fname)


fldreg = cutl.calc_regmean(fld,lat,lon,region=region)

