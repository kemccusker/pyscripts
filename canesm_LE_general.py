""" canesm_LE_general.py

"""
import pandas as pd
import cccmaplots as cplt
import datetime as datetime
from netCDF4 import num2date
import sys
import loadLE as le

timesel1='1979-01-01,1989-12-31'
timesel2='2002-01-01,2012-12-31'

conv=1
region='eurasiamori'
sea='DJF'


field='ts'; ncfield='ts'; comp='Amon'
ftype='climo' # or 'fullts'


fdict = {'field': field, 'ncfield': ncfield, 'comp': comp}

histcdat=le.load_LEdata(fdict,'historical',timesel=timesel1, rettype='ndarray',conv=conv,ftype=ftype)
(numens,ntime,nlatlon) = histcdat.shape
histpdat=le.load_LEdata(fdict,'historical',timesel=timesel2, rettype='ndarray',conv=conv,ftype=ftype)

histnatcdat=le.load_LEdata(fdict,'historicalNat',timesel=timesel1, rettype='ndarray',conv=conv,ftype=ftype)
#(numens,ntime,nlatlon) = histnatcdat.shape
histnatpdat=le.load_LEdata(fdict,'historicalNat',timesel=timesel2, rettype='ndarray',conv=conv,ftype=ftype)


flist=le.build_filenames(fdict,'historical',ftype='fullts')
numens=len(flist)
fname = flist[0]
lat=cnc.getNCvar(fname,'lat')
nlat=len(lat)
lon=cnc.getNCvar(fname,'lon')
nlon=len(lon)

histcdat = histcdat.reshape((numens,ntime,nlat,nlon))
histpdat = histpdat.reshape((numens,ntime,nlat,nlon))

histnatcdat = histnatcdat.reshape((numens,ntime,nlat,nlon))
histnatpdat = histnatpdat.reshape((numens,ntime,nlat,nlon))

# Now, seasonal average and regional average here
nyr = ntime/12.-1 # @@@@@@ hack for laptop. data probably not correct @@@@ actually needed on linux too...

#histcsea = np.zeros((numens,nyr,nlat,lon)
histreg = np.zeros((numens,nyr))
histnatreg = np.zeros((numens,nyr))

for nii in range(0,numens):

    histcsea = cutl.seasonalize_monthlyts(histcdat[nii,...],season=sea)
    histpsea = cutl.seasonalize_monthlyts(histpdat[nii,...],season=sea)
    ctl = cutl.calc_regmean(histcsea,lat,lon,region=region)
    pert = cutl.calc_regmean(histpsea,lat,lon,region=region)
    histreg[nii,...] = pert - ctl

    histnatcsea = cutl.seasonalize_monthlyts(histnatcdat[nii,...],season=sea)
    histnatpsea = cutl.seasonalize_monthlyts(histnatpdat[nii,...],season=sea)
    ctlnat = cutl.calc_regmean(histnatcsea,lat,lon,region=region)
    pertnat = cutl.calc_regmean(histnatpsea,lat,lon,region=region)
    histnatreg[nii,...] = pertnat - ctlnat

plt.figure()
plt.hist(histreg.mean(axis=1),alpha=0.5)
plt.hist(histnatreg.mean(axis=1),color='.5')
