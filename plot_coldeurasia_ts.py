import cccmaplots as cplt
import cccmaNC as cnc
import cccmautils as cutl


cnc=reload(cnc)
cutl=reload(cutl)

plt.close('all')

region='eurasia'
sea='JF'

#basepath = '/Volumes/MyPassport2TB/DATA/OBSERVATIONS/'
basepath = '/ra40/data/ncs/reanalyses/'
file = 'gistemp1200_ERSST.nc'
# base is 1951-1980

fname = basepath + file

styr=1979

timesel=str(styr) + '-01-01,2014-12-31'
fld = cnc.getNCvar(fname,'tempanomaly',timesel=timesel)
lat = cnc.getNCvar(fname,'lat')
lon = cnc.getNCvar(fname,'lon')

dates = cnc.get_NCdates(fname)


fldreg = cutl.calc_regmean(fld,lat,lon,region=region)
fldgm = cutl.global_mean_areawgted3d(fld,lat,lon)


# ####### plotting #######
xx=np.arange(0,len(fldreg))

fig,ax = plt.subplots(1,1)
ax.plot(xx,fldreg)
ax.set_xticks(np.arange(6,len(fldreg),120)) # ticks for each decade, in July
ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,10))
ax.set_title(region + ' (monthly anom)')

fig,ax = plt.subplots(1,1)
ax.plot(xx,fldgm)
ax.set_xticks(np.arange(6,len(fldreg),120)) # ticks for each decade, in July
ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,10))
ax.set_title('Global avg (monthly anom)')

fldrstd=fldreg.std()
fldrm=fldreg.mean()

fig,ax = plt.subplots(1,1)
ax.plot(xx,fldreg)
ax.axhline(fldrm+fldrstd,linewidth=.5,color='k')
ax.axhline(fldrm-fldrstd,linewidth=.5,color='k')    
ax.set_xticks(np.arange(6,len(fldreg),120)) # ticks for each decade, in July
ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,10))
ax.set_title(region + ' (monthly anom)')

# ######### Sea mean ###

fld = cnc.getNCvar(fname,'tempanomaly',timesel=timesel, seas=sea)
fldreg = cutl.calc_regmean(fld,lat,lon,region=region)

xx=np.arange(0,len(fldreg))

fldrstd=fldreg.std()
fldrm=fldreg.mean()

fig,ax = plt.subplots(1,1)
ax.plot(xx,fldreg)
ax.axhline(fldrm+fldrstd,linewidth=.5,color='k')
ax.axhline(fldrm-fldrstd,linewidth=.5,color='k')    
ax.set_xticks(np.arange(0,len(fldreg),5)) # ticks for each decade, in July
#ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,5))
ax.set_title(region + ' (' + sea + ' anom)')


nlat=len(lat)
nlon=len(lon)

fldrs = fld.reshape((len(xx),nlat*nlon))
# TREND
slope,intercept = np.polyfit(xx,fldrs,1) 

slope=slope.reshape((nlat,nlon))

fig,ax=plt.subplots(1,1)
cplt.kemmap(slope,lat,lon,title='trend',axis=ax,cmin=-.2,cmax=.2,cmap='blue2red_20',drawgrid=True)
