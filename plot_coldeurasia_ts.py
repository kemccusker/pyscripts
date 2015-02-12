import cccmaplots as cplt
import cccmaNC as cnc
import cccmautils as cutl


cnc=reload(cnc)
cutl=reload(cutl)

printtofile=True
plt.close('all')

region='eurasia'
sea='ANN'

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
if printtofile:
    fig.savefig(field + '_' + region + '_monanom_' + str(styr) + '-2014_timeseries.pdf')


fig,ax = plt.subplots(1,1)
ax.plot(xx,fldgm)
ax.set_xticks(np.arange(6,len(fldreg),120)) # ticks for each decade, in July
ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,10))
ax.set_title('Global avg (monthly anom)')
if printtofile:
    fig.savefig(field + '_globavg_monanom_' + str(styr) + '-2014_timeseries.pdf')

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
if printtofile: # with standard deviation lines
    fig.savefig(field + '_' + region + '_monanom_' + str(styr) + '-2014_timeseries2.pdf')


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
if printtofile:
    fig.savefig(field + '_' + region + '_' + sea + 'anom_' + str(styr) + '-2014_timeseries2.pdf')


nlat=len(lat)
nlon=len(lon)

fldrs = fld.reshape((len(xx),nlat*nlon))
# ############# TREND
slope,intercept = np.polyfit(xx,fldrs,1) # seasonal trend

slope=slope.reshape((nlat,nlon))
gmtr = cutl.global_mean_areawgted(slope*len(xx),lat,lon)

ttl = '1979-2014 ' + sea + ' trend: gm= ' + '$%.2f$'%(gmtr) + ' per ' + str(len(xx)) + 'yr'


ptype='sq'
fig,ax=plt.subplots(1,1)
cplt.kemmap(slope*len(xx),lat,lon,title=ttl,axis=ax,cmin=-2,cmax=2,
            cmap='blue2red_20',drawgrid=True,type=ptype)
if printtofile:
    fig.savefig(field + 'trend_' + str(styr) + '-2014_' + sea + '_' + ptype + '.pdf')



########  TREND SINCE 2000
timesel2='2000-01-01,2014-12-31'
fld2 = cnc.getNCvar(fname,'tempanomaly',timesel=timesel2, seas=sea)
fldreg2 = cutl.calc_regmean(fld2,lat,lon,region=region)

xx2=np.arange(0,len(fldreg2))
fldrs2 = fld2.reshape((len(xx2),nlat*nlon))



slope,intercept = np.polyfit(xx2,fldrs2,1) # seasonal trend

slope=slope.reshape((nlat,nlon))

gmtr = cutl.global_mean_areawgted(slope*len(xx2),lat,lon)

ttl = '2000-14 ' + sea + ' trend: gm= ' + '$%.2f$'%(gmtr) + ' per ' + str(len(xx2)) + 'yr'

ptype='sq'
fig,ax=plt.subplots(1,1)
cplt.kemmap(slope*len(xx2),lat,lon,title=ttl,axis=ax,cmin=-3,cmax=3,
            cmap='blue2red_20',drawgrid=True,type=ptype)

if printtofile:
    fig.savefig(field + 'trend_2000-2014_' + sea + '_' + ptype + '.pdf')
