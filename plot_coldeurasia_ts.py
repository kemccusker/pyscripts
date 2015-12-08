import cccmaplots as cplt
import cccmaNC as cnc
import cccmautils as cutl
import constants as con
import cccmacmaps as ccm
import platform as platform

cnc=reload(cnc)
cutl=reload(cutl)

euranom=False
nhanom=False # if both False and temptimeseries=True, do eur - nh
glob=False # if True, then write_temptimeseries will write global mean data
bksanom=False # for Z500/circ index: if True and euranom=False, write out bksmori Z500 timeseries

write_temptimeseries=False # write the temperature timeseries to netcdf (GISTEMP)
write_temptimeseriesera=False # write erainterim's ST instead
write_siatimeseries=False # write the sea ice area timeseries to netcdf
write_circtimeseries=True # write circulation index time series. if anom flags all False, do sicregion-region

write_tempmap=False
write_simsicmap=False
write_nsidcsicmap=False
write_simsatmap=False; casename='E4' # E4 cold, E1 warm
write_simz500map=False


printtofile=False
plt.close('all')

styr=1979
timesel=str(styr) + '-01-01,2014-12-31'
timesel=str(styr) + '-01-01,2015-07-01'

region='eurasiamori'
sicregion='bksmori'
sea='DJF'
sicsea='DJF'
field='st'
units='$^\circ$ C'

plat = platform.system()   
if plat == 'Darwin':  # means I'm on my mac
    obasepath = '/Volumes/MyPassport2TB/DATA/OBSERVATIONS/'
    obasepath2=obasepath

else:  # on linux workstation in Vic
    obasepath = '/HOME/rkm/work/BCs/' # obs basepath (boundary conditions)
    obasepath2 = '/HOME/rkm/work/DATA/' # other obs, like GISS


#td_nsidc_merged_197811_latest_128_64_sicn_1978111600-2013121612.nc
filesic = obasepath + 'NSIDC/td_bootstrap_197811_latest_128_64_sicn_1978111600-2013121612.nc'
#filesic='/HOME/rkm/work/BCs/NSIDC/nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc'
#filesic='/Volumes/MyPassport2TB/DATA/OBSERVATIONS/nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc'
latsic=cnc.getNCvar(filesic,'lat')
lonsic=cnc.getNCvar(filesic,'lon')


#basepath = '/Volumes/MyPassport2TB/DATA/OBSERVATIONS/'
#basepath = '/HOME/rkm/work/DATA/GISS/'
gisfile = 'GISS/gistemp1200_ERSST.nc'
gisfile = 'GISS/td_giss_tsurf1200_188001_201509_128_64_st_1880011612-2015091600.nc'
gisfile = 'GISS/td_giss_tsurf1200_188001_201509_fill_128_64_st_1880011612-2015091600.nc'
# base is 1951-1980

fname = obasepath2 + gisfile

#satfield='tempanomaly'
satfield='ST' # for the model-prepared file
#fld = cnc.getNCvar(fname,'tempanomaly',timesel=timesel)
fld = cnc.getNCvar(fname,satfield,timesel=timesel)

lat = cnc.getNCvar(fname,'lat')
lon = cnc.getNCvar(fname,'lon')
gistime = cnc.getNCvar(fname,'time',timesel=timesel)

dates = cnc.get_NCdates(fname)

# ADD Simulated obs == NSIDC
filecnsidc,filepnsidc=con.build_filepathpair('NSIDC','st')

fldreg = cutl.calc_regmean(fld,lat,lon,region=region,model=None)
fldgm = cutl.global_mean_areawgted3d(fld,lat,lon,model=None)
fldnh = cutl.calc_regmean(fld,lat,lon,region='nh',model=None)


# ####### plotting timeseries #######
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

fldcsim = cnc.getNCvar(filecnsidc,'ST',seas=sea,timesel='002-01-01,121-12-31')
latsim = cnc.getNCvar(filecnsidc,'lat')
lonsim = cnc.getNCvar(filecnsidc,'lon')

fldpsim = cnc.getNCvar(filepnsidc,'ST',seas=sea,timesel='002-01-01,121-12-31')
tstatmapsim,pvalmapsim=cutl.ttest_ind(fldcsim,fldpsim)

fldregsim = cutl.calc_regmean(fldpsim-fldcsim,latsim,lonsim,region=region,model=None)
plotregsim = fldpsim.mean(axis=0)-fldcsim.mean(axis=0)
regsimm = fldregsim.mean()
regsimstd=fldregsim.std()


fld = cnc.getNCvar(fname,satfield,timesel=timesel, seas=sea)
fldreg = cutl.calc_regmean(fld,lat,lon,region=region,model=None)
gisreg=fldreg # write to file
gisnh = cutl.calc_regmean(fld,lat,lon,region='nh',model=None) # seasonal
gisgm = cutl.calc_regmean(fld,lat,lon,region='gm',model=None)

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

fig,ax = plt.subplots(1,1)
ax.plot(xx,gisnh)
ax.set_xticks(np.arange(0,len(gisnh),5)) # ticks for each decade, in July
#ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,5))
ax.set_title('NH (' + sea + ' anom)')
if printtofile:
    fig.savefig(field + '_nh_' + sea + 'anom_' + str(styr) + '-2014_timeseries.pdf')

fig,ax = plt.subplots(1,1)
ax.plot(xx,fldreg-gisnh)
ax.set_xticks(np.arange(0,len(gisnh),5)) # ticks for each decade, in July
#ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,5))
ax.set_title(region + '-NH (' + sea + ' anom)')
if printtofile:
    fig.savefig(field + '_' + region + '-nh_' + sea + 'anom_' + str(styr) + '-2014_timeseries.pdf')



plt.figure() # normalize by stddev
plt.plot(gisnh/np.std(gisnh)); 
plt.plot(gisreg/np.std(gisreg),'k')
plt.plot( (gisreg-gisnh)/np.std(gisreg-gisnh),'r')
plt.legend(('NH',region,region+'-NH'),loc='upper left')
plt.title('SAT norm by stdev')



def runmean(input, window=5,axis=0) :
    ret = np.cumsum(input, dtype=float,axis=axis)
    ret[window:,...] = ret[window:,...] - ret[:-window,...]
    return ret[window - 1:,...] / np.float(window)

window=5
xxrun = xx[window/2:-window/2+1]
fig,ax=plt.subplots(1,1)
ax.plot(xxrun, runmean((gisgm)/np.std(gisgm),window=window),'0.1',linewidth=2)
ax.plot(xxrun,runmean(gisnh/np.std(gisnh),window=window),color='0.5',linewidth=2); 
ax.plot(xxrun,runmean(gisreg/np.std(gisreg),window=window),'b',linewidth=2)
ax.plot(xx, (gisgm)/np.std(gisgm),'0.1',linewidth=1,alpha=0.5)
ax.plot(xx,gisnh/np.std(gisnh),color='0.5',linewidth=1,alpha=0.5); 
ax.plot(xx,gisreg/np.std(gisreg),'b',linewidth=1,alpha=0.5)
ax.legend(('GM','NH','Eurasia'),loc='upper left',frameon=False)
ax.set_title(sea + ' Normalized SAT ($\sigma$) + ' + str(window) + ' yr smooth')
ax.set_xticks(np.arange(0,len(gisnh),5)) # ticks for each decade, in July
#ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,5))
if printtofile:
    fig.savefig(field + 'gm_' + sea + 'anom_' + str(styr) + '-2014_run' + str(window)+'timeseries.pdf')


fig,ax=plt.subplots(1,1)
ax.plot(xxrun,runmean(gisnh/np.std(gisnh),window=window),color='orange',linewidth=2); 
ax.plot(xxrun,runmean(gisreg/np.std(gisreg),window=window),'b',linewidth=2)
ax.plot(xxrun, runmean((gisgm)/np.std(gisgm),window=window),'r',linewidth=2)
ax.legend(('NH',region,'GM'),loc='upper left',frameon=False)
ax.set_title(sea + ' Normalized SAT ($\sigma$), ' + str(window) + ' yr smooth')
ax.set_xticks(np.arange(0,len(gisnh),5)) 
#ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,5))


fldsic = cnc.getNCvar(filesic,'SICN',timesel=timesel,seas=sicsea)
fldsicreg = cutl.calc_regtotseaicearea(fldsic,latsic,lonsic,region=sicregion,model=None)#isarea=False
nsidcreg = fldsicreg # to write to file

sicrstd=fldsicreg.std()
sicrm = fldsicreg.mean()

xxs=np.arange(0,len(fldsicreg))
fig,ax = plt.subplots(1,1)
ax.plot(xxs,fldsicreg)
ax.axhline(sicrm+sicrstd,linewidth=.5,color='k')
ax.axhline(sicrm-sicrstd,linewidth=.5,color='k')    
ax.set_xticks(np.arange(0,len(fldsicreg),5)) # ticks for each decade, in July
#ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2011,5))
ax.set_title(sicregion + ' (' + sicsea + ' anom)')
if printtofile:
    fig.savefig( 'sicn_' + sicregion + '_' + sicsea + '_' + str(styr) + '-2011_timeseries2.pdf')



# ######### SAT and SIC timeseries ###
fig,axs=plt.subplots(2,1)
fig.set_size_inches(8,6.5)
ax=axs[0]
ax.plot(xx,fldreg)
ax.axhline(fldrm+fldrstd,linewidth=.5,color='k')
ax.axhline(fldrm-fldrstd,linewidth=.5,color='k')    
ax.set_xticks(np.arange(0,len(fldreg),5)) # ticks for each decade, in July
#ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2014,5))
ax.set_title(sea + ' Eurasia SAT anomaly ($^\circ$C)')# + region)
#ax.set_title(sea)

ax=axs[1]
ax.plot(xxs,fldsicreg/1e12)
ax.axhline((sicrm+sicrstd)/1e12,linewidth=.5,color='k')
ax.axhline((sicrm-sicrstd)/1e12,linewidth=.5,color='k')    
ax.set_xticks(np.arange(0,len(fldsicreg),5)) # ticks for each decade, in July
#ax.minorticks_on()
ax.set_xticklabels(np.arange(styr,2011,5))
ax.set_title(sicsea + ' Barents-Kara SIC (e$^{12}$)')# + sicregion)
ax.set_xlabel('years')
if printtofile:
    fig.savefig(field + '_' + region + 'anom_sicn_' + sicregion + 
                '_' + sicsea + sea + '_' + str(styr) + '-2014_2011_timeseries2.pdf')

siccol=ccm.get_linecolor('yelloworange')

import matplotlib.lines as mlines

sat=mlines.Line2D([],[],color='k',linewidth=2)
sic=mlines.Line2D([],[],color=siccol,linewidth=2)

# #### SAT and SIC timeseries on same fig (twinx)
fig,axl=plt.subplots(1,1)
fig.set_size_inches(8,3)
satlin = axl.plot(xx,fldreg,'k',linewidth=2,label='SAT')
axr=axl.twinx()
siclin = axr.plot(xxs,fldsicreg/1e12,color=siccol,linewidth=2,label='SIC')
axr.set_ylabel(r'Barents-Kara SIC (e$^{12}$)')
axl.set_ylabel(r'$\Delta$ Eurasian SAT ($^\circ$C)')
axl.set_xticklabels(np.arange(styr,2014,5))
axl.set_title('Dec-Jan-Feb')
hl,ll = axl.get_legend_handles_labels()
hr,lr = axr.get_legend_handles_labels()

axl.legend((sat,sic),('SAT','SIC'),frameon=False)
if printtofile:
    fig.savefig(field + '_' + region + 'anom_sicn_' + sicregion + 
                '_' + sicsea + sea + '_' + str(styr) + '-2013_timeseries_onepanel.pdf')





nlat=len(lat)
nlon=len(lon)

fldrs = fld.reshape((len(xx),nlat*nlon))
# ############# TREND
slope,intercept = np.polyfit(xx,fldrs,1) # seasonal trend

slope=slope.reshape((nlat,nlon))
gmtr = cutl.global_mean_areawgted(slope*len(xx),lat,lon,model=None)

# 1D trend:
fldregfull = cutl.calc_regmean(fld,lat,lon,region=region,model=None)

mm, bb, rval, pval, std_err = sp.stats.linregress(xx,fldregfull)
print region + ' avg and ' + sea + ' avg trend (' + timesel + ') pval: ' + str(pval)


ttl = '1979-2014 ' + sea + ' trend: gm= ' + '$%.2f$'%(gmtr) + ' per ' + str(len(xx)) + 'yr'


ptype='sq'
fig,ax=plt.subplots(1,1)
cplt.kemmap(slope*len(xx),lat,lon,title=ttl,axis=ax,cmin=-2,cmax=2,
            cmap='blue2red_20',drawgrid=True,ptype=ptype)
if printtofile:
    fig.savefig(field + 'trend_' + str(styr) + '-2014_' + sea + '_' + ptype + '.pdf')

# sea ice 1D trend
mm, bb, rval, pval, std_err = sp.stats.linregress(xxs,fldsicreg)
print sicregion + ' avg and ' + sicsea + ' avg trend (' + timesel + ') pval: ' + str(pval)



# ######### MOVING TREND ##############
import numpy.ma as ma

mm=np.zeros(fldregfull.shape)
pval=np.zeros(fldregfull.shape)
for xxi in xx:
    print 'start yr: ' + str(1979+xxi)
    print 'trend length: ' + str(xx[xxi:].shape)
    mm[xxi], bb, rval, pval[xxi], std_err = sp.stats.linregress(xx[xxi:],fldregfull[xxi:])
    
mmsic=np.zeros(fldsicreg.shape)
pvalsic=np.zeros(fldsicreg.shape)
for xxi in xxs:
    print 'start yr: ' + str(1979+xxi)
    print 'trend length: ' + str(xxs[xxi:].shape)
    mmsic[xxi], bb, rval, pvalsic[xxi], std_err = sp.stats.linregress(xxs[xxi:],fldsicreg[xxi:])


mmma = ma.masked_where(pval>0.1,mm)
mmsicma = ma.masked_where(pvalsic>0.1,mmsic)

yrs=np.arange(1979,1979+len(xx))    
plt.figure()
plt.plot(yrs,mm,'k',linewidth=1)
plt.plot(yrs,mmma,'k',linewidth=4)


fig,axl=plt.subplots(1,1)
fig.set_size_inches(8,3)
satlin = axl.plot(xx,mm,'k',linewidth=1,label='SAT')
axl.plot(xx,mmma,'k',linewidth=4)
axr=axl.twinx()
axl.axhline(y=0,color='k',linestyle='--')
axr.axhline(y=0,color=siccol,linestyle='--')
axl.legend((sat,sic),('SAT','SIC'),frameon=False,loc='upper left')
siclin = axr.plot(xxs,mmsic/1e12,color=siccol,linewidth=1,label='SIC')
axr.plot(xxs,mmsicma/1e12,color=siccol,linewidth=4)
axr.set_ylabel(r'Barents-Kara SIC (e$^{12}$/yr)')
axl.set_ylabel(r'Eurasian SAT ($^\circ$C/yr)')
axl.set_xticklabels(np.arange(styr,2014,5))
axl.set_title('Dec-Jan-Feb Moving Start Yr')
hl,ll = axl.get_legend_handles_labels()
hr,lr = axr.get_legend_handles_labels()

if printtofile:
    fig.savefig(field + '_' + region + '_sicn_' + sicregion + 
                '_' + sicsea + sea + '_1979on_movingstyrtrend_onepanel.pdf')



########  TREND SINCE 2000
timesel2='2000-01-01,2014-12-31'
fld2 = cnc.getNCvar(fname,satfield,timesel=timesel2, seas=sea)
fldreg2 = cutl.calc_regmean(fld2,lat,lon,region=region,model=None)


xx2=np.arange(0,len(fldreg2))
fldrs2 = fld2.reshape((len(xx2),nlat*nlon))
# 1D trend:
mm, bb, rval, pval, std_err = sp.stats.linregress(xx2,fldreg2)
print region + ' avg and ' + sea + ' avg trend (' + timesel2 + ') pval: ' + str(pval)

slope,intercept = np.polyfit(xx2,fldrs2,1) # seasonal trend

slope=slope.reshape((nlat,nlon))

gmtr = cutl.global_mean_areawgted(slope*len(xx2),lat,lon,model=None)

ttl = '2000-14 ' + sea + ' trend: gm= ' + '$%.2f$'%(gmtr) + ' per ' + str(len(xx2)) + 'yr'

ptype='sq'
fig,ax=plt.subplots(1,1)
cplt.kemmap(slope*len(xx2),lat,lon,title=ttl,axis=ax,cmin=-3,cmax=3,
            cmap='blue2red_20',drawgrid=True,ptype=ptype)
if printtofile:
    fig.savefig(field + 'trend_2000-2014_' + sea + '_' + ptype + '.pdf')


#  # SIC
sic2=cnc.getNCvar(filesic,'SICN',timesel=timesel2,seas=sicsea)
sic2rm = cutl.calc_regtotseaicearea(sic2,latsic,lonsic,region=sicregion,model=None)
xxs2=np.arange(0,len(sic2rm))

mm, bb, rval, pval, std_err = sp.stats.linregress(xxs2,sic2rm)
print sicregion + ' avg and ' + sicsea + ' avg trend (' + timesel2 + ') pval: ' + str(pval)



######## EPOCH differences
cmin=-1; cmax=1 # -1,1 match simulated

print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'

timesel1='1979-01-01,1989-12-31' # @@ adding in another yr to get last DJF
timesel2='2002-01-01,2012-12-31'
print timesel1,timesel2
fld1 = cnc.getNCvar(fname,satfield,timesel=timesel1, seas=sea)
fld2 = cnc.getNCvar(fname,satfield,timesel=timesel2, seas=sea)
tstatmap,pvalmap = cutl.ttest_ind(fld1,fld2)

fld1regm= cutl.calc_regmean(fld1,lat,lon,region,model=None)
fld2regm = cutl.calc_regmean(fld2,lat,lon,region,model=None)
# SIGNIFICANCE b/w two means:
tstat,pval = cutl.ttest_ind(fld2regm,fld1regm)
print 'EPOCH diffs for ' + region + ' ' + sea + ' avg tstat,pval: '+ str(tstat) +',' + str(pval)
print '  --- mean val: ' + str(fld2regm.mean(axis=0)-fld1regm.mean(axis=0))
print ' (expect fld2 > fld1); can use 1-sided pval (pval/2) and tstat>0' # (see below)
print ' tstat < 0 so fld2!>fld1, but pval is still only 0.3, so not signif'
print ''
fld1nhm= cutl.calc_regmean(fld1,lat,lon,'nh',model=None)
fld2nhm = cutl.calc_regmean(fld2,lat,lon,'nh',model=None)
tstat,pval = cutl.ttest_ind(fld2nhm,fld1nhm)
print 'EPOCH diffs for nh ' + sea + ' avg tstat,pval: '+ str(tstat) +',' + str(pval)
print '  --- mean val: ' + str(fld2nhm.mean(axis=0)-fld1nhm.mean(axis=0))
print ' (expect fld2 > fld1); can use 1-sided pval (pval/2) and tstat>0' # (see below)

tstat,pval = cutl.ttest_ind(fld2regm-fld2nhm,fld1regm-fld1nhm)
print 'EPOCH diffs for ' + region + ' ' + sea + '-NH: avg tstat,pval: ' + str(tstat) +',' +str(pval)
print '  --- mean val: ' + str((fld2regm-fld2nhm).mean(axis=0)-(fld1regm-fld1nhm).mean(axis=0))

"""
http://stackoverflow.com/questions/15984221/how-to-perform-two-sample-one-tailed-t-test-with-numpy-scipy

because the one-sided tests can be backed out from the two-sided tests. 
(With symmetric distributions one-sided p-value is just half of the two-sided pvalue)...

It goes on to say that scipy always gives the test statistic as signed. 
This means that given p and t values from a two-tailed test, you would 
reject the null hypothesis of a greater-than test when p/2 < alpha and t > 0, 
and of a less-than test when p/2 < alpha and t < 0.

"""

fldregtm = fld2regm.mean(axis=0) - \
         fld1regm.mean(axis=0)
fldregnhtm = (fld2regm-fld2nhm).mean(axis=0) -\
             (fld1regm-fld1nhm).mean(axis=0)

plotfld = fld2.mean(axis=0) - fld1.mean(axis=0)
gismap=plotfld # for write to file

# Now get observed SIC (NSIDC)
#basepath='/HOME/rkm/work/BCs/NSIDC/'
#file1='nsidc_bt_128x64_1978m11_2011m12_sic_1979-1989climo.nc'
#file2='nsidc_bt_128x64_1978m11_2011m12_sic_2002-2011climo.nc'
#filets='nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc'

sic1=cnc.getNCvar(filesic,'SICN',timesel=timesel1,seas=sicsea)
sic2=cnc.getNCvar(filesic,'SICN',timesel=timesel2,seas=sicsea)
sia1regm = cutl.calc_regtotseaicearea(sic1,latsic,lonsic,region=sicregion,model=None)
sia2regm = cutl.calc_regtotseaicearea(sic2,latsic,lonsic,region=sicregion,model=None)

sic1regm= cutl.calc_regmean(sic1,latsic,lonsic,sicregion,model=None)
sic2regm = cutl.calc_regmean(sic2,latsic,lonsic,sicregion,model=None)
# SIGNIFICANCE b/w two means:
tstat,pval = cutl.ttest_ind(sic1regm,sic2regm)
print 'EPOCH diffs for SIC ' + sicregion + ' ' + sicsea + ' avg tstat,pval: '+ str(tstat) +',' + str(pval)
print '  --- mean val: ' + str(sic2regm.mean(axis=0)-sic1regm.mean(axis=0))

# SIGNIFICANCE b/w two means:
tstat,pval = cutl.ttest_ind(sia1regm,sia2regm)
print 'EPOCH diffs for SIA ' + sicregion + ' ' + sicsea + ' avg tstat,pval: '+ str(tstat) +',' + str(pval)
print '  --- mean val: ' + str(sia2regm.mean(axis=0)-sia1regm.mean(axis=0))



cminsic=-.3; cmaxsic=.21
incr = (cmaxsic-cminsic) / (20.)
contssic = np.arange(cminsic,cmaxsic+incr,incr)
contssic = contssic[::2]
#contssic = np.arange(-.3,.21,.08) # sic contours

plotsic = sic2.mean(axis=0) - sic1.mean(axis=0)

lons, lats = np.meshgrid(lonsic,latsic)

# Observed temp anomaly with Observed (NSIDC) sic anomaly =====================
ptype='eabkslamb' # 'eabksstere'
fig,ax=plt.subplots(1,1)
bm,pc=cplt.kemmap(plotfld,lat,lon,ptype=ptype,axis=ax,cmin=cmin,cmax=cmax,cmap='blue2red_20')
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
cplt.addtsigm(bm,pvalmap,lat,lon,siglevel=0.1)
if printtofile:
    fig.savefig(field + '_' + sea+ 'anomsig_sicncont_' + sicsea + '_' + ptype + '_giss1979-89_2002-12.pdf')

fig,ax=plt.subplots(1,1)
bm,pc=cplt.kemmap(plotfld,lat,lon,ptype='sq',axis=ax,cmin=cmin,cmax=cmax,cmap='blue2red_20')
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
cplt.addtsigm(bm,pvalmap,lat,lon,siglevel=0.1)


# Do both climo SIE contours instead === this one 
#   is no good as contours are on top of e/o in DJF
contssic=[0.15, 0.15]
fig,ax=plt.subplots(1,1)
bm,pc=cplt.kemmap(plotfld,lat,lon,ptype=ptype,axis=ax,cmin=cmin,cmax=cmax,cmap='blue2red_20')
bm.contour(lons,lats,sic1.mean(axis=0),levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
bm.contour(lons,lats,sic2.mean(axis=0),levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='--')



# ### observed temp anomaly w/ NSIDC sea ice AND simulated observed temp anomaly w/ sea ice =====
contssic = np.arange(cminsic,cmaxsic+incr,incr)
contssic = contssic[::2]

fig,axs=plt.subplots(2,1)
fig.set_size_inches((4,8.5))
fig.subplots_adjust(hspace=.02,wspace=.02)
ax=axs[0]
bm,pc=cplt.kemmap(plotfld,lat,lon,ptype=ptype,axis=ax,cmin=cmin,cmax=cmax,cmap='blue2red_20')
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')

ax=axs[1]
bm,pc=cplt.kemmap(plotregsim,latsim,lonsim,ptype=ptype,axis=ax,cmin=cmin,cmax=cmax,cmap='blue2red_20')
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
if printtofile:
    fig.savefig(field + '_' + sea + 'anom_sicncont_' + sicsea + '_' + ptype + '_gissnsidcsim1979-89_2002-12.pdf')


printtofile=False
# # this one is 1 row
fig,axs=plt.subplots(1,2)
fig.set_size_inches((8,4))
fig.subplots_adjust(hspace=.02,wspace=.02)
ax=axs[0]
bm,pc=cplt.kemmap(plotfld,lat,lon,ptype=ptype,axis=ax,cmin=cmin,cmax=cmax,
                  cmap='blue2red_20',lcol='0.3',suppcb=True)
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
cplt.addtsigm(bm,pvalmap,lat,lon,siglevel=0.1)

ax.set_title('a. Observed SAT change')

ax=axs[1]
bm,pc=cplt.kemmap(plotregsim,latsim,lonsim,ptype=ptype,axis=ax,cmin=cmin,cmax=cmax,
                  cmap='blue2red_20',lcol='0.3',suppcb=True)
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
cplt.addtsigm(bm,pvalmapsim,latsim,lonsim,siglevel=0.1)

ax.set_title('b. Simulated SAT change')
cbar_ax = fig.add_axes([.25,0.07, 0.5, .04])
cbar=fig.colorbar(pc,cax=cbar_ax, orientation='horizontal',format='%.1f') 

cbar_ax.tick_params(labelsize=15)
cbar.set_ticks(np.arange(-1,1.2,0.2))
#cbar.set_label('$\Delta$ SAT')
cbar.ax.set_ylabel(units,rotation=0,labelpad=20)
cbar.ax.yaxis.set_label_position('right')
cbar_ax.set_xticklabels(('-1.0','','','','','0','','','','','1.0'))

if printtofile:
    fig.savefig(field + '_' + sea + 'anomSIG_sicncont_' + sicsea + \
                '_' + ptype + '_gissnsidcsim1979-89_2002-12_horiz.png')



## same as above but NO SIC CONTOURS
printtofile=False
# # this one is 1 row
fig,axs=plt.subplots(1,2)
fig.set_size_inches((8,4))
fig.subplots_adjust(hspace=.02,wspace=.02)
ax=axs[0]
bm,pc=cplt.kemmap(plotfld,lat,lon,ptype=ptype,axis=ax,cmin=cmin,cmax=cmax,
                  cmap='blue2red_20',lcol='0.3',suppcb=True)
ax.set_title('a. Observed SAT change')

ax=axs[1]
bm,pc=cplt.kemmap(plotregsim,latsim,lonsim,ptype=ptype,axis=ax,cmin=cmin,cmax=cmax,
                  cmap='blue2red_20',lcol='0.3',suppcb=True)
ax.set_title('b. Simulated SAT change')
cbar_ax = fig.add_axes([.25,0.07, 0.5, .04])
cbar=fig.colorbar(pc,cax=cbar_ax, orientation='horizontal',format='%.1f') 

cbar_ax.tick_params(labelsize=15)
cbar.set_ticks(np.arange(-1,1.2,0.2))
#cbar.set_label('$\Delta$ SAT')
cbar.ax.set_ylabel(units,rotation=0,labelpad=20)
cbar.ax.yaxis.set_label_position('right')
cbar_ax.set_xticklabels(('-1.0','','','','','0','','','','','1.0'))

if printtofile:
    fig.savefig(field + '_' + sea + 'anom_' + \
                '_' + ptype + '_gissnsidcsim1979-89_2002-12_horiz.pdf')




################### WRITE NETCDF #####################
if write_temptimeseries:

    testdata=False

    # get the data again to be safe (since script is a mess):
    #gisfile = '/HOME/rkm/work/DATA/GISS/td_giss_tsurf1200_188001_201509_fill_128_64_st_1880011612-2015091600.nc'
    #gisfile = '/HOME/rkm/work/DATA/GISS/td_giss_tsurf1200_188001_201509_128_64_st_1880011612-2015091600.nc'
    # Use the data obtained directy from GISS
    gisfile = '/HOME/rkm/work/DATA/GISS/gistemp1200_ERSSTv4.nc'
    gisst= cnc.getNCvar(gisfile,'tempanomaly',timesel='1979-01-01,2015-07-01',seas=sea) # 'tempanomaly'
    latgis=cnc.getNCvar(gisfile,'lat')
    longis=cnc.getNCvar(gisfile,'lon')
    print 'write_temptimeseries, region= ' + region
    gisreg = cutl.calc_regmean(gisst,latgis,longis,region,model=None)
    gisgm = cutl.calc_regmean(gisst,latgis,longis,'gm',model=None)
    gisnh = cutl.calc_regmean(gisst,latgis,longis,'nh',model=None)


    if testdata:
        gisfileor = '/HOME/rkm/work/DATA/GISS/gistemp1200_ERSST.nc'
        gisfilev4 = '/HOME/rkm/work/DATA/GISS/gistemp1200_ERSSTv4.nc'
        gisor = cnc.getNCvar(gisfileor,'tempanomaly',timesel='1979-01-01,2015-07-01',seas=sea) 
        lator = cnc.getNCvar(gisfileor,'lat')
        lonor = cnc.getNCvar(gisfileor,'lon')
        gisorv4 = cnc.getNCvar(gisfilev4,'tempanomaly',timesel='1979-01-01,2015-07-01',seas=sea) 

        gisnhv3 = cutl.calc_regmean(gisor,lator,lonor,'nh',model=None)
        gisnhv4 = cutl.calc_regmean(gisorv4,lator,lonor,'nh',model=None)

        plt.figure()
        plt.plot(gisnhv3,'r')
        plt.plot(gisnhv4,'k')
        #plt.plot(gisnh,'b')
        plt.legend(('NH v3.2', 'NH v4'))#,'NH model grid'))

        gfnhornc='/HOME/rkm/pyscripts/netcdf/giss_DJF_nh_1979-2013_timeseries.nc'
        gnhor=cnc.getNCvar(gfnhornc,'tempanomaly')
        plt.plot(gnhor,'g')
        plt.legend(('NH v3.2', 'NH v4','NH model grid','NH 1979-2013 processed data'))

        gfnhnc='/HOME/rkm/pyscripts/netcdf/giss_DJF_nh_1979-2014_timeseries.nc'
        gnh=cnc.getNCvar(gfnhnc,'tempanomaly')
        plt.plot(gnh,'y')
        plt.legend(('NH v3.2', 'NH v4','NH model grid','NH 1979-2013 processed data','NH 1979-2014 processed data'))

        gfnhnc2='/HOME/rkm/pyscripts/giss_DJF_nh_1979-2014_timeseries.nc'
        gnh2=cnc.getNCvar(gfnhnc2,'tempanomaly')
        plt.plot(gnh2,'cyan',linestyle='--')
        plt.legend(('NH v3.2', 'NH v4','NH model grid','NH 1979-2013 processed data',
                    'NH 1979-2014 processed data','NH 1979-2014 processed data fix?'))

    from netCDF4 import Dataset

    if euranom:
        outfile='giss_' + sea + '_' + region + '_1979-2014_timeseries.nc'
    elif nhanom:
        outfile='giss_' + sea + '_nh_1979-2014_timeseries.nc'
    elif glob:
        outfile='giss_' + sea + '_gm_1979-2014_timeseries.nc'
    else:
        outfile='giss_' + sea + '_' + region + '-nh_1979-2014_timeseries.nc'

    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)

    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    outfld = outnc.createVariable('tempanomaly','f4',('time',),fill_value=1.0e38)       

    # add attributes to variables
    outfld.units = 'K'
    if euranom:
        outfld.long_name = 'Surface temperature anomaly from 1951-1980, regional avg: ' + region + ', seasonal avg: ' + sea
    elif nhanom:
        outfld.long_name = 'Surface temperature anomaly from 1951-1980, regional avg: NH, seasonal avg: ' + sea
    elif glob: 
        outfld.long_name = 'Surface temperature anomaly from 1951-1980, regional avg: global mean, seasonal avg: ' + sea
    else: # eurasia - NH
        outfld.long_name = 'Surface temperature anomaly regional avg: ' + region + '-NH avg, seasonal avg: ' + sea

    # don't add scale factor back. already taken into account

    outtimes.long_name = 'time'
    outtimes.units = 'days since 1800-01-01 00:00:00'
    outtimes.calendar = '365_day'
    
    # global attributes
    import time

    if euranom:
        #outnc.title = 'original file: gistemp1200_ERSST.nc. Regional avg: ' + region + ', Seasonal avg: ' + sea
        outnc.title = 'original file: gistemp1200_ERSSTv4.nc. Regional avg: ' + region + ', Seasonal avg: ' + sea
    elif nhanom:
        #outnc.title = 'original file: gistemp1200_ERSST.nc. Regional avg: NH, Seasonal avg: ' + sea
        outnc.title = 'original file: gistemp1200_ERSSTv4.nc. Regional avg: NH, Seasonal avg: ' + sea
    elif glob:
        #outnc.title = 'original file: gistemp1200_ERSST.nc. Regional avg: global mean, Seasonal avg: ' + sea
        outnc.title = 'original file: gistemp1200_ERSSTv4.nc. Regional avg: global mean, Seasonal avg: ' + sea
    else:
        outnc.title = 'original file: gistemp1200_ERSST.nc. Regional avg: ' + region + '-NH avg, Seasonal avg: ' + sea
        outnc.title = 'original file: gistemp1200_ERSSTv4.nc. Regional avg: ' + region + '-NH avg, Seasonal avg: ' + sea

    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    #outtimes[:] = gistime[11:-12:12] # get all ~Dec except last one
    #print 'gistime size ' + str(gistime[11:-12:12].shape)
    outtimes[:] = gistime[11::12] # get all ~Dec 
    print 'gistime size ' + str(gistime[11::12].shape)
    #outlats[:] = bclat
    #outlons[:] = bclon
    if euranom:
        outfld[:] = gisreg
    elif nhanom:
        outfld[:] = gisnh
    elif glob:
        outfld[:] = gisgm
    else:
        outfld[:] = gisreg-gisnh

    print 'gisreg size ' + str(gisreg.shape)

    outnc.close()


if write_circtimeseries:



    graveraint= 9.80665 # m/s2 (different from Canadian models)

    #erafile = '/HOME/rkm/work/DATA/ERAINT/td_era_int_197901_201412_gp_128_64_phi500_1979011612-2014121612.nc'
    erafile = '/HOME/rkm/work/DATA/ERAINT/td_era_int_197901_201507_gp_128_64_phi500_1979011612-2015071612.nc'
    eraz500= cnc.getNCvar(erafile,'PHI',timesel='1979-01-01,2015-07-01',seas=sea)/graveraint
    latera=cnc.getNCvar(erafile,'lat')
    lonera=cnc.getNCvar(erafile,'lon')
    erareg = cutl.calc_regmean(eraz500,latera,lonera,sicregion,model=None)
    eraregeur = cutl.calc_regmean(eraz500,latera,lonera,'eurasiamori',model=None)

    eratime = cnc.getNCvar(erafile,'time',timesel=timesel)

    #eraz5002= cnc.getNCvar(erafile2,'PHI',timesel='1979-01-01,2015-07-01',seas=sea)/graveraint
    ##eraz5002 = cutl.seasonalize_monthlyts(eraz5002,season=sea)# doesn't seem to get last season...check func@@
    #latera2=cnc.getNCvar(erafile2,'lat')
    #lonera2=cnc.getNCvar(erafile2,'lon')
    #erareg2 = cutl.calc_regmean(eraz5002,latera2,lonera2,sicregion,model=None)


    xx=np.arange(0,len(erareg))
    #xx2=np.arange(0,len(erareg2))

    window=5
    xxrun = xx[window/2:-window/2+1]
    #xxrun2 = xx2[window/2:-window/2+1]
    fig,ax=plt.subplots(1,1)
    ax.plot(xxrun, runmean(erareg,window=window),'0.1',linewidth=2)
    #ax.plot(xxrun2, runmean(erareg2,window=window),'0.1',linewidth=2)
    ax.plot(xx, (erareg),'0.1',linewidth=1,alpha=0.5)
    #ax.plot(xx2, (erareg2),'g',linewidth=1,alpha=0.5)

    #ax.legend(('GM','NH','Eurasia'),loc='upper left',frameon=False)
    #ax.set_title(sea + ' Normalized SAT ($\sigma$) + ' + str(window) + ' yr smooth')
    ax.set_xticks(np.arange(0,len(erareg),5)) # ticks for each decade, in July
    #ax.minorticks_on()
    ax.set_xticklabels(np.arange(styr,2015,5))
    #if printtofile:
    #    fig.savefig(field + 'gm_' + sea + 'anom_' + str(styr) + '-2014_run' + str(window)+'timeseries.pdf')

    fig,ax=plt.subplots(1,1)
    ax.plot(xxrun, runmean(erareg-eraregeur,window=window),'0.1',linewidth=2)
    ax.plot(xx, (erareg-eraregeur),'0.1',linewidth=1,alpha=0.5)
    ax.legend(('BKS-EUR',),loc='upper left',frameon=False)
    ax.set_title(sea + ' BKS Z500-EUR Z500 circ index (m) + ' + str(window) + ' yr smooth')
    ax.set_xticks(np.arange(0,len(erareg),5)) # ticks for each decade, in July
    #ax.minorticks_on()
    ax.set_xticklabels(np.arange(styr,2015,5))


    fig,ax=plt.subplots(1,1)
    ax.plot(xxrun, runmean( (erareg-erareg.mean())/np.std(erareg),window=window),'0.5',linewidth=2)
    ax.plot(xxrun, runmean( ((erareg-eraregeur)-(erareg-eraregeur).mean())/np.std(erareg-eraregeur),
                           window=window),'0.1',linewidth=2)
    ax.plot(xxrun, runmean( (erareg-erareg.mean())/np.std(erareg)-(eraregeur-eraregeur.mean())/np.std(eraregeur),
                            window=window),'r',linewidth=2)
    ax.plot(xx, (erareg-erareg.mean())/np.std(erareg),'0.5',linewidth=1,alpha=0.5)
    ax.plot(xx, ((erareg-eraregeur)-(erareg-eraregeur).mean())/np.std(erareg-eraregeur),
            '0.1',linewidth=1,alpha=0.5)

    ax.plot(xx, (erareg-erareg.mean())/np.std(erareg)-(eraregeur-eraregeur.mean())/np.std(eraregeur),
            'r',linewidth=1,alpha=0.5)

    ax.legend(('BKS','BKS-EUR','BKS standardized - EUR standardized'),loc='upper left',frameon=False)
    ax.set_title(sea + ' Normalized Z500 index ($\sigma$) + ' + str(window) + ' yr smooth')
    ax.set_xticks(np.arange(0,len(erareg),5)) # ticks for each decade, in July
    #ax.minorticks_on()
    ax.set_xticklabels(np.arange(styr,2015,5))



    from netCDF4 import Dataset

    if euranom:
        outfile='eraz500_' + sea + '_' + region + '_1979-2014_timeseries.nc'
    elif bksanom:
        outfile='eraz500_' + sea + '_' + sicregion + '_1979-2014_timeseries.nc'
    #elif glob:
    #    outfile='eraz500_' + sea + '_gm_1979-2014_timeseries.nc'
    else:
        #outfile='eraz500_' + sea + '_' + sicregion + '-' + region + '_1979-2014_timeseries.nc'
        outfile='eraz500_' + sea + '_stdized' + sicregion + '-' + region + '_1979-2014_timeseries.nc'

    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)


    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    outfld = outnc.createVariable('z500','f4',('time',),fill_value=1.0e38)       

    # add attributes to variables
    outfld.units = 'm'
    if euranom:
        outfld.long_name = 'Z500 regional avg: ' + region + ', seasonal avg: ' + sea
    elif bksanom:
        outfld.long_name = 'Z500 regional avg: ' + sicregion + ', seasonal avg: ' + sea
    #elif glob: 
    #    outfld.long_name = 'Z500 regional avg: global mean, seasonal avg: ' + sea
    else:
        #outfld.long_name = 'Z500 regional avg: ' + sicregion + '-' + region + ' avg, seasonal avg: ' + sea
        outfld.long_name = 'standardized Z500 regional avg: ' + sicregion + '-' + region + ' avg, seasonal avg: ' + sea


    outtimes.long_name = 'time'
    outtimes.units = 'days since 1850-1-1'
    outtimes.calendar = '365_day'

    
    # global attributes
    import time

    if euranom:
        outnc.title = 'original file: td_era_int_197901_201507_gp_128_64_phi500_1979011612-2015071612.nc. Regional avg: ' + region + ', Seasonal avg: ' + sea
    elif bksanom:
        outnc.title = 'original file: td_era_int_197901_201507_gp_128_64_phi500_1979011612-2015071612.nc. Regional avg: ' + sicregion + ', Seasonal avg: ' + sea
    #elif glob:
    #    outnc.title = 'original file: td_era_int_197901_201507_gp_128_64_phi500_1979011612-2015071612.nc. Regional avg: global mean, Seasonal avg: ' + sea
    else:
        outnc.title = 'original file: td_era_int_197901_201507_gp_128_64_phi500_1979011612-2015071612.nc. Regional avg: ' + sicregion + '-' + region + ' avg, Seasonal avg: ' + sea

    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    #outtimes[:] = eratime[11:-12:12] # get all ~Dec except last one
    outtimes[:] = eratime[11::12] # get all ~Dec 
    print 'eratime size ' + str(eratime[11::12].shape) 
    #outlats[:] = bclat
    #outlons[:] = bclon
    if euranom:
        outfld[:] = eraregeur
    elif bksanom:
        outfld[:] = erareg
    #elif glob:
    #    outfld[:] = gisgm
    else:
        #outfld[:] = erareg-eraregeur
        # NORMALIZE FIRST
        outfld[:] = (erareg-erareg[:10].mean())/erareg.std() - (eraregeur-eraregeur[:10].mean())/eraregeur.std()

    print 'erareg size ' + str(erareg.shape)

    outnc.close()


if write_temptimeseriesera:


    erafile = '/HOME/rkm/work/DATA/ERAINT/td_era_int_197901_201507_gp_128_64_st_1979011612-2015071612.nc'
    erast= cnc.getNCvar(erafile,'ST',timesel='1979-01-01,2015-07-01',seas=sea)
    latera=cnc.getNCvar(erafile,'lat')
    lonera=cnc.getNCvar(erafile,'lon')
    erareg = cutl.calc_regmean(erast,latera,lonera,region,model=None)
    eragm = cutl.calc_regmean(erast,latera,lonera,'gm',model=None)
    eranh = cutl.calc_regmean(erast,latera,lonera,'nh',model=None)

    eratime = cnc.getNCvar(erafile,'time',timesel=timesel)

    xx=np.arange(0,len(erareg))

    window=5
    xxrun = xx[window/2:-window/2+1]
    fig,ax=plt.subplots(1,1)
    ax.plot(xxrun, runmean(erareg,window=window),'0.1',linewidth=2)
    ax.plot(xx, (erareg),'0.1',linewidth=1,alpha=0.5)

    #ax.legend(('GM','NH','Eurasia'),loc='upper left',frameon=False)
    #ax.set_title(sea + ' Normalized SAT ($\sigma$) + ' + str(window) + ' yr smooth')
    ax.set_xticks(np.arange(0,len(erareg),5)) 
    ax.minorticks_on()
    ax.set_xticklabels(np.arange(styr,2015,5))
    #if printtofile:
    #    fig.savefig(field + 'gm_' + sea + 'anom_' + str(styr) + '-2014_run' + str(window)+'timeseries.pdf')



    from netCDF4 import Dataset

    if euranom:
        outfile='erast_' + sea + '_' + region + '_1979-2014_timeseries.nc'
    elif nhanom:
        outfile='erast_' + sea + '_NH_1979-2014_timeseries.nc'
    elif glob:
        outfile='erast_' + sea + '_gm_1979-2014_timeseries.nc'
    else:
        outfile='erast_' + sea + '_' + region + '-NH_1979-2014_timeseries.nc'

    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)


    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    outfld = outnc.createVariable('ST','f4',('time',),fill_value=1.0e38)       

    # add attributes to variables
    outfld.units = 'K'
    if euranom:
        outfld.long_name = 'ST regional avg: ' + region + ', seasonal avg: ' + sea
    elif nhanom:
        outfld.long_name = 'ST regional avg: NH, seasonal avg: ' + sea
    elif glob: 
        outfld.long_name = 'ST regional avg: global mean, seasonal avg: ' + sea
    else: # eurasia - NH
        outfld.long_name = 'ST regional avg: ' + region + '-NH avg, seasonal avg: ' + sea


    outtimes.long_name = 'time'
    outtimes.units = 'days since 1850-1-1'
    outtimes.calendar = '365_day'

    
    # global attributes
    import time

    if euranom:
        outnc.title = 'original file: td_era_int_197901_201507_gp_128_64_st_1979011612-2015071612.nc. Regional avg: ' + region + ', Seasonal avg: ' + sea
    elif nhanom:
        outnc.title = 'original file: td_era_int_197901_201507_gp_128_64_st_1979011612-2015071612.nc. Regional avg: NH, Seasonal avg: ' + sea
    elif glob:
        outnc.title = 'original file: td_era_int_197901_201507_gp_128_64_st_1979011612-2015071612.nc. Regional avg: global mean, Seasonal avg: ' + sea
    else:
        outnc.title = 'original file: td_era_int_197901_201507_gp_128_64_st_1979011612-2015071612.nc. Regional avg: ' + region + '-NH avg, Seasonal avg: ' + sea

    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    #outtimes[:] = eratime[11:-12:12] # get all ~Dec except last one
    outtimes[:] = eratime[11::12] # get all ~Dec 
    print 'eratime size ' + str(eratime[11::12].shape) 
    #outlats[:] = bclat
    #outlons[:] = bclon
    if euranom:
        outfld[:] = erareg
    elif nhanom:
        outfld[:] = eranh
    elif glob:
        outfld[:] = eragm
    else:
        outfld[:] = erareg-eranh

    print 'erareg size ' + str(erareg.shape)

    outnc.close()




if write_siatimeseries:
    from netCDF4 import Dataset

    outfile='nsidcbt_' + sicsea + '_' + sicregion + '_1979-2013_timeseries.nc'
    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)

    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    outfld = outnc.createVariable('sia','f4',('time',),fill_value=1.0e38)

    # add attributes to variables
    outfld.units = 'm2'
    outfld.long_name = 'Sea ice area, regional avg: ' + sicregion + ', seasonal avg: ' + sicsea

    outtimes.long_name = 'time'
    outtimes.units = 'days since 1800-01-01 00:00:00'
    outtimes.calendar = '365_day'

    # global attributes
    import time

    #outnc.title = 'original file: nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc. Regional avg: ' + sicregion + ', Seasonal avg: ' + sicsea
    outnc.title = 'original file: td_bootstrap_197811_latest_128_64_sicn_1978111600-2013121612.nc Regional avg: ' + sicregion + ', Seasonal avg: ' + sicsea
    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    outtimes[:] = gistime[11:-12:12] # get all ~Dec except last 1 for nsidc (?)
    outfld[:] = nsidcreg

    outnc.close()
   
if write_tempmap:

    from netCDF4 import Dataset

    outfile='giss_' + sea + '_2002-12_minus_1979-89_map.nc'
    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)
    outlat = outnc.createDimension('lat',len(lat))
    outlon = outnc.createDimension('lon',len(lon))

    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    outlats = outnc.createVariable('lat','d',('lat',))
    outlons = outnc.createVariable('lon','d',('lon',))
    outfld = outnc.createVariable('tempanomaly','f4',('time','lat','lon',),fill_value=1.0e38)

    # add attributes to variables
    outfld.units = 'K'
    outfld.long_name = 'Surface temperature anomaly from 1951-1980, seasonal avg: ' + sea + ', 2002-12 minus 1979-89'
    # don't add scale factor back. already taken into account

    outtimes.long_name = 'time'
    outtimes.units = 'days since 1800-01-01 00:00:00'
    outtimes.calendar = '365_day'

    outlats.units = 'degrees_north'
    outlats.long_name = 'Latitude'
    outlats.standard_name = 'latitude'

    outlons.units = 'degrees_east'
    outlons.long_name = 'Longitude'
    outlons.standard_name = 'longitude'
    
    # global attributes
    import time

    outnc.title = 'original file: gistemp1200_ERSST.nc. , Seasonal avg: ' + sea + ', 2002-12 minus 1979-89'

    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    outtimes[:] = gistime[11] # one time (filler) for map
    outlats[:] = lat
    outlons[:] = lon
    outfld[:] = np.expand_dims(gismap,axis=0)
    #outfld[:] = gismap
    outnc.close()


if write_simsicmap:
    # mean sea ice conc BC anom
    
    fnamec,fnamep=con.build_filepathpair('ENSE','sicn')
    lat=cnc.getNCvar(fnamec,'lat')
    lon=cnc.getNCvar(fnamep,'lon')
    ctl=cnc.getNCvar(fnamec,'SICN',seas=sea).mean(axis=0)
    pt=cnc.getNCvar(fnamep,'SICN',seas=sea).mean(axis=0)
    sicdiff=pt-ctl

    # test
    plt.figure()
    cplt.kemmap(pt-ctl,lat,lon,ptype='nh',cmap='red2blue_w20',cmin=-.2,cmax=.2)


    from netCDF4 import Dataset

    outfile='ENSEsim_SICN_' + sea + '_2002-12_minus_1979-89_map.nc'
    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)
    outlat = outnc.createDimension('lat',len(lat))
    outlon = outnc.createDimension('lon',len(lon))

    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    outlats = outnc.createVariable('lat','d',('lat',))
    outlons = outnc.createVariable('lon','d',('lon',))
    outfld = outnc.createVariable('SICN','f4',('time','lat','lon',),fill_value=1.0e38)

    # add attributes to variables
    outfld.units = 'frac'
    outfld.long_name = 'Sea ice concentration anomaly, seasonal avg: ' + sea + ', 2002-12 minus 1979-89'

    outtimes.long_name = 'time'
    outtimes.units = 'days since 1800-01-01 00:00:00'
    outtimes.calendar = '365_day'

    outlats.units = 'degrees_north'
    outlats.long_name = 'Latitude'
    outlats.standard_name = 'latitude'

    outlons.units = 'degrees_east'
    outlons.long_name = 'Longitude'
    outlons.standard_name = 'longitude'
    
    # global attributes
    import time

    outnc.title = 'original files: kemctl1ense_sicn_001-121_ts.nc, kem1pert2ense_sicn_001-121_ts.nc. , Seasonal avg: ' + sea + ', 2002-12 minus 1979-89'

    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    outtimes[:] = gistime[11] # one time (filler) for map
    outlats[:] = lat
    outlons[:] = lon
    outfld[:] = np.expand_dims(sicdiff,axis=0)
    outnc.close()


if write_nsidcsicmap:
    # mean sea ice conc BC anom
    
    fnamec,fnamep=con.build_filepathpair('NSIDC','sicn')
    lat=cnc.getNCvar(fnamec,'lat')
    lon=cnc.getNCvar(fnamep,'lon')
    ctl=cnc.getNCvar(fnamec,'SICN',seas=sea).mean(axis=0)
    pt=cnc.getNCvar(fnamep,'SICN',seas=sea).mean(axis=0)
    sicdiff=pt-ctl

    # test
    plt.figure()
    cplt.kemmap(pt-ctl,lat,lon,ptype='nh',cmap='red2blue_w20',cmin=-.2,cmax=.2)


    from netCDF4 import Dataset

    outfile='NSIDCsim_SICN_' + sea + '_2002-11_minus_1979-89_map.nc'
    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)
    outlat = outnc.createDimension('lat',len(lat))
    outlon = outnc.createDimension('lon',len(lon))

    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    outlats = outnc.createVariable('lat','d',('lat',))
    outlons = outnc.createVariable('lon','d',('lon',))
    outfld = outnc.createVariable('SICN','f4',('time','lat','lon',),fill_value=1.0e38)

    # add attributes to variables
    outfld.units = 'frac'
    outfld.long_name = 'Sea ice concentration anomaly, seasonal avg: ' + sea + ', 2002-11 minus 1979-89'

    outtimes.long_name = 'time'
    outtimes.units = 'days since 1800-01-01 00:00:00'
    outtimes.calendar = '365_day'

    outlats.units = 'degrees_north'
    outlats.long_name = 'Latitude'
    outlats.standard_name = 'latitude'

    outlons.units = 'degrees_east'
    outlons.long_name = 'Longitude'
    outlons.standard_name = 'longitude'
    
    # global attributes
    import time

    outnc.title = 'original files: kemnsidcctl_sicn_001-121_ts.nc, kemnsidcpert_sicn_001-121_ts.nc. , Seasonal avg: ' + sea + ', 2002-11 minus 1979-89'

    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    outtimes[:] = gistime[11] # one time (filler) for map
    outlats[:] = lat
    outlons[:] = lon
    outfld[:] = np.expand_dims(sicdiff,axis=0)
    outnc.close()



if write_simsatmap:
    # SAT anom for cold or warm cases in the Average SIC forcing ensemble
    casename='E1'
    fnamec,fnamep=con.build_filepathpair(casename,'st')
    lat=cnc.getNCvar(fnamec,'lat')
    lon=cnc.getNCvar(fnamep,'lon')
    ctl=cnc.getNCvar(fnamec,'ST',seas=sea).mean(axis=0)
    pt=cnc.getNCvar(fnamep,'ST',seas=sea).mean(axis=0)
    stdiff=pt-ctl

    # test
    plt.figure()
    cplt.kemmap(pt-ctl,lat,lon,ptype='eabksstere',cmap='blue2red_20',cmin=-1,cmax=1)


    from netCDF4 import Dataset

    outfile=casename + 'sim_ST_' + sea + '_2002-12_minus_1979-89_map.nc'
    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)
    outlat = outnc.createDimension('lat',len(lat))
    outlon = outnc.createDimension('lon',len(lon))

    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    outlats = outnc.createVariable('lat','d',('lat',))
    outlons = outnc.createVariable('lon','d',('lon',))
    outfld = outnc.createVariable('ST','f4',('time','lat','lon',),fill_value=1.0e38)

    # add attributes to variables
    outfld.units = 'deg C'
    outfld.long_name = 'Surface air temperature anomaly, seasonal avg: ' + sea + ', 2002-12 minus 1979-89'

    outtimes.long_name = 'time'
    outtimes.units = 'days since 1800-01-01 00:00:00'
    outtimes.calendar = '365_day'

    outlats.units = 'degrees_north'
    outlats.long_name = 'Latitude'
    outlats.standard_name = 'latitude'

    outlons.units = 'degrees_east'
    outlons.long_name = 'Longitude'
    outlons.standard_name = 'longitude'
    
    # global attributes
    import time

    outnc.title = 'original files: kemctl1' + casename.lower() + '_st_001-121_ts.nc, kem1pert2' + casename.lower() + '_st_001-121_ts.nc. , Seasonal avg: ' + sea + ', 2002-12 minus 1979-89'

    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    outtimes[:] = gistime[11] # one time (filler) for map
    outlats[:] = lat
    outlons[:] = lon
    outfld[:] = np.expand_dims(stdiff,axis=0)
    outnc.close()

if write_simz500map:
    # Z500 anom for cold or warm cases in the Average SIC forcing ensemble
    casename='E1'
    fnamec,fnamep=con.build_filepathpair(casename,'gz50000')
    lat=cnc.getNCvar(fnamec,'lat')
    lon=cnc.getNCvar(fnamep,'lon')
    ctl=cnc.getNCvar(fnamec,'PHI',seas=sea).mean(axis=0)*1/con.get_g()
    pt=cnc.getNCvar(fnamep,'PHI',seas=sea).mean(axis=0)*1/con.get_g()
    gzdiff=pt-ctl

    # test
    plt.figure()
    cplt.kemmap(pt-ctl,lat,lon,ptype='eabksstere',cmap='blue2red_20',cmin=-10,cmax=10)


    from netCDF4 import Dataset

    outfile=casename + 'sim_Z500_' + sea + '_2002-12_minus_1979-89_map.nc'
    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)
    outlat = outnc.createDimension('lat',len(lat))
    outlon = outnc.createDimension('lon',len(lon))

    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    outlats = outnc.createVariable('lat','d',('lat',))
    outlons = outnc.createVariable('lon','d',('lon',))
    outfld = outnc.createVariable('GZ','f4',('time','lat','lon',),fill_value=1.0e38)

    # add attributes to variables
    outfld.units = 'm'
    outfld.long_name = 'Geopotential height at 500 hPa anomaly, seasonal avg: ' + sea + ', 2002-12 minus 1979-89'

    outtimes.long_name = 'time'
    outtimes.units = 'days since 1800-01-01 00:00:00'
    outtimes.calendar = '365_day'

    outlats.units = 'degrees_north'
    outlats.long_name = 'Latitude'
    outlats.standard_name = 'latitude'

    outlons.units = 'degrees_east'
    outlons.long_name = 'Longitude'
    outlons.standard_name = 'longitude'
    
    # global attributes
    import time

    outnc.title = 'original files: kemctl1' + casename.lower() + '_gz50000_001-121_ts.nc, kem1pert2' + casename.lower() + '_gz50000_001-121_ts.nc. , Seasonal avg: ' + sea + ', 2002-12 minus 1979-89'

    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    outtimes[:] = gistime[11] # one time (filler) for map
    outlats[:] = lat
    outlons[:] = lon
    outfld[:] = np.expand_dims(gzdiff,axis=0)
    outnc.close()
