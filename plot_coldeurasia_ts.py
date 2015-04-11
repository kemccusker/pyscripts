import cccmaplots as cplt
import cccmaNC as cnc
import cccmautils as cutl
import constants as con
import cccmacmaps as ccm
import platform as platform

cnc=reload(cnc)
cutl=reload(cutl)

write_temptimeseries=False # write the temperature timeseries to netcdf
euranom=False
nhanom=False # if both False and temptimeseries=True, do eur - nh

write_siatimeseries=True # write the sea ice area timeseries to netcdf
write_tempmap=False
write_simsicmap=False
write_nsidcsicmap=False

printtofile=False
plt.close('all')

styr=1979
timesel=str(styr) + '-01-01,2014-12-31'

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
# base is 1951-1980

fname = obasepath2 + gisfile

fld = cnc.getNCvar(fname,'tempanomaly',timesel=timesel)
lat = cnc.getNCvar(fname,'lat')
lon = cnc.getNCvar(fname,'lon')
gistime = cnc.getNCvar(fname,'time',timesel=timesel)

dates = cnc.get_NCdates(fname)

# ADD Simulated obs == NSIDC
filecnsidc,filepnsidc=con.build_filepathpair('NSIDC','st')

fldreg = cutl.calc_regmean(fld,lat,lon,region=region)
fldgm = cutl.global_mean_areawgted3d(fld,lat,lon)
fldnh = cutl.calc_regmean(fld,lat,lon,region='nh')


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

fldregsim = cutl.calc_regmean(fldpsim-fldcsim,latsim,lonsim,region=region)
plotregsim = fldpsim.mean(axis=0)-fldcsim.mean(axis=0)
regsimm = fldregsim.mean()
regsimstd=fldregsim.std()


fld = cnc.getNCvar(fname,'tempanomaly',timesel=timesel, seas=sea)
fldreg = cutl.calc_regmean(fld,lat,lon,region=region)
gisreg=fldreg # write to file
gisnh = cutl.calc_regmean(fld,lat,lon,region='nh') # seasonal


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


fldsic = cnc.getNCvar(filesic,'SICN',timesel=timesel,seas=sicsea)
fldsicreg = cutl.calc_regtotseaicearea(fldsic,latsic,lonsic,region=sicregion)#isarea=False
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
gmtr = cutl.global_mean_areawgted(slope*len(xx),lat,lon)

# 1D trend:
fldregfull = cutl.calc_regmean(fld,lat,lon,region=region)

mm, bb, rval, pval, std_err = sp.stats.linregress(xx,fldregfull)
print region + ' avg and ' + sea + ' avg trend (' + timesel + ') pval: ' + str(pval)


ttl = '1979-2014 ' + sea + ' trend: gm= ' + '$%.2f$'%(gmtr) + ' per ' + str(len(xx)) + 'yr'


ptype='sq'
fig,ax=plt.subplots(1,1)
cplt.kemmap(slope*len(xx),lat,lon,title=ttl,axis=ax,cmin=-2,cmax=2,
            cmap='blue2red_20',drawgrid=True,type=ptype)
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
fld2 = cnc.getNCvar(fname,'tempanomaly',timesel=timesel2, seas=sea)
fldreg2 = cutl.calc_regmean(fld2,lat,lon,region=region)


xx2=np.arange(0,len(fldreg2))
fldrs2 = fld2.reshape((len(xx2),nlat*nlon))
# 1D trend:
mm, bb, rval, pval, std_err = sp.stats.linregress(xx2,fldreg2)
print region + ' avg and ' + sea + ' avg trend (' + timesel2 + ') pval: ' + str(pval)

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


#  # SIC
sic2=cnc.getNCvar(filesic,'SICN',timesel=timesel2,seas=sicsea)
sic2rm = cutl.calc_regtotseaicearea(sic2,latsic,lonsic,region=sicregion)
xxs2=np.arange(0,len(sic2rm))

mm, bb, rval, pval, std_err = sp.stats.linregress(xxs2,sic2rm)
print sicregion + ' avg and ' + sicsea + ' avg trend (' + timesel2 + ') pval: ' + str(pval)


######## EPOCH differences
cmin=-1; cmax=1 # -1,1 match simulated

print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'

timesel1='1979-01-01,1989-12-31'
timesel2='2002-01-01,2012-12-31'
fld1 = cnc.getNCvar(fname,'tempanomaly',timesel=timesel1, seas=sea)
fld2 = cnc.getNCvar(fname,'tempanomaly',timesel=timesel2, seas=sea)
tstatmap,pvalmap = cutl.ttest_ind(fld1,fld2)

fld1regm= cutl.calc_regmean(fld1,lat,lon,region)
fld2regm = cutl.calc_regmean(fld2,lat,lon,region)
# SIGNIFICANCE b/w two means:
tstat,pval = cutl.ttest_ind(fld1regm,fld2regm)
print 'EPOCH diffs for ' + region + ' ' + sea + ' avg pval: ' + str(pval)

fldregtm = fld1regm.mean(axis=0) - \
         fld2regm.mean(axis=0)
plotfld = fld2.mean(axis=0) - fld1.mean(axis=0)
gismap=plotfld # for write to file

# Now get observed SIC (NSIDC)
#basepath='/HOME/rkm/work/BCs/NSIDC/'
#file1='nsidc_bt_128x64_1978m11_2011m12_sic_1979-1989climo.nc'
#file2='nsidc_bt_128x64_1978m11_2011m12_sic_2002-2011climo.nc'
#filets='nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc'

sic1=cnc.getNCvar(filesic,'SICN',timesel=timesel1,seas=sicsea)
sic2=cnc.getNCvar(filesic,'SICN',timesel=timesel2,seas=sicsea)
sic1regm= cutl.calc_regmean(sic1,latsic,lonsic,sicregion)
sic2regm = cutl.calc_regmean(sic2,latsic,lonsic,sicregion)
# SIGNIFICANCE b/w two means:
tstat,pval = cutl.ttest_ind(sic1regm,sic2regm)
print 'EPOCH diffs for ' + sicregion + ' ' + sicsea + ' avg pval: ' + str(pval)


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
bm,pc=cplt.kemmap(plotfld,lat,lon,type=ptype,axis=ax,cmin=cmin,cmax=cmax,cmap='blue2red_20')
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
cplt.addtsigm(bm,pvalmap,lat,lon,siglevel=0.1)
if printtofile:
    fig.savefig(field + '_' + sea+ 'anomsig_sicncont_' + sicsea + '_' + ptype + '_giss1979-89_2002-12.pdf')


# Do both climo SIE contours instead === this one 
#   is no good as contours are on top of e/o in DJF
contssic=[0.15, 0.15]
fig,ax=plt.subplots(1,1)
bm,pc=cplt.kemmap(plotfld,lat,lon,type=ptype,axis=ax,cmin=cmin,cmax=cmax,cmap='blue2red_20')
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
bm,pc=cplt.kemmap(plotfld,lat,lon,type=ptype,axis=ax,cmin=cmin,cmax=cmax,cmap='blue2red_20')
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')

ax=axs[1]
bm,pc=cplt.kemmap(plotregsim,latsim,lonsim,type=ptype,axis=ax,cmin=cmin,cmax=cmax,cmap='blue2red_20')
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
if printtofile:
    fig.savefig(field + '_' + sea + 'anom_sicncont_' + sicsea + '_' + ptype + '_gissnsidcsim1979-89_2002-12.pdf')


printtofile=True
# # this one is 1 row
fig,axs=plt.subplots(1,2)
fig.set_size_inches((8,4))
fig.subplots_adjust(hspace=.02,wspace=.02)
ax=axs[0]
bm,pc=cplt.kemmap(plotfld,lat,lon,type=ptype,axis=ax,cmin=cmin,cmax=cmax,
                  cmap='blue2red_20',lcol='0.3',suppcb=True)
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
cplt.addtsigm(bm,pvalmap,lat,lon,siglevel=0.1)

ax.set_title('a. Observed SAT change')

ax=axs[1]
bm,pc=cplt.kemmap(plotregsim,latsim,lonsim,type=ptype,axis=ax,cmin=cmin,cmax=cmax,
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
bm,pc=cplt.kemmap(plotfld,lat,lon,type=ptype,axis=ax,cmin=cmin,cmax=cmax,
                  cmap='blue2red_20',lcol='0.3',suppcb=True)
ax.set_title('a. Observed SAT change')

ax=axs[1]
bm,pc=cplt.kemmap(plotregsim,latsim,lonsim,type=ptype,axis=ax,cmin=cmin,cmax=cmax,
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

    from netCDF4 import Dataset

    if euranom:
        outfile='giss_' + sea + '_' + region + '_1979-2013_timeseries.nc'
    elif nhanom:
        outfile='giss_' + sea + '_nh_1979-2013_timeseries.nc'
    else:
        outfile='giss_' + sea + '_' + region + '-nh_1979-2013_timeseries.nc'

    outnc = Dataset(outfile,'w',format='NETCDF3_CLASSIC')

    # create the dimensions
    outtime = outnc.createDimension('time', None)
    #outlat = outnc.createDimension('lat',len(lat))
    #outlon = outnc.createDimension('lon',len(lon))

    # create variables
    outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
    #outlats = outnc.createVariable('lat','d',('lat',))
    #outlons = outnc.createVariable('lon','d',('lon',))
    outfld = outnc.createVariable('tempanomaly','f4',('time',),fill_value=1.0e38)       

    # add attributes to variables
    outfld.units = 'K'
    if euranom:
        outfld.long_name = 'Surface temperature anomaly from 1951-1980, regional avg: ' + region + ', seasonal avg: ' + sea
    elif nhanom:
        outfld.long_name = 'Surface temperature anomaly from 1951-1980, regional avg: NH, seasonal avg: ' + sea
    else: # eurasia - NH
        outfld.long_name = 'Surface temperature anomaly regional avg: ' + region + '-NH avg, seasonal avg: ' + sea

    # don't add scale factor back. already taken into account

    outtimes.long_name = 'time'
    outtimes.units = 'days since 1800-01-01 00:00:00'
    outtimes.calendar = '365_day'

    #outlats.units = 'degrees_north'
    #outlats.long_name = 'Latitude'
    #outlats.standard_name = 'latitude'

    #outlons.units = 'degrees_east'
    #outlons.long_name = 'Longitude'
    #outlons.standard_name = 'longitude'
    
    # global attributes
    import time

    if euranom:
        outnc.title = 'original file: gistemp1200_ERSST.nc. Regional avg: ' + region + ', Seasonal avg: ' + sea
    elif nhanom:
        outnc.title = 'original file: gistemp1200_ERSST.nc. Regional avg: NH, Seasonal avg: ' + sea
    else:
        outnc.title = 'original file: gistemp1200_ERSST.nc. Regional avg: ' + region + '-NH avg, Seasonal avg: ' + sea

    outnc.creation_date = time.ctime(time.time())
    outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

    # set the data to the variables: important to have [:]!
    outtimes[:] = gistime[11:-12:12] # get all ~Dec except last one
    print 'gistime size ' + str(gistime[11:-12:12].shape)
    #outlats[:] = bclat
    #outlons[:] = bclon
    if euranom:
        outfld[:] = gisreg
    elif nhanom:
        outfld[:] = gisnh
    else:
        outfld[:] = gisreg-gisnh

    print 'gisreg size ' + str(gisreg.shape)

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
    cplt.kemmap(pt-ctl,lat,lon,type='nh',cmap='red2blue_w20',cmin=-.2,cmax=.2)


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
    cplt.kemmap(pt-ctl,lat,lon,type='nh',cmap='red2blue_w20',cmin=-.2,cmax=.2)


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
