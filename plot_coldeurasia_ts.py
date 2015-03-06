import cccmaplots as cplt
import cccmaNC as cnc
import cccmautils as cutl
import constants as con

cnc=reload(cnc)
cutl=reload(cutl)

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

filesic='/HOME/rkm/work/BCs/NSIDC/nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc'
#filesic='/Volumes/MyPassport2TB/DATA/OBSERVATIONS/nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc'
latsic=cnc.getNCvar(filesic,'lat')
lonsic=cnc.getNCvar(filesic,'lon')


#basepath = '/Volumes/MyPassport2TB/DATA/OBSERVATIONS/'
#basepath = '/raid/ra40/data/ncs/reanalyses/'
basepath = '/HOME/rkm/work/DATA/GISS/'
file = 'gistemp1200_ERSST.nc'
# base is 1951-1980

fname = basepath + file

fld = cnc.getNCvar(fname,'tempanomaly',timesel=timesel)
lat = cnc.getNCvar(fname,'lat')
lon = cnc.getNCvar(fname,'lon')

dates = cnc.get_NCdates(fname)

# ADD Simulated obs == NSIDC
filecnsidc,filepnsidc=con.build_filepathpair('NSIDC','st')

fldreg = cutl.calc_regmean(fld,lat,lon,region=region)
fldgm = cutl.global_mean_areawgted3d(fld,lat,lon)


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
fldregsim = cutl.calc_regmean(fldpsim-fldcsim,latsim,lonsim,region=region)
plotregsim = fldpsim.mean(axis=0)-fldcsim.mean(axis=0)
regsimm = fldregsim.mean()
regsimstd=fldregsim.std()


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

fldsic = cnc.getNCvar(filesic,'SICN',timesel=timesel,seas=sicsea)
fldsicreg = cutl.calc_regtotseaicearea(fldsic,latsic,lonsic,region=sicregion)#isarea=False

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
                '_' + sicsea + sea + '_' + str(styr) + '-2014_2011_timeseries_onepanel.pdf')





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


######## EPOCH differences
cmin=-1; cmax=1 # -1,1 match simulated

print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'

timesel1='1979-01-01,1989-12-31'
timesel2='2002-01-01,2012-12-31'
fld1 = cnc.getNCvar(fname,'tempanomaly',timesel=timesel1, seas=sea)
fld2 = cnc.getNCvar(fname,'tempanomaly',timesel=timesel2, seas=sea)

fldregtm = cutl.calc_regmean(fld2,lat,lon,region).mean(axis=0) - \
         cutl.calc_regmean(fld1,lat,lon,region).mean(axis=0)
plotfld = fld2.mean(axis=0) - fld1.mean(axis=0)


# Now get observed SIC (NSIDC)
#basepath='/HOME/rkm/work/BCs/NSIDC/'
#file1='nsidc_bt_128x64_1978m11_2011m12_sic_1979-1989climo.nc'
#file2='nsidc_bt_128x64_1978m11_2011m12_sic_2002-2011climo.nc'
#filets='nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc'

sic1=cnc.getNCvar(filesic,'SICN',timesel=timesel1,seas=sicsea)
sic2=cnc.getNCvar(filesic,'SICN',timesel=timesel2,seas=sicsea)

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
if printtofile:
    fig.savefig(field + '_' + sea+ 'anom_sicncont_' + sicsea + '_' + ptype + '_giss1979-89_2002-12.pdf')


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
ax.set_title('a. Observed SAT change')

ax=axs[1]
bm,pc=cplt.kemmap(plotregsim,latsim,lonsim,type=ptype,axis=ax,cmin=cmin,cmax=cmax,
                  cmap='blue2red_20',lcol='0.3',suppcb=True)
bm.contour(lons,lats,plotsic,levels=contssic,
           colors='w',linewidths=2,latlon=True,linestyles='-')
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
    fig.savefig(field + '_' + sea + 'anom_sicncont_' + sicsea + \
                '_' + ptype + '_gissnsidcsim1979-89_2002-12_horiz.pdf')


## same as above but NO SIC CONTOURS
printtofile=True
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
