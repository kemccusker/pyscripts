
import cccmautils as cutl
import pandas as pd
import scipy as sp
import scipy.stats

printtofile=False
plt.close('all')

# @@ Look at Southern Ocean and Southern Jet


basepath='/ra40/data/kem/CanSISE/CanESM2/LE/'
sims=('historical-r1','historical-r2','historical-r3','historical-r4','historical-r5')
simsnat=('historicalNat-r1','historicalNat-r2','historicalNat-r3','historicalNat-r4','historicalNat-r5')
ensnum=10

# CMIP variable names
field='sic'; comp='OImon'; cmin=-3e11; cmax=3e11;cmap='red2blue_w20'
#field='tas'; comp='Amon'; cmin=-1; cmax=1; cmap='blue2red_w20'
reg2=True # if sea ice, include southern hem? if other, include regional mean?
region= 'etroppac'

sea='ANN'
timesel=None
timesel2 = '2002-01-01,2012-12-31'
superii=0
regstr=''


if reg2:
    if field=='sic':
        regstr='sh'
    else:
        regstr=region
    
years=np.arange(1950,2021)
if sea=='DJF':
    years=years[:-1]

simnhdt={}
simshdt={}
simgmdt={}
simrmdt={}

simtrnddt={}
alltrnddt={}
for sim in sims:

    ensnhdt={} # northern hem
    ensshdt={} # southern hem
    ensgmdt={} # global mean
    ensrmdt={} # regional mean

    ensseldt={}
    enstrnddt={}
    for eii in np.arange(1,ensnum+1):

        fname=basepath + sim + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
               sim + '_r' + str(eii) + 'i1p1_195001-202012.nc'

        if superii==0: # just get once
            lat=cnc.getNCvar(fname,'lat')
            lon=cnc.getNCvar(fname,'lon')
            
        fld = cnc.getNCvar(fname,field,timesel=timesel,seas=sea)
        # select a given time period and keep maps
        fldsel = cnc.getNCvar(fname,field,timesel=timesel2,seas=sea)
        nt = fldsel.shape[0]
        nlon=fldsel.shape[2]
        nlat=fldsel.shape[1]
        fldre = fldsel.reshape((nt,nlon*nlat))

        xx=np.arange(0,nt)
        
        #ensseldt[eii] = fldsel

        if field=='sic':
            ensnhdt[eii],ensshdt[eii] = cutl.calc_totseaicearea(fld/100.,lat,lon,isarea=False)
            fldsel=cutl.calc_seaicearea(fldsel,lat,lon)
            fldre=fldsel.reshape((nt,nlon*nlat))
        else:
            ensgmdt[eii] = cutl.global_mean_areawgted3d(fld,lat,lon)
            ensrmdt[eii] = cutl.calc_regmean(fld,lat,lon,region)

    
        #slope[eii], intercept, r_value, p_value, std_err = sp.stats.linregress(xx,dat) # not good for 3d data?
        # this is just the second timesel (ie 2002-2012)
        slope,intercept = np.polyfit(xx,fldre,1) # supposedly can do w/ higher dims?
        
        enstrnddt[eii] = slope #.reshape((nlat,nlon)) # reshape later @@

        # also save all trends into one dictionary (don't differentiate by seed/base run)
        alltrnddt[superii] = slope
        
        superii+=1

    if field=='sic':
        simnhdt[sim] = pd.DataFrame(ensnhdt,index=years)
        simshdt[sim] = pd.DataFrame(ensshdt,index=years)
    else:
        simgmdt[sim] = pd.DataFrame(ensgmdt,index=years)
        simrmdt[sim] = pd.DataFrame(ensrmdt,index=years)

    

    # calc time mean and trend of ensseldt. Then each ens member is a map of trends
    tmpdf = pd.DataFrame(enstrnddt)#,index=years)
    
    #simseldt[sim] = trenddf
    simtrnddt[sim] = tmpdf

alltrnddf = pd.DataFrame(alltrnddt)



# # first, the timeseries

if field=='sic':
    plotdt = simnhdt # NH ice
    plotdt2 = simshdt # SH ice
else:
    plotdt = simgmdt # global mean
    plotdt2 = simrmdt # regional mean


col=0.0
fig,axs=plt.subplots(2,1,sharey=True)
ax=axs[0]
for sim in sims:
    df=plotdt[sim]
    ax.plot(years,df-df.mean(axis=0),color=str(col))
    if reg2: #field=='sic':
        df=plotdt2[sim]
        ax.plot(years,df-df.mean(axis=0),color=str(col),linestyle='--')
    col+=0.2
ax.set_title(field + ' ' + str(sea) + ' ' + regstr + ' (dashed)')

col=0.0
ax=axs[1]
for sim in sims:
    df=plotdt[sim]
    ax.plot(years,df.mean(axis=1)-df.mean(axis=1).mean(axis=0),color=str(col))
    if reg2: #field=='sic':
        df=plotdt2[sim]
        ax.plot(years,df.mean(axis=1)-df.mean(axis=1).mean(axis=0),color=str(col),linestyle='--')
    col+=0.2
ax.set_xlabel('mean of groups')

if printtofile:
    fig.savefig(field + 'LE_remtm' + str(sea) + regstr + '_timeseries.pdf') # time mean removed


# # Maps 
# plot mean trend map for each "group"
for sim in sims:
    trnds = simtrnddt[sim]
    groupm = trnds.mean(axis=1)
    groupm = groupm.reshape((nlat,nlon))

    plt.figure()
    mh = cplt.kemmap(groupm*nt,lat,lon,cmin=cmin,cmax=cmax,title=sim,cmap=cmap)


# mean of all runs
trndmean = alltrnddf.mean(axis=1).reshape((nlat,nlon))
plt.figure()
mh = cplt.kemmap(trndmean*nt,lat,lon,cmin=cmin,cmax=cmax,title='LE mean trend ' + timesel2,cmap=cmap)

# @@ now combine with cmip5 ensemble, and CESM if possible..
    

# Use files that already have SIE or SIA calculated ======= 1/6/2015
# ==================== plot NH sea ice area ======================== 
field='sianh'; comp='OImon'; 
season='SO'
# if season is set to just September (9), I get an error:
# ValueError: Big-endian buffer not supported on little-endian compiler
#   I suspect that the hadisst file is the opposite endian to what my python install is on:
#  http://pandas.pydata.org/pandas-docs/dev/gotchas.html#byte-ordering-issues
# This may not come up for 'ND' average because the cnc module already does some
# averaging and stores data in a new structure. Whereas just reading in a month reads
# directly from the file and doesn't do any averaging (??)
#  Try the fix in the link if want an individual month.

import cccmacmaps as ccm

subyrs=np.arange(1979,1990)
subyrs2=np.arange(2002,2013)
printtofile=True
# ======= CANESM LE (TOT) ===
superii=0
allflddt={}
for sim in sims:

    for eii in np.arange(1,ensnum+1):

        fname=basepath + sim + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
               sim + '_r' + str(eii) + 'i1p1_195001-202012.nc'

        allflddt[superii] = cnc.getNCvar(fname,field,seas=season)

        superii+=1

allflddf=pd.DataFrame(allflddt,index=years)

# ======= CANESM LE (NAT) ===

superii=0
allnatdt={}
for sim in simsnat:

    for eii in np.arange(1,ensnum+1):

        fname=basepath + sim + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
               sim + '_r' + str(eii) + 'i1p1_195001-202012.nc'

        allnatdt[superii] = cnc.getNCvar(fname,field,seas=season)

        superii+=1

allnatdf=pd.DataFrame(allnatdt,index=years)


# =============== OBS:
basedir='/home/rkm/work/BCs/'
subdir='HadISST/'
hadfile=basedir + subdir + 'hadisst1.1_bc_128_64_1870_2013m03_' + field + '_1870010100-2013030100.nc'
hadsel='1950-01-01,2012-12-31'
hadyrs=np.arange(1950,2013)
hadfld=cnc.getNCvar(hadfile,field,timesel=hadsel,seas=season)
haddf = pd.Series(hadfld,index=hadyrs)
pasthad = haddf.loc[subyrs]
preshad = haddf.loc[subyrs2]
(hadtstat,hadpv) = cutl.ttest_ind(preshad,pasthad)
diffhad=preshad.mean(axis=0) - pasthad.mean(axis=0)

subdir='NSIDC/'
nsidcfile = basedir + subdir + 'nsidc_bt_128x64_1978m11_2011m12_' + field + '_1978111600-2011121612.nc'
nsidcsel='1979-01-01,2011-12-31'
nsidcyrs=np.arange(1979,2012)
nsidcfld=cnc.getNCvar(nsidcfile,field,timesel=nsidcsel,seas=season)
nsidcdf=pd.Series(nsidcfld,index=nsidcyrs)
pastnsidc = nsidcdf.loc[subyrs]
presnsidc = nsidcdf.loc[subyrs2[:-1]]
(nsidctstat,nsidcpv) = cutl.ttest_ind(presnsidc,pastnsidc)
diffnsidc=presnsidc.mean(axis=0) - pastnsidc.mean(axis=0)


# =============== original five:
casename='historicalrcp45'
ensnum=5
basedir='/home/rkm/work/DATA/CanESM2/' + casename
origsel='1950-01-01,2012-12-31'
origyrs = np.arange(1950,2013)

origdt={}
for eii in np.arange(1,ensnum+1):

    fname=basedir + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
           casename + '_r' + str(eii) + 'i1p1_185001-201212.nc'

    origdt[eii] = cnc.getNCvar(fname,field,timesel=origsel,seas=season)

origdf=pd.DataFrame(origdt,index=origyrs)
pastorig = origdf.loc[subyrs]
presorig = origdf.loc[subyrs2]
(origtstat,origpv) = cutl.ttest_ind(presorig,pastorig,axis=0)
difforig = presorig.mean(axis=0) - pastorig.mean(axis=0)

# =============== original five: NAT
casename='historicalNat'
basedir='/home/rkm/work/DATA/CanESM2/' + casename
orignatdt={}
for eii in np.arange(1,ensnum+1):

    fname=basedir + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
           casename + '_r' + str(eii) + 'i1p1_185001-201212.nc'

    orignatdt[eii] = cnc.getNCvar(fname,field,timesel=origsel,seas=season)

orignatdf=pd.DataFrame(orignatdt,index=origyrs)
pastnatorig = orignatdf.loc[subyrs]
presnatorig = orignatdf.loc[subyrs2]
(natorigtstat,natorigpv) = cutl.ttest_ind(presnatorig,pastnatorig,axis=0)
diffnatorig = presnatorig.mean(axis=0) - pastnatorig.mean(axis=0)


# TIMESERIES
firebrick=ccm.get_linecolor('firebrick')

fig,axs=plt.subplots(1,1)
axs.plot(years,allnatdf,color='0.5',alpha=.5)
axs.plot(years,allflddf,color='r',alpha=.3)
axs.plot(origyrs,orignatdf,color='k',linewidth=2)
axs.plot(origyrs,origdf,color='brown',linewidth=2)
axs.plot(hadyrs,haddf,color='b',linewidth=2)
axs.plot(nsidcyrs,nsidcdf,color='g',linewidth=2)
axs.set_title(field + ' ' + str(season))



# NOW DO DIFF PDFS
# this will work too: allflddf.query('index>=1979 & index<1990')
past=allflddf.loc[subyrs]
pres=allflddf.loc[subyrs2]
pastnat=allnatdf.loc[subyrs]
presnat=allnatdf.loc[subyrs2]

alldiff =pres.mean(axis=0)- past.mean(axis=0)
allnatdiff =presnat.mean(axis=0)- pastnat.mean(axis=0)

fig,ax = plt.subplots(1,1)
allnatdiff.hist(color='0.5',alpha=0.5)
alldiff.hist(color='r',alpha=0.5)
for eii in range(1,6):
    ax.axvline(difforig[eii],color=ccm.get_linecolor('firebrick'),linewidth=3)
    ax.axvline(diffnatorig[eii],color='k',linewidth=3)
ax.axvline(diffhad,color='b',linewidth=3)
ax.axvline(diffnsidc,color='g',linewidth=3)
ax.set_title(field + ' anom ' + str(season))


from scipy.stats import norm

# ======== PDFS
anormed=norm.fit(alldiff)
amean=anormed[0]
asd=anormed[1]
#Generate X points
axlims = [-4*asd+amean, 4*asd+amean] # large limits
axx = np.linspace(axlims[0],axlims[1],500)
#Get Y points via Normal PDF with fitted parameters
apdf_fitted = norm.pdf(axx,loc=amean,scale=asd)

# NAT
nnormed=norm.fit(allnatdiff)
nmean=nnormed[0]
nsd=nnormed[1]
#Generate X points
nxlims = [-4*nsd+nmean, 4*nsd+nmean] # large limits
nxx = np.linspace(nxlims[0],nxlims[1],500)
#Get Y points via Normal PDF with fitted parameters
npdf_fitted = norm.pdf(nxx,loc=nmean,scale=nsd)

# normed means: integral of the histogram will sum to 1 @@@ don't think it's right
fig,ax=plt.subplots(1,1)
allnatdiff.hist(normed=True,color='0.5',alpha=0.5)#,histtype='stepfilled')
alldiff.hist(normed=True,color=firebrick,alpha=0.5)#,histtype='stepfilled')
#plt.title(field + ' anom ' + str(season))
ax.plot(axx,apdf_fitted,color=firebrick,linewidth=2)
ax.plot(nxx,npdf_fitted,color='k',linewidth=2)
axylims=ax.get_ylim()
ax.axvline(color='k',linestyle='--')
for eii in range(1,6):
    ax.plot(difforig[eii],axylims[1],marker='|',markersize=20,color=firebrick) 
    ax.plot(diffnatorig[eii],axylims[1],marker='|',markersize=20,color='k')
    ax.axvline(x=difforig[eii],ymin=.97,ymax=1,color=firebrick,linewidth=2)
    ax.axvline(x=diffnatorig[eii],ymin=.97,ymax=1,color='k',linewidth=2)
hh=ax.axvline(diffhad,color='b',linewidth=3)
nn=ax.axvline(diffnsidc,color='g',linewidth=3)
ax.legend((hh,nn),('HadISST','NSIDC'),loc='upper left')
ax.set_ylabel('Density')
ax.set_xlabel('$\Delta$ Sea Ice Area (m$^2$)')
ax.set_xlabel('$\Delta$ Sea Ice Area (millions of km$^2$)')
xt=ax.get_xticks()
ax.set_xticklabels(xt/np.float(1e12))
ax.grid('off')
if printtofile:
    fig.savefig(field + 'diff_PDFHIST_CanESMLE_TOTNAT_' + str(season) + '.pdf')


fig,ax=plt.subplots(1,1)
allnatdiff.hist(normed=True,color='0.5',alpha=0.5)#,histtype='stepfilled')
alldiff.hist(normed=True,color=firebrick,alpha=0.5)#,histtype='stepfilled')
#plt.title(field + ' anom ' + str(season))
ax.plot(axx,apdf_fitted,color=firebrick,linewidth=2)
ax.plot(nxx,npdf_fitted,color='k',linewidth=2)
axylims=ax.get_ylim()
ax.set_ylim((axylims[0],axylims[1]+(axylims[1]*.1)))
axylims=ax.get_ylim()
ax.axvline(color='k',linestyle='--')
for eii in range(1,6):
    #ax.plot(difforig[eii],axylims[1],marker='|',markersize=20,color=firebrick) 
    #ax.plot(diffnatorig[eii],axylims[1],marker='|',markersize=20,color='k')
    ax.axvline(x=difforig[eii],ymin=.97,ymax=1,color=firebrick,linewidth=2)
    ax.axvline(x=diffnatorig[eii],ymin=.97,ymax=1,color='k',linewidth=2)
hh=ax.axvline(diffhad,color='b',linewidth=3)
nn=ax.axvline(diffnsidc,color='g',linewidth=3)
#ax.legend((hh,nn),('HadISST','NSIDC'),loc='upper left')
ax.annotate('NSIDC',xy=(diffnsidc-(diffnsidc*.05),axylims[1]-(axylims[1]*.07)),xycoords='data')#,rotation=45)
ax.annotate('HadISST',xy=(diffhad-(diffhad*.05),axylims[1]-(axylims[1]*.07)),xycoords='data')#,rotation=45)
ax.annotate('Historical',xy=(.36,.57),xytext=(.2,.7), xycoords='figure fraction',
            arrowprops=dict(facecolor=firebrick, edgecolor='None',shrink=0.05))
ax.annotate('HistoricalNat',xy=(.7,.5),xytext=(.75,.6), xycoords='figure fraction',
            arrowprops=dict(facecolor='k', edgecolor='None',shrink=0.05))
ax.annotate('Sea Ice Loss', xy=(0,1.02), xytext=(.53,1.02),xycoords='axes fraction',
            arrowprops=dict(facecolor='k', edgecolor='None',shrink=0.05),verticalalignment='center')
ax.annotate('Sea Ice Gain', xy=(1,1.02), xytext=(.74,1.02),xycoords='axes fraction',
            arrowprops=dict(facecolor='k', edgecolor='None',shrink=0.05),verticalalignment='center')

ax.set_ylabel('Density')
ax.set_xlabel('$\Delta$ Sea Ice Area (m$^2$)')
ax.set_xlabel('$\Delta$ Sea Ice Area (millions of km$^2$)')
xt=ax.get_xticks()
ax.set_xticklabels(xt/np.float(1e12))
ax.grid('off')
if printtofile:
    fig.savefig(field + 'diff_PDFHIST_CanESMLE_TOTNAT_' + str(season) + '_annotate.pdf')

import matplotlib.patches as mpatches
import matplotlib.lines as mlines
#http://matplotlib.org/users/legend_guide.html#proxy-legend-handles

# # # PAPER ############
ms = 7 # markersize
deffs='none' # default fillstyle
mew=1.5

fig,ax=plt.subplots(1,1)
allnatdiff.hist(normed=True,color='0.5',alpha=0.5)#,histtype='stepfilled')
alldiff.hist(normed=True,color=firebrick,alpha=0.5)#,histtype='stepfilled')
# use proxy artist for histogram legend entries
#tot = mpatches.Patch(color=firebrick,alpha=0.5)#, label='The red data')
#nat = mpatches.Patch(color='0.5',alpha=0.5)#, label='The red data')
tot=mlines.Line2D([],[],color=firebrick,linewidth=2)
nat=mlines.Line2D([],[],color='k',linewidth=2)

ax.plot(axx,apdf_fitted,color=firebrick,linewidth=2)
ax.plot(nxx,npdf_fitted,color='k',linewidth=2)
axylims=ax.get_ylim()
ax.set_ylim((axylims[0],axylims[1]+(axylims[1]*.1)))
axylims=ax.get_ylim()
ax.axvline(color='k',linestyle='--')
for eii in range(1,6):
    #ax.plot(difforig[eii],axylims[1],marker='|',markersize=20,color=firebrick) 
    #ax.plot(diffnatorig[eii],axylims[1],marker='|',markersize=20,color='k')
    # move these to ON the pdf fit.
    # search for difforig[eii] nearest to axx, then plot apdf_fitted[axx] as marker, filled or not based on sig.
    plotx = difforig[eii]
    idx = cutl.find_nearest(axx,plotx)
    ploty = apdf_fitted[idx]
    if origpv[eii-1]<=siglevel: #significant, fill the marker
        fs='full'
    else:
        fs=deffs # defaultfs
        
    ax.plot(plotx,ploty,marker='o',color=firebrick,mec=firebrick,fillstyle=fs,mew=mew,markersize=ms)
    #ax.axvline(x=difforig[eii],ymin=.97,ymax=1,color=firebrick,linewidth=2)

    plotx = diffnatorig[eii]
    idx = cutl.find_nearest(nxx,plotx)
    ploty = npdf_fitted[idx]
    if natorigpv[eii-1]<=siglevel: #significant, fill the marker
        fs='full'
    else:
        fs=deffs # defaultfs
    ax.plot(plotx,ploty,marker='o',color='k',mec='k',fillstyle=fs,mew=mew,markersize=ms)
    #ax.axvline(x=diffnatorig[eii],ymin=.97,ymax=1,color='k',linewidth=2)

#nn=ax.axvline(diffnsidc,color='g',linewidth=3)
plotx=diffnsidc
idx=cutl.find_nearest(axx,plotx)
ploty=apdf_fitted[idx]
if nsidcpv<siglevel: #significant, fill the marker
    fs='full'
else:
    fs=deffs # defaultfs
nn=mlines.Line2D([],[],color='g',linestyle='none',marker='o',mec='g',mew=mew,markersize=ms,fillstyle=fs)
ax.plot(plotx,ploty,marker='o',color='g',mec='g',fillstyle=fs,mew=mew,markersize=ms)
#nn = ax.axvline(x=diffnsidc,ymin=.97,ymax=1,color='g',linewidth=3)

plotx=diffhad
idx=cutl.find_nearest(axx,plotx)
ploty=apdf_fitted[idx]
if hadpv<siglevel: #significant, fill the marker
    fs='full'
else:
    fs=deffs # defaultfs
hh=mlines.Line2D([],[],color='b',linestyle='none',marker='o',mec='b',mew=mew,markersize=ms,fillstyle=fs)
ax.plot(plotx,ploty,marker='o',color='b',mec='b',fillstyle=fs,mew=mew,markersize=ms)
#hh = ax.axvline(x=diffhad,ymin=.97,ymax=1,color='b',linewidth=3)

ax.legend((nn,hh,tot,nat),('NSIDC','HadISST','Historical','HistoricalNat'),loc='upper left',frameon=False)
## ax.annotate('NSIDC',xy=(diffnsidc-(diffnsidc*.05),axylims[1]-(axylims[1]*.07)),xycoords='data')#,rotation=45)

## tlaby = 0.94
## txt = ax.annotate('SIA Loss',xy=(0.3, tlaby),  xycoords='axes fraction')#,bbox=txtbox) # just text
## laby=0.98 # 1.03
## ann = ax.annotate('', xy=(0, laby),  xycoords='axes fraction',
##                   xytext=(0.7, laby), textcoords='axes fraction',
##                   arrowprops=dict(arrowstyle='<->',connectionstyle='arc3',facecolor='k',
##                                   edgecolor='k',patchB=txt,patchA=txt),
##                   verticalalignment='center')

## ax.annotate('SIA Gain',xy=(.8,tlaby), xycoords='axes fraction')
## ann = ax.annotate('', xy=(0.72, laby),  xycoords='axes fraction',
##                   xytext=(1, laby), textcoords='axes fraction',
##                   arrowprops=dict(arrowstyle='<->',connectionstyle='arc3',facecolor='k',
##                                   edgecolor='k',patchB=txt,patchA=txt),
##                   verticalalignment='center')

ax.set_ylabel('Density')
ax.set_yticklabels('')
ax.set_xlabel('$\Delta$ Sea Ice Area (millions of km$^2$)')
xt=ax.get_xticks()
ax.set_xticklabels(xt/np.float(1e12))
ax.grid('off')
if season == 'ND':
    ax.set_title('Arctic sea-ice area change (Nov-Dec)')
else:
    ax.set_title('Arctic sea-ice area change (' + season + ')')
if printtofile:
    fig.savefig(field + 'diff_PDFHIST_CanESMLE_TOTNAT_' + str(season) + '_paper2d.pdf')



#printtofile=True
fig,ax=plt.subplots(1,1)
plt.title(field + ' anom ' + str(season))
tot=ax.plot(axx,apdf_fitted,color=ccm.get_linecolor('firebrick'),linewidth=2)
nat=ax.plot(nxx,npdf_fitted,color='k',linewidth=2)
ax.fill_between(axx,apdf_fitted,color='r',alpha=0.5)
ax.fill_between(nxx,npdf_fitted,color='k',alpha=0.5)
#ax.grid()
axylims=ax.get_ylim()

for eii in range(1,6):
    ax.plot(difforig[eii],axylims[1],marker='v',markersize=10,color=firebrick) 
    #ax.plot(diffnatorig[eii],axylims[1],marker='v',markersize=10,color='k') # want original NAT runs?
hh=ax.axvline(diffhad,color='b',linewidth=3)
nn=ax.axvline(diffnsidc,color='g',linewidth=3)
ax.set_ylabel('frequency')
#ax.legend((tot,nat,hh,nn),('TOT LE','NAT LE','HadISST','NSIDC'))
ax.legend((hh,nn),('HadISST','NSIDC'))
if printtofile:
    fig.savefig(field + 'diff_PDF_CanESMLE_TOTNAT_' + str(season) + '.pdf')
