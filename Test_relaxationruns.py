# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>
# @@@@@ FROM NOTEBOOK 9/15/2014

# <headingcell level=3>

# Testing relaxation runs: SICN

# <rawcell>

# Figure titles and text are below the figures.
# Scroll to see monthly sea-ice climatologies as maps, and comparisons between the runs.
# (more runs will be added...)

# <rawcell>

# Goal: check out the test relaxation run data (coupled model)
# 9/4/2014

# <codecell>

import cccmaplots as cplt
import cccmautils as cutl
import cccmaNC as cnc
import constants as con
import cccmacmaps as ccm
import numpy.ma as ma
import pandas as pd
from matplotlib import gridspec
import matplotlib.font_manager as fm

#%matplotlib inline

cutl = reload(cutl)

plt.close('all')
#plt.ion()

printtofile=True
showmaps=False  # show all the monthly map comparisons?


#  set field here ==============
field = 'sic'
# ==============================

model = 'CanESM2'
bp=con.get_basepath()
basepath=bp['basepath'] + model + '/'; subdir=bp['subdir']


timeperc='2921_2930'
timeper2x='2451_2460'

timesel = '2922-01-01,2930-12-31' # skip first year
timesel2x = '2452-01-01,2460-12-31' # skip first year

sims = ('iga','gregory_2xco2','kel11','kel09','kel14','kel15','kel17','kel18','kel20','kel24',
        'kel10','kel12','kel16','kel19','kel21','kel25') # 2:10, nudge to 2xco2. 10:, nudge to preI
timeperdt={'iga':timeperc,
           'gregory_2xco2': timeper2x,
           'kel11': timeperc, 'kel09': timeperc, 'kel14': timeperc, 'kel15': timeperc,
           'kel17': timeperc, 'kel18': timeperc, 'kel20': timeperc,'kel24': timeperc,
           'kel10': timeper2x,'kel12': timeper2x,'kel16': timeper2x,'kel19': timeper2x,
           'kel21': timeper2x, 'kel25': timeper2x}
timeseldt={'iga':timesel,
           'gregory_2xco2': timesel2x,
           'kel11': timesel, 'kel09': timesel, 'kel14': timesel, 'kel15': timesel,
           'kel17': timesel, 'kel18': timesel, 'kel20': timesel, 'kel24': timesel,
           'kel10': timesel2x,'kel12': timesel2x,'kel16': timesel2x,'kel19': timesel2x,
           'kel21': timesel2x,'kel25': timesel2x}


""" The run names to consider from the first set are

kel09, kel11, kel14, kel15, kel17, kel18, kel20,* kel24 -- timeper '2921_2930'

The run names from the second set are

kel10, kel12, kel16, kel19, kel21, *kel25 -- timeper '2451_2460'

'controls'
iga -- 2921_2930
gregory_2xco2 -- 2451_2460
"""

ncfield=field.upper()


if field == 'sicn':
    cminc=0; cmaxc=1       # for climo
    cmin=-.75; cmax=.75    # for 2xco2 vs prei
    cminn=-.25; cmaxn=.25  # for nudge vs original
    ylims = (0,.08)
    conv=1
elif field == 'sic':
    cminc=0; cmaxc=2.5     # for climo
    cmin=-1; cmax=1    # for 2xco2 vs prei
    cminn=-.25; cmaxn=.25  # for nudge vs original
    ylims = (0,.15)
    conv=1/913.

fontP = fm.FontProperties()
fontP.set_size('small')

seasons = ('SON','DJF','MAM','JJA')
flddt = {}

for sii,skey in enumerate(sims):
    fname = basepath + skey + subdir + skey + '_' + field + '_' + timeperdt[skey] + '_ts.nc'

    if sii==0: # only get once
        time = cnc.getNCvar(fname,'time',timesel=timeseldt[skey])
        ntime = len(time)
        lat = cnc.getNCvar(fname,'lat')
        lon = cnc.getNCvar(fname,'lon')

    flddt[skey] = cnc.getNCvar(fname,ncfield,timesel=timeseldt[skey])*conv
    
# <headingcell level=3>

# SICN: Below are sea ice concentration comparisons for various combinations
# (figure titles are below figures)

# <codecell>

if showmaps:
    fldc = flddt['iga']

    # control (iga) climo
    fig1 = cplt.map_allmonths(fldc,lat,lon,type='nh',
                              cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,title='preI')
    if printtofile:
        fig1.savefig(field + '_iga_allmos_nh.pdf')

    for skey in sims[2:10]: # just the Group I sims, climo and diff
        fig2 = cplt.map_allmonths(flddt[skey],lat,lon,type='nh',
                                 cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,title='2xco2 nudge (' + skey + ')')
        fig3 = cplt.map_allmonths(flddt[skey]-fldc,lat,lon,type='nh',
                                 cmap='red2blue_w20',cmin=cmin,cmax=cmax,lmask=1,
                                  title='2xco2 nudge (' + skey + ') - preI')
        if printtofile:
            fig2.savefig(field + '_' + skey + '_allmos_nh.pdf')
            fig3.savefig(field + 'diff_' + skey + '_v_iga_allmos_nh.pdf')
    
    # <headingcell level=3>

    # Compare the original PreI run with nudging to 2xco2 ice in preI climate

    # <codecell>

    fldc2x = flddt['gregory_2xco2']

    # control (2xco2) climo
    fig1 = cplt.map_allmonths(fldc2x,lat,lon,type='nh',
                              cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,title='2xco2')
    if printtofile:
        fig1.savefig(field + '_gregory_2xco2_allmos_nh.pdf')
        
    for skey in sims[10:]: # just the Group II sims, climo and diff
        fig2 = cplt.map_allmonths(flddt[skey],lat,lon,type='nh',
                                 cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,title='preI nudge (' + skey + ')')
        fig3 = cplt.map_allmonths(flddt[skey]-fldc2x,lat,lon,type='nh',
                                 cmap='red2blue_w20',cmin=cmin,cmax=cmax,lmask=1,
                                  title='preI nudge (' + skey + ') - 2xco2')
        if printtofile:
            fig2.savefig(field + '_' + skey + '_allmos_nh.pdf')
            fig3.savefig(field + 'diff_' + skey + '_v_gregory_2xco2_allmos_nh.pdf')
            

# <headingcell level=3>

# Compare the original 2xco2 run with nudging to preI ice in 2xco2 climate

# <codecell>

# # Compare 2xco2 with the 2xco2 nudge: RMSE

# calculate RMSE 
rmsedt={}; rmseclimdt = {}; annrmseclimdt = {}; climrmsedt={}
region='polcap60'

fldclimdt = {}; flddiffdt = {}; flddiffclimdt = {}

# get controls
fldc2x = flddt['gregory_2xco2']
fldc2xclim,std = cutl.climatologize(fldc2x) # climo mean (2xco2)
fldclimdt['gregory_2xco2'] = fldc2xclim

fldc = flddt['iga']
fldcclim,std = cutl.climatologize(fldc) # climo mean (iga)
fldclimdt['iga'] = fldcclim

for skey in sims[2:10]: # Group I
    fld = flddt[skey]
    # timeseries
    flddiffdt[skey] = fldc2x-fld
    # climos
    fldclimdt[skey],std = cutl.climatologize(fld)  
    flddiffclimdt[skey] = fldc2xclim-fldclimdt[skey]
    
    rmse = np.sqrt(np.square(flddiffdt[skey]))
    rmseclim = np.sqrt(np.square(flddiffclimdt[skey]))
                           
    rmsedt[skey] = cutl.calc_regmean(rmse,lat,lon,region)
    rmseclimdt[skey] = cutl.calc_regmean(rmseclim,lat,lon,region)

    annrmseclimdt[skey] = cutl.annualize_monthlyts(rmseclimdt[skey])
    
    climrmsedt[skey],rmsestd = cutl.climatologize(rmsedt[skey])

for skey in sims[10:]: # Group II
    fld = flddt[skey]
    # timeseries
    flddiffdt[skey] = fldc-fld 
    # climos
    fldclimdt[skey],std = cutl.climatologize(fld)  
    flddiffclimdt[skey] = fldcclim-fldclimdt[skey] 
    
    rmse = np.sqrt(np.square(flddiffdt[skey]))
    rmseclim = np.sqrt(np.square(flddiffclimdt[skey]))
                           
    rmsedt[skey] = cutl.calc_regmean(rmse,lat,lon,region)
    rmseclimdt[skey] = cutl.calc_regmean(rmseclim,lat,lon,region)

    annrmseclimdt[skey] = cutl.annualize_monthlyts(rmseclimdt[skey])
    
    climrmsedt[skey],rmsestd = cutl.climatologize(rmsedt[skey])

# <codecell>

from matplotlib import gridspec
colors = ('b','g','m','y','c','k','r',ccm.get_linecolor('mediumpurple4'))
cdict = {'iga': '0.5', 'gregory_2xco2': '0.3',
         'kel11': 'b', 'kel09': 'g', 'kel14': 'm', 'kel15': 'y', 'kel17': 'c',
         'kel18': 'k', 'kel20': 'r', 'kel24': ccm.get_linecolor('limegreen'),
         'kel10': 'b', 'kel12': 'g', 'kel16': 'm', 'kel19': 'y', 'kel21': 'c',
         'kel25': ccm.get_linecolor('limegreen')}

mons=np.arange(1,13)

# RMSE figure: Group I
fig = plt.figure(figsize=(14,4))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2]) 

ax=plt.subplot(gs[0])
for sii,skey in enumerate(sims[2:10]):
    ax.plot(mons,rmseclimdt[skey],color=cdict[skey],linewidth=2) # RMSE of climatologies

for sii,skey in enumerate(sims[2:10]):
    ax.plot(mons,climrmsedt[skey],color=cdict[skey],linestyle='--',linewidth=2)# climatology of the RMSE in time

legstr=sims[2:10] + ('time avg RSME',)
ax.legend(legstr,loc='upper left', fancybox=True,
          prop=fontP,framealpha=0.5,ncol=2)
ax.set_title(region + ': RMSE of climo (solid)')
ax.set_ylim(ylims)
ax.set_xlim((1,12))

ax=plt.subplot(gs[1])
for sii,skey in enumerate(sims[2:10]):
    ax.plot(rmsedt[skey],color=cdict[skey],linewidth=2)# RMSE in time (over 9 years)

ax.legend(sims[2:10],loc='upper left',fancybox=True,
          prop=fontP,framealpha=0.5, ncol=2)
ax.set_title(region + ': RMSE in time')
ax.set_ylim(ylims)
ax.set_xlim((0,110))

if printtofile:
    fig.savefig(field + '_2xnudge_RMSE_climo_ts.pdf')


print 'RMSE over polar cap (>60N) by month'
rmseclimdf = pd.DataFrame(rmseclimdt)
g1sims=sims[2:10]
g2sims=sims[10:]

g1rmseclim=rmseclimdf.loc[:,g1sims]
g2rmseclim=rmseclimdf.loc[:,g2sims]

print g1rmseclim

#print 'ANN mean RMSE'

annrmseclim = rmseclimdf.mean(axis=0)
#print annrmseclim

print '\nANN mean RMSE (day-weighted)'
#print annrmseclimdt
annrmseclimdf = pd.DataFrame(annrmseclimdt)
g1annrmseclim = annrmseclimdf.loc[:,g1sims]
g2annrmseclim = annrmseclimdf.loc[:,g2sims]
print g1annrmseclim

# <rawcell>

# The left plot, solid lines show the monthly RMSE (averaged >60N) where the runs were first averaged into climos, and then the RMSE was calculated. The dashed lines show the RMSE climo (where the RMSE was calculated in time -- see right panel -- and then averaged into a climo). The above table of monthly RMSE values and the following annual means were calculated from the first method (RMSE of the climos). 'kel17' appears to be the winner with the smallest annual mean polar cap (>60N) average RMSE compared to the original state (gregory run). 

# <headingcell level=3>

# RMSE averaged >60N for each relaxation run (preI nudged to 2xco2 ice) compared to gregory_2xco2

# <codecell>

# calc total sea ice area/volume, plot it

mons=np.arange(1,13)

fldpltdt = {}; fldplttsdt = {}; fldpltclimdt = {}
if field=='sicn':

    for skey in sims:
        fldpltclimdt[skey],sh = cutl.calc_totseaicearea(fldclimdt[skey],lat,lon)
        fldplttsdt[skey],sh = cutl.calc_totseaicearea(flddt[skey],lat,lon)

    tstr = 'SIA (millions of km^2)'
    divby = np.float(1e12) # convert sea ice area to millions of km2

elif field=='sic':
    # need sea ice area to calc volume
    field='sicn'; ncfield=field.upper(); conv=1
    
    
    for skey in sims:
        fname = basepath + skey + subdir + skey + '_' + field + '_' + timeperdt[skey] + '_ts.nc'
        flds = cnc.getNCvar(fname,ncfield,timesel=timeseldt[skey])*conv # sea ice conc
        fldplttsdt[skey],sh=cutl.calc_totseaicevol(flddt[skey],flds,lat,lon) # time series of ice vol

        # climo the ice conc
        fldsclimo,std = cutl.climatologize(flds)
        # calc climo ice vol
        fldpltclimdt[skey],sh = cutl.calc_totseaicevol(fldclimdt[skey],fldsclimo,lat,lon) # climo ice vol

    field='sic'; ncfield=field.upper(); conv=1/913.
    tstr = 'SIV (m^3)'
    divby=1

# these comparisions ask how good is the nudging to 2xco2 ice (in a prei climate)
fldc2xplt=fldpltclimdt['gregory_2xco2']
fldc2xpltts = fldplttsdt['gregory_2xco2']
fldcplt = fldpltclimdt['iga']
fldcpltts = fldplttsdt['iga']

# as a baseline: what is the difference b/w 2xco2 and preindustrial ice?
diff2xpiplt = fldc2xplt-fldcplt

# plot seasonal cycle and timeseries of SIC or SICN: Group I
fig = plt.figure(figsize=(12,5))
gs = gridspec.GridSpec(3,3)
ax=plt.subplot(gs[0,0])

ax.plot(mons,fldc2xplt/divby,color='.5',linewidth=2) # control 2xco2 is gray
for sii,skey in enumerate(sims[2:10]):
    ax.plot(mons,(fldpltclimdt[skey])/divby,color=cdict[skey],linestyle='--',linewidth=2)

ax.set_xlim((1,12))
ax.set_title('NH ' + tstr)

ax=plt.subplot(gs[0,1])
for sii,skey in enumerate(sims[2:10]):
    ax.plot(mons,(fldc2xplt-fldpltclimdt[skey])/divby,color=cdict[skey],linewidth=2)

ax.set_xlim((1,12))
ax.set_title('NH diff: gregory-2xco2 nudges')

ax=plt.subplot(gs[0,2])
ax.plot(mons,diff2xpiplt/divby,color='.5',linewidth=3)
ax.set_xlim((1,12))
ax.set_title('NH diff: gregory-iga (for ref)')

# plot timeseries 
ax=plt.subplot(gs[1,:]) # plot DIFF values through time.
for sii,skey in enumerate(sims[2:10]):
    ax.plot((fldc2xpltts-fldplttsdt[skey])/divby,color=cdict[skey],linewidth=2)
ax.set_xlim((0,110))
ax.set_ylabel('anomaly')

ax=plt.subplot(gs[2,:]) # plot CLIMO values through time.
ax.plot(fldc2xpltts/divby,color='.5',linewidth=3)
for sii,skey in enumerate(sims[2:10]):
    ax.plot(fldplttsdt[skey]/divby,color=cdict[skey],linestyle='--',linewidth=2)
ax.set_xlim((0,110))
ax.set_xlabel('time (months)')
ax.set_ylabel('climo')

if printtofile:
    fig.savefig(field + 'diffclimo_2xnudge_seacyc_nh.pdf')


# <codecell>


if showmaps:
    # # Compare 2xco2 with the 2xco2 nudge: maps
    for skey in sims[2:10]:
        fig = cplt.map_allmonths(fldc2x-flddt[skey],lat,lon,type='nh',
                                 cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                                 title='2xco2 - 2xco2 nudge (' + skey + ')')
        if printtofile:
            fig.savefig(field + 'diff_gregory_2xco2_v_' + skey + '_allmos_nh.pdf')


# <rawcell>

# # Blue here means that the original run (2xco2 simulation) has more ice than the nudged run. 
# # This means that the nudging actually was too effective, such that more ice 
# # is melted than should be. This is mostly true in Nov and Dec. Otherwise it
# # appears the nudging is pretty good, on the whole.

# <headingcell level=3>

# Compare the original 2xco2 run with nudging to 2xco2 ice in preI climate (e.g. if it was perfect, there would be no anomalies)

# <codecell>

# # Compare preI with the preI nudge: RMSE

# <codecell>

# RMSE figure: Group II
fig = plt.figure(figsize=(14,4))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2]) 

ax=plt.subplot(gs[0])
for sii,skey in enumerate(sims[10:]):
    ax.plot(mons,rmseclimdt[skey],color=cdict[skey],linewidth=2)  # RMSE of climatologies
for sii,skey in enumerate(sims[10:]):
    ax.plot(mons,climrmsedt[skey],color=cdict[skey],linestyle='--',linewidth=2)# climatology of the RMSE in time

legstr=sims[10:] + ('time avg RMSE',)
ax.legend(legstr,loc='lower left',
          fancybox=True, prop=fontP,framealpha=0.5,ncol=2)
ax.set_title(region + ': RMSE of climo (solid)')

ax.set_xlim((1,12))
ax.set_ylim(ylims)

ax=plt.subplot(gs[1])
for sii,skey in enumerate(sims[10:]):
    ax.plot(rmsedt[skey],color=cdict[skey],linewidth=2)# RMSE in time (over 9 years)
    
ax.legend(sims[10:],
           loc='lower left',fancybox=True, prop=fontP,framealpha=0.5, ncol=2)
ax.set_title(region + ': RMSE in time')
ax.set_ylim(ylims)
ax.set_xlim((0,110))

if printtofile:
    fig.savefig(field + '_preinudge_RMSE_climo_ts.pdf')

print 'RMSE over polar cap (>60N) by month'
## rmseclimdfb = pd.DataFrame(rmseclimdictb)
## print rmseclimdfb
print g2rmseclim
## print 'ANN mean RMSE'

## annrmseclimb = rmseclimdfb.mean(axis=0)
## print annrmseclimb

print '\nANN mean RMSE (day-weighted)'
print g2annrmseclim

## print annrmseclimdictb

# <rawcell>

# Same as 'In [67]', but for comparing  preI nudge to original preI run (iga). Winner appears to be 'kel10'

# <headingcell level=3>

# RMSE averaged >60N for each relaxation run (2xco2 nudged to preI ice) compared to iga

# <codecell>

# seasonal cycle and timeseries of SIA or SIV: Group II

mons=np.arange(1,13)

fig = plt.figure(figsize=(12,5))
gs = gridspec.GridSpec(3,3)
ax=plt.subplot(gs[0,0])

ax.plot(mons,fldcplt/divby,color='.5',linewidth=2) # control prei is gray
for sii,skey in enumerate(sims[10:]):
    ax.plot(mons,fldpltclimdt[skey]/divby,color=cdict[skey],linestyle='--',linewidth=2)
    
ax.set_xlim((1,12))
ax.set_title('NH ' + tstr)


ax = plt.subplot(gs[0,1])
for sii,skey in enumerate(sims[10:]):
    ax.plot(mons,(fldcplt-fldpltclimdt[skey])/divby,color=cdict[skey],linewidth=2)
    
ax.set_xlim((1,12))
ax.set_title('NH diff: iga-prei nudges')

ax=plt.subplot(gs[0,2])
ax.plot(mons,diff2xpiplt/divby,color='.5',linewidth=3)
ax.set_xlim((1,12))
ax.set_title('NH diff: gregory-iga (for ref)')

ax=plt.subplot(gs[1,:]) # plot DIFF values through time.
for sii,skey in enumerate(sims[10:]):
    ax.plot((fldcpltts-fldplttsdt[skey])/divby,color=cdict[skey],linewidth=2)
ax.set_xlim((0,110))
ax.set_ylabel('anomaly')

ax=plt.subplot(gs[2,:]) # plot CLIMO values through time.
ax.plot(fldcpltts/divby,color='.5',linewidth=3)
for sii,skey in enumerate(sims[10:]):
    ax.plot(fldplttsdt[skey]/divby,color=cdict[skey],linestyle='--',linewidth=2)
ax.set_xlim((0,110))
ax.set_xlabel('time (months)')
ax.set_ylabel('climo')

if printtofile:
    fig.savefig(field + 'diffclimo_preinudge_seacyc_nh.pdf')


# <codecell>


if showmaps:
    # # Compare preI with the preI nudge
    for skey in sims[10:]:
        fig = cplt.map_allmonths(fldc-flddt[skey],lat,lon,type='nh',
                                 cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                                 title='preI - preI nudge (' + skey + ')')
        if printtofile:
            fig.savefig(field + 'diff_iga_v_' + skey + '_allmos_nh.pdf')
        

# <rawcell>

# # Blue means the original run (preindustrial control) has more ice than the nudge run.
# # This means that the nudge to preI is not strong enough, such that not enough ice 
# # builds back up as much as it should. By eye, this is especially true at the end 
# # of the growth season, including Apr, and in the marginal ice zone. 
# # The nudging is better just after minimum ice (in Nov-Dec)

# <headingcell level=3>

# Compare the original preI run with nudging to preI ice in 2xco2 climate (e.g. if it was perfect, there would be no anomalies)

# <codecell>

## # # # Compare relaxation runs to the AGCM runs # # # # # # # 
## #   The BCs were taken from CanESM historical ensemble mean

## # @@ probably don't need this in python script (just in notebook for ref)
## model = 'CanAM4'
## bp=con.get_basepath()
## basepath=bp['basepath'] + model + '/'; subdir=bp['subdir']

## casenamec='kemctl1' # 1979-89 climo BC
## casenamep='kem1pert2' # 2002-12 climo BC
## timeper='001-121'
## timesel = '0002-01-01,0121-12-31' 


## #field = 'sicn' # same field as above
## #ncfield=field.upper()


## fnamec = basepath + casenamec + subdir + casenamec + '_' + field + '_' + timeper + '_ts.nc'
## fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timeper + '_ts.nc'

## fldppres = cnc.getNCvar(fnamep,ncfield)*conv # present day
## fldchist = cnc.getNCvar(fnamec,ncfield)*conv # historical


## # <codecell>
## if showmaps:
##     fig = cplt.map_allmonths(fldppres,lat,lon,type='nh',
##                              cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,title='2002-12')
##     fig = cplt.map_allmonths(fldchist,lat,lon,type='nh',
##                              cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,title='1979-89')
##     fig = cplt.map_allmonths(fldppres-fldchist,lat,lon,type='nh',
##                              cmap='red2blue_w20',cmin=cmin,cmax=cmax,lmask=1,title='2002-2012 - 1979-89')

# <headingcell level=3>

# Our current AGCM runs with boundary conditions taken from CanESM2 historical ensemble mean for the given years. For comparison

