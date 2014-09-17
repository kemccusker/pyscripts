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
plt.ion()

printtofile=True
showmaps=False  # show all the monthly map comparisons?


#  set field here ==============
field = 'sicn'




model = 'CanESM2'
bp=con.get_basepath()
basepath=bp['basepath'] + model + '/'; subdir=bp['subdir']

casenamec='iga'; timeperc='2921_2930'
casenamec2x='gregory_2xco2'; timeper2x='2451_2460'

casenamepc='kel11'
casenamepc09='kel09'
casenamepc14='kel14'
casenamepc15='kel15'
casenamepc17='kel17'
casenamepc18='kel18'
casenamepc20='kel20'
timesel = '2922-01-01,2930-12-31' # skip first year

casenamep2x='kel10'
casenamep2x12='kel12'
casenamep2x16='kel16'
casenamep2x19='kel19'
casenamep2x21='kel21'

timesel2x = '2452-01-01,2460-12-31' # skip first year


""" The run names to consider from the first set are

kel09, kel11, kel14, kel15, kel17, kel18, kel20 -- timeper '2921_2930'

The run names from the second set are

kel10, kel12, kel16, kel19, kel21 -- timeper '2451_2460'

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

fnamec = basepath + casenamec + subdir + casenamec + '_' + field + '_' + timeperc + '_ts.nc'
fnamepc = basepath + casenamepc + subdir + casenamepc + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
fnamepc09 = basepath + casenamepc09 + subdir + casenamepc09 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
fnamepc14 = basepath + casenamepc14 + subdir + casenamepc14 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
fnamepc15 = basepath + casenamepc15 + subdir + casenamepc15 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
fnamepc17 = basepath + casenamepc17 + subdir + casenamepc17 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
fnamepc18 = basepath + casenamepc18 + subdir + casenamepc18 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
fnamepc20 = basepath + casenamepc20 + subdir + casenamepc20 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl


fnamec2x = basepath + casenamec2x + subdir + casenamec2x + '_' + field + '_' + timeper2x + '_ts.nc' 
fnamep2x = basepath + casenamep2x + subdir + casenamep2x + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x
fnamep2x12 = basepath + casenamep2x12 + subdir + casenamep2x12 + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x
fnamep2x16 = basepath + casenamep2x16 + subdir + casenamep2x16 + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x
fnamep2x19 = basepath + casenamep2x19 + subdir + casenamep2x19 + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x
fnamep2x21 = basepath + casenamep2x21 + subdir + casenamep2x21 + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x


time = cnc.getNCvar(fnamepc,'time',timesel=timesel)
ntime = len(time)

lat = cnc.getNCvar(fnamepc,'lat')
lon = cnc.getNCvar(fnamepc,'lon')

fldc = cnc.getNCvar(fnamec,ncfield,timesel=timesel)*conv # "control" = iga = preindustrial
fldpc = cnc.getNCvar(fnamepc,ncfield,timesel=timesel)*conv # "pert" matched to this control = 2xco2 nudge
fldpc09 = cnc.getNCvar(fnamepc09,ncfield,timesel=timesel)*conv # "pert" kel09 matched to this control = 2xco2 nudge
fldpc14 = cnc.getNCvar(fnamepc14,ncfield,timesel=timesel)*conv # "pert" kel14 matched to this control = 2xco2 nudge
fldpc15 = cnc.getNCvar(fnamepc15,ncfield,timesel=timesel)*conv # "pert" kel15 matched to this control = 2xco2 nudge
fldpc17 = cnc.getNCvar(fnamepc17,ncfield,timesel=timesel)*conv # "pert" kel17 matched to this control = 2xco2 nudge
fldpc18 = cnc.getNCvar(fnamepc18,ncfield,timesel=timesel)*conv # "pert" kel18 matched to this control = 2xco2 nudge
fldpc20 = cnc.getNCvar(fnamepc20,ncfield,timesel=timesel)*conv # "pert" kel20 matched to this control = 2xco2 nudge


fldc2x = cnc.getNCvar(fnamec2x,ncfield,timesel=timesel2x)*conv # "2xco2 control" = gregory_2xco2
fldp2x = cnc.getNCvar(fnamep2x,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge
fldp2x12 = cnc.getNCvar(fnamep2x12,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge
fldp2x16 = cnc.getNCvar(fnamep2x16,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge
fldp2x19 = cnc.getNCvar(fnamep2x19,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge
fldp2x21 = cnc.getNCvar(fnamep2x21,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge

# <headingcell level=3>

# SICN: Below are sea ice concentration comparisons for various combinations
# (figure titles are below figures)

# <codecell>

if showmaps:
    fig1 = cplt.map_allmonths(fldc,lat,lon,type='nh',
                             cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,title='preI')    
    # kel11
    fig2 = cplt.map_allmonths(fldpc,lat,lon,type='nh',
                             cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,title='2xco2 nudge (' + casenamepc + ')')
    fig3 = cplt.map_allmonths(fldpc-fldc,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cmin,cmax=cmax,lmask=1,
                              title='2xco2 nudge (' + casenamepc + ') - preI')
    if printtofile:
        fig1.savefig(field + '_' + casenamec + '_allmos_nh.pdf')
        fig2.savefig(field + '_' + casenamepc + '_allmos_nh.pdf')
        fig3.savefig(field + 'diff_' + casenamepc + '_v_' + casenamec + '_allmos_nh.pdf')
    

    # kel09 
    fig4 = cplt.map_allmonths(fldpc09,lat,lon,type='nh',
                             cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,
                             title='2xco2 nudge (' + casenamepc09 + ')')

    fig5 = cplt.map_allmonths(fldpc09-fldc,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cmin,cmax=cmax,lmask=1,
                             title='2xco2 nudge (' + casenamepc09 + ') - preI')    
    # kel14
    fig6 = cplt.map_allmonths(fldpc14,lat,lon,type='nh',
                             cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,
                             title='2xco2 nudge (' + casenamepc14 + ')')

    fig7 = cplt.map_allmonths(fldpc14-fldc,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cmin,cmax=cmax,lmask=1,
                             title='2xco2 nudge (' + casenamepc14 + ') - preI')
    if printtofile:
        fig4.savefig(field + '_' + casenamepc09 + '_allmos_nh.pdf')
        fig5.savefig(field + 'diff_' + casenamepc09 + '_v_' + casenamec + '_allmos_nh.pdf')
        fig6.savefig(field + '_' + casenamepc14 + '_allmos_nh.pdf')
        fig7.savefig(field + 'diff_' + casenamepc14 + '_v_' + casenamec + '_allmos_nh.pdf')

    # <headingcell level=3>

    # Compare the original PreI run with nudging to 2xco2 ice in preI climate

    # <codecell>

    fig = cplt.map_allmonths(fldc2x,lat,lon,type='nh',
                             cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,title='2xco2')
    if printtofile:
        fig.savefig(field + '_' + casenamec2x + '_allmos_nh.pdf')
    fig = cplt.map_allmonths(fldp2x,lat,lon,type='nh',
                             cmap='blue2blue_bw10',cmin=cminc,cmax=cmaxc,lmask=1,
                             title='preI nudge (' + casenamep2x + ')')
    if printtofile:
        fig.savefig(field + '_' + casenamep2x + '_allmos_nh.pdf')
    fig = cplt.map_allmonths(fldp2x-fldc2x,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cmin,cmax=cmax,lmask=1,
                             title='preI nudge (' + casenamep2x + ')- 2xco2')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamep2x + '_v_' + casenamec2x + '_allmos_nh.pdf')

# <headingcell level=3>

# Compare the original 2xco2 run with nudging to preI ice in 2xco2 climate

# <codecell>

# # Compare 2xco2 with the 2xco2 nudge: RMSE

# calculate RMSE 
rmsedict={}; rmseclimdict = {}; annrmseclimdict = {}
region='polcap60'

fldc2xclim,fldc2xstd = cutl.climatologize(fldc2x) # climo mean
fldpcclim,fldpcstd = cutl.climatologize(fldpc) # climo mean, kel11
fldpc09clim,fldpc09std = cutl.climatologize(fldpc09) # climo mean, kel09
fldpc14clim,fldpc14std = cutl.climatologize(fldpc14) # climo mean, kel14
fldpc15clim,fldpc15std = cutl.climatologize(fldpc15) # climo mean
fldpc17clim,fldpc17std = cutl.climatologize(fldpc17) # climo mean
fldpc18clim,fldpc18std = cutl.climatologize(fldpc18) # climo mean
fldpc20clim,fldpc20std = cutl.climatologize(fldpc20) # climo mean


diff11 = fldc2x-fldpc
diff09 = fldc2x-fldpc09
diff14 = fldc2x-fldpc14
diff15 = fldc2x-fldpc15
diff17 = fldc2x-fldpc17
diff18 = fldc2x-fldpc18
diff20 = fldc2x-fldpc20


diff11clim = fldc2xclim-fldpcclim
diff09clim = fldc2xclim-fldpc09clim
diff14clim = fldc2xclim-fldpc14clim
diff15clim = fldc2xclim-fldpc15clim
diff17clim = fldc2xclim-fldpc17clim
diff18clim = fldc2xclim-fldpc18clim
diff20clim = fldc2xclim-fldpc20clim

rmse11 = np.sqrt(np.square(diff11))
rmse09 = np.sqrt(np.square(diff09))
rmse14 = np.sqrt(np.square(diff14))
rmse15 = np.sqrt(np.square(diff15))
rmse17 = np.sqrt(np.square(diff17))
rmse18 = np.sqrt(np.square(diff18))
rmse20 = np.sqrt(np.square(diff20))


rmse11clim = np.sqrt(np.square(diff11clim))
rmse09clim = np.sqrt(np.square(diff09clim))
rmse14clim = np.sqrt(np.square(diff14clim))
rmse15clim = np.sqrt(np.square(diff15clim))
rmse17clim = np.sqrt(np.square(diff17clim))
rmse18clim = np.sqrt(np.square(diff18clim))
rmse20clim = np.sqrt(np.square(diff20clim))

rmsedict[casenamepc] = cutl.calc_regmean(rmse11,lat,lon,region)
rmsedict[casenamepc09] = cutl.calc_regmean(rmse09,lat,lon,region)
rmsedict[casenamepc14] = cutl.calc_regmean(rmse14,lat,lon,region)
rmsedict[casenamepc15] = cutl.calc_regmean(rmse15,lat,lon,region)
rmsedict[casenamepc17] = cutl.calc_regmean(rmse17,lat,lon,region)
rmsedict[casenamepc18] = cutl.calc_regmean(rmse18,lat,lon,region)
rmsedict[casenamepc20] = cutl.calc_regmean(rmse20,lat,lon,region)


rmseclimdict[casenamepc] = cutl.calc_regmean(rmse11clim,lat,lon,region)
rmseclimdict[casenamepc09] = cutl.calc_regmean(rmse09clim,lat,lon,region)
rmseclimdict[casenamepc14] = cutl.calc_regmean(rmse14clim,lat,lon,region)
rmseclimdict[casenamepc15] = cutl.calc_regmean(rmse15clim,lat,lon,region)
rmseclimdict[casenamepc17] = cutl.calc_regmean(rmse17clim,lat,lon,region)
rmseclimdict[casenamepc18] = cutl.calc_regmean(rmse18clim,lat,lon,region)
rmseclimdict[casenamepc20] = cutl.calc_regmean(rmse20clim,lat,lon,region)


annrmseclimdict[casenamepc] = cutl.annualize_monthlyts(rmseclimdict[casenamepc])
annrmseclimdict[casenamepc09] = cutl.annualize_monthlyts(rmseclimdict[casenamepc09])
annrmseclimdict[casenamepc14] = cutl.annualize_monthlyts(rmseclimdict[casenamepc14])
annrmseclimdict[casenamepc15] = cutl.annualize_monthlyts(rmseclimdict[casenamepc15])
annrmseclimdict[casenamepc17] = cutl.annualize_monthlyts(rmseclimdict[casenamepc17])
annrmseclimdict[casenamepc18] = cutl.annualize_monthlyts(rmseclimdict[casenamepc18])
annrmseclimdict[casenamepc20] = cutl.annualize_monthlyts(rmseclimdict[casenamepc20])


climrmse,rmsestd = cutl.climatologize(rmsedict['kel11'])
climrmse09,rmse09std = cutl.climatologize(rmsedict['kel09'])
climrmse14,rmse14std = cutl.climatologize(rmsedict['kel14'])
climrmse15,rmse15std = cutl.climatologize(rmsedict['kel15'])
climrmse17,rmse17std = cutl.climatologize(rmsedict['kel17'])
climrmse18,rmse18std = cutl.climatologize(rmsedict['kel18'])
climrmse20,rmse20std = cutl.climatologize(rmsedict['kel20'])

# <codecell>

from matplotlib import gridspec

mons=np.arange(1,13)

fig = plt.figure(figsize=(14,4))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2]) 


ax=plt.subplot(gs[0])
ax.plot(mons,rmseclimdict['kel11'],'b',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdict['kel09'],'g',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdict['kel14'],'m',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdict['kel15'],'y',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdict['kel17'],'c',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdict['kel18'],'k',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdict['kel20'],'r',linewidth=2) # RMSE of climatologies

ax.plot(mons,climrmse,'b--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse09,'g--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse14,'m--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse15,'y--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse17,'c--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse18,'k--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse20,'r--',linewidth=2) # climatology of the RMSE in time


ax.legend(('kel11',
           'kel09','kel14',
           'kel15','kel17',
           'kel18', 'kel20',
           'time avg RMSE'),loc='upper left',
          fancybox=True, prop=fontP,framealpha=0.5,ncol=2)
ax.set_title(region + ': RMSE of climo (solid)')

ax.set_ylim(ylims)
ax.set_xlim((1,12))

ax=plt.subplot(gs[1])
#ax=axs[1]
ax.plot(rmsedict['kel11'],'b',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedict['kel09'],'g',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedict['kel14'],'m',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedict['kel15'],'y',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedict['kel17'],'c',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedict['kel18'],'k',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedict['kel20'],'r',linewidth=2) # RMSE in time (over 9 years)

ax.legend(('kel11','kel09', 
            'kel14', 'kel15',
            'kel17', 'kel18','kel20'),
           loc='upper left',fancybox=True, prop=fontP,framealpha=0.5, ncol=2)
ax.set_title(region + ': RMSE in time')
ax.set_ylim(ylims)
ax.set_xlim((0,110))

if printtofile:
    fig.savefig('2xnudge_RMSE_climo_ts.pdf')

print 'RMSE over polar cap (>60N) by month'
rmseclimdf = pd.DataFrame(rmseclimdict)
print rmseclimdf

print 'ANN mean RMSE'

annrmseclim = rmseclimdf.mean(axis=0)
print annrmseclim

print '\nANN mean RMSE (day-weighted)'
print annrmseclimdict

# <rawcell>

# The left plot, solid lines show the monthly RMSE (averaged >60N) where the runs were first averaged into climos, and then the RMSE was calculated. The dashed lines show the RMSE climo (where the RMSE was calculated in time -- see right panel -- and then averaged into a climo). The above table of monthly RMSE values and the following annual means were calculated from the first method (RMSE of the climos). 'kel17' appears to be the winner with the smallest annual mean polar cap (>60N) average RMSE compared to the original state (gregory run). 

# <headingcell level=3>

# RMSE averaged >60N for each relaxation run (preI nudged to 2xco2 ice) compared to gregory_2xco2

# <codecell>

# calc total sea ice area/volume, plot it

mons=np.arange(1,13)

fldcclim,fldcstd = cutl.climatologize(fldc) # climo mean (prei)

if field=='sicn':
    fldc2xplt,sh=cutl.calc_totseaicearea(fldc2xclim,lat,lon) # 2xco2 control
    fldcplt,sh=cutl.calc_totseaicearea(fldcclim,lat,lon) # preindustrial control

    # these perturbed runs are attempting to match the above 2xco2 control.
    # they are labelled 'pcxx' to denote that they are perts, 'p', 
    #    to be compared with standard control, 'c' (preindustrial run, 'iga')
    fldpcplt,sh=cutl.calc_totseaicearea(fldpcclim,lat,lon)
    fldpc09plt,sh=cutl.calc_totseaicearea(fldpc09clim,lat,lon)
    fldpc14plt,sh=cutl.calc_totseaicearea(fldpc14clim,lat,lon)
    fldpc15plt,sh=cutl.calc_totseaicearea(fldpc15clim,lat,lon)
    fldpc17plt,sh=cutl.calc_totseaicearea(fldpc17clim,lat,lon)
    fldpc18plt,sh=cutl.calc_totseaicearea(fldpc18clim,lat,lon)
    fldpc20plt,sh=cutl.calc_totseaicearea(fldpc20clim,lat,lon)

    tstr = 'SIA (millions of km^2)'
    divby = np.float(1e12) # convert sea ice area to millions of km2

elif field=='sic':
    # need sea ice area
    field='sicn'; ncfield=field.upper(); conv=1
    fnamec = basepath + casenamec + subdir + casenamec + '_' + field + '_' + timeperc + '_ts.nc'
    fnamepc = basepath + casenamepc + subdir + casenamepc + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
    fnamepc09 = basepath + casenamepc09 + subdir + casenamepc09 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
    fnamepc14 = basepath + casenamepc14 + subdir + casenamepc14 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
    fnamepc15 = basepath + casenamepc15 + subdir + casenamepc15 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
    fnamepc17 = basepath + casenamepc17 + subdir + casenamepc17 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
    fnamepc18 = basepath + casenamepc18 + subdir + casenamepc18 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl
    fnamepc20 = basepath + casenamepc20 + subdir + casenamepc20 + '_' + field + '_' + timeperc + '_ts.nc' # match with ctl

    # note that sea ice conc is read in here for use in calc below (nudging to prei ice)
    
    fnamec2x = basepath + casenamec2x + subdir + casenamec2x + '_' + field + '_' + timeper2x + '_ts.nc' 
    fnamep2x = basepath + casenamep2x + subdir + casenamep2x + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x
    fnamep2x12 = basepath + casenamep2x12 + subdir + casenamep2x12 + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x
    fnamep2x16 = basepath + casenamep2x16 + subdir + casenamep2x16 + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x
    fnamep2x19 = basepath + casenamep2x19 + subdir + casenamep2x19 + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x
    fnamep2x21 = basepath + casenamep2x21 + subdir + casenamep2x21 + '_' + field + '_' + timeper2x + '_ts.nc' # match with 2x

    fldcs = cnc.getNCvar(fnamec,ncfield,timesel=timesel)*conv # "control" = iga = preindustrial
    fldpcs = cnc.getNCvar(fnamepc,ncfield,timesel=timesel)*conv # "pert" matched to this control = 2xco2 nudge
    fldpc09s = cnc.getNCvar(fnamepc09,ncfield,timesel=timesel)*conv # "pert" kel09 matched to this control = 2xco2 nudge
    fldpc14s = cnc.getNCvar(fnamepc14,ncfield,timesel=timesel)*conv # "pert" kel14 matched to this control = 2xco2 nudge
    fldpc15s = cnc.getNCvar(fnamepc15,ncfield,timesel=timesel)*conv # "pert" kel15 matched to this control = 2xco2 nudge
    fldpc17s = cnc.getNCvar(fnamepc17,ncfield,timesel=timesel)*conv # "pert" kel17 matched to this control = 2xco2 nudge
    fldpc18s = cnc.getNCvar(fnamepc18,ncfield,timesel=timesel)*conv # "pert" kel18 matched to this control = 2xco2 nudge
    fldpc20s = cnc.getNCvar(fnamepc20,ncfield,timesel=timesel)*conv # "pert" kel20 matched to this control = 2xco2 nudge

    fldc2xs = cnc.getNCvar(fnamec2x,ncfield,timesel=timesel2x)*conv # "2xco2 control" = gregory_2xco2
    fldp2xs = cnc.getNCvar(fnamep2x,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge
    fldp2x12s = cnc.getNCvar(fnamep2x12,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge
    fldp2x16s = cnc.getNCvar(fnamep2x16,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge
    fldp2x19s = cnc.getNCvar(fnamep2x19,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge
    fldp2x21s = cnc.getNCvar(fnamep2x21,ncfield,timesel=timesel2x)*conv # "pert" matched to this control = preI nudge

    # ==============================
    # =============================
    # @@@@@@@@@@ add timeseries calcs here @@@@@@@@@@@@@@@@
    # @@actually want climo, but tile it.#fldc2xpltts,sh = cutl.calc_totseaicevol(fldc2x,fldc2xs,lat,lon) # 2xco2 control timeseries
    
    fldpcpltts,sh=cutl.calc_totseaicevol(fldpc,fldpcs,lat,lon) # timeseries
    fldpc09pltts,sh=cutl.calc_totseaicevol(fldpc09,fldpc09s,lat,lon)
    # @@@@ add more....


    fldc2xsclim,fldc2xsstd = cutl.climatologize(fldc2xs) # climo mean, gregory
    fldcsclim,fldcsstd = cutl.climatologize(fldcs) # climo mean, iga
    fldpcsclim,fldpcsstd = cutl.climatologize(fldpcs) # climo mean, kel11
    fldpc09sclim,fldpc09sstd = cutl.climatologize(fldpc09s) # climo mean, kel09
    fldpc14sclim,fldpc14sstd = cutl.climatologize(fldpc14s) # climo mean, kel14
    fldpc15sclim,fldpc15sstd = cutl.climatologize(fldpc15s) # climo mean
    fldpc17sclim,fldpc17sstd = cutl.climatologize(fldpc17s) # climo mean
    fldpc18sclim,fldpc18sstd = cutl.climatologize(fldpc18s) # climo mean
    fldpc20sclim,fldpc20sstd = cutl.climatologize(fldpc20s) # climo mean
    
    fldc2xplt,sh=cutl.calc_totseaicevol(fldc2xclim,fldc2xsclim,lat,lon) # 2xco2 control
    fldcplt,sh=cutl.calc_totseaicevol(fldcclim,fldcsclim,lat,lon) # preindustrial control

    # these perturbed runs are attempting to match the above 2xco2 control.
    # they are labelled 'pcxx' to denote that they are perts, 'p', 
    #    to be compared with standard control, 'c' (preindustrial run, 'iga')
    fldpcplt,sh=cutl.calc_totseaicevol(fldpcclim,fldpcsclim,lat,lon)
    fldpc09plt,sh=cutl.calc_totseaicevol(fldpc09clim,fldpc09sclim,lat,lon)
    fldpc14plt,sh=cutl.calc_totseaicevol(fldpc14clim,fldpc14sclim,lat,lon)
    fldpc15plt,sh=cutl.calc_totseaicevol(fldpc15clim,fldpc15sclim,lat,lon)
    fldpc17plt,sh=cutl.calc_totseaicevol(fldpc17clim,fldpc17sclim,lat,lon)
    fldpc18plt,sh=cutl.calc_totseaicevol(fldpc18clim,fldpc18sclim,lat,lon)
    fldpc20plt,sh=cutl.calc_totseaicevol(fldpc20clim,fldpc20sclim,lat,lon)

    field='sic'; ncfield=field.upper(); conv=1/913.
    tstr = 'SIV (m^3)'
    divby=1

# these comparisions ask how good is the nudging to 2xco2 ice (in a prei climate)
diff11plt = fldc2xplt-fldpcplt
diff09plt = fldc2xplt-fldpc09plt
diff14plt = fldc2xplt-fldpc14plt
diff15plt = fldc2xplt-fldpc15plt
diff17plt = fldc2xplt-fldpc17plt
diff18plt = fldc2xplt-fldpc18plt
diff20plt = fldc2xplt-fldpc20plt

# as a baseline: what is the difference b/w 2xco2 and preindustrial ice?
diff2xpiplt = fldc2xplt-fldcplt

fig = plt.figure(figsize=(12,4))
#fig,axs = plt.subplots(1,3)
gs = gridspec.GridSpec(2,3)
ax=plt.subplot(gs[0,0])

#fig.set_size_inches(12,4)
#ax=axs[0]

ax.plot(mons,fldc2xplt/divby,color='.5',linewidth=2) # control 2xco2 is gray
ax.plot(mons,fldpcplt/divby,'b--',linewidth=2) # blue, green, magenta, yellow,cyan, black
ax.plot(mons,fldpc09plt/divby,'g--',linewidth=2)
ax.plot(mons,fldpc14plt/divby,'m--',linewidth=2)
ax.plot(mons,fldpc15plt/divby,'y--',linewidth=2)
ax.plot(mons,fldpc17plt/divby,'c--',linewidth=2)
ax.plot(mons,fldpc18plt/divby,'k--',linewidth=2)
ax.plot(mons,fldpc20plt/divby,'r--',linewidth=2)
ax.set_xlim((1,12))
ax.set_title('NH ' + tstr)


#ax = axs[1]
ax=plt.subplot(gs[0,1])
ax.plot(mons,diff11plt/divby,'b',linewidth=2)
ax.plot(mons,diff09plt/divby,'g',linewidth=2)
ax.plot(mons,diff14plt/divby,'m',linewidth=2)
ax.plot(mons,diff15plt/divby,'y',linewidth=2)
ax.plot(mons,diff17plt/divby,'c',linewidth=2)
ax.plot(mons,diff18plt/divby,'k',linewidth=2)
ax.plot(mons,diff20plt/divby,'r',linewidth=2)
ax.set_xlim((1,12))
ax.set_title('NH diff: gregory-2xco2 nudges')

#ax=axs[2]
ax=plt.subplot(gs[0,2])
ax.plot(mons,diff2xpiplt/divby,color='.5',linewidth=3)
ax.set_xlim((1,12))
ax.set_title('NH diff: gregory-iga (for ref)')

# =======================================
# @@@@@@@@@@@@@ plot timeseries here @@@@@@@@@
ax=plt.subplot(gs[1,:]) # plot values through time.
## ax.plot(diff11pltts/divby,'b',linewidth=2)
## ax.plot(diff09pltts/divby,'g',linewidth=2)
## ax.plot(diff14pltts/divby,'m',linewidth=2)
## ax.plot(diff15pltts/divby,'y',linewidth=2)
## ax.plot(diff17pltts/divby,'c',linewidth=2)
## ax.plot(diff18pltts/divby,'k',linewidth=2)
## ax.plot(diff20pltts/divby,'r',linewidth=2)
    
if printtofile:
    fig.savefig(field + 'diffclimo_2xnudge_seacyc_nh.pdf')


# <codecell>

# # Compare 2xco2 with the 2xco2 nudge: maps

if showmaps:
    # kel11
    fig = cplt.map_allmonths(fldc2x-fldpc,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='2xco2 - 2xco2 nudge (' + casenamepc + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec2x + '_v_' + casenamepc + '_allmos_nh.pdf')

    # kel09
    fig = cplt.map_allmonths(fldc2x-fldpc09,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='2xco2 - 2xco2 nudge (' + casenamepc09 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec2x + '_v_' + casenamepc09 + '_allmos_nh.pdf')

    # kel14
    fig = cplt.map_allmonths(fldc2x-fldpc14,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='2xco2 - 2xco2 nudge (' + casenamepc14 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec2x + '_v_' + casenamepc14 + '_allmos_nh.pdf')

    # kel15
    fig = cplt.map_allmonths(fldc2x-fldpc15,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='2xco2 - 2xco2 nudge (' + casenamepc15 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec2x + '_v_' + casenamepc15 + '_allmos_nh.pdf')

    # kel17
    fig = cplt.map_allmonths(fldc2x-fldpc17,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='2xco2 - 2xco2 nudge (' + casenamepc17 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec2x + '_v_' + casenamepc17 + '_allmos_nh.pdf')

    # kel18
    fig = cplt.map_allmonths(fldc2x-fldpc18,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='2xco2 - 2xco2 nudge (' + casenamepc18 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec2x + '_v_' + casenamepc18 + '_allmos_nh.pdf')

    # kel20
    fig = cplt.map_allmonths(fldc2x-fldpc20,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='2xco2 - 2xco2 nudge (' + casenamepc20 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec2x + '_v_' + casenamepc20 + '_allmos_nh.pdf')

# <rawcell>

# # Blue here means that the original run (2xco2 simulation) has more ice than the nudged run. 
# # This means that the nudging actually was too effective, such that more ice 
# # is melted than should be. This is mostly true in Nov and Dec. Otherwise it
# # appears the nudging is pretty good, on the whole.

# <headingcell level=3>

# Compare the original 2xco2 run with nudging to 2xco2 ice in preI climate (e.g. if it was perfect, there would be no anomalies)

# <codecell>

# # Compare preI with the preI nudge: RMSE

# calculate RMSE 
rmsedictb={}; rmseclimdictb = {}; annrmseclimdictb = {}
region='polcap60'

fldcclim,fldcstd = cutl.climatologize(fldc) # climo mean
fldp2xclim,fldp2xcstd = cutl.climatologize(fldp2x) # climo mean, kel10
fldp2x12clim,fldp2x12std = cutl.climatologize(fldp2x12) # climo mean, kel12
fldp2x16clim,fldp2x16std = cutl.climatologize(fldp2x16) # climo mean, kel16
fldp2x19clim,fldp2x19std = cutl.climatologize(fldp2x19) # climo mean, kel19
fldp2x21clim,fldp2x21std = cutl.climatologize(fldp2x21) # climo mean, kel21


diff10 = fldc-fldp2x
diff12 = fldc-fldp2x12
diff16 = fldc-fldp2x16
diff19 = fldc-fldp2x19
diff21 = fldc-fldp2x21

diff10clim = fldcclim-fldp2xclim
diff12clim = fldcclim-fldp2x12clim
diff16clim = fldcclim-fldp2x16clim
diff19clim = fldcclim-fldp2x19clim
diff21clim = fldcclim-fldp2x21clim


rmse10 = np.sqrt(np.square(diff10))
rmse12 = np.sqrt(np.square(diff12))
rmse16 = np.sqrt(np.square(diff16))
rmse19 = np.sqrt(np.square(diff19))
rmse21 = np.sqrt(np.square(diff21))


rmse10clim = np.sqrt(np.square(diff10clim))
rmse12clim = np.sqrt(np.square(diff12clim))
rmse16clim = np.sqrt(np.square(diff16clim))
rmse19clim = np.sqrt(np.square(diff19clim))
rmse21clim = np.sqrt(np.square(diff21clim))


rmsedictb[casenamep2x] = cutl.calc_regmean(rmse10,lat,lon,region)
rmsedictb[casenamep2x12] = cutl.calc_regmean(rmse12,lat,lon,region)
rmsedictb[casenamep2x16] = cutl.calc_regmean(rmse16,lat,lon,region)
rmsedictb[casenamep2x19] = cutl.calc_regmean(rmse19,lat,lon,region)
rmsedictb[casenamep2x21] = cutl.calc_regmean(rmse21,lat,lon,region)

rmseclimdictb[casenamep2x] = cutl.calc_regmean(rmse10clim,lat,lon,region)
rmseclimdictb[casenamep2x12] = cutl.calc_regmean(rmse12clim,lat,lon,region)
rmseclimdictb[casenamep2x16] = cutl.calc_regmean(rmse16clim,lat,lon,region)
rmseclimdictb[casenamep2x19] = cutl.calc_regmean(rmse19clim,lat,lon,region)
rmseclimdictb[casenamep2x21] = cutl.calc_regmean(rmse21clim,lat,lon,region)


annrmseclimdictb[casenamep2x] = cutl.annualize_monthlyts(rmseclimdictb[casenamep2x])
annrmseclimdictb[casenamep2x12] = cutl.annualize_monthlyts(rmseclimdictb[casenamep2x12])
annrmseclimdictb[casenamep2x16] = cutl.annualize_monthlyts(rmseclimdictb[casenamep2x16])
annrmseclimdictb[casenamep2x19] = cutl.annualize_monthlyts(rmseclimdictb[casenamep2x19])
annrmseclimdictb[casenamep2x21] = cutl.annualize_monthlyts(rmseclimdictb[casenamep2x21])


climrmse10,rmse10std = cutl.climatologize(rmsedictb['kel10'])
climrmse12,rmse12std = cutl.climatologize(rmsedictb['kel12'])
climrmse16,rmse16std = cutl.climatologize(rmsedictb['kel16'])
climrmse19,rmse19std = cutl.climatologize(rmsedictb['kel19'])
climrmse21,rmse21std = cutl.climatologize(rmsedictb['kel21'])


# <codecell>


fig = plt.figure(figsize=(14,4))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2]) 

ax=plt.subplot(gs[0])
ax.plot(mons,rmseclimdictb['kel10'],'b',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdictb['kel12'],'g',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdictb['kel16'],'m',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdictb['kel19'],'y',linewidth=2) # RMSE of climatologies
ax.plot(mons,rmseclimdictb['kel21'],'c',linewidth=2) # RMSE of climatologies

ax.plot(mons,climrmse10,'b--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse12,'g--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse16,'m--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse19,'y--',linewidth=2) # climatology of the RMSE in time
ax.plot(mons,climrmse21,'c--',linewidth=2) # climatology of the RMSE in time

ax.legend(('kel10',
           'kel12','kel16',
           'kel19','kel21',
           'time avg RMSE'),loc='lower left',
          fancybox=True, prop=fontP,framealpha=0.5,ncol=2)
ax.set_title(region + ': RMSE of climo (solid)')

ax.set_xlim((1,12))
ax.set_ylim(ylims)

ax=plt.subplot(gs[1])
#ax=axs[1]
ax.plot(rmsedictb['kel10'],'b',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedictb['kel12'],'g',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedictb['kel16'],'m',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedictb['kel19'],'y',linewidth=2) # RMSE in time (over 9 years)
ax.plot(rmsedictb['kel21'],'c',linewidth=2) # RMSE in time (over 9 years)

ax.legend(('kel10','kel12', 
            'kel16', 'kel19','kel21'),
           loc='lower left',fancybox=True, prop=fontP,framealpha=0.5, ncol=2)
ax.set_title(region + ': RMSE in time')
ax.set_ylim(ylims)
ax.set_xlim((0,110))

if printtofile:
    fig.savefig('preinudge_RMSE_climo_ts.pdf')

print 'RMSE over polar cap (>60N) by month'
rmseclimdfb = pd.DataFrame(rmseclimdictb)
print rmseclimdfb

print 'ANN mean RMSE'

annrmseclimb = rmseclimdfb.mean(axis=0)
print annrmseclimb

print '\nANN mean RMSE (day-weighted)'
print annrmseclimdictb

# <rawcell>

# Same as 'In [67]', but for comparing  preI nudge to original preI run (iga). Winner appears to be 'kel10'

# <headingcell level=3>

# RMSE averaged >60N for each relaxation run (2xco2 nudged to preI ice) compared to iga

# <codecell>

# calc total sea ice area, plot it
#diff10clim = fldcclim-fldp2xclim
#diff12clim = fldcclim-fldp2x12clim
#diff16clim = fldcclim-fldp2x16clim
#diff19clim = fldcclim-fldp2x19clim

# @@@ calc sea ice thickness here
mons=np.arange(1,13)

if field=='sicn':
    
    fldc2xplt,sh=cutl.calc_totseaicearea(fldc2xclim,lat,lon) # 2xco2 control
    fldcplt,sh=cutl.calc_totseaicearea(fldcclim,lat,lon) # preindustrial control

    # these perturbed runs are attempting to match the above prei control.
    # they are labelled 'p2xnn' to denote that they are perts, 'p', 
    #    to be compared with 2xco2 control, '2x' ('gregory_2xco2')
    fldp2xplt,sh=cutl.calc_totseaicearea(fldp2xclim,lat,lon)
    fldp2x12plt,sh=cutl.calc_totseaicearea(fldp2x12clim,lat,lon)
    fldp2x16plt,sh=cutl.calc_totseaicearea(fldp2x16clim,lat,lon)
    fldp2x19plt,sh=cutl.calc_totseaicearea(fldp2x19clim,lat,lon)
    fldp2x21plt,sh=cutl.calc_totseaicearea(fldp2x21clim,lat,lon)

    tstr = 'SIA (millions of km^2)'
    divby = np.float(1e12)

elif field=='sic':

    # Note that sea ice conc was read in and climo'd above,
    # but recalc'ing just to be safe. 
    
    fldc2xsclim,fldc2xsstd = cutl.climatologize(fldc2xs) # climo mean, gregory
    fldcsclim,fldcsstd = cutl.climatologize(fldcs) # climo mean, iga
    
    fldp2xsclim,fldp2xsstd = cutl.climatologize(fldp2xs) # climo mean, kel10    
    fldp2x12sclim,fldp2x12sstd = cutl.climatologize(fldp2x12s) # climo mean, kel12
    fldp2x16sclim,fldp2x16sstd = cutl.climatologize(fldp2x16s) # climo mean, kel16
    fldp2x19sclim,fldp2x19sstd = cutl.climatologize(fldp2x19s) # climo mean, kel19
    fldp2x21sclim,fldp2x21sstd = cutl.climatologize(fldp2x21s) # climo mean, kel21
    
    fldc2xplt,sh=cutl.calc_totseaicevol(fldc2xclim,fldc2xsclim,lat,lon) # 2xco2 control
    fldcplt,sh=cutl.calc_totseaicevol(fldcclim,fldcsclim,lat,lon) # preindustrial control

    # these perturbed runs are attempting to match the above prei control.
    # they are labelled 'p2xnn' to denote that they are perts, 'p', 
    #    to be compared with 2xco2 control, '2x' ('gregory_2xco2')
    fldp2xplt,sh=cutl.calc_totseaicevol(fldp2xclim,fldp2xsclim,lat,lon)
    fldp2x12plt,sh=cutl.calc_totseaicevol(fldp2x12clim,fldp2x12sclim,lat,lon)
    fldp2x16plt,sh=cutl.calc_totseaicevol(fldp2x16clim,fldp2x16sclim,lat,lon)
    fldp2x19plt,sh=cutl.calc_totseaicevol(fldp2x19clim,fldp2x19sclim,lat,lon)
    fldp2x21plt,sh=cutl.calc_totseaicevol(fldp2x21clim,fldp2x21sclim,lat,lon)

    # put field back 
    field='sic'; ncfield=field.upper(); conv=1/913.
    tstr = 'SIV (m^3)'
    divby=1

# these comparisions ask how good is the nudging to prei ice (in a 2xco2 climate)
diff10plt = fldcplt-fldp2xplt
diff12plt = fldcplt-fldp2x12plt
diff16plt = fldcplt-fldp2x16plt
diff19plt = fldcplt-fldp2x19plt
diff21plt = fldcplt-fldp2x21plt

# as a baseline: what is the difference b/w 2xco2 and preindustrial ice?
diff2xpiplt = fldc2xplt-fldcplt


fig,axs = plt.subplots(1,3)
fig.set_size_inches(12,4)
ax=axs[0]

ax.plot(mons,fldcplt/divby,color='.5',linewidth=2) # control prei is gray
ax.plot(mons,fldp2xplt/divby,'b--',linewidth=2) # blue, green, magenta, yellow,cyan, black
ax.plot(mons,fldp2x12plt/divby,'g--',linewidth=2)
ax.plot(mons,fldp2x16plt/divby,'m--',linewidth=2)
ax.plot(mons,fldp2x19plt/divby,'y--',linewidth=2)
ax.plot(mons,fldp2x21plt/divby,'c--',linewidth=2)
ax.set_xlim((1,12))
ax.set_title('NH ' + tstr)


ax = axs[1]
ax.plot(mons,diff10plt/divby,'b',linewidth=2)
ax.plot(mons,diff12plt/divby,'g',linewidth=2)
ax.plot(mons,diff16plt/divby,'m',linewidth=2)
ax.plot(mons,diff19plt/divby,'y',linewidth=2)
ax.plot(mons,diff21plt/divby,'c',linewidth=2)
ax.set_xlim((1,12))
ax.set_title('NH diff: iga-prei nudges')

ax=axs[2]
ax.plot(mons,diff2xpiplt/divby,color='.5',linewidth=3)
ax.set_xlim((1,12))
ax.set_title('NH diff: gregory-iga (for ref)')

if printtofile:
    fig.savefig(field + 'diffclimo_preinudge_seacyc_nh.pdf')



# <codecell>

# # Compare preI with the preI nudge

if showmaps:
    # kel10
    fig = cplt.map_allmonths(fldc-fldp2x,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='preI - preI nudge (' + casenamep2x + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec + '_v_' + casenamep2x + '_allmos_nh.pdf')

    # kel12
    fig = cplt.map_allmonths(fldc-fldp2x12,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='preI - preI nudge (' + casenamep2x12 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec + '_v_' + casenamep2x12 + '_allmos_nh.pdf')

    # kel16
    fig = cplt.map_allmonths(fldc-fldp2x16,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='preI - preI nudge (' + casenamep2x16 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec + '_v_' + casenamep2x16 + '_allmos_nh.pdf')

    # kel19
    fig = cplt.map_allmonths(fldc-fldp2x19,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='preI - preI nudge (' + casenamep2x19 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec + '_v_' + casenamep2x19 + '_allmos_nh.pdf')
    
    # kel21
    fig = cplt.map_allmonths(fldc-fldp2x21,lat,lon,type='nh',
                             cmap='red2blue_w20',cmin=cminn,cmax=cmaxn,lmask=1,
                             title='preI - preI nudge (' + casenamep2x21 + ')')
    if printtofile:
        fig.savefig(field + 'diff_' + casenamec + '_v_' + casenamep2x21 + '_allmos_nh.pdf')


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

