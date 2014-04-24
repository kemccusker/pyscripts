# plot_canam4sims2.py
#    2/5/2014
#    start looking at the CanAM4 sea ice simulations
#
#    2/11/2014: try calling my own functions
#    

import numpy as np # for array handling
import scipy as sp # scientific python
import matplotlib.pyplot as plt # for basic plotting
import matplotlib.cm as cm
from subprocess import call # for doing system calls - not really needed
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap # for maps
import datetime as datetime
import matplotlib.colors as col
import platform as platform
import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import matplotlib.font_manager as fm


# while I'm still creating these modules, have to reload to get changes
cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)

plt.close("all")
plt.ion()

printtofile=1

model = 'CanAM4'

# # # ########### set Simulations #############
# Control run
casename = 'kemctl1'
timstr = '001-061'
styr = 2             # skip year 1
enyr = 61 

# Pert run
casenamep1 = 'kem1pert1'  # 2002-2012 sic and sit
casenamep2 = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
casenamep3 = 'kem1pert3'  # 2002-2012 sic, adjusted sst. control sit
timstrp = '001-061'  
styrp = 2             # skip year 1
enyrp = 61

cmap = ''
cmapclimo = 'Spectral_r'

# # # ######## set Field info ###############
## field = 'gt'
## units = 'C'
## conv = 1  # no conversion
## cmin = -2 # for anomaly plots
## cmax = 2  # for anomaly plots
## cminstd = -0.3 # for std anomaly plots
## cmaxstd = 0.3  # for std anomaly plots
## anngmymin = 14
## anngmymax = 14.5
## monstdymin = -0.8
## monstdymax = 0.8

## field = 'st'
## anngmymin = 13.7
## anngmymax = 14.1
## annpmymin = -11
## annpmymax = -8.5
## N=10 # for colormap when I get them working@@

## field = 'pmsl'
## units = 'hPa' # pretty sure hpa @@double check
## conv = 1
## cmin = -.5 # for anomaly plots
## cmax = .5  # for anomaly plots
## cminstd = -0.5 # for std anomaly plots
## cmaxstd = 0.5  # for std anomaly plots
## cminm=-2 # for monthly maps
## cmaxm=2  # for monthly maps
## anngmymin = 1011.37 # annual mean timeseries
## anngmymax = 1011.45
## annpmymin = 1009.0
## annpmymax = 1015.0
## monstdymin = -0.1  # monthly mean standardized timeseries
## monstdymax = 0.1


field = 'pcp'
units = 'mm/day' # original: kg m-2 s-1
conv = 86400  # convert from kg m-2 s-1 to mm/day
cmin = -.4 # for anomaly plots
cmax = .4  # for anomaly plots
cminstd = -0.2 # for std anomaly plots
cmaxstd = 0.2  # for std anomaly plots
cminclimo=0
cmaxclimo=15
cminm=-1 # for monthly maps
cmaxm=1  # for monthly maps
anngmymin = 2.6 # annual mean timeseries
anngmymax = 3
annpmymin = 1
annpmymax = 1.5
monstdymin = -0.1  # monthly mean standardized timeseries
monstdymax = 0.1
cmap = 'PuOr'
cmapclimo = 'YlGnBu'

# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'

fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
fnamep1 = basepath + casenamep1 + subdir + casenamep1 + '_' + field + '_' + timstrp + '_ts.nc'
fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp + '_ts.nc'
fnamep3 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp + '_ts.nc'

# Get the data
ncfilec = Dataset(fnamec,'r') # control
ncfilep1 = Dataset(fnamep1,'r') # pert1
ncfilep2 = Dataset(fnamep2,'r') # pert2
ncfilep3 = Dataset(fnamep3,'r') # pert3

lat = ncfilec.variables['lat'][:]
lon = ncfilec.variables['lon'][:]

fldc = ncfilec.variables[field.upper()][(styr-1)*12:(enyr*12+1),:,:]*conv # time start year to end
fldp1 = ncfilep1.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
fldp2 = ncfilep2.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
fldp3 = ncfilep3.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end


# # # ##################### Do calculations #################
# annual time-series (3d)
anntsc = cutl.annualize_monthlyts(fldc)
anntsp1 = cutl.annualize_monthlyts(fldp1)
anntsp2 = cutl.annualize_monthlyts(fldp2)
anntsp3 = cutl.annualize_monthlyts(fldp3)

# annual global mean time-series
anngmc = cutl.global_mean_areawgted3d(anntsc,lat,lon)
anngmp1 = cutl.global_mean_areawgted3d(anntsp1,lat,lon)
anngmp2 = cutl.global_mean_areawgted3d(anntsp2,lat,lon)
anngmp3 = cutl.global_mean_areawgted3d(anntsp3,lat,lon)

# annual polar (>=60N) mean time-series
annpmc = cutl.polar_mean_areawgted3d(anntsc,lat,lon)
annpmp1 = cutl.polar_mean_areawgted3d(anntsp1,lat,lon)
annpmp2 = cutl.polar_mean_areawgted3d(anntsp2,lat,lon)
annpmp3 = cutl.polar_mean_areawgted3d(anntsp3,lat,lon)

# global mean annual mean time mean
anngmtmc = np.mean(anngmc)
anngmtmp1 = np.mean(anngmp1)
anngmtmp2 = np.mean(anngmp2)
anngmtmp3 = np.mean(anngmp3)


# monthly global mean time-series
gmc = cutl.global_mean_areawgted3d(fldc,lat,lon)
gmp1 = cutl.global_mean_areawgted3d(fldp1,lat,lon)
gmp2 = cutl.global_mean_areawgted3d(fldp2,lat,lon)
gmp3 = cutl.global_mean_areawgted3d(fldp3,lat,lon)

# monthly polar (>=60N) mean time-series
pmc = cutl.polar_mean_areawgted3d(fldc,lat,lon)
pmp1 = cutl.polar_mean_areawgted3d(fldp1,lat,lon)
pmp2 = cutl.polar_mean_areawgted3d(fldp2,lat,lon)
pmp3 = cutl.polar_mean_areawgted3d(fldp3,lat,lon)

# monthly climatology and standard dev (3d)
(climoc, stdc) = cutl.climatologize3d(fldc)
(climop1,stdp1) = cutl.climatologize3d(fldp1)
(climop2,stdp2) = cutl.climatologize3d(fldp2)
(climop3,stdp3) = cutl.climatologize3d(fldp3)

# global mean monthly climatology
climogmc = cutl.global_mean_areawgted3d(climoc,lat,lon)
climogmp1 = cutl.global_mean_areawgted3d(climop1,lat,lon)
climogmp2 = cutl.global_mean_areawgted3d(climop2,lat,lon)
climogmp3 = cutl.global_mean_areawgted3d(climop3,lat,lon)

# polar (>=60N) mean monthly climatology
climopmc = cutl.polar_mean_areawgted3d(climoc,lat,lon)
climopmp1 = cutl.polar_mean_areawgted3d(climop1,lat,lon)
climopmp2 = cutl.polar_mean_areawgted3d(climop2,lat,lon)
climopmp3 = cutl.polar_mean_areawgted3d(climop3,lat,lon)


# global mean monthly standard deviation (in time)
stdgmc = cutl.global_mean_areawgted3d(stdc,lat,lon)
stdgmp1 = cutl.global_mean_areawgted3d(stdp1,lat,lon)
stdgmp2 = cutl.global_mean_areawgted3d(stdp2,lat,lon)
stdgmp3 = cutl.global_mean_areawgted3d(stdp3,lat,lon)

# polar (>=60N) mean monthly standard deviation (in time)
stdpmc = cutl.polar_mean_areawgted3d(stdc,lat,lon)
stdpmp1 = cutl.polar_mean_areawgted3d(stdp1,lat,lon)
stdpmp2 = cutl.polar_mean_areawgted3d(stdp2,lat,lon)
stdpmp3 = cutl.polar_mean_areawgted3d(stdp3,lat,lon)


# remove seasonal cycle from monthly timeseries
#    repeat climo for each year of timeseries
climogmcts = np.squeeze(np.tile(climogmc,(1,fldc.shape[0]/12)))
climogmp1ts = np.squeeze(np.tile(climogmp1,(1,fldp1.shape[0]/12)))
climogmp2ts = np.squeeze(np.tile(climogmp2,(1,fldp2.shape[0]/12)))
climogmp3ts = np.squeeze(np.tile(climogmp3,(1,fldp3.shape[0]/12)))
gmcrem = gmc-climogmcts
gmp1rem = gmp1-climogmp1ts
gmp2rem = gmc-climogmp2ts
gmp3rem = gmc-climogmp3ts

#    POLAR (>=60N)
climopmcts = np.squeeze(np.tile(climopmc,(1,fldc.shape[0]/12)))
climopmp1ts = np.squeeze(np.tile(climopmp1,(1,fldp1.shape[0]/12)))
climopmp2ts = np.squeeze(np.tile(climopmp2,(1,fldp2.shape[0]/12)))
climopmp3ts = np.squeeze(np.tile(climopmp3,(1,fldp3.shape[0]/12)))
pmcrem = pmc-climopmcts
pmp1rem = pmp1-climopmp1ts
pmp2rem = pmc-climopmp2ts
pmp3rem = pmc-climopmp3ts

# normalize each month's anomaly by that month's sigma
stdgmcts = np.squeeze(np.tile(stdgmc,(1,fldc.shape[0]/12)))
stdgmp1ts = np.squeeze(np.tile(stdgmp1,(1,fldp1.shape[0]/12)))
stdgmp2ts = np.squeeze(np.tile(stdgmp2,(1,fldp2.shape[0]/12)))
stdgmp3ts = np.squeeze(np.tile(stdgmp3,(1,fldp3.shape[0]/12)))
gmcremstd = gmcrem/stdgmcts
gmp1remstd = gmp1rem/stdgmp1ts
gmp2remstd = gmp2rem/stdgmp2ts
gmp3remstd = gmp3rem/stdgmp3ts

#    POLAR (>=60N)
stdpmcts = np.squeeze(np.tile(stdpmc,(1,fldc.shape[0]/12)))
stdpmp1ts = np.squeeze(np.tile(stdpmp1,(1,fldp1.shape[0]/12)))
stdpmp2ts = np.squeeze(np.tile(stdpmp2,(1,fldp2.shape[0]/12)))
stdpmp3ts = np.squeeze(np.tile(stdpmp3,(1,fldp3.shape[0]/12)))
pmcremstd = pmcrem/stdpmcts
pmp1remstd = pmp1rem/stdpmp1ts
pmp2remstd = pmp2rem/stdpmp2ts
pmp3remstd = pmp3rem/stdpmp3ts


# # # other things to calculate @@
#     land vs ocean mean
#     pdf
#     standard dev





# # # ################## Make Figures #########################

#   Maps

# control ANN climo (whole timeseries)
#    GLOBE

plotfld = np.mean(anntsc,axis=0)

plt.figure()
cplt.kemmap(plotfld,lat,lon,type='sq',cmap=cmapclimo,\
            cmin=cminclimo,cmax=cmaxclimo,title=field + " ANN mean " + casename)
if printtofile:
    plt.savefig(field + '_ANNclimo_' + casename + '.pdf')
    plt.savefig(field + '_ANNclimo_' + casename + '.png')

#    POLE
plt.figure()
cplt.kemmap(plotfld,lat,lon,type='nh',cmap=cmapclimo,\
            cmin=cminclimo,cmax=cmaxclimo,title=field + " ANN mean " + casename)
if printtofile:
    plt.savefig(field + '_ANNclimo_' + casename + '_nh.pdf')
    plt.savefig(field + '_ANNclimo_' + casename + '_nh.png')

#   annual mean diffs, all sims
fig, ax = plt.subplots(1,3)
fig.set_size_inches(12,4)
fig.subplots_adjust(hspace=0.1)

plotfld =  np.mean(anntsp1-anntsc,axis=0)
plt.subplot(1,3,1)
cplt.kemmap(plotfld,lat,lon,type='nh',cmin=cmin,cmax=cmax,title=casenamep1,cmap=cmap)
plotfld =  np.mean(anntsp2-anntsc,axis=0)
plt.subplot(1,3,2)
cplt.kemmap(plotfld,lat,lon,type='nh',cmin=cmin,cmax=cmax,title=casenamep2,cmap=cmap)
plotfld =  np.mean(anntsp3-anntsc,axis=0)
plt.subplot(1,3,3)
cplt.kemmap(plotfld,lat,lon,type='nh',cmin=cmin,cmax=cmax,title=casenamep3,cmap=cmap)

#cbar_ax = fig.add_axes([.95,.15,.02,.7])
#fig.colorbar(pc,cax=cbar_ax) # @@ pc is pcolormesh which is inside the kemmap function
fig.suptitle('annual mean ' + field + ' diffs')

if printtofile:
    plt.savefig(field + '_ANNdiff_allperts_nh.pdf')
    plt.savefig(field + '_ANNdiff_allperts_nh.png')

#   annual sigma diffs, all sims
fig, ax = plt.subplots(1,3)
fig.set_size_inches(12,4)
fig.subplots_adjust(hspace=0.1)

plotfld =  np.std(anntsp1,axis=0)-np.std(anntsc,axis=0)
plt.subplot(1,3,1)
cplt.kemmap(plotfld,lat,lon,type='nh',cmin=cminstd,cmax=cmaxstd,title=casenamep1,cmap=cmap)
plotfld =  np.std(anntsp2,axis=0)-np.std(anntsc,axis=0)
plt.subplot(1,3,2)
cplt.kemmap(plotfld,lat,lon,type='nh',cmin=cminstd,cmax=cmaxstd,title=casenamep2,cmap=cmap)
plotfld =  np.std(anntsp3,axis=0)-np.std(anntsc,axis=0)
plt.subplot(1,3,3)
cplt.kemmap(plotfld,lat,lon,type='nh',cmin=cminstd,cmax=cmaxstd,title=casenamep3,cmap=cmap)

#cbar_ax = fig.add_axes([.95,.15,.02,.7])
#fig.colorbar(pc,cax=cbar_ax) # @@ pc is pcolormesh which is inside the kemmap function
fig.suptitle('annual mean ' + field + ' diffs in standard dev')

if printtofile:
    plt.savefig(field + 'std_ANNdiff_allperts_nh.pdf')
    plt.savefig(field + 'std_ANNdiff_allperts_nh.png')

#   ALL MONTHS maps
#printtofile=1
fig, ax = plt.subplots(2,6)
fig.set_size_inches(12,6)
fig.subplots_adjust(hspace=0.1)

months=con.get_mon()
for midx in range(len(months)):
    plotfld = climop1[midx,:,:] - climoc[midx,:,:]

    plt.subplot(2,6,midx+1)
    cplt.kemmap(plotfld,lat,lon,type='nh',cmin=cminm,cmax=cmaxm,title=months[midx],cmap=cmap)

if printtofile:
     plt.savefig(field + '_diff_' + casenamep1 + '_v_' + casename + 'allmonths_nh.pdf')
     plt.savefig(field + '_diff_' + casenamep1 + '_v_' + casename + 'allmonths_nh.png')

#   ALL MONTHS maps
fig, ax = plt.subplots(2,6)
fig.set_size_inches(12,6)
fig.subplots_adjust(hspace=0.1)

for midx in range(len(months)):
    plotfld = climop2[midx,:,:] - climoc[midx,:,:]

    plt.subplot(2,6,midx+1)
    cplt.kemmap(plotfld,lat,lon,type='nh',cmin=cminm,cmax=cmaxm,title=months[midx],cmap=cmap)

if printtofile:
     plt.savefig(field + '_diff_' + casenamep2 + '_v_' + casename + 'allmonths_nh.pdf')
     plt.savefig(field + '_diff_' + casenamep2 + '_v_' + casename + 'allmonths_nh.png')


#   ALL MONTHS maps
fig, ax = plt.subplots(2,6)
fig.set_size_inches(12,6)
fig.subplots_adjust(hspace=0.1)

for midx in range(len(months)):
    plotfld = climop3[midx,:,:] - climoc[midx,:,:]

    plt.subplot(2,6,midx+1)
    cplt.kemmap(plotfld,lat,lon,type='nh',cmin=cminm,cmax=cmaxm,title=months[midx],cmap=cmap)

if printtofile:
     plt.savefig(field + '_diff_' + casenamep3 + '_v_' + casename + 'allmonths_nh.pdf')
     plt.savefig(field + '_diff_' + casenamep3 + '_v_' + casename + 'allmonths_nh.png')


##################################
#   Line plots: timeseries
#printtofile=0
fontP = fm.FontProperties()
fontP.set_size('xx-small')

years = np.arange(styr,enyr+1)
plt.figure()
plt.plot(years,anngmc,'k')
plt.plot(years,anngmp1,'g')
plt.plot(years,anngmp2,'r')
plt.plot(years,anngmp3,'b')
plt.title('annual, global mean ' + field)
plt.ylim(anngmymin,anngmymax)
plt.xlim(styr,enyr)
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
if printtofile:
    plt.savefig(field + '_ANN_timeseries_allsims.pdf')
    plt.savefig(field + '_ANN_timeseries_allsims.png')

#     POLAR
plt.figure()
plt.plot(years,annpmc,'k')
plt.plot(years,annpmp1,'g')
plt.plot(years,annpmp2,'r')
plt.plot(years,annpmp3,'b')
plt.title('annual, polar mean (>=60N) ' + field)
plt.ylim(annpmymin,annpmymax)
plt.xlim(styr,enyr)
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
if printtofile:
    plt.savefig(field + '_ANN_timeseries_allsims_np.pdf')
    plt.savefig(field + '_ANN_timeseries_allsims_np.png')

    
fig = plt.figure()
fig.set_size_inches(12,4)
plt.plot(gmcrem,'k')
plt.plot(gmp1rem,'g')
plt.plot(gmp2rem,'r')
plt.plot(gmp3rem,'b')
plt.title('monthly ' + field + ' (monthly climo removed)')
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
if printtofile:
    plt.savefig(field + '_monthly_timeseries_allsims_climorem.pdf')
    plt.savefig(field + '_monthly_timeseries_allsims_climorem.png')

#    POLAR
fig = plt.figure()
fig.set_size_inches(12,4)
plt.plot(pmcrem,'k')
plt.plot(pmp1rem,'g')
plt.plot(pmp2rem,'r')
plt.plot(pmp3rem,'b')
plt.title('Polar (>=60N) mean monthly ' + field + ' (monthly climo removed)')
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
if printtofile:
    plt.savefig(field + '_monthly_timeseries_allsims_climorem_np.pdf')
    plt.savefig(field + '_monthly_timeseries_allsims_climorem_np.png')

fig = plt.figure()
fig.set_size_inches(12,4)
plt.plot(gmcremstd,'k')
plt.plot(gmp1remstd,'g')
plt.plot(gmp2remstd,'r')
plt.plot(gmp3remstd,'b')
plt.title('standardized monthly ' + field + ' (monthly climo removed)')
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
plt.ylim(monstdymin,monstdymax)
if printtofile:
    plt.savefig(field + '_monthly_timeseries_allsims_stdized.pdf')
    plt.savefig(field + '_monthly_timeseries_allsims_stdized.png')

#    POLAR
fig = plt.figure()
fig.set_size_inches(12,4)
plt.plot(pmcremstd,'k')
plt.plot(pmp1remstd,'g')
plt.plot(pmp2remstd,'r')
plt.plot(pmp3remstd,'b')
plt.title('Polar (>=60N) mean standardized monthly ' + field + ' (monthly climo removed)')
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
# @@plt.ylim(monstdymin,monstdymax)
if printtofile:
    plt.savefig(field + '_monthly_timeseries_allsims_stdized_np.pdf')
    plt.savefig(field + '_monthly_timeseries_allsims_stdized_np.png')


############################################
#   Line plots: seasonal cycle
mons = range(1,13)
plt.figure()
plt.plot(mons,climogmc,'k')
plt.plot(mons,climogmp1,'g')
plt.plot(mons,climogmp2,'r')
plt.plot(mons,climogmp3,'b')
plt.title(field + ' global mean seasonal climo')
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
plt.xlim(1,12)
if printtofile:
    plt.savefig(field + '_seacycle_allsims.pdf')
    plt.savefig(field + '_seacycle_allsims.png')

#    POLAR
plt.figure()
plt.plot(mons,climopmc,'k')
plt.plot(mons,climopmp1,'g')
plt.plot(mons,climopmp2,'r')
plt.plot(mons,climopmp3,'b')
plt.title(field + ' polar (>=60N) mean seasonal climo')
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
plt.xlim(1,12)
if printtofile:
    plt.savefig(field + '_seacycle_allsims_np.pdf')
    plt.savefig(field + '_seacycle_allsims_np.png')

plt.figure()
plt.plot(mons,climogmp1-climogmc,'g')
plt.plot(mons,climogmp2-climogmc,'r')
plt.plot(mons,climogmp3-climogmc,'b')
plt.title(field + ' global mean seasonal diff from ' + casename)
plt.legend((casenamep1,casenamep2,casenamep3),loc=9,prop=fontP)
plt.xlim(1,12)
if printtofile:
    plt.savefig(field + '_seacyclediff_allpert.pdf')
    plt.savefig(field + '_seacyclediff_allpert.png')

#   POLAR
plt.figure()
plt.plot(mons,climopmp1-climopmc,'g')
plt.plot(mons,climopmp2-climopmc,'r')
plt.plot(mons,climopmp3-climopmc,'b')
plt.title(field + ' polar (>=60N) mean seasonal diff from ' + casename)
plt.legend((casenamep1,casenamep2,casenamep3),loc=9,prop=fontP)
plt.xlim(1,12)
if printtofile:
    plt.savefig(field + '_seacyclediff_allpert_np.pdf')
    plt.savefig(field + '_seacyclediff_allpert_np.png')
    

plt.figure()
plt.plot(mons,stdgmc,'k')
plt.plot(mons,stdgmp1,'g')
plt.plot(mons,stdgmp2,'r')
plt.plot(mons,stdgmp3,'b')
plt.title(field + ' global mean seasonal standard dev')
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
plt.xlim(1,12)
if printtofile:
    plt.savefig(field + 'std_seacycle_allsims.pdf')
    plt.savefig(field + 'std_seacycle_allsims.png')

#    POLAR
plt.figure()
plt.plot(mons,stdpmc,'k')
plt.plot(mons,stdpmp1,'g')
plt.plot(mons,stdpmp2,'r')
plt.plot(mons,stdpmp3,'b')
plt.title(field + ' polar (>=60N) mean seasonal standard dev')
plt.legend((casename,casenamep1,casenamep2,casenamep3),prop=fontP)
plt.xlim(1,12)
if printtofile:
    plt.savefig(field + 'std_seacycle_allsims_np.pdf')
    plt.savefig(field + 'std_seacycle_allsims_np.png')


plt.figure()
plt.plot(mons,stdgmp1-stdgmc,'g')
plt.plot(mons,stdgmp2-stdgmc,'r')
plt.plot(mons,stdgmp3-stdgmc,'b')
plt.title(field + ' global mean seasonal standard dev')
plt.legend((casenamep1,casenamep2,casenamep3),prop=fontP)
plt.xlim(1,12)
if printtofile:
    plt.savefig(field + 'std_seacyclediff_allpert.pdf')
    plt.savefig(field + 'std_seacyclediff_allpert.png')

#    POLAR
plt.figure()
plt.plot(mons,stdpmp1-stdpmc,'g')
plt.plot(mons,stdpmp2-stdpmc,'r')
plt.plot(mons,stdpmp3-stdpmc,'b')
plt.title(field + ' polar (>=60N) mean seasonal standard dev')
plt.legend((casenamep1,casenamep2,casenamep3),prop=fontP)
plt.xlim(1,12)
if printtofile:
    plt.savefig(field + 'std_seacyclediff_allpert_np.pdf')
    plt.savefig(field + 'std_seacyclediff_allpert_np.png')

