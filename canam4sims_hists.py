# canam4sims_stats.py
#    2/14/2014: taken from plot_canam4sims2.py: 
#               calculate & plot statistical properties of the runs
#

import numpy as np # for array handling
import scipy as sp # scientific python
import scipy.stats
import matplotlib.pyplot as plt # for basic plotting
import matplotlib.cm as cm
from subprocess import call # for doing system calls - not really needed
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap # for maps
import datetime as datetime
import matplotlib.colors as col
import platform as platform
import cccmaplots as cplt    # my module
import constants as con      # my module
import cccmautils as cutl    # my module
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
timstr2 = '062-111'
styr = 2             # skip year 1
enyr = 61 

# Pert run
casenamep1 = 'kem1pert1'  # 2002-2012 sic and sit
casenamep2 = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
casenamep3 = 'kem1pert3'  # 2002-2012 sic, adjusted sst. control sit
timstrp = '001-061'
timstrp2 = '062-111'
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

## field = 'st'
## xlims = -4,4
## xlimsn = -5,5

field = 'pmsl'
units = 'hPa' # pretty sure hpa @@double check
conv = 1
cmin = -.5 # for anomaly plots
cmax = .5  # for anomaly plots
cminstd = -0.5 # for std anomaly plots
cmaxstd = 0.5  # for std anomaly plots
cminm=-2 # for monthly maps
cmaxm=2  # for monthly maps


## field = 'pcp'
## units = 'mm/day' # original: kg m-2 s-1
## conv = 86400  # convert from kg m-2 s-1 to mm/day
## cmin = -.4 # for anomaly plots
## cmax = .4  # for anomaly plots
## cminstd = -0.2 # for std anomaly plots
## cmaxstd = 0.2  # for std anomaly plots


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
fnamec2 = basepath + casename + subdir + casename + '_' + field + '_' + timstr2 + '_ts.nc'
fnamep12 = basepath + casenamep1 + subdir + casenamep1 + '_' + field + '_' + timstrp2 + '_ts.nc'
fnamep22 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp2 + '_ts.nc'
fnamep32 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp2 + '_ts.nc'


# Get the data
ncfilec = Dataset(fnamec,'r') # control
ncfilep1 = Dataset(fnamep1,'r') # pert1
ncfilep2 = Dataset(fnamep2,'r') # pert2
ncfilep3 = Dataset(fnamep3,'r') # pert3
#   time period 2
ncfilec2 = Dataset(fnamec2,'r') # control
ncfilep12 = Dataset(fnamep12,'r') # pert1
ncfilep22 = Dataset(fnamep22,'r') # pert2
ncfilep32 = Dataset(fnamep32,'r') # pert3

lat = ncfilec.variables['lat'][:]
lon = ncfilec.variables['lon'][:]

fldc = ncfilec.variables[field.upper()][(styr-1)*12:(enyr*12+1),:,:]*conv # time start year to end
fldc = np.append(fldc, ncfilec2.variables[field.upper()][...]*conv, axis=0)

fldp1 = ncfilep1.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
fldp1 = np.append(fldp1,ncfilep12.variables[field.upper()][...]*conv, axis=0)

fldp2 = ncfilep2.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
fldp2 = np.append(fldp2,ncfilep22.variables[field.upper()][...]*conv, axis=0)

fldp3 = ncfilep3.variables[field.upper()][(styrp-1)*12:(enyrp*12+1),:,:]*conv # time start year to end
fldp3 = np.append(fldp3,ncfilep32.variables[field.upper()][...]*conv, axis=0)


# annual time-series (3d)
anntsc = cutl.annualize_monthlyts(fldc)
anntsp1 = cutl.annualize_monthlyts(fldp1)
anntsp2 = cutl.annualize_monthlyts(fldp2)
anntsp3 = cutl.annualize_monthlyts(fldp3)


# calc histograms (or pdfs...is there a function?): scipy.stats
#   http://docs.scipy.org/doc/scipy/reference/tutorial/stats.html#t-test-and-ks-test
# calc ttest / statistical significance
# calc minimum ensemble size: how does this vary with N
# calc "time to statistical significance" ...same thing?

casenamep = casenamep2
timeavg = 'ANN'
plotfld = anntsp2 - anntsc
cellwgts = cutl.get_cellwgts(lat,lon,repeat=plotfld.shape)
plotfld = plotfld.flatten()

# Histogram: using pyplot hist()
## fig1 = plt.figure()
## plt.hist(plotfld,bins=100,normed=True) # cumulative=True
## # example of transparency
## # plt.hist(uniform_numbers, bins=20, histtype='stepfilled', \
## #       normed=True, color='r', alpha=0.5, label='Uniform')
## #plt.title("Gaussian Histogram")
## plt.xlabel(field)
## plt.ylabel("Frequency")
## plt.show()

# Histogram: using Numpy and bar plot
location = "global"

(hist,binedges) = np.histogram(plotfld,bins=100,weights=cellwgts.flatten(),density=True)
widths = np.diff(binedges)

fig2 = plt.figure()
plt.bar(binedges[0:-1],hist/np.sum(hist)*100,width=widths,color='0.75',alpha=0.5)
plt.xlabel(field + " (" + units + ")")
plt.ylabel("% of total")
plt.title(timeavg + " " + casenamep + " - " + casename + " (" + location + ")")
#plt.xlim(xlims)
# @@ should put all perts on here and print to file...


# North Pole (>=60N)
northbnd = 60
location = 'gt' + str(northbnd) + 'N'

plotfld = anntsp2 - anntsc
cellwgts = cutl.get_cellwgts(lat,lon,repeat=plotfld.shape)
plotfld = plotfld[:,lat>=northbnd,:]
cellwgts = cellwgts[:,lat>=northbnd,:]# can't do this b/c weights won't add to 1
plotfld = plotfld.flatten()

(histc,binedgesc) = np.histogram(anntsc[:,lat>=northbnd,:].flatten(),bins=100,density=True)
(histp,binedgesp) = np.histogram(anntsp1[:,lat>=northbnd,:].flatten(),bins=100,density=True)

plt.figure()
plt.bar(binedgesc[0:-1],histc / np.sum(histc) *100,width=np.diff(binedgesc),color='b',alpha=0.5)
plt.bar(binedgesp[0:-1],histp / np.sum(histp) *100,width=np.diff(binedgesp),color='0.75',alpha=0.5)
plt.xlabel(field + " (" + units + ")")
plt.ylabel("% of total")
plt.title(timeavg + " (" + location + ")")
plt.legend((casename,casenamep))
#plt.xlim(xlims)
if printtofile:
    plt.savefig(field + 'PDF_' + timeavg + '_' + location + '_' + casenamep + '_' + casename + '.pdf')
    plt.savefig(field + 'PDF_' + timeavg + '_' + location + '_' + casenamep + '_' + casename + '.png')



histccum = np.zeros(histc.shape)
histpcum = np.zeros(histp.shape)
for el in np.arange(1,len(histc)+1):
    histccum[el-1] = np.sum(histc[0:el])
    histpcum[el-1] = np.sum(histp[0:el])

plt.figure()
plt.plot(binedgesc[0:-1],histccum/np.sum(histc),'k')
plt.plot(binedgesp[0:-1],histpcum/np.sum(histp),color='.5')
plt.xlabel(field + " (" + units + ")")
plt.ylabel("Fraction of total")
plt.title("CDF: " + timeavg + " (" + location + ")")
plt.legend((casename,casenamep),loc='lower right')
if printtofile:
    plt.savefig(field + 'CDF_' + timeavg + '_' + location + '_' + casenamep + '_' + casename + '.pdf')
    plt.savefig(field + 'CDF_' + timeavg + '_' + location + '_' + casenamep + '_' + casename + '.png')

plt.figure()
plt.plot(binedgesc[0:-1],((histpcum/np.sum(histp))-(histccum/np.sum(histc)))*100,'k')
plt.xlabel(field + " (" + units + ")")
plt.ylabel("difference in CDF frac of total (%)")
plt.title("CDF: " + timeavg + " " + location  + " (" + casenamep + '-' + casename + ")")
if printtofile:
    plt.savefig(field + 'CDFdiff_' + timeavg + '_' + location + '_' + casenamep + '_v_' + casename + '.pdf')
    plt.savefig(field + 'CDFdiff_' + timeavg + '_' + location + '_' + casenamep + '_v_' + casename + '.png')


# ######### Put all perts on one hist
plotfld2 = anntsp2 - anntsc
cellwgts = cutl.get_cellwgts(lat,lon,repeat=plotfld2.shape)
plotfld2 = plotfld2[:,lat>=northbnd,:]
cellwgts = cellwgts[:,lat>=northbnd,:]
plotfld2 = plotfld2.flatten()

plotfld1 = anntsp1 - anntsc
plotfld1 = plotfld1[:,lat>=northbnd,:].flatten()
plotfld3 = anntsp3 - anntsc
plotfld3 = plotfld3[:,lat>=northbnd,:].flatten()


(hist2,binedges2) = np.histogram(plotfld2,bins=100,weights=cellwgts.flatten(),density=True)
widths2 = np.diff(binedges2)
(hist1,binedges1) = np.histogram(plotfld1,bins=100,weights=cellwgts.flatten(),density=True)
widths1 = np.diff(binedges1)
(hist3,binedges3) = np.histogram(plotfld3,bins=100,weights=cellwgts.flatten(),density=True)
widths3 = np.diff(binedges3)

fig4 = plt.figure()
plt.bar(binedges1[0:-1],hist1/np.sum(hist1)*100,width=widths1,color='b',alpha=0.5,ec='none')
plt.bar(binedges2[0:-1],hist2/np.sum(hist2)*100,width=widths2,color='0.3',alpha=0.5,ec='none')
plt.bar(binedges3[0:-1],hist3/np.sum(hist3)*100,width=widths3,color='0.7',alpha=0.5,ec='none')
plt.xlabel(field + " (" + units + ")")
plt.ylabel("% of total")
plt.title(timeavg + " All Perts - " + casename + " (" + location + ")")
plt.legend((casenamep1,casenamep2,casenamep3))
#plt.xlim(xlimsn)
plt.show()
if printtofile:
    plt.savefig(field + 'PDF_' + timeavg + '_' + location + '_allsims.pdf')
    plt.savefig(field + 'PDF_' + timeavg + '_' + location + '_allsims.png')



## # Kernel-Density Estimation (KDE)
## #     @@ how to area-weight?
## kernc = sp.stats.gaussian_kde(plotfld)

## fig = plt.figure()
## ax = fig.add_subplot(111)

## ax.plot(plotfld, np.zeros(plotfld.shape), 'b+', ms=20)  # rug plot
## x_eval = np.linspace(np.min(plotfld),np.max(plotfld), num=200)
## ax.plot(x_eval, kernc(x_eval), 'k-', label="Scott's Rule")
## plt.title("Note, unweighted!")
## plt.show()


