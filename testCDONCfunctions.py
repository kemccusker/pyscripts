# @@ test getNCvar() function with CDO bindings implemented
#
import numpy as np # for array handling
import scipy as sp # scientific python
import matplotlib.pyplot as plt # for basic plotting
from netCDF4 import Dataset
import cccmautils as cutl
import cccmaNC as cnc
import cccmacmaps
import platform as platform
import constants as con

plt.close('all')

cutl = reload(cutl)
cnc = reload(cnc)

model = 'CanAM4'
ftype = 'ts'    # timeseries
#ftype = 'climo'  # 12-month climatology
timeavg = 'ANN'
 
# # # ########### set Simulations #############
# Control run
casename = 'kemctl1'
timstr = '001-061'
timstr2 = '062-111'
#timstr = '001-111'
styr = 2             # skip year 1
enyr = 61 

# Pert run
casenamep1 = 'kem1pert1'  # 2002-2012 sic and sit
casenamep2 = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
casenamep3 = 'kem1pert3'  # 2002-2012 sic, adjusted sst. control sit
timstrp = '001-061'
#timstrp = '001-111'
styrp = 2             # skip year 1
enyrp = 61

#cmap = 'RdBu_r'
cmap = 'blue2red_20'
cmapclimo = 'Spectral_r'
conv = 1   # conversion factor to convert units, etc

# # # ######## set Field info ###############
field = 't'
ncfield = 'TEMP'
units = 'K' # @@
cminc = 200
cmaxc = 310
cmin = -.5
cmax = .5
cminm = -.8 # for monthly
cmaxm = .8  # for monthly
cminsc = -2.5 # as Screen et al paper
cmaxsc = 2.5  # as Screen et al paper



# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'

fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_' + ftype + '.nc'
fnamec2 = basepath + casename + subdir + casename + '_' + field + '_' + timstr2 + '_' + ftype + '.nc'



################################ VANILLA ######################################

## ncfilec = Dataset(fnamec,'r') # control
## #fldc = ncfilec.variables[ncfield][(styr-1)*12:(enyr*12+1),:,:,:]*conv
## fldc = ncfilec.variables[ncfield][...]*conv
## ncfilec.close()

## # super-most vanilla
## fldctmvan = np.average(fldc,0)
## fldctmzmvan = np.average(fldctmvan,2)

## # vanilla but remove extra lon
## fldctmvanNOLON = fldctmvan[:,:,0:-1]
## fldctmzmvanNOLON = np.average(fldctmvanNOLON,2)

## # vanilla but time-weight
## fldczmvan = np.average(fldc,3)
## fldczmvanann = cutl.annualize_monthlyts(fldczmvan)
## fldczmvananntm = np.average(fldczmvanann,0)

## # vanilla but remove extra lon and time-weight
## fldcNOLON = fldc[:,:,:,0:-1]
## fldczmvanNOLON = np.average(fldcNOLON,3)
## fldczmvanannNOLON = cutl.annualize_monthlyts(fldczmvanNOLON)
## fldczmvananntmNOLON = np.average(fldczmvanannNOLON,0)

################################# PY/NC/CDO ###################################
lat = cnc.getNCvar(fnamec,'lat')
lev = cnc.getNCvar(fnamec,'plev')
nlev = len(lev)
nlat = len(lat)


fldczmnew = np.append(cnc.getNCvar(fnamec,ncfield,calc='zm')*conv,
                      cnc.getNCvar(fnamec2,ncfield,calc='zm')*conv,
                      axis=0)
fldczmnewtm = np.mean(cutl.annualize_monthlyts(fldczmnew),0)
fldczmnewtm2= np.mean(fldczmnew,0) # not time-weighted

fldczmnewANN = np.append(cnc.getNCvar(fnamec,ncfield,calc='zm',seas='ANN')*conv,
                         cnc.getNCvar(fnamec2,ncfield,calc='zm',seas='ANN')*conv,
                         axis=0)
fldczmnewANNtm = np.mean(fldczmnewANN,0)

fldczmnewdjf = np.append(cnc.getNCvar(fnamec,ncfield,calc='zm',seas='DJF')*conv,
                         cnc.getNCvar(fnamec2,ncfield,calc='zm',seas='DJF')*conv,
                         axis=0)
fldczmnewdjftm = np.mean(fldczmnewdjf,0)

fldczmnewjja = np.append(cnc.getNCvar(fnamec,ncfield,calc='zm',seas='JJA')*conv,
                         cnc.getNCvar(fnamec2,ncfield,calc='zm',seas='JJA')*conv,
                         axis=0)
fldczmnewjjatm = np.mean(fldczmnewjja,0)


### test climatologize3d()

fldczmnewclimo,std = cutl.climatologize3d(fldczmnew)



################################# PY/NC ########################################
lat = cnc.getNCvar(fnamec,'lat')
lev = cnc.getNCvar(fnamec,'plev')

#fldczm = cnc.getNCvar_old(fnamec,ncfield,calc='zm',timechunk=((styr-1)*12,enyr*12+1) )*conv
fldczm = cnc.getNCvar_old(fnamec,ncfield,calc='zm')*conv
fldczmtm = np.mean(cutl.annualize_monthlyts(fldczm),0)
fldczmtm2 = np.average(fldczm,0) # not month-weighted


lats,levs = np.meshgrid(lat,lev)



plotfld = fldczmnewANNtm - fldczmnewtm
# check zm and seasonal avg done in cccmaNC module vs just zm done in module
# these are exactly the same! good.
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('AllModule - justZMmod (weighted later) ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()


plotfld = fldczmnewtm - fldczmnewtm2
# check just zm in module (weighted after) vs zm in module (not-weighted)
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('justZMmod (weighted later) - justZMmod (non-weighted) ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()





ftype = 'climo'

fnamec = basepath + casename + subdir + casename + '_' + field + '_001-111_' + ftype + '.nc'

#fldcclim = cnc.getNCvar(fnamec,ncfield)*conv
#print fldcclim.shape

fldczmclim = cnc.getNCvar(fnamec,ncfield,calc="zm")*conv
print fldczmclim.shape

# weight these
wgts = con.get_monweights() # 12
# reshape to 12 x lev x lat
wgts = np.tile(wgts,(nlev,nlat,1)) # make a weights array for each grid pt
# put the 3rd dim (2) in first spot, 1st dim (0) in second spot, 2nd dim (1) in third spot
# this effectively shifts the dimensions back to time,dim2,dim3
wgts = np.transpose(wgts,(2,0,1))
print wgts.shape

fldczmclimANN = np.average(fldczmclim,axis=0,weights=wgts) # weighted annual average

# test against: fldczmnewANNtm

plotfld = fldczmclimANN - fldczmnewANNtm
# the same, good.
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('from 12-mo climo (wgted) - from timeseries (annualized/wgted) ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()

wgtsseas = wgts[[0,1,11],:,:]
fldczmclimdjf = np.average(fldczmclim[[0,1,11],:,:],axis=0,weights=wgtsseas)

plotfld = fldczmclimdjf - fldczmnewdjftm
# DJF definitely shows some differences here. It seems that is is because
# the seasonalize() method keeps the "season" together (ie shortens the DJF
# timeseries by a yr compared with annual timeseries) whereas the CDO climo
# does not worry about this.
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('from 12-mo climo DJF (wgted) - from timeseries DJF (seasonalized/wgted) ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()

wgtsseas = wgts[[5,6,7],:,:]
fldczmclimjja = np.average(fldczmclim[[5,6,7],:,:],axis=0,weights=wgtsseas)

plotfld = fldczmclimjja - fldczmnewjjatm
# the same, good.
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('from 12-mo climo JJA (wgted) - from timeseries JJA (seasonalized/wgted) ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()



# Now test my climatologize3d() method (fldczmnewclimo)

plotfld = fldczmclim[11,:,:] - fldczmnewclimo[11,:,:]
# these all look good (the comparison is the same)
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('CDO climo zm - python climo zm ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()



"""
# checking now whether extra longitude makes a noticeable difference. it does.
plotfld = fldctmzmvan - fldctmzmvanNOLON

fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('Most vanilla - Most vanilla NOLON' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()


plotfld = fldczmvananntmNOLON - fldctmzmvan

fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('Most vanilla - Most vanilla (tim-wgted NOLON)' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()

"""
