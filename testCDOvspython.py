# @@@@ did not get far with this. continue tomorrow. use new NC reading funcs
#
import numpy as np # for array handling
import scipy as sp # scientific python
import matplotlib.pyplot as plt # for basic plotting
from netCDF4 import Dataset
import cccmautils as cutl
import cccmaNC as cnc
import cdo as cdo; cdo = cdo.Cdo()
import os
import cccmacmaps

# cdo bindings leave a lot of temp files
os.system('rm -rf /tmp/cdoPy*')
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
fnamep1 = basepath + casenamep1 + subdir + casenamep1 + '_' + field + '_' + timstrp + '_' + ftype + '.nc'
fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp + '_' + ftype + '.nc'
fnamep3 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp + '_' + ftype + '.nc'

################################ VANILLA ######################################

ncfilec = Dataset(fnamec,'r') # control
#fldc = ncfilec.variables[ncfield][(styr-1)*12:(enyr*12+1),:,:,:]*conv
fldc = ncfilec.variables[ncfield][...]*conv
ncfilec.close()

# super-most vanilla
fldctmvan = np.average(fldc,0)
fldctmzmvan = np.average(fldctmvan,2)

# vanilla but remove extra lon
fldctmvanNOLON = fldctmvan[:,:,0:-1]
fldctmzmvanNOLON = np.average(fldctmvanNOLON,2)

# vanilla but time-weight
fldczmvan = np.average(fldc,3)
fldczmvanann = cutl.annualize_monthlyts(fldczmvan)
fldczmvananntm = np.average(fldczmvanann,0)

# vanilla but remove extra lon and time-weight
fldcNOLON = fldc[:,:,:,0:-1]
fldczmvanNOLON = np.average(fldcNOLON,3)
fldczmvanannNOLON = cutl.annualize_monthlyts(fldczmvanNOLON)
fldczmvananntmNOLON = np.average(fldczmvanannNOLON,0)

################################# PY/NC ########################################
lat = cnc.getNCvar(fnamec,'lat')
lev = cnc.getNCvar(fnamec,'plev')

#fldczm = cnc.getNCvar(fnamec,ncfield,calc='zm',timechunk=((styr-1)*12,enyr*12+1) )*conv
fldczm = cnc.getNCvar(fnamec,ncfield,calc='zm')*conv
#fldpzm = cnc.getNCvar(fnamep2,ncfield,calc='zm',timechunk=((styr-1)*12,enyr*12+1) )*conv
fldczmtm = np.mean(cutl.annualize_monthlyts(fldczm),0)
fldpzmtm = np.mean(cutl.annualize_monthlyts(fldpzm),0)
fldczmtm2 = np.average(fldczm,0)
fldpzmtm2 = np.average(fldpzm,0) # not month-weighted

################################# CDO ###########################################
## @@ test out cdo bindings and compare
# xm_fgco2 = cdo.zonmean( input = ifile ,returnMaArray ='fgco2' )
#   OR: select out specific dates, time-mean then zonal mean
# xm_fgco2 = cdo.zonmean( input = cdo.timmean( input = cdo.seldate('1990-01-01,2005-12-31', input=ifile ) ) ,returnMaArray ='fgco2')

# this works: fldcselcdo = cdo.seldate('0002-01-01,0061-12-31', input = fnamec, returnArray = ncfield )

#fldctmcdo = np.squeeze(cdo.timmean(input = cdo.seldate('0002-01-01,0061-12-31', input = fnamec), returnArray = ncfield ))

# @@@@ move these to functions so can call os.system()
## fldczmcdo = np.squeeze(
##     cdo.zonmean( input =
##                  cdo.timmean(input =
##                              cdo.seldate('0002-01-01,0061-12-31', input = fnamec ) ),
##                  returnMaArray  = ncfield))

fldczmcdo = np.squeeze(
    cdo.zonmean( input =
                 cdo.timmean(input = fnamec),returnMaArray  = ncfield))
os.system('rm -rf /tmp/cdoPy*')

## fldpzmcdo = np.squeeze(
##     cdo.zonmean( input =
##                  cdo.timmean(input =
##                              cdo.seldate('0002-01-01,0061-12-31', input = fnamep2 ) ),
##                  returnMaArray  = ncfield))  



lats,levs = np.meshgrid(lat,lev)



plotfld = fldczmtm - fldczmcdo

fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('Py/NC - CDO ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()


plotfld = fldczmtm2 - fldczmcdo

fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('Py/NC2 (unwgted) - CDO ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()


plotfld = fldctmzmvan - fldczmcdo
# this is the closest comparison. they are effectively the same.
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('Most vanilla - CDO ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()


plotfld = fldczmvananntm - fldczmcdo
# time-weighting appears to matter quite a lot, unless annualize() is wrong...
# there are still other differences here
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('Most vanilla (time-wgted) - CDO ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()


plotfld = fldczmvananntm - fldctmzmvan
# time-weighting appears to matter quite a lot, unless annualize() is wrong...
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('Most vanilla (time-wgted) - Most vanilla ' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()


plotfld = fldczmvananntm - fldczmtm
# differences here indicate that the "timechunk" isn't working as I think in getNCvar() --> nope, it's the extra longitude!
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(lats,levs/100,plotfld,cmap= plt.cm.get_cmap(cmap),shading='gouraud',vmin=-.12,vmax=.12)
ax.set_xlim(-90,90)
ax.set_ylim(10,1000)
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_yticks([1000,800, 500, 300, 100, 10])
ax.set_yticklabels((1000,800,500,300,100,10))
ax.set_title('Most vanilla (time-wgted) - Py/NC (tim-wgted)' + ncfield + ' (' + units + ')')
cbar = plt.colorbar()

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
# Get the data
ncfilec = Dataset(fnamec,'r') # control
ncfilep1 = Dataset(fnamep1,'r') # pert1
ncfilep2 = Dataset(fnamep2,'r') # pert2
ncfilep3 = Dataset(fnamep3,'r') # pert3

lat = ncfilec.variables['lat'][:]
lon = ncfilec.variables['lon'][:]
lev = ncfilec.variables['plev'][:]

if ftype == 'ts':
    fldc = ncfilec.variables[ncfield][(styr-1)*12:(enyr*12+1),:,:,:]*conv # time start year to end
    #fldp1 = ncfilep1.variables[ncfield][(styrp-1)*12:(enyrp*12+1),:,:,:]*conv # time start year to end
    fldp2 = ncfilep2.variables[ncfield][(styrp-1)*12:(enyrp*12+1),:,:,:]*conv # time start year to end
    #fldp3 = ncfilep3.variables[ncfield][(styrp-1)*12:(enyrp*12+1),:,:,:]*conv # time start year to end
else:
    fldc = ncfilec.variables[ncfield][...]*conv 
    fldp1 = ncfilep1.variables[ncfield][...]*conv
    fldp2 = ncfilep2.variables[ncfield][...]*conv # time start year to end
    fldp3 = ncfilep3.variables[ncfield][...]*conv # time start year to end

    
nt,nlev,nlat,nlon = fldc.shape # if nt == 12 then it's a climo

#  ############# set which perturbation run! ######
casenamep = casenamep2
fldp = fldp2
#  ###############################################

# take zonal mean first
# for zonal mean, remove the wraparound lon at the end
fldc = np.squeeze(fldc[:,:,:,0:-1])
fldp = np.squeeze(fldp[:,:,:,0:-1])
fldczm = np.mean(fldc,3)
fldpzm = np.mean(fldp,3)

# seasonalize
seasfldczm = cutl.seasonalize_monthlyts(fldczm,timeavg)
seasfldpzm = cutl.seasonalize_monthlyts(fldpzm,timeavg)
"""
