"""
  compare_canamtocanesm.py
    4/30/2014: compare the CanAM4 runs with historical sea ice
               boundary conditions to the full CanESM2 historical
               simulation: subtracting CanESM-CanAM should give an
               idea of what part of the signal is due to sea ice
               and what is due to the forcing
               Use canam4sims_stats2.py as guide.
"""
import scipy.stats
import matplotlib.cm as cm
import datetime as datetime
import matplotlib.colors as col
import platform as platform
import constants as con      # my module
import cccmautils as cutl    # my module
import matplotlib.font_manager as fm
#import copy

# starting with climos. not sure what to do about stats.
# what is the right way to do this:
#    (historicalrcp45-historical) - (kem1pert2-kemctl1)
#          OR
#    historicalrcp45-kem1pert2
#    ?
# how similar are the control time periods?? historical v kemctl1


plt.close("all")
plt.ion()

printtofile=0
plotann=1    # annual average
plotallmos=1 # each month separately
bimos=0 # averages every 2 mos (JF, MA, MJ, JA, SO, ND) @@ add
seasonal=1 # averages seasons (DJF, MAM, JJA, SON)


sigtype = 'cont' # significance: 'cont' or 'hatch' which is default

seasons = 'DJF','MAM','JJA','SON'

amodel = 'CanAM4'
cmodel = 'CanESM2'

# # # ########### set Simulations #############
# # CANESM: coupled
ccasename = 'historical'
ctimstr = '1979-1989'
ccasenamep = 'historicalrcp45'
ctimstrp = '2002-2012'

# # CANAM: atmosphere
acasename = 'kemctl1'
atimstr = '001-111'
atimstrp = '001-111'
acasenamep = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
# same timstrs as coupled

comp = 'Amon'


# # # ######## set Field info (CanAM name) ###################
# st, sic, gt, pmsl, pcp, hfl, hfs, turb, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)
afield = 'st'

if afield == 'st':
    cfield = 'tas' # coupled (CMIP) field name
    units = 'K'
    conv = 1  # no conversion
    cmin = -2; cmax = 2  # for anomaly plots
    cminp=-.5; cmaxp=.5 # for when pert is 'ctl'
    cminm = -3; cmaxm = 3   # monthly
    ## print 'small clim!'
    ## cmin = -1; cmax = 1  # for anomaly plots
    ## cminm = -1.5; cmaxm = 1.5   # monthly
    
    cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
    cmap = 'blue2red_w20'
elif afield == 'pmsl':
    cfield = 'psl'
    units = 'hPa' # pretty sure hpa @@double check
    conv = 1
    cmin = -1; cmax = 1  # for anomaly plots
    cminm=-2; cmaxm=2  # for monthly maps
    cminp=cmin; cmaxp=cmax # for when pert is 'ctl'
    cminmp=cminm; cmaxmp=cmaxm
    cmap = 'blue2red_20'

else:
    print 'No settings for ' + afield


# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
    subdir = '/'
else:  # on linux workstation in Vic
    abasepath = '/home/rkm/work/DATA/' + amodel + '/'
    asubdir = '/ts/'
    cbasepath = '/home/rkm/work/DATA/' + cmodel + '/'

afnamec = abasepath + acasename + asubdir + acasename + '_' + afield + '_' + atimstr + '_climo.nc'
afnamep = abasepath + acasenamep + asubdir + acasenamep + '_' + afield + '_' + atimstrp + '_climo.nc'
cfnamec = cbasepath + ccasename + '/' + cfield + '/' + cfield + '_' + comp + '_' + cmodel +\
          '_' + ccasename + '_ens_' + ctimstr + 'climo.nc'
cfnamep = cbasepath + ccasenamep + '/' + cfield + '/' + cfield + '_' + comp + '_' + cmodel +\
          '_' + ccasenamep + '_ens_' + ctimstrp + 'climo.nc'

lat = cnc.getNCvar(afnamec,'lat')
lon = cnc.getNCvar(afnamec,'lon')

afldc = cnc.getNCvar(afnamec,afield.upper())*conv
afldp = cnc.getNCvar(afnamep,afield.upper())*conv
cfldc = cnc.getNCvar(cfnamec,cfield)*conv
cfldp = cnc.getNCvar(cfnamep,cfield)*conv
cfldc = np.dstack((cfldc,cfldc[...,0])) # add wraparound lon
cfldp = np.dstack((cfldp,cfldp[...,0])) # add wraparound lon

mon = np.arange(1,afldc.shape[0]+1)

cplt.map_allmonths(afldp-afldc,lat,lon,title='CanAM4',cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh')

cplt.map_allmonths(cfldp-cfldc,lat,lon,title='CanESM2',cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh')

cplt.map_allmonths((cfldp-cfldc)-(afldp-afldc),lat,lon,title='CanESM2diff-CanAM4diff',
                   cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh')

cplt.map_allmonths(cfldp-afldp,lat,lon,title='CanESM2pert-CanAM4pert',
                   cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh') # think this is the wrong way to look at it


# zonal means
afldczm = np.mean(afldc[...,:-1],axis=2)
afldpzm = np.mean(afldp[...,:-1],axis=2)
cfldczm = np.mean(cfldc[...,:-1],axis=2)
cfldpzm = np.mean(cfldp[...,:-1],axis=2)

# Yay, also works.
#afldczm2 = cnc.getNCvar(afnamec,afield.upper(),calc='zm')*conv
#afldpzm2 = cnc.getNCvar(afnamep,afield.upper(),calc='zm')*conv


lats,mons = np.meshgrid(lat,mon)

cmlen=float( plt.cm.get_cmap(cmap).N) 
incr = (cmaxm-cminm) / (cmlen)
conts = np.arange(cminm,cmaxm+incr,incr)



plotfld = afldpzm-afldczm

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
pc = plt.contourf(mons,lats,plotfld,
                  cmap=plt.cm.get_cmap(cmap),levels=conts,
                  vmin=cminm,vmax=cmaxm,extend='both')
ax1.set_title('CanAM4')
ax1.set_xlabel('months')
ax1.set_ylabel('lat')
cbar = fig1.colorbar(pc)


plotfld = cfldpzm-cfldczm

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
pc = plt.contourf(mons,lats,plotfld,
                  cmap=plt.cm.get_cmap(cmap),levels=conts,
                  vmin=cminm,vmax=cmaxm,extend='both')
ax2.set_title('CanESM2')
ax2.set_xlabel('months')
ax2.set_ylabel('lat')
cbar = fig2.colorbar(pc)


plotfld = (cfldpzm-cfldczm)-(afldpzm-afldczm)

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
pc = plt.contourf(mons,lats,plotfld,
                  cmap=plt.cm.get_cmap(cmap),levels=conts,
                  vmin=cminm,vmax=cmaxm,extend='both')
ax3.set_title('CanESM2diff-CanAM4diff')
ax3.set_xlabel('months')
ax3.set_ylabel('lat')
cbar = fig3.colorbar(pc)



