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
#import cccmaplots as cplt


# starting with climos. not sure what to do about stats.
# what is the right way to do this:
#    (historicalrcp45-historical) - (kem1pert2-kemctl1)
#          OR
#    historicalrcp45-kem1pert2
#    ? I think the first way for sure. Care about the forcing diffs.
# how similar are the control time periods?? historical v kemctl1


plt.close("all")
plt.ion()

printtofile=1
plotallmos=0 # make allmonth figs of CanESM/CanAM and their diff
seasonal=0 # averages seasons (DJF, MAM, JJA, SON)
monxlat=0
canesmens=1; ensnum=5 # look at ens members separately

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
# st, sicn, gt, pmsl, pcp, hfl, hfs, turb, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)
afield = 'sicn'

if afield == 'st':
    cfield = 'tas' # coupled (CMIP) field name
    units = 'K'
    aconv = 1; cconv=aconv  # no conversion
    cmin = -2; cmax = 2  # for anomaly plots
    cminp=-.5; cmaxp=.5 # for when pert is 'ctl'
    cminm = -3; cmaxm = 3   # monthly
    ## print 'small clim!'
    ## cmin = -1; cmax = 1  # for anomaly plots
    ## cminm = -1.5; cmaxm = 1.5   # monthly
    
    cminmp = -1; cmaxmp = 1 # for when pert is 'ctl'
    cmap = 'blue2red_w20'
elif afield == 'sicn':
    cfield = 'sic'
    comp = 'OImon'
    units = 'frac'
    aconv = 1; cconv=1/100. # make fraction
    cmin = -.15; cmax = .15  # for anomaly plots
    cminp=-.10; cmaxp=.10 # for when pert is 'ctl'
    cminm = -.15; cmaxm = .15   # monthly
    cmap = 'red2blue_w20'
elif afield == 'pmsl':
    cfield = 'psl'
    units = 'hPa' # pretty sure hpa @@double check
    aconv = 1; cconv=1 # @@@ double check
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

afldc = cnc.getNCvar(afnamec,afield.upper())*aconv
afldp = cnc.getNCvar(afnamep,afield.upper())*aconv
cfldc = cnc.getNCvar(cfnamec,cfield)*cconv
cfldp = cnc.getNCvar(cfnamep,cfield)*cconv
cfldc = np.dstack((cfldc,cfldc[...,0])) # add wraparound lon
cfldp = np.dstack((cfldp,cfldp[...,0])) # add wraparound lon


if plotallmos:
    if afield == 'sicn':
        print "Doesn't really make sense to plotallmos with SICN b/c by def they are the same"
        
    cplt.map_allmonths(afldp-afldc,lat,lon,title='CanAM4',cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh')
    if printtofile:
        plt.savefig('CanAM4_' + afield + '_' + acasenamep + '_v_' + acasename +
                     '_allmos_nh.pdf')

    cplt.map_allmonths(cfldp-cfldc,lat,lon,title='CanESM2',cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh')
    if printtofile:
        plt.savefig('CanESM2_' + cfield + '_' + ccasename + '_' + ctimstrp + '-' + ctimstr +
                     '_allmos_nh.pdf')

    cplt.map_allmonths((cfldp-cfldc)-(afldp-afldc),lat,lon,title='CanESM2diff-CanAM4diff',
                       cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh')
    if printtofile:
        plt.savefig('CanESM2-CanAM4_' + afield + '_' + ccasename + '_' + ctimstrp + '-' + ctimstr +
                     '_allmos_nh.pdf')



if seasonal:
    if afield == 'sicn':
        print "Doesn't really make sense to plot seasonal with SICN b/c by def they are the same"
    # seasonal figs
    cplt.map_allseas(afldp-afldc,lat,lon,title='CanAM4',cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',climo=1)
    if printtofile:
        plt.savefig('CanAM4_' + afield + '_' + acasenamep + '_v_' + acasename +
                     '_seas_nh.pdf')

    cplt.map_allseas(cfldp-cfldc,lat,lon,title='CanESM2',cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',climo=1)
    if printtofile:
        plt.savefig('CanESM2_' + cfield + '_' + ccasename + '_' + ctimstrp + '-' + ctimstr +
                     '_seas_nh.pdf')

    cplt.map_allseas((cfldp-cfldc)-(afldp-afldc),lat,lon,title='CanESM2diff-CanAM4diff',
                       cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',climo=1)
    if printtofile:
        plt.savefig('CanESM2-CanAM4_' + afield + '_' + ccasename + '_' + ctimstrp + '-' + ctimstr +
                     '_seas_nh.pdf')


if monxlat:

    # zonal means
    afldczm = np.mean(afldc[...,:-1],axis=2)
    afldpzm = np.mean(afldp[...,:-1],axis=2)
    cfldczm = np.mean(cfldc[...,:-1],axis=2)
    cfldpzm = np.mean(cfldp[...,:-1],axis=2)

    # Yay, also works.
    #afldczm2 = cnc.getNCvar(afnamec,afield.upper(),calc='zm')*aconv
    #afldpzm2 = cnc.getNCvar(afnamep,afield.upper(),calc='zm')*aconv

    mon = np.arange(1,afldc.shape[0]+1)

    # month by latitude figs
    lats,mons = np.meshgrid(lat,mon)

    cmlen=float( plt.cm.get_cmap(cmap).N) 
    incr = (cmaxm-cminm) / (cmlen)
    conts = np.arange(cminm,cmaxm+incr,incr)


    fig1,ax1 = plt.subplots(1,3)
    fig1.set_size_inches(10,6)
    plotfld = cfldpzm-cfldczm

    pc = ax1[0].contourf(mons,lats,plotfld,
                      cmap=plt.cm.get_cmap(cmap),levels=conts,
                      vmin=cminm,vmax=cmaxm,extend='both')
    ax1[0].set_title('Full historical')
    ax1[0].set_xlabel('months')
    ax1[0].set_ylabel('lat')

    plotfld = afldpzm-afldczm

    pc = ax1[1].contourf(mons,lats,plotfld,
                      cmap=plt.cm.get_cmap(cmap),levels=conts,
                      vmin=cminm,vmax=cmaxm,extend='both')
    ax1[1].set_title('Sea-ice only')
    ax1[1].set_xlabel('months')

    plotfld = (cfldpzm-cfldczm)-(afldpzm-afldczm)

    pc = ax1[2].contourf(mons,lats,plotfld,
                      cmap=plt.cm.get_cmap(cmap),levels=conts,
                      vmin=cminm,vmax=cmaxm,extend='both')
    ax1[2].set_title('Residual')
    ax1[2].set_xlabel('months')
    cbar_ax = fig1.add_axes([.91,.15, .02,.7])
    fig1.colorbar(pc,cax=cbar_ax)

    if printtofile:
        fig1.savefig('CanESM2_CanAM4_andDIFF_' + afield + '_' + ccasename + '_' +
                     ctimstrp + '-' + ctimstr + '_monxlat.pdf')


    # Northern HEM only
    fig1,ax1 = plt.subplots(1,3)
    fig1.set_size_inches(10,4)
    plotfld = cfldpzm-cfldczm

    pc = ax1[0].contourf(mons,lats,plotfld,
                      cmap=plt.cm.get_cmap(cmap),levels=conts,
                      vmin=cminm,vmax=cmaxm,extend='both')
    ax1[0].set_title('Full historical')
    ax1[0].set_xlabel('months')
    ax1[0].set_ylabel('lat')
    ax1[0].set_ylim(0,90)

    plotfld = afldpzm-afldczm

    pc = ax1[1].contourf(mons,lats,plotfld,
                      cmap=plt.cm.get_cmap(cmap),levels=conts,
                      vmin=cminm,vmax=cmaxm,extend='both')
    ax1[1].set_title('Sea-ice only')
    ax1[1].set_xlabel('months')
    ax1[1].set_ylim(0,90)

    plotfld = (cfldpzm-cfldczm)-(afldpzm-afldczm)

    pc = ax1[2].contourf(mons,lats,plotfld,
                      cmap=plt.cm.get_cmap(cmap),levels=conts,
                      vmin=cminm,vmax=cmaxm,extend='both')
    ax1[2].set_title('Residual')
    ax1[2].set_xlabel('months')
    ax1[2].set_ylim(0,90)
    cbar_ax = fig1.add_axes([.91,.15, .02,.7])
    fig1.colorbar(pc,cax=cbar_ax)

    if printtofile:
        fig1.savefig('CanESM2_CanAM4_andDIFF_' + afield + '_' + ccasename + '_' +
                     ctimstrp + '-' + ctimstr + '_monxlat_nh.pdf')


if canesmens:
    # plot each ensemble member individually
    # do they resemble hadisst runs at all? does one?

    for eii in range(1,ensnum+1): # five ens members
        
        cfnamec = cbasepath + ccasename + '/' + cfield + '/' + cfield +\
                  '_' + comp + '_' + cmodel + '_' + ccasename +\
                  '_r' + str(eii) + 'i1p1_' + ctimstr + 'climo.nc'
        cfnamep = cbasepath + ccasenamep + '/' + cfield + '/' + cfield +\
                  '_' + comp + '_' + cmodel + '_' + ccasenamep +\
                  '_r' + str(eii) + 'i1p1_' + ctimstrp + 'climo.nc'


        cfldc = cnc.getNCvar(cfnamec,cfield)*cconv
        cfldp = cnc.getNCvar(cfnamep,cfield)*cconv
        cfldc = np.dstack((cfldc,cfldc[...,0])) # add wraparound lon
        cfldp = np.dstack((cfldp,cfldp[...,0])) # add wraparound lon

        fig = cplt.map_allmonths(cfldp-cfldc,lat,lon,
                                 title='CanESM r' + str(eii) + ' ' + cfield,
                                 cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh')
        if printtofile:
            fig.savefig('CanESMr' + str(eii) + '_' + cfield + '_' +
                        ccasename + '_' + ctimstrp + '-' + ctimstr +
                        '_allmos_nh.pdf')


        fig = cplt.map_allseas(cfldp-cfldc,lat,lon,title='CanESM r' + str(eii) + ' ' + cfield,
                       cmin=cminm,cmax=cmaxm,cmap=cmap,type='nh',climo=1)
        if printtofile:
            fig.savefig('CanESMr' + str(eii) + '_' + cfield + '_' +
                        ccasename + '_' + ctimstrp + '-' + ctimstr +
                        '_seas_nh.pdf')


        


        

    


