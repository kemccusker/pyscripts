"""
  plot_timeseries.py: taken from compare_canamtocanesm.py 5/6/2014
                      Want CanESM2 individual ens members, plus mean,
                      plus obs for sea ice
                      
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



plt.close("all")
plt.ion()

printtofile=False

seasons = 'DJF','MAM','JJA','SON'

amodel = 'CanAM4'
cmodel = 'CanESM2'

# # # ########### set Simulations #############
# # CANESM: coupled
ccasename = 'historical'
ctimstr = '1979-1989'
ccasenamep = 'historicalrcp45'
ctimstrp = '2002-2012'
ctimeper = '185001-201212'

# # CANAM: atmosphere
acasename = 'kemctl1'
atimstr = '001-111'
atimstrp = '001-111'
acasenamep = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
# same timstrs as coupled

comp = 'Amon'
ensnum=5

# # # ######## set Field info (CanAM name) ###################
# st, sicn, gt, pmsl, pcp, hfl, hfs, turb, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)
afield = 'st'

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
    basepath = '/home/rkm/work/DATA/' + cmodel + '/'
    basepath2 = '/home/rkm/work/BCs/'

#afnamec = abasepath + acasename + asubdir + acasename + '_' + afield + '_' + atimstr + '_climo.nc'
#afnamep = abasepath + acasenamep + asubdir + acasenamep + '_' + afield + '_' + atimstrp + '_climo.nc'
# sic_OImon_CanESM2_historicalrcp45_r1i1p1_185001-201212.nc

## if ccasename=='historical':
##     ccasename='historicalrcp45'
    
# loop through ens members:
cfldcall = np.zeros((ensnum,12)) # a climo
cfldpall = np.zeros((ensnum,163*12)) # hard-coded monthly timeseries 1850-2012
for eii in range(1,ensnum+1): # five ens members
        
        cfnamec = basepath + ccasename + '/' + cfield + '/' + cfield +\
                  '_' + comp + '_' + cmodel + '_' + ccasename +\
                  '_r' + str(eii) + 'i1p1_' + ctimstr + 'climo.nc' # to take anomaly from
        cfnamep = basepath + ccasenamep + '/' + cfield + '/' + cfield +\
                  '_' + comp + '_' + cmodel + '_' + ccasenamep +\
                  '_r' + str(eii) + 'i1p1_' + ctimeper + '.nc'  # sic, 128x64


        cfldc = cnc.getNCvar(cfnamec,cfield)*cconv # ctl time period avg
        cfldp = cnc.getNCvar(cfnamep,cfield)*cconv # full timeseries
        lat = cnc.getNCvar(cfnamec,'lat')
        lon = cnc.getNCvar(cfnamec,'lon')
        
        ## cfldc = np.dstack((cfldc,cfldc[...,0])) # add wraparound lon
        ## cfldp = np.dstack((cfldp,cfldp[...,0])) # add wraparound lon

        # if sea ice frac: calc area and save it
        if cfield == 'sic':

            # calc sea ice area
            # mult fraction by grid cell area & sum
            areas = cutl.calc_cellareas(lat,lon)
            careasp = np.tile(areas,(cfldp.shape[0],1,1)) # need one per time
            careasc = np.tile(areas,(cfldc.shape[0],1,1)) # need one per month (climo)
            cfldc = cfldc*careasc
            cfldp = cfldp*careasp

            cfldcall[eii-1,:] = np.sum(np.sum(cfldc[:,lat>0,:],2),1) # NH total ice area
            cfldpall[eii-1,:] = np.sum(np.sum(cfldp[:,lat>0,:],2),1) # NH total ice area
        else:
            latlim = 60 # 60N
            # just do polar mean for now
            print 'doing polar mean of ' + cfield
            cfldcall[eii-1,:] = cutl.polar_mean_areawgted3d(cfldc,lat,lon,latlim=latlim)
            cfldpall[eii-1,:] = cutl.polar_mean_areawgted3d(cfldp,lat,lon,latlim=latlim)
            
        
# Now get obs for sea ice
if cfield=='sic':

    fhadsicc = basepath2 + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sicn_' +\
              ctimstr + 'climo.nc' #SICN, 129x64 CLIMO
    fhadsicp = basepath2 + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sicn_1870010100-2013030100.nc'

    fnsidcsicc = basepath2 + 'NSIDC/nsidc_bt_128x64_1978m11_2011m12_sicn_' + ctimstr + 'climo.nc'
    fnsidcsicp = basepath2 + 'NSIDC/nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc' #SICN, 129x64
    # nsidc_bt_128x64_1978m11_2011m12_sicn_1978111600-2011121612.nc
    

    hadsicc = cnc.getNCvar(fhadsicc,'SICN') # climo
    hadsicp = cnc.getNCvar(fhadsicp,'SICN',timesel='1979-01-01,2012-12-31') # timeseries @@ note will not work on mac

    nsidcsicc = cnc.getNCvar(fnsidcsicc,'SICN') # climo
    nsidcsicp = cnc.getNCvar(fnsidcsicp,'SICN',timesel='1979-01-01,2012-12-31') # timeseries @@ note will not work on mac


    hadsicc = hadsicc[...,:-1]
    hadsicp = hadsicp[...,:-1]
    nsidcsicc = nsidcsicc[...,:-1]
    nsidcsicp = nsidcsicp[...,:-1]

    hareasp = np.tile(areas,(hadsicp.shape[0],1,1)) # need one per time
    hareasc = careasc # need one per month (climo)
    hfldc = hadsicc*hareasc
    hfldp = hadsicp*hareasp
    hfldc = np.sum(np.sum(hfldc[:,lat>0,:],2),1) # NH total ice area
    hfldp = np.sum(np.sum(hfldp[:,lat>0,:],2),1) # NH total ice area

    nareasp = np.tile(areas,(nsidcsicp.shape[0],1,1))
    nareasc = careasc
    nfldc = nsidcsicc*nareasc
    nfldp = nsidcsicp*nareasp
    nfldc = np.sum(np.sum(nfldc[:,lat>0,:],2),1)
    nfldp = np.sum(np.sum(nfldp[:,lat>0,:],2),1)



years = np.arange(1850,2013)
hyears = np.arange(1979,2013)
nyears = np.arange(1979,2012)

darkolivegreen1 = np.array([202, 255, 112])/255 # terrible
darkolivegreen3 = np.array([162, 205, 90])/255.
darkseagreen = np.array([143, 188, 143])/255.
darkseagreen4 = np.array([105, 139, 105])/255.
dodgerblue = np.array([30, 144, 255])/255. 
orangered4 = np.array([139, 37, 0])/255.

#  ANNUAL
plt.figure()
plt.plot(years,cutl.annualize_monthlyts(cfldpall[0,:]),'0.25')
plt.plot(years,cutl.annualize_monthlyts(cfldpall[1,:]),'0.4')
plt.plot(years,cutl.annualize_monthlyts(cfldpall[2,:]),'0.55')
plt.plot(years,cutl.annualize_monthlyts(cfldpall[3,:]),'0.7')
plt.plot(years,cutl.annualize_monthlyts(cfldpall[4,:]),'0.85')
plt.plot(years,cutl.annualize_monthlyts( np.mean(cfldpall,axis=0) ),color='k',linewidth=2)
if afield=='sicn':
    plt.plot(hyears,cutl.annualize_monthlyts(hfldp),color=orangered4,linewidth=2)
    plt.plot(nyears,cutl.annualize_monthlyts(nfldp),color=dodgerblue,linewidth=2)      
plt.xlim((1850,2012))
plt.title('ANN NH SIA')
plt.grid()
if printtofile:
    plt.savefig('CanESMens_OBS_' + cfield + '_ANN_timeseries.pdf')


plt.figure()
plt.plot(years,cutl.annualize_monthlyts(cfldpall[0,:])-
         cutl.annualize_monthlyts(cfldcall[0,:]),'0.25')
plt.plot(years,cutl.annualize_monthlyts(cfldpall[1,:])-
         cutl.annualize_monthlyts(cfldcall[1,:]),'0.4')
plt.plot(years,cutl.annualize_monthlyts(cfldpall[2,:])-
         cutl.annualize_monthlyts(cfldcall[2,:]),'0.55')
plt.plot(years,cutl.annualize_monthlyts(cfldpall[3,:])-
         cutl.annualize_monthlyts(cfldcall[3,:]),'0.7')
plt.plot(years,cutl.annualize_monthlyts(cfldpall[4,:])-
         cutl.annualize_monthlyts(cfldcall[4,:]),'0.85')
plt.plot(years,cutl.annualize_monthlyts(np.mean(cfldpall,axis=0) )-
         cutl.annualize_monthlyts( np.mean(cfldcall,axis=0) ),
         color='k',linewidth=2)
if afield=='sicn':
    plt.plot(hyears,cutl.annualize_monthlyts(hfldp)-
         cutl.annualize_monthlyts(hfldc),color=orangered4,linewidth=2)
    plt.plot(nyears,cutl.annualize_monthlyts(nfldp)-
         cutl.annualize_monthlyts(nfldc),color=dodgerblue,linewidth=2)      
plt.xlim((1850,2012))
plt.title('ANN NH SIA anom from 1979-89')
plt.grid()
if printtofile:
    plt.savefig('CanESMens_OBSanom' + ctimstr + '_' + cfield + '_ANN_timeseries.pdf')



#   MINIMUM (SEPT)
plt.figure()
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[0,:],mo=9),'0.25')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[1,:],mo=9),'0.4')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[2,:],mo=9),'0.55')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[3,:],mo=9),'0.7')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[4,:],mo=9),'0.85')
plt.plot(years,cutl.seasonalize_monthlyts( np.mean(cfldpall,axis=0),mo=9 ),color='k',linewidth=2)
if afield=='sicn':
    plt.plot(hyears,cutl.seasonalize_monthlyts(hfldp,mo=9),color=orangered4,linewidth=2)
    plt.plot(nyears,cutl.seasonalize_monthlyts(nfldp,mo=9),color=dodgerblue,linewidth=2)      
plt.xlim((1850,2012))
plt.title('September NH SIA')
plt.grid()
if printtofile:
    plt.savefig('CanESMens_OBS_' + cfield + '_Sep_timeseries.pdf')


plt.figure()
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[0,:],mo=9)-
         cutl.seasonalize_monthlyts(cfldcall[0,:],mo=9,climo=1),'0.25')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[1,:],mo=9)-
         cutl.seasonalize_monthlyts(cfldcall[1,:],mo=9,climo=1),'0.4')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[2,:],mo=9)-
         cutl.seasonalize_monthlyts(cfldcall[2,:],mo=9,climo=1),'0.55')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[3,:],mo=9)-
         cutl.seasonalize_monthlyts(cfldcall[3,:],mo=9,climo=1),'0.7')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[4,:],mo=9)-
         cutl.seasonalize_monthlyts(cfldcall[4,:],mo=9,climo=1),'0.85')
plt.plot(years,cutl.seasonalize_monthlyts(np.mean(cfldpall,axis=0),mo=9 )-
         cutl.seasonalize_monthlyts( np.mean(cfldcall,axis=0),mo=9,climo=1 ),
         color='k',linewidth=2)
if afield=='sicn':
    plt.plot(hyears,cutl.seasonalize_monthlyts(hfldp,mo=9)-
         cutl.seasonalize_monthlyts(hfldc,mo=9,climo=1),color=orangered4,linewidth=2)
    plt.plot(nyears,cutl.seasonalize_monthlyts(nfldp,mo=9)-
         cutl.seasonalize_monthlyts(nfldc,mo=9,climo=1),color=dodgerblue,linewidth=2)      
plt.xlim((1850,2012))
plt.title('September NH SIA anom from 1979-89')
plt.grid()
if printtofile:
    plt.savefig('CanESMens_OBSanom' + ctimstr + '_' + cfield + '_Sep_timeseries.pdf')



#   MAXIMUM (MAR)
plt.figure()
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[0,:],mo=3),'0.25')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[1,:],mo=3),'0.4')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[2,:],mo=3),'0.55')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[3,:],mo=3),'0.7')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[4,:],mo=3),'0.85')
plt.plot(years,cutl.seasonalize_monthlyts( np.mean(cfldpall,axis=0),mo=3 ),color='k',linewidth=2)
if afield=='sicn':
    plt.plot(hyears,cutl.seasonalize_monthlyts(hfldp,mo=3),color=orangered4,linewidth=2)
    plt.plot(nyears,cutl.seasonalize_monthlyts(nfldp,mo=3),color=dodgerblue,linewidth=2)      
plt.xlim((1850,2012))
plt.title('March NH SIA')
plt.grid()
if printtofile:
    plt.savefig('CanESMens_OBS_' + cfield + '_Mar_timeseries.pdf')


plt.figure()
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[0,:],mo=3)-
         cutl.seasonalize_monthlyts(cfldcall[0,:],mo=3,climo=1),'0.25')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[1,:],mo=3)-
         cutl.seasonalize_monthlyts(cfldcall[1,:],mo=3,climo=1),'0.4')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[2,:],mo=3)-
         cutl.seasonalize_monthlyts(cfldcall[2,:],mo=3,climo=1),'0.55')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[3,:],mo=3)-
         cutl.seasonalize_monthlyts(cfldcall[3,:],mo=3,climo=1),'0.7')
plt.plot(years,cutl.seasonalize_monthlyts(cfldpall[4,:],mo=3)-
         cutl.seasonalize_monthlyts(cfldcall[4,:],mo=3,climo=1),'0.85')
plt.plot(years,cutl.seasonalize_monthlyts(np.mean(cfldpall,axis=0),mo=3 )-
         cutl.seasonalize_monthlyts( np.mean(cfldcall,axis=0),mo=3,climo=1 ),
         color='k',linewidth=2)
if afield=='sicn':
    plt.plot(hyears,cutl.seasonalize_monthlyts(hfldp,mo=3)-
         cutl.seasonalize_monthlyts(hfldc,mo=3,climo=1),color=orangered4,linewidth=2)
    plt.plot(nyears,cutl.seasonalize_monthlyts(nfldp,mo=3)-
         cutl.seasonalize_monthlyts(nfldc,mo=3,climo=1),color=dodgerblue,linewidth=2)      
plt.xlim((1850,2012))
plt.title('March NH SIA anom from 1979-89')
plt.grid()
if printtofile:
    plt.savefig('CanESMens_OBSanom' + ctimstr + '_' + cfield + '_Mar_timeseries.pdf')



        


        

    


