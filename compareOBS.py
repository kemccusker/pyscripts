""" compareOBS.py
    3/7/2014
    compare HadISST, Hurrell blended hadISST/NOAA, NSIDC (and CanESM2)
"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import platform as platform
import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import cccmaNC as cnc
import cccmacmaps

cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)
cnc = reload(cnc)

plt.close("all")
plt.ion()

printtofile=1
testsic=0

mtype='nh'  # map type (nh,sh,sq)
plotallmos=0
plotseacyc=0
plotseacycmag=1 # only for SST

had=0
hurr=0
canesm=1 # assume pert2 if "doBCs"
#else NSIDC
# @@ something is way wrong with CanESM2 SST files and/or figures @@
# @@ OR is is HadISST? basically the comparison for ssts north of 60N is horrible.
#    there is no seasonal cycle in HadISST but there is one that looks ok in CanESM2.
#    what's the seasonal cycle in actual observed ssts? it's only about 3-4 degrees
#    so something is wrong w/ my CanESM calc. Potentially it's just related to masking
#    and/or polar average function.

sic=0
sicn=0
#else SST

doBCs=1  # use the actual BC files @@so far only hadisst sst
latlim=50

flipmask=0 # flip the landmask (for HURRELL)

deni = 913 # density of ice
bcstr = ''

if sicn:
    # SICN caxis info
    cminm=-.25; cmaxm=.25
    cminm=-.15; cmaxm=.15 # @@
    cmap = 'red2blue_w20'
elif sic:
    cminm=-.5*deni
    cmaxm = .5*deni
    cmap = 'red2blue_w20'
else: # sst
    cminm=-2
    cmaxm=2
    cmap = 'blue2red_w20'

# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    #basepath = '/Users/kelly/CCCma/CanSISE/DATA/' # @@@
    basepath = '/Volumes/MyPassport1TB/DATA/'
    basepath2 = '/Users/kelly/CCCma/CanSISE/BoundaryConditionFiles/'
    temppath = '/Users/kelly/CCCma/CanSISE/matlab/'
    
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/BCs/'
    basepath2 = basepath
    temppath = basepath

timeper='1979-1989'
timeperp='2002-2012' # CanESM2 (or HadISST)
timeperp2='2002-2011' # Hadisst, Hurrell and NSIDC bootstrap

if sicn:
    # @@ check what the units are: are they all fraction?
    fhadsic = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sicn_' + timeper + 'climo.nc' #SICN, 129x64
    fhadsicts = basepath + 'HadISST/SICN_BC_HadISST_' + timeper + '_0000120100-0125120100.nc' # test matlab generated BC file
    
    fhurrsic = basepath + 'HURRELL/MODEL_ICE.T42_' + timeper + 'climo.nc' #SEAICE (%), 128x64
    fnsidcsic = basepath + 'NSIDC/nsidc_bt_128x64_1978m11_2011m12_sicn_' + timeper + 'climo.nc' #SICN, 129x64

    fhadsicp = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sicn_' + timeperp2 + 'climo.nc' #SICN, 129x64
    fhadsicpts = basepath + 'HadISST/SICN_BC_HadISST_' + timeperp2 + '_0000120100-0125120100.nc' # test matlab generated BC file
    fhurrsicp = basepath + 'HURRELL/MODEL_ICE.T42_' + timeperp2 + 'climo.nc' #SEAICE (%), 128x64
    fnsidcsicp = basepath + 'NSIDC/nsidc_bt_128x64_1978m11_2011m12_sicn_' + timeperp2 + 'climo.nc' #SICN, 129x64

    hadsicc = cnc.getNCvar(fhadsic,'SICN')
    hadsiccts = cnc.getNCvar(fhadsicts,'SICN')
    hadsiccclimo,junk = cutl.climatologize3d(hadsiccts)
    
    hurrsicc = cnc.getNCvar(fhurrsic,'SEAICE')/100
    hurrtimes = cnc.getNCvar(fhurrsic,'time')
    #hurrsicc.resize(hadsicc.shape) # can't resize like this. ??
    #np.append(hurrsicc,hurrsicc[:,:,0])#,axis=2) # could not get append to work
    #hurrsicc[:,:,len(lon)-1] = hurrsicc[:,:,0] # add a wraparound. nope.
    hurrsicc = np.flipud(hurrsicc)  # the lats are flipped compared to hadisst and nsidc
    nsidcsicc = cnc.getNCvar(fnsidcsic,'SICN')

    hadsicp = cnc.getNCvar(fhadsicp,'SICN')
    hadsicpts = cnc.getNCvar(fhadsicpts,'SICN')
    hadsicpclimo,junk = cutl.climatologize3d(hadsicpts)
    
    hurrsicp = cnc.getNCvar(fhurrsicp,'SEAICE')/100
    #hurrsicc.resize(hadsicc.shape) # other datasets have extra lon.
    #hurrsicp[:,:,len(lon)-1] = hurrsicp[:,:,0] # add a wraparound
    hurrsicp = np.flipud(hurrsicp)  # the lats are flipped compared to hadisst and nsidc
    nsidcsicp = cnc.getNCvar(fnsidcsicp,'SICN')
    nsidcsicc = nsidcsicc[:,:,0:-1]
    nsidcsicp = nsidcsicp[:,:,0:-1]

    nsidcsicd = nsidcsicp-nsidcsicc

    # add CanESM2 boundary conditions too
    fcansic = basepath + 'CanESM2/SICN_BC_CanESM2_historical' + timeper + '_climo.nc'
    fcansicp = basepath + 'CanESM2/SICN_BC_CanESM2_historical' + timeperp + '_climo.nc'
    
    cansicc = cnc.getNCvar(fcansic,'SICN')
    cansicp = cnc.getNCvar(fcansicp,'SICN')
    cansicc = cansicc[:,:,0:-1]
    cansicp = cansicp[:,:,0:-1]
    cansicd = cansicp - cansicc
    
elif sic:
    
    fhadsic = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sic_' + timeper + 'climo.nc' #SIC, 129x64
    fhadsicp = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sic_' + timeperp2 + 'climo.nc' #SIC, 129x64
    fhadsicts = temppath + 'HadISST/SIC_BC_HadISST_' + timeper + '_1870010100-2011020100.nc' # test matlab generated BC file
    fhadsicpts = temppath + 'HadISST/SIC_BC_HadISST_' + timeperp2 + '_1870010100-2011020100.nc' # test matlab generated BC file

    fhadsicorig = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sic_1870010100-2013030100.nc'
    fhadsicnorig = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sicn_1870010100-2013030100.nc'

    hadsicc = cnc.getNCvar(fhadsic,'SIC')
    hadsicp = cnc.getNCvar(fhadsicp,'SIC')
    hadsiccts = cnc.getNCvar(fhadsicts,'SIC')
    hadsiccclimo,junk = cutl.climatologize3d(hadsiccts)
    hadsicpts = cnc.getNCvar(fhadsicpts,'SIC')
    hadsicpclimo,junk = cutl.climatologize3d(hadsicpts)

    # also need sicn for averaging
    fhadsicn = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sicn_' + timeper + 'climo.nc' #SICN, 129x64
    fhadsicnp = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sicn_' + timeperp2 + 'climo.nc' #SICN, 129x64

    hadsicnc = cnc.getNCvar(fhadsicn,'SICN')
    hadsicnp = cnc.getNCvar(fhadsicnp,'SICN')


    # add CanESM2 boundary conditions too
    fcansic = basepath2 + 'CanESM2/SIC_BC_CanESM2_historical' + timeper + '_climo.nc'
    fcansicp = basepath2 + 'CanESM2/SIC_BC_CanESM2_historical' + timeperp + '_climo.nc'
    
    cansicc = cnc.getNCvar(fcansic,'SIC')
    cansicp = cnc.getNCvar(fcansicp,'SIC')
    cansicc = cansicc[:,:,0:-1]
    cansicp = cansicp[:,:,0:-1]
    cansicd = cansicp - cansicc

    # add CanESM2 boundary conditions too
    fcansicn = basepath2 + 'CanESM2/SICN_BC_CanESM2_historical' + timeper + '_climo.nc'
    fcansicnp = basepath2 + 'CanESM2/SICN_BC_CanESM2_historical' + timeperp + '_climo.nc'
    
    cansicnc = cnc.getNCvar(fcansicn,'SICN')
    cansicnp = cnc.getNCvar(fcansicnp,'SICN')
    cansicnc = cansicnc[:,:,0:-1]
    cansicnp = cansicnp[:,:,0:-1]
    cansicnd = cansicnp - cansicnc

else:
    # SST (var names are sic still...)

    if doBCs:
        bcstr='BC'
        
        # use the actual BC files for the simulations
        fhadsic = basepath + 'HadISST/hadisst_kemhadctl_128x64_0001_0125_gt.nc'
        fhadsicp = basepath + 'HadISST/hadisst_kemhadpert_128x64_0001_0125_gt.nc'
        # or GTadjusted_BC_HadISST_2002-2011_0000120100-0125120100_abs10thresh.nc

        hadsicc,hadsiccstd = cutl.climatologize(cnc.getNCvar(fhadsic,'GT')) # I think, K?
        hadsicp,hadsicpstd = cutl.climatologize(cnc.getNCvar(fhadsicp,'GT'))
        
        # @@ Hurrell files are same as not doBCs
        fhurrsic = basepath + 'HURRELL/MODEL_SST.T42_' + timeper + 'climo.nc' #SST, degC, 128x64
        fhurrsicp = basepath + 'HURRELL/MODEL_SST.T42_' + timeperp2 + 'climo.nc' #SST, degC, 128x64
        
    else:
        fhadsic = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_gt_' + timeper + 'climo.nc' #GT, K, 129x64
        fhurrsic = basepath + 'HURRELL/MODEL_SST.T42_' + timeper + 'climo.nc' #SST, degC, 128x64
        fhadsicp = basepath + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_gt_' + timeperp + 'climo.nc' #GT, K, 129x64
        fhurrsicp = basepath + 'HURRELL/MODEL_SST.T42_' + timeperp2 + 'climo.nc' #SST, degC, 128x64

        hadsicc = cnc.getNCvar(fhadsic,'GT') # I think, K?
        hadsicp = cnc.getNCvar(fhadsicp,'GT')

    # CanESM2 uses BC files regardless of doBCs flag
    fcansic = basepath + 'CanESM2/GT_BC_CanESM2_historical1979-1989_1870010100-2011020100.nc'
    fcansicp = basepath + 'CanESM2/GTadjusted_BC_CanESM2_historical2002-2012_1870010100-2011020100_abs10thresh.nc'
#    fcansicp = basepath + 'CanESM2/canesm2_kem1pert2_128x64_0001_0120_gt_0000120100-0120120100.nc' # adjusted ssts, shifted one month from above file
    

    cansicc,cansiccstd = cutl.climatologize(cnc.getNCvar(fcansic,'GT'))
    cansicp,cansicpstd = cutl.climatologize(cnc.getNCvar(fcansicp,'GT'))
    
    hurrsicc = cnc.getNCvar(fhurrsic,'SST')+273
    #hurrsicc.resize(hadsicc.shape) # can't resize like this. ??
    #np.append(hurrsicc,hurrsicc[:,:,0])#,axis=2) # could not get append to work
    #hurrsicc[:,:,len(lon)-1] = hurrsicc[:,:,0] # add a wraparound. nope.
    #hurrsicc = np.flipud(hurrsicc)  # the lats are flipped compared to hadisst and nsidc

    
    hurrsicp = cnc.getNCvar(fhurrsicp,'SST')+273
    #hurrsicc.resize(hadsicc.shape) # other datasets have extra lon.
    #hurrsicp[:,:,len(lon)-1] = hurrsicp[:,:,0] # add a wraparound
    #hurrsicp = np.flipud(hurrsicp)  # the lats are flipped compared to hadisst and nsidc

    # MASK out land & inland lakes
    lsmask=con.get_t63landmask()
    lsmask = np.tile(lsmask,(12,1,1))
    ## if flipmask:
    ##     lsmask=np.flipud(lsmask)            
    hadsicc = ma.masked_where(lsmask!=0,hadsicc) # 0 is ocean
    hadsicp = ma.masked_where(lsmask!=0,hadsicp) # 0 is ocean
    hadsicc = hadsicc[:,:,0:-1]
    hadsicp = hadsicp[:,:,0:-1]
    
    #hurrsicc = ma.masked_where(lsmask!=0,hurrsicc) # 0 is ocean # @@ BOGUS
    cansicc = ma.masked_where(lsmask!=0,cansicc) # 0 is ocean
    cansicp = ma.masked_where(lsmask!=0,cansicp) # 0 is ocean
    cansicc = cansicc[:,:,0:-1]
    cansicp = cansicp[:,:,0:-1]
    
    hadsicd = hadsicp-hadsicc
    hurrsicd = hurrsicp-hurrsicc
    cansicd = cansicp-cansicc

lat = cnc.getNCvar(fhadsic,'lat') # just adjust Hurrell data, which is flipped
lon = cnc.getNCvar(fhadsic,'lon') # same for all datasets
lon = lon[0:-1]

if testsic==0:
    


    if had:
        hadsicc = hadsicc[:,:,0:-1]
        hadsicp = hadsicp[:,:,0:-1]
        #hadsiccclimo = hadsiccclimo[:,:,0:-1] # @@ shouldn't need anymore? already tested
        #hadsicpclimo = hadsicpclimo[:,:,0:-1]

        hadsicd = hadsicp-hadsicc
        #hadsicclimod = hadsicpclimo-hadsiccclimo

        plotflds=hadsicd
        dset = 'HadISST'
        title = dset + ' (' + timeperp2 + ')-(' + timeper + ')'
        savstr=timeperp2 + 'min' + timeper
    elif hurr:
        lathurr = cnc.getNCvar(fhurrsic,'lat')
        hurrsicd = hurrsicp-hurrsicc

        plotflds=hurrsicd
        lat=lathurr
        flipmask=1
        dset = 'HURRELL'
        title = dset + ' (' + timeperp2 + ')-(' + timeper + ')'
        savstr=timeperp2 + 'min' + timeper
    elif canesm:
        dset = 'CanESM2'
        
        if plotseacycmag:
            # get ground cover, GC
            #testcase = 'kemctl1'; plotflds = cansicc # ?? @@
            #testcase = 'kem1pert2'; plotflds = cansicp
            testcase = 'kem1pert1'; plotflds = cansicc # same SSTs as control, will use the GC from this run
            
            filen = '/home/rkm/work/DATA/CanAM4/' + testcase + '/ts/' + testcase + '_gc_001-111_ts.nc'
            gcts = cnc.getNCvar(filen,'GC',timesel='0002-01-01,0111-12-31')
            filen = '/home/rkm/work/DATA/CanAM4/' + testcase + '/ts/' + testcase + '_gc_001-111_climo.nc'
            gcclimo = cnc.getNCvar(filen,'GC') # -1 land/sea ice>.15, 0 sea
            # to test for sea (0), have to use np.isclose() with tolerance 1e-9
            # because the numbers are actually very close to zero but not exactly.
            gctsoneyr = gcts[0:12,:,:-1]
            
            cminscm = -40; cmaxscm = 40
            
            title = dset + ' ' + testcase + ' SEACYC Magnitude (' + timeper + ')'
            savstr=timeper + 'seacycmag'
        else:
            plotflds=cansicd   
            title = dset + ' (' + timeperp + ')-(' + timeper + ')'
            savstr=timeperp + 'min' + timeper             
    else:
        #nsidc
        if sicn==0:
            print 'No SST or SIC dataset for NSIDC!'
            exit()

        plotflds=nsidcsicd
        dset = 'NSIDCbootstrap'
        title = ' (' + timeperp2 + ')-(' + timeper + ')'
        savstr=timeperp2 + 'min' + timeper


    cmlen = float(15)

    if plotseacycmag: # only for CanESM2 SST
        maxfld = np.max(plotflds,axis=0)
        maxidx = np.argmax(plotflds,axis=0)
        minfld = np.min(plotflds,axis=0)
        minfld2 = np.min(ma.masked_where(gctsoneyr==1,plotflds),axis=0)
        minidx = np.argmin(plotflds,axis=0)

        janfld = plotflds[0,:,:]
        julfld = plotflds[6,:,:]
        seacycmon = julfld - janfld

        seacycmag = maxfld - minfld
        # PLOT SEASONAL CYCLE MAG as a map
        figa,axs = plt.subplots(1,2)
        bm,pc = cplt.kemmap(seacycmon,lat,lon,cmin=cminscm,cmax=cmaxscm,cmap=cmap,
                                type=mtype,axis=axs[0],suppcb=1,lmask=1,flipmask=flipmask)
        axs[0].set_title('Jul - Jan')

        bm,pc = cplt.kemmap(seacycmag,lat,lon,cmin=cminscm,cmax=cmaxscm,cmap=cmap,
                                type=mtype,axis=axs[1],suppcb=1,lmask=1,flipmask=flipmask)
        axs[1].set_title('Max - Min')
        
        cbar_ax = figa.add_axes([.91,.25,.02,.5])
        figa.colorbar(pc,cax=cbar_ax)

        lons, lats = np.meshgrid(lon,lat)
        ## # PLOT MINIMUM GT: shouldn't get below freezing in the sea (GC=0)
        ## figb,axs = plt.subplots(1,1)
        ## bm,pc = cplt.kemmap(minfld2,lat,lon,cmin=238,cmax=308,cmap='blue2red_20',
        ##                         type=mtype,axis=axs,lmask=1,flipmask=flipmask)
        ## bm.contour(lons,lats,minfld2,levels=[273, 273],colors='k',linewidths='2',latlon=True)
        ## axs.set_title('Min GT')
        
        plotfldsm = ma.masked_where(np.logical_or(gctsoneyr==-1,gctsoneyr==1),plotflds)

        # PLOT MAP FOR ALL MONTHS
        months = con.get_mon()
        midx=0

        fig,axx = plt.subplots(2,6)
        fig.set_size_inches(12,6)
        fig.subplots_adjust(hspace=.05,wspace=.05)

        # These figures look good: where colder than 271.2, the surface isn't SEA (0), for kemctl1
        for ax in axx.flat:

            plotfld = plotfldsm[midx,:,:];
            bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=236.2,cmax=306.2,cmap='blue2red_20',
                                type=mtype,axis=ax,suppcb=1,lmask=1,flipmask=flipmask)
            bm.contour(lons,lats,plotfld,levels=[271.2, 271.2],colors='k',linewidths='2',latlon=True)
            ax.set_title(months[midx])

            midx=midx+1

        cbar_ax = fig.add_axes([.91,.25,.02,.5])
        fig.colorbar(pc,cax=cbar_ax)
        plt.suptitle(testcase + ' SST masked')
        if printtofile:
            fig.savefig('SSTfreezecont_' + testcase + '_allmos_nh.png')
        

        cellareas = cutl.calc_cellareas(lat,lon)
        
        cellareasall = ma.zeros((12,) + cellareas.shape)
        polcellareas = ma.zeros((12,) + cellareas.shape)
        polplotfldsm = ma.zeros((12,) + cellareas.shape)
        totalarea = np.zeros((12,))
        plotfldsmNth = np.zeros((12,))
        for moidx in range(0,12):
            cellareasall[moidx,:,:] = ma.masked_where(np.logical_or(gctsoneyr[moidx,:,:]==-1,gctsoneyr[moidx,:,:]==1),cellareas)
            polcellareas[moidx,:,:] = ma.masked_where(lats>=latlim,cellareasall[moidx,:,:])
            polplotfldsm[moidx,:,:] = ma.masked_where(lats>=latlim,plotfldsm[moidx,:,:])
            totalarea[moidx] = np.sum(polplotfldsm[moidx,:,:])

            wgts = polcellareas[moidx,:,:]/totalarea[moidx]
            plotfldsmNth[moidx] = np.average(polplotfldsm[moidx,:,:],weights=wgts)
            
        fig = plt.figure()
        plt.plot(np.arange(1,13),plotfldsmNth,'k');
        plt.title(testcase + ' Average SST North of ' + str(latlim) + 'N where GC=0')
        plt.ylim((270,286))
        plt.xlim((1,12))
        if printtofile:
            fig.savefig('SST' + str(latlim) + 'Nseamask_' + testcase + '_seascyc.pdf')
                
    if plotallmos:

        incr = (cmaxm-cminm) / cmlen
        conts = np.arange(cminm,cmaxm+incr,incr)

        months = con.get_mon()
        midx=0

        fig,axx = plt.subplots(2,6)
        fig.set_size_inches(12,6)
        fig.subplots_adjust(hspace=.05,wspace=.05)

        for ax in axx.flat:

            plotfld = plotflds[midx,:,:];

            bm,pc = cplt.kemmap(plotfld,lat,lon,cmin=cminm,cmax=cmaxm,cmap=cmap,
                                type=mtype,axis=ax,suppcb=1,lmask=1,flipmask=flipmask)
            ax.set_title(months[midx])

            midx=midx+1

        cbar_ax = fig.add_axes([.91,.25,.02,.5])
        fig.colorbar(pc,cax=cbar_ax)
        plt.suptitle(title)

        if printtofile:
            if sicn:
                fig.savefig(dset + '_sicn' + bcstr + '_' + savstr + '_allmos_' + mtype + 'smclim.pdf')
            elif sic:
                fig.savefig(dset + '_sic' + bcstr + '_' + savstr + '_allmos_' + mtype + '.pdf')
            else:
                fig.savefig(dset + '_sst'  + bcstr + '_' + savstr + '_allmos_' + mtype + '.pdf')

    if plotseacyc:

        if sicn:
            # calc sea ice area
            # mult fraction by grid cell area & sum
            areas = cutl.calc_cellareas(lat,lon)
            ## plt.figure()
            ## plt.pcolor(areas)
            ## plt.colorbar()

            areas = np.tile(areas,(12,1,1))
            #hadsicc = ma.masked_outside(hadsicc,0,1) #@@ need this? No.
            #hadsicp = ma.masked_outside(hadsicp,0,1)

            hadsiac = hadsicc*areas
            hadsiap = hadsicp*areas
            hadsiacclimo = hadsiccclimo*areas
            hadsiapclimo = hadsicpclimo*areas
            hurrsiac = hurrsicc*areas
            hurrsiap = hurrsicp*areas
            nsidcsiac = nsidcsicc*areas
            nsidcsiap = nsidcsicp*areas
            cansiac = cansicc*areas
            cansiap = cansicp*areas

            ## plt.figure()
            ## plt.pcolor(sia[0,:,:])
            ## plt.colorbar()

            ## plt.figure()
            ## plt.pcolor(hadsicc[0,:,:])
            ## plt.colorbar()

            hadsiacnh = np.sum(np.sum(hadsiac[:,lat>0,:],2),1)
            hadsiapnh = np.sum(np.sum(hadsiap[:,lat>0,:],2),1)
            ## Climos are exactly the same as above! good
            ## hadsiacclimonh = np.sum(np.sum(hadsiacclimo[:,lat>0,:],2),1)
            ## hadsiapclimonh = np.sum(np.sum(hadsiapclimo[:,lat>0,:],2),1)
            hurrsiacnh = np.sum(np.sum(hurrsiac[:,lat>0,:],2),1)
            hurrsiapnh = np.sum(np.sum(hurrsiap[:,lat>0,:],2),1)
            nsidcsiacnh = np.sum(np.sum(nsidcsiac[:,lat>0,:],2),1)
            nsidcsiapnh = np.sum(np.sum(nsidcsiap[:,lat>0,:],2),1)
            cansiacnh = np.sum(np.sum(cansiac[:,lat>0,:],2),1)
            cansiapnh = np.sum(np.sum(cansiap[:,lat>0,:],2),1)

            ## siash = np.sum(np.sum(sia[:,lat<0,:],2),1)
            ## siag = np.sum(np.sum(sia,2),1)

            fig = plt.figure()
            plt.plot(hadsiacnh,'k'); plt.plot(hadsiapnh,'k--')
            #plt.plot(hurrsiacnh,'b'); plt.plot(hurrsiapnh,'b--') # SCREWY
            plt.plot(nsidcsiacnh,'r'); plt.plot(nsidcsiapnh,'r--')
            plt.plot(cansiacnh,'g'); plt.plot(cansiapnh,'g--')
            plt.plot(hadsiacclimonh,'c'); plt.plot(hadsiapclimonh,'c--')
            plt.legend(('HadISST 1979-89','HadISST 2002-11',
                        'NSIDCbt 1979-89','NSIDCbt 2002-11',
                        'CanESM2 1979-89','CanESM2 2002-12'),'lower left')
            plt.title('Arctic SICN')
            if printtofile:
                fig.savefig('SICNNH_seascyc_OBS.pdf')

            fig = plt.figure();
            plt.plot(hadsiapnh-hadsiacnh,'k')
            plt.plot(nsidcsiapnh-nsidcsiacnh,'r')
            plt.plot(cansiapnh-cansiacnh,'g')
            plt.legend(('HadISST DIFF','NSIDCbt DIFF','CanESM2 DIFF'),'lower left')
            plt.title('difference in Arctic SICN (2002-11/12 vs 1979-89)')
            if printtofile:
                fig.savefig('SICNNHDIFF_seascyc_OBS.pdf')

            fig = plt.figure();
            plt.plot((hadsiapnh-hadsiacnh)/hadsiacnh*100,'k')
            plt.plot((nsidcsiapnh-nsidcsiacnh)/nsidcsiacnh*100,'r')
            plt.plot((cansiapnh-cansiacnh)/cansiacnh*100,'g')
            plt.legend(('HadISST DIFF','NSIDCbt DIFF','CanESM2 DIFF'),'lower left')
            plt.title('% difference in Arctic SICN (2002-11/12 vs 1979-89)')
            if printtofile:
                fig.savefig('SICNNHDIFFpct_seascyc_OBS.pdf')    
        elif sic:

            # calculate average NH sea ice thickness for each month
            areas = cutl.calc_cellareas(lat,lon)
            areas = np.tile(areas,(12,1,1))

            # only NH sic
            hadsiccnh = hadsicc[:,lat>0,:]
            hadsicpnh = hadsicp[:,lat>0,:]

            # need sea ice area too
            hadsicnc = hadsicnc[:,:,0:-1]
            hadsicnp = hadsicnp[:,:,0:-1]
            hadsiac = hadsicnc*areas
            hadsiap = hadsicnp*areas
            # only NH
            hadsiacnh = hadsiac[:,lat>0,:]
            hadsiapnh = hadsiap[:,lat>0,:]  
            # total sea ice area for area weighting          
            hadtotsiacnh = np.sum(np.sum(hadsiacnh,2),1) # total for each month
            hadtotsiapnh = np.sum(np.sum(hadsiapnh,2),1)
            # tile totalsia to shape of grid
            hadtotsiacnh = np.tile(hadtotsiacnh,(hadsiacnh.shape[1],hadsiacnh.shape[2],1))
            hadtotsiacnh = np.transpose(hadtotsiacnh,(2,0,1))
            hadtotsiapnh = np.tile(hadtotsiapnh,(hadsiapnh.shape[1],hadsiapnh.shape[2],1))
            hadtotsiapnh = np.transpose(hadtotsiapnh,(2,0,1))

            # average thickness
            hadsiccnhavg = np.sum(np.sum(hadsiccnh/913*(hadsiacnh/hadtotsiacnh),2),1)
            hadsicpnhavg = np.sum(np.sum(hadsicpnh/913*(hadsiapnh/hadtotsiapnh),2),1)

            # total volume: thickness * SIA
            hadsivcnh = np.sum(np.sum(hadsiccnh/913 * hadsiacnh,2),1)
            hadsivpnh = np.sum(np.sum(hadsicpnh/913 * hadsiapnh,2),1)

            # do CanESM
            cansiac = cansicnc*areas
            cansiap = cansicnp*areas
            # only NH
            cansiccnh = cansicc[:,lat>0,:]
            cansicpnh = cansicp[:,lat>0,:]
            cansiacnh = cansiac[:,lat>0,:]
            cansiapnh = cansiap[:,lat>0,:]  
            # total sea ice area for area weighting          
            cantotsiacnh = np.sum(np.sum(cansiacnh,2),1) # total for each month
            cantotsiapnh = np.sum(np.sum(cansiapnh,2),1)
            # tile totalsia to shape of grid
            cantotsiacnh = np.tile(cantotsiacnh,(cansiacnh.shape[1],cansiacnh.shape[2],1))
            cantotsiacnh = np.transpose(cantotsiacnh,(2,0,1))
            cantotsiapnh = np.tile(cantotsiapnh,(cansiapnh.shape[1],cansiapnh.shape[2],1))
            cantotsiapnh = np.transpose(cantotsiapnh,(2,0,1))

            # average thickness
            cansiccnhavg = np.sum(np.sum(cansiccnh/913*(cansiacnh/cantotsiacnh),2),1)
            cansicpnhavg = np.sum(np.sum(cansicpnh/913*(cansiapnh/cantotsiapnh),2),1)

            # total volume: thickness * SIA
            cansivcnh = np.sum(np.sum(cansiccnh/913 * cansiacnh,2),1)
            cansivpnh = np.sum(np.sum(cansicpnh/913 * cansiapnh,2),1)

            fig = plt.figure()
            plt.plot(hadsiccnhavg,'k'); plt.plot(hadsicpnhavg,'k--')
            plt.plot(cansiccnhavg,'g'); plt.plot(cansicpnhavg,'g--')
            plt.legend(('HadISST 1979-89','HadISST 2002-11',
                        'CanESM2 1979-89','CanESM2 2002-12'),'lower left')
            plt.title('Arctic SIC')
            if printtofile:
                fig.savefig('SICNH_seascyc_OBS.pdf')

            fig = plt.figure()
            plt.plot(hadsivcnh,'k'); plt.plot(hadsivpnh,'k--')
            plt.plot(cansivcnh,'g'); plt.plot(cansivpnh,'g--')
            plt.legend(('HadISST 1979-89','HadISST 2002-11',
                        'CanESM2 1979-89','CanESM2 2002-12'),'lower left')
            plt.title('Total Arctic Sea Ice Volume')
            if printtofile:
                fig.savefig('SICNHTotVol_seascyc_OBS.pdf')    

        else:

            # calc average SST near sea-ice (where exactly?)
            hadsstc60N = cutl.polar_mean_areawgted3d(hadsicc,lat,lon,latlim=latlim)
            hadsstp60N = cutl.polar_mean_areawgted3d(hadsicp,lat,lon,latlim=latlim)
            hurrsstc60N = cutl.polar_mean_areawgted3d(hurrsicc,lat,lon,latlim=latlim)
            hurrsstp60N = cutl.polar_mean_areawgted3d(hurrsicp,lat,lon,latlim=latlim)
            cansstc60N = cutl.polar_mean_areawgted3d(cansicc,lat,lon,latlim=latlim)
            cansstp60N = cutl.polar_mean_areawgted3d(cansicp,lat,lon,latlim=latlim)

            fig = plt.figure()
            plt.plot(hadsstc60N,'k'); plt.plot(hadsstp60N,'k--')
            #plt.plot(hurrsstc60N,'b'); plt.plot(hurrsstp60N,'b--') # SCREWY
            plt.plot(cansstc60N,'g'); plt.plot(cansstp60N,'g--')
            plt.legend(('HadISST 1979-89','HadISST 2002-11',
                        'CanESM2 1979-89','CanESM2 2002-12'),'lower center')
            plt.title('Average SST North of ' + str(latlim) + 'N')
            if printtofile:
                fig.savefig('SST' + str(latlim) + 'N_seascyc_OBS.pdf')

            fig = plt.figure()
            plt.plot(hadsstp60N-hadsstc60N,'k');
            #plt.plot(hurrsstp60N-hurrsstc60N,'b');#SCREWY
            plt.plot(cansstp60N-cansstc60N,'g');
            plt.legend(('HadISST DIFF','CanESM2 DIFF'),'upper left')
            plt.title('Average difference in SST North of ' + str(latlim) + 'N')
            if printtofile:
                fig.savefig('SST' + str(latlim) + 'NDIFF_seascyc_OBS.pdf')

if testsic:

    lat = cnc.getNCvar(fhadsicorig,'lat')
    lon = cnc.getNCvar(fhadsicorig,'lon')
    sic = cnc.getNCvar(fhadsicorig,'SIC')
    sicn = cnc.getNCvar(fhadsicnorig,'SICN')

    sicANN = cutl.annualize_monthlyts(sic)
    sicANN=sicANN[:,:,0:-1]
    sicANNnh = sicANN[:,lat>0,:]

    sicnANN = cutl.annualize_monthlyts(sicn)
    sicnANN=sicnANN[:,:,0:-1]
    sicnANNnh = sicnANN[:,lat>0,:]

    areas = cutl.calc_cellareas(lat,lon)
    areas = np.tile(areas,(sicANN.shape[0],1,1))
    areasnh = areas[:,lat>0,0:-1]
    totalareanh=np.sum(np.sum(areasnh,2),1)

    siaANNnh = sicnANNnh*areasnh
    totalsianh = np.sum(np.sum(siaANNnh,2),1)
    totalsianh = np.tile(totalsianh,(sicANNnh.shape[1],sicANNnh.shape[2],1))
    totalsianh = np.transpose(totalsianh,(2,0,1))

    sicANNnhavg = np.sum(np.sum(sicANNnh/913*(areasnh/totalareanh[0]),2),1)
    sicANNnhavg2 = np.sum(np.sum(sicANNnh/913*(siaANNnh/totalsianh),2),1)

    years = range(1870,2013)
    plt.figure()
    #plt.plot(sicANNnhavg)
    plt.plot(years,sicANNnhavg2,'k')
    plt.title('HadISST NH average SIC (m)')
    plt.xlim(1870,2012)
    if printtofile:
        plt.savefig('HadISST_SICnh_timeseries.pdf')
    

                          
