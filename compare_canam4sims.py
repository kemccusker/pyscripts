""" compare_canam4sims.py

    Summarize all simulations. E.g. seasonal cycles, area averages, zonal means
    Start with 2D fields
    3/29/2014

"""
import numpy as np
import numpy.ma as ma
import scipy as sp # scientific python
import scipy.stats
import matplotlib.pyplot as plt
import platform as platform
import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import cccmaNC as cnc
import matplotlib.font_manager as fm

cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)
cnc = reload(cnc)

plt.close("all")
plt.ion()

printtofile=1

plotseacyc=0
plotzonmean=1
plotbar=0 # this just does fluxes

model='CanAM4'
seasons = 'DJF','MAM','JJA','SON'

# @@ should calculate fluxes only NOT on land...
# @@ should calculate snow depth and fraction only ON land...

fontP = fm.FontProperties()
fontP.set_size('small')
legloc = 'upper left'

fluxes = 'hfl','hfs','flg' #,'fsg' Leave off fsg b/c it's the absorption the ground gets. doesn't really say how the atmos changes

# st, gt, pmsl, pcp, hfl, hfs, turb, flg, fsg, fn, pcpn, zn, su, sv (@@later ufs,vfs)
field = 'su'

latlim=70
siglevel=0.05
avgtype = 'gt' + str(latlim) + 'N'


# # # ########### set Simulations #############
# Control run
casename = 'kemctl1'; ctlstr=''
timstr = '001-111'

# Pert run
casenamep1 = 'kem1pert1'  # 2002-2012 sic and sit
casenamep2 = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst
casenamep3 = 'kem1pert3'  # 2002-2012 sic, adjusted sst. control sit
timstrp = '001-111'

#### SET DIFFERENT CONTROL
casename = casenamep3; ctlstr='_ctl' + casename;
print 'CONTROL IS ' + casename


##### set up field info  ############
if field == 'st':
    units = 'K'
    conv = 1  # no conversion
    legloc = 'upper left'
    lims=-.5,3.5
    limsp=-.6,.6
elif field == 'gt':
    units = 'K'
    conv = 1  # no conversion
    legloc = 'upper left'
    lims=-1,3.5
    limsp=-.6,.6
elif field == 'pmsl':
    units = 'hPa' # pretty sure hpa @@double check
    conv = 1
    legloc='upper right'
    lims=-1,.8
    limsp=-1,1
elif field == 'pcp':
    pct=1; units = '%'
#    units = 'mm/day' # original: kg m-2 s-1
    conv = 86400  # convert from kg m-2 s-1 to mm/day
    lims=-.1,.1
    limsp=lims
elif field == 'hfl': # sfc upward LH flux
    units = 'W/m2'
    conv = 1
    lims=-1,2.5
    limsp=-1,1
    print 'HFL: positive atmosphere gains heat'
    print '@@ should hfl be masked to not include land?'
elif field == 'hfs': # sfc upward SH flux
    units = 'W/m2'
    conv = 1
    lims=-2,5
    limsp=-1,1
    print 'HFS: positive atmosphere gains heat'
    print '@@ should hfs be masked to not include land?'
elif field == 'turb': # combine hfl and hfs
    units = 'W/m2'
    conv = 1
    lims=-2.5,6.5
    limsp = -1.5,1.5 # for when a pert is the 'ctl'
    print 'TURB: positive atmosphere gains heat'
    print '@@ should turb be masked to not include land?'
elif field == 'net': # net of all sfc fluxes
    print " 'net ' not yet implemented! @@"
    
elif field == 'flg': # net downward LW at the sfc. Positive down?
    units = 'W/m2'
    conv = 1
    lims=-3.5,1
    limsp=-1.5,1.5
    print 'FLG: net downward LW. positive down. Negative means atmos gains heat'
elif field == 'fsg': # net (absorbed) solar downard at sfc
    units = 'W/m2'
    conv = 1
    lims=-1,7
    limsp=-1.5,1.5
    print 'FSG: net (absorbed) solar down. positive down. positive means sfc gains/absorbs more heat'
elif field == 'fn': # snow fraction
    units = '%'
    conv=100
    legloc='lower left'
    avgtype = 'gt' + str(latlim) + 'Nmasked'
    lims=-.5,.3
    limsp=-.5,.5
elif field == 'pcpn': # snowfall rate (water equivalent, kg/m2/s)
    #pct = 1; units='%'
    units = 'mm/day'
    conv = 86400  # convert from kg m-2 s-1 to mm/day (I think it's same as pcp) @@
    lims=-.03,.07
    limsp=-.03,.04
elif field == 'zn': # snow depth (m)
#    pct=1; units='%'
    units = 'cm'
    conv = 100; # convert to cm
    avgtype = 'gt' + str(latlim) + 'Nmasked'
    legloc='upper right'
    lims=-.3,.5
    limsp=-.3,.5
elif field == 'su':
    units = 'm/s'
    conv = 1;
    lims = -.3,.4
    limsp=-.4,.4
elif field == 'sv':
    units = 'm/s'
    conv = 1;
    legloc='lower right'
else:
    print 'No settings for ' + field

plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'


##### set up file info and read data #######
if field == 'turb':
    field='hfl'; fieldb='hfs'
    fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
    fnamep1 = basepath + casenamep1 + subdir + casenamep1 + '_' + field + '_' + timstrp + '_ts.nc'
    fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp + '_ts.nc'
    fnamep3 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp + '_ts.nc'
    fnamecb = basepath + casename + subdir + casename + '_' + fieldb + '_' + timstr + '_ts.nc'
    fnamep1b = basepath + casenamep1 + subdir + casenamep1 + '_' + fieldb + '_' + timstrp + '_ts.nc'
    fnamep2b = basepath + casenamep2 + subdir + casenamep2 + '_' + fieldb + '_' + timstrp + '_ts.nc'
    fnamep3b = basepath + casenamep3 + subdir + casenamep3 + '_' + fieldb + '_' + timstrp + '_ts.nc'

    fldc = cnc.getNCvar(fnamec,field.upper(),timesel='0002-01-01,0111-12-31')*conv + cnc.getNCvar(fnamecb,fieldb.upper(),timesel='0002-01-01,0111-12-31')*conv
    fldp1 = cnc.getNCvar(fnamep1,field.upper(),timesel='0002-01-01,0111-12-31')*conv + cnc.getNCvar(fnamep1b,fieldb.upper(),timesel='0002-01-01,0111-12-31')*conv
    fldp2 = cnc.getNCvar(fnamep2,field.upper(),timesel='0002-01-01,0111-12-31')*conv + cnc.getNCvar(fnamep2b,fieldb.upper(),timesel='0002-01-01,0111-12-31')*conv
    fldp3 = cnc.getNCvar(fnamep3,field.upper(),timesel='0002-01-01,0111-12-31')*conv + cnc.getNCvar(fnamep3b,fieldb.upper(),timesel='0002-01-01,0111-12-31')*conv

    field='turb'
else:
    fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
    fnamep1 = basepath + casenamep1 + subdir + casenamep1 + '_' + field + '_' + timstrp + '_ts.nc'
    fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp + '_ts.nc'
    fnamep3 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp + '_ts.nc'

    fldc = cnc.getNCvar(fnamec,field.upper(),timesel='0002-01-01,0111-12-31')*conv
    fldp1 = cnc.getNCvar(fnamep1,field.upper(),timesel='0002-01-01,0111-12-31')*conv
    fldp2 = cnc.getNCvar(fnamep2,field.upper(),timesel='0002-01-01,0111-12-31')*conv
    fldp3 = cnc.getNCvar(fnamep3,field.upper(),timesel='0002-01-01,0111-12-31')*conv


lat = cnc.getNCvar(fnamec,'lat')
lon = cnc.getNCvar(fnamec,'lon')

fldcclim,fldcstd = cutl.climatologize3d(fldc)
fldp1clim,fldp1std = cutl.climatologize3d(fldp1)
fldp2clim,fldp2std = cutl.climatologize3d(fldp2)
fldp3clim,fldp3std = cutl.climatologize3d(fldp3)

if avgtype == 'gt' + str(latlim) + 'Nmasked':
    print avgtype + '!!!'
    fldc = ma.masked_where(fldc==0,fldc)
    fldp1 = ma.masked_where(fldp1==0,fldp1)
    fldp2 = ma.masked_where(fldp2==0,fldp2)
    fldp3 = ma.masked_where(fldp3==0,fldp3)

fldcpm = cutl.polar_mean_areawgted3d(fldc,lat,lon,latlim=latlim)
fldp1pm = cutl.polar_mean_areawgted3d(fldp1,lat,lon,latlim=latlim)
fldp2pm = cutl.polar_mean_areawgted3d(fldp2,lat,lon,latlim=latlim)
fldp3pm = cutl.polar_mean_areawgted3d(fldp3,lat,lon,latlim=latlim)

# calculate polar mean statistics

tstatp1 = np.arange(1,13,dtype=float)*0
pvalp1 = np.arange(1,13,dtype=float)*0
tstatp2 = np.arange(1,13,dtype=float)*0
pvalp2 = np.arange(1,13,dtype=float)*0
tstatp3 = np.arange(1,13,dtype=float)*0
pvalp3 = np.arange(1,13,dtype=float)*0

for midx in range(0,12):
    
    tstatp1[midx],pvalp1[midx] = sp.stats.ttest_ind(fldp1pm[midx::12],fldcpm[midx::12])
    tstatp2[midx],pvalp2[midx] = sp.stats.ttest_ind(fldp2pm[midx::12],fldcpm[midx::12])
    tstatp3[midx],pvalp3[midx] = sp.stats.ttest_ind(fldp3pm[midx::12],fldcpm[midx::12])



mons = range(1,13)


if plotseacyc:
    
    plt.figure() 
    plt.plot(mons,pvalp1,'b')
    plt.plot(mons,pvalp2,'r')
    plt.plot(mons,pvalp3,'g') 
    plt.plot((1,12),(siglevel,siglevel),'k')
    plt.xlim((1,12))
    plt.title('pval')
    
    # first try just the polar mean
    fldcclimpm = cutl.polar_mean_areawgted3d(fldcclim,lat,lon,latlim=latlim)
    fldp1climpm = cutl.polar_mean_areawgted3d(fldp1clim,lat,lon,latlim=latlim) 
    fldp2climpm = cutl.polar_mean_areawgted3d(fldp2clim,lat,lon,latlim=latlim) 
    fldp3climpm = cutl.polar_mean_areawgted3d(fldp3clim,lat,lon,latlim=latlim)

    # get stddev of the already polar averaged data!
    tmpc,fldcstdpm = cutl.climatologize(fldcpm)
    tmpp1,fldp1stdpm = cutl.climatologize(fldp1pm)
    tmpp2,fldp2stdpm = cutl.climatologize(fldp2pm)
    tmpp3,fldp3stdpm = cutl.climatologize(fldp3pm)
    
    ## fldcstdpm = cutl.polar_mean_areawgted3d(fldcstd,lat,lon,latlim=latlim) # this is wrong metric
    ## fldp1stdpm = cutl.polar_mean_areawgted3d(fldp1std,lat,lon,latlim=latlim) 
    ## fldp2stdpm = cutl.polar_mean_areawgted3d(fldp2std,lat,lon,latlim=latlim) 
    ## fldp3stdpm = cutl.polar_mean_areawgted3d(fldp3std,lat,lon,latlim=latlim)


    # mean values
    fig = plt.figure()
    plt.plot(mons,fldcclimpm,'k',linewidth=2)
    plt.plot(mons,fldp1climpm,'b',linewidth=2)
    plt.plot(mons,fldp2climpm,'r',linewidth=2)
    plt.plot(mons,fldp3climpm,'g',linewidth=2)
    plt.plot(mons,fldcclimpm+fldcstdpm,color='.5')
    plt.plot(mons,fldcclimpm-fldcstdpm,color='.5')
    plt.xlim((1,12))
    plt.xlabel('Month')
    plt.ylabel(field)
    plt.title('Polar mean ' + '(' + avgtype + ')')
    plt.legend(('CTRL','PERT1','PERT2','PERT3','CTRL 1sigma'),legloc,prop=fontP)
    if printtofile:
        fig.savefig(field + '_seascyc_' + avgtype + '_allsims' + ctlstr + '.pdf')

    # mean standard deviations
    fig = plt.figure()
    plt.plot(mons,fldcstdpm,'k')
    plt.plot(mons,fldp1stdpm,'b')
    plt.plot(mons,fldp2stdpm,'r')
    plt.plot(mons,fldp3stdpm,'g')
    plt.xlim((1,12))
    plt.xlabel('Month')
    plt.ylabel(field)
    plt.title('Polar mean  ' + '(' + avgtype + ') standard dev')
    if printtofile:
        if casename == 'kemctl1': # only need to save one since it's not a diff plot
            fig.savefig(field + 'STDDEV_seascyc_' + avgtype + '_allsim.pdf')

    # differences
    fig = plt.figure()
    plt.plot(mons,fldp1climpm-fldcclimpm,'b',linewidth=2)
    plt.plot(mons,fldp2climpm-fldcclimpm,'r',linewidth=2)
    plt.plot(mons,fldp3climpm-fldcclimpm,'g',linewidth=2)
    plt.plot(ma.masked_where(pvalp1>siglevel,mons),ma.masked_where(
        pvalp1>siglevel,fldp1climpm-fldcclimpm),linestyle='none',
             marker='s',color='b',mec='none',markersize=8)
    plt.plot(ma.masked_where(pvalp2>siglevel,mons),ma.masked_where(
        pvalp2>siglevel,fldp2climpm-fldcclimpm),linestyle='none',
             marker='s',color='r',mec='none',markersize=8)
    plt.plot(ma.masked_where(pvalp3>siglevel,mons),ma.masked_where(
        pvalp3>siglevel,fldp3climpm-fldcclimpm),linestyle='none',
             marker='s',color='g',mec='none',markersize=8)
    plt.plot((1,12),(0,0),'k')
    plt.xlim((1,12))
    plt.xlabel('Month')
    plt.ylabel(field)
    plt.title('Polar mean diffs ' + '(' + avgtype + ')')
    if casename=='kem1pert3' and field=='st':
        legloc='upper right'
    elif (casename=='kem1pert1' or casename=='kem1pert2') and field=='turb':
        legloc='lower left'
    plt.legend(('PERT1','PERT2','PERT3','Sig at 0.05'),legloc,prop=fontP)
    if printtofile:
        fig.savefig(field + 'DIFF_seascyc_' + avgtype + '_allsims' + ctlstr + '.pdf')

    fig = plt.figure()
    plt.plot(mons,fldp1stdpm-fldcstdpm,'b')
    plt.plot(mons,fldp2stdpm-fldcstdpm,'r')
    plt.plot(mons,fldp3stdpm-fldcstdpm,'g')
    plt.plot((1,12),(0,0),'k')
    plt.xlim((1,12))
    plt.xlabel('Month')
    plt.ylabel(field)
    plt.title('Polar mean standard dev diffs')


if plotzonmean:

    if casename != 'kemctl1':
        lims = limsp
        
    fldczmsea = np.zeros((len(seasons),len(lat)))
    fldp1zmsea = np.zeros((len(seasons),len(lat)))
    fldp2zmsea = np.zeros((len(seasons),len(lat)))
    fldp3zmsea = np.zeros((len(seasons),len(lat)))

    tstatp1zmsea = np.zeros((len(seasons),len(lat)))
    pvalp1zmsea = np.zeros((len(seasons),len(lat)))
    tstatp2zmsea = np.zeros((len(seasons),len(lat)))
    pvalp2zmsea = np.zeros((len(seasons),len(lat)))
    tstatp3zmsea = np.zeros((len(seasons),len(lat)))
    pvalp3zmsea = np.zeros((len(seasons),len(lat)))
    
    # seasonalize first, then zonal mean
    for sidx in range(0,len(seasons)):
        fldczm  = np.mean(cutl.seasonalize_monthlyts(fldc,season=seasons[sidx]), axis=2) # zm with time dim
        fldp1zm = np.mean(cutl.seasonalize_monthlyts(fldp1,season=seasons[sidx]), axis=2)
        fldp2zm = np.mean(cutl.seasonalize_monthlyts(fldp2,season=seasons[sidx]), axis=2)
        fldp3zm = np.mean(cutl.seasonalize_monthlyts(fldp3,season=seasons[sidx]), axis=2)

        tstatp1zm,pvalp1zm = sp.stats.ttest_ind(fldp1zm,fldczm,axis=0)
        tstatp2zm,pvalp2zm = sp.stats.ttest_ind(fldp2zm,fldczm,axis=0)
        tstatp3zm,pvalp3zm = sp.stats.ttest_ind(fldp3zm,fldczm,axis=0)

        # now get the seasonal average (season x lat)
        fldczmsea[sidx,:] = np.mean(fldczm,axis=0)
        fldp1zmsea[sidx,:] = np.mean(fldp1zm,axis=0)
        fldp2zmsea[sidx,:] = np.mean(fldp2zm,axis=0)
        fldp3zmsea[sidx,:] = np.mean(fldp3zm,axis=0)

        pvalp1zmsea[sidx,:] = pvalp1zm
        pvalp2zmsea[sidx,:] = pvalp2zm
        pvalp3zmsea[sidx,:] = pvalp3zm
        

    plt.figure()
    plt.plot(lat,fldczmsea[0,:],'k')
    plt.plot(lat,fldczmsea[1,:],'g')
    plt.plot(lat,fldczmsea[2,:],'r')
    plt.plot(lat,fldczmsea[3,:],'b')
    plt.xlim(-90,90)
    plt.xlabel('lat')

    fldp1diff = fldp1zmsea-fldczmsea
    fldp1diffmsk = ma.masked_where(pvalp1zmsea>siglevel,fldp1diff)
    fldp2diff = fldp2zmsea-fldczmsea
    fldp2diffmsk = ma.masked_where(pvalp2zmsea>siglevel,fldp2diff)
    fldp3diff = fldp3zmsea-fldczmsea
    fldp3diffmsk = ma.masked_where(pvalp3zmsea>siglevel,fldp3diff)

    # One case, all seasons
    plt.figure()
    plt.plot(lat,fldp1diff[0,:],'k')
    plt.plot(lat,fldp1diff[1,:],'g')
    plt.plot(lat,fldp1diff[2,:],'r')
    plt.plot(lat,fldp1diff[3,:],'b')
    plt.plot(lat,fldp1diffmsk[0,:],linestyle='none',marker='s',color='k')
    plt.plot(lat,fldp1diffmsk[1,:],linestyle='none',marker='s',color='g')
    plt.plot(lat,fldp1diffmsk[2,:],linestyle='none',marker='s',color='r')
    plt.plot(lat,fldp1diffmsk[3,:],linestyle='none',marker='s',color='b')
    plt.xlim(-90,90)
    plt.xlabel('lat')
    plt.title(field + ' ' + casenamep1 + ' - ' + casename)


    # all cases, all seasons
    plt.figure()
    plt.plot(lat,fldp1diff[0,:],'k')
    plt.plot(lat,fldp2diff[0,:],'k',linestyle='--')
    plt.plot(lat,fldp3diff[0,:],'k',linestyle=':')
    plt.plot(lat,fldp1diffmsk[0,:],linestyle='none',marker='s',mec='none',color='k')
    plt.plot(lat,fldp2diffmsk[0,:],linestyle='none',marker='s',mec='none',color='k')
    plt.plot(lat,fldp3diffmsk[0,:],linestyle='none',marker='s',mec='none',color='k')

    plt.plot(lat,fldp1diff[1,:],'g')
    plt.plot(lat,fldp2diff[1,:],'g',linestyle='--')
    plt.plot(lat,fldp3diff[1,:],'g',linestyle=':')
    plt.plot(lat,fldp1diffmsk[1,:],linestyle='none',marker='s',mec='none',color='g')
    plt.plot(lat,fldp2diffmsk[1,:],linestyle='none',marker='s',mec='none',color='g')
    plt.plot(lat,fldp3diffmsk[1,:],linestyle='none',marker='s',mec='none',color='g')

    plt.plot(lat,fldp1diff[2,:],'r')
    plt.plot(lat,fldp2diff[2,:],'r',linestyle='--')
    plt.plot(lat,fldp3diff[2,:],'r',linestyle=':')
    plt.plot(lat,fldp1diffmsk[2,:],linestyle='none',marker='s',mec='none',color='r')
    plt.plot(lat,fldp2diffmsk[2,:],linestyle='none',marker='s',mec='none',color='r')
    plt.plot(lat,fldp3diffmsk[2,:],linestyle='none',marker='s',mec='none',color='r')

    plt.plot(lat,fldp1diff[3,:],'b')
    plt.plot(lat,fldp2diff[3,:],'b',linestyle='--')
    plt.plot(lat,fldp3diff[3,:],'b',linestyle=':')
    plt.plot(lat,fldp1diffmsk[3,:],linestyle='none',marker='s',mec='none',color='b')
    plt.plot(lat,fldp2diffmsk[3,:],linestyle='none',marker='s',mec='none',color='b')
    plt.plot(lat,fldp3diffmsk[3,:],linestyle='none',marker='s',mec='none',color='b')
    
    plt.xlim(-90,90)
    plt.xlabel('lat')
    plt.title(field)
    if printtofile:
        plt.savefig(field + 'ZMDIFF_allsims_allseas' + ctlstr + '.pdf')


    # now do each season as a separate subplot
    # all cases, all seasons
    # this isn't great b/c can't see the individual sims
    # @@ might be better to do sims in their own panels w/ all seasons
    seacolors = 'k','g','r','b'
    
    sidx=0
    fig,axs = plt.subplots(1,4) 
    fig.set_size_inches(12,3)
    fig.subplots_adjust(hspace=.15,wspace=.2)

    for ax in axs.flat:
        ax.plot(fldp1diff[sidx,:],lat,seacolors[sidx])
        ax.plot(fldp2diff[sidx,:],lat,seacolors[sidx],linestyle='--')
        ax.plot(fldp3diff[sidx,:],lat,seacolors[sidx],linestyle=':')
        ax.plot(fldp1diffmsk[sidx,:],lat,linestyle='none',marker='s',mec='none',color=seacolors[sidx])
        ax.plot(fldp2diffmsk[sidx,:],lat,linestyle='none',marker='s',mec='none',color=seacolors[sidx])
        ax.plot(fldp3diffmsk[sidx,:],lat,linestyle='none',marker='s',mec='none',color=seacolors[sidx])
        ax.set_ylim(-90,90)
        ax.set_ylabel('lat')
        ax.set_xlim(lims)
        ax.set_title(seasons[sidx])
        sidx=sidx+1

    fig.suptitle(field)
    if printtofile:
        fig.savefig(field + 'ZMDIFF_allsims_allseas_sbplt' + ctlstr + '.pdf')


    # Divided by sims =======================
    #printtofile=0
    sidx=0
    fig,axs = plt.subplots(3,1) 
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    ax=axs[0]
    # shade where significant
    for sidx in range(0,len(seasons)):
        ax.plot(lat,fldp1diff[sidx,:],seacolors[sidx])
        ax.plot(lat,fldp1diffmsk[sidx,:],linestyle='none',marker='s',mec='none',color=seacolors[sidx])
        ax.fill_between(lat,lims[1],lims[0],where=~fldp1diffmsk[sidx,:].mask,facecolor=seacolors[sidx],alpha=0.2)
    ax.axhline(y=0,color='.5')
    ax.set_xlim(-90,90)
    #ax.set_xlabel('lat')
    ax.set_ylim(lims)
    ax.set_title(casenamep1)

    sidx=0
    ax=axs[1]
    for sidx in range(0,len(seasons)):
        ax.plot(lat,fldp2diff[sidx,:],seacolors[sidx])
        ax.plot(lat,fldp2diffmsk[sidx,:],linestyle='none',marker='s',mec='none',color=seacolors[sidx])
        ax.fill_between(lat,lims[1],lims[0],where=~fldp2diffmsk[sidx,:].mask,facecolor=seacolors[sidx],alpha=0.2)
    ax.axhline(y=0,color='.5')
    ax.set_xlim(-90,90)
    #ax.set_xlabel('lat')
    ax.set_ylabel(field + ' (' + units + ')')
    ax.set_ylim(lims)
    ax.set_title(casenamep2)

    sidx=0
    ax=axs[2]
    for sidx in range(0,len(seasons)):
        ax.plot(lat,fldp3diff[sidx,:],seacolors[sidx])
        ax.plot(lat,fldp3diffmsk[sidx,:],linestyle='none',marker='s',mec='none',color=seacolors[sidx])
        ax.fill_between(lat,lims[1],lims[0],where=~fldp3diffmsk[sidx,:].mask,facecolor=seacolors[sidx],alpha=0.2)
    ax.axhline(y=0,color='.5')
    ax.set_xlim(-90,90)
    ax.set_xlabel('lat')
    ax.set_ylim(lims)
    ax.set_title(casenamep3)

    if printtofile:
        fig.savefig(field + 'ZMDIFF_allsims_allseas_shade' + ctlstr + '.pdf')


    # Divided by sims =======================
    #    JUST NH and JUST FALL/WINTER
    #printtofile=0
    sidx=0
    fig,axs = plt.subplots(3,1) 
    fig.set_size_inches(6,12)
    fig.subplots_adjust(hspace=.2,wspace=.2)
    ax=axs[0]
    # shade where significant
    for sidx in 0,3:
        ax.plot(lat,fldp1diff[sidx,:],seacolors[sidx])
        ax.plot(lat,fldp1diffmsk[sidx,:],linestyle='none',marker='s',mec='none',color=seacolors[sidx])
        ax.fill_between(lat,lims[1],lims[0],where=~fldp1diffmsk[sidx,:].mask,facecolor=seacolors[sidx],alpha=0.2)
    ax.axhline(y=0,color='.5')
    ax.set_xlim(0,90)
    #ax.set_xlabel('lat')
    ax.set_ylim(lims)
    ax.set_title(casenamep1)

    sidx=0
    ax=axs[1]
    for sidx in 0,3:
        ax.plot(lat,fldp2diff[sidx,:],seacolors[sidx])
        ax.plot(lat,fldp2diffmsk[sidx,:],linestyle='none',marker='s',mec='none',color=seacolors[sidx])
        ax.fill_between(lat,lims[1],lims[0],where=~fldp2diffmsk[sidx,:].mask,facecolor=seacolors[sidx],alpha=0.2)
    ax.axhline(y=0,color='.5')
    ax.set_xlim(0,90)
    #ax.set_xlabel('lat')
    ax.set_ylabel(field + ' (' + units + ')')
    ax.set_ylim(lims)
    ax.set_title(casenamep2)

    sidx=0
    ax=axs[2]
    for sidx in 0,3:
        ax.plot(lat,fldp3diff[sidx,:],seacolors[sidx])
        ax.plot(lat,fldp3diffmsk[sidx,:],linestyle='none',marker='s',mec='none',color=seacolors[sidx])
        ax.fill_between(lat,lims[1],lims[0],where=~fldp3diffmsk[sidx,:].mask,facecolor=seacolors[sidx],alpha=0.2)
    ax.axhline(y=0,color='.5')
    ax.set_xlim(0,90)
    ax.set_xlabel('lat')
    ax.set_ylim(lims)
    ax.set_title(casenamep3)

    if printtofile:
        fig.savefig(field + 'ZMDIFF_allsims_AUTWIN_shade_NH' + ctlstr + '.pdf')


# ===============================================================
if plotbar:
    # bar graphs of averaged values for all sims.
    # most useful for Fluxes?

    avgtype = 'gt' + str(latlim) + 'N' # ?@@ which is best? HM, it is exactly the same as 'masked'. I guess it's becasue there are actually values on land

    fldcdict = {}
    fldp1dict = {}
    fldp2dict = {}
    fldp3dict = {}
    fldp1pvdict = {}
    fldp2pvdict = {}
    fldp3pvdict = {}
    fldcstddict = {}
    fldp1stddict = {}
    fldp2stddict = {}
    fldp3stddict = {}
    
    for field in fluxes:
        ctlstr=''
        # make PERT1,2,3 the 'control' by uncommenting next line
        #casename = 'kem1pert3'; ctlstr='_ctl' + casename
        
        fnamec = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'
        fnamep1 = basepath + casenamep1 + subdir + casenamep1 + '_' + field + '_' + timstrp + '_ts.nc'
        fnamep2 = basepath + casenamep2 + subdir + casenamep2 + '_' + field + '_' + timstrp + '_ts.nc'
        fnamep3 = basepath + casenamep3 + subdir + casenamep3 + '_' + field + '_' + timstrp + '_ts.nc'

        fldc = cnc.getNCvar(fnamec,field.upper(),timesel='0002-01-01,0111-12-31')
        fldp1 = cnc.getNCvar(fnamep1,field.upper(),timesel='0002-01-01,0111-12-31')
        fldp2 = cnc.getNCvar(fnamep2,field.upper(),timesel='0002-01-01,0111-12-31')
        fldp3 = cnc.getNCvar(fnamep3,field.upper(),timesel='0002-01-01,0111-12-31')

        if avgtype == 'gt' + str(latlim) + 'Nmasked':
            # @@ switch to land mask for fluxes??
            print avgtype + '!!!'
            fldc = ma.masked_where(fldc==0,fldc)
            fldp1 = ma.masked_where(fldp1==0,fldp1)
            fldp2 = ma.masked_where(fldp2==0,fldp2)
            fldp3 = ma.masked_where(fldp3==0,fldp3)

        ## fldcpm = cutl.polar_mean_areawgted3d(fldc,lat,lon,latlim=latlim)
        ## fldp1pm = cutl.polar_mean_areawgted3d(fldp1,lat,lon,latlim=latlim)
        ## fldp2pm = cutl.polar_mean_areawgted3d(fldp2,lat,lon,latlim=latlim)
        ## fldp3pm = cutl.polar_mean_areawgted3d(fldp3,lat,lon,latlim=latlim)

        # now for each season, calc mean and std dev
        fldcpmseas = {}
        fldp1pmseas = {}
        fldp2pmseas = {}
        fldp3pmseas = {}
        fldp1pmseaspval = {}
        fldp2pmseaspval = {}
        fldp3pmseaspval = {}
        fldcpmseasstd = {}
        fldp1pmseasstd = {}
        fldp2pmseasstd = {}
        fldp3pmseasstd = {}
        for seas in seasons:

            fldcpmseasts = cutl.polar_mean_areawgted3d(cutl.seasonalize_monthlyts(fldc,season=seas),lat,lon,latlim=latlim)
            fldp1pmseasts = cutl.polar_mean_areawgted3d(cutl.seasonalize_monthlyts(fldp1,season=seas),lat,lon,latlim=latlim)
            fldp2pmseasts = cutl.polar_mean_areawgted3d(cutl.seasonalize_monthlyts(fldp2,season=seas),lat,lon,latlim=latlim)
            fldp3pmseasts = cutl.polar_mean_areawgted3d(cutl.seasonalize_monthlyts(fldp3,season=seas),lat,lon,latlim=latlim)

            ttestp1,pvalp1 = sp.stats.ttest_ind(fldp1pmseasts,fldcpmseasts,axis=0)
            ttestp2,pvalp2 = sp.stats.ttest_ind(fldp2pmseasts,fldcpmseasts,axis=0)
            ttestp3,pvalp3 = sp.stats.ttest_ind(fldp3pmseasts,fldcpmseasts,axis=0)            
            
            fldp1pmseaspval[seas] = pvalp1
            fldp2pmseaspval[seas] = pvalp2
            fldp3pmseaspval[seas] = pvalp3

            fldcpmseasstd[seas] = np.std(fldcpmseasts,axis=0)
            fldp1pmseasstd[seas] = np.std(fldp1pmseasts,axis=0)
            fldp2pmseasstd[seas] = np.std(fldp2pmseasts,axis=0)
            fldp3pmseasstd[seas] = np.std(fldp3pmseasts,axis=0)

            fldcpmseas[seas] = np.mean(fldcpmseasts,axis=0)
            fldp1pmseas[seas] = np.mean(fldp1pmseasts,axis=0)
            fldp2pmseas[seas] = np.mean(fldp2pmseasts,axis=0)
            fldp3pmseas[seas] = np.mean(fldp3pmseasts,axis=0)

        
        # set the calc'd values to the dictionary: tree is field --> season
        fldcdict[field] = fldcpmseas
        fldp1dict[field] = fldp1pmseas
        fldp2dict[field] = fldp2pmseas
        fldp3dict[field] = fldp3pmseas

        fldp1pvdict[field] = fldp1pmseaspval
        fldp2pvdict[field] = fldp2pmseaspval
        fldp3pvdict[field] = fldp3pmseaspval

        fldcstddict[field] = fldcpmseasstd
        fldp1stddict[field] = fldp1pmseasstd
        fldp2stddict[field] = fldp2pmseasstd
        fldp3stddict[field] = fldp3pmseasstd

    coldict = {'DJF':'k','MAM':'g','JJA':'r','SON':'b','hfl':'b','hfs':'r','flg':'g','fsg':'o'}
    
    # Now do the figures ====================================
    
    # for each season, print all fluxes
    printtofile=1
    rects1 = {}; rects2 = {}; rects3 = {}
    wi=.2
    N=12 # 4 seasons
    ind = np.arange(0,N,3)
    incr=0
    fig,ax = plt.subplots()
    fig.set_size_inches(9,5)
    for flux in fluxes:
        vals1 = []; vals2 = []; vals3 = []
        pvals1 = []; pvals2 = []; pvals3 = []
        print flux
        for sea in seasons:
            print sea
            if flux=='flg':
                inv=-1 # use this to invert FLG so + is atmos gains heat
            else:
                inv=1
            vals1.append((fldp1dict[flux][sea] - fldcdict[flux][sea])*inv)
            pvals1.append(fldp1pvdict[flux][sea])
            vals2.append((fldp2dict[flux][sea] - fldcdict[flux][sea])*inv)
            pvals2.append(fldp2pvdict[flux][sea])
            vals3.append((fldp3dict[flux][sea] - fldcdict[flux][sea])*inv)
            pvals3.append(fldp3pvdict[flux][sea])
            
        print pvals1    
        sigvals1 = ma.masked_where(np.array(pvals1)>siglevel,np.array(vals1))
        sigvals2 = ma.masked_where(np.array(pvals2)>siglevel,np.array(vals2))
        sigvals3 = ma.masked_where(np.array(pvals3)>siglevel,np.array(vals3))
        
        rects1[flux] = ax.bar(ind+incr,vals1, width=wi,color=coldict[flux],alpha=.4)
        rects2[flux] = ax.bar(ind+1+incr,vals2, width=wi,color=coldict[flux],alpha=.4,hatch='/')
        rects3[flux] = ax.bar(ind+2+incr,vals3, width=wi,color=coldict[flux],alpha=.4,hatch='.')
        # add significance marks
        ax.plot(ind+incr+wi/2,sigvals1,color='.5',linestyle='none',marker='s')
        ax.plot(ind+1+incr+wi/2,sigvals2,color='.5',linestyle='none',marker='s')
        ax.plot(ind+2+incr+wi/2,sigvals3,color='.5',linestyle='none',marker='s')
        incr=incr+wi
        
    ax.set_xticks(ind+wi*6)
    ax.set_xticklabels((seasons))
    ax.axvline(x=ind[0]+wi*14+.1,color='.5')
    ax.axvline(x=ind[1]+wi*14+.1,color='.5')
    ax.axvline(x=ind[2]+wi*14+.1,color='.5')
    ax.legend( (rects1['hfl'][0],rects1['hfs'][0],rects1['flg'][0]),(fluxes),'upper left',prop=fontP)
    ax.axhline(y=0,color='.5')
    if casename=='kemctl1':
        ax.text(ind[0]+wi/2,-.1,casenamep1,rotation=45,fontsize=8)
        ax.text(ind[0]+wi*5,-.1,casenamep2,rotation=45,fontsize=8)
        ax.text(ind[0]+wi*9,-.1,casenamep3,rotation=45,fontsize=8)
        ax.plot(ind[0]+1+wi*2, 2.85,marker='s',color='.5')
        ax.text(ind[0]+1+wi*2+.12, 2.8,'sig at 95%',fontsize=10)
    ax.set_title('Surface Fluxes (' + avgtype + '): + atmosphere gains heat')
    if printtofile:
        fig.savefig('fluxesDIFF_' + avgtype + '_barbyseas' + ctlstr + '.pdf')
