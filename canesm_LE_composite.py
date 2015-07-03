""" taken from canesm_LE_regress.py: July 2, 2015
"""


import loadLE as le
import cccmaplots as cplt
import cccmautils as cutl
import canesm_LE_general as leg
import loadmodeldata as lmd
import pandas as pd
import constants as con

printtofile=True

local=True


timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'
timeselall = '1979-01-01,2012-12-31'

ftype='fullts' # 'fullclimo' or 'climo' or 'fullts'
# this is not implemented here. but would do regressions in time on LE ensemble avg
ensmean=False 

seasp='DJF' # season of spatial field
sear='DJF' # season of regional avgs


# spatial field1 in color
leconvsp=1
fieldsp='tas'; ncfieldsp='tas'; compsp='Amon'; 
cminsp=-3; cmaxsp=3 # for colors
cminspbig=-3; cmaxspbig=3 # for colors

cminspa=-1; cmaxspa=1 # for colors for AGCM

# spatial field2 in contours
leconvsp2=1
fieldsp2='zg50000.00'; ncfieldsp2='zg'; compsp2='Amon'
cminsp2=-30; cmaxsp2=30 # to calc contour interval
cminsp2big=-30; cmaxsp2big=30 # to calc contour interval

cminsp2a=-10; cmaxsp2a=10 # to calc contour interval for AGCM


leconvr=leconvr2=1
# regional avg field 1
#leconvr=-1 # this way, sea ice loss is linked with positive changes elsewhere
fieldr='sic'; ncfieldr='sic'; compr='OImon'; regionr='bksmori'

# regional avg field 2
#leconvr2=-1 # so cooling=high heights
fieldr2='tas'; ncfieldr2='tas'; compr2='Amon'; regionr2='eurasiamori';



fdictsp = {'field': fieldsp, 'ncfield': ncfieldsp, 'comp': compsp}
fdictsp2 = {'field': fieldsp2, 'ncfield': ncfieldsp2, 'comp': compsp2}
fdictr = {'field': fieldr+regionr, 'ncfield': ncfieldr, 'comp': compr}
fdictr2 = {'field': fieldr2+regionr2, 'ncfield': ncfieldr2, 'comp': compr2}

lat=le.get_lat(local=local)
lon=le.get_lon(local=local)
nlat=len(lat); nlon=len(lon)



def load_field(fdict,casename,timesel,seas,ftype='fullts',conv=1,local=False,verb=False):
    
    """ 
        returns [numens x space.flat] or [numens]

    """

    ledat = le.load_LEdata(fdict,casename,timesel=timesel, 
                           rettype='ndarray',conv=conv,ftype=ftype,local=local,verb=verb)

    print '@@@ ledat.shape ' + str(ledat.shape) # why does the 3d data have space flattened already...

    # time needs to be first dimension
    try:
        if ledat.ndim==2:
            ledat = ledat.T
        elif ledat.ndim==3:
            ledat = np.transpose(ledat,(1,0,2))
        else:
            print 'Loaded data is not 2 or 3 dimensions. Do not understand.'
            raise Exception
    except:
        raise

    lesea = cutl.seasonalize_monthlyts(ledat,season=seas).mean(axis=0)  # numens x space.flat

    return lesea


def load_agcmfield(field,sims,seas,conv=1,region=None,subsampyrs=11,styears=None):
    """ loads subsampled agcm data from specified simulations

              number of total samples will be determined by length of all sim data
                      and numyrs (e.g. (ntime / numyrs)*numsims)

              returns subsamp
                      styears
    """
    threed=False

    simconv=1
    if field=='tas': simfield='st'; simncfield='ST'
    elif field=='zg50000.00': simfield='gz50000'; simncfield='PHI'; simconv=1/con.get_g()
    elif field=='sia': simfield='sicn'; simncfield='SICN'; print '@@ danger, sia actually sicn average'
    elif field=='sic': simfield='sicn'; simncfield='SICN'
    else: print 'cannot addsims for ' + field

    
    if region!=None:
        simflddf = pd.DataFrame(lmd.loaddata((simfield,),sims,ncfields=(simncfield,), timefreq=seas, 
                                             region=region))*simconv
    else:
        # assume threed b/c no region given.
        threed=True
        simflddf = lmd.loaddata((simfield,),sims,ncfields=(simncfield,), timefreq=seas, 
                                region=region,rettype='ndarray')*simconv


    subsamp,styearsss = leg.subsamp_sims(simflddf,numyrs=subsampyrs,styears=styears,threed=threed)

    return subsamp,styearsss




#casenames=('historical','historicalNat','historicalMisc')
casenames=('historical',)

for casename in casenames:
    # # SPATIAL DATA
    print 'SPATIAL'

    lecseasp = load_field(fdictsp,casename,timeselc,seasp,ftype=ftype,conv=leconvsp)
    lepseasp = load_field(fdictsp,casename,timeselp,seasp,ftype=ftype,conv=leconvsp)
    leseasp = lepseasp-lecseasp

    lecseasp2 = load_field(fdictsp2,casename,timeselc,seasp,ftype=ftype,conv=leconvsp2)
    lepseasp2 = load_field(fdictsp2,casename,timeselp,seasp,ftype=ftype,conv=leconvsp2)
    leseasp2 = lepseasp2-lecseasp2


    # # 1D DATA
    print '1D'
    lecsear = load_field(fdictr,casename,timeselc,sear,ftype=ftype,conv=leconvr)
    lepsear = load_field(fdictr,casename,timeselp,sear,ftype=ftype,conv=leconvr)
    lesear = (lepsear-lecsear)#/(lepsear-lecsear).std() # don't normalize for the composites

    lecsear2 = load_field(fdictr2,casename,timeselc,sear,ftype=ftype,conv=leconvr2)
    lepsear2 = load_field(fdictr2,casename,timeselp,sear,ftype=ftype,conv=leconvr2)
    lesear2 = (lepsear2-lecsear2)#/(lepsear2-lecsear2).std()


    if casename!='historical': # really just has to be not the first loop thru
        tmp = np.vstack((tmp,leseasp))
        tmp2 = np.vstack((tmp2,leseasp2))

        tmpr = np.hstack((tmpr,lesear))
        tmpr2 = np.hstack((tmpr2,lesear2))
    else:
        tmp = leseasp
        tmp2 = leseasp2
        tmpr = lesear
        tmpr2 = lesear2

leseasp=tmp
leseasp2=tmp2
lesear=tmpr
lesear2=tmpr2

numens = lesear.shape[0]

# =========== composite on BKS SIC ===============
# calc composites: choose top and bottom 10 BKS SIC
nn=10
highidx = (-lesear).argsort()[:nn]
lowidx = lesear.argsort()[:nn]

meanlesp = leseasp.mean(axis=0).reshape((nlat,nlon))
highlesp = leseasp[highidx,...].mean(axis=0).reshape((nlat,nlon))
lowlesp = leseasp[lowidx,...].mean(axis=0).reshape((nlat,nlon))
# low vs high
(lesptstat, lesppval) = cutl.ttest_ind(leseasp[lowidx,...].reshape((nn,nlat,nlon)), 
                                       leseasp[highidx,...].reshape((nn,nlat,nlon)))
# low vs mean
(lelosptstat, lelosppval) = cutl.ttest_ind(leseasp[lowidx,...].reshape((nn,nlat,nlon)), 
                                           leseasp.reshape((numens,nlat,nlon)))
# high vs mean
(lehisptstat, lehisppval) = cutl.ttest_ind(leseasp[highidx,...].reshape((nn,nlat,nlon)), 
                                           leseasp.reshape((numens,nlat,nlon)))


meanlesp2 = leseasp2.mean(axis=0).reshape((nlat,nlon))
highlesp2 = leseasp2[highidx,...].mean(axis=0).reshape((nlat,nlon))
lowlesp2 = leseasp2[lowidx,...].mean(axis=0).reshape((nlat,nlon))
# low vs high
(lesp2tstat, lesp2pval) = cutl.ttest_ind(leseasp2[lowidx,...].reshape((nn,nlat,nlon)), 
                                       leseasp2[highidx,...].reshape((nn,nlat,nlon)))
# low vs mean
(lelosp2tstat, lelosp2pval) = cutl.ttest_ind(leseasp2[lowidx,...].reshape((nn,nlat,nlon)), 
                                           leseasp2.reshape((numens,nlat,nlon)))
# high vs mean
(lehisp2tstat, lehisp2pval) = cutl.ttest_ind(leseasp2[highidx,...].reshape((nn,nlat,nlon)), 
                                           leseasp2.reshape((numens,nlat,nlon)))



# =========== composite on EUR SAT ===============
# calc composites: choose top and bottom 10 Eur SAT
high2idx = (-lesear2).argsort()[:nn]
low2idx = lesear2.argsort()[:nn]

high2lesp = leseasp[high2idx,...].mean(axis=0).reshape((nlat,nlon))
low2lesp = leseasp[low2idx,...].mean(axis=0).reshape((nlat,nlon))
# low vs high
(le2sptstat, le2sppval) = cutl.ttest_ind(leseasp[low2idx,...].reshape((nn,nlat,nlon)), 
                                         leseasp[high2idx,...].reshape((nn,nlat,nlon)))
# low vs mean
(lelo2sptstat, lelo2sppval) = cutl.ttest_ind(leseasp[low2idx,...].reshape((nn,nlat,nlon)), 
                                             leseasp.reshape((numens,nlat,nlon)))
# high vs mean
(lehi2sptstat, lehi2sppval) = cutl.ttest_ind(leseasp[high2idx,...].reshape((nn,nlat,nlon)), 
                                             leseasp.reshape((numens,nlat,nlon)))


high2lesp2 = leseasp2[high2idx,...].mean(axis=0).reshape((nlat,nlon))
low2lesp2 = leseasp2[low2idx,...].mean(axis=0).reshape((nlat,nlon))
# low vs high
(le2sp2tstat, le2sp2pval) = cutl.ttest_ind(leseasp2[low2idx,...].reshape((nn,nlat,nlon)), 
                                           leseasp2[high2idx,...].reshape((nn,nlat,nlon)))
# low vs mean
(lelo2sp2tstat, lelo2sp2pval) = cutl.ttest_ind(leseasp2[low2idx,...].reshape((nn,nlat,nlon)), 
                                               leseasp2.reshape((numens,nlat,nlon)))
# high vs mean
(lehi2sp2tstat, lehi2sp2pval) = cutl.ttest_ind(leseasp2[high2idx,...].reshape((nn,nlat,nlon)), 
                                               leseasp2.reshape((numens,nlat,nlon)))



# === AGCM ==========
alat=con.get_t63lat(); alon=con.get_t63lon()
alons, alats = np.meshgrid(alon,alat)

simsE=('E1','E2','E3','E4','E5'); 
sims=('R1','R2','R3','R4','R5');

asseasp,styears = load_agcmfield(fieldsp,sims,seasp) # already anomalies
asseasp2,styears = load_agcmfield(fieldsp2,sims,seasp,styears=styears) # already anomalies

assear,styears = load_agcmfield(fieldr,sims,sear,styears=styears,region=regionr) # already anomalies
assear2,styears = load_agcmfield(fieldr2,sims,sear,styears=styears,region=regionr2) # already anomalies
assear=assear # / assear.std()
assear2=assear2 #/ assear2.std()

# ===  E sims:
aesseasp,styearse = load_agcmfield(fieldsp,simsE,seasp) 
aesseasp2,styearse = load_agcmfield(fieldsp2,simsE,seasp,styears=styearse) 
aessear2,styearse = load_agcmfield(fieldr2,simsE,sear,styears=styearse,region=regionr2) 
aessear2=aessear2 #/ aessear2.std()


(anens,anlat,anlon)=asseasp.shape
rshape=(anens,anlat*anlon)


# =========== composite on BKS SIC (AGCM) ===============
# calc composites: choose top and bottom BKS SIC (10 subsamps each)

ahighidx = (-assear).argsort()[:nn]
alowidx = assear.argsort()[:nn]

ameansp = asseasp.mean(axis=0)
ahighsp = asseasp[ahighidx,...].mean(axis=0)
alowsp = asseasp[alowidx,...].mean(axis=0)
# low vs high
(asptstat, asppval) = cutl.ttest_ind(asseasp[alowidx,...], 
                                     asseasp[ahighidx,...])
# low vs mean
(alosptstat, alosppval) = cutl.ttest_ind(asseasp[alowidx,...], 
                                         asseasp)
# low vs mean
(ahisptstat, ahisppval) = cutl.ttest_ind(asseasp[ahighidx,...], 
                                         asseasp)

ameansp2 = asseasp2.mean(axis=0)
ahighsp2 = asseasp2[ahighidx,...].mean(axis=0)
alowsp2 = asseasp2[alowidx,...].mean(axis=0)
# low vs high
(asp2tstat, asp2pval) = cutl.ttest_ind(asseasp2[alowidx,...], 
                                       asseasp2[ahighidx,...])
# low vs mean
(alosp2tstat, alosp2pval) = cutl.ttest_ind(asseasp2[alowidx,...], 
                                           asseasp2)
# low vs mean
(ahisp2tstat, ahisp2pval) = cutl.ttest_ind(asseasp2[ahighidx,...], 
                                           asseasp2)


# =========== composite on EUR SAT (AGCM) ===============
# calc composites: choose top and bottom 10 Eur SAT 

ahigh2idx = (-assear2).argsort()[:nn]
alow2idx = assear2.argsort()[:nn]

ahigh2sp = asseasp[ahigh2idx,...].mean(axis=0)
alow2sp = asseasp[alow2idx,...].mean(axis=0)
# low vs high
(a2sptstat, a2sppval) = cutl.ttest_ind(asseasp[alow2idx,...], 
                                     asseasp[ahigh2idx,...])
# low vs mean
(alo2sptstat, alo2sppval) = cutl.ttest_ind(asseasp[alow2idx,...], 
                                           asseasp)
# high vs mean
(ahi2sptstat, ahi2sppval) = cutl.ttest_ind(asseasp[ahigh2idx,...], 
                                           asseasp)


ahigh2sp2 = asseasp2[ahigh2idx,...].mean(axis=0)
alow2sp2 = asseasp2[alow2idx,...].mean(axis=0)
# low vs high
(a2sp2tstat, a2sp2pval) = cutl.ttest_ind(asseasp2[alow2idx,...], 
                                         asseasp2[ahigh2idx,...])
# low vs mean
(alo2sp2tstat, alo2sp2pval) = cutl.ttest_ind(asseasp2[alow2idx,...], 
                                             asseasp2)
# high vs mean
(ahi2sp2tstat, ahi2sp2pval) = cutl.ttest_ind(asseasp2[ahigh2idx,...], 
                                             asseasp2)



# E sims ----
aehigh2idx = (-aessear2).argsort()[:nn]
aelow2idx = aessear2.argsort()[:nn]

aemeansp = aesseasp.mean(axis=0)
aehigh2sp = aesseasp[aehigh2idx,...].mean(axis=0)
aelow2sp = aesseasp[aelow2idx,...].mean(axis=0)
# low vs high
(ae2sptstat, ae2sppval) = cutl.ttest_ind(aesseasp[aelow2idx,...], 
                                         aesseasp[aehigh2idx,...])
# low vs mean
(aelo2sptstat, aelo2sppval) = cutl.ttest_ind(aesseasp[aelow2idx,...], 
                                             aesseasp)
# high vs mean
(aehi2sptstat, aehi2sppval) = cutl.ttest_ind(aesseasp[aehigh2idx,...], 
                                             aesseasp)


aemeansp2 = aesseasp2.mean(axis=0)
aehigh2sp2 = aesseasp2[aehigh2idx,...].mean(axis=0)
aelow2sp2 = aesseasp2[aelow2idx,...].mean(axis=0)
# low vs high
(ae2sp2tstat, ae2sp2pval) = cutl.ttest_ind(aesseasp2[aelow2idx,...], 
                                           aesseasp2[aehigh2idx,...])
# low vs mean
(aelo2sp2tstat, aelo2sp2pval) = cutl.ttest_ind(aesseasp2[aelow2idx,...], 
                                               aesseasp2)
# high vs mean
(aehi2sp2tstat, aehi2sp2pval) = cutl.ttest_ind(aesseasp2[aehigh2idx,...], 
                                               aesseasp2)



# AGCM RESIDUAL pvals (var vs const)
(aressptstat,aressppval) = cutl.ttest_ind(asseasp[alow2idx,...]- asseasp[ahigh2idx,...],
                                          aesseasp[aelow2idx,...]- aesseasp[aehigh2idx,...])
(aressp2tstat,aressp2pval) = cutl.ttest_ind(asseasp2[alow2idx,...]- asseasp2[ahigh2idx,...],
                                          aesseasp2[aelow2idx,...]- aesseasp2[aehigh2idx,...])

# ==============================================
# ====================== FIGURES ===============

lons, lats = np.meshgrid(lon,lat)
cmlen=15.
incr = (cmaxsp2-cminsp2) / (cmlen)
conts = np.arange(cminsp2,cmaxsp2+incr,incr)

incrbig = (cmaxsp2big-cminsp2big) / (cmlen)
contsbig = np.arange(cminsp2big,cmaxsp2big+incrbig,incrbig)


ttl1=seasp + ': High ' + sear + ' BKS SIC'
ttl2=seasp + ': Low ' + sear + ' BKS SIC' 

#ttl1=ttl2=''



# ========== 4 panel fig =========
addsig=True



if addsig:
    prstr='sig'
else:
    prstr=''

fig,axs=plt.subplots(2,3)
fig.set_size_inches(10,8)
fig.subplots_adjust(wspace=0.05,hspace=0.03)

ax=axs[0,0]
bm,pc=cplt.kemmap(highlesp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl1,suppcb=True,
                  panellab='a',lcol='0.2')
bm.contour(lons,lats,highlesp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)
ax.set_ylabel('CGCM')

ax=axs[0,1]
bm,pc=cplt.kemmap(lowlesp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl2,suppcb=True,
                  panellab='b',lcol='0.2')
bm.contour(lons,lats,lowlesp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[0,2]
bm,pc=cplt.kemmap(meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title='Mean Anom',suppcb=True,
                  panellab='c',lcol='0.2')
bm.contour(lons,lats,meanlesp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)
#cplt.add_colorbar(fig,pc,orientation='horizontal',pos=[.25,.5, .5,.03])

ax=axs[1,0]
bm,pc=cplt.kemmap(highlesp-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='High Anom',suppcb=True,
                  panellab='d',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lehisppval,lat,lon)
bm.contour(lons,lats,highlesp2-meanlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lehisp2pval,lat,lon,type='cont',color='g')
ax.set_ylabel('CGCM Diff')

ax=axs[1,1]
bm,pc=cplt.kemmap(lowlesp-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low Anom',suppcb=True,
                  panellab='e',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lelosppval,lat,lon)
bm.contour(lons,lats,lowlesp2-meanlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lelosp2pval,lat,lon,type='cont',color='g')

ax=axs[1,2]
bm,pc=cplt.kemmap(lowlesp-highlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low-High',suppcb=True,
                  panellab='f',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lesppval, lat, lon)
bm.contour(lons,lats,lowlesp2-highlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2pval,lat,lon,type='cont',color='g')

cplt.add_colorbar(fig,pc,orientation='horizontal')
if printtofile:
    if addsig:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                '_compon_' + fieldr+regionr+ sear + '_CGCM' + prstr + '.png', dpi=400)
    else:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr+regionr+ sear + '_CGCM' + prstr + '.pdf') 

# =========== AGCM ====

fig,axs=plt.subplots(2,3)
fig.set_size_inches(10,8)
fig.subplots_adjust(wspace=0.05,hspace=0.03)

ax=axs[0,0]
bm,pc=cplt.kemmap(ahighsp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl1,suppcb=True,
                  panellab='a',lcol='0.2')
bm.contour(alons,alats,ahighsp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)
ax.set_ylabel('AGCM')

ax=axs[0,1]
bm,pc=cplt.kemmap(alowsp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl2,suppcb=True,
                  panellab='b',lcol='0.2')
bm.contour(alons,alats,alowsp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[0,2]
bm,pc=cplt.kemmap(ameansp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title='Mean Anom',suppcb=True,
                  panellab='c',lcol='0.2')
bm.contour(alons,alats,ameansp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)

#cplt.add_colorbar(fig,pc,orientation='horizontal',pos=[.25,.5, .5,.03])

ax=axs[1,0]
bm,pc=cplt.kemmap(ahighsp-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='High Anom',suppcb=True,
                  panellab='d',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,ahisppval,alat,alon)
bm.contour(alons,alats,ahighsp2-ameansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,ahisp2pval,alat,alon,type='cont',color='g')

ax.set_ylabel('AGCM Diff')

ax=axs[1,1]
bm,pc=cplt.kemmap(alowsp-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low Anom',suppcb=True,
                  panellab='e',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,alosppval,alat,alon)
bm.contour(alons,alats,alowsp2-ameansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,alosp2pval,alat,alon,type='cont',color='g')

ax=axs[1,2]
bm,pc=cplt.kemmap(alowsp-ahighsp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low-High',suppcb=True,
                  panellab='f',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,asppval, alat, alon)
bm.contour(alons,alats,alowsp2-ahighsp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,asp2pval,alat,alon,type='cont',color='g')

cplt.add_colorbar(fig,pc,orientation='horizontal')
if printtofile:
    if addsig:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr+regionr+ sear + '_AGCM' + prstr + '.png',dpi=400) 
    else:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr+regionr+ sear + '_AGCM' + prstr + '.pdf') 




# CGCM and AGCM

fig,axs=plt.subplots(2,3)
fig.set_size_inches(10,8)
fig.subplots_adjust(wspace=0.05,hspace=0.03)

ax=axs[0,0]
bm,pc=cplt.kemmap(highlesp-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='High Anom',suppcb=True,
                  panellab='a',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lehisppval,lat,lon)
bm.contour(lons,lats,highlesp2-meanlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lehisp2pval,lat,lon,type='cont',color='g')
ax.set_ylabel('CGCM Diff')

ax=axs[0,1]
bm,pc=cplt.kemmap(lowlesp-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low Anom',suppcb=True,
                  panellab='b',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lelosppval,lat,lon)
bm.contour(lons,lats,lowlesp2-meanlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lelosp2pval,lat,lon,type='cont',color='g')

ax=axs[0,2]
bm,pc=cplt.kemmap(lowlesp-highlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low-High',suppcb=True,
                  panellab='c',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lesppval, lat, lon)
bm.contour(lons,lats,lowlesp2-highlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2pval,lat,lon,type='cont',color='g')

ax=axs[1,0]
bm,pc=cplt.kemmap(ahighsp-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='High Anom',suppcb=True,
                  panellab='d',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,ahisppval,alat,alon)
bm.contour(alons,alats,ahighsp2-ameansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,ahisp2pval,alat,alon,type='cont',color='g')
ax.set_ylabel('AGCM Diff')

ax=axs[1,1]
bm,pc=cplt.kemmap(alowsp-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low Anom',suppcb=True,
                  panellab='e',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,alosppval,alat,alon)
bm.contour(alons,alats,alowsp2-ameansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,alosp2pval,alat,alon,type='cont',color='g')

ax=axs[1,2]
bm,pc=cplt.kemmap(alowsp-ahighsp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low-High',suppcb=True,
                  panellab='f',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,asppval, alat, alon)
bm.contour(alons,alats,alowsp2-ahighsp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,asp2pval,alat,alon,type='cont',color='g')

cplt.add_colorbar(fig,pc,orientation='horizontal')

if printtofile:
    if addsig:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr+regionr+ sear + '_CGCMAGCM' + prstr + '.png',dpi=400) 
    else:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr+regionr+ sear + '_CGCMAGCM' + prstr + '.pdf') 








# ====================================================
# # # 4 panel Fig: composite on SAT

ttl12=seasp + ': High ' + sear + ' Eur SAT'
ttl22=seasp + ': Low ' + sear + ' Eur SAT' 


fig,axs=plt.subplots(2,3)
fig.set_size_inches(10,8)
fig.subplots_adjust(wspace=0.05,hspace=0.03)

ax=axs[0,0]
bm,pc=cplt.kemmap(high2lesp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl12,suppcb=True,
                  panellab='a',lcol='0.2')
bm.contour(lons,lats,high2lesp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)
ax.set_ylabel('CGCM')

ax=axs[0,1]
bm,pc=cplt.kemmap(low2lesp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl22,suppcb=True,
                  panellab='b',lcol='0.2')
bm.contour(lons,lats,low2lesp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[0,2]
bm,pc=cplt.kemmap(meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title='Mean Anom',suppcb=True,
                  panellab='c',lcol='0.2')
bm.contour(lons,lats,meanlesp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[1,0]
bm,pc=cplt.kemmap(high2lesp-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='High Anom',suppcb=True,
                  panellab='d',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lehi2sppval,lat,lon)
bm.contour(lons,lats,high2lesp2-meanlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lehi2sp2pval,lat,lon,type='cont',color='g')
ax.set_ylabel('CGCM Diff')

ax=axs[1,1]
bm,pc=cplt.kemmap(low2lesp-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low Anom',suppcb=True,
                  panellab='e',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lelo2sppval,lat,lon)
bm.contour(lons,lats,low2lesp2-meanlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lelo2sp2pval,lat,lon,type='cont',color='g')

ax=axs[1,2]
bm,pc=cplt.kemmap(low2lesp-high2lesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low-High',suppcb=True,
                  panellab='f',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,le2sppval, lat, lon)
bm.contour(lons,lats,low2lesp2-high2lesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,le2sp2pval,lat,lon,type='cont',color='g')

cplt.add_colorbar(fig,pc,orientation='horizontal')
if printtofile:
    if addsig:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_CGCM' + prstr + '.png',dpi=400) 
    else:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_CGCM' + prstr + '.pdf') 



# ======= AGCM
fig,axs=plt.subplots(2,3)
fig.set_size_inches(10,8)
fig.subplots_adjust(wspace=0.05,hspace=0.03)

ax=axs[0,0]
bm,pc=cplt.kemmap(ahigh2sp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl12,suppcb=True,
                  panellab='a',lcol='0.2')
bm.contour(alons,alats,ahigh2sp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)
ax.set_ylabel('AGCM')

ax=axs[0,1]
bm,pc=cplt.kemmap(alow2sp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl22,suppcb=True,
                  panellab='b',lcol='0.2')
bm.contour(alons,alats,alow2sp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[0,2]
bm,pc=cplt.kemmap(ameansp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title='Mean Anom',suppcb=True,
                  panellab='c',lcol='0.2')
bm.contour(alons,alats,ameansp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[1,0]
bm,pc=cplt.kemmap(ahigh2sp-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='High Anom',suppcb=True,
                  panellab='d',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,ahi2sppval,alat,alon)
bm.contour(alons,alats,ahigh2sp2-ameansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,ahi2sp2pval,alat,alon,type='cont',color='g')
ax.set_ylabel('AGCM Diff')

ax=axs[1,1]
bm,pc=cplt.kemmap(alow2sp-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low Anom',suppcb=True,
                  panellab='e',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,alo2sppval,alat,alon)
bm.contour(alons,alats,alow2sp2-ameansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,alo2sp2pval,alat,alon,type='cont',color='g')

ax=axs[1,2]
bm,pc=cplt.kemmap(alow2sp-ahigh2sp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low-High',suppcb=True,
                  panellab='f',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,a2sppval, alat, alon)
bm.contour(alons,alats,alow2sp2-ahigh2sp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,a2sp2pval,alat,alon,type='cont',color='g')

cplt.add_colorbar(fig,pc,orientation='horizontal')
if printtofile:
    if addsig:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_AGCM' + prstr + '.png',dpi=400) 
    else:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_AGCM' + prstr + '.pdf') 


# ==== AGCM E SIMS !

fig,axs=plt.subplots(2,3)
fig.set_size_inches(10,8)
fig.subplots_adjust(wspace=0.05,hspace=0.03)

ax=axs[0,0]
bm,pc=cplt.kemmap(aehigh2sp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl12,suppcb=True,
                  panellab='a',lcol='0.2')
bm.contour(alons,alats,aehigh2sp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)
ax.set_ylabel('AGCM')

ax=axs[0,1]
bm,pc=cplt.kemmap(aelow2sp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title=ttl22,suppcb=True,
                  panellab='b',lcol='0.2')
bm.contour(alons,alats,aelow2sp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[0,2]
bm,pc=cplt.kemmap(aemeansp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                  title='Mean Anom',suppcb=True,
                  panellab='c',lcol='0.2')
bm.contour(alons,alats,aemeansp2,levels=contsbig,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[1,0]
bm,pc=cplt.kemmap(aehigh2sp-aemeansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='High Anom',suppcb=True,
                  panellab='d',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,aehi2sppval,alat,alon)
bm.contour(alons,alats,aehigh2sp2-aemeansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,aehi2sp2pval,alat,alon,type='cont',color='g')
ax.set_ylabel('AGCM Diff')

ax=axs[1,1]
bm,pc=cplt.kemmap(aelow2sp-aemeansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low Anom',suppcb=True,
                  panellab='e',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,aelo2sppval,alat,alon)
bm.contour(alons,alats,aelow2sp2-aemeansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,aelo2sp2pval,alat,alon,type='cont',color='g')

ax=axs[1,2]
bm,pc=cplt.kemmap(aelow2sp-aehigh2sp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low-High',suppcb=True,
                  panellab='f',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,ae2sppval, alat, alon)
bm.contour(alons,alats,aelow2sp2-aehigh2sp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,ae2sp2pval,alat,alon,type='cont',color='g')

cplt.add_colorbar(fig,pc,orientation='horizontal')
if printtofile:
    if addsig:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_AGCM_Esims' + prstr + '.png',dpi=400) 
    else:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_AGCM_Esims' + prstr + '.pdf') 


# compare Var SIC and cons SIC (E)
fig,axs=plt.subplots(1,3)
fig.set_size_inches(10,4)
fig.subplots_adjust(wspace=0.05,hspace=0.03)

ax=axs[0]
bm,pc=cplt.kemmap((alow2sp-ahigh2sp),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Var (Low-High)',suppcb=True,
                  panellab='a',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,a2sppval, alat, alon)
bm.contour(alons,alats,(alow2sp2-ahigh2sp2),levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,a2sp2pval,alat,alon,type='cont',color='g')

ax=axs[1]
bm,pc=cplt.kemmap((aelow2sp-aehigh2sp),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Const (Low-High)',suppcb=True,
                  panellab='b',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,ae2sppval, alat, alon)
bm.contour(alons,alats,(aelow2sp2-aehigh2sp2),levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,ae2sp2pval,alat,alon,type='cont',color='g')

ax=axs[2] 
bm,pc=cplt.kemmap((alow2sp-ahigh2sp)-(aelow2sp-aehigh2sp),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Var-Const (Low-High)',suppcb=True,
                  panellab='c',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,aressppval,alat,alon)
bm.contour(alons,alats,(alow2sp2-ahigh2sp2)-(aelow2sp2-aehigh2sp2),levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,aressp2pval,alat,alon,type='cont',color='g')

cplt.add_colorbar(fig,pc,orientation='horizontal')
if printtofile:
    if addsig:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_AGCMR-AGCME' + prstr + '.png',dpi=400) 
    else:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_AGCMR-AGCME' + prstr + '.pdf') 





# CGCM and AGCM--------------
#   composite on Eurasian SAT

fig,axs=plt.subplots(2,3)
fig.set_size_inches(10,8)
fig.subplots_adjust(wspace=0.05,hspace=0.03)

ax=axs[0,0]
bm,pc=cplt.kemmap(high2lesp-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='High Anom',suppcb=True,
                  panellab='a',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lehi2sppval,lat,lon)
bm.contour(lons,lats,high2lesp2-meanlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lehi2sp2pval,lat,lon,type='cont',color='g')

ax.set_ylabel('CGCM Diff')
ax=axs[0,1]
bm,pc=cplt.kemmap(low2lesp-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low Anom',suppcb=True,
                  panellab='b',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lelo2sppval,lat,lon)
bm.contour(lons,lats,low2lesp2-meanlesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lelo2sp2pval,lat,lon,type='cont',color='g')

ax=axs[0,2]
bm,pc=cplt.kemmap(low2lesp-high2lesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low-High',suppcb=True,
                  panellab='c',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,le2sppval, lat, lon)
bm.contour(lons,lats,low2lesp2-high2lesp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,le2sp2pval,lat,lon,type='cont',color='g')

ax=axs[1,0]
bm,pc=cplt.kemmap(ahigh2sp-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='High Anom',suppcb=True,
                  panellab='d',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,ahi2sppval,alat,alon)
bm.contour(alons,alats,ahigh2sp2-ameansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,ahi2sp2pval,alat,alon,type='cont',color='g')
ax.set_ylabel('AGCM Diff')

ax=axs[1,1]
bm,pc=cplt.kemmap(alow2sp-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low Anom',suppcb=True,
                  panellab='e',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,alo2sppval,alat,alon)
bm.contour(alons,alats,alow2sp2-ameansp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,alo2sp2pval,alat,alon,type='cont',color='g')

ax=axs[1,2]
bm,pc=cplt.kemmap(alow2sp-ahigh2sp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='Low-High',suppcb=True,
                  panellab='f',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,a2sppval, alat, alon)
bm.contour(alons,alats,alow2sp2-ahigh2sp2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,a2sp2pval,alat,alon,type='cont',color='g')

cplt.add_colorbar(fig,pc,orientation='horizontal')

if printtofile:
    if addsig:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_CGCMAGCM' + prstr + '.png',dpi=400) 
    else:
        fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_CGCMAGCM' + prstr + '.pdf') 
