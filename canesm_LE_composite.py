""" taken from canesm_LE_regress.py: July 2, 2015
"""


import loadLE as le
import cccmaplots as cplt
import cccmautils as cutl
import canesm_LE_general as leg
import loadmodeldata as lmd
import pandas as pd
import constants as con
import loadCanESM2data as lcd
import matplotlib.lines as mlines


printtofile=True

dataloaded=False
dofigures=False
local=True
addsig=True
compagcm=False # compare R sims and E sims

timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'
timeselall = '1979-01-01,2012-12-31'

ftype='fullts' # 'fullclimo' or 'climo' or 'fullts'
# this is not implemented here. but would do regressions in time on LE ensemble avg
ensmean=False 

seasp='DJF' # season of spatial field
sear='DJF' # season of regional avgs
diffttl1=diffttl2=diffttl3='Low-High' # for comp on BKS SIC and Eur SAT (otherwise High-Low)
diffmult1=diffmult2=diffmult3=1 # if High-Low then need to mult by -1

# spatial field1 in color
leconvsp=1
fieldsp='tas'; ncfieldsp='tas'; compsp='Amon'; 
cminsp=-3; cmaxsp=3 # for colors
cminspbig=-3; cmaxspbig=3 # for colors

cminspa=-1; cmaxspa=1 # for colors for AGCM
cminice=-10; cmaxice=10 # colors for ice

# spatial field2 in contours
leconvsp2=1
fieldsp2='zg50000.00'; ncfieldsp2='zg'; compsp2='Amon'
cminsp2=-30; cmaxsp2=30 # to calc contour interval
cminsp2big=-30; cmaxsp2big=30 # to calc contour interval

cminsp2a=-10; cmaxsp2a=10 # to calc contour interval for AGCM


# COMPOSITE ON THESE VALUES
leconvr=leconvr2=leconvr3=1

# regional avg field 1
fieldr='sic'; ncfieldr='sic'; compr='OImon'; regionr='bksmori'; r1str='BKS SIC'; #leconvr=-1 # sea ice loss is linked with positive changes elsewhere

# regional avg field 2
# cooling=high heights
fieldr2='tas'; ncfieldr2='tas'; compr2='Amon'; regionr2='eurasiamori'; r2str='Eur SAT'; #leconvr2=-1

# regional avg field 3
fieldr3='zg50000.00'; ncfieldr3='zg'; compr3='Amon'; regionr3='bksmori'; r3str='BKS Z500'; diffttl3='High-Low'; diffmult3=-1 
#leconvr=-1; #leconvr2=-1; #both conv -1 to get figs to show low-high equal to high heights and cold continent.

sttl1='Comp on ' + r1str
sttl2='Comp on ' + r2str
sttl3='Comp on ' + r3str

fdictsp = {'field': fieldsp, 'ncfield': ncfieldsp, 'comp': compsp}
fdictsp2 = {'field': fieldsp2, 'ncfield': ncfieldsp2, 'comp': compsp2}
fdictice = {'field': 'sic', 'ncfield': 'sic', 'comp':'OImon'}
fdictr = {'field': fieldr+regionr, 'ncfield': ncfieldr, 'comp': compr}
fdictr2 = {'field': fieldr2+regionr2, 'ncfield': ncfieldr2, 'comp': compr2}
fdictr3 = {'field': fieldr3+regionr3,'ncfield': ncfieldr3, 'comp': compr3}


lat=le.get_lat(local=local)
lon=le.get_lon(local=local)
nlat=len(lat); nlon=len(lon)

if addsig:
    prstr='sig'
else:
    prstr=''


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

def load_canesmfield(fdict,casename,seas,conv=1,subsampyrs=11,numsamp=50,
                     styear=None,anomyears=None,local=False,verb=False,addcyc=True):
    """ loads subsampled CGCM data from specified simulation 
             (assumes a long control, e.g. piControl)

             number of total chunks will be determined by length of sim data
                      and numyrs (ntime / numyrs). First the simulation is chunked into
                      numyrs segments. Then 2 segments at a time are randomly chosen
                      (at least a decade apart) to generate anomalies. This is done 
                      numsamp times (e.g. 50)

              returns subsamp  (subsampled anomalies)
                      styears  (start index of chunking of long run)
                      anomyears (indices of anomaly differences)
    """

    pidat = lcd.load_data(fdict,casename,local=local,conv=conv,verb=verb)
    piseadat = cutl.seasonalize_monthlyts(pidat,season=seas)
    # the data must be seasonalized before using this func.
    pisea,styear,anomyears = leg.subsamp_anom_pi(piseadat, numyrs=subsampyrs,numsamp=numsamp,
                                                 styear=styear,anomyears=anomyears)
    if addcyc:
        st=pisea[...,-1]
        pisea=np.dstack((pisea,st[...,None]))
    return pisea,styear,anomyears


def load_agcmfield(field,sims,seas,conv=1,region=None,subsampyrs=11,styears=None):
    """ loads subsampled agcm data from specified simulations

              number of total samples will be determined by length of all sim data
                      and numyrs (e.g. (ntime / numyrs)*numsims)

              returns subsamp
                      styears
    """
    threed=False

    simconv=conv
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


def composite_lefield(casenames, fdictsp, fdictr, loaddictsp, loaddictr,verb=False,local=False):
    """ composite fdictsp on fdictr 

           returns (hi spatial comp, lo spatial comp, mean spatial, hipval, lopval, lo v hi pval, hi idx, lo idx)
    """

    for cii,casename in enumerate(casenames):
        # loaddict should have: timeselp, timeselc, sea, ftype, leconv
        #     for each field being loaded
        # SPATIAL FIELD
        lecseasp = load_field(fdictsp, casename, loaddictsp['timeselc'], loaddictsp['sea'],
                              ftype=loaddictsp['ftype'], conv=loaddictsp['conv'],verb=verb,local=local)
        lepseasp = load_field(fdictsp, casename, loaddictsp['timeselp'], loaddictsp['sea'],
                              ftype=loaddictsp['ftype'], conv=loaddictsp['conv'],verb=verb,local=local)
        leseasp = lepseasp-lecseasp

        # REGIONAL AVG FIELD: 1D
        lecsear = load_field(fdictr, casename, loaddictr['timeselc'], loaddictr['sea'],
                              ftype=loaddictr['ftype'], conv=loaddictr['conv'],verb=verb,local=local)
        lepsear = load_field(fdictr, casename, loaddictr['timeselp'], loaddictr['sea'],
                              ftype=loaddictr['ftype'], conv=loaddictr['conv'],verb=verb,local=local)
        lesear = lepsear-lecsear

        if cii != 0: 
            tmp = np.vstack((tmp,leseasp))
            tmpr = np.hstack((tmpr,lesear))
        else:
            tmp = leseasp
            tmpr = lesear

    leseasp=tmp
    lesear=tmpr

    return do_composite(lesear,leseasp,rshape=(nlat,nlon),verb=verb,addcyc=True)

# end composite_lefield()

def do_composite(rfld,spfld,nn=10,rshape=None,verb=False,addcyc=False):
    """ 

          returns (hi spatial comp, lo spatial comp, mean spatial, hipval, lopval, lo v hi pval, hi idx, lo idx)
    """
   
    numens = rfld.shape[0] # numens or numsamp
    highidx = (-rfld).argsort()[:nn]
    lowidx = rfld.argsort()[:nn]

    if rshape==None:
        # this should have the effect of keeping shape the same
        # when calling reshape (since it's not necessary)
        rshape=spfld.shape[1:] # lat x lon dim

    rshapenn=(nn,)+rshape
    rshapeens=(numens,)+rshape

    meanlesp = spfld.mean(axis=0).reshape(rshape)
    highlesp = spfld[highidx,...].mean(axis=0).reshape(rshape)
    lowlesp = spfld[lowidx,...].mean(axis=0).reshape(rshape)

    if verb:
        #print 'highidx: ' + str(highidx)
        print 'highidx vals: ' + str(rfld[highidx])
        #print 'lowidx vals: ' + str(lowidx)
        print 'lowidx vals: ' + str(rfld[lowidx])
    
    # low vs high


    (lesptstat, lesppval) = cutl.ttest_ind(spfld[lowidx,...].reshape(rshapenn), 
                                           spfld[highidx,...].reshape(rshapenn))
    # low vs mean
    (lelosptstat, lelosppval) = cutl.ttest_ind(spfld[lowidx,...].reshape(rshapenn), 
                                               spfld.reshape(rshapeens))
    # high vs mean
    (lehisptstat, lehisppval) = cutl.ttest_ind(spfld[highidx,...].reshape(rshapenn), 
                                               spfld.reshape(rshapeens))

    if addcyc:
        st=meanlesp[...,-1]
        meanlesp=np.hstack((meanlesp,st[...,None]))
        st=highlesp[...,-1]
        highlesp=np.hstack((highlesp,st[...,None]))
        st=lowlesp[...,-1]
        lowlesp=np.hstack((lowlesp,st[...,None]))
        st=lesppval[...,-1]
        lesppval=np.hstack((lesppval,st[...,None]))
        st=lehisppval[...,-1]
        lehisppval=np.hstack((lehisppval,st[...,None]))
        st=lelosppval[...,-1]
        lelosppval=np.hstack((lelosppval,st[...,None]))

    compdt = {'highsp': highlesp, 'lowsp': lowlesp, 'meansp': meanlesp,
              'hisppval': lehisppval, 'losppval': lelosppval, 'sppval': lesppval,
              'highidx': highidx, 'lowidx': lowidx}

    return compdt

# end do_composite()



if not dataloaded:
    nn=10 #@@
    #casenames=('historical','historicalNat','historicalMisc')
    casenames=('historical',)
    loaddictsp={'timeselc':timeselc, 'timeselp': timeselp, 
                'sea': seasp, 'ftype': ftype, 'conv':leconvsp}
    loaddictsp2={'timeselc':timeselc, 'timeselp': timeselp, 
                 'sea': seasp, 'ftype': ftype, 'conv':leconvsp2}
    loaddictice={'timeselc':timeselc, 'timeselp': timeselp, 
                'sea': seasp, 'ftype': ftype, 'conv':leconvsp}
    loaddictr={'timeselc':timeselc, 'timeselp': timeselp, 
               'sea': sear, 'ftype': ftype, 'conv':leconvr}
    loaddictr2={'timeselc':timeselc, 'timeselp': timeselp, 
                'sea': sear, 'ftype': ftype, 'conv':leconvr2}
    loaddictr3={'timeselc':timeselc, 'timeselp': timeselp, 
                'sea': sear, 'ftype': ftype, 'conv':leconvr3}

    print fdictsp
    print fdictr
    # "high" for regional BKS SIC is 'high ice loss' (conv=-1)
    # "high" for regional Eur SAT is 'warm temp anomaly' (conv=1)
    # if conv=1, "high" for regional BKS Z500 is 'high height anomaly'

    # comp on region 1
    print 'CGCM comp on BKS SIC'
    # highlesp, lowlesp, meanlesp, lehisppval, lelosppval, lesppval, highidx, lowidx
    lespr1dt = composite_lefield(casenames, 
                                 fdictsp, fdictr, 
                                 loaddictsp, loaddictr,
                                 verb=True,local=local)
    # highlesp2, lowlesp2, meanlesp2, lehisp2pval, lelosp2pval, lesp2pval, highidx, lowidx
    lesp2r1dt = composite_lefield(casenames, 
                                  fdictsp2, fdictr, 
                                  loaddictsp2, loaddictr,
                                  verb=True,local=local)
    # highleice, lowleice, meanleice, lehiicepval, leloicepval, leicepval, highidx, lowidx
    leicer1dt = composite_lefield(casenames, 
                                  fdictice, fdictr, 
                                  loaddictice, loaddictr,
                                  verb=True,local=local)
    # comp on region 2
    print 'CGCM comp on Eur SAT'
    # highlespr2, lowlespr2, meanlespr2, lehispr2pval, lelospr2pval, lespr2pval, highidx2, lowidx2
    lespr2dt = composite_lefield(casenames, 
                                 fdictsp, fdictr2, 
                                 loaddictsp, loaddictr2,
                                 verb=True,local=local)
    # highlesp2r2, lowlesp2r2, meanlesp2r2, lehisp2r2pval, lelosp2r2pval, lesp2r2pval, highidx2, lowidx2
    lesp2r2dt = composite_lefield(casenames, 
                                  fdictsp2, fdictr2, 
                                  loaddictsp2, loaddictr2,
                                  verb=True,local=local)
    # highleicer2, lowleicer2, meanleicer2, lehiicer2pval, leloicer2pval, leicer2pval, highidx2, lowidx2
    leicer2dt = composite_lefield(casenames, 
                                  fdictice, fdictr2, 
                                  loaddictice, loaddictr2,
                                  verb=True,local=local)

    # comp on region 3
    print 'CGCM comp on BKS Z500'
    lespr3dt = composite_lefield(casenames, 
                                 fdictsp, fdictr3, 
                                 loaddictsp, loaddictr3,
                                 verb=True,local=local)
    lesp2r3dt = composite_lefield(casenames, 
                                  fdictsp2, fdictr3, 
                                  loaddictsp2, loaddictr3,
                                  verb=True,local=local)
    leicer3dt = composite_lefield(casenames, 
                                  fdictice, fdictr3, 
                                  loaddictice, loaddictr3,
                                  verb=True,local=local)


    # ========== PRE-IND ==============

    piseasp,styear,anomyears = load_canesmfield(fdictsp,'piControl',seasp,conv=leconvsp,local=local)
    piseasp2,styear,anomyears = load_canesmfield(fdictsp2,'piControl',seasp,conv=leconvsp2,
                                                 local=local,styear=styear,anomyears=anomyears)
    piseaspice,styear,anomyears = load_canesmfield(fdictice,'piControl',sear,conv=1,
                                                 local=local,styear=styear,anomyears=anomyears)

    pisear,styear,anomyears = load_canesmfield(fdictr,'piControl',sear,conv=leconvr,
                                               local=local,styear=styear,anomyears=anomyears,addcyc=False)
    pisear2,styear,anomyears = load_canesmfield(fdictr2,'piControl',sear,conv=leconvr2,
                                                local=local,styear=styear,anomyears=anomyears,addcyc=False)
    pisear3,styear,anomyears = load_canesmfield(fdictr3,'piControl',sear,conv=leconvr3,
                                                local=local,styear=styear,anomyears=anomyears,addcyc=False)

    print 'PI comp on BKS SIC'
    # pihighsp,pilowsp,pimeansp,pihisppval,pilosppval,pisppval,pihighidx,pilowidx
    pispr1dt = do_composite(pisear,piseasp,verb=True)
    # pihighsp2,pilowsp2,pimeansp2,pihisp2pval,pilosp2pval,pisp2pval,pihighidx,pilowidx
    pisp2r1dt = do_composite(pisear,piseasp2,verb=True)
    # pihighice,pilowice,pimeanice,pihiicepval,piloicepval,piicepval,pihighidx,pilowidx
    piicer1dt = do_composite(pisear,piseaspice,verb=True)

    print 'PI comp on Eur SAT'
    # pihighspr2,pilowspr2,pimeanspr2,pihispr2pval,pilospr2pval,pispr2pval,pihighidxr2,pilowidxr2
    pispr2dt = do_composite(pisear2,piseasp,verb=True)
    # pihighsp2r2,pilowsp2r2,pimeansp2r2,pihisp2r2pval,pilosp2r2pval,pisp2r2pval,pihighidxr2,pilowidxr2
    pisp2r2dt = do_composite(pisear2,piseasp2,verb=True)
    # pihighicer2,pilowicer2,pimeanicer2,pihiicer2pval,piloicer2pval,piicer2pval,pihighidxr2,pilowidxr2
    piicer2dt = do_composite(pisear2,piseaspice,verb=True)

    print 'PI comp on BKS Z500'
    pispr3dt = do_composite(pisear3,piseasp,verb=True)
    pisp2r3dt = do_composite(pisear3,piseasp2,verb=True)
    piicer3dt = do_composite(pisear3,piseaspice,verb=True)

    # === AGCM ==========
    alat=con.get_t63lat(); alon=con.get_t63lon()
    alons, alats = np.meshgrid(alon,alat)

    simsE=('E1','E2','E3','E4','E5'); 
    sims=('R1','R2','R3','R4','R5');

    asseasp,styears = load_agcmfield(fieldsp,sims,seasp) # already anomalies
    asseasp2,styears = load_agcmfield(fieldsp2,sims,seasp,styears=styears) # already anomalies
    asseaice,styears = load_agcmfield('sic',sims,seasp,styears=styears,conv=100) # want the SICN pattern

    assear,styears = load_agcmfield(fieldr,sims,sear,styears=styears,region=regionr,conv=leconvr) # already anomalies
    assear2,styears = load_agcmfield(fieldr2,sims,sear,styears=styears,region=regionr2,conv=leconvr2) # already anomalies
    assear3,styears = load_agcmfield(fieldr3,sims,sear,styears=styears,region=regionr3,conv=leconvr3) # already anomalies
    assear=assear # / assear.std()
    assear2=assear2 #/ assear2.std()
    assear3=assear3

    (anens,anlat,anlon)=asseasp.shape
    rshape=(anens,anlat*anlon)


    # =========== composite on BKS SIC (AGCM) ===============
    # calc composites: choose top and bottom BKS SIC (10 subsamps each)

    print 'AGCM comp on BKS SIC'
    # ahighsp,alowsp,ameansp,ahisppval,alosppval,asppval,ahighidx,alowidx
    aspr1dt= do_composite(assear,asseasp,verb=True)
    # ahighsp2,alowsp2,ameansp2,ahisp2pval,alosp2pval,asp2pval,ahighidx,alowidx
    asp2r1dt = do_composite(assear,asseasp2,verb=True)
    # ahighice,alowice,ameanice,ahiicepval,aloicepval,aicepval,ahighidx,alowidx
    aicer1dt = do_composite(assear,asseaice,verb=True)

    print 'AGCM comp on Eur SAT'
    # ahighspr2,alowspr2,ameanspr2,ahispr2pval,alospr2pval,aspr2pval,ahighidxr2,alowidxr2
    aspr2dt = do_composite(assear2,asseasp,verb=True)
    # ahighsp2r2,alowsp2r2,ameansp2r2,ahisp2r2pval,alosp2r2pval,asp2r2pval,ahighidxr2,alowidxr2
    asp2r2dt = do_composite(assear2,asseasp2,verb=True)
    # ahighicer2,alowicer2,ameanicer2,ahiicer2pval,aloicer2pval,aicer2pval,ahighidxr2,alowidxr2
    aicer2dt = do_composite(assear2,asseaice,verb=True)

    print 'AGCM comp on BKS Z500'
    aspr3dt = do_composite(assear3,asseasp,verb=True)
    asp2r3dt = do_composite(assear3,asseasp2,verb=True)
    aicer3dt = do_composite(assear3,asseaice,verb=True)

    if compagcm: # if compare AGCM ensembles

        # ===  E sims:
        aesseasp,styearse = load_agcmfield(fieldsp,simsE,seasp) 
        aesseasp2,styearse = load_agcmfield(fieldsp2,simsE,seasp,styears=styearse) 
        aessear2,styearse = load_agcmfield(fieldr2,simsE,sear,styears=styearse,region=regionr2) 
        aessear2=aessear2 #/ aessear2.std()

        # E sims ----
        aehighspr2,aelowspr2,aemeanspr2,aehispr2pval,aelospr2pval,aespr2pval,aehighidxr2,aelowidxr2= do_composite(aessear2,aesseasp,verb=True)
        aehighsp2r2,aelowsp2r2,aemeansp2r2,aehisp2r2pval,aelosp2r2pval,aesp2r2pval,aehighidxr2,aelowidxr2= do_composite(aessear2,aesseasp2,verb=True)

        # AGCM RESIDUAL pvals (var vs const)
        (aressptstat,aressppval) = cutl.ttest_ind(asseasp[alow2idx,...]- asseasp[ahigh2idx,...],
                                                  aesseasp[aelow2idx,...]- aesseasp[aehigh2idx,...])
        (aressp2tstat,aressp2pval) = cutl.ttest_ind(asseasp2[alow2idx,...]- asseasp2[ahigh2idx,...],
                                                  aesseasp2[aelow2idx,...]- aesseasp2[aehigh2idx,...])





# For each composite (r1, r2, r3), compute the BKS SIC average:
icedt = {'bkssic': {'PI':piicer1dt,'CGCM':leicer1dt,'AGCM':aicer1dt},
         'eursat': {'PI':piicer2dt,'CGCM':leicer2dt,'AGCM':aicer2dt},
         'bksz500': {'PI':piicer3dt,'CGCM':leicer3dt,'AGCM':aicer3dt} }

# For each composite (r1,r2,r3), compute Eur SAT average:
spdt = {'bkssic': {'PI':pispr1dt,'CGCM':lespr1dt,'AGCM':aspr1dt},
        'eursat': {'PI':pispr2dt,'CGCM':lespr2dt,'AGCM':aspr2dt},
        'bksz500': {'PI':pispr3dt,'CGCM':lespr3dt,'AGCM':aspr3dt} }

# For each composite (r1,r2,r3), compute BKS Z500 average: (prob don't need)
sp2dt = {'bkssic': {'PI':pisp2r1dt,'CGCM':lesp2r1dt,'AGCM':asp2r1dt},
         'eursat': {'PI':pisp2r2dt,'CGCM':lesp2r2dt,'AGCM':asp2r2dt},
         'bksz500': {'PI':pisp2r3dt,'CGCM':lesp2r3dt,'AGCM':asp2r3dt} }



#  Here calculate the BKS SIC associated with each composite
allregimdt={}; allregitdt={}
for iii,ikey in enumerate(icedt.keys()):
    regmdt={}
    regtotdt={}
    dts=icedt[ikey]
    for dkey in dts.keys():
        dt=dts[dkey]
        diff = dt['lowsp']-dt['highsp']
        regtot = cutl.calc_regtotseaicearea(diff[:,:-1],lat,lon,'bksmori')
        regm = cutl.calc_regmean(diff[:,:-1],lat,lon,'bksmori')

        regmdt[dkey]=regm
        regtotdt[dkey]=regtot # not sure which one i want

    allregimdt[ikey]=regmdt # i for ice
    allregitdt[ikey]=regtotdt

allregimdf=pd.DataFrame(allregimdt)
allregitdf=pd.DataFrame(allregitdt)

#  Here calculate the Eur SAT associated with each composite
allregspmdt={}
for iii,ikey in enumerate(spdt.keys()):# for each composite region
    regmdt={}
    regtotdt={}
    dts=spdt[ikey] 
    for dkey in dts.keys():# for each ensemble
        dt=dts[dkey]
        diff = dt['lowsp']-dt['highsp']
        regm = cutl.calc_regmean(diff[:,:-1],lat,lon,'eurasiamori')

        regmdt[dkey]=regm

    allregspmdt[ikey]=regmdt # sp for spatial 1 (SAT)

allregspmdf=pd.DataFrame(allregspmdt)


lons, lats = np.meshgrid(lon,lat)
cmlen=15.
incr = (cmaxsp2-cminsp2) / (cmlen)
conts = np.arange(cminsp2,cmaxsp2+incr,incr)

incrbig = (cmaxsp2big-cminsp2big) / (cmlen)
contsbig = np.arange(cminsp2big,cmaxsp2big+incrbig,incrbig)

stxx=1
xx=np.arange(stxx,stxx+len(allregimdf.keys()))
regs=('bkssic','eursat','bksz500')
multfacs=(1,1,-1) # multiply the comp on z500 by -1 to get high over bks
enss=('PI','CGCM','AGCM')
mkrs=('s','o','d')
clrs=('r','b','0.7')

# prepare legend entries
lgs=(mlines.Line2D([],[],color=clrs[0],linestyle='none',marker=mkrs[0]),
     mlines.Line2D([],[],color=clrs[1],linestyle='none',marker=mkrs[1]),
     mlines.Line2D([],[],color=clrs[2],linestyle='none',marker=mkrs[2])) 
lgstrs=('BKS SIC','Eur SAT','BKS Z500')          

#fig,axs=plt.subplots(1,2)
#fig.set_size_inches((10,6))
fig,axs=plt.subplots(2,2) # maps on top
fig.set_size_inches((10,10))
ax=axs[0,0] # CGCM comp on BKS SIC
bm,pc=cplt.kemmap(lespr1dt['lowsp']-lespr1dt['highsp'],alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title=sttl1,suppcb=False,
                  panellab='a',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lespr1dt['sppval'], alat, alon)
bm.contour(alons,alats,lesp2r1dt['lowsp']-lesp2r1dt['highsp'],levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2r1dt['sppval'],alat,alon,type='cont',color='g')

ax=axs[0,1] # CGCM comp on Eur SAT
bm,pc=cplt.kemmap(lespr2dt['lowsp']-lespr2dt['highsp'],alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=sttl2,suppcb=False,
                      panellab='b',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lespr2dt['sppval'], alat, alon)
bm.contour(alons,alats,lesp2r2dt['lowsp']-lesp2r2dt['highsp'],levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2r2dt['sppval'],alat,alon,type='cont',color='g')

ax=axs[1,0]
for rii,reg in enumerate(regs):
    for eii,ens in enumerate(enss):
        mult=multfacs[rii]
        ax.plot(xx[eii],mult*allregimdf[reg][ens],marker=mkrs[rii],
                color=clrs[rii],linestyle='none',markersize=10)
ax.set_ylabel('$\Delta$ ' + r1str)
ax.axhline(y=0,color='k',linestyle='--')
ax.set_xlim((stxx-0.5,len(enss)+1.5))
ax.legend(lgs,lgstrs,loc='lower right')
ax.set_xticks(xx)
ax.set_xticklabels(enss)
ax.annotate('c',xy=(0.01,1.01),
            xycoords='axes fraction',fontsize=16,fontweight='bold')

ax=axs[1,1]
for rii,reg in enumerate(regs):
    for eii,ens in enumerate(enss):
        mult=multfacs[rii]
        ax.plot(xx[eii],mult*allregspmdf[reg][ens],marker=mkrs[rii],
                color=clrs[rii],linestyle='none',markersize=10)
ax.set_ylabel('$\Delta$ ' + r2str)
ax.axhline(y=0,color='k',linestyle='--')
ax.set_xlim((stxx-0.5,len(enss)+0.5))
ax.set_xticks(xx)
ax.set_xticklabels(enss)
ax.annotate('d',xy=(0.01,1.01),
            xycoords='axes fraction',fontsize=16,fontweight='bold')
if addsig:
    prstr='sig'
else:
    prstr=''
if printtofile:
    fig.savefig('Figure4_draft_4panel_CGCMmaps_compavgs' + prstr + '.pdf')
    fig.savefig('Figure4_draft_4panel_CGCMmaps_compavgs' + prstr + '.png',dpi=400)


# TEST FIGURE: All 3 comps
fig,axs=plt.subplots(1,3)
fig.set_size_inches(10,4)
ax=axs[0] # CGCM comp on BKS SIC
bm,pc=cplt.kemmap(lespr1dt['lowsp']-lespr1dt['highsp'],alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title=sttl1,suppcb=True,
                  panellab='a',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lespr1dt['sppval'], alat, alon)
bm.contour(alons,alats,lesp2r1dt['lowsp']-lesp2r1dt['highsp'],levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2r1dt['sppval'],alat,alon,type='cont',color='g')

ax=axs[1] # CGCM comp on Eur SAT
bm,pc=cplt.kemmap(lespr2dt['lowsp']-lespr2dt['highsp'],alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=sttl2,suppcb=True,
                      panellab='b',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lespr2dt['sppval'], alat, alon)
bm.contour(alons,alats,lesp2r2dt['lowsp']-lesp2r2dt['highsp'],levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2r2dt['sppval'],alat,alon,type='cont',color='g')

ax=axs[2]
bm,pc=cplt.kemmap(-1*(lespr3dt['lowsp']-lespr3dt['highsp']),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=sttl3,suppcb=True,
                      panellab='c',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lespr3dt['sppval'], alat, alon)
bm.contour(alons,alats,-1*(lesp2r3dt['lowsp']-lesp2r3dt['highsp']),levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2r3dt['sppval'],alat,alon,type='cont',color='g')
if printtofile:
    fig.savefig('CGCM_all3comps_maps' + prstr + '.pdf')
    fig.savefig('CGCM_all3comps_maps' + prstr + '.png',dpi=400)



# for plotting the ice maps below
rplotdts = (piicer1dt,leicer1dt,aicer1dt)
diffmult=diffmult3; diffttl=diffttl3
sttl=sttl1
panstr='a'; panttl='BKS SIC'


#rplotdts = (piicer2dt,leicer2dt,aicer2dt)
#diffmult=diffmult2; diffttl=diffttl2
#sttl=sttl2
#panstr='b'; panttl='Eur SAT'

# plot ice maps:
ylabs=('PI','CGCM','AGCM')
allplabs = (('a','b','c'),
            ('d','e','f'),
            ('g','h','i'))

cmapice='red2blue_w20'
fig,axs=plt.subplots(3,3)
fig.set_size_inches(10,10)
fig.subplots_adjust(wspace=0.05,hspace=0.03)

for aii in np.arange(0,3): 

    axidx=aii # loop rows
    dt=rplotdts[axidx]
    plabs=allplabs[axidx]

    ax = axs[axidx,0]
    bm,pc=cplt.kemmap(dt['highsp']-dt['meansp'],alat,alon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab=plabs[0],lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,dt['hisppval'],alat,alon)
    ax.set_ylabel(ylabs[axidx])

    ax = axs[axidx,1]
    bm,pc=cplt.kemmap(dt['lowsp']-dt['meansp'],alat,alon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab=plabs[1],lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,dt['losppval'],alat,alon)

    ax = axs[axidx,2]
    bm,pc=cplt.kemmap(diffmult*(dt['lowsp']-dt['highsp']),alat,alon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab=plabs[2],lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,dt['sppval'],alat,alon)

    plt.suptitle(sttl)








# ==============================================
# ====================== FIGURES ===============

if dofigures:
    # FIRST PLOT THE HIGH AND LOW SEA ICE ANOMS
    #printtofile=True
    cmapice='red2blue_w20'
    fig,axs=plt.subplots(3,3)
    fig.set_size_inches(10,10)
    fig.subplots_adjust(wspace=0.05,hspace=0.03)

    ax=axs[0,0]
    bm,pc=cplt.kemmap(pihighice-pimeanice,lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='a',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pihiicepval,lat,lon)
    ax.set_ylabel('PI')
    ax=axs[0,1]
    bm,pc=cplt.kemmap(pilowice-pimeanice,lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='b',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,piloicepval,lat,lon)
    ax=axs[0,2]
    bm,pc=cplt.kemmap(diffmult*(pilowice-pihighice),lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='c',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,piicepval,lat,lon)
    ax=axs[1,0]
    bm,pc=cplt.kemmap(highleice-meanleice,lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehiicepval,lat,lon)
    ax.set_ylabel('CGCM')
    ax=axs[1,1]
    bm,pc=cplt.kemmap(lowleice-meanleice,lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,leloicepval,lat,lon)
    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(lowleice-highleice),lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm, leicepval,lat,lon)
    ax=axs[2,0]
    bm,pc=cplt.kemmap(ahighice-ameanice,alat,alon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='g',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahiicepval,alat,alon)
    ax.set_ylabel('AGCM')
    ax=axs[2,1]
    bm,pc=cplt.kemmap(alowice-ameanice,alat,alon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='h',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aloicepval,alat,alon)
    ax=axs[2,2]
    bm,pc=cplt.kemmap(diffmult*(alowice-ahighice),alat,alon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='i',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aicepval,alat,alon)
    cplt.add_colorbar(fig,pc,orientation='horizontal')

    plt.suptitle(sttl1)

    if printtofile:
        fig.savefig('sic_' + sear + \
                    '_compon_' + fieldr+regionr+ sear + '_PICGCMAGCM' + prstr + '.pdf') 


    fig,axs=plt.subplots(3,3)
    fig.set_size_inches(10,10)
    fig.subplots_adjust(wspace=0.05,hspace=0.03)

    ax=axs[0,0]
    bm,pc=cplt.kemmap(pihighicer2-pimeanice,lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='a',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pihiicer2pval,lat,lon)
    ax.set_ylabel('PI')
    ax=axs[0,1]
    bm,pc=cplt.kemmap(pilowicer2-pimeanice,lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='b',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,piloicer2pval,lat,lon)
    ax=axs[0,2]
    bm,pc=cplt.kemmap(diffmult*(pilowicer2-pihighicer2),lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='c',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,piicer2pval,lat,lon)
    ax=axs[1,0]
    bm,pc=cplt.kemmap(highleicer2-meanleice,lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehiicer2pval,lat,lon)
    ax.set_ylabel('CGCM')
    ax=axs[1,1]
    bm,pc=cplt.kemmap(lowleicer2-meanleice,lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,leloicer2pval,lat,lon)
    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(lowleicer2-highleicer2),lat,lon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,leicer2pval,lat,lon)
    ax=axs[2,0]
    bm,pc=cplt.kemmap(ahighicer2-ameanice,alat,alon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='g',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahiicer2pval,alat,alon)
    ax.set_ylabel('AGCM')
    ax=axs[2,1]
    bm,pc=cplt.kemmap(alowicer2-ameanice,alat,alon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='h',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aloicer2pval,alat,alon)
    ax=axs[2,2]
    bm,pc=cplt.kemmap(diffmult*(alowicer2-ahighicer2),alat,alon,type='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='i',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aicer2pval,alat,alon)
    cplt.add_colorbar(fig,pc,orientation='horizontal')

    plt.suptitle(sttl2)

    if printtofile:
        fig.savefig('sic_' + sear + \
                    '_compon_' + fieldr2+regionr2+ sear + '_PICGCMAGCM' + prstr + '.pdf') 







    ttl1=seasp + ': High ' + sear + r1str
    ttl2=seasp + ': Low ' + sear + r1str

    #ttl1=ttl2=''



    # ========== 4 panel fig =========

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
    bm,pc=cplt.kemmap(diffmult*(lowlesp-highlesp),lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lesppval, lat, lon)
    bm.contour(lons,lats,diffmult*(lowlesp2-highlesp2),levels=conts,
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

    # =========== PI  ===========
    #printtofile=True
    fig,axs=plt.subplots(2,3)
    fig.set_size_inches(10,8)
    fig.subplots_adjust(wspace=0.05,hspace=0.03)

    ax=axs[0,0]
    bm,pc=cplt.kemmap(pihighsp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl1,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(lons,lats,pihighsp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('PreIndustrial')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(pilowsp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl2,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(lons,lats,pilowsp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(pimeansp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(lons,lats,pimeansp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    #cplt.add_colorbar(fig,pc,orientation='horizontal',pos=[.25,.5, .5,.03])

    ax=axs[1,0]
    bm,pc=cplt.kemmap(pihighsp-pimeansp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pihisppval,lat,lon)
    bm.contour(lons,lats,pihighsp2-pimeansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pihisp2pval,lat,lon,type='cont',color='g')
    ax.set_ylabel('PreI Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(pilowsp-pimeansp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pilosppval,lat,lon)
    bm.contour(lons,lats,pilowsp2-pimeansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pilosp2pval,lat,lon,type='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(pilowsp-pihighsp),lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pisppval, lat, lon)
    bm.contour(lons,lats,diffmult*(pilowsp2-pihighsp2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pisp2pval,lat,lon,type='cont',color='g')

    cplt.add_colorbar(fig,pc,orientation='horizontal')
    if printtofile:
        if addsig:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr+regionr+ sear + '_PI' + prstr + '.png', dpi=400)
        else:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr+regionr+ sear + '_PI' + prstr + '.pdf') 

    #printtofile=False

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
    bm,pc=cplt.kemmap(diffmult*(alowsp-ahighsp),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,asppval, alat, alon)
    bm.contour(alons,alats,diffmult*(alowsp2-ahighsp2),levels=conts,
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
    bm,pc=cplt.kemmap(diffmult*(lowlesp-highlesp),lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
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
    bm,pc=cplt.kemmap(diffmult*(alowsp-ahighsp),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,asppval, alat, alon)
    bm.contour(alons,alats,diffmult*(alowsp2-ahighsp2),levels=conts,
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

    ttl12=seasp + ': High ' + sear + ' ' + r2str 
    ttl22=seasp + ': Low ' + sear + ' ' + r2str 


    fig,axs=plt.subplots(2,3)
    fig.set_size_inches(10,8)
    fig.subplots_adjust(wspace=0.05,hspace=0.03)

    ax=axs[0,0]
    bm,pc=cplt.kemmap(highlespr2,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl12,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(lons,lats,highlesp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('CGCM')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(lowlespr2,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl22,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(lons,lats,lowlesp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(lons,lats,meanlesp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[1,0]
    bm,pc=cplt.kemmap(highlespr2-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehispr2pval,lat,lon)
    bm.contour(lons,lats,highlesp2r2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lehisp2r2pval,lat,lon,type='cont',color='g')
    ax.set_ylabel('CGCM Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(lowlespr2-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lelospr2pval,lat,lon)
    bm.contour(lons,lats,lowlesp2r2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lelosp2r2pval,lat,lon,type='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(lowlespr2-highlespr2),lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lespr2pval, lat, lon)
    bm.contour(lons,lats,diffmult*(lowlesp2r2-highlesp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lesp2r2pval,lat,lon,type='cont',color='g')

    cplt.add_colorbar(fig,pc,orientation='horizontal')
    if printtofile:
        if addsig:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr2+regionr2+ sear + '_CGCM' + prstr + '.png',dpi=400) 
        else:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr2+regionr2+ sear + '_CGCM' + prstr + '.pdf') 

    # =========== PI  ===========

    fig,axs=plt.subplots(2,3)
    fig.set_size_inches(10,8)
    fig.subplots_adjust(wspace=0.05,hspace=0.03)

    ax=axs[0,0]
    bm,pc=cplt.kemmap(pihighspr2,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl12,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(lons,lats,pihighsp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('PreIndustrial')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(pilowspr2,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl22,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(lons,lats,pilowsp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(pimeansp,lat,lon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(lons,lats,pimeansp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    #cplt.add_colorbar(fig,pc,orientation='horizontal',pos=[.25,.5, .5,.03])

    ax=axs[1,0]
    bm,pc=cplt.kemmap(pihighspr2-pimeansp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pihispr2pval,lat,lon)
    bm.contour(lons,lats,pihighsp2r2-pimeansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pihisp2r2pval,lat,lon,type='cont',color='g')
    ax.set_ylabel('PreI Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(pilowspr2-pimeansp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pilospr2pval,lat,lon)
    bm.contour(lons,lats,pilowsp2r2-pimeansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pilosp2r2pval,lat,lon,type='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(pilowspr2-pihighspr2),lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pispr2pval, lat, lon)
    bm.contour(lons,lats,diffmult*(pilowsp2r2-pihighsp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pisp2r2pval,lat,lon,type='cont',color='g')

    cplt.add_colorbar(fig,pc,orientation='horizontal')
    if printtofile:
        if addsig:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                    '_compon_' + fieldr2+regionr2+ sear + '_PI' + prstr + '.png', dpi=400)
        else:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr2+regionr2+ sear + '_PI' + prstr + '.pdf') 




    # ======= AGCM
    fig,axs=plt.subplots(2,3)
    fig.set_size_inches(10,8)
    fig.subplots_adjust(wspace=0.05,hspace=0.03)

    ax=axs[0,0]
    bm,pc=cplt.kemmap(ahighspr2,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl12,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(alons,alats,ahighsp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('AGCM')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(alowspr2,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl22,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(alons,alats,alowsp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(ameansp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(alons,alats,ameansp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[1,0]
    bm,pc=cplt.kemmap(ahighspr2-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahispr2pval,alat,alon)
    bm.contour(alons,alats,ahighsp2r2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,ahisp2r2pval,alat,alon,type='cont',color='g')
    ax.set_ylabel('AGCM Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(alowspr2-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,alospr2pval,alat,alon)
    bm.contour(alons,alats,alowsp2r2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,alosp2r2pval,alat,alon,type='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(alowspr2-ahighspr2),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aspr2pval, alat, alon)
    bm.contour(alons,alats,diffmult*(alowsp2r2-ahighsp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,asp2r2pval,alat,alon,type='cont',color='g')

    cplt.add_colorbar(fig,pc,orientation='horizontal')
    if printtofile:
        if addsig:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr2+regionr2+ sear + '_AGCM' + prstr + '.png',dpi=400) 
        else:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr2+regionr2+ sear + '_AGCM' + prstr + '.pdf') 


    # ==== AGCM E SIMS !
    if compagcm:
        fig,axs=plt.subplots(2,3)
        fig.set_size_inches(10,8)
        fig.subplots_adjust(wspace=0.05,hspace=0.03)

        ax=axs[0,0]
        bm,pc=cplt.kemmap(aehighspr2,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                          title=ttl12,suppcb=True,
                          panellab='a',lcol='0.2')
        bm.contour(alons,alats,aehighsp2r2,levels=contsbig,
                   colors='0.5',linewidths=1,latlon=True)
        ax.set_ylabel('AGCM')

        ax=axs[0,1]
        bm,pc=cplt.kemmap(aelowspr2,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                          title=ttl22,suppcb=True,
                          panellab='b',lcol='0.2')
        bm.contour(alons,alats,aelowsp2r2,levels=contsbig,
                   colors='0.5',linewidths=1,latlon=True)

        ax=axs[0,2]
        bm,pc=cplt.kemmap(aemeansp,alat,alon,type='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                          title='Mean Anom',suppcb=True,
                          panellab='c',lcol='0.2')
        bm.contour(alons,alats,aemeansp2,levels=contsbig,
                   colors='0.5',linewidths=1,latlon=True)

        ax=axs[1,0]
        bm,pc=cplt.kemmap(aehighspr2-aemeansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='High Anom',suppcb=True,
                          panellab='d',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aehispr2pval,alat,alon)
        bm.contour(alons,alats,aehighsp2r2-aemeansp2,levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,aehisp2r2pval,alat,alon,type='cont',color='g')
        ax.set_ylabel('AGCM Diff')

        ax=axs[1,1]
        bm,pc=cplt.kemmap(aelowspr2-aemeansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='Low Anom',suppcb=True,
                          panellab='e',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aelospr2pval,alat,alon)
        bm.contour(alons,alats,aelowsp2r2-aemeansp2,levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,aelosp2r2pval,alat,alon,type='cont',color='g')

        ax=axs[1,2]
        bm,pc=cplt.kemmap(diffmult*(aelowspr2-aehighspr2),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title=diffttl,suppcb=True,
                          panellab='f',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aespr2pval, alat, alon)
        bm.contour(alons,alats,diffmult*(aelowsp2r2-aehighsp2r2),levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,aesp2r2pval,alat,alon,type='cont',color='g')

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
        bm,pc=cplt.kemmap(diffmult*(alowspr2-ahighspr2),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='Var (' + diffttl + ')',suppcb=True,
                          panellab='a',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aspr2pval, alat, alon)
        bm.contour(alons,alats,diffmult*(alowsp2r2-ahighsp2r2),levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,asp2r2pval,alat,alon,type='cont',color='g')

        ax=axs[1]
        bm,pc=cplt.kemmap(diffmult*(aelowspr2-aehighspr2),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='Const (' + diffttl + ')',suppcb=True,
                          panellab='b',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aespr2pval, alat, alon)
        bm.contour(alons,alats,diffmult*(aelowsp2r2-aehighsp2r2),levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,aesp2r2pval,alat,alon,type='cont',color='g')

        ax=axs[2] 
        bm,pc=cplt.kemmap((alowspr2-ahighspr2)-(aelowspr2-aehighspr2),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='Var-Const (' + diffttl + ')',suppcb=True,
                          panellab='c',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aressppval,alat,alon)
        bm.contour(alons,alats,(alowsp2r2-ahighsp2r2)-(aelowsp2r2-aehighsp2r2),levels=conts,
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



    # CGCM and PI-------------
    #   Remove PI signal from CGCM

    """fig,axs=plt.subplots(2,3)
    fig.set_size_inches(10,8)
    fig.subplots_adjust(wspace=0.05,hspace=0.03)

    ax=axs[0,0]
    #-- PI comp on SIC:
    bm,pc=cplt.kemmap(lowpisp-highpisp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='PI Low-High',suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pisppval, lat, lon)
    bm.contour(lons,lats,lowpisp2-highpisp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pisp2pval,lat,lon,type='cont',color='g')
    ax.set_ylabel(sttl1)
    #-- CGCM comp on SIC:
    ax=axs[0,1]
    bm,pc=cplt.kemmap(lowlesp-highlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='CGCM Low-High',suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lesppval, lat, lon)
    bm.contour(lons,lats,lowlesp2-highlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lesp2pval,lat,lon,type='cont',color='g')
    #-- CGCM-PI comp on SIC:
    ax=axs[0,2]
    junk,cgpipval = cutl.ttest_ind(meanlesp[lowidx,...].reshape((nn,nlat,nlon))-\
                                   meanlesp[highidx,...].reshape((nn,nlat,nlon)),
                                   piseasp[pilowidx,...]-piseasp[pihighidx,...])
    junk,cgpipval2 = cutl.ttest_ind(meanlesp2[lowidx,...].reshape((nn,nlat,nlon))-\
                                    meanlesp2[highidx,...].reshape((nn,nlat,nlon)),
                                    piseasp2[pilowidx,...]-piseasp2[pihighidx,...])
    junk,cgpi2pval = cutl.ttest_ind(meanlesp[low2idx,...].reshape((nn,nlat,nlon))-\
                                    meanlesp[high2idx,...].reshape((nn,nlat,nlon)),
                                    piseasp[pilow2idx,...]-piseasp[pihigh2idx,...])
    junk,cgpi2pval2 = cutl.ttest_ind(meanlesp2[low2idx,...].reshape((nn,nlat,nlon))-\
                                     meanlesp2[high2idx,...].reshape((nn,nlat,nlon)),
                                     piseasp2[pilow2idx,...]-piseasp2[pihigh2idx,...])

    bm,pc=cplt.kemmap((lowlesp-highlesp)-(lowpisp-highpisp),lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='CGCM-PI (Low-High)',suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,cgpipval, lat, lon)
    bm.contour(lons,lats,(lowlesp2-highlesp2)-(lowpisp2-highpisp2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,cgpipval2,lat,lon,type='cont',color='g')

    ax=axs[1,0]
    # PI comp on SAT:
    bm,pc=cplt.kemmap(low2pisp-high2pisp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='PI Low-High',suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pi2sppval, lat, lon)
    bm.contour(lons,lats,low2pisp2-high2pisp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pi2sp2pval,lat,lon,type='cont',color='g')
    ax.set_ylabel(sttl2)
    ax=axs[1,1]
    # CGCM comp on SAT:
    bm,pc=cplt.kemmap(low2lesp-high2lesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='CGCM Low-High',suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,le2sppval, lat, lon)
    bm.contour(lons,lats,low2lesp2-high2lesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,le2sp2pval,lat,lon,type='cont',color='g')
    #-- CGCM-PI comp on SAT:
    ax=axs[1,2]
    bm,pc=cplt.kemmap((low2lesp-high2lesp)-(low2pisp-high2pisp),lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='CGCM-PI (Low-High)',suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,cgpi2pval, lat, lon)
    bm.contour(lons,lats,(low2lesp2-high2lesp2)-(low2pisp2-high2pisp2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,cgpi2pval2,lat,lon,type='cont',color='g')

    cplt.add_colorbar(fig,pc,orientation='horizontal')
    if printtofile:
        if addsig:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_'+ fieldr+regionr + '_' +fieldr2+regionr2+ sear + '_CGCM-PI' + prstr + '.png',dpi=400) 
        else:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr+regionr+ '_' +fieldr2+regionr2+ sear + '_CGCM-PI' + prstr + '.pdf') 
    """

    # CGCM and AGCM--------------
    #   composite on Eurasian SAT

    fig,axs=plt.subplots(2,3)
    fig.set_size_inches(10,8)
    fig.subplots_adjust(wspace=0.05,hspace=0.03)

    ax=axs[0,0]
    bm,pc=cplt.kemmap(highlespr2-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='a',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehispr2pval,lat,lon)
    bm.contour(lons,lats,highlesp2r2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lehisp2r2pval,lat,lon,type='cont',color='g')

    ax.set_ylabel('CGCM Diff')
    ax=axs[0,1]
    bm,pc=cplt.kemmap(lowlespr2-meanlesp,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='b',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lelospr2pval,lat,lon)
    bm.contour(lons,lats,lowlesp2r2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lelosp2r2pval,lat,lon,type='cont',color='g')

    ax=axs[0,2]
    bm,pc=cplt.kemmap(diffmult*(lowlespr2-highlespr2),lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='c',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lespr2pval, lat, lon)
    bm.contour(lons,lats,diffmult*(lowlesp2r2-highlesp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lesp2r2pval,lat,lon,type='cont',color='g')

    ax=axs[1,0]
    bm,pc=cplt.kemmap(ahighspr2-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahispr2pval,alat,alon)
    bm.contour(alons,alats,ahighsp2r2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,ahisp2r2pval,alat,alon,type='cont',color='g')
    ax.set_ylabel('AGCM Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(alowspr2-ameansp,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,alospr2pval,alat,alon)
    bm.contour(alons,alats,alowsp2r2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,alosp2r2pval,alat,alon,type='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(alowspr2-ahighspr2),alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aspr2pval, alat, alon)
    bm.contour(alons,alats,diffmult*(alowsp2r2-ahighsp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,asp2r2pval,alat,alon,type='cont',color='g')

    cplt.add_colorbar(fig,pc,orientation='horizontal')

    if printtofile:
        if addsig:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr2+regionr2+ sear + '_CGCMAGCM' + prstr + '.png',dpi=400) 
        else:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr2+regionr2+ sear + '_CGCMAGCM' + prstr + '.pdf') 
