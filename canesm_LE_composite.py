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
from matplotlib import gridspec
import matplotlib.font_manager as fm
import scipy.io as sio
import datetime as datetime
import string as string

printtofile=False

#dataloaded=True
loadmat=False; 
when='14:51:28.762886'; styearsR = [ 8.,  7.,  2.,  8.,  8.] # variable SIC styears
styearsE=[ 4.,  1.,  7.,  3.,  1.]; styearsN=[1.] #when for these: 17:01:16.908687
styearPI = 0 # PI styear
#<coldeur for PI> when='14:43:36.586252' 
anomyearsPI = [[79, 26],
        [58, 17],
        [79, 41],
        [28, 24],
        [37, 80],
        [25, 48],
        [30, 49],
        [33, 85],
        [51, 43],
        [82,  6],
        [62, 34],
        [ 3, 17],
        [56, 63],
        [67,  4],
        [29, 73],
        [74,  0],
        [28,  8],
        [46, 56],
        [14, 76],
        [72, 37],
        [88,  4],
        [31, 56],
        [31,  4],
        [40,  0],
        [20, 49],
        [21, 81],
        [68, 56],
        [77, 37],
        [18,  1],
        [82, 26],
        [78, 55],
        [22, 47],
        [78, 16],
        [54, 76],
        [ 5, 47],
        [58, 25],
        [49, 20],
        [64, 36],
        [34,  2],
        [62, 18],
        [89,  3],
        [25, 10],
        [10, 30],
        [84, 55],
        [88, 18],
        [87, 50],
        [68, 26],
        [ 5, 67],
        [77, 13],
        [74, 61]]




saveascii=False
savemat=True
dofigures=False
local=True
addsig=False
compagcm=True # compare R sims and E sims

cisiglevel=0.05
siglevel=0.05

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
fieldsp='tas'; ncfieldsp='tas'; compsp='Amon'; spkey='sat'
cminsp=-3; cmaxsp=3 # for colors
cminspbig=-3; cmaxspbig=3 # for colors

cminspa=-1; cmaxspa=1 # for colors for AGCM
cminice=-10; cmaxice=10 # colors for ice

# spatial field2 in contours
leconvsp2=1
fieldsp2='zg50000.00'; ncfieldsp2='zg'; compsp2='Amon'; sp2key='z500'
cminsp2=-30; cmaxsp2=30 # to calc contour interval
cminsp2big=-30; cmaxsp2big=30 # to calc contour interval

cminsp2a=-10; cmaxsp2a=10 # to calc contour interval for AGCM


# COMPOSITE ON THESE VALUES
leconvr=leconvr2=leconvr3=1

# regional avg field 1
fieldr='sic'; ncfieldr='sic'; compr='OImon'; regionr='bksmori'; 
r1str='BKS SIC'; r1strlong='Barents/Kara sea ice concentration'; r1units='%'; r1key='bkssic'


#fieldr='turb'; ncfieldr='turb'; compr='Amon'; regionr='bksmori'; #@@@
#r1str='BKS Turb'; r1strlong='Barents/Kara turbulent heat flux'; r1units='W/m2'; r1key='bksturb'



# regional avg field 2
# cooling=high heights
fieldr2='tas'; ncfieldr2='tas'; compr2='Amon'; regionr2='eurasiamori'; 
r2str='Eur SAT'; r2strlong='Eurasian surface air temperature'; r2units='$^\circ$C'; r2key='eursat' #leconvr2=-1

# regional avg field 3
fieldr3='zg50000.00'; ncfieldr3='zg'; compr3='Amon'; 
regionr3='bksmori'; r3str='BKS Z500'; r3key='bksz500'
diffttl3='High-Low'; diffmult3=-1; r3units='m'
#leconvr=-1; #leconvr2=-1; #both conv -1 to get figs to show low-high equal to high heights and cold continent.

#fieldr3='turb'; ncfieldr3='turb'; compr3='Amon'; regionr3='bksmori'; #@@@
#r3str='BKS Turb'; r3strlong='Barents/Kara turbulent heat flux'; r3units='W/m2'; r3key='bksturb'


sttl1='Comp on ' + r1str
sttl2='Comp on ' + r2str
sttl3='Comp on ' + r3str

fdictsp = {'field': fieldsp, 'ncfield': ncfieldsp, 'comp': compsp}
fdictsp2 = {'field': fieldsp2, 'ncfield': ncfieldsp2, 'comp': compsp2}
if fieldr=='turb':
    fdictice = {'field': 'turb', 'ncfield': 'turb', 'comp': 'Amon'} #@@@ testing. misleading names.
else:
    fdictice = {'field': 'sic', 'ncfield': 'sic', 'comp':'OImon'}

fdictr = {'field': fieldr+regionr, 'ncfield': ncfieldr, 'comp': compr}
fdictr2 = {'field': fieldr2+regionr2, 'ncfield': ncfieldr2, 'comp': compr2}
fdictr3 = {'field': fieldr3+regionr3,'ncfield': ncfieldr3, 'comp': compr3}

# used for file loading & figures
#regions=('bkssic','eursat','bksz500')
regions=('eursat',r1key,r3key) # SWAP order
#fields=('ice','sat','z500')
fields=('sat','ice','z500') # SWAP order


lat=le.get_lat(local=local)
lon=le.get_lon(local=local)
nlat=len(lat); nlon=len(lon)
alat=con.get_t63lat(); alon=con.get_t63lon()
alons, alats = np.meshgrid(alon,alat)

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
    elif field=='turb': simfield='turb'; simncfield='turb'; # the sim var names are placeholders
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
        print 'highidx vals: ' + str(rfld[highidx])
        print 'lowidx vals: ' + str(rfld[lowidx])

    # low vs high
    splowt=spfld[lowidx,...].reshape(rshapenn) # keep sample dim
    sphight=spfld[highidx,...].reshape(rshapenn) # keep sample dim
    (lesptstat, lesppval) = cutl.ttest_ind(splowt, sphight)

    # low vs mean
    (lelosptstat, lelosppval) = cutl.ttest_ind(splowt, 
                                               spfld.reshape(rshapeens))
    # high vs mean
    (lehisptstat, lehisppval) = cutl.ttest_ind(sphight, 
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
        st=sphight[...,-1]
        sphight=np.dstack((sphight,st[...,None]))
        st=splowt[...,-1]
        splowt=np.dstack((splowt,st[...,None]))

    compdt = {'highsp': highlesp, 'lowsp': lowlesp, 'meansp': meanlesp,
              'hisppval': lehisppval, 'losppval': lelosppval, 'sppval': lesppval,
              'highidx': highidx, 'lowidx': lowidx, 'highspt': sphight, 'lowspt': splowt}

    return compdt

# end do_composite()


if loadmat:

    leallrdt={r1key:dict.fromkeys(fields),
              'eursat':dict.fromkeys(fields),
              r3key:dict.fromkeys(fields)}    

    matbase='pymatfiles/LE_composites_'
    for rkey in regions:
        for fkey in fields:
            matname = matbase + rkey + '_' + fkey + '_' + when + '.mat'
            leallrdt[rkey][fkey]=sio.loadmat(matname,squeeze_me=True)

    piallrdt={r1key:dict.fromkeys(fields),
              'eursat':dict.fromkeys(fields),
              r3key:dict.fromkeys(fields)}    

    matbase='pymatfiles/PI_composites_'
    for rkey in regions:
        for fkey in fields:
            matname = matbase + rkey + '_' + fkey + '_' + when + '.mat'
            piallrdt[rkey][fkey]=sio.loadmat(matname,squeeze_me=True)

    aallrdt={r1key:dict.fromkeys(fields),
              'eursat':dict.fromkeys(fields),
              r3key:dict.fromkeys(fields)}    

    matbase='pymatfiles/AGCMRsims_composites_'
    for rkey in regions:
        for fkey in fields:
            matname = matbase + rkey + '_' + fkey + '_' + when + '.mat'
            aallrdt[rkey][fkey]=sio.loadmat(matname,squeeze_me=True)

    aeallrdt={r1key:dict.fromkeys(fields),
              'eursat':dict.fromkeys(fields),
              r3key:dict.fromkeys(fields)}    

    matbase='pymatfiles/AGCMEsims_composites_'
    for rkey in regions:
        for fkey in fields:
            matname = matbase + rkey + '_' + fkey + '_' + when + '.mat'
            aeallrdt[rkey][fkey]=sio.loadmat(matname,squeeze_me=True)


#if not dataloaded:
else:
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

    now = str(datetime.datetime.now().time())
    now=when

    ler1flds={'sat':lespr1dt,'z500':lesp2r1dt,'ice':leicer1dt}
    ler2flds={'sat':lespr2dt,'z500':lesp2r2dt,'ice':leicer2dt}
    ler3flds={'sat':lespr3dt,'z500':lesp2r3dt,'ice':leicer3dt}

    leallrdt={r1key:ler1flds,'eursat':ler2flds,r3key:ler3flds}    

    if savemat:
        matbase='pymatfiles/LE_composites_'
        for rkey in regions:
            for fkey in fields:
                savedt=leallrdt[rkey][fkey]

                matname = matbase + rkey + '_' + fkey + '_' + now + '.mat'
                sio.savemat(matname,savedt)
 

    # ========== PRE-IND ==============

    piseasp,styear,anomyears = load_canesmfield(fdictsp,'piControl',seasp,conv=leconvsp,
                                                local=local,styear=styearPI,anomyears=anomyearsPI)
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


    pir1flds={'sat':pispr1dt,'z500':pisp2r1dt,'ice':piicer1dt}
    pir2flds={'sat':pispr2dt,'z500':pisp2r2dt,'ice':piicer2dt}
    pir3flds={'sat':pispr3dt,'z500':pisp2r3dt,'ice':piicer3dt}

    piallrdt={r1key:pir1flds,'eursat':pir2flds,r3key:pir3flds}    

    if savemat:
        matbase='pymatfiles/PI_composites_'
        for rkey in regions:
            for fkey in fields:
                savedt=piallrdt[rkey][fkey]

                matname = matbase + rkey + '_' + fkey + '_' + now + '.mat'
                sio.savemat(matname,savedt)
        sio.savemat(matbase+'anomyears_'+now+'.mat',{'anomyears':anomyears})
        sio.savemat(matbase+'styear_'+now+'.mat',{'styear':styear})


    # === AGCM ==========

    simsE=('E1','E2','E3','E4','E5'); 
    sims=('R1','R2','R3','R4','R5');
    simsO=('NSIDC',)

    asseasp,styearsr = load_agcmfield(fieldsp,sims,seasp,styears=styearsR) # already anomalies
    asseasp2,styearsr = load_agcmfield(fieldsp2,sims,seasp,styears=styearsr) # already anomalies
    asseaice,styearsr = load_agcmfield('sic',sims,seasp,styears=styearsr,conv=100) # want the SICN pattern

    assear,styearsr = load_agcmfield(fieldr,sims,sear,styears=styearsr,region=regionr,conv=leconvr) # already anomalies
    assear2,styearsr = load_agcmfield(fieldr2,sims,sear,styears=styearsr,region=regionr2,conv=leconvr2) # already anomalies
    assear3,styearsr = load_agcmfield(fieldr3,sims,sear,styears=styearsr,region=regionr3,conv=leconvr3) # already anomalies
    assear=assear # / assear.std()
    assear2=assear2 #/ assear2.std()
    assear3=assear3

    (anens,anlat,anlon)=asseasp.shape
    rshape=(anens,anlat*anlon)


    print 'AGCMR comp on BKS SIC'
    aspr1dt= do_composite(assear,asseasp,verb=True)
    asp2r1dt = do_composite(assear,asseasp2,verb=True)
    aicer1dt = do_composite(assear,asseaice,verb=True)

    print 'AGCMR comp on Eur SAT'
    aspr2dt = do_composite(assear2,asseasp,verb=True)
    asp2r2dt = do_composite(assear2,asseasp2,verb=True)
    aicer2dt = do_composite(assear2,asseaice,verb=True)

    print 'AGCMR comp on BKS Z500'
    aspr3dt = do_composite(assear3,asseasp,verb=True)
    asp2r3dt = do_composite(assear3,asseasp2,verb=True)
    aicer3dt = do_composite(assear3,asseaice,verb=True)

    ar1flds={'sat':aspr1dt,'z500':asp2r1dt,'ice':aicer1dt}
    ar2flds={'sat':aspr2dt,'z500':asp2r2dt,'ice':aicer2dt}
    ar3flds={'sat':aspr3dt,'z500':asp2r3dt,'ice':aicer3dt}

    aallrdt={r1key:ar1flds,'eursat':ar2flds,r3key:ar3flds}    

    if savemat:
        matbase='pymatfiles/AGCMRsims_composites_'
        for rkey in regions:
            for fkey in fields:
                savedt=aallrdt[rkey][fkey]

                matname = matbase + rkey + '_' + fkey + '_' + now + '.mat'
                sio.savemat(matname,savedt)
        sio.savemat(matbase+'styears_'+now+'.mat',{'styears':styearsr})


    if compagcm: # if compare AGCM ensembles

        # ===  E sims:
        aesseasp,styearse = load_agcmfield(fieldsp,simsE,seasp,styears=styearsE) # already anomalies
        aesseasp2,styearse = load_agcmfield(fieldsp2,simsE,seasp,styears=styearse) # already anomalies
        aesseaice,styearse = load_agcmfield('sic',simsE,seasp,styears=styearse,conv=100) # want the SICN pattern

        aessear,styearse = load_agcmfield(fieldr,simsE,sear,styears=styearse,region=regionr,conv=leconvr) # already anomalies
        aessear2,styearse = load_agcmfield(fieldr2,simsE,sear,styears=styearse,region=regionr2,conv=leconvr2) # already anomalies
        aessear3,styearse = load_agcmfield(fieldr3,simsE,sear,styears=styearse,region=regionr3,conv=leconvr3) # already anomalies
        aessear=aessear # / assear.std()
        aessear2=aessear2 #/ assear2.std()
        aessear3=aessear3


        print 'AGCME comp on BKS SIC'
        aespr1dt= do_composite(aessear,aesseasp,verb=True)
        aesp2r1dt = do_composite(aessear,aesseasp2,verb=True)
        aeicer1dt = do_composite(aessear,aesseaice,verb=True)

        print 'AGCME comp on Eur SAT'
        aespr2dt = do_composite(aessear2,aesseasp,verb=True)
        aesp2r2dt = do_composite(aessear2,aesseasp2,verb=True)
        aeicer2dt = do_composite(aessear2,aesseaice,verb=True)

        print 'AGCME comp on BKS Z500'
        aespr3dt = do_composite(aessear3,aesseasp,verb=True)
        aesp2r3dt = do_composite(aessear3,aesseasp2,verb=True)
        aeicer3dt = do_composite(aessear3,aesseaice,verb=True)

        aer1flds={'sat':aespr1dt,'z500':aesp2r1dt,'ice':aeicer1dt}
        aer2flds={'sat':aespr2dt,'z500':aesp2r2dt,'ice':aeicer2dt}
        aer3flds={'sat':aespr3dt,'z500':aesp2r3dt,'ice':aeicer3dt}

        aeallrdt={r1key:aer1flds,'eursat':aer2flds,r3key:aer3flds}    

        if savemat:
            matbase='pymatfiles/AGCMEsims_composites_'
            for rkey in regions:
                for fkey in fields:
                    savedt=aeallrdt[rkey][fkey]

                    matname = matbase + rkey + '_' + fkey + '_' + now + '.mat'
                    sio.savemat(matname,savedt)
            sio.savemat(matbase+'styears_'+now+'.mat',{'styears':styearse})
        

        # AGCM RESIDUAL pvals (var vs const)
        #(aressptstat,aressppval) = cutl.ttest_ind(asseasp[alow2idx,...]- asseasp[ahigh2idx,...],
        #                                          aesseasp[aelow2idx,...]- aesseasp[aehigh2idx,...])
        #(aressp2tstat,aressp2pval) = cutl.ttest_ind(asseasp2[alow2idx,...]- asseasp2[ahigh2idx,...],
        #                                          aesseasp2[aelow2idx,...]- aesseasp2[aehigh2idx,...])


        
        # ===  Obs sims (NSIDC):
        aosseasp,styearso = load_agcmfield(fieldsp,simsO,seasp,styears=styearsN) # already anomalies
        aosseasp2,styearso = load_agcmfield(fieldsp2,simsO,seasp,styears=styearso) # already anomalies
        aosseaice,styearso = load_agcmfield('sic',simsO,seasp,styears=styearso,conv=100) # want the SICN pattern

        aossear,styearso = load_agcmfield(fieldr,simsO,sear,styears=styearso,region=regionr,conv=leconvr) # already anomalies
        aossear2,styearso = load_agcmfield(fieldr2,simsO,sear,styears=styearso,region=regionr2,conv=leconvr2) # already anomalies
        aossear3,styearso = load_agcmfield(fieldr3,simsO,sear,styears=styearso,region=regionr3,conv=leconvr3) # already anomalies
        aossear=aossear # / assear.std()
        aossear2=aossear2 #/ assear2.std()
        aossear3=aossear3

        if savemat:
            matbase='pymatfiles/AGCMOsims_composites_' # @@ not really a composite. just save startyr
            sio.savemat(matbase+'styears_'+now+'.mat',{'styears':styearso})

    if saveascii:
        # write ascii file for john:
        agcmout1=np.hstack((assear, aossear, aessear))*100 # column 1 (BKS SIC: variable SIC, NSIDC SIC, mean SIC)
        agcmout2=np.hstack((assear2, aossear2, aessear2)) # column 2 (Eur SAT)
        agcmout3=np.hstack((assear3, aossear3, aessear3)) # column 3 (BKS Z500)
        agcmout = np.vstack((agcmout1,agcmout2,agcmout3)).T
        np.savetxt('agcmout3.ascii',agcmout,delimiter='\t')

        # 
        #cgcmout1=np.hstack((lefldcnd,lefldncnd,lefldmcnd)) # (ALL forcing, NAT, AERO/MISC)
        #cgcmout2=np.hstack((lefld2,lefld2n,lefld2m))
        #cgcmout3=np.hstack((lefld1,lefld1n,lefld1m))
        #cgcmout = np.vstack((cgcmout1,cgcmout2,cgcmout3)).T
        #np.savetxt('cgcmout.ascii',cgcmout,delimiter='\t')        




# For each composite (r1, r2, r3), compute the BKS SIC average:
allcasedt = {'Preindustrial':piallrdt, 'CGCM': leallrdt, 'AGCM':aallrdt}
allcasedt = {'CGCM': leallrdt, 'AGCM': aallrdt, 'AGCM_fixed': aeallrdt}

#  Here calculate the BKS SIC associated with each composite
allregimdt={};allregimlodt={};allregimhidt={}; allregitdt={}; allregimcidt={}; allregitcidt={}
allregimlocidt={};allregimhicidt={};
fkey='ice'
print '@@@@@@@@@@@@@ calculating CI for ICE, all composites'
for ckey in allcasedt.keys():
    print '=ENS ' + ckey
    regmdt={}; regmlodt={}; regmhidt={}; regtotdt={}; regmcidt={}; regtcidt={}
    regmlocidt={}; regmhicidt={}
    allrdt=allcasedt[ckey]
    for rkey in regions:

        print '  REGION ' + rkey
        
        dt=allrdt[rkey][fkey]

        diff = dt['lowspt']-dt['highspt'] # with sample dim
        regt = cutl.calc_regtotseaicearea(diff[...,:-1],lat,lon,'bksmori')
        regm = cutl.calc_regmean(diff[...,:-1],lat,lon,'bksmori')

        # Confidence interval: 2.5-97.5% interval on the mean
        regmci = sp.stats.t.interval(1-cisiglevel,len(regm)-1,
                                     loc=regm.mean(axis=0),
                                     scale=regm.std(axis=0)/np.sqrt(len(regm)))
        regtci = sp.stats.t.interval(1-cisiglevel,len(regt)-1,
                                     loc=regt.mean(axis=0),
                                     scale=regt.std(axis=0)/np.sqrt(len(regt)))
        print '   DIFF: ' + str(regm.mean(axis=0)) + ', CI: ' + str(regmci) # @@

        # @@ CONFIDENCE interval on the low! (for John to compare)
        regmlo = cutl.calc_regmean(dt['lowspt'][...,:-1],lat,lon,'bksmori')
        regmhi = cutl.calc_regmean(dt['highspt'][...,:-1],lat,lon,'bksmori')
        regmloci = sp.stats.t.interval(1-cisiglevel,len(regmlo)-1,
                                       loc=regmlo.mean(axis=0),
                                       scale=regmlo.std(axis=0)/np.sqrt(len(regmlo)))
        regmhici = sp.stats.t.interval(1-cisiglevel,len(regmhi)-1,
                                       loc=regmhi.mean(axis=0),
                                       scale=regmhi.std(axis=0)/np.sqrt(len(regmhi)))
        
        
        print '   LO: ' + str(regmlo.mean(axis=0)) + ', CI: ' + str(regmloci) # @@
        #print '   LO vals: ' + str(regmlo)
        print '   HI: ' + str(regmhi.mean(axis=0)) + ', CI: ' + str(regmhici) # @@
        junk,regmdiffpval = cutl.ttest_ind(regmlo,regmhi)
        print '   LO v HI pval: ' + str(regmdiffpval)

        regmdt[rkey]=regm.mean(axis=0)
        regmlodt[rkey]=regmlo.mean(axis=0)
        regmhidt[rkey]=regmhi.mean(axis=0)
        regtotdt[rkey]=regt.mean(axis=0) # not sure which one i want
        regmcidt[rkey]=regmci
        regmlocidt[rkey]=regmloci
        regmhicidt[rkey]=regmhici
        regtcidt[rkey]=regtci

    allregimdt[ckey]=regmdt # i for ice
    allregimlodt[ckey]=regmlodt # i for ice
    allregimhidt[ckey]=regmhidt # i for ice
    allregitdt[ckey]=regtotdt
    allregimcidt[ckey]=regmcidt
    allregimlocidt[ckey]=regmlocidt
    allregimhicidt[ckey]=regmhicidt
    allregitcidt[ckey]=regtcidt

allregimdf=pd.DataFrame(allregimdt)
allregitdf=pd.DataFrame(allregitdt)
allregimcidf=pd.DataFrame(allregimcidt)
allregitcidf=pd.DataFrame(allregitcidt)

#  Here calculate the Eur SAT associated with each composite
allregspmdt={}; allregspmcidt={};

fkey='sat'
print '@@@@@@@@@@@@@ calculating CI for SAT, all composites'
for ckey in allcasedt.keys():
    print '=ENS ' + ckey

    regmdt={}; regtotdt={}; regmcidt={}; regtcidt={}
    allrdt=allcasedt[ckey]
    for rkey in regions:
        print '  REGION ' + rkey
        dt=allrdt[rkey][fkey]

        diff = dt['lowspt']-dt['highspt']
        regm = cutl.calc_regmean(diff[...,:-1],lat,lon,'eurasiamori')

        # Confidence interval: 2.5-97.5% interval on the mean
        regmci = sp.stats.t.interval(1-cisiglevel,len(regm)-1,
                                     loc=regm.mean(axis=0),
                                     scale=regm.std(axis=0)/np.sqrt(len(regm)))
        print '   DIFF: ' + str(regm.mean(axis=0)) + ', CI: ' + str(regmci) # @@

        # @@ CONFIDENCE interval on the low! (for John to compare)
        regmlo = cutl.calc_regmean(dt['lowspt'][...,:-1],lat,lon,'eurasiamori')
        regmhi = cutl.calc_regmean(dt['highspt'][...,:-1],lat,lon,'eurasiamori')
        regmloci = sp.stats.t.interval(1-cisiglevel,len(regmlo)-1,
                                       loc=regmlo.mean(axis=0),
                                       scale=regmlo.std(axis=0)/np.sqrt(len(regmlo)))
        regmhici = sp.stats.t.interval(1-cisiglevel,len(regmhi)-1,
                                       loc=regmhi.mean(axis=0),
                                       scale=regmhi.std(axis=0)/np.sqrt(len(regmhi)))
        print '   LO: ' + str(regmlo.mean(axis=0)) + ', CI: ' + str(regmloci) # @@
        #print '   LO vals: ' + str(regmlo)
        print '   HI: ' + str(regmhi.mean(axis=0)) + ', CI: ' + str(regmhici) # @@
        junk,regmdiffpval = cutl.ttest_ind(regmlo,regmhi)
        print '   LO v HI pval: ' + str(regmdiffpval)

        regmdt[rkey]=regm.mean(axis=0)
        regmcidt[rkey]=regmci

    allregspmdt[ckey]=regmdt # sp for spatial 1 (SAT)
    allregspmcidt[ckey]=regmcidt

allregspmdf=pd.DataFrame(allregspmdt)
allregspmcidf=pd.DataFrame(allregspmcidt)


#  Here calculate the BKS Z500 associated with each composite
allregsp2mdt={}; allregsp2mcidt={};allregsp2mlodt={}; allregsp2mlocidt={};
allregsp2mhidt={}; allregsp2mhicidt={};
fkey='z500'
print '@@@@@@@@@@@@@ calculating CI for Z500, all composites'
for ckey in allcasedt.keys():
    print '=ENS ' + ckey # @@

    regmdt={};regmlodt={};regmhidt={}; regtotdt={}; regmcidt={}; regmlocidt={};regmhicidt={}
    allrdt=allcasedt[ckey]
    for rkey in regions:
        
        print '  REGION ' + rkey # @@

        dt=allrdt[rkey][fkey]

        diff = dt['lowspt']-dt['highspt']
        regm = cutl.calc_regmean(diff[...,:-1],lat,lon,'bksmori')
        
        # Confidence interval: 2.5-97.5% interval on the mean
        regmci = sp.stats.t.interval(1-cisiglevel,len(regm)-1,
                                     loc=regm.mean(axis=0),
                                     scale=regm.std(axis=0)/np.sqrt(len(regm)))
        print '   DIFF: ' + str(regm.mean(axis=0)) + ', CI: ' + str(regmci) # @@

        # @@ CONFIDENCE interval on the low! (for John to compare)
        regmlo = cutl.calc_regmean(dt['lowspt'][...,:-1],lat,lon,'bksmori')
        regmhi = cutl.calc_regmean(dt['highspt'][...,:-1],lat,lon,'bksmori')
        regmloci = sp.stats.t.interval(1-cisiglevel,len(regmlo)-1,
                                       loc=regmlo.mean(axis=0),
                                       scale=regmlo.std(axis=0)/np.sqrt(len(regmlo)))
        regmhici = sp.stats.t.interval(1-cisiglevel,len(regmhi)-1,
                                       loc=regmhi.mean(axis=0),
                                       scale=regmhi.std(axis=0)/np.sqrt(len(regmhi)))
        print '   LO: ' + str(regmlo.mean(axis=0)) + ', CI: ' + str(regmloci) # @@
        #print '   LO vals: ' + str(regmlo)
        print '   HI: ' + str(regmhi.mean(axis=0)) + ', CI: ' + str(regmhici) # @@
        junk,regmdiffpval = cutl.ttest_ind(regmlo,regmhi)
        print '   LO v HI pval: ' + str(regmdiffpval)

        regmdt[rkey]=regm.mean(axis=0)
        regmcidt[rkey]=regmci
        regmlodt[rkey]=regmlo.mean(axis=0)
        regmlocidt[rkey]=regmloci
        regmhidt[rkey]=regmhi.mean(axis=0)
        regmhicidt[rkey]=regmhici

    allregsp2mdt[ckey]=regmdt # sp for spatial 1 (SAT)
    allregsp2mcidt[ckey]=regmcidt
    allregsp2mlodt[ckey]=regmlodt # sp for spatial 1 (SAT)
    allregsp2mlocidt[ckey]=regmlocidt
    allregsp2mhidt[ckey]=regmhidt # sp for spatial 1 (SAT)
    allregsp2mhicidt[ckey]=regmhicidt

allregsp2mdf=pd.DataFrame(allregsp2mdt)
allregsp2mcidf=pd.DataFrame(allregsp2mcidt)



fontP = fm.FontProperties()
fontP.set_size('small')

lons, lats = np.meshgrid(lon,lat)
cmlen=15.
incr = (cmaxsp2-cminsp2) / (cmlen)
conts = np.arange(cminsp2,cmaxsp2+incr,incr)

incrbig = (cmaxsp2big-cminsp2big) / (cmlen)
contsbig = np.arange(cminsp2big,cmaxsp2big+incrbig,incrbig)

stxx=1
xx=np.arange(stxx,(stxx+len(allregimdf.keys())-1)*1.2,1.2)
#regs=('bkssic','eursat','bksz500')
#regs=('eursat','bkssic','bksz500') # SWAP order. @@ Actually, use regions instead

multfacs=(1,1,-1) # multiply the comp on z500 by -1 to get high over bks
enss=('CGCM','Preindustrial','AGCM')
enss=('CGCM','AGCM','AGCM_fixed')
mkrs=('s','o','^')
clrs=('r','b','0.4')
#clrs=('.4','.4','.4')


# prepare legend entries
lgs=(mlines.Line2D([],[],color=clrs[0],linewidth=2),
     mlines.Line2D([],[],color=clrs[1],linewidth=2),
     mlines.Line2D([],[],color=clrs[2],linewidth=2)) #,linestyle='none',marker=mkrs[0]),
     #mlines.Line2D([],[],color=clrs[1],linestyle='none',marker=mkrs[1]),
     #mlines.Line2D([],[],color=clrs[2],linestyle='none',marker=mkrs[2])) 
#lgstrs=('BKS SIC','Eur SAT','BKS Z500')  
lgstrs=(r2str,r1str,r3str)  # SWAP order
        

lespr1dt=leallrdt[r1key][spkey]
lesp2r1dt=leallrdt[r1key][sp2key]
lespr2dt=leallrdt[r2key][spkey]
lesp2r2dt=leallrdt[r2key][sp2key]
lespr3dt=leallrdt[r3key][spkey]
lesp2r3dt=leallrdt[r3key][sp2key]


fig=plt.figure(figsize=(12,9))
gs1 = gridspec.GridSpec(1, 6)
gs1.update(top=0.97,bottom=0.55,left=0.08,right=0.92,wspace=0.04)
ax=plt.subplot(gs1[0,0:2]) 

# CGCM comp on BKS SIC
bm,pc=cplt.kemmap(lespr1dt['lowsp']-lespr1dt['highsp'],alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='',suppcb=True,
                  panellab='a',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lespr1dt['sppval'], alat, alon)
bm.contour(alons,alats,lesp2r1dt['lowsp']-lesp2r1dt['highsp'],levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2r1dt['sppval'],alat,alon,sigtype='cont',color='g')

ax=plt.subplot(gs1[0,2:4]) # CGCM comp on Eur SAT
bm,pc=cplt.kemmap(lespr2dt['lowsp']-lespr2dt['highsp'],alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='',suppcb=True,
                      panellab='b',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lespr2dt['sppval'], alat, alon)
bm.contour(alons,alats,lesp2r2dt['lowsp']-lesp2r2dt['highsp'],levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2r2dt['sppval'],alat,alon,sigtype='cont',color='g')

ax=plt.subplot(gs1[0,4:])
bm,pc=cplt.kemmap(-1*(lespr3dt['lowsp']-lespr3dt['highsp']),alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='',suppcb=True,
                      panellab='c',lcol='0.2')
if addsig:
    cplt.addtsigm(bm,lespr3dt['sppval'], alat, alon)
bm.contour(alons,alats,-1*(lesp2r3dt['lowsp']-lesp2r3dt['highsp']),levels=conts,
           colors='0.5',linewidths=1,latlon=True)
if addsig:
    cplt.addtsigm(bm,lesp2r3dt['sppval'],alat,alon,sigtype='cont',color='g')

gs2 = gridspec.GridSpec(1, 6)
gs2.update(top=0.44,left=0.2,right=0.8,wspace=0.2)
ax=plt.subplot(gs2[0,0:3])
xincr=0.2
stagger=-xincr
for rii,reg in enumerate(regions):
    print reg
    
    for eii,ens in enumerate(enss):
        print '  ' + ens
        
        mult=multfacs[rii]
        xpos=xx[eii]+stagger
        ax.plot(xpos,mult*allregimdf[ens][reg],marker=mkrs[rii],
                color=clrs[rii],mec=clrs[rii],linestyle='none',markersize=10)
        ci=allregimcidf[ens][reg]
        #print mult, (xpos,xpos), ci
        ax.plot((xpos,xpos),mult*np.array(ci),marker='_',
                mew=2,markersize=10,linewidth=2,color=clrs[rii])
        ax.axvspan(xx[eii]-xincr-0.1,xx[eii]+xincr+0.1,color='0.8',alpha=0.5)
    
    stagger+=0.2

ax.set_ylabel('Change in ' + r1str + ' (' + r1units + ')')
ax.axhline(y=0,color='k',linestyle='--')
ax.set_xlim((stxx-0.7,xx[-1]+0.7))
ax.legend(lgs,lgstrs,loc='lower right',fancybox=True,frameon=False,prop=fontP)
ax.set_xticks(xx)
ax.set_xticklabels(enss,rotation=25)
ax.annotate('d',xy=(-0.02,1.04),
            xycoords='axes fraction',fontsize=16,fontweight='bold')

ax=plt.subplot(gs2[0,3:])
stagger=-0.2
for rii,reg in enumerate(regions):
    print reg
    
    for eii,ens in enumerate(enss):
        print '  ' + ens
        mult=multfacs[rii]
        xpos=xx[eii]+stagger
        ax.plot(xpos,mult*allregspmdf[ens][reg],marker=mkrs[rii],
                color=clrs[rii],mec=clrs[rii],linestyle='none',markersize=10)
        ci=allregspmcidf[ens][reg]
        #print mult,(xpos), ci
        ax.plot((xpos,xpos),mult*np.array(ci),
                marker='_',mew=2,markersize=10,linewidth=2,color=clrs[rii])
        ax.axvspan(xx[eii]-xincr-0.1,xx[eii]+xincr+0.1,color='0.8',alpha=0.5)

    stagger+=0.2

ax.set_ylabel('Change in ' + r2str + ' (' + r2units +')')
ax.yaxis.set_label_position('right')
ax.yaxis.tick_right()
ax.axhline(y=0,color='k',linestyle='--')
ax.set_xlim((stxx-0.7,xx[-1]+0.7))
ax.set_xticks(xx)
ax.set_xticklabels(enss,rotation=25)
ax.annotate('e',xy=(-0.02,1.04),
            xycoords='axes fraction',fontsize=16,fontweight='bold')

cbar_ax,cbar = cplt.add_colorbar(fig,pc,orientation='horizontal',
                                 pos=[.25,.53, .5,.03],
                                 label='Change in surface air temperature, $\Delta$SAT (' + r2units +')')
#cbar.add_lines(pc)
if addsig:
    prstr='sig'
else:
    prstr=''
if printtofile:
    fig.savefig('Figure4_draft_4panel_CGCMmaps_compavgs' + prstr + '2.pdf')
    fig.savefig('Figure4_draft_4panel_CGCMmaps_compavgs' + prstr + '2.png',dpi=400)



# VERSION 2: SWAP THINGS AROUND: ==================
# allflddt is organized by spatial field
allflddt={'ice':allregimdf, 'sat':allregspmdf, 'z500':allregsp2mdf}
allcidt={'ice':allregimcidf, 'sat':allregspmcidf, 'z500':allregsp2mcidf}

#fylims=((-12,4),(-2.5,1),(-15,75))
fylims=((-3,1),(-12,4),(-15,75)) # SWAP
fyticklabs=(('','-2','','-1','','0','','1'),
            ('-12','','-8','','-4','','0','','4'),
           ('','0','20','40','60')) # SWAP
mkrs2=('o','^','s')

#gs2 = gridspec.GridSpec(1, 9)
#gs2.update(top=0.44,left=0.02,right=0.98,wspace=0.2)


# loop "region" composited on
#     loop spatial field (ice,sat,z500)
#          loop ensemble

xx=np.arange(0,3); stxx=0
xincr=0.2

plabs = list(string.ascii_lowercase)
#reglabs=('BKS SIC','Eur SAT','BKS Z500')
#reglabunits=('%','$^\circ$C','m')
reglabs=(r2str,r1str,r3str) # SWAP
reglabunits=('$^\circ$C','%','m') # SWAP



print '---- plotting v2 of summary ----'
pii=0
fig = plt.figure(figsize=(14,9.5))
gs1 = gridspec.GridSpec(1, 3)
gs1.update(top=0.97,bottom=0.58,left=0.06,right=0.95,wspace=0.07)

cii=0; ckey='CGCM'; pii=0
for rii,rkey in enumerate(regions):
    mult=multfacs[rii]
    print '  ' + rkey + ' (mult=' + str(mult) + ')'
    ax=plt.subplot(gs1[cii,rii])
    dt=allcasedt[ckey][rkey]

    bm,pc=cplt.kemmap(mult*(dt[spkey]['lowsp']-dt[spkey]['highsp']),alat,alon,
                      ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='',suppcb=True,
                      panellab=plabs[pii],lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,dt[spkey]['sppval'], alat, alon)
    bm.contour(alons,alats,mult*(dt[sp2key]['lowsp']-dt[sp2key]['highsp']),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,dt[sp2key]['sppval'],alat,alon,sigtype='cont',color='g')

    if ax.is_first_row(): ax.set_title(reglabs[rii])
    pii+=1

outer_grid = gridspec.GridSpec(1, 3, wspace=0.25)
outer_grid.update(top=0.44,bottom=0.08,left=0.06,right=0.95)

for rii,rkey in enumerate(regions):

    print rkey

    inner_grid = gridspec.GridSpecFromSubplotSpec(1, 3,
                                                  subplot_spec=outer_grid[rii], 
                                                  wspace=0.46)    
    mult=multfacs[rii]
    for fii,fkey in enumerate(fields):
        print '  ' + fkey

        flddf=allflddt[fkey]
        cidf=allcidt[fkey]

        stagger=-xincr

        ax = plt.Subplot(fig, inner_grid[fii])
        #if fii==1: ax.set_title(reglabs[rii])#'Comp on ' + rkey)

        #ax=plt.subplot(gs2[0,rii+fii])
        for eii,ens in enumerate(enss):

            print '    ' + str(flddf[ens][rkey]) + ' (mult=' + str(mult) + ')'
            xpos=xx[eii]+stagger
            ax.plot(xpos,mult*flddf[ens][rkey],marker=mkrs2[eii],
                    color=clrs[eii],mec=clrs[eii],linestyle='none',markersize=10)
            ci=cidf[ens][rkey]
            ax.plot((xpos,xpos),mult*np.array(ci),marker='_',
                    mew=2,markersize=10,linewidth=2,color=clrs[eii])
            ax.set_xlim((stxx-2,xx[-1]+2))
            ax.axvspan(stxx-1,xx[-1]+1,color='0.8',alpha=0.5)
            ax.axhline(y=0,linestyle='--',color='k')
            stagger+=0.2

        ax.set_xticks(xx)
        ax.set_xticklabels(('',))
        ax.set_ylim(fylims[fii])
        ax.set_yticklabels(fyticklabs[fii])
        ax.set_xlabel('$\Delta$' +reglabs[fii],rotation=25)# fkey)
        ax.set_ylabel('('+reglabunits[fii]+')')
        ax.tick_params(bottom='off',top='off')
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(right='off')
        
        fig.add_subplot(ax)
        if ax.is_first_col():
            ax.annotate(plabs[pii],xy=(-0.02,1.04),
                        xycoords='axes fraction',fontsize=16,fontweight='bold')
            pii+=1

cbar_ax,cbar = cplt.add_colorbar(fig,pc,orientation='horizontal',
                                 pos=[.25,.53, .5,.03],
                                 label='Change in surface air temperature, $\Delta$SAT (' + r2units +')')
if addsig:
    prstr='sig'
else:
    prstr=''
if printtofile:
    fig.savefig('test_fig4_scramble' + prstr +'.pdf')
    fig.savefig('test_fig4_scramble' + prstr +'.png',dpi=400)
    


# VERSION 3 of SUMMARY ====================

lgstrs=(enss)          
mkrs=('o','^','s')


print '---- plotting v3 of summary ----'

fig=plt.figure(figsize=(12,9))
gs1 = gridspec.GridSpec(1, 3)
gs1.update(top=0.97,bottom=0.55,left=0.06,right=0.94,wspace=0.07)

cii=0; ckey='CGCM'; pii=0
for rii,rkey in enumerate(regions):
    mult=multfacs[rii]
    #print '  ' + rkey + ' (mult=' + str(mult) + ')'
    ax=plt.subplot(gs1[cii,rii])
    dt=allcasedt[ckey][rkey]

    bm,pc=cplt.kemmap(mult*(dt[spkey]['lowsp']-dt[spkey]['highsp']),alat,alon,
                      ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='',suppcb=True,
                      panellab=plabs[pii],lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,dt[spkey]['sppval'], alat, alon)
    bm.contour(alons,alats,mult*(dt[sp2key]['lowsp']-dt[sp2key]['highsp']),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,dt[sp2key]['sppval'],alat,alon,sigtype='cont',color='g')

    if ax.is_first_row(): ax.set_title(reglabs[rii])
    pii+=1

gs2 = gridspec.GridSpec(1, 3)
gs2.update(top=0.44,left=0.08,right=0.92,wspace=0.28)
for fii,fkey in enumerate(fields):

    print 'fkey: ' + fkey
    ax=plt.subplot(gs2[0,fii])
    allregmdf=allflddt[fkey]
    allregmcidf=allcidt[fkey]
    xincr=0.2
    stagger=-xincr

    for eii,ens in enumerate(enss):
        print '  ens: ' + ens

        for rii,reg in enumerate(regions):
            print '    reg: ' + reg

            mult=multfacs[rii]
            xpos=xx[rii]+stagger
            ax.plot(xpos,mult*allregmdf[ens][reg],marker=mkrs[eii],
                    color=clrs[eii],mec=clrs[eii],linestyle='none',markersize=10)
            ci=allregmcidf[ens][reg]
            #print mult, (xpos,xpos), ci
            ax.plot((xpos,xpos),mult*np.array(ci),marker='_',
                    mew=2,markersize=10,linewidth=2,color=clrs[eii])
            ax.axvspan(xx[rii]-xincr-0.1,xx[rii]+xincr+0.1,color='0.8',alpha=0.5)

        stagger+=0.2

    ax.set_ylabel('Change in ' + reglabs[fii] + ' (' + reglabunits[fii] + ')')
    ax.axhline(y=0,color='k',linestyle='--')
    ax.set_xlim((stxx-0.7,xx[-1]+0.7))
    ax.set_ylim(fylims[fii])
    #if fii==0:
    if ax.is_first_col():
        ax.legend(lgs,lgstrs,loc='lower right',fancybox=True,frameon=False)#,prop=fontP)
    ax.set_xticks(xx)
    ax.set_xticklabels(reglabs,rotation=25)
    ax.annotate(plabs[pii],xy=(-0.02,1.04),
                xycoords='axes fraction',fontsize=16,fontweight='bold')
    pii+=1

cbar_ax,cbar = cplt.add_colorbar(fig,pc,orientation='horizontal',
                                 pos=[.25,.53, .5,.03],
                                 label='Change in surface air temperature, $\Delta$SAT (' + r2units +')')
#cbar.add_lines(pc)
if addsig:
    prstr='sig'
else:
    prstr=''
if printtofile:
    fig.savefig('Figure4_draft_6panel_CGCMmaps_compavgs' + prstr + 'swap2.pdf')
    fig.savefig('Figure4_draft_6panel_CGCMmaps_compavgs' + prstr + 'swap2.png',dpi=400)



print '---- plotting v4 of summary ----'
xx=xx[:-1]

fig=plt.figure(figsize=(12,9))
gs1 = gridspec.GridSpec(1, 3)
gs1.update(top=0.97,bottom=0.55,left=0.06,right=0.94,wspace=0.07)

cii=0; ckey='CGCM'; pii=0
for rii,rkey in enumerate(regions):
    mult=multfacs[rii]
    print '  ' + rkey + ' (mult=' + str(mult) + ')'
    ax=plt.subplot(gs1[cii,rii])
    dt=allcasedt[ckey][rkey]

    bm,pc=cplt.kemmap(mult*(dt[spkey]['lowsp']-dt[spkey]['highsp']),alat,alon,
                      ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='',suppcb=True,
                      panellab=plabs[pii],lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,dt[spkey]['sppval'], alat, alon)
    bm.contour(alons,alats,mult*(dt[sp2key]['lowsp']-dt[sp2key]['highsp']),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,dt[sp2key]['sppval'],alat,alon,sigtype='cont',color='g')

    if ax.is_first_row(): ax.set_title(reglabs[rii])
    pii+=1

gs2 = gridspec.GridSpec(1, 3)
gs2.update(top=0.44,left=0.08,right=0.92,wspace=0.28)
for fii,fkey in enumerate(fields):

    print 'fkey: ' + fkey
    rlabbools = np.ones(len(reglabs),np.bool)
    rlabbools[fii]=0
                
    ax=plt.subplot(gs2[0,fii])
    allregmdf=allflddt[fkey]
    allregmcidf=allcidt[fkey]
    xincr=0.2
    stagger=-xincr

    for eii,ens in enumerate(enss):
        print '  ens: ' + ens
        xii=0
        for rii,reg in enumerate(regions):
            print '    reg: ' + reg
            if (fkey in reg) or (fkey=='ice' and reg==r1key) or (fkey=='z500' and reg==r3key):
                print '    --SKIP'
            elif (ens=='AGCM_fixed' and fkey=='ice') or (ens=='AGCM_fixed' and reg=='bkssic'):
                print '    --SKIP (AGCM_fixed ice)'
                xii+=1
            else:
                mult=multfacs[rii]
                xpos=xx[xii]+stagger
                ax.plot(xpos,mult*allregmdf[ens][reg],marker=mkrs[eii],
                        color=clrs[eii],mec=clrs[eii],linestyle='none',markersize=10)
                ci=allregmcidf[ens][reg]
                #print mult, (xpos,xpos), ci
                ax.plot((xpos,xpos),mult*np.array(ci),marker='_',
                        mew=2,markersize=10,linewidth=2,color=clrs[eii])
                if len(enss)==3:
                    ax.axvspan(xx[xii]-xincr-0.1,xx[xii]+xincr+0.1,color='0.8',alpha=0.5)
                elif len(enss)==2:
                    ax.axvspan(xx[xii]-xincr-0.1,xx[xii]+0.1,color='0.8',alpha=0.5)
                xii+=1
                
        stagger+=0.2

    ax.set_ylabel('Change in ' + reglabs[fii] + ' (' + reglabunits[fii] + ')')
    ax.axhline(y=0,color='k',linestyle='--')
    if len(enss)==3:
        ax.set_xlim((stxx-0.7,xx[-1]+0.7))
        ax.set_xticks(xx)
    elif len(enss)==2:
        ax.set_xlim((stxx-0.7,xx[-1]+0.5))
        ax.set_xticks(xx-0.1)
    if fii==0:
        ax.legend(lgs,lgstrs,loc='lower left',fancybox=True,frameon=False)#,prop=fontP)
    
    ax.set_xticklabels(np.array(reglabs)[rlabbools],rotation=25)
    ax.annotate(plabs[pii],xy=(-0.02,1.04),
                xycoords='axes fraction',fontsize=16,fontweight='bold')
    pii+=1

cbar_ax,cbar = cplt.add_colorbar(fig,pc,orientation='horizontal',
                                 pos=[.25,.53, .5,.03],
                                 label='Change in surface air temperature, $\Delta$SAT (' + r2units +')')
#cbar.add_lines(pc)
if addsig:
    prstr='sig'
else:
    prstr=''
if printtofile:
    fig.savefig('Figure4_draft_6panel_CGCMmaps_compavgs' + prstr + 'swapv4.pdf')
    fig.savefig('Figure4_draft_6panel_CGCMmaps_compavgs' + prstr + 'swapv4.png',dpi=400)






# ======== plot all composite maps:
#allcasedt = {'Preindustrial':piallrdt, 'CGCM': leallrdt, 'AGCM':aallrdt}

# loop case
#    loop regions

pii=0

fig=plt.figure(figsize=(12,9))
gs1 = gridspec.GridSpec(3,3)
gs1.update(left=0.08,right=0.92,wspace=0.02,hspace=0.05)

for cii,ckey in enumerate(allcasedt.keys()):
    print ckey
    for rii,rkey in enumerate(regions):
        mult=multfacs[rii]
        print '  ' + rkey + ' (mult=' + str(mult) + ')'
        ax=plt.subplot(gs1[cii,rii])
        dt=allcasedt[ckey][rkey]

        bm,pc=cplt.kemmap(mult*(dt[spkey]['lowsp']-dt[spkey]['highsp']),alat,alon,
                          ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='',suppcb=True,
                          panellab=plabs[pii],lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,dt[spkey]['sppval'], alat, alon)
        bm.contour(alons,alats,mult*(dt[sp2key]['lowsp']-dt[sp2key]['highsp']),levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,dt[sp2key]['sppval'],alat,alon,sigtype='cont',color='g')

        if ax.is_first_row(): ax.set_title(reglabs[rii])
        if ax.is_first_col(): ax.set_ylabel(ckey)
        pii+=1

if addsig:
    prstr='sig'
else:
    prstr=''
if printtofile:
    fig.savefig('componall_allens_maps' + prstr +'.pdf')
    fig.savefig('componall_allens_maps' + prstr +'.png',dpi=400)




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
ylabs=('Preindustrial','CGCM','AGCM')
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
    iplabs=allplabs[axidx]

    ax = axs[axidx,0]
    bm,pc=cplt.kemmap(dt['highsp']-dt['meansp'],alat,alon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab=iplabs[0],lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,dt['hisppval'],alat,alon)
    ax.set_ylabel(ylabs[axidx])

    ax = axs[axidx,1]
    bm,pc=cplt.kemmap(dt['lowsp']-dt['meansp'],alat,alon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab=iplabs[1],lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,dt['losppval'],alat,alon)

    ax = axs[axidx,2]
    bm,pc=cplt.kemmap(diffmult*(dt['lowsp']-dt['highsp']),alat,alon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab=iplabs[2],lcol='0.2')
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
    bm,pc=cplt.kemmap(pihighice-pimeanice,lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='a',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pihiicepval,lat,lon)
    ax.set_ylabel('PI')
    ax=axs[0,1]
    bm,pc=cplt.kemmap(pilowice-pimeanice,lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='b',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,piloicepval,lat,lon)
    ax=axs[0,2]
    bm,pc=cplt.kemmap(diffmult*(pilowice-pihighice),lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='c',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,piicepval,lat,lon)
    ax=axs[1,0]
    bm,pc=cplt.kemmap(highleice-meanleice,lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehiicepval,lat,lon)
    ax.set_ylabel('CGCM')
    ax=axs[1,1]
    bm,pc=cplt.kemmap(lowleice-meanleice,lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,leloicepval,lat,lon)
    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(lowleice-highleice),lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm, leicepval,lat,lon)
    ax=axs[2,0]
    bm,pc=cplt.kemmap(ahighice-ameanice,alat,alon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='g',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahiicepval,alat,alon)
    ax.set_ylabel('AGCM')
    ax=axs[2,1]
    bm,pc=cplt.kemmap(alowice-ameanice,alat,alon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='h',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aloicepval,alat,alon)
    ax=axs[2,2]
    bm,pc=cplt.kemmap(diffmult*(alowice-ahighice),alat,alon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
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
    bm,pc=cplt.kemmap(pihighicer2-pimeanice,lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='a',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pihiicer2pval,lat,lon)
    ax.set_ylabel('PI')
    ax=axs[0,1]
    bm,pc=cplt.kemmap(pilowicer2-pimeanice,lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='b',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,piloicer2pval,lat,lon)
    ax=axs[0,2]
    bm,pc=cplt.kemmap(diffmult*(pilowicer2-pihighicer2),lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='c',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,piicer2pval,lat,lon)
    ax=axs[1,0]
    bm,pc=cplt.kemmap(highleicer2-meanleice,lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehiicer2pval,lat,lon)
    ax.set_ylabel('CGCM')
    ax=axs[1,1]
    bm,pc=cplt.kemmap(lowleicer2-meanleice,lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,leloicer2pval,lat,lon)
    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(lowleicer2-highleicer2),lat,lon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title=diffttl,suppcb=True,cmap=cmapice,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,leicer2pval,lat,lon)
    ax=axs[2,0]
    bm,pc=cplt.kemmap(ahighicer2-ameanice,alat,alon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='High anom',suppcb=True,cmap=cmapice,
                      panellab='g',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahiicer2pval,alat,alon)
    ax.set_ylabel('AGCM')
    ax=axs[2,1]
    bm,pc=cplt.kemmap(alowicer2-ameanice,alat,alon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
                      title='Low anom',suppcb=True,cmap=cmapice,
                      panellab='h',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aloicer2pval,alat,alon)
    ax=axs[2,2]
    bm,pc=cplt.kemmap(diffmult*(alowicer2-ahighicer2),alat,alon,ptype='nheur',axis=ax,cmin=cminice,cmax=cmaxice,
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
    bm,pc=cplt.kemmap(highlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl1,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(lons,lats,highlesp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('CGCM')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(lowlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl2,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(lons,lats,lowlesp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(lons,lats,meanlesp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    #cplt.add_colorbar(fig,pc,orientation='horizontal',pos=[.25,.5, .5,.03])

    ax=axs[1,0]
    bm,pc=cplt.kemmap(highlesp-meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehisppval,lat,lon)
    bm.contour(lons,lats,highlesp2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lehisp2pval,lat,lon,sigtype='cont',color='g')
    ax.set_ylabel('CGCM Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(lowlesp-meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lelosppval,lat,lon)
    bm.contour(lons,lats,lowlesp2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lelosp2pval,lat,lon,sigtype='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(lowlesp-highlesp),lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lesppval, lat, lon)
    bm.contour(lons,lats,diffmult*(lowlesp2-highlesp2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lesp2pval,lat,lon,sigtype='cont',color='g')

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
    bm,pc=cplt.kemmap(pihighsp,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl1,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(lons,lats,pihighsp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('PreIndustrial')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(pilowsp,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl2,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(lons,lats,pilowsp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(pimeansp,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(lons,lats,pimeansp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    #cplt.add_colorbar(fig,pc,orientation='horizontal',pos=[.25,.5, .5,.03])

    ax=axs[1,0]
    bm,pc=cplt.kemmap(pihighsp-pimeansp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pihisppval,lat,lon)
    bm.contour(lons,lats,pihighsp2-pimeansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pihisp2pval,lat,lon,sigtype='cont',color='g')
    ax.set_ylabel('PreI Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(pilowsp-pimeansp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pilosppval,lat,lon)
    bm.contour(lons,lats,pilowsp2-pimeansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pilosp2pval,lat,lon,sigtype='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(pilowsp-pihighsp),lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pisppval, lat, lon)
    bm.contour(lons,lats,diffmult*(pilowsp2-pihighsp2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pisp2pval,lat,lon,sigtype='cont',color='g')

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
    bm,pc=cplt.kemmap(ahighsp,alat,alon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl1,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(alons,alats,ahighsp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('AGCM')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(alowsp,alat,alon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl2,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(alons,alats,alowsp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(alons,alats,ameansp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    #cplt.add_colorbar(fig,pc,orientation='horizontal',pos=[.25,.5, .5,.03])

    ax=axs[1,0]
    bm,pc=cplt.kemmap(ahighsp-ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahisppval,alat,alon)
    bm.contour(alons,alats,ahighsp2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,ahisp2pval,alat,alon,sigtype='cont',color='g')

    ax.set_ylabel('AGCM Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(alowsp-ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,alosppval,alat,alon)
    bm.contour(alons,alats,alowsp2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,alosp2pval,alat,alon,sigtype='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(alowsp-ahighsp),alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,asppval, alat, alon)
    bm.contour(alons,alats,diffmult*(alowsp2-ahighsp2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,asp2pval,alat,alon,sigtype='cont',color='g')

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
    bm,pc=cplt.kemmap(highlesp-meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='a',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehisppval,lat,lon)
    bm.contour(lons,lats,highlesp2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lehisp2pval,lat,lon,sigtype='cont',color='g')
    ax.set_ylabel('CGCM Diff')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(lowlesp-meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='b',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lelosppval,lat,lon)
    bm.contour(lons,lats,lowlesp2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lelosp2pval,lat,lon,sigtype='cont',color='g')

    ax=axs[0,2]
    bm,pc=cplt.kemmap(diffmult*(lowlesp-highlesp),lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='c',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lesppval, lat, lon)
    bm.contour(lons,lats,lowlesp2-highlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lesp2pval,lat,lon,sigtype='cont',color='g')

    ax=axs[1,0]
    bm,pc=cplt.kemmap(ahighsp-ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahisppval,alat,alon)
    bm.contour(alons,alats,ahighsp2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,ahisp2pval,alat,alon,sigtype='cont',color='g')
    ax.set_ylabel('AGCM Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(alowsp-ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,alosppval,alat,alon)
    bm.contour(alons,alats,alowsp2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,alosp2pval,alat,alon,sigtype='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(alowsp-ahighsp),alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,asppval, alat, alon)
    bm.contour(alons,alats,diffmult*(alowsp2-ahighsp2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,asp2pval,alat,alon,sigtype='cont',color='g')

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
    bm,pc=cplt.kemmap(highlespr2,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl12,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(lons,lats,highlesp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('CGCM')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(lowlespr2,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl22,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(lons,lats,lowlesp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(lons,lats,meanlesp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[1,0]
    bm,pc=cplt.kemmap(highlespr2-meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehispr2pval,lat,lon)
    bm.contour(lons,lats,highlesp2r2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lehisp2r2pval,lat,lon,sigtype='cont',color='g')
    ax.set_ylabel('CGCM Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(lowlespr2-meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lelospr2pval,lat,lon)
    bm.contour(lons,lats,lowlesp2r2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lelosp2r2pval,lat,lon,sigtype='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(lowlespr2-highlespr2),lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lespr2pval, lat, lon)
    bm.contour(lons,lats,diffmult*(lowlesp2r2-highlesp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lesp2r2pval,lat,lon,sigtype='cont',color='g')

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
    bm,pc=cplt.kemmap(pihighspr2,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl12,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(lons,lats,pihighsp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('PreIndustrial')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(pilowspr2,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl22,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(lons,lats,pilowsp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(pimeansp,lat,lon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(lons,lats,pimeansp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    #cplt.add_colorbar(fig,pc,orientation='horizontal',pos=[.25,.5, .5,.03])

    ax=axs[1,0]
    bm,pc=cplt.kemmap(pihighspr2-pimeansp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pihispr2pval,lat,lon)
    bm.contour(lons,lats,pihighsp2r2-pimeansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pihisp2r2pval,lat,lon,sigtype='cont',color='g')
    ax.set_ylabel('PreI Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(pilowspr2-pimeansp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pilospr2pval,lat,lon)
    bm.contour(lons,lats,pilowsp2r2-pimeansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pilosp2r2pval,lat,lon,sigtype='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(pilowspr2-pihighspr2),lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,pispr2pval, lat, lon)
    bm.contour(lons,lats,diffmult*(pilowsp2r2-pihighsp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,pisp2r2pval,lat,lon,sigtype='cont',color='g')

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
    bm,pc=cplt.kemmap(ahighspr2,alat,alon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl12,suppcb=True,
                      panellab='a',lcol='0.2')
    bm.contour(alons,alats,ahighsp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)
    ax.set_ylabel('AGCM')

    ax=axs[0,1]
    bm,pc=cplt.kemmap(alowspr2,alat,alon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title=ttl22,suppcb=True,
                      panellab='b',lcol='0.2')
    bm.contour(alons,alats,alowsp2r2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[0,2]
    bm,pc=cplt.kemmap(ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                      title='Mean Anom',suppcb=True,
                      panellab='c',lcol='0.2')
    bm.contour(alons,alats,ameansp2,levels=contsbig,
               colors='0.5',linewidths=1,latlon=True)

    ax=axs[1,0]
    bm,pc=cplt.kemmap(ahighspr2-ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahispr2pval,alat,alon)
    bm.contour(alons,alats,ahighsp2r2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,ahisp2r2pval,alat,alon,sigtype='cont',color='g')
    ax.set_ylabel('AGCM Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(alowspr2-ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,alospr2pval,alat,alon)
    bm.contour(alons,alats,alowsp2r2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,alosp2r2pval,alat,alon,sigtype='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(alowspr2-ahighspr2),alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aspr2pval, alat, alon)
    bm.contour(alons,alats,diffmult*(alowsp2r2-ahighsp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,asp2r2pval,alat,alon,sigtype='cont',color='g')

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
        bm,pc=cplt.kemmap(aehighspr2,alat,alon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                          title=ttl12,suppcb=True,
                          panellab='a',lcol='0.2')
        bm.contour(alons,alats,aehighsp2r2,levels=contsbig,
                   colors='0.5',linewidths=1,latlon=True)
        ax.set_ylabel('AGCM')

        ax=axs[0,1]
        bm,pc=cplt.kemmap(aelowspr2,alat,alon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                          title=ttl22,suppcb=True,
                          panellab='b',lcol='0.2')
        bm.contour(alons,alats,aelowsp2r2,levels=contsbig,
                   colors='0.5',linewidths=1,latlon=True)

        ax=axs[0,2]
        bm,pc=cplt.kemmap(aemeansp,alat,alon,ptype='nheur',axis=ax,cmin=cminspbig,cmax=cmaxspbig,
                          title='Mean Anom',suppcb=True,
                          panellab='c',lcol='0.2')
        bm.contour(alons,alats,aemeansp2,levels=contsbig,
                   colors='0.5',linewidths=1,latlon=True)

        ax=axs[1,0]
        bm,pc=cplt.kemmap(aehighspr2-aemeansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='High Anom',suppcb=True,
                          panellab='d',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aehispr2pval,alat,alon)
        bm.contour(alons,alats,aehighsp2r2-aemeansp2,levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,aehisp2r2pval,alat,alon,sigtype='cont',color='g')
        ax.set_ylabel('AGCM Diff')

        ax=axs[1,1]
        bm,pc=cplt.kemmap(aelowspr2-aemeansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='Low Anom',suppcb=True,
                          panellab='e',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aelospr2pval,alat,alon)
        bm.contour(alons,alats,aelowsp2r2-aemeansp2,levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,aelosp2r2pval,alat,alon,sigtype='cont',color='g')

        ax=axs[1,2]
        bm,pc=cplt.kemmap(diffmult*(aelowspr2-aehighspr2),alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title=diffttl,suppcb=True,
                          panellab='f',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aespr2pval, alat, alon)
        bm.contour(alons,alats,diffmult*(aelowsp2r2-aehighsp2r2),levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,aesp2r2pval,alat,alon,sigtype='cont',color='g')

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
        bm,pc=cplt.kemmap(diffmult*(alowspr2-ahighspr2),alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='Var (' + diffttl + ')',suppcb=True,
                          panellab='a',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aspr2pval, alat, alon)
        bm.contour(alons,alats,diffmult*(alowsp2r2-ahighsp2r2),levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,asp2r2pval,alat,alon,sigtype='cont',color='g')

        ax=axs[1]
        bm,pc=cplt.kemmap(diffmult*(aelowspr2-aehighspr2),alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='Const (' + diffttl + ')',suppcb=True,
                          panellab='b',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aespr2pval, alat, alon)
        bm.contour(alons,alats,diffmult*(aelowsp2r2-aehighsp2r2),levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,aesp2r2pval,alat,alon,sigtype='cont',color='g')

        ax=axs[2] 
        bm,pc=cplt.kemmap((alowspr2-ahighspr2)-(aelowspr2-aehighspr2),alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                          title='Var-Const (' + diffttl + ')',suppcb=True,
                          panellab='c',lcol='0.2')
        if addsig:
            cplt.addtsigm(bm,aressppval,alat,alon)
        bm.contour(alons,alats,(alowsp2r2-ahighsp2r2)-(aelowsp2r2-aehighsp2r2),levels=conts,
                   colors='0.5',linewidths=1,latlon=True)
        if addsig:
            cplt.addtsigm(bm,aressp2pval,alat,alon,sigtype='cont',color='g')

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
        cplt.addtsigm(bm,pisp2pval,lat,lon,sigtype='cont',color='g')
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
        cplt.addtsigm(bm,lesp2pval,lat,lon,sigtype='cont',color='g')
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
    bm,pc=cplt.kemmap(highlespr2-meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='a',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lehispr2pval,lat,lon)
    bm.contour(lons,lats,highlesp2r2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lehisp2r2pval,lat,lon,sigtype='cont',color='g')

    ax.set_ylabel('CGCM Diff')
    ax=axs[0,1]
    bm,pc=cplt.kemmap(lowlespr2-meanlesp,lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='b',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lelospr2pval,lat,lon)
    bm.contour(lons,lats,lowlesp2r2-meanlesp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lelosp2r2pval,lat,lon,sigtype='cont',color='g')

    ax=axs[0,2]
    bm,pc=cplt.kemmap(diffmult*(lowlespr2-highlespr2),lat,lon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='c',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,lespr2pval, lat, lon)
    bm.contour(lons,lats,diffmult*(lowlesp2r2-highlesp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,lesp2r2pval,lat,lon,sigtype='cont',color='g')

    ax=axs[1,0]
    bm,pc=cplt.kemmap(ahighspr2-ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='High Anom',suppcb=True,
                      panellab='d',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,ahispr2pval,alat,alon)
    bm.contour(alons,alats,ahighsp2r2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,ahisp2r2pval,alat,alon,sigtype='cont',color='g')
    ax.set_ylabel('AGCM Diff')

    ax=axs[1,1]
    bm,pc=cplt.kemmap(alowspr2-ameansp,alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title='Low Anom',suppcb=True,
                      panellab='e',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,alospr2pval,alat,alon)
    bm.contour(alons,alats,alowsp2r2-ameansp2,levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,alosp2r2pval,alat,alon,sigtype='cont',color='g')

    ax=axs[1,2]
    bm,pc=cplt.kemmap(diffmult*(alowspr2-ahighspr2),alat,alon,ptype='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                      title=diffttl,suppcb=True,
                      panellab='f',lcol='0.2')
    if addsig:
        cplt.addtsigm(bm,aspr2pval, alat, alon)
    bm.contour(alons,alats,diffmult*(alowsp2r2-ahighsp2r2),levels=conts,
               colors='0.5',linewidths=1,latlon=True)
    if addsig:
        cplt.addtsigm(bm,asp2r2pval,alat,alon,sigtype='cont',color='g')

    cplt.add_colorbar(fig,pc,orientation='horizontal')

    if printtofile:
        if addsig:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr2+regionr2+ sear + '_CGCMAGCM' + prstr + '.png',dpi=400) 
        else:
            fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                        '_compon_' + fieldr2+regionr2+ sear + '_CGCMAGCM' + prstr + '.pdf') 
