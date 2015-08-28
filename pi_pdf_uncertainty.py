import canesm_LE_general as leg
import pandas as pd
import loadCanESM2data as lcd
import cccmautils as cutl
import loadmodeldata as lmd
import loadLE as le
import cccmacmaps as ccm
import matplotlib.lines as mlines
import collections as coll
import constants as con


printtofile=False
local=True
nresamp=100 # number of times to resample PI (50-members each time)

fieldr='tas'; ncfieldr='tas'; compr='Amon'; regionr='eurasiamori'; 
rstr='Eur SAT'; rstrlong='Eurasian SAT'; runits='$^\circ$C'; rkey='eursat'
convr=1
aconvr=1 # for agcm sims
ylev=0.3 # for mean dots

fieldr='zg50000.00'; ncfieldr='zg'; compr='Amon'; 
regionr='gt60n'; rstr='POL Z500'; rstrlong='Polar Z500'; runits='m'; rkey='gt60nz500'
#regionr='bksmori'; rstr='BKS Z500'; rstrlong='Barents/Kara Z500'; runits='m'; rkey='bksz500'
convr=1
aconvr=1/con.get_g() # for agcm sims
ylev=0.005 # for mean dots


sear='DJF'
timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'


picol = '0.8'
acol=ccm.get_linecolor('paperblue')

when='14:51:28.762886'; styearsR = [ 8.,  7.,  2.,  8.,  8.] # variable SIC styears
styearsE=[ 4.,  1.,  7.,  3.,  1.]; styearsN=[1.] #when for these: 17:01:16.908687


def load_pifield(fdict,seas,conv=1,subsampyrs=11,numsamp=50,
                 styear=None,anomyears=None,local=False,detrend=True,
                 verb=False,addcyc=False):
    """ TAKEN FROM load_canesmfield() in canesm_LE_composite.py

        loads subsampled CGCM data from specified simulation 
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
    casename='piControl'

    pidat = lcd.load_data(fdict,casename,local=local,conv=conv,verb=verb)
    piseadat = cutl.seasonalize_monthlyts(pidat,season=seas)

    if detrend:
        piseadat = cutl.detrend(piseadat,axis=0)

    # the data must be seasonalized before using this func.
    pisea,styear,anomyears = leg.subsamp_anom_pi(piseadat, numyrs=subsampyrs,numsamp=numsamp,
                                                 styear=styear,anomyears=anomyears)
    if addcyc:
        st=pisea[...,-1]
        pisea=np.dstack((pisea,st[...,None]))

    return pisea,styear,anomyears


def load_agcmfield(field,sims,seas,conv=1,region=None,subsampyrs=11,styears=None):
    """  TAKEN FROM load_agcmfield() in canesm_LE_composite.py

        loads subsampled agcm data from specified simulations

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


def load_field(fdict,casename,timesel,seas,ftype='fullts',conv=1,local=False,verb=False):
    
    """ TAKEN FROM load_field() in canesm_le_composite.py

        returns [numens x space.flat] or [numens]

    """

    ledat = le.load_LEdata(fdict,casename,timesel=timesel, 
                           rettype='ndarray',conv=conv,ftype=ftype,local=local,verb=verb)

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





def calc_pdf(dat, cisiglevel=0.05):

    pipdf_fitted,pimean,pisd,pixx = cutl.calc_normfit(dat)
    pidof = len(dat)-1
    pistder = pisd / np.sqrt(pidof+1)
    pici = sp.stats.t.interval(1-cisiglevel, pidof, loc=pimean, scale=pistder)
    picif = sp.stats.t.interval(1-cisiglevel, pidof, loc=pimean, scale=pisd)

    return pipdf_fitted, pimean, pixx


def test_significance(dat1, dat2, siglevel=0.05,verb=False):

    tstat, pval = sp.stats.ttest_ind(dat1,dat2)
    lstat, lpval = sp.stats.levene(dat1,dat2)
    if verb:
        print '  Mean1: ' + str(dat1.mean()) + ', Mean2: ' + str(dat2.mean())
        print '  TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
        if pval<=siglevel:
             print '   The ensemble means are significantly different (' + str(1-siglevel) + ')'
        print '  LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
        if lpval<=siglevel:
             print '   The ensemble variances are significantly different (' + str(1-siglevel) + ')'

    return pval, lpval


# ########### Processing #############
fdictr = {'field': fieldr+regionr, 'ncfield': ncfieldr, 'comp': compr}

piseardt={}; styeardt={}; anomyearsdt={}
pipdfdt={}; pimeandt={}; pixxdt={}

for ii in range(0,nresamp):

    pisear,styear,anomyears = load_pifield(fdictr,sear,conv=convr,
                                           local=local,addcyc=False)
    piseardt[ii] = pisear
    styeardt[ii] = styear
    anomyearsdt[ii] = anomyears

    pipdf, pimean, pixx = calc_pdf(pisear)

    pipdfdt[ii] = pipdf
    pimeandt[ii] = pimean
    pixxdt[ii] = pixx

piseardf=pd.DataFrame(piseardt)
pipdfdf=pd.DataFrame(pipdfdt)
pixxdf=pd.DataFrame(pixxdt)
pimeandf=pd.Series(pimeandt)
piflatpdf, piflatmean, piflatxx = calc_pdf(piseardf.values.flatten())

# Load AGCM sims
simsR=('R1','R2','R3','R4','R5')
arsear,styearsr = load_agcmfield(fieldr,simsR,sear,region=regionr,
                                 conv=aconvr,styears=styearsR)
arpdf, armean, arxx = calc_pdf(arsear)

simsE=('E1','E2','E3','E4','E5')
aesear,styearse = load_agcmfield(fieldr,simsE,sear,region=regionr,
                                 conv=aconvr,styears=styearsE)
aepdf, aemean, aexx = calc_pdf(aesear)


# Load LE (historical)
lecsear = load_field(fdictr, 'historical', timeselc, sear,
                     ftype='fullts', conv=convr,local=local)
lepsear = load_field(fdictr, 'historical', timeselp, sear,
                     ftype='fullts', conv=convr,local=local)
lesear = lepsear-lecsear
lepdf, lemean, lexx = calc_pdf(lesear)

print '====R SIMS v E SIMS'
test_significance(arsear, aesear, verb=True)
print ''
print '====R SIMS v PI avg subsamp'
test_significance(arsear,piseardf.mean(axis=1), verb=True)
print '====E SIMS v PI avg subsamp'
test_significance(aesear,piseardf.mean(axis=1), verb=True)
print '====LE v PI avg subsamp'
test_significance(lesear,piseardf.mean(axis=1), verb=True)
print ''

print '====R SIMS v PI flatten (50*' + str(nresamp) + ')'
test_significance(arsear,piseardf.values.flatten(), verb=True)
print '====E SIMS v PI flatten (50*' + str(nresamp) + ')'
test_significance(aesear,piseardf.values.flatten(), verb=True)
print '====LE v PI flatten (50*' + str(nresamp) + ')'
test_significance(lesear,piseardf.values.flatten(), verb=True)





# ########### PLOTTING ###################
legdt=coll.OrderedDict()
legdt['Preindustrial ' + str(nresamp) + ' 50-member samples'] = mlines.Line2D([],[],
                                                                              color=picol)
legdt['Preindustrial mean of ' + str(nresamp)] =  mlines.Line2D([],[],
                                                                color='k',linewidth=2)
legdt['Preindustrial mean of ' + str(nresamp*50)] =  mlines.Line2D([],[],
                                                                   color='k',linewidth=2,
                                                                   linestyle='--')
legdt['AGCM_variable'] = mlines.Line2D([],[],color=acol,linewidth=2)
legdt['AGCM_fixed'] = mlines.Line2D([],[],color=acol,linewidth=2,
                                    linestyle='--')
legdt['CGCM']=mlines.Line2D([],[],color='0.7',linewidth=2)


fig,ax=plt.subplots(1,1)
fig.set_size_inches(9,8)
ax.plot(pixxdf.values, pipdfdf.values, color=picol,alpha=0.5)
ax.plot(pixxdf.mean(axis=1), pipdfdf.mean(axis=1), color='k', linewidth=2)
ax.plot(pimeandf.values, np.ones((pixxdf.values.shape[1]))*ylev,  
         linestyle='none',marker='.', color=picol)
ax.plot(pimeandf.mean(), ylev, linestyle='none',marker='.',color='k')
ax.plot(piflatxx, piflatpdf, color='k', linestyle='--',linewidth=2)

ax.plot(arxx, arpdf, color=acol,linewidth=3)
ax.plot(armean, ylev, linestyle='none',marker='s', color=acol,alpha=0.5)
ax.plot(aexx, aepdf, color=acol,linewidth=3,linestyle='--')
ax.plot(aemean, ylev, linestyle='none',marker='o', color=acol,alpha=0.5)

ax.plot(lexx, lepdf, color='0.7',linewidth=2)
ax.plot(lemean, ylev, linestyle='none', marker='^', color='0.7')

#ax.set_xlabel('Full PI sample mean: ' + str(pimeandf.mean()))
ax.set_ylabel('Density')
ax.set_xlabel('Change in ' + rstrlong + ' (' + runits + ')')

ax.legend((legdt.values()), (legdt.keys()),loc='upper left',
          fancybox=True,framealpha=0.5,frameon=False)
if printtofile:
    fig.savefig(fieldr + regionr + '_' + sear + '_' +\
                'pi' + str(nresamp) + 'resamp_AGCMer_LE' +'.pdf')



# SAVE SOME OUTPUT
"""
******* This is for resamp= 1000 and nsamp=50 ******
************ EUR SAT

1. this first output tests against the solid black curve of PI,
  the average of all the sampled pdfs (avg of 1000 pdfs)

R SIMS v PI avg subsamp
==== testing input 1 vs input 2 ====
TSTAT: 0.949425595214 PVAL: 0.344739700599
LSTAT: 78.1725745675 PVAL: 3.93980817057e-14
The ensemble variances are significantly different (0.95)

E SIMS v PI avg subsamp
==== testing input 1 vs input 2 ====
TSTAT: 0.0684858870774 PVAL: 0.945538396429
LSTAT: 77.2108084012 PVAL: 5.16896427751e-14
The ensemble variances are significantly different (0.95)

2. the second output tests against the dashed black curve of PI,
  which is the full distribution of all 50,000 anomalies

R SIMS v PI flatten (50*1000)
==== testing input 1 vs input 2 ====
TSTAT: 1.10295407114 PVAL: 0.270052414469
LSTAT: 1.86544267023 PVAL: 0.172004087696

E SIMS v PI flatten (50*1000)
==== testing input 1 vs input 2 ====
TSTAT: 0.0742550921524 PVAL: 0.940807706367
LSTAT: 0.31666140677 PVAL: 0.573623479656



******* This is for resamp= 10000 and nsamp=50 ******

R SIMS v E SIMS
==== testing input 1 vs input 2 ====
TSTAT: 0.647565524008 PVAL: 0.518780140966
LSTAT: 0.243251494042 PVAL: 0.622970911577
R SIMS v PI avg subsamp
==== testing input 1 vs input 2 ====
TSTAT: 0.959612321708 PVAL: 0.339611949009
LSTAT: 80.6482833253 PVAL: 1.9720870637e-14
The ensemble variances are significantly different (0.95)
E SIMS v PI avg subsamp
==== testing input 1 vs input 2 ====
TSTAT: 0.0791961039654 PVAL: 0.93703813417
LSTAT: 79.8613461701 PVAL: 2.45467442474e-14
The ensemble variances are significantly different (0.95)


R SIMS v PI flatten (50*10000)
==== testing input 1 vs input 2 ====
TSTAT: 1.11050316472 PVAL: 0.266782799126
LSTAT: 1.73038241799 PVAL: 0.188362669839
E SIMS v PI flatten (50*10000)
==== testing input 1 vs input 2 ====
TSTAT: 0.0855281292414 PVAL: 0.931841567739
LSTAT: 0.264296082404 PVAL: 0.607184159031



****************** BKS Z500 *****************
******* This is for resamp= 100 and nsamp=50 ******

====R SIMS v E SIMS
  Mean1: 3.34624519573, Mean2: 4.08476507306
  TSTAT: -0.203966108758 PVAL: 0.838802946149
  LSTAT: 1.30148791711 PVAL: 0.256721855393

====R SIMS v PI avg subsamp
  Mean1: 3.34624519573, Mean2: -0.362034874386
  TSTAT: 1.3039819856 PVAL: 0.195294313572
  LSTAT: 52.4262093868 PVAL: 1.01754909141e-10
   The ensemble variances are significantly different (0.95)
====E SIMS v PI avg subsamp
  Mean1: 4.08476507306, Mean2: -0.362034874386
  TSTAT: 1.95193582793 PVAL: 0.0538000016487
  LSTAT: 60.1696440145 PVAL: 8.34303753016e-12
   The ensemble variances are significantly different (0.95)
====LE v PI avg subsamp
  Mean1: 20.9342897461, Mean2: -0.362034874386
  TSTAT: 9.1311157635 PVAL: 9.30117844593e-15
   The ensemble means are significantly different (0.95)
  LSTAT: 58.1727027078 PVAL: 1.57054953536e-11
   The ensemble variances are significantly different (0.95)

====R SIMS v PI flatten (50*100)
  Mean1: 3.34624519573, Mean2: -0.362034874386
  TSTAT: 1.38785266275 PVAL: 0.165243167144
  LSTAT: 0.025539764451 PVAL: 0.873035838143
====E SIMS v PI flatten (50*100)
  Mean1: 4.08476507306, Mean2: -0.362034874386
  TSTAT: 1.66757010635 PVAL: 0.095463123568
  LSTAT: 2.19780283627 PVAL: 0.138270032351
====LE v PI flatten (50*100)
  Mean1: 20.9342897461, Mean2: -0.362034874386
  TSTAT: 7.98484450166 PVAL: 1.72862927032e-15
   The ensemble means are significantly different (0.95)
  LSTAT: 1.51092221294 PVAL: 0.219055711487





****************** Polar (>60N) Z500 *****************
******* This is for resamp= 100 and nsamp=50 ******

====R SIMS v E SIMS
  Mean1: 4.34329642586, Mean2: 4.05086208953
  TSTAT: 0.124864173526 PVAL: 0.900886856615
  LSTAT: 1.23105032791 PVAL: 0.269918102382

====R SIMS v PI avg subsamp
  Mean1: 4.34329642586, Mean2: -0.292373281517
  TSTAT: 2.60305597301 PVAL: 0.0106742279736
   The ensemble means are significantly different (0.95)
  LSTAT: 65.771539419 PVAL: 1.47728131318e-12
   The ensemble variances are significantly different (0.95)
====E SIMS v PI avg subsamp
  Mean1: 4.05086208953, Mean2: -0.292373281517
  TSTAT: 2.82554660165 PVAL: 0.00572033369574
   The ensemble means are significantly different (0.95)
  LSTAT: 52.9011902239 PVAL: 8.69480715577e-11
   The ensemble variances are significantly different (0.95)
====LE v PI avg subsamp
  Mean1: 22.6614946398, Mean2: -0.292373281517
  TSTAT: 17.1831659206 PVAL: 2.49304911358e-31
   The ensemble means are significantly different (0.95)
  LSTAT: 75.1006068889 PVAL: 9.42967223238e-14
   The ensemble variances are significantly different (0.95)

====R SIMS v PI flatten (50*100)
  Mean1: 4.34329642586, Mean2: -0.292373281517
  TSTAT: 2.9180817548 PVAL: 0.0035374693358
   The ensemble means are significantly different (0.95)
  LSTAT: 0.893780887261 PVAL: 0.344500024223
====E SIMS v PI flatten (50*100)
  Mean1: 4.05086208953, Mean2: -0.292373281517
  TSTAT: 2.73830354628 PVAL: 0.00619736600903
   The ensemble means are significantly different (0.95)
  LSTAT: 0.572568927471 PVAL: 0.44927520812
====LE v PI flatten (50*100)
  Mean1: 22.6614946398, Mean2: -0.292373281517
  TSTAT: 14.4881870482 PVAL: 1.22760725958e-46
   The ensemble means are significantly different (0.95)
  LSTAT: 1.69412618173 PVAL: 0.193117335082


"""
