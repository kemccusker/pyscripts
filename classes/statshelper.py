""" statshelper.py
       11/24/2014
       first function: pattern correlate within ensemble
"""
import numpy as np
import scipy as sp
import scipy.stats
import pandas as pd
import constants as con
import cccmautils as cutl
import cccmaNC as cnc

con = reload(con)
cutl = reload(cutl)
cnc = reload(cnc)

def pattcorr_withinensemble(ename,fdict,latlim=60,timesel='0002-01-01,0121-12-31'):
    """ pattcorr_withinensemble(ename,field,latlim=60)
            pattern corr each member of ensemble with each other one

            return pctable, pctablesea (DataFrames)
    """
    # @@ need diffdict

    field=fdict['field']
    ncfield=fdict['ncfield']
    conv=fdict['conv']
    
    seasons=('SON','DJF','MAM','JJA')
    
    
    if ename=='ANT':
        ename='histIC'
    elif ename=='TOT':
        ename='histBC'

    enssims = con.build_ensemblesims(ename)
    ensnum=len(enssims)

    print 'ENSSIMS: ' # @@@
    print enssims # @@
    
    # generate weights for the pattern corr
    lat = con.get_t63lat()
    lon = con.get_t63lon()

    areas = cutl.calc_cellareas(lat,lon)
    areas = areas[lat>latlim,:]
    weights = areas / np.sum(np.sum(areas,axis=1),axis=0)

    # ========= create diffdict first =====
    diffdict={}
    seadiffdict={}
    for skey in enssims:
        fnamec,fnamep = con.build_filepathpair(skey,field)

        # monthly calc
        fldc = cnc.getNCvar(fnamec,ncfield,timesel=timesel)*conv
        fldp = cnc.getNCvar(fnamep,ncfield,timesel=timesel)*conv
        fldd = fldp-fldc

        # Monthly
        flddclimo,flddstd = cutl.climatologize(fldd) # climo first (don't need to do for BCs technically)
        #flddcclimo,flddcstd = cutl.climatologize(flddc) # climo first. baseline diff data

        diffdict[skey] = flddclimo
        print skey + ' ' + str(flddclimo.shape) # @@@

        # Seasonal
        flddsea = np.zeros((4,len(lat),len(lon)))

        for seaii,sea in enumerate(seasons):
            fldcsea = np.mean(cnc.getNCvar(fnamec,ncfield,timesel=timesel,seas=sea)*conv,axis=0)
            fldpsea = np.mean(cnc.getNCvar(fnamep,ncfield,timesel=timesel,seas=sea)*conv,axis=0)
            flddsea[seaii,...] = fldpsea-fldcsea

        seadiffdict[skey] = flddsea

    # ======= Now do pattern corrs within ensemble ====

    # ======= copied from Canam4_BCpatterncorr-Copy0.py ===========

    outterdict= dict.fromkeys(enssims)

    for skey1 in enssims:

        outfld = diffdict[skey1]

        innerdict = dict.fromkeys(enssims)

        for skey2 in enssims:
            #print skey1 + ' compared to ' + skey2

            infld = diffdict[skey2]

            # for each month, compute pattern corr
            pc = np.zeros((12))
            for mii,mon in enumerate(con.get_mon()):
                tmp = np.squeeze(infld[mii,lat>latlim,...])
                tmpcmp = np.squeeze(outfld[mii,lat>latlim,...])
                pc[mii] = cutl.pattcorr(tmp.flatten()*weights.flatten(),
                                        tmpcmp.flatten()*weights.flatten())

            innerdict[skey2] = pc

        outterdict[skey1] = innerdict

    pctable = pd.DataFrame(outterdict) # 5x5

    # seasonal 
    outterdictsea= dict.fromkeys(enssims)

    for skey1 in enssims:

        outfld = seadiffdict[skey1]

        innerdictsea = dict.fromkeys(enssims)

        for skey2 in enssims:    

            #print skey1 + ' compared to ' + skey2

            infld = seadiffdict[skey2]

            # for each season, compute pattern corr
            pcsea = np.zeros((4))
            for seaii,sea in enumerate(seasons):
                tmp = np.squeeze(infld[seaii,lat>latlim,...])
                tmpcmp = np.squeeze(outfld[seaii,lat>latlim,...])
                pcsea[seaii] = cutl.pattcorr(tmp.flatten()*weights.flatten(),
                                        tmpcmp.flatten()*weights.flatten())

            innerdictsea[skey2] = pcsea

        outterdictsea[skey1] = innerdictsea

    pctablesea = pd.DataFrame(outterdictsea) # 5x5

    return pctable, pctablesea

    # ======= end copy Canam4_BCpatterncorr-Copy0.py ===========


def calc_ensemblestats(datablob,ensnames, seas=None,siglevel=0.05):
    """         
           datablob: 'diff' should be scalar regional mean anomalies
           ensnames: a tuple of ensemble names in the datablob.
                     must match ensname in sims dictionary. e.g. histIC, histBC
           seas: 2 or 3 month seasons

           @@ add return vals
    """
    
    diffdt = datablob['diff']
    sims=diffdt.keys()
    simspdt = con.get_simpairsdict()

    if seas==None:
        seasons=('SON','DJF','MAM','JJA')
    else:
        seasons=seas
        
    # Need all members of ensemble, plus mean to do ttest between
    # ensemble means. This ttest is NOT in time, but across
    # ensemble members. Thus, n=5 for TOT and ANT.

    # return tstat, pval, stddev # all across ensemble, not in time.

   
    allensdt={}; allensmdt={};
            
    for ensname in ensnames:
        ensdt={}; ensmdt={}
        print ensname
        for skey in sims: # for each simulation check if it's in the ensemble

            if simspdt[skey]['pert']['ensname']==ensname:
                # create an ensemble dict
                ensdt[skey] = datablob['diff'][skey]
            if simspdt[skey]['pert']['ensname']==ensname+'mean':
                ensmdt[skey] = datablob['diff'][skey]

        allensdt[ensname] = ensdt # dict of ens -> dict of sims in ens --> data
        allensmdt[ensname] = ensmdt # just the ensemble mean
        
    # end loop through ens
    ensdf1 = pd.DataFrame(allensdt[ensnames[0]])
    ensdf2 = pd.DataFrame(allensdt[ensnames[1]])
    ensm1 = allensmdt[ensnames[0]] # need DataFrame for this? should just be a val for each season
    ensm2 = allensmdt[ensnames[1]]
                    
    for sea in seasons:
         print sea
         
         e1 = ensdf1.loc[sea]
         e2 = ensdf2.loc[sea]
         print sea
         print ensnames[0] + ' MEAN: ' + str(e1.mean()) + ' STD: ' + str(e1.std())
         print ensnames[1] + ' MEAN: ' + str(e2.mean()) + ' STD: ' + str(e2.std())

         tstat, pval = sp.stats.ttest_ind(e1,e2)
         print 'TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
         if pval<=siglevel:
             print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'

         fstat, fpval = sp.stats.f_oneway(e1,e2)
         print 'FSTAT: ' + str(fstat) + ' PVAL: ' + str(fpval)
         if fpval<=siglevel:
             print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'

         lstat, lpval = sp.stats.levene(e1,e2)
         print 'LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
         if lpval<=siglevel:
             print 'The ensemble variances are significantly different (' + str(1-siglevel) + ')'


    print '@@@@ not done, still need to add return vals'



def calc_runstats(datablob,sims, seas=None,siglevel=0.05):
    """         
           datablob: 'ctl','pert','diff' should be timeseries of regional means
           sims: a tuple of sim names in the datablob. (prob unnecessary)
           seas: 2 or 3 month seasons

           @@ add return vals

    """
    ctldt = datablob['ctl']
    pertdt = datablob['pert']
    
    diffdt = datablob['diff']
    simspdt = con.get_simpairsdict()

    if seas==None:
        seasons=('SON','DJF','MAM','JJA')
    else:
        seasons=seas

    allstats={}; alltpval={}; allfpval={}
    
    for sim in sims:
        ctlsim=ctldt[sim]
        pertsim=pertdt[sim]
        seatpval={}; seafpval={}
        
        for sea in seasons:
            thestats={}
            print sim + ' ' + sea
            ctl=ctlsim[sea]
            pert=pertsim[sea]

            ## thestats['ctlmean']=ctl.mean()
            ## thestats['pertmean']=pert.mean()
            ## thestats['meandiff'] = pert.mean()-ctl.mean()
            ## thestats['ctlstd']=ctl.std()
            ## thestats['pertstd']=pert.std()
            ## thestats['stddiff']=pert.std()-ctl.std()
            
            print '    DIFF: ' + str(pert.mean()-ctl.mean()) + ', CTL mean: ' + str(ctl.mean()) + ', PERT mean: ' + str(pert.mean())
            tstat, tpval = sp.stats.ttest_ind(pert,ctl) # add autocorr @@
            print '    TSTAT: ' + str(tstat) + ' PVAL: ' + str(tpval)
            if tpval<=siglevel:
                print '  **The ensemble means are significantly different (' + str(1-siglevel) + ')'
            #thestats['tpval'] = pval
            
            #fstat, fpval = sp.stats.f_oneway(pert,ctl) # same as t statistic
            #print '    FSTAT: ' + str(fstat) + ' PVAL: ' + str(fpval)
            #if fpval<=siglevel:
            #    print '  **The ensemble means are significantly different (' + str(1-siglevel) + ')'

            print '    CTL std: ' + str(ctl.std()) + ', PERT std: ' + str(pert.std())
            lstat, lpval = sp.stats.levene(pert,ctl)
            print '    LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
            if lpval<=siglevel:
                print '  **The ensemble variances are significantly different (' + str(1-siglevel) + ')'
            #thestats['fpval'] = lpval

            seatpval[sea]=tpval
            seafpval[sea]=lpval
            
        alltpval[sim]=seatpval
        allfpval[sim]=seafpval

    allstats['tpval']=alltpval
    allstats['fpval']=allfpval
    
    return allstats
    # @@@@ add return vals



    
def pattcorr_ensemble(ename, field, latlim=60):
    # @@@@@@@@@@@@ is this fully implemented? Don't think so. 12/2/14
    
    if ename=='ANT':
        ename='HistIC'
    elif ename=='TOT':
        ename='HistBC'

    enssims = con.build_ensemblesims(ename)
    ensnum=len(enssims)

    # ======= copied from Canam4_BCpatterncorr-Copy0.py ===========

    #ensnum=5
    diffdict = {}
    pcmeandict = {} # fldp-fldc pattern corr compared to mean BC
    pchaddict = {} # fldp-fldc pattern corr compared to hadisst
    seadiffdict = {} # seasonal mean
    pcseameandict = {}
    pcsea2meandict = {} # to test the other pattern corr calc
    pcsea2pvalmeandict = {} # to test the other pattern corr calc

    # generate weights for the pattern corr
    lat = con.get_t63lat()
    lon = con.get_t63lon()

    areas = cutl.calc_cellareas(lat,lon)
    areas = areas[lat>latlim,:]
    weights = areas / np.sum(np.sum(areas,axis=1),axis=0)

    #for eii in range(1,ensnum+1):
    for skey in enssims:

        #skey = etype + str(eii)
        #casenamec = bcasenamec + skey
        #casenamep = bcasenamep + skey
        #fnamec = basepath + casenamec+ subdir + casenamec + '_' + field + '_001-121_ts.nc'
        #fnamep = basepath + casenamep+ subdir + casenamep + '_' + field + '_001-121_ts.nc'
        fnamec,fnamep = con.build_filepathpair(skey,field)

        # monthly calc
        fldc = cnc.getNCvar(fnamec,ncfield,timesel=timesel)*conv
        fldp = cnc.getNCvar(fnamep,ncfield,timesel=timesel)*conv
        fldd = fldp-fldc

        # take the pattern correlation
        flddclimo,flddstd = cutl.climatologize(fldd) # climo first (don't need to do for BCs technically)
        flddcclimo,flddcstd = cutl.climatologize(flddc) # climo first. baseline diff data

        diffdict[skey] = flddclimo

        # for each month, compute pattern corr
        pc = np.zeros((12))
        for mii,mon in enumerate(con.get_mon()):
            tmp = np.squeeze(flddclimo[mii,lat>latlim,...])
            tmpcmp = np.squeeze(flddcclimo[mii,lat>latlim,...])
            pc[mii] = cutl.pattcorr(tmp.flatten()*weights.flatten(),tmpcmp.flatten()*weights.flatten())

        pcmeandict[skey] = pc # monthly

        # seasonal calc    
        fldcsea = np.zeros((4,len(lat),len(lon)))
        fldpsea = np.zeros((4,len(lat),len(lon)))
        flddsea = np.zeros((4,len(lat),len(lon)))
        pcsea = np.zeros((4))
        pcsea2 = np.zeros((4)) # test pattcorr_pearson() @@
        pcsea2pval = np.zeros((4)) # test pattcorr_pearson()

        for seaii,sea in enumerate(seasons):
            fldcsea[seaii,...] = np.mean(cnc.getNCvar(fnamec,ncfield,timesel=timesel,seas=sea)*conv,axis=0)
            fldpsea[seaii,...] = np.mean(cnc.getNCvar(fnamep,ncfield,timesel=timesel,seas=sea)*conv,axis=0)
            flddsea[seaii,...] = fldpsea[seaii,...]-fldcsea[seaii,...]

            tmp = np.squeeze(flddsea[seaii,lat>latlim,...])
            tmpcmp = np.squeeze(flddcsea[seaii,lat>latlim,...])
            pcsea[seaii] = cutl.pattcorr(tmp.flatten()*weights.flatten(),
                                         tmpcmp.flatten()*weights.flatten())
            pcsea2[seaii],pcsea2pval[seaii] = cutl.pattcorr_pearson(tmp.flatten()*weights.flatten(),
                                                                    tmpcmp.flatten()*weights.flatten())


        seadiffdict[skey] = flddsea        
        pcseameandict[skey] = pcsea
        pcsea2meandict[skey] = pcsea2
        pcsea2pvalmeandict[skey] = pcsea2pval

    # ======= end copy from Canam4_BCpatterncorr-Copy0.py ===========
