""" statshelper.py
       11/24/2014
       first function: pattern correlate within ensemble
"""
import numpy as np
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


def pattcorr_ensemble(ename, field, latlim=60):

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
