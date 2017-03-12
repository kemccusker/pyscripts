"""
    rlx_utils.py
    10/24/2016

    This class is for "nudging" or relaxation runs.

"""



def calc_linpatterncorr(ncdt,lin='one',northof=0,conv=1,vert=False):
    
    """ Computes the pattern correlation b/w the 'SUM' and 'FULL'
           for given 'lin'. Input data is seasonal avg timeseries of spatial data.
           
           computes the correlation for area north of 'northof' value.
           
           returns rval,pval
           
    """
    
    pico = ncdt['prei2xco2iceb']
    pipi = ncdt['preipreiice']
    copi = ncdt['2xco2preiice']
    coco = ncdt['2xco22xco2ice']

    if lin=='one':
        ice = (np.mean(pico,axis=0) - np.mean(pipi,axis=0))*conv
        co2 = (np.mean(copi,axis=0) - np.mean(pipi,axis=0))*conv
        ctl = np.mean(pipi,axis=0)*conv        
        suff=''
        suff1='cold'; suff2='hi'
    elif lin=='two':
        ice = (np.mean(coco,axis=0) - np.mean(copi,axis=0))*conv
        co2 = (np.mean(coco,axis=0) - np.mean(pico,axis=0))*conv
        ctl = np.mean(coco,axis=0)*conv
        suff='2'
        suff1='warm'; suff2='lo'

    combo = ice+co2
    full = (np.mean(coco,axis=0) - np.mean(pipi,axis=0))*conv
    
    if vert: # lat in 2nd dim
        rval,pval = cutl.pattcorr_pearson(combo[:,lat>northof].flatten(),
                                          full[:,lat>northof].flatten())
    else: # lat in 1st dim
        rval,pval = cutl.pattcorr_pearson(combo[lat>northof,:].flatten(),
                                          full[lat>northof,:].flatten())
    #print cutl.pattcorr(combo[lat>northof,:].flatten(),full[lat>northof,:].flatten())
    
    return rval,pval

def calc_monthly_linpatterncorr(ncdt,lin='one',northof=0,conv=1,vert=False):
    
    """ Computes the pattern correlation b/w the 'SUM' and 'FULL'
           for given 'lin'. Input data is climatology of spatial data.
           
           computes the correlation for area north of 'northof' value for e/ month.
           
           returns rval,pval
           
    """
    
    pico = ncdt['prei2xco2iceb']
    pipi = ncdt['preipreiice']
    copi = ncdt['2xco2preiice']
    coco = ncdt['2xco22xco2ice']

    if lin=='one':
        ice = (pico - pipi)*conv
        co2 = (copi - pipi)*conv
        ctl = pipi*conv        
        suff=''
        suff1='cold'; suff2='hi'
    elif lin=='two':
        ice = (coco - copi)*conv
        co2 = (coco - pico)*conv
        ctl = coco*conv
        suff='2'
        suff1='warm'; suff2='lo'

    combo = ice+co2
    full = (coco - pipi)*conv
    
    #print lat>northof
    rvals=np.zeros((12,))
    pvals=np.zeros_like(rvals)
    
    for moii in range(0,12):
        if vert: # lat in 2nd dim
            rvals[moii],pvals[moii] = cutl.pattcorr_pearson(combo[moii,:,lat>northof].flatten(),
                                                            full[moii,:,lat>northof].flatten())
        else: # lat in 1st dim
            rvals[moii],pvals[moii] = cutl.pattcorr_pearson(combo[moii,lat>northof,:].flatten(),
                                                            full[moii,lat>northof,:].flatten())
    
    return rvals,pvals
    
def calc_monthly_rmse(ncdt,lin='one',norm=True,northof=0,conv=1,vert=False):
    
    """ Computes the average RMSE b/w the 'SUM' and 'FULL'
           for given 'lin'. Input data is climatology of spatial data.
           
           computes the RMSE for area north of 'northof' value for e/ month.
           
           if norm=True: normalize by spatial variance
           
           returns rmse
           
    """
    
    pico = ncdt['prei2xco2iceb']
    pipi = ncdt['preipreiice']
    copi = ncdt['2xco2preiice']
    coco = ncdt['2xco22xco2ice']

    if lin=='one':
        ice = (pico - pipi)*conv
        co2 = (copi - pipi)*conv
        ctl = pipi*conv        
        suff=''
        suff1='cold'; suff2='hi'
    elif lin=='two':
        ice = (coco - copi)*conv
        co2 = (coco - pico)*conv
        ctl = coco*conv
        suff='2'
        suff1='warm'; suff2='lo'

    combo = ice+co2
    full = (coco - pipi)*conv
    
    rmse=np.zeros((12,))
    spvar = np.zeros_like(rmse)
    for moii in range(0,12):
        if vert: # lat in 2nd dim
            print 'how to weight vertical coord? dp?'
            #rmse[moii] = cutl.calc_areawgted_rmse(combo[moii,:,lat>northof].flatten(),
            #                                                full[moii,:,lat>northof].flatten())
        else: # lat in 1st dim
            #print 'testing area weighted var func @@@'
            spvar[moii] = cutl.calc_areawgted_var(full[moii,lat>northof,:-1],
                                            lat[lat>northof],lon[:-1],model='CanESM2')
            
            rmse[moii] = cutl.calc_areawgted_rmse(full[moii,lat>northof,:-1],
                                                  combo[moii,lat>northof,:-1],
                                                  lat[lat>northof],lon[:-1],model='CanESM2')
            if norm:
                rmse[moii] = rmse[moii]/spvar[moii]
            
    print 'spvar ' +str(spvar) # @@
            
    return rmse

