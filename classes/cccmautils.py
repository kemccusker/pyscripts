"""
cccmautils.py
Put my utility functions in this module.
   e.g. area-weighted global mean, etc

   2/11/2014

"""

import numpy as np
import constants as con
import collections
import numpy.ma as ma
import scipy as sp
import scipy.stats
from scipy.stats import norm
import copy
from scipy import signal # ???

con = reload(con)

def find_nearest(array,value):
    """ find_nearest(array,value): returns index of array
        element with nearest value to value
    """
    return (np.abs(array-value)).argmin()

def updatedict(dd,ud):
    """ update a nested dictionary without overwriting everything
        http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
        dd: dictionary to update
        ud: the update dictionary
        6/5/2014: note that this didn't solve my problem: was trying to update inner keys before outter keys
                  Switched the order and didn't need this function anymore anyway.
    """
    for key, val in ud.iteritems():
        if isinstance(val, collections.Mapping):
            #print 'recurse'
            rec = updatedict(dd.get(key, {}), val)
            dd[key] = rec
        else:
            dd[key] = ud[key]
    return dd

def pattcorr(x,y):
    """ pattcorr(x,y)
    
           Return the pattern correlation coefficient between x and y (index 0,1 of matrix)
              x and y must be flattened arrays (of a 2D map)
              
           This function is to test the results from np.corrcoef() and ma.corrcoef()
              and gives the same result. @@ what about scipy.stats.pearsonr() @@ ?
               @@ appears to give same result @@
    """

    Cxy = np.cov(x,y)
    Cxx = np.cov(x)
    Cyy = np.cov(y)

    Pxy = Cxy / np.sqrt(Cxx*Cyy)

    return Pxy[0,1]

def pattcorr_pearson(x,y):
    """ pattcorr_pearson(x,y)
         
            Return the correlation coefficient and the p-value of the correlation.
            This function calls scipy.stats.pearsonr(x,y), which removes the mean field
              before computing the correlation.
              Doc: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.pearsonr.html
              Source: https://github.com/scipy/scipy/blob/v0.14.0/scipy/stats/stats.py#L2392
    """
    r, pval = sp.stats.pearsonr(x,y)

    return r, pval

def calc_mincorr(N):
    """ calc_mincorr(N):
              calculate the minimum correlation that would be significant given a sample size, N.
              ASSUMES significance level of 0.05 (95%)
    """

    # original formula (http://vassarstats.net/textbook/ch4apx.html): 
    #    t = corr / np.sqrt((1-corr*corr)/(N-2))
    
    # Solve for corr.
    # trying to figure out what minimum r-value (correlation) is significant
    # given sample size, N
    # Tried plugging 1290 into online calc to see when t > 1.9719 and r was 
    #  about 0.055 (tiny, as expected). http://vassarstats.net/textbook/ch4apx.html
    #  which matches my own test calc.

    tcrit = 1.9719 # 0.05

    denom = ((N-2)/np.power(tcrit,2)) + 1
    rmin = np.sqrt( 1/denom)

    return rmin


def calc_seaicearea(input,lat,lon, model='CanESM2'):
    """ Calculate sea ice area from sea ice concentration
        input: 2D or greater array of sea ice concentration.
               If > 2D, time must be first dimension
        lat: array of coordinate lats
        lon: array of coordinate lons

        if model=='CanESM2' (default), return areacella from file.
        if model==None, use calc_cellareas() function which uses lat/lon

        returns: MASKED array of sea ice area with same
                 shape as input. The array is masked such
                 that land is masked out
        @@@ needs testing
    """

    ishape = input.shape
    ndims=len(ishape)

    if np.mod(input.shape[-1],2) != 0: # if lon is odd, leave cyclic lon in landmask. #@@? keep extra lon tho?
        remcyc=False
    else:
        remcyc=True
        
    if ndims>2: # first dim is time
        if model=='CanESM2':
            areas = con.get_t63cellareas(repeat=ishape)
        elif model==None:
            areas = calc_cellareas(lat,lon,repeat=ishape)
        else:
            print 'model setting incorrect. Either ''CanESM2'' or None' # @@
            return -1

        lmask = con.get_t63landmask(repeat=ishape,remcyclic=remcyc) # @@ note assuming T63 here...
    else:
        if model=='CanESM2':
            areas = con.get_t63cellareas()
        elif model==None:
            areas = calc_cellareas(lat,lon)
        else:
            print 'model setting incorrect. Either ''CanESM2'' or None' # @@
            return -1

        lmask = con.get_t63landmask(remcyclic=remcyc)
    # in case input is not a MaskedArray already
    input = ma.MaskedArray(input) # @@@@ will this preserve a mask if one alread exists?
    areas = ma.masked_where(input.mask, areas)

    sia = input*areas
    sia = ma.masked_where(lmask==-1,sia) # mask where land
    
    return sia

def mask_t63land(input):
    """ masks out land (do not consider land) in input field
             using the canesm t63 landmask
             2/25/2016

        return masked input, input.mask

        @@@ check if input already has a mask?
    """

    ishape = input.shape

    nrep = ishape[0:-2] # leave off last 2 dims (lat, lon)
    nrep = nrep + (1,1) # I think this should work for 2D or >2D

    if np.mod(input.shape[-1],2) != 0: # if lon is odd, leave cyclic lon in landmask. 
        remcyc=False
    else:
        remcyc=True

    lmask = con.get_t63landmask(repeat=nrep,remcyclic=remcyc) 

    input = ma.masked_where(lmask==-1,input) # mask where land

    return input,input.mask

def mask_t63ocean(input):
    """  masks out ocean (do not consider ocean) in input field
             using the canesm t63 landmask

    """

    ishape = input.shape

    nrep = ishape[0:-2] # leave off last 2 dims (lat, lon)
    nrep = nrep + (1,1) # I think this should work for 2D or >2D

    if np.mod(input.shape[-1],2) != 0: # if lon is odd, leave cyclic lon in landmask. 
        remcyc=False
    else:
        remcyc=True

    lmask = con.get_t63landmask(repeat=nrep,remcyclic=remcyc) 

    input = ma.masked_where(lmask==0,input) # mask where ocean


def calc_totseaicearea(fld,lat,lon,isarea=False,model='CanESM2'):
    """ calculate total sea ice area
           returns nh,sh
           input is expected to be time x lat x lon SEA ICE CONC
              unless isarea=True

        @@ needs testing 9/9/2014
    """
    if isarea:
        sia=fld 
    else:
        sia = calc_seaicearea(fld,lat,lon,model=model)

    # Changed to ma.sum() on 12/10/14
    nh = ma.sum(ma.sum(sia[...,lat>0,:],axis=-1),axis=-1)
    sh = ma.sum(ma.sum(sia[...,lat<0,:],axis=-1),axis=-1)

    return nh,sh

def calc_seaiceextent(fld, lat, lon=None, model='CanESM2'):
    """ 
         fld should be time x lat x lon SEA ICE CONC [fraction]

         model: 'CanESM2' uses grid cell area from file
                None uses calc_cellareas(lat,lon)
    
         returns nhsie, shsie
    """
    ishape = fld.shape

    if model=='CanESM2':
        areas=con.get_t63cellareas(repeat=ishape)
    elif model==None:
        if lon==None:
            print 'need lon if model==None!!' # @@
            return -1

        if np.mod(len(lon),2)!=0:
            # remove extra lon before calc
            lon=lon[:-1]
            fld=fld[...,:-1]
            ishape=fld.shape

        areas=calc_cellareas(lat,lon, repeat=ishape)
    else:
        print 'model setting incorrect. Either ''CanESM2'' or None' # @@
        return -1

    tmp = copy.copy(fld)
    tmp[fld>=0.15] = 1
    tmp[fld<0.15] = 0

    nh = tmp[...,lat>0,:]*areas[...,lat>0,:]
    nhsie = ma.sum(ma.sum(nh,axis=-1),axis=-1)

    #nhsie=(tmp[...,lat>0,:]*areas[...,lat>0,:]).sum(axis=-1).sum(axis=-1)
    shsie = (tmp[...,lat<0,:]*areas[...,lat<0,:]).sum(axis=-1).sum(axis=-1)

    return nhsie,shsie


def calc_meanseaicethick(fld,sia,lat,lon):
    """ calc_meanseaicethick(input,sia,lat,lon)
                 input: 2D or greater array of sea ice thickness
                 need sea ice area? @@
    """
    
    print "Not finished! @@"

def calc_totseaicevol(fld,sic,lat,lon,repeat=False,isarea=False,model=None):
    """ calculate total ice volume
           returns nh, sh
           input is expected to be time x lat x lon in units of thickness (not mass)
           sic is fractional sea ice concentration
           if sic must be repeated to match input, set
              repeat=True (@@ not yet tested)
           set isarea=True if sic input is actually sea ice area,
              otherwise sia will be calc'd from incoming sic

       @@@@ NOTE for CMIP sit, this function is incorrect.
            do not need sea ice area to calc SIT.
            See calc_totseaicevol_cmip5() 10/1/2015

    """

    if repeat:
        nrep = repeat[0:-2] # leave off last 2 dims (lat, lon)
        nrep = nrep + (1,1)
        sic = np.tile(sic,nrep)
        

    if isarea==False:
        sia=calc_seaicearea(sic,lat,lon,model=model)
    else:
        sia=sic
    
    gridvoln = fld[:,lat>0,...]*sia[:,lat>0,...] # north
    gridvols = fld[:,lat<0,...]*sia[:,lat<0,...] # south
    # Changed to ma.sum() on 12/10/14
    nh = ma.sum(ma.sum(gridvoln,2),1)
    sh = ma.sum(ma.sum(gridvols,2),1)
    
    return nh,sh
  
def calc_totseaicevol_cmip5(fld,lat,lon, model='CanESM2'):
    """ calculate total ice volume
           returns nh, sh
           input is expected to be time x lat x lon in units of 
              thickness over ocean portion of grid cell

          For cmip5 input, just need to multiply by grid cell area to
            get volume
    """

    ishape = fld.shape

    if model=='CanESM2':
        areas=con.get_t63cellareas(repeat=ishape)
    elif model==None:
        areas=calc_cellareas(lat,lon, repeat=ishape)
    else:
        print 'model setting incorrect. Either ''CanESM2'' or None' # @@
        return -1
    
    gridvoln = fld[:,lat>0,...]*areas[:,lat>0,...] # north
    gridvols = fld[:,lat<0,...]*areas[:,lat<0,...] # south
    # Changed to ma.sum() on 12/10/14
    nh = ma.sum(ma.sum(gridvoln,2),1)
    sh = ma.sum(ma.sum(gridvols,2),1)
    
    return nh,sh

def rem_seasonalcycle(input, remclimo=None):
    """ 
        Remove seasonal cycle from timeseries.
        
            input: the data. time must be first dimension and start with Jan!
            
            remclimo: if given, remove this climatology, otherwise
                      a climatology will be computed over entire input and removed.
                      The shape must be the same as input.shape[1:]
                      
    """
    
    print input.shape
    
    ntime = input.shape[0]
    nyr = ntime/12
    rem = np.mod(ntime,12) # remaining months not making up a full year
    
    if remclimo==None:
        
        remclimo,std = cutl.climatologize(input)
        
        
    otherdims=input.shape[1:]
    oshape = tuple(np.ones(len(otherdims)))
    
    print 'nyr, otherdims, oshape ' + str(nyr),str(otherdims),str(oshape)
    print 'remclimo.shape ' + str(remclimo.shape)
    
    #if len(otherdims)==1:
    #    tshape = (nyr,)
    #else:
    tshape = (nyr,)+oshape
        
    print 'tshape ' + str(tshape)
    
    remclimot=np.tile(remclimo,tshape)

    # remove the monthly climo
    inrem = input-remclimot
        
    return inrem  
    
def global_mean_areawgted3d(fld, lat, lon, model='CanESM2'):
    """
        if model=='CanESM2': use areacella from file
        if model==None: use calc_cellareas() using lat/lon
    """

    earthrad = con.get_earthrad()
    totalarea = 4*np.pi*earthrad**2

    if model=='CanESM2':
        cellareas = con.get_t63cellareas()
    elif model==None:
        cellareas = calc_cellareas(lat,lon)
    else:
        print 'model setting incorrect. Either ''CanESM2'' or None' # @@
        return -1

    nt = fld.shape[0]
    wgts = cellareas/np.float(totalarea)
    
    gm = np.zeros(nt)    
    for tidx in range(0,nt):
        gm[tidx] = np.average(fld[tidx,:,:],weights=wgts)
        
    return gm
   
def global_mean_areawgted(fld, lat, lon,model='CanESM2'):
    """
        if model=='CanESM2': use areacella from file
        if model==None: use calc_cellareas() using lat/lon
    """

    earthrad = con.get_earthrad()
    totalarea = 4*np.pi*earthrad**2

    if model=='CanESM2':
        cellareas = con.get_t63cellareas()
    elif model==None:
        cellareas = calc_cellareas(lat,lon)
    else:
        print 'model setting incorrect. Either ''CanESM2'' or None' # @@
        return -1

    #nt = fld.shape[0]
    wgts = cellareas/np.float(totalarea)
    
    gm = ma.average(fld,weights=wgts)
        
    return gm

def polar_mean_areawgted3d(fld,lat,lon,latlim=60,hem='nh',cellareas=None,includenan=False, model='CanESM2'):
    """ Pass in cellareas if you want some of the cells masked, e.g. if masking ocean or land.
        This ONLY works if the mask is the same for each time in timeseries.
            includenan: whether or not to consider NaN values. if True, a mean of a NaN is NaN.
                        Default is False

        if model=='CanESM2': use areacella from file
        if model==None: use calc_cellareas() using lat/lon
    """

    if cellareas == None:
        if model=='CanESM2':
            cellareas = con.get_t63cellareas()
        elif model==None:
            cellareas = calc_cellareas(lat,lon)
        else:
            print 'model setting incorrect. Either ''CanESM2'' or None' # @@
            return -1

    # else cellareas are provided
        
    # @@ what does this do if just an ocean field, etc?

    if hem=='nh':
        polcellareas = cellareas[lat>=latlim,:]
        polinput = fld[:,lat>=latlim,:]
    else:
        polcellareas = cellareas[lat<=latlim,:]
        polinput = fld[:,lat<=latlim,:]

    totalarea=float(np.sum(np.sum(polcellareas)))

    nt = fld.shape[0]
    wgts = polcellareas/totalarea
    #print wgts #@@
    
    pm = np.zeros(nt)    
    for tidx in range(0,nt):
        if includenan:
            pm[tidx] = np.average(polinput[tidx,:,:],weights=wgts)
        else:
            polinput=ma.masked_where(np.logical_or(polinput==np.nan, polinput==np.inf),polinput)
            pm[tidx] = ma.average(polinput[tidx,:,:],weights=wgts)
        
    return pm

def annualize_monthlyts(input, includenan=0,verb=False):
    """ It is expected that the input variable is 3 dimension, with time being first
        But on 5/6/2014 I tried to generalize to any dimension...needs to be tested.
    """

    # time must be first dimension
    # figure out nyears (truncate if only some months in last year)
    dims = input.shape
    if verb:
        print 'dims ' + str(dims)

    nt = input.shape[0] # number of times
    nyrs = nt/12        # I think this will automatically return an integer

    #ndim2 = input.shape[1]
    #ndim3 = input.shape[2]
    
    wgts = con.get_monweights()

    if len(dims)>1:
        diml = list(dims[1:]) # this one is for wgts array
        diml.append(1) # insert size nyrs into first dimension
        diml = tuple(diml)
        if verb:
            print 'diml ' + str(diml)

        diml2 = list(dims[1:]) # this one is for initializing the timeseries array
        diml2.insert(0,nyrs) # insert size nyrs into first dimension
        if verb:
            print 'diml2 ' + str(diml2) 

            print 'wgts.shape ' + str(wgts.shape)

        wgts = np.tile(wgts,diml) # make a weights array for each grid pt
        if verb:
            print 'wgts.shape after tile ' + str(wgts.shape)

        if len(dims)==3:
            # put the 3rd dim (2) in first spot, 1st dim (0) in second spot, 2nd dim (1) in third spot
            # this effectively shifts the dimensions back to time,dim2,dim3
            wgts = np.transpose(wgts,(2,0,1))
        elif len(dims)==2:
            wgts = wgts.T #np.transpose(wgts,(2,0,1))
    else:
        diml2 = nyrs,

    #yidx=0
    counter=0
    annts = ma.zeros(diml2)
    
    # Loop through 12-month chunks
    for yidx in range(0, nyrs):
        annts[yidx,...] = ma.average(input[counter:counter+12,...],axis=0,weights=wgts)
        counter = counter+12

    return annts

def seasonalize_monthlyts(input,season=None,includenan=0,mo=0,climo=0,verb=False):
    """ suggest using seasonalize() instead. it will not chop off remainder months
          if doing a seasonal avg that spans a year transition.
    """

    if season == None and mo==0:
        print 'Must specify either season or mo! Months indexed starting from 1'
        return

    dims = input.shape
    nt = input.shape[0]
    nyrs = nt/12 # should be integer
    incr=3 # default increment is 3 (for 3-month "seasons")
    
    wgts = con.get_monweights()

    if includenan == 1:
        print 'this flag not yet implemented'
        return
    
    if mo != 0:
        # simply pick out the month and return ts
        incr=0
        seasts = input[(mo-1)::12,...] # increments of 12 months
        return seasts
    elif season=='ANN':
        return annualize_monthlyts(input,includenan,verb=verb)
    elif season=='DJF':
        # do weighted average for requested season
        
        if climo==1:
            subwgts = wgts[[0,1,11]]
            seasts = ma.average(input[[0,1,11],...],weights=subwgts,axis=0)
            return seasts

        subwgts = wgts[[0,1,11]] 
        nyrs=nyrs-1
        start=11 # start with December of first year (skip that year's JF)
    elif season=='JJA':
        subwgts = wgts[5:8] # doesn't include index 8
        start=5
    elif season=='NDJ':
        if climo==1:
            subwgts = wgts[[0,10,11]]
            seasts = ma.average(input[[0,10,11],...],weights=subwgts,axis=0)
            return seasts
        
        subwgts = wgts[[0,10,11]]
        start=10 # start with Nov of first year (skip that year's J)
        nyrs=nyrs-1
    elif season=='MAM':
        subwgts = wgts[2:5]
        start=2
    elif season=='SON':
        subwgts = wgts[8:11]
        start=8
    elif season=='ND':
        subwgts = wgts[10:]
        start=10
        incr=2
    elif season=='JF':
        subwgts = wgts[0:2]
        start=0
        incr=2
    elif season=='SO':
        subwgts = wgts[8:10]
        start=8
        incr=2
    else:
        print 'seasonalize_monthlyts(): Season ' + season + ' not supported!'
        return None

    if len(dims)>1:
        diml = list(dims[1:]) # this one is for subwgts array
        diml.append(1) # insert size nyrs into first dimension
        diml = tuple(diml)
        
        diml2 = list(dims[1:]) # this one is for initializing the timeseries array
        diml2.insert(0,nyrs) # insert size nyrs into first dimension

        subwgts = np.tile(subwgts,diml)
        # put the 3rd dim (2) in first spot, 1st dim (0) in second spot, 2nd dim (1) in third spot
        # this effectively shifts the dimensions back to time,dim2,dim3
        if len(dims)==2:
            subwgts = np.transpose(subwgts,(1,0))
        else: # assume 3 dims
            subwgts = np.transpose(subwgts,(2,0,1))
    
    else:
        diml2 = nyrs,


    seasts = ma.zeros(diml2)
    #print subwgts.shape
    
    for yridx in range(0,nyrs):
        subsamp = range(start+yridx*12,start+yridx*12+incr)
        seasts[yridx,...] = ma.average(input[subsamp,...],weights=subwgts,axis=0)        
    # @@ implement includenan? If set = 1, then the mean will not ignore a NaN
    #    such that if a NaN is present, the mean is NaN
    
    return seasts



def seasonalize(input,season=None,includenan=0,mo=0,climo=0,verb=False):
    """ new version of seasonalize_monthlyts()

             first dim must be time. first time index is expected to be Jan.
             mo is 1-based! e.g. 1 is Jan, etc.
    """

    if season == None and mo==0:
        print 'No season or mo specified! Returning all months. Note month starts with 1.'
        return input

    dims = input.shape
    nt = input.shape[0]
    nyrs = nt/12 # should be integer
    incr=3 # default increment is 3 (for 3-month "seasons")
    
    wgts = con.get_monweights()
    dpm = con.get_dayspermon()

    if includenan == 1:
        print 'this flag not yet implemented'
        return
    
    if mo != 0:
        # simply pick out the month and return ts
        incr=0
        seasts = input[(mo-1)::12,...] # increments of 12 months
        return seasts
    elif season=='ANN':
        return annualize_monthlyts(input,includenan,verb=verb)
    elif season=='DJF':
        # do weighted average for requested season
        indices=[0,1,11] # don't worry, these indices only for weights
        if np.mod(input.shape[0],12)>=2:
            # then we have Jan-Feb of last year
            nyrs=nt/12
        else:
            nyrs=nt/12-1
            
        subwgts = dpm[indices]/np.float(np.sum(dpm[indices]))

        if climo==1:
            #subwgts = wgts[indices]
            seasts = ma.average(input[indices,...],weights=subwgts,axis=0)
            return seasts

        #subwgts = wgts[indices] 
        
        #nyrs=nyrs-1
        start=11 # start with December of first year (skip that year's JF)
    elif season=='JJA':
        indices = [5,6,7]
        subwgts = dpm[indices]/np.float(np.sum(dpm[indices]))
        #subwgts = wgts[5:8] # doesn't include index 8
        start=5
    elif season=='NDJ':
        indices = [0,10,11]
        if np.mod(input.shape[0],12)>=1:
            # then we have Jan of last year
            nyrs=nt/12
        else:
            nyrs=nt/12-1
            
        subwgts = dpm[indices]/np.float(np.sum(dpm[indices]))

        if climo==1:
            #subwgts = wgts[[0,10,11]]
            seasts = ma.average(input[indices,...],weights=subwgts,axis=0)
            return seasts
        
        #subwgts = wgts[[0,10,11]]
        start=10 # start with Nov of first year (skip that year's J)
        #nyrs=nyrs-1
    elif season=='MAM':
        indices=[2,3,4]
        subwgts = dpm[indices]/np.float(np.sum(dpm[indices]))
        #subwgts = wgts[2:5]
        start=2
    elif season=='SON':
        indices=[8,9,10]
        subwgts = dpm[indices]/np.float(np.sum(dpm[indices]))
        #subwgts = wgts[8:11]
        start=8
    elif season=='ND':
        indices=[10,11]
        subwgts = dpm[indices]/np.float(np.sum(dpm[indices]))
        #subwgts = wgts[10:]
        start=10
        incr=2
    elif season=='JF':
        indices=[0,1]
        subwgts = dpm[indices]/np.float(np.sum(dpm[indices]))
        #subwgts = wgts[0:2]
        start=0
        incr=2
    elif season=='SO':
        indices = [8,9]
        subwgts = dpm[indices]/np.float(np.sum(dpm[indices]))
        #subwgts = wgts[8:10]
        start=8
        incr=2
    else:
        print 'seasonalize(): Season ' + season + ' not supported!'
        return None

    if len(dims)>1:
        diml = list(dims[1:]) # this one is for subwgts array
        diml.append(1) # insert size nyrs into first dimension
        diml = tuple(diml)
        
        diml2 = list(dims[1:]) # this one is for initializing the timeseries array
        diml2.insert(0,nyrs) # insert size nyrs into first dimension

        subwgts = np.tile(subwgts,diml)
        # put the 3rd dim (2) in first spot, 1st dim (0) in second spot, 2nd dim (1) in third spot
        # this effectively shifts the dimensions back to time,dim2,dim3
        if len(dims)==2:
            subwgts = np.transpose(subwgts,(1,0))
        else: # assume 3 dims
            subwgts = np.transpose(subwgts,(2,0,1))
    
    else:
        diml2 = nyrs,


    seasts = ma.zeros(diml2)
    #print subwgts.shape
    
    for yridx in range(0,nyrs):
        subsamp = range(start+yridx*12,start+yridx*12+incr)
        seasts[yridx,...] = ma.average(input[subsamp,...],weights=subwgts,axis=0)        
    # @@ implement includenan? If set = 1, then the mean will not ignore a NaN
    #    such that if a NaN is present, the mean is NaN
    
    return seasts

        
def calc_cellareas(lat,lon, repeat=None):

    """ assumes longitudes are evenly spaced and lat and lon cover whole globe
        DOES NOT work for fraction of globe

        repeat should be the shape of the field that
        we want to multiply cellareas by. Assume the last
        two numbers are lat and lon sizes
        If None, no repeating necessary
        """

    nlat = lat.shape[0]
    nlon = lon.shape[0]
    earthrad = con.get_earthrad() #6.37122e6
    # sfc area of Earth (from google): 510,072,000 km2
    deg2rad = np.pi/180.0
    
    if lat[0] > 0:
        # this means the lat array is flipped from what is expected
        lattmp = np.copy(lat)
        lattmp = np.flipud(lattmp)
        lat_bounds = np.array([(lattmp[1:] + lattmp[0:-1])/2.])
    else:
        lat_bounds = np.array([(lat[1:] + lat[0:-1])/2.])
        
#    lat_bounds = np.array([(lat[1:nlat] + lat[0:(nlat-1)]) / 2])
    #print lat_bounds # @@

    lat_bounds = np.insert(lat_bounds,0,np.float32(-90))
    lat_bounds = np.append(lat_bounds,np.float32(90))

    del_phi = np.sin(lat_bounds[1:]*deg2rad) - np.sin(lat_bounds[0:-1]*deg2rad)
    #print del_phi # @@
    
    # convert latcircle radians to latcircle areas and
    # divide by number of longitudes to get area of each cell
    cellareas = del_phi*2*np.pi*earthrad**2 / float(nlon)
    
    # need to repmat this into number of longitudes
    cellareas = np.tile(cellareas,(nlon,1))
    cellareas = np.transpose(cellareas,(1,0))

    if repeat != None:
        nrep = repeat[0:-2] # leave off last 2 dims (lat, lon)
        nrep = nrep + (1,1)
        cellareas = np.tile(cellareas,nrep)
        
    return cellareas
                     
def get_cellwgts(lat,lon,repeat=None, model='CanESM2'):

    """ basically, repeat should be the shape of the field that
    #   we want to weight with cell weights. Assume the last
    #   two numbers are lat and lon sizes
    #   If None, no repeating necessary

        if model=='CanESM2': use areacella from file
        if model==None: use calc_cellareas() using lat/lon
    """
    
    earthrad = con.get_earthrad()
    totalarea = 4*np.pi*earthrad**2
    if model=='CanESM2':
        print 'totalarea from earthrad= ' + str(totalarea) # @@
        print 'total area from sum(grid cell areas)= ' + str(con.get_t63cellareas().sum()) # @@
        print 'diff= ' + str(totalarea-con.get_t63cellareas().sum())

        wgts = con.get_t63cellareas()/totalarea
    elif model==None:
        wgts = calc_cellareas(lat,lon)/totalarea
    else:
        print 'model setting incorrect. Either ''CanESM2'' or None' # @@
        return -1

    if repeat != None:
        nrep = repeat[0:-2] # leave off last 2 dims (lat, lon)
        nrep = nrep + (1,1)
        wgts = np.tile(wgts,nrep)

    return wgts

def climatologize3d(input):
    """ Assumes input is 3D
           time x lat x lon or time x lev x lat
    """
    (nt, ndim2, ndim3) = input.shape

    skip = np.mod(nt,12)
    # chop off the last months
    nt = nt - skip

    climo = np.zeros((12,ndim2,ndim3))
    stddev = np.zeros((12,ndim2,ndim3))

    for monidx in range(0,12):
        mon = np.arange(monidx,nt,12) # steps of 12

        climo[monidx,:,:] = np.mean(input[mon,:,:],axis=0)
        stddev[monidx,:,:] = np.std(input[mon,:,:],axis=0)

    return (climo,stddev)

def climatologize(input):
    """ Does not care how many dims! Time MUST be first though.
    """

    dims = input.shape
    nt = dims[0]

    skip = np.mod(nt,12)
    # chop off the last months
    nt = nt - skip

    if len(dims)==1:
        climo = np.zeros((12))
        stddev = np.zeros((12))
    else:
        diml = list(dims[1:])
        diml.insert(0,12) # insert size 12 into first dimension
        diml = tuple(diml)
        climo = np.zeros(diml)
        stddev = np.zeros(diml)

    for monidx in range(0,12):
        mon = np.arange(monidx,nt,12) # steps of 12

        climo[monidx,...] = np.mean(input[mon,...],axis=0)
        stddev[monidx,...] = np.std(input[mon,...],axis=0)

    return (climo,stddev)

def pooledstd(xx,yy):
    """
      pooledstd(xx,yy): determine the pooled standard deviation of the given data.
                        xx is the 'perturbation' and yy is the 'control'
                        ASSUMES time dimension is axis=0 (first)

                        returns: sp
    """
    #  sp = sqrt(  sum1_to_n[ (xi-xbar)^2 ] + sum1_to_m[ (yi-ybar)^2 ] / (n+m-2)  )
    ## tcritbig = 1.9719
    ## tcritsmall = 1.9803
    ## tcritsmallw = 1.9806

    N = xx.shape[0]
    M = yy.shape[0]
    dof = N+M-2

    xbar = np.average(xx,0)
    ybar = np.average(yy,0)
    xsum = np.sum((xx-xbar)**2,0)
    ysum = np.sum((yy-ybar)**2,0)

    sp = np.sqrt( (xsum + ysum) / dof )
    return sp

"""
   calc_Nmin(xx,yy,tcrit=1.9719):
   
       calculate minimum ensemble size for given data.
       xx is the 'perturbation' and yy is the 'control'
       ASSUMES time dimension is axis=0 (first)
       ASSUMES significance: 0.05 (95%)
       Default critical tvalue is for dof 200: 1.9719
           Others: (half timeseries) dof 118: 1.9803
                   (half timeseries, winter) dof 116: 1.9806
                                
       returns: Nmin
"""
def calc_Nmin(xx,yy,tcrit=1.9719):

    xbar = np.average(xx,0)
    ybar = np.average(yy,0)
    sp = pooledstd(xx,yy)
    
    Nmin = 2*(tcrit**2)*(sp/(xbar-ybar))**2
    return Nmin


def mask_region(fld,lat,lon,region,limsdict=None):
    """ mask_region(fld, lat,lon,limsdict):
                 Mask the input data with the given region either defined
                    already in regiondict, or overridden with limsdict.
                    If limsdict is to be used, set region='other'
                 The data will be masked everywhere BUT the region!

                 fld: lat x lon array of data
                      time x lat x lon OR lev x lat x lon also accepted (not tested)

                 lat: array of lats (coordinates)
                 lon: array of lons (coordinates)
                 region: named region (defined in constants.py)
                 limsdict: a dictionary of 'latlim' array and 'lonlim' array
                           (ie from a region dictionary, or user-defined)

                 Returns: Tuple of masked field, mask
    """

    # @@ if given lons from -180 to 180, have to modify this:
    if np.any(lon<0):
        tmplon=lon
        # will convert the negative lons to positive lons >180
        # The first bit on RHS gets diff from 180
        tmplon[lon<0] = (lon[lon<0] + 180) + 180 
        lon = tmplon

    lons,lats = np.meshgrid(lon,lat)    
        
    # create mask
    if region=='other':
        pass # use given limsdict
    else:
        limsdict = con.get_regionlims(region)
        
    latlims = limsdict['latlims']
    lonlims = limsdict['lonlims']    

    reglatsbool = np.logical_and(lat>latlims[0],lat<latlims[1])
    reglonsbool = np.logical_and(lon>lonlims[0],lon<lonlims[1])

    # create mask of everything but the region of interest
    regmask = np.logical_or( 
                            np.logical_or(lats<latlims[0],lats>latlims[1]), 
                            np.logical_or(lons<lonlims[0],lons>lonlims[1]))
    if fld.ndim>2:
        ndim1 = fld.shape[0]
        regmaskt = np.tile(regmask,(ndim1,1,1))
    else:
        regmaskt = regmask

    fld = ma.masked_where(regmaskt,fld)

    return fld, regmaskt

def calc_regmean(fld,lat,lon,region,limsdict=None, model='CanESM2',alsomask=None):
    """ calc_regmean(fld, lat,lon,region,limsdict=None):
                 Mask the input data with the given region either defined
                    already in regiondict, or overridden with limsdict.
                    If limsdict is to be used, set region='other'
                 The data will be masked everywhere BUT the region!

                 fld: lat x lon array of data
                      time x lat x lon OR lev x lat x lon also accepted (not tested)

                 lat: array of lats (coordinates)
                 lon: array of lons (coordinates)
                 region: named region (defined in constants.py)
                 limsdict: a dictionary of 'latlim' array and 'lonlim' array
                           (ie from a region dictionary, or user-defined)

                 model: 'CanESM2', gets grid cell areas from file.
                        None, uses calc_cellareas(lat,lon)

                 alsomask: None, 'land', or 'ocean'. In addition to region mask, 
                            mask out land or ocean (mask out mean do NOT consider it).

                 Returns: Regional mean (or series of regional means with length ndim1)
    """

    # check if fld has extra lon, if so, remove it before doing reg avg
    if np.mod(fld.shape[-1],2) != 0:
        fld=fld[...,:-1]
        print 'calc_regmean() removing extra lon. fld new shape: ' + str(fld.shape)
    if np.mod(lon.shape[0],2) != 0:
        lon=lon[:-1]
        #print 'calc_regmean() removing extra lon. lon new shape: ' + str(lon.shape)

    if region=='gm':
        print 'Global average!'
        if fld.ndim>2: # @@ hack. just make global mean func better
            fldreg = global_mean_areawgted3d(fld,lat,lon,model=model)
        else:
            fldreg = global_mean_areawgted(fld,lat,lon,model=model)
    else:
        fldm,regmask = mask_region(fld,lat,lon,region,limsdict)
        if alsomask != None:
            if alsomask=='land':
                fldm,regmask = mask_t63land(fldm)
            elif alsomask=='ocean':
                fldm,regmask = mask_t63ocean(fldm)
            else:
                print 'alsomask = land or ocean! @@@'

        # calculate area-weights
        if model=='CanESM2':
            areas = con.get_t63cellareas()
        elif model==None:
            areas = calc_cellareas(lat,lon)
        else:
            print 'model setting incorrect. Either ''CanESM2'' or None' # @@
            return -1

        if regmask.ndim>2:
            areasm = ma.masked_where(np.squeeze(regmask[0,...]),areas)
        else:
            areasm = ma.masked_where(regmask,areas)
        weightsm = areasm / ma.sum(ma.sum(areasm,axis=1),axis=0) # weights masked

        if fld.ndim>2:
            ndim1 = fld.shape[0]
            weightsmt = np.tile(weightsm,(ndim1,1,1)) # weights masked tiled
        else:
            weightsmt = weightsm

        tmp = ma.masked_where(regmask,fldm)
        if tmp.ndim==3:
            tmpreg = ma.sum(ma.sum(tmp*weightsmt,axis=2),axis=1)
        elif tmp.ndim==2:
            tmpreg = ma.sum(tmp*weightsmt)

        fldreg = tmpreg # should be ndim1 of regional mean (or just one regional mean)

    return fldreg


def calc_regtotseaicearea(fld,lat,lon,region,limsdict=None,isarea=False,model='CanESM2',alsomask=None):
    """ calc_regtotseaicearea(fld, lat,lon,region,limsdict=None,isarea=False):
                 Mask the input data with the given region either defined
                    already in regiondict, or overridden with limsdict.
                    If limsdict is to be used, set region='other'
                 The data will be masked everywhere BUT the region!

                 fld: lat x lon array of SEA ICE CONC (as fraction; unless isarea=True)
                      time x lat x lon OR lev x lat x lon also accepted (not tested)

                 lat: array of lats (coordinates)
                 lon: array of lons (coordinates)
                 region: named region (defined in constants.py)
                 limsdict: a dictionary of 'latlim' array and 'lonlim' array
                           (ie from a region dictionary, or user-defined)
                 isarea: specifies whether the incoming data is already converted
                         to sea ice area from concentration

                 model: 'CanESM2', gets grid cell areas from file.
                        None, uses calc_cellareas(lat,lon)

                 alsomask: None, 'land', or 'ocean'. In addition to region mask, 
                            mask out land or ocean (mask out means do NOT consider it).

                 Returns: Regional total SIA (or series of regional totals with length ndim1)
    """
    if isarea:
        sia=fld
    else: # note that this function also masks out land
        sia = calc_seaicearea(fld,lat,lon,model=model)
    
    fldm,regmask = mask_region(sia,lat,lon,region,limsdict)

    ## # calculate area-weights
    ## areas = calc_cellareas(lat,lon)
    ## if regmask.ndim>2:
    ##     areasm = ma.masked_where(np.squeeze(regmask[0,...]),areas)
    ## else:
    ##     areasm = ma.masked_where(regmask,areas)
    ## weightsm = areasm / np.sum(np.sum(areasm,axis=1),axis=0) # weights masked

    ## if fld.ndim>2:
    ##     ndim1 = fld.shape[0]
    ##     weightsmt = np.tile(weightsm,(ndim1,1,1)) # weights masked tiled
    ## else:
    ##     weightsmt = weightsm

    tmp = ma.masked_where(regmask,fldm) # am I masking out twice? does it matter?
    if alsomask != None:
        if alsomask=='land':
            tmp,regmask = mask_t63land(tmp)
        elif alsomask=='ocean':
            tmp,regmask = mask_t63ocean(tmp)
        else:
            print 'alsomask = land or ocean! @@@'

    #tmpreg = np.sum(np.sum(tmp,axis=2),axis=1)
    if tmp.ndim==3:
        tmpreg=ma.sum(ma.sum(tmp,axis=2),axis=1) # @@ was this a bug?
    elif tmp.ndim==2:
        tmpreg=ma.sum(tmp)

    fldreg = tmpreg # should be ndim1 of regional mean (or just one regional mean)

    return fldreg


def calc_monthlystd(fld):
    """ calc_monthlystd(fld):
               Given a monthly timeseries (fld) of size time [x lat x lon],
               return a std deviation for each month such that the return dimension
               is 12 x [shape[1:]].

        Note, can also get this info from climatologize(), which returns
          (climo,stddev)
               
    """
    dims = fld.shape
    nt = dims[0]

    if len(dims)==1:
        stddev = np.zeros((12))
    else:
        diml = list(dims[1:])
        diml.insert(0,12) # insert size 12 into first dimension
        diml = tuple(diml)
        stddev = np.zeros(diml)
        
    mo = np.arange(0,12)

    for moidx in mo:

        stddev[moidx] = np.std(fld[moidx::12,...],axis=0)

    return stddev

def calc_monthlytstat(input1,input2):
    """ calc_monthlytstat(input1,input2):
            Testing input1 against input2 (think of input2 as a control)

            return tuple of tstat,pval with shape 12 [x shape[1:]]
    """

    dims = input1.shape
    nt = dims[0]

    if len(dims)==1:
        tstat = np.zeros((12))
        pval = np.zeros((12))
    else:
        diml = list(dims[1:])
        diml.insert(0,12) # insert size 12 into first dimension
        diml = tuple(diml)
        tstat = np.zeros(diml)
        pval = np.zeros(diml)
        
    mo = np.arange(0,12)

    for moidx in mo:
        tstat[moidx,...],pval[moidx,...] =sp.stats.ttest_ind(input1[moidx::12,...],
                                                             input2[moidx::12,...],
                                                             axis=0)

    return (tstat,pval)

def calc_monthlysigarea(input1,input2,siglevel=0.05,latlim=60,region=None, model='CanESM2'):
    """ calc_monthlysigarea(input1,input2,siglevel=0.05,latlim=60,region=None)
              If 'region' is set, it supercedes latlim!!

              model: 'CanESM2' gets grid cell area from file
                     None uses calc_cellareas(lat,lon)

              returns: sigarea AS PERCENT OF DOMAIN
    """
    
    tstat,pval = calc_monthlytstat(input1,input2) # climo of lat x lon
    
    lat = con.get_t63lat()
    lon = con.get_t63lon()
    lons,lats = np.meshgrid(lon,lat)
    
    if model=='CanESM2':
        cellareas=con.get_t63cellareas(repeat=pval.shape)
    elif model==None:
        cellareas = calc_cellareas(lat,lon,repeat=pval.shape)
    else:
        print 'model setting incorrect. Either ''CanESM2'' or None' # @@
        return -1

    totmask = cellareas # for computing total area
    
    amask = ma.masked_where(pval>siglevel,cellareas) # mask out non-sig cells
    

    # now need to mask out regions of globe we don't care about
    if region != None:
        # region supercedes latlim!
        amask = mask_region(amask,lat,lon,region)
        totmask = mask_region(totmask,lat,lon,region)
    
    elif latlim != None:
        regmask=lats>latlim        
        regmask = np.tile(regmask,(12,1,1))

        amask = ma.masked_where(regmask,amask)
        totmask = ma.masked_where(regmask,totmask)

    mo = np.arange(0,12)
    sigarea= np.zeros(mo.shape)
    totarea= np.zeros(sigarea.shape)
    for moidx in mo:
        sigarea[moidx] = np.sum(amask[moidx,...])
        totarea[moidx] = np.sum(totmask[moidx,...])

    #print totarea
    #print sigarea
    return sigarea/(totarea)*100 # as a percent of total area
    
    

def calc_pvals(pert,ctl,axis=0,center='mean',verb=True,effdof=False,siglevel=0.05):
    """ calc_pvals(pert,ctl, verb=True,siglevel=0.05):
               Calculate whether dataset means and variances are
                 significantly different from each other.
                 
               pert and ctl are the data that will be tested.
               Ttest will be tested over axis=0 by default.
               center='mean' is for f-test (levene())
                  from the doc: 'center : {'mean', 'median', 'trimmed'}, optional
                                 Which function of the data to use in the test.
                                 The default is 'median'. Three variations of
                                 Levene's test are possible. The possibilities
                                 and their recommended usages are:
                                   
              * 'median' : Recommended for skewed (non-normal) distributions>
              * 'mean' : Recommended for symmetric, moderate-tailed distributions.
              * 'trimmed' : Recommended for heavy-tailed distributions. ' 
               
               verb: verbose or not. =True will output values and
                     whether or not they are significant based on
                     the input siglevel
               effdof: if True, use autocorrelation to calc effective degrees of
                       freedom in each timeseries (see calc_effectiveDOF() )
                       for use in cutl.ttest_ind()

               returns: (tstat,tpval,fstat,fpval)
    """

    if effdof:
        tstat,tpval = ttest_ind(pert,ctl,axis=axis,effdof=effdof)
    else:
        tstat, tpval = sp.stats.ttest_ind(pert,ctl,axis=axis) 
    lstat, lpval = sp.stats.levene(pert,ctl,center=center)

    if verb:
        print '    DIFF: ' + str(pert.mean()-ctl.mean()) + ', CTL mean: ' + str(ctl.mean()) + ', PERT mean: ' + str(pert.mean())
        print '    TSTAT: ' + str(tstat) + ' PVAL: ' + str(tpval)
        if tpval<=siglevel:
            print '  **The ensemble means are significantly different (' + str(1-siglevel) + ')'
        print '    DIFF: ' + str(pert.std()-ctl.std()) + ', CTL std: ' + str(ctl.std()) + ', PERT std: ' + str(pert.std())
        print '    LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
        if lpval<=siglevel:
            print '  **The ensemble variances are significantly different (' + str(1-siglevel) + ')'

    return (tstat,tpval,lstat,lpval)




def autocorr(input,axis=None):
    """ autocorr(input)
              1/11/2015
              
              auto correlate input with itself. The mean will be removed first.
              return all auto correlations from time 0 to time n

              @@ add ability to handle multiple dimensions
    """

    inp = input-input.mean()
    
    icorr = np.correlate(inp, inp, mode='full')
    # from numpy.convolve(a,v,mode='full')
    # mode : {'full', 'valid','same'}, optional
    #   'full': This returns the convolution at each point of
    #          overlap, with an output shape of (N+M-1,). At the end-points of the
    #          convolution, the signals do not overlap completely, and boundary effects may be seen.
    #   'same': Mode same returns output of length max(M, N). Boundary effects are still visible.
    #  'valid': Mode valid returns output of length max(M, N) - min(M, N) + 1. The convolution
    #          product is only given for points where the signals overlap completely. Values
    #          outside the signal boundary have no effect.

    icorr = icorr / icorr[icorr.argmax()]
    icorrh = icorr[icorr.size / 2:] # half of the correlations (it is symmetric)
    icorrhmax = icorrh.max()

    print 'icorrhmax: ' + str(icorrhmax)
    #print 'icorrh: ' + str(icorrh)
    #print 'icorr: ' + str(icorr)

    return icorrh
    

def lag1_autocorr(input,axis=None):
    """ lag1_autocorr(input)
              1/12/2015
              auto correlate input with itself.
              return lag-1 autocorrelation

              @@ add ability to handle multiple dimensions
    """

    icorrs = autocorr(input)
    return icorrs[1]
    

def calc_effectiveDOF(input,axis=None):
    """ def calc_effectiveDOF(input):
          1/12/2015
          
          r1 is lag-one autocorrelation of input, then n_eff:
          
          n_eff = n * (1-r1) / (1+r1)  [Santer el al. 2000]

          return n_eff

          @@ add ability to handle multiple dimensions
          
    """

    n = input.size
    print 'n: ' + str(n)
    
    r1 = lag1_autocorr(input)
    print 'r1: ' + str(r1)
    r1=np.abs(r1)
    frac=(1-r1) / (1+r1)
    #print 'frac: ' + str(frac)
    
    n_eff = n * (1-r1) / (1+r1)

    print 'n_eff: ' + str(n_eff)
    
    return n_eff


def ttest_ind(input1,input2,axis=0,effdof=False,equal_var=True):
    """ def ttest_ind(input1,input2,axis=0,effdof=False,equal_var=True):

             compute the t-statistic and pvalue for the 2 populations.
                 Assume equal variance. This is exactly the same as
                 calling scipy.stats.ttest_ind() if effdof=False.

                 if effdof = True, calculate effective degrees of freedom
                                 in the two populations by using the lag-1
                                 autocorrelation. See calc_effectiveDOF().

         @@@as of 1/12/2015, multiple dimensions is only handled if effdof=False
              returns (tstat, pval)
    """

    # http://en.wikipedia.org/wiki/Student%27s_t-test#Independent_two-sample_t-test
    #        Equal or unequal sample sizes, equal variance
    #
    #   test statistic
    # t = (Xbar1 - Xbar2) / ( (s_X1X2) * sqrt(1/n1 + 1/n2) )
    #   pooled standard dev
    # s_X1X2 = sqrt( ( (n1-1)(s_X1)^2 + (n2-1)*(s_X2)^2 ) / (n1+n2 - 2) )
    
    if effdof==False:

        (tstat,pval) = sp.stats.ttest_ind(input1,input2,axis=axis,equal_var=equal_var)

    else:

        # @@@ perhaps change this to anomaly timeseries, and difference from 0?
        neff1 = calc_effectiveDOF(input1)
        neff2 = calc_effectiveDOF(input2)
        
        Xbar1 = input1.mean(axis=axis)
        Xbar2 = input2.mean(axis=axis)
        s1=input1.std(axis=axis)
        s2=input2.std(axis=axis)
        s12 = np.sqrt( ( (neff1 - 1)*(s1**2) + (neff2-1)*(s2**2) ) / (neff1+neff2 - 2) )

        tstat = (Xbar1 - Xbar2) /  (s12*np.sqrt( 1/neff1 + 1/neff2) )
        pval = sp.stats.t.sf(np.abs(tstat),neff1+neff2-2)*2

        # http://docs.scipy.org/doc/scipy/reference/tutorial/stats.html
        # >>> tt = (sm-m)/np.sqrt(sv/float(n))  # t-statistic for mean
        # >>> pval = stats.t.sf(np.abs(tt), n-1)*2  # two-sided pvalue = Prob(abs(t)>tt)
        
        
    return (tstat,pval)


def calc_normfit(input,verb=True):
    """ use scipy.stats.norm.fit() to calc a pdf mean and sigma

            input: 1d array of data

                   default # of return vals is 500

            returns: fittedpdf, mean, sigma, x values
    """

    # Now calc the pdf associated with the hist
    inpdf=norm.fit(input)
    mn=inpdf[0]
    sgm=inpdf[1]

    #Generate X points
    xlims = [-4*sgm+mn, 4*sgm+mn] # large limits
    xx = np.linspace(xlims[0],xlims[1],500)
    #Get Y points via Normal PDF with fitted parameters
    pdf_fitted = norm.pdf(xx,loc=mn,scale=sgm)

    if verb:
        print '===== mean ' + str(mn) + ' sigma ' + str(sgm)

    return pdf_fitted, mn, sgm, xx


def calc_kernel(input):
    """ use scipy.stats.gaussian_kde() to calc a pdf 

            input: 1d array of data
            returns pdf, x values
    """
    
    kernel = sp.stats.gaussian_kde(input)
    pdf,mn,sgm,xx= calc_normfit(input) # just to get x values
    
    return kernel(xx), xx

def regress(input1,input2):
    """ Calculate linear regression b/w two variables (x, y)
           Over axis=0!
           Will only do 1-D or 2-D
             or 1-D and 2-D as long as input2 is the one that's 2-D (2/17/16)

        Uses scipy.stats.linregress(input1,input2)

        returns mm,bb,rval,pval

    """

    if input1.ndim==1 and input2.ndim==1:
        mm, bb, rval, pval, std_err = sp.stats.linregress(input1,input2)

    else:
        mm=np.zeros(input2.shape[1])
        bb=np.zeros(input2.shape[1])
        rval=np.zeros(input2.shape[1])
        pval=np.zeros(input2.shape[1])
        if input2.ndim>1 and input1.ndim==1:
            for ii in np.arange(input2.shape[1]):
                mm[ii],bb[ii],rval[ii],pval[ii],_ = sp.stats.linregress(input1,input2[:,ii])

        else:
            for ii in np.arange(input1.shape[1]):
                mm[ii],bb[ii],rval[ii],pval[ii],_ = sp.stats.linregress(input1[:,ii],input2[:,ii])

    # How to plot: onex=np.linspace(axxlims[0],axxlims[1])    
    #              ax.plot(onex,mm*onex + bb, color='k',linewidth=2)

    
    return mm,bb,rval,pval

def trend(input):
    """  Calculate linear trend on axis=0
        use scipy.stats.linregress() or np.polyfit with order=1

        for 1-D
        returns slope, y-intercept, r-value, p-value, standard error

        for >=2-D
        returns slope, y-intercept
    """
    shape = input.shape

    xx = np.arange(0,shape[0])

    #if input.ndim==1:
    #    mm, bb, rval, pval, std_err = sp.stats.linregress(xx,input)
    #    return mm,bb,rval,pval,std_err

    if input.ndim<=2:
        mm,bb = np.polyfit(xx,input,1) 
        return mm,bb
    elif input.ndim==3:
        # polyfit will accept multiple dims I think?
        tmp = input.reshape((len(xx),shape[1]*shape[2]))
        mm,bb = np.polyfit(xx,tmp,1) 
        mm=mm.reshape((shape[1],shape[2])) # matrix of trends
        bb=bb.reshape((shape[1],shape[2]))
        return mm,bb
    else:
        # THROW EXCEPTION
        print 'trend() is not implemented for dim>3'
        return None


def trend_monthly(input):
    """ calls trend() for each month.
           input is assumed to be a monthly timeseries (axis=0).
           If timeseries partially goes into next year (ie. is not divisible by 12),
           extra months will be ignored.

        returns slope, y-intercept
    """
    shape = input.shape
    
    
    if input.ndim>1:
        initshape=(12,)+shape[1:]
       
    elif input.ndim==1:
        initshape=(12,)
        rval=np.zeros(initshape)
        pval=np.zeros(initshape)
        std_err=np.zeros(initshape)

    mm=np.zeros(initshape)
    bb=np.zeros(initshape)

    
    for moidx in np.arange(0,12): # loop through 12 months

            tmp = input[moidx::12,...]
            #if input.ndim==1:
            #    mm[moidx],bb[moidx],rval[moidx],pval[moidx],std_err[moidx]=trend(tmp)
            #else:
            #    mm[moidx,...],bb[moidx,...] = trend(tmp)
            mm[moidx,...],bb[moidx,...] = trend(tmp)
    
    #if input.ndim==1:
    #    return mm,bb,rval,pval,std_err
    #else:
    #    return mm,bb
    return mm,bb

def detrend(input, axis=0, dttype='linear'):
    """ use scipy.signal.detrend() to linearly detrend the 
           data. This function adds the mean back to data.

    """
    return sp.signal.detrend(input,axis=axis,type=dttype) + input.mean(axis=0)


""" Trying to figure out correct way to calc 95% confidence interval
https://github.com/scipy/scipy/blob/v0.13.3/scipy/stats/distributions.py#L788-L819

def interval(self, alpha, *args, **kwds):
        ""
        Confidence interval with equal areas around the median.
        Parameters
        ----------
        alpha : array_like of float
            Probability that an rv will be drawn from the returned range.
            Each value should be in the range [0, 1].
        arg1, arg2, ... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            location parameter, Default is 0.
        scale : array_like, optional
            scale parameter, Default is 1.
        Returns
        -------
        a, b : ndarray of float
            end-points of range that contain ``100 * alpha %`` of the rv's possible
            values.
        ""
        alpha = asarray(alpha)
        if any((alpha > 1) | (alpha < 0)):
            raise ValueError('alpha must be between 0 and 1 inclusive')
        q1 = (1.0-alpha)/2
        q2 = (1.0+alpha)/2
        a = self.ppf(q1, *args, **kwds)
        b = self.ppf(q2, *args, **kwds)
        return a, b


def ppf(self,q,*args,**kwds):
        ""
        Percent point function (inverse of cdf) at q of the given RV.
        Parameters
        ----------
        q : array_like
            lower tail probability
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information)
        loc : array_like, optional
            location parameter (default=0)
        scale : array_like, optional
            scale parameter (default=1)
        Returns
        -------
        x : array_like
            quantile corresponding to the lower tail probability q.
        ""
        args, loc, scale = self._parse_args(*args, **kwds)
        q, loc, scale = map(asarray,(q, loc, scale))
        args = tuple(map(asarray, args))
        cond0 = self._argcheck(*args) & (scale > 0) & (loc == loc)
        cond1 = (0 < q) & (q < 1)
        cond2 = cond0 & (q == 0)
        cond3 = cond0 & (q == 1)
        cond = cond0 & cond1
        output = valarray(shape(cond), value=self.badvalue)

        lower_bound = self.a * scale + loc
        upper_bound = self.b * scale + loc
        place(output, cond2, argsreduce(cond2, lower_bound)[0])
        place(output, cond3, argsreduce(cond3, upper_bound)[0])

        if any(cond):  # call only if at least 1 entry
            goodargs = argsreduce(cond, *((q,)+args+(scale,loc)))
            scale, loc, goodargs = goodargs[-2], goodargs[-1], goodargs[:-2]
            place(output, cond, self._ppf(*goodargs) * scale + loc)
        if output.ndim == 0:
            return output[()]
        return output



## autocorrelation

http://stackoverflow.com/questions/12269834/is-there-any-numpy-autocorrellation-function-with-standardized-output

In[1]: n = 1000
In[2]: x = randn(n)
In[3]: xc = correlate(x, x, mode='full')
In[4]: xc /= xc[xc.argmax()]
In[5]: xchalf = xc[xc.size / 2:]
In[6]: xchalf_max = xchalf.max()
In[7]: print xchalf_max
Out[1]: 1.0

"""


    
