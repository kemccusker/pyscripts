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


def calc_seaicearea(input,lat,lon):
    """ Calculate sea ice area from sea ice concentration
        input: 2D or greater array of sea ice concentration.
               If > 2D, time must be first dimension
        lat: array of coordinate lats
        lon: array of coordinate lons

        returns: MASKED array of sea ice area with same
                 shape as input. The array is masked such
                 that land is masked out
        @@@ needs testing
    """

    ishape = input.shape
    ndims=len(ishape)

    if ndims>2: # first dim is tim
        areas = calc_cellareas(lat,lon,repeat=ishape)
        lmask = con.get_t63landmask(repeat=ishape) # @@ note assuming T63 here...
    else:
        areas = calc_cellareas(lat,lon)
        lmask = con.get_t63landmask()

    sia = input*areas
    sia = ma.masked_where(lmask==-1,sia) # mask where land
    
    return sia

def calc_totseaicearea(input,lat,lon):
    """ calculate total sea ice area
           returns nh,sh
           input is expected to be time x lat x lon

        @@ needs testing 9/9/2014
    """

    sia = calc_seaicearea(input,lat,lon)
    nh = np.sum(np.sum(sia[:,lat>0,:],axis=2),axis=1)
    sh = np.sum(np.sum(sia[:,lat<0,:],axis=2),axis=1)

    return nh,sh

def calc_meanseaicethick(input,lat,lon):
    """ calc_meanseaicethick(input,lat,lon)
                 input: 2D or greater array of sea ice thickness
                 
    """
    
    print "Not finished!"
    
def global_mean_areawgted3d(input, lat, lon):

    earthrad = con.get_earthrad()
    totalarea = 4*np.pi*earthrad**2

    cellareas = calc_cellareas(lat,lon)

    nt = input.shape[0]
    wgts = cellareas/totalarea
    
    gm = np.zeros(nt)    
    for tidx in range(0,nt):
        gm[tidx] = np.average(input[tidx,:,:],weights=wgts)
        
    return gm
   

def polar_mean_areawgted3d(input,lat,lon,latlim=60,hem='nh',cellareas=None,includenan=False):
    """ Pass in cellareas if you want some of the cells masked, e.g. if masking ocean or land.
        This ONLY works if the mask is the same for each time in timeseries.
            includenan: whether or not to consider NaN values. if True, a mean of a NaN is NaN.
                        Default is False
    """

    if cellareas == None:
        cellareas = calc_cellareas(lat,lon)
        # else cellareas are provided
        
    # @@ what does this do if just an ocean field, etc?

    if hem=='nh':
        polcellareas = cellareas[lat>=latlim,:]
        polinput = input[:,lat>=latlim,:]
    else:
        polcellareas = cellareas[lat<=latlim,:]
        polinput = input[:,lat<=latlim,:]

    totalarea=float(np.sum(np.sum(polcellareas)))

    nt = input.shape[0]
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

def annualize_monthlyts(input, includenan=0):
    """ It is expected that the input variable is 3 dimension, with time being first
        But on 5/6/2014 I tried to generalize to any dimension...needs to be tested.
    """

    # time must be first dimension
    # figure out nyears (truncate if only some months in last year)
    dims = input.shape
    nt = input.shape[0] # number of times
    nyrs = nt/12        # I think this will automatically return an integer

    #ndim2 = input.shape[1]
    #ndim3 = input.shape[2]
    
    wgts = con.get_monweights()

    if len(dims)>1:
        diml = list(dims[1:]) # this one is for wgts array
        diml.append(1) # insert size nyrs into first dimension
        diml = tuple(diml)
        
        diml2 = list(dims[1:]) # this one is for initializing the timeseries array
        diml2.insert(0,nyrs) # insert size nyrs into first dimension

        wgts = np.tile(wgts,diml) # make a weights array for each grid pt
        # put the 3rd dim (2) in first spot, 1st dim (0) in second spot, 2nd dim (1) in third spot
        # this effectively shifts the dimensions back to time,dim2,dim3
        wgts = np.transpose(wgts,(2,0,1))
    else:
        diml2 = nyrs,

    #yidx=0
    counter=0
    annts = np.zeros(diml2)
    
    # Loop through 12-month chunks
    for yidx in range(0, nyrs):
        annts[yidx,...] = np.average(input[counter:counter+12,...],axis=0,weights=wgts)
        counter = counter+12

    return annts

def seasonalize_monthlyts(input,season=None,includenan=0,mo=0,climo=0):

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
        return annualize_monthlyts(input,includenan)
    elif season=='DJF':
        # do weighted average for requested season
        
        if climo==1:
            subwgts = wgts[[0,1,11]]
            seasts = np.average(input[[0,1,11],...],weights=subwgts,axis=0)
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
            seasts = np.average(input[[0,10,11],...],weights=subwgts,axis=0)
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
    else:
        print "Season not supported!"
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
        subwgts = np.transpose(subwgts,(2,0,1))
    
    else:
        diml2 = nyrs,


    seasts = np.zeros(diml2)
    #print subwgts.shape
    
    for yridx in range(0,nyrs):
        subsamp = range(start+yridx*12,start+yridx*12+incr)
        seasts[yridx,...] = np.average(input[subsamp,...],weights=subwgts,axis=0)        
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
    lat_bounds = np.insert(lat_bounds,0,-90)
    lat_bounds = np.append(lat_bounds,90)

    del_phi = np.sin(lat_bounds[1:]*deg2rad) - np.sin(lat_bounds[0:-1]*deg2rad)
    
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
                     
def get_cellwgts(lat,lon,repeat=None):

    """ basically, repeat should be the shape of the field that
    #   we want to weight with cell weights. Assume the last
    #   two numbers are lat and lon sizes
    #   If None, no repeating necessary
    """
    
    earthrad = con.get_earthrad()
    totalarea = 4*np.pi*earthrad**2
    wgts = calc_cellareas(lat,lon)/totalarea
    
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

def calc_regmean(fld,lat,lon,region,limsdict=None):
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

                 Returns: Regional mean (or series of regional means with length ndim1)
    """

    
    fldm,regmask = mask_region(fld,lat,lon,region,limsdict)

    # calculate area-weights
    areas = calc_cellareas(lat,lon)
    if regmask.ndim>2:
        areasm = ma.masked_where(np.squeeze(regmask[0,...]),areas)
    else:
        areasm = ma.masked_where(regmask,areas)
    weightsm = areasm / np.sum(np.sum(areasm,axis=1),axis=0) # weights masked

    if fld.ndim>2:
        ndim1 = fld.shape[0]
        weightsmt = np.tile(weightsm,(ndim1,1,1)) # weights masked tiled
    else:
        weightsmt = weightsm

    tmp = ma.masked_where(regmask,fldm)
    tmpreg = np.sum(np.sum(tmp*weightsmt,axis=2),axis=1)
    fldreg = tmpreg # should be ndim1 of regional mean (or just one regional mean)

    return fldreg


def calc_monthlystd(input):
    """ calc_monthlystd(input):
               Given a monthly timeseries (input) of size time [x lat x lon],
               return a std deviation for each month such that the return dimension
               is 12 x [shape[1:]].

        Note, can also get this info from climatologize(), which returns
          (climo,stddev)
               
    """
    dims = input.shape
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

        stddev[mo] = np.std(input[moidx::12,...],axis=0)

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
