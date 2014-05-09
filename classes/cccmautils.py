"""
cccmautils.py
Put my utility functions in this module.
   e.g. area-weighted global mean, etc

   2/11/2014

"""

import numpy as np
import constants as con

con = reload(con)

def find_nearest(array,value):
    """ find_nearest(array,value): returns index of array
        element with nearest value to value
    """
    return (np.abs(array-value)).argmin()


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
   

def polar_mean_areawgted3d(input,lat,lon,latlim=60,hem='nh',cellareas=None):
    """ Pass in cellareas if you want some of the cells masked, e.g. if masking ocean or land.
        This ONLY works if the mask is the same for each time in timeseries.
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

    totalarea=np.sum(np.sum(polcellareas))

    nt = input.shape[0]
    wgts = polcellareas/totalarea
    
    pm = np.zeros(nt)    
    for tidx in range(0,nt):
        pm[tidx] = np.average(polinput[tidx,:,:],weights=wgts)
        
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
