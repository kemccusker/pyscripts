import xarray as xr
import cccmautils as cutl
import numpy as np

''' This class contains utility functions to operate on xarray Dataset or DataArray
'''




def get_areaweights_ds(ds,model='CanESM2'):
    
    awgts=cutl.calc_areawgts(lat=ds.lat.data,lon=ds.lon.data,model='CanESM2')
    area_weights = xr.DataArray(awgts, coords=[ds.lat,ds.lon],name='area_weight')

    # why is the assertion failed? 
    #  The sum is not within tolerance level for some reason,
    #  even thought weights are produced by dividing by the sum
    #np.testing.assert_allclose(area_weights.sum(),1)
    
    return area_weights

def calc_areaweight_mean_ds(ds,dim=None, mask=None, model='CanESM2'):

    """ Calculates area-weighted mean of the input Dataset (or DataArray?? test)
        If the input ds is a subregion, this should still work. (tested with polcap60)
    
        mask: ['land' | 'ocean'] (mask='land' means mask OUT land) 
              Lakes are not included in land or ocean.
              @@this is hard-coded to CanESM2. @

        return area-weighted values in a Dataset (or DataArray?? test)
    """
    aweights = get_areaweights_ds(ds)

    if mask!=None:
        import constants as con
        
        if model=='CanESM2':
            lmask = xr.DataArray(con.get_t63landmask()[:,:-1],
                                     dims=('lat','lon'),
                                     coords=[ds.lat,ds.lon],name='lmask')
            #print lmask
        else:
            raise Exception('No landmask defined for model ' + model)
        
        if mask=='land': mval = -1;
        elif mask=='ocean': mval=0

        #awgtsub = aweights.values
        #awgtsub[lmask!=mval] = np.nan
        #awgtsub = awgtsub / np.float(np.nansum(awgtsub))# normalize weights
        
        #aweights.values = awgtsub

        #print ds.where(lmask!=mval).values.shape
        
        #ds.where(lmask!=mval)
        #dssub = ds.values
        #dssub[lmask!=mval] = np.nan
        #ds.values = dssub
        
    if dim!=None:
        if mask!=None:
            aw=aweights.where(lmask!=mval) / aweights.where(lmask!=mval).sum()
            return  (aw * ds.where(lmask!=mval)).sum(dim=dim)
        else:
            return (aweights * ds).sum(dim=dim)
    else:
        if mask!=None:
            aw=aweights.where(lmask!=mval) / aweights.where(lmask!=mval).sum()
            return (aw * ds.where(lmask!=mval)).sum()
        else:
            return (aweights * ds).sum()


    
# weighted seasonal / annual average
# code from (modified): http://xarray.pydata.org/en/stable/examples/monthly-means.html
def get_dpm(time):#, calendar='standard'):
    """
    return a array of days per month corresponding to the months provided in `months`
    """
    month_length = np.zeros(len(time), dtype=np.int)

    cal_days = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]#dpm[calendar]

    if len(month_length) == 12:
        #just a climo
        month_length = cal_days[1:]
    else:
        for i, (month, year) in enumerate(zip(time.month, time.year)):
            month_length[i] = cal_days[month]
            #if leap_year(year, calendar=calendar):
            #    month_length[i] += 1
    return month_length

def calc_seasonal_mean_ds(ds):
    
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = xr.DataArray(get_dpm(ds.time.to_index()),
                                coords=[ds.time], name='month_length')

    # Calculate the month weights by grouping by 'time.season'.
    # Conversion to float type ('astype(float)') only necessary for Python 2.x
    mweights = month_length.groupby('time.season') / month_length.astype(float).groupby('time.season').sum()

    # Test that the sum of the weights for each season is 1.0
    np.testing.assert_allclose(mweights.groupby('time.season').sum().values, np.ones(4))

    # Calculate the weighted average
    return (ds * mweights).groupby('time.season').sum(dim='time')
    

def calc_annual_mean_ds(ds):
    """ This only works on monthly or greater time resolution 
        (e.g. won't work on monthly climo b/c 'time' is no longer a dim) -->
           Added functionality to handle both...
    """
    # Make a DataArray with the number of days in each month, size = len(time)
    try:        
        month_length = xr.DataArray(get_dpm(ds.time.to_index()),
                                    coords=[ds.time], name='month_length')
        # Calculate the month weights by grouping by 'time.year'.
        # Conversion to float type ('astype(float)') only necessary for Python 2.x
        mweights = month_length.groupby('time.year') / month_length.astype(float).groupby('time.year').sum()

        nyear = len(ds.groupby('time.year').mean().year)
    
        # Test that the sum of the weights for each year is 1.0
        np.testing.assert_allclose(mweights.groupby('time.year').sum().values, 
                                       np.ones(nyear))

        # Calculate the weighted average
        ann = (ds * mweights).groupby('time.year').sum(dim='time')

    except AttributeError:
        # no attribute 'time' --> already a climo (grouped by month)
        month_length = xr.DataArray(get_dpm(ds.month.to_index()),
                                    coords=[ds.month], name='month_length')
        mweights = month_length / month_length.astype(float).sum()
        np.testing.assert_approx_equal(mweights.sum().values, 1)
        ann = (ds * mweights).sum(dim='month')
        
    
    return ann
    

def calc_monthly_mean_ds(ds):
    
    return ds.groupby('time.month').mean(dim='time')


def calc_totseaicearea(ds,isarea=False): #fld,lat,lon,isarea=False,model='CanESM2'):
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

        
    awgts=cutl.calc_areawgts(lat=ds.lat.data,lon=ds.lon.data,model='CanESM2')
    area_weights = xr.DataArray(awgts, coords=[ds.lat,ds.lon],name='area_weight')
    
    # Changed to ma.sum() on 12/10/14
    nh = ma.sum(ma.sum(sia[...,lat>0,:],axis=-1),axis=-1)
    sh = ma.sum(ma.sum(sia[...,lat<0,:],axis=-1),axis=-1)

    return nh,sh

