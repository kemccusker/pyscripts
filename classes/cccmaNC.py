""" cccmaNC.py
    3/3/2014
    functions to open and read data from Netcdf, with
       option to perform basic operations (some use CDO bindings)
       (seasonal mean, zonal mean)
"""
import numpy as np # for array handling
import scipy as sp # scientific python
from netCDF4 import Dataset, num2date
import cccmautils as cutl
import cdo as cdo #; cdo = cdo.Cdo()
import os
import platform as platform
import datetime as datetime
import matplotlib.dates as dates



def openNC(filename):

    ncfile = Dataset(filename,'r')
    return ncfile


def get_NCtimeunits(filename,field='time'):

    ncfile = openNC(filename)
    time_nc=ncfile.variables[field]

    return time_nc.units

def get_NCtimecal(filename,field='time'):

    ncfile = openNC(filename)
    time_nc=ncfile.variables[field]

    if 'calendar' not in time_nc.ncattrs():
        print 'NO CALENDAR in time variable @@@. Default to 365_day'
        return '365_day'
    else:
        return time_nc.calendar

def NCtime2date(thetimes,timeunit,timecal='365_day'):

    thedates = num2date(thetimes,timeunit,timecal)

    # @@@ need to convert each of these to np.datetime64 if want to convert to pandas DatetimeIndex
    
    return thedates

def get_NCdates(filename,field='time',timesel=None):
    """ will get the field ('time' is default) and convert to datetimes

    """

    thetimes = getNCvar(filename,field=field,timesel=timesel)
    theunit = get_NCtimeunits(filename,field=field)
    thecal = get_NCtimecal(filename,field=field)

    thedates = NCtime2date(thetimes,theunit,thecal)

    return thedates


def getNCvar(filename,field,timesel=None,levsel=None,monsel=None,seas=None,calc=None,remlon=1,sqz=True,att=False):
    """ gets a variable from netcdf file.
        Time is assumed to be the 1st dimension, Lon is assumed to be the last.
        If any calculations are requested to be performed on the data, the user
        needs to make sure that the requested operations can be performed (b/c
        some of the other functions only handle certain # of dims, etc. My bad.)

        filename: full path to file
        field: NC variable to read in
        timesel: comma-delim string of date range in 'YYYY-MM-DD' fmt to select (a la CDO)
        levsel: select level in Pa (e.g. 50000 for 500hPa)
        monsel: select month from timeseries
        seas: seasonally (annually) average or return climatology {climo|ANN|DJF|JJA|NDJ|MAM|SON}
        calc: zm (zonal mean),
        remlon: removes extra wrap-around longitude for zonal mean.
                default is 1, remove it
        sqz:  squeeze the data if 'getting all data'. Default True.
                   Trying to avoid situation where need singular dims and squeeze them out
                   (e.g. MOC variable in CCSM4)
        att: if True, include the full netcdf variable (including attributes). Most useful for time.
        
        returns fld
        """
# IF ON MAC: CDO bindings don't work yet, use old function 3/25/2014 #########

    plat = platform.system()

    if plat == 'Darwin':  # means I'm on my mac
        # Call old func
        if calc != None:
            print '@@ not sure calc will work on mac. calc=' + calc

        if levsel!=None:
            plev= np.array([100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, 12500,
                   15000, 17500, 20000, 22500, 25000, 30000, 35000, 40000, 45000, 50000,
                   55000, 60000, 65000, 70000, 75000, 77500, 80000, 82500, 85000, 87500,
                   90000, 92500, 95000, 97500, 100000])

            level=cutl.find_nearest(plev,levsel)
        else:
            level=None        

        if timesel == '0002-01-01,0061-12-31':
            print 'hard-coded skipping of first year of 61-yr chunk @@'
            fld = getNCvar_old(filename,field,seas=seas, monsel=monsel,timechunk=(12,),level=level,calc=calc,sqz=sqz)
        else:
            # if timesel=='0002-01-01,0121-12-31' then just don't set timechunk because 
            #     files on the mac are already selected to skip first year, and they reside
            #     in the 'timsel' subdirectory. Check for that?
            if timesel=='0002-01-01,0121-12-31':
                if 'timsel/' not in filename:
                    print 'On mac, use files in timsel/ subdirectory! @@ NEEDS TESTING'

            fld = getNCvar_old(filename,field,seas=seas,monsel=monsel,level=level,calc=calc,sqz=sqz,timesel=timesel) # doesn't work with all arguments yet @@

        
        return fld

    else:  # on linux workstation in Vic

        ncfile = openNC(filename)
        ndims = len(ncfile.dimensions)
        ncvar = ncfile.variables[field]
        #print ncvar # @@@@@

        
        #### READ VARIABLE FROM NC FILE ########
        if timesel == None and calc == None:

            if levsel !=None:
                if monsel != None:
                    fld = np.squeeze(cdo.sellevel(levsel,input = cdo.selmon(monsel,input = filename),returnMaArray = field))
                else:
                    fld = np.squeeze(cdo.sellevel(levsel,input = filename, returnMaArray = field))
                os.system('rm -rf /tmp/cdoPy*')
            else:
                if monsel != None:
                    #print 'timesel==None and calc==None and monsel !=None'
                    fld = np.squeeze(cdo.selmon(monsel,input = filename, returnMaArray = field))
                    #print fld.shape
                    os.system('rm -rf /tmp/cdoPy*')
                else: # get everything
                    if sqz:
                        print field + ': squeezing data upon read all' # @@@
                        # for most situations, this is what we want. @@@@
                        fld=np.squeeze(ncfile.variables[field][...])
                    else:
                        fld = ncfile.variables[field][...]

        elif timesel != None and calc == 'zm':
            # have to remove the lon before zonal mean, which means have to separate the
            # select dates and zm. thus can't use CDO for zm (unless can pass it data instead of a file?)

            #fld = np.squeeze(cdo.zonmean( input = cdo.seldate(timesel,input = filename), returnMaArray = field))
            print 'assuming T42(63) 64x128 resolution for zonal mean'
            if levsel != None:
                if monsel != None:
                    fld = np.squeeze(cdo.seldate(timesel,input = cdo.zonmean( input = 
                                                 cdo.selindexbox(1,128,1,64,input =
                                                                 cdo.sellevel(levsel,input =
                                                                              cdo.selmon(monsel, input = filename)))),
                                                 returnMaArray = field))# @@@@
                else:
                    fld = np.squeeze(cdo.seldate(timesel,input = cdo.zonmean( input = 
                                                 cdo.selindexbox(1,128,1,64,input =
                                                                 cdo.sellevel(levsel,input = filename))),
                                                 returnMaArray = field))
            else:
                if monsel != None:
                    fld = np.squeeze(cdo.seldate(timesel,input =
                                                 cdo.zonmean( input =
                                                              cdo.selindexbox(1,128,1,64,input =
                                                                              cdo.selmon(monsel,input = filename))),
                                                              returnMaArray = field))
                else:
                    fld = np.squeeze(cdo.seldate(timesel,input =
                                                 cdo.zonmean(input =
                                                             cdo.selindexbox(1,128,1,64,input = filename)),
                                      returnMaArray = field))
            os.system('rm -rf /tmp/cdoPy*')

            ## if remlon:
            ##     # remove extra lon
            ##     fld = np.squeeze(fld[...,0:-1])

            ## lastdimidx = ndims-1
            ## fld = np.mean(fld,lastdimidx)  
            

        elif timesel != None and calc != None:
            if levsel != None and monsel == None:
                fld = np.squeeze(cdo.seldate(timesel,input =
                                             cdo.sellevel(levsel,input = filename),
                                             returnMaArray = field))
            elif levsel != None and monsel != None:
                fld = np.squeeze(
                    cdo.seldate(timesel,input =
                                cdo.sellevel(levsel,input =
                                             cdo.selmon(monsel,input = filename)),
                                returnMaArray = field))
            elif levsel == None and monsel != None:
                fld = np.squeeze(cdo.seldate(timesel,input =
                                  cdo.selmon(monsel,input = filename),
                                  returnMaArray = field))
            else: # levsel and monsel are both None
                fld = np.squeeze(cdo.seldate(timesel,input = filename, returnMaArray = field))
            os.system('rm -rf /tmp/cdoPy*')
            print "only calc='zm' is implemented now. Returning only selected date range/level/month."

        elif timesel != None:
            if levsel != None and monsel == None:
                fld = np.squeeze(cdo.seldate(timesel,input = cdo.sellevel(levsel,input = filename),returnMaArray = field))
            elif levsel != None and monsel != None:
                fld = np.squeeze(
                    cdo.seldate(timesel,input =
                                cdo.sellevel(levsel,input =
                                             cdo.selmon(monsel,input = filename)),
                                returnMaArray = field))
            elif levsel == None and monsel != None:
                fld = np.squeeze(cdo.seldate(timesel,input =
                                  cdo.selmon(monsel,input = filename),
                                  returnMaArray = field))
            else: # levsel and monsel are both None
                fld = np.squeeze(cdo.seldate(timesel,input = filename, returnMaArray = field))
                
            os.system('rm -rf /tmp/cdoPy*')
            
        elif calc == 'zm': # and timesel must be None
            print 'assuming T42(63) 64x128 resolution for zonal mean'
            
            if levsel != None and monsel == None:
                fld = np.squeeze(cdo.sellevel(levsel,input =
                                              cdo.zonmean(input =
                                                          cdo.selindexbox(1,128,1,64,input =
                                                                          filename)),
                                              returnMaArray = field))

            elif levsel != None and monsel != None:
                fld = np.squeeze(
                    cdo.sellevel(levsel,input =
                                 cdo.zonmean(input =
                                             cdo.selindexbox(1,128,1,64,input =
                                                             cdo.selmon(monsel,input =
                                                                        filename))),
                                 returnMaArray = field))

            elif levsel == None and monsel != None:
                fld = np.squeeze(cdo.zonmean(input =
                                             cdo.selindexbox(1,128,1,64,input =
                                                             cdo.selmon(monsel,input =
                                                                        filename)),
                                                 returnMaArray = field))
               
            else: # get all data
                fld = np.squeeze(cdo.zonmean(input =
                                             cdo.selindexbox(1,128,1,64,input =
                                                             filename),
                                             returnMaArray = field))
                
                #print '@@ getting memory errors here...try using CDO to select appropriate lons for the zm calc'
                #fld = ncfile.variables[field][...] # have to get field before removing lon

            os.system('rm -rf /tmp/cdoPy*')
            ## if remlon:
            ##     # remove extra lon
            ##     if ndims==4:
            ##         fld = np.squeeze(fld[:,:,:,0:-1])
            ##     elif ndims==3:
            ##         fld = np.squeeze(fld[:,:,0:-1])
            ##     else: # shouldn't really get here, not expecting 2D (time x lon?)
            ##         fld = np.squeeze(fld[:,0:-1])
            ## lastdimidx = ndims-1
            ## fld = np.mean(fld,lastdimidx)  
            

        else:
            print "huh? timesel and calc combo doesn't make sense"


        

        ####### TIME AVERAGE the VARIABLE ##########
        # fld has to be 3d by the time it is passed to func
        #  (time,lev,lat) or (time,lat,lon)
        if seas != None:
            #print 'getNCvar(): seas!=None: fld.shape: ' + str(fld.shape) # @@@
            
            ## if fld.ndim != 3:
            ##     ## if 1 in fld.shape:
            ##     ##     fld=fld.squeeze() # attempting to deal with spurious dims of 1 @@@
            ##     ##     if fld.ndim != 3:
            ##     ##         print 'data must be 3 dimensional to seasonalize()'
            ##     ##         return
            ##     ## else:
            ##     print 'data must be 3 dimensional to seasonalize()'
            ##     return
            
            if monsel != None:
                print "Can't do seasonal average when monsel != None"
                return
            elif seas == 'climo':
                fld,stddev = cutl.climatologize(fld)
            elif type(seas) == int: # @@ does this work?
                #elif seas not in ('ANN','DJF','JJA','MAM','SON','NDJ'):
                # means seas is an int value for a month
                
                #fld = cutl.seasonalize_monthlyts(fld,mo=seas)
                fld = cutl.seasonalize(fld,mo=seas)
            else:
                #print 'seasonalizing'
                #fld = cutl.seasonalize_monthlyts(fld,season=seas)
                fld = cutl.seasonalize(fld,season=seas)
                #print fld.shape


        # Apply any scaling and offsetting needed:
        try:
            var_offset = ncvar.add_offset
        except:
            var_offset = 0
            
        try:
            var_scale = ncvar.scale_factor
            print 'var_scale ' + str(var_scale)
        except:
            var_scale = 1
            
        
        fld = fld*var_scale + var_offset

        ncfile.close()
        return fld



def getNCvar_old(filename,field,timechunk=None,monsel=None,level=None,seas=None,calc=None,remlon=1,sqz=True,timesel=None):
    """ gets a variable from netcdf file.
        Time is assumed to be the 1st dimension, Lon is assumed to be the last.
        If any calculations are requested to be performed on the data, the user
        needs to make sure that the requested operations can be performed (b/c
        some of the other functions only handle certain # of dims, etc. My bad.)

        filename: full path to file
        field: NC variable to read in
        timechunk: tuple of start,stop
        timesel: same format as CDO. will be parsed into startyear, startmon, startday and endyr,enmon,enday
        monsel: index of month to choose (1=all Jans, 2=all Febs, etc. really meant for CDO bindings)
        level: index of plev to select
        seas: seasonally (annually) average {ANN|DJF|JJA|NDJ|MAM|SON}
        calc: zm (zonal mean),
        remlon: removes extra wrap-around longitude for zonal mean.
                default is 1, remove it
        sqz:  squeeze the data if 'getting all data'. Default True.
                   Trying to avoid situation where need singular dims and squeeze them out
                   (e.g. MOC variable in CCSM4)
                   
        returns fld
        """
#    chunk = None
    
    ncfile = openNC(filename)
    ndims = len(ncfile.dimensions)
    #print 'getNCvar_old()'
    
    #### READ VARIABLE FROM NC FILE ########
    if (timechunk == None) and (timesel == None):
        if level != None:
            fld = ncfile.variables[field][:,level,...]
        else:
            fld = ncfile.variables[field][...]

        # @@@shit, caused problems with MOC.
        if sqz:
            fld=fld.squeeze() # remove spurious dimensions of 1
    else:
        if timesel == None:
            print 'timesel is None' # @@@
        else:
            #print 'timesel ' + timesel

            # parse timesel:
            # year-month-day
            (timselst,timselen)=timesel.split(',')       
            (styear,stmon,stday) = timselst.split('-')
            (enyear,enmon,enday) = timselen.split('-')
            #print styear + ' to ' + enyear # @@@@@@

            # ############################
            # from Joe's util.py (originally Phil Austin cookbook)
            #
            # first convert to netcdftime datetime objects
            #
            time_nc=ncfile.variables['time']
            the_times=time_nc[...]

            #try:

            the_dates=num2date(the_times,time_nc.units,get_NCtimecal(filename))
            #print the_dates # @@@@@@
            #
            # netCDF4 bug(?) means that netcdftime objects can't be compared/sorted
            # so convert to python datetime objects
            #

            py_dates=[datetime.date(*item.timetuple()[:3]) for item in the_dates] # This only deals in Year/Month/Day.
            py_dates=np.array(py_dates)

            # Convert start and stop dates to python datetimes
            # Presently assume monthly data!
            date1 = datetime.date( int(styear), int(stmon), 1 ) #Starts Jan 1
            date2 = datetime.date( int(enyear), int(enmon), 31 )  #Ends Dec 31  
            #print str(date1) + ' and ' + str(date2) # @@@@@@

            #
            # 
            if styear is None:
                start_index=0
            else:
                start_index=find_index(py_dates,date1)[0]

            if enyear is not None:
                stop_index=find_index(py_dates,date2)[0]
                if stop_index == 0:
                    # The stop index is 0 since you asked for a date beyond that in the 
                    # file. So print a warning and set it to the last date in the file
                    print 'You are asking for a date beyond that in the file.'
                    print 'I will take the last time index in the file and move on.'
                    stop_index = len(the_times)

            else:
                stop_index=None

            time_slice=slice(start_index,stop_index+1) # b/c it's monthly data, want inclusize last month@@
            the_times=py_dates[time_slice] # may not need this...@@
            #print 'the_times: ' + str(the_times) # @@@@@@@@@@@@@
            #  @@@@ end testing new time functionality

        
# ############################
        #if (timesel!= None):# and (len(timechunk)==1): # start time until end

#            firsttime = ncfile.variables['time'][timechunk[0],...] #@@@
#            firstdate = datetime.date(1850, 1, 1) + datetime.timedelta(int(firsttime))
#            print firsttime # @@@
#            print firstdate # @@@

#            if level != None:
#                fld = ncfile.variables[field][time_slice,level,...]
#            else:
#                fld = ncfile.variables[field][time_slice,...]
#        else:
        if level != None:
            fld = ncfile.variables[field][time_slice,level,...]
        else:
            fld = ncfile.variables[field][time_slice,...]

        """ old functionality...@@
        if (timechunk!= None) and (len(timechunk)==1): # start time until end
            print timechunk # @@
            firsttime = ncfile.variables['time'][timechunk[0],...] #@@@
            firstdate = datetime.date(1850, 1, 1) + datetime.timedelta(int(firsttime))
            print firsttime # @@@
            print firstdate # @@@

            if level != None:
                fld = ncfile.variables[field][timechunk[0]:,level,...]
            else:
                fld = ncfile.variables[field][timechunk[0]:,...]
        else:
            if level != None:
                fld = ncfile.variables[field][timechunk[0]:timechunk[1],level,...]
            else:
                fld = ncfile.variables[field][timechunk[0]:timechunk[1],...]
        """
        # @@@shit, caused problems with MOC.
        if sqz:
            fld=fld.squeeze() # remove spurious dimensions of 1

        ## #print 'chunking time. ndims= ' + str(ndims) + ' styr,enyr: ' + str(timechunk[0]) + ',' + str(timechunk[1])
        ## if ndims==4: # must be a better way!@@
        ##     if len(timechunk)==1: # start time until end
        ##         fld = ncfile.variables[field][timechunk[0]:,:,:,:]
        ##     else:
        ##         fld = ncfile.variables[field][timechunk[0]:timechunk[1],:,:,:]
        ## elif ndims==3:
        ##     if len(timechunk)==1: # start time until end
        ##         fld = ncfile.variables[field][timechunk[0]:,:,:]
        ##     else:
        ##         fld = ncfile.variables[field][timechunk[0]:timechunk[1],:,:]
        ## else: # shouldn't get here
        ##     if len(timechunk)==1: # start time until end
        ##         fld = ncfile.variables[field][timechunk[0]:,:]
        ##     else:
        ##         fld = ncfile.variables[field][timechunk[0]:timechunk[1],:]
        #print fld.shape

    ncfile.close()
    
    ###### CALCULATIONS on the VARIABLE ########
    if calc == 'zm':
        #print 'zonal mean'
        # lon is assumed to be last dimension
        if remlon:
            # remove extra lon
            fld = np.squeeze(fld[...,0:-1]) # @@
            ## if ndims==4:
            ##     fld = np.squeeze(fld[:,:,:,0:-1])
            ## elif ndims==3:
            ##     fld = np.squeeze(fld[:,:,0:-1])
            ## else: # shouldn't really get here
            ##     fld = np.squeeze(fld[:,0:-1])

        print 'removed lon. ' + str(fld.shape)
        lastdimidx = ndims-1
        fld = np.mean(fld,lastdimidx)
        #print fld.shape
    elif calc == None:
        pass
    else:
        print 'only zm is implemented so far'
        return

    ####### TIME AVERAGE the VARIABLE ##########
    # fld has to be 3d by the time it is passed to func
    #  (time,lev,lat) or (time,lat,lon)
    if seas != None:
        if fld.ndim != 3:
            ## if 1 in fld.shape:
            ##     fld=fld.squeeze() # attempting to deal with spurious dims of 1 @@@
            ##     if fld.ndim != 3:
            ##         print 'data must be 3 dimensional to seasonalize()'
            ##         return
            ## else:
            print 'data must be 3 dimensional to seasonalize()'
            return
        if monsel != None:
            print "Can't do seasonal average when monsel != None"
            return
        else:
            #print 'seasonalizing ' + seas # @@@
            #fld = cutl.seasonalize_monthlyts(fld,seas)
            fld = cutl.seasonalize(fld,seas)
            #print fld.shape

    if monsel != None:
        # send back only selected month
        # monsel is 1-based! ie 1=Jan, 2=Feb. But the index is zero-based
        #   so shift by 1:
        monsel = monsel-1
        fld = fld[monsel::12,...]

        
    return fld


def find_index(vec_vals,target):
    """
    from Joe util.py, originally Phil Austin
    added 2/4/2015

    returns the first index of vec_vals that contains the value
    closest to target.

    Parameters
    ----------

    vec_vals: list or 1-d array
    target:   list 1-d array or scalar


    Returns
    -------

    list of len(target) containing the index idx such that
    vec_vals[idx] is closest to each item in target

    Example
    -------

    >>> lons=[110,115,120,125,130,135,140]
    >>> find_index(lons,[115.4,134.9])
    [1, 5]

    """
    target=np.atleast_1d(target)  #turn scalar into iterable, no-op if already array
    vec_vals=np.array(vec_vals)  #turn list into ndarray or no-op if already array
    index_list=[]
    for item in target:
        first_index=np.argmin(np.abs(vec_vals - item))
        index_list.append(first_index)
    
    return index_list 

if __name__ != '__main__':
 
    plat = platform.system()
    if plat == 'Darwin':  # means I'm on my mac
        pass
    else:
        cdo = cdo.Cdo()


if __name__ == "__main__":
    # unit test
    import platform as platform

    plat = platform.system()

    if plat == 'Darwin':  # means I'm on my mac
        basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
        subdir = '/'
    else:  # on linux workstation in Vic
        basepath = '/home/rkm/work/DATA/' + model + '/'
        subdir = '/ts/'

    casename = 'kemctl1'
    timstr = '001-111'

    #### test 2D var with time
    field = 'st'
    ncfield = 'ST'
    fname = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'

    fldall = getNCvar(fname,ncfield)
    print 'got all ' + ncfield + ' data'
    print fldall.shape

# @@@ FIGURE OUT HOW TO SELECT BASED ON DATE/TIME
##     ## styrp = 2             # skip year 1
##     ## enyrp = 61
##     ## # @@ wrong: chunk = ( (styrp-1)*12:(enyrp*12+1),:,: )
##     ## fldchunk = getNCvar(fname,ncfield,chunk=chunk)
##     ## print 'got chunked ' + ncfield + ' data'
##     ## print fldchunk.shape


## # try using actual netcdf times
## ## double time(time) ;
## ##                 time:standard_name = "time" ;
## ##                 time:long_name = "time" ;
## ##                 time:units = "days since 1850-01-01 00:00:00" ;
## ##                 time:calendar = "365_day" ;
## ## First 4 times: -652604.499976852, -652575, -652545.499976852, -652515, 
##     from netCDF4 import num2date, date2num
##     import datetime as datetime

##     #times = getNCvar(fname,'time')
## tunits="days since 1850-01-01 00:00:00"
## tcal = "365_day"
## dates=num2date(times,tunits,tcal)
##  gives: array([   1-01-16 12:00:01,    1-02-15 00:00:00,    1-03-16 12:00:01,
##           1-04-16 00:00:00,    1-05-16 12:00:01,    1-06-16 00:00:00,
##           1-07-16 12:00:01,    1-08-16 12:00:01,    1-09-16 00:00:00,
##           1-10-16 12:00:01,    1-11-16 00:00:00,    1-12-16 12:00:01,
##           2-01-16 12:00:01,    2-02-15 00:00:00,    2-03-16 12:00:01, ......
## ETC
## http://unidata.github.io/netcdf4-python/


    fldzm = getNCvar(fname,ncfield,calc='zm')
    print 'got zonal mean ' + ncfield + ' data'
    print fldzm.shape

    fldseas = getNCvar(fname,ncfield,seas='ANN')
    print 'got ANN ' + ncfield + ' data'
    print fldseas.shape

    #fldseaszm = getNCvar(fname,ncfield,seas='ANN',calc='zm')
    #print 'got ANN zonal mean ' + ncfield + ' data'
    #print fldseaszm.shape                     


    #### test 3D var with time
    field = 't'
    ncfield = 'TEMP'
    timstr = '001-061'
    
    fname = basepath + casename + subdir + casename + '_' + field + '_' + timstr + '_ts.nc'

    ## fldall = getNCvar(fname,ncfield)
    ## print 'got all ' + ncfield + ' data'
    ## print fldall.shape

    fldchunk = getNCvar(fname,ncfield,timechunk=(12,12+20*12),calc='zm')
    print 'got zonal mean time chunk 2-22 ' + ncfield + ' data'
    print fldchunk.shape
    print fldchunk.shape[0]/12.
                         
    fldzm = getNCvar(fname,ncfield,calc='zm')
    print 'got zonal mean ' + ncfield + ' data'
    print fldzm.shape

    fldseaszm = getNCvar(fname,ncfield,seas='ANN',calc='zm')
    print 'got ANN zonal mean ' + ncfield + ' data'
    print fldseaszm.shape                         

    # this should fail
    fldseas = getNCvar(fname,ncfield,seas='ANN')
    print 'got ANN ' + ncfield + ' data'
    print fldseas.shape
