""" cccmaNC.py
    3/3/2014
    functions to open and read data from Netcdf, with
       option to perform basic operations (some use CDO bindings)
       (seasonal mean, zonal mean)
"""
import numpy as np # for array handling
import scipy as sp # scientific python
from netCDF4 import Dataset
import cccmautils as cutl
import cdo as cdo #; cdo = cdo.Cdo()
import os
import platform as platform


def openNC(filename):

    ncfile = Dataset(filename,'r')
    return ncfile


def getNCvar(filename,field,timesel=None,levsel=None,monsel=None,seas=None,calc=None,remlon=1):
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

        if timesel == '0002-01-01,0111-12-31' or timesel == '0002-01-01,0061-12-31' or timesel=='0002-01-01,0121-12-31':
            print 'hard-coded skipping of first year @@'
            fld = getNCvar_old(filename,field,seas=seas, monsel=monsel,timechunk=(12,),level=level,calc=calc)
        else:
            fld = getNCvar_old(filename,field,seas=seas,monsel=monsel,level=level,calc=calc) # doesn't work with all arguments yet @@
        return fld

    else:  # on linux workstation in Vic

        ncfile = openNC(filename)
        ndims = len(ncfile.dimensions)

        #### READ VARIABLE FROM NC FILE ########
        if timesel == None and calc == None:

            if levsel !=None:
                if monsel != None:
                    fld = np.squeeze(cdo.sellevel(levsel,input = cdo.selmon(monsel,input = filename),returnArray = field))
                else:
                    fld = np.squeeze(cdo.sellevel(levsel,input = filename, returnArray = field))
                os.system('rm -rf /tmp/cdoPy*')
            else:
                if monsel != None:
                    #print 'timesel==None and calc==None and monsel !=None'
                    fld = np.squeeze(cdo.selmon(monsel,input = filename, returnArray = field))
                    #print fld.shape
                    os.system('rm -rf /tmp/cdoPy*')
                else: # get everything
                    fld = ncfile.variables[field][...]

        elif timesel != None and calc == 'zm':
            # have to remove the lon before zonal mean, which means have to separate the
            # select dates and zm. thus can't use CDO for zm (unless can pass it data instead of a file?)

            #fld = np.squeeze(cdo.zonmean( input = cdo.seldate(timesel,input = filename), returnArray = field))
            print 'assuming T42(63) 64x128 resolution for zonal mean'
            if levsel != None:
                if monsel != None:
                    fld = np.squeeze(cdo.seldate(timesel,input = cdo.zonmean( input = 
                                                 cdo.selindexbox(1,128,1,64,input =
                                                                 cdo.sellevel(levsel,input =
                                                                              cdo.selmon(monsel, input = filename)))),
                                                 returnArray = field))# @@@@
                else:
                    fld = np.squeeze(cdo.seldate(timesel,input = cdo.zonmean( input = 
                                                 cdo.selindexbox(1,128,1,64,input =
                                                                 cdo.sellevel(levsel,input = filename))),
                                                 returnArray = field))
            else:
                if monsel != None:
                    fld = np.squeeze(cdo.seldate(timesel,input =
                                                 cdo.zonmean( input =
                                                              cdo.selindexbox(1,128,1,64,input =
                                                                              cdo.selmon(monsel,input = filename))),
                                                              returnArray = field))
                else:
                    fld = np.squeeze(cdo.seldate(timesel,input =
                                                 cdo.zonmean(input =
                                                             cdo.selindexbox(1,128,1,64,input = filename)),
                                      returnArray = field))
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
                                             returnArray = field))
            elif levsel != None and monsel != None:
                fld = np.squeeze(
                    cdo.seldate(timesel,input =
                                cdo.sellevel(levsel,input =
                                             cdo.selmon(monsel,input = filename)),
                                returnArray = field))
            elif levsel == None and monsel != None:
                fld = np.squeeze(cdo.seldate(timesel,input =
                                  cdo.selmon(monsel,input = filename),
                                  returnArray = field))
            else: # levsel and monsel are both None
                fld = np.squeeze(cdo.seldate(timesel,input = filename, returnArray = field))
            os.system('rm -rf /tmp/cdoPy*')
            print "only calc='zm' is implemented now. Returning only selected date range/level/month."

        elif timesel != None:
            if levsel != None and monsel == None:
                fld = np.squeeze(cdo.seldate(timesel,input = cdo.sellevel(levsel,input = filename),returnArray = field))
            elif levsel != None and monsel != None:
                fld = np.squeeze(
                    cdo.seldate(timesel,input =
                                cdo.sellevel(levsel,input =
                                             cdo.selmon(monsel,input = filename)),
                                returnArray = field))
            elif levsel == None and monsel != None:
                fld = np.squeeze(cdo.seldate(timesel,input =
                                  cdo.selmon(monsel,input = filename),
                                  returnArray = field))
            else: # levsel and monsel are both None
                fld = np.squeeze(cdo.seldate(timesel,input = filename, returnArray = field))
                
            os.system('rm -rf /tmp/cdoPy*')
            
        elif calc == 'zm': # and timesel must be None
            print 'assuming T42(63) 64x128 resolution for zonal mean'
            
            if levsel != None and monsel == None:
                fld = np.squeeze(cdo.sellevel(levsel,input =
                                              cdo.zonmean(input =
                                                          cdo.selindexbox(1,128,1,64,input =
                                                                          filename)),
                                              returnArray = field))

            elif levsel != None and monsel != None:
                fld = np.squeeze(
                    cdo.sellevel(levsel,input =
                                 cdo.zonmean(input =
                                             cdo.selindexbox(1,128,1,64,input =
                                                             cdo.selmon(monsel,input =
                                                                        filename))),
                                 returnArray = field))

            elif levsel == None and monsel != None:
                fld = np.squeeze(cdo.zonmean(input =
                                             cdo.selindexbox(1,128,1,64,input =
                                                             cdo.selmon(monsel,input =
                                                                        filename)),
                                                 returnArray = field))
               
            else: # get all data
                fld = np.squeeze(cdo.zonmean(input =
                                             cdo.selindexbox(1,128,1,64,input =
                                                             filename),
                                             returnArray = field))
                
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


        ncfile.close()

        ####### TIME AVERAGE the VARIABLE ##########
        # fld has to be 3d by the time it is passed to func
        #  (time,lev,lat) or (time,lat,lon)
        if seas != None:
            if fld.ndim != 3:
                print 'data must be 3 dimensional to seasonalize()'
                return
            elif monsel != None:
                print "Can't do seasonal average when monsel != None"
                return
            elif seas == 'climo':
                fld,stddev = cutl.climatologize3d(fld)
            elif seas not in ('ANN','DJF','JJA','MAM','SON','NDJ'):
                # means seas is an int value for a month
                fld = cutl.seasonalize_monthlyts(fld,mo=seas)
            else:
                #print 'seasonalizing'
                fld = cutl.seasonalize_monthlyts(fld,seas)
                #print fld.shape

        return fld



def getNCvar_old(filename,field,timechunk=None,monsel=None,level=None,seas=None,calc=None,remlon=1):
    """ gets a variable from netcdf file.
        Time is assumed to be the 1st dimension, Lon is assumed to be the last.
        If any calculations are requested to be performed on the data, the user
        needs to make sure that the requested operations can be performed (b/c
        some of the other functions only handle certain # of dims, etc. My bad.)

        filename: full path to file
        field: NC variable to read in
        timechunk: tuple of start,stop
        monsel: index of month to choose (1=all Jans, 2=all Febs, etc. really meant for CDO bindings)
        level: index of plev to select
        seas: seasonally (annually) average {ANN|DJF|JJA|NDJ|MAM|SON}
        calc: zm (zonal mean),
        remlon: removes extra wrap-around longitude for zonal mean.
                default is 1, remove it
        
        returns fld
        """
#    chunk = None
    
    ncfile = openNC(filename)
    ndims = len(ncfile.dimensions)
    
    #### READ VARIABLE FROM NC FILE ########
    if timechunk == None:
        if level != None:
            fld = ncfile.variables[field][:,level,...]
        else:
            fld = ncfile.variables[field][...]
    else:
        if len(timechunk)==1: # start time until end
            if level != None:
                fld = ncfile.variables[field][timechunk[0]:,level,...]
            else:
                fld = ncfile.variables[field][timechunk[0]:,...]
        else:
            if level != None:
                fld = ncfile.variables[field][timechunk[0]:timechunk[1],level,...]
            else:
                fld = ncfile.variables[field][timechunk[0]:timechunk[1],...]
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
            print 'data must be 3 dimensional to seasonalize()'
            return
        if monsel != None:
            print "Can't do seasonal average when monsel != None"
            return
        else:
            #print 'seasonalizing'
            fld = cutl.seasonalize_monthlyts(fld,seas)
            #print fld.shape

    if monsel != None:
        # send back only selected month
        # monsel is 1-based! ie 1=Jan, 2=Feb. But the index is zero-based
        #   so shift by 1:
        monsel = monsel-1
        fld = fld[monsel::12,...]

        
    return fld


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
