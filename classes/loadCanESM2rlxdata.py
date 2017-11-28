
import cccmaNC as cnc
import cccmautils as cutl
import constants as con

#basedir='/HOME/rkm/work/DATA/CanESM2/'
bp=con.get_basepath(runtype='CMIP',model='CanESM2')
basedir = bp['basepath']
fldsfx=''


    
casenames = ('preipreiice', 'prei2xco2iceb','pi2xco2ipulse', '2xco22xco2ice','2xco2preiice')
shortnames = {'preipreiice': '1C_1I',
             'prei2xco2iceb': '1C_2I',
             'pi2xco2ipulse': '1C_2Ipulse',
             '2xco22xco2ice': '2C_2I',
             '2xco2preiice': '2C_1I'}



def get_timeinfo(last='last100',includeyr1=False):

    if last=='last100':
        timeper1="3022-3121"; timeper2='2552-2651'; # end of 200 yr runs
        timeper3='2652-2751'; timeper4='3122-3221' # end of 300 yr runs
        timesel1='3022-01-01,3121-12-31'; timesel2='2552-01-01,2651-12-31'; 
        timesel3='2652-01-01,2751-12-31'; timesel4="3122-01-01,3221-12-31"; 

        timepertmp='2921-3021'; 
        timeseltmp='2922-01-01,3021-12-31';

        suff='last100yrof200' # for figure names
    elif last=='last200':
        timeper1='2922-3121'; timeper2='2452-2651'; 
        timesel1='2922-01-01,3121-12-31'; timesel2='2452-01-01,2651-12-31'; 

        timeper1='2922-3121'; timeper2='2452-2651'; # end of 200 yr runs
        timeper3='2552-2751'; timeper4='3022-3221' # end of 300 yr runs
        timesel1='2922-01-01,3121-12-31'; timesel2='2452-01-01,2651-12-31'; 
        timesel3='2552-01-01,2751-12-31'; timesel4='3022-01-01,3221-12-31'; 

        timepertmp='2921-3021'; 
        timeseltmp='2922-01-01,3021-12-31'; 

        #fldsfx='NH'
        #suff='last200yrof300' # for figure names

    elif last=='pen100': # penultimate 100 years
        if includeyr1:
            timeper1='2921-3021'; timeper2='2451-2551'; # same as first100
            timeper3='2552-2651'; timeper4='3022-3121'
            timesel1='2921-01-01,3021-12-31'; timesel2='2451-01-01,2551-12-31'
            timesel3='2552-01-01,2651-12-31'; timesel4='3022-01-01,3121-12-31'

            #suff='first101yr' # for figure names
        else:
            timeper1='2921-3021'; timeper2='2451-2551'; 
            timeper3='2552-2651'; timeper4='3022-3121'
            timesel1='2922-01-01,3021-12-31'; timesel2='2452-01-01,2551-12-31'; 
            timesel3='2552-01-01,2651-12-31'; timesel4='3022-01-01,3121-12-31'

        timepertmp=timeper1; #timepertmp2=timeper2
        timeseltmp=timesel1; #timeseltmp2=timesel2
        
    elif last=='first100': # first 100 years
        if includeyr1:
            timeper1='2921-3021'; timeper2='2451-2551'; timeper3=timeper2
            timesel1='2921-01-01,3021-12-31'; timesel2='2451-01-01,2551-12-31'; timesel3=timesel2

            #suff='first101yr' # for figure names
        else:
            timeper1='2921-3021'; timeper2='2451-2551'; timeper3=timeper2
            timesel1='2922-01-01,3021-12-31'; timesel2='2452-01-01,2551-12-31'; timesel3=timesel2

        timepertmp=timeper1; #timepertmp2=timeper2
        timeseltmp=timesel1; #timeseltmp2=timesel2

            #suff='first100yr' # for figure names

    timepers = {'preipreiice': timeper1, 'prei2xco2iceb': timeper4, 
                'pi2xco2ipulse': timepertmp, '2xco22xco2ice': timeper2,
                '2xco2preiice': timeper3} # for file names
    timesels = {'preipreiice': timesel1, 'prei2xco2iceb': timesel4, 
                'pi2xco2ipulse': timeseltmp, '2xco22xco2ice': timesel2,
                '2xco2preiice': timesel3} # for .nc selection

    return timepers,timesels


def load_nclatlon(field, last='last100',includeyr1=False,verb=False,local=False):
    """ return lat,lon from given field nc file. Uses preipreiice casename
                  
          if local=True, just gets lat/lon from DJF mean files (full globe...). ignores 'field', 'last'
    """
    
    timepers,timesels=get_timeinfo(last,includeyr1)
    
    casename = 'preipreiice'
    if local:
        #bd='/Users/kelly/DropboxWork/UWashCanSISE/DATA/relaxationruns/' #limited fields here
        bd='/Users/kelly/DATA/DataDisk/' # limited but not shared so there are more data
        fname=bd+'preipreiice/preipreiice_st_2922-3121_DJF_mean.nc' 
    else:
        fname= basedir + casename +'/ts/' + casename + '_' + field + '_' + timepers[casename] + '_ts.nc'
        
    if verb:
        print fname

    tmplat=cnc.getNCvar(fname,'lat')
    if verb:
        print 'lat.shape ' + str(tmplat.shape) 
    tmplon=cnc.getNCvar(fname,'lon')
    if verb:
        print 'lon.shape ' + str(tmplon.shape) 

    return tmplat,tmplon

def load_nclev(field,verb=False):

    timepers,timesels=get_timeinfo()
    
    casename = 'preipreiice'
    fname= basedir + casename +'/ts/' + casename + '_' + field + '_' + timepers[casename] + '_ts.nc'
    if verb:
        print fname

    lev=cnc.getNCvar(fname,'plev')    

    return lev

def load_ncfield(field, ncfield, zonal=True,conv=1,mo=0,season=None,last='last100',
                     includeyr1=False,verb=False,local=False,seacyc=False,pulse=False):
    """ 
        If zonal=True, compute zonal mean.
           zonal can also be a region, defined in regiondict in constants.py
        If local, use local data (which is limited)
        If seacycle, return the 12 month climo average, which only exists 
                1. locally (as of Mar 24 2017) and 2. is only avail for the last200

        returns dictionaries: full field w/ time, zonal mean w/ time
    """


    timepers,timesels = get_timeinfo(last,includeyr1)

    if local:
        #bd='/Users/kelly/DropboxWork/UWashCanSISE/DATA/relaxationruns/' #limited fields here
        bd='/Users/kelly/DATA/DataDisk/' # limited but not shared so there are more data
        subdir='/'
    else:
        bd=basedir
        subdir='/ts/'
        
    seasonalizedt = {'mo':mo, 'season': season}
    
    ncflddt={}; ncfldtmdt={}; ncfldzmdt={}
    for casename in casenames:
        
        if pulse==False and ('pulse' in casename): continue
            
        if seacyc:
            if last!='last200':
                raise Exception('load_ncfield(): if seacyc=True, last must = last200')
            if local!=True:
                raise Exception('load_ncfield(): if seacyc=True, local must be True (as of Mar 24 2017)')
                
            fname = bd + casename+subdir+casename + '_' + field + '_' + timepers[casename] + '_seacycle.nc'
        else:
            fname= bd + casename +subdir + casename + '_' + field + '_' + timepers[casename] + '_ts.nc'
        timesel=timesels[casename]
        if verb:
            print fname

        try:
            fld=cnc.getNCvar(fname,ncfield,timesel=timesel)*conv # remlon not operational??? why commented out??
            if verb:
                print fname + ', fld.shape ' + str(fld.shape) # @@@
        except IOError as rte: # previously this was a RuntimeError. Not sure what changed. 10/14/2016
            print rte
            if 'No such file or directory' in rte.args:
                # Here the problem is most likely that pico2ipulse filename fails b/c 
                #   there is no NH version of it (for when last200=True). Remove the NH and try again.
                #   @@ Of course this could be a problem too b/c now we have full globe for just this run.
                if field[-2:] != 'NH':
                    print 'not an NH versus global problem...File really not found'
                    raise
                    
                fname= basedir + casename +'/ts/' + casename + '_' + field[:-2] + '_' + timepers[casename] + '_ts.nc'
                print 'Try again: ' + fname
                fld=cnc.getNCvar(fname,ncfield,timesel=timesel)*conv 
                tmplat=cnc.getNCvar(fname,'lat')
                if last=='last200': # how to tell which spot lat is?
                    if zonal: # then it is likely time x lat x lon
                        #print 'zonal is true. fld dims? ' + str(fld.shape)
                        fld=fld[...,tmplat>0,:]
                    else: # likely time x lev x lat
                        fld=fld[...,tmplat>0]
            else:
                raise

        # just get the full timeseries, and zonal mean if asked.
        if 1: #timeseries: # do not time average (actually save time average anyway)
            if (season==None and mo==0):
                pass
            else:
                fld = cutl.seasonalize(fld,**seasonalizedt)
            #fldtm = fld
            # zonal average
            if zonal==True:
                fldzm = fld[...,:-1].mean(axis=2) # remove extra lon
            elif zonal in con.get_regiondict().keys():
                print 'lat,lon shapes ' + str(lat.shape),str(lon.shape)
                tempo, _ = cutl.mask_region(fld,lat,lon,region=zonal)
                fldzm = tempo.mean(axis=2)                
            else:
                fldzm = None 
            #fldtm = fldtmzm.mean(axis=0)
#        else:
#            # seasonal & time mean
#            fldtm = np.mean(cutl.seasonalize(fld,**seasonalizedt),axis=0) 
#            # zonal average
#            if zonal:
#                fldtmzm = fldtm[...,:-1].mean(axis=1)
#            elif zonal in con.get_regiondict().keys():                
#                tempo, _ = cutl.mask_region(fldtm,lat,lon,region=zonal)
#                fldtmzm = tempo.mean(axis=1)        
#            else:
#                fldtmzm = fldtm

        print 'fld.shape: ' + str(fld.shape)
        if fldzm!=None:
            print 'fldzm.shape: ' + str(fldzm.shape)

        ncflddt[casename] = fld # keep all dims
        #ncfldtmdt[casename] = fldtm # time mean
        ncfldzmdt[casename] = fldzm # zonal mean  (w/ time dim)

        # try this: add coords # @@@@@@@@@@
        lat,lon = load_nclatlon(field,last=last,includeyr1=includeyr1,verb=verb,local=local)
        ncflddt['lat'] = lat; ncflddt['lon'] = lon;
        try:
            lev = load_nclev(field,last=last,includeyr1=includeyr1,verb=verb)
            ncflddt['lev'] = lev
        except:
            print "No lev coord. Leave it"
        
        ncfldzmdt['lat'] = lat
        
    return ncflddt, ncfldzmdt
