""" createallmonthBCsclimo.py
    4/16/2014
    taken from createallmonthsBCsclimo.m, which was on macbook pro

    
"""
import copy
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import platform as platform
import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import cccmaNC as cnc
import cccmacmaps as ccm

cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)
ccm = reload(ccm)
cnc = reload(cnc)

plt.close("all")
plt.ion()

# # # ########## Read NC data ###############
plat = platform.system()

months = con.get_mon()

writefiles = True
skipcontrolwrite=False # skip writing control BC files
testplots = True
doobs=True # 7/8/2014. if True, do obs data
dset='NSIDC'; dtype='bootstrap' # or HadISST,''

field = 'sit' # 'ts', 'sic', 'sit' # these are the model output CMIP names

adjustsst=True # for field='ts'
threshtype = 'abs' # when field='ts' and adjustsst=1. absolute or relative threshold for sea ice conc change
thresh = 10  #when field='ts' and adjustsst=1. threshold for where to change SST

controlsit=0 # for field='sit'. set to 1 if want to keep the control sit.
usesictest=1 # for field='sit' and controlsit=1. use pert sea ice conc to check where sit=0
ensmembers = 0; eindex=4 # make BC files for specified CanESM ensemble member

casenamec = 'historical'; timeperc = '1979-1989'
#casename = 'rcp85';
#timeper = '2022-2032'
#timeper = '2042-2052'
casename = 'historicalrcp45'; timeper='2002-2012' # test against matlab-generated BCs


deni = 913 # kg/m3 density of ice


if field == 'ts': # Surface Temperature (K)
    units = 'K';
    clims=[236.2, 306.2]
    if casename == 'rcp85':
        climsdiff=[-4, 4]
    else:
        climsdiff=[-2, 2]
    cmap='kem_w20'; cmapclimo='blue2red_20'; # actually only 18 colors...
    bcfield='GT'
    bcdescrip='ground temperature (K)'
    bcunits='K'
    
elif field == 'sic': # Sea Ice Fraction (%)
    timeres = 'OImon'
    lmaskflag='lmaskon'
    units = '%'
    clims=[0, 1]
    climsdiff=[-.15, .15]
    cmap='red2blue_w20'
    cmapclimo='Spectral' #'sicclimo_12' # @@ don't think i have in python yet
    bcfield='SICN'
    bcdescrip='sea ice fraction'
    bcunits='frac'# '%' # I think this is an error, should be frac. changed 7/8/2014 but not in bc files up through kem1pert2r4ct
        
elif field == 'sit':
    timeres='OImon'
    lmaskflag='lmaskon'
    units = 'm'
    clims=[-2.5*deni, 2.5*deni] # for now set clims like this until have a better cmap
    climsdiff=[-.5*deni, .5*deni] # @@@
    cmap='red2blue_w20'
    cmapclimo= 'Spectral'# @@ don't have yet? 'sitclimo_12'
    bcfield='SIC'
    bcdescrip='sea ice thickness*density' 
    bcunits='kg/m2' 
else:
    'No such field: ' + field

# set up filenames
if doobs:
    # do observational data
    if plat == 'Darwin':  # means I'm on my mac
        basepath = '/Volumes/MyPassport1TB/DATA/CanSISE/'
        basepath2 = '/Users/kelly/CCCma/CanSISE/BoundaryConditionFiles/'
    else:  # on linux workstation in Vic
        basepath = '/home/rkm/work/BCs/'
        basepath2 = '/home/rkm/work/BCs/'

    timeper = '2002-2011'
    # ########### from matlab
    if dset=='HadISST': # grid 64 x 129
        if field=='ts':
            prfield='gt' 
            #ncfield='GT'
            bcfield='GT'
            bcdescrip='ground temperature (K)'
            bcunits='K'
            units='K' 
            #climsdiff=[-1, 1]
            climsdiff=[-2, 2]
            clims=[236.2, 306.2]
            lmaskflag='lmaskon'
            #conv=1
            cmap='kem_w20'; cmapclimo='blue2red_20'; # actually only 18 colors...
            
        elif field=='sic':
            prfield='sicn'
            #ncfield='SICN'
            bcfield='SICN'
            bcdescrip='sea ice fraction (frac)'
            bcunits='frac'
            units='frac'
            climsdiff=[-.15, .15]
            clims=[0, 1]
            lmaskflag='lmaskon'
            #conv=1
            cmap='red2blue_w20'
            cmapclimo='Spectral' #'sicclimo_12' # @@ don't think i have in python yet
        elif field=='sit':
            prfield='sic'
            lmaskflag='lmaskon'
            units = 'm'
            clims=[-2.5*deni, 2.5*deni] 
            climsdiff=[-.5*deni, .5*deni]
            cmap='red2blue_w20'
            cmapclimo= 'Spectral'# @@ don't have yet? 'sitclimo_12'
            bcfield='SIC'
            bcdescrip='sea ice thickness*density' 
            bcunits='kg/m2'

        filename = basepath + dset + '/hadisst1.1_bc_128_64_1870_2013m03_' + prfield + '_' + timeper + 'climo.nc'
        filenamec = basepath + dset + '/hadisst1.1_bc_128_64_1870_2013m03_' + prfield + '_' + timeperc + 'climo.nc'
        sicfnamep = basepath + dset + '/hadisst1.1_bc_128_64_1870_2013m03_sicn_' + timeper + 'climo.nc'
        sicfnamec = basepath + dset + '/hadisst1.1_bc_128_64_1870_2013m03_sicn_' + timeperc + 'climo.nc'

    elif dset == 'NSIDC': # grid 64 x 129

        if field=='ts':
            # @@ READ HadISST sst
            prfield='gt' 
            #ncfield='GT'
            bcfield='GT'
            bcdescrip='ground temperature (K)'
            bcunits='K'
            units='K' 
            #climsdiff=[-1, 1]
            climsdiff=[-2, 2]
            clims=[236.2, 306.2]
            lmaskflag='lmaskon'
            #conv=1
            cmap='kem_w20'; cmapclimo='blue2red_20'; # actually only 18 colors...
            
        elif field=='sic':
            prfield='sicn'
            #ncfield='SICN'
            bcfield='SICN'
            bcdescrip='sea ice fraction (frac)'
            bcunits='frac'
            units='frac'
            climsdiff=[-.15, .15]
            clims=[0, 1]
            lmaskflag='lmaskon'
            #conv=1
            cmap='red2blue_w20'
            cmapclimo='Spectral' #'sicclimo_12' # @@ don't think i have in python yet
            fbase='nsidc_bt_128x64_1978m11_2011m12'
            
        elif field=='sit':
            prfield='sic'
            lmaskflag='lmaskon'
            units = 'm'
            clims=[-2.5*deni, 2.5*deni] 
            climsdiff=[-.5*deni, .5*deni]
            cmap='red2blue_w20'
            cmapclimo= 'Spectral'# @@ don't have yet? 'sitclimo_12'
            bcfield='SIC'
            bcdescrip='sea ice thickness*density' 
            bcunits='kg/m2'

        filename = basepath + dset + '/nsidc_bt_128x64_1978m11_2011m12_' + prfield + '_' + timeper + 'climo.nc'
        filenamec = basepath + dset + '/nsidc_bt_128x64_1978m11_2011m12_' + prfield + '_' + timeperc + 'climo.nc'
        sicfnamep = basepath + dset + '/nsidc_bt_128x64_1978m11_2011m12_sicn_' + timeper + 'climo.nc'
        sicfnamec = basepath + dset + '/nsidc_bt_128x64_1978m11_2011m12_sicn_' + timeperc + 'climo.nc'

    # ############ end from matlab
# end if doobs
else:
    model = 'CanESM2'
    timeres = 'Amon'

    # Read in CanESM2 data, "pert" first
    
    if plat == 'Darwin':  # means I'm on my mac
        #basepath = '/Users/kelly/CCCma/CanSISE/DATA/' # @@@
        basepath = '/Volumes/MyPassport1TB/DATA/CanSISE/'
        basepath2 = '/Users/kelly/CCCma/CanSISE/BoundaryConditionFiles/'

    else:  # on linux workstation in Vic
        basepath = '/home/rkm/work/DATA/'
        basepath2 = '/home/rkm/work/BCs/'

    enstype = 'ens'
    if ensmembers:
        enstype = 'r' + str(eindex) + 'i1p1'

    filename = basepath + model + '/' + casename + '/' + field + '/' + \
               field + '_' + timeres + '_' + model + '_' + casename + \
               '_' + enstype + '_' + timeper + 'climo.nc'
    # control time period is in the historical run
    filenamec = basepath + model + '/' + casenamec + '/' + field + '/' + \
                field + '_' + timeres + '_' + model + '_' + casenamec + \
                '_' + enstype + '_' + timeperc + 'climo.nc'


landmask = con.get_t63landmask()
if doobs:
    pass # do not remove wraparound lon
else:
    landmask=landmask[...,:-1]
landmask = np.tile(landmask,(12,1,1))
    
print filename

timefldp = cnc.getNCvar(filename,'time')
lat = cnc.getNCvar(filename,'lat')
lon = cnc.getNCvar(filename,'lon')

if doobs:
    fldp = cnc.getNCvar(filename,bcfield)
    fldc = cnc.getNCvar(filenamec,bcfield)
else:
    fldp = cnc.getNCvar(filename,field)
    fldc = cnc.getNCvar(filenamec,field) # we should only need this to set SSTs in near-future BCs

print filenamec


if field == 'ts' and adjustsst:
    # assume adjustsst will always be true now
    print 'TS and adjustsst'

    # also need to get sea-ice concentration to know where to set pert SSTs
    if doobs:
        # make fraction into percent to match the model output
        sicp = cnc.getNCvar(sicfnamep,'SICN')*100 # field is SICN for HadISST and NSIDC
        sicc = cnc.getNCvar(sicfnamec,'SICN')*100
    else:        
        filename = basepath + model + '/' + casename + '/sic/' + \
                   'sic_OImon_' + model + '_' + casename + \
                   '_' + enstype + '_' + timeper + 'climo.nc'
        sicp = cnc.getNCvar(filename,'sic')
        print filename

        filename = basepath + model + '/' + casenamec + '/sic/' + \
                   'sic_OImon_' + model + '_' + casenamec + \
                   '_' + enstype + '_' + timeperc + 'climo.nc'
        print filename

        sicc = cnc.getNCvar(filename,'sic')

    
    sicd = sicp - sicc
    sicpctd = (sicp-sicc)/sicc*100

    fldsave = copy.copy(fldp) # save original pert for comparison & calcs
    fldp = copy.copy(fldc) # start off with control
    # Now where sea ice change > thresh, set the pert SST
    if threshtype == 'abs':
        fldp[sicd<=-thresh] = fldsave[sicd<=-thresh] # @@ will this work in python?
    else:
        fldp[sicpctd<=-thresh] = fldsave[sicpctd<=-thresh]
        
    # Now, also need to check that any open water is not < 271.2 @@
    #print "Not done: check no open water is < 271.2"
    # ideally don't want to change these vals over land though. landmask==-1
    fldptmp=copy.copy(fldp)
    fldptmp=ma.masked_where(np.logical_and(landmask != -1,fldptmp<271.2),fldptmp) # first mask out cold temps where not land    
    fldp[np.logical_and(sicp<=15,fldptmp.mask)] = 271.2 # where SIC is less than 15% (ie "no ice") and not land and it's too cold, set to freezing
    
    #badsst = ma.masked_where(np.logical_and(sicp<=.15,fldp<=271.2),fldp)
    #cplt.map_allmonths(badsst,lat,lon,cmin=238,cmax=308,cmap='blue2red_20',type='nh',lmask=1)
    
elif field == 'ts' and ~adjustsst:
    # perturbation is just control with a check for open water < 271.2K
    print 'TS and ~adjustsst'
    
    fldp = copy.copy(fldc)

    if doobs:
        # make fraction into percent to match the model output
        sicp = cnc.getNCvar(sicfnamep,'SICN')*100 # field is SICN for HadISST and NSIDC
    else:
        # also need to get pert sea-ice concentration to check open water temp against SIConc
        filename = basepath + model + '/' + casename + '/sic/' + \
                   'sic_OImon_' + model + '_' + casename + \
                   '_' + enstype + '_' + timeper + 'climo.nc'
        print filename
        sicp = cnc.getNCvar(filename,'sic')
    
    #print "Not done: check no open water is < 271.2"
    fldptmp=copy.copy(fldp)
    # mask anywhere that is not land and below icefreeze
    fldptmp=ma.masked_where(np.logical_and(landmask!=-1,fldptmp<271.2),fldptmp)
    # anywhere that is not land, below icefreeze AND has sea-ice conc<.15, set to 271.2
#    flddum = fldp[np.logical_and(sicp<=15,fldptmp.mask)]
    fldp[np.logical_and(sicp<=15,fldptmp.mask)] = 271.2
    
    #badsstnoadj = ma.masked_where(np.logical_and(sicp<=.15,fldp<=271.2),fldp)
    #cplt.map_allmonths(badsst,lat,lon,cmin=238,cmax=308,cmap='blue2red_20',type='nh',lmask=1)
    fldsave=copy.copy(fldc)
    
elif field == 'sic':

    if doobs:
        pass # already fractional
    else:
        fldc = fldc/100 # should be fractional
        fldp = fldp/100
    fldsave=copy.copy(fldp) # dummy

elif field == 'sit' and controlsit:
    # keep control sit

    if doobs:
        pass # already a mass
    else:
        fldc = fldc*deni
        fldp = fldp*deni
    
    if usesictest==1: #use pert sea ice concentration as the check
        # this is the one we use for kem1pert3
        
        # need to get pert sea-ice concentration to check thickness against SIConc
        if doobs:
            # make fraction into percent to match the model output
            sicp = cnc.getNCvar(sicfnamep,'SICN')*100 # field is SICN for HadISST and NSIDC
        else:        
            filename = basepath + model + '/' + casename + '/sic/' + \
                       'sic_OImon_' + model + '_' + casename + \
                       '_' + enstype + '_' + timeper + 'climo.nc'
            print filename
            sicp = cnc.getNCvar(filename,'sic')

        fldtmp = copy.copy(fldc) # start with control thickness
        # @@@ should this actually be, anywhere pert sic is <0.15? since the model apparently considers it open water if conc<0.15
        fldtmp[sicp<0.01] = 0 # anywhere pert sic (concentration) is 0, make thickness 0
        fldsave=copy.copy(fldp) # dummy, save orig pert
        fldp = copy.copy(fldtmp)
        print fldtmp.shape # @@
    else:
        #use pert sea ice thickness as the check
        print 'Keep control SIT with usesittest is not implemented. Do not need it anymore?' #@@
    
elif field == 'sit':
    if doobs:
        pass # already a mass
    else:
        fldc = fldc*deni
        fldp = fldp*deni
    fldsave = copy.copy(fldp) # dummy

# put SH back to control since only want to mod NH!
fldp[:,lat<=0,:] = fldc[:,lat<=0,:] # not even sure I need to do this, already started w/ control
# add wraparound lon @@
if doobs:
    pass # already has wraparound lon
else:
    fldc = np.dstack((fldc,fldc[...,0]))
    fldp = np.dstack((fldp,fldp[...,0]))
    fldsave = np.dstack((fldsave,fldsave[...,0]))


# # end setting up the BC data itself ####
#  Now on to dates
# for now just get original dates
if plat == 'Darwin':
    file = basepath + 'CanSISE/CanAM4/amip/uninflatedBCs/temp_gt_1870010100-2011020100.nc'
else:
    file = basepath2 + 'temp_gt_1870010100-2011020100.nc'


timeoutall = cnc.getNCvar(file,'time')
datestr = '1870010100-2010120100'

repfld = len(timeoutall)/12 # want integer division
rem = np.mod(len(timeoutall),12)
# shorten timeoutall by remainder
timeoutall = timeoutall[:-rem]

fldcall = np.tile(fldc,(repfld,1,1))
fldpall = np.tile(fldp,(repfld,1,1))

# in matlab code, I add on the remainder as well (ie 2 months or something).
#   @@ necessary here? it's harder to add this in python. why do I need it?

bclat = cnc.getNCvar(file,'lat')
bclon = cnc.getNCvar(file,'lon')

#timepers = '1979-1989', '2022-2032', # '2042-2052'
outflds = {}
outflds[timeperc] = fldcall
outflds[timeper] = fldpall
casenames = {}
casenames[timeperc] = casenamec
casenames[timeper] = casename

############### TEST PLOTS ###############################
if testplots:
    ## This was to test the selmon operation in CDO
    ## fldtest = np.zeros((12,len(lat),len(lon)))
    ## midx=1
    ## for mon in months:
    ##     fldtest[midx-1,...] = cnc.getNCvar(filenamec,field,monsel=midx)# climo file
    ##     midx=midx+1
        
    ## fig1 = cplt.map_allmonths(fldtest,lat,lon,title='CDO selmon',cmin=clims[0],cmax=clims[1],type='nh',climo=1)
                       
    if field == 'sic' or field == 'sit':
        plotfldp = ma.masked_where(fldp<=0,fldp)
        plotfldc = ma.masked_where(fldc<=0,fldc)
    else:
        plotfldc = fldc
        plotfldp = fldp
    pparams = dict(type='nh',cmin=clims[0],cmax=clims[1],cmap=cmapclimo)

    if field == 'ts':
        pparams['conts'] = 271.2
        
    figp = cplt.map_allmonths(plotfldp,bclat,bclon, title='pert ' + field,**pparams)
    figc = cplt.map_allmonths(plotfldc,bclat,bclon, title='ctrl ' + field,**pparams)
    figd = cplt.map_allmonths(fldp-fldc,bclat,bclon, title='diff',type='nh',cmin=climsdiff[0],cmax=climsdiff[1],cmap=cmap)
    if field == 'ts':
        figd2 = cplt.map_allmonths(fldp-fldsave,bclat,bclon, title='pert - originalpert',type='nh',cmin=climsdiff[0],cmax=climsdiff[1],cmap=cmap)
        figd3 = cplt.map_allmonths(fldsave-fldc,bclat,bclon, title='originalpert - control',type='nh',cmin=climsdiff[0],cmax=climsdiff[1],cmap=cmap)

    figcsq = cplt.map_allmonths(plotfldc,bclat,bclon, title='ctrl ' + field,type='sq',cmin=clims[0],cmax=clims[1],cmap=cmapclimo)
    figdsq = cplt.map_allmonths(plotfldp-plotfldc,bclat,bclon, title='diff',type='sq',cmin=climsdiff[0],cmax=climsdiff[1],cmap=cmap)
    
############### Write the NetCDF file ####################
from netCDF4 import Dataset

if writefiles:

    for thetimeper in outflds.keys():

        if thetimeper==timeperc and skipcontrolwrite:
            pass
        else:

            print thetimeper

            fldout = outflds[thetimeper]
            casenameout = casenames[thetimeper]
            if ensmembers:
                casenameout = casenameout + 'r' + str(eindex)

            print fldout.shape

            # note that timeper is the perturbed time period
            # casenameout is the current casenames element
            
            # If pert time period and adjusted SST:
            if thetimeper==timeper and field == 'ts' and adjustsst:
                if doobs:
                    outfile=bcfield + 'adjusted_BC_' + dset + dtype + '_' + thetimeper + '_' +\
                             datestr + '_' + threshtype + str(thresh) + 'thresh.nc'
                    nctitle = 'Boundary condition dataset generated from ' + dset + dtype +\
                              ' ' + thetimeper + ' adjusted with ' +  outflds.keys()[1] +\
                              ' where SIC ' + threshtype + ' change is >=10%'
                    if dset=='NSIDC':
                        nctitle = nctitle + '. Original SSTs from HadISST1.1'
                else:
                    outfile= bcfield + 'adjusted_BC_CanESM2_' + casenameout + '_' + thetimeper + '_' +\
                         datestr + '_' + threshtype + str(thresh) + 'thresh.nc'
                    
                    nctitle = 'Boundary condition dataset generated from CanESM2 ' +\
                                  casenamec + ' ' + outflds.keys()[0] + ', adjusted with ' +\
                              casename + ' ' + outflds.keys()[1] + ' where SIC ' + threshtype + ' change is >=10%'
                
            # If pert time period and unadjusted SST:
            elif thetimeper==timeper and field == 'ts' and ~adjustsst:
                if ensmembers:
                    outfile = bcfield + 'frzchk' + casename + 'r' + str(eindex) + timeper + '_BC_CanESM2_' +\
                              casenamec + 'r' + str(eindex) + '_' + timeperc + '_' + datestr + '.nc'
                else:
                    if doobs:
                        print 'not planning on needing to make these! @@ No filename'
                        # really should throw exception or something @@
                    else:
                        outfile = bcfield + 'frzchk' + casename + timeper + '_BC_CanESM2_' +\
                              casenamec + '_' + timeperc + '_' + datestr + '.nc'
                # the data is the control, but open water temps are checked for below freezing where no pert ice

                nctitle = 'Boundary condition dataset generated from CanESM2 ' +\
                              casenamec + ' ' + outflds.keys()[0] + ', but open water <271.2K in ' +\
                              casename + ' ' + outflds.keys()[1] + ' is set to 271.2 (where SIC<15)'
            # If pert time period and control SIT:
            elif thetimeper==timeper and field == 'sit' and controlsit:
                
                if usesictest:
                    teststr='usesictest'                   
                else:
                    teststr='usesittest'

                if ensmembers:
                    outfile=bcfield + '_BC_CanESM2_' + casenamec + 'r' + str(eindex) + '_' + timeperc + 'forpert_' +\
                             datestr + '_' + teststr + '.nc'
                else:
                    if doobs:
                        print 'not planning on needing to make these! @@ No filename!'
                        # really should throw exception or something @@
                    else:
                        outfile=bcfield + '_BC_CanESM2_' + casenamec + '_' + timeperc + 'forpert_' +\
                             datestr + '_' + teststr + '.nc'

                nctitle = 'Boundary condition dataset generated from CanESM2 ' +\
                              casenamec + ' ' + outflds.keys()[0] +\
                              ', except set to zero where no ice in ' + casename + ' ' + outflds.keys()[1]
                
            # All other cases have a "regular" outfile name
            else:
                if doobs:
                    outfile = bcfield + '_BC_' + dset + dtype + '_' + thetimeper + '_' + datestr + '.nc'
                    nctitle = 'Boundary condition dataset generated from ' + dset + dtype +\
                              ' ' + thetimeper + ' climo'
                    if dset == 'NSIDC' and bcfield=='SIC':
                        nctitle = nctitle + '. Original data from HadISST thickness (originally from old Can model output)'
                else:
                    outfile=bcfield + '_BC_CanESM2_' + casenameout + '_' + thetimeper + '_' + datestr + '.nc'
                    nctitle = 'Boundary condition dataset generated from CanESM2 ' +\
                              casenameout + ' ' + thetimeper + ' climo'


            outnc = Dataset(outfile,'w')

            # create the dimensions
            #    not sure why these need python variables, they are never used@@
            #    I guess they are the keys in the dict of Dimensions...
            outtime = outnc.createDimension('time', None) # len(timeoutall)) 
            outlat = outnc.createDimension('lat',len(bclat))
            outlon = outnc.createDimension('lon',len(bclon))

            # create variables
            outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
            outlats = outnc.createVariable('lat','d',('lat',))
            outlons = outnc.createVariable('lon','d',('lon',))
            outfld = outnc.createVariable(bcfield,'f4',('time','lat','lon',),fill_value=1.0e38)

            # add attributes to variables
            #outfld._FillValue = np.float(1.0e38) # do in variable creation step
            outfld.units = bcunits
            outfld.long_name = bcdescrip

            outtimes.long_name = 'time'
            outtimes.units = 'days since 1850-1-1'
            outtimes.calendar = '365_day'

            outlats.units = 'degrees_north'
            outlats.long_name = 'Latitude'
            outlats.standard_name = 'latitude'

            outlons.units = 'degrees_east'
            outlons.long_name = 'Longitude'
            outlons.standard_name = 'longitude'

            # global attributes
            import time

            outnc.title = nctitle
                
            if ensmembers:
                outnc.title = outnc.title + '. Ensemble member r' + str(eindex) + 'i1p1'
                
            outnc.creation_date = time.ctime(time.time())
            outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

            # set the data to the variables: important to have [:]!
            outtimes[:] = timeoutall
            outlats[:] = bclat
            outlons[:] = bclon
            outfld[:] = fldout

            outnc.close()

