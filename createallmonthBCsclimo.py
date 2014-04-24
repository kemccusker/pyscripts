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

writefiles = 1
skipcontrolwrite=1 # skip writing control BC files
testplots = 0

model = 'CanESM2'
timeres = 'Amon'

field = 'ts' # 'ts', 'sic', 'sit' # these are the model output CMIP names
adjustsst=1
threshtype = 'abs'
thresh = 10
controlsit=0 # set to 1 if want to keep the control sit. not implemented 4/22/14


casenamec = 'historical'; timeperc = '1979-1989'
casename = 'rcp85';
#timeper = '2022-2032'
timeper = '2042-2052'
#casename = 'historicalrcp45'; timeper='2002-2012' # test against matlab-generated BCs


deni = 913 # kg/m3 density of ice


if field == 'ts': # Surface Temperature (K)
    units = 'K';
    clims=[238, 308]
    if casename == 'rcp85':
        climsdiff=[-4, 4]
    else:
        climsdiff=[-2, 2]
    cmap='kem_w20'; cmapclimo=cmap; # actually only 18 colors...
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
    bcunits='%'
        
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


# Read in CanESM2 data, "pert" first

if plat == 'Darwin':  # means I'm on my mac
    #basepath = '/Users/kelly/CCCma/CanSISE/DATA/' # @@@
    basepath = '/Volumes/MyPassport1TB/DATA/CanSISE/'
    basepath2 = '/Users/kelly/CCCma/CanSISE/BoundaryConditionFiles/'
    temppath = '/Users/kelly/CCCma/CanSISE/matlab/'
    
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/'
    basepath2 = '/home/rkm/work/BCs/'
    temppath = basepath

landmask = con.get_t63landmask()
landmask=landmask[...,:-1]
landmask = np.tile(landmask,(12,1,1))


filename = basepath + model + '/' + casename + '/' + field + '/' + \
           field + '_' + timeres + '_' + model + '_' + casename + \
           '_ens_' + timeper + 'climo.nc'
print filename

timefldp = cnc.getNCvar(filename,'time')
lat = cnc.getNCvar(filename,'lat')
lon = cnc.getNCvar(filename,'lon')

fldp = cnc.getNCvar(filename,field)

# control time period is in the historical run
filenamec =  basepath + model + '/' + casenamec + '/' + field + '/' + \
            field + '_' + timeres + '_' + model + '_' + casenamec + \
            '_ens_' + timeperc + 'climo.nc'
print filenamec

fldc = cnc.getNCvar(filenamec,field) # we should only need this to set SSTs in near-future BCs

if field == 'ts' and adjustsst:
    # assume adjustsst will always be true now
    print 'TS and adjustsst'

    # also need to get sea-ice concentration to know where to set pert SSTs
    filename = basepath + model + '/' + casename + '/sic/' + \
               'sic_OImon_' + model + '_' + casename + \
               '_ens_' + timeper + 'climo.nc'
    sicp = cnc.getNCvar(filename,'sic')
    print filename
    
    filename = basepath + model + '/' + casenamec + '/sic/' + \
               'sic_OImon_' + model + '_' + casenamec + \
               '_ens_' + timeperc + 'climo.nc'
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
    fldptmp=ma.masked_where(np.logical_and(landmask != -1,fldptmp<271.2),fldptmp) # first mask out land.    
    fldp[np.logical_and(sicp<=.15,fldptmp.mask)] = 271.2
    
    #badsst = ma.masked_where(np.logical_and(sicp<=.15,fldp<=271.2),fldp)
    #cplt.map_allmonths(badsst,lat,lon,cmin=238,cmax=308,cmap='blue2red_20',type='nh',lmask=1)
    
elif field == 'ts' and ~adjustsst:
    # perturbation is just control with a check for open water < 271.2K
    print 'TS and ~adjustsst'
    
    fldp = copy.copy(fldc)

    # also need to get sea-ice concentration to check open water temp against pert SIConc
    filename = basepath + model + '/' + casename + '/sic/' + \
               'sic_OImon_' + model + '_' + casename + \
               '_ens_' + timeper + 'climo.nc'
    print filename
    sicp = cnc.getNCvar(filename,'sic')
    
    #print "Not done: check no open water is < 271.2"
    fldptmp=copy.copy(fldp)
    fldptmp=ma.masked_where(np.logical_and(landmask!=-1,fldptmp<271.2),fldptmp) # first mask out land.
    fldp[np.logical_and(sicp<=.15,fldptmp.mask)] = 271.2
    
    #badsstnoadj = ma.masked_where(np.logical_and(sicp<=.15,fldp<=271.2),fldp)
    #cplt.map_allmonths(badsst,lat,lon,cmin=238,cmax=308,cmap='blue2red_20',type='nh',lmask=1)
    fldsave=copy.copy(fldc)
    
elif field == 'sic':
    
    fldc = fldc/100 # should be fractional
    fldp = fldp/100
    fldsave=copy.copy(fldp) # dummy

elif field == 'sit' and controlsit:
    # keep control sit
    print 'Keeping control SIT is not implemented. Do not need it anymore?' #@@

elif field == 'sit':
    fldc = fldc*deni
    fldp = fldp*deni
    fldsave = copy.copy(fldp) # dummy

# put SH back to control since only want to mod NH!
fldp[:,lat<=0,:] = fldc[:,lat<=0,:] # not even sure I need to do this, already started w/ control
# add wraparound lon @@
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
        
    figp = cplt.map_allmonths(plotfldp,bclat,bclon, title='pert ' + field,type='nh',cmin=clims[0],cmax=clims[1],cmap=cmapclimo)
    figc = cplt.map_allmonths(plotfldc,bclat,bclon, title='ctrl ' + field,type='nh',cmin=clims[0],cmax=clims[1],cmap=cmapclimo)
    figd = cplt.map_allmonths(plotfldp-plotfldc,bclat,bclon, title='diff',type='nh',cmin=climsdiff[0],cmax=climsdiff[1],cmap=cmap)
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

            print fldout.shape


            if thetimeper==timeper and field == 'ts' and adjustsst: 
                outfile= bcfield + 'adjusted_BC_CanESM2_' + casenameout + '_' + thetimeper + '_' +\
                         datestr + '_' + threshtype + str(thresh) + 'thresh.nc'
            elif thetimeper==timeper and field == 'ts' and ~adjustsst:
                outfile = bcfield + 'freezechkfor_' + casename + timeper + '_BC_CanESM2_' +\
                          casenamec + '_' + timeperc + '_' + datestr + '.nc'           
                # the data is the control, but open water temps are checked for below freezing where no pert ice
            elif thetimeper==timeper and field == 'sit' and controlsit:
                if usesictest:
                    outfile=bcfield + '_BC_CanESM2_' + casenamec + '_' + timeperc + 'forpert_' +\
                             datestr + '_usesictest.nc'
                else:
                    outfile=bcfield + '_BC_CanESM2_' + casenamec + '_' + timeperc + 'forpert_' +\
                             datestr + '_usesittest.nc'

            else:
                outfile=bcfield + '_BC_CanESM2_' + casenameout + '_' + thetimeper + '_' + datestr + '.nc'


            outnc = Dataset(outfile,'w')

            # create the dimensions
            #    not sure why these need python variables, they are never used@@
            #    I guess they are the keys in the dict of Dimensions...
            outtime = outnc.createDimension('time',None) 
            outlat = outnc.createDimension('lat',len(bclat))
            outlon = outnc.createDimension('lon',len(bclon))

            # create variables
            outtimes = outnc.createVariable('time','f8',('time',))
            outlats = outnc.createVariable('lat','f4',('lat',))
            outlons = outnc.createVariable('lon','f4',('lon',))
            outfld = outnc.createVariable(bcfield,'f4',('time','lat','lon',))

            # add attributes to variables
            outfld.FillValue_ = np.float(1.0e38)
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
            if thetimeper == timeper and field == 'ts' and adjustsst:
                outnc.title = 'Boundary condition dataset generated from CanESM2 ' +\
                              casenamec + ' ' + outflds.keys()[0] + ', adjusted with ' +\
                              casename + ' ' + outflds.keys()[1] + ' where SIC ' + threshtype + ' change is >=10%'
            elif thetimeper == timeper and field == 'ts' and ~adjustsst:
                outnc.title = 'Boundary condition dataset generated from CanESM2 ' +\
                              casenamec + ' ' + outflds.keys()[0] + ', but open water <271.2K in ' +\
                              casename + ' ' + outflds.keys()[1] + ' is set to 271.2 (where SIC<.15)'
            elif thetimeper == timeper and field == 'sit' and controlsit:
                outnc.title = 'Boundary condition dataset generated from CanESM2 ' +\
                              casenamec + ' ' + outflds.keys()[0] +\
                              ', except set to zero where no ice in ' + casename + ' ' + outflds.keys()[1]
            else:
                outnc.title = 'Boundary condition dataset generated from CanESM2 ' +\
                              casenameout + ' ' + thetimeper + ' climo'

            outnc.creation_date = time.ctime(time.time())
            outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

            # set the data to the variables: important to have [:]!
            outtimes[:] = timeoutall
            outlats[:] = bclat
            outlons[:] = bclon
            outfld[:] = fldout

            outnc.close()
