""" createNSIDCclimos.py
    July 10, 2014

    Need to generate and SST and SIT dataset that is consistent
      with NSIDC sea ice concentration. We are using the climos
      from HadISST (and CCCma for thickness but used w/ HadISST).
"""

import copy as cp
import numpy as np
import numpy.ma as ma
import platform as platform
import cccmaplots as cplt
import cccmaNC as cnc



plt.close("all")
plt.ion()

testplots=True
testsicplots=False # just the sea-ice concentration plots
writefiles=True
dossts=False # else do SIT

# # # ########## Read NC data ###############
plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Volumes/MyPassport1TB/DATA/CanSISE/'
    basepath2 = '/Users/kelly/CCCma/CanSISE/BoundaryConditionFiles/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/'
    basepath2 = '/home/rkm/work/BCs/'
    
# read in HadISST SST & Thickness (and SIC for comparison)
# read in NSIDC SIC

timeper='2002-2011'
timeperc='1979-1989'

# NSIDC first
dset = 'NSIDC'

nsicfnamep = basepath2 + dset + '/nsidc_bt_128x64_1978m11_2011m12_sicn_' + timeper + 'climo.nc'
nsicfnamec = basepath2 + dset + '/nsidc_bt_128x64_1978m11_2011m12_sicn_' + timeperc + 'climo.nc'

nsicp = cnc.getNCvar(nsicfnamep,'SICN')
nsicc = cnc.getNCvar(nsicfnamec,'SICN')


# HadISST
dset='HadISST'
sstfnamep = basepath2 + dset + '/hadisst1.1_bc_128_64_1870_2013m03_gt_' + timeper + 'climo.nc'
sstfnamec = basepath2 + dset + '/hadisst1.1_bc_128_64_1870_2013m03_gt_' + timeperc + 'climo.nc'
sitfnamep = basepath2 + dset + '/hadisst1.1_bc_128_64_1870_2013m03_sic_' + timeper + 'climo.nc'
sitfnamec = basepath2 + dset + '/hadisst1.1_bc_128_64_1870_2013m03_sic_' + timeperc + 'climo.nc'

hsicfnamep = basepath2 + dset + '/hadisst1.1_bc_128_64_1870_2013m03_sicn_' + timeper + 'climo.nc'
hsicfnamec = basepath2 + dset + '/hadisst1.1_bc_128_64_1870_2013m03_sicn_' + timeperc + 'climo.nc'

hsicp = cnc.getNCvar(hsicfnamep,'SICN')
hsicc = cnc.getNCvar(hsicfnamec,'SICN')

cmap='red2blue_w20'
if testsicplots:
    # compare sic
    #   looks like NSIDC has more ice in the control time period and ~less in the pert. so, greater trend
    cplt.map_allmonths(nsicp-hsicp,lat,lon,cmin=-.1,cmax=.1,cmap=cmap,type='nh',lmask=1,title='pert sic nsidc-hadisst',latlim=45)
    cplt.map_allmonths(nsicc-hsicc,lat,lon,cmin=-.1,cmax=.1,cmap=cmap,type='nh',lmask=1,title='ctl sic nsidc-hadisst',latlim=45)

    cplt.map_allmonths((nsicp-nsicc),lat,lon,cmin=-.1,cmax=.1,cmap=cmap,type='nh',lmask=1,title='pert-ctl sic nsidc',latlim=45)
    cplt.map_allmonths((hsicp-hsicc),lat,lon,cmin=-.1,cmax=.1,cmap=cmap,type='nh',lmask=1,title='pert-ctl sic hadisst',latlim=45)

    cplt.map_allmonths((nsicp-nsicc) - (hsicp-hsicc),lat,lon,cmin=-.1,cmax=.1,cmap=cmap,type='nh',lmask=1,title='pert-ctl sic nsidc-hadisst',latlim=45)


if dossts:
    # ########### do SST  ####################################
    hsstp = cnc.getNCvar(sstfnamep,'GT')
    hsstc = cnc.getNCvar(sstfnamec,'GT')
    lat = cnc.getNCvar(sstfnamec,'lat')
    lon = cnc.getNCvar(sstfnamec,'lon')
    timefld = cnc.getNCvar(sstfnamec,'time')

    #    basically take these hadisst sst's and mod for NSIDC

    # 1. if NSIDC SIC >=0.15 (ice) and HadISST SIC < 0.15 (no ice)
    #      SST should be set to freezing
    # 2. if NSIDC SIC >=0.15 (ice) and HadISST SIC >= 0.15 (ice)
    #      use HadISST SST as is
    # 3. if NSIDC SIC < 0.15 (no ice) and HadISST SIC >= 0.15 (ice)
    #      SST should be set to freezing if not already (easiest solution. will underestimate SST)
    # 4. if NSIDC SIC < 0.15 (no ice) and HadISST SIC < 0.15 (no ice)
    #      use HadISST SST as is
    #
    # So, only change the SST field for case 1. and 3. Set to freezing for both

    nsstp = cp.copy(hsstp)
    nsstc = cp.copy(hsstc)

    # case 1.
    nsstp[np.logical_and(nsicp>=0.15, hsicp<0.15)] = 271.2
    nsstc[np.logical_and(nsicc>=0.15, hsicc<0.15)] = 271.2

    # case 3.
    nsstp[np.logical_and(nsicp<0.15, hsicp>=0.15)] = 271.2
    nsstc[np.logical_and(nsicc<0.15, hsicc>=0.15)] = 271.2

    if testplots:
        # compare original hadisst to new nsidcsst: each has very small regions of adjusted SST.
        cplt.map_allmonths(nsstp-hsstp,lat,lon,cmin=-.1,cmax=.1,cmap=cmap,type='nh',lmask=1,title='pert sst nsidc-hadisst',latlim=45)
        cplt.map_allmonths(nsstc-hsstc,lat,lon,cmin=-.1,cmax=.1,cmap=cmap,type='nh',lmask=1,title='ctl sst nsidc-hadisst',latlim=45)

    # set up SST data for writing out
    outflds = {}
    outflds[timeperc] = nsstc
    outflds[timeper] = nsstp
    bcfield = 'GT'
    bcunits = 'K'
    bcdescrip = 'Ground Temperature'
    printfld = 'gt'

else: # do sea ice thickness

    deni = 913. # kg/m3
    
    hsitp = cnc.getNCvar(sitfnamep,'SIC')
    hsitc = cnc.getNCvar(sitfnamec,'SIC')
    lat = cnc.getNCvar(sitfnamec,'lat')
    lon = cnc.getNCvar(sitfnamec,'lon')
    timefld = cnc.getNCvar(sitfnamec,'time')

    """ Email from Mike Lazare 6/24/2014 entitled Re: sea ice thickness boundary condition
    
           In general, a point is considered water only the sea-ice 
    concentration reaches 0.15 (15%). For concentrations less than this, any 
    sea-ice mass
    is ignored essentially, even if it is non-zero. Typically, a 
    concentration of 0.15 maps to a sea-ice mass of 45 Kg/m2 (5 cms of ice), 
    so it is likely
    that there will be non-zero mass present. The BC forcing data typically 
    tries to do consistency checks for mass vs concentration, so when
    the concentration exceeds 0.15 the associated mass will typically be > 
    45 Kg/m2, and the model will interpolate between the previous value
    (less than 45 Kg/m2) and the target (value >45 Kg/m2).

       All this is a roundabout (and somewhat rambling) confirmation of your 
    hypothesis that it ignores the sea ice thickness. The ice present does
    not contribute to albedo until the concentration exceeds 0.15; 
    thereafter both snow/ice and open water (ie leads) factor in.

    """

    #    basically take these hadisst thickness and mod for NSIDC
    
    #    Based on Mike Lazare's email, probably don't have to care too much about the
    #      tolerance for checking SIC b/c the model considers the grid cell to
    #      be open water only if <0.15 (and it doesn't count for albedo either)
    #    So actually just switch to 0.15, similar to SST tests above.

    # 1. if NSIDC SIC >0 (ice) and HadISST SIC = 0 (no ice)
    #      SIT should be what?? 2m? No, 45kg/m3 (~5cm thick) based on the original
    #        way the hadisst thickness was created from some model climatology
    #        script is: to_gt_sic_sicn_job
    # 2. if NSIDC SIC >0 (ice) and HadISST SIC > 0 (ice)
    #      use HadISST SIT as is
    # 3. if NSIDC SIC = 0 (no ice) and HadISST SIC > 0 (ice)
    #      SIT should be set to zero 
    # 4. if NSIDC SIC = 0 (no ice) and HadISST SIC = 0 (no ice)
    #      use HadISST SIT as is (0)
    #
    # So, only change the SIT field for case 1. and 3.
    # Should I use a tolerance for the zero check? yes, but what. @@

    nsitp = cp.copy(hsitp)
    nsitc = cp.copy(hsitc)

    # case 1.
    nsitp[np.logical_and(nsicp>=0.15, hsicp<0.15)] = 45 # ~5cm thickness
    nsitc[np.logical_and(nsicc>=0.15, hsicc<0.15)] = 45 # ~5cm thickness
    
    # case 3.
    nsitp[np.logical_and(nsicp<0.15, hsicp>=0.15)] = 0
    nsitc[np.logical_and(nsicc<0.15, hsicc>=0.15)] = 0

    if testplots:

        # compare original hadisst SIT to new nsidcsit:
        #   if blue, nsidc has greater thickness than hadisst, which is where was set to 2m.
        cplt.map_allmonths(nsitp-hsitp,lat,lon,cmin=-.1,cmax=.1,cmap=cmap,type='nh',lmask=1,title='pert sit nsidc-hadisst',latlim=45)
        cplt.map_allmonths(nsitc-hsitc,lat,lon,cmin=-.1,cmax=.1,cmap=cmap,type='nh',lmask=1,title='ctl sit nsidc-hadisst',latlim=45)


    # set up SIT data for writing out
    outflds = {}
    outflds[timeperc] = nsitc
    outflds[timeper] = nsitp
    bcfield = 'SIC'
    bcunits = 'kg/m2'
    bcdescrip = 'sea ice thickness*density (kg/m2)'
    printfld = 'sic'



############### Write the NetCDF file ####################
from netCDF4 import Dataset

if writefiles:

    for thetimeper in outflds.keys():

        print thetimeper

        fldout = outflds[thetimeper]
        # for some reason the long filename doesn't work, but a short filename does.
        # some talk about that here, but I don't understand why?
        # https://code.google.com/p/netcdf4-python/issues/detail?id=141
        
        #outfile = 'nsidc_bt_128x64_1978m11_2011m12_' + printfld + '_' + thetimeper + 'climo.nc'
        outfile = 'nsidc_' + printfld + '_' + thetimeper + '.nc'

        outnc = Dataset(outfile,'w')

        # create the dimensions
        #    not sure why these need python variables, they are never used@@
        #    I guess they are the keys in the dict of Dimensions...
        outtime = outnc.createDimension('time', None)
        outlat = outnc.createDimension('lat',len(lat))
        outlon = outnc.createDimension('lon',len(lon))

        # create variables
        outtimes = outnc.createVariable('time','f8',('time',)) # f8 and d are the same dtype
        outlats = outnc.createVariable('lat','d',('lat',))
        outlons = outnc.createVariable('lon','d',('lon',))
        outfld = outnc.createVariable(bcfield,'f4',('time','lat','lon',),fill_value=1.0e38)

        # add attributes to variables
        outfld.units = bcunits
        outfld.long_name = bcdescrip
        outfld.grid_type = 'gaussian'
        outfld.missing_value = 1.0e38

        outtimes.long_name = 'time'
        outtimes.units = 'days since 1850-01-01 00:00:00'
        outtimes.calendar = '365_day'
        outtimes.standard_name = 'time'

        outlats.units = 'degrees_north'
        outlats.long_name = 'Latitude'
        outlats.standard_name = 'latitude'
        outlats.axis = 'Y'

        outlons.units = 'degrees_east'
        outlons.long_name = 'Longitude'
        outlons.standard_name = 'longitude'
        outlons.axis = 'X'

        # global attributes
        import time
        outnc.title = 'original ' + bcfield + ' from hadisst1.1_bc_128_64_1870_2013m03_' + \
                      printfld + '_' + thetimeper + 'climo.nc, adjusted for NSIDC sea ice conc'

        outnc.creation_date = time.ctime(time.time())
        outnc.created_by = 'Kelly E. McCusker, CCCma / U. of Victoria'

        # set the data to the variables: important to have [:]!
        outtimes[:] = timefld
        outlats[:] = lat
        outlons[:] = lon
        outfld[:] = fldout

        outnc.close()

