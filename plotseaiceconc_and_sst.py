import numpy as np
import numpy.ma as ma
import scipy as sp # scientific python
import scipy.stats
import matplotlib.pyplot as plt
import platform as platform
import cccmaplots as cplt
import constants as con
import cccmautils as cutl
import cccmaNC as cnc
import matplotlib.font_manager as fm

cplt = reload(cplt)
con = reload(con)
cutl = reload(cutl)
cnc = reload(cnc)

plt.close("all")
plt.ion()

printtofile=1

casename = 'kemctl1'
timstr = '001-111'
timesel = '0002-01-01,0111-12-31'
casenamep2 = 'kem1pert2'  # 2002-2012 sic, sit, adjusted sst

model = 'CanAM4'


plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/'
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'


fields = 'sicn'
unitss='%'
convs=100
cmincs = -100; cmaxcs = 100 # climo
cmins=-15; cmaxs=15
cminms=-30; cmaxms=30 # monthly anom
cmaps = 'red2blue_w20'
cmaps = 'bone'

fieldg = 'gt'
unitsg='C'
convg=1
#cmincg = 238; cmaxcg = 308
cmincg = -35; cmaxcg = 35
cmincnhg = -45; cmaxcnhg = 45 # polar NH
cming=-2; cmaxg=2
cminmg=-3; cmaxmg=3
cmapg = 'blue2red_20'


fnamecs = basepath + casename + subdir + casename + '_' + fields + '_' + timstr + '_ts.nc'
fnamep2s = basepath + casenamep2 + subdir + casenamep2 + '_' + fields + '_' + timstr + '_ts.nc'

fldcs = cnc.getNCvar(fnamecs,fields.upper(),timesel='0002-01-01,0111-12-31')*convs
fldp2s = cnc.getNCvar(fnamep2s,fields.upper(),timesel='0002-01-01,0111-12-31')*convs

fnamecg = basepath + casename + subdir + casename + '_' + fieldg + '_' + timstr + '_ts.nc'
fnamep2g = basepath + casenamep2 + subdir + casenamep2 + '_' + fieldg + '_' + timstr + '_ts.nc'

fldcg = cnc.getNCvar(fnamecg,fieldg.upper(),timesel='0002-01-01,0111-12-31')*convg
fldp2g = cnc.getNCvar(fnamep2g,fieldg.upper(),timesel='0002-01-01,0111-12-31')*convg

lat = cnc.getNCvar(fnamecs,'lat')
lon = cnc.getNCvar(fnamecs,'lon')

fldcsclim,fldcsstd = cutl.climatologize3d(fldcs)
fldpsclim,fldpsstd = cutl.climatologize3d(fldp2s)

fldcgclim,fldcgstd = cutl.climatologize3d(fldcg)
fldpgclim,fldpgstd = cutl.climatologize3d(fldp2g)

print fldcgclim.shape


lons,lats = np.meshgrid(lon,lat)


plotfldg = fldcgclim[2,:,:] # test Jan sst
plotflds = fldcsclim[2,:,:] # test Jan sicn
plotflds = ma.masked_where(plotflds<=0,plotflds)

#fig = plt.figure()
#cplt.kemmap(plotfldg,lat,lon,cmin=cmincg,cmax=cmaxcg,cmap=cmapg)

fig2 = plt.figure()
bm,pc = cplt.kemmap(plotfldg,lat,lon,cmin=cmincnhg,cmax=cmaxcnhg,cmap=cmapg,type='nh')
#cf = bm.contour(lons,lats,plotflds,colors='k',latlon=True)
cplt.kemmap(plotflds,lat,lon,cmin=cmincs,cmax=cmaxcs,cmap=cmaps,type='nh',suppcb=1)
cf2 = bm.contour(lons,lats,plotflds,colors='k',latlon=True,linewidths='2',levels=[15, 15])


#fig3 = plt.figure()
#cplt.kemmap(plotflds,lat,lon,cmin=cmincs,cmax=cmaxcs,cmap=cmaps,type='nh')


months=con.get_mon()
title = casename + ' ' + fields + ' and ' + fieldg

fig3, spax = plt.subplots(2,6)
fig3.set_size_inches(12,6)
fig3.subplots_adjust(hspace=0,wspace=0)

midx=0
for ax in spax.flat:

    plotfldg = fldcgclim[midx,:,:]
    plotflds = fldcsclim[midx,:,:]
    plotflds = ma.masked_where(plotflds<=0,plotflds)

    bm,pc = cplt.kemmap(plotfldg,lat,lon,cmin=cmincnhg,cmax=cmaxcnhg,cmap=cmapg,type='nh',suppcb=1,axis=ax)
    #cf = bm.contour(lons,lats,plotflds,colors='k',latlon=True)
    cplt.kemmap(plotflds,lat,lon,cmin=cmincs,cmax=cmaxcs,cmap=cmaps,type='nh',suppcb=1,axis=ax)
    cf2 = bm.contour(lons,lats,plotflds,colors='k',latlon=True,linewidths='2',levels=[15, 15])
    ax.set_title(months[midx])
    midx = midx+1

cbar_ax = fig3.add_axes([.91,.25, .02,.5])
fig3.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
plt.suptitle(title)

if printtofile: # @@ pdf/eps/ps save doesn't work (errors), why?
    fig3.savefig(fields + '_' + fieldg + '_' + casename + '_allmos_nh.png')
