from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from eofs.standard import Eof
import constants as con
import cccmautils as cutl


# example here: http://ajdawson.github.io/eofs/examples/nao_standard.html
printtofile=True

docorr=True # EOF as correlation (otherwise covariance)
enum=1 # which EOF to retrieve and plot


# As Mori et al: 20-90N, 0-180E DJF SAT anomalies

sim='NSIDC'
field='gz50000'; ncfield='PHI'
#field='st'; ncfield='ST'
#field='pmsl'; ncfield='PMSL'
sea='DJF'

subset=True # rather than calc eof on whole globe, do a sub region
substr='morieof' # @@ make sure to also change lims below when changing this.
#substr='nh20N'

type='nh' # plot type
latlims=[20,90]
lonlims=[0,180]
#lonlims=[0,360]
limsdict = {'latlims':latlims, 'lonlims': lonlims}


if docorr: # show EOF as correlation. otherwise, covariance
    cmin=-1; cmax=1
else:
    cmin=''; cmax=''

fnamec,fnamep=con.build_filepathpair(sim,field)

fldc=cnc.getNCvar(fnamec,ncfield,seas=sea)
fldp=cnc.getNCvar(fnamep,ncfield,seas=sea)

lon=cnc.getNCvar(fnamec,'lon')
lat=cnc.getNCvar(fnamec,'lat')

fldca=np.squeeze(fldc-np.mean(fldc,axis=0)) # remove time mean
fldpa=np.squeeze(fldp-np.mean(fldp,axis=0))




if subset:
    msk1 = con.get_t63regionmask('other',limsdict=limsdict)
    fldca,msk=cutl.mask_region(fldca,lat,lon,'other',limsdict=limsdict)
    fldpa,msk=cutl.mask_region(fldpa,lat,lon,'other',limsdict=limsdict)
    
    oshape=fldca.shape # original shape
    print fldca.shape
    
    fldca1=fldca[~fldca.mask] # @@@ why do I have to invert it. so dumb.
    #fldca2=fldca[msk]
    fldpa1=fldpa[~fldpa.mask]

    nshape=fldca1.shape # new shape
    print fldca1.shape
    #print fldca2.shape
    
    # only keep region of interest?
    lato=lat
    lono=lon
    lat=lat[np.logical_and(lat>=latlims[0],lat<=latlims[1])]
    lon=lon[np.logical_and(lon>=lonlims[0],lon<=lonlims[1])]

    fldca1 = fldca1.reshape( (oshape[0], len(lat),len(lon))  )  #fldcash[0]/oshape[0])) # keep time dim
    fldca=fldca1
    fldpa1 = fldpa1.reshape( (oshape[0], len(lat),len(lon))  ) 
    fldpa=fldpa1

coslat = np.cos(np.deg2rad(lat)).clip(0., 1.) # why square root of cos lat? (better for EOF for some reason.)
wgts = np.sqrt(coslat)[..., np.newaxis]



solverc=Eof(fldca,weights=wgts)
if docorr:
    eof1c=solverc.eofsAsCorrelation(neofs=enum)
    eof1c=eof1c[enum-1,...]
else:
    eof1c=solverc.eofsAsCovariance(neofs=enum)
    eof1c=eof1c[enum-1,...]

eof1c=eof1c.squeeze()
eigsc=solverc.eigenvalues()
vexpc=eigsc[enum-1]/eigsc.sum() *100 # percent variance explained

fig,axs=plt.subplots(1,4)
fig.set_size_inches(12,5)
ax=axs[0]
cplt.kemmap(eof1c,lat,lon,type=type, title=sim + ' control EOF' + str(enum),axis=ax,cmin=cmin,cmax=cmax)
ax.set_ylabel(str(np.round(vexpc)))

solverp=Eof(fldpa,weights=wgts)
if docorr:
    eof1p=solverp.eofsAsCorrelation(neofs=enum)
    eof1p=eof1p[enum-1,...]
else:
    eof1p=solverp.eofsAsCovariance(neofs=enum)
    eof1p=eof1p[enum-1,...]

eof1p=eof1p.squeeze()
eigsp=solverp.eigenvalues()
vexpp=eigsp[enum-1]/eigsp.sum() *100 # percent variance explained

ax=axs[1]
cplt.kemmap(eof1p,lat,lon,type=type, title='pert EOF',axis=ax,cmin=cmin,cmax=cmax)
ax.set_ylabel(str(np.round(vexpp)))

ax=axs[2]
cplt.kemmap(eof1p-eof1c,lat,lon,type=type,title='diff of EOF',axis=ax,cmin=cmin,cmax=cmax)


# Now because I have two 'equilibrium' simulations that I want to difference,
#   do I calc the eofs first for the pert and control, and then difference?
#   Or do I calc the difference, remove the time mean, and calc eof? @@@@
# Also, what domain should I include in the eof calc? this is the whole globe I think.

# try taking the eof of the difference. Note it is completely different.
# This version seems more arbitrary to me, because we are differencing two
# timeseries where each year could theoretically be in any order, they just
# happen to be in this order....

# Then again, the sign of the EOF is supposedly arbitrary, so differencing 2 EOFS
#   may be meaningless. For now, calc difference timeseries as pert-mean(control). Then, remove the
#   time mean of that, and calc EOF

# @@@ Also, I want to actually plot the EOFs in physical units. I forget how -- check notes
#     from objective analysis


#fldd=fldp-fldc
fldd=fldp-np.mean(fldc,axis=0) # this gets rid of the arbitrary timeseries-ness of it...results are very similar but not exact
fldda = np.squeeze(fldd-np.mean(fldd,axis=0))

if subset:
    fldda,msk=cutl.mask_region(fldda,lato,lono,'other',limsdict=limsdict)
    fldda1=fldda[~fldda.mask]
    fldda1 = fldda1.reshape( (oshape[0], len(lat),len(lon))  ) 
    fldda=fldda1




solverd=Eof(fldda,weights=wgts)
if docorr:
    eof1d=solverd.eofsAsCorrelation(neofs=enum)
    eof1d=eof1d[enum-1,...]
else:
    eof1d=solverd.eofsAsCovariance(neofs=enum)
    eof1d=eof1d[enum-1,...]

eigsd=solverd.eigenvalues()
vexpd=eigsd[enum-1]/eigsd.sum() *100 # percent variance explained


ax=axs[3]
cplt.kemmap(eof1d,lat,lon,type=type,title='Eof of diff',axis=ax,cmin=cmin,cmax=cmax)
ax.set_ylabel(str(np.round(vexpd)))

if printtofile:
    if docorr:
        fig.savefig(field + '_globEOF' + str(enum) + 'ascorr_' + sim + '_' + sea + '_' + substr + '_' + type + '2.pdf')
    else:
        fig.savefig(field + '_globEOF' + str(enum) + 'ascovar_' + sim + '_' + sea + '_' +substr + '_' + type + '2.pdf')
    

