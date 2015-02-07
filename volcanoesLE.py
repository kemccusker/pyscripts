
import pandas as pd
import cccmaplots as cplt
import datetime as datetime
from netCDF4 import num2date

import loadLE as le

le=reload(le)

plt.close('all')


#field='taZM'; ncfield='ta'; cmin=-2; cmax=2; cminz=-0.5; cmaxz=0.5
field='uaZM'; ncfield='ua'; cmin=-2; cmax=2; cminz=-0.5; cmaxz=0.5

printtofile=True

fdict = {'field': field, 'ncfield': ncfield, 'comp': 'Amon'}
offset=5
duration=9

# Agung: Feb. 1963
# Chichon: Apr 1982
# Pinatubo: Jul 1991

# These time selections start with month of eruption
## timselagung ='1963-02-01,1965-12-31'
## timselchichon = '1982-04-01,1984-12-31'
## timselpinatubo = '1991-07-01,1993-12-31'

# These time selections start with beginning of year of eruption
#  (so can easily remove seasonal cycle)
timselagung ='1963-01-01,1965-12-31'; agstidx=1 # index into time series for eruption month
timselchichon = '1982-01-01,1984-12-31'; chstidx=3
timselpinatubo = '1991-01-01,1993-12-31'; pistidx=6


# Need to remove each run's annual cycle: get seasonal cycle climo
natclimo = le.load_LEdata(fdict,'historicalNat',ftype='1950-2020_climo',rettype='ndarray')
histclimo = le.load_LEdata(fdict,'historical',ftype='1950-2020_climo',rettype='ndarray')

# Get months surrounding each volcano
agungdat=le.load_LEdata(fdict,'historicalNat',timesel=timselagung, rettype='ndarray')
chichondat=le.load_LEdata(fdict,'historicalNat',timesel=timselchichon, rettype='ndarray')
pinatubodat=le.load_LEdata(fdict,'historicalNat',timesel=timselpinatubo, rettype='ndarray')

flist=le.build_filenames(fdict,'historicalNat',ftype='fullts')
numens=len(flist)

fname = flist[0]

lev = cnc.getNCvar(fname,'plev')
lat = cnc.getNCvar(fname,'lat')
nlev=len(lev)
nlat=len(lat)

atime = cnc.getNCvar(fname,'time',timesel=timselagung)
adates = cnc.get_NCdates(fname,timesel=timselagung)


ctime = cnc.getNCvar(fname,'time',timselchichon)
ptime = cnc.getNCvar(fname,'time',timselpinatubo)

ntime=len(atime)
nyr=ntime/12.

# eventually also do this w/ the historical ensemble (not just historicalNat)
#   for wind only

#======= Agung
## anatrem = np.zeros(agungdat.shape)

## for nii in range(0,numens):
    
##     # tile the climo
##     climo = np.squeeze(natclimo[nii,...])
##     ag = np.squeeze(agungdat[nii,...])
    
##     nspace = climo.shape[1]
##     climot=np.tile(climo,(nyr,1))

##     anatrem[nii,...] = ag-climot

## anatrem = anatrem.reshape((numens,ntime,nlev,nlat))

## # Now select the 9 months starting 5 months after eruption
## agslice=slice(agstidx+5,agstidx+14)

## acomptime=anatrem.mean(axis=0)
## acomp = acomptime[agslice,...].mean(axis=0)


def volcano_avg(volcdat,remclimo, eridx, coords,offset=5, duration=9):
    """ volcano_comp(volcdat,remclimo):
            For each ensemble member, remove its climo sea cycle,
            then select months to average.
            
            ermon: eruption index into timeseries
            Default offset months is 5 months after eruption
            Default duration is 9 months

            returns the average of volcano time period over ensemble
    """
    volcrem = np.zeros(volcdat.shape)
    (numens,ntime,nlevlat) = volcdat.shape
    n1=coords[0]
    n2=coords[1]
    

    for nii in range(0,numens):

        # tile the climo
        climo = np.squeeze(remclimo[nii,...])
        vo = np.squeeze(volcdat[nii,...])

        nspace = climo.shape[1]
        climot=np.tile(climo,(nyr,1))

        volcrem[nii,...] = vo-climot

    volcrem = volcrem.reshape((numens,ntime,n1,n2))

    # Now select the 9 months starting 5 months after eruption
    voslice=slice(eridx+offset,eridx+offset+duration)

    vcomptime=volcrem.mean(axis=0)
    vcomp = vcomptime[voslice,...].mean(axis=0)

    return vcomp


acomp = volcano_avg(agungdat,natclimo, agstidx,(nlev,nlat))
aclim = natclimo.mean(axis=0)
aclim = aclim[agstidx+offset:agstidx+offset+duration,...].mean(axis=0)
aclim = np.reshape(aclim,(nlev,nlat))

ccomp = volcano_avg(chichondat,natclimo, chstidx,(nlev,nlat))
cclim = natclimo.mean(axis=0)
cclim = cclim[chstidx+offset:chstidx+offset+duration,...].mean(axis=0)
cclim = np.reshape(aclim,(nlev,nlat))

pcomp = volcano_avg(pinatubodat,natclimo, pistidx,(nlev,nlat))
pclim = natclimo.mean(axis=0)
pclim = pclim[pistidx+offset:pistidx+offset+duration,...].mean(axis=0)
pclim = np.reshape(aclim,(nlev,nlat))


## #======= Chichon
## cnatrem = np.zeros(chichondat.shape)

## for nii in range(0,numens):
    
##     # tile the climo
##     climo = np.squeeze(natclimo[nii,...])
##     ch = np.squeeze(chichondat[nii,...])
    
##     nspace = climo.shape[1]
##     climot=np.tile(climo,(nyr,1))

##     cnatrem[nii,...] = ch-climot

## cnatrem = cnatrem.reshape((numens,ntime,nlev,nlat))

## # Now select the 9 months starting 5 months after eruption
## chslice=slice(chstidx+5,chstidx+14)

## ccomptime=cnatrem.mean(axis=0)
## ccomp = ccomptime[chslice,...].mean(axis=0)



## #======= Chichon
## pnatrem = np.zeros(pinatubodat.shape)

## for nii in range(0,numens):
    
##     # tile the climo
##     climo = np.squeeze(natclimo[nii,...])
##     pi = np.squeeze(pinatubodat[nii,...])
    
##     nspace = climo.shape[1]
##     climot=np.tile(climo,(nyr,1))

##     pnatrem[nii,...] = pi-climot

## pnatrem = pnatrem.reshape((numens,ntime,nlev,nlat))

## # Now select the 9 months starting 5 months after eruption
## pislice=slice(pistidx+5,pistidx+14)

## pcomptime=pnatrem.mean(axis=0)
## pcomp = pcomptime[pislice,...].mean(axis=0)


# #######plotting #####################

levlim=10

levs,lats=np.meshgrid(lat,lev/100.)

fig,axs=plt.subplots(1,3)
fig.set_size_inches(13,6)

ax=axs[0]
cplt.vert_plot(acomp,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='agung',suppcb=True,
               levlim=levlim,addcontlines=True,cmap='blue2red_20')
ax.contour(levs,lats,aclim,colors='0.5',linewidths=2)

ax=axs[1]
cplt.vert_plot(ccomp,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='chichon',suppcb=True,
               suppylab=True,levlim=levlim,addcontlines=True,cmap='blue2red_20')
ax.contour(levs,lats,cclim,colors='0.5',linewidths=2)

ax=axs[2]
vp = cplt.vert_plot(pcomp,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='pinatubo',suppcb=True,
                    suppylab=True,levlim=levlim,addcontlines=True,cmap='blue2red_20')
ax.contour(levs,lats,pclim,colors='0.5',linewidths=2)

cbar_ax = fig.add_axes([.25,0.07, 0.5, .02])
cbor='horizontal'
fig.colorbar(vp,cax=cbar_ax, orientation=cbor) 
if printtofile:
    fig.savefig('agung_chichon_pinatubo_' + field + '_climocont.pdf')


fig,ax=plt.subplots(1,1)
#fig.set_size_inches(12,6)
cplt.vert_plot((acomp+ccomp+pcomp)/3.,lev,lat,axis=ax,cmin=cmin,cmax=cmax,
               title='volc composite',levlim=levlim,addcontlines=True,cmap='blue2red_20')
ax.contour(levs,lats,(aclim+cclim+pclim)/3.,colors='0.5',linewidths=2)

if printtofile:
    fig.savefig('volc_composite_' + field + '_climocont.pdf')


# zoom in

fig,ax=plt.subplots(1,1)
#fig.set_size_inches(12,6)
cplt.vert_plot((acomp+ccomp+pcomp)/3.,lev,lat,axis=ax,cmin=cminz,cmax=cmaxz,
               title='volc composite',levlim=100,addcontlines=True,cmap='blue2red_20')
ax.contour(levs,lats,(aclim+cclim+pclim)/3.,colors='0.5',linewidths=2)

ax.set_xlim(-30,30)
if printtofile:
    fig.savefig('volc_composite_' + field + '_zoom_climocont.pdf')

"""
#adt=le.load_LEdata(fdict,'historicalNat',timesel=timselagung)
#apan=pd.Panel(adt,major_axis=adates)
#amon = apan.groupby(apan.major_axis.month) # this returns a pandas.core.groupby.PanelGroupBy object. @@ what to do w/ it??

# @@ Note that these are currently not weighted by month:
agungsr=pd.Series(agungdt)
agmean = agungsr.mean(axis=0) # avg over ensemble, so it's months by height by last
# John's avg is 9 months starting with 5th month after eruption
agcomp = np.mean(agmean[4:13,:,:],axis=0) # agung composite

chichonsr=pd.Series(chichondt)
chmean = chichonsr.mean(axis=0) # avg over ensemble, so it's months by height by last
# John's avg is 9 months starting with 5th month after eruption
chcomp = np.mean(chmean[4:13,:,:],axis=0) # chichon composite

pinatubosr=pd.Series(pinatubodt)
pimean = pinatubosr.mean(axis=0) # avg over ensemble, so it's months by height by last
# John's avg is 9 months starting with 5th month after eruption
picomp = np.mean(pimean[4:13,:,:],axis=0) # pinatubo composite

plt.figure(); cplt.vert_plot(agcomp,lev,lat)
plt.figure(); cplt.vert_plot(chcomp,lev,lat)
plt.figure(); cplt.vert_plot(picomp,lev,lat)


levlim=10
cmin=-4; cmax=4

fig,axs=plt.subplots(1,3)
fig.set_size_inches(12,6)
ax=axs[0]
cplt.vert_plot(agcomp-bamean,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='agung',suppcb=True,levlim=levlim)
ax.set_xlim(-30,30)
ax=axs[1]
cplt.vert_plot(chcomp-bamean,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='chichon',suppcb=True,levlim=levlim)
ax.set_xlim(-30,30)
ax=axs[2]
vp = cplt.vert_plot(picomp-bamean,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='pinatubo',suppcb=True,levlim=levlim)
ax.set_xlim(-30,30)
# @@@ add full composite

cbar_ax = fig.add_axes([.25,0.07, 0.5, .02])
cbor='horizontal'
fig.colorbar(vp,cax=cbar_ax, orientation=cbor) 


"""
