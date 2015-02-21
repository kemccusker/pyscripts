
import pandas as pd
import cccmaplots as cplt
import datetime as datetime
from netCDF4 import num2date
import sys
import loadLE as le

le=reload(le)

plt.close('all')
conv=1
# @@@ create a gif animation: convert -delay 10 -loop 0 *.png volcevo.gif

# regular zonal mean anomaly lims, zoom lims, and climo lims
cmint=-0.5; cmaxt=0.5 # height x time plots. good for ua, ta

#field='taZM'; ncfield='ta'; cmin=-2; cmax=2; cminz=-0.5; cmaxz=0.5; cminc=-60+273; cmaxc=40+273
#field='taZMpac'; ncfield='ta'; cmin=-2; cmax=2; cminz=-0.5; cmaxz=0.5; cminc=-60+273; cmaxc=40+273
#field='uaZM'; ncfield='ua'; cmin=-2; cmax=2; cminz=-0.5; cmaxz=0.5; cminc=-40; cmaxc=40
#field='uaZMpac'; ncfield='ua'; cmin=-2; cmax=2; cminz=-1; cmaxz=1; cminc=-40; cmaxc=40
field='zmpsi'; ncfield='ZMPSI'; cmin=-.05; cmax=.05; cminz=-.05; cmaxz=.05; cminc=-1; cmaxc=1; cmint=-0.1; cmaxt=0.1; conv=1e-11

printtofile=True
eps=sys.float_info.epsilon

fdict = {'field': field, 'ncfield': ncfield, 'comp': 'Amon'}
offset=5
duration=3
troplims=(-10,10)

if field=='zmpsi':
    divby=5.0 # to calc increments in climo contours
elif field=='uaZM':
    divby=10.0
else:
    divby=15.0

levels=np.arange(cminc,cmaxc+eps,(cmaxc-cminc)/divby) # 11 for uaZM
levelsn=np.arange(cminc,0,(0-cminc)/divby)
levelsp=np.arange(0,cmaxc+eps,(cmaxc-0)/divby)
levelsp=levelsp[1:] # skip 0

# Agung: Feb. 1963
# Chichon: Apr 1982
# Pinatubo: Jun 1991

# These time selections start with month of eruption
## timselagung ='1963-02-01,1965-12-31'
## timselchichon = '1982-04-01,1984-12-31'
## timselpinatubo = '1991-06-01,1993-12-31'

# These time selections start with beginning of year of eruption
#  (so can easily remove seasonal cycle)
# index into time series for eruption month, and index of eruption month (in a year)
timselagung ='1962-01-01,1968-12-31'; agstidx=13; ageridx=1 
timselchichon = '1981-01-01,1987-12-31'; chstidx=15; cheridx=3
timselpinatubo = '1990-01-01,1996-12-31'; pistidx=17; pieridx=5


# Need to remove each run's annual cycle: get seasonal cycle climo
natclimo = le.load_LEdata(fdict,'historicalNat',ftype='1950-2020_climo',rettype='ndarray',conv=conv)
#@@@histclimo = le.load_LEdata(fdict,'historical',ftype='1950-2020_climo',rettype='ndarray',conv=conv)

# Get months surrounding each volcano
agungdat=le.load_LEdata(fdict,'historicalNat',timesel=timselagung, rettype='ndarray',conv=conv)
chichondat=le.load_LEdata(fdict,'historicalNat',timesel=timselchichon, rettype='ndarray',conv=conv)
pinatubodat=le.load_LEdata(fdict,'historicalNat',timesel=timselpinatubo, rettype='ndarray',conv=conv)

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


def volcano_avg(volcdat,remclimo, eridx, coords,offset=offset, duration=duration):
    """ volcano_comp(volcdat,remclimo):
            For each ensemble member, remove its climo sea cycle,
            then select months to average.
            
            eridx: eruption index into timeseries
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

        # remove the monthly climo
        volcrem[nii,...] = vo-climot

    volcrem = volcrem.reshape((numens,ntime,n1,n2))

    # Now select the 9 months starting 5 months after eruption
    voslice=slice(eridx+offset,eridx+offset+duration)

    vcomptime=volcrem.mean(axis=0)
    vcomp = vcomptime[voslice,...].mean(axis=0)

    return vcomp


def volc_regevolution(volcdat,remclimo, eridx, coords,latlims):
    """ average over given latitudes and keep time dimension
        latlims: tuple of southern and northern limits for averaging

        note: this just works with zonal-average data
        
    """
    deg2rad=2*np.pi/180.0

    troplats = np.logical_and(lat>troplims[0],lat<troplims[1])
    tot=np.cos(lat[troplats]*deg2rad).sum()
    tropwgts=np.cos(lat[troplats]*deg2rad)/tot
    tropwgtst = np.tile(tropwgts,(nlev,1))

    # same as volcano_avg()
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

        # remove the monthly climo
        volcrem[nii,...] = vo-climot

    volcrem = volcrem.reshape((numens,ntime,n1,n2))

    # slice surrounding the eruption only
    # one year before, 5 years after
    voslice=slice(eridx-12,eridx+60)

    vcomptime=volcrem.mean(axis=0) # average over ensemble
    ertime=vcomptime.shape[0]
    #print ertime
    
    print tropwgtst.shape    
    #tropwgtst=np.tile(tropwgtst,(ertime,1,1))
    #print tropwgtst.shape
    
    vcomp = vcomptime[voslice,:,troplats]
    tropwgtst=np.tile(tropwgtst,(vcomp.shape[0],1,1))
    
    print vcomp.shape
    
    vcomp=(vcomp*tropwgtst).sum(axis=2)

    return vcomp


def plot_volcevolution(compev,lev,time,cmap='blue2red_20',cmin=None,cmax=None,axis=None,
                       addcontlines=False,levlim=100,suppylab=False,suppcb=False,title=''):
    """ plot_volcevolution(compev,lev,time,cmap='blue2red_20',cmin=None,cmax=None,axis=None,
                       addcontlines=False,levlim=100,suppylab=False,suppcb=False,title=''):
             
            Plot time x height figure
    """
    fld=compev
    
    #fig,ax=plt.subplots(1,1)
    times,levs = np.meshgrid(time,lev/100.)

    incmap = plt.cm.get_cmap(cmap)
    
    cmlen=float(incmap.N) # or: from __future__ import division

    if cmlen>200:
        cmlen = float(15) # ie. if using a built in colormap that is continuous, want smaller # of colors

        if cmin =='':
            # parameters for plot
            pparams = dict(cmap=incmap)
        
    else:
        
        incr = (cmax-cmin) /cmlen
        conts = np.arange(cmin,cmax+incr,incr)
        pparams = dict(cmap=incmap,levels=conts,vmin=cmin,vmax=cmax,extend='both')

    if axis != None:
        # what to do w/ axis?
        ax=axis
    else:
        ax=plt.gca()

    if cmin=='':
        cf = ax.contourf(times,levs,fld,**pparams) #cmlen autolevels@@removed
    else:
        cf = ax.contourf(times,levs,fld,**pparams)

    if addcontlines:
        if cmin=='':
            ax.contour(times,levs,fld,colors='.3') #cmlen autolevels@@removed
        else:
            ax.contour(times,levs,fld,levels=conts,colors='.3') # may have to make a dict for levels key

    ax.set_yscale('log')
    ax.set_yticks([1000,800, 500, 300, 100, 10])

    if suppylab==False:
        ax.set_yticklabels((1000,800,500,300,100,10))
    else:
        ax.set_yticklabels('')           

    ax.set_ylim(levlim,1000)
    ax.invert_yaxis()
    if suppylab==False:
        ax.set_ylabel('Pressure (hPa)')

    ax.set_title(title)

    if suppcb==False:
        #cbar = ax.colorbar(cf)
        
        #cbar_ax = plt.add_axes([.91,.15, .02,.7])
        plt.colorbar(cf) #pc,cax=cbar_ax)

    return cf

# ======== plot height with time ========

acompev = volc_regevolution(agungdat,natclimo,agstidx,(nlev,nlat),troplims)
times=np.arange(0,acompev.shape[0])
fig,ax=plt.subplots(1,1)
plot_volcevolution(acompev.T,lev,times,cmin=cmint,cmax=cmaxt,levlim=10,title='Agung',axis=ax)
ax.axvline(x=12,color='k',linestyle='--',linewidth=1.5)

ccompev = volc_regevolution(chichondat,natclimo,chstidx,(nlev,nlat),troplims)
fig,ax=plt.subplots(1,1)
plot_volcevolution(ccompev.T,lev,times,cmin=cmint,cmax=cmaxt,levlim=10,title='Chichon',axis=ax)
ax.axvline(x=12,color='k',linestyle='--',linewidth=1.5)

pcompev = volc_regevolution(pinatubodat,natclimo,pistidx,(nlev,nlat),troplims)
fig,ax=plt.subplots(1,1)
plot_volcevolution(pcompev.T,lev,times,cmin=cmint,cmax=cmaxt,levlim=10,title='Pinatubo',axis=ax)
ax.axvline(x=12,color='k',linestyle='--',linewidth=1.5)

fig,ax=plt.subplots(1,1)
fig.set_size_inches(12,4)
plot_volcevolution((acompev.T+ccompev.T+pcompev.T)/3.,lev,times,
                   cmin=cmint,cmax=cmaxt,levlim=10,title='Composite',axis=ax)
ax.axvline(x=12,color='k',linestyle='--',linewidth=1.5)
if printtofile:
    fig.savefig('volc_compositeintime_latlims' + str(troplims[0]) +'-' + str(troplims[1]) + '_' + field + '_' + '.pdf')

fig,ax=plt.subplots(1,1)
fig.set_size_inches(12,4)
plot_volcevolution((acompev.T+ccompev.T+pcompev.T)/3.,lev,times,
                   cmin=cmint,cmax=cmaxt,levlim=100,title='Composite',axis=ax)
ax.axvline(x=12,color='k',linestyle='--',linewidth=1.5)
if printtofile:
    fig.savefig('volc_compositeintime_latlims' + str(troplims[0]) +'-' + str(troplims[1]) + '_' + field + '_' + 'zoom.pdf')



acomp = volcano_avg(agungdat,natclimo, agstidx,(nlev,nlat))
aclim = natclimo.mean(axis=0) # average over ensemble
aclim = aclim[ageridx+offset:ageridx+offset+duration,...].mean(axis=0) # average over months
aclim = np.reshape(aclim,(nlev,nlat))


ccomp = volcano_avg(chichondat,natclimo, chstidx,(nlev,nlat))
cclim = natclimo.mean(axis=0)
cclim = cclim[cheridx+offset:cheridx+offset+duration,...].mean(axis=0)
cclim = np.reshape(cclim,(nlev,nlat))

pcomp = volcano_avg(pinatubodat,natclimo, pistidx,(nlev,nlat))
pclim = natclimo.mean(axis=0)
pclim = pclim[pieridx+offset:pieridx+offset+duration,...].mean(axis=0)
pclim = np.reshape(pclim,(nlev,nlat))


# #######plotting #####################


levlim=10
suff='off' + str(offset) + 'dur' + str(duration) + 'mo'

levs,lats=np.meshgrid(lat,lev/100.)

# WITH CLIMO CONTOURS
fig,axs=plt.subplots(1,3)
fig.set_size_inches(13,6)

ax=axs[0]
cplt.vert_plot(acomp,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='Agung (Feb 1963)',suppcb=True,
               levlim=levlim,cmap='blue2red_20')
#ax.contour(levs,lats,aclim,colors='0.6',linewidths=2,levels=levels)
ax.contour(levs,lats,aclim,colors='0.6',linewidths=1.5,levels=levelsp)
ax.contour(levs,lats,aclim,colors='0.6',linewidths=1.5,levels=levelsn,linestyles='--')
if field=='zmpsi':
    ax.set_yticks([1000,800,500,300,100])
    ax.set_yticklabels((1000,800,500,300,100))

ax=axs[1]
cplt.vert_plot(ccomp,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='Chichon (Apr 1982)',suppcb=True,
               suppylab=True,levlim=levlim,cmap='blue2red_20')
#ax.contour(levs,lats,cclim,colors='0.6',linewidths=2,levels=levels)
ax.contour(levs,lats,cclim,colors='0.6',linewidths=1.5,levels=levelsp)
ax.contour(levs,lats,cclim,colors='0.6',linewidths=1.5,levels=levelsn,linestyles='--')

ax=axs[2]
vp = cplt.vert_plot(pcomp,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='Pinatubo (Jun 1991)',suppcb=True,
                    suppylab=True,levlim=levlim,cmap='blue2red_20')
#ax.contour(levs,lats,pclim,colors='0.6',linewidths=2,levels=levels)
ax.contour(levs,lats,pclim,colors='0.6',linewidths=1.5,levels=levelsp)
ax.contour(levs,lats,pclim,colors='0.6',linewidths=1.5,levels=levelsn,linestyles='--')

cbar_ax = fig.add_axes([.25,0.07, 0.5, .02])
cbor='horizontal'
fig.colorbar(vp,cax=cbar_ax, orientation=cbor) 
if printtofile:
    fig.savefig('agung_chichon_pinatubo_' + field + '_' + suff + '_climocont.pdf')


# NO CLIMO CONTOURS
fig,axs=plt.subplots(1,3)
fig.set_size_inches(13,6)
ax=axs[0]
cplt.vert_plot(acomp,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='Agung (Feb 1963)',suppcb=True,
               levlim=levlim,addcontlines=True,cmap='blue2red_w20')

ax=axs[1]
cplt.vert_plot(ccomp,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='Chichon (Apr 1982)',suppcb=True,
               suppylab=True,levlim=levlim,addcontlines=True,cmap='blue2red_w20')

ax=axs[2]
vp = cplt.vert_plot(pcomp,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='Pinatubo (Jun 1991)',suppcb=True,
                    suppylab=True,levlim=levlim,addcontlines=True,cmap='blue2red_w20')

cbar_ax = fig.add_axes([.25,0.07, 0.5, .02])
cbor='horizontal'
fig.colorbar(vp,cax=cbar_ax, orientation=cbor) 
if printtofile:
    fig.savefig('agung_chichon_pinatubo_' + field + '_' + suff +'.pdf')

    
# COMPOSITE ================
#   WITH CONTOURS
fig,ax=plt.subplots(1,1)
cplt.vert_plot((acomp+ccomp+pcomp)/3.,lev,lat,axis=ax,cmin=cmin,cmax=cmax,
               title='volc composite',levlim=levlim,cmap='blue2red_20')
#ax.contour(levs,lats,(aclim+cclim+pclim)/3.,colors='0.6',linewidths=2,levels=levels)
ax.contour(levs,lats,(aclim+cclim+pclim)/3.,colors='0.6',linewidths=1.5,levels=levelsp)
ax.contour(levs,lats,(aclim+cclim+pclim)/3.,colors='0.6',linewidths=1.5,levels=levelsn,linestyles='--')

if printtofile:
    fig.savefig('volc_composite_' + field + '_' + suff + '_climocont.pdf')

#   NO CONTOURS
fig,ax=plt.subplots(1,1)
#fig.set_size_inches(12,6)
cplt.vert_plot((acomp+ccomp+pcomp)/3.,lev,lat,axis=ax,cmin=cmin,cmax=cmax,
               title='volc composite',addcontlines=True,levlim=levlim,cmap='blue2red_w20')
if printtofile:
    fig.savefig('volc_composite_' + field + '_' + suff  + '.pdf')


# COMPOSITE zoom in ==========
#   WITH CONTOURS
fig,ax=plt.subplots(1,1)
#fig.set_size_inches(12,6)
cplt.vert_plot((acomp+ccomp+pcomp)/3.,lev,lat,axis=ax,cmin=cminz,cmax=cmaxz,
               title='volc composite',levlim=100,cmap='blue2red_20')
#ax.contour(levs,lats,(aclim+cclim+pclim)/3.,colors='0.6',linewidths=2,levels=levels)
ax.contour(levs,lats,(aclim+cclim+pclim)/3.,colors='0.6',linewidths=1.5,levels=levelsp)
ax.contour(levs,lats,(aclim+cclim+pclim)/3.,colors='0.6',linewidths=1.5,levels=levelsn,linestyles='--')

ax.set_xlim(-30,30)
if printtofile:
    fig.savefig('volc_composite_' + field + '_' + suff + '_zoom_climocont.pdf')

#   NO CONTOURS
fig,ax=plt.subplots(1,1)
#fig.set_size_inches(12,6)
cplt.vert_plot((acomp+ccomp+pcomp)/3.,lev,lat,axis=ax,cmin=cminz,cmax=cmaxz,
               title='volc composite',levlim=100,addcontlines=True,cmap='blue2red_w20')

ax.set_xlim(-30,30)
if printtofile:
    fig.savefig('volc_composite_' + field + '_' + suff + '_zoom.pdf')


# =================== MMC notes:
#mmc=cnc.getNCvar(fname,'ZMPSI',seas='ANN')
#lat=cnc.getNCvar(fname,'lat')
#lev=cnc.getNCvar(fname,'plev')
#lats,levs = np.meshgrid(lat,lev)

#fig,ax=plt.subplots(1,1) 
#cnt = ax.contour(lats,levs/100,mmc.mean(axis=0))
#ax.clabel(cnt) # get values: cnt.cvalues. need to format
#ax.invert_yaxis()
#ax.set_yscale('log')





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

"""
