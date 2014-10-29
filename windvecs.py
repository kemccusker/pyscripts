
# @@@@ Double check the wind vectors are correct
#   with the projection. Over the pacific, the zonal wind direction
#   seems reversed when compared to SLP (ie the winds go clockwise
#   around a low! They are correct on the Atlantic side though, where
#   the winds go clockwise around a high.
#  Also, when plot robinson, it all looks right. So, potentially
#  half the globe has uvel flipped?   10/29/14
from mpl_toolkits.basemap import Basemap, shiftgrid


printtofile=True

plt.close('all')

ufield='su'; uncfield='SU'
vfield='sv'; vncfield='SV'

level=50000 # for when threed is True
threed=False

sim = 'ENSE'
sea = 'DJF' #1 #'DJF'


if threed:
    uncfield='U'
    vncfield='V'
    ufield='u' + str(level)
    vfield='v' + str(level)
    level=str(level/100) # for printing figs
else:
    level=''
    
fuc,fup=con.build_filepathpair(sim,ufield)

lat=cnc.getNCvar(fuc,'lat')
lon=cnc.getNCvar(fuc,'lon')

uc=cnc.getNCvar(fuc,uncfield,timesel='0002-01-01,0121-12-31',seas=sea)
up=cnc.getNCvar(fup,uncfield,timesel='0002-01-01,0121-12-31',seas=sea)

ucm=np.mean(uc,axis=0)
upm=np.mean(up,axis=0)
ud=upm-ucm

fvc,fvp=con.build_filepathpair(sim,vfield)

vc=cnc.getNCvar(fvc,vncfield,timesel='0002-01-01,0121-12-31',seas=sea)
vp=cnc.getNCvar(fvp,vncfield,timesel='0002-01-01,0121-12-31',seas=sea)

vcm=np.mean(vc,axis=0)
vpm=np.mean(vp,axis=0)
vd=vpm-vcm

#  From: http://matplotlib.org/basemap/users/examples.html
# shift grid so it goes from -180 to 180 (instead of 0 to 360
# in longitude).  Otherwise, interpolation is messed up.
ugrid,newlons = shiftgrid(180.,ud,lon,start=False)
vgrid,newlons = shiftgrid(180.,vd,lon,start=False)

spud=ud[...,::2]
spvd=vd[...,::2]


lons,lats = np.meshgrid(lon,lat)
splons=lons[...,::2] # sparse lons
splats=lats[...,::2] # sparse lons


# plot this one with all vectors for reference
mag=np.sqrt((upm-ucm)**2+(vpm-vcm)**2)
plt.figure()
bm,cf = cplt.kemmap(mag,lat,lon,cmin=-1,cmax=1,type='nh')
# transform vectors to projection grid.
uproj,vproj,xx,yy = bm.transform_vector(ugrid,vgrid,newlons,lat,46,54,returnxy=True,masked=True)
qu=bm.quiver(xx,yy,uproj,vproj)#,scale=700)
#qu=bm.quiver(lons,lats,ud,vd,latlon=True)
plt.quiverkey(qu,0.1,0.1,1,'1 m/s',labelpos='W')

ftc,ftp=con.build_filepathpair(sim,'st')

tc=cnc.getNCvar(ftc,'ST',timesel='0002-01-01,0121-12-31',seas=sea)
tp=cnc.getNCvar(ftp,'ST',timesel='0002-01-01,0121-12-31',seas=sea)

tcm=np.mean(tc,axis=0)
tpm=np.mean(tp,axis=0)

td=tpm-tcm

plt.figure()
bm,cf = cplt.kemmap(td,lat,lon,cmin=-1.5,cmax=1.5,type='nh')
# transform vectors to projection grid.
uproj,vproj,xx,yy = bm.transform_vector(ugrid,vgrid,newlons,lat,46,54,returnxy=True,masked=True)
qu=bm.quiver(xx,yy,uproj,vproj)#,scale=700)
#qu=bm.quiver(splons,splats,spud,spvd,latlon=True)
plt.quiverkey(qu,0.1,0.1,1,'1 m/s',labelpos='W')
if printtofile:
    plt.savefig('st_windvecs' + level + '_' + str(sea) + '_' + sim + '.png')


fslpc,fslpp=con.build_filepathpair(sim,'pmsl')


slpc=cnc.getNCvar(fslpc,'PMSL',timesel='0002-01-01,0121-12-31',seas=sea)
slpp=cnc.getNCvar(fslpp,'PMSL',timesel='0002-01-01,0121-12-31',seas=sea)

slpcm=np.mean(slpc,axis=0)
slppm=np.mean(slpp,axis=0)
slpd=slppm-slpcm

plt.figure()
bm,cf = cplt.kemmap(slpd,lat,lon,cmin=-2,cmax=2,type='nh',cmap='blue2red_20')
# transform vectors to projection grid.
uproj,vproj,xx,yy = bm.transform_vector(ugrid,vgrid,newlons,lat,46,54,returnxy=True,masked=True)
qu=bm.quiver(xx,yy,uproj,vproj)#,scale=700)
#qu=bm.quiver(splons,splats,spud,spvd,latlon=True)
plt.quiverkey(qu,0.1,0.1,1,'1 m/s',labelpos='W')
if printtofile:
    plt.savefig('pmsl_windvecs' + level + '_' + str(sea) + '_' + sim + '.png')



fuc,fup=con.build_filepathpair(sim,ufield)
fvc,fvp=con.build_filepathpair(sim,vfield)

monstr=con.get_mon()
months=[12,1,2,3]

fig,axs=plt.subplots(2,2) # D, J, F, M
fig.set_size_inches(10,8)
for aii,ax in enumerate(axs.flat):

    sea=months[aii]
    #print sea

    uc=cnc.getNCvar(fuc,uncfield,timesel='0002-01-01,0121-12-31',seas=sea)
    up=cnc.getNCvar(fup,uncfield,timesel='0002-01-01,0121-12-31',seas=sea)

    ucm=np.mean(uc,axis=0)
    upm=np.mean(up,axis=0)
    ud=upm-ucm

    vc=cnc.getNCvar(fvc,vncfield,timesel='0002-01-01,0121-12-31',seas=sea)
    vp=cnc.getNCvar(fvp,vncfield,timesel='0002-01-01,0121-12-31',seas=sea)

    vcm=np.mean(vc,axis=0)
    vpm=np.mean(vp,axis=0)
    vd=vpm-vcm

    #  From: http://matplotlib.org/basemap/users/examples.html
    # shift grid so it goes from -180 to 180 (instead of 0 to 360
    # in longitude).  Otherwise, interpolation is messed up.
    ugrid,newlons = shiftgrid(180.,ud,lon,start=False)
    vgrid,newlons = shiftgrid(180.,vd,lon,start=False)

    slpc=cnc.getNCvar(fslpc,'PMSL',timesel='0002-01-01,0121-12-31',seas=sea)
    slpp=cnc.getNCvar(fslpp,'PMSL',timesel='0002-01-01,0121-12-31',seas=sea)

    slpcm=np.mean(slpc,axis=0)
    slppm=np.mean(slpp,axis=0)
    slpd=slppm-slpcm
    tstat,pval = sp.stats.ttest_ind(slpp,slpc,axis=0)

    bm,cf = cplt.kemmap(slpd,lat,lon,cmin=-2,cmax=2,type='nh',cmap='blue2red_20',axis=ax,suppcb=True)
    cplt.addtsigm(bm,pval,lat,lon,type='cont')
    # transform vectors to projection grid.
    uproj,vproj,xx,yy = bm.transform_vector(ugrid,vgrid,newlons,lat,31,31,returnxy=True,masked=True)
    qu=bm.quiver(xx,yy,uproj,vproj)#,scale=700)
    #qu=bm.quiver(splons,splats,spud,spvd,latlon=True)
    #if aii==1: # not sure it's the same ref vec each plot...
    ax.quiverkey(qu,0.1,0.1,1,'1 m/s',labelpos='W')
        
    ax.set_title(monstr[sea-1])
cbar_ax = fig.add_axes([.91,.25, .02,.5])
fig.colorbar(cf,cax=cbar_ax)
if printtofile:
    plt.savefig('pmsl_windvecs' + level + '_DJFMsbplt_' + sim + '.png')



fig,axs=plt.subplots(2,2) # D, J, F, M
fig.set_size_inches(10,8)
for aii,ax in enumerate(axs.flat):

    sea=months[aii]
    #print sea

    uc=cnc.getNCvar(fuc,uncfield,timesel='0002-01-01,0121-12-31',seas=sea)
    up=cnc.getNCvar(fup,uncfield,timesel='0002-01-01,0121-12-31',seas=sea)

    ucm=np.mean(uc,axis=0)
    upm=np.mean(up,axis=0)
    ud=upm-ucm

    vc=cnc.getNCvar(fvc,vncfield,timesel='0002-01-01,0121-12-31',seas=sea)
    vp=cnc.getNCvar(fvp,vncfield,timesel='0002-01-01,0121-12-31',seas=sea)

    vcm=np.mean(vc,axis=0)
    vpm=np.mean(vp,axis=0)
    vd=vpm-vcm

    #  From: http://matplotlib.org/basemap/users/examples.html
    # shift grid so it goes from -180 to 180 (instead of 0 to 360
    # in longitude).  Otherwise, interpolation is messed up.
    ugrid,newlons = shiftgrid(180.,ud,lon,start=False)
    vgrid,newlons = shiftgrid(180.,vd,lon,start=False)

    tc=cnc.getNCvar(ftc,'ST',timesel='0002-01-01,0121-12-31',seas=sea)
    tp=cnc.getNCvar(ftp,'ST',timesel='0002-01-01,0121-12-31',seas=sea)

    tcm=np.mean(tc,axis=0)
    tpm=np.mean(tp,axis=0)
    td=tpm-tcm
    tstat,pval = sp.stats.ttest_ind(tp,tc,axis=0)

    bm,cf = cplt.kemmap(td,lat,lon,cmin=-1.5,cmax=1.5,type='nh',axis=ax,suppcb=True)
    cplt.addtsigm(bm,pval,lat,lon,type='cont')
    # transform vectors to projection grid.
    uproj,vproj,xx,yy = bm.transform_vector(ugrid,vgrid,newlons,lat,31,31,returnxy=True,masked=True)
    qu=bm.quiver(xx,yy,uproj,vproj)
    #qu=bm.quiver(splons,splats,spud,spvd,latlon=True)
    #if aii==1: # not sure it's the same ref vec each plot...
    ax.quiverkey(qu,0.1,0.1,1,'1 m/s',labelpos='W')
        
    ax.set_title(monstr[sea-1])
cbar_ax = fig.add_axes([.91,.25, .02,.5])
fig.colorbar(cf,cax=cbar_ax)
if printtofile:
    plt.savefig('st_windvecs' + level + '_DJFMsbplt_' + sim + '.png')






"""# ===================================
#  try: http://matplotlib.org/basemap/users/examples.html

plt.figure()
bm,cf = cplt.kemmap(slpd,lat,lon,cmin=-2,cmax=2,type='nh',cmap='blue2red_20')

## # plot wind vectors on projection grid.
## # first, shift grid so it goes from -180 to 180 (instead of 0 to 360
## # in longitude).  Otherwise, interpolation is messed up.
## ugrid,newlons = shiftgrid(180.,ud,lon,start=False)
## vgrid,newlons = shiftgrid(180.,vd,lon,start=False)
# transform vectors to projection grid.
uproj,vproj,xx,yy = bm.transform_vector(ugrid,vgrid,newlons,lat,46,54,returnxy=True,masked=True)

qu=bm.quiver(xx,yy,uproj,vproj)#,scale=700)
plt.quiverkey(qu,0.1,0.1,1,'1 m/s',labelpos='W')


# ====================================== """


## plt.figure()
## bm,cf = cplt.kemmap(slpd,lat,lon,cmin=-2,cmax=2,cmap='blue2red_20')
## qu=bm.quiver(splons,splats,spud,spvd,latlon=True)
## plt.quiverkey(qu,0.1,0.1,1,'1 m/s',labelpos='W')
