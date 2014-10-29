

printtofile=True

plt.close('all')


sim = 'R3'
sea = 'DJF'

fuc,fup=con.build_filepathpair(sim,'su')

lat=cnc.getNCvar(fuc,'lat')
lon=cnc.getNCvar(fuc,'lon')

uc=cnc.getNCvar(fuc,'SU',timesel='0002-01-01,0121-12-31',seas=sea)
up=cnc.getNCvar(fup,'SU',timesel='0002-01-01,0121-12-31',seas=sea)

ucm=np.mean(uc,axis=0)
upm=np.mean(up,axis=0)
ud=upm-ucm

fvc,fvp=con.build_filepathpair(sim,'sv')

vc=cnc.getNCvar(fvc,'SV',timesel='0002-01-01,0121-12-31',seas=sea)
vp=cnc.getNCvar(fvp,'SV',timesel='0002-01-01,0121-12-31',seas=sea)

vcm=np.mean(vc,axis=0)
vpm=np.mean(vp,axis=0)
vd=vpm-vcm

lons,lats = np.meshgrid(lon,lat)

mag=np.sqrt((upm-ucm)**2+(vpm-vcm)**2)
plt.figure()
bm,cf = cplt.kemmap(mag,lat,lon,cmin=-1,cmax=1,type='nh')
bm.quiver(lons,lats,ud,vd,latlon=True)


ftc,ftp=con.build_filepathpair(sim,'st')

tc=cnc.getNCvar(ftc,'ST',timesel='0002-01-01,0121-12-31',seas=sea)
tp=cnc.getNCvar(ftp,'ST',timesel='0002-01-01,0121-12-31',seas=sea)

tcm=np.mean(tc,axis=0)
tpm=np.mean(tp,axis=0)

td=tpm-tcm

plt.figure()
bm,cf = cplt.kemmap(td,lat,lon,cmin=-1.5,cmax=1.5,type='nh')
bm.quiver(lons,lats,ud,vd,latlon=True)
if printtofile:
    plt.savefig('st_windvecs_' + sea + '_' + sim + '.png')


fslpc,fslpp=con.build_filepathpair(sim,'pmsl')


slpc=cnc.getNCvar(fslpc,'PMSL',timesel='0002-01-01,0121-12-31',seas=sea)
slpp=cnc.getNCvar(fslpp,'PMSL',timesel='0002-01-01,0121-12-31',seas=sea)

slpcm=np.mean(slpc,axis=0)
slppm=np.mean(slpp,axis=0)
slpd=slppm-slpcm

plt.figure()
bm,cf = cplt.kemmap(slpd,lat,lon,cmin=-2,cmax=2,type='nh',cmap='blue2red_20')
bm.quiver(lons,lats,ud,vd,latlon=True)
if printtofile:
    plt.savefig('pmsl_windvecs_' + sea + '_' + sim + '.png')

    # @@@ add reference vector
