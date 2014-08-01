"""
   plot_fluxbar.py
   5/8/2014
   Plot bar graph of arctic fluxes, plus seacycle of SST/SIC
"""
import constants as con
import cccmautils as cutl
import numpy.ma as ma

plt.close("all")
plt.ion()

con = reload(con)
cutl = reload(cutl)


printtofile=1


casenamep = con.casenamep2()
casenamec = con.casenamec()

casenamepb = con.casenameph()
casenamecb = con.casenamech()

timeselb = '0002-01-01,0121-12-31' # don't need anymore. timesels are the same.
timstrb = '001-121'
timstrpb = timstrb

timesel = '0002-01-01,0121-12-31'
timstr = '001-111'
timstrp = timstr

if casenamec == 'kemhadctl':
    timesel = '0002-01-01,0121-12-31'
    timstr = '001-121'
    timstrp = timstr

model = 'CanAM4'

plat = platform.system()

if plat == 'Darwin':  # means I'm on my mac
    basepath = '/Users/kelly/CCCma/CanSISE/RUNS/' # @@ needs updating
    subdir = '/'
else:  # on linux workstation in Vic
    basepath = '/home/rkm/work/DATA/' + model + '/'
    subdir = '/ts/'


# get sea ice concentration for averaging purposes
field = 'sicn'
fnamec = basepath + casenamec + subdir + casenamec + '_' + field + '_' + timstr + '_ts.nc'
sicnc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)
sicnc,siccstd = cutl.climatologize(sicnc) # use to mask the flux region

fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstr + '_ts.nc'
sicnp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)
sicnp,sicpstd = cutl.climatologize(sicnp) # use to mask the flux region

# get second set of sims:
fnamecb = basepath + casenamecb + subdir + casenamecb + '_' + field + '_' + timstrb + '_ts.nc'
sicncb = cnc.getNCvar(fnamecb,field.upper(),timesel=timeselb)
sicncb,siccstdb = cutl.climatologize(sicncb) # use to mask the flux region

fnamepb = basepath + casenamepb + subdir + casenamepb + '_' + field + '_' + timstrb + '_ts.nc'
sicnpb = cnc.getNCvar(fnamepb,field.upper(),timesel=timeselb)
sicnpb,sicpbstd = cutl.climatologize(sicnpb) # use to mask the flux region


lat = cnc.getNCvar(fnamec,'lat')
lon = cnc.getNCvar(fnamec,'lon')

latlim = 40

fluxes = 'hfl','hfs','flg' # LH, SH, LWdown
fluxdict = {}
fluxdictb = {}

# prepare the flux data by averaging, etc 
for field in fluxes:

     fnamec = basepath + casenamec + subdir + casenamec + '_' + field + '_' + timstr + '_ts.nc'
     fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstrp + '_ts.nc'

     fldc,cstd = cutl.climatologize(cnc.getNCvar(fnamec,field.upper(),timesel=timesel))
     fldp,pstd = cutl.climatologize(cnc.getNCvar(fnamep,field.upper(),timesel=timesel))

     fldc = ma.masked_where(sicnc<.10,fldc) # masking out non-ice region (as P&M 2014 JClim)
     fldp = ma.masked_where(sicnc<.10,fldp)

     
     # now calc area-weighted mean
     cellareas = cutl.calc_cellareas(lat,lon,repeat=fldc.shape)   
     cellareas = ma.masked_where(sicnc<.10,cellareas)
     totarea = np.sum(np.sum(cellareas[:,lat>latlim,:],axis=2),axis=1)
     totarea = np.tile(totarea,(cellareas.shape[1],cellareas.shape[2],1))
     totarea = np.transpose(totarea,(2,0,1))
     #cellwgts = cellareas/totarea

     fldcam = np.sum(np.sum(fldc[:,lat>latlim,:]*(cellareas[:,lat>latlim,:]/totarea[:,lat>latlim,:]),axis=2),axis=1)
     fldpam = np.sum(np.sum(fldp[:,lat>latlim,:]*(cellareas[:,lat>latlim,:]/totarea[:,lat>latlim,:]),axis=2),axis=1)
     
     fluxdict[field] = fldpam - fldcam


     fnamecb = basepath + casenamecb + subdir + casenamecb + '_' + field + '_' + timstrb + '_ts.nc'
     fnamepb = basepath + casenamepb + subdir + casenamepb + '_' + field + '_' + timstrpb + '_ts.nc'

     fldcb,cbstd = cutl.climatologize(cnc.getNCvar(fnamecb,field.upper(),timesel=timeselb))
     fldpb,pbstd = cutl.climatologize(cnc.getNCvar(fnamepb,field.upper(),timesel=timeselb))

     fldcb = ma.masked_where(sicncb<.10,fldcb)
     fldpb = ma.masked_where(sicncb<.10,fldpb)

      # now calc area-weighted mean
     cellareasb = cutl.calc_cellareas(lat,lon,repeat=fldcb.shape)   
     cellareasb = ma.masked_where(sicncb<.10,cellareasb)
     totareab = np.sum(np.sum(cellareasb[:,lat>latlim,:],axis=2),axis=1)
     totareab = np.tile(totareab,(cellareasb.shape[1],cellareasb.shape[2],1))
     totareab = np.transpose(totareab,(2,0,1))
     #cellwgtsb = cellareasb/totareab

     fldcbam = np.sum(np.sum(fldcb[:,lat>latlim,:]*(cellareasb[:,lat>latlim,:]/totareab[:,lat>latlim,:]),axis=2),axis=1)
     fldpbam = np.sum(np.sum(fldpb[:,lat>latlim,:]*(cellareasb[:,lat>latlim,:]/totareab[:,lat>latlim,:]),axis=2),axis=1)
     
     fluxdictb[field] = fldpbam - fldcbam


# now calc SIA and avg SST
areas = cutl.calc_cellareas(lat,lon,repeat=sicnc.shape)
siac = np.sum( np.sum(sicnc[:,lat>latlim,:]*areas[:,lat>latlim,:],2),1)
siap = np.sum(np.sum(sicnp[:,lat>latlim,:]*areas[:,lat>latlim,:],2),1)

areasb = cutl.calc_cellareas(lat,lon,repeat=sicncb.shape)
siacb = np.sum( np.sum(sicncb[:,lat>latlim,:]*areasb[:,lat>latlim,:],2),1)
siapb = np.sum(np.sum(sicnpb[:,lat>latlim,:]*areasb[:,lat>latlim,:],2),1)


# all this SST stuff is wrong: difference is tiny b/w pert and ctl 5/8/2014 @@
# @@ prob using the wrong masks...

field = 'gt'
fnamec = basepath + casenamec + subdir + casenamec + '_' + field + '_' + timstr + '_ts.nc'
gtc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)
gtc,gtcstd = cutl.climatologize(gtc)

fnamep = basepath + casenamep + subdir + casenamep + '_' + field + '_' + timstr + '_ts.nc'
gtp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)
gtp,gtstd = cutl.climatologize(gtp)
gtp[sicnp>.15] = 271.2-273 # this is what the model would read
gtc[sicnc>.15] = 271.2-273
gtp = ma.masked_where(sicnc<.10,gtp)
gtc = ma.masked_where(sicnc<.10,gtc)

gtcam = np.sum(np.sum(gtc[:,lat>latlim,:]*(cellareas[:,lat>latlim,:]/totarea[:,lat>latlim,:]),axis=2),axis=1)
gtpam = np.sum(np.sum(gtp[:,lat>latlim,:]*(cellareas[:,lat>latlim,:]/totarea[:,lat>latlim,:]),axis=2),axis=1)


net = fluxdict['hfl']+fluxdict['hfs']+fluxdict['flg']*-1
netb = fluxdictb['hfl']+fluxdictb['hfs']+fluxdictb['flg']*-1

darkolivegreen1 = np.array([202, 255, 112])/255 # terrible
darkolivegreen3 = np.array([162, 205, 90])/255.
darkseagreen = np.array([143, 188, 143])/255.
darkseagreen4 = np.array([105, 139, 105])/255.
dodgerblue = np.array([30, 144, 255])/255. 
orangered4 = np.array([139, 37, 0])/255.



wi=.2
N=12 # 12 months, one var for now
ind = np.arange(1,N+1)
incr=0

fig,ax = plt.subplots()
fig.set_size_inches(9,5)

rects1 = ax.bar(ind+incr,fluxdict['hfl'],width=wi,color='.5')
incr=incr+wi
rects2 = ax.bar(ind+incr,fluxdict['hfs'],width=wi,color='.3')
incr=incr+wi
rects3 = ax.bar(ind+incr,fluxdict['flg']*-1,width=wi,color='k')

incr=incr+wi
rects4 = ax.bar(ind+incr,net,width=wi,color='g')

ax.plot(np.arange(1.5,13.5),((siap-siac)/1e11)*-1,color='0.7',linewidth=3)
ax.plot((.5,13),(0,0),'k')
ax.set_xlim(.5,13)
ax.set_ylim(-1,20)
ax.set_ylabel('W/m2 (bars) and SIA/1e11*-1')
ax.set_xlabel('Month')
ax.set_title(casenamep + ' - ' + casenamec)
#ax.legend( (rects1['hfl'][0],rects1['hfs'][0],rects1['flg'][0]),fluxes,'upper left',prop=fontP)
if printtofile:
    fig.savefig('FluxesDIFF_' + casenamep + '_v_' + casenamec + '_barbymonth_withSIAdiff.pdf')


# just do Net flux
wi=.5
N=12 # 12 months, one var for now
ind = np.arange(.75,N)
incr=0

fig,ax = plt.subplots()
fig.set_size_inches(9,5)

rects = ax.bar(ind+incr,net,width=wi,color='.5')

ax.plot(np.arange(1,13),((siap-siac)/1e11)*-1,color='0.7',linewidth=3)
# SST avg is wrong 5/8/2014 @@ 
#ax.plot(np.arange(1.5,13.5),(gtpam-gtcam)*10,color='b',linewidth=3)
ax.plot((.5,13),(0,0),'k')
ax.set_xlim(.5,13)
ax.set_ylim(-1,20)
ax.set_ylabel('W/m2 (bars) and SIA/1e11*-1')
ax.set_xlabel('Month')
ax.set_title(casenamep + ' - ' + casenamec)
if printtofile:
    fig.savefig('NetFluxDIFF_' + casenamep + '_v_' + casenamec + '_barbymonth_withSIAdiff.pdf')


# BOTH SETS OF SIMULATIONS: just do Net flux
wi=.3
N=12 # 12 months, one var for now
ind = np.arange(.75,N)
incr=0

fig,ax = plt.subplots()
fig.set_size_inches(9,5)

rects = ax.bar(ind+incr,net,width=wi,color='.5')
rectsb = ax.bar(ind+wi+incr,netb,width=wi,color=orangered4,alpha=.7)

ax.plot(np.arange(1,13),((siap-siac)/1e11)*-1,color='0.7',linewidth=3)
ax.plot(np.arange(1,13),((siapb-siacb)/1e11)*-1,color=orangered4,linewidth=3)

# SST avg is wrong 5/8/2014 @@ 
#ax.plot(np.arange(1.5,13.5),(gtpam-gtcam)*10,color='b',linewidth=3)
ax.plot((.5,13),(0,0),'k')
ax.set_xlim(.5,13)
ax.set_ylim(-1,20)
ax.set_ylabel('W/m2 (bars) and SIA/1e11*-1')
ax.set_xlabel('Month')
#ax.set_title(casenamep + ' - ' + casenamec)
if printtofile:
    fig.savefig('NetFluxDIFF_' + casenamep + '_' + casenamepb + '_v_ctrls_barbymonth_withSIAdiff.pdf')
