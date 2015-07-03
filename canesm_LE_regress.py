
import loadLE as le
import cccmaplots as cplt
import cccmautils as cutl
import canesm_LE_general as leg
import loadmodeldata as lmd
import pandas as pd
import constants as con

printtofile=False

local=True


timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'
timeselall = '1979-01-01,2012-12-31'

ftype='fullts' # 'fullclimo' or 'climo' or 'fullts'
# this is not implemented here. but would do regressions in time on LE ensemble avg
ensmean=False 

seasp='DJF' # season of spatial field
sear='SON' #'DJF' # season of regional avgs


# spatial field1 in color
leconvsp=1
fieldsp='tas'; ncfieldsp='tas'; compsp='Amon'; 
cminsp=-1; cmaxsp=1 # for colors
cminspa=-1; cmaxspa=1 # for colors for AGCM

# spatial field2 in contours
leconvsp2=1
fieldsp2='zg50000.00'; ncfieldsp2='zg'; compsp2='Amon'
cminsp2=-10; cmaxsp2=10 # to calc contour interval
cminsp2a=-10; cmaxsp2a=10 # to calc contour interval for AGCM


# regional avg field 1
leconvr=-1 # this way, sea ice loss is linked with positive changes elsewhere
fieldr='sic'; ncfieldr='sic'; compr='OImon'; regionr='bksmori'

# regional avg field 2
leconvr2=-1 # so cooling=high heights
fieldr2='tas'; ncfieldr2='tas'; compr2='Amon'; regionr2='eurasiamori';



fdictsp = {'field': fieldsp, 'ncfield': ncfieldsp, 'comp': compsp}
fdictsp2 = {'field': fieldsp2, 'ncfield': ncfieldsp2, 'comp': compsp2}
fdictr = {'field': fieldr+regionr, 'ncfield': ncfieldr, 'comp': compr}
fdictr2 = {'field': fieldr2+regionr2, 'ncfield': ncfieldr2, 'comp': compr2}

lat=le.get_lat(local=local)
lon=le.get_lon(local=local)
nlat=len(lat); nlon=len(lon)



def load_field(fdict,casename,timesel,seas,ftype='fullts',conv=1,local=False,verb=False):
    
    """ 
        returns [numens x space.flat] or [numens]

    """

    ledat = le.load_LEdata(fdict,casename,timesel=timesel, 
                           rettype='ndarray',conv=conv,ftype=ftype,local=local,verb=verb)

    print '@@@ ledat.shape ' + str(ledat.shape) # why does the 3d data have space flattened already...

    # time needs to be first dimension
    try:
        if ledat.ndim==2:
            ledat = ledat.T
        elif ledat.ndim==3:
            ledat = np.transpose(ledat,(1,0,2))
        else:
            print 'Loaded data is not 2 or 3 dimensions. Do not understand.'
            raise Exception
    except:
        raise

    lesea = cutl.seasonalize_monthlyts(ledat,season=seas).mean(axis=0)  # numens x space.flat

    return lesea


def load_agcmfield(field,sims,seas,conv=1,region=None,subsampyrs=11,styears=None):
    """ loads subsampled agcm data from specified simulations

              number of total samples will be determined by length of all sim data
                      and numyrs (e.g. (ntime / numyrs)*numsims)

              returns subsamp
                      styears
    """
    threed=False

    simconv=1
    if field=='tas': simfield='st'; simncfield='ST'
    elif field=='zg50000.00': simfield='gz50000'; simncfield='PHI'; simconv=1/con.get_g()
    elif field=='sia': simfield='sicn'; simncfield='SICN'; print '@@ danger, sia actually sicn average'
    elif field=='sic': simfield='sicn'; simncfield='SICN'
    else: print 'cannot addsims for ' + field

    
    if region!=None:
        simflddf = pd.DataFrame(lmd.loaddata((simfield,),sims,ncfields=(simncfield,), timefreq=seas, 
                                             region=region))*simconv
    else:
        # assume threed b/c no region given.
        threed=True
        simflddf = lmd.loaddata((simfield,),sims,ncfields=(simncfield,), timefreq=seas, 
                                region=region,rettype='ndarray')*simconv


    subsamp,styearsss = leg.subsamp_sims(simflddf,numyrs=subsampyrs,styears=styears,threed=threed)

    return subsamp,styearsss

def slopemap(inr,insp,dims):
    """   
          inr is 1D [time or numens]
          insp is 2D [time or numens x space.flat]
          dims are a tuple of dims to reshape space to (nlat,nlon) 

          returns slopemap [dims]           
    """ 
    slope,intercept = np.polyfit(inr,insp, 1)
    slopemap = slope.reshape(dims)

    return slopemap



#casenames=('historical','historicalNat','historicalMisc')
casenames=('historical',)

for casename in casenames:
    # # SPATIAL DATA
    print 'SPATIAL'

    lecseasp = load_field(fdictsp,casename,timeselc,seasp,ftype=ftype,conv=leconvsp)
    lepseasp = load_field(fdictsp,casename,timeselp,seasp,ftype=ftype,conv=leconvsp)
    leseasp = lepseasp-lecseasp

    lecseasp2 = load_field(fdictsp2,casename,timeselc,seasp,ftype=ftype,conv=leconvsp2)
    lepseasp2 = load_field(fdictsp2,casename,timeselp,seasp,ftype=ftype,conv=leconvsp2)
    leseasp2 = lepseasp2-lecseasp2


    # # 1D DATA
    print '1D'
    lecsear = load_field(fdictr,casename,timeselc,sear,ftype=ftype,conv=leconvr)
    lepsear = load_field(fdictr,casename,timeselp,sear,ftype=ftype,conv=leconvr)
    lesear = (lepsear-lecsear)/(lepsear-lecsear).std()

    lecsear2 = load_field(fdictr2,casename,timeselc,sear,ftype=ftype,conv=leconvr2)
    lepsear2 = load_field(fdictr2,casename,timeselp,sear,ftype=ftype,conv=leconvr2)
    lesear2 = (lepsear2-lecsear2)/(lepsear2-lecsear2).std()


    if casename!='historical': # really just has to be not the first loop thru
        tmp = np.vstack((tmp,leseasp))
        tmp2 = np.vstack((tmp2,leseasp2))

        tmpr = np.hstack((tmpr,lesear))
        tmpr2 = np.hstack((tmpr2,lesear2))
    else:
        tmp = leseasp
        tmp2 = leseasp2
        tmpr = lesear
        tmpr2 = lesear2

leseasp=tmp
leseasp2=tmp2
lesear=tmpr
lesear2=tmpr2


# calc regression slopes
fldsponfldr = slopemap(lesear,leseasp,(nlat,nlon)) # SAT regress on regSIC
fldsp2onfldr = slopemap(lesear,leseasp2,(nlat,nlon)) # Z500 regress on regSIC

fldsponfldr2 = slopemap(lesear2,leseasp,(nlat,nlon)) # SAT regress on regSAT
fldsp2onfldr2 = slopemap(lesear2,leseasp2,(nlat,nlon)) # Z500 regress on regSAT

# === AGCM ==========

simsE=('E1','E2','E3','E4','E5'); 
sims=('R1','R2','R3','R4','R5');

asseasp,styears = load_agcmfield(fieldsp,sims,seasp) # already anomalies
asseasp2,styears = load_agcmfield(fieldsp2,sims,seasp,styears=styears) # already anomalies

assear,styears = load_agcmfield(fieldr,sims,sear,styears=styears,region=regionr) # already anomalies
assear2,styears = load_agcmfield(fieldr2,sims,sear,styears=styears,region=regionr2) # already anomalies
assear=assear / assear.std()
assear2=assear2 / assear2.std()

# Make panel (d) with E sims:
aesseasp,styearse = load_agcmfield(fieldsp,simsE,seasp) 
aesseasp2,styearse = load_agcmfield(fieldsp2,simsE,seasp,styears=styearse) 
aessear2,styearse = load_agcmfield(fieldr2,simsE,sear,styears=styearse,region=regionr2) 
aessear2=aessear2 / aessear2.std()


(anens,anlat,anlon)=asseasp.shape
rshape=(anens,anlat*anlon)
# calc regression slopes: multiply regional avgs by -1 to get colors/signs right.
#                         this was accounted for in LE data upon loading.
asponfldr = slopemap(assear*-1,asseasp.reshape(rshape),(anlat,anlon)) # SAT regress on regSIC
asp2onfldr = slopemap(assear*-1,asseasp2.reshape(rshape),(anlat,anlon)) # Z500 regress on regSIC
asponfldr2 = slopemap(assear2*-1,asseasp.reshape(rshape),(anlat,anlon)) # SAT regress on regSAT
asp2onfldr2 = slopemap(assear2*-1,asseasp2.reshape(rshape),(anlat,anlon)) # Z500 regress on regSAT

aesponfldr2 = slopemap(aessear2*-1,aesseasp.reshape(rshape),(anlat,anlon)) # SAT regress on regSAT
aesp2onfldr2 = slopemap(aessear2*-1,aesseasp2.reshape(rshape),(anlat,anlon)) # Z500 regress on regSAT



# ====================== FIGURES ===============
printtofile=False

lons, lats = np.meshgrid(lon,lat)
cmlen=15.
incr = (cmaxsp2-cminsp2) / (cmlen)
conts = np.arange(cminsp2,cmaxsp2+incr,incr)

ttl1=seasp + ' regress on ' + sear + ' BKS SIC'
ttl2=seasp + ' regress on ' + sear + ' Eur SAT' 

#ttl1=ttl2=''

fig,axs=plt.subplots(1,2)
fig.set_size_inches(10,5)
fig.subplots_adjust(wspace=0.05)
ax=axs[0]
bm,pc=cplt.kemmap(fldsponfldr,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title=ttl1,suppcb=True,
                  panellab='a',lcol='0.2')
bm.contour(lons,lats,fldsp2onfldr,levels=conts,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[1]
bm,pc=cplt.kemmap(fldsponfldr2,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title=ttl2,suppcb=True,
                  panellab='b',lcol='0.2')

bm.contour(lons,lats,fldsp2onfldr2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)

cplt.add_colorbar(fig,pc,orientation='horizontal')

if printtofile:
    fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                '_regresson_' + fieldr+regionr + '_' + fieldr2 + regionr2 + sear + '.pdf') 


# ========== AGCM

alat=con.get_t63lat(); alon=con.get_t63lon()
alons, alats = np.meshgrid(alon,alat)

cmlen=15.
incra = (cmaxsp2a-cminsp2a) / (cmlen)
contsa = np.arange(cminsp2a,cmaxsp2a+incra,incra)


fig,axs=plt.subplots(1,2)
fig.set_size_inches(10,5)
fig.subplots_adjust(wspace=0.05)
ax=axs[0]
bm,pc=cplt.kemmap(asponfldr,alat,alon,type='nheur',axis=ax,cmin=cminspa,cmax=cmaxspa,
                  title=ttl1,suppcb=True,
                  panellab='a',lcol='0.2')
bm.contour(alons,alats,asp2onfldr,levels=contsa,
           colors='0.5',linewidths=1,latlon=True)

ax=axs[1]
bm,pc=cplt.kemmap(asponfldr2,alat,alon,type='nheur',axis=ax,cmin=cminspa,cmax=cmaxspa,
                  title=ttl2,suppcb=True,
                  panellab='b',lcol='0.2')

bm.contour(alons,alats,asp2onfldr2,levels=contsa,
           colors='0.5',linewidths=1,latlon=True)

cplt.add_colorbar(fig,pc,orientation='horizontal')

if printtofile:
    fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                '_regresson_' + fieldr+regionr + '_' + fieldr2 + regionr2 + sear + '_AGCM.pdf') 


# ========== 4 panel fig =========

fig,axs=plt.subplots(2,2)
fig.set_size_inches(10,10)
fig.subplots_adjust(wspace=0.05,hspace=0.05)
ax=axs[0,0]
bm,pc=cplt.kemmap(fldsponfldr,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title=ttl1,suppcb=True,
                  panellab='a',lcol='0.2')
bm.contour(lons,lats,fldsp2onfldr,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
ax.set_ylabel('CGCM')

ax=axs[0,1]
bm,pc=cplt.kemmap(fldsponfldr2,lat,lon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title=ttl2,suppcb=True,
                  panellab='b',lcol='0.2')

bm.contour(lons,lats,fldsp2onfldr2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)
cplt.add_colorbar(fig,pc,orientation='horizontal')

ax=axs[1,0]
bm,pc=cplt.kemmap(asponfldr,alat,alon,type='nheur',axis=ax,cmin=cminspa,cmax=cmaxspa,
                  title='',suppcb=True,
                  panellab='c',lcol='0.2')
bm.contour(alons,alats,asp2onfldr,levels=contsa,
           colors='0.5',linewidths=1,latlon=True)
ax.set_ylabel('AGCM var ICE')

ax=axs[1,1]
bm,pc=cplt.kemmap(asponfldr2,alat,alon,type='nheur',axis=ax,cmin=cminspa,cmax=cmaxspa,
                  title='',suppcb=True,
                  panellab='d',lcol='0.2')

bm.contour(alons,alats,asp2onfldr2,levels=contsa,
           colors='0.5',linewidths=1,latlon=True)

cplt.add_colorbar(fig,pc,orientation='horizontal')

if printtofile:
    fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                '_regresson_' + fieldr+regionr + '_' + fieldr2 +\
                regionr2 + sear + '_CGCMAGCM.pdf') 

# PANEL D with E sims & residual: --------------
fig,axs=plt.subplots(1,2)
fig.set_size_inches(10,5)
ax=axs[0]
#ax=axs[2,1]
#bbox=ax.get_position()

#ax=fig.add_axes([bbox.x0,bbox.y0-bbox.height,bbox.width,bbox.height])
bm,pc=cplt.kemmap(aesponfldr2,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title=ttl2,suppcb=True,
                  panellab='e',lcol='0.2')
bm.contour(alons,alats,aesp2onfldr2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)

#cplt.add_colorbar(fig,pc,orientation='horizontal')
ax.set_ylabel('AGCM const ICE')


#fig,axs=plt.subplots(1,1)
#fig.set_size_inches(5,5)
ax=axs[1]
bm,pc=cplt.kemmap(asponfldr2-aesponfldr2,alat,alon,type='nheur',axis=ax,cmin=cminsp,cmax=cmaxsp,
                  title='residual: var - const',suppcb=True,
                  panellab='f',lcol='0.2')
bm.contour(alons,alats,asp2onfldr2-aesp2onfldr2,levels=conts,
           colors='0.5',linewidths=1,latlon=True)

cplt.add_colorbar(fig,pc,orientation='horizontal')

if printtofile:
    fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                '_regresson_' + fieldr+regionr + '_' + fieldr2 +\
                regionr2 + sear + '_AGCM_E_R-Eresid.pdf') 
