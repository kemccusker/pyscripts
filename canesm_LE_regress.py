
import loadLE as le
import cccmaplots as cplt

printtofile=False

addobs=True # to scatter plot
addnat=False
addsims=True # add the idealized simulations. only good for DJF polar amp vs eurasia SAT

local=True


timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'
timeselall = '1979-01-01,2012-12-31'

ftype='fullts' # 'fullclimo' or 'climo' or 'fullts'
ensmean=False

seasp='DJF' # season of spatial field
sear='DJF' # season of regional avgs


# spatial field1 in color
leconvsp=1
fieldsp='tas'; ncfieldsp='tas'; compsp='Amon'; 

# spatial field2 in contours
leconvsp2=1
fieldsp2='zg50000.00'; ncfieldsp2='zg'; compsp2='Amon'
cminsp2=-10; cmaxsp2=10



# regional avg field 1
leconvr=-1 # this way, sea ice loss is linked with positive changes elsewhere
fieldr='sic'; ncfieldr='sic'; compr='OImon'; regionr='bksmori'

# regional avg field 2
leconvr2=1
fieldr2='tas'; ncfieldr2='tas'; compr2='Amon'; regionr2='eurasiamori'; leconvr2=-1 # so cooling=high heights


# what are the units of these regressions? @@
if fieldsp=='zg50000.00':
    cmin=-10; cmax=10
    cmin2=-10; cmax2=10
elif fieldsp=='tas':
    cmin=-1; cmax=1
    cmin2=-1; cmax2=1


fdictsp = {'field': fieldsp, 'ncfield': ncfieldsp, 'comp': compsp}
fdictsp2 = {'field': fieldsp2, 'ncfield': ncfieldsp2, 'comp': compsp2}
fdictr = {'field': fieldr+regionr, 'ncfield': ncfieldr, 'comp': compr}
fdictr2 = {'field': fieldr2+regionr2, 'ncfield': ncfieldr2, 'comp': compr2}

casename = 'historical'

lat=le.get_lat(local=local)
lon=le.get_lon(local=local)
nlat=len(lat); nlon=len(lon)

def load_field(fdict,casename,timesel,seas,ftype='fullts',conv=1,local=False,verb=False):
    
    """ 
        returns [numens x space.flat] or [numens]

    """

    ledat = le.load_LEdata(fdict,casename,timesel=timesel, 
                           rettype='ndarray',conv=conv,ftype=ftype,local=local,verb=verb)

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



casenames=('historical','historicalNat','historicalMisc')
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
    lepsear = load_field(fdictr,casename,timeselp,seasp,ftype=ftype,conv=leconvr)
    lesear = (lepsear-lecsear)/(lepsear-lecsear).std()

    lecsear2 = load_field(fdictr2,casename,timeselc,sear,ftype=ftype,conv=leconvr2)
    lepsear2 = load_field(fdictr2,casename,timeselp,seasp,ftype=ftype,conv=leconvr2)
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




# ====================== FIGURES ===============
printtofile=False

lons, lats = np.meshgrid(lon,lat)
cmlen=15.
incr = (cmaxsp2-cminsp2) / (cmlen)
conts = np.arange(cminsp2,cmaxsp2+incr,incr)

ttl1=seasp + ' regressions onto ' + sear + ' ' + fieldr+regionr
ttl2=seasp + ' regressions onto ' + sear + ' ' + fieldr2+regionr2

ttl1=ttl2=''

fig,axs=plt.subplots(1,2)
fig.set_size_inches(10,5)
fig.subplots_adjust(wspace=0.05)
ax=axs[0]
bm,pc=cplt.kemmap(fldsponfldr,lat,lon,type='nheur',axis=ax,cmin=cmin,cmax=cmax,
                  title=ttl1,suppcb=True,
                  panellab='a')
bm.contour(lons,lats,fldsp2onfldr,levels=conts,
           colors='k',linewidths=1,latlon=True)

ax=axs[1]
bm,pc=cplt.kemmap(fldsponfldr2,lat,lon,type='nheur',axis=ax,cmin=cmin,cmax=cmax,
                  title=ttl2,suppcb=True,
                  panellab='b')

bm.contour(lons,lats,fldsp2onfldr2,levels=conts,
           colors='k',linewidths=1,latlon=True)

cplt.add_colorbar(fig,pc,orientation='horizontal')

if printtofile:
    fig.savefig(fieldsp + '_' + fieldsp2 + seasp + \
                '_regresson_' + fieldr+regionr + '_' + fieldr2 + regionr2 + sear + '_allLE.pdf') 
