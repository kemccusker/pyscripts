
import cccmautils as cutl
import pandas as pd
import scipy as sp
import scipy.stats

printtofile=False
plt.close('all')

# @@ Look at Southern Ocean and Southern Jet


basepath='/ra40/data/kem/CanSISE/CanESM2/LE/'
sims=('historical-r1','historical-r2','historical-r3','historical-r4','historical-r5')
ensnum=10

#field='sic'; comp='OImon' # CMIP variable names
field='tas'; comp='Amon'
reg2=True # if sea ice, include southern hem? if other, include regional mean?
region= 'etroppac'

sea='ANN'
timesel=None
timesel2 = '2002-01-01,2012-12-31'
superii=0
regstr=''


if reg2:
    if field=='sic':
        regstr='sh'
    else:
        regstr=region
    
years=np.arange(1950,2021)
if sea=='DJF':
    years=years[:-1]

simnhdt={}
simshdt={}
simgmdt={}
simrmdt={}

simtrnddt={}
for sim in sims:

    ensnhdt={}
    ensshdt={}
    ensgmdt={}
    ensrmdt={}

    ensseldt={}
    enstrnddt={}
    for eii in np.arange(1,ensnum+1):

        fname=basepath + sim + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
               sim + '_r' + str(eii) + 'i1p1_195001-202012.nc'

        if superii==0:
            lat=cnc.getNCvar(fname,'lat')
            lon=cnc.getNCvar(fname,'lon')
            
        fld = cnc.getNCvar(fname,field,timesel=timesel,seas=sea)
        # select a given time period and keep maps
        fldsel = cnc.getNCvar(fname,field,timesel=timesel2,seas=sea)
        nt = fldsel.shape[0]
        nlon=fldsel.shape[2]
        nlat=fldsel.shape[1]
        fldre = fldsel.reshape((nt,nlon*nlat))

        xx=np.arange(0,nt)
        
        #ensseldt[eii] = fldsel

        if field=='sic':
            ensnhdt[eii],ensshdt[eii] = cutl.calc_totseaicearea(fld/100.,lat,lon,isarea=False)
        else:
            #print 'Only sic implemented so far. Do global mean?'
            ensgmdt[eii] = cutl.global_mean_areawgted3d(fld,lat,lon)
            ensrmdt[eii] = cutl.calc_regmean(fld,lat,lon,region)

    
        #slope[eii], intercept, r_value, p_value, std_err = sp.stats.linregress(xx,dat) # not good for 3d data?
        slope,intercept = np.polyfit(xx,fldre,1) # supposedly can do w/ higher dims?
        
        enstrnddt[eii] = slope #.reshape((nlat,nlon)) # reshape later @@
        superii+=1

    if field=='sic':
        simnhdt[sim] = pd.DataFrame(ensnhdt,index=years)
        simshdt[sim] = pd.DataFrame(ensshdt,index=years)
    else:
        simgmdt[sim] = pd.DataFrame(ensgmdt,index=years)
        simrmdt[sim] = pd.DataFrame(ensrmdt,index=years)

    # calc time mean and trend of ensseldt. Then each ens member is a map of trends
    tmpdf = pd.DataFrame(enstrnddt)#,index=years)
    
    #simseldt[sim] = trenddf
    simtrnddt[sim] = tmpdf


if field=='sic':
    plotdt = simnhdt # NH ice
    plotdt2 = simshdt # SH ice
else:
    plotdt = simgmdt # global mean
    plotdt2 = simrmdt # regional mean


col=0.0
fig,axs=plt.subplots(2,1,sharey=True)
ax=axs[0]
for sim in sims:
    df=plotdt[sim]
    ax.plot(years,df-df.mean(axis=0),color=str(col))
    if reg2: #field=='sic':
        df=plotdt2[sim]
        ax.plot(years,df-df.mean(axis=0),color=str(col),linestyle='--')
    col+=0.2
ax.set_title(field + ' ' + str(sea) + ' ' + regstr)

col=0.0
ax=axs[1]
for sim in sims:
    df=plotdt[sim]
    ax.plot(years,df.mean(axis=1)-df.mean(axis=1).mean(axis=0),color=str(col))
    if reg2: #field=='sic':
        df=plotdt2[sim]
        ax.plot(years,df.mean(axis=1)-df.mean(axis=1).mean(axis=0),color=str(col),linestyle='--')
    col+=0.2
ax.set_xlabel('mean of groups')

if printtofile:
    fig.savefig(field + 'LE_remtm' + str(sea) + regstr + '_timeseries.pdf') # time mean removed


# plot mean trend map for each "group"
for sim in sims:
    trnds = simtrnddt[sim]
    groupm = trnds.mean(axis=1)
    groupm = groupm.reshape((nlat,nlon))

    plt.figure()
    mh = cplt.kemmap(groupm*nt,lat,lon,cmin=-1,cmax=1)


   
# @@ now combine with cmip5 ensemble, and CESM if possible..
    
