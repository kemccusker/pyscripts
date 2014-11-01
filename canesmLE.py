

basepath='/ra40/data/kem/CanSISE/CanESM2/LE/'
sims=('historical-r1','historical-r2','historical-r3','historical-r4','historical-r5')
ensnum=10

field='sic'; comp='OImon' # CMIP variable names
sea='ANN'
timesel=None
superii=0

years=np.arange(1950,2021)

simnhdt={}
simshdt={}
for sim in sims:

    ensnhdt={}
    ensshdt={}
    for eii in np.arange(1,ensnum+1):

        fname=basepath + sim + '/' + field + '/' + field + '_' + comp + '_CanESM2_' +\
               sim + '_r' + str(eii) + 'i1p1_195001-202012.nc'

        if superii==0:
            lat=cnc.getNCvar(fname,'lat')
            lon=cnc.getNCvar(fname,'lon')
            
        fld = cnc.getNCvar(fname,field,timesel=timesel,seas=sea)

        if field=='sic':
            ensnhdt[eii],ensshdt[eii] = cutl.calc_totseaicearea(fld/100.,lat,lon,isarea=False)
        else:
            print 'Only sic implemented so far. Do global mean?'

        superii+=1

    simnhdt[sim] = pd.DataFrame(ensnhdt,index=years)
    simshdt[sim] = pd.DataFrame(ensshdt,index=years)


col=0.0
fig,axs=plt.subplots(2,1,sharey=True)
ax=axs[0]
for sim in sims:
    df=simnhdt[sim]
    ax.plot(years,df,color=str(col))
    df=simshdt[sim]
    ax.plot(years,df,color=str(col))
    col+=0.2
ax.set_title(field + ' ' + sea)
col=0.0
ax=axs[1]
for sim in sims:
    df=simnhdt[sim]
    ax.plot(years,df.mean(axis=1),color=str(col))
    df=simshdt[sim]
    ax.plot(years,df.mean(axis=1),color=str(col))
    col+=0.2


# @@ now combine with cmip5 ensemble, and CESM if possible..
    
