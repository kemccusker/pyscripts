"""
  scatter_regions.py
      This script makes scatter plots of regional means. Can do two
      variables, two regions, two seasons. There are figures of just the simulation mean
      values with confidence intervals, and figures with all regional averages
      in time with regression lines. 10/10/2014
"""

import cccmautils as cutl
import constants as con
import cccmaplots as cplt
import cccmacmaps as ccm
import pandas as pd

printtofile=True
plt.close('all')

conv1=1; conv2=1

plotscatter=True

#field1='st'; ncfield1='ST'
#field1='sia'; ncfield1='SICN'
field1='gz50000'; ncfield1='PHI'; conv1=1/con.get_g()
#field1='pmsl'; ncfield1='PMSL'
region1='bksmori' #'bksmori' #'polcap65'
sea1='ND' #'DJF'

field2='st'; ncfield2='ST'
#field2='sia'; ncfield2='SICN'
#field2='pmsl'; ncfield2='PMSL'
#field2='gz50000'; ncfield2='PHI'; conv2=1/con.get_g()
region2= 'eurasia' #'eurasiamori'
sea2='ND' #'DJF'

sims = ('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5','HAD','NSIDC','ENS','ENSE')
TOT = ('R1','R2','R3','R4','R5')
ANT = ('E1','E2','E3','E4','E5')

siglevel=0.05
flddregdt= {}
flddregdt2 = {}
flddregtsdt= {}
flddregtsdt2 = {}
ciregdt = {}
ciregdt2 = {}

ntime=120
# if one of the scatter vars is winter, need to remove a year from the other one
#  if it's not winter too
if sea1 in ('DJF','NDJ') or sea2 in ('DJF','NDJ'):
    ntime=ntime-1
    
nsim=len(TOT)
timesel='0002-01-01,0121-12-31'
alltotcr1=np.zeros(ntime*nsim)
alltotcr2=np.zeros(ntime*nsim)
alltotpr1=np.zeros(ntime*nsim)
alltotpr2=np.zeros(ntime*nsim)
alltotr1=np.zeros(ntime*nsim)
alltotr2=np.zeros(ntime*nsim)

allantcr1=np.zeros(ntime*nsim)
allantcr2=np.zeros(ntime*nsim)
allantpr1=np.zeros(ntime*nsim)
allantpr2=np.zeros(ntime*nsim)
allantr1=np.zeros(ntime*nsim)
allantr2=np.zeros(ntime*nsim)
tallii=0 # index to keep track of time in accumulated TOT ensemble
aallii=0 # index to keep track of time in accumulated ANT ensemble
for sim in sims:

    fnamec,fnamep=con.build_filepathpair(sim,field1)

    fldc=cnc.getNCvar(fnamec,ncfield1,timesel=timesel,seas=sea1)*conv1
    fldp=cnc.getNCvar(fnamep,ncfield1,timesel=timesel,seas=sea1)*conv1


    lat=cnc.getNCvar(fnamec,'lat')
    lon=cnc.getNCvar(fnamec,'lon')

    if field1=='sia':
        fldcreg=cutl.calc_regtotseaicearea(fldc,lat,lon,region1) # isarea=False
        fldpreg=cutl.calc_regtotseaicearea(fldp,lat,lon,region1) # isarea=False
    else:
        fldcreg=cutl.calc_regmean(fldc,lat,lon,region1)
        fldpreg=cutl.calc_regmean(fldp,lat,lon,region1)
        
    flddregdt[sim] = np.mean(fldpreg-fldcreg,axis=0) # time mean
    if sea2 in ('DJF','NDJ') and sea1 not in ('DJF','NDJ'): # have to shorten other timeseries
        flddregtsdt[sim] = fldpreg[:-1,...]-np.mean(fldcreg[:-1,...],axis=0) # anomaly timeseries from ctl mean
    else:
        flddregtsdt[sim] = fldpreg-np.mean(fldcreg,axis=0)

    if field2 != field1:
        fnamec,fnamep=con.build_filepathpair(sim,field2)
        fldc2=cnc.getNCvar(fnamec,ncfield2,timesel=timesel,seas=sea2)*conv2
        fldp2=cnc.getNCvar(fnamep,ncfield2,timesel=timesel,seas=sea2)*conv2
    elif sea2 != sea1:
        fldc2=cnc.getNCvar(fnamec,ncfield2,timesel=timesel,seas=sea2)*conv2
        fldp2=cnc.getNCvar(fnamep,ncfield2,timesel=timesel,seas=sea2)*conv2
    else:
        fldc2=fldc
        fldp2=fldp

    if field2=='sia':
        fldcreg2=cutl.calc_regtotseaicearea(fldc,lat,lon,region2) # isarea=False
        fldpreg2=cutl.calc_regtotseaicearea(fldp,lat,lon,region2) # isarea=False
    else:
        fldcreg2=cutl.calc_regmean(fldc2,lat,lon,region2)
        fldpreg2=cutl.calc_regmean(fldp2,lat,lon,region2)
        
    flddregdt2[sim] = np.mean(fldpreg2-fldcreg2,axis=0)
    if sea2 in ('DJF','NDJ') and sea1 not in ('DJF','NDJ'): # have to shorten other timeseries
        flddregtsdt2[sim] = fldpreg2[:-1,...]-np.mean(fldcreg2[:-1,...],axis=0)
    else:
        flddregtsdt2[sim] = fldpreg2-np.mean(fldcreg2,axis=0) # anomaly timeseries from ctl mean

    # @@@ there are files with concated ensembles now: e.g. sim name kemctl1enscat, diffname ENSCAT
    if sim in TOT:
        #print 'concat tot ens members'
        alltotcr1[tallii:tallii+ntime] = fldcreg # control reg1
        alltotcr2[tallii:tallii+ntime] = fldcreg2 # control reg2
        alltotpr1[tallii:tallii+ntime] = fldpreg # pert reg1
        alltotpr2[tallii:tallii+ntime] = fldpreg2 # pert reg2
        alltotr1[tallii:tallii+ntime] = flddregtsdt[sim] # diff reg1
        alltotr2[tallii:tallii+ntime] = flddregtsdt2[sim] # diff reg2
        tallii+=ntime
        #totser.append(pd.Series(flddregtsdt))
        #totser2.append(pd.Series(flddregtsdt2))
    elif sim in ANT:
        #print 'concat ant ens members'
        allantcr1[aallii:aallii+ntime] = fldcreg
        allantcr2[aallii:aallii+ntime] = fldcreg2
        allantpr1[aallii:aallii+ntime] = fldpreg
        allantpr2[aallii:aallii+ntime] = fldpreg2
        allantr1[aallii:aallii+ntime] = flddregtsdt[sim]
        allantr2[aallii:aallii+ntime] = flddregtsdt2[sim]
        aallii+=ntime
        
    # calculate confidence interval
    # double-check the scale setting
    ciregdt[sim] = sp.stats.t.interval(1-siglevel,len(fldpreg)-1,loc=np.mean(fldpreg,axis=0)-np.mean(fldcreg,axis=0),
                             scale=np.std(fldpreg-fldcreg,axis=0)/np.sqrt(len(fldpreg)))
    ciregdt2[sim] = sp.stats.t.interval(1-siglevel,len(fldpreg2)-1,loc=np.mean(fldpreg2,axis=0)-np.mean(fldcreg2,axis=0),
                             scale=np.std(fldpreg2-fldcreg2,axis=0)/np.sqrt(len(fldpreg2)))


# @@@ add confidence int for allant and alltot (means have to concat control and pert separately as well)
# calculate confidence interval
# double-check the scale setting
ciregalltotr1 = sp.stats.t.interval(1-siglevel,len(alltotpr1)-1,loc=np.mean(alltotpr1,axis=0)-np.mean(alltotcr1,axis=0),
                         scale=np.std(alltotpr1-alltotcr1,axis=0)/np.sqrt(len(alltotpr1)))
ciregalltotr2 = sp.stats.t.interval(1-siglevel,len(alltotpr2)-1,loc=np.mean(alltotpr2,axis=0)-np.mean(alltotcr2,axis=0),
                         scale=np.std(alltotpr2-alltotcr2,axis=0)/np.sqrt(len(alltotpr2)))
ciregallantr1 = sp.stats.t.interval(1-siglevel,len(allantpr1)-1,loc=np.mean(allantpr1,axis=0)-np.mean(allantcr1,axis=0),
                         scale=np.std(allantpr1-allantcr1,axis=0)/np.sqrt(len(allantpr1)))
ciregallantr2 = sp.stats.t.interval(1-siglevel,len(allantpr2)-1,loc=np.mean(allantpr2,axis=0)-np.mean(allantcr2,axis=0),
                         scale=np.std(allantpr2-allantcr2,axis=0)/np.sqrt(len(allantpr2)))


if plotscatter:
    # SCATTER SUPER ENSEMBLES IN TIME
    # super ensembles
    fig,ax= plt.subplots(1,1) #plt.figure();
    plt.scatter(allantr1,allantr2,color='r')
    plt.scatter(alltotr1,alltotr2,color='k')
    antr1std=allantr1.std()
    antr2std=allantr2.std()
    totr1std=alltotr1.std()
    totr2std=alltotr2.std()
    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5)
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5)


    # SCATTER SUPER ENS CTL and PERT separately
    #pltc1=allantcr1; pltc2=allantcr2
    #pltp1=allantpr1; pltp2=allantpr2
    # or as anomaly
    pltc1=allantcr1-allantcr1.mean(); pltc2=allantcr2-allantcr2.mean()
    pltp1=allantpr1-allantcr1.mean(); pltp2=allantpr2-allantcr2.mean()
    
    fig,axs= plt.subplots(1,2,sharex=True,sharey=True) #plt.figure();
    fig.set_size_inches(9,5)
    ax=axs[0]
    ax.scatter(pltc1,pltc2,color='k')
    ax.scatter(pltp1,pltp2,color='b')    
    ax.set_title('ANT')
    ax.legend(('CTL','PERT'),loc='best',fancybox=True,framealpha=0.5)
    ax.set_xlabel(field1 + region1 + str(sea1))
    ax.set_ylabel(field2 + region2 + str(sea2))
    
    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5)
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5)
                 
    onex=np.linspace(axxlims[0],axxlims[1])    
    mm, bb, rvalc1, pval, std_err = sp.stats.linregress(pltc1,pltc2)
    ax.plot(onex,mm*onex + bb, color='k',linewidth=2)
    mm, bb, rvalp1, pval, std_err = sp.stats.linregress(pltp1,pltp2)
    ax.plot(onex,mm*onex + bb, color='b',linewidth=2)
    ## axylims = ax.get_ylim()
    ## axxlims = ax.get_xlim()
    ## ax.annotate('$Rctl$= ' + '$%.2f$'%(rvalc) + ',$Rctl$='+'$%.2f$'%(rvalp),
    ##             xy=(axxlims[0]+.1*axxlims[1], axylims[1]-.1*axylims[1]),
    ##             xycoords='data')
    
    ax=axs[1]
    #pltc1=alltotcr1; pltc2=alltotcr2
    #pltp1=alltotpr1; pltp2=alltotpr2
    # or as anomaly
    pltc1=alltotcr1-alltotcr1.mean(); pltc2=alltotcr2-alltotcr2.mean()
    pltp1=alltotpr1-alltotcr1.mean(); pltp2=alltotpr2-alltotcr2.mean()
    
    ax.scatter(pltc1,pltc2,color='k')
    ax.scatter(pltp1,pltp2,color='b')
    ax.set_title('TOT')
    ax.legend(('CTL','PERT'),loc='best',fancybox=True,framealpha=0.5)
    ax.set_xlabel(field1 + region1 + str(sea1))
    #ax.set_ylabel(field2 + region2 + str(sea2))
    
    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5)
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5)
    onex=np.linspace(axxlims[0],axxlims[1])    
    mm, bb, rvalc, pval, std_err = sp.stats.linregress(pltc1,pltc2)
    ax.plot(onex,mm*onex + bb, color='k',linewidth=2)
    mm, bb, rvalp, pval, std_err = sp.stats.linregress(pltp1,pltp2)
    ax.plot(onex,mm*onex + bb, color='b',linewidth=2)
    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    ax.annotate('$Rctl$= ' + '$%.2f$'%(rvalc) + ',$Rpt$='+ '$%.2f$'%(rvalp),
                xy=(axxlims[0]+.1*axxlims[1], axylims[1]-.1*axylims[1]),
                xycoords='data')
    
    ax=axs[0] # go back and annotate first panel
    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    ax.annotate('$Rctl$= ' + '$%.2f$'%(rvalc1) + ',$Rpt$='+ '$%.2f$'%(rvalp1),
                xy=(axxlims[0]+.1*axxlims[1], axylims[1]-.1*axylims[1]),
                xycoords='data')
    if printtofile:
        fig.savefig('scatterregress_' + field1 + region1 + str(sea1) + '_v_' + field2 + region2 + str(sea2) + '_ANTTOTsuperens_ctlpertsbplt.pdf')
    

    #SCATTER SUPER ENSEMBLE MEANS with conf interval
    fig,ax= plt.subplots(1,1) #plt.figure();
    ant=plt.scatter(allantr1.mean(),allantr2.mean(),color='r',marker='s',s=8**2)
    tot=plt.scatter(alltotr1.mean(),alltotr2.mean(),color='k',marker='s',s=8**2)
    plt.plot(ciregalltotr1,(alltotr2.mean(),alltotr2.mean()),
            color='k',linewidth=2,marker='_',markersize=6)
    plt.plot((alltotr1.mean(),alltotr1.mean()),ciregalltotr2,
            color='k',linewidth=2,marker='_',markersize=6)
    plt.plot(ciregallantr1,(allantr2.mean(),allantr2.mean()),
            color='r',linewidth=2,marker='_',markersize=6)
    plt.plot((allantr1.mean(),allantr1.mean()),ciregallantr2,
            color='r',linewidth=2,marker='_',markersize=6)
    plt.xlabel(str(sea1) + ' ' + field1 + ' ' + region1)
    plt.ylabel(str(sea2) + ' ' + field2 + ' ' + region2)
    plt.legend((ant,tot),('ANT','TOT'),loc='best',fancybox=True,framealpha=0.5)
    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5)
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5)
    if printtofile:
        plt.savefig('scatterCI_' + field1 + region1 + str(sea1) + '_v_' + field2 + region2 + str(sea2) + '_ANTTOTsuperens.pdf')



    ser1 = pd.Series(flddregdt)
    ser2 = pd.Series(flddregdt2)
    serts1=pd.Series(flddregtsdt) 
    serts2=pd.Series(flddregtsdt2)
    serci1 = pd.Series(ciregdt)
    serci2 = pd.Series(ciregdt2)


    df = pd.DataFrame([ser1,ser2],index=(region1,region2))
    dfts = pd.DataFrame([serts1,serts2],index=(region1,region2))
    cidf = pd.DataFrame([serci1,serci2],index=(region1,region2))
    print df


    cd = ccm.get_colordict()

    # SCATTER PLOT OF SIM MEANS
    fig,ax = plt.subplots(1) #plt.figure()
    plt.scatter(df.filter(regex='R').values[0],df.filter(regex='R').values[1],color='k',marker='s',s=8**2)
    plt.scatter(df.filter(regex='E').values[0],df.filter(regex='E').values[1],color='0.5',marker='o',s=8**2)
    plt.scatter(df['HAD'].values[0],df['HAD'].values[1],color=cd['HAD'],marker='s',s=8**2)
    plt.scatter(df['NSIDC'].values[0],df['NSIDC'].values[1],color=cd['NSIDC'],marker='s',s=8**2)
    plt.scatter(df['ENS'].values[0],df['ENS'].values[1],color=cd['ENS'],marker='s',s=10**2)
    plt.scatter(df['ENSE'].values[0],df['ENSE'].values[1],color=cd['ENSE'],marker='o',s=10**2)
    plt.legend(('TOT','ANT','HAD','NSIDC','$\overline{TOT}$','$\overline{ANT}$'),loc='best',fancybox=True,framealpha=0.5)
    plt.xlabel(str(sea1) + ' ' + field1 + ' ' + region1)
    plt.ylabel(str(sea2) + ' ' + field2 + ' ' + region2)
    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5)
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5)

    if printtofile:
        fig.savefig('scatter_' + field1 + region1 + str(sea1) + '_v_' + field2 + region2 + str(sea2) + '.pdf')


    # ################# PAPER ###############
    
    # SCATTER PLOT OF SIM MEANS
    # prepare data for regression line
    #mm, bb, rval, pval, std_err = sp.stats.linregress(dfts[sim].values[0],dfts[sim].values[1])
    #  a way to subset: df.ix[:,'R1':'R5'] -- but it can't get both R and E data.
    #rdf=df.ix[:,'R1':'R5']
    #edf=df.ix[:,'E1':'E5']
    #regressdf = rdf.append(edf) # didn't append like I expected. inserted NaNs in places
    
    # this will work:
    allens=['R1','R2','R3','R4','R5','E1','E2','E3','E4','E5'] # has to be an array!
    sub = df.loc[:,allens]
    subx=sub.ix[0,:]
    suby=sub.ix[1,:]
    mm, bb, rval, pval, std_err = sp.stats.linregress(subx,suby)

    #sube = df.loc[0,['E1','E2','E3','E4','E5']]
    firebrick=ccm.get_linecolor('firebrick')

    printtofile=True
    fig,ax = plt.subplots(1)
    fig.set_size_inches(6,5)
    rr = plt.scatter(df.filter(regex='R').values[0],df.filter(regex='R').values[1],
                     color='0.3',marker='o',s=8**2,alpha=0.7)
    ee = plt.scatter(df.loc[region1,['E1','E2','E3','E4','E5']],df.loc[region2,['E1','E2','E3','E4','E5']],
                     color=firebrick,marker='o',s=8**2,alpha=0.7)
    #plt.scatter(df['HAD'].values[0],df['HAD'].values[1],color=cd['HAD'],marker='s',s=8**2)
    #plt.scatter(df['NSIDC'].values[0],df['NSIDC'].values[1],color=cd['NSIDC'],marker='s',s=8**2)
    #plt.scatter(df['ENS'].values[0],df['ENS'].values[1],color=cd['ENS'],marker='s',s=10**2)
    #plt.scatter(df['ENSE'].values[0],df['ENSE'].values[1],color=cd['ENSE'],marker='o',s=10**2)

    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    onex=np.linspace(axxlims[0],axxlims[1])
    plt.plot(onex,mm*onex + bb, color='0.5',linewidth=1)#,linestyle='--')
    rsq = rval**2
    print 'PAPER --- mean sims regr info: RVAL: ' + str(rval) + ', PVAL: ' + str(pval) + ', R2: ' + str(rsq)
    
    plt.legend((rr,ee),('Individual SIC forcings','Average SIC forcing'),
               loc='best',fancybox=True,framealpha=0.5)#,frameon=False)

    if sea1 == 'ND' and sea2=='ND' and region1=='bksmori' and region2 in ('eurasia','eurasiathin','eurasiamori'):
        # PAPER FIG!
        xlab = '$\Delta$ Nov-Dec Barents-Kara Seas'
        if field1=='pmsl':
            xlab = xlab + ' SLP (hPa)'
        else:
            xlab = xlab + ' Z500 (m)'

        ylab = '$\Delta$ Nov-Dec Eurasian SAT ($^\circ$C)'
        
    else:
        xlab = str(sea1) + ' ' + field1 + ' ' + region1
        ylab = str(sea2) + ' ' + field2 + ' ' + region2
        
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5,linestyle='--')
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5,linestyle='--')

    if printtofile:
        fig.savefig('scatterregress_' + field1 + region1 + str(sea1) + '_v_' +
                    field2 + region2 + str(sea2) + 'wacepap.pdf')

    #printtofile=False

    # #######################################

    # SCATTER PLOT OF SIM MEANS v2 (with conf intervals)
    fig,ax = plt.subplots(1) #plt.figure()
    for sim in ('R1','R2','R3','R4','R5'): #,'ENS','HAD','NSIDC'):
        rs = ax.scatter(df[sim].values[0],df[sim].values[1],color='k',marker='s',s=8**2,alpha=0.5)
        ax.plot(cidf[sim].values[0],
                (df[sim].values[1],df[sim].values[1]),
                color='k',linewidth=2,marker='_',markersize=6)
        ax.plot((df[sim].values[0],df[sim].values[0]),
            cidf[sim].values[1],
            color='k',linewidth=2,marker='_',markersize=6)
    for sim in ('E1','E2','E3','E4','E5'):#,'ENSE','HAD','NSIDC'):
        es = ax.scatter(df[sim].values[0],df[sim].values[1],color='0.5',marker='o',s=8**2,alpha=0.5)
        ax.plot(cidf[sim].values[0],
                (df[sim].values[1],df[sim].values[1]),
                color='0.5',linewidth=2,marker='_',markersize=6)
        ax.plot((df[sim].values[0],df[sim].values[0]),
            cidf[sim].values[1],
            color='0.5',linewidth=2,marker='_',markersize=6)
    for sim in ('HAD',):
        had = ax.scatter(df[sim].values[0],df[sim].values[1],color=cd[sim],marker='s',s=8**2,alpha=0.5)
        ax.plot(cidf[sim].values[0],
                (df[sim].values[1],df[sim].values[1]),
                color=cd[sim],linewidth=2,marker='_',markersize=6)
        ax.plot((df[sim].values[0],df[sim].values[0]),
            cidf[sim].values[1],
            color=cd[sim],linewidth=2,marker='_',markersize=6)
    for sim in ('NSIDC',):
        ns = ax.scatter(df[sim].values[0],df[sim].values[1],color=cd[sim],marker='s',s=8**2,alpha=0.5)
        ax.plot(cidf[sim].values[0],
                (df[sim].values[1],df[sim].values[1]),
                color=cd[sim],linewidth=2,marker='_',markersize=6)
        ax.plot((df[sim].values[0],df[sim].values[0]),
            cidf[sim].values[1],
            color=cd[sim],linewidth=2,marker='_',markersize=6)
    for sim in ('ENSE',): # different marker
        em = ax.scatter(df[sim].values[0],df[sim].values[1],color=cd[sim],marker='o',s=8**2,alpha=0.5)
        ax.plot(cidf[sim].values[0],
                (df[sim].values[1],df[sim].values[1]),
                color=cd[sim],linewidth=2,marker='_',markersize=6)
        ax.plot((df[sim].values[0],df[sim].values[0]),
            cidf[sim].values[1],
            color=cd[sim],linewidth=2,marker='_',markersize=6)
    for sim in ('ENS',):
        rm = ax.scatter(df[sim].values[0],df[sim].values[1],color=cd[sim],marker='s',s=8**2,alpha=0.5)
        ax.plot(cidf[sim].values[0],
                (df[sim].values[1],df[sim].values[1]),
                color=cd[sim],linewidth=2,marker='_',markersize=6)
        ax.plot((df[sim].values[0],df[sim].values[0]),
            cidf[sim].values[1],
            color=cd[sim],linewidth=2,marker='_',markersize=6)
    plt.legend((rs,es,had,ns,rm,em),
               ('TOT','ANT','HAD','NSIDC','$\overline{TOT}$','$\overline{ANT}$'),
               loc='best',fancybox=True,framealpha=0.5)
    plt.xlabel(str(sea1) + ' ' + field1 + ' ' + region1)
    plt.ylabel(str(sea2) + ' ' + field2 + ' ' + region2)
    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5)
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5)
    if printtofile:
        fig.savefig('scatterCI_' + field1 + region1 + str(sea1) + '_v_' + field2 + region2 + str(sea2) + '.pdf')



    ## fig,ax = plt.subplots(1) #plt.figure()
    ## plt.scatter(dfts['R4'].values[0],dfts['R4'].values[1],color=cd['R4'],marker='s',s=3**2,alpha=0.5)
    ## plt.scatter(dfts['R3'].values[0],dfts['R3'].values[1],color=cd['R3'],marker='s',s=3**2,alpha=0.5)
    ## plt.scatter(dfts['ENS'].values[0],dfts['ENS'].values[1],color=cd['ENS'],marker='s',s=3**2,alpha=0.5)
    ## plt.scatter(dfts['ENSE'].values[0],dfts['ENSE'].values[1],color=cd['ENSE'],marker='s',s=3**2,alpha=0.5)

    ## plt.scatter(df['R4'].values[0],df['R4'].values[1],color=cd['R4'],marker='s',s=8**2,edgecolor='k',alpha=0.3)
    ## plt.scatter(df['R3'].values[0],df['R3'].values[1],color=cd['R3'],marker='s',s=8**2,edgecolor='k',alpha=0.3)
    ## plt.scatter(df['ENS'].values[0],df['ENS'].values[1],color=cd['ENS'],marker='s',s=8**2,edgecolor='k',alpha=0.3)
    ## plt.scatter(df['ENSE'].values[0],df['ENSE'].values[1],color=cd['ENSE'],marker='s',s=8**2,edgecolor='k',alpha=0.3)

    ## plt.xlabel(sea1 + ' ' + field1 + ' ' + region1)
    ## plt.ylabel(sea2 + ' ' + field2 + ' ' + region2)
    ## axylims = ax.get_ylim()
    ## if axylims[0]<=0 and axylims[1]>=0:
    ##     ax.axhline(y=0,color='k',linewidth=.5)
    ## axxlims = ax.get_xlim()
    ## if axxlims[0]<=0 and axxlims[1]>=0:
    ##     ax.axvline(x=0,color='k',linewidth=.5)



    # SCATTER PLOT IN TIME AND SIM
    fig,axs = plt.subplots(1,2,sharey=True,sharex=True) #plt.figure()
    fig.set_size_inches(9,5)
    ax=axs[0]
    for sim in ('R1','R2','R3','R4','R5','ENS'):#,'HAD','NSIDC'):
        ax.scatter(dfts[sim].values[0],dfts[sim].values[1],color=cd[sim],marker='s',s=3**2,alpha=0.5)
    for sim in ('R1','R2','R3','R4','R5','ENS'): # means
        ax.scatter(df[sim].values[0],df[sim].values[1],color=cd[sim],marker='s',s=8**2,edgecolor='k',alpha=0.7)

    ax.set_xlabel(str(sea1) + ' ' + field1 + ' ' + region1)
    ax.set_ylabel(str(sea2) + ' ' + field2 + ' ' + region2)
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5)
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5)
    ax.set_title('TOT')

    ax=axs[1]
    for sim in ('E1','E2','E3','E4','E5','ENSE'):#,'HAD','NSIDC'):
        ax.scatter(dfts[sim].values[0],dfts[sim].values[1],color=cd[sim],marker='s',s=3**2,alpha=0.5)
    for sim in ('E1','E2','E3','E4','E5','ENSE'):
        ax.scatter(df[sim].values[0],df[sim].values[1],color=cd[sim],marker='s',s=8**2,edgecolor='k',alpha=0.7)
    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    onex=np.linspace(axxlims[0],axxlims[1])
    ax=axs[0]
    for sim in ('R1','R2','R3','R4','R5','ENS'):#,'HAD','NSIDC'):
        mm, bb, rval, pval, std_err = sp.stats.linregress(dfts[sim].values[0],dfts[sim].values[1])
        #mm, bb = np.polyfit(dfts[sim].values[0],dfts[sim].values[1], 1)
        ax.plot(onex,mm*onex + bb, color=cd[sim],linewidth=2)#,linestyle='--')
        #rr = np.corrcoef(dfts[sim].values[0],dfts[sim].values[1] )[0,1]
        rsq = rval**2
        print 'TOT: R: ' + str(rval) + ', R squared: ' + str(rsq)
        ## val = '$%.2f$'%(rsq)
        ## ax.annotate('$R^2$= ' + val, xy=(axxlims[0]+.1*axxlims[1], axylims[1]-.2*axylims[1]),
        ##             xycoords='data') # upper left?
    ax=axs[1]
    for sim in ('E1','E2','E3','E4','E5','ENSE'):#,'HAD','NSIDC'):
        mm, bb, rval, pval, std_err = sp.stats.linregress(dfts[sim].values[0],dfts[sim].values[1])
        #mm, bb = np.polyfit(dfts[sim].values[0],dfts[sim].values[1], 1)
        ax.plot(onex,mm*onex + bb, color=cd[sim],linewidth=2)#,linestyle='--')
        #rr = np.corrcoef(dfts[sim].values[0],dfts[sim].values[1] )[0,1]
        rsq = rval**2
        print 'ANT: R: ' + str(rval) + ', R squared: ' + str(rsq)

    ax.set_xlabel(str(sea1) + ' ' + field1 + ' ' + region1)
    ax.set_ylabel(str(sea2) + ' ' + field2 + ' ' + region2)
    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5)
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5)
    ax.set_title('ANT')

    if printtofile:
        fig.savefig('scatterregress_' + field1 + region1 + str(sea1) +\
                    '_v_' + field2 + region2 + str(sea2) + '_ANTTOTobssbplt.pdf')

    #flddregdf = pd.DataFrame(flddregdt,i)

    ## plt.figure()
    ## cplt.kemscatter(flddreg,flddreg2)
    ## plt.xlabel(region1)
    ## plt.ylabel(region2)
    ## plt.title(sea)



# ###### HISTOGRAMS #########

print 'HISTOGRAMS' # @@@
print allantr2.shape

## fig,ax=plt.subplots(1,1)
## nna,binsa,patchesa = ax.hist(allantr2,color='.5',alpha=0.5,normed=True,histtype='stepfilled')#,orientation='horizontal')#,cumulative=True)
## nnt, binst, patchest = ax.hist(alltotr2,color='b',alpha=0.5,normed=True,histtype='stepfilled')#,orientation='horizontal')#,cumulative=True)

from scipy.stats import norm

## anormed=norm.fit(allantr2)
## amean=anormed[0]
## asd=anormed[1]
## #Generate X points
## axlims = [-6*asd+amean, 6*asd+amean] # large limits
## axx = np.linspace(axlims[0],axlims[1],500)
## #Get Y points via Normal PDF with fitted parameters
## apdf_fitted = norm.pdf(axx,loc=amean,scale=asd)

## tnormed=norm.fit(alltotr2)
## tmean=tnormed[0]
## tsd=tnormed[1]
## #Generate X points
## txlims = [-6*tsd+tmean, 6*tsd+tmean] # large limits
## txx = np.linspace(txlims[0],txlims[1],500)
## #Get Y points via Normal PDF with fitted parameters
## tpdf_fitted = norm.pdf(txx,loc=tmean,scale=tsd)


## ant=ax.plot(axx,apdf_fitted,'k')
## tot=ax.plot(txx,tpdf_fitted,'b')
## ax.set_xlabel((field2 + region2 + '_' + str(sea2)))
## ax.set_title('Anomalies')
## ax.legend(('ANT','TOT'))
## if printtofile:
##     fig.savefig('histpdf_' + (field2 + region2 + '_' + str(sea2)) + '_ANTTOTanom.pdf')
## print 'ANT: ' + str(anormed)
## print 'TOT: ' + str(tnormed)


def plot_anttotanom_histpdf(antdata,totdata,printtofile=False,label=''):
    
    fig,ax=plt.subplots(1,1)
    nna,binsa,patchesa = ax.hist(antdata,color='.5',alpha=0.5,normed=True,histtype='stepfilled')#,orientation='horizontal')#,cumulative=True)
    nnt, binst, patchest = ax.hist(totdata,color='b',alpha=0.5,normed=True,histtype='stepfilled')#,orientation='horizontal')#,cumulative=True)

    anormed=norm.fit(antdata)
    amean=anormed[0]
    asd=anormed[1]
    #Generate X points
    axlims = [-4*asd+amean, 4*asd+amean] # large limits
    axx = np.linspace(axlims[0],axlims[1],500)
    #Get Y points via Normal PDF with fitted parameters
    apdf_fitted = norm.pdf(axx,loc=amean,scale=asd)

    tnormed=norm.fit(totdata)
    tmean=tnormed[0]
    tsd=tnormed[1]
    #Generate X points
    txlims = [-4*tsd+tmean, 4*tsd+tmean] # large limits
    txx = np.linspace(txlims[0],txlims[1],500)
    #Get Y points via Normal PDF with fitted parameters
    tpdf_fitted = norm.pdf(txx,loc=tmean,scale=tsd)

    ax.plot(axx,apdf_fitted,'k')
    ax.plot(txx,tpdf_fitted,'b')
    ax.set_xlabel(label)
    ax.set_title('Anomalies')
    ax.legend(('ANT','TOT'))
    if printtofile:
        fig.savefig('histpdf_' + label + '_ANTTOTanom.pdf')

    print '\n'
    print label
    print 'ANT (mean anom, anom std): ' + str(anormed)
    print 'TOT (mean anom, anom std): ' + str(tnormed)


def plot_anttotsbplt_histpdf(antdata,totdata,printtofile=False,label=''):

    """ antdata and totdata must be tuples of CTL, PERT !!
    """
    
    allantc=antdata[0]
    allantp=antdata[1]
    alltotc=totdata[0]
    alltotp=totdata[1]
    
    antcfld = allantc-np.mean(allantc)
    antpfld = allantp-np.mean(allantc) # subtract its own mean or control mean??
    print '\n'
    print label
    print 'antcfld mean ' + str(antcfld.mean())
    print 'antpfld mean ' + str(antpfld.mean())
    
    totcfld = alltotc-np.mean(alltotc)
    totpfld = alltotp-np.mean(alltotc)
    
    fig,axs=plt.subplots(1,2) ########### SUBPLOT
    fig.set_size_inches(10,5)
    ax=axs[0]
    nnac,binsac,patchesac = ax.hist(antcfld,color='.5',alpha=0.5,normed=True,histtype='stepfilled')#,orientation='horizontal')#,cumulative=True)
    nnap, binsap, patchesap = ax.hist(antpfld,color='r',alpha=0.5,normed=True,histtype='stepfilled')#,orientation='horizontal')#,cumulative=True)

    acnormed=norm.fit(antcfld)
    amean=acnormed[0]
    asd=acnormed[1]
    #Generate X points
    axlims = [-4*asd+amean, 4*asd+amean] # large limits
    axx = np.linspace(axlims[0],axlims[1],500)
    #Get Y points via Normal PDF with fitted parameters
    acpdf_fitted = norm.pdf(axx,loc=amean,scale=asd)

    apnormed=norm.fit(antpfld)
    apmean=apnormed[0]
    apsd=apnormed[1]
    #Generate X points
    apxlims = [-4*apsd+apmean, 4*apsd+apmean] # large limits
    apxx = np.linspace(apxlims[0],apxlims[1],500)
    #Get Y points via Normal PDF with fitted parameters
    appdf_fitted = norm.pdf(apxx,loc=apmean,scale=apsd)

    ax.plot(axx,acpdf_fitted,'k')
    ax.plot(apxx,appdf_fitted,'r')
    ax.set_xlabel(label)
    ax.set_title('ANT')
    ax.legend(('CTL','PERT'))

    print label
    print 'ANT CTL (mean,std): ' + str(acnormed)
    print 'ANT PERT (mean,std): ' + str(apnormed)

    ax=axs[1]
    nntc, binstc, patchestc = ax.hist(totcfld,color='0.5',alpha=0.5,normed=True,histtype='stepfilled')
    nntp, binstp, patchestp = ax.hist(totpfld,color='r',alpha=0.5,normed=True,histtype='stepfilled')

    tcnormed=norm.fit(totcfld)
    tmean=tcnormed[0]
    tsd=tcnormed[1]
    #Generate X points
    txlims = [-4*tsd+tmean, 4*tsd+tmean] # large limits
    txx = np.linspace(txlims[0],txlims[1],500)
    #Get Y points via Normal PDF with fitted parameters
    tcpdf_fitted = norm.pdf(txx,loc=tmean,scale=tsd)

    tpnormed=norm.fit(totpfld)
    tpmean=tpnormed[0]
    tpsd=tpnormed[1]
    #Generate X points
    tpxlims = [-4*tpsd+tpmean, 4*tpsd+tpmean] # large limits
    tpxx = np.linspace(tpxlims[0],tpxlims[1],500)
    #Get Y points via Normal PDF with fitted parameters
    tppdf_fitted = norm.pdf(tpxx,loc=tpmean,scale=tpsd)

    ax.plot(txx,tcpdf_fitted,'k')
    ax.plot(tpxx,tppdf_fitted,'r')
    ax.set_xlabel(label)
    ax.set_title('TOT')
    ax.legend(('CTL','PERT'))

    print label
    print 'TOT CTL (mean,std): ' + str(tcnormed)
    print 'TOT PERT (mean,std): ' + str(tpnormed)
    if printtofile:
        fig.savefig('histpdf_' + label + '_ANTTOTsbplt.pdf')


printtofile=False
plot_anttotanom_histpdf(allantr2,alltotr2,
                        label=(field2 + region2 + '_' + str(sea2)),printtofile=printtofile)
cutl.calc_pvals(allantr2,alltotr2)
print '-------------------------'

plot_anttotanom_histpdf(allantr1,alltotr1,
                        label=(field1 + region1 + '_' + str(sea1)),printtofile=printtofile)
cutl.calc_pvals(allantr1,alltotr1)


plot_anttotsbplt_histpdf((allantcr2,allantpr2),(alltotcr2,alltotpr2),
                         label=(field2 + region2 + '_' + str(sea2)),printtofile=printtofile)
plot_anttotsbplt_histpdf((allantcr1,allantpr1),(alltotcr1,alltotpr1),
                         label=(field1 + region1 + '_' + str(sea1)),printtofile=printtofile)

