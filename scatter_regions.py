
import cccmautils as cutl
import constants as con
import cccmaplots as cplt
import cccmacmaps as ccm
import pandas as pd

printtofile=True
plt.close('all')

conv1=1; conv2=1

field1='st'; ncfield1='ST'
#field1='sia'; ncfield1='SICN'
#field1='gz50000'; ncfield1='PHI'; conv1=1/con.get_g()
#field1='pmsl'; ncfield1='PMSL'
region1='polcap60' #'bksmori' #'polcap65'
sea1='SON'

#field2='st'; ncfield2='ST'
#field2='sia'; ncfield2='SICN'
field2='pmsl'; ncfield2='PMSL'
#field2='gz50000'; ncfield2='PHI'; conv2=1/con.get_g()
region2= 'polcap60' #'eurasiamori'
sea2='SON'

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

alltotcr1=np.zeros(120*5)
alltotcr2=np.zeros(120*5)
alltotpr1=np.zeros(120*5)
alltotpr2=np.zeros(120*5)
alltotr1=np.zeros(120*5)
alltotr2=np.zeros(120*5)

allantcr1=np.zeros(120*5)
allantcr2=np.zeros(120*5)
allantpr1=np.zeros(120*5)
allantpr2=np.zeros(120*5)
allantr1=np.zeros(120*5)
allantr2=np.zeros(120*5)
tallii=0 # index to keep track of time in accumulated TOT ensemble
aallii=0 # index to keep track of time in accumulated ANT ensemble
for sim in sims:

    fnamec,fnamep=con.build_filepathpair(sim,field1)

    fldc=cnc.getNCvar(fnamec,ncfield1,timesel='0002-01-01,0121-12-31',seas=sea1)*conv1
    fldp=cnc.getNCvar(fnamep,ncfield1,timesel='0002-01-01,0121-12-31',seas=sea1)*conv1


    lat=cnc.getNCvar(fnamec,'lat')
    lon=cnc.getNCvar(fnamec,'lon')

    if field1=='sia':
        fldcreg=cutl.calc_regtotseaicearea(fldc,lat,lon,region1) # isarea=False
        fldpreg=cutl.calc_regtotseaicearea(fldp,lat,lon,region1) # isarea=False
    else:
        fldcreg=cutl.calc_regmean(fldc,lat,lon,region1)
        fldpreg=cutl.calc_regmean(fldp,lat,lon,region1)
        
    flddregdt[sim] = np.mean(fldpreg-fldcreg,axis=0) # time mean
    if sea2 is 'DJF' and sea1 is not 'DJF': # have to shorten other timeseries
        flddregtsdt[sim] = fldpreg[:-1,...]-np.mean(fldcreg[:-1,...],axis=0) # anomaly timeseries from ctl mean
    else:
        flddregtsdt[sim] = fldpreg-np.mean(fldcreg,axis=0)

    if field2 != field1:
        fnamec,fnamep=con.build_filepathpair(sim,field2)
        fldc2=cnc.getNCvar(fnamec,ncfield2,timesel='0002-01-01,0121-12-31',seas=sea2)*conv2
        fldp2=cnc.getNCvar(fnamep,ncfield2,timesel='0002-01-01,0121-12-31',seas=sea2)*conv2
    elif sea2 != sea1:
        fldc2=cnc.getNCvar(fnamec,ncfield2,timesel='0002-01-01,0121-12-31',seas=sea2)*conv2
        fldp2=cnc.getNCvar(fnamep,ncfield2,timesel='0002-01-01,0121-12-31',seas=sea2)*conv2
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
    if sea1 is 'DJF' and sea2 is not 'DJF':
        flddregtsdt2[sim] = fldpreg2[:-1,...]-np.mean(fldcreg2[:-1,...],axis=0)
    else:
        flddregtsdt2[sim] = fldpreg2-np.mean(fldcreg2,axis=0) # anomaly timeseries from ctl mean

    if sim in TOT:
        #print 'concat tot ens members'
        alltotcr1[tallii:tallii+120] = fldcreg
        alltotcr2[tallii:tallii+120] = fldcreg2
        alltotpr1[tallii:tallii+120] = fldpreg
        alltotpr2[tallii:tallii+120] = fldpreg2
        alltotr1[tallii:tallii+120] = flddregtsdt[sim]
        alltotr2[tallii:tallii+120] = flddregtsdt2[sim]
        tallii+=120
        #totser.append(pd.Series(flddregtsdt))
        #totser2.append(pd.Series(flddregtsdt2))
    elif sim in ANT:
        #print 'concat ant ens members'
        allantcr1[aallii:aallii+120] = fldcreg
        allantcr2[aallii:aallii+120] = fldcreg2
        allantpr1[aallii:aallii+120] = fldpreg
        allantpr2[aallii:aallii+120] = fldpreg2
        allantr1[aallii:aallii+120] = flddregtsdt[sim]
        allantr2[aallii:aallii+120] = flddregtsdt2[sim]
        aallii+=120
        
    # calculate confidence interval
    # double-check the scale setting
    ciregdt[sim] = sp.stats.t.interval(1-siglevel,len(fldpreg)-1,loc=np.mean(fldpreg,axis=0)-np.mean(fldcreg,axis=0),
                             scale=np.std(fldpreg,axis=0)/np.sqrt(len(fldpreg)))
    ciregdt2[sim] = sp.stats.t.interval(1-siglevel,len(fldpreg2)-1,loc=np.mean(fldpreg2,axis=0)-np.mean(fldcreg2,axis=0),
                             scale=np.std(fldpreg2,axis=0)/np.sqrt(len(fldpreg2)))


# @@@ add confidence int for allant and alltot (means have to concat control and pert separately as well)
# calculate confidence interval
# double-check the scale setting
ciregalltotr1 = sp.stats.t.interval(1-siglevel,len(alltotpr1)-1,loc=np.mean(alltotpr1,axis=0)-np.mean(alltotcr1,axis=0),
                         scale=np.std(alltotpr1,axis=0)/np.sqrt(len(alltotpr1)))
ciregalltotr2 = sp.stats.t.interval(1-siglevel,len(alltotpr2)-1,loc=np.mean(alltotpr2,axis=0)-np.mean(alltotcr2,axis=0),
                         scale=np.std(alltotpr2,axis=0)/np.sqrt(len(alltotpr2)))
ciregallantr1 = sp.stats.t.interval(1-siglevel,len(allantpr1)-1,loc=np.mean(allantpr1,axis=0)-np.mean(allantcr1,axis=0),
                         scale=np.std(allantpr1,axis=0)/np.sqrt(len(allantpr1)))
ciregallantr2 = sp.stats.t.interval(1-siglevel,len(allantpr2)-1,loc=np.mean(allantpr2,axis=0)-np.mean(allantcr2,axis=0),
                         scale=np.std(allantpr2,axis=0)/np.sqrt(len(allantpr2)))


# super ensembles
plt.figure();
plt.scatter(allantr1,allantr2,color='r')
plt.scatter(alltotr1,alltotr2,color='k')
antr1std=allantr1.std()
antr2std=allantr2.std()
totr1std=alltotr1.std()
totr2std=alltotr2.std()

print 'really should be printing the conf int range on the figure@@'
plt.figure();
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


# SCATTER PLOT OF SIM MEANS v2
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

