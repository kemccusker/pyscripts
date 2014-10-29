
import cccmautils as cutl
import constants as con
import cccmaplots as cplt
import cccmacmaps as ccm
import pandas as pd

printtofile=True
plt.close('all')

conv1=1; conv2=1

#field1='sia'; ncfield1='SICN'
#field1='gz50000'; ncfield1='PHI'; conv1=1/con.get_g()
field1='pmsl'; ncfield1='PMSL'
region1='polcap70' #'polcap65'
sea1='DJF'

field2='st'; ncfield2='ST'
#field2='sia'; ncfield2='SICN'
region2= 'ntham'
sea2='DJF'

sims = ('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5','HAD','NSIDC','ENS','ENSE')

flddregdt= {}
flddregdt2 = {}
flddregtsdt= {}
flddregtsdt2 = {}

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
    if sea2 is 'DJF' and sea1 is not 'DJF':
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

ser1 = pd.Series(flddregdt)
ser2 = pd.Series(flddregdt2)
serts1=pd.Series(flddregtsdt) 
serts2=pd.Series(flddregtsdt2)

df = pd.DataFrame([ser1,ser2],index=(region1,region2))
dfts = pd.DataFrame([serts1,serts2],index=(region1,region2))

print df

cd = ccm.get_colordict()

fig,ax = plt.subplots(1) #plt.figure()
plt.scatter(df.filter(regex='R').values[0],df.filter(regex='R').values[1],color='k',marker='s',s=8**2)
plt.scatter(df.filter(regex='E').values[0],df.filter(regex='E').values[1],color='0.5',marker='*',s=8**2)
plt.scatter(df.filter(regex='HAD').values[0],df.filter(regex='HAD').values[1],color=cd['HAD'],s=8**2)
plt.scatter(df.filter(regex='NSIDC').values[0],df.filter(regex='NSIDC').values[1],color=cd['NSIDC'],s=8**2)
plt.scatter(df.filter(regex='ENS').values[0],df.filter(regex='ENS').values[1],color=cd['ENS'],s=10**2)
plt.scatter(df.filter(regex='ENSE').values[0],df.filter(regex='ENSE').values[1],color=cd['ENSE'],s=10**2)

plt.legend(('TOT','ANT','HAD','NSIDC','$\overline{TOT}$','$\overline{ANT}$'),loc='best',fancybox=True,framealpha=0.5)
plt.xlabel(sea1 + ' ' + field1 + ' ' + region1)
plt.ylabel(sea2 + ' ' + field2 + ' ' + region2)
axylims = ax.get_ylim()
if axylims[0]<=0 and axylims[1]>=0:
    ax.axhline(y=0,color='k',linewidth=.5)
axxlims = ax.get_xlim()
if axxlims[0]<=0 and axxlims[1]>=0:
    ax.axvline(x=0,color='k',linewidth=.5)

if printtofile:
    fig.savefig('scatter_' + field1 + region1 + sea1 + '_v_' + field2 + region2 + sea2 + '.pdf')


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


fig,axs = plt.subplots(1,2,sharey=True,sharex=True) #plt.figure()
fig.set_size_inches(9,5)
ax=axs[0]
for sim in ('R1','R2','R3','R4','R5','ENS'):
    ax.scatter(dfts[sim].values[0],dfts[sim].values[1],color=cd[sim],marker='s',s=3**2,alpha=0.5)
for sim in ('R1','R2','R3','R4','R5','ENS'):
    ax.scatter(df[sim].values[0],df[sim].values[1],color=cd[sim],marker='s',s=8**2,edgecolor='k',alpha=0.7)

ax.set_xlabel(sea1 + ' ' + field1 + ' ' + region1)
ax.set_ylabel(sea2 + ' ' + field2 + ' ' + region2)
axylims = ax.get_ylim()
if axylims[0]<=0 and axylims[1]>=0:
    ax.axhline(y=0,color='k',linewidth=.5)
axxlims = ax.get_xlim()
if axxlims[0]<=0 and axxlims[1]>=0:
    ax.axvline(x=0,color='k',linewidth=.5)
ax.set_title('TOT')

ax=axs[1]

#fig,ax = plt.subplots(1) #plt.figure()
for sim in ('E1','E2','E3','E4','E5','ENSE'):
    ax.scatter(dfts[sim].values[0],dfts[sim].values[1],color=cd[sim],marker='s',s=3**2,alpha=0.5)
for sim in ('E1','E2','E3','E4','E5','ENSE'):
    ax.scatter(df[sim].values[0],df[sim].values[1],color=cd[sim],marker='s',s=8**2,edgecolor='k',alpha=0.7)

ax.set_xlabel(sea1 + ' ' + field1 + ' ' + region1)
ax.set_ylabel(sea2 + ' ' + field2 + ' ' + region2)
axylims = ax.get_ylim()
if axylims[0]<=0 and axylims[1]>=0:
    ax.axhline(y=0,color='k',linewidth=.5)
axxlims = ax.get_xlim()
if axxlims[0]<=0 and axxlims[1]>=0:
    ax.axvline(x=0,color='k',linewidth=.5)
ax.set_title('ANT')

if printtofile:
    fig.savefig('scatter_' + field1 + region1 + sea1 + '_v_' + field2 + region2 + sea2 + '_ANTTOTsbplt.pdf')

#flddregdf = pd.DataFrame(flddregdt,i)

## plt.figure()
## cplt.kemscatter(flddreg,flddreg2)
## plt.xlabel(region1)
## plt.ylabel(region2)
## plt.title(sea)

