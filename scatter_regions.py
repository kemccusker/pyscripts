
import cccmautils as cutl
import constants as con
import cccmaplots as cplt
import cccmacmaps as ccm
import pandas as pd

printtofile=True

field1='sia'; ncfield1='SICN'
region1='bksmori' #'polcap65'
sea1='DJF'

field2='st'; ncfield2='ST'
region2='eurasiamori'
sea2='DJF'

sims = ('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5','HAD','NSIDC','ENS','ENSE')

flddregdt= {}
flddregdt2 = {}

for sim in sims:

    fnamec,fnamep=con.build_filepathpair(sim,field1)

    fldc=cnc.getNCvar(fnamec,ncfield1,timesel='0002-01-01,0121-12-31',seas=sea)
    fldp=cnc.getNCvar(fnamep,ncfield1,timesel='0002-01-01,0121-12-31',seas=sea)


    lat=cnc.getNCvar(fnamec,'lat')
    lon=cnc.getNCvar(fnamec,'lon')

    if field1=='sia':
        fldcreg=cutl.calc_regtotseaicearea(fldc,lat,lon,region1) # isarea=False
        fldpreg=cutl.calc_regtotseaicearea(fldp,lat,lon,region1) # isarea=False
    else:
        fldcreg=cutl.calc_regmean(fldc,lat,lon,region1)
        fldpreg=cutl.calc_regmean(fldp,lat,lon,region1)
        
    flddregdt[sim] = np.mean(fldpreg-fldcreg,axis=0) # time mean

    if field2 != field1:
        fnamec,fnamep=con.build_filepathpair(sim,field2)
        fldc2=cnc.getNCvar(fnamec,ncfield2,timesel='0002-01-01,0121-12-31',seas=sea2)
        fldp2=cnc.getNCvar(fnamep,ncfield2,timesel='0002-01-01,0121-12-31',seas=sea2)
    elif sea2 != sea1:
        fldc2=cnc.getNCvar(fnamec,ncfield2,timesel='0002-01-01,0121-12-31',seas=sea2)
        fldp2=cnc.getNCvar(fnamep,ncfield2,timesel='0002-01-01,0121-12-31',seas=sea2)
    else:
        fldc2=fldc
        fldp2=fldp

    if field2=='sia':
        fldcreg2=cutl.calc_regtotseaicearea(fldc,lat,lon,region1) # isarea=False
        fldpreg2=cutl.calc_regtotseaicearea(fldp,lat,lon,region1) # isarea=False
    else:
        fldcreg2=cutl.calc_regmean(fldc2,lat,lon,region2)
        fldpreg2=cutl.calc_regmean(fldp2,lat,lon,region2)
        
    flddregdt2[sim] = np.mean(fldpreg2-fldcreg2,axis=0)

ser1 = pd.Series(flddregdt)
ser2 = pd.Series(flddregdt2)

df = pd.DataFrame([ser1,ser2],index=(region1,region2))

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
    
#flddregdf = pd.DataFrame(flddregdt,i)

## plt.figure()
## cplt.kemscatter(flddreg,flddreg2)
## plt.xlabel(region1)
## plt.ylabel(region2)
## plt.title(sea)

