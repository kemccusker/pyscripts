
import cccmautils as cutl
import constants as con
import cccmaplots as cplt
import pandas as pd

field='st'
ncfield=field.upper()
sea='DJF'
region1='polcap65'
region2='eurasia'

sims = ('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5')

flddregdt= {}
flddregdt2 = {}

for sim in sims:

    fnamec,fnamep=con.build_filepathpair(sim,field)

    fldc=cnc.getNCvar(fnamec,ncfield,timesel='0002-01-01,0121-12-31',seas=sea)
    fldp=cnc.getNCvar(fnamep,ncfield,timesel='0002-01-01,0121-12-31',seas=sea)


    lat=cnc.getNCvar(fnamec,'lat')
    lon=cnc.getNCvar(fnamec,'lon')

    fldcreg=cutl.calc_regmean(fldc,lat,lon,region1)
    fldpreg=cutl.calc_regmean(fldp,lat,lon,region1)
    flddregdt[sim] = np.mean(fldpreg-fldcreg,axis=0)
    
    fldcreg2=cutl.calc_regmean(fldc,lat,lon,region2)
    fldpreg2=cutl.calc_regmean(fldp,lat,lon,region2)
    flddregdt2[sim] = np.mean(fldpreg2-fldcreg2,axis=0)

ser1 = pd.Series(flddregdt)
ser2 = pd.Series(flddregdt2)

df = pd.DataFrame([ser1,ser2],index=(region1,region2))

print df

plt.figure()
plt.scatter(df.values[0],df.values[1])
plt.scatter(df.filter(regex='E').values[0],df.filter(regex='E').values[1],color='r')

#flddregdf = pd.DataFrame(flddregdt,i)

## plt.figure()
## cplt.kemscatter(flddreg,flddreg2)
## plt.xlabel(region1)
## plt.ylabel(region2)
## plt.title(sea)

