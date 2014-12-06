""" From Neil:
import pandas as pd
df = pd.DataFrame.from_csv('/raid/ra40/data/ncs/for_edhawkins_sic/rcp45/cmip5_rcp45_sie.csv')
df_CanESM2 = df.filter(regex='CanESM2')

# plot September extend
df[df.index.month==9].resample('A').plot(color='k',alpha=0.5, legend=False)

df_CanESM2[df_CanESM2.index.month==9].resample('A').plot(color='g', 
legend=False)

   cmip5seaice.py: 10/7/2014
       calc the range of sea ice extent trends for within model ensembles, etc.
       
"""

import pandas as pd
import constants as con
import collections as coll
import scipy as sp
import scipy.stats
import cccmacmaps as ccm
import cccmautils as cutl
import cmip5meta as cm5

plt.close('all')
plt.ion()

printtofile=False


df = pd.DataFrame.from_csv('/raid/ra40/data/ncs/for_edhawkins_sic/rcp45/cmip5_rcp45_sie.csv')
df_CanESM2 = df.filter(regex="CanESM2")

models = ('ACCESS1-0','ACCESS1-3','BNU-ESM','CCSM4','CESM1-BGC','CESM1-CAM5','CMCC-CMS',
          'CMCC-CM','CNRM-CM5','CSIRO-Mk3-6-0','CanESM2','EC-EARTH','FGOALS-g2','GFDL-CM3',
          'GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR',
          'IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC-ESM-CHEM','MIROC-ESM','MIROC5','MPI-ESM-LR',
          'MPI-ESM-MR','MRI-CGCM3','NorESM1-ME')
styr='1979'
endyr='2012' # inclusive

mons = con.get_mon()
montrenddt = coll.OrderedDict() #dict.fromkeys(mons)
monmndt = coll.OrderedDict()
monmxdt = coll.OrderedDict()
monavgdt = coll.OrderedDict()
hmontrenddt = coll.OrderedDict()
nmontrenddt = coll.OrderedDict()


#mondf = pd.DataFrame() # index will be np.arange
totens=len(df.keys()) # total # of ens members

## # trend calc
## mm, bb = np.polyfit(xx, dat, 1)
## ax.plot(onex,mm*onex + bb, color='k')#,linestyle='--')

## better trend
## slope, intercept, r_value, p_value, std_err = stats.linregress(xi,y)



# ================= OBS ===========
basepath2 = '/home/rkm/work/BCs/'
#fhadsicc = basepath2 + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_sicn_' +\
#         ctimstr + 'climo.nc' #SICN, 129x64 CLIMO
fhad = basepath2 + 'HadISST/hadisst1.1_bc_128_64_1870_2013m03_SIEnhfrac_1870010100-2013030100.nc'

#fnsidcsicc = basepath2 + 'NSIDC/nsidc_bt_128x64_1978m11_2011m12_sicn_' + ctimstr + 'climo.nc'
fnsidc = basepath2 + 'NSIDC/nsidc_bt_128x64_1978m11_2011m12_SIEnhfrac_1978111600-2011121612.nc' #SICN, 129x64

had = cnc.getNCvar(fhad,'SICN',timesel='1979-01-01,2012-12-31')
hadtime = cnc.getNCvar(fhad,'time',timesel='1979-01-01,2012-12-31')
nsidc = cnc.getNCvar(fnsidc,'SICN',timesel='1979-01-01,2012-12-31')

earthrad = con.get_earthrad()
totalarea = 4*np.pi*earthrad**2
print had.shape
print nsidc.shape
had=had*(totalarea/2.) # convert to SIE from hemispheric fraction
nsidc=nsidc*(totalarea/2.) # convert to SIE from hemispheric fraction

# =================================




#superii=0 # super index for all ensemble member trends
superslopes = np.zeros((len(mons),totens))
superkeys = ()#np.zeros((len(mons),totens))

for monii,mo in enumerate(mons):

     superii=0 # super index for all ensemble member trends
     thismonskey=()

     modensdt = dict.fromkeys(models)
     modtrenddt = dict.fromkeys(models)
     mndt = dict.fromkeys(models)
     mxdt = dict.fromkeys(models)
     avgdt = dict.fromkeys(models)

     enumdt={} # number of ensemble members
     for mod in models:

         onemod = df.filter(regex=mod)[styr:endyr] #DataFrame of DataFrames where e/ sub DF is one model's ensemble

         modata = onemod[onemod.index.month==monii+1].resample('A') # resample annually for given month
         
         modensdt[mod] = modata # shape of values is # months x # ens mems
         ensnum = modata.values.shape[-1]
         enumdt[mod] = ensnum
         xx = np.arange(0,modata.values.shape[0])

         slope= np.zeros((ensnum,))
         for eii in np.arange(0,ensnum):
             # calculate trend
             dat = np.squeeze(modata.values[...,eii])
             #print dat.shape
             slope[eii], intercept, r_value, p_value, std_err = sp.stats.linregress(xx,dat)

             superslopes[monii,superii] = slope[eii]
             thismonskey = thismonskey + (mod,)
             #superkeys[monii,superii] = mod
             superii=superii+1
         
         modtrenddt[mod] = slope # trends for e/ ens member in the model group
         mntrnd = np.min(slope)
         mxtrnd = np.max(slope)
         avgtrnd = np.mean(slope)
         mndt[mod] = mntrnd
         mxdt[mod] = mxtrnd
         avgdt[mod] = avgtrnd

     # do obs e/ month
     homon = cutl.seasonalize_monthlyts(had,mo=monii+1)
     xx = np.arange(0,homon.shape[0])
     hslope, intercept, r_value, p_value, std_err = sp.stats.linregress(xx,homon)
     hmontrenddt[mo] = hslope
     
     nomon = cutl.seasonalize_monthlyts(nsidc,mo=monii+1)
     xx = np.arange(0,nomon.shape[0])
     nslope, intercept, r_value, p_value, std_err = sp.stats.linregress(xx,nomon)
     nmontrenddt[mo] = nslope

     superkeys = superkeys + (thismonskey,)
     montrenddt[mo] = modtrenddt
     monmndt[mo]= mndt
     monmxdt[mo]= mxdt
     monavgdt[mo] = avgdt

montrenddf = pd.DataFrame(montrenddt)
mndf = pd.DataFrame(monmndt)
mxdf = pd.DataFrame(monmxdt)
avgdf = pd.DataFrame(monavgdt)

superensdf=pd.DataFrame(data=superslopes,index=mons,columns=superkeys[0])
canesm=superensdf.CanESM2

# plot seasonal trends for all ensemble members
hmontrenddf = pd.DataFrame(hmontrenddt,index=(1,))
nmontrenddf = pd.DataFrame(nmontrenddt,index=(1,))

superensdf.plot(color='k',alpha=0.5,legend=False)
plt.plot(canesm,color='r',linewidth=3)
plt.plot(hmontrenddf.values[0],'green',linewidth=3) # not sure why but this has an array w/in an array....
plt.plot(nmontrenddf.values[0],color=ccm.get_linecolor('dodgerblue'),linewidth=3)# not sure why but this has an array w/in an array....
plt.title('SIE 1979-2012 trends (/yr)')
plt.xlabel('Month')
plt.xlim((0,11))
if printtofile:
    plt.savefig('SIE_CMIP5_allmemstrend_seacycle_1979-2012_wCanESM2.pdf')

#superensdf.hist()


# min / max histograms each month
fig,axs = plt.subplots(3,4)
fig.set_size_inches(12,9)
for aii,ax in enumerate(axs.flat):
    mon = mons[aii]
    ax.hist(mndf[mon],color='.5',alpha=0.5)
    ax.hist(mxdf[mon],color='orange',alpha=0.5)
    ax.axvline(mndf[mon]['CanESM2'],color='k',linewidth=3)
    ax.axvline(mxdf[mon]['CanESM2'],color='r',linewidth=3)
    ax.set_title(mon)
if printtofile:
    fig.savefig('SIE_CMIP5_minmaxtrendhist_allmos_1979-2012_wCanESM2.pdf')

# average trend histograms each month
fig,axs = plt.subplots(3,4)
fig.set_size_inches(12,9)
for aii,ax in enumerate(axs.flat):
    mon = mons[aii]
    ax.hist(avgdf[mon],color='.5',alpha=0.5)
    ax.axvline(avgdf[mon]['CanESM2'],color='k',linewidth=3)
    ax.axvline(hmontrenddt[mon],color='green',linewidth=3)
    ax.axvline(nmontrenddt[mon],color=ccm.get_linecolor('dodgerblue'),linewidth=3)
    ax.set_title(mon)
if printtofile:
    fig.savefig('SIE_CMIP5_avgtrendhist_allmos_1979-2012_wCanESM2.pdf')


# want ALL ens trends in one PDF (per month).
fig,axs = plt.subplots(3,4,sharex=True)
fig.set_size_inches(12,9)
for rii,row in enumerate(superensdf.iterrows()):
    ax = axs.flat[rii]
    mon=row[0] # the name is first in tuple
    rseries = row[1] # pandas Series
    
    ax.hist(rseries,color='.5',alpha=0.5)
    for vline in rseries.CanESM2:
        ax.axvline(vline,color='k',linewidth=3)
    ax.axvline(hmontrenddt[mon],color='green',linewidth=3)
    ax.axvline(nmontrenddt[mon],color=ccm.get_linecolor('dodgerblue'),linewidth=3)
    
    ax.set_title(mon)
if printtofile:
    fig.savefig('SIE_CMIP5_allmemstrendhist_allmos_1979-2012_wCanESM2.pdf')


# CanESM is index 10
cmap='red2blue_w20'

mycmap = plt.cm.get_cmap(cmap)
    
#plt.rc('axes', color_cycle=mycc) # default colorcycle


canidx=10
canxx=np.squeeze(np.ones((len(mons),1))*canidx)
fig,ax = plt.subplots()
fig.set_size_inches(14,3)
#ax.set_color_cycle(mycc)

mons2= ('Nov','Dec','Jan','Feb',
        'Mar','Apr','May',
        'Jun','Jul','Aug','Sep',
        'Oct')
sep=10
mar=4
colors=('skyblue','steelblue3','steelblue4','mediumblue',
        'darkseagreen4','darkseagreen','darkolivegreen3',
        'warm5','warm3','warm2','warm1',
        'chocolate4')


for mii,mon in enumerate(mons2):
    plt.plot(mndf[mon],color=ccm.get_linecolor(colors[mii]),linestyle='None',marker='s')
    plt.plot(mxdf[mon],color=ccm.get_linecolor(colors[mii]),linestyle='None',marker='o')
plt.plot(mndf['Sep'],color='orange')
plt.plot(mxdf['Sep'],color=ccm.get_linecolor(colors[sep]))

#plt.plot(mndf,linestyle='None',marker='s')
#plt.plot(mxdf,'r',linestyle='None',marker='o')
plt.plot(canxx,mndf.loc['CanESM2'],color=ccm.get_linecolor(colors[mii]),linestyle='None',marker='s')
plt.plot(canxx,mxdf.loc['CanESM2'],color=ccm.get_linecolor(colors[mii]),linestyle='None',marker='o')
ax.set_xticks(np.arange(0,29))
ax.set_xlim((-1,30))
ax.set_xticklabels(mxdf['Sep'].keys())
ax.set_title('SIE min (square), max (circle) trends')
fig.autofmt_xdate()
if printtofile:
    fig.savefig('SIE_CMIP5_minmaxtrend_allmodels_1979-2012.pdf')


# Trend RANGE: By model
fig,ax = plt.subplots()
fig.set_size_inches(14,4)

for mii,mon in enumerate(mons2):
    plt.plot(mxdf[mon]-mndf[mon],color=ccm.get_linecolor(colors[mii]),linestyle='None',marker='s')
plt.plot(mxdf['Sep']-mndf['Sep'],color=ccm.get_linecolor(colors[sep]))
plt.plot(mxdf['Mar']-mndf['Mar'],color=ccm.get_linecolor(colors[mar]))

plt.plot(canxx,mxdf.loc['CanESM2']-mndf.loc['CanESM2'],color=ccm.get_linecolor(colors[mii]),linestyle='None',marker='s')
ax.set_xticks(np.arange(0,29))
ax.set_xlim((-1,30))
ax.set_xticklabels(mxdf['Sep'].keys())
ax.set_title('Within model SIE trend range')
fig.autofmt_xdate()
if printtofile:
    fig.savefig('SIE_CMIP5_mnmxtrendrange_allmodels_1979-2012.pdf')


rngdf=mxdf-mndf
rngdf=rngdf[rngdf>0] # get rid of ranges of zero (1 ens member)

# ALL MONTHS hists
rngdf.hist()
plt.title('Within model SIE trend ranges')


# Sep and March hists together
fig,ax = plt.subplots()
rngdf['Sep'].hist(color=ccm.get_linecolor(colors[sep]),alpha=0.5)
rngdf['Mar'].hist(color=ccm.get_linecolor(colors[mar]),alpha=0.5)
ax.axvline(rngdf['Sep']['CanESM2'],color=ccm.get_linecolor(colors[sep]),linewidth=3)
ax.axvline(rngdf['Mar']['CanESM2'],color=ccm.get_linecolor(colors[mar]),linewidth=3)
ax.set_title('Sep/Mar within-CMIP5-model trend range (1979-2012)')


# SEA CYCLE trend RANGE
fig,ax = plt.subplots()
fig.set_size_inches(7,4)

for mii,mod in enumerate(models):
    plt.plot(mxdf.loc[mod]-mndf.loc[mod],color='0.7')#,color=mycmap.colors[mii])#,linestyle='None',marker='s')

plt.plot(mxdf.loc['CanESM2']-mndf.loc['CanESM2'],'k',linewidth=3)#,linestyle='None',marker='s')
ax.set_xticks(np.arange(0,12))
ax.set_xlim((0,11))
ax.set_xticklabels(mons)
ax.set_title('Range in within-model SIE trends (1979-2012)')
if printtofile:
    fig.savefig('SIE_CMIP5_trendrange_seacyc_1979-2012_wCanESM2.pdf')


## for monkey in monmnmxdt.keys():

##     rng = monmnmxdt[monkey][1] - monmnmxdt[monkey][0]
##     print rng

##          #print ensnum

##     modata.plot()
    
##     #mondf[monii] = pd.DataFrame(modensdt) # a dictionary of DataFrames of the model ens (one for each month)

## #mondf = pd.DataFrame(mondt)

## #plt.figure()
## #mondf[8].plot()



###### want cmip5 multi-model mean sea ice concentration maps for
# 2002-12 minus 1979-89   @@@
# will have to concatenate rcp4.5 too

othermodels=('CMCC-CMS','CMCC-CM','EC-EARTH','FGOALS-g2',
             'GFDL-CM2p1','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G',
             'MIROC4h')

field='sic'
comp='OImon'
bp='/rd40/data/CMIP5/historical/' + field + '/'
timeselc='1979-01-01,1989-12-31'
timeselp1='2002-01-01,2005-12-31'
timeselp2='2006-01-01,2012-12-31'

## ctldt = dict.fromkeys(models)
## p1dt = dict.fromkeys(models) # 2002-2005

## for model in models:
##      enum=enumdt[model] # number of ens members
##      ctlmoddt={}; p1moddt={}
##      for eii in range(1,enum+1):
##           if model in (othermodels): # these models have 5 or 10-year files
##                print 'skipping ' + model
##           else: # add filenames to 
##                if model in ('HadCM3',):
##                     fname=bp+model+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
##                       model + '_historical_r' + str(eii)+'i1p1_195912-198411.nc'
##                elif model in ('HadGEM2-AO',):
##                     fname=bp+model+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
##                            model + '_historical_r' + str(eii)+'i1p1_186001-200512.nc'
##                elif model in ('HadGEM2-CC',):
##                     fname=bp+model+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
##                            model + '_historical_r' + str(eii)+'i1p1_195912-200511.nc'
##                elif model in ('HadGEM2-ES',):
##                     fname=bp+model+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
##                            model + '_historical_r' + str(eii)+'i1p1_195912-200512.nc'
##                elif model in ('MIROC5',):
##                     fname=bp+model+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
##                            model + '_historical_r' + str(eii)+'i1p1_185001-201212.nc'
##                elif model in ('MPI-ESM-MR',):
##                     fname=bp+model+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
##                            model + '_historical_r' + str(eii)+'i1p1_190001-199912.nc'
##                elif model in ('MRI-ESM1',):
##                     fname=bp+model+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
##                            model + '_historical_r' + str(eii)+'i1p1_185101-200512.nc'
##                elif model in ('GISS-E2-H',) and eii==6:
##                     fname=bp+model+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
##                            model + '_historical_r' + str(eii)+'i1p1_195101-200012.nc'
##                else:
##                     fname=bp+model+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
##                            model + '_historical_r' + str(eii)+'i1p1_185001-200512.nc'

##                print fname
##                ctlmoddt[eii],junk = cutl.climatologize(cnc.getNCvar(fname,field,timesel=timeselc))
##                p1moddt[eii],junk = cutl.climatologize(cnc.getNCvar(fname,field,timesel=timeselp1))
          
##      ctldt[model]=ctlmoddt
##      p1dt[model]=p1moddt

    

print len(cm5.models)
#bp='/ra40/data/ncs/sic/'


ctldt = dict.fromkeys(models)
p1dt = dict.fromkeys(models) # 2002-2005

for mod in cm5.models:

     print mod
     meta=cm5.allmodeldt[mod]
     if meta==None:
          print mod + ' not yet implemented!' # @@@

     elif type(meta['fyears'])!=str:
          # skip
          print ' skipping models w/ annoying years for now'
     else:
          fyears=meta['fyears']
          enum=meta['numens'] # number of ens members
          ctlmoddt={}; p1moddt={}
          for eii in range(1,enum+1):
               if 'skipens' in meta.keys() and eii in meta['skipens']:
                    print 'skipping ensemble: ' + str(eii)
               else:
                    fname=bp+mod+'/r'+str(eii)+'i1p1/' + field + '_' + comp + '_' +\
                                mod + '_historical_r' + str(eii)+'i1p1_' + fyears +'.nc'
                    print fname
                    ctlmoddt[eii],junk = cutl.climatologize(cnc.getNCvar(fname,field,timesel=timeselc))
                    p1moddt[eii],junk = cutl.climatologize(cnc.getNCvar(fname,field,timesel=timeselp1))
          
     ctldt[mod]=ctlmoddt
     p1dt[mod]=p1moddt
