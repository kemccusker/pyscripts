""" canesm_LE_general.py

"""
import pandas as pd
import cccmaplots as cplt
import cccmautils as cutl
import cccmacmaps as ccm
import datetime as datetime
from netCDF4 import num2date
import sys
import loadLE as le
import constants as con

printtofile=False
dohist=False
doscatter=True

timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'

field1='zg50000.00'; ncfield1='zg'; comp1='Amon'; region1='bksmori'
conv1= 1 
sea1='DJF'

field2='tas'; ncfield2='tas'; comp2='Amon'; region2='eurasiamori'
conv2=1
sea2='DJF'

#field='sia'; ncfield='sianh'; comp='OImon'; region='nh'

ftype='fullts' # 'fullclimo' or 'climo' or 'fullts'


fdict1 = {'field': field1+region1, 'ncfield': ncfield1, 'comp': comp1}
fdict2 = {'field': field2+region2, 'ncfield': ncfield2, 'comp': comp2}

if doscatter:

    cdat1 = le.load_LEdata(fdict1,'historical',timesel=timeselc, rettype='ndarray',conv=conv1,ftype=ftype)
    (numens1,ntime1) = cdat1.shape
    pdat1=le.load_LEdata(fdict1,'historical',timesel=timeselp, rettype='ndarray',conv=conv1,ftype=ftype)
    csea1 = cutl.seasonalize_monthlyts(cdat1.T,season=sea1).T
    psea1 = cutl.seasonalize_monthlyts(pdat1.T,season=sea1).T
    fld1=psea1.mean(axis=1)-csea1.mean(axis=1)

    cdat2 = le.load_LEdata(fdict2,'historical',timesel=timeselc, rettype='ndarray',conv=conv2,ftype=ftype)
    (numens2,ntime2) = cdat1.shape
    pdat2=le.load_LEdata(fdict2,'historical',timesel=timeselp, rettype='ndarray',conv=conv2,ftype=ftype)
    csea2 = cutl.seasonalize_monthlyts(cdat2.T,season=sea2).T
    psea2 = cutl.seasonalize_monthlyts(pdat2.T,season=sea2).T
    fld2=psea2.mean(axis=1)-csea2.mean(axis=1)


    mm, bb, rval, pval, std_err = sp.stats.linregress(fld1,fld2)

    fig,ax=plt.subplots(1,1)
    ax.scatter(fld1,fld2)
    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    onex=np.linspace(axxlims[0],axxlims[1])
    ax.plot(onex,mm*onex + bb, color='0.5',linewidth=1)#,linestyle='--')



if dohist:
    histcdat=le.load_LEdata(fdict,'historical',timesel=timeselc, rettype='ndarray',conv=conv,ftype=ftype)
    (numens,ntime) = histcdat.shape
    histpdat=le.load_LEdata(fdict,'historical',timesel=timeselp, rettype='ndarray',conv=conv,ftype=ftype)

    # Now have 11 years of monthly data. Grab DJF:
    histc = cutl.seasonalize_monthlyts(histcdat.T,season='DJF').T
    histp = cutl.seasonalize_monthlyts(histpdat.T,season='DJF').T

    histnatcdat=le.load_LEdata(fdict,'historicalNat',timesel=timeselc, rettype='ndarray',conv=conv,ftype=ftype)
    #(numens,ntime,nlatlon) = histnatcdat.shape
    histnatpdat=le.load_LEdata(fdict,'historicalNat',timesel=timeselp, rettype='ndarray',conv=conv,ftype=ftype)
    histnatc = cutl.seasonalize_monthlyts(histnatcdat.T,season='DJF').T
    histnatp = cutl.seasonalize_monthlyts(histnatpdat.T,season='DJF').T

    firebrick=ccm.get_linecolor('firebrick')

    plt.figure()
    plt.hist(histp.mean(axis=1)-histc.mean(axis=1),color=firebrick,alpha=0.5)
    plt.hist(histnatp.mean(axis=1)-histnatc.mean(axis=1),color='b',alpha=0.5)
    plt.title('CanESM LE: ' + sea + ' ' + '2002-02 - 1979-89')
    plt.xlabel(field + ' ' + region)
    if printtofile:
        plt.savefig(field+'_' + region + 'historical_historicalNat_2002-12_1979-89_PDF_' + sea + '.pdf')



    natdf=pd.DataFrame(histnatp.mean(axis=1)-histnatc.mean(axis=1))
    histdf=pd.DataFrame(histp.mean(axis=1)-histc.mean(axis=1))

    #fig,ax=plt.subplots(1,1)
    histdf.hist(normed=True,color='0.5',alpha=0.5)#,histtype='stepfilled')
    natdf.hist(normed=True,color=firebrick,alpha=0.5)



""" Got memory errors when trying to read in everything @@
flist=le.build_filenames(fdict,'historical',ftype='fullts')
numens=len(flist)
fname = flist[0]
lat=cnc.getNCvar(fname,'lat')
nlat=len(lat)
lon=cnc.getNCvar(fname,'lon')
nlon=len(lon)

histcdat = histcdat.reshape((numens,ntime,nlat,nlon))
histpdat = histpdat.reshape((numens,ntime,nlat,nlon))

histnatcdat = histnatcdat.reshape((numens,ntime,nlat,nlon))
histnatpdat = histnatpdat.reshape((numens,ntime,nlat,nlon))

# Now, seasonal average and regional average here
nyr = ntime/12.-1 # @@@@@@ hack for laptop. data probably not correct @@@@ actually needed on linux too...

#histcsea = np.zeros((numens,nyr,nlat,lon)
histreg = np.zeros((numens,nyr))
histnatreg = np.zeros((numens,nyr))

for nii in range(0,numens):

    histcsea = cutl.seasonalize_monthlyts(histcdat[nii,...],season=sea)
    histpsea = cutl.seasonalize_monthlyts(histpdat[nii,...],season=sea)
    ctl = cutl.calc_regmean(histcsea,lat,lon,region=region)
    pert = cutl.calc_regmean(histpsea,lat,lon,region=region)
    histreg[nii,...] = pert - ctl

    histnatcsea = cutl.seasonalize_monthlyts(histnatcdat[nii,...],season=sea)
    histnatpsea = cutl.seasonalize_monthlyts(histnatpdat[nii,...],season=sea)
    ctlnat = cutl.calc_regmean(histnatcsea,lat,lon,region=region)
    pertnat = cutl.calc_regmean(histnatpsea,lat,lon,region=region)
    histnatreg[nii,...] = pertnat - ctlnat

plt.figure()
plt.hist(histreg.mean(axis=1),alpha=0.5)
plt.hist(histnatreg.mean(axis=1),color='.5')"""
