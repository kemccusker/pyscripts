""" canesm_LE_general.py

    This script basically reads in files that are already processed
    into area-averages and such.

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

# exception handling works below. add to other clauses @@

printtofile=True
dohist=False
doscatter=True
addobs=True # to scatter plot
addnat=False
addsims=False # add the idealized simulations. only good for DJF polar amp vs eurasia SAT

performop1 = False
op1='div'; region1op='gm' # polar amp: gt60n / gm
performop2 = False
op2='sub'; region2op='nh'

timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'

field1='zg50000.00'; ncfield1='zg'; comp1='Amon'; region1='bksmori'
#field1='sia'; ncfield1='sianh'; comp1='OImon'; region1='nh'
#field1='tas'; ncfield1='tas'; comp1='Amon'; region1='gt60n'
leconv1= 1 
sea1='SON'

field2='tas'; ncfield2='tas'; comp2='Amon'; region2='eurasiamori'
leconv2=1
sea2='DJF'



ftype='fullts' # 'fullclimo' or 'climo' or 'fullts'


fdict1 = {'field': field1+region1, 'ncfield': ncfield1, 'comp': comp1}
fdict2 = {'field': field2+region2, 'ncfield': ncfield2, 'comp': comp2}

if doscatter:

    # historical
    casename='historical'
    lecdat1 = le.load_LEdata(fdict1,casename,timesel=timeselc, rettype='ndarray',conv=leconv1,ftype=ftype)
    (numens1,ntime1) = lecdat1.shape
    lepdat1=le.load_LEdata(fdict1,casename,timesel=timeselp, rettype='ndarray',conv=leconv1,ftype=ftype)
    lecsea1 = cutl.seasonalize_monthlyts(lecdat1.T,season=sea1).T
    lepsea1 = cutl.seasonalize_monthlyts(lepdat1.T,season=sea1).T
    lesea1 = lepsea1 - lecsea1
    
    if performop1:
        try:
            fdict1op = {'field': field1+region1op, 'ncfield': ncfield1, 'comp': comp1}
            # should rename these variables to be 'op' so it's more general
            subc1 = le.load_LEdata(fdict1op,casename,timesel=timeselc, rettype='ndarray',conv=leconv1,ftype=ftype)
            subp1 = le.load_LEdata(fdict1op,casename,timesel=timeselp, rettype='ndarray',conv=leconv1,ftype=ftype)
            sub1 = cutl.seasonalize_monthlyts(subp1.T,season=sea1).T - cutl.seasonalize_monthlyts(subc1.T,season=sea1).T

            if op1=='sub': # subtract
                lefld1 = lesea1.mean(axis=1) - sub1.mean(axis=1)
            elif op1=='div': # divide
                lefld1 = lesea1.mean(axis=1) / sub1.mean(axis=1)
            else:
                print 'operation not supported!'
                raise Exception
        except:
            raise

    else:
        lefld1=lepsea1.mean(axis=1)-lecsea1.mean(axis=1)

    lecdat2 = le.load_LEdata(fdict2,casename,timesel=timeselc, rettype='ndarray',conv=conv2,ftype=ftype)
    (numens2,ntime2) = lecdat1.shape
    lepdat2=le.load_LEdata(fdict2,casename,timesel=timeselp, rettype='ndarray',conv=conv2,ftype=ftype)
    lecsea2 = cutl.seasonalize_monthlyts(lecdat2.T,season=sea2).T
    lepsea2 = cutl.seasonalize_monthlyts(lepdat2.T,season=sea2).T
    lesea2 = lepsea2 - lecsea2
    if performop2:
        if op2=='sub': # subtract
            fdict2sub = {'field': field2+region2op, 'ncfield': ncfield2, 'comp': comp2}
            subc2 = le.load_LEdata(fdict2sub,casename,timesel=timeselc, rettype='ndarray',conv=leconv2,ftype=ftype)
            subp2 = le.load_LEdata(fdict2sub,casename,timesel=timeselp, rettype='ndarray',conv=leconv2,ftype=ftype)
            sub2 = cutl.seasonalize_monthlyts(subp2.T,season=sea2).T - cutl.seasonalize_monthlyts(subc2.T,season=sea2).T
            lefld2 = lesea2.mean(axis=1) - sub2.mean(axis=1)
    else:
        lefld2=lepsea2.mean(axis=1)-lecsea2.mean(axis=1)

    lemm, lebb, lerval, lepval, lestd_err = sp.stats.linregress(lefld1,lefld2)
    

    # historicalNat
    casename2='historicalNat'
    lecdat1n = le.load_LEdata(fdict1,casename2,timesel=timeselc, rettype='ndarray',conv=leconv1,ftype=ftype)
    (numens1,ntime1) = lecdat1n.shape
    lepdat1n=le.load_LEdata(fdict1,casename2,timesel=timeselp, rettype='ndarray',conv=leconv1,ftype=ftype)
    lecsea1n = cutl.seasonalize_monthlyts(lecdat1n.T,season=sea1).T
    lepsea1n = cutl.seasonalize_monthlyts(lepdat1n.T,season=sea1).T
    lesea1n = lepsea1n - lecsea1n
    if performop1:
        fdict1op = {'field': field1+region1op, 'ncfield': ncfield1, 'comp': comp1}
        subc1n = le.load_LEdata(fdict1op,casename2,timesel=timeselc, rettype='ndarray',conv=leconv1,ftype=ftype)
        subp1n = le.load_LEdata(fdict1op,casename2,timesel=timeselp, rettype='ndarray',conv=leconv1,ftype=ftype)
        sub1n = cutl.seasonalize_monthlyts(subp1n.T,season=sea1).T - cutl.seasonalize_monthlyts(subc1n.T,season=sea1).T

        if op1=='sub': # subtract
            lefld1n = lesea1n.mean(axis=1) - sub1n.mean(axis=1)
        elif op1=='div': # divide
            lefld1n = lesea1n.mean(axis=1) / sub1n.mean(axis=1)
    else:
        lefld1n=lepsea1n.mean(axis=1)-lecsea1n.mean(axis=1)

    lecdat2n = le.load_LEdata(fdict2,casename2,timesel=timeselc, rettype='ndarray',conv=conv2,ftype=ftype)
    (numens2,ntime2) = lecdat1n.shape
    lepdat2n=le.load_LEdata(fdict2,casename2,timesel=timeselp, rettype='ndarray',conv=conv2,ftype=ftype)
    lecsea2n = cutl.seasonalize_monthlyts(lecdat2n.T,season=sea2).T
    lepsea2n = cutl.seasonalize_monthlyts(lepdat2n.T,season=sea2).T
    lesea2n = lepsea2n - lecsea2n
    if performop2:
        if op2=='sub': # subtract
            fdict2sub = {'field': field2+region2op, 'ncfield': ncfield2, 'comp': comp2}
            subc2n = le.load_LEdata(fdict2sub,casename2,timesel=timeselc, rettype='ndarray',conv=leconv2,ftype=ftype)
            subp2n = le.load_LEdata(fdict2sub,casename2,timesel=timeselp, rettype='ndarray',conv=leconv2,ftype=ftype)
            sub2n = cutl.seasonalize_monthlyts(subp2n.T,season=sea2).T - cutl.seasonalize_monthlyts(subc2n.T,season=sea2).T
            lefld2n = lesea2n.mean(axis=1) - sub2n.mean(axis=1)
    else:
        lefld2n=lepsea2n.mean(axis=1)-lecsea2n.mean(axis=1)

    lemmn, lebbn, lervaln, lepvaln, lestd_errn = sp.stats.linregress(lefld1n,lefld2n)

    if addobs: # DATA MUST EXIST:hard coded for TAS and Z500 @@
        graveraint= 9.80665 # m/s2 (different from Canadian models)

        if field1 == 'zg50000.00':
            erafile = '/HOME/rkm/work/DATA/ERAINT/td_era_int_197901_201412_gp_128_64_phi500_1979011612-2014121612.nc'
            eraz500c= cnc.getNCvar(erafile,'PHI',timesel=timeselc,seas=sea1)/graveraint
            eraz500p= cnc.getNCvar(erafile,'PHI',timesel=timeselp,seas=sea1)/graveraint
            latera=cnc.getNCvar(erafile,'lat')
            lonera=cnc.getNCvar(erafile,'lon')
            obsreg1 = cutl.calc_regmean(eraz500p-eraz500c,latera,lonera,region1)
            obsreg1=obsreg1.mean()
        elif field1 == 'tas':
            gisfile = '/HOME/rkm/work/DATA/GISS/gistemp1200_ERSST.nc'
            latgis=cnc.getNCvar(gisfile,'lat')
            longis=cnc.getNCvar(gisfile,'lon')
            gissatc= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselc,seas=sea2) 
            gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselp,seas=sea2)
            obsreg1 =  cutl.calc_regmean(gissatp-gissatc,latgis,longis,region1)
            if performop1:
                opfld1 = cutl.calc_regmean(gissatp-gissatc,latgis,longis,region1op)

                if op1=='sub': # subtract
                    obsreg1 = obsreg1.mean() - opfld1.mean()
                elif op1=='div': # divide
                    obsreg1 = obsreg1.mean() / opfld1.mean()
            else:
                obsreg1=obsreg1.mean()

        if field2 == 'tas':
            gisfile = '/HOME/rkm/work/DATA/GISS/gistemp1200_ERSST.nc'
            latgis=cnc.getNCvar(gisfile,'lat')
            longis=cnc.getNCvar(gisfile,'lon')
            gissatc= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselc,seas=sea2) 
            gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselp,seas=sea2)
            
            if sea1==sea2=='DJF' and performop2==True and op2=='sub' and region2=='eurasiamori' and region2op=='nh':
                # a file already exists for this
                gisfile='/HOME/rkm/work/DATA/GISS/giss_DJF_eurasiamori-nh_1979-2013_timeseries.nc'
                gissatc= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselc) 
                gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselp)
                obsreg2=gissatp.mean()-gissatc.mean()
            else:
                obsreg2 =  cutl.calc_regmean(gissatp-gissatc,latgis,longis,region2)
                obsreg2=obsreg2.mean()
            
    fig,ax=plt.subplots(1,1)
    ledat=plt.scatter(lefld1,lefld2,color=ccm.get_linecolor('darkolivegreen3'),marker='*',s=5**2,alpha=0.5)
    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    onex=np.linspace(axxlims[0],axxlims[1])
    ax.plot(onex,lemm*onex + lebb, color=ccm.get_linecolor('darkolivegreen3'),linewidth=1)

    obstr=''
    leghnds=(ledat,)
    legstrs=(casename + ' LE',)

    if addnat:
        ledatn=plt.scatter(lefld1n,lefld2n,color=ccm.get_linecolor('steelblue3'),marker='*',s=5**2,alpha=0.5)
        axylims = ax.get_ylim()
        axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,lemmn*onex + lebbn, color=ccm.get_linecolor('steelblue3'),linewidth=1)
        leghnds=leghnds + (ledatn,)
        legstrs=legstrs + (casename2 + ' LE',)
    if addobs:
        obs=plt.scatter(obsreg1,obsreg2,color='blue',marker='s',s=8**2)
        leghnds=leghnds + (obs,)
        legstrs=legstrs + ('Observations',)
        obstr='obs'

    if addsims and sea1==sea2=='DJF' and performop1==True and op1=='div' and region1=='gt60n' and region1op=='gm':
        # calculated using, for example:
        # lmd.loaddata(('st',),('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5'),
        #              timefreq='DJF',meantype='time',region='eurasiamori')
        # polaramp is gt60n / gm

        # these are DJF values for E1-5 and R1-5
        simpolamp = [[  9.98685376],
                     [ 11.2359957 ],
                     [ 13.62247628],
                     [ 12.03047158],
                     [ 12.47854705],
                     [ 10.8695116 ],
                     [ 13.04175298],
                     [ 10.06997661],
                     [ 10.7259977 ],
                     [ 10.69236807]]
        simeur = [[ 0.27966023],
                  [-0.03198602],
                  [ 0.12460408],
                  [-0.24664229],
                  [-0.03082406],
                  [ 0.12214981],
                  [-0.10484801],
                  [-0.00995005],
                  [ 0.35847141],
                  [ 0.27739474]]
        sims = ax.scatter(simpolamp,simeur,color='0.3',marker='o',s=6**2,alpha=0.7)
        leghnds=leghnds + (sims,)
        legstrs=legstrs + ('Modelled SIC forcing',)
        obstr=obstr+'sims'

    fontP = fm.FontProperties()
    fontP.set_size('small')
    plt.legend(leghnds,legstrs,
               loc='best',fancybox=True,framealpha=0.5,prop=fontP)#,frameon=False) 

    plt.annotate(casename + ' R= ' + '$%.2f$'%(lerval) + ', p='+ '$%.2f$'%(lepval),
                xy=(axxlims[0]+.05*axxlims[1], axylims[1]-.1*axylims[1]),
                xycoords='data')
    if addnat:
        plt.annotate(casename2 + ' R= ' + '$%.2f$'%(lervaln) + ', p='+ '$%.2f$'%(lepvaln),
                    xy=(axxlims[0]+.05*axxlims[1], axylims[1]-.2*axylims[1]),
                    xycoords='data')

    opstr1 = opstr2 = ''
    if performop1:
        opstr1= op1 + '_' + region1op
    plt.xlabel(sea1 + ' ' + field1 + ' ' + region1 + ' ' + opstr1)
    if performop2:
        opstr2= op2 + '_' + region2op
    plt.ylabel(sea2 + ' ' + field2 + ' ' + region2 + ' ' + opstr2)

    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5,linestyle='--')
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5,linestyle='--')

    if printtofile:
        if addnat:
            plt.savefig('scatterregress_' + field1 + region1 + opstr1 + str(sea1) + '_v_' + 
                        field2 + region2 +opstr2 + str(sea2) + '_' + casename + '_' + 
                        casename2 + '_2002-12_1979-89' + obstr + '.pdf')
        else:
            plt.savefig('scatterregress_' + field1 + region1 + opstr1 + str(sea1) + '_v_' + 
                        field2 + region2 +opstr2 + str(sea2) + '_' + casename + 
                        '_2002-12_1979-89' + obstr + '.pdf')

if dohist:
    histcdat=le.load_LEdata(fdict2,'historical',timesel=timeselc, rettype='ndarray',conv=conv,ftype=ftype)
    (numens,ntime) = histcdat.shape
    histpdat=le.load_LEdata(fdict2,'historical',timesel=timeselp, rettype='ndarray',conv=conv,ftype=ftype)

    # Now have 11 years of monthly data. Grab DJF:
    histc = cutl.seasonalize_monthlyts(histcdat.T,season='DJF').T
    histp = cutl.seasonalize_monthlyts(histpdat.T,season='DJF').T

    histnatcdat=le.load_LEdata(fdict2,'historicalNat',timesel=timeselc, rettype='ndarray',conv=conv,ftype=ftype)
    #(numens,ntime,nlatlon) = histnatcdat.shape
    histnatpdat=le.load_LEdata(fdict2,'historicalNat',timesel=timeselp, rettype='ndarray',conv=conv,ftype=ftype)
    histnatc = cutl.seasonalize_monthlyts(histnatcdat.T,season='DJF').T
    histnatp = cutl.seasonalize_monthlyts(histnatpdat.T,season='DJF').T

    firebrick=ccm.get_linecolor('firebrick')
    darkolivegreen3=ccm.get_linecolor('darkolivegreen3')
    steelblue3=ccm.get_linecolor('steelblue3')

    plt.figure()
    plt.hist(histp.mean(axis=1)-histc.mean(axis=1),color=darkolivegreen3,alpha=0.5)
    plt.hist(histnatp.mean(axis=1)-histnatc.mean(axis=1),color=steelblue3,alpha=0.5)
    plt.title('CanESM LE: ' + sea + ' ' + '2002-02 - 1979-89')
    plt.xlabel(field2 + ' ' + region2)
    if printtofile:
        plt.savefig(field2+'_' + region2 + 'historical_historicalNat_2002-12_1979-89_PDF_' + sea + '.pdf')



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
