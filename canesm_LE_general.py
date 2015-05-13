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
import matplotlib.font_manager as fm
import loadmodeldata as lmd

# exception handling works below. add to other clauses @@

printtofile=True
dohist=False
doregress=True
doscatter=False
addobs=False # to scatter plot
addnat=True
addsims=True # add the idealized simulations. only good for DJF polar amp vs eurasia SAT

performop1 = False
#op1='div'; region1op='gm' # polar amp: gt60n / gm
op1='sub'; region1op='deeptrop' # pole-eq temp gradient: gt60n - deeptrop (or trop)
performop2 = True
op2='sub'; region2op='nh'

timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'

#field1='zg50000.00'; ncfield1='zg'; comp1='Amon'; region1='bksmori'
#field1='sia'; ncfield1='sianh'; comp1='OImon'; region1='nh'
field1='sic'; ncfield1='sic'; comp1='OImon'; region1='bksmori' # @@ a hack. prefer SIA
#field1='tas'; ncfield1='tas'; comp1='Amon'; region1='gt60n' #region1='bksmori'
leconv1= 1 
sea1='DJF'

field2='tas'; ncfield2='tas'; comp2='Amon'; region2='eurasiamori'
#field2='zg50000.00'; ncfield2='zg'; comp2='Amon'; region2='bksmori'
leconv2=1
sea2='DJF'

xlab=ylab=None
#xlab = '$\Delta$ DJF Barents-Kara Seas Z500 (m)'
#ylab = '$\Delta$ DJF SAT(Eurasia) - SAT(NH) ($^\circ$C)'
#ylab = '$\Delta$ DJF Eurasian SAT ($^\circ$C)'

simconv1=simconv2=1
if field1=='tas': simfield1='st'; simncfield1='ST'
elif field1=='zg50000.00': simfield1='gz50000'; simncfield1='PHI'; simconv1=1/con.get_g()
elif field1=='sia': simfield1='sicn'; simncfield1='SICN'; print '@@ danger, sia actually sicn average'
else: print 'cannot addsims for ' + field1; addsims=False

if field2=='tas': simfield2='st'; simncfield2='ST'
elif field2=='zg50000.00': simfield2='gz50000'; simncfield2='PHI'; simconv2=1/con.get_g()
elif field2=='sia': simfield2='sicn'; simncfield2='SICN'; print '@@ danger, sia actually sicn average'
else: print 'cannot addsims for ' + field2; addsims=False


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

    lecdat2 = le.load_LEdata(fdict2,casename,timesel=timeselc, rettype='ndarray',conv=leconv2,ftype=ftype)
    (numens2,ntime2) = lecdat1.shape
    lepdat2=le.load_LEdata(fdict2,casename,timesel=timeselp, rettype='ndarray',conv=leconv2,ftype=ftype)
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

    lecdat2n = le.load_LEdata(fdict2,casename2,timesel=timeselc, rettype='ndarray',conv=leconv2,ftype=ftype)
    (numens2,ntime2) = lecdat1n.shape
    lepdat2n=le.load_LEdata(fdict2,casename2,timesel=timeselp, rettype='ndarray',conv=leconv2,ftype=ftype)
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
        else:
            print 'this field not supported with addobs @@: ' + field1; addobs=False

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
        else:
            print 'this field not supported with addobs @@: ' + field2; addobs=False
            
    fig,ax=plt.subplots(1,1)
    ledat=plt.scatter(lefld1,lefld2,color=ccm.get_linecolor('darkolivegreen3'),marker='*',s=5**2,alpha=0.5)
    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    onex=np.linspace(axxlims[0],axxlims[1])
    ax.plot(onex,lemm*onex + lebb, color=ccm.get_linecolor('darkolivegreen3'),linewidth=1)
    

    obstr=''
    leghnds=(ledat,)
    legstrs=(casename + ' LE',)

    print '-------- LE slope, rval, pval: ' + str(lemm),str(lerval),str(lepval)
    if addnat:
        ledatn=plt.scatter(lefld1n,lefld2n,color=ccm.get_linecolor('steelblue3'),marker='*',s=5**2,alpha=0.5)
        axylims = ax.get_ylim()
        axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,lemmn*onex + lebbn, color=ccm.get_linecolor('steelblue3'),linewidth=1)
        leghnds=leghnds + (ledatn,)
        legstrs=legstrs + (casename2 + ' LE',)
        print '-------- LENat slope, rval, pval: ' + str(lemmn),str(lervaln),str(lepvaln)
    if addobs:
        obs=plt.scatter(obsreg1,obsreg2,color='blue',marker='s',s=8**2)
        leghnds=leghnds + (obs,)
        legstrs=legstrs + ('Observations',)
        obstr='obs'

    """if addsims and sea1==sea2=='DJF' and performop1==True and op1=='div' and region1=='gt60n' and region1op=='gm':
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
        simsh = ax.scatter(simpolamp,simeur,color='0.3',marker='o',s=6**2,alpha=0.7)
        leghnds=leghnds + (simsh,)
        legstrs=legstrs + ('Modelled SIC forcing',)
        obstr=obstr+'sims'"""
    if addsims:
        sims=('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5')
        flddf1 = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea1, 
                                           meantype='time',region=region1))*simconv1

        if performop1:
            flddf1op = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea1, 
                                                 meantype='time',region=region1op))*simconv1

            if op1=='sub': # subtract
                flddf1 = flddf1 - flddf1op
            elif op1=='div': # divide
                flddf1 = flddf1 / flddf1op

        flddf2 = pd.DataFrame(lmd.loaddata((simfield2,),sims,ncfields=(simncfield2,), timefreq=sea2, 
                                           meantype='time',region=region2))*simconv2

        if performop2:
            flddf2op = pd.DataFrame(lmd.loaddata((simfield2,),sims,ncfields=(simncfield2,), timefreq=sea2, 
                                                 meantype='time',region=region2op))*simconv2

            if op2=='sub': # subtract
                flddf2 = flddf2 - flddf2op
            elif op2=='div': # divide
                flddf2 = flddf2 / flddf2op

        simmm, simbb, simrval, simpval, simstd_err = sp.stats.linregress(np.squeeze(flddf1.values),np.squeeze(flddf2.values))
        print '-------- SIMS slope, rval, pval: ' + str(simmm),str(simrval),str(simpval)

        simsh = ax.scatter(flddf1,flddf2,color='0.3',marker='o',s=6**2,alpha=0.7)
        #axylims = ax.get_ylim()
        #axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,simmm*onex + simbb, color='0.3',linewidth=1)
        leghnds=leghnds + (simsh,)
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
    if xlab==None:
        plt.xlabel(sea1 + ' ' + field1 + ' ' + region1 + ' ' + opstr1)
    else:
        plt.xlabel(xlab)
    if performop2:
        opstr2= op2 + '_' + region2op

    if ylab==None:
        plt.ylabel(sea2 + ' ' + field2 + ' ' + region2 + ' ' + opstr2)
    else:
        plt.ylabel(ylab)

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
    conv=leconv1 # just assume we are doing variable 1
    sea=sea1

    histcdat=le.load_LEdata(fdict2,'historical',timesel=timeselc, rettype='ndarray',conv=conv,ftype=ftype)
    (numens,ntime) = histcdat.shape
    histpdat=le.load_LEdata(fdict2,'historical',timesel=timeselp, rettype='ndarray',conv=conv,ftype=ftype)

    # Now have 11 years of monthly data. Grab DJF:
    histc = cutl.seasonalize_monthlyts(histcdat.T,season=sea1).T
    histp = cutl.seasonalize_monthlyts(histpdat.T,season=sea1).T

    histnatcdat=le.load_LEdata(fdict2,'historicalNat',timesel=timeselc, rettype='ndarray',conv=conv,ftype=ftype)
    #(numens,ntime,nlatlon) = histnatcdat.shape
    histnatpdat=le.load_LEdata(fdict2,'historicalNat',timesel=timeselp, rettype='ndarray',conv=conv,ftype=ftype)
    histnatc = cutl.seasonalize_monthlyts(histnatcdat.T,season=sea1).T
    histnatp = cutl.seasonalize_monthlyts(histnatpdat.T,season=sea1).T

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

if doregress:
    donorm=True

    seasp=sea1 # season of spatial field
    sear=sea2 # season of regional avgs


    # spatial field
    leconvsp=1
    fieldsp='zg50000.00'; ncfieldsp='zg'; compsp='Amon'; 
    #fieldsp='tas'; ncfieldsp='tas'; compsp='Amon'; 

    # regional avg field 1
    leconvr=-1 # this way, sea ice loss is linked with positive changes elsewhere
    fieldr='sic'; ncfieldr='sic'; compr='OImon'; regionr='bksmori'

    # regional avg field 2
    leconvr2=1
    fieldr2='tas'; ncfieldr2='tas'; compr2='Amon'; regionr2='eurasiamori'; leconvr2=-1 # so cooling=high heights
    #fieldr2='zg50000.00'; ncfieldr2='zg'; compr2='Amon'; regionr2='bksmori'
    

    # what are the units of these regressions? @@
    if fieldsp=='zg50000.00':
        # these are probably : m/% and m/C
        # try normalizing the 1D field by its sigma so that plots will be m/SD
        cmin=-3; cmax=3 # I think m/%
        cmin2=-16; cmax2=16 # m/C
        normstr=''
        if donorm:
            cmin=-10; cmax=10
            cmin2=-10; cmax2=10
            normstr='norm'
    elif fieldsp=='tas':
        # these are probably : C/% and C/m
        # try normalizing the 1D field by its sigma so that plots will be C/SD
        cmin=-.5; cmax=.5 # I think C/%
        cmin2=-.1; cmax2=.1 # C/m
        normstr=''
        if donorm:
            cmin=-.7; cmax=.7
            cmin2=-.7; cmax2=.7
            normstr='norm'


    fdictsp = {'field': fieldsp, 'ncfield': ncfieldsp, 'comp': compsp}
    fdictr = {'field': fieldr+regionr, 'ncfield': ncfieldr, 'comp': compr}
    fdictr2 = {'field': fieldr2+regionr2, 'ncfield': ncfieldr2, 'comp': compr2}

    casename = 'historical'

    lat=le.get_lat()
    lon=le.get_lon()
    nlat=len(lat); nlon=len(lon)

    # LOAD SPATIAL DATA
    lecdatsp = le.load_LEdata(fdictsp,casename,timesel=timeselc, rettype='ndarray',conv=leconvsp,ftype=ftype)
    (numensp,ntimesp,nspacesp) = lecdatsp.shape
    lepdatsp=le.load_LEdata(fdictsp,casename,timesel=timeselp, rettype='ndarray',conv=leconvsp,ftype=ftype)
    # time needs to be first dimension
    lecdatsp = np.transpose(lecdatsp,(1,0,2))
    lepdatsp = np.transpose(lepdatsp,(1,0,2))
    lecseasp = cutl.seasonalize_monthlyts(lecdatsp,season=seasp).mean(axis=0) # average the 11 seasonal avgs 
    lepseasp = cutl.seasonalize_monthlyts(lepdatsp,season=seasp).mean(axis=0) # average the 11 seasonal avgs 
    leseasp = lepseasp - lecseasp # numens x space.flat
    
    # LOAD 1D DATA
    lecdatr = le.load_LEdata(fdictr,casename,timesel=timeselc, rettype='ndarray',conv=leconvr,ftype=ftype)
    (numenr,ntimer) = lecdatr.shape
    lepdatr=le.load_LEdata(fdictr,casename,timesel=timeselp, rettype='ndarray',conv=leconvr,ftype=ftype)
    lecsear = cutl.seasonalize_monthlyts(lecdatr.T,season=sear).mean(axis=0)
    lepsear = cutl.seasonalize_monthlyts(lepdatr.T,season=sear).mean(axis=0)
    lesear = lepsear - lecsear # numens
    if donorm:
        lesear=lesear / lesear.std()

    #mm, bb, rval, pval, std_err = sp.stats.linregress(lesear,leseasp) # errors. why? @@
    slope,intercept = np.polyfit(lesear,leseasp,1)
    bkszg=slope.reshape((nlat,nlon))

    # LOAD 1D DATA (2)
    lecdatr2 = le.load_LEdata(fdictr2,casename,timesel=timeselc, rettype='ndarray',conv=leconvr2,ftype=ftype)
    (numenr2,ntimer2) = lecdatr2.shape
    lepdatr2=le.load_LEdata(fdictr2,casename,timesel=timeselp, rettype='ndarray',conv=leconvr2,ftype=ftype)
    lecsear2 = cutl.seasonalize_monthlyts(lecdatr2.T,season=sear).mean(axis=0)
    lepsear2 = cutl.seasonalize_monthlyts(lepdatr2.T,season=sear).mean(axis=0)
    lesear2 = lepsear2 - lecsear2 # numens
    if donorm:
        lesear2=lesear2 / lesear2.std()

    #mm, bb, rval, pval, std_err = sp.stats.linregress(lesear,leseasp) # errors. why? @@
    slope,intercept = np.polyfit(lesear2,leseasp,1)
    eurzg=slope.reshape((nlat,nlon))




    fig,axs=plt.subplots(1,2)
    fig.set_size_inches(10,5)
    ax=axs[0]#
    cplt.kemmap(bkszg,lat,lon,type='nh',axis=ax,cmin=cmin,cmax=cmax,
                title=sear + ' ' +fieldr+regionr + ' regressed onto ' + seasp + ' ' + fieldsp)

    ax=axs[1] #
    cplt.kemmap(eurzg,lat,lon,type='nh',axis=ax, cmin=cmin2,cmax=cmax2,
                title=sear + ' ' +fieldr2+regionr2 + ' regressed onto ' + seasp + ' ' + fieldsp)

    if printtofile:
        fig.savefig(fieldr +regionr + '_' + fieldr2 + regionr2 + sear \
                    + '_regresson_' + fieldsp + seasp + normstr + '.pdf') 


    """    #Got memory errors when trying to read in everything @@
           #But one LE works
    flist=le.build_filenames(fdict,casename,ftype='fullts')
    numens=len(flist)
    fname = flist[0]
    lat=cnc.getNCvar(fname,'lat')
    nlat=len(lat)
    lon=cnc.getNCvar(fname,'lon')
    nlon=len(lon)

    histcdat = histcdat.reshape((numens,ntime,nlat,nlon))
    histpdat = histpdat.reshape((numens,ntime,nlat,nlon))

    #histnatcdat = histnatcdat.reshape((numens,ntime,nlat,nlon))
    #histnatpdat = histnatpdat.reshape((numens,ntime,nlat,nlon))

    # Now, seasonal average and regional average here
    nyr = ntime/12.-1 # @@@@@@ hack for laptop. data probably not correct @@@@ actually needed on linux too...

    #histcsea = np.zeros((numens,nyr,nlat,lon)
    histreg = np.zeros((numens,nyr))
    #histnatreg = np.zeros((numens,nyr))

    for nii in range(0,numens):

        histcsea = cutl.seasonalize_monthlyts(histcdat[nii,...],season=sea)
        histpsea = cutl.seasonalize_monthlyts(histpdat[nii,...],season=sea)
        ctl = cutl.calc_regmean(histcsea,lat,lon,region=region)
        pert = cutl.calc_regmean(histpsea,lat,lon,region=region)
        histreg[nii,...] = pert - ctl

        #histnatcsea = cutl.seasonalize_monthlyts(histnatcdat[nii,...],season=sea)
        #histnatpsea = cutl.seasonalize_monthlyts(histnatpdat[nii,...],season=sea)
        #ctlnat = cutl.calc_regmean(histnatcsea,lat,lon,region=region)
        #pertnat = cutl.calc_regmean(histnatpsea,lat,lon,region=region)
        #histnatreg[nii,...] = pertnat - ctlnat

    plt.figure()
    plt.hist(histreg.mean(axis=1),alpha=0.5)
    #plt.hist(histnatreg.mean(axis=1),color='.5')"""
