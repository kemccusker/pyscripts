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
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from scipy.stats import norm
import random
import numpy.ma as ma

# exception handling works below. add to other clauses @@

printtofile=False
dohist=False
doregress=False
doscatter=True
dolongtermavg=False

addobs=True # to scatter plot
addnat=False
addsims=True # add the idealized simulations. only good for DJF polar amp vs eurasia SAT

local=False

conditional=True # plot scatter conditional on 3rd var

performop1 = False
#op1='div'; region1op='gm' # polar amp: gt60n / gm
op1='sub'; region1op='deeptrop' # pole-eq temp gradient: gt60n - deeptrop (or trop)
performop2 = False
op2='sub'; region2op='nh'

timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'
timeselall = '1979-01-01,2012-12-31'

# x field
field1='zg50000.00'; ncfield1='zg'; comp1='Amon'; region1='bksmori'
#field1='sia'; ncfield1='sianh'; comp1='OImon'; region1='nh'
#field1='sic'; ncfield1='sic'; comp1='OImon'; region1='bksmori' # @@ a hack. prefer SIA
#field1='tas'; ncfield1='tas'; comp1='Amon'; region1='eurasiamori' #region1='bksmori'
leconv1= 1 
sea1='DJF'

# y field
field2='tas'; ncfield2='tas'; comp2='Amon'; region2='eurasiamori'
#field2='zg50000.00'; ncfield2='zg'; comp2='Amon'; region2='bksmori'
leconv2=1
sea2='DJF'

fieldcnd = 'sic'; ncfieldcnd='sic'; compcnd='OImon'; regioncnd='bksmori'
leconvcnd=1
seacnd=sea1
cmincnd=-9; cmaxcnd=3 # for 4 color colorbar

xlab=ylab=None
if field2=='tas' and region2=='eurasiamori' and performop2==True and op2=='sub':
    ylab = '$\Delta$ DJF SAT(Eurasia) - SAT(NH) ($^\circ$C)'
elif field2=='tas' and region2=='eurasiamori' and performop2==False:
    ylab = '$\Delta$ DJF Eurasian SAT ($^\circ$C)'
if field1=='zg50000.00' and region1=='bksmori':
    xlab = '$\Delta$ DJF Barents-Kara Seas Z500 (m)'


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
fdictcnd = {'field': fieldcnd+regioncnd, 'ncfield': ncfieldcnd, 'comp': compcnd}



def subsamp_sims(simsdf,numyrs=11):
    """ select 11 year segments from given simsdf (data should be anoms)
             simsdf.keys(): sims
             simsdf.index(): time index

             number of total samples will be determined by length of all sim data
                      and numyrs (e.g. (ntime / numyrs)*numsims)

             returns ndarray of subsample averages
    """
    #print simsdf

    ntime,numsims=simsdf.values.shape
    samp = ntime/numyrs
    allsii=0 # keep track of all sims and all subsamps
    subsampavg=np.zeros(samp*numsims)

    print 'sample each of ' + str(numsims) + ' sims ' + str(samp) + ' times'

    for nii,sim in enumerate(simsdf.keys()):
        
        vals=simsdf[sim].values     
        
        # random index to start looping, since we have a remainder when ntime/numyrs
        startyr = np.random.randint(np.mod(ntime,numyrs))
        print 'start ' + str(startyr)
        for sii in np.arange(startyr,ntime-numyrs,numyrs):

            #print 'sii ' + str(sii) + ', allsii ' + str(allsii)
            #subsampavg[allsii] = vals[sii+sii*numyrs:sii+sii*numyrs+numyrs].mean()
            subsampavg[allsii] = vals[sii:sii+numyrs].mean()
            allsii+=1

    return subsampavg



if doscatter:
    printtofile=True

    shorttermsims=True
    ymin=-2; ymax=3 # for Eurasia SAT v BKS Z500

    # historical
    casename='historical'
    lecdat1 = le.load_LEdata(fdict1,casename,timesel=timeselc, rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
    (numens1,ntime1) = lecdat1.shape
    lepdat1=le.load_LEdata(fdict1,casename,timesel=timeselp, rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
    lecsea1 = cutl.seasonalize_monthlyts(lecdat1.T,season=sea1).T
    lepsea1 = cutl.seasonalize_monthlyts(lepdat1.T,season=sea1).T
    lesea1 = lepsea1 - lecsea1
    
    if performop1:
        try:
            fdict1op = {'field': field1+region1op, 'ncfield': ncfield1, 'comp': comp1}
            # should rename these variables to be 'op' so it's more general
            subc1 = le.load_LEdata(fdict1op,casename,timesel=timeselc, rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
            subp1 = le.load_LEdata(fdict1op,casename,timesel=timeselp, rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
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

    lecdat2 = le.load_LEdata(fdict2,casename,timesel=timeselc, rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
    (numens2,ntime2) = lecdat2.shape
    lepdat2=le.load_LEdata(fdict2,casename,timesel=timeselp, rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
    lecsea2 = cutl.seasonalize_monthlyts(lecdat2.T,season=sea2).T
    lepsea2 = cutl.seasonalize_monthlyts(lepdat2.T,season=sea2).T
    lesea2 = lepsea2 - lecsea2
    if performop2:
        if op2=='sub': # subtract
            fdict2sub = {'field': field2+region2op, 'ncfield': ncfield2, 'comp': comp2}
            subc2 = le.load_LEdata(fdict2sub,casename,timesel=timeselc, rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
            subp2 = le.load_LEdata(fdict2sub,casename,timesel=timeselp, rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
            sub2 = cutl.seasonalize_monthlyts(subp2.T,season=sea2).T - cutl.seasonalize_monthlyts(subc2.T,season=sea2).T
            lefld2 = lesea2.mean(axis=1) - sub2.mean(axis=1)
    else:
        lefld2=lepsea2.mean(axis=1)-lecsea2.mean(axis=1)

    lemm, lebb, lerval, lepval, lestd_err = sp.stats.linregress(lefld1,lefld2)

    if conditional:
        lecdatcnd = le.load_LEdata(fdictcnd,casename,timesel=timeselc, 
                                   rettype='ndarray',conv=leconvcnd,ftype=ftype,local=local)
        (numenscnd,ntimecnd) = lecdatcnd.shape
        lepdatcnd=le.load_LEdata(fdictcnd,casename,timesel=timeselp, 
                                 rettype='ndarray',conv=leconvcnd,ftype=ftype,local=local)
        lecseacnd = cutl.seasonalize_monthlyts(lecdatcnd.T,season=seacnd).T
        lepseacnd = cutl.seasonalize_monthlyts(lepdatcnd.T,season=seacnd).T
        leseacnd = lepseacnd - lecseacnd
        lefldcnd=lepseacnd.mean(axis=1)-lecseacnd.mean(axis=1)

    if addnat:
        # historicalNat
        casename2='historicalNat'
        lecdat1n = le.load_LEdata(fdict1,casename2,timesel=timeselc, 
                                  rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
        (numens1,ntime1) = lecdat1n.shape
        lepdat1n=le.load_LEdata(fdict1,casename2,timesel=timeselp, 
                                rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
        lecsea1n = cutl.seasonalize_monthlyts(lecdat1n.T,season=sea1).T
        lepsea1n = cutl.seasonalize_monthlyts(lepdat1n.T,season=sea1).T
        lesea1n = lepsea1n - lecsea1n
        if performop1:
            fdict1op = {'field': field1+region1op, 'ncfield': ncfield1, 'comp': comp1}
            subc1n = le.load_LEdata(fdict1op,casename2,timesel=timeselc, 
                                    rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
            subp1n = le.load_LEdata(fdict1op,casename2,timesel=timeselp, 
                                    rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
            sub1n = cutl.seasonalize_monthlyts(subp1n.T,season=sea1).T -\
                    cutl.seasonalize_monthlyts(subc1n.T,season=sea1).T

            if op1=='sub': # subtract
                lefld1n = lesea1n.mean(axis=1) - sub1n.mean(axis=1)
            elif op1=='div': # divide
                lefld1n = lesea1n.mean(axis=1) / sub1n.mean(axis=1)
        else:
            lefld1n=lepsea1n.mean(axis=1)-lecsea1n.mean(axis=1)

        lecdat2n = le.load_LEdata(fdict2,casename2,timesel=timeselc, 
                                  rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
        (numens2,ntime2) = lecdat1n.shape
        lepdat2n=le.load_LEdata(fdict2,casename2,timesel=timeselp, 
                                rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
        lecsea2n = cutl.seasonalize_monthlyts(lecdat2n.T,season=sea2).T
        lepsea2n = cutl.seasonalize_monthlyts(lepdat2n.T,season=sea2).T
        lesea2n = lepsea2n - lecsea2n
        if performop2:
            if op2=='sub': # subtract
                fdict2sub = {'field': field2+region2op, 'ncfield': ncfield2, 'comp': comp2}
                subc2n = le.load_LEdata(fdict2sub,casename2,timesel=timeselc, 
                                        rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
                subp2n = le.load_LEdata(fdict2sub,casename2,timesel=timeselp, 
                                        rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
                sub2n = cutl.seasonalize_monthlyts(subp2n.T,season=sea2).T -\
                        cutl.seasonalize_monthlyts(subc2n.T,season=sea2).T
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
            gissatc= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselc,seas=sea1) 
            gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselp,seas=sea1)
            obsreg1 =  cutl.calc_regmean(gissatp-gissatc,latgis,longis,region1)
            if performop1:
                opfld1 = cutl.calc_regmean(gissatp-gissatc,latgis,longis,region1op)

                if op1=='sub': # subtract
                    obsreg1 = obsreg1.mean() - opfld1.mean()
                elif op1=='div': # divide
                    obsreg1 = obsreg1.mean() / opfld1.mean()
            else:
                obsreg1=obsreg1.mean()
        elif field1 == 'sic':
            nsidcfile = '/HOME/rkm/work/BCs/NSIDC/td_bootstrap_197811_latest_128_64_sicn_1978111600-2013121612.nc'
            latns = cnc.getNCvar(nsidcfile,'lat')
            lonns = cnc.getNCvar(nsidcfile,'lon')
            nssicc= cnc.getNCvar(nsidcfile,'SICN',timesel=timeselc,seas=sea1)*100
            nssicp= cnc.getNCvar(nsidcfile,'SICN',timesel=timeselp,seas=sea1)*100
            obsreg1 =  cutl.calc_regmean(nssicp-nssicc,latns,lonns,region1)
            if performop1:
                print '@@ no performop1 for field1=sic in addobs'
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


    if addsims:
        # 120-yr averages
        sims=('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5')
        # why did this work before? had to change to Series (jun 11 2015)
        flddf1 = pd.Series(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea1, 
                                           meantype='time',region=region1))*simconv1

        if performop1:
            flddf1op = pd.Series(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea1, 
                                                 meantype='time',region=region1op))*simconv1

            if op1=='sub': # subtract
                if region1op=='nh':
                    flddf1 = flddf1 - flddf1op.mean()
                else:
                    flddf1 = flddf1 - flddf1op
            elif op1=='div': # divide
                flddf1 = flddf1 / flddf1op
        

        flddf2 = pd.Series(lmd.loaddata((simfield2,),sims,ncfields=(simncfield2,), timefreq=sea2, 
                                           meantype='time',region=region2))*simconv2

        if performop2:
            flddf2op = pd.Series(lmd.loaddata((simfield2,),sims,ncfields=(simncfield2,), timefreq=sea2, 
                                                 meantype='time',region=region2op))*simconv2

            if op2=='sub': # subtract
                if region2op=='nh':
                    flddf2 = flddf2 - flddf2op.mean()
                else:
                    flddf2 = flddf2 - flddf2op
            elif op2=='div': # divide
                flddf2 = flddf2 / flddf2op

        simmm, simbb, simrval, simpval, simstd_err = sp.stats.linregress(np.squeeze(flddf1.values),np.squeeze(flddf2.values))
        print '-------- SIMS slope, rval, pval: ' + str(simmm),str(simrval),str(simpval)

    if shorttermsims:
        simsss=('E1','E2','E3','E4','E5'); simsstr='sbEonly' # sub smp
        #simsss=('R1','R2','R3','R4','R5'); simsstr='sbRonly'
        fldssdf1 = pd.DataFrame(lmd.loaddata((simfield1,),simsss,ncfields=(simncfield1,), timefreq=sea1, 
                                             region=region1))*simconv1 # "simss" subsample sims

        if performop1:
            fldssdf1op = pd.Series(lmd.loaddata((simfield1,),simsss,ncfields=(simncfield1,), timefreq=sea1, 
                                                 meantype='time',region=region1op))*simconv1

            if op1=='sub': # subtract
                if region1op=='nh':
                    fldssdf1 = fldssdf1 - fldssdf1op.mean() # want to subtract ens mean nh
                else:
                    fldssdf1 = fldssdf1 - fldssdf1op # not sure if want ens mean here...
            elif op1=='div': # divide
                fldssdf1 = fldssdf1 / fldssdf1op
        
        fldssdf1 = subsamp_sims(fldssdf1,numyrs=11)

        fldssdf2 = pd.DataFrame(lmd.loaddata((simfield2,),simsss,ncfields=(simncfield2,), timefreq=sea2, 
                                        region=region2))*simconv2

        if performop2:
            fldssdf2op = pd.Series(lmd.loaddata((simfield2,),simsss,ncfields=(simncfield2,), timefreq=sea2, 
                                                 meantype='time',region=region2op))*simconv2

            if op2=='sub': # subtract
                if region2op=='nh':
                    fldssdf2 = fldssdf2 - fldssdf2op.mean() # want to subtract ens mean nh
                else:
                    fldssdf2 = fldssdf2 - fldssdf2op # not sure if want ens mean here...
            elif op2=='div': # divide
                fldssdf2 = fldssdf2 / fldssdf2op
        
        fldssdf2 = subsamp_sims(fldssdf2,numyrs=11)

        simssmm, simssbb, simssrval, simsspval, simssstd_err = sp.stats.linregress(fldssdf1,
                                                                                   fldssdf2)
        print '-------- (subsamp) SIMS slope, rval, pval: ' + str(simssmm),str(simssrval),str(simsspval)


        

    # Example of conditional scatter: 
    # http://stackoverflow.com/questions/12965075/matplotlib-scatter-plot-colour-as-function-of-third-variable
    """x = np.linspace(0, 20, 100)
    y = np.sin(x)
    z = x + 20 * y

    scaled_z = (z - z.min()) / z.ptp()
    colors = plt.cm.coolwarm(scaled_z)

    plt.scatter(x, y, marker='+', edgecolors=colors, s=150, linewidths=4)
    plt.show()"""

    # ========== FIGURE ===========

    fig,ax=plt.subplots(1,1)
    fig.set_size_inches(8,5)
    if conditional:
        cmap=plt.cm.afmhot #autumn # Greens_r #OrRd_r # coolwarm_r
        #markerdt={'marker': 'o', 'linewidth':0}
        scaled_cnd = (lefldcnd - lefldcnd.min()) / lefldcnd.ptp()
        cndcolors = plt.cm.afmhot(scaled_cnd,alpha=0.7) #edgecolors=cndcolors

        b2r20cm = plt.cm.get_cmap('blue2red_20')
        #testcm = plt.cm.autumn.from_list('testcm',plt.cm.autumn(scaled_cnd,alpha=0.7),N=4)
        testclrs = np.flipud(b2r20cm.colors[[7,11,14,-1],:])
        #testcm = b2r20cm.from_list('testcm',testclrs,N=4)
        testcm = matplotlib.colors.ListedColormap(testclrs, name='testcm')

        ledatlg=mlines.Line2D([],[],color=cndcolors[len(cndcolors)/2+5],
                              markerfacecolor=cndcolors[len(cndcolors)/2+5],marker='o',
                              linestyle='none',mew=0,markersize=10)

        ledat=ax.scatter(lefld1,lefld2,c=lefldcnd,cmap=testcm,vmin=cmincnd,
                         vmax=cmaxcnd,s=90,alpha=0.7,marker='o', linewidth=0)
        ledat.set_clim([cmincnd,cmaxcnd])

        cbar_ax = fig.add_axes([.91,.25, .02,.5]) 
        cb=fig.colorbar(ledat,cax=cbar_ax)#,orientation='horizontal')
        cb.set_label('$\Delta$ BKS SIC (%)')
        lcolor='brown'
        leghnds=(ledatlg,)
    else:
        ledat=ax.scatter(lefld1,lefld2,color=ccm.get_linecolor('darkolivegreen3'),
                          marker='*',s=100,alpha=0.5)
        lcolor=ccm.get_linecolor('darkolivegreen3')
        leghnds=(ledat,)

    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    onex=np.linspace(axxlims[0],axxlims[1])
    ax.plot(onex,lemm*onex + lebb, color=lcolor,linewidth=1)
    
    obstr=''
    
    legstrs=(casename + ' LE',)

    print '-------- LE slope, rval, pval: ' + str(lemm),str(lerval),str(lepval)
    if addnat:
        ledatn=ax.scatter(lefld1n,lefld2n,color=ccm.get_linecolor('steelblue3'),marker='*',s=5**2,alpha=0.5)
        axylims = ax.get_ylim()
        axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,lemmn*onex + lebbn, color=ccm.get_linecolor('steelblue3'),linewidth=1)
        leghnds=leghnds + (ledatn,)
        legstrs=legstrs + (casename2 + ' LE',)
        print '-------- LENat slope, rval, pval: ' + str(lemmn),str(lervaln),str(lepvaln)
    if addobs:
        obs=ax.scatter(obsreg1,obsreg2,color='blue',marker='s',s=60)
        leghnds=leghnds + (obs,)
        legstrs=legstrs + ('Observations',)
        obstr='obs'

    if addsims:
        simsh = ax.scatter(flddf1,flddf2,color='0.3',marker='o',s=120,alpha=0.7)
        #axylims = ax.get_ylim()
        #axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,simmm*onex + simbb, color='0.3',linewidth=1)
        leghnds=leghnds + (simsh,)
        legstrs=legstrs + ('Modelled SIC forcing',)
        obstr=obstr+'sims'

    if shorttermsims:
        simssh = ax.scatter(fldssdf1,fldssdf2,color='0.7',marker='o',s=60,alpha=0.7)
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,simssmm*onex + simssbb, color='0.7',linewidth=1)
        leghnds=leghnds + (simssh,)
        legstrs=legstrs + ('11-yr Modelled SIC forcing',)
        obstr=obstr+'sims'

    fontP = fm.FontProperties()
    fontP.set_size('small')
    ax.legend(leghnds,legstrs,
               loc='best',fancybox=True,framealpha=0.5,prop=fontP)#,frameon=False) 
    ax.set_ylim((ymin,ymax))
    ax.annotate(casename + ' R= ' + '$%.2f$'%(lerval) + ', p='+ '$%.2f$'%(lepval),
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
        ax.set_xlabel(sea1 + ' ' + field1 + ' ' + region1 + ' ' + opstr1)
    else:
        ax.set_xlabel(xlab)
    if performop2:
        opstr2= op2 + '_' + region2op

    if ylab==None:
        ax.set_ylabel(sea2 + ' ' + field2 + ' ' + region2 + ' ' + opstr2)
    else:
        ax.set_ylabel(ylab)

    axylims = ax.get_ylim()
    if axylims[0]<=0 and axylims[1]>=0:
        ax.axhline(y=0,color='k',linewidth=.5,linestyle='--')
    axxlims = ax.get_xlim()
    if axxlims[0]<=0 and axxlims[1]>=0:
        ax.axvline(x=0,color='k',linewidth=.5,linestyle='--')

    if printtofile:
        if conditional:
            cnd='cnd'
        else: cnd=''
        if addnat:
            plt.savefig('scatterregress_' + field1 + region1 + opstr1 + str(sea1) + '_v_' + 
                        field2 + region2 +opstr2 + str(sea2) + '_' + casename + '_' + 
                        casename2 + '_2002-12_1979-89' + obstr + cnd + '.pdf')
        else:
            plt.savefig('scatterregress_' + field1 + region1 + opstr1 + str(sea1) + '_v_' + 
                        field2 + region2 +opstr2 + str(sea2) + '_' + casename + 
                        '_2002-12_1979-89' + obstr + cnd + '.pdf')

if dohist:
    conv=leconv1 # just assume we are doing variable 1
    sea=sea1

    histcdat=le.load_LEdata(fdict2,'historical',timesel=timeselc, rettype='ndarray',conv=conv,ftype=ftype,local=local)
    (numens,ntime) = histcdat.shape
    histpdat=le.load_LEdata(fdict2,'historical',timesel=timeselp, rettype='ndarray',conv=conv,ftype=ftype,local=local)

    # Now have 11 years of monthly data. Grab DJF:
    histc = cutl.seasonalize_monthlyts(histcdat.T,season=sea1).T
    histp = cutl.seasonalize_monthlyts(histpdat.T,season=sea1).T

    histnatcdat=le.load_LEdata(fdict2,'historicalNat',timesel=timeselc, rettype='ndarray',conv=conv,ftype=ftype,local=local)
    #(numens,ntime,nlatlon) = histnatcdat.shape
    histnatpdat=le.load_LEdata(fdict2,'historicalNat',timesel=timeselp, rettype='ndarray',conv=conv,ftype=ftype,local=local)
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
    addfld=True # add contours of second field
    ensmean=True

    seasp=sea1 # season of spatial field
    sear=sea2 # season of regional avgs


    # spatial field1 in color
    leconvsp=1
    #fieldsp='zg50000.00'; ncfieldsp='zg'; compsp='Amon'; 
    fieldsp='tas'; ncfieldsp='tas'; compsp='Amon'; 

    # spatial field2 in contours
    leconvsp2=1
    fieldsp2='zg50000.00'; ncfieldsp2='zg'; compsp2='Amon'
    cminsp2=-10; cmaxsp2=10

    # regional avg field 1
    leconvr=-1 # this way, sea ice loss is linked with positive changes elsewhere
    #fieldr='sic'; ncfieldr='sic'; compr='OImon'; regionr='bksmori'
    fieldr='tas'; ncfieldr='tas'; compr='Amon'; regionr='eurasiamori'; leconvr=-1 # not sure want opp?

    # regional avg field 2
    leconvr2=1
    #fieldr2='tas'; ncfieldr2='tas'; compr2='Amon'; regionr2='eurasiamori'; leconvr2=-1 # so cooling=high heights
    fieldr2='zg50000.00'; ncfieldr2='zg'; compr2='Amon'; regionr2='bksmori'
    

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
    fdictsp2 = {'field': fieldsp2, 'ncfield': ncfieldsp2, 'comp': compsp2}
    fdictr = {'field': fieldr+regionr, 'ncfield': ncfieldr, 'comp': compr}
    fdictr2 = {'field': fieldr2+regionr2, 'ncfield': ncfieldr2, 'comp': compr2}

    casename = 'historical'

    lat=le.get_lat(local=local)
    lon=le.get_lon(local=local)
    nlat=len(lat); nlon=len(lon)

    
    if ensmean: # then do a regression in time, over the ensemble mean
         # LOAD SPATIAL DATA (contours)
        # need to instead get the full time 1979-2012 
        lealldatsp = le.load_LEdata(fdictsp,casename,timesel=timeselall, rettype='ndarray',conv=leconvsp,ftype='ensmean',local=local)
        (numensp,ntimesp,nspacesp) = lealldatsp.shape
        # time needs to be first dimension
        lealldatsp = np.transpose(lealldatsp,(1,0,2))
        leallseasp = cutl.seasonalize_monthlyts(lealldatsp,season=seasp).mean(axis=1) # average the ensemble 
        leseasp = leallseasp # time x space.flat

        # need to instead get the full time 1979-2012 
        lealldatsp2 = le.load_LEdata(fdictsp2,casename,timesel=timeselall, rettype='ndarray',conv=leconvsp2,ftype='ensmean',local=local)
        (numensp,ntimesp,nspacesp) = lealldatsp2.shape
        # time needs to be first dimension
        lealldatsp2 = np.transpose(lealldatsp2,(1,0,2))
        leallseasp2 = cutl.seasonalize_monthlyts(lealldatsp2,season=seasp).mean(axis=1) # average the ensemble 
        leseasp2 = leallseasp2 # time x space.flat

        # LOAD 1D DATA
        lealldatr = le.load_LEdata(fdictr,casename,timesel=timeselall, rettype='ndarray',conv=leconvr,ftype='ensmean',local=local)
        (numenr,ntimer) = lealldatr.shape
        leallsear = cutl.seasonalize_monthlyts(lealldatr.T,season=sear).mean(axis=1)
        lesear = leallsear # ntime
        print '@@@@@@@@@@@@@@@@@@@@@@@@ ' + str(lesear.shape)
        if donorm:
            lesear=lesear / lesear.std()

    else:
        # LOAD SPATIAL DATA (contours)
        lecdatsp = le.load_LEdata(fdictsp,casename,timesel=timeselc, rettype='ndarray',conv=leconvsp,ftype=ftype,local=local)
        (numensp,ntimesp,nspacesp) = lecdatsp.shape
        lepdatsp=le.load_LEdata(fdictsp,casename,timesel=timeselp, rettype='ndarray',conv=leconvsp,ftype=ftype,local=local)
        # time needs to be first dimension
        lecdatsp = np.transpose(lecdatsp,(1,0,2))
        lepdatsp = np.transpose(lepdatsp,(1,0,2))
   
        lecseasp = cutl.seasonalize_monthlyts(lecdatsp,season=seasp).mean(axis=0) # average the 11 seasonal avgs 
        lepseasp = cutl.seasonalize_monthlyts(lepdatsp,season=seasp).mean(axis=0) # average the 11 seasonal avgs 
        leseasp = lepseasp - lecseasp # numens x space.flat

        # LOAD SPATIAL DATA 2 (colors)
        lecdatsp2 = le.load_LEdata(fdictsp2,casename,timesel=timeselc, rettype='ndarray',conv=leconvsp2,ftype=ftype,local=local)
        (numensp,ntimesp,nspacesp) = lecdatsp2.shape
        lepdatsp2=le.load_LEdata(fdictsp2,casename,timesel=timeselp, rettype='ndarray',conv=leconvsp2,ftype=ftype,local=local)
        # time needs to be first dimension
        lecdatsp2 = np.transpose(lecdatsp2,(1,0,2))
        lepdatsp2 = np.transpose(lepdatsp2,(1,0,2))
        lecseasp2 = cutl.seasonalize_monthlyts(lecdatsp2,season=seasp).mean(axis=0) # average the 11 seasonal avgs 
        lepseasp2 = cutl.seasonalize_monthlyts(lepdatsp2,season=seasp).mean(axis=0) # average the 11 seasonal avgs 
        leseasp2 = lepseasp2 - lecseasp2 # numens x space.flat

        # LOAD 1D DATA
        lecdatr = le.load_LEdata(fdictr,casename,timesel=timeselc, rettype='ndarray',conv=leconvr,ftype=ftype,local=local)
        (numenr,ntimer) = lecdatr.shape
        lepdatr=le.load_LEdata(fdictr,casename,timesel=timeselp, rettype='ndarray',conv=leconvr,ftype=ftype,local=local)
        lecsear = cutl.seasonalize_monthlyts(lecdatr.T,season=sear).mean(axis=0) # this is right b/c of the .T
        lepsear = cutl.seasonalize_monthlyts(lepdatr.T,season=sear).mean(axis=0)
        lesear = lepsear - lecsear # numens
        if donorm:
            lesear=lesear / lesear.std()
    


    #mm, bb, rval, pval, std_err = sp.stats.linregress(lesear,leseasp) # errors. why? @@
    slope,intercept = np.polyfit(lesear,leseasp,1)
    bkssat=slope.reshape((nlat,nlon)) # SAT regress on SIC
    slope,intercept = np.polyfit(lesear,leseasp2,1)
    bkszg=slope.reshape((nlat,nlon)) # Z500 regress on SIC

    if ensmean:
        # LOAD 1D DATA
        lealldatr2 = le.load_LEdata(fdictr2,casename,timesel=timeselall, rettype='ndarray',conv=leconvr2,ftype='ensmean',local=local)
        (numenr2,ntimer2) = lealldatr2.shape
        leallsear2 = cutl.seasonalize_monthlyts(lealldatr2.T,season=sear).mean(axis=1)
        lesear2 = leallsear2 # ntime
        print '@@@@@@@@@@@@@@@@@@@@@@@@ ' + str(lesear2.shape)
        if donorm:
            lesear2=lesear2 / lesear2.std()

        print '@@@@ figure out regression in time....check units etc'
    else:
        # LOAD 1D DATA (2)
        lecdatr2 = le.load_LEdata(fdictr2,casename,timesel=timeselc, rettype='ndarray',conv=leconvr2,ftype=ftype,local=local)
        (numenr2,ntimer2) = lecdatr2.shape
        lepdatr2=le.load_LEdata(fdictr2,casename,timesel=timeselp, rettype='ndarray',conv=leconvr2,ftype=ftype,local=local)
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
                title=seasp + ' ' + fieldsp + ' regressed onto ' + sear + ' ' +fieldr+regionr )

    ax=axs[1] #
    cplt.kemmap(eurzg,lat,lon,type='nh',axis=ax, cmin=cmin2,cmax=cmax2,
                title= seasp + ' ' + fieldsp + ' regressed onto ' + sear + ' ' +fieldr2+regionr2)

    if printtofile:
        fig.savefig(fieldr +regionr + '_' + fieldsp + seasp + \
                    + '_regresson_' + fieldr2 + regionr2 + sear + normstr + '.pdf') 

    # @@ create a figure with regression contours on top of other regression:
    #   e.g. z500 regressed onto BKS SIC contours on SAT regressed onto BKS SIC map

    printtofile=False
    lons, lats = np.meshgrid(lon,lat)
    cmlen=15.
    incr = (cmaxsp2-cminsp2) / (cmlen)
    conts = np.arange(cminsp2,cmaxsp2+incr,incr)


    fig,ax=plt.subplots(1,1)
    fig.set_size_inches(5,5)
    bm,pc=cplt.kemmap(bkssat,lat,lon,type='nh',axis=ax,cmin=cmin,cmax=cmax,
                title=seasp + ' regressions onto ' + sear + ' ' + fieldr+regionr)
    bm.contour(lons,lats,bkszg,levels=conts,
               colors='k',linewidths=1,latlon=True)
    if printtofile:
        fig.savefig(fieldsp + seasp + '_' + fieldsp2 + seasp \
                    + '_regresson_' + fieldr + regionr + sear + normstr + '.pdf') 


    tmp=np.zeros(bkssat.shape)
    fig,ax=plt.subplots(1,1)
    fig.set_size_inches(5,5)
    bm,pc=cplt.kemmap(tmp,lat,lon,type='nh',axis=ax,cmin=cmin,cmax=cmax,
                      title=seasp + ' regressions onto ' + sear + ' ' + fieldr+regionr,suppcb=True)
    
    bm.contour(lons,lats,bkszg,levels=conts,
               colors='k',linewidths=1,latlon=True)
    if printtofile:
        fig.savefig(fieldsp2 + seasp \
                    + '_regresson_' + fieldr + regionr + sear + normstr + '.pdf') 


def sample120yravg(lesea,numsamp,nummems=11,allowreps=True):
    """ choose nummems random ens members, and do it numsamp
        times.

        nummems: number of members to choose. default 11 (for 11mem x 11yr = 122yr)
        allowreps: default True. allow repeat ensemble indices to be chosen in different
                   numsamp chunks (never dupes within each nummems selection)

        returns ltavg, ltsigma (average & std of each nummem set. len numsamp)
    """

    # sample 11 ensemble members: 11 members x 11 years = ~120 yrs to equal AGCM sims
    # (do this numsamp times with diff combos of 11).
    savesel=np.zeros((numsamp,nummems))
    ltavg=np.zeros((numsamp)) # 'lt' = longterm avg
    ltsigma=np.zeros((numsamp))

    lelist=list(lesea)

    for ii in np.arange(0,numsamp):
        
        #sel = np.random.randint(0,len(lesea),size=nummems) # out of the number of ens members
        
        sel = random.sample(np.arange(0,len(lelist)),nummems) # no duplicates!
        print len(lelist) 
        print sel
        savesel[ii]=sel

        #c=[ a[i] for i in b] # how to select multiple elements of list w/ indices
        tmp = np.array([lelist[ind] for ind in sel]) # select members
        if allowreps==False:
            # now delete those members
            sel.sort(reverse=True) # if remove the indices starting from end, shouldn't fail.
            for ind in sel: 
                print 'deleting index ' + str(ind)
                del lelist[ind] 

        ltavg[ii] = tmp.mean()
        ltsigma[ii]= tmp.std()

    return ltavg,ltsigma


if dolongtermavg:

    longtermLE=False # else, subsample AGCM sims

    printtofile=False

    numsamp=100 # how many times to sample 11 ens members

    addraw=False # add the decadal diffs?
    subnh=False # @@@@ when subtracting, prob should subtract the mean NH temp, not individual runs
    substr='' # for filename
    ttlstr=' ' 
    prstr=''
    option2=True # otherwise, original method of selecting 11 members numsamp times. (longtermLE=True)
    if option2:
        numsamp=50

    sea='DJF'
    siglevel=0.1

    leconv=1
    field='tas'
    ncfield='tas'
    comp='Amon'
    region='eurasiamori'
    
    simconv1=1
    if field=='tas': simfield1='st'; simncfield1='ST'
    elif field=='zg50000.00': simfield1='gz50000'; simncfield1='PHI'; simconv1=1/con.get_g()
    elif field=='sia': simfield1='sicn'; simncfield1='SICN'; print '@@ danger, sia actually sicn average'
    else: print 'cannot addsims for ' + field;

    fdict = {'field': field+region, 'ncfield': ncfield, 'comp': comp}

    casename = 'historical'

    mew=1.5; ms=7
    firebrick=ccm.get_linecolor('firebrick')
    hcol=ccm.get_linecolor('darkolivegreen3')
    hcolline=ccm.get_linecolor('darkolivegreen3')#'darkseagreen4')
    ltcol=ccm.get_linecolor('darkseagreen4')
    natcol=ccm.get_linecolor('steelblue3')
    natcolline=ccm.get_linecolor('steelblue3')#4')

    lat=le.get_lat(local=local)
    lon=le.get_lon(local=local)
    nlat=len(lat); nlon=len(lon)

    # LOAD 1D DATA
    lecdat = le.load_LEdata(fdict,casename,timesel=timeselc, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
    (numen,ntime) = lecdat.shape
    lepdat=le.load_LEdata(fdict,casename,timesel=timeselp, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
    lecsea = cutl.seasonalize_monthlyts(lecdat.T,season=sea).mean(axis=0)
    lepsea = cutl.seasonalize_monthlyts(lepdat.T,season=sea).mean(axis=0)
    lesea = lepsea - lecsea # numens

    if subnh:
        fdictsub = {'field': field+'nh', 'ncfield': ncfield, 'comp': comp}
        subc = le.load_LEdata(fdictsub,casename,timesel=timeselc, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
        subp = le.load_LEdata(fdictsub,casename,timesel=timeselp, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
        sub = cutl.seasonalize_monthlyts(subp.T,season=sea).T - cutl.seasonalize_monthlyts(subc.T,season=sea).T
        # @@@ change this to subtract the MEAN hemispheric anom from all ens members
        #lesea = lesea - sub.mean(axis=1)
        lesea = lesea - sub.mean(axis=1).mean()

    # RAW anomalies (decadal diffs)
    # calc the pdf associated with the hist
    rawpdf=norm.fit(lesea)
    rawmean=rawpdf[0]
    rawsd=rawpdf[1]
    #Generate X points
    rawxlims = [-4*rawsd+rawmean, 4*rawsd+rawmean] # large limits
    rawxx = np.linspace(rawxlims[0],rawxlims[1],500)
    #Get Y points via Normal PDF with fitted parameters
    rawpdf_fitted = norm.pdf(rawxx,loc=rawmean,scale=rawsd)


    if longtermLE:
        # second option: choose 4 random, non-overlapping groups of 11 members x 11-yrs
        #                do this numsamp times. Compare with sigma of 5 AGCM sims
        if option2:
            ltavg2=np.zeros(numsamp)
            ltsigma2=np.zeros(numsamp)
            for ns in np.arange(0,numsamp):
                # returns 4 avgs and sigmas across e/ of the 12 selected members
                avgtmp,sigmajunk = sample120yravg(lesea,4,nummems=11,allowreps=False)
                ltavg2[ns]=avgtmp.mean()
                ltsigma2[ns]=avgtmp.std()

            ltavg=ltavg2
            ltsigma=ltsigma2
            prstr='opt2'        
        else: # I don't think this makes as much sense for estimating sigma convergence? @@

            # sample 11 ensemble members: 11 members x 11 years = ~120 yrs to equal AGCM sims
            # (do this numsamp times with diff combos of 11).        
            ltavg,ltsigma = sample120yravg(lesea,numsamp) # option one.

            # the sigma returned here isn't really appropriate for comparison with the AGCM sims.
            # This is because it is the sigma across each 12 member selection
            # instead, now should randomly select 3-5 elements from ltavg and compute 
            #   sigma across those small sets. For non-repeating selections of 4, 
            #   can only run 12 times if original numsamp is 50.
            nummems2=4
            ltavg2,ltsigma2= sample120yravg(ltavg,numsamp/nummems2,nummems=nummems2,allowreps=False)

        # Now calc the pdf associated with the hist
        ltpdf=norm.fit(ltavg)
        ltmean=ltpdf[0]
        ltsd=ltpdf[1]
        #Generate X points
        ltxlims = [-4*ltsd+ltmean, 4*ltsd+ltmean] # large limits
        ltxx = np.linspace(ltxlims[0],ltxlims[1],500)
        #Get Y points via Normal PDF with fitted parameters
        ltpdf_fitted = norm.pdf(ltxx,loc=ltmean,scale=ltsd)


    else: # not longtermLE (short term sims instead)

        # here we want to subsample 11-year segments from the sims
        sims=('E1','E2','E3','E4','E5'); simsstr='sbEonly' # sub
        #sims=('R1','R2','R3','R4','R5'); simsstr='sbRonly'
        simflddf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                             region=region))*simconv1

        

        if subnh:
            simsubdf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                                 filetype='diff',region='nh'))*simconv1

            # want to subtract the mean hemispheric avg anomaly
            simflddf = simflddf - simsubdf.mean(axis=1).mean() # average over sims and then time (should be a scalar)


        subsamp = subsamp_sims(simflddf,numyrs=11)
        plotsims = subsamp

        # Now calc the pdf associated with the hist
        sspdf=norm.fit(plotsims) # ss = subsamp
        ssmean=sspdf[0]
        sssd=sspdf[1]
        #Generate X points
        ssxlims = [-4*sssd+ssmean, 4*sssd+ssmean] # large limits
        ssxx = np.linspace(ssxlims[0],ssxlims[1],500)
        #Get Y points via Normal PDF with fitted parameters
        sspdf_fitted = norm.pdf(ssxx,loc=ssmean,scale=sssd)

        print '===== AGCM mean ' + str(ssmean) + ' sigma ' + str(sssd)


    # ========= add sims (120-yr averages)
    sims=('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5','NSIDC'); #simsstr=''
    #simsRN=('R1','R2','R3','R4','R5','NSIDC')
    simsR=('R1','R2','R3','R4','R5')
    #sims=simsRN; simsstr='Ronly'
    simflddf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                         meantype='time',region=region),index=sims)*simconv1


    # --- for estimating sigma in Individual SIC forcing ensemble only:
    simflddfr = pd.DataFrame(lmd.loaddata((simfield1,),simsR,ncfields=(simncfield1,), timefreq=sea, 
                                          meantype='time',region=region),index=simsR)*simconv1
    simval=simflddfr.values[0,:]
    simstds = [simval[1:].std(), simval[[0,2,3,4]].std(), simval[[0, 1, 3, 4]].std(), 
               simval[[0, 1, 2, 4]].std(), simval[0:-1].std()] # hack

    simsigma=simval.std()


    if subnh:
        simsubcdf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                             filetype='ctl',region='nh'))*simconv1
        simsubpdf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                             filetype='pert',region='nh'))*simconv1
        simfldcdf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                       filetype='ctl',region=region))*simconv1
        simfldpdf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                       filetype='pert',region=region))*simconv1
        simfldcdf=simfldcdf-simsubcdf
        simfldpdf=simfldpdf-simsubpdf
        (tstat,pvals)=cutl.ttest_ind(simfldpdf.values,simfldcdf.values,axis=0) 
        print pvals.shape

        simflddf = simfldpdf.mean(axis=0)-simfldcdf.mean(axis=0)
        simfldm = ma.masked_where(pvals>siglevel,simflddf.values)
        # have to calc pvals 

    else:
        simpvdf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                            filetype='pval',region=region),index=sims)
        simfldm = ma.masked_where(simpvdf.values>siglevel,simflddf.values)


    # add obs
    if field == 'tas':
        gisfile = '/HOME/rkm/work/DATA/GISS/gistemp1200_ERSST.nc'
        latgis=cnc.getNCvar(gisfile,'lat')
        longis=cnc.getNCvar(gisfile,'lon')
        gissatc= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselc,seas=sea) 
        gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselp,seas=sea)
        obsregc=cutl.calc_regmean(gissatc,latgis,longis,region)
        obsregp=cutl.calc_regmean(gissatp,latgis,longis,region)
        (tstat,obspv) = cutl.ttest_ind(obsregp,obsregc)
        obsreg = obsregp.mean()-obsregc.mean()

        if subnh:
            obssubc=cutl.calc_regmean(gissatc,latgis,longis,'nh')
            obssubp=cutl.calc_regmean(gissatp,latgis,longis,'nh')
            (tstat,obspv) = cutl.ttest_ind(obsregp-obssubp,obsregc-obssubc)
            obsreg = (obsregp-obssubp).mean() - (obsregc-obssubc).mean()

            substr='_subnh'
            ttlstr='-NH '



    # ========= make the figures =========
    if longtermLE:
        if option2:
            # how does sigma converge with n?
            fig,ax=plt.subplots(1,1)
            ax.hist(ltsigma,color=ltcol,alpha=0.5)
            ax.axvline(simsigma,color='0.3',linewidth=3)
            for ll in np.arange(0,5):
                ax.axvline(simstds[ll],color='0.3',linewidth=1)
            ax.set_title('$\sigma$ across sets of 4 Estimated 120-year Eurasian' + ttlstr + 'SAT change (DJF)')
            ax.set_ylabel('Density')
            ltlg=mpatches.Patch(color=ltcol,alpha=0.5)
            simlg=mlines.Line2D([],[],color='0.3',linewidth=3) # all five
            simslg=mlines.Line2D([],[],color='0.3',linewidth=1)
            ax.set_xlabel('n=' + str(numsamp))
            ax.legend((ltlg,simlg,simslg),('120-yr avg Historical','5 SIC forcings', 'combos of 4 SIC forcings'), 
                      loc='upper left',frameon=False)
            if printtofile:
                fig.savefig(field + region + substr + 'SIGMA_' + sea + '_LEsims' +\
                            simsstr + 'obs_est120yravg_hist' + prstr + '.pdf')
        else:
            fig,ax=plt.subplots(1,1)
            ax.hist(ltsigma,color=ltcol,alpha=0.5)
            ax.axvline(simsigma,color='0.3',linewidth=3)
            ax.set_title('$\sigma$ across sets of 4 Estimated 120-year Eurasian' + ttlstr + 'SAT change (DJF)')
            ax.set_ylabel('Density')
            ltlg=mpatches.Patch(color=ltcol,alpha=0.5)
            simlg=mlines.Line2D([],[],color='0.3',linewidth=3)
            ax.set_xlabel('n=' + str(numsamp))
            ax.legend((ltlg,simlg),('120-yr avg Historical','SIC forcings'), 
                      loc='upper right',frameon=False)
            if printtofile:
                fig.savefig(field + region + substr + 'SIGMA_' + sea + '_LEsims' +\
                            simsstr + 'obs_est120yravg_hist' + prstr + '.pdf')


        ltlg=mlines.Line2D([],[],color=ltcol,linewidth=2)
        simlg=mlines.Line2D([],[],color='0.3',linestyle='none',marker='o')
        nsidclg=mlines.Line2D([],[],color='g',linestyle='none',marker='o')
        obslg=mlines.Line2D([],[],color='b',linestyle='none',marker='s')

        fig,ax=plt.subplots(1,1)
        ax.hist(ltavg,normed=True,color=ltcol,alpha=0.5)
        ax.plot(ltxx,ltpdf_fitted,color=ltcol,linewidth=2)
        ax.plot(simflddf.values, np.zeros(len(simflddf)),color='0.3',
                 marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms)
        ax.plot(simfldm,np.zeros(len(simflddf)),color='0.3',
                 marker='o',linestyle='',fillstyle='full')
        if subnh: # HACK @@
            ax.plot(simflddf['NSIDC'],0,color='g',
                    marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms,mec='g')
        else:
            ax.plot(simflddf['NSIDC'].values[0],0,color='g',
                    marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms,mec='g')
        ax.plot(obsreg,0,color='b',marker='s',fillstyle='full',markersize=ms)# blue for actual obs. green for simulated obs
        ylims=ax.get_ylim()
        ax.set_ylim(-0.1,ylims[1])
        ax.set_yticklabels('')
        ax.axhline(y=0,color='k')
        ax.axvline(x=0,color='k',linestyle='--')
        ax.set_ylabel('Density')
        ax.set_title('Estimated 120-year Eurasian' + ttlstr + 'SAT change (DJF)')
        ax.set_xlabel('$\Delta$ SAT ($^\circ$C); n=' + str(numsamp))
        ax.legend((ltlg,simlg,nsidclg,obslg),('120-yr avg Historical','SIC forcings','NSIDC SIC forcing','Observations'),
                  loc='upper left',frameon=False)


        if printtofile:
            fig.savefig(field + region + substr + '_' + sea + '_LEsims' + simsstr + 'obs_est120yravg_hist' + prstr + '.pdf')

        if addraw:
            # add RAW for comparison
            ax.hist(lesea,normed=True,color=hcol,alpha=0.5)
            ax.plot(rawxx,rawpdf_fitted,color=hcolline,linewidth=2)
            rawlg=mlines.Line2D([],[],color=hcolline,linewidth=2)
            ax.legend((ltlg,rawlg,simlg,obslg),
                      ('120-yr avg Historical','Historical','SIC forcings','Observations'),
                      loc='upper left',frameon=False)
            if printtofile:
                fig.savefig(field + region + substr + '_' + sea + '_LEsims' + simsstr + 'obs_est120yravgraw_hist' + prstr + '.pdf')

    else: # not longtermLE

        sslg=mlines.Line2D([],[],color=ltcol,linewidth=2) # subsamp sims
        simlg=mlines.Line2D([],[],color='0.3',linestyle='none',marker='o')
        nsidclg=mlines.Line2D([],[],color='g',linestyle='none',marker='o')
        obslg=mlines.Line2D([],[],color='b',linestyle='none',marker='s')

        fig,ax=plt.subplots(1,1)
        #ax.hist(ltavg,normed=True,color=ltcol,alpha=0.5)
        #ax.plot(ltxx,ltpdf_fitted,color=ltcol,linewidth=2)
        ax.hist(plotsims,normed=True,color=ltcol,alpha=0.5)
        ax.plot(ssxx,sspdf_fitted,color=ltcol,linewidth=2)
        ax.plot(simflddf.values, np.zeros(len(simflddf)),color='0.3',
                 marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms)
        ax.plot(simfldm,np.zeros(len(simflddf)),color='0.3',
                 marker='o',linestyle='',fillstyle='full')
        if subnh: # HACK @@
            ax.plot(simflddf['NSIDC'],0,color='g',
                    marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms,mec='g')
        else:
            ax.plot(simflddf['NSIDC'].values[0],0,color='g',
                    marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms,mec='g')
        ax.plot(obsreg,0,color='b',marker='s',fillstyle='full',markersize=ms)# blue for actual obs. green for simulated obs
        ylims=ax.get_ylim()
        ax.set_ylim(-0.1,ylims[1])
        ax.set_yticklabels('')
        ax.axhline(y=0,color='k')
        ax.axvline(x=0,color='k',linestyle='--')
        ax.set_ylabel('Density')
        ax.set_title('11-year epochs Eurasian' + ttlstr + 'SAT change (DJF)')
        ax.set_xlabel('$\Delta$ SAT ($^\circ$C); LE avg= $%.2f$'%(rawmean) + ' AGCM avg= $%.2f$'%(ssmean) )

        # add RAW 
        ax.hist(lesea,normed=True,color=hcol,alpha=0.5)
        ax.plot(rawxx,rawpdf_fitted,color=hcolline,linewidth=2)
        rawlg=mlines.Line2D([],[],color=hcolline,linewidth=2)
        ax.legend((rawlg,sslg,simlg,obslg),
                  ('Historical LE','11-yr subsample SIC forcings','120-yr average SIC forcings','Observations'),
                  loc='upper left',frameon=False)
        if printtofile:
            now = str(datetime.datetime.now().time())
            fig.savefig(field + region + substr + '_' + sea + '_LEsims' + simsstr +\
                        'obs_11yrsubsmpavgraw_hist' + prstr + now + '.pdf')



    # ==== estimate sigma uncertainty ====
    # Need sia too
    leconvsia=1
    fieldsia='sia'
    ncfieldsia='sia'
    compsia='OImon'
    regionsia='nh'

    fdictsia = {'field': fieldsia+regionsia, 'ncfield': ncfieldsia+regionsia, 'comp': compsia}

    lecsia = le.load_LEdata(fdictsia,casename,timesel=timeselc, rettype='ndarray',conv=leconvsia,ftype=ftype,local=local)
    (numensi,ntimesi) = lecsia.shape
    lepsia=le.load_LEdata(fdictsia,casename,timesel=timeselp, rettype='ndarray',conv=leconvsia,ftype=ftype,local=local)
    lecsisea = cutl.seasonalize_monthlyts(lecsia.T,season=sea).mean(axis=0)
    lepsisea = cutl.seasonalize_monthlyts(lepsia.T,season=sea).mean(axis=0)
    lesisea = lepsisea - lecsisea # numens

    lesistd = lesisea.std() # compare to std of option 1 histogram?

    # add sims
    simsR=('R1','R2','R3','R4','R5')
    # make sia: @@@@@@@@@@@@@@@@@@@@@
    """# for estimating sigma:
    simflddfr = pd.DataFrame(lmd.loaddata((simfield1,),simsR,ncfields=(simncfield1,), timefreq=sea, 
                                          meantype='time',region=region),index=simsR)*simconv1
    simval=simflddfr.values[0,:]
    simstds = [simval[1:].std(), simval[[0,2,3,4]].std(), simval[[0, 1, 3, 4]].std(), 
               simval[[0, 1, 2, 4]].std(), simval[0:-1].std()]

    simsigma=simval.std()"""
