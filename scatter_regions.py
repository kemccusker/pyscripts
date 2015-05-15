"""
  scatter_regions.py
      This script makes scatter plots of regional means. Can do two
      variables, two regions, two seasons. There are figures of just the simulation mean
      values with confidence intervals, and figures with all regional averages
      in time with regression lines. 10/10/2014

      1/26/15: try to make into a function that uses datablobs.
      
"""

import cccmautils as cutl
import constants as con
import cccmaplots as cplt
import cccmacmaps as ccm
import pandas as pd
import numpy.ma as ma
import matplotlib.font_manager as fm

printtofile=False
plt.close('all')

conv1=1; conv2=1
graveraint= 9.80665 # m/s2 (ERA-Int, different from Canadian models)

plotscatter=True
addobs=True
addle=False

# field1 is x
#field1='st'; ncfield1='ST'
#field1='sia'; ncfield1='SICN'
field1='gz50000'; ncfield1='PHI'; conv1=1/con.get_g()
#field1='pmsl'; ncfield1='PMSL'
#field1='net'; 

region1='bksmori' #'polcap65'
sea1='DJF' #'ND' #'DJF'

# field2 is y
field2='st'; ncfield2='ST'
#field2='sia'; ncfield2='SICN'
#field2='pmsl'; ncfield2='PMSL'
#field2='gz50000'; ncfield2='PHI'; conv2=1/con.get_g()
region2= 'eurasiamori' #'eurasiamori'
sea2='DJF' #'ND' #'DJF'

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
if (sea1 in ('DJF','NDJ')) or (sea2 in ('DJF','NDJ')):
    ntime=ntime-1
#if (sea1 in ('DJF','NDJ') and sea2 not in ('DJF','NDJ')) or (sea2 in ('DJF','NDJ') and sea1 not in ('DJF','NDJ')):
#    ntime=ntime-1
    
nsim=len(TOT)
#if sea2 in ('DJF','NDJ') and sea1 not in ('DJF','NDJ'): # have to shorten other timeseries
#    inittime=ntime*nsim-1
#else:
inittime=ntime*nsim

timesel='0002-01-01,0121-12-31'
alltotcr1=np.zeros(inittime)
alltotcr2=np.zeros(inittime)
alltotpr1=np.zeros(inittime)
alltotpr2=np.zeros(inittime)
alltotr1=np.zeros(inittime)
alltotr2=np.zeros(inittime)

allantcr1=np.zeros(inittime)
allantcr2=np.zeros(inittime)
allantpr1=np.zeros(inittime)
allantpr2=np.zeros(inittime)
allantr1=np.zeros(inittime)
allantr2=np.zeros(inittime)
tallii=0 # index to keep track of time in accumulated TOT ensemble
aallii=0 # index to keep track of time in accumulated ANT ensemble
for sim in sims:

    

    if field1 in ('turb','net'):
        
        fielda='hfl'; fieldb='hfs'
        fnamec,fnamepa=con.build_filepathpair(sim,fielda)
        fnamecb,fnamepb=con.build_filepathpair(sim,fieldb)
        fldc = cnc.getNCvar(fnamec,fielda.upper(),timesel=timesel,seas=sea1) +\
               cnc.getNCvar(fnamecb,fieldb.upper(),timesel=timesel,seas=sea1)
        fldp = cnc.getNCvar(fnamepa,fielda.upper(),timesel=timesel,seas=sea1) +\
               cnc.getNCvar(fnamepb,fieldb.upper(),timesel=timesel,seas=sea1)

        if field1=='net':
            fieldb='flg'
            conv=-1
            fnamecb,fnamepb=con.build_filepathpair(sim,fieldb)

            fldc = fldc + cnc.getNCvar(fnamecb,fieldb.upper(),
                                       timesel=timesel,seas=sea1)*conv
            fldp = fldp + cnc.getNCvar(fnamepb,fieldb.upper(),
                                       timesel=timesel,seas=sea1)*conv
            conv=1
    else:
        fnamec,fnamep=con.build_filepathpair(sim,field1)
        fldc=cnc.getNCvar(fnamec,ncfield1,timesel=timesel,seas=sea1)*conv1
        fldp=cnc.getNCvar(fnamep,ncfield1,timesel=timesel,seas=sea1)*conv1


    lat=cnc.getNCvar(fnamec,'lat')
    lon=cnc.getNCvar(fnamec,'lon')

    if field1=='sia':
        fldcreg=cutl.calc_regtotseaicearea(fldc,lat,lon,region1) # isarea=False
        fldpreg=cutl.calc_regtotseaicearea(fldp,lat,lon,region1) # isarea=False
    else:
        # IF FLUX, MASK OUT OPEN-WATER:
        #print 'Implement masking out of open water in control for FLUXES @@@'
        if field1 in ('hfl','hfs','flg','fsg','turb','net'):
            print 'field1, ' + field1 + ' is a flux. Mask out non-ice in control' # @@
            # mask out regions that are not ice in the control (as P&M 2014 JClim)
            sicfname,junk=con.build_filepathpair(sim,'sicn')
            sicnc = cnc.getNCvar(sicfname,'SICN',timesel=timesel,seas=sea1)

            fldc = ma.masked_where(sicnc<.10,fldc) 
            fldp = ma.masked_where(sicnc<.10,fldp)
            
        fldcreg=cutl.calc_regmean(fldc,lat,lon,region1)
        fldpreg=cutl.calc_regmean(fldp,lat,lon,region1)
        
    flddregdt[sim] = np.mean(fldpreg-fldcreg,axis=0) # time mean
    if sea2 in ('DJF','NDJ') and sea1 not in ('DJF','NDJ'): # have to shorten other timeseries
        flddregtsdt[sim] = fldpreg[:-1,...]-np.mean(fldcreg[:-1,...],axis=0) # anomaly timeseries from ctl mean
    else:
        flddregtsdt[sim] = fldpreg-np.mean(fldcreg,axis=0)

    if (field2 != field1) or (sea2 != sea1):
        print field2 + ' != ' + field1 + ' or ' + str(sea2) + ' != ' + str(sea1)
        
        if field2 in ('turb','net'):

            fielda='hfl'; fieldb='hfs'
            fnamec,fnamepa=con.build_filepathpair(sim,fielda)
            fnamecb,fnamepb=con.build_filepathpair(sim,fieldb)
            fldc2 = cnc.getNCvar(fnamec,fielda.upper(),timesel=timesel,seas=sea2) +\
                   cnc.getNCvar(fnamecb,fieldb.upper(),timesel=timesel,seas=sea2)
            fldp2 = cnc.getNCvar(fnamepa,fielda.upper(),timesel=timesel,seas=sea2) +\
                   cnc.getNCvar(fnamepb,fieldb.upper(),timesel=timesel,seas=sea2)

            if field2=='net':
                fieldb='flg'
                conv=-1
                fnamecb,fnamepb=con.build_filepathpair(sim,fieldb)

                fldc2 = fldc2 + cnc.getNCvar(fnamecb,fieldb.upper(),
                                           timesel=timesel,seas=sea2)*conv
                fldp2 = fldp2 + cnc.getNCvar(fnamepb,fieldb.upper(),
                                           timesel=timesel,seas=sea2)*conv
                conv=1
        else:
            fnamec,fnamep=con.build_filepathpair(sim,field2)
            fldc2=cnc.getNCvar(fnamec,ncfield2,timesel=timesel,seas=sea2)*conv2
            fldp2=cnc.getNCvar(fnamep,ncfield2,timesel=timesel,seas=sea2)*conv2
            
    ## elif sea2 != sea1:
    ##     fldc2=cnc.getNCvar(fnamec,ncfield2,timesel=timesel,seas=sea2)*conv2
    ##     fldp2=cnc.getNCvar(fnamep,ncfield2,timesel=timesel,seas=sea2)*conv2
    else:
        fldc2=fldc
        fldp2=fldp

    if field2=='sia':
        fldcreg2=cutl.calc_regtotseaicearea(fldc,lat,lon,region2) # isarea=False
        fldpreg2=cutl.calc_regtotseaicearea(fldp,lat,lon,region2) # isarea=False
    else:
        if field2 in ('hfl','hfs','flg','fsg','turb','net'):
            # mask out regions that are not ice in the control (as P&M 2014 JClim)
            print 'field2, ' + field2 + ' is a flux. Mask out non-ice in control' # @@
            sicfname,junk=con.build_filepathpair(sim,'sicn')
            sicnc2 = cnc.getNCvar(sicfname,'SICN',timesel=timesel,seas=sea2)

            fldc2 = ma.masked_where(sicnc2<.10,fldc2)
            fldp2 = ma.masked_where(sicnc2<.10,fldp2)
            
        fldcreg2=cutl.calc_regmean(fldc2,lat,lon,region2)
        fldpreg2=cutl.calc_regmean(fldp2,lat,lon,region2)
        
    flddregdt2[sim] = np.mean(fldpreg2-fldcreg2,axis=0)
#    if sea2 in ('DJF','NDJ') and sea1 not in ('DJF','NDJ'): # have to shorten other timeseries
#        flddregtsdt2[sim] = fldpreg2[:-1,...]-np.mean(fldcreg2[:-1,...],axis=0)
#    else:
    flddregtsdt2[sim] = fldpreg2-np.mean(fldcreg2,axis=0) # anomaly timeseries from ctl mean

    # @@@ there are files with concated ensembles now: e.g. sim name kemctl1enscat, diffname ENSCAT
    if sim in TOT:
        #print 'concat tot ens members'
        alltotcr1[tallii:tallii+ntime] = fldcreg[:ntime] # control reg1
        alltotcr2[tallii:tallii+ntime] = fldcreg2[:ntime] # control reg2
        alltotpr1[tallii:tallii+ntime] = fldpreg[:ntime] # pert reg1
        alltotpr2[tallii:tallii+ntime] = fldpreg2[:ntime] # pert reg2
        alltotr1[tallii:tallii+ntime] = flddregtsdt[sim] # diff reg1
        alltotr2[tallii:tallii+ntime] = flddregtsdt2[sim] # diff reg2
        tallii+=ntime
        #totser.append(pd.Series(flddregtsdt))
        #totser2.append(pd.Series(flddregtsdt2))
    elif sim in ANT:
        #print 'concat ant ens members'
        allantcr1[aallii:aallii+ntime] = fldcreg[:ntime] # to deal w/ SON vs DJF for example
        allantcr2[aallii:aallii+ntime] = fldcreg2[:ntime]
        allantpr1[aallii:aallii+ntime] = fldpreg[:ntime]
        allantpr2[aallii:aallii+ntime] = fldpreg2[:ntime]
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
    
    fontP = fm.FontProperties()
    fontP.set_size('small')

    # SCATTER PLOT OF SIM MEANS
    # prepare data for regression line
    #mm, bb, rval, pval, std_err = sp.stats.linregress(dfts[sim].values[0],dfts[sim].values[1])
    #  a way to subset: df.ix[:,'R1':'R5'] -- but it can't get both R and E data.
    #rdf=df.ix[:,'R1':'R5']
    #edf=df.ix[:,'E1':'E5']
    #regressdf = rdf.append(edf) # didn't append like I expected. inserted NaNs in places
    

    # this will work:
    allens=['R1','R2','R3','R4','R5','E1','E2','E3','E4','E5'] # has to be an array!
    rens=['R1','R2','R3','R4','R5']
    eens=['E1','E2','E3','E4','E5'] 
    sub = df.loc[:,allens]
    subx=sub.ix[0,:]
    suby=sub.ix[1,:]
    mm, bb, rval, pval, std_err = sp.stats.linregress(subx,suby)

    #sube = df.loc[0,['E1','E2','E3','E4','E5']]
    firebrick=ccm.get_linecolor('firebrick')

    printtofile=False
    fig,ax = plt.subplots(1)
    fig.set_size_inches(6,5)
    rr = plt.scatter(df.filter(regex='R').values[0],df.filter(regex='R').values[1],
                     color='0.3',marker='o',s=8**2,alpha=0.7)
    ee = plt.scatter(sub.ix[0,5:],sub.ix[1,5:],
                     color=firebrick,marker='o',s=8**2,alpha=0.7)
    #plt.scatter(df['HAD'].values[0],df['HAD'].values[1],color=cd['HAD'],marker='s',s=8**2)
    #ns=plt.scatter(df['NSIDC'].values[0],df['NSIDC'].values[1],color=cd['NSIDC'],marker='o',s=8**2,alpha=0.7) # @@@@ add back?
    ns=plt.scatter(df['NSIDC'].values[0],df['NSIDC'].values[1],color='green',marker='o',s=8**2,alpha=0.7)
    #plt.scatter(df['ENS'].values[0],df['ENS'].values[1],color=cd['ENS'],marker='s',s=10**2)
    #plt.scatter(df['ENSE'].values[0],df['ENSE'].values[1],color=cd['ENSE'],marker='o',s=10**2)

    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    onex=np.linspace(axxlims[0],axxlims[1])
    plt.plot(onex,mm*onex + bb, color='0.5',linewidth=1)#,linestyle='--')
    rsq = rval**2
    print 'PAPER --- mean sims regr info: SLOPE: ' + str(mm) + ', RVAL: ' + str(rval) + ', PVAL: ' + str(pval) + ', R2: ' + str(rsq)
    

    if sea1 in ('ND','DJF') and sea2 in ('ND','DJF') and region1=='bksmori' and region2 in ('eurasia','eurasiathin','eurasiamori') and field1=='gz50000' and field2=='st':
        # PAPER FIG!
        if sea1 == sea2 == 'DJF':
            monlab = 'Dec-Jan-Feb'
        else:
            monlab = 'Nov-Dec'

        xlab = '$\Delta$ ' + monlab + ' Barents-Kara Seas'
        if field1=='pmsl':
            xlab = xlab + ' SLP (hPa)'
        else:
            xlab = xlab + ' Z500 (m)'

        ylab = '$\Delta$ ' + monlab + ' Eurasian SAT ($^\circ$C)'
        
        # ### ADD ACTUAL OBS HERE:
        # GISS (SAT) and ERA-INT (Z500)
        erafile = '/HOME/rkm/work/DATA/ERAINT/td_era_int_197901_201412_gp_128_64_phi500_1979011612-2014121612.nc'
        gisfile = '/HOME/rkm/work/DATA/GISS/gistemp1200_ERSST.nc'
        eraz500c= cnc.getNCvar(erafile,ncfield1,timesel='1979-01-01,1989-12-31',seas=sea1)/graveraint
        eraz500p= cnc.getNCvar(erafile,ncfield1,timesel='2002-01-01,2012-12-31',seas=sea1)/graveraint
        erareg = cutl.calc_regmean(eraz500p-eraz500c,lat,lon,region1)
        gissatc= cnc.getNCvar(gisfile,'tempanomaly',timesel='1979-01-01,1989-12-31',seas=sea2) # @@ need to deal w/ scale factor? OH IT HAPPENS IN MY FUNCTION
        gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel='2002-01-01,2012-12-31',seas=sea2)
        latgis=cnc.getNCvar(gisfile,'lat')
        longis=cnc.getNCvar(gisfile,'lon')
        gisreg =  cutl.calc_regmean(gissatp-gissatc,latgis,longis,region2)

        if addobs:
            obs=plt.scatter(erareg.mean(),gisreg.mean(),color='blue',marker='s',s=8**2)
        
    else:
        addobs=False
        print 'Cannot add obs for this field combo' #@@
        xlab = str(sea1) + ' ' + field1 + ' ' + region1
        ylab = str(sea2) + ' ' + field2 + ' ' + region2

    if addobs:
        probs='obs' # print obs string
        plt.legend((rr,ee,ns,obs),('Individual SIC forcings',
                                   'Average SIC forcing','NSIDC SIC forcing','Observations'),
                   loc='best',fancybox=True,framealpha=0.5,prop=fontP)#,frameon=False)  
    else:
        probs=''
        plt.legend((rr,ee,ns),('Individual SIC forcings',
                                   'Average SIC forcing','NSIDC SIC forcing'),
                   loc='best',fancybox=True,framealpha=0.5,prop=fontP)#,frameon=False) 
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
                    field2 + region2 + str(sea2) + 'wacepap' + probs + '.pdf')
        fig.savefig('scatterregress_' + field1 + region1 + str(sea1) + '_v_' +
                    field2 + region2 + str(sea2) + 'wacepap' + probs + '.eps')

    # ### ADD LE DATA #####
    if addle:
        printtofile=True
        import loadLE as le

        timeselc='1979-01-01,1989-12-31'
        timeselp='2002-01-01,2012-12-31'

        lefield1='zg50000.00'; lencfield1='zg'; comp1='Amon'; leconv1=1 #region1='bksmori'
        lefield2='tas'; lencfield2='tas'; comp2='Amon'; #region2='eurasiamori'
        fdict1 = {'field': lefield1+region1, 'ncfield': lencfield1, 'comp': comp1}
        fdict2 = {'field': lefield2+region2, 'ncfield': lencfield2, 'comp': comp2}


        ftype='fullts' # 'fullclimo' or 'climo' or 'fullts'
        lefdict1 = {'field': lefield1+region1, 'ncfield': lencfield1, 'comp': comp1}
        lefdict2 = {'field': lefield2+region2, 'ncfield': lencfield2, 'comp': comp2}

        # historical
        lecdat1 = le.load_LEdata(fdict1,'historical',timesel=timeselc, rettype='ndarray',conv=leconv1,ftype=ftype)
        (numens1,ntime1) = lecdat1.shape
        lepdat1=le.load_LEdata(fdict1,'historical',timesel=timeselp, rettype='ndarray',conv=leconv1,ftype=ftype)
        lecsea1 = cutl.seasonalize_monthlyts(lecdat1.T,season=sea1).T
        lepsea1 = cutl.seasonalize_monthlyts(lepdat1.T,season=sea1).T
        lefld1=lepsea1.mean(axis=1)-lecsea1.mean(axis=1)

        lecdat2 = le.load_LEdata(fdict2,'historical',timesel=timeselc, rettype='ndarray',conv=conv2,ftype=ftype)
        (numens2,ntime2) = lecdat1.shape
        lepdat2=le.load_LEdata(fdict2,'historical',timesel=timeselp, rettype='ndarray',conv=conv2,ftype=ftype)
        lecsea2 = cutl.seasonalize_monthlyts(lecdat2.T,season=sea2).T
        lepsea2 = cutl.seasonalize_monthlyts(lepdat2.T,season=sea2).T
        lefld2=lepsea2.mean(axis=1)-lecsea2.mean(axis=1)

        obs=plt.scatter(erareg.mean(),gisreg.mean(),color='blue',marker='s',s=8**2)

        ledat=plt.scatter(lefld1,lefld2,color=ccm.get_linecolor('darkolivegreen3'),marker='*',s=5**2,alpha=0.5)
        lemm, lebb, lerval, lepval, lestd_err = sp.stats.linregress(lefld1,lefld2)
        axylims = ax.get_ylim()
        axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,lemm*onex + lebb, color=ccm.get_linecolor('darkolivegreen3'),linewidth=1)

        # historicalNat
        lecdat1n = le.load_LEdata(fdict1,'historicalNat',timesel=timeselc, rettype='ndarray',conv=leconv1,ftype=ftype)
        (numens1,ntime1) = lecdat1n.shape
        lepdat1n=le.load_LEdata(fdict1,'historicalNat',timesel=timeselp, rettype='ndarray',conv=leconv1,ftype=ftype)
        lecsea1n = cutl.seasonalize_monthlyts(lecdat1n.T,season=sea1).T
        lepsea1n = cutl.seasonalize_monthlyts(lepdat1n.T,season=sea1).T
        lefld1n=lepsea1n.mean(axis=1)-lecsea1n.mean(axis=1)

        lecdat2n = le.load_LEdata(fdict2,'historicalNat',timesel=timeselc, rettype='ndarray',conv=conv2,ftype=ftype)
        (numens2,ntime2) = lecdat1n.shape
        lepdat2n=le.load_LEdata(fdict2,'historicalNat',timesel=timeselp, rettype='ndarray',conv=conv2,ftype=ftype)
        lecsea2n = cutl.seasonalize_monthlyts(lecdat2n.T,season=sea2).T
        lepsea2n = cutl.seasonalize_monthlyts(lepdat2n.T,season=sea2).T
        lefld2n=lepsea2n.mean(axis=1)-lecsea2n.mean(axis=1)

        ledatn=plt.scatter(lefld1n,lefld2n,color=ccm.get_linecolor('steelblue3'),marker='*',s=5**2,alpha=0.5)
        lemmn, lebbn, lervaln, lepvaln, lestd_errn = sp.stats.linregress(lefld1n,lefld2n)
        axylims = ax.get_ylim()
        axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,lemmn*onex + lebbn, color=ccm.get_linecolor('steelblue3'),linewidth=1)

        plt.legend((rr,ee,ns,obs,ledat,ledatn),('Individual SIC forcings','Average SIC forcing',
                                                'NSIDC SIC forcing','Observations',
                                                'historical LE', 'historicalNat LE'),
                   loc='best',fancybox=True,framealpha=0.5,prop=fontP)#,frameon=False)    

        if printtofile:
            fig.savefig('scatterregress_' + field1 + region1 + str(sea1) + '_v_' +
                        field2 + region2 + str(sea2) + 'wacepapobs_LE.pdf')
            #fig.savefig('scatterregress_' + field1 + region1 + str(sea1) + '_v_' +
            #            field2 + region2 + str(sea2) + 'wacepapobs_LE.eps')


    # ##### JUST SIMS ONE COLOR #########
    fig,ax = plt.subplots(1)
    fig.set_size_inches(6,5)
    rr = plt.scatter(df.filter(regex='R').values[0],df.filter(regex='R').values[1],
                     color='0.3',marker='o',s=8**2,alpha=0.7)
    ee = plt.scatter(sub.ix[0,5:],sub.ix[1,5:],
                     color='0.3',marker='o',s=8**2,alpha=0.7)
    #plt.scatter(df['HAD'].values[0],df['HAD'].values[1],color=cd['HAD'],marker='s',s=8**2)
    #ns=plt.scatter(df['NSIDC'].values[0],df['NSIDC'].values[1],color=cd['NSIDC'],marker='o',s=8**2,alpha=0.7) # @@@@ add back?
    ns=plt.scatter(df['NSIDC'].values[0],df['NSIDC'].values[1],color='green',marker='o',s=8**2,alpha=0.7)
    #plt.scatter(df['ENS'].values[0],df['ENS'].values[1],color=cd['ENS'],marker='s',s=10**2)
    #plt.scatter(df['ENSE'].values[0],df['ENSE'].values[1],color=cd['ENSE'],marker='o',s=10**2)

    axylims = ax.get_ylim()
    axxlims = ax.get_xlim()
    onex=np.linspace(axxlims[0],axxlims[1])
    plt.plot(onex,mm*onex + bb, color='0.5',linewidth=1)#,linestyle='--')
    if addobs:
        obs=plt.scatter(erareg.mean(),gisreg.mean(),color='blue',marker='s',s=8**2)

    if addobs:
        probs='obs' # print obs string
        plt.legend((ee,ns,obs),('Modelled SIC forcings',
                                'NSIDC SIC forcing','Observations'),
                   loc='best',fancybox=True,framealpha=0.5,prop=fontP)#,frameon=False)  
    else:
        probs=''
        plt.legend((ee,ns),('Modelled SIC forcings','NSIDC SIC forcing'),
                   loc='best',fancybox=True,framealpha=0.5,prop=fontP)#,frameon=False) 
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
                    field2 + region2 + str(sea2) + 'wacepap' + probs + '_onecol.pdf')
        fig.savefig('scatterregress_' + field1 + region1 + str(sea1) + '_v_' +
                    field2 + region2 + str(sea2) + 'wacepap' + probs + '_onecol.eps')

    if addle:
        # ### JUST LE DATA #####
        fig,ax = plt.subplots(1)
        fig.set_size_inches(6,5)

        obs=ax.scatter(erareg.mean(),gisreg.mean(),color='blue',marker='s',s=8**2)

        ledat=ax.scatter(lefld1,lefld2,color=ccm.get_linecolor('darkolivegreen3'),marker='*',s=5**2,alpha=0.5)
        lemm, lebb, lerval, lepval, lestd_err = sp.stats.linregress(lefld1,lefld2)
        axylims = ax.get_ylim()
        axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,lemm*onex + lebb, color=ccm.get_linecolor('darkolivegreen3'),linewidth=1)

        ledatn=ax.scatter(lefld1n,lefld2n,color=ccm.get_linecolor('steelblue3'),marker='*',s=5**2,alpha=0.5)
        lemmn, lebbn, lervaln, lepvaln, lestd_errn = sp.stats.linregress(lefld1n,lefld2n)
        axylims = ax.get_ylim()
        axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,lemmn*onex + lebbn, color=ccm.get_linecolor('steelblue3'),linewidth=1)

        plt.legend((obs,ledat,ledatn),('Observations',
                                       'historical LE', 'historicalNat LE'),
                   loc='best',fancybox=True,framealpha=0.5,prop=fontP)#,frameon=False)    

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
                        field2 + region2 + str(sea2) + '_obs_LE.pdf')

        # ### JUST LE HISTORICAL DATA #####
        fig,ax = plt.subplots(1)
        fig.set_size_inches(6,5)

        obs=plt.scatter(erareg.mean(),gisreg.mean(),color='blue',marker='s',s=8**2)

        ledat=plt.scatter(lefld1,lefld2,color=ccm.get_linecolor('darkolivegreen3'),marker='*',s=5**2,alpha=0.5)
        lemm, lebb, lerval, lepval, lestd_err = sp.stats.linregress(lefld1,lefld2)
        axylims = ax.get_ylim()
        axxlims = ax.get_xlim()
        onex=np.linspace(axxlims[0],axxlims[1])
        ax.plot(onex,lemm*onex + lebb, color=ccm.get_linecolor('darkolivegreen3'),linewidth=1)

        plt.legend((obs,ledat),('Observations','historical LE'),
                   loc='best',fancybox=True,framealpha=0.5,prop=fontP)#,frameon=False)    

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
                        field2 + region2 + str(sea2) + '_obs_historicalLE.pdf')


    printtofile=False

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
        print 'TOT: R: ' + str(rval) + ', R squared: ' + str(rsq) + ', SLOPE: ' + str(mm) + ', PVAL: ' + str(pval)
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
        print 'ANT: R: ' + str(rval) + ', R squared: ' + str(rsq) + ', SLOPE: ' + str(mm) + ', PVAL: ' + str(pval)

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
print 'ANT: ' + field2, region2
cutl.calc_pvals(allantpr2,allantcr2)
print 'TOT: ' + field2, region2
cutl.calc_pvals(alltotpr2,alltotcr2)

plot_anttotsbplt_histpdf((allantcr1,allantpr1),(alltotcr1,alltotpr1),
                         label=(field1 + region1 + '_' + str(sea1)),printtofile=printtofile)
print 'ANT: ' + field1, region1
cutl.calc_pvals(allantpr1,allantcr1)
print 'TOT: ' + field1, region1
cutl.calc_pvals(alltotpr1,alltotcr1)
