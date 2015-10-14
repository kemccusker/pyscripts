""" canesm_LE_general.py

    This script basically reads in files that are already processed
    into area-averages and such.

"""
import matplotlib.pyplot as plt 
import numpy as np
import scipy as sp
import pandas as pd
import cccmaNC as cnc
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
import matplotlib.colors as mcolors
from scipy.stats import norm
import random
import numpy.ma as ma
import scipy.io as sio
import loadCanESM2data as lcd
import corrstats as corrstats

# exception handling works below. add to other clauses @@
cutl=reload(cutl)

local=True

conditional=True # plot scatter conditional on 3rd var

styearsR = [ 8.,  7.,  2.,  8.,  8.] # variable SIC styears
styearsE=[ 4.,  1.,  7.,  3.,  1.]; # mean SIC styears
styearsN=[1.] # NSIDC sim
# composite when for 2: '14:51:28.762886'
# composite when for ss and ns: 17:01:16.908687
styearPI = 2
# PI anomyrs  when: '14:51:28.762886sic'. mean eursat=-0.0026
anomyearsPI =[[74,  8],
       [57,  5],
       [50, 53],
       [51, 42],
       [ 7, 81],
       [15,  9],
       [61, 57],
       [48, 21],
       [49, 24],
       [65, 24],
       [33, 66],
       [71, 36],
       [29, 43],
       [48, 68],
       [18, 67],
       [ 5, 21],
       [17,  1],
       [40, 16],
       [79, 18],
       [10, 82],
       [ 0, 19],
       [70, 36],
       [82, 67],
       [14,  7],
       [57, 13],
       [18, 34],
       [59, 29],
       [20,  6],
       [12,  6],
       [58, 68],
       [82, 49],
       [66,  2],
       [ 0, 10],
       [28, 59],
       [63, 73],
       [28, 83],
       [51, 54],
       [19,  2],
       [26, 47],
       [52, 66],
       [27, 10],
       [37, 48],
       [38, 64],
       [33, 36],
       [42, 85],
       [28, 24],
       [79, 51],
       [80, 10],
       [44, 57],
       [32, 41]]

"""styearPI = 0 # PI styear
# PI anomyrs  when: '14:51:28.762886'
anomyearsPI = [[79, 26],
        [58, 17],
        [79, 41],
        [28, 24],
        [37, 80],
        [25, 48],
        [30, 49],
        [33, 85],
        [51, 43],
        [82,  6],
        [62, 34],
        [ 3, 17],
        [56, 63],
        [67,  4],
        [29, 73],
        [74,  0],
        [28,  8],
        [46, 56],
        [14, 76],
        [72, 37],
        [88,  4],
        [31, 56],
        [31,  4],
        [40,  0],
        [20, 49],
        [21, 81],
        [68, 56],
        [77, 37],
        [18,  1],
        [82, 26],
        [78, 55],
        [22, 47],
        [78, 16],
        [54, 76],
        [ 5, 47],
        [58, 25],
        [49, 20],
        [64, 36],
        [34,  2],
        [62, 18],
        [89,  3],
        [25, 10],
        [10, 30],
        [84, 55],
        [88, 18],
        [87, 50],
        [68, 26],
        [ 5, 67],
        [77, 13],
        [74, 61]]"""



performop1 = True
#op1='div'; region1op='gm' # polar amp: gt60n / gm
#op1='sub'; region1op='deeptrop' # pole-eq temp gradient: gt60n - deeptrop (or trop)
#op1='sub'; region1op='less' # bks sea ice minus less (laptev/east siberian)
#op1='sub'; region1op='eurasiamori' # z500 grad b/w bks and eurasia
#op1='sub'; region1op='eurasiagrad2' # z500 grad b/w bks and eurasia
#op1='sub'; region1op='eurasiamorie' # z500 grad b/w bks and eurasia
op1='sub'; region1op='eurasiasemed' # z500 grad b/w bks and eurasia

performop2 = False
op2='div'; region2op='gm' # polar amp: gt60n / gm
#op2='sub'; region2op='nh'
#op2='sub'; region2op='eurasiamori' # sub z500

timeselc='1979-01-01,1989-12-31'
timeselp='2002-01-01,2012-12-31'
timeselall = '1979-01-01,2012-12-31'

# x field
#field1='turb'; ncfield1='turb'; comp1='Amon'; region1='bksmori'
field1='zg50000.00'; ncfield1='zg'; comp1='Amon'; region1='bksmori'
#field1='sia'; ncfield1='sianh'; comp1='OImon'; region1='nh'
#field1='sic'; ncfield1='sic'; comp1='OImon'; region1='bksmori'; #'gt60n' #region1='bksmori' # @@ a hack. prefer SIA
#field1='tas'; ncfield1='tas'; comp1='Amon'; region1='eurasiamori' #region1='bksmori'
leconv1= 1 
sea1='DJF'

# y field
field2='tas'; ncfield2='tas'; comp2='Amon'; region2='eurasiamori' #'eurasiathicke'; 
#field2='tas'; ncfield2='tas'; comp2='Amon'; region2='gt60n'
#field2='zg50000.00'; ncfield2='zg'; comp2='Amon'; region2='bksmori'
#field2='turb'; ncfield2='turb'; comp2='Amon'; region2='bksmori'
leconv2=1
sea2='DJF'

fieldcnd = 'sic'; ncfieldcnd='sic'; compcnd='OImon'; regioncnd='bksmori'
leconvcnd=1
seacnd=sea1
cmincnd=-12; cmaxcnd=3 # for 5 color colorbar
cmapcnd='blue2red_20'

#fieldcnd = 'tas'; ncfieldcnd='tas'; compcnd='Amon'; regioncnd='gm' #regioncnd='eurasiamori'
#leconvcnd=1
#seacnd=sea1
#if regioncnd=='eurasiamori':
#    cmincnd=-2; cmaxcnd=3 # for 5 color colorbar 
#elif regioncnd=='gm':
#    cmincnd=0.5; cmaxcnd=1
#cmapcnd='blue2red_20'

cisiglevel=0.05
siglevel=0.05
simsR=('R1','R2','R3','R4','R5')
simsE=('E1','E2','E3','E4','E5')

xlab=ylab=None
if field2=='tas' and region2=='eurasiamori' and performop2==True and op2=='sub':
    ylab = '$\Delta$ DJF SAT(Eurasia) - SAT(NH) ($^\circ$C)'
elif field2=='tas' and region2=='eurasiamori' and performop2==False:
    ylab = '$\Delta$ DJF Eurasian SAT ($^\circ$C)'
if field1=='zg50000.00' and region1=='bksmori' and performop1==False:
    xlab = '$\Delta$ DJF Barents-Kara Seas Z500 (m)'


simconv1=simconv2=1
if field1=='tas': simfield1='st'; simncfield1='ST'
elif field1=='zg50000.00': simfield1='gz50000'; simncfield1='PHI'; simconv1=1/con.get_g()
elif field1=='sia': simfield1='sicn'; simncfield1='SICN'; print '@@ danger, sia actually sicn average'
elif field1=='sic': simfield1='sicn'; simncfield1='SICN'; simconv1=100
elif field1=='turb': simfield1='turb'; simncfield1='turb'; # the sim var names are placeholders
else: print 'cannot addsims for ' + field1; addsims=False

if field2=='tas': simfield2='st'; simncfield2='ST'
elif field2=='zg50000.00': simfield2='gz50000'; simncfield2='PHI'; simconv2=1/con.get_g()
elif field2=='sia': simfield2='sicn'; simncfield2='SICN'; print '@@ danger, sia actually sicn average'
elif field2=='sic': simfield2='sicn'; simncfield2='SICN'; simconv1=100
elif field2=='turb': simfield2='turb'; simncfield2='turb'; # the sim var names are placeholders
else: print 'cannot addsims for ' + field2; addsims=False

if fieldcnd=='sic': simfieldcnd='sicn'; simncfieldcnd='SICN'; simconvcnd=100
elif fieldcnd=='tas': simfieldcnd='st'; simncfieldcnd='ST'; simconvcnd=1
else: print 'cannot addsims for conditional: ' + fieldcnd; addsims=False


ftype='fullts' # 'fullclimo' or 'climo' or 'fullts'


fdict1 = {'field': field1+region1, 'ncfield': ncfield1, 'comp': comp1}
fdict2 = {'field': field2+region2, 'ncfield': ncfield2, 'comp': comp2}
fdictcnd = {'field': fieldcnd+regioncnd, 'ncfield': ncfieldcnd, 'comp': compcnd}



def subsamp_sims(simsdf,numyrs=11,styears=None,threed=False):
    """ select 11 year segments from given simsdf (data should be anoms)
             simsdf.keys(): sims
             simsdf.index(): time index

             if threed = True: the incoming data (simsdf) should be an ndarray
                     with dims [nsims x ntime x nlat x nlon]. Return data will be
                     nsamp x nlat x nlon

             number of total samples will be determined by length of all sim data
                      and numyrs (e.g. (ntime / numyrs)*numsims)

             returns ndarray of subsample averages, and startyrs
    """
    #print simsdf

    
    if threed:
        (numsims,ntime,nlat,nlon)=simsdf.shape
        initshape=(nlat,nlon)
        keys=np.arange(0,numsims)
    else:
        ntime,numsims=simsdf.values.shape
        initshape=()
        keys=simsdf.keys()

    samp = ntime/numyrs
    allsii=0 # keep track of all sims and all subsamps
    initshape=(samp*numsims,)+initshape

    subsampavg=np.zeros(initshape)

    print 'sample each of ' + str(numsims) + ' sims ' + str(samp) + ' times'
    savstyears = np.zeros((numsims))

    
    for nii,sim in enumerate(keys):
        
        if threed:
            vals = simsdf[nii,...]
        else:
            print sim
            vals=simsdf[sim].values     
        
        if styears == None:
            # random index to start looping, since we have a remainder when ntime/numyrs
            startyr = np.random.randint(np.mod(ntime,numyrs))
        else:
            # start years were passed in: use them
            startyr = styears[nii]

        print 'start ' + str(startyr)
        savstyears[nii] = startyr
        for sii in np.arange(startyr,ntime-numyrs,numyrs):

            #print 'sii ' + str(sii) + ', allsii ' + str(allsii)
            #subsampavg[allsii] = vals[sii+sii*numyrs:sii+sii*numyrs+numyrs].mean()
            subsampavg[allsii,...] = vals[sii:sii+numyrs,...].mean(axis=0)
            allsii+=1

    return subsampavg,savstyears


def subsamp_anom_pi(pidat, numsamp=50, numyrs=11,styear=None,anomyears=None,verb=False):
    """ subsample piControl and produce numsamp anomalies of numyrs-length periods

                 Make sure the data is seasonalized before passing into function.

                 returns subsampanom, styear,anomyears
    """

    ndim=pidat.ndim
    if ndim==1:
        ntime = pidat.shape[0]
        initshape=()
    elif ndim==3:
        (ntime,nlat,nlon)=pidat.shape
        initshape=(nlat,nlon)
    else:
        # what is this mysterious shape?
        (ntime,nspace)=pidat.shape
        

    samp = ntime/numyrs
    allsii=0 # keep track of all subsamps
    anomshape = (numsamp,)+initshape 
    initshape=(samp,)+initshape

    subsampavg=np.zeros(initshape)
    subsampanom=np.zeros(anomshape)

    # chunk up non-overlapping time periods
    # then take anomalies
    if styear == None:
        # random index to start looping, since we have a remainder when ntime/numyrs
        startyr = np.random.randint(np.mod(ntime,numyrs))
    else:
        # start years were passed in: use them
        startyr = styear
    
    if verb:
        print 'start ' + str(startyr)

    for sii in np.arange(startyr,ntime-numyrs,numyrs):

        subsampavg[allsii,...] = pidat[sii:sii+numyrs,...].mean(axis=0)
        allsii+=1

    # now select 2 random time periods to generate anomalies
    # make sure they are at least a decade apart
    anomyrs=[]
    for ii in np.arange(0,numsamp):

        keepgoing=True

        while keepgoing:
            if anomyears == None:
                sel = random.sample(np.arange(0,samp),2) # 2 time periods, no dupes
            else:
                sel = anomyears[ii]

            if sel[0] in np.arange(sel[1]-1,sel[1]+2):
                # should not happen with anomyrs that are passed in.
                if verb:
                    print 'Bad anom time period indices. Keepgoing. ' + str(sel)
            else:
                if verb:
                    print 'anom time period indices: ' + str(sel)
                subsampanom[ii,...] = subsampavg[sel[1],...] - subsampavg[sel[0],...]
                anomyrs.append(sel)
                keepgoing=False
                

    return subsampanom, startyr, anomyrs


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


def plot_shorttermpdf(fig,ax,field,region,xxdat,pdfdat,histdat,meandat,cidat,
                      cifdat,pointdat,startyearsdat=None,addsims=False,addnat=True,addmisc=True,addpi=False,
                      combagcm=False,addanno=True,plab=None,pversion=''):
    """ startyearsdat is an unnecessary input, but is in 'retdict' from calc_shorttermpdf()
    """

    fsz=18

    mew=1.5; ms=9
    firebrick=ccm.get_linecolor('firebrick')
    hcol='0.7'#ccm.get_linecolor('darkolivegreen3')
    hcolline='0.7'#ccm.get_linecolor('darkolivegreen3')#'darkseagreen4')

    #hcol=ccm.get_linecolor('niceblue')
    #hcolline=ccm.get_linecolor('niceblue')
    ltcol='0.5'#ccm.get_linecolor('darkseagreen4')
    ltcol=ccm.get_linecolor('paperblue') #'niceblue2') #'steelblue4')
    natcol=ccm.get_linecolor('steelblue3')
    natcolline=ccm.get_linecolor('steelblue3')#4')
    miscol=ccm.get_linecolor('orange3') #'deepskyblue')
    miscolline=ccm.get_linecolor('orange3') #'deepskyblue')
    obscol='r' #'k' #ccm.get_linecolor('midnightblue') #'b'
    picol='purple'

    sslg=mlines.Line2D([],[],color=ltcol,linewidth=2) # subsamp sims
    ss2lg=mlines.Line2D([],[],color='k',linewidth=2) # subsamp sims2 (variable bc)
    ss2lg2=mlines.Line2D([],[],color='k',linewidth=2,mew=2,marker='|')
    simlg=mlines.Line2D([],[],color='0.3',linestyle='none',marker='o')
    nsidclg=mlines.Line2D([],[],color='g',linestyle='none',marker='o')
    nsidclg2=mlines.Line2D([],[],color='g',linewidth=2,mew=2,marker='|')
    #obslg=mlines.Line2D([],[],color=obscol,linestyle='none',marker='s')
    obslg=mlines.Line2D([],[],color=obscol,linewidth=2)
    rawlg=mlines.Line2D([],[],color=hcolline,linewidth=2)
    rawnlg=mlines.Line2D([],[],color=natcolline,linewidth=2)
    rawmlg=mlines.Line2D([],[],color=miscolline,linewidth=2)
    pilg=mlines.Line2D([],[],color=picol,linewidth=2)

    prstr =''
    #legh = (sslg,rawlg)
    #legstr = ('AGCM Ice','CGCM All')
    # actually add AGCM to end
    legh = (rawlg,)
    legstr = ('CGCM',)

    skey='obs'
    obsreg= np.squeeze(pointdat[skey])
    skey='orig'
    or5sea= np.squeeze(pointdat[skey])
    skey='AGCM120'
    simflddf={}
    simflddf['NSIDC']= np.squeeze(pointdat[skey])

    if not combagcm:
        skey='AGCM2'
    else:
        skey='BOTHAGCM'

    ss2xx = np.squeeze(xxdat[skey]); ss2pdf_fitted=np.squeeze(pdfdat[skey]); 
    ss2mean= np.squeeze(meandat[skey]); ss2ci= np.squeeze(cidat[skey]); ss2cif= np.squeeze(cifdat[skey])
    ss2hist=histdat[skey]
    if not combagcm:
        skey='AGCM'
        ssxx =  np.squeeze(xxdat[skey]); sspdf_fitted= np.squeeze(pdfdat[skey]); 
        ssmean= np.squeeze(meandat[skey]); ssci= np.squeeze(cidat[skey]); sscif= np.squeeze(cifdat[skey])
        
    skey='histle'
    lesea= histdat[skey]; rawxx=xxdat[skey]; rawpdf_fitted=pdfdat[skey]; 
    rawmean=meandat[skey]; rawci=cidat[skey]; rawcif=cifdat[skey]
    if addnat:
        skey='histnat'
        lensea= histdat[skey]; rawnxx=xxdat[skey]; rawnpdf_fitted=pdfdat[skey]; 
        rawnmean=meandat[skey]; rawnci=cidat[skey]; rawncif=cifdat[skey]
    if addmisc:
        skey='histmisc'
        lemsea= histdat[skey]; rawmxx=xxdat[skey]; rawmpdf_fitted=pdfdat[skey]; 
        rawmmean=meandat[skey]; rawmci=cidat[skey]; rawmcif=cifdat[skey]
    if addpi:
        skey='pi'
        pisea= histdat[skey]; pixx=xxdat[skey]; pipdf_fitted=pdfdat[skey]; 
        pimean=meandat[skey]; pici=cidat[skey]; picif=cifdat[skey]

    if field not in ('sia','sic') and pversion != 'c':
        
        ax.hist(ss2hist,normed=True,color=ltcol,alpha=0.4,histtype='stepfilled')#@@ maybe?
        ax.plot(ss2xx,ss2pdf_fitted,color=ltcol,linewidth=3)

        #ax.plot(ss2xx,ss2pdf_fitted,color='k',linewidth=2) # variable boundary forcings

        #ax.hist(osubsamp,normed=True,color='g',alpha=0.3)
        #ax.plot(ossxx,osspdf_fitted,color='g',linewidth=2)
        #prstr=prstr+'addRsimsnsidcsim_'
        #legstr=legstr+('AGCM NSIDC','AGCM var ICE')
        #legh=legh+(nsidclg2,ss2lg)

        #prstr=prstr+'addRsims_'
        #legstr=legstr+('AGCM var ICE',)
        #legh=legh+(ss2lg,)
    if addsims:
        prstr=prstr+'120yrsims_'
        simsy=-0.02
        for es in simsE:
            ax.plot(simflddf[es],simsy,color=ltcol,
                    marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms)
            ax.plot(simfldm[es],simsy,color=ltcol,alpha=0.6,
                    marker='o',linestyle='',fillstyle='full')
        for rs in simsR:
            ax.plot(simflddf[rs],simsy,color='k',
                    marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms)
            ax.plot(simfldm[rs],simsy,color='k',alpha=0.6,
                    marker='o',linestyle='',fillstyle='full') 

        if subnh: # HACK @@
            ax.plot(simflddf['NSIDC'],simsy,color='g',
                    marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms,mec='g')
        else:
            ax.plot(simflddf['NSIDC'],simsy,color='g',
                    marker='o',linestyle='',fillstyle='none',mew=mew,markersize=ms,mec='g')
            ax.plot(simfldm['NSIDC'],simsy,color='g',alpha=0.6,
                    marker='o',linestyle='',fillstyle='full')

    # add RAW LE
    #print 'plotting LE histogram.............' + str(lesea) #@@
    # @@@ NOTES about np.histogram:
    # counts,bins=np.histogram(agcm2,density=True)
    #  The counts only add up to one if bin width is 1.
    #  So np.sum(counts*np.diff(bins)) = 1.0
    #  If want better y-axis ticks for "density" (now not plotted at all)
    #    can do counts*np.diff(bins) as y-ticklabels (I'm guessing) Aug 20 2015 @@

    # ALSO if want better binedges (for AGCM eursat anyway) where no bins
    # cross over zero (encompass zero): (this example assumes 10 bins)
    # binedges=np.arange(-1*max(np.abs(agcm2)), max(np.abs(agcm2))+2*max(np.abs(agcm2))/10., 2*max(np.abs(agcm2))/10.)

    if pversion not in ('d',):
        ax.hist(lesea,normed=True,color=hcol,alpha=0.3,histtype='stepfilled') 
    ax.plot(rawxx,rawpdf_fitted,color=hcolline,linewidth=3)

    fs='full'
    for eii in range(0,5): # orig five
        # search for difforig[eii] nearest to axx, 
        # then plot apdf_fitted[axx] as marker, filled or not based on sig.
        plotx = or5sea[eii]
        idx = cutl.find_nearest(rawxx,plotx)
        ploty = rawpdf_fitted[idx]
        ax.plot(plotx,ploty,marker='o',color=hcolline,mec='k',fillstyle=fs,mew=1,markersize=ms)
        # @@@ check significance?


    if field!='sia':
        nsx= simflddf['NSIDC']
        idx=cutl.find_nearest(ss2xx,nsx)
        nsy=ss2pdf_fitted[idx]
        #ax.plot(nsx,nsy,marker='o',color='g',mec='g',fillstyle=fs,mew=mew,markersize=ms)
        #if addanno:
        #    ax.annotate('Obs AGCM',xy=(nsx,nsy),xytext=(10,20),
        #                textcoords='offset points',fontsize=fsz,
        #                arrowprops=dict(arrowstyle='-',connectionstyle='arc3',
        #                                facecolor='k', edgecolor='k')) # label the marker on the pdf curve

    ylims=ax.get_ylim()
    print '@@@ YLIMS: ' + str(ylims) # @@@@
    if field=='sia':
        yytop = 1.9e-12
        boxywi= 0.002
        yincr=boxywi+0.001
        ax.set_ylim(0,2e-12)
        obsanny=1.8e-12
        offsetpts=(15,0)
    elif field=='tas':
        yytop = 0.95
        boxywi = 0.01
        yincr=boxywi+0.02
        #ax.set_ylim(-0.05,ylims[1])
        obsanny=0.8
        offsetpts=(-50,0)
    elif field == 'sic':
        yytop = 0.19
        boxywi = 0.002
        yincr=boxywi+0.001
        obsanny=0.14
        offsetpts=(15,0)

    if pversion=='c' and field=='tas':
        pass
    else:
        ax.set_ylabel('Density',fontsize=fsz)
    ax.set_yticklabels('')
    ax.axhline(y=0,color='k')
    #ax.axvline(x=0,color='k',linestyle='--')
    
    if pversion in ('c','d'):
        # add uncertainty boxes
        ax.plot((rawmean,rawmean),(yytop-boxywi,yytop+boxywi),linewidth=2,color=hcolline)
        ax.add_patch(mpatches.Rectangle((rawcif[0],yytop-boxywi),
                                        rawcif[1]-rawcif[0],boxywi*2,ec=hcol,fc='white',linewidth=2))
        ax.add_patch(mpatches.Rectangle((rawci[0],yytop-boxywi),
                                        rawci[1]-rawci[0],boxywi*2,ec=hcol,fc=hcol,linewidth=2,alpha=0.5))
        yytop-=yincr
        if pversion in ('d',) and field in ('tas',):
            ax.plot((ss2mean,ss2mean),(yytop-boxywi,yytop+boxywi),linewidth=2,color=ltcol)
            ax.add_patch(mpatches.Rectangle((ss2cif[0],yytop-boxywi),
                                            ss2cif[1]-ss2cif[0],boxywi*2,ec=ltcol,fc='white',linewidth=2))
            ax.add_patch(mpatches.Rectangle((ss2ci[0],yytop-boxywi),
                                            ss2ci[1]-ss2ci[0],boxywi*2,ec=ltcol,fc=ltcol,linewidth=2,alpha=0.5))

        
    if not combagcm:
        xlab = '$\Delta$ SAT ($^\circ$C); AGCM= $%.2f$'%(ssmean)+ ' CGCM= $%.2f$'%(rawmean)
    else:
        xlab = '$\Delta$ SAT ($^\circ$C); BOTHAGCM= $%.2f$'%(ss2mean)+ ' CGCM= $%.2f$'%(rawmean)

    #obslg=ax.axvline(obsreg,color=obscol,linewidth=2) # vert line for obs
    #if field=='sia' and region=='nh':
    if pversion in ('d',) and field=='tas': 
        pass
    else:
        #plotx=obsreg
        #idx=cutl.find_nearest(rawxx,plotx)
        #ploty=rawpdf_fitted[idx]
        #ax.plot(plotx,ploty,marker='o',color=obscol,mec=obscol,fillstyle=fs,mew=mew,markersize=ms)
        #ax.annotate('Obs',xy=(plotx,ploty),xytext=(-28,2),
        #            textcoords='offset points') # label the marker on the pdf curve
        ax.axvline(x=obsreg,color=obscol,linewidth=2)
        if addanno:
            ax.annotate('Obs',xy=(obsreg,obsanny),xytext=offsetpts,
                        textcoords='offset points',
                        arrowprops=dict(arrowstyle='-',connectionstyle='arc3',
                                        facecolor='k', edgecolor='k'),fontsize=fsz) # label vertical obs line

    al =0.9; histtype='step'
    #al =0.3; histtype='stepfilled'
    if addnat:
        #ax.hist(lensea,normed=True,color=natcol,alpha=al,histtype=histtype)
        ax.plot(rawnxx,rawnpdf_fitted,color=natcolline,linewidth=2,alpha=al)
        legh = legh + (rawnlg,)
        legstr = legstr + ('CGCM Nat',)
        prstr = prstr+'nat'
        xlab = xlab + ' Nat=$%.2f$'%(rawnmean) 
    if addmisc:
        #ax.hist(lemsea,normed=True,color=miscol,alpha=al,histtype=histtype)
        ax.plot(rawmxx,rawmpdf_fitted,color=miscolline,linewidth=2,alpha=al)
        legh = legh + (rawmlg,)
        legstr = legstr + ('CGCM Aero',)
        prstr = prstr+'misc'
        xlab = xlab + ' Aero=$%.2f$'%(rawmmean) 
    if addpi:
        if field=='sia':
            # for the paper: add PI hist too, to sia
            ax.hist(pisea,normed=True,color=picol,alpha=0.5,histtype='stepfilled')
        
        ax.plot(pixx,pipdf_fitted,color=picol,linewidth=2,alpha=al)
        legh = legh + (pilg,)
        legstr = legstr + ('PreIndustrial',)
        prstr = prstr+'PI'
        xlab = xlab + ' PI=$%.2f$'%(pimean) 

    xlab = xlab + ' Obs=$%.2f$'%(obsreg) 
    #legh=legh+(obslg,)
    #if field=='sia':
    #    legstr=legstr+('NSIDC Obs',)
    #else:
    #    legstr=legstr+('Obs',)

    if field=='sia':
        ax.set_xlim((-2.5e12,1.6e12))
        units='millions of km$^2$'
        xlab = 'Change in sea ice area, $\Delta$SIC (' + units + ')' # for paper.
        xtlabs = ax.get_xticks()/1e12
        ax.set_xticklabels(xtlabs)
        axesloc=[.71, .69, .18, .2]
    elif field=='tas':
        if pversion in ('d',):
            ax.set_xlim((-2.5,4))
        else:
            ax.set_xlim((-2.5,5))
        units = '$^\circ$C'
        xlab = 'Change in Eurasian surface air temperature, $\Delta$SAT ('+units +')' # for paper.
        axesloc=[.71, .69, .18, .2]
    elif field=='sic':
        units='%'
        xlab = 'Change in Barents/Kara sea ice concentration, $\Delta$SIC (' + units + ')'
        axesloc=[.71, .69, .18, .2]

    #axesloc=[.69, .69, .29, .3] # for inset
    axesloc=[.7, .79, .29, .2] # for inset. smaller version

    """ http://matplotlib.org/examples/pylab_examples/axes_demo.html
    a = axes([.65, .6, .2, .2], axisbg='y')

    """
    if pversion not in ('c','d'):
        # #  INSET!
        # http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
        pos=ax.get_position()
        wid=pos.width
        height=pos.height
        inax_pos = ax.transAxes.transform(axesloc[0:2])
        transFigure = fig.transFigure.inverted()
        infig_pos = transFigure.transform(inax_pos)
        xpos=infig_pos[0]
        ypos=infig_pos[1]
        wid*=axesloc[2] # not sure about this?
        height*=axesloc[3] # and this?
        print 'add_axes: ' + str([xpos,ypos,wid,height])
        inax = fig.add_axes([xpos,ypos,wid,height])

        #inax = fig.add_axes(axesloc)
        top=2.5
        yy=top
        if field!='sia':
            yy=top=3.5
            boxwi=0.25
        else:
            boxwi=0.1
        #inax.plot(rawmean, yy, linestyle='none',marker='s',mec=hcolline,color=hcolline)
        inax.plot((rawmean,rawmean),(yy-boxwi,yy+boxwi),linewidth=2,color=hcolline)
        # from matplotlib.patches import Rectangle
        # someX, someY = 2, 3
        # currentAxis = plt.gca()
        # currentAxis.add_patch(Rectangle((someX - .5, someY - .5), 1, 1, facecolor="grey"))
        inax.add_patch(mpatches.Rectangle((rawcif[0],yy-boxwi),
                                          rawcif[1]-rawcif[0],boxwi*2,ec=hcolline,fc='white',linewidth=2))
        inax.add_patch(mpatches.Rectangle((rawci[0],yy-boxwi),
                                          rawci[1]-rawci[0],boxwi*2,ec=hcolline,fc=hcolline,linewidth=2,alpha=0.5))
        #inax.plot(rawci, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=hcolline,color=hcolline)
        #inax.plot(rawcif, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=hcolline,color=hcolline)
        yy=yy-1
        if addnat:
            inax.plot(rawnmean, yy, linestyle='none',marker='s',mec=natcolline,color=natcolline)
            inax.plot(rawnci, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=natcolline,color=natcolline)
            inax.plot(rawncif, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=natcolline,color=natcolline)
            yy=yy-1
        if addmisc:
            inax.plot(rawmmean, yy, linestyle='none',marker='s',mec=miscolline,color=miscolline)
            inax.plot(rawmci, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=miscolline,color=miscolline)   
            inax.plot(rawmcif, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=miscolline,color=miscolline)   
            yy=yy-1
        if addpi:
            inax.plot((pimean,pimean),(yy-boxwi,yy+boxwi),linewidth=2,color=picol)
            inax.add_patch(mpatches.Rectangle((picif[0],yy-boxwi),
                                              picif[1]-picif[0],boxwi*2,ec=picol,fc='white',linewidth=2))
            inax.add_patch(mpatches.Rectangle((pici[0],yy-boxwi),
                                              pici[1]-pici[0],boxwi*2,ec=picol,fc=picol,linewidth=2,alpha=0.5))
            #inax.plot(pimean, yy, linestyle='none',marker='s',mec=picol,color=picol)
            #inax.plot(pici, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=picol,color=picol)   
            #inax.plot(picif, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=picol,color=picol)   
            yy=yy-1
        if field not in ('sia','sic'):
            inax.plot((ss2mean,ss2mean),(yy-boxwi,yy+boxwi),linewidth=2,color=ltcol)
            inax.add_patch(mpatches.Rectangle((ss2cif[0],yy-boxwi),
                                              ss2cif[1]-ss2cif[0],boxwi*2,ec=ltcol,fc='white',linewidth=2))
            inax.add_patch(mpatches.Rectangle((ss2ci[0],yy-boxwi),
                                              ss2ci[1]-ss2ci[0],boxwi*2,ec=ltcol,fc=ltcol,linewidth=2,alpha=0.5))
            #inax.plot(ss2mean, yy, linestyle='none',marker='s',mec=ltcol,color=ltcol)
            #inax.plot(ss2ci, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=ltcol,color=ltcol)
            #inax.plot(ss2cif, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec=ltcol,color=ltcol)
            yy=yy-1
        ##if field!='sia':
        ##    inax.plot(ss2mean, yy, linestyle='none',marker='s',mec='k',color='k')
        ##    inax.plot(ss2ci, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec='k',color='k')
        ##    inax.plot(ss2cif, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec='k',color='k')
        ##    yy=yy-1
        ##    inax.plot(ossmean, yy, linestyle='none',marker='s',mec='g',color='g')
        ##    inax.plot(ossci, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec='g',color='g')
        ##    inax.plot(osscif, (yy,yy), linewidth=2,marker='|',markersize=6,mew=2,mec='g',color='g')

        #inax.axvline(x=0,linestyle='--',color='k')
        inax.axvline(x=obsreg,color=obscol,linewidth=2,alpha=0.7)
        inax.set_yticklabels('')
        inax.set_yticks([])
        inax.set_ylim(yy+0.5,top+.5)
        if field=='tas' and region=='eurasiamori':
            inax.set_xlim(-1.4,1.8) # taseurasiamori limits
            xticks=inax.get_xticks()
            inax.set_xticklabels(('','-1','','0','','1','','2'))
            #legh=legh+(ss2lg2,nsidclg2)
            #legstr=legstr+('AGCM var ICE','AGCM NSIDC ICE')
        elif field=='sia' and region=='nh':
            inax.set_xlim(-1.4e12,0.1e12) 
            inax.set_xticklabels(('',-1.2,'',-0.8,'',-0.4,'',0))
        elif field=='sic' and region=='bksmori':
            inax.set_xlim(-12,4)
            inax.set_xticklabels((-12,'',-8,'',-4,'',0,'',4))
        inax.set_xlabel('Mean change\n ('+units+')',fontsize=fsz-3)
        # add 'Obs' annotation
        if field=='sia':
            #tlabx=0.5;tlabx2=0.49; tlaby=tlaby2=0.7
            #labx=0.41; laby=0.56
            laby=1.8 # data coord
            offsetpts=(20,15)
        elif field=='tas':
            # NOT WORKING?
            #tlabx=0.02; tlaby=0.85 # text position
            #labx=0.03; laby=0.84; # start of arrow
            #tlabx2=0.25; tlaby2=0.69 # head of arrow
            laby=2.65 # data coord
            offsetpts=(-40,20)
        elif field=='sic':
            laby=1.8
            offsetpts=(20,15)

        if addanno:
            inax.annotate('Obs',xy=(obsreg,laby),xytext=offsetpts,
                          textcoords='offset points',
                          arrowprops=dict(arrowstyle='-',connectionstyle='arc3',
                                          facecolor='k', edgecolor='k'),fontsize=fsz-2) # label vertical obs line
        inax.xaxis.set_ticks_position('bottom')

    if addsims:
        legstr=legstr+('120-yr AGCM',)
        legh=legh+(simlg,)

    fontP = fm.FontProperties()
    fontP.set_size(fsz-2)
        
    # add AGCM Ice to legend
    if field=='tas':
        legloc=(0.71,0.47)
        if pversion in ('c','d'): # no inset
            legloc=(0.65,0.8) #0.71,0.8)
        if pversion not in ('c',):
            legstr=legstr+('AGCM',)# variable',)
            legh=legh+(sslg,)
        if pversion not in ('d',):
            legstr=legstr+('GIStemp',)
            legh=legh+(obslg,)
        ax.legend(legh,legstr, loc=legloc,frameon=False,prop=fontP)
    elif field in ('sia','sic'):
        legstr=legstr+('NSIDC',)
        legh=legh+(obslg,)
        legloc=(0.69,0.5)
        if pversion=='c':
            legloc=(0.69,0.8)
        ax.legend(legh,legstr, loc=legloc,frameon=False,prop=fontP)

    ax.set_xlabel(xlab,fontsize=fsz)
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_tick_params(direction='out')

    # panel labels
    if plab != None:
        ax.annotate(plab,xy=(0.01,1.01),xycoords='axes fraction',fontsize=fsz,fontweight='bold')

    return ax,prstr
# end plot_shorttermpdf() function -------------------------------------

def calc_shorttermpdf(fdict,field,region,sea,timesel,leconv=1,subnh=False,combagcm=False,comblenat=False,
                      addnat=True,addmisc=True,addsims=False,appendorig=False,addpi=False,verb=True):
    """
         returns retdict, simsstr
            retdict = {'xxdat': xxdat, 'pdfdat': pdfdat, 'meandat': meandat,
                       'cidat': cidat, 'cifdat': cifdat, 'histdat': histdat,
                       'pointdat':pointdat, 'startyearsdat': startyearsdat }
    
    """
    #pinumsamp=100; print 'pinumsamp= 100!'
    pinumsamp=50

    timeselc=timesel[0]; timeselp=timesel[1]
    ncfield=fdict['ncfield']
    comp=fdict['comp']

    simconv1=1
    if field=='tas': simfield1='st'; simncfield1='ST'
    elif field=='zg50000.00': simfield1='gz50000'; simncfield1='PHI'; simconv1=1/con.get_g()
    elif field=='sia': simfield1='sicn'; simncfield1='SICN'; print '@@ danger, sia actually sicn average'
    elif field=='sic': simfield1='sicn'; simncfield1='SICN';
    else: print 'cannot addsims for ' + field;

    casename = 'historical'

    # LOAD 1D DATA

    lecdat = le.load_LEdata(fdict,casename,timesel=timeselc, 
                            rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
    (numen,ntime) = lecdat.shape
    lepdat=le.load_LEdata(fdict,casename,timesel=timeselp, 
                          rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
    lecsea = cutl.seasonalize_monthlyts(lecdat.T,season=sea,verb=True).mean(axis=0)
    lepsea = cutl.seasonalize_monthlyts(lepdat.T,season=sea).mean(axis=0)
    lesea = lepsea - lecsea # numens

    # load original five in case:
    or5c=le.load_originalfive(fdict, casename,timesel=timeselc, 
                              rettype='ndarray',conv=leconv,ftype=ftype,verb=verb)
    or5p=le.load_originalfive(fdict, casename,timesel=timeselp, 
                              rettype='ndarray',conv=leconv,ftype=ftype,verb=verb)
    or5csea=cutl.seasonalize_monthlyts(or5c.T,season=sea,verb=True).mean(axis=0)
    or5psea=cutl.seasonalize_monthlyts(or5p.T,season=sea,verb=True).mean(axis=0)
    or5sea = or5psea-or5csea

    if subnh:
        fdictsub = {'field': field+'nh', 'ncfield': ncfield, 'comp': comp}
        subc = le.load_LEdata(fdictsub,casename,timesel=timeselc, 
                              rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
        subp = le.load_LEdata(fdictsub,casename,timesel=timeselp, 
                              rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
        sub = cutl.seasonalize_monthlyts(subp.T,season=sea).T - cutl.seasonalize_monthlyts(subc.T,season=sea).T
        # @@@ change this to subtract the MEAN hemispheric anom from all ens members
        #lesea = lesea - sub.mean(axis=1)
        lesea = lesea - sub.mean(axis=1).mean()

        # load original five in case:
        subor5c=le.load_originalfive(fdictsub, casename,timesel=timeselc, 
                                     rettype='ndarray',conv=leconv,ftype=ftype,verb=verb)
        subor5p=le.load_originalfive(fdictsub, casename,timesel=timeselp, 
                                     rettype='ndarray',conv=leconv,ftype=ftype,verb=verb)
        subor5=cutl.seasonalize_monthlyts(or5csub.T,season=sea,verb=True).T - \
                    cutl.seasonalize_monthlyts(or5psub.T,season=sea,verb=True).T
        or5sea = or5sea - subor5.mean(axis=1).mean() # @@@@ should remove the total mean of all 55 in this case?

    if appendorig:
        lesea = np.hstack((lesea,or5sea))

    # RAW anomalies (decadal diffs)
    # calc the pdf associated with the hist
    rawpdf_fitted,rawmean,rawsd,rawxx = cutl.calc_normfit(lesea) # to get mean
    #rawpdf_fitted,rawxx = cutl.calc_kernel(lesea)
    rawdf = len(lesea)-1
    rawstder = rawsd / np.sqrt(rawdf+1)
    rawci = sp.stats.t.interval(1-cisiglevel, rawdf, loc=rawmean, scale=rawstder)
    rawcif = sp.stats.t.interval(1-cisiglevel, rawdf, loc=rawmean, scale=rawsd)

    if addnat or comblenat:
        lencdat = le.load_LEdata(fdict,'historicalNat',timesel=timeselc, 
                                 rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
        (numenn,ntimen) = lencdat.shape
        lenpdat=le.load_LEdata(fdict,'historicalNat',timesel=timeselp, 
                               rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
        lencsea = cutl.seasonalize_monthlyts(lencdat.T,season=sea).mean(axis=0)
        lenpsea = cutl.seasonalize_monthlyts(lenpdat.T,season=sea).mean(axis=0)
        lensea = lenpsea - lencsea # numens

        if subnh:
            fdictsub = {'field': field+'nh', 'ncfield': ncfield, 'comp': comp}
            nsubc = le.load_LEdata(fdictsub,'historicalNat',timesel=timeselc, 
                                   rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
            nsubp = le.load_LEdata(fdictsub,'historicalNat',timesel=timeselp, 
                                   rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
            nsub = cutl.seasonalize_monthlyts(nsubp.T,season=sea).T - cutl.seasonalize_monthlyts(nsubc.T,season=sea).T
            # @@@ change this to subtract the MEAN hemispheric anom from all ens members
            #lesea = lesea - sub.mean(axis=1)
            lensea = lensea - nsub.mean(axis=1).mean()

        if comblenat:
            # combine them and then do the fit
            lemean=lesea.mean()
            lenmean=lensea.mean()
            # remove nat mean and add in historical mean
            lensea=lensea-lenmean+lemean
            lesea = np.hstack((lesea,lensea))
            print '@@@@@@@@@@@@@@@@@@' + str(lesea.shape)
            
            rawpdf_fitted,rawmean,rawsd,rawxx = cutl.calc_normfit(lesea) # to get mean
            #rawpdf_fitted,rawxx = cutl.calc_kernel(lesea)
            rawdf = len(lesea)-1
            rawstder = rawsd / np.sqrt(rawdf+1)
            rawci = sp.stats.t.interval(1-cisiglevel, rawdf, loc=rawmean, scale=rawstder)
            rawcif = sp.stats.t.interval(1-cisiglevel, rawdf, loc=rawmean, scale=rawsd)
        else:
            rawnpdf_fitted,rawnmean,rawnsd,rawnxx = cutl.calc_normfit(lensea)
            #rawnpdf_fitted,rawnxx = cutl.calc_kernel(lensea)
            rawndf = len(lensea)-1
            rawnstder = rawnsd / np.sqrt(rawndf+1)
            rawnci = sp.stats.t.interval(1-cisiglevel, rawndf, loc=rawnmean, scale=rawnstder)
            rawncif = sp.stats.t.interval(1-cisiglevel, rawndf, loc=rawnmean, scale=rawnsd)

    if addmisc:
        lemcdat = le.load_LEdata(fdict,'historicalMisc',timesel=timeselc, 
                                 rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
        (numenm,ntimem) = lemcdat.shape
        lempdat=le.load_LEdata(fdict,'historicalMisc',timesel=timeselp, 
                               rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
        lemcsea = cutl.seasonalize_monthlyts(lemcdat.T,season=sea).mean(axis=0)
        lempsea = cutl.seasonalize_monthlyts(lempdat.T,season=sea).mean(axis=0)
        lemsea = lempsea - lemcsea # numens

        if subnh:
            fdictsub = {'field': field+'nh', 'ncfield': ncfield, 'comp': comp}
            msubc = le.load_LEdata(fdictsub,'historicalMisc',timesel=timeselc, 
                                   rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
            msubp = le.load_LEdata(fdictsub,'historicalMisc',timesel=timeselp, 
                                   rettype='ndarray',conv=leconv,ftype=ftype,local=local,verb=verb)
            msub = cutl.seasonalize_monthlyts(msubp.T,season=sea).T - cutl.seasonalize_monthlyts(msubc.T,season=sea).T
            # @@@ change this to subtract the MEAN hemispheric anom from all ens members
            #lesea = lesea - sub.mean(axis=1)
            lemsea = lemsea - msub.mean(axis=1).mean()

        rawmpdf_fitted,rawmmean,rawmsd,rawmxx = cutl.calc_normfit(lemsea)
        #rawmpdf_fitted,rawmxx = cutl.calc_kernel(lemsea)
        rawmdf = len(lemsea)-1
        rawmstder = rawmsd / np.sqrt(rawmdf+1)
        rawmci = sp.stats.t.interval(1-cisiglevel, rawmdf, loc=rawmmean, scale=rawmstder)
        rawmcif = sp.stats.t.interval(1-cisiglevel, rawmdf, loc=rawmmean, scale=rawmsd)

    if addpi:
        pidat = lcd.load_data(fdict,'piControl',local=local,conv=leconv,detrend=False,verb=verb)
        piseadat = cutl.seasonalize_monthlyts(pidat,season=sea)
        piseadat = cutl.detrend(piseadat,axis=0)

        # this data needs to be put into anomalies
        # for now, set years to none
        styear=styearPI; anomyears=anomyearsPI
        #styear=None; anomyears=None #@@@@@
        #print '@@@@@@@@@ new PI subsampling @@@@@@@@'
        # the data must be seasonalized before using this func.
        pisea,styear,anomyears = subsamp_anom_pi(piseadat, numyrs=11,numsamp=pinumsamp,
                                                 styear=styear,anomyears=anomyears)

        if subnh:
            fdictsub = {'field': field+'nh', 'ncfield': ncfield, 'comp': comp}
            pidatsub = lcd.load_data(fdictsub,'piControl',local=local,conv=leconv,verb=verb)
            piseadatsub = cutl.seasonalize_monthlyts(pidatsub,season=sea)
            piseadatsub = cutl.detrend(piseadatsub,axis=0)
            piseasub,styear,anomyears = subsamp_anom_pi(piseadatsub, numyrs=11,numsamp=pinumsamp,
                                                        styear=styear,anomyears=anomyears)
            pisea = pisea - piseasub.mean(axis=0) # ??

        pipdf_fitted,pimean,pisd,pixx = cutl.calc_normfit(pisea)
        pidf = len(pisea)-1
        pistder = pisd / np.sqrt(pidf+1)
        pici = sp.stats.t.interval(1-cisiglevel, pidf, loc=pimean, scale=pistder)
        picif = sp.stats.t.interval(1-cisiglevel, pidf, loc=pimean, scale=pisd)

        # save styear and anomyears for mat
        styearpi = styear; anomyearspi = anomyears

    """if longtermLE:
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

        ltpdf_fitted,ltmean,ltsd,ltxx = cutl.calc_normfit(ltavg)"""


    #else: # not longtermLE (short term sims instead)

    # here we want to subsample 11-year segments from the sims
    
    sims=('E1','E2','E3','E4','E5'); simsstr='sbEonly' # sub
    #sims=('R1','R2','R3','R4','R5'); simsstr='sbRonly'
    if combagcm:
        sims=('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5'); simsstr='sbERboth'
    simflddf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                         region=region))*simconv1
    osimflddf = pd.DataFrame(lmd.loaddata((simfield1,),('NSIDC',),ncfields=(simncfield1,), timefreq=sea, 
                                         region=region))*simconv1
    if subnh:
        simsubdf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                             filetype='diff',region='nh'))*simconv1
        # want to subtract the mean hemispheric avg anomaly
        simflddf = simflddf - simsubdf.mean(axis=1).mean() # average over sims and then time (should be a scalar)
        osimsubdf = pd.DataFrame(lmd.loaddata((simfield1,),('NSIDC',),ncfields=(simncfield1,), timefreq=sea, 
                                             filetype='diff',region='nh'))*simconv1
        # want to subtract the mean hemispheric avg anomaly
        osimflddf = osimflddf - osimsubdf.mean(axis=1).mean() # average over sims and then time (should be a scalar)

    subsamp,styearsss = subsamp_sims(simflddf,numyrs=11,styears=styearsE)
    plotsims = subsamp
    print '==== AGCM '
    sspdf_fitted,ssmean,sssd,ssxx = cutl.calc_normfit(plotsims)
    #sspdf_fitted,ssxx = cutl.calc_kernel(plotsims)
    ssdf = len(subsamp)-1
    ssstder = sssd / np.sqrt(ssdf+1)

    # mean and conf-int: @@@@@
    # ci = sp.stats.t.interval(1-cisiglevel, df, loc= meananom, scale=stder)
    ssci = sp.stats.t.interval(1-cisiglevel, ssdf, loc=ssmean, scale=ssstder) # this is mean value's 5-95%
    sscif = sp.stats.t.interval(1-cisiglevel, ssdf, loc=ssmean, scale=sssd) # full range. don't divide by n. this is pdf 5-95%

    osubsamp,styearsns = subsamp_sims(osimflddf,numyrs=11,styears=styearsN)
    print '==== NSIDC AGCM '
    osspdf_fitted,ossmean,osssd,ossxx = cutl.calc_normfit(osubsamp)
    ossdf = len(osubsamp)-1
    ossstder = osssd / np.sqrt(ossdf+1)
    ossci = sp.stats.t.interval(1-cisiglevel, ossdf, loc=ossmean, scale=ossstder)
    osscif = sp.stats.t.interval(1-cisiglevel, ossdf, loc=ossmean, scale=osssd)

    if not combagcm:
        # do each ens separately, for testing stats (maybe not plotting)
        sims2=('R1','R2','R3','R4','R5'); #simsstr='sbRonly'
        sim2flddf = pd.DataFrame(lmd.loaddata((simfield1,),sims2,ncfields=(simncfield1,), timefreq=sea, 
                                             region=region))*simconv1
        if subnh:
            sim2subdf = pd.DataFrame(lmd.loaddata((simfield1,),sims2,ncfields=(simncfield1,), timefreq=sea, 
                                                 filetype='diff',region='nh'))*simconv1

            # want to subtract the mean hemispheric avg anomaly
            sim2flddf = sim2flddf - sim2subdf.mean(axis=1).mean() # average over sims and then time (should be a scalar)


        subsamp2,styears2 = subsamp_sims(sim2flddf,numyrs=11,styears=styearsR)
        plotsims2 = subsamp2
        print '==== AGCM2 '
        ss2pdf_fitted,ss2mean,ss2sd,ss2xx = cutl.calc_normfit(plotsims2)
        #ss2pdf_fitted,ss2xx = cutl.calc_kernel(plotsims2)
        ss2df = len(subsamp2)-1
        ss2stder = ss2sd / np.sqrt(ss2df+1)
        ss2ci = sp.stats.t.interval(1-cisiglevel, ss2df, loc=ss2mean, scale=ss2stder)
        ss2cif = sp.stats.t.interval(1-cisiglevel, ss2df, loc=ss2mean, scale=ss2sd)

        tstat, pval = sp.stats.ttest_ind(plotsims,plotsims2)
        lstat, lpval = sp.stats.levene(plotsims,plotsims2)
        print '==== testing ANT vs TOT ===='
        print 'TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
        if pval<=siglevel:
             print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'
        print 'LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
        if lpval<=siglevel:
             print 'The ensemble variances are significantly different (' + str(1-siglevel) + ')'


    tstat, pval = sp.stats.ttest_ind(plotsims,lesea)
    lstat, lpval = sp.stats.levene(plotsims,lesea)
    print '==== testing ANT vs ALL LE ===='
    print 'TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
    if pval<=siglevel:
         print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'
    print 'LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
    if lpval<=siglevel:
         print 'The ensemble variances are significantly different (' + str(1-siglevel) + ')'

    if addnat:
        tstat, pval = sp.stats.ttest_ind(lesea,lensea)
        lstat, lpval = sp.stats.levene(lesea,lensea)
        print '==== testing LE vs Nat LE ===='
        print 'TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
        if pval<=siglevel:
             print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'
        print 'LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
        if lpval<=siglevel:
             print 'The ensemble variances are significantly different (' + str(1-siglevel) + ')'

    if addmisc:
        tstat, pval = sp.stats.ttest_ind(lesea,lemsea)
        lstat, lpval = sp.stats.levene(lesea,lemsea)
        print '==== testing LE vs Aero LE ===='
        print 'TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
        if pval<=siglevel:
             print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'
        print 'LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
        if lpval<=siglevel:
             print 'The ensemble variances are significantly different (' + str(1-siglevel) + ')'

    if addsims:
        # ========= add sims (120-yr averages)
        sims=('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5','NSIDC'); #simsstr=''
        #simsRN=('R1','R2','R3','R4','R5','NSIDC')
        #simsR=('R1','R2','R3','R4','R5')
        #simsE=('E1','E2','E3','E4','E5')
        #sims=simsRN; simsstr='Ronly'
        #simflddf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
        #                                     meantype='time',region=region),index=sims)*simconv1
        simflddf = pd.Series(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                             meantype='time',region=region))*simconv1

        # --- for estimating sigma in Individual SIC forcing ensemble only:
        #simflddfr = pd.DataFrame(lmd.loaddata((simfield1,),simsR,ncfields=(simncfield1,), timefreq=sea, 
        #                                      meantype='time',region=region),index=simsR)*simconv1
        simflddfr = pd.Series(lmd.loaddata((simfield1,),simsR,ncfields=(simncfield1,), timefreq=sea, 
                                              meantype='time',region=region))*simconv1
        simval=simflddfr.values
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
            simfldm =simflddf.mask(simpvdf>siglevel)# ma.masked_where(pvals>siglevel,simflddf.values)
            # have to calc pvals 

        else:
            simpvdf = pd.Series(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                                filetype='pval',region=region))#,index=sims)
            simfldm = simflddf.mask(simpvdf>siglevel)
            #simfldm = ma.masked_where(simpvdf.values>siglevel,simflddf.values)

    else: # just add NSIDC
        sims=('NSIDC',)
        simflddf = pd.Series(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                             meantype='time',region=region))*simconv1
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
            simfldm =simflddf.mask(simpvdf>siglevel)# ma.masked_where(pvals>siglevel,simflddf.values)

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
            ttlstr=ttlstr+'-NH '
    elif field in ('sia','sic'):
        nsidcfile = '/HOME/rkm/work/BCs/NSIDC/td_bootstrap_197811_latest_128_64_sicn_1978111600-2013121612.nc'
        latns=cnc.getNCvar(nsidcfile,'lat')
        lonns=cnc.getNCvar(nsidcfile,'lon')
        nsidcsel='1979-01-01,2011-12-31'
        nsidcyrs=np.arange(1979,2012)
        if sea=='DJF':
            nsidcyrs=nsidcyrs[:-1]
        nsidcfldc=cnc.getNCvar(nsidcfile,'SICN',timesel=timeselc,seas=sea)
        nsidcfldp=cnc.getNCvar(nsidcfile,'SICN',timesel=timeselp,seas=sea)
        if field == 'sia':
            obsreg=cutl.calc_regtotseaicearea(nsidcfldp,latns,lonns,region=region,isarea=False).mean() -\
                    cutl.calc_regtotseaicearea(nsidcfldc,latns,lonns,region=region,isarea=False).mean()
        else:
            obsreg=cutl.calc_regmean(nsidcfldp,latns,lonns,region=region).mean() -\
                    cutl.calc_regmean(nsidcfldc,latns,lonns,region=region).mean()
            obsreg = obsreg*100 # convert to %

    xxdat={}; pdfdat={}; meandat={}; cidat={}; cifdat={}; histdat={}; pointdat={}; startyearsdat={}
    # save data to call plot function
    if not combagcm:
        skey='AGCM2'
        histdat[skey]=plotsims2; xxdat[skey]=ss2xx; pdfdat[skey]=ss2pdf_fitted
        meandat[skey]=ss2mean; cidat[skey]=ss2ci; cifdat[skey]=ss2cif
        startyearsdat[skey]=styears2
        skey='AGCM'
    else:
        skey='BOTHAGCM'
    histdat[skey]=plotsims; xxdat[skey]=ssxx; pdfdat[skey]=sspdf_fitted
    meandat[skey]=ssmean; cidat[skey]=ssci; cifdat[skey]=sscif; startyearsdat[skey]=styearsss

    skey='histle'
    histdat[skey]=lesea; xxdat[skey]=rawxx; pdfdat[skey]=rawpdf_fitted
    meandat[skey]=rawmean; cidat[skey]=rawci; cifdat[skey]=rawcif; #startyearsdat[skey]=@@placeholder
    if addnat:
        skey='histnat'
        histdat[skey]=lensea; xxdat[skey]=rawnxx; pdfdat[skey]=rawnpdf_fitted
        meandat[skey]=rawnmean; cidat[skey]=rawnci; cifdat[skey]=rawncif; #startyearsdat[skey]=@@placeholder
    if addmisc:
        skey='histmisc'
        histdat[skey]=lemsea; xxdat[skey]=rawmxx; pdfdat[skey]=rawmpdf_fitted
        meandat[skey]=rawmmean; cidat[skey]=rawmci; cifdat[skey]=rawmcif; #startyearsdat[skey]=@@placeholder
    if addpi:
        skey='pi'
        histdat[skey]=pisea; xxdat[skey]=pixx; pdfdat[skey]=pipdf_fitted
        meandat[skey]=pimean; cidat[skey]=pici; cifdat[skey]=picif; 
        startyearsdat[skey]=styearpi; startyearsdat[skey+'anomyrs'] = anomyearspi

    pointdat['obs']=obsreg
    pointdat['orig']=or5sea
    pointdat['AGCM120']=simflddf['NSIDC'] # hack: to save to mat file has to be an array only

    retdict = {'xxdat': xxdat, 'pdfdat': pdfdat, 'meandat': meandat,
               'cidat': cidat, 'cifdat': cifdat, 'histdat': histdat,
               'pointdat':pointdat, 'startyearsdat': startyearsdat}

    return retdict, simsstr

# end calc_shorttermpdf() function ------------------------------------------------------------------------


# ================================================================
# ==================== main() ====================================
def main(dowhat=None,addobs=True,addsims=False,addnat=False,
         addmisc=False,addpi=False,verb=False,combagcm=False,
         comblenat=False,pversion='',local=True,printtofile=False):


    """ dowhat options are: doscatter, dohist, doregress, dolongtermavg
             default is doscatter
    """
    doscatter=dohist=doregress=dolongtermavg=False

    #addobs=True # @@@ why didn't global var work?
    #addsims=True

    if dowhat==None or dowhat=='doscatter':
        doscatter=True
    elif dowhat=='dohist':
        dohist=True
    elif dowhat=='doregress':
        doregress=True
    elif dowhat=='dolongtermavg':
        dolongtermavg=True

    if doscatter:
        #printtofile=True
        #addnat=True
        #addmisc=True
        savedat=True # save data to ascii for John
        #pinumsamp=100; print 'pinumsamp= 100!!!' # @@
        pinumsamp=50
        shorttermsims=True
        if field2=='tas':
            if region2=='eurasiamori':
                ymin=-2; ymax=3 # for Eurasia SAT v BKS Z500
            elif performop2 and region2op=='gm':
                ymin=1; ymax=5
        elif field2=='turb':
            ymin=-15; ymax=25
        elif field2=='zg50000.00':
            ymin=-60; ymax=80

        # historical
        casename='historical'
        lecdat1 = le.load_LEdata(fdict1,casename,timesel=timeselc, 
                                 rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
        (numens1,ntime1) = lecdat1.shape
        lepdat1=le.load_LEdata(fdict1,casename,timesel=timeselp, 
                               rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
        lecsea1 = cutl.seasonalize_monthlyts(lecdat1.T,season=sea1).T
        lepsea1 = cutl.seasonalize_monthlyts(lepdat1.T,season=sea1).T
        lesea1 = lepsea1 - lecsea1

        if performop1:
            try:
                fdict1op = {'field': field1+region1op, 'ncfield': ncfield1, 'comp': comp1}
                # should rename these variables to be 'op' so it's more general
                subc1 = le.load_LEdata(fdict1op,casename,timesel=timeselc, 
                                       rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
                subp1 = le.load_LEdata(fdict1op,casename,timesel=timeselp, 
                                       rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
                sub1 = cutl.seasonalize_monthlyts(subp1.T,season=sea1).T - \
                       cutl.seasonalize_monthlyts(subc1.T,season=sea1).T

                
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

        lecdat2 = le.load_LEdata(fdict2,casename,timesel=timeselc, 
                                 rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
        (numens2,ntime2) = lecdat2.shape
        lepdat2=le.load_LEdata(fdict2,casename,timesel=timeselp, 
                               rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
        lecsea2 = cutl.seasonalize_monthlyts(lecdat2.T,season=sea2).T
        lepsea2 = cutl.seasonalize_monthlyts(lepdat2.T,season=sea2).T
        lesea2 = lepsea2 - lecsea2
        if performop2:
            
            fdict2sub = {'field': field2+region2op, 'ncfield': ncfield2, 'comp': comp2}
            subc2 = le.load_LEdata(fdict2sub,casename,timesel=timeselc, 
                                   rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
            subp2 = le.load_LEdata(fdict2sub,casename,timesel=timeselp, 
                                   rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
            sub2 = cutl.seasonalize_monthlyts(subp2.T,season=sea2).T -\
                   cutl.seasonalize_monthlyts(subc2.T,season=sea2).T

            print '######### gm tas: ' + str(sub2.mean(axis=1))

            if op2=='sub': # subtract
                lefld2 = lesea2.mean(axis=1) - sub2.mean(axis=1)
            elif op2=='div': # divide
                lefld2 = lesea2.mean(axis=1) / sub2.mean(axis=1)
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
                
                fdict2sub = {'field': field2+region2op, 'ncfield': ncfield2, 'comp': comp2}
                subc2n = le.load_LEdata(fdict2sub,casename2,timesel=timeselc, 
                                        rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
                subp2n = le.load_LEdata(fdict2sub,casename2,timesel=timeselp, 
                                        rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
                sub2n = cutl.seasonalize_monthlyts(subp2n.T,season=sea2).T -\
                        cutl.seasonalize_monthlyts(subc2n.T,season=sea2).T

                if op2=='sub': # subtract
                    lefld2n = lesea2n.mean(axis=1) - sub2n.mean(axis=1)
                elif op2=='div': # divide
                    lefld2n = lesea2n.mean(axis=1) / sub2n.mean(axis=1)
                    
            else:
                lefld2n=lepsea2n.mean(axis=1)-lecsea2n.mean(axis=1)

            lemmn, lebbn, lervaln, lepvaln, lestd_errn = sp.stats.linregress(lefld1n,lefld2n)

            if conditional:
                lecdatncnd = le.load_LEdata(fdictcnd,casename2,timesel=timeselc, 
                                           rettype='ndarray',conv=leconvcnd,ftype=ftype,local=local)
                (numensncnd,ntimencnd) = lecdatncnd.shape
                lepdatncnd=le.load_LEdata(fdictcnd,casename2,timesel=timeselp, 
                                         rettype='ndarray',conv=leconvcnd,ftype=ftype,local=local)
                lecseancnd = cutl.seasonalize_monthlyts(lecdatncnd.T,season=seacnd).T
                lepseancnd = cutl.seasonalize_monthlyts(lepdatncnd.T,season=seacnd).T
                leseancnd = lepseancnd - lecseancnd
                lefldncnd=lepseancnd.mean(axis=1)-lecseancnd.mean(axis=1)


        if addmisc:
            # historicalMisc
            casename3='historicalMisc'
            lecdat1m = le.load_LEdata(fdict1,casename3,timesel=timeselc, 
                                      rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
            (numens1m,ntime1m) = lecdat1m.shape
            lepdat1m=le.load_LEdata(fdict1,casename3,timesel=timeselp, 
                                    rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
            lecsea1m = cutl.seasonalize_monthlyts(lecdat1m.T,season=sea1).T
            lepsea1m = cutl.seasonalize_monthlyts(lepdat1m.T,season=sea1).T
            lesea1m = lepsea1m - lecsea1m
            if performop1:
                fdict1op = {'field': field1+region1op, 'ncfield': ncfield1, 'comp': comp1}
                subc1m = le.load_LEdata(fdict1op,casename3,timesel=timeselc, 
                                        rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
                subp1m = le.load_LEdata(fdict1op,casename3,timesel=timeselp, 
                                        rettype='ndarray',conv=leconv1,ftype=ftype,local=local)
                sub1m = cutl.seasonalize_monthlyts(subp1m.T,season=sea1).T -\
                        cutl.seasonalize_monthlyts(subc1m.T,season=sea1).T

                if op1=='sub': # subtract
                    lefld1m = lesea1m.mean(axis=1) - sub1m.mean(axis=1)
                elif op1=='div': # divide
                    lefld1m = lesea1m.mean(axis=1) / sub1m.mean(axis=1)
            else:
                lefld1m=lepsea1m.mean(axis=1)-lecsea1m.mean(axis=1)

            lecdat2m = le.load_LEdata(fdict2,casename3,timesel=timeselc, 
                                      rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
            (numens2m,ntime2m) = lecdat1m.shape
            lepdat2m=le.load_LEdata(fdict2,casename3,timesel=timeselp, 
                                    rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
            lecsea2m = cutl.seasonalize_monthlyts(lecdat2m.T,season=sea2).T
            lepsea2m = cutl.seasonalize_monthlyts(lepdat2m.T,season=sea2).T
            lesea2m = lepsea2m - lecsea2m
            if performop2:
                
                fdict2sub = {'field': field2+region2op, 'ncfield': ncfield2, 'comp': comp2}
                subc2m = le.load_LEdata(fdict2sub,casename3,timesel=timeselc, 
                                        rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
                subp2m = le.load_LEdata(fdict2sub,casename3,timesel=timeselp, 
                                        rettype='ndarray',conv=leconv2,ftype=ftype,local=local)
                sub2m = cutl.seasonalize_monthlyts(subp2m.T,season=sea2).T -\
                        cutl.seasonalize_monthlyts(subc2m.T,season=sea2).T

                if op2=='sub': # subtract
                    lefld2m = lesea2m.mean(axis=1) - sub2m.mean(axis=1)
                elif op2=='div': # divide
                    lefld2m = lesea2m.mean(axis=1) / sub2m.mean(axis=1)
            else:
                lefld2m=lepsea2m.mean(axis=1)-lecsea2m.mean(axis=1)

            lemmm, lebbm, lervalm, lepvalm, lestd_errm = sp.stats.linregress(lefld1m,lefld2m)

            if conditional:
                lecdatmcnd = le.load_LEdata(fdictcnd,casename3,timesel=timeselc, 
                                           rettype='ndarray',conv=leconvcnd,ftype=ftype,local=local)
                (numensmcnd,ntimemcnd) = lecdatmcnd.shape
                lepdatmcnd=le.load_LEdata(fdictcnd,casename3,timesel=timeselp, 
                                         rettype='ndarray',conv=leconvcnd,ftype=ftype,local=local)
                lecseamcnd = cutl.seasonalize_monthlyts(lecdatmcnd.T,season=seacnd).T
                lepseamcnd = cutl.seasonalize_monthlyts(lepdatmcnd.T,season=seacnd).T
                leseamcnd = lepseamcnd - lecseamcnd
                lefldmcnd=lepseamcnd.mean(axis=1)-lecseamcnd.mean(axis=1)

        if addpi:
            # piControl
            pidat1 = lcd.load_data(fdict1,'piControl',conv=leconv1,local=local)
            ntime1pi = pidat1.shape

            piseadat1 = cutl.seasonalize_monthlyts(pidat1,season=sea1)
            piseadat1 = cutl.detrend(piseadat1,axis=0)
            styear=styearPI; anomyears=anomyearsPI
            # the data must be seasonalized before using this func.
            pisea1,styear1,anomyears1 = subsamp_anom_pi(piseadat1, numyrs=11,numsamp=pinumsamp,
                                                        styear=styear,anomyears=anomyears)

            if verb:
                print 'pidat1.shape, piseadat1.shape ' + str(pidat1.shape) + ',' + str(piseadat1.shape) # @@@
                print 'pisea1.shape ' + str(pisea1.shape) # @@@

            if performop1:
                fdict1op = {'field': field1+region1op, 'ncfield': ncfield1, 'comp': comp1}
                pisub1 = lcd.load_data(fdict1op,'piControl',conv=leconv1,local=local)
                pisubsea1 = cutl.seasonalize_monthlyts(pisub1,season=sea1)
                pisubsea1 = cutl.detrend(pisubsea1,axis=0)
                pisubsea1,styear1,anomyears1 = subsamp_anom_pi(pisubsea1, numyrs=11,numsamp=pinumsamp,
                                                               styear=styear1,anomyears=anomyears1)
                

                if op1=='sub': # subtract
                    pifld1 = pisea1.mean(axis=1) - pisubsea1.mean(axis=1)
                elif op1=='div': # divide
                    pifld1 = pisea1.mean(axis=1) / pisubsea1.mean(axis=1)
            else:
                pifld1 = pisea1

            pidat2 = lcd.load_data(fdict2,'piControl',conv=leconv2,local=local)
            ntime2pi = pidat2.shape

            piseadat2 = cutl.seasonalize_monthlyts(pidat2,season=sea1)
            piseadat2 = cutl.detrend(piseadat2,axis=0)
            # the data must be seasonalized before using this func.
            pisea2,styear1,anomyears1 = subsamp_anom_pi(piseadat2, numyrs=11,numsamp=pinumsamp,
                                                        styear=styear1,anomyears=anomyears1)
            if performop2:
                fdict2op = {'field': field2+region2op, 'ncfield': ncfield2, 'comp': comp2}
                pisub2 = lcd.load_data(fdict2op,'piControl',conv=leconv2,local=local)
                pisubsea2 = cutl.seasonalize_monthlyts(pisub2,season=sea2)
                pisubsea2 = cutl.detrend(pisubsea2,axis=0)
                pisubsea2,styear1,anomyears1 = subsamp_anom_pi(pisubsea2, numyrs=11,numsamp=pinumsamp,
                                                               styear=styear1,anomyears=anomyears1)

                if op2=='sub': # subtract
                    pifld2 = pisea2 - pisubsea2
                elif op2=='div': # divide
                    pifld2 = pisea2 / pisubsea2
            else:
                pifld2 = pisea2
                
            if verb:
                print 'pidat2.shape, piseadat2.shape ' + str(pidat2.shape) + ',' + str(piseadat2.shape) # @@@
                print 'pisea2.shape ' + str(pisea2.shape) # @@@

            pimm, pibb, pirval, pipval, pistd_err = sp.stats.linregress(pifld1,pifld2)


        if addobs: # DATA MUST EXIST:hard coded for TAS and Z500 @@
            graveraint= 9.80665 # m/s2 (different from Canadian models)

            if field1 == 'zg50000.00':
                erafile = '/HOME/rkm/work/DATA/ERAINT/td_era_int_197901_201412_gp_128_64_phi500_1979011612-2014121612.nc'
                eraz500c= cnc.getNCvar(erafile,'PHI',timesel=timeselc,seas=sea1)/graveraint
                eraz500p= cnc.getNCvar(erafile,'PHI',timesel=timeselp,seas=sea1)/graveraint
                latera=cnc.getNCvar(erafile,'lat')
                lonera=cnc.getNCvar(erafile,'lon')
                obsreg1 = cutl.calc_regmean(eraz500p-eraz500c,latera,lonera,region1,model=None)
                if performop1:
                    opfld1 = cutl.calc_regmean(eraz500p-eraz500c,latera,lonera,region1op,model=None)

                    if op1=='sub': # subtract
                        obsreg1 = obsreg1.mean() - opfld1.mean()
                    elif op1=='div': # divide
                        obsreg1 = obsreg1.mean() / opfld1.mean()
                else:
                    obsreg1=obsreg1.mean()
            elif field1 == 'tas':
                gisfile = '/HOME/rkm/work/DATA/GISS/gistemp1200_ERSST.nc'
                latgis=cnc.getNCvar(gisfile,'lat')
                longis=cnc.getNCvar(gisfile,'lon')
                gissatc= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselc,seas=sea1) 
                gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselp,seas=sea1)
                obsreg1 =  cutl.calc_regmean(gissatp-gissatc,latgis,longis,region1,model=None)
                if performop1:
                    opfld1 = cutl.calc_regmean(gissatp-gissatc,latgis,longis,region1op,model=None)

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
                print 'OBS SIC: SELECT UP THROUGH 2011 ONLY @@@@@'
                nssicc= cnc.getNCvar(nsidcfile,'SICN',timesel=timeselc,seas=sea1)*100
                nssicp= cnc.getNCvar(nsidcfile,'SICN',timesel='2002-01-01,2011-12-31',seas=sea1)*100
                obsreg1 =  cutl.calc_regmean(nssicp-nssicc.mean(axis=0),latns,lonns,region1,model=None)
                if performop1:
                    opfld1 = cutl.calc_regmean(nssicp-nssicc.mean(axis=0),latns,lonns,region1op,model=None)

                    if op1=='sub': # subtract
                        print '========== OBS SIC: ' + str(obsreg1.mean()) + '-' + str(opfld1.mean()) +\
                            ' = ' + str(obsreg1.mean() - opfld1.mean())
                        obsreg1 = obsreg1.mean() - opfld1.mean()
                    elif op1=='div': # divide
                        obsreg1 = obsreg1.mean() / opfld1.mean()
                else:
                    obsreg1=obsreg1.mean()
                print 'obsreg1 (sic): ' + str(obsreg1)
            else:
                print 'this field not supported with addobs @@: ' + field1; addobs=False

            if field2 == 'tas':
                gisfile = '/HOME/rkm/work/DATA/GISS/gistemp1200_ERSST.nc'
                latgis=cnc.getNCvar(gisfile,'lat')
                longis=cnc.getNCvar(gisfile,'lon')
                gissatc= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselc,seas=sea2) 
                if field1=='sic':
                    print 'OBS TAS: SELECT UP THROUGH 2011 ONLY to match SIC @@@@@'
                    gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel='2002-01-01,2011-12-31',seas=sea2)
                else:
                    gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselp,seas=sea2)

                if sea1==sea2=='DJF' and performop2==True and op2=='sub' and region2=='eurasiamori' and region2op=='nh':
                    # a file already exists for this
                    gisfile='/HOME/rkm/work/DATA/GISS/giss_DJF_eurasiamori-nh_1979-2013_timeseries.nc'
                    gissatc= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselc) 
                    gissatp= cnc.getNCvar(gisfile,'tempanomaly',timesel=timeselp)
                    obsreg2=gissatp.mean()-gissatc.mean()
                else:
                    obsreg2 =  cutl.calc_regmean(gissatp-gissatc.mean(axis=0),latgis,longis,region2,model=None)
                    if performop2:
                        opfld2 = cutl.calc_regmean(gissatp-gissatc.mean(axis=0),latgis,longis,region2op,model=None)

                        if op2=='sub': # subtract
                            obsreg2 = obsreg2.mean() - opfld2.mean()
                        elif op2=='div': # divide
                            obsreg2 = obsreg2.mean() / opfld2.mean()
                    else:
                        obsreg2=obsreg2.mean()
                    
            else:
                print 'this field not supported with addobs @@: ' + field2; addobs=False


        if addsims:
            # 120-yr averages
            sims=('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5')
            #sims=('R1','R2','R3','R4','R5')
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

            fldssdf1,styrs = subsamp_sims(fldssdf1,numyrs=11,styears=styearsE)

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

            fldssdf2,styrs = subsamp_sims(fldssdf2,numyrs=11,styears=styrs)

            simssmm, simssbb, simssrval, simsspval, simssstd_err = sp.stats.linregress(fldssdf1,
                                                                                       fldssdf2)
            print '-------- (subsamp) SIMS ' + simsstr + ' slope, rval, pval: ' + str(simssmm),str(simssrval),str(simsspval)

            # SIC for output file
            if conditional:
                fldssdfcnd = pd.DataFrame(lmd.loaddata((simfieldcnd,),simsss,ncfields=(simncfieldcnd,), timefreq=seacnd, 
                                                     region=region1))*simconvcnd # "simss" subsample sims
                fldssdfcnd,styrs = subsamp_sims(fldssdfcnd,numyrs=11,styears=styrs)



            simsssr=('R1','R2','R3','R4','R5'); simsstr=simsstr+'sbRonly'
            fldssrdf1 = pd.DataFrame(lmd.loaddata((simfield1,),simsssr,ncfields=(simncfield1,), timefreq=sea1, 
                                                 region=region1))*simconv1 # "simss" subsample sims

            if performop1:
                fldssrdf1op = pd.Series(lmd.loaddata((simfield1,),simsssr,ncfields=(simncfield1,), timefreq=sea1, 
                                                     meantype='time',region=region1op))*simconv1

                if op1=='sub': # subtract
                    if region1op=='nh':
                        fldssrdf1 = fldssrdf1 - fldssrdf1op.mean() # want to subtract ens mean nh
                    else:
                        fldssrdf1 = fldssrdf1 - fldssrdf1op # not sure if want ens mean here...
                elif op1=='div': # divide
                    fldssrdf1 = fldssrdf1 / fldssrdf1op

            fldssrdf1,styrsr = subsamp_sims(fldssrdf1,numyrs=11,styears=styearsR)

            fldssrdf2 = pd.DataFrame(lmd.loaddata((simfield2,),simsssr,ncfields=(simncfield2,), timefreq=sea2, 
                                            region=region2))*simconv2

            if performop2:
                fldssrdf2op = pd.Series(lmd.loaddata((simfield2,),simsssr,ncfields=(simncfield2,), timefreq=sea2, 
                                                     meantype='time',region=region2op))*simconv2

                if op2=='sub': # subtract
                    if region2op=='nh':
                        fldssrdf2 = fldssrdf2 - fldssrdf2op.mean() # want to subtract ens mean nh
                    else:
                        fldssrdf2 = fldssrdf2 - fldssrdf2op # not sure if want ens mean here...
                elif op2=='div': # divide
                    fldssrdf2 = fldssrdf2 / fldssrdf2op

            fldssrdf2,styrsr = subsamp_sims(fldssrdf2,numyrs=11,styears=styrsr)

            simssrmm, simssrbb, simssrrval, simssrpval, simssrstd_err = sp.stats.linregress(fldssrdf1,
                                                                                       fldssrdf2)
            print '-------- (subsamp) SIMS sbRonly slope, rval, pval: ' + str(simssrmm),str(simssrrval),str(simssrpval)

            # test correlations b/w R and E sims
            zre, pvre = corrstats.independent_corr(simssrrval, simssrval, len(fldssdf1), n2=len(fldssrdf1))
            print '-------- (subsamp) SIMS R vs E diff b/w correlations zscore, pval : ' + str(zre),str(pvre)

            # test correlations b/w R sims and LE
            zrle, pvrle = corrstats.independent_corr(simssrrval, lerval, len(fldssrdf1), n2=len(lefld1))
            print '-------- (subsamp) SIMS R vs LE diff b/w correlations zscore, pval : ' + str(zrle),str(pvrle)
            # test correlations b/w E sims and LE
            zele, pvele = corrstats.independent_corr(simssrval, lerval, len(fldssdf1), n2=len(lefld1))
            print '-------- (subsamp) SIMS E vs LE diff b/w correlations zscore, pval : ' + str(zele),str(pvele)

            if addpi:
                # test correlations b/w E sims and PI
                zepi, pvepi = corrstats.independent_corr(simssrval, pirval, len(fldssdf1), n2=len(pifld1))
                print '-------- (subsamp) SIMS E vs PI diff b/w correlations zscore, pval : ' + str(zepi),str(pvepi)
                
                # test correlations b/w R sims and PI
                zrpi, pvrpi = corrstats.independent_corr(simssrrval, pirval, len(fldssrdf1), n2=len(pifld1))
                print '-------- (subsamp) SIMS R vs PI diff b/w correlations zscore, pval : ' + str(zrpi),str(pvrpi)

                # test correlations b/w LE sims and PI
                zlepi, pvlepi = corrstats.independent_corr(lerval, pirval, len(lefld1), n2=len(pifld1))
                print '-------- (subsamp) SIMS LE vs PI diff b/w correlations zscore, pval : ' + str(zlepi),str(pvlepi)


            if conditional:
                fldssrdfcnd = pd.DataFrame(lmd.loaddata((simfieldcnd,),simsssr,ncfields=(simncfieldcnd,), timefreq=seacnd, 
                                                     region=region1))*simconvcnd # "simss" subsample sims
                fldssrdfcnd,styrsr = subsamp_sims(fldssrdfcnd,numyrs=11,styears=styrsr)


            simssso=('NSIDC',); simsstr=simsstr+'sbNSIDC'
            fldssodf1 = pd.DataFrame(lmd.loaddata((simfield1,),simssso,ncfields=(simncfield1,), timefreq=sea1, 
                                                 region=region1))*simconv1 # "simss" subsample sims

            if performop1:
                fldssodf1op = pd.Series(lmd.loaddata((simfield1,),simssso,ncfields=(simncfield1,), timefreq=sea1, 
                                                     meantype='time',region=region1op))*simconv1

                if op1=='sub': # subtract
                    if region1op=='nh':
                        fldssodf1 = fldssodf1 - fldssodf1op.mean() # want to subtract ens mean nh
                    else:
                        #print '========== NSIDC SIC: ' + str(fldssodf1.values.mean()) + '-' + str(fldssodf1op.values.mean()) +\
                        #    ' = ' + str(fldssodf1.values.mean() - fldssodf1op.values.mean())
                        fldssodf1 = fldssodf1 - fldssodf1op # not sure if want ens mean here...
                elif op1=='div': # divide
                    fldssodf1 = fldssodf1 / fldssodf1op

            fldssodf1,styrso = subsamp_sims(fldssodf1,numyrs=11)
            

            fldssodf2 = pd.DataFrame(lmd.loaddata((simfield2,),simssso,ncfields=(simncfield2,), timefreq=sea2, 
                                            region=region2))*simconv2

            if performop2:
                fldssodf2op = pd.Series(lmd.loaddata((simfield2,),simssso,ncfields=(simncfield2,), timefreq=sea2, 
                                                     meantype='time',region=region2op))*simconv2

                if op2=='sub': # subtract
                    if region2op=='nh':
                        fldssodf2 = fldssodf2 - fldssodf2op.mean() # want to subtract ens mean nh
                    else:
                        fldssodf2 = fldssodf2 - fldssodf2op # not sure if want ens mean here...
                elif op2=='div': # divide
                    fldssodf2 = fldssodf2 / fldssodf2op

            fldssodf2,styrso = subsamp_sims(fldssodf2,numyrs=11,styears=styearsN)

            simssomm, simssobb, simssorval, simssopval, simssostd_err = sp.stats.linregress(fldssodf1,
                                                                                       fldssodf2)
            print '-------- (subsamp) SIMS NSIDC slope, rval, pval: ' + str(simssomm),str(simssorval),str(simssopval)

            if conditional:
                fldssodfcnd = pd.DataFrame(lmd.loaddata((simfieldcnd,),simssso,ncfields=(simncfieldcnd,), timefreq=seacnd, 
                                                     region=region1))*simconvcnd # "simss" subsample sims
                fldssodfcnd,styrso = subsamp_sims(fldssodfcnd,numyrs=11,styears=styrso)



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
        fig.set_size_inches(10,5.5)
        if conditional:
            #cmap=plt.cm.afmhot #autumn # Greens_r #OrRd_r # coolwarm_r
            #markerdt={'marker': 'o', 'linewidth':0}
            #scaled_cnd = (lefldcnd - lefldcnd.min()) / lefldcnd.ptp()
            #cndcolors = plt.cm.afmhot(scaled_cnd,alpha=0.7) #edgecolors=cndcolors

            if fieldcnd=='sic':
                b2r20cm = plt.cm.get_cmap(cmapcnd) #'blue2red_20')
                testclrs = np.flipud(b2r20cm.colors[[6,11,14,-3,-1],:]) 
                clabcnd='$\Delta$ BKS SIC (%)'
            elif fieldcnd=='tas':
                b2r20cm = plt.cm.get_cmap(cmapcnd)
                testclrs = b2r20cm.colors[[3,6,11,14,-3],:]
                clabcnd='$\Delta$ Eur SAT ($^\circ$C)'

            #testcm = plt.cm.autumn.from_list('testcm',plt.cm.autumn(scaled_cnd,alpha=0.7),N=4)
            #testclrs = np.flipud(b2r20cm.colors[[7,11,14,-1],:]) # 4 color clims -9 to 3
            ### testclrs = np.flipud(b2r20cm.colors[[6,11,14,-3,-1],:]) 
            #testclrs = np.flipud(b2r20cm.colors[[6,11,11,-3,-1],:]) 
            #testclrs = np.flipud(b2r20cm.colors[[6,11,13,13,-1],:]) # 5 color clims -12 to 3
            #testcm = b2r20cm.from_list('testcm',testclrs,N=4)
            testcm = mcolors.ListedColormap(testclrs, name='testcm')

            ledatlg=mlines.Line2D([],[],color=testclrs[2],
                                  markerfacecolor=testclrs[2],marker='o',
                                  linestyle='none',mew=0,markersize=8)

            ledat=ax.scatter(lefld1,lefld2,c=lefldcnd,cmap=testcm,vmin=cmincnd,
                             vmax=cmaxcnd,s=60,alpha=0.9,marker='o', linewidth=0)
            ledat.set_clim([cmincnd,cmaxcnd])

            cbar_ax = fig.add_axes([.91,.25, .02,.5]) 
            cb=fig.colorbar(ledat,cax=cbar_ax)#,orientation='horizontal')
            cb.set_label(clabcnd)
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
        if addmisc:
            ledatm=ax.scatter(lefld1m,lefld2m,color=ccm.get_linecolor('orange3'),marker='*',s=5**2,alpha=0.5)
            axylims = ax.get_ylim()
            axxlims = ax.get_xlim()
            onex=np.linspace(axxlims[0],axxlims[1])
            ax.plot(onex,lemmm*onex + lebbm, color=ccm.get_linecolor('orange3'),linewidth=1)
            leghnds=leghnds + (ledatm,)
            legstrs=legstrs + (casename3 + ' LE',)
            print '-------- LEMisc slope, rval, pval: ' + str(lemmm),str(lervalm),str(lepvalm)   

        if addpi:
            pidathndl=ax.scatter(pifld1,pifld2,color='purple',marker='*',s=5**2,alpha=0.5)
            axylims = ax.get_ylim()
            axxlims = ax.get_xlim()
            onex=np.linspace(axxlims[0],axxlims[1])
            ax.plot(onex,pimm*onex + pibb, color='purple',linewidth=1)
            leghnds=leghnds + (pidathndl,)
            legstrs=legstrs + ('PreIndustrial',)
            print '-------- Preindustrial slope, rval, pval: ' + str(pimm),str(pirval),str(pipval)   
        
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
            simssh = ax.scatter(fldssdf1,fldssdf2,color='0.7',marker='o',s=40,alpha=0.7)
            onex=np.linspace(axxlims[0],axxlims[1])
            ax.plot(onex,simssmm*onex + simssbb, color='0.7',linewidth=1)
            # R sims too
            ax.scatter(fldssrdf1,fldssrdf2,color='cyan',marker='o',s=40,alpha=0.7)
            onex=np.linspace(axxlims[0],axxlims[1])
            ax.plot(onex,simssrmm*onex + simssrbb, color='cyan',linewidth=1)

            leghnds=leghnds + (simssh,)
            legstrs=legstrs + ('11-yr Modelled SIC forcing',)
            obstr=obstr+'sims'

            # NSIDC sim too
            simssho = ax.scatter(fldssodf1,fldssodf2,color='g',marker='o',s=50,alpha=0.7)
            onex=np.linspace(axxlims[0],axxlims[1])
            ax.plot(onex,simssomm*onex + simssobb, color='g',linewidth=1)

            leghnds=leghnds + (simssho,)
            legstrs=legstrs + ('11-yr NSIDC SIC forcing',)
            obstr=obstr+'nsidcsim'

        fontP = fm.FontProperties()
        fontP.set_size('small')
        ax.legend(leghnds,legstrs,
                   loc='best',fancybox=True,framealpha=0.5,prop=fontP,frameon=False) 
        if field2 in ('tas','turb','zg50000.00') and performop2==False:
            ax.set_ylim((ymin,ymax))
            
        """ax.annotate(casename + ' R= ' + '$%.2f$'%(lerval) + ', p='+ '$%.2f$'%(lepval),
                    xy=(axxlims[0]+.05*axxlims[1], axylims[1]-.1*axylims[1]),
                    xycoords='data')
        if addnat:
            plt.annotate(casename2 + ' R= ' + '$%.2f$'%(lervaln) + ', p='+ '$%.2f$'%(lepvaln),
                        xy=(axxlims[0]+.05*axxlims[1], axylims[1]-.2*axylims[1]),
                        xycoords='data')"""

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
            now = str(datetime.datetime.now().time())
            if conditional:
                cnd='cnd'
            else: cnd=''
            if addnat:
                plt.savefig('scatterregress_' + field1 + region1 + opstr1 + str(sea1) + '_v_' + 
                            field2 + region2 +opstr2 + str(sea2) + '_' + casename + '_' + 
                            casename2 + '_2002-12_1979-89' + obstr + cnd + now+ '.pdf')
            else:
                plt.savefig('scatterregress_' + field1 + region1 + opstr1 + str(sea1) + '_v_' + 
                            field2 + region2 +opstr2 + str(sea2) + '_' + casename + 
                            '_2002-12_1979-89' + obstr + cnd + now+ '.pdf')


        if savedat:

            # plot the sic data:
            plt.figure()
            plt.plot(fldssrdfcnd,color='cyan',linestyle='none',marker='o')
            plt.plot(fldssdfcnd,color='0.5',linestyle='none',marker='o')
            plt.plot(fldssodfcnd,color='g',linestyle='none',marker='o')
            plt.plot(lefldcnd,color='blue',linestyle='none',marker='o',alpha=0.5)
            plt.title('SIC data')
            plt.legend(('AGCM var','AGCM mean','AGCM NSIDC','Hist LE'),loc='best')

            # write ascii file for john:
            print 'OBS DATA: ' + str(obsreg1),str(obsreg2)
            agcmout1=np.hstack((fldssrdfcnd,fldssodfcnd,fldssdfcnd)) # column 1
            agcmout2=np.hstack((fldssrdf2,fldssodf2,fldssdf2)) # column 2
            agcmout3=np.hstack((fldssrdf1,fldssodf1,fldssdf1)) # column 3
            agcmout = np.vstack((agcmout1,agcmout2,agcmout3)).T
            np.savetxt('agcmout3.ascii',agcmout,delimiter='\t')

            cgcmout1=np.hstack((lefldcnd,lefldncnd,lefldmcnd))
            cgcmout2=np.hstack((lefld2,lefld2n,lefld2m))
            cgcmout3=np.hstack((lefld1,lefld1n,lefld1m))
            cgcmout = np.vstack((cgcmout1,cgcmout2,cgcmout3)).T
            np.savetxt('cgcmout3.ascii',cgcmout,delimiter='\t')

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
        ensmean=False

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
        bkssat=slope.reshape((nlat,nlon)) # SAT regress on SIC (fieldr actually)
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




        """fig,axs=plt.subplots(1,2)
        fig.set_size_inches(10,5)
        ax=axs[0]#
        cplt.kemmap(bkszg,lat,lon,ptype='nheur',axis=ax,cmin=cmin,cmax=cmax,
                    title=seasp + ' ' + fieldsp + ' regressed onto ' + sear + ' ' +fieldr+regionr )

        ax=axs[1] #
        cplt.kemmap(eurzg,lat,lon,ptype='nheur',axis=ax, cmin=cmin2,cmax=cmax2,
                    title= seasp + ' ' + fieldsp + ' regressed onto ' + sear + ' ' +fieldr2+regionr2)

        if printtofile:
            fig.savefig(fieldr +regionr + '_' + fieldsp + seasp + \
                        + '_regresson_' + fieldr2 + regionr2 + sear + normstr + '.pdf') """

        # @@ create a figure with regression contours on top of other regression:
        #   e.g. z500 regressed onto BKS SIC contours on SAT regressed onto BKS SIC map

        #printtofile=False
        lons, lats = np.meshgrid(lon,lat)
        cmlen=15.
        incr = (cmaxsp2-cminsp2) / (cmlen)
        conts = np.arange(cminsp2,cmaxsp2+incr,incr)


        fig,ax=plt.subplots(1,1)
        fig.set_size_inches(5,5)
        bm,pc=cplt.kemmap(bkssat,lat,lon,ptype='nheur',axis=ax,cmin=cmin,cmax=cmax,
                    title=seasp + ' regressions onto ' + sear + ' ' + fieldr+regionr)
        bm.contour(lons,lats,bkszg,levels=conts,
                   colors='k',linewidths=1,latlon=True)
        if printtofile:
            fig.savefig(fieldsp + seasp + '_' + fieldsp2 + seasp \
                        + '_regresson_' + fieldr + regionr + sear + normstr + '.pdf') 


        """tmp=np.zeros(bkssat.shape)
        fig,ax=plt.subplots(1,1)
        fig.set_size_inches(5,5)
        bm,pc=cplt.kemmap(tmp,lat,lon,ptype='nheur',axis=ax,cmin=cmin,cmax=cmax,
                          title=seasp + ' regressions onto ' + sear + ' ' + fieldr+regionr,suppcb=True)

        bm.contour(lons,lats,bkszg,levels=conts,
                   colors='k',linewidths=1,latlon=True)
        if printtofile:
            fig.savefig(fieldsp2 + seasp \
                        + '_regresson_' + fieldr + regionr + sear + normstr + '.pdf') """


        return bkssat,bkszg


    if dolongtermavg:

        extra=False # timeseries, testing kernel density estimates

        longtermLE=False # else, subsample AGCM sims
        #addsims=False # 120-yr avg AGCM sims
        #addnat=True
        #addmisc=True
        appendorig=False # append original historical runs? (for now 6/22/15)
        addorig=True # just add orig points to distribution curve (not histogram)
        #addpi=True # add preindustrial control

        #printtofile=True

        numsamp=100 # how many times to sample 11 ens members, for longtermLE=True
        addraw=False # add the decadal diffs?  for longtermLE=True
        option2=True # otherwise, original method of selecting 11 members numsamp times. (longtermLE=True)
        if option2:
            numsamp=50

        subnh=False # @@@@ when subtracting, prob should subtract the mean NH temp, not individual runs
        substr='' # for filename
        prstr=''

        sea='DJF'
        #siglevel=0.1
        #cisiglevel=0.05

        leconv=1
        field='tas'
        ncfield='tas'
        comp='Amon'
        region='eurasiamori'
        #region='eurasiathicke'
        #region='gt60n'

        #field='sia'; ncfield='sianh'; comp='OImon'; region='nh';

        if region=='eurasiamori' or region=='eurasiathicke': ttlstr='Eurasian'
        elif region=='gm': ttlstr='Global'
        elif region=='gt60': ttlstr='Polar (>60N)'
        else: ttlstr=''

        mew=1.5; ms=7
        firebrick=ccm.get_linecolor('firebrick')
        hcol='0.7'#ccm.get_linecolor('darkolivegreen3')
        hcolline='0.7'#ccm.get_linecolor('darkolivegreen3')#'darkseagreen4')

        ltcol='0.5'#ccm.get_linecolor('darkseagreen4')

        natcol=ccm.get_linecolor('steelblue3')
        natcolline=ccm.get_linecolor('steelblue3')#4')
        miscol=ccm.get_linecolor('orange3') #'deepskyblue')
        miscolline=ccm.get_linecolor('orange3') #'deepskyblue')
        obscol='k' #ccm.get_linecolor('midnightblue') #'b'

        lat=le.get_lat(local=local)
        lon=le.get_lon(local=local)
        nlat=len(lat); nlon=len(lon)

        fdict = {'field': field+region, 'ncfield': ncfield, 'comp': comp}

        # ===== put in function =====
        """simconv1=1
        if field=='tas': simfield1='st'; simncfield1='ST'
        elif field=='zg50000.00': simfield1='gz50000'; simncfield1='PHI'; simconv1=1/con.get_g()
        elif field=='sia': simfield1='sicn'; simncfield1='SICN'; print '@@ danger, sia actually sicn average'
        else: print 'cannot addsims for ' + field;

        casename = 'historical'


        # LOAD 1D DATA
        lecdat = le.load_LEdata(fdict,casename,timesel=timeselc, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
        (numen,ntime) = lecdat.shape
        lepdat=le.load_LEdata(fdict,casename,timesel=timeselp, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
        lecsea = cutl.seasonalize_monthlyts(lecdat.T,season=sea,verb=True).mean(axis=0)
        lepsea = cutl.seasonalize_monthlyts(lepdat.T,season=sea).mean(axis=0)
        lesea = lepsea - lecsea # numens

        # load original five in case:
        or5c=le.load_originalfive(fdict, casename,timesel=timeselc, rettype='ndarray',conv=leconv,ftype=ftype)
        or5p=le.load_originalfive(fdict, casename,timesel=timeselp, rettype='ndarray',conv=leconv,ftype=ftype)
        or5csea=cutl.seasonalize_monthlyts(or5c.T,season=sea,verb=True).mean(axis=0)
        or5psea=cutl.seasonalize_monthlyts(or5p.T,season=sea,verb=True).mean(axis=0)
        or5sea = or5psea-or5csea

        if subnh:
            fdictsub = {'field': field+'nh', 'ncfield': ncfield, 'comp': comp}
            subc = le.load_LEdata(fdictsub,casename,timesel=timeselc, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
            subp = le.load_LEdata(fdictsub,casename,timesel=timeselp, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
            sub = cutl.seasonalize_monthlyts(subp.T,season=sea).T - cutl.seasonalize_monthlyts(subc.T,season=sea).T
            # @@@ change this to subtract the MEAN hemispheric anom from all ens members
            #lesea = lesea - sub.mean(axis=1)
            lesea = lesea - sub.mean(axis=1).mean()

            # load original five in case:
            subor5c=le.load_originalfive(fdictsub, casename,timesel=timeselc, rettype='ndarray',conv=leconv,ftype=ftype)
            subor5p=le.load_originalfive(fdictsub, casename,timesel=timeselp, rettype='ndarray',conv=leconv,ftype=ftype)
            subor5=cutl.seasonalize_monthlyts(or5csub.T,season=sea,verb=True).T - \
                        cutl.seasonalize_monthlyts(or5psub.T,season=sea,verb=True).T
            or5sea = or5sea - subor5.mean(axis=1).mean() # @@@@ should remove the total mean of all 55 in this case?

        if appendorig:
            lesea = np.hstack((lesea,or5sea))

        # RAW anomalies (decadal diffs)
        # calc the pdf associated with the hist
        rawpdf_fitted,rawmean,rawsd,rawxx = cutl.calc_normfit(lesea) # to get mean
        #rawpdf_fitted,rawxx = cutl.calc_kernel(lesea)
        rawdf = len(lesea)-1
        rawstder = rawsd / np.sqrt(rawdf+1)
        rawci = sp.stats.t.interval(1-cisiglevel, rawdf, loc=rawmean, scale=rawstder)
        rawcif = sp.stats.t.interval(1-cisiglevel, rawdf, loc=rawmean, scale=rawsd)

        if addnat:
            lencdat = le.load_LEdata(fdict,'historicalNat',timesel=timeselc, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
            (numenn,ntimen) = lencdat.shape
            lenpdat=le.load_LEdata(fdict,'historicalNat',timesel=timeselp, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
            lencsea = cutl.seasonalize_monthlyts(lencdat.T,season=sea).mean(axis=0)
            lenpsea = cutl.seasonalize_monthlyts(lenpdat.T,season=sea).mean(axis=0)
            lensea = lenpsea - lencsea # numens

            if subnh:
                fdictsub = {'field': field+'nh', 'ncfield': ncfield, 'comp': comp}
                nsubc = le.load_LEdata(fdictsub,'historicalNat',timesel=timeselc, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
                nsubp = le.load_LEdata(fdictsub,'historicalNat',timesel=timeselp, rettype='ndarray',conv=leconv,ftype=ftype,local=local)
                nsub = cutl.seasonalize_monthlyts(nsubp.T,season=sea).T - cutl.seasonalize_monthlyts(nsubc.T,season=sea).T
                # @@@ change this to subtract the MEAN hemispheric anom from all ens members
                #lesea = lesea - sub.mean(axis=1)
                lensea = lensea - nsub.mean(axis=1).mean()

            rawnpdf_fitted,rawnmean,rawnsd,rawnxx = cutl.calc_normfit(lensea)
            #rawnpdf_fitted,rawnxx = cutl.calc_kernel(lensea)
            rawndf = len(lensea)-1
            rawnstder = rawnsd / np.sqrt(rawndf+1)
            rawnci = sp.stats.t.interval(1-cisiglevel, rawndf, loc=rawnmean, scale=rawnstder)
            rawncif = sp.stats.t.interval(1-cisiglevel, rawndf, loc=rawnmean, scale=rawnsd)

        if addmisc:
            lemcdat = le.load_LEdata(fdict,'historicalMisc',timesel=timeselc, 
                                     rettype='ndarray',conv=leconv,ftype=ftype,local=local)
            (numenm,ntimem) = lemcdat.shape
            lempdat=le.load_LEdata(fdict,'historicalMisc',timesel=timeselp, 
                                   rettype='ndarray',conv=leconv,ftype=ftype,local=local)
            lemcsea = cutl.seasonalize_monthlyts(lemcdat.T,season=sea).mean(axis=0)
            lempsea = cutl.seasonalize_monthlyts(lempdat.T,season=sea).mean(axis=0)
            lemsea = lempsea - lemcsea # numens

            if subnh:
                fdictsub = {'field': field+'nh', 'ncfield': ncfield, 'comp': comp}
                msubc = le.load_LEdata(fdictsub,'historicalMisc',timesel=timeselc, 
                                       rettype='ndarray',conv=leconv,ftype=ftype,local=local)
                msubp = le.load_LEdata(fdictsub,'historicalMisc',timesel=timeselp, 
                                       rettype='ndarray',conv=leconv,ftype=ftype,local=local)
                msub = cutl.seasonalize_monthlyts(msubp.T,season=sea).T - cutl.seasonalize_monthlyts(msubc.T,season=sea).T
                # @@@ change this to subtract the MEAN hemispheric anom from all ens members
                #lesea = lesea - sub.mean(axis=1)
                lemsea = lemsea - msub.mean(axis=1).mean()

            rawmpdf_fitted,rawmmean,rawmsd,rawmxx = cutl.calc_normfit(lemsea)
            #rawmpdf_fitted,rawmxx = cutl.calc_kernel(lemsea)
            rawmdf = len(lemsea)-1
            rawmstder = rawmsd / np.sqrt(rawmdf+1)
            rawmci = sp.stats.t.interval(1-cisiglevel, rawmdf, loc=rawmmean, scale=rawmstder)
            rawmcif = sp.stats.t.interval(1-cisiglevel, rawmdf, loc=rawmmean, scale=rawmsd)

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

            ltpdf_fitted,ltmean,ltsd,ltxx = cutl.calc_normfit(ltavg)


        else: # not longtermLE (short term sims instead)

            # here we want to subsample 11-year segments from the sims
            sims=('E1','E2','E3','E4','E5'); simsstr='sbEonly' # sub
            #sims=('R1','R2','R3','R4','R5'); simsstr='sbRonly'
            #sims=('R1','R2','R3','R4','R5','E1','E2','E3','E4','E5'); simsstr='sbERboth'
            simflddf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                                 region=region))*simconv1
            osimflddf = pd.DataFrame(lmd.loaddata((simfield1,),('NSIDC',),ncfields=(simncfield1,), timefreq=sea, 
                                                 region=region))*simconv1
            if subnh:
                simsubdf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                                     filetype='diff',region='nh'))*simconv1
                # want to subtract the mean hemispheric avg anomaly
                simflddf = simflddf - simsubdf.mean(axis=1).mean() # average over sims and then time (should be a scalar)
                osimsubdf = pd.DataFrame(lmd.loaddata((simfield1,),('NSIDC',),ncfields=(simncfield1,), timefreq=sea, 
                                                     filetype='diff',region='nh'))*simconv1
                # want to subtract the mean hemispheric avg anomaly
                osimflddf = osimflddf - osimsubdf.mean(axis=1).mean() # average over sims and then time (should be a scalar)

            subsamp,styearsss = subsamp_sims(simflddf,numyrs=11)
            plotsims = subsamp
            print '==== AGCM '
            sspdf_fitted,ssmean,sssd,ssxx = cutl.calc_normfit(plotsims)
            #sspdf_fitted,ssxx = cutl.calc_kernel(plotsims)
            ssdf = len(subsamp)-1
            ssstder = sssd / np.sqrt(ssdf+1)

            # mean and conf-int: @@@@@
            # ci = sp.stats.t.interval(1-cisiglevel, df, loc= meananom, scale=stder)
            ssci = sp.stats.t.interval(1-cisiglevel, ssdf, loc=ssmean, scale=ssstder) # this is mean value's 5-95%
            sscif = sp.stats.t.interval(1-cisiglevel, ssdf, loc=ssmean, scale=sssd) # full range. don't divide by n. this is pdf 5-95%

            osubsamp,styearsns = subsamp_sims(osimflddf,numyrs=11)
            print '==== NSIDC AGCM '
            osspdf_fitted,ossmean,osssd,ossxx = cutl.calc_normfit(osubsamp)
            ossdf = len(osubsamp)-1
            ossstder = osssd / np.sqrt(ossdf+1)
            ossci = sp.stats.t.interval(1-cisiglevel, ossdf, loc=ossmean, scale=ossstder)
            osscif = sp.stats.t.interval(1-cisiglevel, ossdf, loc=ossmean, scale=osssd)

            # do each ens separately, for testing stats (maybe not plotting)
            sims2=('R1','R2','R3','R4','R5'); #simsstr='sbRonly'
            sim2flddf = pd.DataFrame(lmd.loaddata((simfield1,),sims2,ncfields=(simncfield1,), timefreq=sea, 
                                                 region=region))*simconv1
            if subnh:
                sim2subdf = pd.DataFrame(lmd.loaddata((simfield1,),sims2,ncfields=(simncfield1,), timefreq=sea, 
                                                     filetype='diff',region='nh'))*simconv1

                # want to subtract the mean hemispheric avg anomaly
                sim2flddf = sim2flddf - sim2subdf.mean(axis=1).mean() # average over sims and then time (should be a scalar)


            subsamp2,styears2 = subsamp_sims(sim2flddf,numyrs=11)
            plotsims2 = subsamp2
            print '==== AGCM2 '
            ss2pdf_fitted,ss2mean,ss2sd,ss2xx = cutl.calc_normfit(plotsims2)
            #ss2pdf_fitted,ss2xx = cutl.calc_kernel(plotsims2)
            ss2df = len(subsamp2)-1
            ss2stder = ss2sd / np.sqrt(ss2df+1)
            ss2ci = sp.stats.t.interval(1-cisiglevel, ss2df, loc=ss2mean, scale=ss2stder)
            ss2cif = sp.stats.t.interval(1-cisiglevel, ss2df, loc=ss2mean, scale=ss2sd)

            tstat, pval = sp.stats.ttest_ind(plotsims,plotsims2)
            lstat, lpval = sp.stats.levene(plotsims,plotsims2)
            print '==== testing ANT vs TOT ===='
            print 'TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
            if pval<=siglevel:
                 print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'
            print 'LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
            if lpval<=siglevel:
                 print 'The ensemble variances are significantly different (' + str(1-siglevel) + ')'


            tstat, pval = sp.stats.ttest_ind(plotsims,lesea)
            lstat, lpval = sp.stats.levene(plotsims,lesea)
            print '==== testing ANT vs ALL LE ===='
            print 'TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
            if pval<=siglevel:
                 print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'
            print 'LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
            if lpval<=siglevel:
                 print 'The ensemble variances are significantly different (' + str(1-siglevel) + ')'

            tstat, pval = sp.stats.ttest_ind(lesea,lensea)
            lstat, lpval = sp.stats.levene(lesea,lensea)
            print '==== testing LE vs Nat LE ===='
            print 'TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
            if pval<=siglevel:
                 print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'
            print 'LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
            if lpval<=siglevel:
                 print 'The ensemble variances are significantly different (' + str(1-siglevel) + ')'

            tstat, pval = sp.stats.ttest_ind(lesea,lemsea)
            lstat, lpval = sp.stats.levene(lesea,lemsea)
            print '==== testing LE vs Aero LE ===='
            print 'TSTAT: ' + str(tstat) + ' PVAL: ' + str(pval)
            if pval<=siglevel:
                 print 'The ensemble means are significantly different (' + str(1-siglevel) + ')'
            print 'LSTAT: ' + str(lstat) + ' PVAL: ' + str(lpval)
            if lpval<=siglevel:
                 print 'The ensemble variances are significantly different (' + str(1-siglevel) + ')'

        if addsims:
            # ========= add sims (120-yr averages)
            sims=('E1','E2','E3','E4','E5','R1','R2','R3','R4','R5','NSIDC'); #simsstr=''
            #simsRN=('R1','R2','R3','R4','R5','NSIDC')
            simsR=('R1','R2','R3','R4','R5')
            simsE=('E1','E2','E3','E4','E5')
            #sims=simsRN; simsstr='Ronly'
            #simflddf = pd.DataFrame(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
            #                                     meantype='time',region=region),index=sims)*simconv1
            simflddf = pd.Series(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                                 meantype='time',region=region))*simconv1

            # --- for estimating sigma in Individual SIC forcing ensemble only:
            #simflddfr = pd.DataFrame(lmd.loaddata((simfield1,),simsR,ncfields=(simncfield1,), timefreq=sea, 
            #                                      meantype='time',region=region),index=simsR)*simconv1
            simflddfr = pd.Series(lmd.loaddata((simfield1,),simsR,ncfields=(simncfield1,), timefreq=sea, 
                                                  meantype='time',region=region))*simconv1
            simval=simflddfr.values
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
                simfldm =simflddf.mask(simpvdf>siglevel)# ma.masked_where(pvals>siglevel,simflddf.values)
                # have to calc pvals 

            else:
                simpvdf = pd.Series(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                                    filetype='pval',region=region))#,index=sims)
                simfldm = simflddf.mask(simpvdf>siglevel)
                #simfldm = ma.masked_where(simpvdf.values>siglevel,simflddf.values)

        else: # just add NSIDC
            sims=('NSIDC',)
            simflddf = pd.Series(lmd.loaddata((simfield1,),sims,ncfields=(simncfield1,), timefreq=sea, 
                                                 meantype='time',region=region))*simconv1
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
                simfldm =simflddf.mask(simpvdf>siglevel)# ma.masked_where(pvals>siglevel,simflddf.values)

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
                ttlstr=ttlstr+'-NH '
        elif field == 'sia':
            nsidcfile = '/HOME/rkm/work/BCs/NSIDC/td_bootstrap_197811_latest_128_64_sicn_1978111600-2013121612.nc'
            latns=cnc.getNCvar(nsidcfile,'lat')
            lonns=cnc.getNCvar(nsidcfile,'lon')
            nsidcsel='1979-01-01,2011-12-31'
            nsidcyrs=np.arange(1979,2012)
            if sea=='DJF':
                nsidcyrs=nsidcyrs[:-1]
            nsidcfldc=cnc.getNCvar(nsidcfile,'SICN',timesel=timeselc,seas=sea)
            nsidcfldp=cnc.getNCvar(nsidcfile,'SICN',timesel=timeselp,seas=sea)
            obsreg=cutl.calc_regtotseaicearea(nsidcfldp,latns,lonns,region='nh',isarea=False).mean() -\
                      cutl.calc_regtotseaicearea(nsidcfldc,latns,lonns,region='nh',isarea=False).mean()"""


        # ========= make the figures =========
        if longtermLE:
            if option2:
                # how does sigma converge with n?
                fig,ax=plt.subplots(1,1)
                ax.hist(ltsigma,color=ltcol,alpha=0.5)
                ax.axvline(simsigma,color='0.3',linewidth=3)
                for ll in np.arange(0,5):
                    ax.axvline(simstds[ll],color='0.3',linewidth=1)
                ax.set_title('$\sigma$ across sets of 4 Estimated 120-year ' + ttlstr + 'SAT change (DJF)')
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
                ax.set_title('$\sigma$ across sets of 4 Estimated 120-year ' + ttlstr + 'SAT change (DJF)')
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
            ax.set_title('Estimated 120-year ' + ttlstr + 'SAT change (' + sea +')')
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

        else: # not longtermLE figs ===============================================

            loadmat=False
            savemat=False
            vertical=False # vertical plot 
            #combagcm=True
            #comblenat=True
            keys=('xxdat', 'histdat', 'meandat', 'pdfdat', 'cidat', 'pointdat', 'cifdat')
            simsstr=''

            #fielda='sia'; ncfielda='sianh'; compa='OImon'; regiona='nh';
            fielda='sic'; ncfielda='sic'; compa='OImon'; regiona='bksmori';
            #fielda='zg50000.00'; ncfielda='zg'; compa='Amon'; regiona='bksmori'

            fdicta = {'field': fielda+regiona, 'ncfield': ncfielda, 'comp': compa}

            if loadmat:
                if fielda=='sic':
                    when='14:51:28.762886sic'; 
                else:
                    when='14:51:28.762886'; 


                #when='18:19:46.394326' #'16:06:52.938245' # # choose which set of files to load
                print 'loadmat!! when=' + when # @@
                retdict = dict.fromkeys(keys)
                retdicta = dict.fromkeys(keys)

                #@@@@ simsstr +\ # remove from filename!
                matbase = 'pymatfiles/' + field + region + substr + '_' + sea + '_LEsims' +\
                          'obs_11yrsubsmpavgraw_histinset' # hack in filename. 

                for key in keys:
                    matname =  matbase + '_' + key +'_' + when + '.mat'
                    print matname
                    retdict[key] = sio.loadmat(matname, squeeze_me=True)

                    matname =  matbase + '_' + key +'a_' + when + '.mat'
                    retdicta[key] = sio.loadmat(matname,squeeze_me=True)


            else:
                retdicta,simsstra = calc_shorttermpdf(fdicta,fielda,regiona,
                                                      sea,(timeselc,timeselp),
                                                      verb=verb,addpi=addpi,
                                                      addnat=addnat,addmisc=addmisc,
                                                      combagcm=combagcm,comblenat=comblenat)

                retdict,simsstr = calc_shorttermpdf(fdict,field,region,
                                                    sea,(timeselc,timeselp),
                                                    verb=verb,addpi=addpi,
                                                    addnat=addnat,addmisc=addmisc,
                                                    combagcm=combagcm,comblenat=comblenat)

            if vertical:
                fig,axs=plt.subplots(2,1)
                fig.set_size_inches(10,18)
                fig.subplots_adjust(hspace=.15)
            else:
                fig,axs=plt.subplots(1,2)
                fig.set_size_inches(18,8)
                fig.subplots_adjust(wspace=.07)

            ax2=axs[0]
            ax2,prstr=plot_shorttermpdf(fig,ax2,fielda,regiona,addpi=addpi,
                                        addnat=addnat,addmisc=addmisc,
                                        combagcm=combagcm,plab='a',pversion=pversion[0],**retdicta)
            ax=axs[1]
            ax,prstr=plot_shorttermpdf(fig,ax,field,region,addpi=addpi,
                                       addnat=addnat,addmisc=addmisc,
                                       combagcm=combagcm,plab='b',pversion=pversion[0],**retdict)

            if loadmat:
                now = when # have filename match the data
            else:
                now = when='14:51:28.762886sic';#@@@@@str(datetime.datetime.now().time())

            if printtofile:
                fig.savefig(field + region + substr + '_' + sea + '_LEsims' + simsstr +\
                            'obs_11yrsubsmpavgraw_histinset' + prstr + now + '_' + pversion[0] +'.pdf')# _b has AGCM hist and not CGCM All@@


            # plot them separately!
            fig,ax=plt.subplots(1,1)
            fig.set_size_inches(9,8)
            ax,prstr=plot_shorttermpdf(fig,ax,fielda,regiona,addpi=addpi,
                                       addnat=addnat,addmisc=addmisc,addanno=False,
                                       combagcm=combagcm,**retdicta)
            if printtofile:
                fig.savefig(fielda + regiona + substr + '_' + sea + '_' +\
                            'obs_11yrsubsmpavgraw_histinset' + prstr + now + '_' + pversion[0] +'.pdf')


            fig,ax=plt.subplots(1,1)
            fig.set_size_inches(9,8)
            ax,prstr=plot_shorttermpdf(fig,ax,field,region,addpi=addpi,
                                       addnat=addnat,addmisc=addmisc,addanno=False,
                                       combagcm=combagcm,pversion=pversion[1],**retdict)
            if printtofile:
                fig.savefig(field + region + substr + '_' + sea + '_' +\
                            'obs_11yrsubsmpavgraw_histinset' + prstr + now + '_' + pversion[1] + '.pdf')


            if savemat:

                matbase = 'pymatfiles/' + field + region + substr + '_' + sea + '_LEsims'+\
                          'obs_11yrsubsmpavgraw_histinset' # + simsstr +\

                for key in retdict.keys():
                    matname =  matbase + '_' + key +'_' + now + '.mat'
                    sio.savemat(matname, retdict[key])

                    matname =  matbase + '_' + key +'a_' + now + '.mat'
                    sio.savemat(matname, retdicta[key])



        if extra: # ===================================================
            # TESTING kernel shape
            kernel=sp.stats.gaussian_kde(lesea)
            kernel2=sp.stats.gaussian_kde(lesea,bw_method='silverman')
            kernelf1=sp.stats.gaussian_kde(lesea,bw_method=1)
            kernelfp7=sp.stats.gaussian_kde(lesea,bw_method=.7)
            kernelfp35=sp.stats.gaussian_kde(lesea,bw_method=.35)
            kernelfp25=sp.stats.gaussian_kde(lesea,bw_method=.25)

            plt.figure()
            plt.plot(rawxx,kernel(rawxx),'b--')
            plt.plot(rawxx,kernel2(rawxx),'r--')
            plt.plot(rawxx,kernelf1(rawxx),'g--')
            plt.plot(rawxx,kernelfp7(rawxx),'c--')
            plt.plot(rawxx,kernelfp35(rawxx),'m--')
            plt.plot(rawxx,kernelfp25(rawxx),color='orange',linestyle='--')
            plt.hist(lesea,normed=True,color=hcol,alpha=0.3,histtype='stepfilled')
            plt.legend(('scott (scale=0.457)','silverman (scale=0.484)',
                        'scale=1','scale=0.7','scale=0.35','scale=0.25','orig'),
                        loc='best')



            # TEST OUT TIMESERIES OF LE VS OBS
            testfld='tas'; testreg='eurasiamori'; testsea='DJF'
            les = ('historical','historicalNat','historicalMisc')
            lecols = {'historical':hcol,'historicalNat':natcol,'historicalMisc':miscol}

            gmdt={'field':testfld+testreg,'ncfield': ncfield,'comp':comp}

            legmanom={}
            for lename in les:

                lebase = le.load_LEdata(gmdt,lename,rettype='ndarray',
                                        conv=leconv,timesel='1951-01-01,1980-12-31',ftype=ftype,local=local)
                #lebase = cutl.annualize_monthlyts(lebase.T).mean(axis=0)
                lebase = cutl.seasonalize_monthlyts(lebase.T,season=testsea).mean(axis=0)
                legm = le.load_LEdata(gmdt,lename,rettype='ndarray',
                                      conv=leconv,ftype=ftype,local=local)
                legm=cutl.seasonalize_monthlyts(legm.T,season=testsea)
                legmanom[lename]=legm-lebase

            gisgm=cnc.getNCvar(gisfile,'tempanomaly',timesel='1950-01-01,2020-12-31',seas=testsea)
            gisgm=cutl.calc_regmean(gisgm,latgis,longis,testreg)
            xxle=np.arange(1950,1950+len(legm))
            xxgis=np.arange(1950,1950+len(gisgm))

            leleg=mlines.Line2D([],[],color='0.5')
            gisleg=mlines.Line2D([],[],color='black',linewidth=2)

            legs=[]
            legsstr=[]
            plt.figure()
            for lename in les:
                legs.append(mlines.Line2D([],[],color=lecols[lename]))
                legsstr.append(lename)
                plt.plot(xxle,legmanom[lename],color=lecols[lename],alpha=0.5)

            legs.append(gisleg)
            legsstr.append('GIStemp1200')

            plt.plot(xxgis,gisgm,color='black',linewidth=2)
            plt.title(testreg + ' ' + testsea +' SAT anom (base 1951-80)')
            plt.xlabel('years')
            plt.legend((leleg,gisleg),('Historical LE','GIStemp1200'),loc='best',frameon=False)
            plt.legend(legs,legsstr,loc='best',frameon=False)
            if printtofile:
                plt.savefig(testfld + testreg + '_anom_1951-80base_HistoricalLENatMisc_GIStemp_timeseries.pdf')

            # Which ens member matches observations the best?
            # Test with Historical LE:
            lehist = legmanom[lename]
            lehist=lehist[:-7,:] # to match obs
            ccoef=np.zeros(lehist.shape[1])
            for ii in np.arange(0,lehist.shape[1]):
                # for each ens, calc corr coef
                ccoef[ii] = np.corrcoef(lehist[:,ii],gisgm)[0,1]

            # @@@@ ideally would like to plot the "max correlation" ens
            # member in the timeseries. Also ideally, if testing eurasiamori,
            # check only since the year 2000 or something.


            # OR pick out the 5 coldest Eurasia's and 5 warmest: what is the diff?


        """# ==== estimate sigma uncertainty ====
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
        # make sia: @@@@@@@@@@@@@@@@@@@@@"""
        """# for estimating sigma:
        simflddfr = pd.DataFrame(lmd.loaddata((simfield1,),simsR,ncfields=(simncfield1,), timefreq=sea, 
                                              meantype='time',region=region),index=simsR)*simconv1
        simval=simflddfr.values[0,:]
        simstds = [simval[1:].std(), simval[[0,2,3,4]].std(), simval[[0, 1, 3, 4]].std(), 
                   simval[[0, 1, 2, 4]].std(), simval[0:-1].std()]

        simsigma=simval.std()"""

# ===== end main() =======


# this means, if we are running this module as main script
#    (not importing from another), register and show
if __name__ == "__main__":
    main(dowhat='doscatter',addobs=True,addnat=True,
         addmisc=True,addsims=True,addpi=True) 
