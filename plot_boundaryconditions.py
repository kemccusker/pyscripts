import cccmautils as cutl
import constants as con
import loadmodeldata as lmd
import matplotlib.lines as mlines

plt.close('all')
masksicn=True # use present day sicn to mask SST
printtofile=False

deni = 913. # kg/m3 density of ice
timesel='0011-01-01,0011-12-31' # don't have to load much b/c constant
#djf=[0,1,11]
sea='DJF'
fields=('sicn','gt','sic')


def load_BCs(fields, season=None, timesel='0011-01-01,0011-12-31'):

    simcbase='kemctl1'
    simpbase='kem1pert2'
    obsbase = 'kemnsidc'

    flddat={}; fldcdat={}; fldpdat={}; sicnpdat={}
    for fii,field in enumerate(fields):
        fnamesc={}; fnamesp={}
        flddiff={}; fldcavg={}; fldpavg={}

        print field
        # Loop through sim data:
        # load all 5 sim BCs
        for ii in range(0,6):

            # do obs data at ii=5
            if ii==5:
                simname=obsbase+'ctl'
                fnamec = con.getBCfilenames(field,sim=simname)
                print fnamec
                simname=obsbase+'pert'
                fnamep = con.getBCfilenames(field,sim=simname)
            else:
                simname=simcbase+'r'+str(ii+1)
                fnamec = con.getBCfilenames(field,sim=simname)
                print fnamec
                simname=simpbase+'r'+str(ii+1)
                fnamep = con.getBCfilenames(field,sim=simname)

            fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)
            fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)

            if field=='sicn':
                # save pert data to mask gt
                if season==None:
                    sicnpdat[ii] = fldp
                else:
                    sicnpdat[ii] = cutl.seasonalize_monthlyts(fldp,season=sea,climo=1)

            if season==None:
                flddiff[ii] = fldp-fldc
                fldcavg[ii] = fldc
                fldpavg[ii] = fldp
            else:
                flddiff[ii] = cutl.seasonalize_monthlyts(fldp-fldc,season=sea,climo=1) 
                fldcavg[ii] = cutl.seasonalize_monthlyts(fldc,season=sea,climo=1) 
                fldpavg[ii] = cutl.seasonalize_monthlyts(fldp,season=sea,climo=1) 

            fnamesc[ii]=fnamec
            fnamesp[ii]=fnamep

        flddat[field] = flddiff
        fldcdat[field] = fldcavg
        fldpdat[field] = fldpavg

    return flddat, fldcdat, fldpdat, sicnpdat



# === info for plot:

gt={}; sicn={}; sic={}
gt['cmin'] = -2.; gt['cmax'] = 2.; gt['cmap']='blue2red_w20'
#sicn['cmin'] = -30; sicn['cmax'] = 30; sicn['cmap']='red2blue_w24sic'
sicn['cmin'] = -25; sicn['cmax'] = 25; sicn['cmap']='red2blue_w20sic'
sic['cmin'] = -0.5; sic['cmax'] = 0.5; sic['cmap']='red2blue_w20'

meta = {'gt': gt, 'sicn': sicn, 'sic': sic}
ylabs = {0: 'OBS', 1: 'E1', 2: 'E2', 3: 'E3', 4: 'E4', 5:'E5'} 
convdt = {'gt': 1, 'sicn': 100, 'sic': 1/np.float(deni)}


# === load one year for seasonal cycle:
#flddat,fldcdat,fldpdat,sicnpdat = load_BCs(fields, season=None, timesel=timesel)




# ============== PLOT ==============
printtofile=True
lat = con.get_t63lat() #cnc.getNCvar(fnamec,'lat')
lon = con.get_t63lon() #cnc.getNCvar(fnamep,'lon')

fsz=16
fliptohoriz=True
latlim=57
pparams = {'lcol':'0.9', 'coastres': 'c', 'coastwidth': 0.5, 
           'area_thresh':70000,'lmask':True}

fields=('sicn','sic','gt')
fields=('sicn',)
tlabs = {'gt': '$\Delta$ SST ($^\circ$C)', 'sicn': '$\Delta$ SIC (%)', 'sic': '$\Delta$ SIT (m)'}

if fliptohoriz:
    psuff='_horizor'+str(len(fields))+'b' # horizontal orientation
    # rows are fields and cols are sims
    if len(fields)==1:
        fig,axs = plt.subplots(2,3)
        fig.set_size_inches((6,len(fields)*2+2))
        fig.subplots_adjust(wspace=0.01,hspace=0.3)
    else:
        fig,axs = plt.subplots(len(fields),6)
        fig.set_size_inches((10,len(fields)*2))
        fig.subplots_adjust(wspace=0.05,hspace=0.1)
else:
    psuff='_vertor'+str(len(fields))+'b' # vertical orientation
    # cols are fields and rows are sims
    if len(fields)==1:
        fig,axs = plt.subplots(3,2)
        fig.set_size_inches((len(fields)*2+1,6))
        fig.subplots_adjust(wspace=0.05,hspace=0.05)
    else:
        fig,axs = plt.subplots(6,len(fields))
        fig.set_size_inches((len(fields)*2+1,14))
        fig.subplots_adjust(wspace=0.1,hspace=0.05)

for fii,field in enumerate(fields):

    metap = meta[field]
    print field, metap
    conv = convdt[field]

    flddiff=flddat[field]

    for ii in range(0,6): # 6 sims

        if len(fields)==1:
            ax=axs.flatten()[ii]

        else:
            if fliptohoriz:
                ax=axs[fii,ii]
            else:
                ax=axs[ii,fii]

        plotfld=flddiff[ii]*conv
        plotfld = cutl.seasonalize_monthlyts(plotfld,season=sea, climo=1)

        if field=='gt' and masksicn:
            tmp = sicnpdat[ii]
            tmp = cutl.seasonalize_monthlyts(tmp,season=sea, climo=1)
            plotfld[tmp>=0.15] = 0 # mask where model sees grid cell as sea ice (not SST)

        bm,pc = cplt.kemmap(plotfld, lat, lon, ptype='nheur',
                            latlim=latlim, axis=ax, suppcb=True, 
                            round=False,lcol=pparams['lcol'], coastres= 'c', 
                            coastwidth= 0.5, area_thresh=70000,
                            lmask=True,**metap)

        if len(fields)==1:
            ax.set_title(ylabs[ii],fontsize=fsz)
        else:
            if fliptohoriz:
                if ax.is_first_col(): ax.set_ylabel(tlabs[field])
                if ax.is_first_row(): ax.set_title(ylabs[ii])
            else:
                if ax.is_first_row(): ax.set_title(tlabs[field])
                if ax.is_first_col(): ax.set_ylabel(ylabs[ii])

    axposx = ax.get_position().get_points()[0][0] # x position of last ax
    axposy = ax.get_position().get_points()[0][1] # y position of last ax
    axwi = ax.get_position().width # width of last ax
    axht = ax.get_position().height # height of last ax
    if fliptohoriz:
        if len(fields)==1:
            pass
            cpos=[.91,.15, .02,.7]
        else:
            cpos=[0.91, axposy,.02,axht]
        cbar_ax,cbar = cplt.add_colorbar(fig,pc,orientation='vertical',pos=cpos,label='(%)')
    else:
        if len(fields)==1:
            cpos=[0.12,0.07,0.8,.02]
        else:
            cpos=[axposx,.07, axwi,.02]
        cbar_ax,cbar = cplt.add_colorbar(fig,pc,orientation='horizontal',pos=cpos)
    if field=='sicn':
        cbar.set_ticks([metap['cmin'],-20,-15,-10,-5,0,5,10,15,20,metap['cmax']])
        cbar.set_ticklabels(('',-20,'',-10,'',0,'',10,'',20,''))
    else:
        cbar.set_ticks([metap['cmin'],0,metap['cmax']])
        cbar.set_ticklabels((metap['cmin'],0,metap['cmax']))

if printtofile:
    fig.savefig('suppfig_boundaryconditions' + psuff+'.pdf',bbox_inches='tight')
    fig.savefig('suppfig_boundaryconditions' + psuff+'.eps',bbox_inches='tight',format='eps',dpi=600)


#tst=flddat['sicn'][1]
#cplt.map_allmonths(tst,lat,lon,ptype='nheur',
#                   cmin=-.15,cmax=.15,cmap='red2blue_w20')


printtofile=False


# ======== plot seasonal cycle

simslg=mlines.Line2D([],[],color='0.5',linewidth=2) 
osimlg=mlines.Line2D([],[],color='r',linewidth=2) 


sie=False  # else it's sea ice area
fields=('sicn','sic')
xx=range(0,12)
xxticks=np.arange(0,12)
xxlabs=('Jan','','Mar','','May','','Jul','','Sep','','Nov','')
tlabs = {'gt': 'SST ($^\circ$C)', 
         'sicn': 'SIE (millions of km$^2$)', 
         'sic': 'SIV (thousands of km$^3$)'}


fig,axs = plt.subplots(2,1)
fig.set_size_inches((5.2,5))
fig.subplots_adjust(wspace=0.1,hspace=0.1)

# climos!
fig2,axs2 = plt.subplots(2,1)
fig2.set_size_inches((5.2,5))
fig2.subplots_adjust(wspace=0.1,hspace=0.1)

savedjf={}

for fii,field in enumerate(fields):

    metap = meta[field]
    print field, metap
    conv = convdt[field]

    flddiff=flddat[field]
    fldcavg=fldcdat[field]
    fldpavg=fldpdat[field]

    ax=axs[fii]
    ax2=axs2[fii]

    for ii in range(0,6):

        simfld=flddiff[ii]
        simcfld=fldcavg[ii]
        simpfld=fldpavg[ii]

        if field=='sicn':

            if sie:
                pstr='SIEnh'
                plotfldc,_ = cutl.calc_seaiceextent(simcfld,lat,lon,model=None)
                plotfldp,_ = cutl.calc_seaiceextent(simpfld,lat,lon,model=None)
            else: # sia
                pstr='SIAnh'
                tlabs['sicn'] = 'SIA (millions of km$^2$)'
                plotfldc,_ = cutl.calc_totseaicearea(simcfld,lat,lon,model=None)
                plotfldp,_ = cutl.calc_totseaicearea(simpfld,lat,lon,model=None)

            plotfld = plotfldp-plotfldc
            #plotfld,_ = cutl.calc_totseaicearea(simfld,lat,lon,model=None)
        elif field=='sic':
            tmpc=fldcdat['sicn'][ii]
            tmpp=fldpdat['sicn'][ii]
            plotfldc,_ = cutl.calc_totseaicevol(simcfld*conv,tmpc,lat,lon)
            plotfldp,_ = cutl.calc_totseaicevol(simpfld*conv,tmpp,lat,lon)
            plotfld = plotfldp-plotfldc
            #print conv
            #print 'plotfldc,plotfldp' + str(plotfldc)+',' + str(plotfldp)
        else:
            print '@@ field ' + field + ' not recognized' 
            break
        
        if ii==5:
            clr='r'
        else:
            clr='0.5'
        ax.plot(xx,plotfld, color=clr, linewidth=2)
        if field=='sicn':
            print ii, plotfld[[0,1,11]] # @@@@
            print '   MEAN: ' + str(plotfld[[0,1,11]].mean())
            savedjf[ii] = plotfld[[0,1,11]].mean()

        ax2.plot(xx,plotfldc,color=clr,linewidth=2)
        ax2.plot(xx,plotfldp,color=clr,linewidth=2,linestyle='--')
        
    yticks=ax.get_yticks()
    ax.set_yticklabels(yticks/1e12)
    ax.set_xticks(xxticks)
    if ax.is_last_row():
        ax.set_xticklabels(xxlabs)
    else:
        ax.legend((simslg,osimlg),('Model BCs','Observational BC'),loc='lower left',
                  frameon=False)
        ax.set_xticklabels('')
    ax.set_xlim(0,11)
    ax.set_ylabel('$\Delta$ ' + tlabs[field])

    yticks=ax2.get_yticks()
    ax2.set_yticklabels(yticks/1e12)
    ax2.set_xticks(xxticks)
    if ax2.is_last_row():
        ax2.set_xticklabels(xxlabs)
    else:
        ax2.set_xticklabels('')
    ax2.set_xlim(0,11)
    ax2.set_ylabel(tlabs[field])

if printtofile:
    fig.savefig(pstr+'_SITnh_diff_seacycle_BCs.pdf',bbox_inches='tight')
    fig.savefig(pstr+'_SITnh_diff_seacycle_BCs.eps',bbox_inches='tight',format='eps',dpi=600)
    fig2.savefig(pstr+'_SITnh_climo_seacycle_BCs.pdf',bbox_inches='tight')
    fig2.savefig(pstr+'_SITnh_climo_seacycle_BCs.eps',bbox_inches='tight',format='eps',dpi=600)

#tst=flddat['sicn'][1]
#cplt.map_allmonths(tst,lat,lon,ptype='nheur',
#                   cmin=-.15,cmax=.15,cmap='red2blue_w20')



# @@@@@@@@@@@@@ also, plot regions for paper: supp fig (save code)
cplt.plot_regions(('bksmori','eurasiamori'),colors=('k','blue'),ptype='nheur'); 
if printtofile:
    plt.savefig('suppfig_defineregions.pdf',bbox_inches='tight')
    plt.savefig('suppfig_defineregions.eps',bbox_inches='tight',format='eps',dpi=600)
