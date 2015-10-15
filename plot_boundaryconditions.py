import cccmautils as cutl
import constants as con
import loadmodeldata as lmd


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

    flddat={}; sicnpdat={}
    for fii,field in enumerate(fields):
        fnamesc={}; fnamesp={}
        flddiff={}

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
            else:
                flddiff[ii] = cutl.seasonalize_monthlyts(fldp-fldc,season=sea,climo=1) 
            fnamesc[ii]=fnamec
            fnamesp[ii]=fnamep

        flddat[field] = flddiff


    return flddat, sicnpdat



# === info for plot:

gt={}; sicn={}; sic={}
gt['cmin'] = -2.; gt['cmax'] = 2.; gt['cmap']='blue2red_w20'
sicn['cmin'] = -15; sicn['cmax'] = 15; sicn['cmap']='red2blue_w20'
sic['cmin'] = -0.5; sic['cmax'] = 0.5; sic['cmap']='red2blue_w20'

meta = {'gt': gt, 'sicn': sicn, 'sic': sic}
tlabs = {'gt': '$\Delta$ SST ($^\circ$C)', 'sicn': '$\Delta$ SIC (%)', 'sic': '$\Delta$ SIT (m)'}
ylabs = {5: 'OBS', 0: 'E1', 1: 'E2', 2: 'E3', 3: 'E4', 4:'E5'} 
convdt = {'gt': 1, 'sicn': 100, 'sic': 1/np.float(deni)}


# === load one year for seasonal cycle:
flddat, sicnpdat = load_BCs(fields, season=None, timesel=timesel)


# ============== PLOT ==============
lat = con.get_t63lat() #cnc.getNCvar(fnamec,'lat')
lon = con.get_t63lon() #cnc.getNCvar(fnamep,'lon')


fields=('sicn','sic','gt')

fig,axs = plt.subplots(6,3)
fig.set_size_inches((7,14))
fig.subplots_adjust(wspace=0.1,hspace=0.05)

for fii,field in enumerate(fields):

    metap = meta[field]
    print field, metap
    conv = convdt[field]

    flddiff=flddat[field]

    for ii in range(0,6):

        ax=axs[ii,fii]

        plotfld=flddiff[ii]*conv
        plotfld = cutl.seasonalize_monthlyts(plotfld,season=sea, climo=1)

        if field=='gt' and masksicn:
            tmp = sicnpdat[ii]
            tmp = cutl.seasonalize_monthlyts(tmp,season=sea, climo=1)
            plotfld[tmp>=0.15] = 0 # mask where model sees grid cell as sea ice (not SST)

        bm,pc = cplt.kemmap(plotfld, lat, lon, ptype='nheur',
                            latlim=45, axis=ax, suppcb=True, 
                            round=False,lmask=True,**metap)

        if ax.is_first_row(): ax.set_title(tlabs[field])
        if ax.is_first_col(): ax.set_ylabel(ylabs[ii])

    axposx = ax.get_position().get_points()[0][0] # x position of last row
    axwi = ax.get_position().width # width of last row
    
    cpos=[axposx,.07, axwi,.02]
    cbar_ax,cbar = cplt.add_colorbar(fig,pc,orientation='horizontal',pos=cpos)
    cbar.set_ticks([metap['cmin'],0,metap['cmax']])
    cbar.set_ticklabels((metap['cmin'],0,metap['cmax']))

if printtofile:
    fig.savefig('suppfig_boundaryconditions.pdf')
    fig.savefig('suppfig_boundaryconditions.eps',format='eps',dpi=600)


#tst=flddat['sicn'][1]
#cplt.map_allmonths(tst,lat,lon,ptype='nheur',
#                   cmin=-.15,cmax=.15,cmap='red2blue_w20')

# ======== plot seasonal cycle

fields=('sicn','sic')
xx=range(0,12)
xxticks=np.arange(0,12)
xxlabs=('Jan','','Mar','','May','','Jul','','Sep','','Nov','')
tlabs = {'gt': '$\Delta$ SST ($^\circ$C)', 
         'sicn': '$\Delta$ SIA (millions of km$^2$)', 
         'sic': '$\Delta$ SIT (millions of m$^3$)'}


fig,axs = plt.subplots(2,1)
fig.set_size_inches((7,7))
fig.subplots_adjust(wspace=0.1,hspace=0.1)

for fii,field in enumerate(fields):

    metap = meta[field]
    print field, metap
    conv = convdt[field]

    flddiff=flddat[field]

    ax=axs[fii]

    for ii in range(0,6):

        simfld=flddiff[ii]

        if field=='sicn':
            # Can't do extent here b/c only saved the difference of SICN.
            #  doesn't make sense to calc SIE @@@@
            #plotfld,_ = cutl.calc_seaiceextent(simfld,lat,lon,model=None)
            plotfld,_ = cutl.calc_totseaicearea(simfld,lat,lon,model=None)
        elif field=='sic':
            tmp=flddat['sicn'][ii]
            # have to do CTL and PERT separately to get this right. @@@@@
            plotfld,_ = cutl.calc_totseaicevol(simfld*conv,tmp,lat,lon)
        else:
            print '@@ field ' + field + ' not recognized' 
            break
        
        if ii==5:
            clr='r'
        else:
            clr='0.5'
        ax.plot(xx,plotfld, color=clr, linewidth=2)
        
    yticks=ax.get_yticks()
    #ax.set_yticklabels(yticks/1e12)
    #ax.set_xticks(xxticks)
    if ax.is_last_row():
        ax.set_xticklabels(xxlabs)
    else:
        ax.set_xticklabels('')
    ax.set_xlim(0,11)
    ax.set_ylabel(tlabs[field])


tst=flddat['sicn'][1]
cplt.map_allmonths(tst,lat,lon,ptype='nheur',
                   cmin=-.15,cmax=.15,cmap='red2blue_w20')



# @@@@@@@@@@@@@ also, plot regions for paper: supp fig (save code)
cplt.plot_regions(('bksmori','eurasiamori'),colors=('k','blue'),ptype='nheur'); 
if printtofile:
    plt.savefig('suppfig_defineregions.pdf')
    plt.savefig('suppfig_defineregions.eps',format='eps',dpi=600)
