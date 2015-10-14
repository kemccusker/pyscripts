import cccmautils as cutl
import constants as con
import loadmodeldata as lmd

deni = 913. # kg/m3 density of ice

field='gt'
timesel='0011-01-01,0011-12-31' # don't have to load much b/c constant

djf=[0,1,11]
sea='DJF'

# from CanAM4_BCmaps (notebook)
if field=='sicn':
    cmap='red2blue_w20'
    cmin=-.15
    cmax=.15
    #cmin=-.30
    #cmax=.30
elif field =='sic':
    #plotfld = plotfld/deni
    cmap='red2blue_w20'
    cmin=-.5
    cmax=.5
elif field == 'gt':
    cmap='blue2red_w20'
    cmin=-2.5
    cmax=2.5
    cmin=-2
    cmax=2
    # need pert SICN too to check if concentration < 0.15
    #sicnp = cnc.getNCvar(con.getBCfilenames('sicn',casenamep),'SICN')
    #sicnp = sicnp[1:,...]
    #plotfld = fldp-fldc
    #plotfld[sicnp>=0.15] = 0 # this is better. @@ but does it make sense when sea ice grows??
    #plotfld = ma.masked_where(sicnp>=0.15,plotfld)

gt={}; sicn={}; sic={}
gt['cmin'] = -2.; gt['cmax'] = 2.; gt['cmap']='blue2red_w20'
sicn['cmin'] = -0.15; sicn['cmax'] = 0.15; sicn['cmap']='red2blue_w20'
sic['cmin'] = -0.5*deni; sic['cmax'] = 0.5*deni; sic['cmap']='red2blue_w20'

meta = {'gt': gt, 'sicn': sicn, 'sic': sic}
tlabs = {'gt': 'SST ($^\circ$C)', 'sicn': 'SIC (%)', 'sic': 'SIT (m)'}
ylabs = {0: 'OBS', 1: 'E1', 2: 'E2', 3: 'E3', 4: 'E4', 5:'E5'} 

fields=('sicn','gt','sic')

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

        # do obs data at ii=0
        if ii==0:
            simname=obsbase+'ctl'
            fnamec = con.getBCfilenames(field,sim=simname)
            print fnamec
            simname=obsbase+'pert'
            fnamep = con.getBCfilenames(field,sim=simname)
        else:
            simname=simcbase+'r'+str(ii)
            fnamec = con.getBCfilenames(field,sim=simname)
            print fnamec
            simname=simpbase+'r'+str(ii)
            fnamep = con.getBCfilenames(field,sim=simname)

        fldc = cnc.getNCvar(fnamec,field.upper(),timesel=timesel)
        fldp = cnc.getNCvar(fnamep,field.upper(),timesel=timesel)

        if field=='sicn':
            # save pert data to mask gt
            sicnpdat[ii] = cutl.seasonalize_monthlyts(fldp,season=sea,climo=1)

        flddiff[ii] = cutl.seasonalize_monthlyts(fldp-fldc,season=sea,climo=1) 
        fnamesc[ii]=fnamec
        fnamesp[ii]=fnamep

    flddat[field] = flddiff



# ============== PLOT ==============
lat = cnc.getNCvar(fnamec,'lat')
lon = cnc.getNCvar(fnamep,'lon')

masksicn=True
printtofile=True

fields=('sicn','sic','gt')

fig,axs = plt.subplots(6,3)
fig.set_size_inches((7,14))
for fii,field in enumerate(fields):

    metap = meta[field]
    print field, metap

    flddiff=flddat[field]

    for ii in range(0,6):

        ax=axs[ii,fii]

        plotfld=flddiff[ii]

        if field=='gt' and masksicn:
            tmp = sicnpdat[ii]
            plotfld[tmp>=0.15] = 0 # mask where model sees grid cell as sea ice (not SST)

        bm,pc = cplt.kemmap(plotfld, lat, lon, ptype='nheur',
                            latlim=45, axis=ax, suppcb=True, 
                            round=False,lmask=True,**metap)

        if ax.is_first_row(): ax.set_title(tlabs[field])
        if ax.is_first_col(): ax.set_ylabel(ylabs[ii])

    # @@@@@@@@ add colorbars!

    #cplt.add_colorbar(fig,pc,orientation='horizontal')
if printtofile:
    fig.savefig('suppfig_boundaryconditions.pdf')
    fig.savefig('suppfig_boundaryconditions.eps',format='eps',dpi=600)



# @@@@@@@@@@@@@ also, plot regions for paper: supp fig (save code)
cplt.plot_regions(('bksmori','eurasiamori'),colors=('k','blue'),ptype='nheur'); 
