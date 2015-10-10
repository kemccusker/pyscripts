

import loadmodeldata as lmd

deni = 913. # kg/m3 density of ice

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
    sicnp = cnc.getNCvar(con.getBCfilenames('sicn',casenamep),'SICN')
    sicnp = sicnp[1:,...]
    plotfld = fldp-fldc
    plotfld[sicnp>=0.15] = 0 # this is better. @@ but does it make sense when sea ice grows??
    #plotfld = ma.masked_where(sicnp>=0.15,plotfld)



# load obs BC here:


# now load all 5 sim BCs
field='gt'
for ii in range(1,6):

    simname='kemctlr'+str(ii)
    fnamesc[ii] = con.getBCfilenames(field,sim=simname)
    

    simname='kempert2r'+str(ii)
    fnamesp[ii] = con.getBCfilenames(field,sim=simname)
