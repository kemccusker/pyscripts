# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Create a plot of streamfunction in southern ocean -- geotimescales project.<br>
# Overlay it on mean and anomalous temperature. Then calculate heat transport of anomalous W working on the mean Temp.<br>
# <br>
# <ul>
# <li>Figure 1: Climatological Eulerian MOC (1970-1999 mean)</li>
# <li>Figure 2: Sulfate engineering anomalous Eulerian MOC (2045-2054 compared to climo)</li>
# <li>Figure 3: GHG removal anomalous Eulerian MOC (2045-2054 compared to climo)</li>
# <li>Figure 4: Climatological TEMP (1970-1999 mean) with sulfate engineering anom eulerian MOC (contours)</li>
# </ul>
# <br>
# For figures 1-4: + MOC has solid lines (cw), - MOC has dashed lines (ccw). no zero contour.
# <br><br>
# <ul>
# <li>Figure 5: Anomalous TEMP (shading) with anomalous eulerian MOC (contours) for sulfate eng (left), ghg removal (right)</li>
# <li>Figure 6: Same as Fig. 5 but anomalous eddy-induced MOC</li>
# <li>Figure 7: Same as Fig. 5 but anomalous total MOC (eulerian + eddy-induced)</li>
# <li>Figure 8: Same as Fig. 7 but zoomed closer in on SH region</li>
# <li>Figure 9: Same as Fig. 5 but anomalous submeso-scale eddy MOC</li>
# </ul>
# <br>
# For figures 5-9: + MOC has thick solid lines (cw), - MOC has thin solid lines (ccw). no zero contour.
# <br><br>
# <ul>
# <li>Figure 10a: components of heat trans calc for sulfates run: climo dT/dz; anomalous WVEL+WISOP; anom WVEL; anom WISOP</li>
# <li>Figure 10b: Heating rate (wprime*dTbar/dz; K/s) at each level for sulfate engineering (left) and GHG removal (right)</li>
# <li>Figure 11: same as Fig. 10 but in K/day and zoomed in</li>
# <li>Figure 12: Same as Fig. 10, but multiply by layer thickness to convert to W/m2 in each grid cell. Overlaid with climo TEMP contours (gray) and anomalous eulerian MOC (green)</li>
# </ul>
# <br>
# For figures 10-12: having trouble removing the spurious data in lower left. It's probably related to the bathymetry and not appropriately handling missing data.
# <br><br>
# <em>Figure headings are below the figures</em>

# <codecell>

#%matplotlib inline 

import cccmaplots as cplt
import cccmaNC as cnc

plt.close('all')

printtofile=False

basepath = '/Users/kelly/School/DATA/'

casenamec = 'b40.20th.track1.1deg.006'
casenamep = 'geo2035ensavg' # or 'rcp8_5GHGrem1850'

filenamec = basepath + casenamec + '/' + casenamec + '.pop.ANN.1970-1999.nc'
filenamep = basepath + casenamep + '/' + casenamep + '.pop.ANN.2045-2054.nc'

print filenamec
print filenamep

pig=True # do pine island glacier region # @@@ testing

# From ocean_temps_vert2_SOforthesis.m @@@
# indices in the lon direction for the various SO basins
# probably need to shift by 1 for python (0-based), also I don't think this indexing works in python
#iwed = [294:320 1:30];
#iross = [183:227];
#ipig = [250:276]; 
# 80W to 120W ? http://www.awi.de/fileadmin/user_upload/News/Press_Releases/2013/3._Quartal/Pine_Island_Glacier/PIG_map_Bodentopograpfie_beschriftet_p.jpg
#if pig: 
#    gridslicelonw = 80;
#    gridslicelone = 120;

# <codecell>

latauxgrid = cnc.getNCvar(filenamec,'lat_aux_grid',sqz=False)
transreg = cnc.getNCvar(filenamec, 'transport_regions',sqz=False)
moccomp = cnc.getNCvar(filenamec, 'moc_components',sqz=False)
mocz = cnc.getNCvar(filenamec,'moc_z',sqz=False)


""" float MOC(time, transport_reg, moc_comp, moc_z, lat_aux_grid) ;
                MOC:long_name = "Meridional Overturning Circulation" ;
                MOC:units = "Sverdrups" ;
                MOC:coordinates = "lat_aux_grid moc_z moc_components transport_region time" ;
                MOC:missing_value = 9.96921e+36f ;
                
                
    transport_regions =
  "Global Ocean - Marginal Seas",
  "Atlantic Ocean + Mediterranean Sea + Labrador Sea + GIN Sea + Arctic Ocean + Hudson Bay" 
  
    moc_components =
  "Eulerian Mean",
  "Eddy-Induced (bolus)",
  "Submeso" ;
   
     float moc_z(moc_z) ;
                moc_z:long_name = "depth from surface to top of layer" ;
                moc_z:units = "centimeters" ;
                moc_z:positive = "down" ;
                moc_z:valid_min = 0.f ;
                moc_z:valid_max = 549999.1f ;
                
    moc_z = 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 
    11000, 12000, 13000, 14000, 15000, 16000, 17019.68, 18076.13, 19182.12, 
    20349.93, 21592.34, 22923.31, 24358.45, 25915.58, 27615.26, 29481.47, 
    31542.37, 33831.23, 36387.47, 39258.05, 42498.89, 46176.66, 50370.69, 
    55174.91, 60699.67, 67072.86, 74439.8, 82960.7, 92804.35, 104136.8, 
    117104, 131809.4, 148290.1, 166499.2, 186301.4, 207487.4, 229803.9, 
    252990.4, 276809.8, 301067.1, 325613.8, 350344.9, 375189.2, 400101.2, 
    425052.5, 450026.1, 475012, 500004.7, 525000.9, 549999.1 ;
"""

totmocc = cnc.getNCvar(filenamec,'MOC',sqz=False)
totmocp = cnc.getNCvar(filenamep,'MOC',sqz=False)

print totmocc.shape 
mocc=totmocc[0,0,0,...]
print mocc.shape

print totmocp.shape
mocp=totmocp[0,0,0,...]

# <codecell>

lats,levs = np.meshgrid(latauxgrid,mocz/100.)

contsp = np.arange(0,36,3)#[0,2,4,6,8,10,12,14,16,18,20]
contsn = np.arange(-10,-.5,.5)

fig = plt.figure()
ax = fig.add_subplot(111)
CS1 = plt.contour(lats,levs,mocc,contsp,\
            colors='k',linestyles='solid')
plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CS1 = plt.contour(lats,levs,mocc,contsn,\
            colors='k',linestyles='dashed')
plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)

ax.set_ylim((0,3000))
ax.invert_yaxis()
ax.set_xlim((-80,-40))
ax.set_title(casenamec + ' Eul MOC')

# <headingcell level=3>

# Fig. 1: Climo Eulerian MOC (1970-1999)

# <codecell>

contspd = np.arange(0,5,.2)
contsnd = np.arange(-5,-.2,.3)

fig = plt.figure()
ax = fig.add_subplot(111)
CS1 = plt.contour(lats,levs,mocp-mocc,contspd,\
            colors='k',linestyles='solid')
plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CS1 = plt.contour(lats,levs,mocp-mocc,contsnd,\
            colors='k',linestyles='dashed')
plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
#CS1 = plt.contour(lats,levs,mocp-mocc,[0,0],\
#            colors='r',linestyles='solid')

ax.set_ylim((0,3000))
ax.invert_yaxis()
ax.set_xlim((-80,-40))
ax.set_title(casenamep + ' - ' + casenamec + ' Eul MOC anomaly')

# <headingcell level=3>

# Fig. 2: Sulfate engineering (geo2035ensavg): Anomalous Eulerian MOC

# <codecell>

casenamep2 = 'rcp8_5GHGrem1850' # or 'rcp8_5GHGrem1850'
filenamep2 = basepath + casenamep2 + '/' + casenamep2 + '.pop.ANN.2045-2054.nc'

print filenamep
totmocp2 = cnc.getNCvar(filenamep2,'MOC',sqz=False)
mocp2=totmocp2[0,0,0,...]

contspd = np.arange(.2,5,.2)
contsnd = np.arange(-5,-.2,.3)


fig = plt.figure()
ax = fig.add_subplot(111)
CS1 = plt.contour(lats,levs,mocp2-mocc,contspd,\
            colors='k',linestyles='solid')
plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CS1 = plt.contour(lats,levs,mocp2-mocc,contsnd,\
            colors='k',linestyles='dashed')
plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
#CS1 = plt.contour(lats,levs,mocp-mocc,[0,0],\
#            colors='r',linestyles='solid')

ax.set_ylim((0,3000))
ax.invert_yaxis()
ax.set_xlim((-80,-40))
ax.set_title(casenamep2 + ' - ' + casenamec + ' Eul MOC anomaly')

# <headingcell level=3>

# Fig. 3: GHG removal (rcp8_5GHGrem1850): Anomalous Eulerian MOC

# <codecell>

kmt = cnc.getNCvar(filenamec,'KMT')


# Now add TEMP
tempc = np.squeeze(cnc.getNCvar(filenamec,'TEMP'))
tempp = np.squeeze(cnc.getNCvar(filenamep,'TEMP'))
tempp2 = np.squeeze(cnc.getNCvar(filenamep2,'TEMP'))
zt = cnc.getNCvar(filenamec, 'z_t')
tlat = cnc.getNCvar(filenamec,'TLAT')
tlon = cnc.getNCvar(filenamec,'TLONG')
tarea=cnc.getNCvar(filenamec,'TAREA')


import numpy.ma as ma


print tempc.shape

for lii,zz in enumerate(zt):
    # first mask out levels below sea floor
    tempc[lii,...] = ma.masked_where(kmt <= lii,tempc[lii,...])
    tempp[lii,...] = ma.masked_where(kmt <= lii,tempp[lii,...])
    tempp2[lii,...] = ma.masked_where(kmt <= lii,tempp2[lii,...])
    
if pig:  # 80W to 120W. Or 280 to 230
    printtofile=False

    """ @@@@ need to weight the grid cells by area, even for zonal means. @@@@
    piglon = range(230,281) # lon indices
    piglat = range(0,187) # lat indices
    kmtpig=kmt[piglat,:]
    kmtpig=kmtpig[:,piglon]
    tareapig=tarea[piglat,:]
    tareapig=tareapig[:,piglon]

    kmtpigt = np.tile(kmtpig,(fld.shape[0],1,1))

    tareapigt = np.tile(tareapig,(len(zt),1,1))
    tareat = np.tile(tarea,(len(zt),1,1))
    """

    lonlims = [230,280]; region = 'PIG'; strlims='80W-120W'
    #lonlims = [230,260]; region='PIG2'; strlims='100W-120W'
    
    rmaskout = np.logical_and(tlon>lonlims[0], tlon<lonlims[1]) # region mask! this masks OUT the region itself
    rmask = np.logical_or(tlon<=lonlims[0],tlon>=lonlims[1]) # use this one for averaging. keep only the region
    testmask = tempc[0,...]
    testmask = ma.masked_where(rmask,testmask)
    
    bm = cplt.kemmap(testmask,tlat[:,0],tlon[0,:],title='region mask',type='sh')
    if printtofile:
        plt.savefig(region + '_' + strlims + '_map.pdf')

    # tile the mask
    rmask = np.tile(rmask,(len(zt),1,1))
    print rmask.shape

    tempcreg = ma.masked_where(rmask,tempc)
    temppreg = ma.masked_where(rmask,tempp)
    tempp2reg = ma.masked_where(rmask,tempp2)

    # do I need to also take ocean fraction into account?!
    tempcreg=np.squeeze(np.mean(tempcreg,axis=2))
    temppreg=np.squeeze(np.mean(temppreg,axis=2))
    tempp2reg=np.squeeze(np.mean(tempp2reg,axis=2))

    tlats,zlevs = np.meshgrid(np.squeeze(tlat[:,1]),zt/100.)
    cmap='jet'
    cmin=-2; cmax=8

    #cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
    cmlen=float(30)
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    CF1 = plt.contourf(tlats,zlevs,tempcreg,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

    ax.set_ylim((0,800))
    ax.invert_yaxis()
    ax.set_xlim((-75,-50))
    ax.set_title(region + ' climo TEMP')
    cbar = fig.colorbar(CF1)


# # ====================== paper ====================
# # PIG zonal mean TEMP
    #cmap='blue2red_w20'
    cmap='blue2red_w20'  # @@@
    cmin=-.5; cmax=.5
    

    printtofile= True
    ylim=1000
    #cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
    cmlen=float(20)
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    #fig = plt.figure()
    fig,axs=plt.subplots(1,2,sharey=True)
    ax=axs[0]
    fig.set_size_inches(14,3)
    #ax = fig.add_subplot(121)
    CF1 = ax.contourf(tlats,zlevs,temppreg-tempcreg,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

    # contours for MOC anomaly
    contspd = np.arange(.1,5,.3)
    contsnd = np.arange(-5,-.1,.3)

    ax.set_ylim((0,ylim))
    ax.invert_yaxis()
    ax.set_xlim((-75,-50))
    ax.set_yticks(np.arange(0,900,100))
    ax.set_yticklabels([0,'',200,'',400,'',600,'',800],fontsize=18)
    ax.set_xticks(np.arange(-75,-45,5))
    ax.set_xticklabels(['75$^\circ$S', '70$^\circ$S', '65$^\circ$S', \
                        '60$^\circ$S', '55$^\circ$S', '50$^\circ$S'],fontsize=18)

    #ax.set_title(casenamep + ' ' + region + ' anom TEMP')
    ax.set_title('PIG Sulf',fontsize=18)
    ax.set_ylabel('Depth (m)',fontsize=18)
    ax.axvline(x=-65,linestyle='--',color='k') # @@@ the vert line is to show the area averaged in the VHT plots
    #cbar = fig.colorbar(CF1)

    #ax2 = fig.add_subplot(122)
    ax2=axs[1]
    CF2 = ax2.contourf(tlats,zlevs,tempp2reg-tempcreg,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

    ax2.set_ylim((0,ylim))
    ax2.invert_yaxis()
    ax2.set_xlim((-75,-50))
    ax2.set_yticks(np.arange(0,900,100))
    ax2.set_yticklabels([0,'',200,'',400,'',600,'',800],fontsize=18)
    ax2.set_xticks(np.arange(-75,-45,5))
    ax2.set_xticklabels(['75$^\circ$S', '70$^\circ$S', '65$^\circ$S', \
                        '60$^\circ$S', '55$^\circ$S', '50$^\circ$S'],fontsize=18)
    #ax2.set_title(casenamep2 + ' ' + region + ' anom TEMP')
    ax2.set_title('PIG GHGrem',fontsize=18)
    ax2.axvline(x=-65,linestyle='--',color='k') # @@@ the vert line is to show the area averaged in the VHT plots
    #cbar = fig.colorbar(CF2)
    cbar_ax = fig.add_axes([.91,.15, .02,.7])
    fig.colorbar(CF2,cax=cbar_ax)

    if printtofile:
        fig.savefig('TEMPanom_subplotSHzm_' + region + '_c_ylim' + str(ylim) + '.pdf')


    # ===== TEST figure showing the difference b/w Sulf and GHGrem differences

    printtofile=False

    cmin=-1.2; cmax=1.2
    cmap='blue2red_20'

    cmlen=float(20)
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    fig = plt.figure()
    fig.set_size_inches(7,3)
    ax = fig.add_subplot(111)
    CF1 = plt.contourf(tlats,zlevs,(temppreg-tempcreg)-(tempp2reg-tempcreg),
                       cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

    # contours for MOC anomaly
    contspd = np.arange(.1,5,.3)
    contsnd = np.arange(-5,-.1,.3)

    ax.set_ylim((0,ylim))
    ax.invert_yaxis()
    ax.set_xlim((-75,-50))
    ax.set_title(casenamep + '-' + casenamep2 + ' ' + region + ' anom TEMP')
    cbar = fig.colorbar(CF1)
    if printtofile:
        fig.savefig('TEMPanom_Sulf-GHGrem_' + region + '.pdf')

    # end if pig region
    # #########################


tempc=np.squeeze(np.mean(tempc,axis=2))
tempp=np.squeeze(np.mean(tempp,axis=2))
tempp2=np.squeeze(np.mean(tempp2,axis=2))

print tempc.shape

# <codecell>

tlats,zlevs = np.meshgrid(np.squeeze(tlat[:,1]),zt/100.)
cmap='jet'
cmin=-2; cmax=8

#cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
cmlen=float(30)
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

fig = plt.figure()
ax = fig.add_subplot(111)
CF1 = plt.contourf(tlats,zlevs,tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')


contspd = np.arange(.2,5,.2)
contsnd = np.arange(-5,-.2,.3)


CS1 = plt.contour(lats,levs,mocp-mocc,contspd,\
            colors='k',linestyles='solid')
plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CS1 = plt.contour(lats,levs,mocp-mocc,contsnd,\
            colors='k',linestyles='dashed')
plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)


ax.set_ylim((0,3000))
ax.invert_yaxis()
ax.set_xlim((-80,-40))
ax.set_title('climo TEMP & ' + casenamep + ' anom Eul MOC')
cbar = fig.colorbar(CF1)

# <headingcell level=3>

# Fig. 4: Climo TEMP (shading) with sulfate enginering anomalous Eulerian MOC (contours)

# <codecell>

# EULERIAN MOC


cmap='blue2red_w20'
cmin=-.5; cmax=.5

printtofile=False

moctype='eul'

#cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
cmlen=float(20)
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

fig = plt.figure()
fig.set_size_inches(14,3)
ax = fig.add_subplot(121)
CF1 = plt.contourf(tlats,zlevs,tempp-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

# contours for MOC anomaly
contspd = np.arange(.1,5,.3)
contsnd = np.arange(-5,-.1,.3)


CS1 = plt.contour(lats,levs,mocp-mocc,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CS1 = plt.contour(lats,levs,mocp-mocc,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)

ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_xlim((-80,-40))
ax.set_title(casenamep + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF1)



ax2 = fig.add_subplot(122)
CF2 = plt.contourf(tlats,zlevs,tempp2-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')


CS2 = plt.contour(lats,levs,mocp2-mocc,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)
CS2 = plt.contour(lats,levs,mocp2-mocc,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)

ax2.set_ylim((0,2000))
ax2.invert_yaxis()
ax2.set_xlim((-80,-40))
ax2.set_title(casenamep2 + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF2)
if printtofile:
    fig.savefig('MOC' + moctype + 'anom_climoTEMP_subplotSHzm.pdf')

# <headingcell level=3>

# Fig. 5: Anomalous TEMP (shading) with anomalous Eulerian MOC in contours (left: sulfate, right: ghgrem)

# <codecell>

# EDDY-INDUCED MOC

cmap='blue2red_w20'
cmin=-.5; cmax=.5

#cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
cmlen=float(20)
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

plotmoc = totmocp[0,0,1,...]-totmocc[0,0,1,...] # eddy-induced
plotmoc2 = totmocp2[0,0,1,...]-totmocc[0,0,1,...] # eddy-induced

moctype='eddy'

fig = plt.figure()
fig.set_size_inches(14,3)
ax = fig.add_subplot(121)
CF1 = plt.contourf(tlats,zlevs,tempp-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

# contours for MOC anomaly
contspd = np.arange(.2,5,.3)
contsnd = np.arange(-5,-.2,.3)

CS1 = plt.contour(lats,levs,plotmoc,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CS1 = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)

ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_xlim((-80,-40))
ax.set_title(casenamep + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF1)



ax2 = fig.add_subplot(122)
CF2 = plt.contourf(tlats,zlevs,tempp2-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')


CS2 = plt.contour(lats,levs,plotmoc2,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)
CS2 = plt.contour(lats,levs,plotmoc2,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)

ax2.set_ylim((0,2000))
ax2.invert_yaxis()
ax2.set_xlim((-80,-40))
ax2.set_title(casenamep2 + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF2)

if printtofile:
    fig.savefig('MOC' + moctype + 'anom_climoTEMP_subplotSHzm.pdf')

# <headingcell level=3>

# Fig. 6: Anomalous TEMP (shading) with anomalous Eddy-induced MOC in contours (left: sulfate, right: ghgrem)

# <codecell>

# TOTAL (EULERIAN + EDDY-INDUCED) MOC

cmap='blue2red_w20'
cmin=-.5; cmax=.5

printtofile=False

#cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
cmlen=float(20)
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

plotmoc = (totmocp[0,0,1,...]+totmocp[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced
plotmoc2 = (totmocp2[0,0,1,...]+totmocp2[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced

moctype='eul+edd'

fig = plt.figure()
fig.set_size_inches(14,3)
ax = fig.add_subplot(121)
CF1 = plt.contourf(tlats,zlevs,tempp-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

# contours for MOC anomaly
contspd = np.arange(.1,5,.3)
contsnd = np.arange(-5,-.1,.3)

CS1 = plt.contour(lats,levs,plotmoc,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CS1 = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)

ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_xlim((-80,-40))
ax.set_title(casenamep + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF1)



ax2 = fig.add_subplot(122)
CF2 = plt.contourf(tlats,zlevs,tempp2-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')


CS2 = plt.contour(lats,levs,plotmoc2,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)
CS2 = plt.contour(lats,levs,plotmoc2,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)

ax2.set_ylim((0,2000))
ax2.invert_yaxis()
ax2.set_xlim((-80,-40))
ax2.set_title(casenamep2 + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF2)

if printtofile:
    fig.savefig('MOC' + moctype + 'anom_climoTEMP_subplotSHzm.pdf')

# <headingcell level=3>

# Fig. 7: Same as above, but total anomalous MOC (Eulerian + eddy-induced)

# <codecell>

# ZOOM: TOTAL (EULERIAN + EDDY-INDUCED) MOC

cmap='blue2red_w20'
cmin=-.5; cmax=.5
printtofile=False

#cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
cmlen=float(20)
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

plotmoc = (totmocp[0,0,1,...]+totmocp[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced
plotmoc2 = (totmocp2[0,0,1,...]+totmocp2[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced

moctype='eul+edd'

fig = plt.figure()
fig.set_size_inches(14,3)
ax = fig.add_subplot(121)
CF1 = plt.contourf(tlats,zlevs,tempp-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

# contours for MOC anomaly
#contspd = np.arange(.01,5,.1)
contspd = [.05,.1,.2,.3,.4,.5,1,1.5]

#contsnd = np.arange(-5,-.01,.1)
contsnd = [-1.5,-1,-.5,-.4,-.3,-.2,-.1,-.05]

CS1 = plt.contour(lats,levs,plotmoc,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CS1 = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)

ax.set_ylim((0,1000))
ax.invert_yaxis()
ax.set_xlim((-80,-55))
ax.set_title(casenamep + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF1)



ax2 = fig.add_subplot(122)
CF2 = plt.contourf(tlats,zlevs,tempp2-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')


CS2 = plt.contour(lats,levs,plotmoc2,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)
CS2 = plt.contour(lats,levs,plotmoc2,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)

ax2.set_ylim((0,1000))
ax2.invert_yaxis()
ax2.set_xlim((-80,-55))
ax2.set_title(casenamep2 + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF2)

if printtofile:
    fig.savefig('MOC' + moctype + 'anom_climoTEMP_subplotSHzmZOOM.pdf')

# <headingcell level=3>

# Fig. 8: Same as Fig. 7 (total MOC) but zoomed in

# <codecell>

# SUBMESO-SCALE EDDY

cmap='blue2red_w20'
cmin=-.5; cmax=.5
printtofile=False


#cmlen=float( plt.cm.get_cmap(cmap).N) # or: from __future__ import division
cmlen=float(20)
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

plotmoc = totmocp[0,0,2,...]-totmocc[0,0,2,...] # submeso
plotmoc2 = totmocp2[0,0,2,...]-totmocc[0,0,2,...] # submeso

moctype='submeso'

fig = plt.figure()
fig.set_size_inches(14,3)
ax = fig.add_subplot(121)
CF1 = plt.contourf(tlats,zlevs,tempp-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

# contours for MOC anomaly
contspd = np.arange(.1,5,.2)
contsnd = np.arange(-5,-.1,.2)

CS1 = plt.contour(lats,levs,plotmoc,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CS1 = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)

ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_xlim((-80,-40))
ax.set_title(casenamep + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF1)



ax2 = fig.add_subplot(122)
CF2 = plt.contourf(tlats,zlevs,tempp2-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')


CS2 = plt.contour(lats,levs,plotmoc2,contspd,\
            colors='k',linestyles='solid',linewidths=2)
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)
CS2 = plt.contour(lats,levs,plotmoc2,contsnd,\
            colors='k',linestyles='solid')
#plt.clabel(CS2,fmt = '%2.1f',inline=1,fontsize=10)

ax2.set_ylim((0,2000))
ax2.invert_yaxis()
ax2.set_xlim((-80,-40))
ax2.set_title(casenamep2 + ' anom TEMP (shade) & ' + moctype + ' MOC (cont)')
cbar = fig.colorbar(CF2)

if printtofile:
    fig.savefig('MOC' + moctype + 'anom_climoTEMP_subplotSHzm.pdf')

# <headingcell level=3>

# Fig. 9: same as Figs 5,6, etc but submeso-scale eddy MOC

# <headingcell level=3>

# Next, Calculate wprime*dTbar/dz at all levels

# <codecell>

print filenamec
print filenamep
print filenamep2

kmt = cnc.getNCvar(filenamec,'KMT')

wvc = np.squeeze(cnc.getNCvar(filenamec,'WVEL'))/100. # convert to m/s
wic = np.squeeze(cnc.getNCvar(filenamec,'WISOP'))/100.

wvp = np.squeeze(cnc.getNCvar(filenamep,'WVEL'))/100.
wip = np.squeeze(cnc.getNCvar(filenamep,'WISOP'))/100.

wvp2 = np.squeeze(cnc.getNCvar(filenamep2,'WVEL'))/100.
wip2 = np.squeeze(cnc.getNCvar(filenamep2,'WISOP'))/100.

# for each level, calc wprime*dTvar/dz
wprime = (wvp+wip)-(wvc+wic)
wprime2= (wvp2+wip2)-(wvc+wic)
wvelprime = wvp-wvc
wisopprime = wip-wic
wvelprime2 = wvp2-wvc
wisopprime2 = wip2-wic

# also wbar*dTprime/dz
wbar = wvc+wic
wisopbar = wic
wvelbar = wvc

print wprime.shape

# use index-of deepest grid cell to determine if the cell should be included in zonal mean
#   for each level: 
#       if kmt<=level index:
#            mask the cell
# aka: for each level, mask where kmt<=level index
import numpy.ma as ma
import scipy.stats

for lii,zz in enumerate(zt):
    # first mask out levels below sea floor
    wprime[lii,...] = ma.masked_where(kmt <= lii,wprime[lii,...])
    wprime2[lii,...] = ma.masked_where(kmt <= lii,wprime2[lii,...])
    wvelprime[lii,...] = ma.masked_where(kmt <= lii,wvelprime[lii,...])
    wvelprime2[lii,...] = ma.masked_where(kmt <= lii,wvelprime2[lii,...])
    wisopprime[lii,...] = ma.masked_where(kmt <= lii,wisopprime[lii,...])
    wisopprime2[lii,...] = ma.masked_where(kmt <= lii,wisopprime2[lii,...])
    wbar[lii,...] = ma.masked_where(kmt <= lii,wbar[lii,...])
    wvelbar[lii,...] = ma.masked_where(kmt <= lii,wvelbar[lii,...])
    wisopbar[lii,...] = ma.masked_where(kmt <= lii,wisopbar[lii,...])

#plt.figure()
#plt.pcolor(wprime[10,...],vmin=-2e-6,vmax=2e-6)
#plt.colorbar()

# Now take zonal mean    
#wprime = np.squeeze(sp.stats.nanmean(wprime,axis=2)) # give wtrans of Inf and 0 in most spots?
#wprime2 = np.squeeze(sp.stats.nanmean(wprime2,axis=2))

if pig:
    wprimereg = ma.masked_where(rmask,wprime)
    wprime2reg = ma.masked_where(rmask,wprime2)
    wvelprimereg = ma.masked_where(rmask,wvelprime)
    wvelprime2reg = ma.masked_where(rmask,wvelprime2)
    wisopprimereg = ma.masked_where(rmask,wisopprime)
    wisopprime2reg = ma.masked_where(rmask,wisopprime2)

    wbarreg = ma.masked_where(rmask,wbar)
    wvelbarreg = ma.masked_where(rmask,wvelbar)
    wisopbarreg = ma.masked_where(rmask,wisopbar)
    
    wprimereg = np.squeeze(np.mean(wprimereg,axis=2))
    wprime2reg = np.squeeze(np.mean(wprime2reg,axis=2))
    wvelprimereg = np.squeeze(np.mean(wvelprimereg,axis=2))
    wvelprime2reg = np.squeeze(np.mean(wvelprime2reg,axis=2))
    wisopprimereg = np.squeeze(np.mean(wisopprimereg,axis=2))
    wisopprime2reg = np.squeeze(np.mean(wisopprime2reg,axis=2))

    wbarreg = np.squeeze(np.mean(wbarreg,axis=2))
    wvelbarreg = np.squeeze(np.mean(wvelbarreg,axis=2))
    wisopbarreg = np.squeeze(np.mean(wisopbarreg,axis=2))


    tbarreg = tempcreg# already zonal meaned MEAN T
    dtbarreg = np.diff(tbarreg,axis=0) # delta of MEAN T with height

    tprimereg = temppreg-tempcreg
    tprime2reg = tempp2reg-tempcreg
    dtprimereg = np.diff(tprimereg, axis=0) # delta of ANOM T with height
    dtprime2reg = np.diff(tprime2reg, axis=0)

    # thickness of each layer
    dzt = np.diff(zt/100.) # convert to m
    wtransreg = np.zeros((len(dzt),wprimereg.shape[1]))
    wtrans2reg = np.zeros((len(dzt),wprime2reg.shape[1]))
    wtranswvreg = np.zeros((len(dzt),wprimereg.shape[1]))
    wtranswv2reg = np.zeros((len(dzt),wprime2reg.shape[1]))
    wtranswireg = np.zeros((len(dzt),wprimereg.shape[1]))
    wtranswi2reg = np.zeros((len(dzt),wprime2reg.shape[1]))

    wbartransreg = np.zeros((len(dzt),wbarreg.shape[1])) # mean w time anom dT
    wbartrans2reg = np.zeros((len(dzt),wbarreg.shape[1])) # mean w time anom dT
    wbartranswvreg = np.zeros((len(dzt),wbarreg.shape[1]))
    wbartranswv2reg = np.zeros((len(dzt),wbarreg.shape[1]))
    wbartranswireg = np.zeros((len(dzt),wbarreg.shape[1]))
    wbartranswi2reg = np.zeros((len(dzt),wbarreg.shape[1]))

    # calc heat transport (heating rate) for each level
    for lii,dz in enumerate(dzt):
        #print 'ind: ' + str(lii) + ', dz: ' + str(dz)
        #  W prime * (dTbar / dz)
        wtransreg[lii,...] = wprimereg[lii,...]*(dtbarreg[lii,...]/dz)
        wtrans2reg[lii,...] = wprime2reg[lii,...]*(dtbarreg[lii,...]/dz)

        wtranswvreg[lii,...] = wvelprimereg[lii,...]*(dtbarreg[lii,...]/dz) # WVEL only
        wtranswv2reg[lii,...] = wvelprime2reg[lii,...]*(dtbarreg[lii,...]/dz)

        wtranswireg[lii,...] = wisopprimereg[lii,...]*(dtbarreg[lii,...]/dz) # WISOP only
        wtranswi2reg[lii,...] = wisopprime2reg[lii,...]*(dtbarreg[lii,...]/dz)

        # W bar * (dTprime / dz)
        wbartransreg[lii,...] = wbarreg[lii,...]*(dtprimereg[lii,...]/dz)
        wbartrans2reg[lii,...] = wbarreg[lii,...]*(dtprime2reg[lii,...]/dz)

        wbartranswvreg[lii,...] = wvelbarreg[lii,...]*(dtprimereg[lii,...]/dz)
        wbartranswv2reg[lii,...] = wvelbarreg[lii,...]*(dtprime2reg[lii,...]/dz)

        wbartranswireg[lii,...] = wisopbarreg[lii,...]*(dtprimereg[lii,...]/dz)
        wbartranswi2reg[lii,...] = wisopbarreg[lii,...]*(dtprime2reg[lii,...]/dz)


# # 
# #   PIG zonal mean vertical heat trans and velocity subplot
    printtofile=False
    meshlats,meshdz = np.meshgrid(np.squeeze(tlat[:,1]),zt[1:]/100.)

    cmlen=float(20)
    ylims = (0,500)
    xlims = (-77,-40)

    fig = plt.figure()
    fig.set_size_inches(14,8)

    cmax=.8; cmin=-.8
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(241)
    plt.contourf(meshlats,meshdz,dtbarreg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both') 
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title(region + ' dTbar/dz')

    cmax=1e-6; cmin=-1e-6
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(242)
    plt.contourf(tlats,zlevs,wprimereg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title('sulfates: ' + region + ' wprime')

    cmax=3e-7; cmin=-3e-7
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(243)
    plt.contourf(tlats,zlevs,wvelprimereg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title('sulfates: ' + region + ' wvelprime')

    cmax=3e-7; cmin=-3e-7
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(244)
    plt.contourf(tlats,zlevs,wisopprimereg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title('sulfates: ' + region + ' wisopprime')

    # do GHGrem =================
    cmax=1e-6; cmin=-1e-6
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(246)
    plt.contourf(tlats,zlevs,wprime2reg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title('GHGrem: ' + region + ' wprime')

    cmax=3e-7; cmin=-3e-7
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(247)
    plt.contourf(tlats,zlevs,wvelprime2reg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title('GHGrem: ' + region + ' wvelprime')

    cmax=3e-7; cmin=-3e-7
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(248)
    plt.contourf(tlats,zlevs,wisopprime2reg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title('GHGrem: ' + region + ' wisopprime')

    if printtofile:
        fig.savefig('dTbardz_Wcomps_Sulf_GHGrem_subplotSHzm_' + region + '.pdf')

    # THIS PLOT IS W bar * dTprime ! ==================
# #   PIG zonal mean vertical heat trans and velocity subplot
    printtofile=False
    meshlats,meshdz = np.meshgrid(np.squeeze(tlat[:,1]),zt[1:]/100.)

    cmlen=float(20)
    ylims = (0,500)
    xlims = (-77,-40)

    fig = plt.figure()
    fig.set_size_inches(14,8)

    cmax=.1; cmin=-.1
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(241)
    plt.contourf(meshlats,meshdz,dtprimereg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both') 
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title('Sulf: ' + region + ' dTprime/dz')

    cmax=.1; cmin=-.1
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(245) # second row
    plt.contourf(meshlats,meshdz,dtprime2reg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both') 
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title('GHGrem: ' + region + ' dTprime/dz')

    cmax=3e-6; cmin=-3e-6
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(242)
    plt.contourf(tlats,zlevs,wbarreg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title(region + ' wbar')

    cmax=3e-6; cmin=-3e-6
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(243)
    plt.contourf(tlats,zlevs,wvelbarreg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title(region + ' wvelbar')

    cmax=3e-6; cmin=-3e-6
    incr = (cmax-cmin) / (cmlen)
    conts = np.arange(cmin,cmax+incr,incr)

    ax = fig.add_subplot(244)
    plt.contourf(tlats,zlevs,wisopbarreg,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
    plt.colorbar()
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    ax.invert_yaxis()
    ax.set_title(region + ' wisopbar')

    if printtofile:
        fig.savefig('dTprimedz_Wbarcomps_Sulf_GHGrem_subplotSHzm_' + region + '.pdf')

##### end if pig region


wprime = np.squeeze(np.mean(wprime,axis=2))
wprime2 = np.squeeze(np.mean(wprime2,axis=2))
wvelprime = np.squeeze(np.mean(wvelprime,axis=2))
wvelprime2 = np.squeeze(np.mean(wvelprime2,axis=2))
wisopprime = np.squeeze(np.mean(wisopprime,axis=2))
wisopprime2 = np.squeeze(np.mean(wisopprime2,axis=2))

print wprime.shape

tbar = tempc# already zonal meaned 
tmptbar=np.flipud(tbar)
print tmptbar.shape
tmpdtbar = np.diff(tmptbar,axis=0) # lower level minus upper level
tmpdtbar = np.flipud(tmpdtbar) # positive T gradient sfc to seafloor

dtbarsave = np.diff(tbar,axis=0)
#dtbar=tmpdtbar
dtbar=dtbarsave # negative T gradient sfc to seafloor globally -- original way, correct way

print dtbar.shape

# thickness of each layer
dzt = np.diff(zt/100.) # convert to m
wtrans = np.zeros((len(dzt),wprime.shape[1]))
wtrans2 = np.zeros((len(dzt),wprime2.shape[1]))
wtranswv = np.zeros((len(dzt),wprime.shape[1]))
wtranswv2 = np.zeros((len(dzt),wprime2.shape[1]))
wtranswi = np.zeros((len(dzt),wprime.shape[1]))
wtranswi2 = np.zeros((len(dzt),wprime2.shape[1]))


print dzt.shape
# calc heat transport (heating rate) for each level
for lii,dz in enumerate(dzt):
    #print 'ind: ' + str(lii) + ', dz: ' + str(dz)
    
    wtrans[lii,...] = wprime[lii,...]*(dtbar[lii,...]/dz)
    wtrans2[lii,...] = wprime2[lii,...]*(dtbar[lii,...]/dz)
    
    wtranswv[lii,...] = wvelprime[lii,...]*(dtbar[lii,...]/dz) # WVEL only
    wtranswv2[lii,...] = wvelprime2[lii,...]*(dtbar[lii,...]/dz)
    
    wtranswi[lii,...] = wisopprime[lii,...]*(dtbar[lii,...]/dz) # WISOP only
    wtranswi2[lii,...] = wisopprime2[lii,...]*(dtbar[lii,...]/dz)

    
print wtrans.shape
#print dzt

# <codecell>

printtofile=False # True

meshlats,meshdz = np.meshgrid(np.squeeze(tlat[:,1]),zt[1:]/100.)

cmlen=float(20)
ylims = (0,500)
xlims = (-77,-40)

fig = plt.figure()
fig.set_size_inches(14,8)

"""ax = fig.add_subplot(151)
plt.pcolor(meshlats,meshdz,tmpdtbar,cmap='blue2red_20')
plt.colorbar()
ax.set_ylim((0,3000))
ax.set_xlim((-80, -40))
ax.invert_yaxis()"""

cmax=.8; cmin=-.8
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

ax = fig.add_subplot(241)
#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
plt.contourf(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
ax.set_title('dTbar/dz')

cmax=1e-6; cmin=-1e-6
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

ax = fig.add_subplot(242)
#plt.pcolor(tlats,zlevs,wprime,cmap='blue2red_20',vmin=-1e-6,vmax=1e-6)
plt.contourf(tlats,zlevs,wprime,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
ax.set_title('sulfates: wprime')

cmax=3e-7; cmin=-3e-7
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

ax = fig.add_subplot(243)
#plt.pcolor(tlats,zlevs,wvelprime,cmap='blue2red_20',vmin=-3e-7,vmax=3e-7)
plt.contourf(tlats,zlevs,wvelprime,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
ax.set_title('sulfates: wvelprime')

cmax=3e-7; cmin=-3e-7
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

ax = fig.add_subplot(244)
#plt.pcolor(tlats,zlevs,wisopprime,cmap='blue2red_20',vmin=-3e-7,vmax=3e-7)
plt.contourf(tlats,zlevs,wisopprime,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
ax.set_title('sulfates: wisopprime')

# do GHGrem =================
cmax=1e-6; cmin=-1e-6
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

ax = fig.add_subplot(246)
plt.contourf(tlats,zlevs,wprime2,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
ax.set_title('GHGrem: wprime')

cmax=3e-7; cmin=-3e-7
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

ax = fig.add_subplot(247)
plt.contourf(tlats,zlevs,wvelprime2,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
ax.set_title('GHGrem: wvelprime')

cmax=3e-7; cmin=-3e-7
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

ax = fig.add_subplot(248)
plt.contourf(tlats,zlevs,wisopprime2,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
ax.set_title('GHGrem: wisopprime')

if printtofile:
    fig.savefig('dTbardz_Wcomps_Sulf_GHGrem_subplotSHzm.pdf')

# <headingcell level=3>

# Figure 10a: components of heat trans calc for sulfates/GHGrem runs: climo dT/dz; anomalous WVEL+WISOP; anom WVEL; anom WISOP

# <codecell>

# w*dTbar/dz

meshlats,meshdz = np.meshgrid(np.squeeze(tlat[:,1]),zt[1:]/100.)
cmap='blue2red_w20'
cmin=-2e-9; cmax=2e-9

cmlen=float(20)
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)



fig = plt.figure()
fig.set_size_inches(14,3)
ax = fig.add_subplot(121)
CF1 = plt.contourf(meshlats,meshdz,wtrans,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

ax.set_ylim((0,1000))
ax.invert_yaxis()
ax.set_xlim((-80,-40))
cbar = fig.colorbar(CF1)
ax.set_title(casenamep)


ax2 = fig.add_subplot(122)
CF2 = plt.contourf(meshlats,meshdz,wtrans2,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

ax2.set_ylim((0,1000))
ax2.invert_yaxis()
ax2.set_xlim((-80,-40))
cbar = fig.colorbar(CF2)
ax2.set_title(casenamep2)

# <headingcell level=3>

# Fig. 10b: This should be heating rate (K/s) at each level for sulfate engineering (left) and GHG removal (right)

# <codecell>

# w dTbar/dz ZOOM
s2day = 60*60*24
print s2day

cmin=-2e-4; cmax=2e-4
cmlen=float(20)
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)


fig = plt.figure()
fig.set_size_inches(14,3)
ax = fig.add_subplot(121)
CF1 = plt.contourf(meshlats,meshdz,wtrans*s2day,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

ax.set_ylim((0,1000))
ax.invert_yaxis()
ax.set_xlim((-80,-60))
ax.set_title(casenamep)
cbar = fig.colorbar(CF1)


ax2 = fig.add_subplot(122)
CF2 = plt.contourf(meshlats,meshdz,wtrans2*s2day,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

ax2.set_ylim((0,1000))
ax2.invert_yaxis()
ax2.set_xlim((-80,-60))
ax2.set_title(casenamep2)
cbar = fig.colorbar(CF2)


# <headingcell level=3>

# Fig. 11: same as Fig. 10 but in K/day and zoomed in a little bit more (having trouble getting rid of spurious data in lower left)

# <codecell>

printtofile=False #True

rho_sw=cnc.getNCvar(filenamec,'rho_sw')
cp_sw = cnc.getNCvar(filenamec,'cp_sw')
rhocp = 1e-1*cp_sw*rho_sw # [J/K/m^3]

# tile the layer thickness on T grid
dzttile = np.tile(dzt,(wtrans.shape[1],1))
dzttile = np.transpose(dzttile)

xlims=(-77,-50)
ylims=(0,800)

# for vertical heat trans W/m2
cmin=-.2; cmax=.2
cmlen=float(20)
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)


fig = plt.figure()
fig.set_size_inches(14,3)
ax = fig.add_subplot(121)
CF1 = plt.contourf(meshlats,meshdz,wtrans*rhocp*dzttile,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')

# ---- add TEMP contours
tcmin=-2; tcmax=8
cmlen=float(15)
incr = (tcmax-tcmin) / (cmlen)
tconts = np.arange(tcmin,tcmax+incr,incr)

CS1 = plt.contour(tlats,zlevs,tempc,vmin=tcmin,vmax=tcmax,levels=tconts,colors='0.7',linewidths=2)

# ---- Add MOC contours
# contours for MOC anomaly
contspd = np.arange(.1,5,.5) # pos
contsnd = np.arange(-5,-.1,.5) # neg
# total MOC
plotmoc = (totmocp[0,0,1,...]+totmocp[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced
plotmoc2 = (totmocp2[0,0,1,...]+totmocp2[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced


CSm = plt.contour(lats,levs,plotmoc,contspd,\
            colors='g',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='g',linestyles='solid')

ax.set_ylim(ylims)
ax.invert_yaxis()
ax.set_xlim(xlims)
ax.set_title(casenamep)
cbar = fig.colorbar(CF1)


ax2 = fig.add_subplot(122)
CF2 = plt.contourf(meshlats,meshdz,wtrans2*rhocp*dzttile,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
CS2 = plt.contour(tlats,zlevs,tempc,vmin=tcmin,vmax=tcmax,levels=tconts,colors='0.7',linewidths=2)

CSm2 = plt.contour(lats,levs,plotmoc2,contspd,\
            colors='g',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm2 = plt.contour(lats,levs,plotmoc2,contsnd,\
            colors='g',linestyles='solid')

ax2.set_ylim(ylims)
ax2.invert_yaxis()
ax2.set_xlim(xlims)
ax2.set_title(casenamep2)
cbar = fig.colorbar(CF2)

if printtofile:
    fig.savefig('MOCeul+eddanom_climoTEMPcont_vertheattrans_subplotSH_ylim' + str(ylims[1]) + '.pdf')

# <headingcell level=3>

# Fig. 12: Same as 10, but multiply by layer thickness to convert to W/m2 in each grid cell. Overlaid with climo TEMP contours (gray) and anomalous eulerian MOC (green)

# <codecell>

printtofile=False #True

# plot MOC contours over dTbar/dz
rho_sw=cnc.getNCvar(filenamec,'rho_sw')
cp_sw = cnc.getNCvar(filenamec,'cp_sw')
rhocp = 1e-1*cp_sw*rho_sw # [J/K/m^3]

# tile the layer thickness on T grid
dzttile = np.tile(dzt,(wtrans.shape[1],1))
dzttile = np.transpose(dzttile)

xlims=(-77,-50)
ylims=(0,800)

cmlen=float(20)


fig = plt.figure()
fig.set_size_inches(14,3)
ax = fig.add_subplot(121)
# dTbar/dz here ==============

cmax=.8; cmin=-.8
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
plt.contourf(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
#ax.set_title('dTbar/dz')

# ---- Add MOC contours
# contours for MOC anomaly
contspd = np.arange(.1,5,.5) # pos
contsnd = np.arange(-5,-.1,.5) # neg
# total MOC
plotmoc = (totmocp[0,0,1,...]+totmocp[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced
plotmoc2 = (totmocp2[0,0,1,...]+totmocp2[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced


CSm = plt.contour(lats,levs,plotmoc,contspd,\
            colors='g',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='g',linestyles='solid')

ax.set_ylim(ylims)
ax.invert_yaxis()
ax.set_xlim(xlims)
ax.set_title(casenamep)


ax2 = fig.add_subplot(122)
cmax=.8; cmin=-.8
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
plt.contourf(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
#ax.set_title('dTbar/dz')

CSm2 = plt.contour(lats,levs,plotmoc2,contspd,\
            colors='g',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm2 = plt.contour(lats,levs,plotmoc2,contsnd,\
            colors='g',linestyles='solid')

ax2.set_ylim(ylims)
ax2.invert_yaxis()
ax2.set_xlim(xlims)
ax2.set_title(casenamep2)

if printtofile:
    fig.savefig('MOCeul+eddanom_dTbardz_subplotSH_ylim' + str(ylims[1]) + '.pdf')

# <codecell>

# # ===================== paper ======
# #  Zonal mean TEMP plus total MOC contour

printtofile=False

# plot MOC contours over T anomaly
rho_sw=cnc.getNCvar(filenamec,'rho_sw')
cp_sw = cnc.getNCvar(filenamec,'cp_sw')
rhocp = 1e-1*cp_sw*rho_sw # [J/K/m^3]
ylim=1000
xlims=(-77,-50)
ylims=(0,ylim)

cmlen=float(20)

cmap='blue2red_w20' # @@@

fig = plt.figure()
#fig.set_size_inches(10,3)
fig.set_size_inches(14,3) # to match PIG region
ax = fig.add_subplot(121)
# dTbar/dz here ==============

cmax=.5; cmin=-.5
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
plt.contourf(tlats,zlevs,tempp-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
#ax.set_title('dTbar/dz')

# ---- Add MOC contours
# contours for MOC anomaly
contspd = np.arange(.1,5,.15) # pos
contsnd = np.arange(-5,-.1,.15) # neg
# total MOC
plotmoc = (totmocp[0,0,1,...]+totmocp[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced
plotmoc2 = (totmocp2[0,0,1,...]+totmocp2[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced


CSm = plt.contour(lats,levs,plotmoc,contspd,\
            colors='.3',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='.3',linestyles='solid')

ax.set_ylim(ylims)
ax.invert_yaxis()
ax.set_xlim(xlims)
ax.set_title(casenamep)


ax2 = fig.add_subplot(122)
cmax=.5; cmin=-.5
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
plt.contourf(tlats,zlevs,tempp2-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
#ax.set_title('dTbar/dz')

CSm2 = plt.contour(lats,levs,plotmoc2,contspd,\
            colors='.3',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm2 = plt.contour(lats,levs,plotmoc2,contsnd,\
            colors='.3',linestyles='solid')

ax2.set_ylim(ylims)
ax2.invert_yaxis()
ax2.set_xlim(xlims)
ax2.set_title(casenamep2)

if printtofile:
    fig.savefig('MOCeul+eddanom_TEMPanom_subplotSH_ylim' + str(ylims[1]) + 'xlim' + str(xlims[1]) + '_smcont.pdf')
# ==================== end paper fig ==========

# # ===================== paper ======
# #  Zonal mean TEMP --- NO MOC

printtofile=True

# plot MOC contours over T anomaly
rho_sw=cnc.getNCvar(filenamec,'rho_sw')
cp_sw = cnc.getNCvar(filenamec,'cp_sw')
rhocp = 1e-1*cp_sw*rho_sw # [J/K/m^3]

ylim=1000
xlims=(-77,-50)
ylims=(0,ylim)

cmlen=float(20)

cmap='blue2red_w20' # @@@

#fig = plt.figure()

fig,axs = plt.subplots(1,2,sharey=True)
ax=axs[0]
#fig.set_size_inches(10,3)
fig.set_size_inches(14,3) # to match PIG region
#ax = fig.add_subplot(121)
# dTbar/dz here ==============

cmax=.5; cmin=-.5
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
CF = ax.contourf(tlats,zlevs,tempp-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
#cbar = fig.colorbar(CF)

"""# ---- Add MOC contours
# contours for MOC anomaly
contspd = np.arange(.1,5,.15) # pos
contsnd = np.arange(-5,-.1,.15) # neg
# total MOC
plotmoc = (totmocp[0,0,1,...]+totmocp[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced
plotmoc2 = (totmocp2[0,0,1,...]+totmocp2[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced


CSm = plt.contour(lats,levs,plotmoc,contspd,\
            colors='.3',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='.3',linestyles='solid')
"""
ax.set_ylim(ylims)
ax.invert_yaxis()
ax.set_yticks(np.arange(0,900,100))
ax.set_yticklabels([0,'',200,'',400,'',600,'',800],fontsize=18)
ax.set_xlim(xlims)
ax.set_xticks(np.arange(-75,-45,5))
ax.set_xticklabels(['75$^\circ$S', '70$^\circ$S', '65$^\circ$S', \
                    '60$^\circ$S', '55$^\circ$S', '50$^\circ$S'],fontsize=18)
#ax.set_title(casenamep)
ax.set_title('Sulf',fontsize=18)
ax.set_ylabel('Depth (m)',fontsize=18)

#ax2 = fig.add_subplot(122)
ax2=axs[1]
cmax=.5; cmin=-.5
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
CF2 = ax2.contourf(tlats,zlevs,tempp2-tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
#cbar = fig.colorbar(CF2)

#CSm2 = plt.contour(lats,levs,plotmoc2,contspd,\
#            colors='.3',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
#CSm2 = plt.contour(lats,levs,plotmoc2,contsnd,\
#            colors='.3',linestyles='solid')

ax2.set_ylim(ylims)
ax2.invert_yaxis()
ax2.set_yticks(np.arange(0,900,100))
ax2.set_yticklabels([0,'',200,'',400,'',600,'',800],fontsize=18)
ax2.set_xlim(xlims)
ax2.set_xticks(np.arange(-75,-45,5))
ax2.set_xticklabels(['75$^\circ$S', '70$^\circ$S', '65$^\circ$S', \
                    '60$^\circ$S', '55$^\circ$S', '50$^\circ$S'],fontsize=18)
#ax2.set_title(casenamep2)
ax2.set_title('GHGrem',fontsize=18)
cbar_ax = fig.add_axes([.91,.15, .02,.7])
fig.colorbar(CF,cax=cbar_ax)

if printtofile:
    fig.savefig('TEMPanom_subplotSH_ylim' + str(ylims[1]) + 'xlim' + str(xlims[1]) + '_c.pdf')
# ==================== end paper fig ==========
# =============================


# ====== Test Figure ===================
# ======= Climo Temp with total MOC anom =========
# <codecell>

printtofile=False

# plot MOC contours over T
rho_sw=cnc.getNCvar(filenamec,'rho_sw')
cp_sw = cnc.getNCvar(filenamec,'cp_sw')
rhocp = 1e-1*cp_sw*rho_sw # [J/K/m^3]


xlims=(-77,-50)
ylims=(0,800)

cmlen=float(20)

fig = plt.figure()
fig.set_size_inches(10,3)
ax = fig.add_subplot(121)
# dTbar/dz here ==============

cmax=5; cmin=-2
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)
cmap='jet'

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
plt.contourf(tlats,zlevs,tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
#ax.set_title('dTbar/dz')

# ---- Add MOC contours
# contours for MOC anomaly
contspd = np.arange(.1,5,.2) # pos
contsnd = np.arange(-5,-.1,.2) # neg

# total MOC
plotmoc = (totmocp[0,0,1,...]+totmocp[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced
plotmoc2 = (totmocp2[0,0,1,...]+totmocp2[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced


CSm = plt.contour(lats,levs,plotmoc,contspd,\
            colors='.3',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='.3',linestyles='solid')

ax.set_ylim(ylims)
ax.invert_yaxis()
ax.set_xlim(xlims)
ax.set_title(casenamep)


ax2 = fig.add_subplot(122)
cmax=5; cmin=-2
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
plt.contourf(tlats,zlevs,tempc,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
#ax.set_title('dTbar/dz')

CSm2 = plt.contour(lats,levs,plotmoc2,contspd,\
            colors='.3',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm2 = plt.contour(lats,levs,plotmoc2,contsnd,\
            colors='.3',linestyles='solid')

ax2.set_ylim(ylims)
ax2.invert_yaxis()
ax2.set_xlim(xlims)
ax2.set_title(casenamep2)

if printtofile:
    fig.savefig('MOCeul+eddyanom_TEMPclimo_subplotSH_ylim' + str(ylims[1]) + '_2smcont.pdf')


# === test fig 2: dTbar / dz with total MOC contours ===

printtofile=False

# plot MOC contours over T
rho_sw=cnc.getNCvar(filenamec,'rho_sw')
cp_sw = cnc.getNCvar(filenamec,'cp_sw')
rhocp = 1e-1*cp_sw*rho_sw # [J/K/m^3]


xlims=(-77,-50)
ylims=(0,800)

cmlen=float(20)

fig = plt.figure()
fig.set_size_inches(10,3)
ax = fig.add_subplot(121)
# dTbar/dz here ==============

cmax=0.8; cmin=-0.8
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)
cmap='blue2red_20'

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
plt.contourf(meshlats,meshdz,dtbar,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both')
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
#ax.set_title('dTbar/dz')

# ---- Add MOC contours
# contours for MOC anomaly
contspd = np.arange(.1,5,.2) # pos
contsnd = np.arange(-5,-.1,.2) # neg

# total MOC
plotmoc = (totmocp[0,0,1,...]+totmocp[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced
plotmoc2 = (totmocp2[0,0,1,...]+totmocp2[0,0,0,...])-(totmocc[0,0,1,...]+totmocc[0,0,0,...]) # Eulerian+eddy-induced


CSm = plt.contour(lats,levs,plotmoc,contspd,\
            colors='.3',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm = plt.contour(lats,levs,plotmoc,contsnd,\
            colors='.3',linestyles='solid')

ax.set_ylim(ylims)
ax.invert_yaxis()
ax.set_xlim(xlims)
ax.set_title(casenamep)

ax2 = fig.add_subplot(122)
cmax=0.8; cmin=-0.8
incr = (cmax-cmin) / (cmlen)
conts = np.arange(cmin,cmax+incr,incr)

#plt.pcolor(meshlats,meshdz,dtbarsave,cmap='blue2red_20',vmin=-1,vmax=1)
plt.contourf(meshlats,meshdz,dtbar,cmap=cmap,vmin=cmin,vmax=cmax,levels=conts,extend='both') # $$$
plt.colorbar()
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.invert_yaxis()
#ax.set_title('dTbar/dz')

CSm2 = plt.contour(lats,levs,plotmoc2,contspd,\
            colors='.3',linestyles='solid',linewidths=2)
#plt.clabel(CS1,fmt = '%2.1f',inline=1,fontsize=10)
CSm2 = plt.contour(lats,levs,plotmoc2,contsnd,\
            colors='.3',linestyles='solid')

ax2.set_ylim(ylims)
ax2.invert_yaxis()
ax2.set_xlim(xlims)
ax2.set_title(casenamep2)

if printtofile:
    fig.savefig('MOCeul+eddyanom_dTbar_subplotSH_ylim' + str(ylims[1]) + '_2smcont.pdf')



# ============ end test fig ============


# <codecell>

import cccmacmaps as ccm

printtofile=False #True

sec2yr = s2day*365

mediumblue = ccm.get_linecolor('mediumblue') # Sulf
dodgerblue = ccm.get_linecolor('dodgerblue') # GHGrem
darkolivegreen3 = ccm.get_linecolor('darkolivegreen3') # GHGrem
firebrick = ccm.get_linecolor('firebrick') # RCP8.5

                  
onelat = tlat[:,1]
totw = wtrans*rhocp*dzttile
totw2 = wtrans2*rhocp*dzttile
totwv = wtranswv*rhocp*dzttile
totwv2 = wtranswv2*rhocp*dzttile
totwi = wtranswi*rhocp*dzttile
totwi2 = wtranswi2*rhocp*dzttile

if pig:
    totwreg = wtransreg*rhocp*dzttile
    totw2reg = wtrans2reg*rhocp*dzttile
    totwvreg = wtranswvreg*rhocp*dzttile
    totwv2reg = wtranswv2reg*rhocp*dzttile
    totwireg = wtranswireg*rhocp*dzttile
    totwi2reg = wtranswi2reg*rhocp*dzttile

    totwbarreg = wbartransreg*rhocp*dzttile
    totwbar2reg = wbartrans2reg*rhocp*dzttile
    totwbarwvreg = wbartranswvreg*rhocp*dzttile
    totwbarwv2reg = wbartranswv2reg*rhocp*dzttile
    totwbarwireg = wbartranswireg*rhocp*dzttile
    totwbarwi2reg = wbartranswi2reg*rhocp*dzttile


"""totwL450 = np.squeeze(totw[30,...]) # 408m depth using zlevs
totw2L450 = np.squeeze(totw2[30,...])
totwvL450 = np.squeeze(totwv[30,...]) 
totwv2L450 = np.squeeze(totwv2[30,...]) 
totwiL450 = np.squeeze(totwi[30,...]) 
totwi2L450 = np.squeeze(totwi2[30,...]) """

ylim=500
#dep=19; pigvhylims=[-.8,.5] # 197m depth using zlevs (in meters already)
#dep=15 # 155m depth using zlevs
dep=28; pigvhylims=[-.2,.2] # 351.09347534179688 depth using zlevs

totwL1 = np.squeeze(wtrans[dep,...]) 
totw2L1 = np.squeeze(wtrans2[dep,...])
totwvL1 = np.squeeze(wtranswv[dep,...]) 
totwv2L1 = np.squeeze(wtranswv2[dep,...]) 
totwiL1 = np.squeeze(wtranswi[dep,...]) 
totwi2L1 = np.squeeze(wtranswi2[dep,...])

"""dep=22 #
totwL450 = np.squeeze(wtrans[dep,...]) # 236m depth using zlevs
totw2L450 = np.squeeze(wtrans2[dep,...])
totwvL450 = np.squeeze(wtranswv[dep,...]) 
totwv2L450 = np.squeeze(wtranswv2[dep,...]) 
totwiL450 = np.squeeze(wtranswi[dep,...]) 
totwi2L450 = np.squeeze(wtranswi2[dep,...])

dep=26 # 
totwL200 = np.squeeze(totw[dep,...]) # 305m depth using zlevs
totw2L200 = np.squeeze(totw2[dep,...])
totwvL200 = np.squeeze(totwv[dep,...]) 
totwv2L200 = np.squeeze(totwv2[dep,...]) 
totwiL200 = np.squeeze(totwi[dep,...]) 
totwi2L200 = np.squeeze(totwi2[dep,...]) """

# trying to avoid cells too far south 
Nlim=-65   # 65 to 74 is good for pig
Slim=-74

tareay = tarea[:,0] # because we are dealing w/ SH only, doesn't matter what lon we choose
# now tile tareay for each depth
tareayt = np.tile(tareay,(len(zt)-1,1)) # for transport (or when a dt or dz is involved)
# create weights:
totareay = np.sum(tareayt[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
totareayt = np.tile(totareay,(len(tareayt[0,np.logical_and(onelat<= Nlim,onelat>Slim)]),1))
totareayt = np.transpose(totareayt,(1,0))
wgts = tareayt[:,np.logical_and(onelat<= Nlim,onelat>Slim)] / totareayt

tareaytw = np.tile(tareay,(len(zt),1)) # for all 60 levels
totareayw = np.sum(tareaytw[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
totareaytw = np.tile(totareayw,(len(tareaytw[0,np.logical_and(onelat<= Nlim,onelat>Slim)]),1))
totareaytw = np.transpose(totareaytw,(1,0))
wgtsw = tareaytw[:,np.logical_and(onelat<= Nlim,onelat>Slim)] / totareaytw


# try area-weighting:
#  vertical heat trans
totw = np.average(totw[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
totw2 = np.average(totw2[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
totwv = np.average(totwv[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
totwv2 = np.average(totwv2[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
totwi = np.average(totwi[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
totwi2 = np.average(totwi2[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)

#  vertical velocity
wavg=np.average(wprime[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
wavg2=np.average(wprime2[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
wvavg=np.average(wvelprime[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
wvavg2=np.average(wvelprime2[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
wiavg=np.average(wisopprime[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
wiavg2=np.average(wisopprime2[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)

if pig:
# # =================== paper =======
# #  PIG vertical heat trans and velocity

    printtofile=False

    #ylim=450
    ylim=1000
    totwL1reg = np.squeeze(wtransreg[dep,...]) 
    totw2L1reg = np.squeeze(wtrans2reg[dep,...])
    totwvL1reg = np.squeeze(wtranswvreg[dep,...]) 
    totwv2L1reg = np.squeeze(wtranswv2reg[dep,...]) 
    totwiL1reg = np.squeeze(wtranswireg[dep,...]) 
    totwi2L1reg = np.squeeze(wtranswi2reg[dep,...])

    # meridional average with depth: vertical heat transport
    totwreg = np.average(totwreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totw2reg = np.average(totw2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwvreg = np.average(totwvreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwv2reg = np.average(totwv2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwireg = np.average(totwireg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwi2reg = np.average(totwi2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)

    # meridional average WBAR*dTPRIME with depth: vertical heat transport
    totwbarreg = np.average(totwbarreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwbar2reg = np.average(totwbar2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwbarwvreg = np.average(totwbarwvreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwbarwv2reg = np.average(totwbarwv2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwbarwireg = np.average(totwbarwireg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    totwbarwi2reg = np.average(totwbarwi2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)

    # meridional average w PRIME with depth
    wavgreg=np.average(wprimereg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    wavg2reg=np.average(wprime2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    wvavgreg=np.average(wvelprimereg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    wvavg2reg=np.average(wvelprime2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    wiavgreg=np.average(wisopprimereg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    wiavg2reg=np.average(wisopprime2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)

    # meridional average w BAR with depth
    wbaravgreg = np.average(wbarreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    wbarwvavgreg = np.average(wvelbarreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)
    wbarwiavgreg = np.average(wisopbarreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgtsw)

    
    # climo dTbar, averaged meridionally
    dTbaravgreg=np.average(dtbarreg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)

    # anom dTprime, averaged meridionally
    dTprimeavgreg = np.average(dtprimereg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)
    dTprime2avgreg = np.average(dtprime2reg[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1,weights=wgts)


    fig2 = plt.figure()
    fig2.set_size_inches(14,4)
    ax = fig2.add_subplot(131)

    plt.plot(totwreg,zt[1:]/100.,color=mediumblue,linewidth=3)
    plt.plot(totw2reg,zt[1:]/100.,color=darkolivegreen3,linewidth=3) # GHGrem
    plt.plot(totwvreg,zt[1:]/100.,color=mediumblue,linestyle='--')
    plt.plot(totwv2reg,zt[1:]/100.,color=darkolivegreen3,linestyle='--')
    plt.plot(totwireg,zt[1:]/100.,color=mediumblue)
    plt.plot(totwi2reg,zt[1:]/100.,color=darkolivegreen3)

    plt.plot([0,0],[0,1000],'k')
    plt.legend(('sulfate','ghgrem'),loc='best')
    plt.ylim((0,ylim))
    plt.xlim(-.15,.15)
    #plt.xticks([-0.15,-0.10,-0.05, 0, 0.05, 0.1, 0.15])
    #ax.set_xticklabels([-0.15,'',-0.05, 0, .05,'', 0.15], fontsize=16)
    plt.yticks([0,100,200,300,400])
    ax.set_yticklabels([0,100,200,300,400],fontsize=16)
    #plt.title(region + ' Avg vert heat trans (W/m2) ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
    plt.title('VHT ($W/m^2$)')
    ax.invert_yaxis()

    ax2=fig2.add_subplot(132)
    plt.plot(wavgreg*sec2yr,zt/100.,color=mediumblue,linewidth=3)
    plt.plot(wavg2reg*sec2yr,zt/100.,color=darkolivegreen3,linewidth=3)
    plt.plot(wvavgreg*sec2yr,zt/100.,color=mediumblue,linestyle='--')
    plt.plot(wvavg2reg*sec2yr,zt/100.,color=darkolivegreen3,linestyle='--')
    plt.plot(wiavgreg*sec2yr,zt/100.,color=mediumblue)
    plt.plot(wiavg2reg*sec2yr,zt/100.,color=darkolivegreen3)

    plt.plot([0,0],[0,1000],'k')
    plt.ylim((0,ylim))
    plt.xlim(-10,8)
    plt.xticks([-10,-8,-6,-4,-2, 0, 2, 4, 6, 8, 10])
    ax2.set_xticklabels([-10,'',-6,'',-2,0,2, '',6,'',10], fontsize=16)
    plt.yticks([0,100,200,300,400])
    ax2.set_yticklabels([0,100,200,300,400],fontsize=16)

    plt.title('Vert velocity (m/yr) ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
    ax2.invert_yaxis()

    plottotw = totwL1reg; plottotw2 = totw2L1reg
    plottotwv = totwvL1reg; plottotwv2 = totwv2L1reg
    plottotwi = totwiL1reg; plottotwi2 = totwi2L1reg


    ax3=fig2.add_subplot(133)
    plt.plot(onelat,plottotw*sec2yr,color=mediumblue,linewidth=3)
    plt.plot(onelat,plottotwv*sec2yr,color=mediumblue,linestyle='--')
    plt.plot(onelat,plottotwi*sec2yr,color=mediumblue)

    plt.plot(onelat,plottotw2*sec2yr,color=darkolivegreen3,linewidth=3)
    plt.plot(onelat,plottotwv2*sec2yr,color=darkolivegreen3,linestyle='--')
    plt.plot(onelat,plottotwi2*sec2yr,color=darkolivegreen3)

    plt.plot((Slim,Nlim),(0,0),'k')

    ax3.set_title(str(zlevs[dep,1]) + 'm heat trans')
    ax3.set_xlabel('Lat')
    ax3.set_xlim((Slim,Nlim))
    #ax3.set_ylim((-.8,.5))
    ax3.set_ylim(pigvhylims)

    if printtofile:
        fig2.savefig('vertheattrans_wvels_' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 
                     'S_lev' + str(np.round(zlevs[dep,1])) + '_ylim' + str(ylim) + '_' + region + '_b.pdf')

    printtofile=True

    # VERSION 2 has dTbar instead of vert heat through a layer ======= PAPER
    #fig2 = plt.figure()
    fig2,axs = plt.subplots(1,3,sharey=True)
    fig2.set_size_inches(14,4)
    ax = axs[0] #fig2.add_subplot(131,sharey=True)

    ax.plot(-1*totwreg,zt[1:]/100.,color=mediumblue,linewidth=4)
    ax.plot(-1*totw2reg,zt[1:]/100.,color=darkolivegreen3,linewidth=4) # GHGrem
    ax.plot(-1*totwvreg,zt[1:]/100.,color=mediumblue,linewidth=2,linestyle='--')
    ax.plot(-1*totwv2reg,zt[1:]/100.,color=darkolivegreen3,linewidth=2,linestyle='--')
    ax.plot(-1*totwireg,zt[1:]/100.,color=mediumblue,linewidth=2)#,linestyle=':')
    ax.plot(-1*totwi2reg,zt[1:]/100.,color=darkolivegreen3,linewidth=2)#3,linestyle=':')

    yticks=np.arange(0,ylim,100)
    ax.plot([0,0],[0,1000],'k')
    ax.legend(('Sulf','GHGrem'),loc='best',fancybox=True,framealpha=0.5)
    ax.set_ylim((0,ylim))
    ax.set_xlim(-.15,.15)
    ax.set_xticks([-0.15,-0.10,-0.05, 0, 0.05, 0.1, 0.15])
    ax.set_xticklabels([-0.15,'',-0.05, 0, .05,'', 0.15], fontsize=18)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks,fontsize=18)
    #ax.set_title('VHT (W/m$^2$)',fontsize=18)
    ax.set_title('w$^{\prime}$*d$\overline{T}$/dz (W/m$^2$)',fontsize=18)
    ax.set_ylabel('Depth (m)',fontsize=18)
    #plt.title(region + ' Avg vert heat trans (W/m2) ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
    ax.invert_yaxis()

    ax2=axs[1] #fig2.add_subplot(132,sharey=True)
    ax2.plot(wavgreg*sec2yr,zt/100.,color=mediumblue,linewidth=4)
    ax2.plot(wavg2reg*sec2yr,zt/100.,color=darkolivegreen3,linewidth=4)
    ax2.plot(wvavgreg*sec2yr,zt/100.,color=mediumblue,linewidth=2,linestyle='--')
    ax2.plot(wvavg2reg*sec2yr,zt/100.,color=darkolivegreen3,linewidth=2,linestyle='--')
    ax2.plot(wiavgreg*sec2yr,zt/100.,color=mediumblue,linewidth=2)#,linestyle='-.')
    ax2.plot(wiavg2reg*sec2yr,zt/100.,color=darkolivegreen3,linewidth=2)#,linestyle='-.')

    ax2.plot([0,0],[0,1000],'k')
    ax2.set_ylim((0,ylim))
    ax2.set_xlim(-10,8)
    ax2.set_xticks([-10,-8,-6,-4,-2, 0, 2, 4, 6, 8, 10])
    ax2.set_xticklabels([-10,'',-6,'',-2,0,2, '',6,'',10], fontsize=18)
    #ax2.set_yticks([0,100,200,300,400])
    #ax2.set_yticklabels([0,100,200,300,400],fontsize=18)
    #plt.title('Vert velocity (m/yr) ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
    ax2.set_title('w$^{\prime}$ (m/yr)',fontsize=18)
    ax2.invert_yaxis()

    ax3=axs[2] #fig2.add_subplot(133)
#    # POSITIVE means getting warmer with depth
#    plt.plot(dTbaravgreg,zt[1:]/100.,color='k',linewidth=3)
    # NEGATIVE means getting warmer with depth (if mult by -1)
    ax3.plot(-1*dTbaravgreg,zt[1:]/100.,color='k',linewidth=3) 

    ax3.plot([0,0],[0,1000],'k')
    ax3.set_ylim((0,ylim))
    #plt.xlim(-10,8)
    ax3.set_xticks([-0.3,-0.2,-0.1,0,0.1]) # for when dT/dz is neg for warm with depth
    ax3.set_xticklabels([-0.3,-0.2,-0.1,0,0.1], fontsize=18)
    #ax3.set_xticks([-0.1,0,0.1,0.2,0.3]) # for when dT/dz is pos for warm wiht depth
    #ax3.set_xticklabels([-0.1,0,0.1,0.2,0.3], fontsize=18)

    #ax3.set_yticks([0,100,200,300,400])
    #ax3.set_yticklabels([0,100,200,300,400],fontsize=18)

    #ax3.set_title('dTbar ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
    ax3.set_title('d$\overline{T}$/dz ($^\circ$C/m)',fontsize=18)
    ax3.invert_yaxis()

    if printtofile:
        fig2.savefig('vertheattrans_wvels_dTbarnegwrmwithdepth_' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 
                     'S_ylim' + str(ylim) + '_' + region + '_b.pdf')


    # VERSION 3 is wbar * dTprime
    printtofile=False

    fig2 = plt.figure()
    fig2.set_size_inches(14,4)
    ax = fig2.add_subplot(131)

    plt.plot(totwbarreg,zt[1:]/100.,color=mediumblue,linewidth=2)
    plt.plot(totwbar2reg,zt[1:]/100.,color=dodgerblue,linewidth=2) # GHGrem
    plt.plot(totwbarwvreg,zt[1:]/100.,color=mediumblue,linestyle='--')
    plt.plot(totwbarwv2reg,zt[1:]/100.,color=dodgerblue,linestyle='--')
    plt.plot(totwbarwireg,zt[1:]/100.,color=mediumblue)
    plt.plot(totwbarwi2reg,zt[1:]/100.,color=dodgerblue)

    plt.plot([0,0],[0,1000],'k')
    plt.legend(('sulfate','ghgrem'),loc='best')
    plt.ylim((0,ylim))
    plt.xlim(-15,15)
    plt.title(region + ' Avg vert heat trans wbar*dTprime (W/m2) ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
    ax.invert_yaxis()

    ax2=fig2.add_subplot(132)
    plt.plot(wbaravgreg*sec2yr,zt/100.,color='k',linewidth=2)
    plt.plot(wbarwvavgreg*sec2yr,zt/100.,color='k',linestyle='--')
    plt.plot(wbarwiavgreg*sec2yr,zt/100.,color='k')

    plt.plot([0,0],[0,1000],'k')
    plt.ylim((0,ylim))
    #@@plt.xlim(-10,8)
    plt.title('Climo w (m/yr) ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
    ax2.invert_yaxis()

    ax3=fig2.add_subplot(133)
    # POSITIVE means getting warmer with depth
    plt.plot(dTprimeavgreg,zt[1:]/100.,color=mediumblue,linewidth=2)
    plt.plot(dTprime2avgreg,zt[1:]/100.,color=dodgerblue,linewidth=2)

    plt.plot([0,0],[0,1000],'k')
    plt.ylim((0,ylim))
    #plt.xlim(-10,8)
    plt.title('dTprime ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
    ax3.invert_yaxis()

    if printtofile:
        fig2.savefig('vertheattrans_wbarvels_dTprime_' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 
                     'S_ylim' + str(ylim) + '_' + region + '.pdf')


""" comment out non-area weighted @@
totw = np.mean(totw[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
totw2 = np.mean(totw2[:,np.logical_and(onelat<=Nlim,onelat>Slim)],axis=1)
totwv = np.mean(totwv[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
totwv2 = np.mean(totwv2[:,np.logical_and(onelat<=Nlim,onelat>Slim)],axis=1)
totwi = np.mean(totwi[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
totwi2 = np.mean(totwi2[:,np.logical_and(onelat<=Nlim,onelat>Slim)],axis=1)

wavg=np.mean(wprime[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
wavg2=np.mean(wprime2[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
wvavg=np.mean(wvelprime[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
wvavg2=np.mean(wvelprime2[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
wiavg=np.mean(wisopprime[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
wiavg2=np.mean(wisopprime2[:,np.logical_and(onelat<= Nlim,onelat>Slim)],axis=1)
"""

# # === full zonal vertical heat trans and velocity ====

printtofile=False # True

fig2 = plt.figure()
fig2.set_size_inches(14,4)
ax = fig2.add_subplot(131)

plt.plot(totw,zt[1:]/100.,color=mediumblue,linewidth=2)
plt.plot(totw2,zt[1:]/100.,color=dodgerblue,linewidth=2) # GHGrem
plt.plot(totwv,zt[1:]/100.,color=mediumblue,linestyle='--')
plt.plot(totwv2,zt[1:]/100.,color=dodgerblue,linestyle='--')
plt.plot(totwi,zt[1:]/100.,color=mediumblue)
plt.plot(totwi2,zt[1:]/100.,color=dodgerblue)

plt.plot([0,0],[0,1000],'k')
plt.legend(('sulfate','ghgrem'))
plt.ylim((0,ylim))
plt.xlim(-.15,.15)
plt.title('Avg vert heat trans (W/m2) ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
ax.invert_yaxis()


ax2=fig2.add_subplot(132)
plt.plot(wavg*sec2yr,zt/100.,color=mediumblue,linewidth=2)
plt.plot(wavg2*sec2yr,zt/100.,color=dodgerblue,linewidth=2)
plt.plot(wvavg*sec2yr,zt/100.,color=mediumblue,linestyle='--')
plt.plot(wvavg2*sec2yr,zt/100.,color=dodgerblue,linestyle='--')
plt.plot(wiavg*sec2yr,zt/100.,color=mediumblue)
plt.plot(wiavg2*sec2yr,zt/100.,color=dodgerblue,)

plt.plot([0,0],[0,1000],'k')
#plt.legend(('sulfate','ghgrem'))
plt.ylim((0,ylim))
plt.xlim(-5,5)
plt.title('Vert velocity (m/yr) ' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 'S at e/ lev')
ax2.invert_yaxis()

plottotw = totwL1; plottotw2 = totw2L1
plottotwv = totwvL1; plottotwv2 = totwv2L1
plottotwi = totwiL1; plottotwi2 = totwi2L1


ax3=fig2.add_subplot(133)
plt.plot(onelat,plottotw*sec2yr,color=mediumblue,linewidth=2)
plt.plot(onelat,plottotwv*sec2yr,color=mediumblue,linestyle='--')
plt.plot(onelat,plottotwi*sec2yr,color=mediumblue)

plt.plot(onelat,plottotw2*sec2yr,color=dodgerblue,linewidth=2)
plt.plot(onelat,plottotwv2*sec2yr,color=dodgerblue,linestyle='--')
plt.plot(onelat,plottotwi2*sec2yr,color=dodgerblue)

ax3.set_title(str(zlevs[dep,1]) + 'm heat trans')
ax3.set_xlabel('Lat')
ax3.set_xlim((Slim,Nlim))
ax3.set_ylim((-.1,.1))

if printtofile:
    fig2.savefig('vertheattrans_wvels_' + str(np.abs(Slim)) + 'S-' + str(np.abs(Nlim)) + 
                 'S_lev' + str(zlevs[dep,1]) + '_ylim' + str(ylim) + '.pdf')

# <headingcell level=3>

# Fig. 13: Average vertical heat trans (W/m2) between given lats, at each depth. Thick line is using WVEL+WISOP anom to calculate, dashed line is using WVEL only, thin solid is using WISOP only

# <codecell>


# <rawcell>

# 
# 
# 
# ====== Not done below ======
# Now want to show zonal mean net surface flux into the ocean

# <codecell>

#include only ocean grid cells
#zmtemp = plotfld.*ofrac./1;
#zmfld = squeeze(mean(zmtemp,2));

shfc = cnc.getNCvar(filenamec,'SHF')
shfp = cnc.getNCvar(filenamep,'SHF')
shfp2 = cnc.getNCvar(filenamep2,'SHF')

tlon = cnc.getNCvar(filenamec,'TLONG')
tlon.shape
tarea=cnc.getNCvar(filenamec,'TAREA')
kmt=cnc.getNCvar(filenamec,'KMT') # k Index of Deepest Grid Cell on T Grid" . index of 0 means land at sfc

