"""
cccmaplots.py
Use this module to house specific plotting functions
    and colormaps.
    e.g. kemmap()

    2/11/2014

"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits as mpltk
import copy
import constants as con
import cccmacmaps as ccm


"""
kemmap(fld, lat, lon, title='', units='', cmap='blue2red_w20', type='sq', cmin='', cmax='',
       axis=None, suppcb=0,lmask=0,flipmask=0,latlim=None)

    Inputs: fld: 2D matrix of data [lat x lon]
            lat: 1D lat array
            lon: 1D lon array
            title: plot title
            units: data units (colorbar label)
            cmap: map colormap. default blue to red, continuous
            type: type of map. default 'sq' for Robinson projection. 
                  else: 'nh', 'sh' orthographic projections
            cmin, cmax: color scale limits
            axis: handle to the axes instance to give to basemap. default None.
            suppcb: suppress colorbar. default = 0 (do not suppress)
            lmask: add land mask (e.g. for ocean-only data). default is no mask (0)
            flipmask: a hack to flip the lat direction of lmask (e.g. for HURRELL data)
            latlim: if polar projection, it will be the equatorward limit (so far only
                    works if do a stereographic proj, not ortho 5/13/2014

    Returns: basemap handle (to add to bm after function call),
             pcolormesh handle (typically to add colorbar after func call)

"""

def kemmap(fld, lat, lon,title='',units='',cmap='blue2red_w20',type='sq',
           cmin='',cmax='',axis=None, suppcb=0,lmask=0,flipmask=0,latlim=None):

    if cmap =='' or cmap==None:
        cmap='blue2red_w20'
        
    incmap = plt.cm.get_cmap(cmap)

    
    # default Basemap dictionary
    if type == 'sq':
        mapparams = dict(projection='robin',lon_0=180,lat_0=0, resolution='c')
    elif type == 'nh':
        if latlim != None:
            mapparams = dict(projection='npstere',boundinglat=latlim,lon_0=0,resolution='c')
        else:
            # try mill, hammer, merc
            mapparams = dict(projection='ortho',lon_0=0.,lat_0=90.,\
                             resolution='c') #llcrnrlon='-180',llcrnrlat='45',urcrnrlon='180',urcrnrlat='90'
        # I thought the above corner limits would work to zoom on NH but I'm getting
        # AttributeError: 'Basemap' object has no attribute '_height'
        # 5/12/14 -- don't know why. same goes for lat_0=0.
    elif type == 'sh':
        if latlim != None:
            mapparams = dict(projection='spstere',boundinglat=latlim,lon_0=0,resolution='c')
        else:
            mapparams = dict(projection='ortho',lon_0=0.,lat_0=-90., resolution='c')
            # same error if add: llcrnrlon='-180',llcrnrlat='-90',urcrnrlon='180',urcrnrlat='-45'
    else:
        print "Incorrect map type. Choose sq,nh,sh"
        exit
        
    # default pcolormesh dictionary
    if cmin =='':
        #pcparams = dict(shading='gouraud',latlon=True,cmap=incmap)
        pcparams = dict(latlon=True,cmap=incmap)
    else:
        cmlen=float(incmap.N) # or: from __future__ import division

        if cmlen>200:
            cmlen = float(20) # ie. if using a built in colormap that is continuous, want smaller # of colors
        
        incr = (cmax-cmin) /cmlen
        conts = np.arange(cmin,cmax+incr,incr)

        #pcparams = dict(shading='gouraud',latlon=True,cmap=incmap,vmin=cmin,vmax=cmax)
        pcparams = dict(latlon=True,cmap=incmap,levels=conts,vmin=cmin,vmax=cmax,extend='both')

    if axis != None: # if an axis is given, add to dict for basemap
        mapparams['ax'] = axis

    """ m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution='l',\
    llcrnrx=0.,llcrnry=0.,urcrnrx=m1.urcrnrx/2.,urcrnry=m1.urcrnry/2.)
    """

    bm = Basemap(**mapparams)
        
    if np.mod(lon.shape,2) == 0:
        # add cyclic lon
        fld,lon = mpltk.basemap.addcyclic(fld,lon)
     
    if lmask==1: # the landmask has an extra lon already
        # add land mask
        lsmask=con.get_t63landmask()
        if flipmask:
            lsmask=np.flipud(lsmask)
            
        fld = ma.masked_where(lsmask!=0,fld) # 0 is ocean
        
    lons, lats = np.meshgrid(lon,lat)
        

#    pc = bm.pcolormesh(lons,lats,fld,**pcparams)
    if cmin=='':
        pc = bm.contourf(lons,lats,fld,20,**pcparams) # 20 auto levels
    else:
        pc = bm.contourf(lons,lats,fld,**pcparams)
    bm.drawcoastlines(color='0.7')
    #bm.drawmapboundary(fill_color='#99ffff')

    if lmask==1:
        bm.drawmapboundary(fill_color='0.7')
        bm.drawcoastlines(color='.3')

    # I think drawlsmask puts the mask on the bottom
    #bm.drawlsmask(land_color='0.7',lsmask=con.get_t63landmask(),ax=axis)
    #bm.drawparallels(np.arange(-90.,120.,30.))
    #bm.drawmeridians(np.arange(0.,360.,30.))
    if axis!=None:
        axis.set_title(title,fontsize=10)
    else:
        plt.title(title,fontsize=10)

    # add colorbar.
    if suppcb == 0:
        cbar = bm.colorbar(pc,location='bottom',pad="5%")
        cbar.set_label(units)

    return bm,pc


def addtsigm(basem, pvals, lat, lon, siglevel=0.05,color='k',type='hatch'):
    """ 
    addtsigm(basem, pvals, lat, lon, siglevel=0.05,color='k',type='hatch')

    Add significance contour to given basemap(!)
    """
    
    # Make pval field binary
    plotfld = copy.copy(pvals) # shallow copy
    plotfld[pvals<=siglevel] = 1.5
    plotfld[pvals>siglevel] = 0

    lons, lats = np.meshgrid(lon,lat)

    if type == 'hatch':
        basem.contourf(lons,lats,plotfld,levels=[1,2],colors='none',hatches='.',latlon=True)
    else:
        basem.contour(lons,lats,plotfld,[0, 1.5],colors=color,linewidths=2,latlon=True)



def addtsig(ploth, pvals, dim1, dim2, siglevel=0.05,color='k',type='hatch',cmap='YlGnBu_r'):
    """ 
    addtsig(ploth, pvals, dim1, dim2, siglevel=0.05,color='k',type='hatch')

    Add significance contour to given pyplot handle.
    type can be 'hatch', 'color', 'cont' (the else case is cont)
    """

    print 'be careful, I am not sure this works properly 4/29/14' #@@ might just be when i screwed w/ dims and tried to plot lat x time?
    
    plotfld = copy.copy(pvals) # shallow copy
    plotfld = ma.masked_where(pvals>siglevel,pvals)
    plotfld[pvals<=siglevel] = 1.5 # significant!
    plotfld[pvals>siglevel] = 0

    # expecting lats,levs OR lats,lons OR times,lats
    dim1s, dim2s = np.meshgrid(dim1,dim2)


    if type == 'hatch':
        pc = ploth.contourf(dim1s,dim2s,plotfld,levels=[1,2],colors='none',hatches='.')
        # ploth.contourf(dim1s,dim2s,plotfld,levels=[1,2],colors='none',hatches='.')
    elif type == 'color':
        pc = ploth.pcolormesh(dim1s,dim2s,plotfld,\
                              cmap= plt.cm.get_cmap(cmap),shading='flat',\
                              vmin=0,vmax=0.05) # gouraud
    else: # 'cont'
        pc = ploth.contour(dim1s,dim2s,plotfld,levels=[0,1.5],colors='k', linewidths=2)

    return pc


#@@def vert_plot(field, lev, lat, title='',units='',cmap='blue2red_w20',cmin='',cmax='',axis=None, suppcb=0):


def map_allmonths(fld, lat, lon,title='',units='',cmap='blue2red_w20',type='sq',
           cmin='',cmax='',axis=None, suppcb=0,lmask=0,climo=0,flipmask=0,
                  pvals = None,sigtype='hatch',conts=None,latlim=None):

    months = con.get_mon()

    midx=0
    fig, spax = plt.subplots(2,6)
    fig.set_size_inches(12,4.5) # changed from 12,6 on 4/24/14. see how it saves
    if latlim!=None:
        fig.subplots_adjust(hspace=0.05,wspace=0.05)
    else:
        fig.subplots_adjust(hspace=0,wspace=0)

    for ax in spax.flat:
        if climo == 0:
            monfld = fld[midx::12,...]
            plotfld = np.mean(monfld,axis=0)
        else:
            plotfld = fld[midx,...]

        #tstat,pval = sp.stats.ttest_ind(monfldp,monfldc,axis=0)
        #sigs[midx,:,:] = ma.masked_where(pval>0.05,pval) 

        
        bm,pc = kemmap(plotfld,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,type=type,\
                     axis=ax,suppcb=1,lmask=lmask,flipmask=flipmask,units=units,latlim=latlim)
        ax.set_title(months[midx])
        if pvals != None:
            cplt.addtsigm(bm,pvals,lat,lon,type=sigtype)
        if conts != None:
            # add specified contour(s)
            lons, lats = np.meshgrid(lon,lat)
            bm.contour(lons,lats,plotfld,[conts],colors='k',linewidths=2,latlon=True)
                

        midx = midx+1

    cbar_ax = fig.add_axes([.91,.25, .02,.5])
    fig.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    plt.suptitle(title)
    
    return fig


def map_allseas(fld, lat, lon,title='',units='',cmap='blue2red_w20',type='sq',
           cmin='',cmax='',axis=None, suppcb=0,lmask=0,climo=0,flipmask=0,
                  pvals = None,sigtype='hatch',conts=None,latlim=None):

    """ default seasons = 'SON','DJF','MAM','JJA'
    """
    import cccmautils as cutl
    
    seasons = 'SON','DJF','MAM','JJA'

    midx=0
    fig,axs = plt.subplots(1,4) 
    fig.set_size_inches(12,3)
    fig.subplots_adjust(hspace=.15,wspace=.05)

    for ax in axs.flat:
        
        #print seasons[midx]
        tmpfld=np.squeeze(cutl.seasonalize_monthlyts(fld,season=seasons[midx],climo=climo))
        #print tmpfld.shape # already 2D... only when climo=1?
        #if seasons[midx]=='DJF':
        #    print tmpfld

        if climo:
            plotfld = tmpfld
        else:
            plotfld = np.mean(tmpfld,axis=0)

        #print plotfld.shape
        
        bm,pc = kemmap(plotfld,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,type=type,\
                       axis=ax,suppcb=1,lmask=lmask,flipmask=flipmask,units=units,latlim=latlim)
        ax.set_title(seasons[midx])
        if pvals != None:
            cplt.addtsigm(bm,pvals,lat,lon,type=sigtype)

        if conts != None:
            # add specified contour(s)
            lons, lats = np.meshgrid(lon,lat)
            bm.contour(lons,lats,plotfld,[conts],colors='k',linewidths=2,latlon=True)

        midx = midx+1

    cbar_ax = fig.add_axes([.91,.25, .02,.5])
    fig.colorbar(pc,cax=cbar_ax)
    plt.suptitle(title)

    return fig
        
