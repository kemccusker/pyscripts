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
import cccmautils as cutl

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
    If type is anything other than 'hatch', it will be a contour
    """
    
    # Make pval field binary
    plotfld = copy.copy(pvals) # shallow copy
    #plotfld=ma.masked_where(pvals,plotfld) # @@ didn't work...
    
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

    #print 'be careful, I am not sure this works properly 4/29/14' #@@ might just be when i screwed w/ dims and tried to plot lat x time?
    # I think it's only a potential problem for the stats through time??
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


def vert_plot(fld, lev, lat, title='',units='',cmap='blue2red_w20',cmin='',cmax='', type=None,
              axis=None, suppcb=False, latlim=None, levlim=None,addcontlines=False,screen=False,
              suppylab=False):
    """ screen=False: this flag tells whether the plot should be after Screen et al. 2013, ClimDyn
                          1000-300hPa, 20-90N. Ignores latlim/levlim
        suppylab=False: suppress the y labels
    """

    lats,levs = np.meshgrid(lat,lev/100.)

    if cmap =='' or cmap==None:
        cmap='blue2red_w20'
        
    incmap = plt.cm.get_cmap(cmap)

    if cmin =='':
        # parameters for plot
        pparams = dict(cmap=incmap)
    else:
        cmlen=float(incmap.N) # or: from __future__ import division

        if cmlen>200:
            cmlen = float(15) # ie. if using a built in colormap that is continuous, want smaller # of colors
        
        incr = (cmax-cmin) /cmlen
        conts = np.arange(cmin,cmax+incr,incr)
        pparams = dict(cmap=incmap,levels=conts,vmin=cmin,vmax=cmax,extend='both')

    if axis != None:
        # what to do w/ axis?
        ax=axis
    else:
        ax=plt.gca()
    
    cf = ax.contourf(lats,levs,fld,**pparams)

    if addcontlines:
        ax.contour(lats,levs,fld,levels=conts,colors='.3') # may have to make a dict for levels key

    ax.set_ylim(levlim,1000)
    ax.invert_yaxis()

    if screen==True:
        ax.set_yticks([900, 700, 500, 300])
        if suppylab==False:
            ax.set_yticklabels((900, 700, 500, 300))
        else:
            ax.set_yticklabels('')
        ax.set_xlim(20,90)
        ax.set_xticks([40, 60, 80])
        ax.set_xticklabels((40, 60, 80))
    else:
        ax.set_yscale('log')
        ax.set_yticks([1000,800, 500, 300, 100, 10])
        if type=='nh':
            ax.set_xlim(latlim,90)
            ax.set_xticks([20, 45, 70])
            ax.set_xticklabels((20, 45, 70))
        elif type=='sh':
            ax.set_xlim(-90,latlim)
            ax.set_xticks([-70,-45,-20])
            ax.set_xticklabels((-70,-45,-20))
        else:
            ax.set_xlim(-90,90)
            ax.set_xticks([-45, 0, 45])
            ax.set_xticklabels((-45, 0, 45))
            
        if suppylab==False:
            ax.set_yticklabels((1000,800,500,300,100,10))
        else:
            ax.set_yticklabels('')

    if suppylab==False:
        ax.set_ylabel('Pressure (hPa)')
    ax.set_xlabel('Latitude')
    ax.set_title(title)

    if suppcb==False:
        cbar = ax.colorbar(cf)
        
        #cbar_ax = fig4.add_axes([.91,.15, .02,.7])
        #fig4.colorbar(pc,cax=cbar_ax)

    return cf

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
            addtsigm(bm,pvals[midx,...],lat,lon,type=sigtype)
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
            addtsigm(bm,pvals,lat,lon,type=sigtype)

        if conts != None:
            # add specified contour(s)
            lons, lats = np.meshgrid(lon,lat)
            bm.contour(lons,lats,plotfld,[conts],colors='k',linewidths=2,latlon=True)

        midx = midx+1

    cbar_ax = fig.add_axes([.91,.25, .02,.5])
    fig.colorbar(pc,cax=cbar_ax)
    plt.suptitle(title)

    return fig
        

def plotvert_allseas(fld, lev, lat,title='',units='',cmap='blue2red_w20',type=None,
                      cmin='',cmax='',axis=None, suppcb=0,climo=0,
                      pvals = None,sigtype='hatch',latlim=None,levlim=None,
                      addcontlines=True,screen=False):

    """ default seasons = 'SON','DJF','MAM','JJA'
    """
    import cccmautils as cutl
    
    seasons = 'SON','DJF','MAM','JJA'

    midx=0
    fig,axs = plt.subplots(1,4) 
    fig.set_size_inches(12,4)
    fig.subplots_adjust(hspace=.15,wspace=.05)

    for axii,ax in enumerate(axs.flat):
        
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

        if axii==0:
            vpparams = dict(cmin=cmin,cmax=cmax,cmap=cmap,type=type,
                            axis=ax,suppcb=True,units=units,
                            latlim=latlim,levlim=levlim,
                            addcontlines=addcontlines,screen=screen)
        else:
            vpparams = dict(cmin=cmin,cmax=cmax,cmap=cmap,type=type,
                            axis=ax,suppcb=True,units=units,
                            latlim=latlim,levlim=levlim,
                            addcontlines=addcontlines,screen=screen,
                            suppylab=True)
            
        cf = vert_plot(plotfld,lev,lat,**vpparams)
        
        ax.set_title(seasons[midx])
        if pvals != None:
            addtsig(cf,pvals,lat,lon,type=sigtype)

        midx = midx+1

    cbar_ax = fig.add_axes([.91,.25, .02,.5])
    fig.colorbar(cf,cax=cbar_ax)
    plt.suptitle(title)

    return fig
        

def plot_region(regname,type='nh',axis=None,latlim=None): 
    """ plot_region(regname,type='nh',axis=None,latlim=None):
                     Given a region name, plot it for reference.
    """

    reglims = con.get_regionlims(regname)
    latlims = reglims['latlims']
    lonlims = reglims['lonlims']

    dummy = con.get_t63landmask()
    lat = con.get_t63lat()
    lon = con.get_t63lon()

    lons,lats = np.meshgrid(lon,lat)

    reglatsbool = np.logical_and(lat>latlims[0],lat<latlims[1])
    reglonsbool = np.logical_and(lon>lonlims[0],lon<lonlims[1])

    # mask everything but the region of interest
    regmask = np.logical_or( 
                            np.logical_or(lats<latlims[0],lats>latlims[1]), 
                            np.logical_or(lons<lonlims[0],lons>lonlims[1]))
    dummym = ma.masked_where(regmask,dummy)

    plt.figure()
    kemmap(dummym,lat,lon,type=type,axis=axis,latlim=latlim,suppcb=1,cmin=-11,cmax=2,cmap='blue2blue_w10')

    

def plot_allregions(type='nh'):
    """ plot_allregions(type='nh'): plot all defined regions
    """

    regdict = con.get_regiondict()
    nreg = len(regdict)
    rem = np.mod(nreg,2)
    rows = 2
    
    if rem==0:
        cols = nreg/rows
    else:
        cols = nreg/rows + rem # rem will always be 1 if dividing by 2

    if cols>8:
        rows = 3
        cols = nreg/rows + np.mod(nreg,rows)/2

    print 'nrows: ' + str(rows) + ' ncols: ' + str(cols) # @@@
    
    lat = con.get_t63lat()
    lon = con.get_t63lon()
    
    fig,spax = plt.subplots(rows,cols)
    
    for aii,ax in enumerate(spax.flat):

        dummy = con.get_t63landmask() # dummy data
        
        regkey = regdict.keys()[aii]
        limsdict = regdict[regkey]
        # mask the dummy data
        dummym,dmask = cutl.mask_region(dummy,lat,lon,limsdict)
                
        kemmap(dummym,lat,lon,type=type,axis=ax,suppcb=1,cmin=-11,cmax=2,cmap='blue2blue_w10')        
        ax.set_title(regkey)
        
        if aii==nreg-1: break # get out of loop if done with regions
        
    
