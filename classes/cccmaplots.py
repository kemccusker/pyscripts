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
from matplotlib.patches import Polygon
import matplotlib.font_manager as fm

"""
kemmap(fld, lat, lon, title='', units='', cmap='blue2red_w20', ptype='sq', cmin='', cmax='',
       axis=None, suppcb=0,lmask=0,flipmask=0,latlim=None)

    Inputs: fld: 2D matrix of data [lat x lon]
            lat: 1D lat array
            lon: 1D lon array
            title: plot title
            units: data units (colorbar label)
            cmap: map colormap. default blue to red, continuous
            ptype: ptype of plot/map. default 'sq' for Robinson projection. 
                  else: 'nh', 'sh' orthographic projections,
                        'eastere' is Eurasia 'stere' projection
            cmin, cmax: color scale limits
            axis: handle to the axes instance to give to basemap. default None.
            suppcb: suppress colorbar. default = 0 (do not suppress)
            lmask: add land mask (e.g. for ocean-only data). default is no mask (0)
            flipmask: a hack to flip the lat direction of lmask (e.g. for HURRELL data)
            latlim: if polar projection, it will be the equatorward limit (so far only
                    works if do a stereographic proj, not ortho 5/13/2014
            drawgrid: if True, draw parallels and meridians
            round: default True. Used if ptype 'nh' or 'sh' and latlim provided
                   otherwise the zoomed in figure will be square
            coastres: resolution of coastline dataset to plot. default 'c'=coarsest. 'l' is low res
            coastwidth: linewidth of coastlines
            area_thresh: threshold of rivers lakes to be drawn in km^2. default is for 'c' resolution
            lcol: color of land contours, default 0.7 (lightish gray)

    Returns: basemap handle (to add to bm after function call),
             pcolormesh handle (typically to add colorbar after func call)

"""

def kemmap(fld, lat, lon,title='',units='',cmap='blue2red_w20',ptype='sq',
           cmin='',cmax='',axis=None, suppcb=0,lmask=0,flipmask=0,latlim=None,drawgrid=False,
           round=True,lcol='0.7',coastres='c',coastwidth=1,area_thresh=10000,panellab=None):
    """ returns bm,pc (Basemap,Pcolor handle)
    """

    
    if cmap =='' or cmap==None:
        cmap='blue2red_w20'
        
    incmap = plt.cm.get_cmap(cmap)

    
    # default Basemap dictionary
    if ptype == 'sq':
        mapparams = dict(projection='robin',lon_0=180,lat_0=0, 
                         resolution=coastres,area_thresh=area_thresh)
    elif ptype == 'nh' or ptype=='nheur':
        if ptype=='nheur':
            lon0 = 90.
        else:
            lon0=0.
        if latlim != None: # try 'round=True' !@@@
            mapparams = dict(projection='npstere',boundinglat=latlim,lon_0=lon0,
                             resolution=coastres,area_thresh=area_thresh)
            if round==True:
                mapparams['round'] = True
            
        else:
            # try mill, hammer, merc
            mapparams = dict(projection='ortho',lon_0=lon0,lat_0=89.5,\
                             resolution=coastres,area_thresh=area_thresh) #llcrnrlon='-180',llcrnrlat='45',urcrnrlon='180',urcrnrlat='90'
        # I thought the above corner limits would work to zoom on NH but I'm getting
        # AttributeError: 'Basemap' object has no attribute '_height'
        # 5/12/14 -- don't know why. same goes for lat_0=0.        
        # 10/13/2015: changed lat_0 to 89.5 from 90. b/c .eps figure file wouldn't save
        
    elif ptype == 'sh':
        if latlim != None: # try 'round=True' !@@@
            mapparams = dict(projection='spstere',boundinglat=latlim,lon_0=0,
                             resolution=coastres,area_thresh=area_thresh)
            if round==True:
                mapparams['round'] = True
        else:
            mapparams = dict(projection='ortho',lon_0=0.,lat_0=-89.5, 
                             resolution=coastres,area_thresh=area_thresh)
            # same error if add: llcrnrlon='-180',llcrnrlat='-90',urcrnrlon='180',urcrnrlat='-45'
    elif ptype == 'eastere': # Eurasia stere projection
        #mapparams = dict(width=2500000,height=2700000,resolution='i',projection='laea',\
        #    lat_ts=62.5,lat_0=62.5,lon_0=77.0)
        mapparams = dict(llcrnrlon=40.,llcrnrlat=10.,urcrnrlon=160.,urcrnrlat=50.,
                         resolution=coastres,area_thresh=area_thresh,projection='stere',lat_0=45.,lon_0=80.)
        #mapparams = dict(width=3000000,height=3000000,resolution='c',projection='laea',\
        #                 lat_0=55.,lon_0=80.)# can't get width/height big enough -- errors
    elif ptype == 'eabksstere': # Eurasia + Barents Kara attempt
        mapparams = dict(llcrnrlon=40.,llcrnrlat=10.,urcrnrlon=175.,urcrnrlat=60.,
                         resolution=coastres,area_thresh=area_thresh,projection='stere',lat_0=45.,lon_0=80.)
    elif ptype == 'ealamb': # Lambert azimuthal equal-area
        mapparams = dict(llcrnrlon=40.,llcrnrlat=10.,urcrnrlon=160.,urcrnrlat=50.,
                         resolution=coastres,area_thresh=area_thresh,projection='laea',lat_0=45.,lon_0=80.)
    elif ptype == 'eabkslamb': # Lambert azimuthal equal-area
        mapparams = dict(llcrnrlon=40.,llcrnrlat=10.,urcrnrlon=175.,urcrnrlat=60.,
                         resolution=coastres,area_thresh=area_thresh,projection='laea',lat_0=45.,lon_0=80.)
    elif ptype == 'nastere': # North America stere projection
        mapparams = dict(llcrnrlon=220.,llcrnrlat=20.,urcrnrlon=320.,urcrnrlat=50.,
                         resolution=coastres,area_thresh=area_thresh,projection='stere',lat_0=45.,lon_0=230.)
    else:
        print "Incorrect map ptype. Choose sq,nh,sh,nastere,eastere,eabksstere,ealamb,eabkslamb"
        return -1
        
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
        pcparams = dict(latlon=True,cmap=incmap,vmin=cmin,vmax=cmax)
        if ptype not in ('eastere','eabksstere','nastere','ealamb','eabkslamb'):
            pcparams['levels']=conts
            pcparams['extend']='both'

    if axis != None: # if an axis is given, add to dict for basemap
        mapparams['ax'] = axis
    
    if ptype in ('eastere','eabksstere','nastere','ealamb','eabkslamb'):
        pcparams['shading']='flat' #'gouraud' # @@ gouraud pdf files fail/crash. flat and interp are same??

    """ m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution='l',\
    llcrnrx=0.,llcrnry=0.,urcrnrx=m1.urcrnrx/2.,urcrnry=m1.urcrnry/2.)
    """

    bm = Basemap(**mapparams)
        
    #if np.mod(lon.shape[0],2) == 0:
    if np.mod(len(lon),2) == 0:
        # add cyclic lon
        fld,lon = mpltk.basemap.addcyclic(fld,lon)
     
    if lmask==1: # the landmask has an extra lon already
        # add land mask
        lsmask=con.get_t63landmask()
        if flipmask:
            lsmask=np.flipud(lsmask)
            
        #@@@fld = ma.masked_where(lsmask!=0,fld) # 0 is ocean
        
    lons, lats = np.meshgrid(lon,lat)
        

#    pc = bm.pcolormesh(lons,lats,fld,**pcparams)
    if cmin=='':
        if ptype in ('eastere','eabksstere','nastere','ealamb','eabkslamb'):
            # for some reason these projections don't plot colors correctly so have to use pcolormesh()@@@
            pc=bm.pcolormesh(lons,lats,fld,**pcparams) 
        else:
            pc = bm.contourf(lons,lats,fld,20,**pcparams) # 20 auto levels
    else:
        if ptype in ('eastere','eabksstere','nastere','ealamb','eabkslamb'):
            pc=bm.pcolormesh(lons,lats,fld,**pcparams)
        else:
            #print '@@@@@@@@@@@ ' + str(pcparams)
            pc = bm.contourf(lons,lats,fld,**pcparams)

    if drawgrid:
        bm.drawparallels(np.arange(-90.,90.,20.),labels=[1,0,0,0])
        bm.drawmeridians(np.arange(0.,360.,30.),labels=[0,0,0,1])


    bm.drawcoastlines(color=lcol,linewidth=coastwidth)
    #bm.drawmapboundary(fill_color='#99ffff')

    if lmask:
        #@@@bm.drawmapboundary(fill_color='0.7')
        if float(lcol)>0.7:
            bm.drawcoastlines(color='.5',linewidth=coastwidth)
        else:
            bm.drawcoastlines(color='.3',linewidth=coastwidth)
        bm.fillcontinents(color=lcol)

    # I think drawlsmask puts the mask on the bottom
    #bm.drawlsmask(land_color='0.7',lsmask=con.get_t63landmask(),ax=axis)
    
    if axis!=None:
        axis.set_title(title,fontsize=10)
        if panellab!=None:
            axis.annotate(panellab,xy=(0.01,1.02),
                          xycoords='axes fraction',fontsize=18,fontweight='bold')
    else:
        plt.title(title,fontsize=10)
        if panellab!=None:
            plt.annotate(panellab,xy=(0.01,1.02),
                         xycoords='axes fraction',fontsize=18,fontweight='bold')

    # add colorbar.
    if suppcb == 0:
        cbar = bm.colorbar(pc,location='bottom',pad="5%")
        cbar.set_label(units)

    return bm,pc


def add_contours(basem, fld, lat, lon, levels=None, colors='0.5',linewidths=1):

    if np.mod(len(lon),2) == 0:
        # add cyclic lon
        fld,lon = mpltk.basemap.addcyclic(fld,lon)

    lons, lats = np.meshgrid(lon,lat)
    basem.contour(lons,lats,fld,levels=levels,
           colors=colors,linewidths=linewidths,latlon=True)


def addtsigm(basem, pvals, lat, lon, siglevel=0.05,color='k',sigtype='hatch'):
    """ 
    addtsigm(basem, pvals, lat, lon, siglevel=0.05,color='k',sigtype='hatch')

    Add significance contour to given basemap(!)
    If sigtype is anything other than 'hatch', it will be a contour
    """
    
    if np.mod(lon.shape,2) == 0:
        # add cyclic lon
        pvals,lon = mpltk.basemap.addcyclic(pvals,lon)

    # Make pval field binary
    plotfld = copy.copy(pvals) # shallow copy
    #plotfld=ma.masked_where(pvals,plotfld) # @@ didn't work...
    
    plotfld[pvals<=siglevel] = 1.5
    plotfld[pvals>siglevel] = 0

    lons, lats = np.meshgrid(lon,lat)

    if sigtype == 'hatch':
        basem.contourf(lons,lats,plotfld,levels=[1,2],colors='none',latlon=True,hatches='o')#hatches='.')
    else:
        basem.contour(lons,lats,plotfld,[0, 1.5],colors=color,linewidths=2,latlon=True)



def addtsig(ploth, pvals, dim1, dim2, siglevel=0.05,color='k',sigtype='hatch',cmap='YlGnBu_r'):
    """ 
    addtsig(ploth, pvals, dim1, dim2, siglevel=0.05,color='k',sigtype='hatch')

    Add significance contour to given pyplot handle.
    sigtype can be 'hatch', 'color', 'cont' (the else case is cont)
    """

    #print 'be careful, I am not sure this works properly 4/29/14' #@@ might just be when i screwed w/ dims and tried to plot lat x time?
    # I think it's only a potential problem for the stats through time??
    plotfld = copy.copy(pvals) # shallow copy
    plotfld = ma.masked_where(pvals>siglevel,pvals)
    plotfld[pvals<=siglevel] = 1.5 # significant!
    plotfld[pvals>siglevel] = 0

    # expecting lats,levs OR lats,lons OR times,lats
    dim1s, dim2s = np.meshgrid(dim1,dim2)


    if sigtype == 'hatch':
        pc = ploth.contourf(dim1s,dim2s,plotfld,levels=[1,2],colors='none',hatches='o')#'.')
        # ploth.contourf(dim1s,dim2s,plotfld,levels=[1,2],colors='none',hatches='.')
    elif sigtype == 'color':
        pc = ploth.pcolormesh(dim1s,dim2s,plotfld,\
                              cmap= plt.cm.get_cmap(cmap),shading='flat',\
                              vmin=0,vmax=0.05) # gouraud
    else: # 'cont'
        pc = ploth.contour(dim1s,dim2s,plotfld,levels=[0,1.5],colors='k', linewidths=2)

    return pc


def vert_plot(fld, lev, lat, title='',units='',cmap='blue2red_w20',cmin='',cmax='', ptype=None,
              axis=None, suppcb=False, latlim=None, levlim=None,addcontlines=False,screen=False,
              suppylab=False,suppxlab=False,hPa=False):
    """ screen=False: this flag tells whether the plot should be after Screen et al. 2013, ClimDyn
                          1000-300hPa, 20-90N. Ignores latlim/levlim
        suppylab=False: suppress the y labels
    """

    if hPa:
        lats,levs = np.meshgrid(lat,lev) # do not divide by 100, already in hPa
    else:
        lats,levs = np.meshgrid(lat,lev/100.)

    if cmap =='' or cmap==None:
        cmap='blue2red_w20'
        
    incmap = plt.cm.get_cmap(cmap)
    
    cmlen=float(incmap.N) # or: from __future__ import division

    if cmlen>200:
        cmlen = float(15) # ie. if using a built in colormap that is continuous, want smaller # of colors

    if cmin =='':
        # parameters for plot
        pparams = dict(cmap=incmap)
        
    else:
        
        incr = (cmax-cmin) /cmlen
        conts = np.arange(cmin,cmax+incr,incr)
        pparams = dict(cmap=incmap,levels=conts,vmin=cmin,vmax=cmax,extend='both')

    if axis != None:
        # what to do w/ axis?
        ax=axis
    else:
        ax=plt.gca()

    if cmin=='':
        cf = ax.contourf(lats,levs,fld,**pparams) #cmlen autolevels@@removed
    else:
        cf = ax.contourf(lats,levs,fld,**pparams)

    if addcontlines:
        if cmin=='':
            ax.contour(lats,levs,fld,colors='.3') #cmlen autolevels@@removed
        else:
            ax.contour(lats,levs,fld,levels=conts,colors='.3') # may have to make a dict for levels key

    #ax.set_ylim(levlim,1000)
    #ax.invert_yaxis()

    if screen==True:
        ax.set_yticks([900, 700, 500, 300])
        if suppylab==False:
            ax.set_yticklabels((900, 700, 500, 300))
        else:
            ax.set_yticklabels('')
        ax.set_xlim(20,90)
        ax.set_xticks([40, 60, 80])
        if suppxlab==False:
            ax.set_xticklabels((40, 60, 80))
        else:
            ax.set_xticklabels('')
    else:
        ax.set_yscale('log')
        ax.set_yticks([1000,800, 500, 300, 100, 10])
        if ptype=='nh':
            ax.set_xlim(latlim,90)
            ax.set_xticks([20, 45, 70])
            if suppxlab:
                ax.set_xticklabels((20, 45, 70))
            else:
                ax.set_xticklabels('')
        elif ptype=='sh':
            ax.set_xlim(-90,latlim)
            ax.set_xticks([-70,-45,-20])
            if suppxlab==False:
                ax.set_xticklabels((-70,-45,-20))
            else:
                ax.set_xticklabels('')
        else:
            ax.set_xlim(-90,90)
            ax.set_xticks([-45, 0, 45])
            if suppxlab==False:
                ax.set_xticklabels((-45, 0, 45))
            else:
                ax.set_xticklabels('')
            
        if suppylab==False:
            ax.set_yticklabels((1000,800,500,300,100,10))
        else:
            ax.set_yticklabels('')

            

    ax.set_ylim(levlim,1000)
    ax.invert_yaxis()
    if suppylab==False:
        ax.set_ylabel('Pressure (hPa)')
    if suppxlab==False:
        ax.set_xlabel('Latitude')
    ax.set_title(title)

    if suppcb==False:
        #cbar = ax.colorbar(cf)
        
        #cbar_ax = plt.add_axes([.91,.15, .02,.7])
        plt.colorbar(cf) #pc,cax=cbar_ax)

    return cf

def map_allmonths(fld, lat, lon,title='',units='',cmap='blue2red_w20',ptype='sq',
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

        
        bm,pc = kemmap(plotfld,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,ptype=ptype,\
                     axis=ax,suppcb=1,lmask=lmask,flipmask=flipmask,units=units,latlim=latlim)
        ax.set_title(months[midx])
        if pvals != None:
            addtsigm(bm,pvals[midx,...],lat,lon,sigtype=sigtype)
        if conts != None:
            # add specified contour(s)
            lons, lats = np.meshgrid(lon,lat)
            bm.contour(lons,lats,plotfld,[conts],colors='k',linewidths=2,latlon=True)
                

        midx = midx+1

    cbar_ax = fig.add_axes([.91,.25, .02,.5])
    fig.colorbar(pc,cax=cbar_ax) # or do bm.colorbar....
    plt.suptitle(title)
    
    return fig

def add_colorbar(fig,phand,orientation='vertical',pos=None,label=None,fontsize=14):
    """
        returns cbar_ax, cbar

    """

    try:
        if pos==None:
            if orientation=='vertical':
                pos = [.91,.25, .02,.5]
            elif orientation=='horizontal':
                pos = [.25,.05, .5,.03]
            else:
                print 'orientation not recognized!'
                raise Exception
    except:
        raise

    cbar_ax = fig.add_axes(pos)
    cparams={'orientation': orientation}

    if label!=None:
        cparams['label']=label
    
        if orientation=='horizontal':
            lab = cbar_ax.xaxis.label
        else:
            lab = cbar_ax.yaxis.label

        font = fm.FontProperties(size=fontsize)
        lab.set_font_properties(font)

    cparams['cax'] = cbar_ax
    cbar = fig.colorbar(phand,**cparams) #cax=cbar_ax,orientation=orientation,label=label)

    return cbar_ax, cbar

def map_allseas(fld, lat, lon,title='',units='',cmap='blue2red_w20',ptype='sq',
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
        
        bm,pc = kemmap(plotfld,lat,lon,cmin=cmin,cmax=cmax,cmap=cmap,ptype=ptype,\
                       axis=ax,suppcb=1,lmask=lmask,flipmask=flipmask,units=units,latlim=latlim)
        ax.set_title(seasons[midx])
        if pvals != None:
            addtsigm(bm,pvals,lat,lon,sigtype=sigtype)

        if conts != None:
            # add specified contour(s)
            lons, lats = np.meshgrid(lon,lat)
            bm.contour(lons,lats,plotfld,[conts],colors='k',linewidths=2,latlon=True)

        midx = midx+1

    cbar_ax = fig.add_axes([.91,.25, .02,.5])
    fig.colorbar(pc,cax=cbar_ax)
    plt.suptitle(title)

    return fig
        

def plotvert_allseas(fld, lev, lat,title='',units='',cmap='blue2red_w20',ptype=None,
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
            vpparams = dict(cmin=cmin,cmax=cmax,cmap=cmap,ptype=ptype,
                            axis=ax,suppcb=True,units=units,
                            latlim=latlim,levlim=levlim,
                            addcontlines=addcontlines,screen=screen)
        else:
            vpparams = dict(cmin=cmin,cmax=cmax,cmap=cmap,ptype=ptype,
                            axis=ax,suppcb=True,units=units,
                            latlim=latlim,levlim=levlim,
                            addcontlines=addcontlines,screen=screen,
                            suppylab=True)
            
        cf = vert_plot(plotfld,lev,lat,**vpparams)
        
        ax.set_title(seasons[midx])
        if pvals != None:
            #addtsig(cf,pvals,lat,lon,sigtype=sigtype) # @@@ untested? no lat... 5/10/15
            addtsig(cf,pvals,lat,lev,sigtype=sigtype) 

        midx = midx+1

    cbar_ax = fig.add_axes([.91,.25, .02,.5])
    fig.colorbar(cf,cax=cbar_ax)
    plt.suptitle(title)

    return fig
        

def plot_region(regname,ptype='nh',axis=None,latlim=None,limsdict=None): 
    """ plot_region(regname,ptype='nh',axis=None,latlim=None,limsdict=None):
                     Given a region name, plot it for reference.

                     latlims is unused right now
                     if passing limsdict (to override regname), set regname='other'
    """

    if regname=='other':
        reglims=limsdict
    else:
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
    kemmap(dummym,lat,lon,ptype=ptype,axis=axis,latlim=latlim,suppcb=1,
           cmin=-11,cmax=2,cmap='blue2blue_w10',drawgrid=True)

    

def plot_allregions(ptype='nh'):
    """ plot_allregions(ptype='nh'): plot all defined regions
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
    fig.set_size_inches(cols*2,rows*3)
    
    for aii,ax in enumerate(spax.flat):

        dummy = con.get_t63landmask() # dummy data
        
        regkey = regdict.keys()[aii]
        limsdict = regdict[regkey]
        # mask the dummy data
        dummym,dmask = cutl.mask_region(dummy,lat,lon,regkey,limsdict=limsdict)
                
        kemmap(dummym,lat,lon,ptype=ptype,axis=ax,suppcb=1,cmin=-11,cmax=2,cmap='blue2blue_w10')        
        ax.set_title(regkey)
        ax.set_xlabel(str(limsdict['latlims']) + ',' +str(limsdict['lonlims']) )
        
        if aii==nreg-1: break # get out of loop if done with regions
   

def plot_regions(regions, ptype='nh', colors='k', latlim=None, drawgrid=False):
    """ plot multiple region boundaries on one map

           regions: tuple of region names

           default color = black
           colors can be a tuple of color specifications. 
                index will wrap if not enough

    """

    dummy = con.get_t63landmask()
    lat = con.get_t63lat()
    lon = con.get_t63lon()

    lons,lats = np.meshgrid(lon,lat)
    dummy = dummy*np.nan

    plt.figure()
    bm,pc = kemmap(dummy,lat,lon,ptype=ptype,latlim=latlim,suppcb=1,
                   drawgrid=drawgrid, lmask=True)
    
    for ii,reg in enumerate(regions):

        if len(colors)==1: clr=colors
        else:
            if ii==len(colors): ii=0 # wrap index
            clr=colors[ii]

        add_regionpolym(reg,bm, ec=clr)

    
def kemscatter(fldx,fldy,weights=None,axis=None,xlims=None,ylims=None,suppregress=False,
               supponetoone=False, marker='.',color='k',grid=True):
    """ fldx, fldy, and weights if given will be flattened.
          suppregress: suppress regression line if True.
          supponetoone: suppress one-to-one line if True.
          
    """

    if axis!=None:
        ax=axis
    else:
        ax=plt.gca()

    if weights!=None:
        scatxx=fldx.flatten()*weights.flatten()
        scatyy=fldy.flatten()*weights.flatten()        
    else:
        scatxx=fldx.flatten()
        scatyy=fldy.flatten()

    ax.scatter(scatxx,scatyy, marker=marker,color=color)


    if xlims==None:
        axxlims = ax.get_xlim()
    else:
        axxlims = xlims
    if ylims==None:
        axylims = ax.get_ylim()
    else:
        axylims = ylims
    
    onex=np.arange(-100,100) # a hack. Need to actuall infer it from data @@
    oney=onex
    if ~supponetoone:
        # one-to-one line
        ax.plot(onex,oney,color='b',linestyle='--')


    if ~suppregress:
        mm, bb = np.polyfit(scatxx, scatyy, 1)
        ax.plot(onex,mm*onex + bb, color=color)
        rr = np.corrcoef(scatxx, scatyy)[0,1]
        #print rr
        rsq = rr**2
        #print 'R squared: ' + str(rsq)
        val = '$%.2f$'%(rsq)
        ax.annotate('$R^2$= ' + val, xy=(axxlims[0]+.1*axxlims[1], axylims[1]-.2*axylims[1]),
                    xycoords='data') # upper left?
        mval = '$%.2f$'%(mm)
        ax.annotate('$m$= ' + mval, xy=(axxlims[0]+.1*axxlims[1], axylims[1]-.4*axylims[1]),
                    xycoords='data')

    # reset the limits
    ax.set_ylim(axylims)
    ax.set_xlim(axxlims)
    ax.grid(grid)


    return ax


def plot_pattcorrs(pcdf, pcdf2=None, rmin=None, axis=None, tftype='seasonal'):
    """ plot_pattcorrs(pcdf, axis):
                          this figure plots range of pattern correlations
                          as bars, with mean patt corr marker, plus % common
                          variation annotated.

             pcdf: DataFrame of pattern correlations
             pcdf2: second set of data to plot
             rmin: minimum r value (pattern correlation) that is significant for the data
             axis: plot axis
             tftype: time frequency type: 'seasonal' or 'monthly'
    """
    import math as math
    
    if tftype=='seasonal':
        numt=4
        seasons=('SON','DJF','MAM','JJA')
    else:
        print 'ONLY SEASONAL IMPLEMENTED @@ 11/25/14'
        numt=12 # monthly
        seasons=con.get_mon()

    # within ensemble pattern correlation combinations
    # number of combinations: n choose k (here k=2)
    #    n! / (k! *(n-k)!)
    ncomb = math.factorial(len(pcdf)) / (2 * math.factorial(len(pcdf)-2) )
    
    # first, create matrix of data to take mean, max, min, etc
    pcstack = np.zeros((ncomb,numt)) # ncomb combos x 4 seasons (or 12 mon)
    stackii=0
    keys=pcdf.keys()
    for eii,skey1 in enumerate(keys):
        for ineii,skey2 in enumerate(keys[eii+1:len(keys)]):
            # stack the pc's so I can take max/min, etc
            pcstack[stackii] = pcdf[skey1][skey2]
            stackii+=1
    mxsea = np.max(pcstack,axis=0)
    mnsea = np.min(pcstack,axis=0)
    avgsea = np.mean(pcstack,axis=0)
    tmpstack=copy.copy(pcstack)
    tmpstack[pcstack<0] = 0 # get rid of negative correlations. they are zero for all intents & purposes in this case
    pcstacksq = np.power(tmpstack,2) # square the patt corr (the ones that are positive)
    avgseasq = np.mean(pcstacksq,axis=0) # average squared patt corr (coef of determination)

    if not pcdf2.empty:
        ncomb2 = math.factorial(len(pcdf2)) / (2 * math.factorial(len(pcdf2)-2) )
        pcstack2 = np.zeros((ncomb2,numt)) # ncomb combos x 4 seasons
        stackii=0
        keys=pcdf2.keys()
        for eii,skey1 in enumerate(keys):
            for ineii,skey2 in enumerate(keys[eii+1:len(keys)]):
                # stack the pc's so I can take max/min, etc
                pcstack2[stackii] = pcdf2[skey1][skey2]           
                stackii+=1   
        mxsea2 = np.max(pcstack2,axis=0)
        mnsea2 = np.min(pcstack2,axis=0)
        avgsea2 = np.mean(pcstack2,axis=0)
        tmpstack2=copy.copy(pcstack2)
        tmpstack2[pcstack2<0] = 0 # get rid of negative correlations. they are zero for all intents & purposes in this case
        pcstacksq2 = np.power(tmpstack2,2) # square the patt corr (the ones that are positive)
        avgseasq2 = np.mean(pcstacksq2,axis=0) # average squared patt corr (coef of determination)

    # add annual mean
    annwgts = np.array((91,90,92,92))/365.
    print annwgts

    #fig,axs = plt.subplots(1,1)
    #fig.set_size_inches(6,4) # more squat

    ax = axis #[0]

    legtpl=()

    wi=0.1 # width of bar
    incr=0.3 # how much to shift in b/w the 2 sets of data
    xxsea2=np.arange(1,6)
    xboxmin=xxsea2-.2

    if rmin!=None:
        ax.axhspan(-1*rmin,rmin,color='orange',alpha=.3) # shade where corr is NOT significant

    fillcol='0.7'
    fillcol2=ccm.get_linecolor('dodgerblue')
    for boxii in range(0,4): # loop through seasons
        # shaded bars/boxes
        ax.fill_between((xboxmin[boxii],xboxmin[boxii]+wi),mnsea[boxii],mxsea[boxii],color=fillcol, alpha=.5)
        if not pcdf2.empty:
            ax.fill_between((xboxmin[boxii]+incr,xboxmin[boxii]+incr+wi),mnsea2[boxii], mxsea2[boxii],color=fillcol2, alpha=.5)

        # markers
        ax.plot(xboxmin[boxii]+wi/2.,avgsea[boxii],color='k',marker='_',linestyle='none',markersize=15)#,alpha=.7)
        if not pcdf2.empty:
            ax.plot(xboxmin[boxii]+wi/2.+incr,avgsea2[boxii],color='b',marker='_',linestyle='none',markersize=15)#,alpha=.7)   

        # mean values (text)
        val = '$%.0f$'%(avgseasq[boxii]*100) # @@ square the corrs before taking mean
        ax.annotate(val+'%', xy=(xboxmin[boxii]+wi/2. -.09, .95),  xycoords='data')
        if not pcdf2.empty:
            val = '$%.0f$'%(avgseasq2[boxii]*100) # @@ square the corrs before taking mean
            ax.annotate(val+'%', xy=(xboxmin[boxii]+wi/2.+incr -.06, .95),  xycoords='data')

    boxii=boxii+1

    # add annual mean -----------
    annwgtst = np.tile(annwgts,(len(pcdf),1)) # tile
    # each ens w/ each other
    # pcstacksq is ncomb x numt
    annwgtst = np.tile(annwgts,(ncomb,1)) # tile
    ann = np.average(pcstack,weights=annwgts,axis=1) # ann mean per patt corr
    annmax = np.max(ann)
    annmin = np.min(ann)
    avgann = np.mean(ann)
    annsq = np.average(pcstacksq,weights=annwgts,axis=1) # ann mean per patt corr
    avgannsq = np.mean(annsq)

    if not pcdf2.empty:
        # add annual mean
        annwgtst = np.tile(annwgts,(len(pcdf2),1)) # tile
        # each ens w/ each other
        # pcstacksq is ncomb x numt
        annwgtst = np.tile(annwgts,(ncomb2,1)) # tile
        ann = np.average(pcstack2,weights=annwgts,axis=1) # ann mean per patt corr
        annmax2 = np.max(ann)
        annmin2 = np.min(ann)
        avgann2 = np.mean(ann)
        annsq2 = np.average(pcstacksq2,weights=annwgts,axis=1) # ann mean per patt corr
        avgannsq2 = np.mean(annsq2)

    # plot annual mean markers
    ax.fill_between((xboxmin[boxii],xboxmin[boxii]+wi),annmin,annmax,color=fillcol,alpha=.5)
    if not pcdf2.empty:
        ax.fill_between((xboxmin[boxii]+incr,xboxmin[boxii]+incr+wi),annmin2,annmax2,color=fillcol2,alpha=.5)

    ax.plot(xboxmin[boxii]+wi/2.,avgann,color='k',marker='_',linestyle='none',markersize=15)#,alpha=.9)
    if not pcdf2.empty:
        ax.plot(xboxmin[boxii]+wi/2.+incr,avgann2,color='b',marker='_',linestyle='none',markersize=15)#,alpha=.9)

    val = '$%.0f$'%(avgannsq*100) # @@ square the corrs before taking mean
    ax.annotate(val+'%', xy=(xboxmin[boxii]+wi/2. -.09, .95),  xycoords='data')
    if not pcdf2.empty:
        val = '$%.0f$'%(avgannsq2*100) # @@ square the corrs before taking mean
        ax.annotate(val+'%', xy=(xboxmin[boxii]+wi/2.+incr -.06, .95),  xycoords='data')

    ax.set_ylabel('Pattern Correlation')
    ax.set_xlim((.5,5.5))
    ax.set_xticks(xxsea2)
    ax.set_xticklabels((seasons)+('ANN',))
    ax.set_ylim((0,1))
    ax.grid(True)

def plot_onecascade(topdata,bottdata,topy,botty,ax=None, mparams=None, lparams=None,color='r',pvalst=None,pvalsb=None,siglevel=0.05):
    """ def plot_onecascade(topdata,bottdata,topy,botty,ax=None, mparams=None, lparams=None):

           mparams: dictionary of key/values for marker
           lparams: dictionary of key/values for line
           pvalst and pvalsb should be pvalues for the data in top and bottom, respectively (for the markers)
    """
    mew=1.5 # marker edge width
    lw=1
    ms=7
    if ax==None:
        ax=plt.gca()
    if mparams ==None:
        mparams = dict(marker='o',markersize=ms,linestyle='none',color=color,mec=color,mew=mew)
    if lparams ==None:
        lparams = dict(color=color,linewidth=lw)
        
    tlen = len(topdata) # this should match bottdata length
    if tlen != len(bottdata):
        print 'topdata length is not equal to bottdata length! ' + str(tlen) + '!=' + str(len(bottdata))
        return

    for tii,td in enumerate(topdata):
        ax.plot(td,topy,fillstyle='none',**mparams) # marker
        # test for significance
        if pvalst != None:
            if pvalst[tii] <= siglevel:
                ax.plot(td,topy,**mparams) # sig marker
            
        #print 'topdata els: ' + str(td) + ', ' + str(topy)
        
        bd=bottdata[tii]
        
        #print 'bottdata els: ' + str(bd) + ', ' + str(botty)
        ploty=botty*np.ones(len(bd))
        ax.plot(bd,ploty,fillstyle='none',**mparams) # all markers
        if pvalsb!=None:
            bdpv = pvalsb[tii]
            print 'bdpv: ' + str(bdpv) # @@@@@
            print 'bd: ' + str(bd)
            pvmask=ma.masked_where(np.array(bdpv)>siglevel,np.array(bd))
            print pvmask # @@@@
            ax.plot(pvmask,ploty,**mparams) # all sig markers

        for b in bd:
            #print 'element of bottdata: ' + str(b)

            ax.plot((b,td), (botty,topy),**lparams) # individual lines

def add_regionpolym(region, bm, limsdict=None,ax=None, fc='none', ec='black', lw=2, alpha=1):
    """ add region box in map coords to a basemap.

            if region=='other', use limsdict 
    """

    if region=='other':
        rlims=limsdict
    else:
        rlims=con.get_regionlims(region)

    lats,lons = corners_to_poly(rlims)
    add_polym(lats,lons,bm,ax=ax,fc=fc,ec=ec,lw=lw,alpha=alpha)


def add_polym( lats, lons, bm, ax=None, fc='red', ec='none', lw=1.0, alpha=0.4):
    """ from:http://stackoverflow.com/questions/12251189/how-to-draw-rectangles-on-a-basemap

        Adds a polygon given lat/lon vertices. Use corners_to_poly() to convert
             limsdict (latlims/lonlims) to lats, lons.

        lats: array of lat vertices (in degrees)
        lons: array of lon vertices matching lats (in degrees)
        bm: Basemap instance
        ax: axis
        fc: facecolor for polygon patch
        ec: edgecolor of patch
        lw= linewidth for edge
        alpha: opacity of patch

    """

    x, y = bm( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, facecolor=fc, alpha=alpha, edgecolor=ec, linewidth=lw )
    if ax==None:
        ax=plt.gca()

    ax.add_patch(poly)


def corners_to_poly(limsdict, resolution=30):
    """ utility function to convert corners in a limsdict (from constants.get_regionlims())
           to vertices for a polygon, by defining each of the 4
           line segments (with num points = resolution). This is useful
           for drawing a region polygon on a basemap using 
           
           limsdict: {'latlims': [southern lat lim, northern lat lim], 
                      'lonlims': [western lon lim, eastern lon lim] }
           resolution: number of points to include in line segment
           
           returns lats, lons (stacked arrays of all lat points and all lon points)
    """
    
    lllat = limsdict['latlims'][0]; ullat = limsdict['latlims'][1] # lower/upper left lats
    lllon = limsdict['lonlims'][0]; ullon = limsdict['lonlims'][0] # lower/upper left lons
    urlat = limsdict['latlims'][1]; lrlat = limsdict['latlims'][0] # upper/lower right lats
    urlon = limsdict['lonlims'][1]; lrlon = limsdict['lonlims'][1] # upper/lower right lons

    # left side line segment
    latsL = np.linspace(lllat,ullat,resolution)
    lonsL = np.linspace(lllon,ullon,resolution)

    # top line segment
    latsT = np.linspace(ullat,urlat,resolution)
    lonsT = np.linspace(ullon,urlon,resolution)

    # right line segment
    latsR = np.linspace(urlat,lrlat,resolution)
    lonsR = np.linspace(urlon,lrlon,resolution)

    # bottom line segment
    latsB = np.linspace(lrlat,lllat, resolution )
    lonsB = np.linspace(lrlon,lllon, resolution )


    lats=np.hstack((latsL,latsT,latsR,latsB))
    lons=np.hstack((lonsL,lonsT,lonsR,lonsB))

    return lats,lons


def add_regressline(input1,input2,ax=None,axxlims=None,color='k',linewidth=1):

    
    if ax==None:
        ax=plt.gca()

    if axxlims==None:
        axxlims = ax.get_xlim()

    mm,bb,rval,pval = cutl.regress(input1,input2)
    print 'regress pval,rval: ' + str(pval),str(rval)

    onex=np.linspace(axxlims[0],axxlims[1])    
    ax.plot(onex,mm*onex + bb, color=color,linewidth=linewidth)
