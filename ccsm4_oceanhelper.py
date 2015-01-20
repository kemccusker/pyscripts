# ccsm4 ocean helper functions

import numpy as np
import numpy.ma as ma


def mask_seafloor(input, kmt, axis=0):
    """ def mask_seafloor(input, kmt):
              default is that input has depth (zt) as first dimension.
              axis is axis of depth dimension.
              Change axis=1 if time is first and depth is second.

              returns input as masked array 
    """

    tmp = ma.zeros(input.shape)
    ndepth = input.shape[axis]

    if axis==1:
        # need to tile kmt
        kmtt = np.tile(kmt,(input.shape[0],1,1))

    for lii in np.arange(0,ndepth): # loop through each depth layer

        # mask out levels below sea floor
        if axis==0:
            tmp[lii,...] = ma.masked_where(kmt <= lii, input[lii,...])
        elif axis==1:
            tmp[:,lii,...] = ma.masked_where(kmtt <= lii, input[:,lii,...])
        else:
            print 'axis must = 1 or 0'

    return tmp


def ocn_zonal_mean(input, kmt, axis=2):
    """ def ocn_zonal_mean(input, kmt, axis=2):
            
            Take zonal mean of ocean data. To do this, must know sea floor depth.
            Default is that longitude is axis=2 (ie variable depth x lat x lon). 
            Change to axis=3 if time dimension is first.

            returns zonal mean of size dim(input) - 1
    """

    if axis==3:
        # there is a time dimension, so depth is axis 1
        tmp = mask_seafloor(input, kmt, axis=1)
    else:
        # no time dimension, depth is axis 0
        tmp = mask_seafloor(input, kmt, axis=0)

    zm = ma.mean(tmp,axis=axis)

    return zm


def ocn_regional_average(input, tlat,tlon, tarea, limsdict):
    """ def ocn_regional_average(input, tlat,tlon, tarea, limsdict):

            horizontal regional average, given lat and lon limits.
            For now this is only good in the southern hem because of 
                the more regular grid.

            limsdict = ('latlims': (slim,nlim), 'lonlims': (wlim,elim))
            @@@ not done yet 1/17/2015
    """
    
    print 'ocn_regional_average() not implemented yet @@@@@'
    latlims=limsdict['latlims']
    slim = latlims[0]; nlim = latlims[1]

    onelat=tlat[:,0]


    # pick sub-region in tarea
    # @@@@ I am not sure what to do about the weird grid in the northern hem. prob
    # @@@@ need to remap first.
    tareasub = tarea[np.logical_and(onelat<=nlim,onelat>slim),:] # sublat x lon
    
    """
    Nlim=-72   # how close to "ice shelf"?
    Slim=-74
    # control: ==============================
    #   pick sublat range from area
    tareapigcsub = tempareapigc[:,np.logical_and(onelon<= Nlim,onelon>Slim),:] # depth x sublat x lon
    #   total area in the sublat x lon box 
    totareapigxyc = ma.sum(ma.sum(tareapigcsub,axis=2),axis=1) # depth
    #   total area in the zonal dir
    totareapigxc = ma.sum(tareapigcsub,axis=2) # depth x sublat
    #   tile sublat/lon box tot area with sublat range
    totareapigxytc = np.tile(totareapigxyc,(tareapigcsub.shape[1],1))
    totareapigxytc = np.transpose(totareapigxytc,(1,0))
    #   Now weights for averaging a zonal mean, with latitude
    wgtspigc = totareapigxc / totareapigxytc
    #   Average meridionally
    avgpigc,ret1 = ma.average(tempfldcpigzm[:,np.logical_and(onelon<= Nlim,onelon>Slim)],
                     axis=1,weights=wgtspigc,returned=True)

    # === and with time dim
    # Pert runs, with time =================
    # Now pick the sublat range from tareapig and create weights
    tareapigsub1 = tempareapig1[:,:,np.logical_and(onelon<= Nlim,onelon>Slim),:] # time x depth x sublat x lon
    tareapigsub2 = tempareapig2[:,:,np.logical_and(onelon<= Nlim,onelon>Slim),:] # time x depth x sublat x lon
    ##  total area in sublat x lon box
    totareapigxy1 = ma.sum(ma.sum(tareapigsub1,axis=2),axis=2) # time x depth 
    totareapigxy2 = ma.sum(ma.sum(tareapigsub2,axis=2),axis=2) # time x depth

    totareapigx1 = ma.sum(tareapigsub1,axis=3) # total area in the lon range for e/ time x depth x lat
    totareapigx2 = ma.sum(tareapigsub2,axis=3) # total area in the lon range for e/ time x depth x lat

    # should be able to take zonal mean no prob
    #   then average in meridional direction with weights totareapigx / totareapigxy per latitude
    # zm already done: tempfldpigzm1,2

    # first tile over sublat to create weights
    totareapigxyt2 = np.tile(totareapigxy2,(tareapigsub2.shape[2],1,1)) 
    totareapigxyt2=np.transpose(totareapigxyt2,(1,2,0)) # time x depth x lat?
    wgtspig2 = totareapigx2 / totareapigxyt2

    totareapigxyt1 = np.tile(totareapigxy1,(tareapigsub1.shape[2],1,1)) 
    totareapigxyt1=np.transpose(totareapigxyt1,(1,2,0)) # time x depth x lat?
    wgtspig1 = totareapigx1 / totareapigxyt1
    print 'wgtspig1.shape: ' + str(wgtspig1.shape)
    print 'wgtspig2.shape: ' + str(wgtspig2.shape)
    print 'tempfldpigzm1.shape: ' + str(tempfldpigzm1.shape)

    # then average meridionally
    avgpig1,ret1 = ma.average(tempfldpigzm1[:,:,np.logical_and(onelon<= Nlim,onelon>Slim)],
                              axis=2,weights=wgtspig1,returned=True)
    avgpig2,ret2 = ma.average(tempfldpigzm2[:,:,np.logical_and(onelon<= Nlim,onelon>Slim)],
                              axis=2,weights=wgtspig2,returned=True)


    """

    
