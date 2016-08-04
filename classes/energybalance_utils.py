"""
   energybalance_utils.py
   July 7, 2016

   Functions to compute energy-balance related diagnostics,
     e.g. poleward heat transport

   Started from funcs originally written in PHT_compare.ipynb
"""
import constants as con
import numpy as np


erad = con.get_earthrad() # m


# ########### modified from Brian Rose's notes:
# http://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture13%20--%20Heat%20transport.html#section5 
def inferred_heat_transport( energy_in, lat_deg ):
    '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
    from scipy import integrate

    lat_rad = np.deg2rad( lat_deg )
    return ( 1E-15 * 2 * np.math.pi * erad**2 * 
            integrate.cumtrapz( np.cos(lat_rad)*energy_in,
            x=lat_rad, initial=0. ) )

def calc_inferred_pht(energy,lat):
    """ 
        energy: energ in, e.g. net top of atmosphere energy imbalance, net surface imbalance
        lat: latitude of grid boxes

        returns total pht (in whatever input units are, prob W)
    """
    from scipy import integrate

    latrad = np.deg2rad(lat)

    #print 'calc_inferred_pht() ' + str(2 * np.math.pi * erad**2 * np.cos(latrad)) #@@@@@
    #if nettoa.ndim >1:
    #    # time is included
    #    #retpht = 2 * np.math.pi * erad**2 * np.cumsum(nettoa*np.cos(latrad),axis=1)
    #    retpht = 2 * np.math.pi * erad**2 * integrate.cumtrapz( np.cos(latrad)*nettoa,
    #        x=latrad,initial=0. )
    #else:
    #    retpht = 2 * np.math.pi * erad**2 * np.cumsum(nettoa*np.cos(latrad))
    
    retpht = 2 * np.math.pi * erad**2 * integrate.cumtrapz( np.cos(latrad)*energy,
                                                            x=latrad,initial=0. )            
    return retpht

#def calc_inferred_ocnpht(netsfc,lat):
    """ SAME AS calc_inferred_pht(). No need for 2 funcs. @@@

        netsfc: net surface energy imbalance

        returns ocean pht (in whatever input units are, prob W)
    """
"""    from scipy import integrate

    latrad = np.deg2rad(lat)

    #if netsfc.ndim >1:
    #    # time is included
    #    retpht = 2 * np.math.pi * erad**2 * np.cumsum(netsfc*np.cos(latrad),axis=1)
    #else:
    #    retpht = 2 * np.math.pi * erad**2 * np.cumsum(netsfc*np.cos(latrad))

    retpht = 2 * np.math.pi * erad**2 * integrate.cumtrapz( np.cos(latrad)*netsfc,
                                                            x=latrad,initial=0. )
            
    return retpht
"""
def calc_inferred_pht_components(nettoa,netsfc,lat):
    """ 
        nettoa: net top of atmosphere energy imbalance
        netsfc: net surface energy imbalance
          dims can be time x lat, or just lat
        
        returns total pht, ocn pht, atmos pht as residual
                (in whatever input units are, prob W)
    """

    totpht = calc_inferred_pht(nettoa,lat)
    ocnpht = calc_inferred_pht(netsfc,lat)
    atmpht = totpht - ocnpht

            
    return totpht,ocnpht,atmpht
