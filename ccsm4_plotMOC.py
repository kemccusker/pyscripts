"""  
    ccsm4_plotMOC.py

       for geotimescales project: plot MOC in southern ocean, on top
               of temperature isotherms climo (or change?)
       7/15/2014
"""

import numpy as np # for array handling
import numpy.ma as ma # masked array
import scipy as sp # scientific python
import matplotlib.pyplot as plt # for basic plotting


plt.close("all")
plt.ion()

casenamec = 'b40.20th.track1.1deg.006'
casenamep = 'geo2035ensavg' # or 'rcp8_5GHGrem1850'

filenamec = basepath + casenamec + '.pop.ANN.1980-1999.nc'
filenamep = basepath + casenamep + '.pop.ANN.2045-2054.nc'

latauxgrid = cnc.getNCvar(filenamec,'lat_aux_grid')
transreg = cnc.getNCvar(filenamec, 'transport_regions')
moccomp = cnc.getNCvar(filenamec, 'moc_components')
mocz = cnc.getNCvar(filenamec,'moc_z')
zt = nc_varget(filenamec, 'z_t')




# @@@@@ see notebook
