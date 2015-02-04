
import pandas as pd
import cccmaplots as cplt
import datetime as datetime
from netCDF4 import num2date

# Agung: Feb. 1963
# Chichon: Apr 1982
# Pinatubo: Jul 1991

# These time selections start with month of eruption
## timselagung ='1963-02-01,1965-12-31'
## timselchichon = '1982-04-01,1985-12-31'
## timselpinatubo = '1991-07-01,1993-12-31'

# These time selections start with beginning of year of eruption
#  (so can easily remove seasonal cycle)
timselagung ='1963-01-01,1965-12-31'
timselchichon = '1982-01-01,1985-12-31'
timselpinatubo = '1991-01-01,1993-12-31'


## timselbase = '1960-01-01,1999-12-31'
## basedt = le.load_LEdata(fdict,'historicalNat',timesel=timselbase,seas=('ANN',))# annual avg fine for now
## basedf=pd.DataFrame(basedt)
## tmp = basedf.loc['ANN']
## bamean = np.mean(tmp.mean(axis=0), axis=0) # ens mean and time mean


# Need to remove each run's annual cycle
natclimo = le.load_LEdata(fdict,'historicalNat',ftype='1950-2020_climo')
histclimo = le.load_LEdata(fdict,'historical',ftype='1950-2020_climo')

agungdt=le.load_LEdata(fdict,'historicalNat',timesel=timselagung)
chichondt=le.load_LEdata(fdict,'historicalNat',timesel=timselchichon)
pinatubodt=le.load_LEdata(fdict,'historicalNat',timesel=timselpinatubo)

fname = agungdt.keys()[0]

lev = cnc.getNCvar(fname,'plev')
lat = cnc.getNCvar(fname,'lat')
atime = cnc.getNCvar(fname,'time',timesel=timselagung)
adates = cnc.get_NCdates(fname,timesel=timselagung)


#atime = pd.date_range(start='1963/1/1',end='1965/12/31',freq='M')


ctime = cnc.getNCvar(fname,'time',timselchichon)
ptime = cnc.getNCvar(fname,'time',timselpinatubo)



adt=le.load_LEdata(fdict,'historicalNat',timesel=timselagung)
apan=pd.Panel(adt,major_axis=adates)
amon = apan.groupby(apan.major_axis.month) # this returns a pandas.core.groupby.PanelGroupBy object. @@ what to do w/ it??

# @@ Note that these are currently not weighted by month:
agungsr=pd.Series(agungdt)
agmean = agungsr.mean(axis=0) # avg over ensemble, so it's months by height by last
# John's avg is 9 months starting with 5th month after eruption
agcomp = np.mean(agmean[4:13,:,:],axis=0) # agung composite

chichonsr=pd.Series(chichondt)
chmean = chichonsr.mean(axis=0) # avg over ensemble, so it's months by height by last
# John's avg is 9 months starting with 5th month after eruption
chcomp = np.mean(chmean[4:13,:,:],axis=0) # chichon composite

pinatubosr=pd.Series(pinatubodt)
pimean = pinatubosr.mean(axis=0) # avg over ensemble, so it's months by height by last
# John's avg is 9 months starting with 5th month after eruption
picomp = np.mean(pimean[4:13,:,:],axis=0) # pinatubo composite

plt.figure(); cplt.vert_plot(agcomp,lev,lat)
plt.figure(); cplt.vert_plot(chcomp,lev,lat)
plt.figure(); cplt.vert_plot(picomp,lev,lat)


levlim=10
cmin=-4; cmax=4

fig,axs=plt.subplots(1,3)
fig.set_size_inches(12,6)
ax=axs[0]
cplt.vert_plot(agcomp-bamean,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='agung',suppcb=True,levlim=levlim)
ax.set_xlim(-30,30)
ax=axs[1]
cplt.vert_plot(chcomp-bamean,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='chichon',suppcb=True,levlim=levlim)
ax.set_xlim(-30,30)
ax=axs[2]
vp = cplt.vert_plot(picomp-bamean,lev,lat,axis=ax,cmin=cmin,cmax=cmax,title='pinatubo',suppcb=True,levlim=levlim)
ax.set_xlim(-30,30)
# @@@ add full composite

cbar_ax = fig.add_axes([.25,0.07, 0.5, .02])
cbor='horizontal'
fig.colorbar(vp,cax=cbar_ax, orientation=cbor) 


