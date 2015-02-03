
import pandas as pd
import cccmaplots as cplt

# Agung: Feb. 1963
# Chichon: Apr 1982
# Pinatubo: Jul 1991

# These time selections start with month of eruption
timselagung ='1963-02-01,1965-12-31'
timselchichon = '1982-04-01,1985-12-31'
timselpinatubo = '1991-07-01,1993-12-31'

timselbase = '1960-01-01,1999-12-31'
basedt = le.load_LEdata(fdict,'historicalNat',timesel=timselbase,seas=('ANN',))# annual avg fine for now
basedf=pd.DataFrame(basedt)
tmp = basedf.loc['ANN']
bamean = np.mean(tmp.mean(axis=0), axis=0) # ens mean and time mean



agungdt=le.load_LEdata(fdict,'historicalNat',timesel=timselagung)
chichondt=le.load_LEdata(fdict,'historicalNat',timesel=timselchichon)
pinatubodt=le.load_LEdata(fdict,'historicalNat',timesel=timselpinatubo)

fname = agungdt.keys()[0]

agungtime = cnc.getNCvar(fname,'time',timesel=timselagung) # not really using this yet.
lev = cnc.getNCvar(fname,'plev')
lat = cnc.getNCvar(fname,'lat')

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


