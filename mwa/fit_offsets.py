#Run Agean on image (see Fornax A example) Crossmatch with GLEAM using topcat 
from numpy import *
import matplotlib.pyplot as plt
from astropy.table import Table
from statsmodels.robust.scale import mad

data = Table.read('cross_match.fits')
#data = Table.read('../CenA-field-m-s-n-t.fits')

orig_ras = data['ra']
orig_decs = data['dec']
up_ras = data['RAJ2000']
up_decs = data['DEJ2000']

offsets = data['Separation']

ra_offs = []
dec_offs = []


for source in xrange(len(orig_ras)):
    if abs(offsets[source]) < 80.0:
    	ra_offs.append(orig_ras[source] - up_ras[source])
    	dec_offs.append(orig_decs[source] - up_decs[source])


fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

ra_offs = array(ra_offs)*3600.0
dec_offs = array(dec_offs)*3600.0


med_ra = median(ra_offs)
mad_ra = mad(ra_offs)
med_dec = median(dec_offs)
mad_dec = mad(dec_offs)

min_x, max_x = -40,100

bins = linspace(min_x, max_x, 30)


ax.hist(ra_offs,histtype='step',linewidth=3.5,color='k',normed=True,bins=bins,label='$\mathrm{RA}$ $\mathrm{offset}$')
ax.hist(dec_offs,histtype='step',linewidth=3.5,color='r',normed=True,bins=bins,label='$\mathrm{Dec}$ $\mathrm{offset}$')

ax.axvline(med_ra,color='k',linestyle='--',linewidth=2.0,label='$\mu_{\mathrm{RA}}=%02d\pm%02d$' %(med_ra,mad_ra))
ax.axvline(med_dec,color='r',linestyle='dashed',linewidth=2.0,label='$\mu_{\mathrm{Dec}}=%02d\pm%02d$' %(med_dec,mad_dec))

ax.set_xlabel("Positional offset (arcsecs)",fontsize=16)
ax.set_ylabel("Normalised Count",fontsize=16)
ax.set_xlim(min_x, max_x)


handles,labels = ax.get_legend_handles_labels()

order_handles = [handles[2],handles[0],handles[3],handles[1]]
order_labels = [labels[2],labels[0],labels[3],labels[1]]


ax.legend(order_handles,order_labels,loc='best')

fig.savefig('initial_offset.png',bbox_inches='tight')
