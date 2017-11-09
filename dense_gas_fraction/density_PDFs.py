# Create mask of SMA footprint, apply to Herschel data,
# derive PDF of enclosed Herschel data.
# Highlight the pixels that overlap with SMA cores
# from the dendrogram core catalog

import os
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import aplpy

sma_path = os.path.expanduser('~/Dropbox/SMA_CMZ/CMZoom_Images/November17_continuum_files/')
herschel_path = os.path.expanduser('~/Dropbox/SMA_CMZ_FITS_files/')
catalog_path = os.path.expanduser('~/Dropbox/SMA_CMZ/prototype_catalog/')
figure_path = os.path.expanduser('~/Dropbox/SMA_CMZ/figures/catalog/')

column_file = os.path.join(herschel_path, 'column_properunits_conv36_source_only.fits')
columnlist = fits.open(column_file)
column = columnlist[0].data


sma_file = os.path.join(sma_path, 'mosaic.fits')
sma_mosaic = fits.open(sma_file)
sma_orig = sma_mosaic[0].data
smamask = np.isfinite(sma_orig)

catalog = Table.read(os.path.join(catalog_path, 'mosaic_Nov2017_Jy_per_Ster.fits_datatab.fits'))

#multiply column map with SMA mask, so only get N(H2) values
# where there is also SMA data
column = column*smamask
#make all 0s nans, otherwise have problems with histograms
column[column==0]=np.nan
columnlist[0].data = column
#columnlist.writeto(path+source+'test.fits',clobber=True)

column_flat=column.flatten()
colmin=np.nanmin(column_flat)
colmax=np.nanmax(column_flat)

# Visually inspect mask on the Herschel file with the mask as a contour
plt.style.use('seaborn-colorblind')
plt.rcParams.update({'font.size': 16}) #set fontsize

fig = aplpy.FITSFigure(column_file)

fig.show_colorscale(cmap='Greys_r',vmin=colmin,vmax=colmax)
fig.show_colorbar()
fig.colorbar.set_width(0.3)
fig.set_theme('publication')
fig.set_tick_labels_format(xformat='ddd.d', yformat='dd.d')

# Plot contour of mask ...
fig.show_contour(smamask,linewidths=[1],color='C0',levels=[0])
plt.savefig(os.path.join(figure_path, 'mask_on_herschelcolumn.pdf'),
            format='pdf', dpi=100, bbox_inches='tight')


# Plot Histogram
plt.rcParams.update({'font.size': 24}) #set fontsize
histogram = plt.figure(1,figsize=(13,6))
ax=histogram.gca()
plt.ylabel('Number of Pixels')
plt.xlabel('Column Density N(H$_2$) [cm$^{-2}$]')

#bins = np.logspace(22.9,23.9,100)
bins = np.logspace(np.log10(colmin),np.log10(colmax),100)
plt.hist(column_flat, bins, color='gray',alpha=0.7, log='True', label='Full cloud (oversampled): '+source)
plt.hist(sourcecol[3], bins,alpha=0.7, log='True', label='Dendrogram leaves: '+source)
#plt.hist(G1602[3], bins, color='C0',alpha=0.7, log='True', label=source)
#plt.hist(sourcecol[3], bins, color='magenta', alpha=0.7, log='True', label=source)


#plt.hist(dendro,bins, color='cyan', alpha=0.7, log='True', label='Dendrogram Sources')

plt.gca().set_xscale("log")
#ax.set_xlim(10**22.9, 10**23.9)
ax.set_xlim(colmin, colmax)
#ax.set_ylim(1,200)
ax.set_ylim(1,5000)
plt.gcf().subplots_adjust(bottom=0.2) #make room for x-axis

legend = plt.legend(loc='upper left', shadow=False, fontsize=18)#'x-large')

# Save figure
# Need to save it as a PDF, otherewise, lose transparency
plt.savefig(sourcecol[2],format='pdf', dpi=100, bbox_inches='tight')

plt.show()
