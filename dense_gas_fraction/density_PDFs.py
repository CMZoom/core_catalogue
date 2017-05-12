# Create mask of SMA footprint, apply to Herschel data, 
# derive PDF of enclosed Herschel data.
# Highlight the pixels that overlap with SMA cores
# from the dendrogram core catalog  

# User, update your path for the location of your files locally.
path='/Users/battersby/Work/cmz/cmzoom_catalog_files/'

# packages
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import aplpy

# Read the Herschel column density file
column_file=path+'column_properunits_conv36_source_only.fits'
columnlist=fits.open(column_file)
column=columnlist[0].data


# Read SMA data file, make mask, multiply by column map and output an array
def get_col_inmask(source):
	
	# Assumes a certain filename convention, edit here as needed
	SMA_file=path+source+'.continuum.clean_regridded.fits'
	mask_file=path+source+'_mask.fits'
	hist_file=path+source+'_histogram.pdf'
	
	# Read the SMA data file
	# Already regridded and projected to have identical file to the
	# target above in CASA, info in README about this
	smalist=fits.open(SMA_file)
	smadat=smalist[0].data

	# Create mask
	mask=smadat
	mask[np.isfinite(mask)]=1
	mask[np.isnan(mask)]=0

	# Write mask to fits file
	# There must be a better / easier way to write a new data file to
	# a FITS file using its header... but I don't know what it is!!
	masklist=columnlist 
	masklist[0].data=mask
	masklist.writeto(mask_file, clobber=True)

	# Apply mask to Herschel data
	column_wmask=column*mask
	column_mask_flat = column_wmask.flatten()
	
	# Save all the things you need together 
	# OMG, definitely a better way to do this!!
	a=(source,mask_file,hist_file,column_mask_flat)
 
	return a

#assumes a filename convention of 'source.continuum.clean.fits'
source='G1.602+0.018'
G1602=get_col_inmask(source)

source='G0.489+0.010'
#G0489=get_col_inmask(source)

# Visually inspect mask on the Herschel file with the mask as a contour
plt.style.use('seaborn-colorblind')
plt.rcParams.update({'font.size': 16}) #set fontsize
fig=aplpy.FITSFigure(column_file)
fig.recenter(0.5,0.0, width=2.8, height=1.0)
#fig.show_colorscale(cmap=cm.inferno, vmin=1e+22,vmax=1e+23)
fig.show_colorscale(cmap='Greys_r',vmin=1e+22,vmax=1e+23)
fig.show_colorbar()
fig.colorbar.set_width(0.3)
fig.set_theme('publication')
fig.set_tick_labels_format(xformat='ddd.d', yformat='dd.d')
# Plot contour of mask ...
#fig.show_contour(G1602[1],linewidths='1',colors='cyan',levels=[0])
fig.show_contour(G1602[1],linewidths='1',color='C0',levels=[0])
#fig.show_contour(G0489[1],linewidths='1',colors='magenta',levels=[0])

plt.show()


# Plot Histogram
plt.rcParams.update({'font.size': 24}) #set fontsize
histogram = plt.figure(1,figsize=(13,6))
ax=histogram.gca()
plt.ylabel('Number of Pixels')
plt.xlabel('Column Density N(H$_2$) [cm$^{-2}$]')

bins = np.logspace(22.5,23.5,100)
plt.hist(G1602[3], bins,alpha=0.7, log='True', label=source)
#plt.hist(G1602[3], bins, color='C0',alpha=0.7, log='True', label=source)
#plt.hist(G0489[3], bins, color='magenta', alpha=0.7, log='True', label=source)


#plt.hist(dendro,bins, color='cyan', alpha=0.7, log='True', label='Dendrogram Sources')

plt.gca().set_xscale("log")
ax.set_xlim(10**22.5, 10**23.5)
ax.set_ylim(1,200)
plt.gcf().subplots_adjust(bottom=0.2) #make room for x-axis

legend = plt.legend(loc='upper right', shadow=False, fontsize=24)#'x-large')

# Save figure
# Need to save it as a PDF, otherewise, lose transparency
plt.savefig(G1602[2],format='pdf', dpi=100, bbox_inches='tight')

plt.show()





####






