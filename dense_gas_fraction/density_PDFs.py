# Create mask of SMA footprint, apply to Herschel data, 
# derive PDF of enclosed Herschel data.
# Highlight the pixels that overlap with SMA cores
# from the dendrogram core catalog  

# User, update your path for the location of your files locally.
path='/Users/battersby/Work/cmz/cmzoom_catalog_files/'
source='G1.602+0.018'

# packages
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import aplpy

# Set up file names 
column_file=path+'column_properunits_conv36_source_only.fits'
#assumes a filename convention of 'source.continuum.clean.fits'
SMA_file=path+source+'.continuum.clean_regridded.fits'
mask_file=path+source+'_mask.fits'
hist_file=path+source+'_histogram.pdf'

# Read the Herschel column density file#, flatten the data for histograms
columnlist=fits.open(column_file)
column=columnlist[0].data
#dat1=hdulist1[0].data
#herschel=dat1.flatten()

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


####
# Visually inspect mask on the Herschel file
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
fig.show_contour(mask_file,linewidths='1',colors='cyan',levels=[0])

plt.show()
####


# Flatten image to 1D for histogram
col=column_wmask.flatten() 

# Plot Histogram
#plot.style.use('seaborn-colorblind')
plt.rcParams.update({'font.size': 24}) #set fontsize
histogram = plt.figure(1,figsize=(13,6))
ax=histogram.gca()
plt.ylabel('Number of Pixels')
plt.xlabel('Column Density N(H$_2$) [cm$^{-2}$]')

bins = np.logspace(22.0,24.0,100)
plt.hist(col, bins, color='gray', alpha=0.7, log='True', label=source)
#plt.hist(dendro,bins, color='cyan', alpha=0.7, log='True', label='Dendrogram Sources')

plt.gca().set_xscale("log")
ax.set_xlim(10**22.0, 10**24.0)
ax.set_ylim(1,2e3)
plt.gcf().subplots_adjust(bottom=0.2) #make room for x-axis

legend = plt.legend(loc='upper right', shadow=False, fontsize=24)#'x-large')

# Save figure
# Need to save it as a PDF, otherewise, lose transparency
plt.savefig(hist_file,format='pdf', dpi=100, bbox_inches='tight')

plt.show()





####






