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

# Read the Herschel column density file, flatten the data for histograms
hdulist1=fits.open(path+'column_properunits_conv36_source_only.fits')
hdu1=hdulist1[0]
dat1=hdulist1[0].data
herschel=dat1.flatten()

# Create mask
# ....


# Apply mask to Herschel data

# Visually inspect mask on the Herschel file
plt.rcParams.update({'font.size': 16}) #set fontsize
fig=aplpy.FITSFigure(path+'column_properunits_conv36_source_only.fits')
fig.recenter(0.5,0.0, width=2.8, height=1.0)
#fig.show_colorscale(cmap=cm.inferno, vmin=1e+22,vmax=1e+23)
fig.show_colorscale(cmap='Greys_r',vmin=1e+22,vmax=1e+23)
fig.show_colorbar()
fig.colorbar.set_width(0.3)

fig.set_theme('publication')
fig.set_tick_labels_format(xformat='ddd.d', yformat='dd.d')

# Plot contour of mask ...
#fig.show_contour(path+'sgrb2.fits',linewidths='1',colors='cyan',levels=[0])

plt.show()



