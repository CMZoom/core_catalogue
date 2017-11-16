import numpy as np
import os
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.table import Column
import matplotlib.pyplot as plt
from matplotlib import rc
import pylab as pl

# Files
catalog_path = os.path.expanduser('~/Dropbox/SMA_CMZ/prototype_catalog/')
catalogfile = catalog_path+'mosaic_Nov2017_Jy_per_Ster_pruned_datab_with_ColumnDensity.fits'
leaftable = catalog_path+'leaf_properties_tab.fits'
figure_path = os.path.expanduser('~/Dropbox/SMA_CMZ/figures/catalog/')
figure_name = figure_path+'leaf_properties_mass_histo.png'

# Read in catalogs
catalog = Table.read(catalogfile)
table = Table.read(leaftable)

#catalog.show_in_browser(jsviewer=True)

distance = 8340. # distance to GC; Reid et al. 2014
arcsec2pc = distance/((360./(2.*np.pi))*60.*60.)

mass = table['mass'].data
physprop = table['mass'].data
xtitle = 'Mass (M$_{\odot}$)'
highmass = (mass > 10.)

# Histogram plot
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
fsize = 14
fweight = 'normal'

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10.0, 7.0))
plt.subplots_adjust(bottom=0.2,left=0.15,top=0.9,hspace=0, wspace=0)
ax = plt.subplot(111)

ax.set_ylabel('N', size=fsize, labelpad=15, rotation=0)
ax.set_xlabel(xtitle, size=fsize, labelpad=10)
ax.tick_params(axis='x', which='major', pad=5, labelsize=fsize)
ax.tick_params(axis='y', which='major', pad=5, labelsize=fsize)
ax.set_xscale('log')

nbins = 20
minval = np.nanmin(physprop)
maxval = np.nanmax(physprop)
bins = np.logspace(np.log10(minval), np.log10(maxval), nbins)

# plot
n, bins, patches = ax.hist(physprop[(np.isfinite(physprop))], bins=bins, facecolor='grey', edgecolor='black',alpha=1.0, label='SMA')#

plt.savefig(figure_name, dpi=600, bbox_inches='tight')
plt.show()
