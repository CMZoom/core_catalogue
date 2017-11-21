import numpy as np
import os
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.table import Column
import matplotlib.pyplot as plt
from matplotlib import rc
import pylab as pl
from astropy import wcs
from astropy.wcs import WCS
import matplotlib.colors as colors
from astrodendro import Dendrogram, pp_catalog
from astrodendro.analysis import PPStatistic

# Files
image_path = os.path.expanduser('~/Dropbox/CMZoom_Data/continuum_images/')
image = 'mosaic_JySr.fits'
#image_j2000 = 'mosaic.fits'
catalog_path = os.path.expanduser('~/Dropbox/CMZoom_Data/prototype_catalog/')
dendrofile = catalog_path+'prototype_dendrogram_nov2017_dv1_dd1_dp10_pp5_pm2_gal.fits'
catalogfile = catalog_path+'prototype_dendrogram_nov2017_datatab_dv1_dd1_dp10_pp5_pm2_gal.fits'
leaftable = catalog_path+'leaf_properties_dv1_dd1_dp10_pp5_pm2_gal.fits'
figure_path = os.path.expanduser('~/Dropbox/SMA_CMZ/figures/catalog/')
figure_name = figure_path+'leaf_map_test_jonny.png'

# Load map
hdu  = fits.open(image_path+image)
hdr  = hdu[0].header
data = hdu[0].data
wcs  = WCS(hdr)

catalog = Table.read(catalogfile)
table = Table.read(leaftable)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

fsize = 14
fweight = 'normal'

fig   = plt.figure(figsize=( 15, 5.0))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)

lon = ax.coords['glon']
lat = ax.coords['glat']
lon.set_axislabel('Galactic Longitude (deg)', minpad=1., size=fsize)
lat.set_axislabel('Galactic Latitude (deg)', minpad=1., size=fsize)
lon.set_ticks(spacing=0.2 * u.deg, color='black', exclude_overlapping=True, size=5, width=1.1)
lat.set_ticks(spacing=0.2 * u.deg, color='black', exclude_overlapping=True, size=5, width=1.1)
lon.set_ticklabel(size=fsize)
lat.set_ticklabel(size=fsize)
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
lon.set_minor_frequency(4)
lat.set_minor_frequency(4)
lon.set_major_formatter('d.d')
lat.set_major_formatter('d.d')

minval = np.nanmin(data)
maxval = np.nanmax(data)

im = ax.imshow(data, interpolation='nearest', cmap='Blues_r', origin='lower', norm=colors.PowerNorm(gamma=1./4.0))

d = Dendrogram.load_from(dendrofile)

leaf_mask = np.zeros(data.shape)
count=0
for leaf in d.leaves:
    foundleaf = (leaf.idx == catalog['index'].data)
    if np.any(foundleaf):
        arr = leaf.get_mask()
        leaf_mask[(arr==True)] = table['index'].data[count]
        count+=1

ax.contour(leaf_mask, transform=ax.get_transform(wcs), cmap='Reds', linewidths=0.25)

plt.savefig(figure_name, dpi=800, bbox_inches='tight')

plt.show()
