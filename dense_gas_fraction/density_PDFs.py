# Create mask of SMA footprint, apply to Herschel data,
# derive PDF of enclosed Herschel data.
# Highlight the pixels that overlap with SMA cores
# from the dendrogram core catalog

import os
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column
from astropy import wcs
from astropy.utils.console import ProgressBar
import reproject
import numpy as np
import matplotlib.pyplot as plt
import aplpy

sma_path = os.path.expanduser('~/Dropbox/SMA_CMZ/CMZoom_Images/November17_continuum_fits/')
herschel_path = os.path.expanduser('~/Dropbox/SMA_CMZ_FITS_files/')
catalog_path = os.path.expanduser('~/Dropbox/SMA_CMZ/prototype_catalog/')
figure_path = os.path.expanduser('~/Dropbox/SMA_CMZ/figures/catalog/')

column_file = os.path.join(herschel_path, 'column_properunits_conv36_source_only.fits')
# column density file handle
column_fh = fits.open(column_file)

temperature_file = os.path.join(herschel_path, 'temp_conv36_source_only.fits')
# temperature density file handle
temperature_fh = fits.open(temperature_file)

assert temperature_fh[0].data.shape == column_fh[0].data.shape


sma_file = os.path.join(sma_path, 'mosaic.fits')
sma_mosaic = fits.open(sma_file)
sma_orig = sma_mosaic[0].data
sma_observed = np.isfinite(sma_orig)

smaobserved_projto_herschel,_ = reproject.reproject_interp((sma_observed,
                                                            wcs.WCS(sma_mosaic[0].header).celestial),
                                                           column_fh[0].header)
smaobserved_projto_herschel[np.isnan(smaobserved_projto_herschel)] = 0
smaobserved_projto_herschel = smaobserved_projto_herschel.astype('bool')

catalog = Table.read(os.path.join(catalog_path, 'mosaic_Nov2017_pruned_v1_datatab.fits'))

colwcs = wcs.WCS(column_fh[0].header)
pix = colwcs.wcs_world2pix(catalog['GLON'], catalog['GLAT'], 0)
column_dens = column_fh[0].data[pix[1].astype('int'), pix[0].astype('int')]
temperature = temperature_fh[0].data[pix[1].astype('int'), pix[0].astype('int')]
catalog.add_column(Column(name='DustTemperature', data=temperature))
catalog.add_column(Column(name='ColumnDensity', data=column_dens))
catalog.write(os.path.join(catalog_path,
                           'mosaic_Nov2017_Jy_per_Ster_pruned_datab_with_ColumnDensity.fits'),
              overwrite=True)

column_masked = column_fh[0].data
column_masked[~smaobserved_projto_herschel] = np.nan

colmin = np.nanmin(column_masked)
colmax = np.nanmax(column_masked)


smasourcemask = np.zeros_like(smaobserved_projto_herschel, dtype='bool')

# for each source, "mask" a region that is within a herschel beam as "this region contains a source"
herschel_beamsize_pixels = (25*u.arcsec / np.mean(wcs.utils.proj_plane_pixel_scales(colwcs)*u.deg)).decompose()

# make a small radial mask
npix = int(np.ceil(herschel_beamsize_pixels))
yy,xx = np.indices([npix*2]*2, dtype='float')
rad = ((yy-npix)**2 + (xx-npix)**2)**0.5
radmask = rad < herschel_beamsize_pixels

pb = ProgressBar(len(catalog))
for row,(cx,cy) in zip(catalog, zip(*pix)):
    cy = int(np.round(cy))
    cx = int(np.round(cx))
    smasourcemask[cy-npix:cy+npix,cx-npix:cx+npix][radmask] = True
    pb.update()



# Visually inspect mask on the Herschel file with the mask as a contour
plt.style.use('seaborn-colorblind')
plt.rcParams.update({'font.size': 16}) #set fontsize

fig = aplpy.FITSFigure(column_file)

fig.show_colorscale(cmap='Greys_r',vmin=colmin,vmax=colmax)
fig.add_colorbar()
fig.colorbar.set_width(0.3)
fig.set_theme('publication')
#fig.set_tick_labels_format(xformat='ddd.d', yformat='dd.d')

# Plot contour of mask ...
fig.show_contour(smaobserved_projto_herschel.astype('int'), linewidths=[1],
                 colors=['b'], levels=[0])

fig.show_contour(smasourcemask.astype('int'), linewidths=[1],
                 colors=['r'], levels=[0])

fig.recenter(0,0,width=2,height=0.75)

plt.savefig(os.path.join(figure_path, 'mask_on_herschelcolumn.pdf'),
            format='pdf', dpi=100, bbox_inches='tight')


# Plot Histogram
plt.rcParams.update({'font.size': 24}) #set fontsize
histfig = plt.figure(2,figsize=(13,6))
ax = histfig.gca()
plt.ylabel('Number of Pixels')
plt.xlabel('Column Density N(H$_2$) [cm$^{-2}$]')

#bins = np.logspace(22.9,23.9,100)
bins = np.logspace(np.log10(colmin),np.log10(colmax),100)
plt.hist(column_masked[np.isfinite(column_masked)], bins,
         color='gray',alpha=0.7, log='True',
         label='Full cloud')
plt.hist(column_masked[smasourcemask & np.isfinite(column_masked)], bins, alpha=0.7, log='True',
         label='SMA Sources')
#plt.hist(G1602[3], bins, color='C0',alpha=0.7, log='True', label=source)
#plt.hist(sourcecol[3], bins, color='magenta', alpha=0.7, log='True', label=source)


#plt.hist(dendro,bins, color='cyan', alpha=0.7, log='True', label='Dendrogram Sources')

ax.set_xscale("log")
#ax.set_xlim(10**22.9, 10**23.9)
ax.set_xlim(colmin, colmax)
#ax.set_ylim(1,200)
ax.set_ylim(1,5000)
plt.gcf().subplots_adjust(bottom=0.2) #make room for x-axis

legend = plt.legend(loc='upper left', shadow=False, fontsize=18)#'x-large')

# Save figure
# Need to save it as a PDF, otherewise, lose transparency
plt.savefig(os.path.join(figure_path, 'column_histograms.pdf'),format='pdf',
            dpi=100, bbox_inches='tight')

plt.show()
