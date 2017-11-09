# Create mask of SMA footprint, apply to Herschel data,
# derive PDF of enclosed Herschel data.
# Highlight the pixels that overlap with SMA cores
# from the dendrogram core catalog

# User, update your path for the location of your files locally.
path='/Users/battersby/Work/cmz/cmzoom_catalog_files/'
# test
# packages
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import aplpy

# Read the Herschel column density file
#column_file=path+'column_properunits_conv36_source_only.fits'
#column_file=path+'column_properunits_conv36_source_only_OVERSAMPLED_with_G0.489+0.010_objmask_SELF_REGRIDDED.fits'
source='G0.489+0.010'
#source='G1.602+0.018'
#source='G359.889-0.093'
sma_list=fits.open(path+source+'.continuum.clean_REGRIDDED_with_G0.489+0.010_objmask_SELF_REGRIDDED.fits')
#sma_list=fits.open(path+source+'.continuum.clean_SELF_REGRIDDED.fits')

column_file=path+'column_properunits_conv36_source_only_OVERSAMPLED_with_'+source+'_objmask_SELF_REGRIDDED.fits'
columnlist=fits.open(column_file)
column=columnlist[0].data


sma_orig=sma_list[0].data
smamask = np.isfinite(sma_orig)

#multiply column map with SMA mask, so only get N(H2) values
# where there is also SMA data
column = column*smamask
#make all 0s nans, otherwise have problems with histograms
column[column==0]=np.nan
columnlist[0].data = column
columnlist.writeto(path+source+'test.fits',clobber=True)

column_flat=column.flatten()
colmin=np.nanmin(column_flat)
colmax=np.nanmax(column_flat)
print('colmin:', colmin)
print('colmax:', colmax)

# Read SMA data file, make mask, multiply by column map and output an array
def get_col_inmask(source):

        # Assumes a certain filename convention, edit here as needed
        #SMA_file=path+source+'.continuum.clean_regridded.fits'
        SMA_file=path+source+'_objmask_SELF_REGRIDDED.fits'
        mask_file=path+source+'_mask.fits'
        hist_file=path+source+'_histogram.pdf'

        # Read the SMA data file
        # Already regridded and projected to have identical file to the
        # target above in CASA, info in README about this
        smalist=fits.open(SMA_file)
        smadat=smalist[0].data

        # Create mask
        mask = np.isfinite(smadat)

        # Write mask to fits file
        masklist = fits.PrimaryHDU(data=mask, header=columnlist[0].header)
        masklist.writeto(mask_file, overwrite=True)

        # Apply mask to Herschel data
        column_wmask = column*mask
        column_mask_flat = column_wmask.flatten()

        # Save all the things you need together
        return source, mask_file, hist_file, column_mask_flat

#assumes a filename convention of 'source.continuum.clean.fits'
#source='G1.602+0.018'
#G1602=get_col_inmask(source)

#source='G0.489+0.010'
#sourcecol=get_col_inmask(source)
sourcecol=get_col_inmask(source)

# Visually inspect mask on the Herschel file with the mask as a contour
plt.style.use('seaborn-colorblind')
plt.rcParams.update({'font.size': 16}) #set fontsize
#fig=aplpy.FITSFigure(column_file)
fig=aplpy.FITSFigure(path+source+'test.fits')
#fig.recenter(0.5,0.0, width=2.8, height=1.0)
#fig.recenter(0.49,0.008, width=0.1, height=0.05)
#fig.show_colorscale(cmap=cm.inferno, vmin=1e+22,vmax=1e+23)
#fig.show_colorscale(cmap='Greys_r',vmin=5e+22,vmax=5e+23)
fig.show_colorscale(cmap='Greys_r',vmin=colmin,vmax=colmax)
fig.show_colorbar()
fig.colorbar.set_width(0.3)
fig.set_theme('publication')
fig.set_tick_labels_format(xformat='ddd.d', yformat='dd.d')
# Plot contour of mask ...
#fig.show_contour(G1602[1],linewidths='1',colors='cyan',levels=[0])
fig.show_contour(sourcecol[1],linewidths='1',color='C0',levels=[0])
#fig.show_contour(sourcecol[1],linewidths='1',colors='magenta',levels=[0])
plt.savefig(path+source+'_cloud_dendro_img.pdf', format='pdf', dpi=100, bbox_inches='tight')

plt.show()

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
