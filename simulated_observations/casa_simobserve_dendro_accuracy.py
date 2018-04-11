#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 14:53:18 2017

@author: hph
"""

import pyfits
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.ndimage
import astrodendro
from astropy.utils.console import ProgressBar
from astropy import wcs
from astropy import units as u
from astropy.table import Column
from astropy.convolution import convolve_fft, Gaussian2DKernel
import warnings
#import scipy.ndimage
errmsgs = np.seterr(all='ignore') # silence warning messages about div-by-zero

path = '/Users/hph/current_fits/simulated_observations/'

def run_dendro(data,header,noise,delta,value,pp,pm,low_noise):
    ##Calculate the dendrogram
    dend = astrodendro.Dendrogram.compute(data,
                                    min_value=value*low_noise,
                                    min_delta=delta*low_noise,
                                    min_npix=25)
    #PixelAreaArcsec = 3600. * abs(header['CDELT1']) * 3600. * abs(header['CDELT2'])
    
    ##Compile Metadata
    Metadata = {}
    Metadata['data_unit'] = u.Jy / u.sr
    Metadata['spatial_scale'] =  PixelAreaArcsec**0.5 * u.arcsec
    
    ##Actually build the catalog from the initial dendrogram
    Catalogue = astrodendro.pp_catalog(dend, Metadata)
    Catalogue.rename_column('_idx', 'index')
    ##Add noise column to catalog
    keys = ['noise', 'is_leaf', 'peak_cont_flux', 'min_cont_flux',
            'mean_cont_flux']
    columns = {k:[] for k in (keys)}
    
    for ii, row in enumerate(ProgressBar(Catalogue)):
        structure = dend[row['index']]
        assert structure.idx == row['index'] == ii
        dend_inds = structure.indices()
        columns['noise'].append(noise[dend_inds].mean())
        columns['is_leaf'].append(structure.is_leaf)
        peakflux = data[dend_inds].max()
        columns['peak_cont_flux'].append(peakflux)
        columns['min_cont_flux'].append(data[dend_inds].min())
        columns['mean_cont_flux'].append(data[dend_inds].mean())
        
    for k in columns:
        if k not in Catalogue.keys():
            Catalogue.add_column(Column(name=k, data=columns[k]))
    
    ##Do the pruning
        noise_object = (resid_image-smresid)**2
        print "starting second fft..."
        noise_sqr = convolve_fft(noise_object,  Gaussian2DKernel(2), psf_pad=False, fft_pad=False, allow_huge=True)
        print "fft complete!"
        noise =noise_sqr**0.5
        noise[~np.isfinite(resid_image)] = np.nan
        resid_file[0].data = noise
                
        #Run the dendrogram to get the mask and leafnumber
        
        leafnum1 = run_dendro(newImage,header,noise,min_delta,min_value,pp,pm,low_noise)
        maskImage1 = pyfits.getdata(path+'temp_mask.fits')
        
        #compile found peaks for analysis
        peakFound1 = 0.
        for i in range(n):
            if maskImage1[yVal[i],xVal[i]] >= 0:
                peakFound1 += 1.
        print "The number of peaks found is " +str(peakFound1)
        fraction1 = np.append(fraction1,float(peakFound1)/float(n))
        
        init_Peak += -stepSize
    arrayList1 = np.zeros((fraction1.size,numRuns),dtype = float)
    arrayList1[:,k] = fraction1
    print "Analysis done!"

avg1 = arrayList1[:,0]


k = 1
while k < numRuns:
    for j in range(len(avg1)):
        avg1[j] = avg1[j] + arrayList1[j,k]
    k += 1

## Something is wrong with the scaling here. Off by a factor of numRuns, but I can't figure out why.
## Probably something silly. Let me know if you see it before I do!
avg1 = avg1/float(numRuns)

peakVals = peakVals/(4*24.7)

print "Plotting..."


plt.close()
plt.plot(peakVals, avg1, label='min delta = x8')
#plt.legend(loc='best')
#plt.title('min value = x'+str(min_value))
plt.xlabel('Intensity of Source (Jy/Sr)')
plt.ylabel('Percent Sources Found')

plt.show()
print "Plotting complete!"

plt.savefig(path+'dendro_accuracy_april10.png')

