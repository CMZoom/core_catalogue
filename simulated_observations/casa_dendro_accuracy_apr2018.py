#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 14:53:18 2017

@author: hph
"""

import pyfits
import matplotlib.pyplot as plt
import numpy as np
import resid_file[0].data = noise
        scipy
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
                                    min_npix=15)
    PixelAreaArcsec = 3600. * abs(header['CDELT1']) * 3600. * abs(header['CDELT2'])
    
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
    #do pruning
    cat_mask = (Catalogue['is_leaf'] &
                (Catalogue['peak_cont_flux']>pp*Catalogue['noise']) &
                (Catalogue['mean_cont_flux']>pm*Catalogue['noise']))
    pruned_ppcat = Catalogue[cat_mask]
    mask = dend.index_map.copy()
    for ii in range(len(Catalogue)):
        if ii not in pruned_ppcat['index']:
            mask[mask == ii] = -1

    pyfits.writeto(path+'temp_mask.fits',mask,header,clobber=True)
    print len(pruned_ppcat['index'])
    return len(pruned_ppcat['index'])



rawImage = pyfits.open(path+'fakecloud_4.fits')[0].data
header = pyfits.open(path+'G1.602+0.018_v0.continuum.nocontsub.robust_plus05.multiscale.1e4niter.HALF_arcsec.autothresh.tclean.fits')[0].header
imagesize = 200


## Set parameters for accuracy probe

## How many total iterations to be averaged together?
numRuns = 1

## Initial dendrogram parameters, scaled by "low_noise"
min_value = 3.
min_delta = 8.

## Pruning parameters, scaled by local noise estimate.
pp = 5.
pm = 2.

## Low noise estimate, made by hand from fake cloud low noise regions
low_noise = 1.*10**6

## The peak flux you want to see in the injected cores
desiredPeakValue = 2.*10**7

## Scaling for max peak, since the flux is smoothed out by Gaussian Kernel 
max_init_Peak = round(desiredPeakValue *4.* 24.7,2)

##The desired minimum peak flux in the injected cores.
desired_Min_Peak = 1.2*10**6
min_Peak = round(desired_Min_Peak * 4.* 24.7,2)

## The desired number of steps between the max and min flux for the injected cores.
numSteps = 4
stepSize = (max_init_Peak - min_Peak)/numSteps

arrayList1 = np.zeros((numSteps,numRuns),dtype = float)


## Begin the iteration
for k in range(numRuns):
    print str(k*100./numRuns) + " percent complete..."
    newImage = np.zeros(rawImage.shape)
    ##initialize coordinates of artificial peaks
    xMax = newImage[0,:].size
    yMax = newImage[:,0].size
    n= 50 ## This n is the number of artificial gaussian peaks to add
    sigma = 4.
    xVal= np.random.randint(low=0., high=imagesize-1., size=n)
    yVal= np.random.randint(low=0., high=imagesize-1., size=n)
    init_Peak = max_init_Peak

    peakVals = np.empty(shape = (0),dtype=float)
    fraction1 = np.empty(shape = (0),dtype=float)
    ##Loop through peak amplitudes
    while init_Peak > min_Peak:
        
        newImage = np.zeros(rawImage.shape)
        peakVals = np.append(peakVals,init_Peak)
        for i in range(n):
            newImage[yVal[i],xVal[i]] += init_Peak
        newImage = scipy.ndimage.filters.gaussian_filter(newImage,4)
        newImage = newImage+rawImage  
        
        ## Casa commands to apply simobserve - simanalyze business
        pyfits.writeto(path+'initial_temp_image.fits',newImage,header,clobber=True)
        simobserve(project='sim_temp',skymodel=path+'initial_temp_image.fits',incenter='230GHz',inwidth='50MHz',antennalist='sma.subcompact.cfg',totaltime='14400s')
        simanalyze(project='sim_temp',analyze=True,showuv=False,showpsf=False,showmodel=False,showdifference=False,niter = 10000,threshold='5mJy')
        exportfits(imagename='/Users/hph/Documents/Python_Scripts/sim_temp/sim_temp.sma.subcompact.noisy.image.flat', fitsimage=path+'final_temp_image.fits',overwrite=True)
        simulatedImage = pyfits.getdata(path+'final_temp_image.fits')
        
        ##generate noisemap from residual
        exportfits(imagename='/Users/hph/Documents/Python_Scripts/sim_temp/sim_temp.sma.subcompact.noisy.residual', fitsimage=path+'temp_residual.fits',overwrite=True)
        resid_file = pyfits.open(path+'temp_residual.fits')
        resid_image = resid_file[0].data[:][:][0][0]
        print "starting first fft..."
        smresid = convolve_fft(np.nan_to_num(resid_image), Gaussian2DKernel(2), psf_pad=False, fft_pad=False, allow_huge=True)
    
        noise_object = (resid_image-smresid)**2
        print "starting second fft..."
        noise_sqr = convolve_fft(noise_object,  Gaussian2DKernel(2), psf_pad=False, fft_pad=False, allow_huge=True)
        print "fft complete!"
        noise =noise_sqr**0.5
        noise[~np.isfinite(resid_image)] = np.nan
                
        ##Run the dendrogram to get the mask and leafnumber
        
        leafnum1 = run_dendro(newImage,header,noise,min_delta,min_value,pp,pm,low_noise)
        maskImage1 = pyfits.getdata(path+'temp_mask.fits')
        
        ##compile found peaks for analysis
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

## Average together the results    
avg1 = arrayList1[:,0]

k = 1
while k < numRuns:
    for j in range(len(avg1)):
        avg1[j] = avg1[j] + arrayList1[j,k]
    k += 1

avg1 = avg1/float(numRuns)

peakVals = peakVals/(4*24.7)


## Build the accuracy plot
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

