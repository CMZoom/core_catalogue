#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Edited on Thurs June 1st 2017
Compilation of Casa commands to regrid one fits file's image onto another. 

to run, from casa enter the command "execfile('casaRegrid.py')" and it should work ok.
Some errors might pop up but as long as they're not severe its probably ok <read: I don't know what to do about them>.

@author: hph
"""
#from astropy.io import fits
import os
try:
    from casa import casac
except ImportError, e:
    print("Failed to load Casa: \n")
    exit(1)

'''
    Edit the following file names accordingly (all files are expected to be in your working directory):
        originalFITS should be the exact name of the FITS file that will be regridded.
        mosaicFITS should be the exact name of the FITS file that originalFITS will be regridded onto.
        outputFITS will be the file name of the FITS file outputted into your working directory. 
'''

templateFITS = 'G0.489+0.010_objmask.fits'
originalFITS = 'G0.489+0.010_objmask.fits'
outputFITS = 'G0.489+0.010_objmask_SELF_REGRIDDED.fits'

importfits(fitsimage=templateFITS,imagename='templateTemp.image',overwrite=True)
importfits(fitsimage=originalFITS,imagename='originalTemp.image',overwrite=True)

#regrid the images to galactic coordinates
imregrid(imagename='templateTemp.image',template='GALACTIC',output='templateGalactic.image',overwrite=True)
imregrid(imagename='originalTemp.image',template='GALACTIC',output='originalGalactic.image',overwrite=True)

#perform the actual regridding of sma (in galactic) onto the mosaic (also in galactic)
imregrid(imagename='originalGalactic.image',template='templateGalactic.image',output='original_regridded.image',overwrite=True)

#export the new FITS image
exportfits(imagename='original_regridded.image',fitsimage=outputFITS,overwrite=True)
#viewer(outputFITS)

#delete image file artifacts to clear space
os.system("rm -R *.image")

print("Regridding complete!")
