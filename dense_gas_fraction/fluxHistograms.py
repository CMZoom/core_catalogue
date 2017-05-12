
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import math



###   Read in FITS files:
    
###   This one is used to make the mask
G1602 = fits.open('/Users/hph/Documents/FITS_Images/G1602_regrid.fits')

###   The following four are masked and then plotted as histograms.
GCMosaic250 = fits.open('/Users/hph/Documents/FITS_Images/gcmosaic_250um.fits')

GCMosaic70 = fits.open('/Users/hph/Documents/FITS_Images/gcmosaic_70um.fits')

GCMosaic24 = fits.open('/Users/hph/Documents/FITS_Images/gcmosaic_24um.fits')

GCMosaic8 = fits.open('/Users/hph/Documents/FITS_Images/gcmosaic_8um.fits')

print "FITS files loaded!"

###   Accesses the FITS images
imageData1 = G1602[0].data
imageData2 = GCMosaic250[0].data
imageData3 = GCMosaic70[0].data
imageData4 = GCMosaic24[0].data
imageData5 = GCMosaic8[0].data

print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

###   function used to build mask:
def buildMask(imageData1):
    maskImage = imageData1[0,0,:,:]
    for pix2 in range(imageData1[0,0,0,:].size):
        for pix1 in range(imageData1[0,0,:,0].size):
            pixel = imageData1[0,0,pix1,pix2]
            if math.isnan(pixel):
                maskImage[pix1,pix2] = 0
            else:
                maskImage[pix1,pix2] = 1
    #plt.imshow(maskImage)
    return maskImage

###   multiply image by mask
def applyMask(maskImage, imageData2):
    return maskImage * imageData2

###   wrapper for summing arrays with NaNs.
def sumPixels(mosaicPixels):
    return np.nansum(mosaicPixels)

###   UNCOMMENT THIS NEXT SECTION TO BUILD A MASK, THEN COMMENT IT
###   AGAIN SO IT RUNS FASTER.
'''
print "Building mask..."
mask1602 = G1602
mask1602[0].data[0,0,:,:] = buildMask(imageData1)
mask1602.writeto('G1602_Mask')
'''

maskFITS = fits.open('G1602_Mask')
mask = maskFITS[0].data[0,0,:,:]

print "Multiplying mosaics by mask..."
GC250_Masked = applyMask(mask,imageData2)
GC70_Masked = applyMask(mask,imageData3)
GC24_Masked = applyMask(mask,imageData4)
GC8_Masked = applyMask(mask,imageData5)


print "The sum of fluxes (in mJy/sr) for the 250um mosaic is:"
print sumPixels(GC250_Masked)
print "The sum of fluxes (in mJy/sr) for the 70um mosaic is:"
print sumPixels(GC70_Masked)
print "The sum of fluxes (in mJy/sr) in for the 24um mosaic is:"
print sumPixels(GC24_Masked)
print "The sum of fluxes (in mJy/sr) for the 8um mosaic is:"
print sumPixels(GC8_Masked)


print "Prepping histograms..."
list1 = filter(lambda k: k != 0 and math.isnan(k) == False, GC250_Masked.flatten())
print "Histogram for first FITS image complete!"
list2 = filter(lambda k: k != 0 and math.isnan(k) == False, GC70_Masked.flatten())
print "Histogram for second FITS image complete!"
list3 = filter(lambda k: k != 0 and math.isnan(k) == False, GC24_Masked.flatten())
print "Histogram for third FITS image complete!"
list4 = filter(lambda k: k != 0 and math.isnan(k) == False, GC8_Masked.flatten())
print "Histogram for fourth FITS image complete!"

fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
n, bins, patches = ax1.hist(list1,bins=255)
ax1.set_title('250um')
ax1.set_xlabel('mJy/sr')
ax1.set_ylabel('# of Occurences')

fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
n, bins, patches = ax2.hist(list2,bins=255)
ax2.set_title('70um')
ax2.set_xlabel('mJy/sr')
ax2.set_ylabel('# of Occurences')

fig3 = plt.figure()
ax3 = fig3.add_subplot(1, 1, 1)
n, bins, patches = ax3.hist(list3,bins=255)
ax3.set_title('24um')
ax3.set_xlabel('mJy/sr')
ax3.set_ylabel('# of Occurences')

fig4 = plt.figure()
ax4 = fig4.add_subplot(1, 1, 1)
n, bins, patches = ax4.hist(list4,bins=255)
ax4.set_title('8um')
ax4.set_xlabel('mJy/sr')
ax4.set_ylabel('# of Occurences')
