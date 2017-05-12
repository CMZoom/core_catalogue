# import modules

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy import wcs
import warnings
warnings.simplefilter('ignore', wcs.FITSFixedWarning)

def get_axis(header):
    """
    Generate & return the velocity axis from the fits header.
    """
    mywcs = wcs.WCS(header)
    specwcs = mywcs.sub([wcs.WCSSUB_SPECTRAL])
    return specwcs.wcs_pix2world(np.arange(header['NAXIS{0}'.format(mywcs.wcs.spec+1)]), 0)

def generate_table(data_cube, axis, indices):
    """
    Generate an astropy table to store the information.
    """
    table = Table(meta={'name': data_cube})
    table['Velocity'] = Column(axis, unit='km/s', description = 'Velocity')

    return table

def get_indices(index, mask):
    """
    Retrieve the indices of the current structure from the mask file.
    """
    indices = np.array(np.where(mask == index))
    return indices

def get_flux(cube, indices):
    """
    Get the spatially-averaged flux over the current structure.
    """
    # Split out the x, y indices individually
    idx = indices[0,:]
    idy = indices[1,:]
    # derive spatially-averaged flux per channel
    flux = np.zeros(len(cube[:,0,0]))
    for i in range(len(cube[:,0,0])):
        flux[i] = np.sum(cube[i, idx, idy])
    flux = flux/np.size(indices[0,:])
    return flux

def make_table(region, maskfile, linelist):
    """
    The main routine for creating the astropy table. Note - this relies on the
    naming convention for the spectral cubes as of 12/05/17 e.g.

    G0.489+0.010.cube.C18O.clean.fits

    The code will generate an output file with the following naming convention

    G0.489+0.010.cube.C18O.clean.table.fits

    """
    # create list of cube names
    cube_prefix = '.cube.'
    suffix = '.clean.fits'
    data_cubes = [region+cube_prefix+item+suffix for item in linelist]

    # read maskfile
    hdu = fits.open(maskfile)
    maskhead = hdu[0].header
    mask = hdu[0].data
    hdu.close()

    # Nans give an annoying error - change to a value which will never get used
    nanids = np.isnan(mask)
    mask[nanids] = -1000.0

    for data_cube in data_cubes:
        # Loop over each molecular line data cube in linelist
        # Open the fits file
        hdu = fits.open(data_cube)
        cubehead = hdu[0].header
        cube = hdu[0].data
        hdu.close()
        cube = np.squeeze(cube)
        # Generate the velocity axis
        axis = get_axis(cubehead)
        axis = axis[0]/1000.0 # km/s

        # Retrieve all the structure ids from the mask
        structids = list(set(mask[mask >= 0.0]))
        # Prime the table to store the information
        table = generate_table(data_cube, axis, structids)
        # Prime the output file
        outfile = data_cube.replace(".fits",".table.fits")

        # Loop through all the structures
        for idx in structids:
            indices = get_indices(idx, mask)
            mean_flux = get_flux(cube, indices)
            # Add column to table
            table[str(int(idx))] = Column(mean_flux, unit='Jy/beam', description = 'Mean Flux')

        # Write the table to the current directory
        table.write(outfile, format='fits', overwrite=True)

    return
