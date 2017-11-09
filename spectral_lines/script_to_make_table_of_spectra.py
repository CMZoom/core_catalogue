# This is an example script detailing how to generate an astropy table containing
# the spatially-averaged spectra extracted from each of the dendrogram structures.
#
# Add line names to the list - the code will generate a table for each file.
#
# It has been written this way so that the code can be looped over.

from spectable import make_table
from astropy.table import Table

# User needs to input the name of the region and the name of the mask
region = 'G0.489+0.010'
mask = './masks/G0.489+0.010_objmask.fits'

# User should input a list of lines - a table will be generated for each line.
linelist = ['C18O']

# Run the code
make_table(region, mask, linelist)

# You can read the table in like this
#t = Table.read("G0.489+0.010.cube.C18O.clean.table.fits")

# If you want to loop over multiple regions then you could do something like
# this...

#region = ['G0.489+0.010', 'G0.489+0.010']
#mask = ['G0.489+0.010_objmask.fits', 'G0.489+0.010_objmask.fits']
#linelist = ['C18O', 'H2CO']

#for i in range(len(region)):
#    sma_region = region[i]
#    sma_mask = mask[i]
#    make_table(sma_region, sma_mask, linelist)
