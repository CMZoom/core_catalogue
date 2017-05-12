# Simple plot of a spectrum. Will plot the spatially-averaged spectra
# corresponding to multiple spectral lines, but for a single core.
# Note - this will also work for a single line.
#
#

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.pyplot import cm

# Import the table containing the spatially-averaged spectra

table_dir = "/Users/Jonathan/Dropbox/Work/CMZ/CMZ_SMA/output_spec/"
region = "G0.489+0.010"
linelist = ['C18O', 'C18O']

table = Table.read(table_dir+region+".cube."+linelist[0]+".clean.table.fits")

# Get a list of all available structure indices
indices = table.colnames
indices.pop(0)
indices = np.array([int(val) for val in indices])
print "Structure indices: ", indices

structure_to_plot = indices[0] # change this number

# Plot information - update y_axis for different structure
x_axis = np.array([val for val in table['Velocity']])
y_axis = np.array([val for val in table[str(structure_to_plot)]])

stddev_y_axis = np.std(y_axis)
miny = np.min(y_axis)-stddev_y_axis
maxy = np.max(y_axis)+stddev_y_axis

# Plotting
fig   = plt.figure(figsize=( 8.0, 8.0))
ax = fig.add_subplot(111)
ax.set_xlabel('Velocity  [km/s]')
ax.set_ylabel('Mean flux [Jy/beam]')
ax.set_xlim([np.min(x_axis),np.max(x_axis)])
plt.plot(x_axis, y_axis, 'k-', drawstyle='steps')

if len(linelist) > 1.0:
    
    # Generate a new colour for each line
    n = len(linelist)
    colour=iter(cm.rainbow(np.linspace(0,1,n)))

    for i in range(len(linelist)-1):
        c=next(colour)
        # Read in the new table
        table = Table.read(table_dir+region+".cube."+linelist[i+1]+".clean.table.fits")
        # Identify the same structure
        y_axis = np.array([val for val in table[str(structure_to_plot)]])
        # Calculate the y range - set the ylimits to their maximum values
        stddev_y_axis = np.std(y_axis)
        _miny = np.min(y_axis)-stddev_y_axis
        _maxy = np.max(y_axis)+stddev_y_axis
        if _miny < miny:
            miny = _miny
        if _maxy > maxy:
            maxy = _maxy

        # Overplot the new line
        plt.plot(x_axis, y_axis, color = c, linestyle='-', drawstyle='steps')

ax.set_ylim([miny, maxy])

plt.show()
