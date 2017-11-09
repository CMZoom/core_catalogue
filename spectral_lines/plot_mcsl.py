# Simple plot. Will plot the spatially-averaged spectra for a single line for
# each core in a particular region.
#
#

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.pyplot import cm

# Import the table containing the spatially-averaged spectra

table_dir = "./tables/"
region = "G0.489+0.010"
linelist = ['C18O']

table = Table.read(table_dir+region+".cube."+linelist[0]+".clean.table.fits")

# Get a list of all available structure indices
indices = table.colnames
indices.pop(0)
indices = [int(val) for val in indices]
print "Structure indices: ", indices

structure_to_plot = indices[0] # change this number

# Plot information - update y_axis for different structure
x_axis = np.array([val for val in table['Velocity']])
y_axis = np.array([val for val in table[str(structure_to_plot)]])

stddev_y_axis = np.std(y_axis[(y_axis != np.nan)])
miny = np.nanmin(y_axis)-stddev_y_axis
maxy = np.nanmax(y_axis)+stddev_y_axis

# Plotting
fig   = plt.figure(figsize=( 8.0, 8.0))
ax = fig.add_subplot(111)
fig.suptitle(linelist[0], fontsize=20)
ax.set_xlabel('Velocity  [km/s]')
ax.set_ylabel('Mean flux [Jy/beam]')
ax.set_xlim([np.min(x_axis),np.max(x_axis)])
plt.plot(x_axis, y_axis, 'k-', drawstyle='steps', label=str(structure_to_plot))

if len(indices) > 1.0:

    # Generate a new colour for each line
    n = len(indices)
    colour=iter(cm.rainbow(np.linspace(0,1,n)))

    for i in range(len(indices)-1):
        c=next(colour)
        structure_to_plot = indices[i+1] # change this number
        # Identify the same structure
        y_axis = np.array([val for val in table[str(structure_to_plot)]])
        # Calculate the y range - set the ylimits to their maximum values
        stddev_y_axis = np.std(y_axis)
        _miny = np.nanmin(y_axis)-stddev_y_axis
        _maxy = np.nanmax(y_axis)+stddev_y_axis
        if _miny < miny:
            miny = _miny
        if _maxy > maxy:
            maxy = _maxy

        # Overplot the new line
        plt.plot(x_axis, y_axis, color = c, linestyle='-', drawstyle='steps', label=str(structure_to_plot))

# Create a legend
legend = ax.legend(loc='upper left', numpoints=1, borderpad=1., markerscale=8,frameon=False, ncol=2)
for label in legend.get_texts():
    label.set_fontsize(10)
for label in legend.get_lines():
    label.set_linewidth(2.5)  # the legend line width

ax.set_ylim([miny, maxy])

plt.show()
