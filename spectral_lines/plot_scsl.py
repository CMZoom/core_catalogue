# Simple plot of a spectrum. Will plot the spatially-averaged spectrum for a
# single core and a single line.
#
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

# Import the table containing the spatially-averaged spectra

table_dir = "/Users/Jonathan/Dropbox/Work/CMZ/CMZ_SMA/output_spec/"
region = "G0.489+0.010"
linelist = ['C18O']

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

# Plotting
fig   = plt.figure(figsize=( 8.0, 8.0))
ax = fig.add_subplot(111)
fig.suptitle(linelist[0], fontsize=20)
ax.set_xlabel('Velocity  [km/s]')
ax.set_ylabel('Mean flux [Jy/beam]')
ax.set_xlim([np.min(x_axis),np.max(x_axis)])
ax.set_ylim([np.min(y_axis)-stddev_y_axis,np.max(y_axis)+stddev_y_axis])
plt.plot(x_axis, y_axis, 'k-', drawstyle='steps', label=str(structure_to_plot))

# Create a legend
legend = ax.legend(loc='upper left', numpoints=1, borderpad=1., markerscale=8,frameon=False, ncol=2)
for label in legend.get_texts():
    label.set_fontsize(10)
for label in legend.get_lines():
    label.set_linewidth(2.5)  # the legend line width

plt.show()
