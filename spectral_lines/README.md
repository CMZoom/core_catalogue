# spectral_lines

Code descriptions:
==================

Brief description of contents:

spectable.py
------------

Use script_to_make_table_of_spectra.py to generate a table. spectable.py takes
the output from the dendrogram analysis (a fits image with the indices of each
identified structure included), takes the pixels corresponding to each of the
structures from a supplied data cube, and generates a spatially-averaged
spectrum over each structure. This information is output as an astropy table in
a .fits file. The 0th column in the table is the velocity and each subsequent
column corresponds to an identified strucuture and gives the mean_flux in each
channel. Note that in script_to_make_table_of_spectra.py I have added some lines
of code that could be used if you want to loop over a number of regions.

plot_scsl.py
------------

Will plot the spatially-averaged spectrum of a given line extracted from a core.

plot_scml.py
------------

Will plot multiple spatially-averaged spectra extracted from a core.

plot_mcsl.py
------------

Will plot spatially-averaged spectra extracted from all cores.
