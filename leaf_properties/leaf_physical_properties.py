import numpy as np
import os
from astropy.io import fits
from astropy imports wcs
from astrodendro import Dendrogram, pp_catalog
from astrodendro.analysis import PPStatistic
from astropy.table import Table
from astropy import units as u
from astropy.table import Column

# Files
catalog_path = os.path.expanduser('~/Dropbox/SMA_CMZ/prototype_catalog/')
dendrofile = catalog_path+'prototype_dendrogram_nov2017.fits'
catalogfile = catalog_path+'mosaic_Nov2017_Jy_per_Ster_pruned_datab_with_ColumnDensity.fits'
outputtable = catalog_path+'leaf_properties_tab.fits'

# Read in catalog
catalog = Table.read(catalogfile)
#catalog.show_in_browser(jsviewer=True)

# Image information and assumptions
distance        = 8340. # distance to GC; Reid et al. 2014
#Temp            = 20.0
Wave            = (3.0e8/226.e9)
Wave0           = 1.3e-3
k0              = 0.899
nu              = 3.e08/Wave
nu0             = 3.e08/Wave0
beta            = 1.75
Kappag2d        = k0*((nu/nu0)**beta)
g2d             = 100.0
Kappa           = Kappag2d / g2d
mu              = 2.8 # express everything in H2

# Constants
G = 6.67408e-11
msun = 1.989e33
mh = 1.6737236e-27
pc2cm = 3.08567758e18
percm2perm = 1.e6
hplanck = 6.63e-34
clight = 2.99792e8
kboltzmann = 1.381e-23
sin1yr = 3.15569e7
arcsec2pc = distance/((360./(2.*np.pi))*60.*60.)

def mass_calc_submm( Wave, Temp, Kappa, Integrated_Flux, Obj_Dist ):

    Obj_Dist = Obj_Dist * pc2cm

    from planck_func import planck_wave

    B = planck_wave( Wave, Temp )

    Mass = (Obj_Dist**2. * Integrated_Flux) / ( Kappa * B )

    Mass = Mass / msun

    return Mass

def column_density(Wave, Temp, Kappa, Flux, mu):

    from planck_func import planck_wave

    B = planck_wave( Wave, Temp )

    N = Flux / (mu * (mh*1.e3) * Kappa * B)

    return N

def number_density_sphere_pc( Mass_sol, Radius_pc, mu ):

    # This subroutine accepts mass in solar masses and radius in pc and calculates the number density.

    Mass = Mass_sol * (msun/1000.0)
    Radius = Radius_pc * (pc2cm/100.0)

    n = Mass / (((4. / 3.)*np.pi) * mu * mh * Radius**3.0)

    # Convert to particles per cubic centimetre

    n = n / percm2perm

    return n

def mass_density_sphere( Mass_sol, Radius_pc ):

    # This subroutine accepts mass in solar masses and radius in pc and calculates the mass density.

    Mass = Mass_sol * (msun/1000.0)
    Radius = Radius_pc * (pc2cm/100.0)

    rho = Mass / (((4. / 3.)*np.pi) * Radius**3.0)

    return rho

def tff_spherical( number_density, mu ):

    # Accepts a number density in units of particles per cubic centimetre and derives the free fall time in yrs

    mass_density = mu * mh * number_density * percm2perm

    tff = np.sqrt( (3. * np.pi) / (32. * G * mass_density) )

    tff = tff / sin1yr # free-fall time in years

    return tff

def planck_wave( Wave, Temp ):

    # Planck function using wavelength

    planck_conv_wave = 1.e-26 * clight / Wave**2.0

    planck = ((2.0*hplanck*clight**2.0)/(Wave**5.0))*(1.0/(np.exp((hplanck*clight)/(Wave*kboltzmann*Temp))-1.0))
    planck = planck/planck_conv_wave

    return planck

# Generate an astropy table to store all of the information
table = Table(meta={'name': 'Leaf Properties'})
headings = ['index', 'mass', 'N', 'Sigma', 'n', 'rho', 'tff']
description = ['', 'g2d=100, beta=1.75', 'beta=1.75, mu=2.8', '', '', '', '']

# Update the table with Herschel-derived temperature
fh = fits.open(os.path.expanduser('~/Dropbox/SMA_CMZ/gcmosaic_temp_conv25.fits'))
hwcs = wcs.WCS(fh[0].header)
for row in catalog:
    # TODO: check this; I'm not sure if x_cen/y_cen are in the appropriate coordinate system
    xpix,ypix = map(int, map(np.round, hwcs.wcs_world2pix(catalog['x_cen'], catalog['y_cen'], 0)[0]))
    row['DustTemperature'] = fh[0].data[ypix, xpix]


# physical Properties

mass = mass_calc_submm(Wave, catalog['DustTemperature'].data, Kappa, catalog['flux_integrated'].data, distance)
N = column_density(Wave, catalog['DustTemperature'].data, Kappa, catalog['peak_cont_flux'].data, mu)
Sigma = N * mh * 1000.0
n = number_density_sphere_pc(mass, catalog['r_eff'].data*arcsec2pc, mu)
rho = mass_density_sphere(mass, catalog['r_eff'].data*arcsec2pc)
tff = tff_spherical(n, mu)

mass = mass * u.Msun
n = n * u.cm**-3
N = N * u.cm**-2
Sigma = Sigma * u.g * u.cm**-2
rho = rho * u.kg * u.m**-3
tff = tff * u.yr

properties = [catalog['index'],mass, N, Sigma, n, rho, tff]

for i in range(len(headings)):
    table[headings[i]] = Column(properties[i], description=description[i])

table.write(outputtable, format='fits', overwrite=True)
#table.write(outputtable.replace(".fits",".ecsv"), format='ascii.ecsv', overwrite=True)
#table = Table.read(outputtable)
#table.show_in_browser(jsviewer=True)
