__author__ = 'swebb'

"""Example demonstration of the two-stream instability using the
multisymplectic electrostatic algorithm. The beam is assumed initially cold,
while the plasma has some RMS velocity spread.
"""

from opal.fields import discrete_fourier_electrostatic as dfe
from opal.interpolaters_depositers import tent_dfes as depinterp
from opal.particles import non_rel_ptcl as ptcls
from opal.boundaries import particle_boundaries
from opal.auxiliary import constants

import numpy as np

dimensions = 2

plasma_density = 1.e13 # 10^13 ptcls/cm^2
beam_density = 2.e13

beam_vel = 1.e5 # cm/sec
plasma_temp = 1.e3 # (cm/sec)^2

plasma_frequency = np.sqrt(4.*np.pi*plasma_density*
                           constants.elementary_charge/constants.electron_mass)
omega_p_res = 10 # Steps per plasma frequency
debye_length = plasma_temp/plasma_frequency

beam_width = 2.*debye_length

dt = 1./(omega_p_res*plasma_frequency)

domain_lengths = np.array([20.*debye_length, 10.*debye_length])
ptcls_per_debye = 10
ptcl_width = np.array([debye_length/2., debye_length/2.])

num_ptcls_plasma = np.product(domain_lengths)*plasma_density
num_macro_plasma = int((np.product(
                        domain_lengths)/debye_length**2)*ptcls_per_debye)
# Use one macro_weight for simplicity
macro_weight = num_ptcls_plasma/num_macro_plasma
num_ptcls_beam = domain_lengths[0]*beam_width*beam_density
num_macro_beam = int(num_ptcls_beam/macro_weight)

delta_k = 2.*np.pi/domain_lengths
n_modes = np.round(domain_lengths/ptcl_width)
for idx in range(np.shape(n_modes)[0]):
    n_modes[idx] = int(n_modes[idx])

####
#
# Set up simulation
#
####

# Periodic boundary conditions
class periodic_boundary:

    def __init__(self, lengths):

        self.lengths = np.array(lengths)


    def apply_boundary(self, particles):

        particles.pos[:] = particles.pos[:] % self.lengths

my_boundary = periodic_boundary(domain_lengths)

sim_parameters = {}

sim_parameters['number of particles'] = num_macro_beam+num_macro_plasma
sim_parameters['charge'] = -constants.elementary_charge
sim_parameters['mass'] = constants.electron_mass

sim_parameters['dimensions'] = dimensions
sim_parameters['dt'] = dt

# Number of particles per macroparticle
sim_parameters['macro weight'] = macro_weight
sim_parameters['particle size'] = ptcl_width

# Field parameters
sim_parameters['n_modes'] = [n_modes]*dimensions# 20 modes/dimension
sim_parameters['delta k'] = delta_k

# Create the depositer/interpolater, particles, and field solvers

the_depinterp = depinterp.tent_dfes(sim_parameters)
the_particles = ptcls.non_rel_ptcl(sim_parameters)
the_fields    = dfe.discrete_fourier_electrostatic(sim_parameters)
the_boundary  = particle_boundaries.particle_boundaries(sim_parameters)
the_boundary.add_boundary(my_boundary)
the_depinterp.add_field(the_fields)