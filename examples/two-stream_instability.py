__author__ = 'swebb'

"""Example demonstration of the two-stream instability using the
multisymplectic electrostatic algorithm. The beam is assumed initially cold,
while the plasma has some RMS velocity spread.

Recall that all units are CGS.
"""

from opal.fields import discrete_fourier_electrostatic as dfe
from opal.interpolaters_depositers import tent_dfes as depinterp
from opal.particles import non_rel_ptcl as ptcls
from opal.boundaries import particle_boundaries
from opal.auxiliary import constants
import random

import numpy as np

from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rc('font',**{'family':'serif','serif':['Palatino'], 'size':16})
mpl.rc('text', usetex=True)

dimensions = 2

plasma_density = 1.e17
beam_density = 1.e17

beam_vel = 1.e5 # cm/sec
plasma_temp = 1.e20 # (cm/sec)^2

plasma_frequency = np.sqrt(4.*np.pi*plasma_density*
                           constants.elementary_charge/constants.electron_mass)
omega_p_res = 10 # Steps per plasma frequency
debye_length = np.sqrt(4.*np.pi*plasma_density*constants.elementary_charge**2/
                       (0.5*constants.electron_mass*plasma_temp))

n_plasma_periods = 20
n_steps = n_plasma_periods*omega_p_res

beam_width = 0.5*debye_length

dt = 1./(omega_p_res*plasma_frequency)

domain_lengths = np.array([10.*debye_length, 10.*debye_length])

ptcls_per_debye = 10
ptcl_width = np.array([2*debye_length/ptcls_per_debye,
                       2*debye_length/ptcls_per_debye])

num_ptcls_plasma = np.product(domain_lengths)*plasma_density
num_macro_plasma = int((np.product(
                        domain_lengths)/debye_length**2)*ptcls_per_debye)
# Use one macro_weight for simplicity
macro_weight = num_ptcls_plasma/num_macro_plasma

num_ptcls_beam = domain_lengths[0]*beam_width*beam_density
num_macro_beam = int(num_ptcls_beam/macro_weight)

delta_k = 2.*np.pi/domain_lengths
modes = np.round(2*domain_lengths/ptcl_width)
n_modes = np.zeros(np.shape(modes), dtype='int')
for idx in range(np.shape(n_modes)[0]):
    n_modes[idx] = int(modes[idx])

plot_potential = True
dump_period = omega_p_res/4

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
sim_parameters['n_modes'] = n_modes
sim_parameters['delta k'] = delta_k

print 'Setting up simulation:'
for key in sim_parameters.keys():
    print key, '=', sim_parameters[key]

# Create the depositer/interpolater, particles, and field solvers

the_depinterp = depinterp.tent_dfes(sim_parameters)
the_particles = ptcls.non_rel_ptcl(sim_parameters)
the_fields    = dfe.discrete_fourier_electrostatic(sim_parameters)
the_boundary  = particle_boundaries.particle_boundaries(sim_parameters)
the_boundary.add_boundary(my_boundary)
the_depinterp.add_field(the_fields)

# Create particles from the parameters
for idx in range(0, num_macro_plasma):
    pos = np.array([random.random()*domain_lengths[0],
                    random.random()*domain_lengths[1]])
    vel = np.array([random.gauss(0., plasma_temp),
                    random.gauss(0., plasma_temp)])
    the_particles.add_particle(pos, vel, macro_weight)

for idx in range(0, num_macro_beam):
    pos = np.array([random.random()*domain_lengths[0],
                    (0.5 - random.random())*beam_width+0.5*domain_lengths[1]])
    vel = np.array([random.gauss(0., plasma_temp),
                    random.gauss(0., plasma_temp)])
    the_particles.add_particle(pos, vel, macro_weight)

####
#
# Run the simulation loop
#
####

x = np.arange(0., domain_lengths[0], 0.1*debye_length)
y = np.arange(0., domain_lengths[1], 0.1*debye_length)

XX, YY = np.meshgrid(x, y)

the_particles.half_move_back()

for step in range(0, n_steps):
    print 'Running time step', step
    the_particles.move()

    the_depinterp.deposit_sources(the_particles.pos,
                                  the_particles.vel,
                                  the_particles.weights)

    if step%dump_period == 0:

        print 'Generating plot'

        if plot_potential:
            #phi = the_fields.get_fields()
            rhotilde = the_depinterp.get_rho()
            kvecs = the_fields.get_kvectors()
            rho = 0.
            # compute the test charge force at 1 cm away from the point charge.
            for idx in range(0, np.shape(kvecs)[0]):
                rho += \
                    rhotilde[idx]*np.exp(1.j*(XX*kvecs[idx,0]+YY*kvecs[idx,1]))

            plt.imshow(rho.real,
                       extent=[0., domain_lengths[0],
                               0., domain_lengths[1]],
                       origin='lower',
                       cmap=mpl.cm.bone_r)
            plt.colorbar()
            pltname = 'rho_%07d.png' % step
            plt.savefig(pltname)

            plt.clf()

    the_particles.accelerate(the_depinterp)

the_particles.half_move_forward()
