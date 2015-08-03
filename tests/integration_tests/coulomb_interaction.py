__author__ = 'swebb'

"""
This is a test to determine whether the two-body Coulomb attraction between
an electron and a positron is correct.
"""
from opal.fields import discrete_fourier_electrostatic as dfe
from opal.interpolaters_depositers import tent_dfes as depinterp
from opal.particles import non_rel_ptcl as ptcls
from opal.boundaries import particle_boundaries

from matplotlib import pyplot as plt
import matplotlib as mpl
#mpl.rc('font',**{'family':'sans-serif','sans-serif':[
#  'Helvetica']})
mpl.rc('font',**{'family':'serif','serif':['Palatino'], 'size':16})
mpl.rc('text', usetex=True)

__author__ = 'swebb'
__email__ = 'swebb@radiasoft.net'

from opal.auxiliary import constants
import numpy as np

# Set all simulation parameters at the top for convenience

dimensions = 2
dt = 1.e-10
nsteps = 1#2*10**4
plot_potential = True
plot_diagnostics = False

# Particle properties
num_particles = 2
macro_weight = 1
num_macro = num_particles/macro_weight

simulation_lengths = np.array([30., 30.])

# Define the periodic boundary conditions

class periodic_boundary:

    def __init__(self, lengths):

        self.lengths = np.array(lengths)


    def apply_boundary(self, particles):

        particles.pos[:] = particles.pos[:] % self.lengths


my_boundary = periodic_boundary(simulation_lengths)

# Field properties
n_modes = 250
delta_k = 2*np.pi/simulation_lengths
macro_size = 0.15

# The params_dictionary for the electrostatic field + particles
sim_parameters = {}

# Simulation calls for one million electrons in a Gaussian ball Coulomb
# exploding over time
sim_parameters['number of particles'] = num_particles
sim_parameters['charge'] = -constants.elementary_charge
sim_parameters['mass'] = constants.electron_mass

sim_parameters['dimensions'] = dimensions
sim_parameters['dt'] = dt

# Number of particles per macroparticle
sim_parameters['macro weight'] = macro_weight
sim_parameters['particle size'] = np.array([macro_size, macro_size, macro_size])

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

weight = []

pos = [0.5*(simulation_lengths[0]-1.), 0.5*simulation_lengths[1]+0.1]
vel = [0., 0.]
weight = 1.
the_particles.add_particle(pos, vel, weight)
pos = [0.5*(simulation_lengths[0]+1.), 0.5*simulation_lengths[1]+0.1]
vel = [0., 0.]
weight = -1.
the_particles.add_particle(pos, vel, weight)
# Run the simulation
the_particles.half_move_back()

for idx in range(0, nsteps):
    the_particles.move()
    the_depinterp.deposit_sources(the_particles.pos,
                                  the_particles.vel,
                                  the_particles.weights)

    acceleration = the_depinterp.compute_forces(the_particles.pos,
                                                the_particles.vel,
                                                the_particles.weights)

    the_answer = np.array([[  6.29176115e+07, -3.91382417e-07],
                           [ -6.29176115e+07, -3.32348075e-07]])

    error = the_answer-acceleration
    #right_value =

    #error_metric

    print acceleration