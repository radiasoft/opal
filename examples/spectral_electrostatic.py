"""
This is an example file for using the OPAL libraries. This particular
example measures the total energy of fields + particles + coupling for a
Coulomb explosion in two dimensions.
"""
from opal.fields import discrete_fourier_electrostatic as dfe
from opal.interpolaters_depositers import tent_dfes as depinterp

__author__ = 'swebb'
__email__ = 'swebb@radiasoft.net'

from opal.auxiliary import constants
import numpy as np

# Set all simulation parameters at the top for convenience

dimensions = 1
dt = 1.e-10
nsteps = 1

# Particle properties
num_particles = 1000000
macro_weight = 500000
num_macro = num_particles/macro_weight

vel_spread = 1. #cm/s
ball_radius = 0.001 #cm, RMS

macro_size = ball_radius/100.

# Field properties
n_modes = 100
delta_k = 0.1/macro_size

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
sim_parameters['delta k'] = 1/sim_parameters['particle size']

# Create the depositer/interpolater, particles, and field solvers

the_depinterp = depinterp.tent_dfes(sim_parameters)
the_particles = ptcls.non_rel_electrostatic(sim_parameters)
the_fields    = dfe.discrete_fourier_electrostatic(sim_parameters)
the_depinterp.add_field(the_fields)

# generate the particle distribution
#for idx in range(0, num_macro):
#    x  = 1.77777#random.gauss(0., ball_radius)
#    y  = random.gauss(0., ball_radius)
#    vx = random.gauss(0., vel_spread)
#    vy = random.gauss(0., vel_spread)
#    pos = [x] #[x, y]
#    vel = [vx] #[vx, vy]
#    the_particles.add_particle(pos, vel)

pos = [0.]
vel = [0.]
the_particles.add_particle(pos, vel)
pos = [1.]
vel = [0.]
the_particles.add_particle(pos, vel)

# Run the simulation
ptcl_history = []
KE = []
t = []
the_particles.half_move_back()
for idx in range(0, nsteps):
    the_particles.move()
    the_particles.deposit(the_depinterp)
    the_particles.accelerate(the_depinterp)
    if idx%10 == 0:
        pos, vel = the_particles.get_particles()
        #plt.scatter(pos[:,0], pos[:,1])
        #plt.show()
the_particles.half_move_forward()

