"""
This is an example file for using the OPAL libraries. This particular
example measures the total energy of fields + particles + coupling for a
Coulomb explosion in two dimensions.
"""
from opal.fields import discrete_fourier_electrostatic as dfe
from opal.interpolaters_depositers import tent_dfes as depinterp
from opal.particles import non_rel_ptcl as ptcls
from opal.boundaries import particle_boundaries

from matplotlib import pyplot as plt
import matplotlib as mpl
#mpl.rc('font',**{'family':'sans-serif','sans-serif':[
#  'Helvetica']})
mpl.rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rc('text', usetex=True)

__author__ = 'swebb'
__email__ = 'swebb@radiasoft.net'

from opal.auxiliary import constants
import numpy as np

# Set all simulation parameters at the top for convenience

dimensions = 2
dt = 1.e-15
nsteps = 10000

# Particle properties
num_particles = 2
macro_weight = 1
num_macro = num_particles/macro_weight

vel_spread = 1. #cm/s

simulation_lengths = np.array([2., 1.])

# Define the periodic boundary conditions

class periodic_boundary:

    def __init__(self, lengths):

        self.lengths = np.array(lengths)


    def apply_boundary(self, particles):

        particles.pos[:] = particles.pos[:] % self.lengths


my_boundary = periodic_boundary(simulation_lengths)

# Field properties
n_modes = 1
delta_k = 2*np.pi/simulation_lengths
macro_size = 0.01


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

# generate the particle distribution
#for idx in range(0, num_macro):
#    x  = 1.77777#random.gauss(0., ball_radius)
#    y  = random.gauss(0., ball_radius)
#    vx = random.gauss(0., vel_spread)
#    vy = random.gauss(0., vel_spread)
#    pos = [x] #[x, y]
#    vel = [vx] #[vx, vy]
#    the_particles.add_particle(pos, vel)

weight = []

#pos = [0., 0.]
#vel = [0.1, 0.]
#weight.append(1.)
#the_particles.add_particle(pos, vel)
pos = [1.5, 0.9]
vel = [-0.2, 0.]
weight.append(1.)
the_particles.add_particle(pos, vel)
pos = [0.5, 0.1]
vel = [0.2, 0.]
weight.append(1.)
the_particles.add_particle(pos, vel)
weights = np.array(weight)
# Run the simulation
ptcl_history = []
KE = []
t = []
x1 = []
y1 = []
x2 = []
y2 = []
the_particles.half_move_back()
for idx in range(0, nsteps):
    the_particles.move()
    the_depinterp.deposit_sources(the_particles.pos,
                                  the_particles.vel,
                                  weights)
    acceleration = the_depinterp.compute_forces(the_particles.pos,
                                                the_particles.vel)
    the_particles.accelerate(the_depinterp)
    if idx%50 == 0:
        pos, vel = the_particles.get_particles()
        x1.append(pos[1, 0])
        y1.append(pos[1, 1])
        x2.append(pos[0, 0])
        y2.append(pos[0, 1])
        t.append(idx*dt)

    the_boundary.apply_boundary(the_particles)

the_particles.half_move_forward()
#plt.scatter(x, y)
plt.plot(t, x2, label=r'$x_2$')
plt.plot(t, x1, label=r'$x_1$')
plt.xlabel('$t$ [sec]')
plt.ylabel('$x(t)$ [cm]')
plt.legend()
plt.tight_layout()
plt.savefig('periodix_x.png')

plt.clf()

plt.plot(t, y2, label=r'$y_2$')
plt.plot(t, y1, label=r'$y_1$')
plt.xlabel('$t$ [sec]')
plt.ylabel('$y(t)$ [cm]')
plt.legend()
plt.tight_layout()
plt.savefig('periodic_y.png')

plt.clf()

plt.scatter(x1, y1, label=r'ptcl 1', c='r')
plt.scatter(x2, y2, label=r'ptcl 2', c='b')
plt.xlabel('$x$ [cm]')
plt.ylabel('$y$ [cm]')
plt.legend()
plt.tight_layout()
plt.savefig('xy_scatter.png')