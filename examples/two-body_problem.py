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

simulation_lengths = np.array([6., 5.])

# Define the periodic boundary conditions

class periodic_boundary:

    def __init__(self, lengths):

        self.lengths = np.array(lengths)


    def apply_boundary(self, particles):

        particles.pos[:] = particles.pos[:] % self.lengths


my_boundary = periodic_boundary(simulation_lengths)

# Field properties
n_modes = 100
delta_k = 2*np.pi/simulation_lengths
macro_size = 0.015


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
pos = [4.1, 3.]
vel = [0., 1.e3]
weight = 1.
the_particles.add_particle(pos, vel, weight)
pos = [3.1, 3.]
vel = [0., -1.e3]
weight = -1.
the_particles.add_particle(pos, vel, weight)
# Run the simulation
ptcl_history = []
KE = []
t = []
x1 = []
y1 = []
x2 = []
y2 = []
vx1 = []
vx2 = []
vy1 = []
vy2 = []
the_particles.half_move_back()
x = np.arange(0., simulation_lengths[0], 0.1)
y = np.arange(0., simulation_lengths[1], 0.1)

XX, YY = np.meshgrid(x, y)
for idx in range(0, nsteps):
    the_particles.move()
    the_depinterp.deposit_sources(the_particles.pos,
                                  the_particles.vel,
                                  the_particles.weights)

    acceleration = the_depinterp.compute_forces(the_particles.pos,
                                                the_particles.vel,
                                                the_particles.weights)

    print 'a =', acceleration, 'cm/sec^2'

    the_particles.accelerate(the_depinterp)

    if plot_potential:
        phi = the_fields.get_fields()
        kvecs = the_fields.get_kvectors()
        potential = 0.
        # compute the test charge force at 1 cm away from the point charge.
        phiat1 = np.zeros(2, dtype='complex')
        for idx in range(0, np.shape(kvecs)[0]):
            potential += phi[idx]*np.exp(1.j*(XX*kvecs[idx,0]+YY*kvecs[idx,1]))
            phiat1 += phi[idx]*\
                      np.exp(1.j*((pos[0]-1.)*kvecs[idx,0]+
                                  pos[1]*kvecs[idx,1]))

        phiat1 *= constants.elementary_charge
        print 'phi at 1 cm =', phiat1, 'ergs'

        potential *= constants.elementary_charge

        #plt.contour(XX, YY, potential.real)

        phiAnalytic = constants.elementary_charge**2/\
                      (np.sqrt((XX-pos[0])**2 + (YY-pos[1])**2))

        plt.contour(XX, YY, phiAnalytic)
        plt.colorbar()

        #plt.contourf(XX, YY, potential.real, cmap=mpl.cm.bone_r)

        plt.imshow(potential.real,
                   extent=[0., simulation_lengths[0],
                           0., simulation_lengths[1]],
                   origin='lower',
                   cmap=mpl.cm.bone_r)
        plt.colorbar()
        plt.show()

        plt.clf()

    if idx%100 == 0:
        pos, vel = the_particles.get_particles()
        x1.append(pos[1, 0])
        y1.append(pos[1, 1])
        x2.append(pos[0, 0])
        y2.append(pos[0, 1])
        vx1.append(vel[0, 0])
        vy1.append(vel[0, 1])
        vx2.append(vel[1, 0])
        vy2.append(vel[1, 1])
        t.append(idx*dt)

    the_boundary.apply_boundary(the_particles)

the_particles.half_move_forward()

if plot_diagnostics:
    plt.plot(t, x2, label=r'$x_2$')
    plt.plot(t, x1, label=r'$x_1$')
    plt.xlabel('$t$ [sec]')
    plt.ylabel('$x(t)$ [cm]')
    plt.legend()
    plt.tight_layout()
    plt.savefig('periodic_x.png')

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
    plt.scatter(x2, y2, label=r'ptcl 2', c='b', alpha=0.5)
    plt.xlabel('$x$ [cm]')
    plt.ylabel('$y$ [cm]')
    plt.legend()
    plt.tight_layout()
    plt.savefig('xy_scatter.png')

    plt.clf()

    plt.scatter(vx1, vy1, label=r'ptcl 1', c='r')
    plt.scatter(vx2, vy2, label=r'ptcl 2', c='b', alpha=0.5)
    plt.xlabel('$v_x$ [cm/sec]')
    plt.ylabel('$v_y$ [cm/sec]')
    plt.legend()
    plt.tight_layout()
    plt.savefig('vxy_scatter.png')