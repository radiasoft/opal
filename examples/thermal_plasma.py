__author__ = 'swebb'
"""
Test of the energy conservation
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

# Plasma parameters
plasma_density = 1.e17
plasma_temp = 1.e14 # (cm/sec)^2

plasma_frequency = np.sqrt(4.*np.pi*plasma_density*
                           constants.elementary_charge/constants.electron_mass)
plasma_wavelength = constants.c/plasma_frequency

# Time step resolutions
omega_p_res = 4 # Steps per plasma frequency
dt = 2.*np.pi/(omega_p_res*plasma_frequency)
n_plasma_oscillations = 1
n_steps = omega_p_res*n_plasma_oscillations

domain_lengths = np.array([10.*plasma_wavelength, 10.*plasma_wavelength])

# Macroparticle properties
ptcl_width = np.array([.5*plasma_wavelength, .5*plasma_wavelength])
ptcls_per_plasma_wavelength = 10
num_ptcls_plasma = np.product(domain_lengths)*plasma_density
num_macro_plasma = int((np.product(
                        domain_lengths)/plasma_wavelength**2)*
                       ptcls_per_plasma_wavelength)
# Use one macro_weight for simplicity
macro_weight = num_ptcls_plasma/num_macro_plasma

delta_k = 2.*np.pi/domain_lengths
modes = np.round(2*domain_lengths/ptcl_width)
n_modes = np.zeros(np.shape(modes), dtype='int')
for idx in range(np.shape(n_modes)[0]):
    n_modes[idx] = int(modes[idx])

dump_period = 100
plot_potential = True

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

sim_parameters['number of particles'] = num_macro_plasma
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
x = []
y = []
for idx in range(0, num_macro_plasma):
    pos = np.array([0.2*domain_lengths[0],
                    0.2*domain_lengths[1]])
    vel = np.array([random.gauss(0., plasma_temp),
                    random.gauss(0., plasma_temp)])
    x.append(pos[0]/domain_lengths[0])
    y.append(pos[1]/domain_lengths[1])
    the_particles.add_particle(pos, vel, macro_weight)

####
#
# Run the simulation loop
#
####

x = np.arange(0., domain_lengths[0], 0.05*domain_lengths[0])
y = np.arange(0., domain_lengths[1], 0.05*domain_lengths[1])
XX, YY = np.meshgrid(x, y)

E = []
ptclE = []
fieldE = []
intE = []
t = []

the_particles.half_move_back()

for step in range(0, n_steps):
    if step%10 == 0:
        print 'Running time step', step

    # Update sequence
    the_particles.move()
    the_depinterp.deposit_sources(the_particles.pos,
                                  the_particles.vel,
                                  the_particles.weights)
    the_particles.accelerate(the_depinterp)

    # Histories
    if step%dump_period == 0:

        the_depinterp.reset()
        the_fields.reset()

        the_depinterp.deposit_sources(the_particles.pos,
                              the_particles.vel,
                              the_particles.weights)

        rhotilde = the_depinterp.get_rho()
        the_fields.compute_fields(rhotilde)
        phitilde = the_fields.get_fields()

        if plot_potential:

            kvecs = the_fields.get_kvectors()
            phi = 0.
            # compute the test charge force at 1 cm away from the point charge.
            for idx in range(0, np.shape(kvecs)[0]):
                phi += \
                    phitilde[idx]*np.exp(1.j*(XX*kvecs[idx,0]+YY*kvecs[idx,1]))

            plt.imshow(phi.real,
                       extent=[0., domain_lengths[0],
                               0., domain_lengths[1]],
                       origin='lower',
                       cmap=mpl.cm.bone_r)
            plt.colorbar()
            pltname = 'phi_%07d.png' % step
            plt.savefig(pltname)

            plt.clf()

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

        rhophi = the_depinterp.compute_energy()
        U = the_fields.compute_energy()
        ke = the_particles.compute_energy()

        fieldE.append(U)
        ptclE.append(ke)
        intE.append(rhophi)
        E.append(ke + U + rhophi)
        t.append(step*dt)

the_particles.half_move_forward()

the_depinterp.reset()
the_fields.reset()

the_depinterp.deposit_sources(the_particles.pos,
                      the_particles.vel,
                      the_particles.weights)

rhotilde = the_depinterp.get_rho()
the_fields.compute_fields(rhotilde)
phitilde = the_fields.get_fields()

rhophi = the_depinterp.compute_energy()
U = the_fields.compute_energy()
ke = the_particles.compute_energy()

fieldE.append(U)
ptclE.append(ke)
intE.append(rhophi)
E.append(ke + U + rhophi)
t.append(step*dt)

print 'E =', E

plt.plot(t, E, t, ptclE, t, fieldE, t, intE)
plt.show()