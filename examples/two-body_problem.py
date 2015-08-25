"""
This is an example file for using the OPAL libraries. This particular
example measures the total energy of fields + particles + coupling for a
Coulomb explosion in two dimensions.
"""
from opal.fields import discrete_fourier_electrostatic as dfe
from opal.interpolaters_depositers import tent_dfes as depinterp
from opal.particles import non_rel_ptcl as ptcls
from opal.boundaries import particle_boundaries

import time

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
import scipy.signal as signal

# Set all simulation parameters at the top for convenience

dimensions = 2
dt = 5.e-9
nsteps = 2*10**6
plot_potential = False
plot_diagnostics = True
dump_step=250

# Particle properties
num_particles = 2
macro_weight = 1
num_macro = num_particles/macro_weight

simulation_lengths = np.array([10., 10.])

# Define the periodic boundary conditions

class periodic_boundary:

    def __init__(self, lengths):

        self.lengths = np.array(lengths)


    def apply_boundary(self, particles):

        particles.pos[:] = particles.pos[:] % self.lengths


my_boundary = periodic_boundary(simulation_lengths)

# Field properties
delta_k = 2*np.pi/simulation_lengths
macro_size = 0.25
n_modes = 2*int(simulation_lengths[0]/macro_size)

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

#pos = [0., 0.]
#vel = [0.1, 0.]
#weight.append(1.)
#the_particles.add_particle(pos, vel)
pos = [0.5*(simulation_lengths[0]+1.), 0.5*simulation_lengths[1]+0.1]
vel = [0., 1.e3]
weight = 1.
the_particles.add_particle(pos, vel, weight)
pos = [0.5*(simulation_lengths[0]-1.), 0.5*simulation_lengths[1]-0.1]
vel = [0., -1.e3]
weight = -2.
the_particles.add_particle(pos, vel, weight)
# Run the simulation

# Set up the histories
ptcl_history = []
E = []
rhophi = []
KE = []
U = []
t = []
x1 = []
y1 = []
x2 = []
y2 = []
vx1 = []
vx2 = []
vy1 = []
vy2 = []
mmntmx = []
mmntmy = []
the_particles.half_move_back()
x = np.arange(0., simulation_lengths[0], 0.025*simulation_lengths[0])
y = np.arange(0., simulation_lengths[1], 0.025*simulation_lengths[1])


the_particles.half_move_forward()

# Compute the fields at the end of the time step
the_depinterp.reset()
the_fields.reset()
the_depinterp.deposit_sources(the_particles.pos,
                      the_particles.vel,
                      the_particles.weights)
rhotilde = the_depinterp.get_rho()
the_fields.compute_fields(rhotilde)
phitilde = the_fields.get_fields()

pos, vel = the_particles.get_particles()
total_momentum = the_particles.compute_momentum()

x1.append(pos[1, 0])
y1.append(pos[1, 1])
x2.append(pos[0, 0])
y2.append(pos[0, 1])
vx1.append(vel[0, 0])
vy1.append(vel[0, 1])
vx2.append(vel[1, 0])
vy2.append(vel[1, 1])
ke = the_particles.compute_energy()
rp = the_depinterp.compute_energy()
uu = the_fields.compute_energy()
U.append(uu)
KE.append(ke)
rhophi.append(rp)
E.append(ke+rp+uu)
mmntmx.append(total_momentum[0])
mmntmy.append(total_momentum[1])
t.append(0.)

the_particles.half_move_back()


t_i = time.time()

XX, YY = np.meshgrid(x, y)
for idx in range(0, nsteps):
    the_particles.move()
    the_depinterp.deposit_sources(the_particles.pos,
                                  the_particles.vel,
                                  the_particles.weights)

    the_particles.accelerate(the_depinterp)

    # Histories
    if plot_potential:

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

            #plt.imshow(phi.real,
            #           extent=[0., simulation_lengths[0],
            #                   0., simulation_lengths[1]],
            #           origin='lower',
            #           cmap=mpl.cm.bone_r)
            #plt.colorbar()

            CS1 = plt.contour(XX, YY, phi, colors='red')
            plt.clabel(CS1)
            plt.colorbar()
            analyticresult = the_particles.weights[0]/\
                             np.sqrt((XX-the_particles.pos[0,0])**2 +
                              (YY-the_particles.pos[0,1])**2)

            analyticresult += the_particles.weights[1]/\
                              np.sqrt((XX-the_particles.pos[1,0])**2 +
                               (YY-the_particles.pos[1,1])**2)

            analyticresult*=constants.elementary_charge

            levels = np.linspace(np.min(analyticresult), 0., 10)

            CS2 = plt.contour(XX, YY, analyticresult, colors='blue',
                           levels=levels)
            plt.clabel(CS2)
            plt.colorbar()
            plt.show()

            plt.clf()

            kvecs = the_fields.get_kvectors()
            rho = 0.
            # compute the test charge force at 1 cm away from the point charge.
            for idx in range(0, np.shape(kvecs)[0]):
                rho += \
                    rhotilde[idx]*np.exp(1.j*(XX*kvecs[idx,0]+YY*kvecs[idx,1]))

            #plt.imshow(rho.real,
            #           extent=[0., simulation_lengths[0],
            #                   0., simulation_lengths[1]],
            #           origin='lower',
            #           cmap=mpl.cm.bone_r)
            #plt.colorbar()
            #plt.show()

            #plt.clf()


    the_boundary.apply_boundary(the_particles)

    if idx%dump_step == 0:

        t_f = time.time()

        if not idx == 0:
            t_left =((t_f - t_i) / idx) * nsteps / 60. -  (t_f - t_i)/60.
            print 'Estimated complete in', t_left, 'min.'

        the_particles.half_move_forward()

        # Compute the fields at the end of the time step
        the_depinterp.reset()
        the_fields.reset()
        the_depinterp.deposit_sources(the_particles.pos,
                              the_particles.vel,
                              the_particles.weights)
        rhotilde = the_depinterp.get_rho()
        the_fields.compute_fields(rhotilde)
        phitilde = the_fields.get_fields()

        pos, vel = the_particles.get_particles()

        x1.append(pos[1, 0])
        y1.append(pos[1, 1])
        x2.append(pos[0, 0])
        y2.append(pos[0, 1])
        vx1.append(vel[0, 0])
        vy1.append(vel[0, 1])
        vx2.append(vel[1, 0])
        vy2.append(vel[1, 1])

        total_momentum = the_particles.compute_momentum()

        ke = the_particles.compute_energy()
        rp = the_depinterp.compute_energy()
        uu = the_fields.compute_energy()
        U.append(uu)
        KE.append(ke)
        rhophi.append(rp)
        E.append(ke+rp+uu)
        mmntmx.append(total_momentum[0])
        mmntmy.append(total_momentum[1])
        t.append((idx+1)*dt)

        the_particles.half_move_back()

the_particles.half_move_forward()

if plot_diagnostics:
    print 'plotting trajectories'
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

    #plt.scatter(x1, y1, label=r'ptcl 1', c='r')
    #plt.scatter(x2, y2, label=r'ptcl 2', c='b', alpha=0.5)
    #plt.xlabel('$x$ [cm]')
    #plt.ylabel('$y$ [cm]')
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig('xy_scatter.png')

    #plt.clf()

    #plt.scatter(vx1, vy1, label=r'ptcl 1', c='r')
    #plt.scatter(vx2, vy2, label=r'ptcl 2', c='b', alpha=0.5)
    #plt.xlabel('$v_x$ [cm/sec]')
    #plt.ylabel('$v_y$ [cm/sec]')
    #plt.legend()
    #plt.tight_layout()
    #plt.savefig('vxy_scatter.png')

    #plt.clf()

    print 'plotting energy'

    E0 = E[0]

    for idx in range(0, len(E)):
        E[idx]/=E0
        KE[idx] /= E0
        rhophi[idx] /= E0
        U[idx] /= E0

    plt.plot(t, E, c='0.5')
    plt.xlabel('$t$ [sec]')
    plt.ylabel(r'$\frac{E}{E_0}$')
    plt.tight_layout()
    plt.savefig('energy.png')

    plt.clf()

    px0 = mmntmx[0]
    py0 = mmntmy[0]

    for idx in range(0, len(mmntmx)):
        mmntmx[idx] -= px0
        mmntmy[idx] -= py0

    plt.plot(t, np.log10(mmntmx), label=r'$\Sigma p_x$')
    plt.plot(t, np.log10(mmntmy), label=r'$\Sigma p_y$')
    plt.xlabel('$t$ [sec]')
    plt.ylabel(r'$\log(\Delta \Sigma p)$')
    plt.legend()
    plt.tight_layout()
    plt.savefig('momentum.png')

    plt.clf()

    plt.plot(t, KE, 'darkolivegreen',
             label=r'$\frac{1}{2} m \mathbf{v}^2$')
    plt.plot(t, rhophi, 'cornflowerblue',
             label=r'$\rho \phi$')
    plt.plot(t, U, 'lightsalmon',
             label=r'$\nabla \varphi \cdot \nabla \varphi$')
    plt.legend()
    plt.xlabel('$t$ [sec]')
    plt.ylabel(r'$\frac{E}{E_0}$')
    #plt.tight_layout()
    plt.savefig('energy_breakdown.png')

    plt.clf()

    # Compute the envelopes and plot the energy breakdown that way

    #
    # Find rhophi envelop
    #

    maxTrhophi = []
    minTrhophi = []
    maxrhophi = []
    minrhophi = []

    maxrhophi.append(rhophi[0])
    minrhophi.append(rhophi[0])
    maxTrhophi.append(t[0])
    minTrhophi.append(t[0])

    rhophimaxima = signal.argrelextrema(np.array(rhophi), np.greater)
    rhophiminima = signal.argrelextrema(np.array(rhophi), np.less)

    for keidx in rhophimaxima[0]:
        maxTrhophi.append(t[keidx])
        maxrhophi.append(rhophi[keidx])

    for keidx in rhophiminima[0]:
        minTrhophi.append(t[keidx])
        minrhophi.append(rhophi[keidx])

    maxrhophi.append(rhophi[-1])
    minrhophi.append(rhophi[-1])
    maxTrhophi.append(t[-1])
    minTrhophi.append(t[-1])

    #
    # Find U envelope
    #

    maxTU = []
    minTU = []
    maxU = []
    minU = []

    maxU.append(U[0])
    minU.append(U[0])
    maxTU.append(t[0])
    minTU.append(t[0])

    KEmaxima = signal.argrelextrema(np.array(U), np.greater)
    KEminima = signal.argrelextrema(np.array(U), np.less)

    for keidx in KEmaxima[0]:
        maxTU.append(t[keidx])
        maxU.append(U[keidx])

    for keidx in KEminima[0]:
        minTU.append(t[keidx])
        minU.append(U[keidx])

    maxU.append(U[-1])
    minU.append(U[-1])
    maxTU.append(t[-1])
    minTU.append(t[-1])

    #
    # Find kinetic energy envelope
    #

    maxTKE = []
    minTKE = []
    maxKE = []
    minKE = []
    maxKE.append(KE[0])
    minKE.append(KE[0])
    maxTKE.append(t[0])
    minTKE.append(t[0])

    KEmaxima = signal.argrelextrema(np.array(KE), np.greater)
    KEminima = signal.argrelextrema(np.array(KE), np.less)

    for keidx in KEmaxima[0]:
        maxTKE.append(t[keidx])
        maxKE.append(KE[keidx])

    for keidx in KEminima[0]:
        minTKE.append(t[keidx])
        minKE.append(KE[keidx])

    maxKE.append(KE[-1])
    minKE.append(KE[-1])
    maxTKE.append(t[-1])
    minTKE.append(t[-1])

    plt.plot(t, KE, c='lightsalmon', alpha=0.25)
    plt.plot(maxTKE, maxKE, c='lightsalmon',
             label=r'$\frac{1}{2} m \mathbf{v}^2$')
    plt.plot(minTKE, minKE, c='lightsalmon')

    plt.plot(t, U, c='cornflowerblue', alpha=0.25)
    plt.plot(maxTU, maxU, c='cornflowerblue',
             label=r'$\nabla \varphi \cdot \nabla \varphi$')
    plt.plot(minTU, minU, c='cornflowerblue')

    plt.plot(t, rhophi, c='darkolivegreen', alpha=0.25)
    plt.plot(maxTrhophi, maxrhophi, c='darkolivegreen',
             label=r'$\rho \varphi$')
    plt.plot(minTrhophi, minrhophi, c='darkolivegreen')

    plt.xlabel(r'$t$ [sec]')
    plt.ylabel(r'$\frac{E}{\Sigma E_0^{(i)}}$')

    plt.legend()
    plt.tight_layout()
    plt.savefig('energy_envelopes.png')

    plt.clf()

    r = np.sqrt((np.array(x1) - np.array(x2))**2 +
                (np.array(y1) - np.array(y2))**2)

    plt.plot(t, r)
    plt.xlabel(r'$t$ [sec]')
    plt.ylabel(r'$r(t)$ [cm]')
    plt.savefig('radius.png')

    plt.clf()

t_f = time.time()

print 'simulation took', (t_f - t_i)/60., 'minutes'