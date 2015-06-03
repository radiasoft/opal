"""Class for moving particles moving at nonrelativistic speeds."""

__name__ = 'non_rel_electrostatic'
__author__ = 'Stephen Webb'

import numpy as np
from matplotlib import pyplot as plt

class non_rel_electrostatic:

    def __init__(self, params_dictionary):
        """Class constructor for non_rel_electrostatic. Sets several
        class-wide parameters such as the charge-to-mass ratio and the time
        step.

        Args:
            params_dictionary (dictionary): The parameters dictionary for
            this particle species.
        """

        nparticles = int(params_dictionary['number of particles']/\
                     params_dictionary['macro weight'])
        self.dt = params_dictionary['dt']
        dims = params_dictionary['dimensions']

        self.charge = \
            params_dictionary['charge']*params_dictionary['macro weight']
        self.mass = \
            params_dictionary['mass']*params_dictionary['macro weight']
        self.dt = params_dictionary['dt']

        self.charge2mass = self.charge/self.mass
        self.dtc2m = self.dt*self.charge2mass

        self.pos = np.zeros([nparticles,dims])
        self.vel = np.zeros([nparticles,dims])
        self.exists = [False] * nparticles


    def add_particle(self, _pos, _vel):
        """ Add a particle with the given position and velocity

        Args:
            _pos (numpy array): a dims-dimensional array of positions
            _vel (numpy array): a dims-dimensional array of velocities
        """
        idx = 0
        while self.exists[idx]:
            idx += 1
        self.pos[idx,:] = _pos[:]
        self.vel[idx,:] = _vel[:]
        self.exists[idx] = True


    def move(self):
        """ Move the particles a full step forward in time."""
        self.pos += self.vel*self.dt


    def half_move_back(self):
        """ Move the particles a half step backward in time. This is used at
        the end of a simulation to complete a second-order
        drift-kick configuration."""
        self.pos -= self.vel*0.5*self.dt


    def half_move_forward(self):
        """ Move the particles a half step forward in time. This is used at
        the beginning of a simulation to start a second-order drift-kick
        configuration."""
        self.pos += self.vel*0.5*self.dt


    def accelerate(self, interpolater):
        """ Accelerate the particle velocity from the interpolated forces.
        Assumes no external magnetic forces.

        Args:
            interpolater (interpolater/depositer): The class which
            interpolates the fields from the solver.
        """
        efield, bfield = interpolater.compute_forces(self.pos, self.vel)
#        plt.quiver(self.pos[:,0], self.pos[:,1], efield[:,0], efield[:,1])
#        plt.show()
        self.vel += efield*self.dtc2m


    def deposit(self, depositer):
        """ Deposit the particle charge and current to a field.

        Args:
            depositer (interpolater/depositer): The class which deposits to
            the desired field
        """
        depositer.deposit_sources(self.pos, self.vel, self.charge)

    def get_particles(self):

        return self.pos, self.vel