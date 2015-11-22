__author__ = 'swebb'

import numpy as np
#from numba import jit

class non_rel_ptcl:


    def __init__(self, pd):

        self.num_ptcls = pd['number of particles']
        self.dimension = pd['dimensions']
        self.mass = pd['mass']
        self.pos = np.zeros((self.num_ptcls, self.dimension))
        self.vel = np.zeros((self.num_ptcls, self.dimension))
        self.weights = np.zeros(self.num_ptcls)
        self.exists = [False]*self.num_ptcls
        self.dt = pd['dt']


    def half_move_back(self):

        self.pos -= 0.5*self.dt * self.vel


    def half_move_forward(self):

        self.pos += 0.5*self.dt * self.vel


    def move(self):

        self.pos += self.dt * self.vel


    def accelerate(self, depinterp):

        self.vel += self.dt * depinterp.compute_forces(self.pos, self.vel,
                                                self.weights)


    def get_particles(self):

        return self.pos, self.vel


    def add_particle(self, position, velocity, weight):

        added_ptcl = False
        for idx in range(0, self.num_ptcls):
            if not self.exists[idx]:
                self.pos[idx,:] = position[:]
                self.vel[idx,:] = velocity[:]
                self.weights[idx] = weight
                self.exists[idx] = True
                added_ptcl = True
                break

        if not added_ptcl:
            np.append(np.zeros(self.dimension), self.pos)
            np.append(np.zeros(self.dimension), self.vel)
            self.pos[-1,:] = position[:]
            self.vel[-1,:] = velocity[:]
            self.weights[idx] = weight
            self.exists.append(True)

            added_ptcl = True


    def compute_energy(self):
        """ Compute the total kinetic energy of the particles

        Returns:
            KE (scalar) total kinetic energy in ergs
        """

        ke = 0.

        for idx in range(0, np.shape(self.vel)[0]):
            ke += 0.5*self.mass*abs(self.weights[idx])*\
                  np.dot(self.vel[idx],self.vel[idx])

        return ke


#    @jit
    def compute_momentum(self):
        """ Compute the total momentum of the particles

        Returns
            P (vector) total momentum
        """

        p = 0.
        for idx in range(0, np.shape(self.vel)[0]):
            p += self.mass * abs(self.weights[idx])*self.vel[idx]

        return p