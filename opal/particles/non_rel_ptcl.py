__author__ = 'swebb'

import numpy as np

class non_rel_ptcl:


    def __init__(self, pd):

        self.num_ptcls = pd['number of particles']
        self.dimension = pd['dimensions']
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

        acceleration = depinterp.compute_forces(self.pos, self.vel,
                                                self.weights)
        self.vel += self.dt * acceleration


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