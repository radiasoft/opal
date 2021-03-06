__author__ = 'swebb'

import numpy as np
#from numba import jit

class tent_dfes:


    def __init__(self, pd):

        self.type = 'spectral'
        self.ptcl_size = pd['particle size']
        self.charge2mass = pd['charge']/pd['mass']
        self.charge = pd['charge']
        self.rho = 0.


    def add_field(self, fields):
        """ Add a field to the interpolater/depositer.

        Args:
            fields (fields class): The corresponding field. Must be spectral
             or else it will raise an exception.
        """

        self.fields = fields
        if not self.fields.get_type() == 'spectral':
            msg = 'Fields' + self.fields.get_name()
            msg += 'is not spectral. You must use matching field and ' \
                   'interpolater_depositer types.'
            raise Exception(msg)
        self.k_modes = self.fields.get_kvectors()

        # This shape function is specific for a tent particle with finite
        # transverse size.
        self.shape_function = np.ones(np.shape(self.k_modes)[0])
        for dim in range (0, np.shape(self.k_modes)[1]):
            for idx in range(0, np.shape(self.k_modes)[0]):
                form_factor = \
                    (2.-2.*np.cos(self.ptcl_size[dim]*
                                  self.k_modes[idx,dim]))/ \
                    (self.ptcl_size[dim]*self.k_modes[idx,dim])**2
                if np.isnan(form_factor):
                    form_factor = 1.
                self.shape_function[idx] *= form_factor


    def compute_forces(self, pos, vel, weights):
        """ Calculates the forces for each particle based on the fields.
        For the discrete Fourier electrostatic class, this is a simple
        electric field acceleration that includes no magnetic rotation.

        Args:
            pos (numpy array): an array of the particle positions
            vel (numpy array): an array of the particle velocities. This is
            the standard API for interps/deps, but is not used for
            electrostatic systems.

        Returns:
            fields (numpy array): a numpy array of the fields computed from
            the fields. This must be multiplied by the charge-to-mass ratio
            of the particles to compute the force.
        """

        phimodes = self.fields.compute_fields(self.rho)
        self.rho = 0.

        efield = -1.j*\
                 np.array(np.matrix(
                     phimodes*self.shape_function*
                     np.exp(1.j*(np.dot(pos, self.k_modes.T))))*\
                 np.matrix(self.k_modes))

        for idx in xrange(0, np.shape(weights)[0]):
            efield[idx] *= weights[idx]/abs(weights[idx])

        # Truncate imaginary part that is lingering after sum due to roundoff
        acceleration = self.charge2mass*(efield.real)

        return acceleration


    def deposit_sources(self, pos, vel, weight):
        """ Calculates the source terms required for the fields and passes
        them off to the fields. Because this class is electrostatic, it only
        computes the charge density.

        Args:
            pos (numpy array): an array of the particle positions
            vel (numpy array): an array of the particle velocities. This is
            the standard API for interps/deps, but is not used for
            electrostatic systems.
        """

        self.rho += self.charge*np.dot(weight,
                        np.exp(-1.j*np.dot(pos, self.k_modes.T))*self.shape_function)

    def get_rho(self):

        return self.rho


    def compute_energy(self):

        rhotilde = np.conj(self.get_rho())
        phitilde = self.fields.get_fields()
        rhophi = np.sum(rhotilde*phitilde)

        return rhophi.real


    def reset(self):
        """ Resets to a pre-step case
        """

        self.rho = 0.