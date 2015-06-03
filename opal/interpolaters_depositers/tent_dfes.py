__author__ = 'swebb'

import numpy as np

from matplotlib import pyplot as plt

class tent_dfes:


    def __init__(self, params_dictionary):

        self.type = 'spectral'
        self.ptcl_size = params_dictionary['particle size']



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
        self.shape_function = 1.
        for dim in range (0, np.shape(self.k_modes)[1]):
            self.shape_function *= \
            (2. - 2.*np.cos(self.ptcl_size[dim]*self.k_modes[:,dim]))/ \
            (self.ptcl_size[dim]*self.k_modes[:,dim])**2


    def compute_forces(self, pos, vel):
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

        efield = np.zeros(np.shape(vel), dtype=np.complex)

        phimodes = self.fields.get_fields()
        coeffs = phimodes*self.shape_function

        for idx in range(0, np.shape(vel)[0]):
            # compute e^-ikx first, then add the coefficients and vectors
            fourier = np.exp(1.j*(np.dot(pos[idx], self.k_modes.T)))
            fourier *= coeffs
            efield[idx] = np.dot(fourier, self.k_modes)

        efield *= 1.j

        bfield = np.zeros(np.shape(efield))

        return efield, bfield


    def deposit_sources(self, pos, vel, charge):
        """ Calculates the source terms required for the fields and passes
        them off to the fields. Because this class is electrostatic, it only
        computes the charge density.

        Args:
            pos (numpy array): an array of the particle positions
            vel (numpy array): an array of the particle velocities. This is
            the standard API for interps/deps, but is not used for
            electrostatic systems.
        """

        fourier = np.exp(1.j*(np.dot(pos,self.k_modes.T)))
        fourier = np.sum(fourier,axis=0)

        rho = self.shape_function*fourier #charge*

        #plt.plot(self.k_modes, np.real(rho)) #, s=1)
        #plt.plot(self.k_modes, np.imag(rho))
        #plt.ylim([min(rho), max(rho)])
        #plt.show()
        #plt.plot(x, density)
        #plt.show()
        self.fields.compute_fields(rho)