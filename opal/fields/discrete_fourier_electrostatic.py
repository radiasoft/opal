__author__ = 'swebb'


import numpy as np


class discrete_fourier_electrostatic:

    def __init__(self, params_dictionary):

        self.type = 'spectral'
        self.n_modes = params_dictionary['n_modes']
        self.dk = params_dictionary['delta k']
        dims = params_dictionary['dimensions']

        # Use linear strides to generate the k-space vectors
        # Positive and negative k-vectors
        for idx in range(dims):
            self.n_modes[idx] *= 2
        num_modes = np.prod(self.n_modes)

        self.k_vectors = np.zeros((num_modes, dims))
        self.k_squared = np.zeros(num_modes)
        kmin = np.zeros(dims)
        for idx in range(0, dims):
            kmin[idx] = -0.5*self.dk[idx]*self.n_modes[idx] + 0.5*self.dk[idx]
        for idx in range(0, num_modes):
            linear_stride = idx
            for dim_idx in range(0, dims):
                this_idx = linear_stride%self.n_modes[dim_idx]
                self.k_vectors[idx,dim_idx] = (this_idx)*self.dk[dim_idx]
                linear_stride -= this_idx
                linear_stride /= self.n_modes[dim_idx]
            self.k_vectors[idx] += kmin
            #self.k_vectors[idx] *= 2.*np.pi
            self.k_squared[idx] = np.dot(self.k_vectors[idx],
                                         self.k_vectors[idx])


    def compute_fields(self, rho):
        """ Adds to the charge distribution for the electrostatic spectral
        code.

        Args:
            rho (np array): Numpy array of the individual Fourier amplitudes
            for the charge density. Should match the indexing properly for
            the k-vectors
        """

        self.phi = -4.*np.pi*rho/self.k_squared
        return self.phi


    def get_fields(self):

        return self.phi


    def get_kvectors(self):
        """Write-safe accessor method to get the k-vector arrays

        Returns:
            k_modes (np array): the
        """

        return self.k_vectors

    def get_type(self):

        return 'spectral'