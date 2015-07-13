__author__ = 'swebb'

failed = False

try:
    from opal.interpolaters_depositers import tent_dfes
    import numpy as np
    from matplotlib import pyplot as plt

    weights = np.array([10., 15.]) #, 5.])


    ptclsize = 0.01
    nmodes = 2
    dimension = 2

    pos = np.zeros((2,dimension))
    vel = np.zeros((2,dimension))
    pos[0,0] = 0.35
    pos[0,1] = 0.
    pos[1,0] = 0.
    pos[1,1] = -0.15

    size_array = [ptclsize]*dimension
    params_dict = {}
    params_dict['n_modes'] = [nmodes]*dimension
    params_dict['delta k'] = [0.1/ptclsize]*dimension
    params_dict['dimensions'] = dimension
    params_dict['particle size'] = [ptclsize]*dimension

    # Fake class to emulate the spectral field solver
    class dummy_fields:

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
                self.k_vectors[idx] *= 2.*np.pi
                self.k_squared[idx] = np.dot(self.k_vectors[idx],
                                             self.k_vectors[idx])

        def get_kvectors(self):

            return self.k_vectors

        def compute_fields(self, rho):
            self.phi = -4.*np.pi*rho/(self.k_squared)

            return self.phi

        def get_fields(self):

            return self.phi

        def get_type(self):

            return 'spectral'

    my_dummy_field = dummy_fields(params_dict)
    my_depositer = tent_dfes.tent_dfes(params_dict)
    my_depositer.add_field(my_dummy_field)
    k_vectors = my_dummy_field.get_kvectors()
    rho = my_depositer.deposit_sources(pos, vel, weights)

    rho_expected = np.array([  8.70576717e-14 +4.30717519j,
                               2.81915491e-14+23.01316777j,
                               2.81915491e-14 +4.60263355j,
                               8.70576717e-14+21.53587596j,
                               5.83559083e-14-23.01316777j,
                               -6.92674107e-15 -4.91835941j,
                               -6.92674107e-15-24.59179705j,
                               5.83559083e-14 -4.60263355j,
                               5.83559083e-14 +4.60263355j,
                               -6.92674107e-15+24.59179705j,
                               -6.92674107e-15 +4.91835941j,
                               5.83559083e-14+23.01316777j,
                               8.70576717e-14-21.53587596j,
                               2.81915491e-14 -4.60263355j,
                               2.81915491e-14-23.01316777j,
                               8.70576717e-14 -4.30717519j])

    testvalue = np.abs(np.dot(rho-rho_expected,rho-rho_expected))/\
                np.abs(np.dot(rho, rho))

    if testvalue > 1.e-16:
        failed = True

    # Test the force calculation

    E, B = my_depositer.compute_forces(pos, vel)

except:
    print 'tent_dfes failed tests'
    raise