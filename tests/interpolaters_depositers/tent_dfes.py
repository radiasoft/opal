__author__ = 'swebb'


try:
    from opal.interpolaters_depositers import tent_dfes
    import numpy as np
    from matplotlib import pyplot as plt

    weights = np.array([10., 15.]) #, 5.])
    pos = np.zeros((2,2))
    vel = np.zeros((2,2))
    pos[0,0] = 0.
    pos[0,1] = 0.
    pos[1,0] = -0.35
    pos[1,1] = -0.15
#    pos[2,0] = 1.
#    pos[2,1] = 0.5

    ptclsize = 0.01
    nmodes = 100

    size_array = [ptclsize]*2
    params_dict = {}
    params_dict['n_modes'] = [nmodes]*2
    params_dict['delta k'] = [0.1/ptclsize]*2
    params_dict['dimensions'] = 2
    params_dict['particle size'] = [ptclsize]*2

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

            return 0.

        def get_type(self):

            return 'spectral'

    my_dummy_field = dummy_fields(params_dict)
    my_depositer = tent_dfes.tent_dfes(params_dict)
    my_depositer.add_field(my_dummy_field)
    k_vectors = my_dummy_field.get_kvectors()
    rho = my_depositer.deposit_sources(pos, vel, weights)

    x = np.linspace(-2., 2.)
    y = np.linspace(-2., 2.)

    XX, YY = np.meshgrid(x, y)

    mydensity = 0.
    total_modes = np.shape(k_vectors)[0]
    for idx in range(0, total_modes):
        phase = XX*k_vectors[idx,0] + YY*k_vectors[idx,1]
        mydensity += rho[idx]*np.exp(-1.j*phase)

    mydensity /= (2*np.pi)**2

    plt.contourf(XX, YY, mydensity.real)
    plt.colorbar()
    plt.show()
    plt.clf()


except:
    raise