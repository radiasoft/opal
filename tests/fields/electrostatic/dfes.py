__author__ = 'swebb'

try:

    from opal.fields import discrete_fourier_electrostatic as dfes
    from matplotlib import pyplot as plt
    import numpy as np

    # test has to emulate the interpolater/depositer regarding shape functions
    # Should create two point charges of different weight at two different
    # locations, with charges 1. and 1.5.
    weights = np.array([10., 15.])
    pos = np.zeros((2,2))
    vel = np.zeros((2,2))
    pos[0,0] = 0.
    pos[0,1] = 0.
    pos[1,0] = -0.35
    pos[1,1] = -0.

    ptclsize = 0.01

    size_array = [ptclsize]*2
    params_dict = {}
    params_dict['n_modes'] = [10]*2
    params_dict['delta k'] = [0.1/ptclsize]*2
    params_dict['dimensions'] = 2

    my_fields = dfes.discrete_fourier_electrostatic(params_dict)

    # Compute the charge density in the Fourier basis
    shape_function = 1.
    k_modes = my_fields.get_kvectors()
    for dim in range(0, np.shape(k_modes)[1]):
        shape_function *= 2.*(1.-np.cos(size_array[dim]*k_modes[:,dim]))/\
                          (size_array[dim]*k_modes[:,dim])**2

    # computes the Fourier components for each particle
    fourier = np.exp(-1.j*np.dot(pos, k_modes.T))
    fourier *= shape_function
    rho = np.dot(weights, fourier)

    kxarray = k_modes[:,0]
    kyarray = k_modes[:,1]

    KX, KY = np.meshgrid(kxarray, kyarray)

    rhotilde = 0.
    myshapefunction = 4.*(1-np.cos(size_array[0]*KX))*(1-np.cos(size_array[
        1])*KY)/(size_array[0]*KX*size_array[1]*KY)**2

    for idx in range(0, np.shape(pos)[0]):
        rhotilde += np.exp(1.j*pos[idx,0]*KX)*\
                    np.exp(1.j*pos[idx,1]*KY)*\
                    myshapefunction

    plt.contourf(KX, KY, rhotilde.real)
    plt.show()
    plt.clf()
    plt.contourf(KX, KY, rhotilde.imag)
    plt.show()
    plt.clf()

    phi = my_fields.compute_fields(rho)

    x = np.linspace(-2., 2.)
    y = np.linspace(-2., 2.)

    XX, YY = np.meshgrid(x, y)

    phiofx = 0.
    for idx in range(0, np.shape(kxarray)[0]):
        phiofx += phi[idx]*\
                  np.exp(1.j*XX*kxarray[idx])*\
                  np.exp(1.j*YY*kyarray[idx])

    rhoofx = 0.
    for idx in range(0, np.shape(kxarray)[0]):
        rhoofx += rho[idx]*\
                  np.exp(1.j*XX*kxarray[idx])*\
                  np.exp(1.j*YY*kyarray[idx])

    #plt.contourf(XX, YY, phiofx)
    #plt.show()

except:

    print 'dfes failed tests'
    raise