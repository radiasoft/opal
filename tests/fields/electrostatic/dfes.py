__author__ = 'swebb'

try:

    from opal.fields import discrete_fourier_electrostatic as dfes
    from matplotlib import pyplot as plt
    import numpy as np

    dimensions = 2
    dk = 0.1
    nmodes = 2
    pd = {}
    pd['dimensions'] = dimensions
    pd['delta k'] = np.array([dk]*dimensions)
    pd['n_modes'] = np.array([nmodes]*dimensions)

    es_solver = dfes.discrete_fourier_electrostatic(pd)

    mykvectors = es_solver.get_kvectors()
    nkvecs = np.shape(mykvectors)[0]
    rho = np.zeros(nkvecs)

    ksquared = np.zeros(nkvecs)

    for idx in range(0, nkvecs):
        ksquared[idx] = np.dot(mykvectors[idx], mykvectors[idx])
        rho[idx] = 1.*idx
        # This is an unphysical rho, but is a simple test
        # that phi computes the right arithmetic

    phi_expected = np.array([  -0.,
                               -12.73239545,
                               -25.46479089,
                               -21.22065908,
                               -50.92958179,
                               -318.30988618,
                               -381.97186342,
                               -89.12676813,
                               -101.85916358,
                               -572.95779513,
                               -636.61977237,
                               -140.05634992,
                               -84.88263632,
                               -165.52114082,
                               -178.25353626,
                               -106.10329539])

    phi = es_solver.compute_fields(rho)

    error = abs(phi - phi_expected)
    for idx in range(0,nkvecs):
        if error[idx] > 1.e-8:
            failed = True


except:

    print 'dfes failed tests'
    raise