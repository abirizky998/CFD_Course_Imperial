import numpy as np
from algorithm_and_tools import tdma

def coeff_init(N):
    return np.zeros((4,N))

def cds(N, gamma, dx, rho, u, T, TW, TE):
    aP, aE, aW, b = coeff_init(N)
    for i in range(1,N-1):
        aE[i] = gamma / dx - 0.5 * rho * u
        aW[i] = gamma / dx + 0.5 * rho * u
        aP[i] = aE[i] + aW[i] + b[i]

    return tdma(N, aP, aE, aW, b, T, TW, TE)

def uds(N, gamma, dx, rho, u, T, TW, TE):
    aP, aE, aW, b = coeff_init(N)
    for i in range(1, N-1):
        aE[i] = gamma / dx + max(-(rho * u), 0)
        aW[i] = gamma / dx + max((rho * u), 0)
        aP[i] = aE[i] + aW[i] + b[i]

    return tdma(N, aP, aE, aW, b, T, TW, TE)

def pds(N, gamma, dx, rho, u, T, TW, TE, Pe):
    aP, aE, aW, b = coeff_init(N)
    for i in range(1, N-1):
        aE[i] = (gamma / dx) * max((1 - 0.1 * abs(Pe))**5,0) + max(-(rho * u), 0)
        aW[i] = (gamma / dx) * max((1 - 0.1 * abs(Pe))**5,0) + max((rho * u), 0)
        aP[i] = aE[i] + aW[i] + b[i]
        
    return tdma(N, aP, aE, aW, b, T, TW, TE)