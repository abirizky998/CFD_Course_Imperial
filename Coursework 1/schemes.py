import numpy as np
from algorithm_and_tools import tdma

# a is current, b is upstream, c is downstream, d is source
def cds(N, gamma, dx, rho, u, temp):
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    for i in range(1,N-1):
        b[i] = gamma / dx - 0.5 * rho * u
        c[i] = gamma / dx + 0.5 * rho * u
        a[i] = b[i] + c[i] + d[i]

    return tdma(N, a, b, c, d, temp)

def uds(N, gamma, dx, rho, u, temp):
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    for i in range(1, N-1):
        b[i] = gamma / dx + max(-(rho * u), 0)
        c[i] = gamma / dx + max((rho * u), 0)
        a[i] = b[i] + c[i] + d[i]

    return tdma(N, a, b, c, d, temp)

def plds(N, gamma, dx, rho, u, temp, Pe):
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    for i in range(1, N-1):
        b[i] = (gamma / dx) * max((1 - 0.1 * abs(Pe))**5,0) + max(-(rho * u), 0)
        c[i] = (gamma / dx) * max((1 - 0.1 * abs(Pe))**5,0) + max((rho * u), 0)
        a[i] = b[i] + c[i] + d[i]
        
    return tdma(N, a, b, c, d, temp)