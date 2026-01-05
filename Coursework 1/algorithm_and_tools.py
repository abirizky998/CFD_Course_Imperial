import numpy as np
import matplotlib.pyplot as plt

# a is current, b is upstream, c is downstream, d is source
def tdma(N, a, b, c, d, T, TW, TE):
    P = np.zeros(N)
    Q = np.zeros(N)

    # west BC, T[0] is known so P[0] = 0 and Q[0] is T[0]
    P[0] = 0
    Q[0] = TW
    # setting up the inside elements forwards
    for i in range(1, N-1): 
        P[i] = b[i] / (a[i] - c[i] * P[i-1])
        Q[i] = (d[i] + c[i] * Q[i-1]) / (a[i] - c[i] * P[i-1])
    # east BC, T[N] is known
    P[N-1] = 0
    Q[N-1] = TE
    for i in range(N-2, -1, -1): # backwards solve
        T[i] = P[i] * T[i+1] + Q[i]

    return T.copy()

def analytical(grid, L, Pe, Phi_0, Phi_L):
    return Phi_0 + ((np.exp(grid * Pe / L) - 1) / (np.exp(Pe) - 1) * (Phi_L-Phi_0))

def error_calc(grid, analytical, numerical):
    return 100 * (sum(abs((numerical - analytical)/analytical)) / grid)

def plotting(grid, analytical_results, title, scheme_results, error_value, Pe_G, Pe_L):
    fig, ax = plt.subplots() # values are rounded for display, not processing
    ax.plot(grid, analytical_results, label=f'Analytical', c='black')
    ax.plot(grid, scheme_results, label=f'Numerical, Error:{round(error_value,2)}%', 
            c='blue', ls='--', marker='o', mfc='red') 

    ax.set_title(f"{title}, Pe Global: {round(Pe_G,3)}, Pe Local: {round(Pe_L,3)}")
    ax.set_xlabel('x / L')
    ax.set_ylabel('Temperature (°C)')
    ax.legend()

    return plt.show()

def error_dx(dx, cds_err, uds_err, pds_err, u, gamma, density):
    fig, ax = plt.subplots()
    ax.plot(dx, cds_err, label='CDS', color='blue', marker='o')
    ax.plot(dx, uds_err, label='UDS', color='green', marker='s')
    ax.plot(dx, pds_err, label='PLDS', color='red', marker='^')

    ax.set_title("Grid-Point Spacing vs. Cell Count\n"
    fr"$\rho = {density}$ kg/m$^3$, $\Gamma_\phi = {gamma}$ kg/ms, $u = {u}$ m/s" + "\n")
    ax.set_xlabel(fr"$\delta$ x = domain length / cell count")
    ax.set_ylabel('Error (%)')
    ax.set_yscale('log')
    ax.secondary_xaxis(
        'top', functions=(lambda x: x * density * u / gamma, lambda x: x / (density * u / gamma))
        ).set_xlabel("Local Peclet Number ($Pe_{local}$)")
    ax.grid(True, which="both", ls="-", alpha=0.2)
    
    ax.legend()
    
    return plt.show()

def error_cell(N, cds_err, uds_err, pds_err, u, gamma, density):
    fig, ax = plt.subplots()
    ax.plot(N, cds_err,label='CDS', color='blue', marker='o')
    ax.plot(N, uds_err,label='UDS', color='green', marker='s')
    ax.plot(N, pds_err,label='PLDS', color='red', marker='^')

    ax.set_title("Grid-Point vs Error Value\n"
    fr"$\rho = {density}$ kg/m$^3$, $\Gamma_\phi = {gamma}$ kg/ms, $u = {u}$ m/s")
    ax.set_xlabel("Grid-Point Count")
    ax.set_ylabel("Error (%)")
    ax.set_yscale('log')
    ax.grid(True, which="both", ls="-", alpha=0.2)
    ax.legend()

    return plt.show()