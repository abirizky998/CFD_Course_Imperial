import numpy as np
import matplotlib.pyplot as plt

def tdma(N, a, b, c, d, T):
    P = np.zeros(N)
    Q = np.zeros(N)

    # west BC, T[0] is known so P[0] = 0 and Q[0] is T[0]
    P[0] = 0
    Q[0] = 100

    # setting up the inside elements
    for i in range(1, N-1): 
        P[i] = b[i] / (a[i] - c[i] * P[i-1])
        Q[i] = (d[i] + c[i] * Q[i-1]) / (a[i] - c[i] *P[i-1])

    # east BC, T[N] is known
    P[N-1] = 0
    Q[N-1] = 20

    for i in range(N-2, -1, -1): # backwards solve
        T[i] = P[i] * T[i+1] + Q[i]

    return T.copy()

def analytical(grid, L, Pe):
    return 100 + ((np.exp(grid * Pe / L) - 1) / (np.exp(Pe) - 1) * (20-100))

def error_calc(grid, analytical, numerical):
    E = 100 * round((sum(abs((numerical - analytical)/analytical)) / grid),4)
    return f'{E} %' 

def plotting(grid, analytical_results, title, scheme_results, error_value, Pe_G, Pe_L):
    fig, ax = plt.subplots()
    ax.plot(grid, analytical_results, label=f'Analytical', c='black')
    ax.plot(grid, scheme_results, label=f'Numerical, Error:{error_value}', 
            c='blue', ls='--', marker='o', mfc='red')

    ax.set_title(f"{title} \n Pe Global: {Pe_G}, Pe Local: {Pe_L}")
    ax.set_xlabel('x / L')
    ax.set_ylabel('Temperature')
    ax.legend()

    return fig