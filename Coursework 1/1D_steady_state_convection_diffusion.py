import numpy as np
from algorithm_and_tools import analytical, plotting
from schemes import *

def solve():
    # Domain setup
    length = 1.0 # m
    cell_count = 10
    delta_x = length / cell_count
    x_grid = np.linspace(0, length, cell_count)

    # Fluid and Flow Properties setup
    gamma_T = 0.2 # diffusivity, between 0.1 and 1.0
    density = 1 # kg/m^3, value between 0.1 and 1.0
    Temp = np.full(cell_count,20.0) # deg C, T field initialisation
    vel = 1.0 # m/s

    # Equations Setup
    Pe_global = density * vel * length / gamma_T
    Pe_local = density * vel * delta_x / gamma_T
    
    print("Global Peclet Number:", Pe_global)
    print("Local Peclet Number:", Pe_local)
    
    T_cds = cds(cell_count, gamma_T, delta_x, density, vel, Temp)
    T_uds = uds(cell_count, gamma_T, delta_x, density, vel, Temp)
    T_plds = plds(cell_count, gamma_T, delta_x, density, vel, Temp, Pe_local)
    T_analytical = analytical(x_grid, length, Pe_global)

    result_lists = [T_cds, T_uds, T_plds]
    result_names = ["Central Differencing Scheme",
                    "Upwind Differencing Scheme",
                    "Power Law Differencing Scheme",]
    figures = []

    for i, j in zip(result_lists, result_names):
        figure = plotting(x_grid, i, T_analytical, j)
        figures.append(figure)

    return figures

solve()