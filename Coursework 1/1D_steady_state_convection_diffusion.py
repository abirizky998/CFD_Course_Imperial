import numpy as np
from algorithm_and_tools import analytical, plotting, error_calc
from schemes import *

def solve(length, cell_count, gamma_T, density, vel):
    # Domain setup
    delta_x = length / cell_count
    x = np.linspace(0, length, cell_count)

    # T field init 
    Temp = np.full(cell_count,20.0) # deg C

    # Equations Setup
    Pe_global = round((density * vel * length / gamma_T),3)
    Pe_local = round((density * vel * delta_x / gamma_T),3)
    T_analytical = analytical(x, length, Pe_global)

    schemes = ["Central Differencing Scheme",
                "Upwind Differencing Scheme",
                "Power Law Differencing Scheme",]
    results = [cds(cell_count, gamma_T, delta_x, density, vel, Temp), 
                uds(cell_count, gamma_T, delta_x, density, vel, Temp), 
                pds(cell_count, gamma_T, delta_x, density, vel, Temp, Pe_local)]
    errors = [error_calc(cell_count, T_analytical, results[0]),
              error_calc(cell_count, T_analytical, results[1]),
              error_calc(cell_count, T_analytical, results[2])]

    figures = []
    for i, j, k in zip(schemes, results, errors):
        figure = plotting(x, T_analytical, i, j, k, Pe_global, Pe_local)
        figures.append(figure)
        
    return (i for i in figures)

solve(length = 1.0, cell_count = 10, gamma_T = 0.1, density = 1.0, vel = 0.5)