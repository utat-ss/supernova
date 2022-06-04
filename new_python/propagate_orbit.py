import numpy as np
from typing import Callable
from propagator_solution import PropagatorSolution
from ode_solver import RK1012_solve
from simplified_physics import J2_gravity
from matplotlib import pyplot as plt

# import 3D axes
from mpl_toolkits.mplot3d import Axes3D

def propagate_orbit(t_span: np.array, y0: np.array,
                    force_model: Callable, solver: Callable,
                    h0: float = 50, etol: float = 1e-6) -> PropagatorSolution:
    '''
    Propagates orbit using a given solver.

    Parameters
    ----------
    t_span: np.array
        start and end of propagation, Julian day
    y0: np.array
        initial state; formatted as [x, y, z, vx, vy, vz] in m and m/s
        coordinate system: ECI J2000
    force_model: Callable
        function for which y' = f(t, y) (ODE)
    solver: Callable
        solver which performs integration
    h0: float
        initial step size, in seconds
    etol: float
        error tolerance for position (m)
    '''

    solution = solver(force_model, y0, t_span, h0, etol)

    return solution

if __name__ == '__main__':
    # Test code
    t_span = np.array([0, 0.1])  # propagate for one day
    y0 = np.array([6928e3, 0, 0, 0, 0.1e3, 7.5e3])  # initial state

    # Propagate using RK1012
    solution = propagate_orbit(t_span, y0, J2_gravity, RK1012_solve, etol=1e-12)
    print(solution)

    # Plot results
    plt.plot(solution.seconds_elapsed, solution.y[0, :], label='x')
    plt.plot(solution.seconds_elapsed, solution.y[1, :], label='y')
    plt.plot(solution.seconds_elapsed, solution.y[2, :], label='z')
    plt.legend()
    plt.show()

    # 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(solution.y[0, :], solution.y[1, :], solution.y[2, :], "o-")
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')

    plt.show()
