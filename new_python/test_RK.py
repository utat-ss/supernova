# Test solvers
from ode_solver import RK1012_solve
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp, OdeSolution


# Simple harmonic oscillator
def f(t, y):
    return [y[1], -y[0]]


if __name__ == "__main__":
    y0 = [1, 0]
    t_span = [0, 10/86400]
    etol = 1e-4
    h0 = 0.1

    sol_1012 = RK1012_solve(f, y0, t_span, etol, h0)

    print(sol_1012)

    plt.plot(sol_1012.t, sol_1012.y[:, 0], 'bo')
    plt.plot(sol_1012.t, sol_1012.y[:, 1], 'ro')

    ref: OdeSolution = solve_ivp(f, (0, 10), y0, rtol=1e-8)

    plt.plot(ref.t, ref.y[0, :], 'b-')
    plt.plot(ref.t, ref.y[1, :], 'r-')

    plt.legend(("RK1012", "RK1012", "scipy.integrate.solve_ivp",
    "scipy.integrate.solve_ivp"), loc="best")

    plt.title("Test Case: Simple Harmonic Oscillator")
    plt.xlabel("Time (s)")
    plt.grid()
    plt.show()
