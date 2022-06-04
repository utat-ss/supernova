# Test solvers
from ode_solver import RK1012_solve, RK1012_fixed_step
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp, OdeSolution
import numpy as np


def f_SHM(t, y):
    # Simple harmonic oscillator
    return [y[1], -y[0]]


def f_test(t, y):
    # Simple function whose solution is known
    # y = t**2 + 1 for y(0) = 1
    return 2*(t**3+t)/y[0]


def f_cos(t, y):
    # simple ODE for which the solution is cos(t)
    return np.cos(t) - y[0] - np.sin(t)


def test_SHM():
    y0 = [1, 0]
    t_span = [0, 10/86400]
    etol = 1e-4
    h0 = 0.1

    sol_1012 = RK1012_solve(f_SHM, y0, t_span, etol, h0)

    print(sol_1012)

    plt.plot(sol_1012.t, sol_1012.y[0, :], 'bo')
    plt.plot(sol_1012.t, sol_1012.y[1, :], 'ro')

    ref: OdeSolution = solve_ivp(f_SHM, (0, 10), y0, rtol=1e-8)

    plt.plot(ref.t, ref.y[0, :], 'b-')
    plt.plot(ref.t, ref.y[1, :], 'r-')

    plt.legend(("RK1012", "RK1012", "scipy RK45",
    "scipy RK45"), loc="best")

    plt.title("Test Case: Simple Harmonic Oscillator")
    plt.xlabel("Time (s)")
    plt.grid()
    plt.show()


def test_analytical():
    y0 = [1]
    t_span = [0, 10/86400]
    etol = 1e-4
    h0 = 0.1

    sol_1012 = RK1012_solve(f_test, y0, t_span, etol, h0)

    print(sol_1012)

    ref: OdeSolution = solve_ivp(f_test, (0, 10), y0, rtol=1e-8)

    plt.plot(sol_1012.t, sol_1012.y[0, :], 'bo')

    plt.plot(ref.t, ref.y[0, :], 'r-')

    analy_t = np.linspace(0, 10, 100)
    analy_y = analy_t**2 + 1

    plt.plot(analy_t, analy_y, 'b--')

    plt.legend(("RK1012", "scipy RK45", "analytical"), loc="best")

    plt.title("Test Case: Analytically Solvable ODE")
    plt.xlabel("Time (s)")
    plt.grid()
    plt.show()


def h_convergence():
    # Test error convergence for h
    hs = 10 / np.arange(4, 30)
    y0 = [1]
    t_span = [0, 100/86400]
    err = []

    for h in hs:
        sol_1012 = RK1012_fixed_step(f_cos, y0, t_span, 1e-4, h)
        err.append(abs(sol_1012.y[0, -1] - np.cos(100)))

    # fit log-log slope
    fit = np.polyfit(np.log(hs), np.log(err), 1)

    # log-log plot of error convergence
    plt.loglog(hs, err)
    plt.title("Error convergence for h")
    plt.xlabel("log(h)")
    plt.ylabel("log(error at final time)")
    plt.grid()
    plt.show()

def test_cos():
    # Test error convergence for h
    hs = [10/3, 2.5, 1, 0.1]
    y0 = [1]
    t_span = [0, 10/86400]

    for h in hs:
        sol_1012 = RK1012_fixed_step(f_cos, y0, t_span, 1e-4, h)
        plt.plot(sol_1012.t, sol_1012.y[0, :], label=f"h={h}")

    plt.xlabel("Time (s)")
    plt.ylabel("y(t)")
    plt.legend()
    plt.title("RK1012 Instability for cos(t) solution")
    plt.grid()
    plt.show()
    

if __name__ == "__main__":
    test_SHM()
    test_analytical()
    h_convergence()
    test_cos()
