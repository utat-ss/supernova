# Test RK1012 solver
from ode_solver import RK_solve
from matplotlib import pyplot as plt


# Simple harmonic oscillator
def f(t, y):
    return [y[1], -y[0]]


if __name__ == "__main__":
    y0 = [1, 0]
    t_span = [0, 10]
    etol = 1e-6
    h0 = 50

    sol = RK_solve(f, y0, t_span, etol, h0)

    plt.plot(sol.t_jd, sol.y[:, 0], 'b-')
    plt.plot(sol.t_jd, sol.y[:, 1], 'r-')
    plt.show()
