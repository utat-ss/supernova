# Generates stability region of some ODE methods.
from matplotlib import pyplot as plt
import numpy as np
from ode_solver import RK1012_step, RK4_step


def test_f(t, y):
        return y


def RK4_stability_region(radius: float):
    xx, yy = np.meshgrid(np.linspace(-radius, radius, 200),
                         np.linspace(-radius, radius, 200))
    h = xx + 1j*yy  # step sizes to test with
    result = np.zeros(h.shape, dtype=complex)

    for i in range(len(h)):
        for j in range(len(h)):
            result[i, j] = RK4_step(test_f, np.array([1]), 0, h[i, j])[0]

    # Plot the stability region
    mask = np.abs(result) > 1
    mask = mask.astype(int)

    plt.contourf(xx, yy, mask, cmap='Blues')


def RK1012_stability_region(radius: float):
    '''
    Generates the stability region for the linear test problem
    y' = lambda * y, where lambda is a complex number.
    '''
    xx, yy = np.meshgrid(np.linspace(-radius, radius, 400),
                         np.linspace(-radius, radius, 400))
    h = xx + 1j*yy  # step sizes to test with
    result = np.zeros(h.shape, dtype=complex)

    for i in range(len(h)):
        for j in range(len(h)):
            result[i, j] = RK1012_step(test_f, np.array([1]), 0, h[i, j])[0][0]

    # Plot the stability region
    mask = np.abs(result) > 1
    mask = mask.astype(int)

    plt.contourf(xx, yy, mask, cmap='Blues')


if __name__ == "__main__":
    # SET dtype=complex on the solvers page to see the stability region

    plt.figure(figsize=(12, 6))

    plt.subplot(121)
    RK1012_stability_region(4)
    plt.title('RK1012')
    plt.xlabel("Re(L*dt)")
    plt.ylabel("Im(L*dt)")
    plt.grid()

    plt.subplot(122)
    RK4_stability_region(4)
    plt.title('RK4')
    plt.xlabel("Re(L*dt)")
    plt.ylabel("Im(L*dt)")
    plt.grid()

    plt.show()