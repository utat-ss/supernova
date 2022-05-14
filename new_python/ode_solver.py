from typing import Callable
import numpy as np
from propagator_solution import PropagatorSolution
from weight_reformatter import A, B, C

# High order RK method for solving IVP


def RK_solve(f: Callable, y0: np.array, t_span: np.array,
             etol: float, h0: float = 50) -> PropagatorSolution:
    '''12th Order Adaptive Runge-Kutta method

    Uses embedded 10th order method to produce error estimate.

    Parameters
    ----------
    f: Callable
        function for which y' = f(t, y)
    y0: np.array
        initial state
    t_span: np.array
        start and end of propagation, Julian day
    etol: float
        error tolerance
    h0: float
        initial step size, in seconds

    Returns
    -------
    PropagatorSolution
        propagated orbit
    '''

    if t_span[0] > t_span[1]:
        raise ValueError('t_span must be increasing')

    # Initialize initial solution array
    n = max((t_span[1]-t_span[0])/h0, 10)

    y = np.zeros((n, len(y0)))
    y[0, :] = y0

    t_seconds = np.zeros(y.shape[0])
    t_jd = np.zeros(y.shape[0])
    t_jd[0] = t_span[0]

    step = 0
    h = h0

    while (t_jd[step] < t_span[1]):
        yi, err = RK1012_step(f, y[step, :], t_jd[step], h)

        # Check error tolerance
        if err < etol:
            # Accept step
            step += 1

            # Expand solution array if needed
            if step >= y.shape[0]:
                y = np.resize(y, (2*y.shape[0], y.shape[1]))
                t_seconds = np.resize(t_seconds, (2*t_seconds.shape[0]))
                t_jd = np.resize(t_jd, (2*t_jd.shape[0]))

            # Update solution array
            y[step, :] = yi
            t_seconds[step] = t_seconds[step-1] + h
            t_jd[step] = t_jd[step-1] + h/86400

        # Else, reject step
        # Adjust stepsize according to 12th order error estimate
        h = 0.9 * h * (etol/err)**(1/12)
        if h > h0 * 100:
            h = h0 * 90  # Failsafe for large step sizes

    # Trim solution array
    y = np.resize(y, (step+1, y.shape[1]))
    t_seconds = np.resize(t_seconds, (step+1))
    t_jd = np.resize(t_jd, (step+1))

    # Convert to datetime64
    t_datetime = np.datetime64(t_seconds * 1e9, 'ns')

    return PropagatorSolution(t_jd, y, t_datetime)


def RK1012_step(f, y, t, h) -> "tuple(np.array, float)":
    '''
    Makes a single step using the RK1012 method
    RK1012 is a 12th order method with embedded 10th order method
    consisting of 25 explicit stages.
    '''
    y_i = np.zeros((25, len(y)))  # Intermediate solutions
    k = np.zeros((25, len(y)))  # Intermediate derivatives

    # RK Stage 0
    y_i[0, :] = y
    k[0, :] = f(t, y)

    for stage in range(1, 25):
        # Generate Input y_i
        for w in range(stage):
            y_i[stage, :] += h * A[stage, w] * k[w, :]

        # Compute intermediate derivatives
        k[stage, :] = f(t + C[stage] * h, y_i[stage, :])

    # Compute error estimate
    # use only the first 3 elements of the vector (position error)
    # error estimate is given by (49/640) *  h  * (Fk[1]-Fk[23])
    err = np.linalg.norm((49/640) * h * (k[1, :] - k[23, :]))
    for w in range(24):
        y += h * B[w] * k[w, :]

    return y, err
