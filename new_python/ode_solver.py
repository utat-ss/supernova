from math import ceil
from typing import Callable
import numpy as np
from time import perf_counter

from propagator_solution import PropagatorSolution

from RK1210_weights import A as RK1012A, B as RK1012B, C as RK1012C
from RK4_weights import A as RK4A, B as RK4B, C as RK4C

# High order RK method for solving IVP


def RK1012_solve(f: Callable, y0: np.array, t_span: np.array,
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
    n = max(ceil((t_span[1]-t_span[0])*86400/h0), 10)

    y = np.zeros((n, len(y0)))
    y[0, :] = y0

    t_seconds = np.zeros(y.shape[0])
    t_jd = np.zeros(y.shape[0])
    t_jd[0] = t_span[0]

    step = 0
    h = h0

    time_start = perf_counter()

    while (t_jd[step] < t_span[1]):
        yi, err = RK1012_step(f, y[step, :], t_jd[step], h)
        
        # Check error tolerance to determine accept/reject step
        if err < etol:
            step += 1  # Accept step

            # Expand solution array if needed
            if step >= y.shape[0]:
                y = np.resize(y, (2*y.shape[0], y.shape[1]))
                t_seconds = np.resize(t_seconds, (2*t_seconds.shape[0]))
                t_jd = np.resize(t_jd, (2*t_jd.shape[0]))

            # Update solution array
            y[step, :] = yi
            t_seconds[step] = t_seconds[step-1] + h
            t_jd[step] = t_jd[step-1] + h/86400

        # Adjust stepsize according to 12th order error estimate
        # prevent division by zero
        h *= (0.9 * (etol/err)**(1/12) if err != 0 else 2.0)

    # Trim solution array
    y = np.resize(y, (step+1, y.shape[1]))
    t_seconds = np.resize(t_seconds, (step+1))
    t_jd = np.resize(t_jd, (step+1))

    elapsed_time = perf_counter() - time_start

    return PropagatorSolution(t_jd, y, elapsed_time)


def RK1012_step(f, y, t, h) -> "tuple(np.array, float)":
    '''
    Makes a single step using the RK1012 method
    RK1012 is a 12th order method with embedded 10th order method
    consisting of 25 explicit stages.
    '''
    k = np.zeros((25, len(y)))  # Intermediate derivatives

    for stage in range(25):
        k[stage, :] = f(t + h * RK1012C[stage],  # Derivatives
                        y + h * RK1012A[stage, :stage] @ k[:stage, :])
        # Reduce dimension of A and k to avoid wasting compute on zero elements

    # Compute error estimate
    err = abs(np.linalg.norm((49/640) * h * (k[1, :] - k[23, :])))
    
    return y + h * RK1012B @ k, err

def RK4_solve(f: Callable, y0: np.array, t_span: np.array,
             etol: float, h: float = 50):
    '''4th Order Classic Runge-Kutta method

    As standard as it gets.

    Parameters
    ----------
    f: Callable
        function for which y' = f(t, y)
    y0: np.array
        initial state
    t_span: np.array
        start and end of propagation, Julian day
    etol: float
        error tolerance, not used
    h0: float
        initial step size, in seconds

    Returns
    -------
    PropagatorSolution
        propagated orbit
    '''
    if t_span[0] > t_span[1]:
        raise ValueError('t_span must be increasing')

    n = ceil((t_span[1]-t_span[0])*86400/h)
    y = np.zeros((n, len(y0)))
    y[0, :] = y0

    t_seconds = np.zeros(y.shape[0])
    t_jd = np.zeros(y.shape[0])
    t_jd[0] = t_span[0]

    for step in range(n-1):
        y[step+1, :] = RK4_step(f, y[step, :], t_jd[step], h)
        t_seconds[step+1] = t_seconds[step] + h
        t_jd[step+1] = t_jd[step] + h/86400

    return (t_seconds, y)

def RK4_step(f, y, t, h) -> np.array:
    '''
    Makes a single step using the RK4 method
    '''
    k = np.zeros((4, len(y))) 
    for stage in range(4):
        k[stage, :] = f(t + h * RK4C[stage],
                        y + h * RK4A[stage, :stage] @ k[:stage, :])
    return y + h * RK4B @ k
 