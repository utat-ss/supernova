from __future__ import annotations

import numpy as np
import numpy.typing as npt

from supernova.backend import ffi
from supernova.backend.lib import orbitPTR, stateFromKeplerian, free


def propagate_orbit(
    solver: str, model: str, t_span: tuple[float], y0: tuple[float], atol: float
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Propagate orbit using given solver and model.

    Parameters
    ----------
    solver : str
        Name of the solver to use, either "RK810" or "RK1012"
        RK810 is a 10th order, 17 stage Runge-Kutta method
        RK1012 is a 12th order, 25 stage Runge-Kutta method
    model : str
        Name of the model to use, either "simplified" or "advanced"
        simplified: only J2 perturbation is considered
        advanced: JGM-3 20th order spherical harmonic model is used
    t_span : tuple[float]
        Time span to propagate the orbit over, in seconds
    y0 : tuple[float]
        Initial state vector, in meters and meters per second (ECI)
        (x, y, z, vx, vy, vz)
    atol : float
        Absolute tolerance for the solver

    Returns
    -------
    t : npt.NDArray[np.float64]
        Time steps, in seconds
    y : npt.NDArray[np.float64]
        State vector, in meters and meters per second (ECI)
    """
    solution = orbitPTR(solver.encode("utf-8"), model.encode("utf-8"), t_span, y0, atol)

    t_buf = ffi.buffer(solution.t, solution.n * 8)
    t = np.frombuffer(t_buf, dtype=np.float64).copy()
    t.shape = (solution.n,)

    y_buf = ffi.buffer(solution.y, solution.n * 8 * 6)
    y = np.frombuffer(y_buf, dtype=np.float64).copy()
    y.shape = (solution.n, 6)

    free(solution.t)
    free(solution.y)

    return t, y


def state_from_keplerian(
    sma: float, ecc: float, inc: float, raan: float, aop: float, m: float
) -> tuple[float]:
    """Calculate state vector from Keplerian elements.

    Parameters
    ----------
    sma : float
        Semi-major axis, in meters
    ecc : float
        Eccentricity
    inc : float
        Inclination, in radians
    raan : float
        Right ascension of the ascending node, in radians
    aop : float
        Argument of periapsis, in radians
    m : float
        Mean anomaly, in radians
    """
    d_ptr = stateFromKeplerian(sma, ecc, inc, raan, aop, m)

    d_buf = ffi.buffer(d_ptr, 6 * 8)
    d = np.frombuffer(d_buf, dtype=np.float64).copy()

    free(d_ptr)
    return tuple(*d)
