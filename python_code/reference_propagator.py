from numba import njit
from poliastro.twobody.propagation import func_twobody, propagate, cowell
from poliastro.core.perturbations import J2_perturbation
from astropy import units as u
from poliastro.bodies import Earth
import numpy as np


#@title Aerodynamic Characteristics

#@markdown C_D: Drag coefficeint -- dim ensionless

#@markdown FRONTAL_AREA: Frontal area of S/C in m^2

#@markdown MASS: Spacecraft mass in kg

@njit
def a_d(t0, state, k, J2, R):
    return J2_perturbation(t0, state, k, J2, R)

# Force model
def f(t0, state, k):
    du_kep = func_twobody(t0, state, k)
    # standard twobody

    ax, ay, az = a_d(
        t0,
        state,
        k,
        R=Earth.R.to(u.km).value,
        J2=Earth.J2.value,
    )
    du_ad = np.array([0, 0, 0, ax, ay, az])

    return du_kep + du_ad


def prop(orb, times, fname):
    rr = propagate(orb, times * u.s, method=cowell, f=f)
    arr = rr.xyz.value.T * 1000
    times = np.reshape(times, (times.shape[0], 1))
    np.savetxt(fname, np.concatenate((times, arr), axis=-1), delimiter=",")
