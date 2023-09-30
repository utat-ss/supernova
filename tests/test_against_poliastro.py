# Poliastro comparison tool to evaluate performance

import time

import numpy as np
from astropy import units as u
from astropy.time import Time, TimeDelta
from poliastro.bodies import Earth
from poliastro.core.perturbations import J2_perturbation
from poliastro.core.propagation import func_twobody
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import CowellPropagator
from poliastro.twobody.sampling import EpochsArray

import supernova.plotter as plotter
from supernova.api import propagate_orbit


# Force model
def f(t0, u_, k):
    du_kep = func_twobody(t0, u_, k)
    ax, ay, az = J2_perturbation(t0, u_, k, J2=Earth.J2.value, R=Earth.R.to(u.km).value)
    du_ad = np.array([0, 0, 0, ax, ay, az])
    return du_kep + du_ad


def prop(orb: Orbit, times):
    rr, vv = orb.to_ephem(
        EpochsArray(
            orb.epoch + TimeDelta(np.array(times) * u.s),
            method=CowellPropagator(f=f),
        )
    ).rv()

    return times, np.hstack((rr.to(u.m).value, vv.to(u.m / u.s).value))


def get_orbit() -> Orbit:
    """
    Get SSO from poliastro
    """
    SMA = 6903 * u.km
    LTAN = 22.5 * u.hourangle
    Eccentricity = 0.003621614 * u.one
    argument_of_perigee = 0 * u.deg
    total_anomaly = 0 * u.deg
    epoch_date = "2023-06-01"
    epoch_time = "00:00"

    # Store time as UTC
    epoch = Time(f"{epoch_date} {epoch_time}", format="iso", scale="utc")

    print("Determining intial orbit from parameters...")
    # Define orbital parameters, along with units
    orb = Orbit.heliosynchronous(
        attractor=Earth,
        a=SMA,
        ecc=Eccentricity,
        raan=LTAN,
        argp=argument_of_perigee,
        nu=total_anomaly,
        epoch=epoch,
    )

    print(orb)
    return orb


def compare_propagators(orbit: Orbit, t_seconds: float = 86400 * 1):
    # Initial conditions
    t_span = [0, t_seconds]
    supernova_y0 = [*(orbit.r.to(u.m).value), *(orbit.v.to(u.m / u.s).value)]

    t_start = time.perf_counter()
    t_supernova, y_supernova = propagate_orbit(
        "RK1012", "simplified", t_span, supernova_y0, 1e-6
    )
    t_C = time.perf_counter() - t_start
    print(f"Supernova completion in {t_C} seconds.")

    t_start = time.perf_counter()
    t_pa, y_pa = prop(orbit, t_supernova)
    t_Poli = time.perf_counter() - t_start
    print(f"Poliastro propagator completion in {t_Poli} seconds.")
    print(f"Supernova is {t_Poli/t_C}x faster than Poliastro")

    # Error analysis
    diff = y_supernova - y_pa
    print(
        f"RMS error between Supernova and Poliastro is {np.sqrt(np.average(diff**2))} metres across {len(t_supernova)} steps."
    )
    print(
        f"Max error between Supernova and Poliastro is {np.max(np.abs(diff))} metres across {len(t_supernova)} steps."
    )
    plotter.plot_error(np.abs(diff))


if __name__ == "__main__":
    orbit = get_orbit()
    # REMEMBER TO RECOMPILE WITH SUPERNOVA DRAG TURNED OFF
    compare_propagators(orbit, 86400 * 50)
