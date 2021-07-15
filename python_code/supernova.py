from ctypes import CDLL, POINTER, c_double, c_char_p
import numpy as np
from poliastro.twobody import Orbit
from poliastro.bodies import Earth
from astropy import units as u
from astropy.time import Time
import time
import plotter
from reference_propagator import prop

'''
Loading C Library
'''
supernova = CDLL(("./python_code/supernova.so"))
# Load library


def wrap_func(lib, funcname, restype, argtypes):
    ''' Referenced from
    https://dbader.org/blog/python-ctypes-tutorial-part-2
    Taken from ESC190
    '''
    func = lib.__getattr__(funcname)
    func.restype = restype
    func.argtypes = argtypes
    return func


orbit = wrap_func(supernova, "orbit", None,
                  [c_char_p, POINTER(c_double), POINTER(c_double),
                   c_char_p, c_double])
# Define function from supernova C in Python


def get_orbit():
    '''
    Get SSO from poliastro
    '''
    SMA = 6878 * u.km
    LTAN = 22.5 * u.hourangle
    Eccentricity = 0 * u.one
    argument_of_perigee = 0 * u.deg
    total_anomaly = 0 * u.deg
    epoch_date = "2023-06-01"
    epoch_time = "00:00"

    # Store time as UTC
    epoch = Time(f"{epoch_date} {epoch_time}", format='iso', scale="utc")

    print("Determining intial orbit from parameters...")
    # Define orbital parameters, along with units
    orb = Orbit.heliosynchronous(attractor=Earth,
                                 a=SMA, ecc=Eccentricity,
                                 ltan=LTAN,
                                 argp=argument_of_perigee,
                                 nu=total_anomaly,
                                 epoch=epoch)

    print(orb)
    return orb


def compare_propagators(days: int):
    # Initial conditions
    py_tSpan = [0, 86400 * days]
    tSpan = (c_double * len(py_tSpan))(*py_tSpan)

    t_start = time.perf_counter()
    orb = get_orbit()
    t_orb = time.perf_counter()-t_start
    print(f"parameter fetch completion in {t_orb} seconds.")
    py_y0 = [*(orb.r.value*1000), *(orb.v.value*1000)]
    y0 = (c_double * len(py_y0))(*py_y0)

    t_start = time.perf_counter()
    orbit("RK810".encode("utf-8"), tSpan, y0, "FINCH.csv".encode("utf-8"), 1e-4)
    t_C = time.perf_counter()-t_start
    print(f"Supernova completion in {t_C} seconds.")

    arr = np.loadtxt("FINCH.csv", delimiter=",")  # Load times
    ref_times = arr[:, 0]

    t_start = time.perf_counter()
    prop(orb, ref_times, "poliastro.csv")
    t_Poli = time.perf_counter()-t_start
    print(f"Poliastro propagator completion in {t_Poli} seconds.")
    print(f"Supernova is {t_Poli/t_C}x faster than Poliastro")

    # Error analysis
    arrPoli = np.loadtxt("poliastro.csv", delimiter=",")
    diff = arr[:, 1:4] - arrPoli[:, 1:]
    print(f"RMS error between Supernova and Poliastro is {np.sqrt(np.average(diff**2))} metres across {len(ref_times)} steps.")
    print(f"Max error between Supernova and Poliastro is {np.max(np.abs(diff))} metres across {len(ref_times)} steps.")

    plotter.plot_trajectory("FINCH.csv", "poliastro.csv")
    plotter.plot_error(np.abs(diff))


def propagate_orbit(days: int):
    # Initial conditions
    py_tSpan = [0, 86400 * days]
    tSpan = (c_double * len(py_tSpan))(*py_tSpan)

    orb = get_orbit()
    py_y0 = [*(orb.r.value*1000), *(orb.v.value*1000)]
    y0 = (c_double * len(py_y0))(*py_y0)

    orbit("RK810".encode("utf-8"), tSpan, y0, "FINCH.csv".encode("utf-8"), 1e-12)

    plotter.plot_trajectory("FINCH.csv")


if __name__ == "__main__":
    compare_propagators(30)
