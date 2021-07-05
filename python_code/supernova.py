from ctypes import CDLL, POINTER, c_double, c_int, c_char_p
import time
from poliastro.twobody import Orbit
from poliastro.bodies import Earth
from astropy import units as u
from astropy.time import Time
import plotter

'''
Loading C Library
'''
supernova = CDLL(("./python_code/supernova.so"))
# Load library


def wrap_func(lib, funcname, restype, argtypes):
    ''' Referenced from
    https://dbader.org/blog/python-ctypes-tutorial-part-2
    '''
    func = lib.__getattr__(funcname)
    func.restype = restype
    func.argtypes = argtypes
    return func


orbit = wrap_func(supernova, "orbit", None,
                  [c_int, POINTER(c_double), c_char_p, c_char_p])
# Define function from supernova C in Python

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
'''
Propagate using supernova
'''
# Initial conditions
py_y0 = [*(orb.r.value*1000), *(orb.v.value*1000)]
print("inputs rr vv:", py_y0)
y0 = (c_double * len(py_y0))(*py_y0)

t = time.perf_counter()
orbit(100000, y0,
      "./c_code/tables/RK8.txt".encode("utf-8"), "./results/RK8.csv".encode("utf-8"))
print(f"Propagation Completed in time: {time.perf_counter() - t}s")

plotter.plot_trajectory("./results/RK8.csv")
