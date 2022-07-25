from ctypes import CDLL, POINTER, c_double, c_int, c_char_p, Structure
import numpy as np
from plotter import plot_from_array, plot_3d_from_array
from celest.satellite import Time, Coordinate, Satellite
from celest.encounter import GroundPosition, windows

'''
Code Responsible for Loading C Library
'''
supernova = CDLL(("./python_code/supernova.so"))
# Load library


class AdaptiveSolution(Structure):
    _fields_ = [("t", POINTER(c_double)),
                ("y", POINTER(c_double)),
                ("n", c_int)]
    '''
    Structure:
    t: Times array of shape (n,)
    y: position/velocity array of shape (n, 6) with units [m, s]
    n: number of steps used
    '''


def wrap_func(lib, funcname, restype, argtypes):
    ''' Referenced from
    https://dbader.org/blog/python-ctypes-tutorial-part-2
    Taken from ESC190
    '''
    func = lib.__getattr__(funcname)
    func.restype = restype
    func.argtypes = argtypes
    return func


orbit = wrap_func(supernova, "orbitPTR", POINTER(AdaptiveSolution),
                  [c_char_p, c_char_p, POINTER(c_double), POINTER(c_double),
                   c_double])
state = wrap_func(supernova, "stateFromKeplerian", POINTER(c_double),
                  [c_double, c_double, c_double, c_double, c_double, c_double])
# Define functions from supernova C in Python


if __name__ == "__main__":
    days = 50

    py_tSpan = [0, 86400 * days]
    tSpan = (c_double * len(py_tSpan))(*py_tSpan)

    state: POINTER(c_double) = state(6378e3+270e3, 0, 97 * np.pi/180, 0, 0, 0)

    py_y0 = [-4749231.102296294, -4975106.82469687, 0.0, -
             719.6538503589323, 686.980716081442, 7561.282735263496]
    y0 = (c_double * len(py_y0))(*py_y0)
    solution: POINTER(AdaptiveSolution) = orbit("RK810".encode(
        "utf-8"), "simplified".encode("utf-8"), tSpan, y0, 1e-6)

    n = solution[0].n  # Dereference pointer
    print(f"Steps taken: {n}")
    t = np.ctypeslib.as_array(solution[0].t, shape=(n,))
    y = np.ctypeslib.as_array(solution[0].y, shape=(n, 6))
    initial_state = np.ctypeslib.as_array(state, shape=(6,))
    plot_from_array(t, y)
    plot_3d_from_array(t, y)

    # Celest stuff
    julian = Time(t/86400)
    position = Coordinate(y[:, :3], "gcrs", t/86400)

    toronto = GroundPosition(latitude=43.6532, longitude=-79.3832)
    southern_manitoba = GroundPosition(latitude=53.7167, longitude=-98.8167)

    # Get ITRS data.
    itrs_x, itrs_y, itrs_z = position.itrs()

    # plot_from_array(t, np.vstack((itrs_x, itrs_y, itrs_z)).T)

    alt, az = position.horizontal(southern_manitoba)

    satellite = Satellite(position=y[:, :3], frame='gcrs', julian=t/86400,
                          offset=2414900)

    satellite.itrs()

    # Generate ground location windows.
    IMG_windows = windows.generate(
        satellite=satellite, location=southern_manitoba, enc="image", ang=30, lighting=1)
    GL_windows = windows.generate(
        satellite=satellite, location=toronto, enc="data_link", ang=10, lighting=0)

    print("hello")
