from ctypes import CDLL, POINTER, c_double, c_int, c_char_p, Structure
from comparison import compare_propagators, get_orbit
import numpy as np

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
                  [c_char_p, POINTER(c_double), POINTER(c_double),
                   c_double])
# Define function from supernova C in Python


if __name__ == "__main__":
    # compare_propagators(orbit, 10)
    days = 10
    # orb = get_orbit()

    py_tSpan = [0, 86400 * days]
    tSpan = (c_double * len(py_tSpan))(*py_tSpan)

    # py_y0 = [*(orb.r.value*1000), *(orb.v.value*1000)]
    py_y0 = [-4749231.102296294, -4975106.82469687, 0.0, -719.6538503589323, 686.980716081442, 7561.282735263496]

    y0 = (c_double * len(py_y0))(*py_y0)

    solution: POINTER(AdaptiveSolution) = orbit("RK810".encode("utf-8"), tSpan, y0, 1e-12)

    n = solution[0].n  # Dereference pointer
    print(f"Steps taken: {n}")
    t = np.ctypeslib.as_array(solution[0].t, shape=(n,))
    y = np.ctypeslib.as_array(solution[0].y, shape=(n, 6))
    print(t)
    print(y)
