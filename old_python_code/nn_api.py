# Neural Network Training API
from ctypes import CDLL, POINTER, c_double, c_void_p
import numpy as np

'''
Code Responsible for Loading C Library
'''
supernova = CDLL(("./python_code/supernova.so"))
# Load library
# Make sure you have a supernova.so binary compiled for your system


def _wrap_func(lib, funcname, restype, argtypes):
    ''' Referenced from https://dbader.org/blog/python-ctypes-tutorial-part-2
    '''
    func = lib.__getattr__(funcname)
    func.restype = restype
    func.argtypes = argtypes
    return func


force_model = _wrap_func(
    supernova, "simplified_perturbations", None,
    [c_double, POINTER(c_double), POINTER(c_double)])

state_from_keplerian = _wrap_func(
    supernova, "stateFromKeplerian", POINTER(c_double),
    [c_double, c_double, c_double, c_double, c_double, c_double])

next_state = _wrap_func(
    supernova, "nextStatefromCurrent", POINTER(c_double),
    [c_double, c_double, c_double, c_double, c_double, c_double, c_double]
)

free = _wrap_func(supernova, "freemem", None, [c_void_p])
# Define functions from supernova C in Python


def get_state_and_acceleration(state: np.ndarray) -> \
        "tuple[np.ndarray, np.ndarray]":
    '''
    state: 6 keplerian elements presented as [SMA, ECC, INC, RAAN, AOP, MA]
    all distances in m
    all angles in radians
    '''
    rv = state_from_keplerian(*state)  # position/velocity vector
    result = (c_double*6)()  # place to store the resulting velocity/acceleration
    force_model(0, rv, result)

    py_rv = np.copy(np.ctypeslib.as_array(rv, shape=(6,)))
    py_result = np.copy(np.ctypeslib.as_array(result, shape=(6,)))
    free(rv)

    return py_rv, py_result


def get_state_and_next_state(state: np.ndarray, timestep) -> \
        "np.ndarray":
    '''
    state: 6 position/velocity state presented as [X, Y, Z, VX, VY, VZ]
    timestep: how many seconds forward to advance orbit (keep it less than 60 and more than 1)
    all distances in m
    '''
    next = next_state(*state, timestep)  # position/velocity vector
    py_next = np.copy(np.ctypeslib.as_array(next, shape=(6,)))

    free(next)

    return py_next


if __name__ == "__main__":
    # Here's a demo using orbit from http://spacesys.utat.ca/confluence/display/FIN/Current+Orbital+Parameters
    position, effect = get_state_and_acceleration([6903e3, 0.004, np.radians(97.5), 0, 0, 0])
    print(f"Input: {position} == current position and velocity")
    print(f"Output: {effect} == resulting change in position and velocity")
    timestep = 10
    next = get_state_and_next_state(position, timestep)
    print(f"Next state in {timestep} seconds: {next} == next position and velocity")

    # Vary the keplerian elements bit by bit to train
    # See https://replit.com/@itchono/Orbit-Visualizer#main.py
    # to see how orbital elements affect the orbit


