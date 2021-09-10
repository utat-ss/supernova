# Neural Network Training API
from ctypes import CDLL, POINTER, c_double
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
# Define functions from supernova C in Python


def get_state_and_acceleration(state: np.ndarray) -> \
        "tuple[np.ndarray, np.ndarray]":
    '''
    state: 6 keplerian elements presented as [SMA, ECC, INC, RAAN, AOP, MA]
    all distances in m
    all angles in radians
    '''
    rv = state_from_keplerian(*state)  # position/velocity vector
    result = (c_double*6)()    # place to store the resulting velocity/acceleration
    force_model(0, rv, result)

    return (np.ctypeslib.as_array(rv, shape=(6,)), np.ctypeslib.as_array(result, shape=(6,)))


if __name__ == "__main__":
    # Here's a demo using orbit from http://spacesys.utat.ca/confluence/display/FIN/Current+Orbital+Parameters
    position, effect = get_state_and_acceleration([6903e3, 0.004, np.radians(97.5), 0, 0, 0])
    print(f"Input: {position} == current position and velocity")
    print(f"Output: {effect} == resulting change in position and velocity")

    # Vary the keplerian elements bit by bit to train
    # See https://replit.com/@itchono/Orbit-Visualizer#main.py
    # to see how orbital elements affect the orbit


