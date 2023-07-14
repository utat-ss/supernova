# Simplified physics model for orbit propagation
import numpy as np

GM = 3.986004418e14  # m^3/s^
J2 = 1.082626925638742e-3  # m^2/s^2
RE = 6378.137e3  # m

def simple_gravity(t: float, y: np.array) -> np.array:
    '''
    Calculates gravitational forces on spacecraft
    using GM/r^2
    '''
    position = y[:3]

    r = np.linalg.norm(position)
    a_g = -GM/r**3 * position  # acceleration due to GM/r^2
    return np.concatenate((y[3:], a_g))

def J2_gravity(t: float, y: np.array) -> np.array:
    '''
    Calculates gravitational forces on spacecraft
    using GM/r^2 and J2 perturbation
    '''
    position = y[:3]

    r = np.linalg.norm(position)
    factor = (3.0/2.0) * GM * RE**2 * J2 / r**5
    z = position[2]

    a_g = -GM/r**3 * position  # acceleration due to GM/r^2
    a_J2 = factor * position * np.array([5*z**2/r**2 - 1, 5*z**2/r**2 - 1, 5*z**2/r**2 - 3])

    return np.concatenate((y[3:], a_J2 + a_g))




    