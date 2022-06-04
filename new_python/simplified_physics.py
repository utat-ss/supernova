# Simplified physics model for orbit propagation
import numpy as np

GM = 3.986004418e14  # m^3/s^
J2 = 1.082626925638742e-3  # m^2/s^2

def J2_gravity(t: float, y: np.array) -> np.array:
    '''
    Calculates gravitational forces on spacecraft
    using GM/r^2 and J2 perturbation
    '''
    position = y[:3]

    r = np.linalg.norm(position)
    r_3 = r**3

    a_g = -GM/r_3 * position  # acceleration due to GM/r^2
    a_J2 = -J2 * r_3 * np.cross(position, position)  # acceleration due to J2

    return np.concatenate((y[3:], a_J2 + a_g))




    