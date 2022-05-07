import numpy as np
from scipy.interpolate import interp1d


# Container class for a propagated orbit

class PropagatorSolution:
    def __init__(self, t: np.array, y: np.ndarray):
        self.r = y[:, 0:3]  # Position, in m
        self.v = y[:, 3:6]  # Velocity, in m/s
        self.jd = t  # Julian day
        self.elapsed_seconds = (self.jd - self.jd[0]) * 86400  # s
        # Convert jd to datetime
        self.t = np.datetime64((self.jd - 2440587.5) * 86400, 's')  # datetime

    def sample(self, t_sample: np.array) -> "tuple(np.array, np.array)":
        '''
        Samples the position and velocity at time range t_sample.
        '''
        # Interpolate position and velocity
        f_r = interp1d(self.jd, self.r, axis=0)
        f_v = interp1d(self.jd, self.v, axis=0)
        r = f_r(t_sample)
        v = f_v(t_sample)

        return r, v

    @property
    def epoch(self):
        '''
        Gives initial time of solution in julian day.
        '''
        return self.jd[0]
    
