import numpy as np
from scipy.interpolate import interp1d


# Container class for a propagated orbit

class PropagatorSolution:
    def __init__(self, t: np.array, y: np.ndarray, solution_time: float):
        self.y = y  # solution array
        self.jd = t  # Julian day
        self.seconds_elapsed = (self.jd - self.jd[0]) * 86400  # s
        self.datetime = ((self.jd - 2440587.5) * 86400).astype("datetime64[s]")  # datetime
        self.solution_time = solution_time

    def sample(self, t_sample: np.array) -> np.array:
        '''
        Samples the position and velocity at time range t_sample.
        '''
        # Interpolate position and velocity
        f_y = interp1d(self.jd, self.y, axis=0)
        y = f_y(t_sample)
        return y

    @property
    def epoch(self):
        '''
        Gives initial time of solution in julian day.
        '''
        return self.jd[0]

    @property
    def t(self):
        '''
        Returns time in format similar to scipy solution
        '''
        return self.seconds_elapsed

    def __len__(self) -> int:
        '''
        Length
        '''
        return len(self.t)
    
    def __repr__(self) -> str:
        return f"PropagatorSolution with {len(self)} steps taken " \
               f"in {self.solution_time:.5f} seconds "\
               f"({self.solution_time/len(self):.5f} sec per step)"
