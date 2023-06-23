# RK4 Weights
import numpy as np

A = np.array([
    [0, 0, 0, 0],
    [1/2, 0, 0, 0],
    [0, 1/2, 0, 0],
    [0, 0, 1, 0]
    ])

C = np.array([0, 1/2, 1/2, 1])

B = np.array([1/6, 1/3, 1/3, 1/6])