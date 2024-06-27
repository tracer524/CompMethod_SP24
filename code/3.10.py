import numpy as np
from scipy import linalg


def jacobi(y: np.array([(), ()])):
    return np.array([(2*y[0], 2*y[1]), (3*y[0]**2, -1)])


def func(y: np.array([(), ()])):
    return np.array([y[0]**2+y[1]**2-1, y[0]**3-y[1]]).T


x = np.array([0.8, 0.6]).T
delta = np.array([1, 1]).T
while delta.max() > 0.001:
    delta = linalg.inv(jacobi(x))@func(x)
    x = x - delta
    print(x)
