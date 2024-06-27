import numpy as np


A = np.array([[1, 2, -2], [1, 1, 1], [2, 2, 1]])
B = np.array([[2, -1, 1], [1, 1, 1], [1, 1, -2]])
D = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
I = np.eye(3)
C = I - np.dot(np.linalg.inv(D), A)
print(C)
print(max(abs(np.linalg.eigvals(C))))
