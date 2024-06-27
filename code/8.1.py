'''
import numpy as np


A = np.array([[-4, 4, 1], [0, 1, 0], [0, 0, 1]])
X = np.array([[1], [1], [4]])
i = 1
temp = np.array([0, 0, 0])
epsilon = 0.0001
while True:
    print(i, "th iteration", sep="")
    i += 1
    Y = A @ (X/np.linalg.norm(X, ord=np.inf))
    print("Y=", Y)
    print("Y/X=", Y/X)
    if np.linalg.norm(Y / X - temp) < epsilon:
        break
    temp = Y / X
    X = Y

print(np.linalg.eigvals(A))
'''

import numpy as np

def power_method(A, tol=1e-9, max_iterations=1000):
    """
    使用幂法求矩阵A的按模最大的特征值及其对应的特征向量，并输出最后5次迭代的情况。

    参数:
    - A: numpy array，代表矩阵A。
    - tol: float, 迭代的容差。
    - max_iterations: int, 最大迭代次数。

    返回:
    - lambda_approx: 按模最大的特征值的近似值。
    - x_approx: 对应的特征向量的近似值。
    """
    n, _ = A.shape
    x = np.random.rand(n)  # 随机初始化一个向量
    x = x / np.linalg.norm(x)  # 归一化这个向量

    lambda_history = []
    x_history = []

    lambda_approx = 0
    for iteration in range(max_iterations):
        Ax = A @ x
        lambda_new = np.linalg.norm(Ax)
        x_new = Ax / lambda_new

        lambda_history.append(lambda_new)
        x_history.append(x_new)

        # 检查收敛性
        if np.abs(lambda_new - lambda_approx) < tol:
            print_last_5_iterations(lambda_history, x_history)
            return lambda_new, x_new

        x = x_new
        lambda_approx = lambda_new

    print("警告: 达到最大迭代次数。可能未收敛。")
    print_last_5_iterations(lambda_history, x_history)
    return lambda_approx, x

def print_last_5_iterations(lambdas, vectors):
    """
    输出最后5次迭代的特征值和特征向量。

    参数:
    - lambdas: 特征值的历史列表。
    - vectors: 特征向量的历史列表。
    """
    last_5 = max(0, len(lambdas) - 5)
    print("最后5次迭代情况：")
    for i in range(last_5, len(lambdas)):
        print(f"Iteration {i+1}: Lambda = {lambdas[i]}, Vector = {vectors[i]}")

# 示例矩阵
A = np.array([[-4, 4, 1], [0, 1, 0], [0, 0, 1]])

# 使用幂法求按模最大的特征值和特征向量
eigenvalue, eigenvector = power_method(A)

print("按模最大的特征值:", eigenvalue)
print("对应的特征向量:", eigenvector)
