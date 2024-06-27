import numpy as np


def inverse_power_method(A, tol=1e-9, max_iterations=1000):
    n, _ = A.shape
    x = np.random.rand(n)
    x = x / np.linalg.norm(x)

    lambda_history = []
    x_history = []

    for _ in range(max_iterations):
        # 解线性系统 Ax = lambda*x，求解x
        x_new = np.linalg.solve(A, x)
        # 重新归一化向量
        x_new = x_new / np.linalg.norm(x_new)

        # 计算Rayleigh商，得到当前迭代的特征值近似
        lambda_new = np.dot(x_new, A @ x_new)

        lambda_history.append(1 / lambda_new)
        x_history.append(x_new)

        # 检查收敛性
        if np.linalg.norm(A @ x_new - lambda_new * x_new) < tol:
            print_last_iterations(lambda_history, x_history, 3)
            return 1 / lambda_new, x_new

        x = x_new

    print("警告: 达到最大迭代次数。可能未收敛。")
    print_last_iterations(lambda_history, x_history, 3)
    return 1 / lambda_new, x_new


def print_last_iterations(lambdas, vectors, num_iterations):
    start_index = max(0, len(lambdas) - num_iterations)
    print(f"最后{num_iterations}次迭代的情况：")
    for i in range(start_index, len(lambdas)):
        print(f"Iteration {i + 1}: Lambda = {lambdas[i]}, Vector = {vectors[i]}")


# 示例矩阵
A = np.array([[5, 3], [-2, 0]])

# 使用反幂法求按模最小的特征值和特征向量
min_eigenvalue, min_eigenvector = inverse_power_method(A)

print("按模最小的特征值:", min_eigenvalue)
print("对应的特征向量:", min_eigenvector)
