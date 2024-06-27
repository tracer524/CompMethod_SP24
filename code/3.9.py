def func(x):
    return (x**3-2)/3


x0 = 1
x1 = 3
temp = 0
epsilon = 0.001
while abs(func(x1)) > epsilon:
    temp = x1-func(x1)*(x1-x0)/(func(x1)-func(x0))
    x0 = x1
    x1 = temp
    print(temp, func(temp))
print(x1)
