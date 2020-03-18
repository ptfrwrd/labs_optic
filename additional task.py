"""
    TO DO: двумерное преобразование Фурье.
    [u,v] = [-4,4]
    sigma ** 2 = 1
    K(x, y, u, v) = exp( - 2 * pi * i * (x*u - y*v)
    f(x, y) = H_2(x) * H_3(x) * exp( (-x **2 - y **2)/ sigma ** 2)

"""


import cmath
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scp
from mpl_toolkits.mplot3d import Axes3D


# полином Эрмита, n - степень
def hermite(n, x):
    H_n= scp.eval_hermite(n, x)
    return H_n


# подсчёт функции
def integral_function(x, sigma, u, n_herm):
    exp_f = cmath.exp(- x ** 2 / sigma ** 2)
    exp_K = cmath.exp(- 1j * 2 * np.pi * x * u)
    n_hermite = hermite(n_herm, x)
    return exp_f * exp_K * n_hermite


# численный подсчёт
def numerical_method(a, b, n, sigma, n_herm, u, v, m):
    points = []
    result_function = 0
    h_x = (b - a) / n
    h_e = (v - u) / m
    for j in range(m):
        u_l = u + j * h_e
        for k in range(n):
            x_k = a + k * h_x
            result_function += integral_function(x_k, sigma, u_l, n_herm) * h_x
        points.append(result_function)
    return points


# итоговый результат
def result(a, b, u, v, n, m, sigma, n_herm_x, n_herm_y):
    x_part = numerical_method(a, b, n, sigma, n_herm_x, u, v, m)
    y_part = numerical_method(a, b, n, sigma, n_herm_y, u, v, m)

    result_conversion = []
    for i in range(len(x_part)):
        result_conversion.append(x_part[i]*y_part[i])
    chart(result_conversion, a, b)
    return result_conversion

# амплитуда преобразования
def amplitude(conversion_points):
    amplitude_conv = []
    for i in range(len(conversion_points)):
        amplitude_conv.append(np.sqrt(conversion_points[i].imag ** 2 +
                                      + conversion_points[i].real ** 2))
    return amplitude_conv

# 3D график
def chart_3D(mass, a, b):
    x = np.linspace(float(a), float(b), len(mass))
    y = []
    z = []
    for i in range(len(mass)):
        y.append(mass[i].real)
        z.append(mass[i].imag)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, label='conversion')
    plt.savefig('график преобразования')
    plt.show()


# построение графика
def chart(conversion_points, a, b):
    x_points = np.linspace(float(a), float(b), 100)
    y_points = amplitude(conversion_points)
    plt.subplot(2, 1, 1)
    plt.plot(x_points, y_points, color='red')
    plt.title("Амплитуда оптического сигнала", fontsize=10)
    plt.xlabel("x", fontsize=10)
    plt.ylabel("A", fontsize=10)
    plt.grid(True)
    y_points.clear()
    plt.subplot(2, 1, 2)
    for i in range(len(conversion_points)):
        temp = np.angle(conversion_points[i])
        y_points.append(temp)
    plt.plot(x_points, y_points,  color='blue')
    plt.xlabel("x", fontsize=10)
    plt.ylabel("phase", fontsize=10)
    plt.title("Фаза оптического сигнала", fontsize=10)
    plt.grid(True)
    plt.savefig('Амплитуда и фаза')
    plt.show()


if __name__ == '__main__':
    a, b = -3, 3
    n, m = 100, 100
    sigma = 1
    n_herm_x, n_herm_y = 2, 3
    u, v = -4, 4
    res = result(a, b, u, v, n, m, sigma, n_herm_x, n_herm_y)
    chart_3D(res, a, b)
    print(res)