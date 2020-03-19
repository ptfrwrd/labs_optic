"""
    TO DO: двумерное преобразование Фурье.
    [u,v] = [-4,4]
    sigma ** 2 = 1
    K(x, y, u, v) = exp( - 2 * pi * i * (x*u + y*v)
    f(x, y) = H_2(x) * H_3(x) * exp( (-x **2 - y **2)/ sigma ** 2)

"""


import cmath
import scipy.special as scp
import numpy as np
import matplotlib.pyplot as plt

T = 1

def rect(x, T):
    return abs(x) < T / 2.0


# полином Эрмита, n - степень
def hermite(n, x):
    H_n= scp.eval_hermite(n, x) * np.exp(-x**2)
    return H_n


# подсчёт функции
def integral_function(x, sigma, u, n_herm):
    exp_f = cmath.exp(- x ** 2 / sigma ** 2)
    exp_K = cmath.exp(- 1j * 2 * np.pi * x * u)
    #n_hermite = hermite(n_herm, x)
    return exp_f * exp_K * rect(x, T)


# численный подсчёт
def numerical_method(a, b, n, sigma, n_herm, u, v, m):
    points = []
    result_function = 0
    h_x = (b - a) / n
    h_e = (v - u) / m
    for j in range(m + 1):
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
    size = len(x_part)
    for i in range(size):
        result_conversion.append([])
        for j in range(size):
            result_conversion[i].append(x_part[i] * y_part[j])
    chart(result_conversion, u, v, size)
    return result_conversion


# амплитуда преобразования
def amplitude(conversion_points, size):
    amplitude_conv = []
    for i in range(size):
        amplitude_conv.append([])
        for j in range(size):
            amplitude_conv[i].append(np.sqrt(conversion_points[i][j].imag ** 2 +
                                      + conversion_points[i][j].real ** 2))
    return amplitude_conv



def phase(conversion_points, size):
    phase_con = []
    for i in range(size):
        phase_con.append([])
        for j in range(size):
            phase_con[i].append(np.angle(conversion_points[i][j]))
    return phase_con


# построение графика
def chart(conversion_points, a, b, size):

    x = np.linspace(float(a), float(b), 1000)
    y = np.linspace(float(a), float(b), 1000)
    amplitude_conv = amplitude(conversion_points, size)
    fig, ax = plt.subplots()
    ax.imshow(amplitude_conv, extent=(a, b, a, b))

    plt.savefig("амплитуда.png")
    plt.show()
    phase_con = phase(conversion_points, size)
    fig, ax = plt.subplots()
    ax.imshow(phase_con, extent=(a, b, a, b))

    plt.savefig("фаза.png")
    plt.show()


if __name__ == '__main__':
    a, b = -3, 3
    n, m = 1000, 1000
    sigma = 1
    n_herm_x, n_herm_y = 2, 3
    u, v = -4, 4
    res = result(a, b, u, v, n, m, sigma, n_herm_x, n_herm_y)
