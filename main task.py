""" K(x,e) = exp(-alfa*x^2*e^2)H_4(x*e)
    [a,b]=[-3,3]
    [p,q]=[-3,3]
    Полином Эрмита
    H_n(x) = (-1)^n*exp^(x^2)*(d^n/dx^n e(-x^2))
    Входное ядро: f(x) = exp(i*betta*x)
    alfa = 1, m, n = 1000
    Рассчётная формула: F(e_l) = sum(k=0, n-1) (K(x_k,e_l)*f_k*h_x; l = 0,...,m
    m  - количество интервалов разбиения
    h_e = (q-p)/m , e_l = p + l*h_e, l =0,...,m
    e_0 = p, e_m = q
    x_k = a + k*h_x
    h_x = (b-a)/n

    1. реализовать численный рассчет
    2. построить график амплитуды и фазы оптического сигнала
    3. построить график амплитуды и фазы преобразования
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scp
from mpl_toolkits.mplot3d import Axes3D

# входное ядро
def f_k(a, b, n, k, betta):
    h_x = (b - a) / n
    x_k = a + k * h_x
    return cmath.exp(1j * betta * x_k)


# полином Эрмита, n - порядок
def hermite(n, x):
    H_n= scp.eval_hermite(n ,x)
    return H_n


# ядро из варианта
def kernel(e, alpha, x, n):
    Kernel = cmath.exp(-alpha * x ** 2 * e ** 2) * hermite(n, x*e)
    return Kernel


# преобразование
def conversion(p, q, m, n, alpha, a, b, betta, n_hermite):
    res = []
    h_e = (q - p) / m
    temp = 0
    for l in range(0, m):
        e_l = p + l * h_e
        temp = series(e_l, alpha, n_hermite, a, b,n, betta)
        res.append(temp)
    return res


def series(e_l, alpha, n_hermite, a, b,n, betta):
    temp = 0
    h_x = (b - a) / n
    x_k = a
    for k in range(0, n):
        x_k += h_x
        temp += kernel(e_l, alpha, x_k, n_hermite) * f_k(a, b, n, k, betta) * h_x
    return temp

# Графики для оптического сигнала
def signal_chart(betta, a, b):
    x_points, y_points = [],[]

    for i in range(-3, 3):
        x_points.append(i)
        y_points.append(1)

    plt.subplot(2, 1, 1)
    plt.plot(x_points, y_points, color='red')
    plt.title("Амплитуда оптического сигнала", fontsize=10)
    plt.xlabel("x", fontsize=10)
    plt.ylabel("A", fontsize=10)
    plt.grid(True)
    x_points.clear()

    plt.subplot(2, 1, 2)
    phase_points = []
    x_p = np.linspace(-3.0, 3.0, 1000)
    for i in x_p:
        temp = np.angle(np.exp(1j * i * betta))
        phase_points.append(temp)
    plt.plot(x_p, phase_points)

    plt.xlabel("x", fontsize=10)
    plt.ylabel("phase", fontsize=10)
    plt.title("Фаза оптического сигнала", fontsize=10)
    plt.grid(True)
    plt.show()


# подсчёт амплитуды для преобразования
def amplitude_conversion(mass):
    abs_mass = []
    for i in range(len(mass)):
        abs_mass.append(np.sqrt(mass[i].imag ** 2 + mass[i].real ** 2))
    print(abs_mass)
    return abs_mass


# графики для преобразования
def conversion_chart(mass, a, b, n):
    x_points = np.linspace(float(a), float(b), n)
    plt.subplot(2, 1, 1)
    plt.plot(x_points, amplitude_conversion(mass), color='red')
    plt.title("Амплитуда преобразования", fontsize=10)
    plt.xlabel("x", fontsize=10)
    plt.ylabel("A", fontsize=10)
    plt.grid(True)
    plt.subplot(2, 1, 2)
    pha = []
    for i in range(len(mass)):
        pha.append(np.angle(mass[i]))
    plt.plot(x_points, pha, color='blue')
    plt.title("Фаза преобразования", fontsize=10)
    plt.xlabel("x", fontsize=10)
    plt.ylabel("phase", fontsize=10)
    plt.grid(True)
    plt.show()


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
    plt.show()

if __name__ == '__main__':
    p, a, q, b = -3, -3, 3, 3
    m, n = 100, 100
    alfa = 1
    betta = 3
    mass = []
    x = np.linspace(-5, 5, 1000)
    h = []
    for i in x:
        h.append(hermite(4,i))
    plt.plot(x, h)
    plt.show()
    n_herm = 4
    mass = conversion(p, q, m, n, alfa, a, b, betta, n_herm)
    print(mass)
    chart_3D(mass, a, b)
    x = np.linspace(float(a), float(b), len(mass))
    signal_chart(betta, a, b)
    conversion_chart(mass, a, b, n)
