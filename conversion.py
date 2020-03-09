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


def f_k(a, b, n, k, betta):
    h_x = (b - a) / n
    x_k = a + k * h_x
    return (1j*betta*x_k)


def hermite(x):
    H_4 = cmath.exp(x ** 2) * diff(x)
    return H_4


def diff(x):
    exp = (4 * x ** 2 - 2) * cmath.exp(-x ** 2)
    return exp


def kernel(e, alfa, x):
    Kernel = cmath.exp(-alfa * x ** 2 * e ** 2) * hermite(x*e)
    return Kernel


def conversion(p, q, m, n, alfa, a, b, betta):
    res = []
    temp = 0
    h_e = (q - p) / m
    for l in range(0, n):
        e_l = p + l * h_e
        for k in range(0, n-1):
            h_x = (b - a) / n
            x_k = a + k * h_x
            temp += kernel(e_l, alfa, x_k) * f_k(a, b, n, k, betta) * h_x
        res.append(temp)
    return res


def interval():
    mass = []
    for i in np.arange(-cmath.pi, cmath.pi):
        mass.append(i)
    return mass


def phase(betta, x_points):
    return betta * x_points


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
    phase_points = []
    x_points.clear()

    # [-pi,pi]
    for x in np.arange(- np.pi,  2 * np.pi, np.pi):
        x_points.append(x)
        phase_points.append(phase(betta, x))

    # [-3pi, pi]
    x_points_left = []
    phase_points_left = []
    for x in np.arange(- 3 * np.pi,  0, np.pi):
        x_points_left.append(x)
        phase_points_left.append(phase(betta, x + 2*np.pi))


    # [pi, 3pi]
    x_points_right = []
    phase_points_right = []
    for x in np.arange(np.pi, 4 * np.pi, np.pi):
        x_points_right.append(x)
        phase_points_right.append(phase(betta, x - 2 * np.pi))

    plt.subplot(2, 1, 2)
    plt.plot(x_points, phase_points, color='blue')
    plt.plot(x_points_left, phase_points_left, color = 'blue')
    plt.plot(x_points_right, phase_points_right, color='blue')
    plt.xlabel("x", fontsize=10)
    plt.ylabel("phase", fontsize=10)
    plt.title("Фаза оптического сигнала", fontsize=10)
    plt.grid(True)
    plt.show()


def amplitude_conversion(mass):
    abs_mass = []
    for i in range(len(mass)):
        abs_mass.append(np.sqrt(mass[i].imag ** 2 + mass[i].real ** 2))
    return abs_mass

def phase_conversion(mass):
    phase_points = []
    for i in range(len(mass)):
        if mass[i].imag < 0:
            phase_points.append(3 * np.pi / 2)
        if mass[i].imag > 0:
            phase_points.append(np.pi / 2)
    return phase_points


def conversion_chart(amplitide, phase, a, b, n):
    x_points = np.linspace(float(a), float(b), n)
    plt.subplot(2, 1, 1)
    plt.plot(x_points, amplitide, color='red')
    plt.title("Амплитуда преобразования", fontsize=10)
    plt.xlabel("x", fontsize=10)
    plt.ylabel("A", fontsize=10)
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(x_points, phase, color='blue')
    plt.title("Фаза преобразования", fontsize=10)
    plt.xlabel("x", fontsize=10)
    plt.ylabel("phase", fontsize=10)
    plt.grid(True)
    plt.show()


if __name__ == '__main__':
    p, a, q, b = -3, -3, 3, 3
    m, n = 1000, 1000
    alfa = 1
    betta = 3
    mass = []
    mass = conversion(p, q, m, n, alfa, a, b, betta)
    x = np.linspace(float(a), float(b), len(mass))

    con_amplitude = amplitude_conversion(mass)
    print(con_amplitude)
    con_phase = phase_conversion(mass)
    for i in range(len(mass)):
        if mass[i].imag < 0:
            print (mass[i], "  ", con_phase[i])
        if mass[i].imag > 0:
            print(mass[i], " ", con_phase[i])

    signal_chart(betta, a , b)
    # график преобразования
    img=[]
    for i in range(n):
        img.append(mass[i].imag)
    plt.plot(x, img)
    plt.show()

    conversion_chart(con_amplitude, con_phase, a, b, n)