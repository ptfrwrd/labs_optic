''' K(x,e) = exp(-alfa*x^2*e^2)H_4(x*e)
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
'''


import cmath
import numpy as np


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

p, a, q, b = -3, -3, 3, 3
m, n = 1000, 1000
alfa = 1
betta = 10
mass = []
mass = conversion(p, q, m, n, alfa, a, b, betta)
print(mass)