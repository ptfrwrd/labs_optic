"""
    TO DO: двумерное преобразование Фурье.
    [u,v] = [-4,4]
    sigma ** 2 = 1
    K(x, y, u, v) = exp( - 2 * pi * i * (x*u + y*v)
    f(x, y) = rect(x) * rect(y) * exp( (-x **2 - y **2)/ sigma ** 2)
"""
import numpy as np
import matplotlib.pyplot as plt

sigma = 1
T = 1


def rect(x):
    return abs(x) < T / 2.0


# temp =  x or y
def f_temp(temp):
    return rect(temp) * np.exp(- temp ** 2 / sigma ** 2)


# temp = x or y, temp_space = u or v
def kernel(temp, temp_space):
    return np.exp(-1j * np.pi * 2 * temp * temp_space)


# 1D integral; [a, b], n - upper series limit, temp_space - u_l or v_l
def finite_integral_1D(a, b, n, temp_space):
    h_temp = (b - a) / n
    temp_k, temp_series = a, 0
    for k in range(n):
        temp_k += h_temp
        temp_series += h_temp * kernel(temp_k, temp_space) * f_temp(temp_k)
    return temp_series


# 1D integral at point: point - temp_space
def finite_integral_at_point(a, b, n, m, u, v):
    integral = []
    h_space = (v - u) / m
    for l in range(m):
        temp_space = u + l * h_space
        integral.append(finite_integral_1D(a, b, n, temp_space))
    return integral


# 2D integral
def finite_integral_2D(a, b, n, m, u, v):
    dx_integral = finite_integral_at_point(a, b, n, m, u, v)
    dy_integral = finite_integral_at_point(a, b, n, m, u, v)
    dx_dy_integral = []
    for i in range(n):
        dx_dy_integral.append([])
        for j in range(m):
            dx_dy_integral[i].append(dx_integral[i] * dy_integral[j])
    return dx_dy_integral


# амплитуда преобразования
def amplitude(conversion_points, size):
    amplitude_conv = []
    for i in range(size):
        amplitude_conv.append([])
        for j in range(size):
            amplitude_conv[i].append(np.sqrt(conversion_points[j][i].imag ** 2 +
                                      + conversion_points[j][i].real ** 2))
    return amplitude_conv



def phase(conversion_points, size):
    phase_con = []
    for i in range(size):
        phase_con.append([])
        for j in range(size):
            phase_con[i].append(np.angle(conversion_points[i][j]))
    return phase_con


# построение графика
def chart(conversion_points, u, v, size):

    x = np.linspace(float(u), float(v), 1000)
    y = np.linspace(float(u), float(v), 1000)
    amplitude_conv = amplitude(conversion_points, size)
    fig, ax = plt.subplots(ncols=1, nrows=2)
    ax[0].imshow(amplitude_conv, extent=(u, v, u, v))
    phase_con = phase(conversion_points, size)
    ax[1].imshow(phase_con, extent=(u, v, u, v))

    plt.savefig("амплитуда и фаза.png")
    plt.show()


if __name__ == '__main__':
    a, b = -3, 3
    u, v = -4, 4
    n, m = 1000, 1000
    integral = finite_integral_2D(a, b, n, m, u, v)
    chart(integral, u, v, n)
