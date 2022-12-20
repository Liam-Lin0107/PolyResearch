import numpy as np
from numpy import pi, sin, cos, log, exp
i = 1j

def characteristic_function(u):
    if process == "GBM.txt":
        cf_value = exp(((x0 / T + r - 0.5 * sigma ** 2) * 1j * u -
                        0.5 * sigma ** 2 * u ** 2) * T)
    else:
        raise "Error"
    return cf_value

def Chi(k, kd_down, kd_up):
    n = len(poly_coef) - 1
    sin_integral = lambda x, v: sin(k * pi * (x - a) / (b - a)) * exp(v * x)
    cos_integral = lambda x, v: cos(k * pi * (x - a) / (b - a)) * exp(v * x)
    chi = 0
    for j in range(1, n + 1):
        aj = poly_coef[j]
        chi += aj * (1 / (1 + (k * pi / (j * (b - a))) ** 2)) * (
                1 / j * (cos_integral(kd_up, j) - cos_integral(kd_down, j)) +
                k * pi / (j ** 2 * (b - a)) * (sin_integral(kd_up, j) - sin_integral(kd_down, j)))
    return chi

def Posi(k, kd_down, kd_up):
    a0 = poly_coef[0]
    sin_integral = lambda x: sin(k * pi * (x - a) / (b - a))
    if k == 0:
        return a0 * (kd_up - kd_down)
    else:
        return a0 * (b - a) / (k * pi) *(sin_integral(kd_up) - sin_integral(kd_down))

def V(k):
    v = 0
    for d in range(D):
        kd_down = K[2 * d]
        kd_up = K[2 * d + 1]
        v += (Chi(k, kd_down, kd_up) + Posi(k , kd_down, kd_up))
    return 2 / (b - a) * v

def MainPricing():
    val = 0
    for k in range(N):
        if k == 0:
            val += 0.5 * exp(-r * T) * np.real(characteristic_function(k * pi / (b - a)) *
                                               exp(-i * k * a * pi / (b - a))) * V(k)
        else:
            val += exp(-r * T) * np.real(characteristic_function(k * pi / (b - a)) *
                                         exp(-i * k * a * pi / (b - a))) * V(k)
    return val

if __name__ == "__main__":
    # ============== Market Data ==============
    S0 = 5
    T = 1
    r = 0.1
    sigma = 0.25
    process = "GBM.txt"

    # ============== Option Setting ==============
    poly_coef = np.array([840, -638, 179, -22, 1]) # store the coefficients of the polynomial
    D = 3
    K = np.log([exp(-9), 4, 5, 6, 7, exp(9)]) # all the positive real solution the k0 can't below 0, kD can't above b

    # ============== Hyperparameter ==============
    N = 1000 # number of frequency for fitting density function
    a, b = -10, 10

    # ============== Preproccessing ==============
    x0 = log(S0)

    # ============== Pricing ==============
    print("Price: ", MainPricing())





