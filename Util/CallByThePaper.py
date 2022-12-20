import numpy as np
from numpy import pi, sin, cos, log, exp

def characteristic_function(u):
    if process == "GBM.txt":
        cf_value = exp((r - 0.5 * sigma ** 2)) * 1j * u * T - \
                   0.5 * sigma ** 2 * u ** 2 * T


    else: 
        return -1

    return cf_value

def Chi(k, c, d):
    cos_part = cos(k * pi * (d - a) / (b - a)) * exp(d) - \
               cos(k * pi * (c - a) / (b - a) * exp(c))
    sin_part = sin(k * pi * (d - a) / (b - a)) * exp(d) - \
               sin(k * pi * (c - a) / (b - a)) * exp(c)
    chi = 1 / (1 + (k * pi / (b - a)) ** 2) * (cos_part +
                                               k * pi / (b - a) * sin_part)
    return chi

def Phi(k, c, d):
    if k == 0:
        phi = d - c
    else:
        phi = sin(k * pi * ((d - a) / (b - a))) - sin(k * pi * (c - a) / (b - a)) * (b - a) / (k * pi)
    return phi

def Vk(k):
    vk = 2 / (b - a) * K * (Chi(k, 0 , b) - Phi(k, 0, b))
    return vk

def pricing():
    val_sum = 0
    for k in range(N):
        if k == 0:
            val_sum += 0.5 * (characteristic_function(k * pi / (b - a)) * exp(-1j * k * pi * a / (b - a))).real * Vk(k)
        else:
            val_sum += (characteristic_function(k * pi / (b - a)) * exp(-1j * k * pi * a / (b - a))).real * Vk(k)
    value = exp(-r * T) * val_sum
    return value

if __name__ == "__main__":
    # ============== Market Data ==============
    S0 = 100
    T = 1
    r = 0.03
    sigma = 0.05
    K = 100
    # ============== Option Setting ==============
    process = "GBM.txt"
    # ============== Hyperparameter ==============
    integral_range = (0, 2000) # based on price <--- need to check!!
    N = 1000 # number of frequency for fitting density function
    # ============== Preproccessing ==============
    x0 = log(S0)
    a, b = log(integral_range)
    # ============== Pricing ==============
    print(pricing())





