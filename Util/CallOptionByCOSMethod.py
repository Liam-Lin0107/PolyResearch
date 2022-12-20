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

def Chi(k):
    chi = 1 / (1 + (k * pi / (b - a)) ** 2) * \
          (cos(k * pi) * exp(b) - cos(k * pi * (log(K) - a) / (b - a)) * K +
           k * pi / (b - a) * (sin(k * pi) * exp(b) - sin(k * pi * (log(K) - a) / (b - a)) * K))
    return chi

def Posi(k):
    if k == 0:
        return K * (b - log(K))
    else:
        return K * (b - a) / (k * pi) * (sin(k * pi) - sin(k * pi * (log(K) - a) / (b - a)))

def V(k):
    return 2 / (b - a) * (Chi(k) - Posi(k))

def MainPrcing():
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
    S0 = 100
    T = 0.1
    r = 0.1
    sigma = 0.25
    K = 80
    # ============== Option Setting ==============
    process = "GBM.txt"
    # ============== Hyperparameter ==============
    N = 160 # number of frequency for fitting density function
    # ============== Preproccessing ==============
    x0 = log(S0)
    a, b = -10, 10
    # ============== Pricing ==============
    print(MainPrcing())





