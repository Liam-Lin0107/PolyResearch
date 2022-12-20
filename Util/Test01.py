import numpy as np
from numpy import pi, sin, cos, log, exp
i = 1j
class PolyByCOSMethod:
    def __init__(self, S0, T, r, sigma, poly_coef, positive_interval, process="BGM", number_of_frequency=256, lower_limit=-10, upper_limit=10):
        # ============== Market Data ==============
        self.x0 = log(S0)
        self.T = T
        self.r = r
        self.sigma = sigma
        self.process = process

        # ============== Option Setting ==============
        self.poly_coef = np.array(poly_coef)  # store the coefficients of the polynomial
        self.positive_interval = np.log(positive_interval)  # all the positive real solution the k0 can't below 0, kD can't above b

        # ============== Hyperparameter ==============
        self.N = number_of_frequency  # number of frequency for fitting density function
        self.a = lower_limit
        self.b = upper_limit

    def characteristic_function(self, u):
        if self.process == "GBM.txt":
            cf_value = exp(((self.x0 / self.T + self.r - 0.5 * self.sigma ** 2) * 1j * u -
                            0.5 * self.sigma ** 2 * u ** 2) * self.T)
        else:
            raise "Error"
        return cf_value

    def Chi(self, k, kd_down, kd_up):
        n = len(self.poly_coef) - 1
        sin_integral = lambda x, v: sin(k * pi * (x - self.a) / (self.b - self.a)) * exp(v * x)
        cos_integral = lambda x, v: cos(k * pi * (x - self.a) / (self.b - self.a)) * exp(v * x)
        chi = 0
        for j in range(1, n + 1):
            aj = self.poly_coef[j]
            chi += aj * (1 / (1 + (k * pi / (j * (self.b - self.a))) ** 2)) * (
                    1 / j * (cos_integral(kd_up, j) - cos_integral(kd_down, j)) +
                    k * pi / (j ** 2 * (self.b - self.a)) * (sin_integral(kd_up, j) - sin_integral(kd_down, j)))
        return chi

    def Posi(self, k, kd_down, kd_up):
        a0 = self.poly_coef[0]
        sin_integral = lambda x: sin(k * pi * (x - self.a) / (self.b - self.a))
        if k == 0:
            return a0 * (kd_up - kd_down)
        else:
            return a0 * (self.b - self.a) / (k * pi) *(sin_integral(kd_up) - sin_integral(kd_down))

    def V(self, k):
        v = 0
        for d in range(int(len(self.poly_coef) / 2)):
            kd_down = self.positive_interval[2 * d]
            kd_up = self.positive_interval[2 * d + 1]
            v += (self.Chi(k, kd_down, kd_up) + self.Posi(k , kd_down, kd_up))
        return 2 / (self.b - self.a) * v


    def Pricing(self):
        val = 0
        for k in range(self.N):
            if k == 0:
                val += 0.5 * exp(-self.r * self.T) * np.real(self.characteristic_function(k * pi / (self.b - self.a)) *
                                                   exp(-i * k * self.a * pi / (self.b - self.a))) * self.V(k)
            else:
                val += exp(-self.r * self.T) * np.real(self.characteristic_function(k * pi / (self.b - self.a)) *
                                             exp(-i * k * self.a * pi / (self.b - self.a))) * self.V(k)
        return val

if __name__ == "__main__":
    polynomial = PolyByCOSMethod(S0=5, T=1, r=0.1, sigma=0.25,
                                 process="GBM.txt",
                                 poly_coef=[840, -638, 179, -22, 1],
                                 positive_interval=[exp(-9), 4, 5, 6, 7, exp(9)],
                                 number_of_frequency=1000,
                                 lower_limit=-10, upper_limit=10)
    print("Price: ", polynomial.Pricing())





