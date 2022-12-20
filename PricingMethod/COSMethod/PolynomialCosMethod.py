import numpy as np
from numpy import pi, sin, cos
from PricingMethod.COSMethod.CharacteristicFunction import *
i = 1j

# 排查流程
# 1. 檢查cf是否有問題
# 2. 檢查積分區間的影響
# 3. 是否有必要進行St/ K 最後手段

class PolynomialCosMethod:
    def __init__(self, S0, T, r, sigma, process, poly_coef, positive_interval, N, lower_limit, upper_limit):
        self.S0 = S0
        self.x0 = log(S0)
        self.T = T
        self.r = r
        self.sigma = sigma
        self.process = process
        self.poly_coef = np.array(poly_coef)  # store the coefficients of the polynomial
        if positive_interval[0] == 0:
            positive_interval[0] = 1e-4
        self.positive_interval = np.log(
            positive_interval)  # all the positive real solution the k0 can't below 0, kD can't above b

        # ============== Hyperparameter ==============
        self.N = N  # number of frequency for fitting density function
        self.lower_limit = log(lower_limit)
        self.upper_limit = log(upper_limit)


    # 特徵函數
    # 為x0 = 0 下的特徵函數，要自己乘上exp(x0)
    # 但我後免在Pricing裡面已經改了故這邊不用乘
    def characteristic_function(self, u):
        cf_value = self.process.getCFValue(u)

        return cf_value


    def Chi(self, k, kd_down, kd_up):
        n = len(self.poly_coef)
        # 因為積分完之後，含有cos與sin的部份太冗長，所以用lambda定義一個小函數
        sin_integral = lambda x, v: sin(k * pi * (x - self.lower_limit) / (self.upper_limit - self.lower_limit)) * exp(v * x)
        cos_integral = lambda x, v: cos(k * pi * (x - self.lower_limit) / (self.upper_limit - self.lower_limit)) * exp(v * x)

        # 計算polynomial的每一項
        chi = 0
        for j in range(1, n):
            aj = self.poly_coef[j]
            chi += aj * (1 / (1 + (k * pi / (j * (self.upper_limit - self.lower_limit))) ** 2)) * (
                    1 / j * (cos_integral(kd_up, j) - cos_integral(kd_down, j)) +
                    k * pi / (j ** 2 * (self.upper_limit - self.lower_limit)) * (sin_integral(kd_up, j) - sin_integral(kd_down, j)))
        return chi

    # K的影響
    def Posi(self, k, kd_down, kd_up):
        a0 = self.poly_coef[0]
        sin_integral = lambda x: sin(k * pi * (x - self.lower_limit) / (self.upper_limit - self.lower_limit))
        if k == 0:
            return a0 * (kd_up - kd_down)
        else:
            return a0 * (self.upper_limit - self.lower_limit) / (k * pi) * (sin_integral(kd_up) - sin_integral(kd_down))

    # 計算不同頻率下的值
    def V(self, k):
        v = 0.0
        for d in range(int(len(self.positive_interval) / 2)):
            kd_down = self.positive_interval[2 * d]
            kd_up = self.positive_interval[2 * d + 1]
            v += (self.Chi(k, kd_down, kd_up) + self.Posi(k , kd_down, kd_up))
        return 2 / (self.upper_limit - self.lower_limit) * v

    # COS model
    def getValue(self):
        # 這邊的cf是假設 lnSt = 0 如果分配不能提出exp(iux)需要改
        val = 0.0
        for k in range(self.N):
            if k == 0:
                val += 0.5 * np.real(self.characteristic_function(k * pi / (self.upper_limit - self.lower_limit)) *
                                                   exp(i * k * pi * (self.x0 - self.lower_limit) / (self.upper_limit - self.lower_limit))) * self.V(k)
            else:
                val += np.real(self.characteristic_function(k * pi / (self.upper_limit - self.lower_limit)) *
                                             exp(i * k * pi * (self.x0 - self.lower_limit) / (self.upper_limit - self.lower_limit))) * self.V(k)
        return exp(-self.r * self.T) * val


if __name__ == "__main__":
    process_Heston = Heston(T=1, r=0.0, sigma=0.0175, mean_reversion_speed=1.5768, long_term_var_mean=0.0398,
                             var_variance_process=0.5751, corr=-0.5711)
    # process_GBM = GBM(T=1, r=0.1, sigma=0.25)
    # process_CGMY = CGMY(r=0.1, C=0, G=5, M=5, T=1, Y=0.5, sigma=0.0175)

    # process_VG = VG(theta=-0.14, sigma=0.12, nu=0.2, T=2)
    polynomialCosMethod = PolynomialCosMethod(S0=100, r=0.0, sigma=0.0175, T=1, process=process_Heston,
                             N=100, lower_limit=0.01, upper_limit=500,
                             poly_coef=[-100, 1], positive_interval=[100, 500])

    print(polynomialCosMethod.getValue())