import numpy as np
from numpy import pi, sin, cos, log, exp, sqrt
i = 1j

# 排查流程
# 1. 檢查cf是否有問題
# 2. 檢查積分區間的影響
# 3. 是否有必要進行St/ K 最後手段

class Polynomial:
    def __init__(self, S0, T, r, sigma, process, poly_coef, positive_interval, N, lower_limit, upper_limit,
                mean_reversion_speed=1.0, long_term_var_mean=0.16, heston_corr=-0.8, var_variance_process=2.0):
        self.S0 = S0
        self.x0 = log(S0)
        self.T = T
        self.r = r
        self.sigma = sigma
        self.process = process

        self.poly_coef = np.array(poly_coef)  # store the coefficients of the polynomial

        self.positive_interval = np.log(
            positive_interval)  # all the positive real solution the k0 can't below 0, kD can't above b

        # ============== Hyperparameter ==============
        self.N = N  # number of frequency for fitting density function
        self.lower_limit = log(lower_limit)
        self.upper_limit = log(upper_limit)

        # =============== Heston Model Parameter =======
        self.mean_reversion_speed = mean_reversion_speed
        self.long_term_var_mean = long_term_var_mean
        self.heston_corr = heston_corr
        self.var_variance_process = var_variance_process

    # 特徵函數
    # 為x0 = 0 下的特徵函數，要自己乘上exp(x0)
    # 但我後免在Pricing裡面已經改了故這邊不用乘
    def characteristic_function(self, u):
        if self.process == "GBM":
            cf_value = exp((self.r - 0.5 * self.sigma ** 2) * self.T * 1j * u -
                            0.5 * self.sigma ** 2 * u ** 2 * self.T)

            return cf_value

        elif self.process == "Heston":
            alpha = (self.mean_reversion_speed - i * self.heston_corr * self.var_variance_process * u) ** 2
            beta = (u ** 2 + i * u) * self.var_variance_process ** 2
            D = sqrt(alpha + beta)
            g_minus = self.mean_reversion_speed - i * self.heston_corr * self.var_variance_process * u - D
            g_plus = self.mean_reversion_speed - i * self.heston_corr * self.var_variance_process * u + D
            G = g_minus / g_plus
            # 不確定這裡的mu是啥 我假設為riskfree rate
            eta = 1 - G * exp(-D * self.T)
            gamma = i * u * self.r * self.T + self.sigma / self.var_variance_process ** 2 * (1 - exp(-D * self.T) / eta) * g_minus
            theta = self.mean_reversion_speed * self.long_term_var_mean / self.var_variance_process ** 2 * (self.T * g_minus - 2 * log(eta / (1 - G)))
            cf_value = exp(gamma) * exp(theta)

            return cf_value

        raise "Error"


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
    def PricingByCOSMethod(self):
        val = 0.0
        for k in range(self.N):
            if k == 0:
                val += 0.5 * np.real(self.characteristic_function(k * pi / (self.upper_limit - self.lower_limit)) *
                                                   exp(i * k * pi * (self.x0 - self.lower_limit) / (self.upper_limit - self.lower_limit))) * self.V(k)
            else:
                val += np.real(self.characteristic_function(k * pi / (self.upper_limit - self.lower_limit)) *
                                             exp(i * k * pi * (self.x0 - self.lower_limit) / (self.upper_limit - self.lower_limit))) * self.V(k)
        return exp(-self.r * self.T) * val

    # 樹模型
    def StockPrice(self, i, j, u, d):
        return self.S0 * (u ** (i - j) * d ** j)

    def Payoff(self, price):
        payoff = 0
        for power, coef in enumerate(self.poly_coef):
            payoff += coef * price ** power
        return max(payoff, 0)

    def PricingByTree(self):
        delta_t = self.T / self.N
        u = exp(self.sigma * (delta_t ** 0.5))
        d = 1 / u
        p = (exp(self.r * delta_t) - d) / (u - d)
        payoff_arr = np.zeros(self.N + 1).T
        for j in range(self.N + 1):
            payoff_arr[j] = self.Payoff(self.StockPrice(self.N, j, u, d))

        for i in reversed(list(range(self.N))):
            for j in range(i + 1):
                payoff_arr[j] = (p * payoff_arr[j] + (1 - p) * payoff_arr[j + 1]) * exp(-self.r * delta_t)
        return payoff_arr[0]


if __name__ == "__main__":
    polynomial2 = Polynomial(S0=100, T=10, r=0.0, sigma=0.0175, process="Heston",
                             N=320000, lower_limit=0.001, upper_limit=1000,
                             poly_coef=[-100, 1], positive_interval=[100, 1000],
                             mean_reversion_speed=1.5768, long_term_var_mean=0.0398,
                             var_variance_process=0.5751, heston_corr=-0.5711, )
    print(polynomial2.PricingByCOSMethod())






