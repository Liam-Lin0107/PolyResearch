import math

import numpy as np
from numpy import log, exp, sqrt
import warnings
from scipy import stats
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt


class GBMByMC:
    def __init__(self, S0, r, sigma, T, poly_coeff,
                 N_line=int(1e6), n=252, N_repeat=20):
        self.S0 = S0
        self.r = r
        self.sigma = sigma
        self.T = T

        self.poly_coeff = poly_coeff

        self.N_line = N_line
        self.n = n
        self.N_repeat = N_repeat
        # 以下用來檢查蒙地卡羅是否有問題(可以不用加上去)
        self.maxST = []
        self.minST = []
    
    # 選擇權的payoff計算
    def _payoff(self, price):
        payoff = 0
        for power, coeff in enumerate(self.poly_coeff):
            payoff += coeff * price ** power
        return max(payoff, 0)
    
    # 計算一次蒙地卡羅的價格
    def getValue(self, repeat):
        dt = self.T / self.n
        # 生成S0的array
        lnSt = np.full(self.N_line, log(self.S0))
        # 每次dt的變化
        for i in range(self.n):
            if i % 20 == 0:
                print(f"GBM: {repeat + 1} round, {((i + 1) / self.n) * 100:.1f}%")

            # antithetic + moment match
            norm_rv = np.random.randn(int(self.N_line / 2))
            norm_rv = np.append(norm_rv, -norm_rv)
            norm_rv = norm_rv / np.std(norm_rv)

            dW = norm_rv * pow(dt, 0.5)
            # 計算dS
            dlnSt = (self.r - pow(self.sigma, 2) / 2) * dt + self.sigma * dW
        
            lnSt = lnSt + dlnSt
        # 因為檢查方便故進ST排序
        ST = np.sort(exp(lnSt))
        # 放入極端值以利後續檢查
        self.maxST.append(np.round(ST[-5:], 2).tolist())
        self.minST.append(np.round(ST[:5], 2).tolist())
        # 將所有股價利用_payoff方法maping成payoff，利用mapping會更快
        payoff = list(map(self._payoff, ST))
        value = np.mean(payoff) * exp(-self.r * self.T)
        return value
    
    # 進行多次蒙地卡羅:為使用者調用的方法，回傳平均值與標準差
    # save_data如果為True則會建立文件存放資料
    def getStatistic(self, file_name="Heston", save_data=False, save_dir=""):
        values = []
        print("GBM simulation starting...")
        for i in range(self.N_repeat):
            values.append(self.getValue(i))

        mean = np.mean(values).item()
        std = np.std(values).item()
        lower_bound = mean - 2 * std
        upper_bound = mean + 2 * std
        if save_data:
            with open(save_dir, "a+") as file:
                max_ST = ',\n\t\t'.join(str(e) for e in self.maxST)
                min_ST = ',\n\t\t'.join(str(e) for e in self.minST)
                s = f"\nPolynomial Coefficient: {self.poly_coeff}\n" \
                    f"Basic:\n" \
                    f"\tS0: {self.S0}, T: {self.T}, sigma: {self.sigma}\n" \
                    f"Result:\n" \
                    f"\tC.V.: ({lower_bound:.8f}, {upper_bound:.8f})\n" \
                    f"\tmean: {mean}\n" \
                    f"\tstd: {std}\n" \
                    f"\tmax: \n\t\t{max_ST}\n" \
                    f"\tmin: \n\t\t{min_ST}\n"

                file.write(s)
        return mean, std


class HestonByMC:
    def __init__(self, S0, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, std_of_var_process,
                 poly_coeff,
                 N_line=int(1e6), n=252, N_repeat=20):
        self.S0 = S0
        self.r = r
        self.sigma = sigma
        self.T = T
        self.mean_reversion_speed = mean_reversion_speed
        self.long_term_var_mean = long_term_var_mean
        self.corr = corr
        self.std_of_var_process = std_of_var_process

        self.poly_coeff = poly_coeff

        self.N_line = N_line
        self.n = n
        self.N_repeat = N_repeat

        self.maxST = []
        self.minST = []

    def _payoff(self, price):
        payoff = 0
        for power, coeff in enumerate(self.poly_coeff):
            payoff += coeff * price ** power
        return max(payoff, 0)

    def getValue(self, repeat):
        dt = self.T / self.n

        lnSt = np.full(self.N_line, log(self.S0))
        Vt = np.full(self.N_line, pow(self.sigma, 2))

        for i in range(self.n):
            if i % 20 == 0:
                print(f"Heston: {repeat + 1} round, {((i + 1) / self.n) * 100:.1f}%")

            # antithetic + moment match
            norm_rv1 = np.random.randn(int(self.N_line / 2))
            norm_rv1 = np.append(norm_rv1, -norm_rv1)
            norm_rv1 = norm_rv1 / np.std(norm_rv1)

            norm_rv2 = np.random.randn(int(self.N_line / 2))
            norm_rv2 = np.append(norm_rv2, -norm_rv2)
            norm_rv2 = norm_rv2 / np.std(norm_rv2)

            dW1t = norm_rv1 * pow(dt, 0.5)
            dW2t = pow(dt, 0.5) * (norm_rv1 * self.corr + norm_rv2 * pow(1 - pow(self.corr, 2), 0.5))

            dlnSt = (self.r - 0.5 * Vt) * dt + sqrt(Vt) * dW1t
            lnSt = lnSt + dlnSt
            Vt = Vt + (self.mean_reversion_speed * (self.long_term_var_mean - Vt) * dt + self.std_of_var_process * sqrt(
                Vt) * dW2t)
            Vt = np.where(Vt > 0, Vt, 0)

        ST = np.sort(exp(lnSt))
        self.maxST.append(np.round(ST[-5:], 2).tolist())
        self.minST.append(np.round(ST[:5], 2).tolist())
        payoff = list(map(self._payoff, ST))
        value = np.mean(payoff) * exp(-self.r * self.T)
        return value

    def getStatistic(self, save_data=False, save_dir=""):
        print("Heston simulation starting...")
        values = []
        for i in range(self.N_repeat):
            values.append(self.getValue(i))

        mean = np.mean(values).item()
        std = np.std(values).item()
        lower_bound = mean - 2 * std
        upper_bound = mean + 2 * std
        if save_data:
            with open(save_dir, "a+") as file:
                max_ST = ',\n\t\t'.join(str(e) for e in self.maxST)
                min_ST = ',\n\t\t'.join(str(e) for e in self.minST)
                s = f"\nPolynomial Coefficient: {self.poly_coeff}\n" \
                    f"Basic:\n" \
                    f"\tS0: {self.S0}, T: {self.T}, sigma: {self.sigma}\n" \
                    f"SV related Coeff:\n" \
                    f"\tmean_reversion: {self.mean_reversion_speed}, long_term_var_mean: {self.long_term_var_mean}, corr: {self.corr}\n" \
                    f"Result:\n" \
                    f"\tC.V.: ({lower_bound:.8f}, {upper_bound:.8f})\n" \
                    f"\tmean: {mean}\n" \
                    f"\tstd: {std}\n" \
                    f"\tmax: \n\t\t{max_ST}\n" \
                    f"\tmin: \n\t\t{min_ST}\n"

                file.write(s)
        return mean, std


class MJDByMC:
    def __init__(self, S0, r, sigma, T, poly_coeff, jump_intensity, jump_mean, jump_var,
                 N_line=int(1e6), n=252, N_repeat=20):
        self.S0 = S0
        self.r = r
        self.sigma = sigma
        self.T = T
        self.poly_coeff = poly_coeff
        self.jump_intensity = jump_intensity
        self.jump_mean = jump_mean
        self.jump_var = jump_var

        self.N_line = N_line
        self.n = n
        self.N_repeat = N_repeat

        self.maxST = []
        self.minST = []
        self.jumps = np.array([0])  # record jumps
        self.jumps_count = np.array([0])  # record jumps count

    def _payoff(self, price):
        payoff = 0
        for power, coeff in enumerate(self.poly_coeff):
            payoff += coeff * price ** power
        return max(payoff, 0)

    def getValue(self, repeat):
        dt = self.T / self.n

        lnSt = np.full(self.N_line, log(self.S0))

        k = exp(self.jump_mean + 0.5 * self.jump_var) - 1  # k = E(Y - 1)
        mu = self.r - self.jump_intensity * k  # risk neutral adjust rate

        for i in range(self.n):
            if i % 20 == 0:
                print(f"MDJ: {repeat + 1} round, {((i + 1) / self.n) * 100:.1f}%")

            # antithetic + moment match
            norm_rv1 = np.random.randn(int(self.N_line / 2))
            norm_rv1 = np.append(norm_rv1, -norm_rv1)
            norm_rv1 = norm_rv1 / np.std(norm_rv1)

            norm_rv2 = np.random.randn(int(self.N_line / 2))
            norm_rv2 = np.append(norm_rv2, -norm_rv2)
            norm_rv2 = norm_rv2 / np.std(norm_rv2)

            dWt = norm_rv1 * pow(dt, 0.5)
            lnYt = norm_rv2 * sqrt(self.jump_var) + self.jump_mean

            jumps = np.random.poisson(self.jump_intensity * dt, self.N_line)
            # deal with jump > 1
            over_jump_ind = np.argwhere(jumps > 1)
            J = lnYt * jumps
            for j in over_jump_ind:
                J[j] = np.sum(np.random.normal(self.jump_mean, sqrt(self.jump_var), size=jumps[j]))

            dlnSt = (mu - 0.5 * pow(self.sigma, 2)) * dt + self.sigma * dWt + J
            lnSt = lnSt + dlnSt

            # compute the jump in order to analysis
            jumps_count = np.unique(jumps, return_counts=True)
            self.jumps = np.union1d(jumps_count[0], self.jumps)
            if len(jumps_count[1]) >= len(self.jumps_count):
                for j in range(len(self.jumps_count)):
                    temp = self.jumps_count.copy()
                    self.jumps_count = jumps_count[1]
                    self.jumps_count[j] += temp[j]
                else:
                    for j in range(len(jumps_count[1])):
                        self.jumps_count[j] += jumps_count[1][j]

        ST = np.sort(exp(lnSt))
        self.maxST.append(np.round(ST[-5:], 2).tolist())
        self.minST.append(np.round(ST[:5], 2).tolist())
        payoff = list(map(self._payoff, ST))
        value = np.mean(payoff) * exp(-self.r * self.T)
        return value

    def getStatistic(self, save_data=False, save_dir=""):
        values = []
        print("MJD simulation starting...")
        for i in range(self.N_repeat):
            values.append(self.getValue(i))
        mean = np.mean(values).item()
        std = np.std(values).item()
        lower_bound = mean - 2 * std
        upper_bound = mean + 2 * std
        if save_data:
            with open(save_dir, "a+") as file:
                max_ST = ',\n\t\t'.join(str(e) for e in self.maxST)
                min_ST = ',\n\t\t'.join(str(e) for e in self.minST)

                jumps = self.jumps.tolist()

                jumps_count = self.jumps_count.tolist()
                s = f"\nPolynomial Coefficient: {self.poly_coeff}\n" \
                    f"Basic:\n" \
                    f"\tS0: {self.S0}, T: {self.T}, sigma: {self.sigma}\n" \
                    f"Jump related Coeff:\n" \
                    f"\tjump_intensity: {self.jump_intensity}, jump_mean: {self.jump_mean}, jump_var: {self.jump_var}\n" \
                    f"Result:\n" \
                    f"\tC.V.: ({lower_bound:.8f}, {upper_bound:.8f})\n" \
                    f"\tmean: {mean}\n" \
                    f"\tstd: {std}\n" \
                    f"\tjumps: \n\t\t{jumps}\n\n\t\t{jumps_count}\n" \
                    f"\tmax: \n\t\t{max_ST}\n" \
                    f"\tmin: \n\t\t{min_ST}\n"
                file.write(s)
        return mean, std


class KJDByMC:
    def __init__(self, S0, r, sigma, T, poly_coeff, jump_intensity, p, eta1, eta2,
                 N_line=int(1e6), n=252, N_repeat=20):
        self.S0 = S0
        self.r = r
        self.sigma = sigma
        self.T = T
        self.poly_coeff = poly_coeff
        self.jump_intensity = jump_intensity
        self.p = p
        self.eta1 = eta1
        self.eta2 = eta2

        self.N_line = N_line
        self.n = n
        self.N_repeat = N_repeat

        self.maxST = []
        self.minST = []
        self.jumps = np.array([0])  # record jumps
        self.jumps_count = np.array([0])  # record jumps count

    def _payoff(self, price):
        payoff = 0
        for power, coeff in enumerate(self.poly_coeff):
            payoff += coeff * price ** power
        return max(payoff, 0)

    def getValue(self, repeat):
        dt = self.T / self.n

        lnSt = np.full(self.N_line, log(self.S0))

        k = self.p * self.eta1 / (self.eta1 - 1) + (1 - self.p) * self.eta2 / (self.eta2 + 1) - 1  # k = E(Y - 1)
        mu = self.r - self.jump_intensity * k  # risk neutral adjust rate

        for i in range(self.n):
            if i % 20 == 0:
                print(f"KJD: {repeat + 1} round, {((i + 1) / self.n) * 100:.1f}%")

            # antithetic + moment match
            norm_rv = np.random.randn(int(self.N_line / 2))
            norm_rv = np.append(norm_rv, -norm_rv)
            norm_rv = norm_rv / np.std(norm_rv)
            dWt = norm_rv * pow(dt, 0.5)

            jumps = np.random.poisson(self.jump_intensity * dt, size=self.N_line)

            exp_rv1 = np.random.exponential(1 / self.eta1, size=self.N_line)
            exp_rv2 = np.random.exponential(1 / self.eta2, size=self.N_line)
            up = np.random.binomial(n=1, p=self.p, size=self.N_line)
            lnYt = up * exp_rv1 + (1 - up) * (-exp_rv2)

            dlnSt = (mu - 0.5 * pow(self.sigma, 2)) * dt + self.sigma * dWt + lnYt * jumps
            lnSt = lnSt + dlnSt

            # compute the jump in order to analysis
            jumps_count = np.unique(jumps, return_counts=True)
            self.jumps = np.union1d(jumps_count[0], self.jumps)
            if len(jumps_count[1]) >= len(self.jumps_count):
                for j in range(len(self.jumps_count)):
                    temp = self.jumps_count.copy()
                    self.jumps_count = jumps_count[1]
                    self.jumps_count[j] += temp[j]
                else:
                    for j in range(len(jumps_count[1])):
                        self.jumps_count[j] += jumps_count[1][j]

        ST = np.sort(exp(lnSt))
        self.maxST.append(np.round(ST[-5:], 2).tolist())
        self.minST.append(np.round(ST[:5], 2).tolist())
        payoff = list(map(self._payoff, ST))
        value = np.mean(payoff) * exp(-self.r * self.T)
        return value

    def getStatistic(self, save_data=False, save_dir=""):
        values = []
        print("KJD simulation starting...")
        for i in range(self.N_repeat):
            values.append(self.getValue(i))

        mean = np.mean(values).item()
        std = np.std(values).item()
        lower_bound = mean - 2 * std
        upper_bound = mean + 2 * std
        if save_data:
            with open(save_dir, "a+") as file:
                max_ST = ',\n\t\t'.join(str(e) for e in self.maxST)
                min_ST = ',\n\t\t'.join(str(e) for e in self.minST)

                jumps = self.jumps.tolist()

                jumps_count = self.jumps_count.tolist()
                s = f"\nPolynomial Coefficient: {self.poly_coeff}\n" \
                    f"Basic:\n" \
                    f"\tS0: {self.S0}, T: {self.T}, sigma: {self.sigma}\n" \
                    f"Jump related Coeff:\n" \
                    f"\tjump_intensity: {self.jump_intensity}, p: {self.p}, eta1: {self.eta1}, eta2: {self.eta2}\n" \
                    f"Result:\n" \
                    f"\tC.V.: ({lower_bound:.8f}, {upper_bound:.8f})\n" \
                    f"\tmean: {mean}\n" \
                    f"\tstd: {std}\n" \
                    f"\tjumps: \n\t\t{jumps}\n\n\t\t{jumps_count}\n" \
                    f"\tmax: \n\t\t{max_ST}\n" \
                    f"\tmin: \n\t\t{min_ST}\n"
                file.write(s)
        return mean, std


class SVJByMC:
    def __init__(self, S0, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, std_of_var_process,
                 poly_coeff, jump_intensity, jump_mean, jump_var,
                 N_line=int(1e6), n=252, N_repeat=20):
        self.S0 = S0
        self.r = r
        self.sigma = sigma
        self.T = T
        self.mean_reversion_speed = mean_reversion_speed
        self.long_term_var_mean = long_term_var_mean
        self.corr = corr
        self.std_of_var_process = std_of_var_process
        self.poly_coeff = poly_coeff
        self.jump_intensity = jump_intensity
        self.jump_mean = jump_mean
        self.jump_var = jump_var

        self.N_line = N_line
        self.n = n
        self.N_repeat = N_repeat

        self.maxST = []
        self.minST = []
        self.jumps = np.array([0])  # record jumps
        self.jumps_count = np.array([0])  # record jumps count

    def _payoff(self, price):
        payoff = 0
        for power, coeff in enumerate(self.poly_coeff):
            payoff += coeff * price ** power
        return max(payoff, 0)

    def getValue(self, repeat):
        dt = self.T / self.n

        lnSt = np.full(self.N_line, log(self.S0))
        Vt = np.full(self.N_line, pow(self.sigma, 2))

        k = exp(self.jump_mean + 0.5 * self.jump_var) - 1  # k = E(Y - 1)
        mu = self.r - self.jump_intensity * k  # risk neutral adjust rate

        for i in range(0, self.n):
            if i % 20 == 0:
                print(f"SVJ: {repeat + 1} round, {((i + 1) / self.n) * 100:.1f}%")

            # antithetic + moment match
            norm_rv1 = np.random.randn(int(self.N_line / 2))
            norm_rv1 = np.append(norm_rv1, -norm_rv1)
            norm_rv1 = norm_rv1 / np.std(norm_rv1)

            norm_rv2 = np.random.randn(int(self.N_line / 2))
            norm_rv2 = np.append(norm_rv2, -norm_rv2)
            norm_rv2 = norm_rv2 / np.std(norm_rv2)

            norm_rv3 = np.random.randn(int(self.N_line / 2))
            norm_rv3 = np.append(norm_rv3, -norm_rv3)
            norm_rv3 = norm_rv3 / np.std(norm_rv3)

            dW1t = norm_rv1 * pow(dt, 0.5)
            dW2t = pow(dt, 0.5) * (norm_rv1 * self.corr + norm_rv2 * pow(1 - pow(self.corr, 2), 0.5))
            lnYt = norm_rv3 * sqrt(self.jump_var) + self.jump_mean

            jumps = np.random.poisson(self.jump_intensity * dt, self.N_line)
            # deal with jump > 1
            over_jump_ind = np.argwhere(jumps)
            J = lnYt * jumps
            for j in over_jump_ind:
                J[j] = np.sum(np.random.normal(self.jump_mean,
                                               sqrt(self.jump_var), size=jumps[j]))

            dlnSt = (mu - 0.5 * Vt) * dt + sqrt(Vt) * dW1t + J
            lnSt = lnSt + dlnSt

            Vt = Vt + (self.mean_reversion_speed * (self.long_term_var_mean - Vt) * dt + self.std_of_var_process * sqrt(
                Vt) * dW2t)
            Vt = np.where(Vt > 0, Vt, 0)

            # compute the jump in order to analysis
            jumps_count = np.unique(jumps, return_counts=True)
            self.jumps = np.union1d(jumps_count[0], self.jumps)
            if len(jumps_count[1]) >= len(self.jumps_count):
                for j in range(len(self.jumps_count)):
                    temp = self.jumps_count.copy()
                    self.jumps_count = jumps_count[1]
                    self.jumps_count[j] += temp[j]
                else:
                    for j in range(len(jumps_count[1])):
                        self.jumps_count[j] += jumps_count[1][j]

        ST = np.sort(exp(lnSt))
        self.maxST.append(np.round(ST[-5:], 2).tolist())
        self.minST.append(np.round(ST[:5], 2).tolist())
        payoff = list(map(self._payoff, ST))
        value = np.mean(payoff) * exp(-self.r * self.T)
        return value

    def getStatistic(self, save_data=False, save_dir=""):
        values = []
        print("SVJ simulation starting...")
        for i in range(self.N_repeat):
            values.append(self.getValue(i))

        mean = np.mean(values).item()
        std = np.std(values).item()
        lower_bound = mean - 2 * std
        upper_bound = mean + 2 * std
        if save_data:
            with open(save_dir, "a+") as file:
                max_ST = ',\n\t\t'.join(str(e) for e in self.maxST)
                min_ST = ',\n\t\t'.join(str(e) for e in self.minST)

                jumps = self.jumps.tolist()
                jumps_count = self.jumps_count.tolist()
                s = f"\nPolynomial Coefficient: {self.poly_coeff}\n" \
                    f"Basic:\n" \
                    f"\tS0: {self.S0}, T: {self.T}, sigma: {self.sigma}\n" \
                    f"SV related Coeff:\n" \
                    f"\tmean_reversion: {self.mean_reversion_speed}, long_term_var_mean: {self.long_term_var_mean}, corr: {self.corr}\n" \
                    f"Jump related Coeff:\n" \
                    f"\tjump_intensity: {self.jump_intensity}, jump_mean: {self.jump_mean}, jump_var: {self.jump_var}\n" \
                    f"Result:\n" \
                    f"\tC.V.: ({lower_bound:.8f}, {upper_bound:.8f})\n" \
                    f"\tmean: {mean}\n" \
                    f"\tstd: {std}\n" \
                    f"\tjumps: \n\t\t{jumps}\n\n\t\t{jumps_count}\n" \
                    f"\tmax: \n\t\t{max_ST}\n" \
                    f"\tmin: \n\t\t{min_ST}\n"
                file.write(s)
        return mean, std


class VGByMC:
    def __init__(self, S0, r, sigma, T, gamma_mean, gamma_var,
                 poly_coeff,
                 N_line=int(1e6), n=252, N_repeat=20):
        self.S0 = S0
        self.r = r
        self.sigma = sigma
        self.T = T

        self.gamma_mean = gamma_mean
        self.gamma_var = gamma_var

        self.poly_coeff = poly_coeff

        self.N_line = N_line
        self.n = n
        self.N_repeat = N_repeat

        self.maxST = []
        self.minST = []

    def _payoff(self, price):
        payoff = 0
        for power, coeff in enumerate(self.poly_coeff):
            payoff += coeff * price ** power
        return max(payoff, 0)

    def Xt_cf(self, u):
        return pow(1 - 1j * u * self.gamma_mean * self.gamma_var + 0.5 * self.sigma ** 2 * self.gamma_var * u ** 2,
                   -self.T / self.gamma_var)

    def getValue(self, repeat):
        # this method from Wikipwdia simulation 1
        # https://en.wikipedia.org/wiki/Variance_gamma_process

        mu = self.r - log(self.Xt_cf(-1j)) / self.T
        Xt = np.zeros(self.N_line)


        # Gamma process rv.
        G = np.random.gamma(shape=self.T / self.gamma_var, scale=self.gamma_var, size=self.N_line)
        # antithetic + moment match
        Z = np.random.randn(int(self.N_line / 2))
        Z = np.append(Z, -Z)
        Z = Z / np.std(Z)
        # random time normal increment
        W = Z * np.sqrt(G)
        # Brownian motion with random time
        Xt = self.gamma_mean * G + self.sigma * W


        ST = np.sort(np.real(self.S0 * exp(mu * self.T + Xt)))
        self.maxST.append(np.round(ST[-5:], 2).tolist())
        self.minST.append(np.round(ST[:5], 2).tolist())
        payoff = list(map(self._payoff, ST))
        value = np.mean(payoff) * exp(-self.r * self.T)
        return value

    def getStatistic(self, save_data=False, save_dir=""):
        print("Variance gamma simulation starting...")
        values = []
        for i in range(self.N_repeat):
            print(f"VG: {i + 1} round.")
            values.append(self.getValue(i))

        mean = np.mean(values).item()
        std = np.std(values).item()
        lower_bound = mean - 2 * std
        upper_bound = mean + 2 * std
        if save_data:
            with open(save_dir, "a+") as file:
                max_ST = ',\n\t\t'.join(str(e) for e in self.maxST)
                min_ST = ',\n\t\t'.join(str(e) for e in self.minST)
                s = f"\nPolynomial Coefficient: {self.poly_coeff}\n" \
                    f"Basic:\n" \
                    f"\tS0: {self.S0}, T: {self.T}, sigma: {self.sigma}\n" \
                    f"VG related Coeff:\n" \
                    f"\tgamma_mean: {self.gamma_mean}, gamma_var: {self.gamma_var}\n" \
                    f"Result:\n" \
                    f"\tC.V.: ({lower_bound:.8f}, {upper_bound:.8f})\n" \
                    f"\tmean: {mean}\n" \
                    f"\tstd: {std}\n" \
                    f"\tmax: \n\t\t{max_ST}\n" \
                    f"\tmin: \n\t\t{min_ST}\n"

                file.write(s)
        return mean, std


class NIGByMC:
    def __init__(self, S0, r, sigma, T, delta, alpha, beta, poly_coeff, N_line=int(1e6), n=252, N_repeat=20):
        self.S0 = S0
        self.r = r
        self.sigma = sigma
        self.T = T

        self.delta = delta
        self.alpha = alpha
        self.beta = beta

        self.poly_coeff = poly_coeff

        self.N_line = N_line
        self.n = n
        self.N_repeat = N_repeat

        self.maxST = []
        self.minST = []

    def _payoff(self, price):
        payoff = 0
        for power, coeff in enumerate(self.poly_coeff):
            payoff += coeff * price ** power
        return max(payoff, 0)

    def Xt_cf(self, u):
        return exp(self.delta * self.T * (pow(self.alpha ** 2 - self.beta ** 2, 0.5) - pow(self.alpha ** 2 - (self.beta + 1j * u) ** 2, 0.5)))

    def getValue(self, repeat):

        mu = self.r - log(self.Xt_cf(-1j)) / self.T
        Xt = stats.norminvgauss.rvs(a=self.alpha * self.delta * self.T, b=self.beta * self.delta * self.T, scale=self.delta * self.T, size=self.N_line)

        ST = np.sort(np.real(self.S0 * exp(mu * self.T + Xt)))
        self.maxST.append(np.round(ST[-5:], 2).tolist())
        self.minST.append(np.round(ST[:5], 2).tolist())
        payoff = list(map(self._payoff, ST))
        value = np.mean(payoff) * exp(-self.r * self.T)
        return value

    def getStatistic(self, save_data=False, save_dir=""):
        print("Normal Inverse Gaussian simulation starting...")
        values = []
        for i in range(self.N_repeat):
            print(f"NIG: {i + 1} round.")
            values.append(self.getValue(i))

        mean = np.mean(values).item()
        std = np.std(values).item()
        lower_bound = mean - 2 * std
        upper_bound = mean + 2 * std
        if save_data:
            with open(save_dir, "a+") as file:
                max_ST = ',\n\t\t'.join(str(e) for e in self.maxST)
                min_ST = ',\n\t\t'.join(str(e) for e in self.minST)
                s = f"\nPolynomial Coefficient: {self.poly_coeff}\n" \
                    f"Basic:\n" \
                    f"\tS0: {self.S0}, T: {self.T}, sigma: {self.sigma}\n" \
                    f"NIG related Coeff:\n" \
                    f"\tdelta: {self.delta}, alpha: {self.alpha}, beta: {self.beta}\n" \
                    f"Result:\n" \
                    f"\tC.V.: ({lower_bound:.8f}, {upper_bound:.8f})\n" \
                    f"\tmean: {mean}\n" \
                    f"\tstd: {std}\n" \
                    f"\tmax: \n\t\t{max_ST}\n" \
                    f"\tmin: \n\t\t{min_ST}\n"

                file.write(s)
        return mean, std

