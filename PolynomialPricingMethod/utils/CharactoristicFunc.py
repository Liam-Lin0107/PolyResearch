from numpy import exp, sqrt, log
from abc import abstractmethod, ABC
from math import gamma

i = 1j
"""
目前需要查詢的部份
1. VG model的參數意義
2. CGMY參數意義
3. 其他process參數意義
"""


class CharacteristicFunction(ABC):
    @abstractmethod
    def getCFValue(self, u):
        pass


class GBM(CharacteristicFunction):
    """
    Black-Scholes Model
    """

    def __init__(self, r, sigma, T):
        self.r = r
        self.sigma = sigma
        self.T = T

    def getCFValue(self, u):
        cf_value = exp((self.r - 0.5 * self.sigma ** 2) * self.T * i * u -
                       0.5 * self.sigma ** 2 * u ** 2 * self.T)

        return cf_value


class Heston(CharacteristicFunction):
    """
    Heston's SV Model
    """

    def __init__(self, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, std_of_var_process):
        self.r = r
        self.sigma = sigma
        self.T = T
        self.mean_reversion_speed = mean_reversion_speed
        self.long_term_var_mean = long_term_var_mean
        self.corr = corr
        self.std_of_var_process = std_of_var_process

    def getCFValue(self, u):
        alpha = (self.mean_reversion_speed - i * self.corr * self.std_of_var_process * u) ** 2
        beta = (u ** 2 + i * u) * self.std_of_var_process ** 2
        b = sqrt(alpha + beta)
        a_minus = self.mean_reversion_speed - i * self.corr * self.std_of_var_process * u - b
        a_plus = self.mean_reversion_speed - i * self.corr * self.std_of_var_process * u + b
        a = a_minus / a_plus
        eta = 1 - a * exp(-b * self.T)
        C = i * u * self.r * self.T + self.mean_reversion_speed * self.long_term_var_mean / pow(self.std_of_var_process,
                                                                                                2) * (
                    self.T * a_minus - 2 * log(eta / (1 - a)))
        D = a_minus / pow(self.std_of_var_process, 2) * ((1 - exp(-b * self.T)) / eta)
        cf_value = exp(C + D * pow(self.sigma, 2))

        return cf_value


class SVJ(CharacteristicFunction):
    """
    Bates's Stochastic Volatility Jump Model
    """

    def __init__(self, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, std_of_var_process, jump_intensity,
                 jump_mean, jump_var):
        self.r = r
        self.sigma = sigma
        self.T = T
        self.mean_reversion_speed = mean_reversion_speed
        self.long_term_var_mean = long_term_var_mean
        self.corr = corr
        self.std_of_var_process = std_of_var_process
        self.jump_intensity = jump_intensity
        self.jump_mean = jump_mean
        self.jump_var = jump_var

    def getCFValue(self, u):
        k = exp(self.jump_mean + 0.5 * self.jump_var) - 1
        mu = self.r - self.jump_intensity * k

        alpha = (self.mean_reversion_speed - i * self.corr * self.std_of_var_process * u) ** 2
        beta = (u ** 2 + i * u) * self.std_of_var_process ** 2
        b = sqrt(alpha + beta)
        a_minus = self.mean_reversion_speed - i * self.corr * self.std_of_var_process * u - b
        a_plus = self.mean_reversion_speed - i * self.corr * self.std_of_var_process * u + b
        a = a_minus / a_plus
        eta = 1 - a * exp(-b * self.T)
        C = i * u * mu * self.T + self.mean_reversion_speed * self.long_term_var_mean / pow(self.std_of_var_process,
                                                                                            2) * (
                    a_minus * self.T - 2 * log(eta / (1 - a)))
        D = a_minus / pow(self.std_of_var_process, 2) * ((1 - exp(-b * self.T)) / eta)

        cf_value = exp(C + D * pow(self.sigma, 2) - self.jump_intensity * self.T * (
                1 - exp(i * u * self.jump_mean - 0.5 * pow(u, 2) * self.jump_var)))
        return cf_value


class MJD(CharacteristicFunction):
    """
    Merton's Normal Jump Diffusion Model
    """

    def __init__(self, r, sigma, T, jump_intensity, jump_mean, jump_var):
        self.r = r
        self.sigma = sigma
        self.T = T
        self.jump_intensity = jump_intensity
        self.jump_mean = jump_mean
        self.jump_var = jump_var

    def getCFValue(self, u):
        k = exp(self.jump_mean + 0.5 * self.jump_var) - 1
        mu = self.r - self.jump_intensity * k
        lnY_cf = exp(i * u * self.jump_mean - 0.5 * u ** 2 * self.jump_var)
        cf_value = exp(i * u * (mu - pow(self.sigma, 2) / 2) * self.T - 0.5 * pow(u * self.sigma,
                                                                                  2) * self.T - self.jump_intensity * self.T * (
                               1 - lnY_cf))
        return cf_value


# class BJD(CharacteristicFunction):
#     """
#     Amin's Bivariate Jump Diffusion Model
#     """
#
#     def __init__(self, p, delta):
#         self.p = p
#         self.delta = delta
#
#     def getCFValue(self, u):
#         cf_value = self.p * exp(i * u * self.delta) + (1 - self.p) * exp(-i * u * self.delta)
#         return cf_value
#
#
class KJD(CharacteristicFunction):
    """
    Kou's Double Exponential Jump Diffusion Model
    """

    def __init__(self, r, sigma, T, jump_intensity, p, eat1, eat2):
        self.r = r
        self.sigma = sigma
        self.T = T
        self.jump_intensity = jump_intensity
        self.p = p
        self.eta1 = eat1
        self.eta2 = eat2

    def getCFValue(self, u):
        k = self.p * self.eta1 / (self.eta1 - 1) + (1 - self.p) * self.eta2 / (self.eta2 + 1) - 1
        mu = self.r - self.jump_intensity * k
        lnY_cf = self.p * self.eta1 / (self.eta1 - i * u) + (1 - self.p) * self.eta2 / (self.eta2 + i * u)
        cf_value = exp(i * u * (mu - pow(self.sigma, 2) / 2) * self.T - 0.5 * pow(u * self.sigma,
                                                                                  2) * self.T - self.jump_intensity * self.T * (
                               1 - lnY_cf))

        return cf_value


class VG(CharacteristicFunction):
    """
    Madan, D.B., and Seneta, E. (1990). The Variance Gamma Model
    """

    def __init__(self, gamma_mean, gamma_var, sigma, T, r):
        self.gamma_mean = gamma_mean
        self.gamma_var = gamma_var
        self.sigma = sigma
        self.r = r
        self.T = T

    def Xt_cf(self, u):
        return pow(1 - i * u * self.gamma_mean * self.gamma_var + 0.5 * self.sigma ** 2 * self.gamma_var * u ** 2,
                   -self.T / self.gamma_var)

    def getCFValue(self, u):
        mu = self.r - log(self.Xt_cf(-i)) / self.T
        return_cf = exp(i * u * mu * self.T) * self.Xt_cf(u)

        return return_cf


class NIG(CharacteristicFunction):
    """
    Barndorff-Nielsen, O.L. (1997) Normal Inverse Gaussian Model
    """

    def __init__(self, r, T, delta, alpha, beta):
        self.r = r
        self.T = T
        self.delta = delta
        self.alpha = alpha
        self.beta = beta

    def Xt_cf(self, u):
        return exp(self.delta * self.T * (pow(self.alpha ** 2 - self.beta ** 2, 0.5) - pow(self.alpha ** 2 - (self.beta + i * u) ** 2, 0.5)))

    def getCFValue(self, u):
        mu = self.r - log(self.Xt_cf(-i)) / self.T
        return_cf = exp(i * u * mu * self.T) * self.Xt_cf(u)

        return return_cf

# class CGMY(CharacteristicFunction):
#     """
#     Carr, P., Geman, H., Madan, D., and Yor, M. (2002) CGMY Model
#     """
#
#     def __init__(self, r, T, sigma, C, G, M, Y):
#         self.r = r
#         self.T = T
#         self.sigma = sigma
#         self.C = C
#         self.G = G
#         self.M = M
#         self.Y = Y
#
#     def getCFValue(self, u):
#         A = i * u * self.r * self.T - 0.5 * u ** 2 * self.sigma ** 2 * self.T
#         B = self.T * self.C * gamma(-self.Y)
#         C = (self.M - i * u) ** self.Y - self.M ** self.Y + (self.G + i * u) ** self.Y - self.G ** self.Y
#         cf_value = exp(A + B * C)
#         return cf_value
