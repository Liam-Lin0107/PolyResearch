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


class MJD(CharacteristicFunction):
    """
    Merton's Normal Jump Diffusion Model
    """
    def __init__(self, gamma, delta):
        self.gamma = gamma
        self.delta = delta

    def getCFValue(self, u):
        cf_value = exp(i * u * self.gamma - 0.5 * u ** 2 * self.delta ** 2)
        return cf_value

class BJD(CharacteristicFunction):
    """
    Amin's Bivariate Jump Diffusion Model
    """
    def __init__(self, p, delta):
        self.p = p
        self.delta = delta

    def getCFValue(self, u):
        cf_value = self.p * exp(i * u * self.delta) + (1 - self.p) * exp(-i * u * self.delta)
        return cf_value

class KJD(CharacteristicFunction):
    """
    Kou's Double Exponential Jump Diffusion Model
    """
    def __init__(self, p, eat1, eat2):
        self.p = p
        self.eta1 = eat1
        self.eta2 = eat2

    def getCFValue(self, u):
        cf_value = self.p * self.eta1 / (self.eta1 - i * u) + (1 - self.p) * self.eta2 / (self.eta2 + i * u)
        return cf_value

class MEM(CharacteristicFunction):
    """
    Kou’s Mixed-Exponential Jump Diffusion Model
    """
    def __init__(self):
        pass

class Heston(CharacteristicFunction):
    """
    Heston's SV Model
    """
    def __init__(self, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process):
        self.r = r
        self.sigma = sigma
        self.T = T
        self.mean_reversion_speed = mean_reversion_speed
        self.long_term_var_mean = long_term_var_mean
        self.corr = corr
        self.var_variance_process = var_variance_process


    def getCFValue(self, u):
        alpha = (self.mean_reversion_speed - i * self.corr * self.var_variance_process * u) ** 2
        beta = (u ** 2 + i * u) * self.var_variance_process ** 2
        D = sqrt(alpha + beta)
        g_minus = self.mean_reversion_speed - i * self.corr * self.var_variance_process * u - D
        g_plus = self.mean_reversion_speed - i * self.corr * self.var_variance_process * u + D
        G = g_minus / g_plus
        # 不確定這裡的mu是啥 我假設為riskfree rate
        eta = 1 - G * exp(-D * self.T)
        epslon = i * u * self.r * self.T + self.sigma / self.var_variance_process ** 2 * (
                1 - exp(-D * self.T) / eta) * g_minus
        theta = self.mean_reversion_speed * self.long_term_var_mean / self.var_variance_process ** 2 * (
                self.T * g_minus - 2 * log(eta / (1 - G)))
        cf_value = exp(epslon + theta)

        return cf_value

class SVJ(CharacteristicFunction):
    """
    Bates's Stochastic Volatility Jump Model
    """
    def __init__(self, T, r, sigma, mean_reversion_speed, long_term_var_mean, corr, var_variance_process, gamma, delta, jump_density):
        self.r = r
        self.sigma = sigma
        self.T = T
        self.mean_reversion_speed = mean_reversion_speed
        self.long_term_var_mean = long_term_var_mean
        self.corr = corr
        self.var_variance_process = var_variance_process
        self.jump_density = jump_density
        self.gamma = gamma
        self.delta = delta

    def getCFValue(self, u):
        k = exp(self.gamma + 0.5 * self.delta ** 2) - 1
        mu = self.r - self.jump_density * k

        alpha = (self.mean_reversion_speed - i * self.corr * self.var_variance_process * u) ** 2
        beta = (u ** 2 + i * u) * self.var_variance_process ** 2
        D = sqrt(alpha + beta)
        g_minus = self.mean_reversion_speed - i * self.corr * self.var_variance_process * u - D
        g_plus = self.mean_reversion_speed - i * self.corr * self.var_variance_process * u + D
        G = g_minus / g_plus
        # 不確定這裡的mu是啥 我假設為riskfree rate
        eta = 1 - G * exp(-D * self.T)
        epslon = i * u * mu * self.T + self.sigma / self.var_variance_process ** 2 * (
                1 - exp(-D * self.T) / eta) * g_minus
        theta = self.mean_reversion_speed * self.long_term_var_mean / self.var_variance_process ** 2 * (
                self.T * g_minus - 2 * log(eta / (1 - G)))

        cf_value = exp(epslon + theta - self.jump_density * self.T * (1 - exp(i * u * self.gamma - 0.5 * u ** 2 * self.delta ** 2)))
        return cf_value


class VG(CharacteristicFunction):
    """
    Madan, D.B., and Seneta, E. (1990). The Variance Gamma Model
    """
    def __init__(self, theta, nu, sigma, T):
        self.theta = theta
        self.nu = nu
        self.sigma= sigma
        self.T = T

    def getCFValue(self, u):
        inner = 1 - i * u * self.theta * self.nu + 0.5 * self.sigma ** 2 * self.nu * u ** 2
        cf_value = pow(inner, -self.T / self.nu)
        return cf_value

class NIG(CharacteristicFunction):
    """
    Barndorff-Nielsen, O.L. (1997) Normal Inverse Gaussian Model
    """
    def __init__(self, delta, T, alpha, beta):
        self.delta = delta
        self.T = T
        self.alpha = alpha
        self.beta = beta

    def getCFValue(self, u):
        inner = self.delta * self.T * (pow(pow(self.alpha, 2) - pow(self.beta, 2), 0.5) - pow(pow(self.alpha, 2) - pow(self.beta + i * u, 2), 2), 0.5)
        cf_value = exp(inner)
        return cf_value


class CGMY(CharacteristicFunction):
    """
    Carr, P., Geman, H., Madan, D., and Yor, M. (2002) CGMY Model
    """
    def __init__(self, r, T, sigma, C, G, M, Y):
        self.r = r
        self.T = T
        self.sigma = sigma
        self.C = C
        self.G = G
        self.M = M
        self.Y = Y

    def getCFValue(self, u):
        A = i * u * self.r * self.T - 0.5 * u ** 2 * self.sigma ** 2 * self.T
        B = self.T * self.C * gamma(-self.Y)
        C = (self.M - i * u) ** self.Y - self.M ** self.Y + (self.G + i * u) ** self.Y - self.G ** self.Y
        cf_value = exp(A + B * C)
        return cf_value