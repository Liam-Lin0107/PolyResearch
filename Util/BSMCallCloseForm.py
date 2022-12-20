import numpy as np
from numpy import log, sqrt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')
class BSMCloseForm:
    def __init__(self, S0, T, r, sigma, K):
        self.S0 = S0
        self.T = T
        self.r = r
        self.sigma = sigma
        self.K = K
    def BSM_call_value(self):
        d1 = (log(self.S0 / self.K) + (self.r + 0.5 * self.sigma ** 2) * self.T) \
            / (self.sigma * sqrt(self.T))
        d2 = (log(self.S0 / self.K) + (self.r - 0.5 * self.sigma ** 2) * self.T) \
            / (self.sigma * sqrt(self.T))
        BS_C = (self.S0 * stats.norm.cdf(d1, 0.0, 1.0) -
                self.K * np.exp(-self.r * self.T) * stats.norm.cdf(d2, 0.0, 1.0))
        return BS_C
if __name__ == "__main__":
    call = BSMCloseForm(S0=5, T=1, r=0.1, sigma=0.25, K=100)
