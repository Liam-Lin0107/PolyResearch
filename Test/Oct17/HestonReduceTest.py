from PolynomialPricingMethod.COSMethod import PolyByCosMethod
from PolynomialPricingMethod.utils.CharactoristicFunc import *
from PolynomialPricingMethod.utils.DensityTools import DensityRecover
from PolynomialPricingMethod.utils.Tools import timeit
from math import inf


# cos method
S0 = 100
T = 1
r = 0.03
sigma = 0.25 # standard deviation
poly_coef = [-100, 1]
positive_interval = [100, inf]
mean_reversion_speed = 1
long_term_var_mean = 0.25 ** 2 # variance
corr = -0.5711
var_variance_process = 1e-10

process = Heston(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                 long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=var_variance_process)
densityRecover = DensityRecover(S0=S0, process_cf=process, poly_coeff=poly_coef, positive_interval=positive_interval, step=0.1)
lower_limit, upper_limit, positive_interval[0], positive_interval[-1] = densityRecover.getBestIntervalPoints()

N = 256
polynomialCosMethod = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process,
                                      poly_coeff=poly_coef, positive_interval=positive_interval,
                                      N=N, lower_limit=lower_limit, upper_limit=upper_limit)
ans = timeit(polynomialCosMethod.getValue)
print(f"N={N}: value: {ans:.6f}")

# # MonteCarlo
# polynomialHestonByMC= PolynomialHestonByMC(S0=S0, T=T, r=r, sigma=sigma, poly_coef=poly_coef,
#                      mean_reversion_speed=mean_reversion_speed, long_term_var_mean=long_term_var_mean,
#                      corr=corr, var_variance_process=var_variance_process)
# print("C.I:", polynomialHestonByMC.getStastic(plot=True, plot_save=True,file_name="MC_Test01"))
# print("")
