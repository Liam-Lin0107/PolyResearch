# import sys
# sys.path.append("D:\\Users\\DZLin\\Oplynomial_Option_Reserch-master")

from PolynomialPricingMethod.COSMethod import PolyByCosMethod
from PolynomialPricingMethod.utils.CharactoristicFunc import *
from PolynomialPricingMethod.utils.DensityTools import DensityRecover
from PolynomialPricingMethod.utils.Tools import timeit

# 左翹
S0 = 100
T = 0.5
r = 0.03
sigma = 0.25
poly_coef = [947100, -30164, 309, -1]
positive_interval = [0, 77, 82, 150]
mean_reversion_speed = 1
long_term_var_mean = 0.25
corr = -0.5711
var_variance_process = 1e-8

process = Heston(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                 long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=var_variance_process)
densityRecover = DensityRecover(S0=S0, process_cf=process, poly_coeff=poly_coef, positive_interval=positive_interval, step=0.05)
lower_limit, upper_limit, positive_interval[0], positive_interval[-1] = densityRecover.getBestIntervalPoints()
print("best integral lower bound:", lower_limit)
print("best integral upper bound:", upper_limit)
print("best interval lower bound:", positive_interval[0])
print("best interval upper bound:", positive_interval[-1])
for N in [64, 128, 256, 512, 1024]:
    polynomialCosMethod = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process,
                                          poly_coeff=poly_coef, positive_interval=positive_interval,
                                          N=N, lower_limit=lower_limit, upper_limit=upper_limit)
    ans = timeit(polynomialCosMethod.getValue)
    print(f"N={N}: value: {ans:.6f}")
densityRecover.plotFittingDetal(plot_save=True, file_name_prefix="Test01")
print("")

# 雙翹
S0 = 5
T = 1
r = 0.1
sigma = 0.25
poly_coef = [840, -638, 179, -22, 1]
positive_interval  = [0, 4, 5, 6, 7, inf]
mean_reversion_speed = 1.5768
long_term_var_mean = 0.0398
corr = -0.5711
var_variance_process = 0.5751

process = Heston(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                 long_term_var_mean=long_term_var_mean, corr=corr, var_variance_process=var_variance_process)
densityRecover = DensityRecover(S0=S0, process=process, poly_coeff=poly_coef, positive_interval=positive_interval, step=0.01)
lower_limit, upper_limit, positive_interval[0], positive_interval[-1] = densityRecover.getBestIntervalPoints()
print("best integral lower bound:", lower_limit)
print("best integral upper bound:", upper_limit)
print("best interval lower bound:", positive_interval[0])
print("best interval upper bound:", positive_interval[-1])
for N in [64, 128, 256, 512, 1024]:
    polynomialCosMethod = PolynomialCosMethod(S0=S0, T=T, r=r, sigma=sigma, process=process,
                        poly_coef=poly_coef,positive_interval=positive_interval,
                        N=N, lower_limit=lower_limit, upper_limit=upper_limit)
    ans = timeit(polynomialCosMethod.getValue)
    print(f"N={N}: value: {ans:.6f}")
densityRecover.plotFittingDetal(plot_save=True, file_name_prefix="Test02")
print("")


# 單邊右翹
S0 = 50
T = 0.5
r = 0.1
sigma = 0.4
poly_coef = [-2725, -100, 1]
positive_interval  = [5*(10+math.sqrt(209)), inf]
mean_reversion_speed = 1.5768
long_term_var_mean = 0.0398
corr = -0.5711
var_variance_process = 0.5751

process = Heston(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                 long_term_var_mean=long_term_var_mean, corr=corr, var_variance_process=var_variance_process)
densityRecover = DensityRecover(S0=S0, process=process, poly_coeff=poly_coef, positive_interval=positive_interval, step=0.01)
lower_limit, upper_limit, positive_interval[0], positive_interval[-1] = densityRecover.getBestIntervalPoints()
print("best integral lower bound:", lower_limit)
print("best integral upper bound:", upper_limit)
print("best interval lower bound:", positive_interval[0])
print("best interval upper bound:", positive_interval[-1])
for N in [64, 128, 256, 512, 1024]:
    polynomialCosMethod = PolynomialCosMethod(S0=S0, T=T, r=r, sigma=sigma, process=process,
                        poly_coef=poly_coef,positive_interval=positive_interval,
                        N=N, lower_limit=lower_limit, upper_limit=upper_limit)
    ans = timeit(polynomialCosMethod.getValue)
    print(f"N={N}: value: {ans:.6f}")
densityRecover.plotFittingDetal(plot_save=True, file_name_prefix="Test03")
print("")