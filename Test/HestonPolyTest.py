from PolynomialPricingMethod.COSMethod import PolyByCosMethod
from PolynomialPricingMethod.utils.CharactoristicFunc import *
from PolynomialPricingMethod.PathSimulationMethod import HestonByMC
from math import exp

print("=" * 10)
print("Test: 1")
print("=" * 10)

S0 = 100
T = 0.5
r = 0.03
sigma = 0.25
poly_coef = [947100, -30164, 309, -1]
positive_interval = [exp(-9), 77, 82, 150]
lower_limit = exp(-9)
upper_limit = exp(400)
N = 10000

mean_reversion_speed = 1.5768
long_term_var_mean = 0.0398
var_variance_process = 0.0751
corr = -0.5711

process = Heston(r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process)
cosMethod = PolyByCosMethod(S0, T, r, sigma, process, poly_coef, positive_interval, N, lower_limit, upper_limit)
print("COSMethod: ", cosMethod.getValue())

polynomialHestonByMC = HestonByMC(S0, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process,
                                  poly_coef, N_line=10000, n=252, N_repeat=20)
a, b = polynomialHestonByMC.getStatistic()
print("MC: ", "[", a, ",", b, "]")

# print("=" * 10)
# print("Test: 1")
# print("=" * 10)
#
# S0 = 100
# T = 0.5
# r = 0.03
# sigma = 0.25
# poly_coef = [947100, -30164, 309, -1]
# positive_interval = [exp(-9), 77, 82, 150]
# lower_limit = exp(-9)
# upper_limit = exp(400)
# N = 10000
#
# mean_reversion_speed=1.5768
# long_term_var_mean=0.0398
# var_variance_process=0.0751
# corr=-0.5711
#
# process = Heston(r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process)
# cosMethod= PolynomialCosMethod(S0, T, r, sigma, process, poly_coef, positive_interval, N, lower_limit, upper_limit)
# print("COSMethod: ",cosMethod.getValue())
#
# polynomialHestonByMC = PolynomialHestonByMC(S0, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process, poly_coef, N_line=10000, n=252, N_repeat=20)
# a, b = polynomialHestonByMC.getStastic()
# print("MC: ","[", a,",", b, "]")
#
#
# print("=" * 10)
# print("Test: 1")
# print("=" * 10)
#
# S0 = 100
# T = 0.5
# r = 0.03
# sigma = 0.25
# poly_coef = [947100, -30164, 309, -1]
# positive_interval = [exp(-9), 77, 82, 150]
# lower_limit = exp(-9)
# upper_limit = exp(400)
# N = 10000
#
# mean_reversion_speed=1.5768
# long_term_var_mean=0.0398
# var_variance_process=0.0751
# corr=-0.5711
#
# process = Heston(r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process)
# cosMethod= PolynomialCosMethod(S0, T, r, sigma, process, poly_coef, positive_interval, N, lower_limit, upper_limit)
# print("COSMethod: ",cosMethod.getValue())
#
# polynomialHestonByMC = PolynomialHestonByMC(S0, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process, poly_coef, N_line=10000, n=252, N_repeat=20)
# a, b = polynomialHestonByMC.getStastic()
# print("MC: ","[", a,",", b, "]")
#
#
# print("=" * 10)
# print("Test: 1")
# print("=" * 10)
#
# S0 = 100
# T = 0.5
# r = 0.03
# sigma = 0.25
# poly_coef = [947100, -30164, 309, -1]
# positive_interval = [exp(-9), 77, 82, 150]
# lower_limit = exp(-9)
# upper_limit = exp(400)
# N = 10000

mean_reversion_speed = 1.5768
long_term_var_mean = 0.0398
var_variance_process = 0.0751
corr = -0.5711

process = Heston(r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process)
cosMethod = PolyByCosMethod(S0, T, r, sigma, process, poly_coef, positive_interval, N, lower_limit, upper_limit)
print("COSMethod: ", cosMethod.getValue())

polynomialHestonByMC = HestonByMC(S0, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process,
                                  poly_coef, N_line=10000, n=252, N_repeat=20)
a, b = polynomialHestonByMC.getStatistic()
print("MC: ", "[", a, ",", b, "]")

print("=" * 10)
print("Test: 1")
print("=" * 10)

S0 = 100
T = 0.5
r = 0.03
sigma = 0.25
poly_coef = [947100, -30164, 309, -1]
positive_interval = [exp(-9), 77, 82, 150]
lower_limit = exp(-9)
upper_limit = exp(400)
N = 10000

mean_reversion_speed = 1.5768
long_term_var_mean = 0.0398
var_variance_process = 0.0751
corr = -0.5711

process = Heston(r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process)
cosMethod = PolyByCosMethod(S0, T, r, sigma, process, poly_coef, positive_interval, N, lower_limit, upper_limit)
print("COSMethod: ", cosMethod.getValue())

polynomialHestonByMC = HestonByMC(S0, r, sigma, T, mean_reversion_speed, long_term_var_mean, corr, var_variance_process,
                                  poly_coef, N_line=10000, n=252, N_repeat=20)
a, b = polynomialHestonByMC.getStatistic()
print("MC: ", "[", a, ",", b, "]")
