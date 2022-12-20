from PolynomialPricingMethod.COSMethod import PolyByCosMethod
from PolynomialPricingMethod.utils.CharactoristicFunc import *
from PolynomialPricingMethod.utils.DensityTools import DensityRecover
from math import inf
import math
from PolynomialPricingMethod.utils.Tools import timeit

# ===================== Test: GBM.txt Call ======================
# S0 = 100
# r = 0.03
# T = 1
# sigma = 0.25
#
# poly_coef = [-100, 1]
# positive_interval = [100, inf]
# process_cf = GBM(r=r, sigma=sigma, T=T)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef,positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[-1] = densityRecover.getIntergralRangeAndInterval()
# call = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                        positive_interval=best_positive_interval, N=256,
#                        lower_limit=lower_limit,
#                        upper_limit=upper_limit)
# print(call.getValue())

# ===================== Test: Heston Call ======================
# S0 = 100
# r = 0.05
# T = 0.5
# sigma = 0.2
#
# std_of_var_process = 0.1
# mean_reversion_speed = 3
# long_term_var_mean = 0.04
# corr = -0.1
#
# poly_coef = [-100, 1]
# positive_interval = [100, inf]
# process_cf = Heston(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
#                     long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef,positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[-1] = densityRecover.getIntergralRangeAndInterval()
#
# call = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                        positive_interval=best_positive_interval, N=256,
#                        lower_limit=lower_limit,
#                        upper_limit=upper_limit)
# print(call.getValue())

# ===================== Test: Heston Polynomial right up  ======================
# S0 = 50
# T = 0.5
# r = 0.1
# sigma = 0.4
#
# std_of_var_process = 0.1
# mean_reversion_speed = 3
# long_term_var_mean = 0.06
# corr = -0.1
#
# poly_coef = [-2725, -100, 1]
# positive_interval  = [5*(10+math.sqrt(209)), inf]
#
# process_cf = Heston(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
#                     long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef,positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[-1] = densityRecover.getIntergralRangeAndInterval()
#
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                        positive_interval=best_positive_interval, N=N,
#                        lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")

# ===================== Test: Heston Polynomial left up  ======================
# S0 = 110
# T = 0.5
# r = 0.03
# sigma = 0.4
#
# std_of_var_process = 0.1
# mean_reversion_speed = 3
# long_term_var_mean = 0.06
# corr = -0.1
#
# poly_coef = [947.1, -30.164, 0.309, -0.001]
# positive_interval = [0, 77, 82, 150]
#
#
# process_cf = Heston(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
#                     long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef,positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[-1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                        positive_interval=best_positive_interval, N=N,
#                        lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")

# ===================== Test: Heston Polynomial both up  ======================
# S0 = 15
# T = 1
# r = 0.1
# sigma = 0.3
#
# std_of_var_process = 0.1
# mean_reversion_speed = 3
# long_term_var_mean = 0.06
# corr = -0.1
#
# poly_coef = [44.235, -39.474, 5.4793, -0.2358, 0.0031]
# positive_interval  = [0, 1.363962, 10.620047, 25.599102, 38.481405, inf]
#
#
# process_cf = Heston(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
#                     long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef,positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[-1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                        positive_interval=best_positive_interval, N=N,
#                        lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")

# # ===================== Test: MJD Call ======================
# S0 = 100
# r = 0.05
# T = 0.5
# sigma = 0.2
#
# # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
# jump_intensity = 140
# jump_mean = 0.01
# jump_var = 0.02 ** 2
#
# poly_coef = [-100, 1]
# positive_interval = [100, inf]
# process_cf = MJD(r=r, sigma=sigma, T=T, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef,
#                                 positive_interval=positive_interval, N=1e5)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# call = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                        positive_interval=best_positive_interval, N=256,
#                        lower_limit=lower_limit,
#                        upper_limit=upper_limit)
# print(call.getValue())
# densityRecover.poltDensity()

# # ===================== Test: MJD Polynomial right up  ======================
# S0 = 50
# T = 0.5
# r = 0.1
# sigma = 0.4
#
# # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
# jump_intensity = 140
# jump_mean = 0.01
# jump_var = 0.02 ** 2
#
# poly_coef = [-2725, -100, 1]
# positive_interval = [5 * (10 + math.sqrt(209)), inf]
#
# process_cf = MJD(r=r, sigma=sigma, T=T, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef, positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                                  positive_interval=best_positive_interval, N=N,
#                                  lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")
#
# # ===================== Test: MJD Polynomial left up  ======================
# S0 = 100
# T = 0.5
# r = 0.03
# sigma = 0.25
#
# # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
# jump_intensity = 140
# jump_mean = 0.01
# jump_var = 0.02 ** 2
#
# poly_coef = [947100, -30164, 309, -1]
# positive_interval = [0, 77, 82, 150]
#
# process_cf = MJD(r=r, sigma=sigma, T=T, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef, positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                                  positive_interval=best_positive_interval, N=N,
#                                  lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")
#
# # ===================== Test: MJD Polynomial both up  ======================
# S0 = 5
# T = 1
# r = 0.1
# sigma = 0.25
#
# # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
# jump_intensity = 140
# jump_mean = 0.01
# jump_var = 0.02 ** 2
#
# poly_coef = [840, -638, 179, -22, 1]
# positive_interval = [0, 4, 5, 6, 7, inf]
#
# process_cf = MJD(r=r, sigma=sigma, T=T, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
#
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef, positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                                  positive_interval=best_positive_interval, N=N,
#                                  lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")




# # ===================== Test: SVJ Call ======================
# S0 = 100
# r = 0.05
# T = 0.5
# sigma = 0.2
#
# std_of_var_process = 0.1
# mean_reversion_speed = 3
# long_term_var_mean = 0.04
# corr = -0.1
#
# jump_intensity = 140
# jump_mean = 0.01
# jump_var = 0.02 ** 2
#
# poly_coeff = [-100, 1]
# positive_interval = [100, inf]
# process_cf = SVJ(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
#                  long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
#                  jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coeff=poly_coeff,
#                                 positive_interval=positive_interval, N=1e5)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# call = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coeff,
#                        positive_interval=best_positive_interval, N=256,
#                        lower_limit=lower_limit,
#                        upper_limit=upper_limit)
# print(call.getValue())
# densityRecover.poltDensity()

# # ===================== Test: SVJ Polynomial right up  ======================
# S0 = 50
# T = 0.5
# r = 0.1
# sigma = 0.4
#
# std_of_var_process = 0.1
# mean_reversion_speed = 3
# long_term_var_mean = 0.04
# corr = -0.1
#
# # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
# jump_intensity = 140
# jump_mean = 0.01
# jump_var = 0.02 ** 2
#
# poly_coef = [-2725, -100, 1]
# positive_interval = [5 * (10 + math.sqrt(209)), inf]
#
# process_cf = SVJ(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
#                  long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
#                  jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef, positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                                  positive_interval=best_positive_interval, N=N,
#                                  lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")
#
# # ===================== Test: SVJ Polynomial left up  ======================
# S0 = 100
# T = 0.5
# r = 0.03
# sigma = 0.25
#
# std_of_var_process = 0.1
# mean_reversion_speed = 3
# long_term_var_mean = 0.04
# corr = -0.1
#
# # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
# jump_intensity = 140
# jump_mean = 0.01
# jump_var = 0.02 ** 2
#
# poly_coef = [947100, -30164, 309, -1]
# positive_interval = [0, 77, 82, 150]
#
# process_cf = SVJ(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
#                  long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
#                  jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef, positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                                  positive_interval=best_positive_interval, N=N,
#                                  lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")

# # ===================== Test: SVJ Polynomial both up  ======================
# S0 = 5
# T = 1
# r = 0.1
# sigma = 0.25
#
# std_of_var_process = 0.1
# mean_reversion_speed = 3
# long_term_var_mean = 0.04
# corr = -0.1
#
# # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
# jump_intensity = 140
# jump_mean = 0.01
# jump_var = 0.02 ** 2
#
# poly_coef = [840, -638, 179, -22, 1]
# positive_interval = [0, 4, 5, 6, 7, inf]
#
# process_cf = SVJ(r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
#                  long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
#                  jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef, positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                                  positive_interval=best_positive_interval, N=N,
#                                  lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")


# ===================== Test: KJD Call  ======================
S0 = 100
r = 0.05
T = 0.5
sigma = 0.16

# 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095  Ans: 9.14732
jump_intensity = 1
p = 0.4
eta1 = 10
eta2 = 5

poly_coeff = [-98, 1]
positive_interval = [98, inf]

process_cf = KJD(r=r, sigma=sigma, T=T, jump_intensity=jump_intensity, p=p, eat1=eta1, eat2=eta2)
best_positive_interval = positive_interval.copy()
densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coeff=poly_coeff, positive_interval=positive_interval)
lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
    -1] = densityRecover.getIntergralRangeAndInterval()

for N in [64, 128, 256, 512, 1024]:
    polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coeff=poly_coeff,
                                 positive_interval=best_positive_interval, N=N,
                                 lower_limit=lower_limit, upper_limit=upper_limit)
    ans = timeit(polynomial.getValue)
    print(f"N={N}: value: {ans:.6f}")


# # ===================== Test: KJD right up  ======================
# S0 = 50
# T = 0.5
# r = 0.1
# sigma = 0.4
#
# # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095  Ans: 9.14732
# jump_intensity = 1
# p = 0.4
# eta1 = 10
# eta2 = 5
#
# poly_coef = [-2725, -100, 1]
# positive_interval = [5 * (10 + math.sqrt(209)), inf]
#
# process_cf = KJD(r=r, sigma=sigma, T=T, jump_intensity=jump_intensity, p=p, eat1=eta1, eat2=eta2)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef, positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                                  positive_interval=best_positive_interval, N=N,
#                                  lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")
#
#
# # ===================== Test: KJD left up  ======================
# S0 = 100
# T = 0.5
# r = 0.03
# sigma = 0.25
#
# # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095  Ans: 9.14732
# jump_intensity = 1
# p = 0.4
# eta1 = 10
# eta2 = 5
#
# poly_coef = [947100, -30164, 309, -1]
# positive_interval = [0, 77, 82, 150]
#
# process_cf = KJD(r=r, sigma=sigma, T=T, jump_intensity=jump_intensity, p=p, eat1=eta1, eat2=eta2)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef, positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                                  positive_interval=best_positive_interval, N=N,
#                                  lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")
#
#
# # ===================== Test: KJD both up  ======================
# S0 = 5
# T = 1
# r = 0.1
# sigma = 0.25
#
# # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095  Ans: 9.14732
# jump_intensity = 1
# p = 0.4
# eta1 = 10
# eta2 = 5
#
# poly_coef = [840, -638, 179, -22, 1]
# positive_interval = [0, 4, 5, 6, 7, inf]
#
# process_cf = KJD(r=r, sigma=sigma, T=T, jump_intensity=jump_intensity, p=p, eat1=eta1, eat2=eta2)
# best_positive_interval = positive_interval.copy()
# densityRecover = DensityRecover(S0=S0, process_cf=process_cf, poly_coef=poly_coef, positive_interval=positive_interval)
# lower_limit, upper_limit, best_positive_interval[0], best_positive_interval[
#     -1] = densityRecover.getIntergralRangeAndInterval()
#
# for N in [64, 128, 256, 512, 1024]:
#     polynomial = PolyByCosMethod(S0=S0, T=T, r=r, sigma=sigma, process_cf=process_cf, poly_coef=poly_coef,
#                                  positive_interval=best_positive_interval, N=N,
#                                  lower_limit=lower_limit, upper_limit=upper_limit)
#     ans = timeit(polynomial.getValue)
#     print(f"N={N}: value: {ans:.6f}")