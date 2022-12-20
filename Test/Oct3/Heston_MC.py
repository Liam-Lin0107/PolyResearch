# import sys
# sys.path.append("D:\\Users\\DZLin\\Oplynomial_Option_Reserch-master")

from PolynomialPricingMethod.PathSimulationMethod import HestonByMC

# 左翹
S0 = 100
T = 0.5
r = 0.03
sigma = 0.25
poly_coef = [947100, -30164, 309, -1]
mean_reversion_speed = 1.5768
long_term_var_mean = 0.0398
corr = -0.5711
var_variance_process = 0.5751


polynomialHestonByMC= HestonByMC(S0=S0, T=T, r=r, sigma=sigma, poly_coef=poly_coef,
                                 mean_reversion_speed=mean_reversion_speed, long_term_var_mean=long_term_var_mean,
                                 corr=corr, std_of_var_process=var_variance_process)
print("C.I:", polynomialHestonByMC.getStatistic(plot=True, plot_save=True, file_name="MC_Test01"))
print("")
#
# # 雙翹
# S0 = 5
# T = 1
# r = 0.1
# sigma = 0.25
# poly_coef = [840, -638, 179, -22, 1]
# mean_reversion_speed = 1
# long_term_var_mean = 0.25
# corr = -0.5711
# var_variance_process = 1e-8
#
# polynomialHestonByMC= PolynomialHestonByMC(S0=S0, T=T, r=r, sigma=sigma, poly_coef=poly_coef,
#                      mean_reversion_speed=mean_reversion_speed, long_term_var_mean=long_term_var_mean,
#                      corr=corr, var_variance_process=var_variance_process)
# print("C.I:", polynomialHestonByMC.getStastic(plot=True, plot_save=True,file_name="MC_Test02"))
# print("")

#
# # 單邊右翹
# S0 = 50
# T = 0.5
# r = 0.1
# sigma = 0.4
# poly_coef = [-2725, -100, 1]
# mean_reversion_speed = 1.5768
# long_term_var_mean = 0.0398
# corr = -0.5711
# var_variance_process = 0.5751
#
# polynomialHestonByMC= PolynomialHestonByMC(S0=S0, T=T, r=r, sigma=sigma, poly_coef=poly_coef,
#                      mean_reversion_speed=mean_reversion_speed, long_term_var_mean=long_term_var_mean,
#                      corr=corr, var_variance_process=var_variance_process)
# print("C.I:", polynomialHestonByMC.getStastic(plot=True, plot_save=True,file_name="MC_Test03"))
# print("")