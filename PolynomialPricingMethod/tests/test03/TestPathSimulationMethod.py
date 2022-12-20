import sys
import time

my_pc = True
if my_pc:
    sys.path.append(
        r"/")  # for my computer
else:
    sys.path.append("D:\\Users\\DZLin\\Oplynomial_Option_Reserch-newFindInterval")  # for 85
import threading
from PolynomialPricingMethod.PathSimulationMethod import *


def GBMTest():
    # ===================== Test: Heston Call ======================
    print("GBM Call 開始模擬...")
    S0 = 100
    r = 0.03
    T = 1
    sigma = 0.25

    poly_coeff = [-100, 1]

    GBM = GBMByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff,
                  N_repeat=20, N_line=100000, n=252)

    print(GBM.getStatistic(save_data=True, save_dir=r"../Data/GBMSimulation"))


def HestonTest():
    S0 = 100
    r = 0.05
    T = 0.5
    sigma = 0.2

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.04
    corr = -0.1

    poly_coeff = [-100, 1]

    heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coeff=poly_coeff, N_repeat=20, N_line=10000, n=252)
    print(heston.getStatistic(save_data=True, save_dir=r"../Data/HestonSimulationTest.txt"))


def MJDTest():
    print("Merton Call 開始模擬...")

    S0 = 100
    r = 0.05
    T = 0.5
    sigma = 0.2

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coeff = [-100, 1]

    MJD = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
                  jump_mean=jump_mean, jump_var=jump_var, N_repeat=20, N_line=10000, n=2000)

    print(MJD.getStatistic(save_data=True, save_dir=r"../Data/MJDSimulationTest.txt"))


def SVJTest():
    # ===================== Test: SVJ Call ======================
    print("SVJ Call 開始模擬...")

    S0 = 100
    r = 0.05
    T = 0.5
    sigma = 0.2

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.04
    corr = -0.1

    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coeff = [-100, 1]

    process = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                      long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                      poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var,
                      N_repeat=20, N_line=10000, n=2000)

    print(process.getStatistic(save_data=True, save_dir="../Data/SVJSimulationTest.txt"))


def KJDTest():
    # ===================== Test: KDJ Call ======================
    print("KDJ Call 開始模擬...")

    S0 = 100
    r = 0.05
    T = 0.5
    sigma = 0.16

    # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095
    jump_intensity = 1
    p = 0.4
    eta1 = 10
    eta2 = 5

    poly_coeff = [-98, 1]

    process = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity, p=p,
                      eta1=eta1, eta2=eta2, N_repeat=20, N_line=10000, n=2000)

    print(process.getStatistic(save_data=True, save_dir="../Data/KJDSimulationTest.txt"))
