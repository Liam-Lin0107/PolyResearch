import sys

my_pc = False
if not my_pc:
    sys.path.append("D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval")  # for 85
import threading
from PolynomialPricingMethod.PathSimulationMethod import *
import time
n = 5000
def SVJTest1():

    # ===================== Test: SVJ Polynomial left up  ======================
    print("SVJ left up 開始模擬...")
    S0 = 110
    T = 0.25
    r = 0.03
    sigma = 0.4

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.04
    corr = -0.1

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coeff = [947.1, -30.164, 0.309, -0.001]

    process = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                      long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                      poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n)

    print(process.getStatistic(save_data=True, save_dir=r"D:\Users\DZLin\Oplynomial_Option_Reserch-newFindInterval\Data\Simulation\SVJ.txt"))
def SVJTest2():
    # ===================== Test: SVJ Polynomial both up  ======================
    print("SVJ both up 開始模擬...")
    S0 = 15
    T = 0.25
    r = 0.1
    sigma = 0.3

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.04
    corr = -0.1

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coeff = [44.235, -39.474, 5.4793, -0.2358, 0.0031]

    process = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                      long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                      poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n)

    print(process.getStatistic(save_data=True, save_dir=r"D:\Users\DZLin\Oplynomial_Option_Reserch-newFindInterval\Data\Simulation\SVJ.txt"))


if __name__ == "__main__":
    

    thread_SVJ1 = threading.Thread(target=SVJTest1, name="SVJ1")
    thread_SVJ2 = threading.Thread(target=SVJTest2, name="SVJ2")
    thread_SVJ1.start()
    thread_SVJ2.start()
   