import sys
import time

my_pc = False
if my_pc:
    sys.path.append(
        r"/")  # for my computer
else:
    sys.path.append("D:\\Users\\DZLin\\Oplynomial_Option_Reserch-newFindInterval")  # for 85
import threading
from PolynomialPricingMethod.PathSimulationMethod import *


def HestonTest():
    # ===================== Test: Heston Call ======================
    print("Heston Call 開始模擬...")
    S0 = 100
    r = 0.05
    T = 0.5
    sigma = 0.2

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.04
    corr = -0.1

    poly_coef = [-100, 1]

    heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coef=poly_coef)
    if not my_pc:
        print(heston.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                           "\\PolynomialPricingMethod\\tests\\Data\\HestonSimulationTest"
                                                           ".txt"))
    else:
        print(heston.getStatistic(save_data=True,
                                  save_dir=r"/PolynomialPricingMethod/tests/Data/HestonSimulationTest.txt"))

    # ===================== Test: Heston Polynomial right up  ======================
    print("Heston right up 開始模擬...")
    S0 = 50
    T = 0.5
    r = 0.1
    sigma = 0.4

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.06
    corr = -0.1

    poly_coef = [-2725, -100, 1]

    heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coef=poly_coef)

    if not my_pc:
        print(heston.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                           "\\PolynomialPricingMethod\\tests\\Data\\HestonSimulationTest"
                                                           ".txt"))
    else:
        print(heston.getStatistic(save_data=True,
                                  save_dir=r"/PolynomialPricingMethod/tests/Data/HestonSimulationTest.txt"))

    # ===================== Test: Heston Polynomial left up  ======================
    print("Heston left up 開始模擬...")
    S0 = 100
    T = 0.5
    r = 0.03
    sigma = 0.25

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.06
    corr = -0.1

    poly_coef = [947100, -30164, 309, -1]
    positive_interval = [0, 77, 82, 150]

    heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coef=poly_coef)

    if not my_pc:
        print(heston.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                           "\\PolynomialPricingMethod\\tests\\Data\\HestonSimulationTest"
                                                           ".txt"))
    else:
        print(heston.getStatistic(save_data=True,
                                  save_dir=r"/PolynomialPricingMethod/tests/Data/HestonSimulationTest.txt"))

    # ===================== Test: Heston Polynomial both up  ======================
    print("Heston both up 開始模擬...")

    S0 = 5
    T = 1
    r = 0.1
    sigma = 0.25

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.06
    corr = -0.1

    poly_coef = [840, -638, 179, -22, 1]

    heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coef=poly_coef)

    if not my_pc:
        print(heston.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                           "\\PolynomialPricingMethod\\tests\\Data\\HestonSimulationTest"
                                                           ".txt"))
    else:
        print(heston.getStatistic(save_data=True,
                                  save_dir=r"/PolynomialPricingMethod/tests/Data/HestonSimulationTest.txt"))


def MertonTest():
    # ===================== Test: Merton Call ======================
    print("Merton Call 開始模擬...")

    S0 = 100
    r = 0.05
    T = 0.5
    sigma = 0.2

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coef = [-100, 1]

    process = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coef=poly_coef, jump_intensity=jump_intensity,
                      jump_mean=jump_mean, jump_var=jump_var)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                          "\\PolynomialPricingMethod\\tests\\Data\\MJDSimulationTest.txt"))
    else:
        print(process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/MJDSimulationTest.txt"))

    # # ===================== Test: Merton Polynomial right up  ======================
    print("Merton right up 開始模擬...")
    S0 = 50
    T = 0.5
    r = 0.1
    sigma = 0.4

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coef = [-2725, -100, 1]

    process = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coef=poly_coef, jump_intensity=jump_intensity,
                      jump_mean=jump_mean, jump_var=jump_var)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                          "\\PolynomialPricingMethod\\tests\\Data\\MJDSimulationTest.txt"))
    else:
        print(
            process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/MJDSimulationTest.txt"))

    # # ===================== Test: Merton Polynomial left up  ======================
    print("Merton left up 開始模擬...")

    S0 = 100
    T = 0.5
    r = 0.03
    sigma = 0.25

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coef = [947100, -30164, 309, -1]

    process = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coef=poly_coef, jump_intensity=jump_intensity,
                      jump_mean=jump_mean, jump_var=jump_var)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                          "\\PolynomialPricingMethod\\tests\\Data\\MJDSimulationTest.txt"))
    else:
        print(
            process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/MJDSimulationTest.txt"))

    # # ===================== Test: Merton Polynomial both up  ======================
    print("Merton both up 開始模擬...")

    S0 = 5
    T = 1
    r = 0.1
    sigma = 0.25

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coef = [840, -638, 179, -22, 1]

    process = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coef=poly_coef, jump_intensity=jump_intensity,
                      jump_mean=jump_mean, jump_var=jump_var)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                          "\\PolynomialPricingMethod\\tests\\Data\\MJDSimulationTest.txt"))
    else:
        print(
            process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/MJDSimulationTest.txt"))


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

    poly_coef = [-98, 1]

    process = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coef=poly_coef, jump_intensity=jump_intensity, p=p, eta1=eta1,
                      eta2=eta2)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                            "\\PolynomialPricingMethod\\tests\\Data\\KJDSimulationTest.txt"))
    else:
        print(process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/KJDSimulationTest.txt"))

    # ===================== Test: KDJ Polynomial right up  ======================
    print("KDJ right up 開始模擬...")

    S0 = 50
    T = 0.5
    r = 0.1
    sigma = 0.4

    # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095
    jump_intensity = 1
    p = 0.4
    eta1 = 10
    eta2 = 5

    poly_coef = [-2725, -100, 1]

    process = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, jump_intensity=jump_intensity,
                      p=p, eta1=eta1, eta2=eta2, poly_coef=poly_coef)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                            "\\PolynomialPricingMethod\\tests\\Data\\KJDSimulationTest.txt"))
    else:
        print(process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/KJDSimulationTest.txt"))

    # ===================== Test: KJD Polynomial left up  ======================
    print("KDJ left up 開始模擬...")

    S0 = 100
    T = 0.5
    r = 0.03
    sigma = 0.25

    # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095
    jump_intensity = 1
    p = 0.4
    eta1 = 10
    eta2 = 5

    poly_coef = [947100, -30164, 309, -1]

    process = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, jump_intensity=jump_intensity,
                      p=p, eta1=eta1, eta2=eta2, poly_coef=poly_coef)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                            "\\PolynomialPricingMethod\\tests\\Data\\KJDSimulationTest.txt"))
    else:
        print(process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/KJDSimulationTest.txt"))

    # ===================== Test: KJD Polynomial both up  ======================
    print("KDJ both up 開始模擬...")

    S0 = 5
    T = 1
    r = 0.1
    sigma = 0.25

    # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095
    jump_intensity = 1
    p = 0.4
    eta1 = 10
    eta2 = 5

    poly_coef = [840, -638, 179, -22, 1]

    process = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, jump_intensity=jump_intensity,
                      p=p, eta1=eta1, eta2=eta2, poly_coef=poly_coef)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                            "\\PolynomialPricingMethod\\tests\\Data\\KJDSimulationTest.txt"))
    else:
        print(process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/KJDSimulationTest.txt"))


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

    poly_coef = [-100, 1]

    process = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                      long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                      poly_coef=poly_coef, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                            "\\PolynomialPricingMethod\\tests\\Data\\SVJSimulationTest.txt"))
    else:
        print(process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/SVJSimulationTest.txt"))

    # ===================== Test: SVJ Polynomial right up  ======================
    print("SVJ right up 開始模擬...")
    S0 = 50
    T = 0.5
    r = 0.1
    sigma = 0.4

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.04
    corr = -0.1

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coef = [-2725, -100, 1]

    process = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                      long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                      poly_coef=poly_coef, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                            "\\PolynomialPricingMethod\\tests\\Data\\SVJSimulationTest.txt"))
    else:
        print(process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/SVJSimulationTest.txt"))
    # ===================== Test: SVJ Polynomial left up  ======================
    print("SVJ left up 開始模擬...")

    S0 = 100
    T = 0.5
    r = 0.03
    sigma = 0.25

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.04
    corr = -0.1

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coef = [947100, -30164, 309, -1]

    process = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                      long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                      poly_coef=poly_coef, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                            "\\PolynomialPricingMethod\\tests\\Data\\SVJSimulationTest.txt"))
    else:
        print(process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/SVJSimulationTest.txt"))

    # ===================== Test: SVJ Polynomial both up  ======================
    print("SVJ both up 開始模擬...")

    S0 = 5
    T = 1
    r = 0.1
    sigma = 0.25

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.04
    corr = -0.1

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coef = [840, -638, 179, -22, 1]

    process = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                      long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                      poly_coef=poly_coef, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var)
    if not my_pc:
        print(process.getStatistic(save_data=True, save_dir="D:\\Users\DZLin\\Oplynomial_Option_Reserch-newFindInterval"
                                                            "\\PolynomialPricingMethod\\tests\\Data\\SVJSimulationTest.txt"))
    else:
        print(process.getStatistic(save_data=True, save_dir=r"/PolynomialPricingMethod/tests/Data/SVJSimulationTest.txt"))


if __name__ == "__main__":
    # thread_heston = threading.Thread(target=HestonTest, name="Heston")
    # thread_MJD = threading.Thread(target=MertonTest, name="MJD")
    # thread_KJD = threading.Thread(target=KJDTest, name="KJD")
    # thread_SVJ = threading.Thread(target=SVJTest, name="SVJ")
    # thread_heston.start()
    # time.sleep(0.2)
    # thread_MJD.start()
    # time.sleep(0.2)
    # thread_KJD.start()
    # time.sleep(0.2)
    # thread_SVJ.start()
    HestonTest()
    MertonTest()
    KJDTest() 
    SVJTest()
