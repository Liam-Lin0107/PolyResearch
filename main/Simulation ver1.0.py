import sys
# 放入文件路徑
sys.path.append(
    "/Users/lindazhong/Downloads/PolyOptionPricingRearch/Program")  # for 85
import threading
# 導入自己寫的蒙地卡羅
from PolynomialPricingMethod.PathSimulationMethod import *
import time
# 設定蒙地卡羅次數
n = 252  # for setting dt
N = int(1e6)  # N_line
N_repeat = 20
save_directory = "./Data/Simulation/"



def GBMTest():
    # ===================== Test: GBM Call  ======================
    print("GBM Call 開始模擬...")
    # 基本資訊
    S0 = 100
    r = 0.05
    T = 0.25
    sigma = 0.2
    # option payoff係數
    poly_coeff = [-100, 1]
    # 建立一個GMBByMC對象並賦值
    GBM = GBMByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat, n=n)
    # 調用GMBByMC的getStatistic方法會回傳平均數與變異數，並print出來
    # 裡面有save_data=True表示會產生一個txt文件存放數據
    print(GBM.getStatistic(save_data=True, save_dir=save_directory +  "GBM.txt"))
    # ===================== Test: GBM Polynomial right up  ======================
    print("GBM right up 開始模擬...")
    S0 = 80
    T = 0.25
    r = 0.1
    sigma = 0.4

    poly_coeff = [-2725, -100, 1]

    GBM = GBMByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat, n=n)
    print(GBM.getStatistic(save_data=True, save_dir=save_directory +  "GBM.txt"))

    # ===================== Test: GBM Polynomial left up  ======================
    print("GBM left up 開始模擬...")
    S0 = 110
    T = 0.25
    r = 0.03
    sigma = 0.4

    poly_coeff = [947.1, -30.164, 0.309, -0.001]

    GBM = GBMByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat, n=n)
    print(GBM.getStatistic(save_data=True, save_dir=save_directory + "GBM.txt"))

    # ===================== Test: GBM Polynomial both up  ======================
    print("GBM both up 開始模擬...")
    S0 = 15
    T = 0.25
    r = 0.1
    sigma = 0.3

    poly_coeff = [44.235, -39.474, 5.4793, -0.2358, 0.0031]

    GBM = GBMByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat, n=n)
    print(GBM.getStatistic(save_data=True, save_dir=save_directory +  "GBM.txt"))

    # # ===================== Test: GBM Polynomial both up[老師版]  ======================
    # print("GBM both up 開始模擬...")
    # S0 = 5
    # T = 1
    # r = 0.1
    # sigma = 0.25
    #
    # poly_coeff = [840, -638, 179, -22, 1]
    #
    # GBM = GBMByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, N_line=N)
    # print(GBM.getStatistic(save_data=True, save_dir=save_directory +  "GBM.txt"))


def HestonTest():
    # ===================== Test: Heston Call ======================
    print("Heston Call 開始模擬...")
    S0 = 100
    r = 0.05
    T = 0.25
    sigma = 0.2

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.04
    corr = -0.1

    poly_coeff = [-100, 1]

    heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat, n=n)

    print(heston.getStatistic(save_data=True, save_dir=save_directory +"Heston.txt"))

    # ===================== Test: Heston Polynomial right up  ======================
    print("Heston right up 開始模擬...")
    S0 = 80
    T = 0.25
    r = 0.1
    sigma = 0.4

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.06
    corr = -0.1

    poly_coeff = [-2725, -100, 1]

    heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat, n=n)

    print(heston.getStatistic(save_data=True, save_dir=save_directory + "Heston.txt"))

    # ===================== Test: Heston Polynomial left up  ======================
    print("Heston left up 開始模擬...")
    S0 = 110
    T = 0.25
    r = 0.03
    sigma = 0.4

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.06
    corr = -0.1

    poly_coeff = [947.1, -30.164, 0.309, -0.001]

    heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat, n=n)

    print(heston.getStatistic(save_data=True, save_dir=save_directory +"Heston.txt"))

    # ===================== Test: Heston Polynomial both up  ======================
    print("Heston both up 開始模擬...")
    S0 = 15
    T = 0.25
    r = 0.1
    sigma = 0.3

    std_of_var_process = 0.1
    mean_reversion_speed = 3
    long_term_var_mean = 0.06
    corr = -0.1

    poly_coeff = [44.235, -39.474, 5.4793, -0.2358, 0.0031]

    heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat, n=n)

    print(heston.getStatistic(save_data=True, save_dir=save_directory +"Heston.txt"))

    # # ===================== Test: Heston Polynomial both up[老師版]  ======================
    # print("Heston both up 開始模擬...")
    # S0 = 5
    # T = 1
    # r = 0.1
    # sigma = 0.25

    # std_of_var_process = 0.1
    # mean_reversion_speed = 3
    # long_term_var_mean = 0.06
    # corr = -0.1

    # poly_coeff = [840, -638, 179, -22, 1]

    # heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
    #                     long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
    #                     poly_coeff=poly_coeff, n=n)

    # print(heston.getStatistic(save_data=True, save_dir=save_directory +"Heston.txt"))


def MertonTest():
    # ===================== Test: Merton Call ======================
    print("Merton Call 開始模擬...")
    S0 = 100
    r = 0.05
    T = 0.25
    sigma = 0.2

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coeff = [-100, 1]

    process = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
                      jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)
    print(process.getStatistic(save_data=True, save_dir=save_directory +"MJD.txt"))

    # ===================== Test: Merton Polynomial right up  ======================
    print("Merton right up 開始模擬...")
    S0 = 80
    T = 0.25
    r = 0.1
    sigma = 0.4

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coeff = [-2725, -100, 1]

    process = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
                      jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)
    print(process.getStatistic(save_data=True, save_dir=save_directory + "MJD.txt"))

    # ===================== Test: Merton Polynomial left up  ======================
    print("Merton left up 開始模擬...")
    S0 = 110
    T = 0.25
    r = 0.03
    sigma = 0.4

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coeff = [947.1, -30.164, 0.309, -0.001]

    process = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
                      jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=30)
    print(process.getStatistic(save_data=True, save_dir=save_directory +"MJD.txt"))

    # ===================== Test: Merton Polynomial both up  ======================
    print("Merton both up 開始模擬...")
    S0 = 15
    T = 0.25
    r = 0.1
    sigma = 0.3

    # 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
    jump_intensity = 140
    jump_mean = 0.01
    jump_var = 0.02 ** 2

    poly_coeff = [44.235, -39.474, 5.4793, -0.2358, 0.0031]

    process = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
                      jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)
    print(process.getStatistic(save_data=True, save_dir=save_directory +"MJD.txt"))

    # # ===================== Test: Merton Polynomial both up[老師版]  ======================
    # print("Merton both up 開始模擬...")
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
    # poly_coeff = [840, -638, 179, -22, 1]
    #
    # process = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
    #                   jump_mean=jump_mean, jump_var=jump_var, n=n)
    # print(process.getStatistic(save_data=True, save_dir=save_directory +"MJD.txt"))


def KJDTest():
    # ===================== Test: KDJ Call ======================
    print("KDJ Call 開始模擬...")
    S0 = 100
    r = 0.05
    T = 0.25
    sigma = 0.2

    # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095
    jump_intensity = 1
    p = 0.4
    eta1 = 10
    eta2 = 5

    poly_coeff = [-100, 1]

    process = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity, p=p,
                      eta1=eta1, eta2=eta2, N_line=N, N_repeat=N_repeat, n=n)
    print(process.getStatistic(save_data=True, save_dir=save_directory + "KJD.txt"))

    # ===================== Test: KDJ Polynomial right up  ======================
    print("KDJ right up 開始模擬...")
    S0 = 80
    T = 0.25
    r = 0.1
    sigma = 0.4

    # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095  Ans: 9.14732
    jump_intensity = 1
    p = 0.4
    eta1 = 10
    eta2 = 5

    poly_coeff = [-2725, -100, 1]

    process = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity, p=p,
                      eta1=eta1, eta2=eta2, N_line=N, N_repeat=N_repeat, n=n)
    print(process.getStatistic(save_data=True, save_dir=save_directory + "KJD.txt"))

    # ===================== Test: KJD Polynomial left up  ======================
    print("KDJ left up 開始模擬...")
    S0 = 110
    T = 0.25
    r = 0.03
    sigma = 0.4

    # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095  Ans: 9.14732
    jump_intensity = 1
    p = 0.4
    eta1 = 10
    eta2 = 5

    poly_coeff = [947.1, -30.164, 0.309, -0.001]

    process = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity, p=p,
                      eta1=eta1, eta2=eta2, N_line=N, N_repeat=N_repeat, n=n)
    print(process.getStatistic(save_data=True, save_dir=save_directory + "KJD.txt"))

    # ===================== Test: KJD Polynomial both up  ======================
    print("KDJ both up 開始模擬...")
    S0 = 15
    T = 0.25
    r = 0.1
    sigma = 0.3

    # 資料來自於於 2002 A Jump-Diffusion Model for Option Pricing p1095  Ans: 9.14732
    jump_intensity = 1
    p = 0.4
    eta1 = 10
    eta2 = 5

    poly_coeff = [44.235, -39.474, 5.4793, -0.2358, 0.0031]

    process = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity, p=p,
                      eta1=eta1, eta2=eta2, N_line=N, N_repeat=N_repeat, n=n)
    print(process.getStatistic(save_data=True, save_dir=save_directory + "KJD.txt"))


def SVJTest():
    # ===================== Test: SVJ Call ======================
    print("SVJ Call 開始模擬...")
    S0 = 100
    r = 0.05
    T = 0.25
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
                      poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "SVJ.txt"))

    # ===================== Test: SVJ Polynomial right up  ======================
    print("SVJ right up 開始模擬...")
    S0 = 80
    T = 0.25
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

    poly_coeff = [-2725, -100, 1]

    process = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                      long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                      poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n,
                      N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "SVJ.txt"))

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
                      poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "SVJ.txt"))

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
                      poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "SVJ.txt"))

    # # ===================== Test: SVJ Polynomial both up[老師版]  ======================
    # print("SVJ both up 開始模擬...")
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
    # poly_coeff = [840, -638, 179, -22, 1]
    # process = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
    #                   long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
    #                   poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n)
    #
    # print(process.getStatistic(save_data=True, save_dir=save_directory + "SVJ.txt"))

def VGTest():
    # ===================== Test: VG Call ======================
    print("VG Call 開始模擬...")
    S0 = 100
    r = 0.05
    T = 0.25
    sigma = 0.2

    gamma_mean = -0.14
    gamma_var = 0.2


    poly_coeff = [-90, 1]

    process = VGByMC(S0=S0, r=r, sigma=sigma, T=T, gamma_mean=gamma_mean, gamma_var=gamma_var,
                      poly_coeff=poly_coeff,
                      N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "VG.txt"))

    # ===================== Test: VG Polynomial right up  ======================
    print("VG right up 開始模擬...")
    S0 = 80
    T = 0.25
    r = 0.1
    sigma = 0.4

    gamma_mean = -0.14
    gamma_var = 0.2

    poly_coeff = [-2725, -100, 1]

    process = VGByMC(S0=S0, r=r, sigma=sigma, T=T, gamma_mean=gamma_mean, gamma_var=gamma_var,
                     poly_coeff=poly_coeff,
                     N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "VG.txt"))

    # ===================== Test: VG Polynomial left up  ======================
    print("VG left up 開始模擬...")
    S0 = 110
    T = 0.25
    r = 0.03
    sigma = 0.4

    gamma_mean = -0.14
    gamma_var = 0.2

    poly_coeff = [947.1, -30.164, 0.309, -0.001]

    process = VGByMC(S0=S0, r=r, sigma=sigma, T=T, gamma_mean=gamma_mean, gamma_var=gamma_var,
                     poly_coeff=poly_coeff,
                     N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "VG.txt"))

    # ===================== Test: VG Polynomial both up  ======================
    print("VG both up 開始模擬...")
    S0 = 15
    T = 0.25
    r = 0.1
    sigma = 0.3

    gamma_mean = -0.14
    gamma_var = 0.2

    poly_coeff = [44.235, -39.474, 5.4793, -0.2358, 0.0031]

    process = VGByMC(S0=S0, r=r, sigma=sigma, T=T, gamma_mean=gamma_mean, gamma_var=gamma_var,
                     poly_coeff=poly_coeff,
                     N_line=N, N_repeat=N_repeat)
    print(process.getStatistic(save_data=True, save_dir=save_directory + "VG.txt"))

def NIGTest():
    # ===================== Test: NIG Call ======================
    print("NIG Call 開始模擬...")
    S0 = 100
    r = 0.05
    T = 0.25
    sigma = 0.2

    delta = 1.326
    alpha = 15.624
    beta = 4.025

    poly_coeff = [-90, 1]

    process = NIGByMC(S0=S0, r=r, sigma=sigma, T=T, delta=delta, alpha=alpha, beta=beta,
                     poly_coeff=poly_coeff,
                     N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "NIG.txt"))

    # ===================== Test: NIG Polynomial right up  ======================
    print("NIG right up 開始模擬...")
    S0 = 80
    T = 0.25
    r = 0.1
    sigma = 0.4

    delta = 1.326
    alpha = 15.624
    beta = 4.025

    poly_coeff = [-2725, -100, 1]

    process = NIGByMC(S0=S0, r=r, sigma=sigma, T=T, delta=delta, alpha=alpha, beta=beta,
                      poly_coeff=poly_coeff,
                      N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "NIG.txt"))

    # ===================== Test: NIG Polynomial left up  ======================
    print("NIG left up 開始模擬...")
    S0 = 110
    T = 0.25
    r = 0.03
    sigma = 0.4

    delta = 1.326
    alpha = 15.624
    beta = 4.025

    poly_coeff = [947.1, -30.164, 0.309, -0.001]

    process = NIGByMC(S0=S0, r=r, sigma=sigma, T=T, delta=delta, alpha=alpha, beta=beta,
                      poly_coeff=poly_coeff,
                      N_line=N, N_repeat=N_repeat)

    print(process.getStatistic(save_data=True, save_dir=save_directory + "NIG.txt"))

    # ===================== Test: NIG Polynomial both up  ======================
    print("NIG both up 開始模擬...")
    S0 = 15
    T = 0.25
    r = 0.1
    sigma = 0.3

    delta = 1.326
    alpha = 15.624
    beta = 4.025

    poly_coeff = [44.235, -39.474, 5.4793, -0.2358, 0.0031]

    process = NIGByMC(S0=S0, r=r, sigma=sigma, T=T, delta=delta, alpha=alpha, beta=beta,
                      poly_coeff=poly_coeff,
                      N_line=N, N_repeat=N_repeat)
    print(process.getStatistic(save_data=True, save_dir=save_directory + "NIG.txt"))

if __name__ == "__main__":
    # # threading_GBM = threading.Thread(target=GBMTest, name="GBM")
    # thread_heston = threading.Thread(target=HestonTest, name="Heston")
    # thread_MJD = threading.Thread(target=MertonTest, name="MJD")
    # # threading_GBM.start()
    # # time.sleep(1)
    # thread_heston.start()
    # # time.sleep(1)
    # thread_MJD.start()
    # # time.sleep(1)
    # # threading_GBM.join()
    # thread_heston.join()
    # thread_MJD.join()
    # thread_KJD = threading.Thread(target=KJDTest, name="KJD")
    # thread_SVJ = threading.Thread(target=SVJTest, name="SVJ")

    # thread_KJD.start()
    # # time.sleep(1)
    # thread_SVJ.start()

    # thread_KJD.join()

    # GBMTest()
    # HestonTest()
    # MertonTest()
    # KJDTest()
    # SVJTest()

    VGTest()
    # NIGTest()