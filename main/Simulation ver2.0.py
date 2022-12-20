import sys
sys.path.append(
    "D:\geffc\Documents\Research\PolyOptionPricingRearch\Program")  # for 85
import threading
from PolynomialPricingMethod.PathSimulationMethod import *
import time

n = 252  # for setting dt
N = int(1e3)  # N_line
N_repeat = 10
save_directory = "./Data/Simulation/"

def callOption():
      print("Call simulation start...")

      S0 = 100
      r = 0.05
      T = 0.25
      sigma = 0.2
      # SV parameters
      std_of_var_process = 0.1
      mean_reversion_speed = 3
      long_term_var_mean = 0.04
      corr = -0.1
      # Normal Jump parameters
      jump_intensity = 140
      jump_mean = 0.01
      jump_var = 0.02 ** 2
      # Double Exponential jump parameters
      jump_intensity_KDJ = 1
      p = 0.4
      eta1 = 10
      eta2 = 5

      poly_coeff = [-100, 1]
    #   gbm = GBMByMC(S0=S0, r=r, sigma=sigma, T=T,
    #               poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat)

      heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat)

      mjd = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
                  jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)

      kjd = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity_KDJ, p=p,
                  eta1=eta1, eta2=eta2, N_line=N, N_repeat=N_repeat)

      svj = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                  long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                  poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)
    #   print(gbm.getStatistic(save_data=True,
    #         save_dir=(save_directory + "GBM.txt")))

      print(heston.getStatistic(save_data=True,
            save_dir=(save_directory + "Heston.txt")))

      print(mjd.getStatistic(save_data=True,
                              save_dir=(save_directory + "MJD.txt")))

      print(kjd.getStatistic(save_data=True,
                              save_dir=(save_directory + "KJD.txt")))

      print(svj.getStatistic(save_data=True,
                              save_dir=(save_directory + "SVJ.txt")))


def rightUpOption():
      print("right up simulation start...")
      S0 = 80
      T = 0.25
      r = 0.1
      sigma = 0.4
      # SV parameters
      std_of_var_process = 0.1
      mean_reversion_speed = 3
      long_term_var_mean = 0.06
      corr = -0.1
      # Normal Jump parameters
      jump_intensity = 140
      jump_mean = 0.01
      jump_var = 0.02 ** 2
      # Double Exponential jump parameters
      jump_intensity_KDJ = 1
      p = 0.4
      eta1 = 10
      eta2 = 5

      poly_coeff = [-2725, -100, 1]
    #   gbm = GBMByMC(S0=S0, r=r, sigma=sigma, T=T,
    #               poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat)

      heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat)

      mjd = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
                  jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)

      kjd = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity_KDJ, p=p,
                  eta1=eta1, eta2=eta2, N_line=N, N_repeat=N_repeat)

      svj = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                  long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                  poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)
    #   print(gbm.getStatistic(save_data=True,
    #         save_dir=(save_directory + "GBM.txt")))

      print(heston.getStatistic(save_data=True,
            save_dir=(save_directory + "Heston.txt")))

      print(mjd.getStatistic(save_data=True,
                              save_dir=(save_directory + "MJD.txt")))

      print(kjd.getStatistic(save_data=True,
                              save_dir=(save_directory + "KJD.txt")))

      print(svj.getStatistic(save_data=True,
                              save_dir=(save_directory + "SVJ.txt")))


def leftUption():
      print("left up simulation start...")

      S0 = 110
      T = 0.25
      r = 0.03
      sigma = 0.4
      # SV parameters
      std_of_var_process = 0.1
      mean_reversion_speed = 3
      long_term_var_mean = 0.06
      corr = -0.1
      # Normal Jump parameters
      jump_intensity = 140
      jump_mean = 0.01
      jump_var = 0.02 ** 2
      # Double Exponential jump parameters
      jump_intensity_KDJ = 1
      p = 0.4
      eta1 = 10
      eta2 = 5

      poly_coeff = [947.1, -30.164, 0.309, -0.001]
    #   gbm = GBMByMC(S0=S0, r=r, sigma=sigma, T=T,
    #               poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat)

      heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat)

      mjd = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
                  jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)

      kjd = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity_KDJ, p=p,
                  eta1=eta1, eta2=eta2, N_line=N, N_repeat=N_repeat)

      svj = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                  long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                  poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)
    #   print(gbm.getStatistic(save_data=True,
    #         save_dir=(save_directory + "GBM.txt")))

      print(heston.getStatistic(save_data=True,
            save_dir=(save_directory + "Heston.txt")))

      print(mjd.getStatistic(save_data=True,
                              save_dir=(save_directory + "MJD.txt")))

      print(kjd.getStatistic(save_data=True,
                              save_dir=(save_directory + "KJD.txt")))

      print(svj.getStatistic(save_data=True,
                              save_dir=(save_directory + "SVJ.txt")))


def bothUpOption():
      print("both up simulation start...")

      S0 = 15
      T = 0.25
      r = 0.1
      sigma = 0.3
      # SV parameters
      std_of_var_process = 0.1
      mean_reversion_speed = 3
      long_term_var_mean = 0.04
      corr = -0.1
      # Normal Jump parameters
      jump_intensity = 140
      jump_mean = 0.01
      jump_var = 0.02 ** 2
      # Double Exponential jump parameters
      jump_intensity_KDJ = 1
      p = 0.4
      eta1 = 10
      eta2 = 5

      poly_coeff = [44.235, -39.474, 5.4793, -0.2358, 0.0031]
    #   gbm = GBMByMC(S0=S0, r=r, sigma=sigma, T=T,
    #               poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat)

      heston = HestonByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                        long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                        poly_coeff=poly_coeff, N_line=N, N_repeat=N_repeat)

      mjd = MJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity,
                  jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)

      kjd = KJDByMC(S0=S0, r=r, sigma=sigma, T=T, poly_coeff=poly_coeff, jump_intensity=jump_intensity_KDJ, p=p,
                  eta1=eta1, eta2=eta2, N_line=N, N_repeat=N_repeat)

      svj = SVJByMC(S0=S0, r=r, sigma=sigma, T=T, mean_reversion_speed=mean_reversion_speed,
                  long_term_var_mean=long_term_var_mean, corr=corr, std_of_var_process=std_of_var_process,
                  poly_coeff=poly_coeff, jump_intensity=jump_intensity, jump_mean=jump_mean, jump_var=jump_var, n=n, N_line=N, N_repeat=N_repeat)
    #   print(gbm.getStatistic(save_data=True,
    #         save_dir=(save_directory + "GBM.txt")))

      print(heston.getStatistic(save_data=True,
            save_dir=(save_directory + "Heston.txt")))

      print(mjd.getStatistic(save_data=True,
                              save_dir=(save_directory + "MJD.txt")))

      print(kjd.getStatistic(save_data=True,
                              save_dir=(save_directory + "KJD.txt")))

      print(svj.getStatistic(save_data=True,
                              save_dir=(save_directory + "SVJ.txt")))


def bothDownOption():
    pass

if __name__ == "__main__":
    callOption()
    rightUpOption()
    leftUption()
    bothUpOption()


