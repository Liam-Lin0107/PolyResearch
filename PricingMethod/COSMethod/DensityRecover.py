import numpy as np
from PricingMethod.COSMethod.PolynomialCosMethod import PolynomialCosMethod
from PricingMethod.COSMethod.CharacteristicFunction import *
from numpy import pi, cos, log
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from math import inf
# warnings.filterwarnings("ignore")
i = 1j

class DensityRecover:
    """
    尋找在這個payoff下使得誤差能達到0.1%水準的切割點
    N為假設多少個頻率疊合能使density recover最好，默認值為2000，注意過程中N都不會改變
    新版設定上界為3000, 下界為0.01
    """
    def __init__(self, S0, process, poly_coef, positive_interval, N=1e6, step=0.05, assume_true_a=1e-10, assume_true_b=2000):
        self.x = log(S0)
        self.process = process
        self.poly_coef = poly_coef
        self.positive_interval = positive_interval
        self.N = int(N)
        self.step = step
        self.assume_true_a = log(assume_true_a)
        self.assume_true_b = log(assume_true_b)


        if len(positive_interval) != 0:
            self.ini_guess_a = log(self.positive_interval[0]) if self.positive_interval[0] != 0 else log(self.positive_interval[1])
            self.ini_guess_b = log(self.positive_interval[-1]) if log(self.positive_interval[-1]) != inf else log(self.positive_interval[-2])

            # 確認a是否正確
            if self.positive_interval[0] == 0:
                if self.ini_guess_a == self.ini_guess_b:  # 若起始a與起始b相同，會產生divided zero error
                    self.ini_guess_a = log(self.positive_interval[1] - 1) # 因為為0的關係所以將a向左移動10元股價

            # 確認b是否正確
            if self.positive_interval[-1] == inf:
                if self.ini_guess_a == self.ini_guess_b: # 若起始a與起始b相同，會產生divided zero error
                    self.ini_guess_b = log(self.positive_interval[-2] + 1) # 因為為inf的關係所以將b向右移動10元股價

        else:
            raise "Need to input one of parameter which is (poly_coef) or (initial_guess_range)"

        self.Fk_Norm_ls = [] # use to store the Cos expension coefficients
        self.Fk_ls = [] # use to store the Cos expension coefficients
        self.right_density_plot = []
        self.left_density_plot = []
        self.error_plot = []
        print("fitting the true density...")
        self.calFkCoeff()

    def calFkCoeff(self):
        for k in range(self.N):
            self.Fk_ls.append(self.Fk(k, isNorm=False))
            self.Fk_Norm_ls.append(self.Fk(k, isNorm=True))

    def Fk(self, k, isNorm):
        """
         是使用lnS0 = S0下的cf
        """
        if isNorm:
            inner = self.process.getCFValue(k * pi / (self.assume_true_b - self.assume_true_a)) * exp(
                -i * k * pi * self.assume_true_a / (self.assume_true_b - self.assume_true_a))
            return 2 / (self.assume_true_b - self.assume_true_a) * np.real(inner)
        else:
            inner = self.process.getCFValue(k * pi / (self.assume_true_b - self.assume_true_a)) * \
                    exp(i * k * pi * (self.x - self.assume_true_a) / (self.assume_true_b - self.assume_true_a))
            return 2 / (self.assume_true_b - self.assume_true_a) * np.real(inner)

    def density(self, x, isNorm = False):
        """
        isNorm表示lnST是否為零，畫圖時要使用isNorm
        """
        if isNorm:
            density = 0
            for k in range(self.N):
                if k == 0:
                    density += 0.5 * self.Fk_Norm_ls[k] * \
                               cos(k * pi * (x - self.assume_true_a) / (self.assume_true_b - self.assume_true_a))
                else:
                    density += self.Fk_Norm_ls[k] * \
                               cos(k * pi * (x - self.assume_true_a) / (self.assume_true_b - self.assume_true_a))
            return density
        else:
            density = 0
            for k in range(self.N):
                if k == 0:
                    density += 0.5 * self.Fk_ls[k] * \
                               cos(k * pi * (x - self.assume_true_a) / (self.assume_true_b - self.assume_true_a))
                else:
                    density += self.Fk_ls[k] * \
                               cos(k * pi * (x - self.assume_true_a) / (self.assume_true_b - self.assume_true_a))
            return density

    def Payoff(self, price):
        payoff = 0
        for power, coef in enumerate(self.poly_coef):
            payoff += coef * price ** power
        return max(payoff, 0)


    def getBestIntervalPoints(self):
        """
        回傳值分別為：
        1. 最佳積分下界
        2. 最佳積分上界
        3. 最佳區間下界
        4. 最佳區間上界
        我在density計算左右兩邊各加減0.01目的是不要讓所求的值在函數邊界
        """

        # 兩側尋找
        error = inf
        if self.positive_interval[0] == 0 and self.positive_interval[-1] == inf:
            print("finding density best fitting points...")
            best_intergral_low = self.ini_guess_a
            best_intergral_high = self.ini_guess_b
            density_low = self.density(best_intergral_low)
            density_high = self.density(best_intergral_high)

            count = 0
            while density_low > 0.01 or density_high > 0.01:
                if density_low > 0.01:
                    best_intergral_low -= self.step
                    if best_intergral_low < self.assume_true_a:
                        raise "need to lower the assume true density lower bound"
                if density_high > 0.01:
                    best_intergral_high += self.step
                    if best_intergral_high > self.assume_true_b:
                        raise "need to higher the assume true density upper bound"
                if density_low > 0.01:
                    density_low = self.density(best_intergral_low)
                if density_high > 0.01:
                    density_high = self.density(best_intergral_high)
                    density_low = density_low if density_low > 0 else 0
                # count += 1
                # print("=" * 50)
                # print("i = ", count)
                # print("=" * 50)
                # print("left density:", density_low)
                # print("right density:", density_high)
                # self.left_density_plot.append(density_low)
                # self.right_density_plot.append(density_high)


            print("searching for both side...")
            best_interval_low = min(self.ini_guess_a, best_intergral_low)
            best_interval_high = max(self.ini_guess_b, best_intergral_high)

            # count = 0
            while error > 1e-4:
                best_interval_low -= self.step
                if best_interval_low < self.assume_true_a:
                    raise "need to lower the assume true density lower bound"
                best_interval_high += self.step
                if best_interval_high > self.assume_true_b:
                    raise "need to higher the assume true density upper bound"


                left_payoff = self.Payoff(exp(best_interval_low))
                left_density = self.density(best_interval_low)

                right_payoff = self.Payoff(exp(best_interval_high))
                right_density = self.density(best_interval_high)


                error = left_payoff * left_density + right_payoff * right_density
                # count += 1
                # print("=" * 50)
                # print("i = ",count)
                # print("=" * 50)
                # print("left_density:", left_density)
                # print("left_payoff:", left_payoff)
                # print("left point", best_interval_low)
                # print("right_density:", right_density)
                # print("right_payoff:", right_payoff)

                self.left_density_plot.append(left_density)
                self.right_density_plot.append(right_density)
                self.error_plot.append(error)


            return exp(best_interval_low), exp(best_interval_high), exp(best_interval_low), exp(best_interval_high)

        # 右側尋找
        elif self.positive_interval[0] != 0 and self.positive_interval[-1] == inf:
            print("finding density best fitting points...")
            best_intergral_low = self.ini_guess_a
            best_intergral_high = self.ini_guess_b
            density_low = self.density(best_intergral_low)
            density_high = self.density(best_intergral_high)

            # count = 0
            while density_low > 0.01 or density_high > 0.01:
                if density_low > 0.01:
                    best_intergral_low -= self.step
                    if best_intergral_low < self.assume_true_a:
                        raise "need to lower the assume true density lower bound"
                if density_high > 0.01:
                    best_intergral_high += self.step
                    if best_intergral_high > self.assume_true_b:
                        raise "need to higher the assume true density upper bound"
                if density_low > 0.01:
                    density_low = self.density(best_intergral_low)
                if density_high > 0.01:
                    density_high = self.density(best_intergral_high)
                # count += 1
                # print("=" * 50)
                # print("i = ", count)
                # print("=" * 50)
                # print("left density:", density_low)
                # print("right density:", density_high)

                self.left_density_plot.append(density_low)
                self.right_density_plot.append(density_high)

            print("right searching....")
            # count = 0
            best_interval_low = min(self.ini_guess_a, best_intergral_low)
            best_interval_high = max(self.ini_guess_b, best_intergral_high)
            while error > 1e-4:
                best_interval_high += self.step
                right_payoff = self.Payoff(exp(best_interval_high))
                right_density = self.density(best_interval_high)
                error = right_density * right_payoff
                # count += 1
                # print("=" * 50)
                # print("i = ", count)
                # print("=" * 50)
                # print("right_density:", right_density)
                # print("right_payoff:", right_payoff)


                self.right_density_plot.append(right_density)
                self.error_plot.append(error)

            return exp(best_intergral_low), exp(best_interval_high), exp(self.ini_guess_a), exp(best_interval_high)

        # 左側尋找
        elif self.positive_interval[0] == 0 and self.positive_interval[-1] != inf:
            print("finding density best fitting points...")
            best_intergral_low = self.ini_guess_a
            best_intergral_high = self.ini_guess_b
            density_low = self.density(best_intergral_low)
            density_high = self.density(best_intergral_high)

            # count = 0
            while density_low > 0.01 or density_high > 0.01:
                if density_low > 0.01:
                    best_intergral_low -= self.step
                    if best_intergral_low < self.assume_true_a:
                        raise "need to lower the assume true density lower bound"
                if density_high > 0.01:
                    best_intergral_high += self.step
                    if best_intergral_high > self.assume_true_b:
                        raise "need to higher the assume true density upper bound"
                if density_low > 0.01:
                    density_low = self.density(best_intergral_low)
                if density_high > 0.01:
                    density_high = self.density(best_intergral_high)
                # count += 1
                # print("=" * 50)
                # print("i = ", count)
                # print("=" * 50)
                # print("left density:", density_low)
                # print("right density:", density_high)

                self.left_density_plot.append(density_low)
                self.right_density_plot.append(density_high)

            print("left searching....")
            best_interval_low = min(self.ini_guess_a, best_intergral_low)
            best_interval_high = max(self.ini_guess_b, best_intergral_high)
            # count = 0
            while error > 1e-4:
                best_interval_low -= self.step

                left_payoff = self.Payoff(exp(best_interval_low))
                left_density = self.density(best_interval_low)

                error = left_density * left_payoff
                # count += 1
                # print("=" * 50)
                # print("i = ", count)
                # print("=" * 50)
                # print("left_density:", left_density)
                # print("left_payoff:", left_payoff)

                self.left_density_plot.append(left_density)
                self.error_plot.append(error)
            return exp(best_interval_low), exp(best_intergral_high), exp(best_interval_low), exp(self.ini_guess_b)

        # 不用尋找
        else:
            print("finding density best fitting points...")
            best_intergral_low = self.ini_guess_a
            best_intergral_high = self.ini_guess_b
            density_low = self.density(best_intergral_low)
            density_high = self.density(best_intergral_high)

            # count = 0
            while density_low > 0.01 or density_high > 0.01:
                if density_low > 0.01:
                    best_intergral_low -= self.step
                    if best_intergral_low < self.assume_true_a:
                        raise "need to lower the assume true density lower bound"
                if density_high > 0.01:
                    best_intergral_high += self.step
                    if best_intergral_high > self.assume_true_a:
                        raise "need to higher the assume true density upper bound"
                if density_low > 0.01:
                    density_low = self.density(best_intergral_low)
                if density_high > 0.01:
                    density_high = self.density(best_intergral_high)
                # count += 1
                # print("=" * 50)
                # print("i = ", count)
                # print("=" * 50)
                # print("left density:", density_low)
                # print("right density:", density_high)

                self.left_density_plot.append(density_low)
                self.right_density_plot.append(density_high)

            print("no searching....")
            return exp(self.ini_guess_a), exp(self.ini_guess_b), exp(self.ini_guess_a), exp(self.ini_guess_b)

    ############# 以下是為了畫圖用的 ######################

    def poltDensity(self, plot_save=False, file_name_prefix=""):
        """
        畫的是lnS0 = 0下的密度函數
        """
        density_ls = []
        x_ls = np.linspace(-1.5, 1.5, 100)
        for x in x_ls:
            density_ls.append(self.density(x, isNorm=True))
        plt.title("Density")
        sns.scatterplot(x=x_ls, y=density_ls)
        if plot_save:
            plt.savefig(f"{self.process.__class__.__name__}_{file_name_prefix}_density.jpg")
        plt.show();

    def plotFittingDetal(self, plot_save=False, file_name_prefix=""):
        self.poltDensity(plot_save, file_name_prefix)
        plt.title("Right Density Curve")
        sns.scatterplot(x=np.arange(len(self.right_density_plot)), y=self.right_density_plot)
        if plot_save:
            plt.savefig(f"{self.process.__class__.__name__}_{file_name_prefix}_{1}.jpg")
        plt.show();
        plt.title("Left Density curve")
        sns.scatterplot(x=np.arange(len(self.left_density_plot)), y=self.left_density_plot)
        if plot_save:
            plt.savefig(f"{self.process.__class__.__name__}_{file_name_prefix}_{2}.jpg")
        plt.show();
        plt.title("Truncate Error Curve")
        sns.scatterplot(x=np.arange(len(self.error_plot)), y=self.error_plot)
        if plot_save:
            plt.savefig(f"{self.process.__class__.__name__}_{file_name_prefix}_{3}.jpg")
        plt.show();



if __name__ == "__main__":
    positive_interval = [90, inf]
    process_GBM = GBM(T=1, r=0.1, sigma=0.25)
    poly_coef = [-90, 1]
    densityRecover = DensityRecover(S0=100, process=process_GBM,
                                    poly_coef=poly_coef, N=1e5,
                                    positive_interval=positive_interval, step=0.05)
    # densityRecover.poltDensity()
    a, b, positive_interval[0], positive_interval[-1] = densityRecover.getBestIntervalPoints()
    print("a=",a, "b=",b, positive_interval[0], positive_interval[-1])

    polynomialCosMethod = PolynomialCosMethod(S0=100, r=0.1, sigma=0.25, T=1, process=process_GBM,
                                              N=1000, lower_limit=a, upper_limit=b,
                                              poly_coef=poly_coef, positive_interval=positive_interval)

    print(polynomialCosMethod.getValue())



    #
    # S0 = 5
    # T = 1
    # r = 0.1
    # sigma = 0.25
    # poly_coef = [840, -638, 179, -22, 1]
    # positive_interval = [0, 4, 5, 6, 7, inf]
    #
    # proccess = GBM(T=T, r=r, sigma=sigma)
    #
    # densityRecover = DensityRecover(S0=S0, process=proccess,
    #                                 poly_coef=poly_coef,
    #                                 positive_interval=positive_interval,
    #                                 step=0.05)
    #
    # a, b, positive_interval[0], positive_interval[-1] = densityRecover.getBestIntervalPoints()
    # print("a=", a, "b=", b, positive_interval[0], positive_interval[-1])
    #
    # polynomialCosMethod = PolynomialCosMethod(S0=S0, r=r, sigma=sigma, T=T, process=proccess,
    #                                           N=10000, lower_limit=a, upper_limit=b,
    #                                           poly_coef=poly_coef, positive_interval=positive_interval)
    #
    # densityRecover.plotFittingDetal()
    # print(polynomialCosMethod.getValue())
    #
