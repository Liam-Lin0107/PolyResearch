import numpy as np
from PolynomialPricingMethod.COSMethod import PolyByCosMethod
from PolynomialPricingMethod.utils.CharactoristicFunc import *
from numpy import pi, cos, log
import matplotlib.pyplot as plt
import seaborn as sns
from math import inf

# warnings.filterwarnings("ignore")
i = 1j


class DensityRecover:
    """
    尋找在這個payoff下使得誤差能達到0.1%水準的切割點
    N為假設多少個頻率疊合能使density recover最好，默認值為1e6，注意過程中N都不會改變
    默認上界：2000
    默認下界：0.0001
    """

    def __init__(self, S0, process_cf, poly_coeff, positive_interval, N=1e5, step=0.05, assume_true_a=1e-4,
                 assume_true_b=2000, debug=False, verbal=True): # 可以設大
        self.x = log(S0)
        self.process = process_cf
        self.poly_coeff = poly_coeff
        self.positive_interval = positive_interval
        self.N = int(N)
        self.step = step
        self.assume_true_a = log(assume_true_a)
        self.assume_true_b = log(assume_true_b)
        self.debug = debug
        self.verbal = verbal  # whetehr need to print the final range and interval

        if len(positive_interval) != 0:
            self.ini_guess_a = log(self.positive_interval[0]) if self.positive_interval[0] != 0 else log(
                self.positive_interval[1])
            self.ini_guess_b = log(self.positive_interval[-1]) if log(self.positive_interval[-1]) != inf else log(
                self.positive_interval[-2])

            # 確認a是否正確
            if self.positive_interval[0] == 0:
                if self.ini_guess_a == self.ini_guess_b:  # 若起始a與起始b相同，會產生divided zero error
                    self.ini_guess_a = log(self.positive_interval[1] - 1)  # 因為為0的關係所以將a向左移動10元股價

            # 確認b是否正確
            if self.positive_interval[-1] == inf:
                if self.ini_guess_a == self.ini_guess_b:  # 若起始a與起始b相同，會產生divided zero error
                    self.ini_guess_b = log(self.positive_interval[-2] + 1)  # 因為為inf的關係所以將b向右移動10元股價

        else:
            raise "Need to input one of parameter which is (poly_coef) or (initial_guess_range)"

        self.Fk_Norm_ls = []  # use to store the Cos expension coefficients
        self.Fk_ls = []  # use to store the Cos expension coefficients
        self.right_density_plot = []
        self.left_density_plot = []
        self.error_plot = []
        print("fitting the true density...")
        self._calFkCoeff()

    def _calFkCoeff(self):
        k = np.arange(self.N)

        self.Fk_ls = np.array(list(map(self._Fk, k, np.full(self.N, False))))
        self.Fk_Norm_ls = np.array(list(map(self._Fk, k, np.full(self.N, True))))  # 用 k>=0 做一個全部都是true的array

    def _Fk(self, k, isNorm):
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

    def _density(self, x, isNorm=False):
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

    def _Payoff(self, price):
        payoff = 0
        for power, coef in enumerate(self.poly_coeff):
            payoff += coef * price ** power
        return max(payoff, 0)

    def _FindBestDensityRange(self):
        """
        用來尋找density < 0.01的兩側位置
        """
        print("find density cutting point...")
        best_intergral_low = self.ini_guess_a  # lnST
        best_intergral_high = self.ini_guess_b  # lnST
        density_low = self._density(best_intergral_low)  # density of lnST
        density_high = self._density(best_intergral_high)  # density of lnST
        # 加入向下趨勢判斷 verson 2.0
        # 用兩個step的是因為太靠近可能COS疊合沒有完全會有震盪
        density_low_minus = self._density(best_intergral_low - 2 * self.step)
        density_high_plus = self._density(best_intergral_high + 2 * self.step)

        if self.debug:
            count = 0
        while density_low > 0.01 or density_high > 0.01 or density_low_minus > density_low or density_high_plus > density_high:
            if density_low > 0.01 or density_low_minus > density_low:
                best_intergral_low -= self.step
                if best_intergral_low < self.assume_true_a:
                    raise "need to lower the assume true density lower bound"
                # update
                density_low = self._density(best_intergral_low)
                density_low_minus = self._density(best_intergral_low - self.step)
            if density_high > 0.01 or density_high_plus > density_high:
                best_intergral_high += self.step
                if best_intergral_high > self.assume_true_b:
                    raise "need to higher the assume true density upper bound"
                # update
                density_high = self._density(best_intergral_high)
                density_high_plus = self._density(best_intergral_high + self.step)

            if self.debug:
                count += 1
                print("=" * 50)
                print("i = ", count)
                print("=" * 50)
                print("left density:", density_low)
                print("right density:", density_high)
            self.left_density_plot.append(density_low)
            self.right_density_plot.append(density_high)
        return best_intergral_low, best_intergral_high

    def getIntergralRangeAndInterval(self):
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

            best_intergral_low, best_intergral_high = self._FindBestDensityRange()

            print("searching for both side...")
            best_interval_low = min(self.ini_guess_a, best_intergral_low)
            best_interval_high = max(self.ini_guess_b, best_intergral_high)
            if self.debug:
                count = 0
            while error > 1e-4:
                best_interval_low -= self.step
                if best_interval_low < self.assume_true_a:
                    raise "need to lower the assume true density lower bound"
                best_interval_high += self.step
                if best_interval_high > self.assume_true_b:
                    raise "need to higher the assume true density upper bound"

                left_payoff = self._Payoff(exp(best_interval_low))
                left_density = self._density(best_interval_low)

                right_payoff = self._Payoff(exp(best_interval_high))
                right_density = self._density(best_interval_high)

                error = left_payoff * left_density + right_payoff * right_density
                if self.debug:
                    count += 1
                    print("=" * 50)
                    print("i = ", count)
                    print("=" * 50)
                    print("left_density:", left_density)
                    print("left_payoff:", left_payoff)
                    print("left point", best_interval_low)
                    print("right_density:", right_density)
                    print("right_payoff:", right_payoff)

                self.left_density_plot.append(left_density)
                self.right_density_plot.append(right_density)
                self.error_plot.append(error)

            if self.verbal:
                print("best integral lower bound:", exp(best_interval_low))
                print("best integral upper bound:", exp(best_interval_high))
                print("best interval lower bound:", exp(best_interval_low))
                print("best interval upper bound:", exp(best_interval_high))
            return exp(best_interval_low), exp(best_interval_high), exp(best_interval_low), exp(best_interval_high)

        # 右側尋找
        elif self.positive_interval[0] != 0 and self.positive_interval[-1] == inf:

            best_intergral_low, best_intergral_high = self._FindBestDensityRange()

            print("right searching....")
            best_interval_high = max(self.ini_guess_b, best_intergral_high)
            if self.debug:
                count = 0
            while error > 1e-4:
                best_interval_high += self.step
                if best_interval_high > self.assume_true_b:
                    raise "need to higher the assume true density upper bound"
                right_payoff = self._Payoff(exp(best_interval_high))
                right_density = self._density(best_interval_high)
                error = right_density * right_payoff

                if self.debug:
                    count += 1
                    print("=" * 50)
                    print("i = ", count)
                    print("=" * 50)
                    print("right_density:", right_density)
                    print("right_payoff:", right_payoff)

                self.right_density_plot.append(right_density)
                self.error_plot.append(error)
            if self.verbal:
                print("best integral lower bound:", exp(best_intergral_low))
                print("best integral upper bound:", exp(best_interval_high))
                print("best interval lower bound:", exp(self.ini_guess_a))
                print("best interval upper bound:", exp(best_interval_high))
            return exp(best_intergral_low), exp(best_interval_high), exp(self.ini_guess_a), exp(best_interval_high)

        # 左側尋找
        elif self.positive_interval[0] == 0 and self.positive_interval[-1] != inf:

            best_intergral_low, best_intergral_high = self._FindBestDensityRange()

            print("left searching....")
            best_interval_low = min(self.ini_guess_a, best_intergral_low)
            if self.debug:
                count = 0
            while error > 1e-4:
                best_interval_low -= self.step
                if best_interval_low < self.assume_true_a:
                    raise "need to lower the assume true density lower bound"
                left_payoff = self._Payoff(exp(best_interval_low))
                left_density = self._density(best_interval_low)

                error = left_density * left_payoff

                if self.debug:
                    count += 1
                    print("=" * 50)
                    print("i = ", count)
                    print("=" * 50)
                    print("left_density:", left_density)
                    print("left_payoff:", left_payoff)

                self.left_density_plot.append(left_density)
                self.error_plot.append(error)

            if self.verbal:
                print("best integral lower bound:", exp(best_interval_low))
                print("best integral upper bound:", exp(best_intergral_high))
                print("best interval lower bound:", exp(best_interval_low))
                print("best interval upper bound:", exp(self.ini_guess_b))
            return exp(best_interval_low), exp(best_intergral_high), exp(best_interval_low), exp(self.ini_guess_b)

        # 不用尋找
        else:

            best_intergral_low, best_intergral_high = self._FindBestDensityRange()
            best_intergral_low = best_intergral_low if best_intergral_low < self.ini_guess_a else self.ini_guess_a
            best_intergral_high = best_intergral_high if best_intergral_high > self.ini_guess_b else self.ini_guess_b
            print("no searching....")

            if self.verbal:
                print("best integral lower bound:", exp(best_intergral_low))
                print("best integral upper bound:", exp(best_intergral_high))
                print("best interval lower bound:", exp(self.ini_guess_a))
                print("best interval upper bound:", exp(self.ini_guess_b))
            return exp(best_intergral_low), exp(best_intergral_high), exp(self.ini_guess_a), exp(self.ini_guess_b)

    ############# 以下是為了畫圖用的 ######################
    def poltDensity(self, plot_save=False, file_name_prefix=""):
        """
        畫的是lnS0 = 0下的密度函數
        """
        density_ls = []
        x_ls = np.linspace(-1.5, 1.5, 100)
        for x in x_ls:
            density_ls.append(self._density(x, isNorm=True))
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
