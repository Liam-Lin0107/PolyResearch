# import sys
# sys.path.append("D:\\Users\\DZLin\\Oplynomial_Option_Reserch-master")

from PolynomialPricingMethod.TreeMethod import PolynomialByBinomialTree


# 左翹
S0 = 100
T = 0.5
r = 0.03
sigma = 0.25
poly_coef = [947100, -30164, 309, -1]

polynomialByBinomialTree = PolynomialByBinomialTree(S0=S0, T=T, r=r, sigma=sigma, poly_coef=poly_coef, N=20000)
ans = polynomialByBinomialTree.getValue()
print(f"N={200000}: value: {ans:.6f}")
print("")



# 雙翹
S0 = 5
T = 1
r = 0.1
sigma = 0.25
poly_coef = [840, -638, 179, -22, 1]

polynomialByBinomialTree = PolynomialByBinomialTree(S0=S0, T=T, r=r, sigma=sigma, poly_coef=poly_coef, N=20000)
ans = polynomialByBinomialTree.getValue()
print(f"N={200000}: value: {ans:.6f}")
print("")

#

# 單邊右翹
S0 = 50
T = 0.5
r = 0.1
sigma = 0.4
poly_coef = [-2725, -100, 1]

polynomialByBinomialTree = PolynomialByBinomialTree(S0=S0, T=T, r=r, sigma=sigma, poly_coef=poly_coef, N=20000)
ans = polynomialByBinomialTree.getValue()
print(f"N={20000}: value: {ans:.6f}")
print("")