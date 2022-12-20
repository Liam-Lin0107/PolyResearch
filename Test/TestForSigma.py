from math import exp

from PolynomialPricingMethod.COSMethod import Polynomial
import numpy as np

for sigma in np.linspace(0.01, 1, 50):
    polynomial1 = Polynomial(S0=10, T=1, r=0.1, sigma=sigma, process="GBM.txt",
                             poly_coef=[-10, 1],
                             positive_interval=[10, exp(9)],
                             N=2000,
                             lower_limit=exp(-9), upper_limit=exp(9))
    COS_ans = polynomial1.PricingByCOSMethod()
    Tree_ans = polynomial1.PricingByTree()

    print("\nTest for sigma =", sigma)
    print("COS: ", COS_ans)
    print("Tree: ", Tree_ans)
