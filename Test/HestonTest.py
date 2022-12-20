# sys.path.append("D:\\Users\\DZLin\\Oplynomial_Option_Reserch-master")

from PolynomialPricingMethod.COSMethod import Polynomial
from PolynomialPricingMethod.utils.Tools import timeit

print("\nT = 1")
print("reference val = 5.785155450")
for N in [23, 64, 96, 128, 160, 192]:
    print("----------------------------")
    polynomial1 = Polynomial(S0=100, T=1, r=0.0, sigma=0.0175, process="Heston",
                             N=N, lower_limit=0.01, upper_limit=350,
                             poly_coef=[-100, 1], positive_interval=[100, 200],
                             mean_reversion_speed=1.5768, long_term_var_mean=0.0398,
                             var_variance_process=0.5751, heston_corr=-0.5711)
    COSAns = timeit(polynomial1.PricingByCOSMethod)
    print(f"N={N}", COSAns, f"error={abs(COSAns - 5.785155450):.6f}")

print("\nT = 10")
print("reference val = 22.318945791")
for N in [32, 63, 96, 128, 160, 192]:
    print("----------------------------")
    polynomial2 = Polynomial(S0=100, T=10, r=0.0, sigma=0.0175, process="Heston",
                             N=N, lower_limit=0.01, upper_limit=1000,
                             poly_coef=[-100, 1], positive_interval=[100, 1000],
                             mean_reversion_speed=1.5768, long_term_var_mean=0.0398,
                             var_variance_process=0.5751, heston_corr=-0.5711, )
    COSAns = timeit(polynomial2.PricingByCOSMethod)
    print(f"N={N}", COSAns, f"error={abs(COSAns - 22.318945791):.6f}")
