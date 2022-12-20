# sys.path.append("D:\\Users\\DZLin\\Oplynomial_Option_Reserch-master")

from PolynomialPricingMethod.COSMethod import Polynomial
from PolynomialPricingMethod.utils.Tools import timeit

# K = 80
print("\nfor strike price = 80")
print("reference val = 20.7799226309")
for N in [16, 32, 63, 128]:
    print("----------------------------")
    polynomial1 = Polynomial(S0=100, r=0.1, T=0.1, sigma=0.25, process="GBM",
                             N=N, lower_limit=50, upper_limit=150,
                             poly_coef=[-80, 1], positive_interval=[80, 150])
    COSAns = timeit(polynomial1.PricingByCOSMethod)
    print(f"N={N}", COSAns, f"error={abs(COSAns - 20.7799226309):.6f}")

# K = 100
print("\nfor strike price = 100")
print("reference val = 3.659968453")
for N in [16, 32, 63, 128]:
    print("----------------------------")
    polynomial2 = Polynomial(S0=100, r=0.1, T=0.1, sigma=0.25, process="GBM",
                             N=N, lower_limit=50, upper_limit=150,
                             poly_coef=[-100, 1], positive_interval=[100, 150])
    COSAns = timeit(polynomial2.PricingByCOSMethod)
    print(f"N={N}", COSAns, f"error={abs(COSAns - 3.659968453):.6f}")

# K = 120
print("\nfor strike price = 120")
print("reference val = 0.044577814")
for N in [16, 32, 63, 128]:
    print("----------------------------")
    polynomial3 = Polynomial(S0=100, r=0.1, T=0.1, sigma=0.25, process="GBM",
                             N=N, lower_limit=50, upper_limit=150,
                             poly_coef=[-120, 1], positive_interval=[120, 150])
    COSAns = timeit(polynomial3.PricingByCOSMethod)
    print(f"N={N}", COSAns, f"error={abs(COSAns - 0.044577814):.6f}")
