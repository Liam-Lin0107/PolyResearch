# sys.path.append("D:\\Users\\DZLin\\Oplynomial_Option_Reserch-master")

from math import exp
from PolynomialPricingMethod.COSMethod import Polynomial
from PolynomialPricingMethod.CloseFormMethod import BSMCall
from PolynomialPricingMethod.utils.Tools import timeit

print("=====Call========")
# ==================Test for S==================
print("Test for S0=5")
polynomial1 = Polynomial(S0=5, T=1, r=0.1, sigma=0.25, process="GBM",
                         poly_coef=[-10, 1],
                         positive_interval=[10, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
bsmCloseForm = BSMCall(S0=5, T=1, r=0.1, sigma=0.25, K=10)
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns, f"error={abs(COSAns - bsmCloseForm.BSM_call_value()):.6f}")
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns, f"error={abs(TreeAns - bsmCloseForm.BSM_call_value()):.6f}")

print("\nTest for S0=100")
polynomial1 = Polynomial(S0=100, T=1, r=0.1, sigma=0.25, process="GBM",
                         poly_coef=[-10, 1],
                         positive_interval=[10, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
bsmCloseForm = BSMCall(S0=100, T=1, r=0.1, sigma=0.25, K=10)
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns, f"error={abs(COSAns - bsmCloseForm.BSM_call_value()):.6f}")
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns, f"error={abs(TreeAns - bsmCloseForm.BSM_call_value()):.6f}")

# ==================Test for T==================
print("\nTest for T=0.1")
polynomial1 = Polynomial(S0=5, T=0.5, r=0.1, sigma=0.25, process="GBM",
                         poly_coef=[-10, 1],
                         positive_interval=[10, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
bsmCloseForm = BSMCall(S0=5, T=0.5, r=0.1, sigma=0.25, K=10)
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns, f"error={abs(COSAns - bsmCloseForm.BSM_call_value()):.6f}")
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns, f"error={abs(TreeAns - bsmCloseForm.BSM_call_value()):.6f}")

print("\nTest for T=50")
polynomial1 = Polynomial(S0=5, T=50, r=0.1, sigma=0.25, process="GBM",
                         poly_coef=[-10, 1],
                         positive_interval=[10, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
bsmCloseForm = BSMCall(S0=5, T=50, r=0.1, sigma=0.25, K=10)
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns, f"error={abs(COSAns - bsmCloseForm.BSM_call_value()):.6f}")
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns, f"error={abs(TreeAns - bsmCloseForm.BSM_call_value()):.6f}")

# ==================Test for r==================
print("\nTest for r=0.01")
polynomial1 = Polynomial(S0=5, T=1, r=0.01, sigma=0.25, process="GBM",
                         poly_coef=[-10, 1],
                         positive_interval=[10, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
bsmCloseForm = BSMCall(S0=5, T=1, r=0.01, sigma=0.25, K=10)
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns, f"error={abs(COSAns - bsmCloseForm.BSM_call_value()):.6f}")
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns, f"error={abs(TreeAns - bsmCloseForm.BSM_call_value()):.6f}")
print("\nTest for r=1")

polynomial1 = Polynomial(S0=5, T=1, r=1, sigma=0.25, process="GBM",
                         poly_coef=[-10, 1],
                         positive_interval=[10, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
bsmCloseForm = BSMCall(S0=5, T=1, r=1, sigma=0.25, K=10)
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns, f"error={abs(COSAns - bsmCloseForm.BSM_call_value()):.6f}")
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns, f"error={abs(TreeAns - bsmCloseForm.BSM_call_value()):.6f}")

# ==================Test for simga==================
print("\nTest for sigma=0.01")
polynomial1 = Polynomial(S0=5, T=1, r=0.1, sigma=0.01, process="GBM",
                         poly_coef=[-10, 1],
                         positive_interval=[10, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
bsmCloseForm = BSMCall(S0=5, T=1, r=0.1, sigma=0.01, K=10)
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns, f"error={abs(COSAns - bsmCloseForm.BSM_call_value()):.6f}")
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns, f"error={abs(TreeAns - bsmCloseForm.BSM_call_value()):.6f}")

print("\nTest for sigma=0.9")
polynomial1 = Polynomial(S0=5, T=1, r=0.1, sigma=0.9, process="GBM",
                         poly_coef=[-10, 1],
                         positive_interval=[10, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
bsmCloseForm = BSMCall(S0=5, T=1, r=0.1, sigma=0.9, K=10)
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns, f"error={abs(COSAns - bsmCloseForm.BSM_call_value()):.6f}")
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns, f"error={abs(TreeAns - bsmCloseForm.BSM_call_value()):.6f}")

print("\n\n")

print("=====Polynomial========")
# ==================Test for S==================
print("Test for S0=1")
polynomial1 = Polynomial(S0=1, T=1, r=0.1, sigma=0.25, process="GBM",
                         poly_coef=[840, -638, 179, -22, 1],
                         positive_interval=[exp(-9), 4, 5, 6, 7, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns)
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns)

print("\nTest for S0=100")
polynomial1 = Polynomial(S0=100, T=1, r=0.1, sigma=0.25, process="GBM",
                         poly_coef=[840, -638, 179, -22, 1],
                         positive_interval=[exp(-9), 4, 5, 6, 7, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns)
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns)

# ==================Test for T==================
print("\nTest for T=0.5")
polynomial1 = Polynomial(S0=1, T=0.5, r=0.1, sigma=0.25, process="GBM",
                         poly_coef=[840, -638, 179, -22, 1],
                         positive_interval=[exp(-9), 4, 5, 6, 7, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns)
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns)

print("\nTest for T=50")
polynomial1 = Polynomial(S0=1, T=50, r=0.1, sigma=0.25, process="GBM",
                         poly_coef=[840, -638, 179, -22, 1],
                         positive_interval=[exp(-9), 4, 5, 6, 7, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns)
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns)

# ==================Test for r==================
print("\nTest for r=0.01")
polynomial1 = Polynomial(S0=1, T=1, r=0.01, sigma=0.25, process="GBM",
                         poly_coef=[840, -638, 179, -22, 1],
                         positive_interval=[exp(-9), 4, 5, 6, 7, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns)
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns)

print("\nTest for r=1")
polynomial1 = Polynomial(S0=1, T=1, r=1, sigma=0.25, process="GBM",
                         poly_coef=[840, -638, 179, -22, 1],
                         positive_interval=[exp(-9), 4, 5, 6, 7, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns)
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns)

# ==================Test for simga==================
print("\nTest for sigma=0.01")
polynomial1 = Polynomial(S0=1, T=1, r=0.1, sigma=0.01, process="GBM",
                         poly_coef=[840, -638, 179, -22, 1],
                         positive_interval=[exp(-9), 4, 5, 6, 7, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns)
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns)

print("\nTest for sigma=0.9")
polynomial1 = Polynomial(S0=1, T=1, r=0.1, sigma=0.9, process="GBM",
                         poly_coef=[840, -638, 179, -22, 1],
                         positive_interval=[exp(-9), 4, 5, 6, 7, exp(9)],
                         N=2000,
                         lower_limit=exp(-9), upper_limit=exp(9))
COSAns = timeit(polynomial1.PricingByCOSMethod)
print("COS: ", COSAns)
TreeAns = timeit(polynomial1.PricingByTree)
print("Tree: ", TreeAns)
