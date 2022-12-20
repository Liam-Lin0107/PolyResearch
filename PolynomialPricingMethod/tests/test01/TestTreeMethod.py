from PolynomialPricingMethod.TreeMethod import PolyByTree
S0 = 100
r = 0.03
T = 1
sigma = 0.25

poly_coef = [-100, 1]

call = PolyByTree(S0=S0, T=T, r=r, sigma=sigma, poly_coeff=poly_coef, N=10000)
print(call.getValue())
