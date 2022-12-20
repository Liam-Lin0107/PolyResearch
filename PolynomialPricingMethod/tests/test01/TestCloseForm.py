from PolynomialPricingMethod.CloseFormMethod import BSMCall

S0 = 100
K = 100
T = 1
sigma = 0.25
r = 0.03
call = BSMCall(S0=S0, T=T, r=r, sigma=sigma, K=K)
print(call.getValue())