from math import exp, log
import numpy as np

def stockPrice(i, j):
    return S0 * (u ** (i - j) * d ** j)

def payoff(price):
    payoff = 0
    for power, coef in enumerate(poly_coef):
        payoff += coef * price ** power
    return max(payoff, 0)
def tree():
    payoff_arr = np.zeros(n + 1).T
    for j in range(n + 1):
        payoff_arr[j] = payoff(stockPrice(n, j))

    for i in reversed(list(range(n))):
        for j in range(i + 1):
            payoff_arr[j] = (p * payoff_arr[j] + (1-p) * payoff_arr[j+1]) * exp(-r*delta_t)
    return payoff_arr[0]



if __name__ == "__main__":
    # ============== Market Data ==============
    S0 = 5
    T = 1
    r = 0.1
    sigma = 0.25

    # ============== Hyperparameter ==============
    n = 1000


    # ============== Tree Parameter ===========
    delta_t = T / n
    u = exp(sigma * (delta_t ** 0.5))
    d = 1/ u
    p = (exp(r * delta_t) - d) / (u - d)

    # ============== Option Setting ==============
    poly_coef = np.array([exp(-9), 77, 82, 159])  # store the coefficients of the polynomial

    # ============== Preproccessing ==============
    x0 = log(S0)

    # ============== Pricing ==============
    print(tree())