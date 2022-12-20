import numpy as np


def merton_jump_paths(S, T, r, sigma, lam, m, v, steps, Npaths):
    size = (steps, Npaths)
    dt = T / steps
    poi_rv = np.multiply(np.random.poisson(lam * dt, size=size),
                         np.random.normal(m, v, size=size)).cumsum(axis=0)
    geo = np.cumsum(((r - sigma ** 2 / 2 - lam * (m + v ** 2 * 0.5)) * dt +
                     sigma * np.sqrt(dt) * np.random.normal(size=size)), axis=0)

    return np.exp(geo + poi_rv) * S


S = 100  # current stock price
T = 0.5  # time to maturity
r = 0.05  # risk free rate
m = 0.01  # meean of jump size
v = 0.02  # standard deviation of jump
lam = 140 # intensity of jump i.e. number of jumps per annum
steps = 252  # time steps
Npaths = 1000000  # number of paths to simulate
sigma = 0.2  # annaul standard deviation , for weiner process

# S0 = 100
# r = 0.05
# T = 0.5
# sigma = 0.2

# 資料來自於於 1998 Pitfalls in Estimating Jump-Diffusion Models p10
# jump_intensity = 140
# jump_mean = 0.01
# jump_var = 0.02 ** 2

poly_coef = [-100, 1]


def _Payoff(price):
    payoff = 0
    for power, coef in enumerate(poly_coef):
        payoff += coef * price ** power
    return max(payoff, 0)


values = []
for i in range(20):
    j = merton_jump_paths(S, T, r, sigma, lam, m, v, steps, Npaths)

    mean = np.mean(list(map(_Payoff, j[:, -1])))
    value = mean * np.exp(-r * T)
    values.append(value)

print(np.mean(values) - 2 * np.std(values), np.mean(values) + 2 * np.std(values))
print(np.mean(values))
