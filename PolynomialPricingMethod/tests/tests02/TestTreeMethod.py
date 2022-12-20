from PolynomialPricingMethod.TreeMethod import PolyByTree
from PolynomialPricingMethod.utils.Tools import timeit
import math


def GBMTest():
    print("GBM")
    # ===================== Test: GBM Call ======================
    S0 = 100
    r = 0.05
    T = 0.5
    sigma = 0.2

    poly_coeff = [-100, 1]

    for N in [640, 1280, 2560, 5120, 10240]:
        call = PolyByTree(S0=S0, T=T, r=r, sigma=sigma, poly_coeff=poly_coeff, N=N)
        ans = timeit(call.getValue)
        print(f"N={N}: value: {ans:.6f}")

    # ===================== Test: GBM Polynomial right up  ======================
    S0 = 50
    T = 0.5
    r = 0.1
    sigma = 0.4

    poly_coeff = [-2725, -100, 1]

    for N in [640, 1280, 2560, 5120, 10240]:
        call = PolyByTree(S0=S0, T=T, r=r, sigma=sigma, poly_coeff=poly_coeff, N=N)
        ans = timeit(call.getValue)
        print(f"N={N}: value: {ans:.6f}")

    # ===================== Test: GBM Polynomial left up  ======================
    S0 = 110
    T = 0.5
    r = 0.03
    sigma = 0.4

    poly_coeff = [947.1, -30.164, 0.309, -0.001]
    for N in [640, 1280, 2560, 5120, 10240]:
        call = PolyByTree(S0=S0, T=T, r=r, sigma=sigma, poly_coeff=poly_coeff, N=N)
        ans = timeit(call.getValue)
        print(f"N={N}: value: {ans:.6f}")

    # ===================== Test: GBM Polynomial both up  ======================
    S0 = 15
    T = 1
    r = 0.1
    sigma = 0.3

    poly_coeff = [44.235, -39.474, 5.4793, -0.2358, 0.0031]
    for N in [640, 1280, 2560, 5120, 10240]:
        call = PolyByTree(S0=S0, T=T, r=r, sigma=sigma, poly_coeff=poly_coeff, N=N)
        ans = timeit(call.getValue)
        print(f"N={N}: value: {ans:.6f}")


if __name__ == "__main__":
    # thread_heston = threading.Thread(target=HestonTest, name="Heston")
    # thread_MJD = threading.Thread(target=MJDTest, name="MJD")
    # thread_KJD = threading.Thread(target=KJDTest, name="KJD")
    # thread_SVJ = threading.Thread(target=SVJTest, name="SVJ")
    # thread_heston.start()
    # time.sleep(0.2)
    # thread_MJD.start()
    # time.sleep(0.2)
    # thread_KJD.start()
    # time.sleep(0.2)
    # thread_SVJ.start()

    GBMTest()
