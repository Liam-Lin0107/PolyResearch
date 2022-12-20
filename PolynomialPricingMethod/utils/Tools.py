import time


def timeit(f):
    start = time.perf_counter()
    result = f()
    end = time.perf_counter()
    print("CPU time  is", f"{end - start:.6f}")
    return result
