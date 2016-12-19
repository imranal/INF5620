import numpy as np
import time
def primes(kmax):
    t0 = time.clock()
    p = np.ones(2000)
    result = []
    if kmax > 2000:
        kmax = 2000
    k = 0
    n = 2
    while k < kmax:
        i = 0
        while i < k and n % p[i] != 0:
            i = i + 1
        if i == k:
            p[k] = n
            k = k + 1
            result.append(n)
        n = n + 1
    t1 = time.clock()

    return result, t1-t0

if __name__ == "__main__":
    result, t = primes(2000)
    print "time taken to finish calculations for 2000 prime numbers : %f seconds"%t
