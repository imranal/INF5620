import time
def primes(int kmax):
    t0 = time.clock()
    cdef int n, k, i
    cdef int p[10000]
    result = []
    if kmax > 10000:
        kmax = 10000
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
