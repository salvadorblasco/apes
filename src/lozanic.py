import math

def comb(n, k):
    if n < k:
        raise ValueError("%d must be greater of equal than %d" % (n,k))
    return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))

def glozanic(n,k,q):
    return int(comb(n,k) + comb(n%q, k%q)*comb(n//q,k//q)) / q
