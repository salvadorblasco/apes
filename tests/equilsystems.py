import sys
sys.path.append('../src/')

import numpy as np


def ccalc(T, B, Kw=1e-14):
    C = np.zeros((T.shape[0], len(B)+3))
    H = np.fromiter(hcalc(T, B, Kw), dtype=np.float)
    L = T[:, 0]/(1+np.sum(B[np.newaxis, :]*H[:, np.newaxis]**np.arange(1, 1+len(B)), axis=1))
    C[:, 0] = L
    C[:, 1] = H
    for i, b in enumerate(B, start=2):
        C[:, i] = b*L*H**(i-1)
    C[:, -1] = Kw/H
    return C


def hcalc(T, B, Kw):
    for t in T:
        coeffs = hcoeffs(t, B, Kw)
        r = np.roots(tuple(reversed(coeffs)))
        h = filter_h(r)
        yield h


def hcoeffs(T, B, Kw):
    coeffs = (len(B) + 3)*[0.0]
    for i in range(len(B)):
        coeffs[i+3] += B[i]
        coeffs[i+2] += T[0]*(i+1)*B[i]
        coeffs[i+2] -= T[1]*B[i]
        coeffs[i+1] -= Kw*B[i]
    coeffs[2] += 1
    coeffs[1] -= T[1]
    coeffs[0] -= Kw
    return coeffs


def filter_h(roots):
    try:
        be_real = [float(x) for x in roots if np.isreal(x)]
    except np.exceptions.ComplexWarning:
        pass

    be_positive = [x for x in be_real if x >= 0.0]
    if len(be_positive) == 1:
        retv = be_positive[0]
    else:
        retv = min(be_positive)
    return float(retv)


def _setup_generic_protonation(TL=0.001, TH=0.010, V0=20.0):
    import libaux
    T0 = np.array([TL, TH])
    buret = np.array([0.0, -0.1])
    N = 20
    V = np.linspace(0.0, 1.2, N)
    return np.array(libaux.build_T_titr2(T0, buret, V0, V))


def generic_monoprotic(logK=9.75, logKw=-13.73):
    P = np.array([[1, 1], [0, -1]])
    B = 10**np.array([logK, logKw])
    T = _setup_generic_protonation(TH=0.05, V0=10.0)
    real_c = ccalc(T, np.array([B[0]]), B[1])
    return real_c, P, B, T


def generic_zeroprotic(logKw = -13.73):
    P = np.array([[-1]])
    B = 10**np.array([logKw])
    T = _setup_generic_protonation(TH=0.06)[:,1]
    Kw = 10**logKw
    real_c = (T + np.sqrt(T**2 + 4*Kw))/2
    return real_c, P, B, T
