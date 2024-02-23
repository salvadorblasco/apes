#!/usr/bin/python3
# -*- coding: utf-8 -*-

# This file should contain the routines dedicated for constant fitting

import numpy as np
import scipy as sp
from scipy.optimize import leastsq

import libeq
import libaux


T0 = np.array([[0.021, 0.01, 0.003, 0.0015, 0.0006, 0.0004, 0.0002, 0.0001, 0.00005]]).T
delta_obs = np.array([
    [9.62,	9.02,	8.99,	8.64,	8.27,	8.52,	8.38,	8.20,	7.94,	7.86],
    [9.68,	9.00,	8.98,	8.64,	8.37,	8.49,	8.37,	8.19,	8.03,	7.86],
    [9.75,	8.99,	8.98,	8.64,	8.47,	8.47,	8.36,	8.17,	8.12,	7.87],
    [9.76,	8.98,	8.97,	8.64,	8.48,	8.46,	8.36,	8.17,	8.14,	7.87],
    [9.76,	8.98,	8.97,	8.64,	8.49,	8.46,	8.36,	8.17,	8.14,	7.87],
    [9.76,	8.98,	8.97,	8.64,	8.49,	8.46,	8.35,	8.17,	8.14,	7.87],
    [9.76,	8.98,	8.97,	8.64,	8.49,	8.46,	8.35,	8.17,	8.15,	7.87],
    [9.76,	8.98,	8.97,	8.64,	8.49,	8.46,	8.36,	8.17,	8.15,	7.87],
    [9.76,	8.98,	8.97,	8.64,	8.49,	8.46,	8.35,	8.17,	8.15,	7.87]])

B0 = np.array([ 1.0 ])
B_mask = np.array([ 1 ])
P = np.array([[ 2 ]])

d0 = np.array([\
    [9.56,	9.10,	9.00,	8.64,	8.19,	8.66,	8.45,	8.27,	7.85,	7.80],
    [9.76,	8.98,	8.97,	8.64,	8.49,	8.46,	8.36,	8.17,	8.15,	7.87]])

d_mask = np.array([\
    [1,	1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0,	0, 0, 0, 0, 0, 0, 0, 0, 0]])

X = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

def numeric_jacobian(f, C, B, P, T, eps=1e-6): 
    """from http://old.nabble.com/calculating-numerical-jacobian-td20506078.html
    Evaluate partial derivatives of f(u) numerically. 

    :note: 
        This routine is currently naive and could be improved. 

    :returns: 
        (*f.shape, *u.shape) array ``df``, where df[i,j] ~= (d f_i / u_j)(u) 
    """ 

    f0 = f(C, B, P, T) # asarray: because of matrices 

    u_shape = C.shape 
    nu = u_shape[0]

    f_shape = f0.shape
    nf = np.prod(f_shape) 

    df = np.empty([nf, nu]) 
    
    for k in range(nu): 
        du = np.zeros(nu) 
        du[k] = max(eps*abs(C.flat[k]), eps) 
        f1 = f(C + np.reshape(du, u_shape), B, P, T)
        df[:,k] = np.reshape((f1 - f0) / eps, [nf]) 

    df.shape = f_shape + u_shape 
    return df 

def jacobian(p, B, B_mask, P, d, d_mask, T, X):
    """
           / ∂δ_1/∂α1 ... ∂δ_1/∂αnd \ 
           | ∂δ_2/∂α1 ... ∂δ_2/∂αnd |
       J = |    ...   ...    ...    |
           \ ∂δ_N/∂α1 ... ∂δ_N/∂αnd /

        
    """

    f = calc_fj(B, P, T)
    E, S = P.shape
    N = T.shape[0]
    n = d.shape[1]

    #nonlocal E, S, N, n

    assert p.ndim == 1

    pidx = 0
    J = np.zeros((N*n, len(p)))

    for j in range(E):
        if B_mask[j] == PARM_REFINE:
            for i1 in range(N):
                for i2 in range(n):
                    row = ((i1+1)*(i2+1)-1)
                    _ = (d[X[k],i2]*f[i1,X[k]]/P[j,X[k]] for k in range(n))
                    J[row,pidx] = -sum(_) / p[pidx]
            pidx = pidx + 1

    bref = pidx
    for j1 in range(n):
        for j2 in range(E+S):
            if d_mask[j2,j1] == PARM_REFINE:
                for i in range(N):
                    row = ((i+1)*(j1+1)-1)
                    J[row,pidx] = f[i,j2]
                pidx = pidx + 1

    return J

def calc_fj(B, P, T):
    "Calculates the fractional populations"

    C = libeq.consol(B.reshape(E,1), P, T)
    assert C.ndim == 2
    assert C.shape == (N, E+S)

    xT = np.empty((N, E+S))
    for i in range(N):
        for j in range(S):
            #if x < S:
            xT[i,j] = C[i,j]/T[i,j]
        for j in range(S, S+E):
            # k is the index of the first non zero element from P
            k = P[j-S,:].nonzero()[0][0]
            xT[i,j] = P[j-S,k]*C[i,j]/T[i,k]
    return xT

print(jacobian(p, B0, B_mask, P, d0, d_mask, T0, X))

#def nmrfit(delta_obs, T, P, B0, B_mask, d0, d_mask, X):
#    """See from C. Frassineti, S. Ghelli, P. Gans, A. Sabatini,
#        M. S. Moruzzi and A. Vacca, Anal. Biochem. 1995, 231, 374-382
#
#    Function to minimise:
#        U = Σ_j w_j(δ_calc - δ_obs)
#        δ_calc = Σ_j(f_j.d_j) where 
#            d_j chemical shift of a particular nucleus in the various species
#                present.
#        f_j = x_j . C_j/T_x        
#            T_x is the total concentration of reagent containing the nucleus
#                under consideration
#            x_j is the stoichiometric coefficient X in the j_th species
#    """
#
#    # N [int] is the number of points in titration
#    # n [int] is the number of nuclei
#    # B0[i=0..E-1] {float} is the initial guess for the equilibrium constant
#    #    in logarithmic units
#    # B_mask[i=0..E-1] {int} is a boolean matrix indicating which constants 
#    #   are to be refined or kept constant.
#    # f[i=0..N-1,j=0..S+E] is the fractional population of species j in point i
#    # d[i=0..n-1, j=0..S+E] shift of the nuclei i in the species j
#    # d_mask[i=0..n-1, j=0..S+E] {int} array of the same shape than d containing
#    #   the information on which d values are known, unknown of refined / unrefined
#    # T[i=0..N-1,j=0..S] Total concentrations of species j in point i 
#    # C[i=0..N-1,j=0..S+E] Concentration of species j in point i
#    # X[i=0..n] Index of the species to which the nucleus i belongs
#
#    # Constant section
#    PARM_CONSTANT = 0
#    PARM_REFINE   = 1
#    PARM_UNKNOWN  = 2
#
#    # Input check section
#    assert isinstance(delta_obs, np.ndarray)
#    N, n = delta_obs.shape
#
#    E, S = check_PB(P, B0)
#
#    assert isinstance(B_mask, np.ndarray)
#    assert B_mask.shape == B0.shape
#
#    assert isinstance(d0, np.ndarray)
#    assert isinstance(d_mask, np.ndarray)
#    assert d0.shape == d_mask.shape
#
#    assert isinstance(T, np.ndarray)
#    assert T.ndim == 2
#    assert T.shape == (N, S)
#
#    assert isinstance(X, np.ndarray)
#    assert X.ndim == 1
#    assert len(X) == n
#
#    # Define nonlocal variables
#    C = np.empty((N, E+S))
#    d = np.empty_like(d0)
#    f = np.empty_like(C)
#
#    def jacobian(p):
#        """
#               / ∂δ_1/∂α1 ... ∂δ_1/∂αnd \ 
#               | ∂δ_2/∂α1 ... ∂δ_2/∂αnd |
#           J = |    ...   ...    ...    |
#               \ ∂δ_N/∂α1 ... ∂δ_N/∂αnd /
#
#            
#        """
#
#        nonlocal f
#        nonlocal d
#        nonlocal E, S, N, n
#
#        assert p.ndim == 1
#
#        pidx = 0
#        J = np.zeros((N*n, len(p)))
#
#        for j in range(E):
#            if B_mask[j] == PARM_REFINE:
#                for i1 in range(N):
#                    for i2 in range(n):
#                        row = ((i1+1)*(i2+1)-1)
#                        _ = (d[X[k],i2]*f[i1,X[k]]/P[j,X[k]] for k in range(n))
#                        J[row,pidx] = -sum(_) / p[pidx]
#                pidx = pidx + 1
#
#        bref = pidx
#        for j1 in range(n):
#            for j2 in range(E+S):
#                if d_mask[j2,j1] == PARM_REFINE:
#                    for i in range(N):
#                        row = ((i+1)*(j1+1)-1)
#                        J[row,pidx] = f[i,j2]
#                    pidx = pidx + 1
#
#        return J
#
#    def calc_fj(B, P, T):
#        "Calculates the fractional populations"
#
#        nonlocal C
#
#        C = libeq.consol(B.reshape(E,1), P, T)
#        assert C.ndim == 2
#        assert C.shape == (N, E+S)
#
#        xT = np.empty((N, E+S))
#        for i in range(N):
#            for j in range(S):
#                #if x < S:
#                xT[i,j] = C[i,j]/T[i,j]
#            for j in range(S, S+E):
#                #else:
#                # k is the index of the first non zero element from P
#                k = P[j-S,:].nonzero()[0][0]
#                xT[i,j] = P[j-S,k]*C[i,j]/T[i,k]
#        return xT
#
#    def delta_calc():
#        nonlocal d, f
#        assert isinstance(f, np.ndarray)
#        assert f.ndim == 2
#        assert f.shape == (N, E+S)
#
#        assert isinstance(d, np.ndarray)
#        assert d.ndim == 2
#        assert d.shape == (E+S, n)
#
#        # f is N × (S+E)
#        # d is (S+E) × n
#        # delta is N × n
#
#        ret = np.dot(f, d)
#        assert isinstance(ret, np.ndarray)
#        assert ret.ndim == 2
#        assert ret.shape == (N, n), ret.shape
#
#        return ret
#
#    def bundle_p(B, _d):
#        # TODO change to work as generator, for efficiency
#        assert B.shape == B_mask.shape
#        assert _d.shape == d_mask.shape
#
#        p = [ B[i] for i in range(len(B)) if B_mask[i] == PARM_REFINE ]
#        p.extend([ _d[i,j] \
#            for i in range(_d.shape[0]) \
#            for j in range(_d.shape[1]) \
#            if d_mask[i,j] == PARM_REFINE ])
#        return np.array(p)
#
#    def unboundle_p(p):
#        # TODO there must be a better way than using a counter 'k'
#        k = 0
#        B = np.empty_like(B0)
#        for i in range(len(B_mask)):
#            if B_mask[i] == PARM_CONSTANT:
#                B[i] = B0[i]
#                continue
#            if B_mask[i] == PARM_REFINE:
#                B[i] = p[k]
#                k = k + 1
#
#        d = np.empty_like(d0)
#
#        for i in range(d_mask.shape[0]):
#            for j in range(d_mask.shape[1]):
#                if d_mask[i,j] == PARM_CONSTANT:
#                    d[i,j] = d0[i,j]
#                    continue
#                if d_mask[i,j] == PARM_REFINE:
#                    d[i,j] = p[k]
#                    k = k + 1
#        
#        assert k == len(p), "k = %d, len(p) = %d" % (k, len(p))
#        return B, d
#
#    def f_obj(p):
#        nonlocal d
#        _B, d = unboundle_p(p)
#        assert isinstance(_B, np.ndarray)
#        assert _B.ndim == 1
#        assert len(_B) == E
#        
#        assert isinstance(d, np.ndarray)
#        assert d.ndim == 2
#        assert d.shape == (E+S, n)
#
#        nonlocal f
#        f = calc_fj(_B, P, T)
#        _res = delta_obs - delta_calc()
#        return _res.flatten()
#
#    p0 = bundle_p(B0, d0)
#    print(leastsq(f_obj, p0)) #, Dfun=jacobian))
#    #bfit, bcov, inf, msg, ier = leastsq(f_obj, p0, Dfun=jacobian)
#    #print(f_obj(p0))
#    #print(jacobian(p0))
#    #return bfit
#
