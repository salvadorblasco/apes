#!/usr/bin/python3

import ctypes
import numpy as np

_eq = ctype.LoadLibrary('./eq.so')
_eq.newton_raphson.argtypes()


def NewtonRaphson(x0, beta, stoichiometry, analyticalc, max_iterations=1000,
                  threshold=1e-10, damping=True, forcer=True, scaling=True,
                  step_limiter=False, zero_offdiag=False, logc=False,
                  log=False, **kwargs):

    _x0 = x0.ctypes

    eq.newton_raphson()

    return np.ctypeslib.as_array(x0)
