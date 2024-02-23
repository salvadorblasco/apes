import itertools

import numpy as np

from libeq.cpluss import cpluss


def damping(x0, beta, stoichiometry, T, damping_factor=5.0):
    r"""Damp concentrations.

    The damping procedure is a variation of the one described by De Robertis
    et al., [#]_
    and the procedure is as follows: If in a certain stage *n*, the ratio

    .. math:: R_k^{(n)} = C_k / C_{k, calc.}^{(n)}

    for the *k*-th component lies outside the range
    :math:`\rho^{-1} < R_k^{(n)} < \rho` (:math:`\rho` is a limit chosen in
    the range :math:`1 < \rho < 10`, then the free concentration  of the
    *k*-th component is damped by the equation

    .. math:: c_{k,damped}^{(n)} = c_k^{(n)}(R_k^q)^{(n)}

    where :math:`q = |p_{ik}|^{-1}_{max}` (i.e. *q* is the reciprocal of
    the largest stoichiometric coefficient of the species containing
    the *k*-th component). This procedure is applied to the component for
    which :math:`\ln{|R_k|}` assumes the maximum value, and is repeated until
    :math:`R_k` for all components satisfies the condiction
    :math:`\rho^{-1} < R_k^{(n)} < \rho`. Then a new Newton-Rhapson
    iteration is performed. (Notes: *concentration* here is the analytical
    concentration *T* and *c* is the free concentrations)

    Parameters:
        x0 (:class:`numpy.ndarray`): The initial guess for the free
            concentrations. It must be an (*N*, *S*)-sized array of
            floats or an *S*-sized array of floats which will be
            converted into (1, *S*)-sized array. This variable contains
            the output in the end.
        beta (:class:`numpy.ndarray`): values of equilibrium constants
        stoichiometry (:class:`numpy.ndarray`): The :term:`stoichiometry array`
        T (:class:`numpy.ndarray`): The :term:`total concentrations array`.
            It must be broadcastable with **x0**.
        damping_factor (float): The threshold that will be applied.

    Returns:
        tuple: The indices of the outliers that did not converge after a
            number of iterations.

    .. [#] *Analytica Chimica Acta* 1986, **191**, 385-398
    """
    n_species = stoichiometry.shape[1]
    mP = 1.0/np.max(np.abs(stoichiometry), axis=0)
    best = x0.size
    max_rep = n_species   # maximum number of replicas when no improvement

    # Compute order of components.
    # TODO there is no need to recalculate this over and over again
    _, c = np.nonzero(stoichiometry)
    o = [c.tolist().count(i) for i in range(n_species)]
    order = [o.index(i) for i in sorted(o)]

    # Compute limits based on damping_factor and sign of T
    # TODO there is no need to recalculate this over and over again
    CLIM1, CLIM2 = _clims(T, damping_factor)

    for i in itertools.cycle(order):
        Tcalc = _damping1(x0, i, beta, stoichiometry, T, damping_factor, mP)
        nout, q = _outliers(CLIM1, CLIM2, Tcalc)

        if q == 0:
            return None

        if q < best:
            best = q
            max_rep += 1
        else:
            max_rep -= 1

        if max_rep == 0:
            outl = np.unique(np.nonzero(~nout)[0])
            return outl


def _damping1(x0, species, beta, stoichiometry, T, damping_factor=5.0,
               mP=None):
    if mP is None:
        mP = 1.0/np.max(np.abs(stoichiometry), axis=0)
    cext = cpluss(x0, beta, stoichiometry)
    Tcalc = x0 + np.dot(cext, stoichiometry)
    R = T[:, species]/Tcalc[:, species]

    # Negative values lead to disaster. Ignore them by making R=1.
    R[R <= 0.0] = 1.0

    w = np.logical_or(R < 1/damping_factor, R > damping_factor)
    x0[np.nonzero(w), species] *= R[w]**mP[species]
    return Tcalc


def _outliers(lowerlim, upperlim, tcalc):
    nout = np.logical_and(tcalc > lowerlim, tcalc < upperlim)
    q = np.count_nonzero(~nout)
    return nout, q


def _clims(T, damping_factor):
    CLIM1 = T/damping_factor
    CLIM2 = T*damping_factor
    w = T < 0
    CLIM2[w] = CLIM1[w]
    CLIM1[w] = T[w]*damping_factor
    return CLIM1, CLIM2
