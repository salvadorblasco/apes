"""FILE: nr.py

Contains hard core numeric recipes.
"""

import numpy as np

import excepts
from libeq.cpluss import cpluss
from libeq.jacobian import jacobian
from libeq.fobj import fobj
from libeq.damping import damping as fdamping


__version__ = '0.2'


def NewtonRaphson(x0, beta, stoichiometry, T, max_iterations=1000,
                  threshold=1e-10, damping=True, forcer=True, scaling=False,
                  step_limiter=False, zero_offdiag=False, logc=False,
                  log=False, panic=True, mask_flag=False, **kwargs):
    r"""Solve the set of equations J·δc = -F using Newton-Raphson's method.

    Given an initial guess **x0** and the parameters of the titration, this
    function uses the Newton-Raphson method to solve the free
    concentrations.  It includes several
    convergence tecniques such as scaling, damping and step limiting. For the
    original source of damping and scaling see De Robertis, *et al.*,
    [#f1]_ and also Del Piero, *et al.* [#f2]_

    Parameters:
        x0 (:class:`numpy.ndarray`): initial guess for iterations. If *logc* is
            True, the natural logarithm of the concentrations is expected to be
            input. Otherwise, the regular concentration are expected.
        beta (:class:`numpy.ndarray`): values of equilibrium constants
        stoichiometry (:class:`numpy.ndarray`): stoichiometric coefficients
        T (:class:`numpy.ndarray`): total concentration values
        scaling (bool, optional): Whether to use or not scaling techniques.
            False by default. See :func:`DRscaling`.
        damping (bool, optional): Whether to use or not damping techniques.
            False by default. See :func:`DRdamping`.
        step_limiter (bool, optional[True]): Whether to use step limiter or
            not.  Step limiter attempts to avoid a step which leads to negative
            concentrations. It searches the increments which will result in
            negative concentrations and replaces the step
            :math:`x_n = x_{n-1} + dx` for
            :math:`x_n = x_{n+1} \cdot e^{dx}`
            which always gives a positive concentration as a result.
            It can jointly be used with the forcer.
        forcer (bool, optional[False]): Whether to use forcer techniques.
            False by default.
        threshold (float): threshold criterium for convergence
        max_iterations (int): maximum number of iterations allowed
        zero_offdiag (bool): If this option is applied, all non-diagonal
            elements of the jacobian are set to zero. This option is useful
            when only an estimate of the free concentrations is wanted.
        max_damps (int):  default, 2. maximum number of dampimgs allowed.
            It applies only if damping = True.
        logc (bool): Fit the logarithmic value of the concentration
        log (bool): The handle for logging
        do_iterations (int, optional): Perform exatly *do_iterations*
            iterations and do not test for convergence.
        panic (bool, optional, default True): If convergence fails, dump to a
            file the data to debug.
    Returns:
        :class:`numpy.ndarray`: An array containing the values of the free
            concentrations.
    Raises:
        :class:`RuntimeError`: if the concentrations cannot be calculated
        :class:`excepts.TooManyIterations`: if the maximum number of iterations
            is reached.

    .. warning:: This function is the very core of the program. Do not change
        the default parameters unless you know what you are doing.

    .. [#f1] *Analytica Chimica Acta* 1986, **191**, 385-398
    .. [#f2] Del Piero, *et al.*: *Annali di Chimica* 2006, 96.
    """

    def _panic_save():
        if panic:
            np.savez_compressed('consol_panic.npz', free_concentration=x0,
                                beta=beta, stoichiometry=stoichiometry,
                                analytc=T)

    if zero_offdiag and scaling:
        raise ValueError("Options scaling and zero_offdiag are not" +
                         "compatible with each other")

    if 'do_iterations' in kwargs:
        do_iterations = kwargs['do_iterations']
        if not isinstance(do_iterations, int):
            raise TypeError('do_iteration must be a positive int.')
        if do_iterations < 1:
            raise ValueError('do_iteration must be a positive int.')
    else:
        do_iterations = None

    iterations = 0
    n_species = stoichiometry.shape[1]

    # copy x0 so that it is not modified outside this function.
    x = np.copy(x0)
    # ------ check input --------
    # libaux.assert_array_dim(2, x0)
    x, T = np.atleast_2d(x, T)

    # ------ main loop ----------
    while True:
        if mask_flag and not logc:
            np.ma.masked_less_equal(x, 0.0, copy=False)

        c0 = cpluss(x, beta, stoichiometry, full=True, logc=logc)

        if logc:
            _c = np.exp(c0)
        else:
            _c = c0

        if mask_flag:
            F = np.ma.masked_invalid(fobj(_c, stoichiometry, T))
            J = np.ma.masked_invalid(jacobian(_c, stoichiometry, logc=logc))
            assert not np.any(np.isnan(J))
            assert not np.any(np.isnan(F))
        else:
            F = fobj(_c, stoichiometry, T)
            J = jacobian(_c, stoichiometry, logc=logc)
            msg = 'could not calculate %s (iteration %d)'
            if np.any(np.isnan(J)):
                _panic_save()
                msg2 = msg % ('jacobian', iterations)
                raise excepts.FailedCalculateConcentrations(msg2, x)
            if np.any(np.isnan(F)):
                msg2 = msg % ('residuals', iterations)
                raise excepts.FailedCalculateConcentrations(msg2, x)

        # TODO This should be deleted when debug is not necessary
        # if not iterations % 50:
        #     print('chisq(it:', iterations, 'n:', T.shape[0], ') = ',
        #           np.sum(F**2), np.max(np.abs(F)))

        if zero_offdiag:
            J *= np.eye(n_species)  # zerom

        if scaling:
            d = DRScaling(J, F)
            dx = np.linalg.solve(J, -F) / np.sqrt(d)
        else:
            dx = np.linalg.solve(J, -F)

        # np.nan_to_num(dx)

        if forcer:
            # TOLX = 1e-5
            # #!
            # tmp = x0/dx
            # tmp[np.isnan(tmp)] = 0.0
            # stpmax = np.max(-tmp, axis=1)
            # # stpmax = np.max(-x0/dx, axis=1)
            # dx *= stpmax[:, np.newaxis]/(1.0+TOLX)
            # slope = -np.sum(F**2, axis=1)
            # g2 = generate_data2(beta, stoichiometry, T, x, dx)
            l, _ = linesearch3(x, dx, beta, stoichiometry, T)
            # # l, _ = linesearch(x0, dx, slope, g2)
            x += l[:, None]*dx

            # ## linesearch1 ##
            # slope = -np.sum(F**2, axis=1)
            # g2 = generate_data2(beta, stoichiometry, T, x, dx)
            # l, _, x = linesearch1(x0, dx, slope, g2)
        else:
            if step_limiter:
                x = limit_step(x, dx)
            else:
                x += dx

        iterations += 1
        if do_iterations:
            if iterations >= do_iterations:
                break
        else:
            # -------- Test for convergence ---------
            # if np.all(np.abs(dx/x0) < threshold):
            if np.all(np.abs(F) < threshold):
                break

        if damping:
            fdamping(x, beta, stoichiometry, T)

        if iterations > max_iterations:
            raise excepts.TooManyIterations('too many iterations', x)

        if mask_flag:
            np.ma.masked_invalid(x, copy=False)
            np.ma.masked_less_equal(x, 0.0, copy=False)
            np.ma.mask_rows(x)

    return x


# def linesearch2(x0, dx, J0, F0, F1, step_limiter=True):
#     r"""2-point parabolic line search along Newton direction.
#
#     Implementation of 2-point parabolic linesearch. This is a variation
#     of the line search for 3 points.[#]_ This line search defines a function
#     :math:`f=\frac12F\cdot F` which is to be minimized. Then we define a
#     parameter λ that 0<λ<1 which is the fraction of the Newton step and
#     then another function *g* which is function of the fractional is defined
#     so that
#
#     .. math:: g(\lambda) = f(x_0 + \lambda \delta x)
#
#     We model :math:`g(\lambda)` as a parabolic function for which we know the
#     values of :math:`g(\lambda=0)` and :math:`g(\lambda=1)` as well as the
#     derivative :math:`g'(\lambda=0)=`
#
#     .. math:: \lambda = \frac{g_0'}{2(g_1-g_0-g_0')}
#
#     .. [#] W. H. Press, S. A. Teukolksy, W. T. Vetterling, Brian P. Flannery,
#        Numerical Recipes in C. The Art of Scientific Computing, Second
#        Edition 1997, pages 384--385.
#     """
#     slope = np.sum(np.sum(F0[..., None] * J0, axis=1)*dx, axis=1)
#     lam = np.full_like(slope, 0.1)
#     who = slope < 0.0
#     g0 = 0.5*np.sum(np.square(F0[who]), axis=1)
#     g1 = 0.5*np.sum(np.square(F1[who]), axis=1)
#     lam[who] = -0.5*slope[who] / (g1-g0-slope[who])
#     lam[lam > 1.0] = 1.0
#     if step_limiter:
#         return limit_step(x0, lam[:, None]*dx)
#     else:
#         return x0 + lam[:, None]*dx


def linesearch3(x0, dx, beta, stoichiometry, T, lmax=None, g0=None, g2=None):
    r"""Three-point parabolic line search.

    This functions implements a 3-point linesearch in the Newton direction.
    This is a variation of the line search for 2 points.[#]_ The function to
    be minimized is the same though but the approach is different and it is
    adapted to the nature of the problem of concentration solving. We define
    a function :math:`f=\frac12F\cdot F` which is to be minimized. Then we
    define a parameter λ that 0<λ<1 which is the fraction of the Newton step
    and then another function *g* which is function of the fractional is
    defined so that

    .. math:: g(\lambda) = f(x_0 + \lambda \delta x)

    We know that negative x₀ values are forbidden, therefore λ might limited
    to values lower than 1. The maximum value allowed for λ is that that makes
    any concentration equal to 0, therefore
    :math:`\lambda_{max} = -x_0/\delta` if :math:`-x_0/\delta<1`

    We model :math:`g(\lambda)` as a parabolic function for which we calculate
    the values for λ=0, λ=½λ(max) and λ=0.99λ(max).

    .. [#] W. H. Press, S. A. Teukolksy, W. T. Vetterling, Brian P. Flannery,
       Numerical Recipes in C. The Art of Scientific Computing, Second Edition
       1997, pages 384--385.
    """
    nerr = np.geterr()
    np.seterr(all='ignore')

    def g(l):
        "Auxiliary function."
        FF = fobj(cpluss(x0 + l[:, np.newaxis]*dx, beta, stoichiometry,
                         full=True),
                  stoichiometry, T)
        return 0.5*np.sum(np.square(FF), axis=1)

    if lmax is None:
        lmax = -x0/dx       # may cause division by 0
        # w = np.logical_or(lmax > 1.0, lmax < 0.0)
        # lmax[w] = 1.0
        lmax[lmax < 0.0] = 1.0
        lmax = np.min(lmax, axis=1)

    if g0 is None:
        g0 = g(np.zeros_like(lmax))

    g1 = g(lmax/2)
    x1 = x0+lmax[:, None]*dx
    x1[x1 < 0.0] = 0.0
    if g2 is None:
        g2 = g(0.99*lmax)

    b = (-g2+4*g1-3*g0)
    a = (g2-g0-b)/lmax
    lmin = -0.5*b/a         # may cause division by 0

    # In the unlikely case where a == 0.0 meaning g0 == g1 == g2 we set the
    # step halfway.
    w = a == 0.0
    if np.any(w):
        lmin[w] = lmax[w]/2

    w = lmin < 0.1*lmax     # Set minimum step as 0.1λ(max)
    lmin[w] = 0.1*lmax[w]

    w = lmin > lmax         # In the unlikely case where λ(min)>λ(max), set
    lmin[w] = 0.95*lmax[w]  # λ(min) close enough to λ(max)

    gmin = g(lmin)

    # Check g(lmin) < g0
    w2 = gmin > g0
    if beta.ndim == 1:
        _beta = beta
    else:
        _beta = beta[w2]
    if np.any(w2):
        lmin[w2], gmin[w2] = linesearch3(x0[w2], dx[w2], _beta, stoichiometry,
                                         T[w2], lmax=lmax[w2]/2, g0=g0[w2],
                                         g2=gmin[w2])
    np.seterr(**nerr)
    return lmin, gmin


def limit_step(x, dx):
    r"""Limit step.

    Given a state (**x**) and a step (**dx**), the next state is expected to be
    :math:`x+dx`. However, in some cases negative values are forbidden. This
    may happen for small values of the state. In the case of *x* being small
    we can approximate :math:`x+1 \simeq e^{-x}`

    Parameters:
        x (:class:`numpy.ndarray`): The state. It must be 1D.
        dx (:class:`numpy.ndarray`): The step. It must have  the same length
            than *x*.

    """
    if len(x) != len(dx):
        raise ValueError('both arguments must have the same size')
    who = (x + dx) > 0.0
    # return np.where(who, x+dx, x*np.exp(dx))
    return np.where(who, x+dx, x + (1 - np.exp(dx)))


def DRScaling(J, F):
    """Apply scaling to both jacobian and objective function.

    Applies scaling technique according to De Robertis, *et al.* [#f1]_
    The scaling technique overcomes the problem of divergence
    when the jacobian (here called G) matrix is near singular.
    "... scaling was applied to matrix *G* and to vector *e*
    (residuals, :math:`e = C_{k, calcd} - C_k`)
    according to the equations
    :math:`g_{kj}^* = g_{kj}(g_{kk}g_{jj})^{-1/2}`
    and
    :math:`e_k^* = e_kg_{kk}^{-1/2}` where :math:`g_{kj}^*`
    and :math:`e_k^*` are the elements of the scaled matrix and vector
    respectively."

    .. math:: J^*_{kj} = J_{kj}(J_{kk}J_{jj})^{-1/2}
    .. math:: F^*_{k} = F_{k}J_{kk}^{-1/2}

    Parameters:
        J (:class:`numpy.ndarray`): jacobian array, which will be modified.
            It can be of any dimensionality provided that the last two are
            of the same size.
        F (:class:`numpy.ndarray`): residuals array, which will be modified.
            If must have one dimmension less than J and the rest of the
            axes be of the same size than J.

    Returns:
        :class:`numpy.ndarray`: The diagonal of the jacobian, :math:`J_{ii}`
            to scale back the result.
    """
    # s = J.shape
    # i = np.indices(s)
    # assert s[-2] == s[-1]
    # d = (J[i[-2] == i[-1]]).reshape(s[:-1])
    d = J[np.arange(J.shape[0])[:, None], np.eye(J.shape[1], dtype=np.bool)]
    J /= np.sqrt(d[..., np.newaxis]*d[..., np.newaxis, :])
    F /= np.sqrt(d)
    return d


# DEPRECATED
def linesearch(x0, dx, slope, gfunc):
    if np.any(slope >= 0.0):
        raise ValueError('Negative slope values.')

    ALF = 1e-4
    TOLX = 1e-7

    fold = gfunc(np.zeros(dx.shape[0]))

    lmin = TOLX / np.max(np.absolute(dx/x0), axis=1)
    assert lmin.shape == (x0.shape[0], )

    lam = np.ones_like(lmin, dtype=np.float)
    lam2 = np.zeros_like(lam)
    f2 = np.zeros_like(fold)
    first = True

    while True:
        # x1 = x0+lam[:, None]*dx
        f = gfunc(lam)

        done_on_x = lam < lmin
        done_on_f = f < (fold+ALF*lam*slope)
        done = np.logical_or(done_on_x, done_on_f)
        if np.all(done):
            break

        if first:
            tmplam = -0.5*slope / (f-fold-slope)
            first = False
        else:
            rhs1 = f-fold-lam*slope
            rhs2 = f2-fold-lam2*slope
            a = (rhs1/lam**2-rhs2/lam2**2)/(lam-lam2)
            b = (-lam2*rhs1/lam**2+lam*rhs2/lam2**2)/(lam-lam2)

            w = conditions(~done, a == 0)
            tmplam[w] = -slope[w]/(2.*b[w])

            disc = b**2-3.*a*slope
            w = conditions(~done, disc < 0.0)
            tmplam[w] = 0.5*lam[w]

            w = conditions(~done, disc >= 0.0, b <= 0.0)
            tmplam[w] = (-b[w]+np.sqrt(disc[w]))/(3.*a[w])

            w = conditions(~done, disc >= 0.0, b > 0.0)
            tmplam[w] = -slope[w]/(b[w]+np.sqrt(disc[w]))

            w = conditions(~done, tmplam > 0.5*lam)
            tmplam[w] = 0.5*lam[w]

        lam2[~done] = lam[~done]
        f2[~done] = f[~done]
        lam[~done] = np.where(tmplam[~done] > 0.1*lam[~done],
                              tmplam[~done], 0.1*lam[~done])

    return lam, f


def generate_data2(beta, stoich, analyticalc, x0, dx):
    # From sandbox/linsrch.py
    def g(lam):
        x1 = x0+lam[:, None]*dx
        Fast = fobj(cpluss(x1, beta, stoich, full=True), stoich, analyticalc)
        gast = 0.5*np.sum(np.square(Fast), axis=1)
        return gast

    return g


def conditions(*args):
    # From sandbox/linsrch.py
    from functools import reduce
    return reduce(np.logical_and, args)


def linesearch1(x0, dx, slope, gfunc):
    if np.any(slope >= 0.0):
        raise ValueError('Negative slope values.')

    ALF = 1e-5
    TOLX = 1e-7

    g0 = gfunc(np.zeros(dx.shape[0]))

    lmin = TOLX / np.max(np.absolute(dx/x0), axis=1)
    assert lmin.shape == (x0.shape[0], )

    lam = np.ones_like(lmin, dtype=np.float)
    lam2 = np.zeros_like(lam)
    f2 = np.zeros_like(lam)
    first = True

    while True:
        g_lam = gfunc(lam)

        done_on_x = lam < lmin
        done_on_f = g_lam < (g0 + ALF*lam*slope)
        done = np.logical_or(done_on_x, done_on_f)
        if np.all(done):
            break

        _lam = lam[~done]
        _slope = slope[~done]
        _g_lam = g_lam[~done]
        _g0 = g0[~done]

        if first:
            tmplam = -0.5*_slope / (_g_lam - _g0 - _slope)
            first = False
        else:
            assert len(lam2[~done]) == len(f2[~done]) == len(_lam) == \
                   len(_slope) == len(_g0) == len(_g_lam)
            tmplam = minimum_poly3(_lam, lam2[~done], _g_lam, f2[~done], _g0,
                                   _slope)

        lam2 = lam[...]
        f2 = g_lam[...]
        import pudb
        pudb.set_trace()
        assert len(lam2[~done]) == len(f2[~done]) == len(_lam) == len(_slope)
        lam[~done] = np.where(tmplam > 0.1*_lam, tmplam, 0.1*_lam)

    assert np.all(np.logical_and(lam <= 1.0, lam >= 0.0))
    x_lam = x0 + lam[:, None]*dx
    return lam, g_lam, x_lam


def minimum_poly3(lam, lam2, f, f2, fold, slope):
    tmplam = np.empty_like(lam)
    rhs1 = f-fold-lam*slope
    rhs2 = f2-fold-lam2*slope
    a = (rhs1/lam**2-rhs2/lam2**2)/(lam-lam2)
    b = (-lam2*rhs1/lam**2+lam*rhs2/lam2**2)/(lam-lam2)

    w = (a == 0)
    tmplam[w] = -slope[w]/(2.*b[w])

    disc = b**2-3.*a*slope
    w = (disc < 0.0)
    tmplam[w] = 0.5*lam[w]

    w = np.logical_and(disc >= 0.0, b <= 0.0)
    tmplam[w] = (-b[w]+np.sqrt(disc[w]))/(3.*a[w])

    w = np.logical_and(disc >= 0.0, b > 0.0)
    tmplam[w] = -slope[w]/(b[w]+np.sqrt(disc[w]))

    w = (tmplam > 0.5*lam)
    tmplam[w] = 0.5*lam[w]
    return tmplam


def __mask(array):
    new_mask = collapse_mask_axis(array <= 0.0)
    np.ma.asarray(array)
    array.mask = new_mask


def collapse_mask_axis(array, axis=0):
    """Compact a multidimmensional mask into 1D mask.

    Given a mask array, it returns a 1D array based on one axis of it where
    :func:`numpy.any` is performed for every element of the axis in question.

    Parameters:
        array (:class:`numpy.ndarray`): Boolean mask array
        axis (int): the axis in which the operation is performed.

    Returns:
        :class:`numpy.ndarray`: the collapsed mask

    .. seealso:: :func:`mask_dim`
    """
    dims = tuple(i for i in range(array.ndim) if i != axis)
    return np.any(array, axis=dims)


def mask_dim(array, mask, axis):
    """Apply 1D-mask to a multidimmensional array along axis.

    Parameters:
        array (:class:`numpy.ndarray`): Boolean mask array
        mask (iterable): Boolean mask array
        axis (int): the axis in which the operation is performed.

    Returns:
        :class:`numpy.ndarray`: the collapsed mask
    """
    if array.shape[axis] != len(mask):
        raise ValueError

    shape1 = tuple(len(mask) if _ == axis else 1 for _ in range(array.ndim))
    shape2 = array.shape
    full_mask = np.broadcast_to(mask.reshape(shape1), shape2)
    return np.ma.array(array, mask=full_mask)
