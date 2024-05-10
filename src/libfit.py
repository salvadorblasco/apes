"""General functions for nonlinear fitting."""

import math

import numpy as np

import consts
import excepts
import libmath
import report


def levenberg_marquardt(bridge, **kwargs):
    r"""Non linear fitting by means of the Levenberg-Marquardt method.

    Parameters:
        x0 (:class:`numpy.ndarray`): initial guess.
        weights (1D-array of floats): containing the values for weighting. It
            must be the same shape and type as *y*.
        func (callable): A function that accepts the values of *x0* and
            return both the residuals vector and the jacobian matrix.
        max_iterations (int, optional): maximum number of iterations allowed
        threshold (float, optional): criteria for convergence
        out_chisq (list, optional): If provided, the successive values for
            χ² in each iteration will be stored.
        verbosity (int, optional): An 0-2 number indicating the level of
            verbosity to be printed. 0 for mute, 1 for normal and 2 for
            pedantic output.
        report (callable, optional): A callable function that accepts the
            values of x0, iteration counter, free concentration values,
            etc., and is called every iteration in order to report on the
            progress of the fitting.
        one_iter (bool, optional): Performs one iterations and returns the
            result.
        quiet_maxits (bool, optional): Prevents this funcyion from throwing
            :class:`excepts.TooManyIterations` and quietly exits and returns
            the result when the maximum number of iterations is reached.

    Returns:
        tuple:
        - :class:`numpy.ndarray`: The refined constants in natural logarithmic
            units
        - :class:`numpy.ndarray`: The free concentrations
        - dict: Extra optional parameters

    Raises:
        ValueError: If invalid parameters are passed.
    """
    # def _report(*kws):
    #     if report is not None:
    #         report(*kws)

    report_buffer = kwargs.get('report', DummyReport())
    one_iter = kwargs.get('one_iter', False)
    threshold = kwargs.pop('threshold', 1e-4)
    max_iterations = kwargs.pop('max_iterations', 100)
    quiet_maxits = kwargs.get('quiet_maxits', False)
    damping = kwargs.pop('damping', 100.0)
    # fcapping = trivial_capping if capping is None else capping 

    n_points, n_vars = bridge.size()
    chisq_hist = []
    sigma_hist = []

    iterations = 1
    weights = bridge.weights()
    W = np.diag(weights)
    # assert W.shape == (n_points, n_points)

    # x = np.copy(x0)

    # # compute χ₂(dx)
    # J, resid = func(x)
    # assert resid.shape == weights.shape
    # chisq = np.sum(resid**2)
    # sigma = fit_sigma(resid, weights, n_points, n_vars)
    # assert isinstance(chisq, float)
    # assert J.ndim == 2 and J.shape == (n_points, n_vars)
    # M = J.T @ W @ J
    # # M = np.dot(np.dot(J.T, W), J)
    # D = np.diag(np.diag(M))

    # assert W.shape == (n_points, n_points)
    chisq = 1e99
    sigma = math.inf
    advance = 0

    # breakpoint()
    for iteration in range(max_iterations):
        # try:
        if iteration:
            # dx = np.linalg.solve(M+damping*D, np.dot(np.dot(J.T, W), resid))
            dx = np.linalg.solve(M+damping*D, J.T @ W @ resid)
        else:
            dx = np.zeros(n_vars)
        # except np.linalg.linalg.LinAlgError:
        #     damping *= 10
        #     continue

        # new_x = x + dx
        #new_x = fcapping(x, dx)
        #J, resid = func(new_x)
        bridge.step_values(dx)
        J, resid = bridge.build_matrices()

        new_chisq = np.sum(resid**2)
        test = (chisq-new_chisq)/chisq

        # print(iteration, dx)
        # print(f'\t {damping:10.4e}  {test:10.4e} {sigma:10.4e}')

        if new_chisq >= chisq:
            damping *= 10
            # print('\tnot decreasing')
        else:
            advance += 1
            # report_buffer.write(report.iteration(new_x, dx))
            print(f"{iteration=:4d}, {damping=:6.2e}, {test=:10.4e}, {dx=}")
            # print('\tdecreasing')
            bridge.accept_values()
            # _report(iterations, x/consts.LOGK, dx/consts.LOGK, chisq)
            # iterations += 1
            damping /= 5
            sigma = fit_sigma(resid, weights, n_points, n_vars)
            # x = new_x
            if one_iter:
                break
            M = J.T @ W @ J
            D = np.diag(np.diag(M))
            chisq = new_chisq
            # chisq_hist.append(chisq)
            # sigma_hist.append(sigma)

            bridge.report_step(iteration=advance, damping=damping, chisq=chisq, sigma=sigma)

            if (test < threshold) and iteration > 2:
                break

    else:
        # if np.all(np.abs(dx)/x < threshold):
        #     break
        ret = {'jacobian': J, 
               'residuals': resid,
               'damping': damping, 'convergence': chisq_hist,
               'iterations': iterations}
        raise excepts.TooManyIterations(msg=("Maximum number of iterations reached"),
                                        last_value=ret)

    ret_extra = {'jacobian': J, 'residuals': resid,
                 'damping': damping, 'convergence': chisq_hist,
                 'sigma': sigma_hist,
                 'iterations': iterations}
    return ret_extra


def simplex(x0, y, fnc, free_conc, weights, **kwargs):
    r"""Non linear fitting by means of the Nelder-Mead method (SIMPLEX).

    See `Numerical Recipes in C, §10.4
    <http://www.nrbook.com/a/bookcpdf/c10-4.pdf>`_
    and `http://www.scholarpedia.org/article/Nelder-Mead_algorithm`_

    Parameters:
        x0 (:class:`numpy.ndarray`): initial guess.
        y (:class:`numpy.ndarray`): the experimental magnitude to be fitted.
        weights (1D-array of floats): containing the values for weighting
        fnc (callable): A function that accepts the values of *x0* as well as
            the free concentrations and return the calculated values for *y*.
        free_conc (callable): A function that accepts *x0* and returns the
            values of the free concentration.
        weights (1D-array of floats): containing the values for weighting. It
            must be the same shape and type as *y*.
        max_iterations (int, optional, default: 20): maximum number of
            iterations allowed.
        verbosity (int, optional): An 0-2 number indicating the level of
            verbosity to be printed. 0 for mute, 1 for normal and 2 for
            pedantic output.
        report (callable, optional): A callable function that accepts the
            values of x0, iteration counter, free concentration values,
            etc., and is called every iteration in order to report on the
            progress of the fitting.
        term_x (float, optional, default=1e-5): Minimum simplex size threshold
            defining the convergence.
        term_f (float, optional, default=1e-7): Minimum objective function
            change to define convergence.

    Returns:
        tuple:
        - :class:`numpy.ndarray`: The refined constants in natural logarithmic
            units
        - :class:`numpy.ndarray`: The free concentrations
        - dict: Extra optional parameters

    Raises:
        ValueError: If invalid parameters are passed.
    """

    # i. Check input parameters
    # -------------------------
    report = kwargs.get('report', None)
    max_iterations = kwargs.get('max_iterations', 20)
    term_x = kwargs.get('term_x', 1e-5)
    term_f = kwargs.get('term_f', 1e-7)

    # ii. initial parameters
    # ----------------------
    Nr = len(x0)    # number of variables to refine
    h = Nr*[0.1]    # the initial steps in each dimmension
    alpha = 1.0     # constant for reflection operation
    beta = 0.5      # constant for contraction operation
    gamma = 2       # constant for expansion operation
    delta = 0.5     # constant for shrinking operation
    iteration = 0
    chisq_hist = []

    def _report(**kws):
        if report is not None:
            report(**kws)

    def _fobj(_x_):
        _concs = free_conc(_x_)
        y_calc = fnc(_x_, _concs)
        _f = np.sum((weights*(y - y_calc))**2)  # compute χ₂(dx)
        return _f, _concs

    # iii. build initial simplex, n+1 points
    # --------------------------------------
    x = [x0]                # simplex coordinates
    _f, _c = _fobj(x0)
    f = [_f]                # simplex values
    concs = [_c]

    for i in range(Nr):
        x_new = np.copy(x[0])
        x_new[i] += h[i]
        x.append(x_new)
        _f, _c = _fobj(x_new)
        f.append(_f)
        concs.append(_c)

    while True:
        iteration += 1

        # iv. Ordering: Determine the indices h, s, l of the worst, second
        #  worst and the best vertex, respectively, in the current working
        #  simplex S
        _h, _s, _l = _hsl(f)

        _report(iteration=iteration, x=x[_l]/consts.LOGK, chisq=f[_l])

        # v. Centroid: Calculate the centroid c of the best side—this is the
        #   one opposite the worst vertex xₕ
        c = _centroid(x[:_h] + x[(_h+1):])

        # vi. Transformation: Compute the new working simplex from the current
        #  one. First, try to replace only the worst vertex xₕ with a better
        #  point by using reflection, expansion or contraction with respect
        #  to the best side. All test points lie on the line defined by xₕ
        #  and c, and at most two of them are computed in one iteration. If
        #  this succeeds, the accepted point becomes the new vertex of the
        #  working simplex. If this fails, shrink the simplex towards the
        #  best vertex xₗ . In this case, n new vertices are computed.

        # vi(a). Reflect: Compute the reflection point xᵣ:=c+α(c−xₕ) and
        #  fᵣ:=f(xᵣ). If fₗ≤fᵣ<f_s , accept xᵣ and terminate the
        #  iteration.
        x_r = c + alpha*(c - x[_h])
        f_r, c_r = _fobj(x_r)

        if f[_l] <= f_r <= f[_s]:
            x_next, c_next, f_next = x_r, c_r, f_r

        elif f_r < f[_l]:
            # vi(b). Expand: If fᵣ<f_l, compute the expansion point
            #  xₑ:=c+γ(xᵣ−c) and fₑ:=f(xₑ) . If fₑ<fᵣ, accept xₑ and
            #  terminate the iteration. Otherwise (if fₑ≥fᵣ), accept xᵣ
            #  and terminate the iteration.
            #  This "greedy minimization" approach includes the better of
            #  the two points xᵣ, xₑ in the new simplex, and the simplex
            #  is expanded only if fₑ<fᵣ<fₗ. It is used in most
            #  implementations, and in theory (Lagarias et al., 1998).
            #  The original Nelder-Mead paper uses “greedy expansion”,
            #  where xₑ is accepted if fₑ<fₗ and fᵣ<fₗ, regardless
            #  of the relationship between fᵣ and fₑ. It may happen that
            #  fᵣ<fₑ, so xᵣ would be a  better new point than xₑ, and
            #  xₑ is still accepted for the new simplex. The working
            #  simplex is kept as large as possible, to
            #  avoid premature termination of iterations, which is sometimes
            #  useful for non-smooth functions.
            x_e = c + gamma*(x_r - c)
            f_e, c_e = _fobj(x_e)
            if f_e < f[_l]:
                x_next, c_next, f_next = x_e, c_e, f_e
            else:
                x_next, c_next, f_next = x_r, c_r, f_r

        elif f_r >= f[_s]:
            # vi(c). Contract: If fᵣ≥fₛ, compute the contraction point x_c by
            #   using the better of the two points xₕ and xᵣ.
            if f_r < f[_h]:
                # Outside: If fₛ≤fᵣ<fₕ, compute x_c:=c+β(xᵣ−c) and
                #  f_c:=f(x_c). If f_c≤fᵣ, accept x_c and terminate the
                #  iteration. Otherwise, perform a shrink transformation.
                x_c = c + beta*(x_r - c)
                f_c, c_c = _fobj(x_c)
                cond = f_c <= f_r
            else:
                # Inside: If f_r≥fₕ, compute x_c:=c+β(xₕ−c) and f_c:=f(x_c).
                #  If f_c<fₕ, accept x_c and terminate the iteration.
                #  Otherwise, perform a shrink transformation.
                x_c = c + beta*(x[_h] - c)
                f_c, c_c = _fobj(x_c)
                cond = f_c < f[_h]

            if cond:
                x_next, c_next, f_next = x_c, c_c, f_c
            else:
                # vi(d). Shrink: Compute n new vertices x_j:=xₗ+δ(x_j−xₗ)
                #  and f_j:=f(x_j) , for j=0,…,n , with j≠l.
                #  The shrink transformation was introduced to prevent the
                #  algorithm from failing in the following case, described
                #  by the quote from the original paper:
                #    A failed contraction is much rarer, but can occur when a
                #    valley is curved and one point of the simplex is much
                #    farther from the valley bottom than the others;
                #    contraction may then cause the reflected point to move
                #    away from the valley bottom instead of towards it.
                #    Further contractions are then useless. The action
                #    proposed contracts the simplex towards the lowest
                #    point, and will eventually bring all points into the
                #    valley.
                for j in range(Nr):
                    if j != _l:
                        x[j] = x[_l] + delta*(x[j]-x[_l])
                        f[j], concs[j] = _fobj(x[j])
                continue

        # vii. Termination test: A practical implementation of the Nelder-Mead
        #   method must include a test that ensures termination in a finite
        #   amount of time. The termination test is often composed of three
        #   different parts: term_x, term_f and fail.
        #   - term_x is the domain convergence or termination test. It becomes
        #   true when the working simplex S is sufficiently small in some sense
        #   (some or all vertices x_j are close enough).
        #   - term_f is the function-value convergence test. It becomes true
        #   when (some or all) function values f_j are close enough in some
        #   sense.
        #   - fail is the no-convergence test. It becomes true if the number
        #   of iterations or function evaluations exceeds some prescribed
        #   maximum allowed value.
        #   The algorithm terminates as soon as at least one of these tests
        #   becomes true.

        x[_h] = x_next      # replace worst point
        f[_h] = f_next   # _fobj(x_next)
        concs[_h] = c_next
        chisq_hist.append(f_next)

        if all(np.linalg.norm(x[i] - x[j]) < term_x
               for i in range(Nr) for j in range(Nr) if i != j):
            break

        if all(ff < term_f for ff in f):
            break

        if iteration > max_iterations:
            break

    _h, _s, _l = _hsl(f)
    x_final = x[_l]
    c_final = concs[_l]

    ret_extra = {'iterations': iteration, 'convergence': chisq_hist}

    return x_final, c_final, ret_extra


def final_params(jacobian, weights, resid):
    error_beta = libmath.error_params(jacobian, weights)
    covar = libmath.covariance(jacobian, weights, resid)
    correl = libmath.correlation_matrix(covar)
    return error_beta, covar, correl


def _centroid(x):
    """Given a list of vectors, return the centroid.

    Parameters:
        x (:class:`numpy.ndarray`): A list of 1D :class:`numpy.ndarray` of
            the same length representing the vector for calculating the
            centroid.

    Returns:
        :class:`numpy.ndarray`: A vector representing the centroid.

    >>> _centroid([np.array([1,2,3,4]),
    ...            np.array([5,6,7,8]),
    ...            np.array([9,10,11,12])])
    array([ 5.,  6.,  7.,  8.])
    """
    return np.sum(np.vstack(x), axis=0) / len(x)


def _hsl(lst):
    """Given a list of floats, return the indices of the worst (biggest
    value), second worst and best (lowest value)

    >>> _hsl([5,3,7,1,0])
    3, 0, 4
    """
    if not len(lst) > 2:
        breakpoint()
    assert len(lst) > 2
    f_idx = [lst.index(ff) for ff in sorted(lst, reverse=True)]
    return f_idx[0], f_idx[1], f_idx[-1]


def fit_sigma(residuals, weights, npoints, nparams):
    return np.sum(weights*residuals**2)/(npoints-nparams)


class DummyReport:
    def write(self, data):
        pass

fitting_functions = {
    consts.METHOD_LM: levenberg_marquardt,
    consts.METHOD_NM: simplex
}
