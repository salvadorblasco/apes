# pylint: disable=invalid-name
# pylint: disable=no-member

"""
Module :mod:`libcal` contains the routines needed for fitting equilibrium
constants and enthalpies from calorimetric data. The main function to invoke
here is :func:`calfit` which handles most of the work. This is the only public
function for this module."""

import numpy as np

import consts
import excepts
import libaux
import libeq
import libio


def calfit(logB0, Bflags, H, Hflags, P, Q, titration, free_concentration, method=consts.METHOD_LM,
           **kwargs):
    r"""This function fits equilibrium constants for a calorimetric titration.

    Parameters:
        logB (:class:`numpy.ndarray`): 1D-array of floats representing the
            initial guess for the constants in the form of log\ :sub:`10`\ .
        Bflags (list): refinement flags for values in logB.
        H (:class:`numpy.ndarray`): 1D-array of floats representing the
            initial guess for the enthalpy values.
        Hflags (list): refinement flags for values in H.
        P (:class:`numpy.ndarray`): The stoichiometric coefficients in the
            usual format.
        Q (list of :class:`numpy.ndarray`): The stoichiometric coefficients
            in the usual format.
        titration (dict): Information about the titrations. Two options are
            accepted: (1) The following key must be provided (a) 'V'
            (:class:`numpy.ndarray`, volume of titre in mL), (b) 'V0' (float,
            the initial volume in mL), (c) 'T0' (:class:`numpy.ndarray`, total
            amount of initial components in mmol), (d) buret
            (:class:`numpy.ndarray`, concentration in the buret for each
            component in mmol/mL) and (e) 'Tflags' (list, optional) or (2)
        method (int): A number indicating the method to be used to mimize
            the objective functions. Allowed values are 0 for the algorithm
            by Levenberg-Marquardt (nonlinear least squares, this is the
            default one) and 1 for the Nelder-Mead algorithm (the simplex
            method).

    Returns:
        tuple: containing the follogin data (1) :class:`numpy.ndarray`: The
            refined constants, (2) the refined enthalpy values, (3) the free
            concentrations and (4) additional fitting data
    """

    if Bflags.count(0) == len(Bflags):
        if 'Tflags' in titration:
            raise NotImplementedError()

        if 'T' in titration:
            T = titration['T']
        else:
            T = libaux.build_T_titr2(titration['T0'], titration['buret'],
                                     titration['V0'], titration['V'])
        C = libeq.consol.consol(10**logB0, P, T, free_concentration)
        N = titration['V'][:, np.newaxis] * C
        refined_enthalpy = calfit_Honly(Q, N, H, Hflags)
        refined_beta = logB0
    else:
        raise NotImplementedError()

    return refined_beta, refined_enthalpy


def calfit_LM1(logB, Bflags, H, Hflags, P, Q, V, T, **kwargs):
    r"""This function fits equilibrium constants for a calorimetric titration
    by means of the Levenberg-Marquardt algoritm for the case when there is
    no refinement flag in the initial state.

    Parameters:
        logB (:class:`numpy.ndarray`): 1D-array of floats representing the
            initial guess for the constants in the form of log\ :sub:`10`\ .
        Bflags (:class:`numpy.ndarray`): refinement flags for values in
            logB0.
        H (:class:`numpy.ndarray`): 1D-array of floats representing the
            initial guess for the enthalpy values.
        Hflags (list): refinement flags for values in H.
        P (:class:`numpy.ndarray`): The stoichiometric coefficients in the
            usual format.
        Q (list of :class:`numpy.ndarray`): The stoichiometric coefficients
            in the usual format.
        T (:class:`numpy.ndarray`): The analytical concentrations for the
            initial components in mmol/mL

    Returns:
        tuple:
        - :class:`numpy.ndarray`: The refined constants in natural logarithmic
            units
        - :class:`numpy.ndarray`: The free concentrations
        - dict: Extra optional parameters

    Raises:
        ValueError: If invalid parameters are passed.
    """

    import functools

    # threshold = libaux.setkwpop(kwargs, 'threshold', 1e-10)
    # qt_handle = libaux.setkwpop(kwargs, 'qt_handle', None)
    # htmlout = libaux.setkwpop(kwargs, 'htmlout', None)
    # quiet_maxits = libaux.setkwpop(kwargs, 'quiet_maxits', False)
    # out_chisq = libaux.setkwpop(kwargs, 'out_chisq', [])
    # weights = libaux.setkwpop(kwargs, 'weights', np.ones(len(Q)))
    # c0 = libaux.setkwpop(kwargs, 'c0', None)
    # verbosity = libaux.setkwpop(kwargs, 'verbosity', 1)
    # max_iterations = libaux.setkwpop(kwargs, 'max_iterations', 20)
    threshold = kwargs.pop('threshold', 1e-10)
    qt_handle = kwargs.pop('qt_handle', None)
    htmlout = kwargs.pop('htmlout', None)
    quiet_maxits = kwargs.pop('quiet_maxits', False)
    out_chisq = kwargs.pop('out_chisq', [])
    weights = kwargs.pop('weights', np.ones(len(Q)))
    c0 = kwargs.pop('c0', None)
    verbosity = kwargs.pop('verbosity', 1)
    max_iterations = kwargs.pop('max_iterations', 20)

    iterations = 1
    var = np.flatnonzero(Bflags)
    libaux.assert_same_len(Q, weights)
    W = np.diag(weights)

    E, S = P.shape
    N = len(Q)
    if len(Bflags) != E:
        if qt_handle is not None:
            qt_handle('Bflags has a wrong number of data')
        raise ValueError('Bflags has a wrong number of data')

    if htmlout is not None:
        htmlout("Starting Levenberg-Marquardt iterations")

    # gather kwargs to pass forward to consol routine
    consol_kwargs = {}
    for k in kwargs.keys():
        if k[:7] == 'consol_':
            consol_kwargs[k[7:]] = kwargs[k]

    C = np.zeros((N, E+S))
    if c0 is not None:
        if c0.shape[1] == S:
            C[:, :S] = c0
        elif c0.shape[1] == E+S:
            C = c0
        else:
            raise ValueError
    else:
        C = libeq.initial_guess(np.exp(logB), P, T)

    fcalcC = functools.partial(libeq.consol, P=P, T=T, x0=C[:, :S],
                               **consol_kwargs)
    C = fcalcC(logB)

    # compute χ₂(dx)
    F = Q - calcQ(V, C[:, S:], H)
    chisq = np.sum(F**2)
    out_chisq.append(chisq)

    if 'damping' in kwargs:
        damping = kwargs.pop('damping')
    else:
        damping = 0.001

    x = libaux.unravel(logB, Bflags)
    J = cal_jac1(P, C, var)
    M = np.dot(np.dot(J.T, W), J)
    D = np.diag(np.diag(M))

    while True:
        try:
            dx = np.linalg.solve(M+damping*D, np.dot(np.dot(J.T, W), F.T))
        except np.linalg.linalg.LinAlgError:
            damping *= 10
            continue

        try:
            C = fcalcC(np.exp(libaux.ravel(logB, x+dx, Bflags)))
        except:
            diagn = {'iteration': iterations, 'x': x, 'dx': dx, 'J': J,
                     'chisq': chisq, 'residuals': F, 'damping': damping}
            raise excepts.FailedCalculateConcentrations(**diagn)
        F = Q - calcQ(V, C[:, S:], H)

        if qt_handle is not None:
            qt_handle("iter %d. chisq = %f" % (iterations, chisq))

        if np.sum(F**2) >= chisq:
            damping *= 10
        else:
            if htmlout:
                htmlout(libio.html_iteration(iterations, x/consts.LOGK,
                                             dx/consts.LOGK, chisq))

            iterations += 1
            damping /= 10
            chisq = np.sum(F**2)
            out_chisq.append(chisq)
            x += dx
            if 'one_iter' in kwargs:
                break
            J = cal_jac1(P, C, var)
            M = np.dot(np.dot(J.T, W), J)
            D = np.diag(np.diag(M))

        if np.all(np.abs(dx)/x < threshold):
            break

        if iterations > max_iterations:
            if htmlout and verbosity:
                htmlout("Maximum number of iterations has been reached." +
                        "Fitting will stop now.")
            if quiet_maxits:
                break
            else:
                raise excepts.TooManyIterations(libaux.ravel(logB,
                                                             x, Bflags))

    error_B = np.diag(np.linalg.inv(M))
    CV = libaux.covariance(J, weights, F)
    D = np.diag(CV)
    nD = len(D)
    CR = CV/np.sqrt(np.dot(D.reshape((nD, 1)), D.reshape((1, nD))))

    if htmlout and verbosity:
        htmlout(libio.html_finalconvg(x/consts.LOGK, error_B/consts.LOGK, CR))

    ret_extra = {'error': libaux.ravel(np.zeros(E), error_B, Bflags),
                 'damping': damping,
                 'covariance': CV,
                 'correlation': CR}

    return libaux.ravel(logB, x, Bflags), C, ret_extra


def calcQ(V, C, H):
    """C is the concentrations of the complexes only.  """
    aux1 = V[:, np.newaxis] * C
    aux2 = aux1[:-1, :] - aux1[1:, :]
    return np.sum(aux2*H, axis=1)


def calcN(C, V):
    """C is the concentrations of the complexes only.  """
    return V[:, np.newaxis] * C


def calfit_Honly(Q, N, H, Hflags):
    r"""Calculates the enthalpy values for the case when the equilibrium
    constants are not to be fit. This case simplifies significantly
    because a linear fitting is enough to find the values.

    .. math::

       \textbf{Q} - \mathbf{N}^0\mathbf{H}^0 = \mathbf{N}^1\mathbf{H}^1 +
       \mathbf{N}^2\mathbf{H}^2 + \ldots + \mathbf{N}^x\mathbf{H}^x

    .. math::

       q_i - \sum_{j=[0]} n_{ij} h_j = \sum_{j=[1]} n_{ij} h_j +
       \sum_k \left(\sum_{j=[k]} n_{ij} \frac{h_j}{h_{[k]}} \right) h_{[k]}
    """

    from numpy import linalg

    # I. define constants
    maxconstr = 6
    # TODO flags should be checked before, not at this point
    libaux.check_refine_flags(Hflags)

    # the keys are the constraint number and the values are the indices of the
    # components to which the constraint apply
    indices = {n: [i for i, j in enumerate(Hflags) if n == j]
               for n in range(max(Hflags))}

    Qret = Q - np.dot(np.take(N, indices[0], axis=1),
                      np.take(H, indices[0]))

    nsize = len(indices[1]) + sum(i in indices for i in range(2, maxconstr))

    # Compute the modified matrix.
    N1 = np.empty(len(Q), nsize)
    N1[:, len(indices[1])] = np.take(N, indices[1], axis=1)

    cols = iter(range(len(indices[1]-nsize, nsize)))
    for k in range(2, maxconstr):
        if k in indices:
            c = next(cols)
            N1[:, c] = np.take(N, indices[k], axis=1) * \
                np.take(H, indices[k]) / H[c]

    Hfit = linalg.lstsq(Qret, N1)[0]
    return libaux.ravel(H, Hfit, Hflags)


def cal_jac1() -> np.ndarray:
    pass


def cal_jac2(H, V, C, P):
    r"""This routine calculates the part of the jacobian that deals with the β

    .. math::

      \frac{\partial Q_i}{\partial\log\beta_j} =
      \sum_k^E H_k \left(
        V_i\frac{\partial c_{i,k+S}}{\partial\log\beta_j} -
        V_{i-1}\frac{\partial c_{i-1,k+S}}{\partial\log\beta_j}
      \right)

    Parameters:
        H (:class:`numpy.ndarray`): The enthalpy values. Must be N-1 sized
        V (:class:`numpy.ndarray`): The total volume of each point in the
            titration (initial plus titre). Must be N-sized.
        C (:class:`numpy.ndarray`): The :term:`free concentrations array`.
        P (:class:`numpy.ndarray`): The :term:`stoichiometry array`.
    """

    S = P.shape[1]
    dcdlB = libeq.extdd_dcdB(C, P) / C
    aux1 = V[:, np.newaxis] * dcdlB[:, S:]
    return np.sum(H[:, np.newaxis] * aux1[:, 1:] - aux1[:, :-1],
                  axis=1)


def cal_jac3():
    pass
