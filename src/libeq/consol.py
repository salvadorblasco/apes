r"""This module contains the routines needed for equilibria solving.

Basically it is all reduced to the solving of the mass balance equation:

.. math:: c_i + \sum_{j=1}^E {p_{ji} c_{j+S}} - T_i = 0

where

.. math:: c_{i+S} = \beta_i \prod_{j=1}^S c_j^{p_{ij}}

The function that does the basic task is the function :func:`consol`
that takes the equilibrium constants, the stoichiometric coefficients and
the starting conditions and returns the free concentrations.

- *N* is the number of experimental points
- *S* is the number of free components
- *E* is the number of equilibria
"""

import itertools

import numpy as np

import libaux
import excepts

from . import LIBEQ_ENGINE

__version__ = '0.5'


def consol(beta, stoichiometry, analytc, initial_values, **kwargs):
    """Solve the equilibrium and calculate the free concentrations.

    The computation of the free concentation is done by using the
    `Newton-Raphson method <https://en.wikipedia.org/wiki/Newton%27s_method>`_
    task that is passed to :func:`NewtonRaphson`

    Parameters:
        beta (:class:`numpy.ndarray`): The equilibrium constants. It must be an
            array of *E* floats or an array of (*N*, *E*) floats.
        stoichiometry (:class:`numpy.ndarray`): stoichiometric coefficients. It
            must be an array of (*E*, *S*) ints.
        analytc (:class:`numpy.ndarray`): total concentrations of every point
            in mmol/mL. It must be an array of (*N*, *S*) floats.
        initial_values (:class:`numpy.ndarray`): initial guess for the first
            point. If they are unknown, :func:`initial_guess` may be called
            first.  It is highly advisable to provide a good initial estimate
            in order to have a
            fast convergence. It must be an (*N*, *S*)-sized or an *S*-sized
            array of floats. If it is an 1D array then the option **stepwise**
            will be set to **True** regardless of the value
            set for this flag by the user. If it is a 2D array and the
            flag **stepwise** is set to **True**, only the first row of **x0**
            will be used. It this parameter is not provided, then the function
            :func:`initial_guess` will be called to provide an initial
            estimate.
        stepwise (:obj:`bool`, optional[False]): If this flag is set, the
            concentrations in the titration curve will be calculated one by
            one taking as initial estimate for one curve the value found for
            the previous one. This results in more reliable values but it is
            significantly slower that refining by the bulk.
        extrapolation (:obj:`int`, optional[0]): Defines the order of the
            polynomial extrapolation to use. If extrapolation=0 no
            extrapolation will be used and the previous value.
            Only used together with stepwise=True,
            otherwise this option is ignored.
        errors (:obj:`str`, optional['silent']): Accepted values are : (1)
            'silent': this will silence TooManyIterations and UnstableIteration
            and will return as a result the value of the last iteration. (2)
            'print': this does the same as 'silent' but prints a warning in the
            standard output. (3) 'raise': raises the exception
        kwargs (dict): Additional arguments that will be passed to
            :func:`NewtonRaphson`
    Returns:
        :class:`numpy.ndarray`: 2D-array containing the free concentrations
            in mol/L
    Notes:
        Here *E* is the number of equilibria, *S* is the number of species and
        *N* is the number of experimental points.

    .. [#] *Talanta* 1996, **43**, 1739-1756
    """
    libaux.assert_BPT_consistency(beta, stoichiometry, analytc)

    caller = {'python': _consol_python, 'fortran': _consol_fortran}
    if LIBEQ_ENGINE not in caller:
        msg = ("Invalid parameter 'engine'."
               "Accepted values: 'python' or 'fortran'.")
        raise ValueError(msg)
    fcall = caller[LIBEQ_ENGINE]

    calc_concentration = fcall(beta, stoichiometry, analytc, initial_values,
                               **kwargs)
    return calc_concentration


def _consol_fortran(beta, stoichiometry, analytc, initial_values, **kwargs):
    """Solve the for the free concentration using the fortran engine."""
    fcall = __select_fortran_call(beta)

    n_points = analytc.shape[0]
    n_equil, n_species = stoichiometry.shape

    f_beta = __forder(beta)
    f_stoich = __forder(stoichiometry)
    f_analc = __forder(analytc)
    f_inivals = np.empty((n_points, n_species+n_equil), order='F')
    f_inivals[:, :n_species] = initial_values[:, :n_species]

    ctrl = 1
    # fcall(f_inivals, f_analc, f_beta, f_stoich, ctrl, n_points, n_equil,
    #       n_species)
    fcall(f_inivals, f_analc, f_beta, f_stoich, ctrl)
    return f_inivals


def _consol_python(beta, stoichiometry, analytc, initial_values,
                   stepwise=False, extrapolation=0, errors='silent',
                   **kwargs):
    """Solve the equilibrium and calculate the free concentrations.

    The computation of the free concentation is done by using the
    `Newton-Raphson method <https://en.wikipedia.org/wiki/Newton%27s_method>`_
    task that is passed to :func:`NewtonRaphson`

    Parameters:
        beta (:class:`numpy.ndarray`): The equilibrium constants. It must be an
            array of *E* floats or an array of (*N*, *E*) floats.
        stoichiometry (:class:`numpy.ndarray`): stoichiometric coefficients. It
            must be an array of (*E*, *S*) ints.
        analytc (:class:`numpy.ndarray`): total concentrations of every point
            in mmol/mL. It must be an array of (*N*, *S*) floats.
        initial_values (:class:`numpy.ndarray`): initial guess for the first
            point. If not provided, it will be estimated. It is highly
            advisable to provide a good initial estimate in order to have a
            fast convergence. It must be an (*N*, *S*)-sized or an *S*-sized
            array of floats. If it is an 1D array then the option **stepwise**
            will be set to **True** regardless of the value
            set for this flag by the user. If it is a 2D array and the
            flag **stepwise** is set to **True**, only the first row of **x0**
            will be used. It this parameter is not provided, then the function
            :func:`initial_guess` will be called to provide an initial
            estimate.
        stepwise (:obj:`bool`, optional[False]): If this flag is set, the
            concentrations in the titration curve will be calculated one by
            one taking as initial estimate for one curve the value found for
            the previous one. This results in more reliable values but it is
            significantly slower that refining by the bulk.
        extrapolation (:obj:`int`, optional[0]): Defines the order of the
            polynomial extrapolation to use. If extrapolation=0 no
            extrapolation will be used and the previous value.
            Only used together with stepwise=True,
            otherwise this option is ignored.
        errors (:obj:`str`, optional['silent']): Accepted values are : (1)
            'silent': this will silence TooManyIterations and UnstableIteration
            and will return as a result the value of the last iteration. (2)
            'print': this does the same as 'silent' but prints a warning in the
            standard output. (3) 'raise': raises the exception
        Cbuffer (:class:`numpy.ndarray`, optional): free concentrations of the
            previous iteration. Used in conjunction with interpolation=
            "linear". Also used as initial guess if *x0* is not provided.
        kwargs (dict): Additional arguments that will be passed to
            :func:`NewtonRaphson`
    Returns:
        :class:`numpy.ndarray`: 2D-array containing the free concentrations
            in mol/L
    Notes:
        Here *E* is the number of equilibria, *S* is the number of species and
        *N* is the number of experimental points.

    .. [#] *Talanta* 1996, **43**, 1739-1756
    """
    from libeq.nr import NewtonRaphson
    from libeq.cpluss import cpluss
    # libaux.assert_array_dim(1, beta)
    # libaux.assert_array_dim(2, T, stoichiometry)
    new_c = None
    c0 = None

    n_data = analytc.shape[0]
    n_equilibria, n_species = stoichiometry.shape

    result = np.zeros((n_data, n_equilibria+n_species), dtype=float)

    if stepwise:
        if beta.ndim == 2:
            iter_beta = np.nditer(beta, flags=['external_loop'])
        else:
            iter_beta = itertools.repeat(beta)

        for row, (t, b) in enumerate(zip(analytc, iter_beta)):
            if initial_values.ndim == 1 and row == 0:
                c0 = initial_values
            elif initial_values.ndim == 1 and row > 0:
                c0 = new_c
                # TODO Consider extrapolation
                # TODO lastn = slice(row-n, row)
                # TODO c0=libmath.extrapoly(np.atleast_2d(t).T, analytc[lastn],
                # TODO                        result[lastn, :n_species])
            elif initial_values.ndim == 2:
                c0 = initial_values[row, :]
            else:
                ValueError("Malformed 'initial_values'")

            try:
                # pr = __start_profile():
                new_c = NewtonRaphson(c0, b, stoichiometry, t, **kwargs)
                # __stop_profile(pr):
            except excepts.TooManyIterations as err:
                # TODO perhaps repeat point with different settings?
                if errors == 'raise':
                    raise err
                if errors == 'print':
                    print("Too many iterations in point %d" % row)
                new_c = err.last_value

            result[row, 0:n_species] = new_c
    else:
        try:
            # pr = __start_profile():
            new_c = NewtonRaphson(initial_values, beta, stoichiometry, analytc,
                                  **kwargs)
        except excepts.TooManyIterations as err:
            # TODO perhaps repeat point with different settings?
            new_c = err.last_value
            raise
        finally:
            # __stop_profile(pr):
            result[:, 0:n_species] = new_c

    result[:, n_species:n_species+n_equilibria] = \
        cpluss(result[:, :n_species], beta, stoichiometry)
    # TODO Check for mistakes in the values
    return result


def simulation(beta, stoichiometry, analytc, initial_values, reference=-1,
               **kwargs):
    r"""Simulate titration.

    This routine should work exactly as :func:`consol()` does, except that
    some of the concentrations are not unknowns but data. This function
    basically recalculates the input parameters to account on that factor and
    then call :func:`consol` to calculate the free concentrations and results
    the result.

    The internals of this function is very simple and just requires a
    rewriting of some equations.

    .. math:: \beta_i^{*} = \beta_i c_k^{p_{ik}}

    and then the species *k* is eliminated from *P* and *T*.

    Parameters:
        beta (:class:`numpy.ndarray`): The :term:`equilibrium constants array`
        stoichiometry (:class:`numpy.ndarray`): The :term:`stoichiometry array`
        analytc (:class:`numpy.ndarray`): The
            :term:`total concentrations array`
        reference (int): index of the independent variable. By default, -1 (the
                last one)
        kwargs: kwargs are not directly used but instead they are passed to
            inner routines. For more info see :func:`consol()` and
            :func:`NewtonRaphson()`.
    Returns:
        :class:`numpy.ndarray`: An array containing the calculated
            concentrations. The number of points returned is the same length
            as **analytc**.
    """
    raise DeprecationWarning
    libaux.assert_BPT_consistency(beta, stoichiometry, analytc)

    new_x = analytc[:, reference]      # new X
    beta_prime = beta[np.newaxis, ...] * \
        new_x[..., np.newaxis]**stoichiometry[:, reference]
    stoich_new = np.delete(stoichiometry, reference, axis=1)
    analc_new = np.delete(analytc, reference, axis=1)
    initial_new = np.delete(initial_values, reference, axis=1)
    concentration = consol(beta_prime, stoich_new, analc_new, initial_new,
                           **kwargs)
    return new_x, concentration


def initial_guess(beta, stoichiometry, analyticalc, **kwargs):
    """Provide an initial guess for the free concentration array.

    This function tries to evaluate an initial guess when no one is
    provided.

    Parameters:
        beta (:class:`numpy.ndarray`): The equilibrium constants array.
        stoichiometry (:class:`numpy.ndarray`): The stoichiometric coefficient
            array
        analyticalc (:class:`numpy.ndarray`): Analytical concentrations array.
        kwargs: Extra arguments are passed to :func:`consol`

    Returns:
        :class:`numpy.ndarray`: An estimate of the free concentrations

    .. note:: If :func:`NewtonRaphson` raises :class:`TooManyIterations` those
        exceptions are not re-raised.
    """
    libaux.assert_BPT_consistency(beta, stoichiometry, analyticalc)
    caller = {'python': _initial_guess_python,
              'fortran': _initial_guess_fortran}
    if LIBEQ_ENGINE not in caller:
        msg = ("Invalid parameter 'engine'."
               "Accepted values: 'python' or 'fortran'.")
        raise ValueError(msg)
    fcall = caller[LIBEQ_ENGINE]
    return fcall(beta, stoichiometry, analyticalc, **kwargs)


def _initial_guess_fortran(beta, stoichiometry, analyticalc, **kwargs):
    fcall = __select_fortran_call(beta)
    n_points = analyticalc.shape[0]
    n_equil, n_species = stoichiometry.shape

    x0 = np.copy(analyticalc)
    x0[x0 < 1e-6] = 1e-6
    # REMOVE. fortran call includes preconditioning
    # from libeq.damping import damping
    # damping(x0, beta, stoichiometry, analyticalc)

    f_beta = __forder(beta)
    f_stoich = __forder(stoichiometry)
    f_analc = __forder(analyticalc)
    f_inivals = np.empty((n_points, n_species+n_equil), order='F')
    f_inivals[:,:n_species] = x0[:,:]
    # bullet = analyticalc[0, :]
    # bullet[bullet <= 0.0] = 1e-6
    # f_inivals[0,:n_species] = bullet

    ctrl = 1
    try:
        fcall(f_inivals, f_analc, f_beta, f_stoich, ctrl)
        # fcall(f_inivals, f_analc, f_beta, f_stoich, ctrl, n_points, n_equil, n_species)
    except:
        print("Error catched")
    # fcall(f_inivals, f_analc, f_beta, f_stoich, ctrl)
    return f_inivals


def _initial_guess_python(beta, stoichiometry, analyticalc, **kwargs):
    """Pure Python version of :func:`initial_guess`.

    Parameters:
        beta (:class:`numpy.ndarray`): The equilibrium constants array.
        stoichiometry (:class:`numpy.ndarray`): The stoichiometric coefficient
            array
        analyticalc (:class:`numpy.ndarray`): Analytical concentrations array.
        kwargs: Extra arguments are passed to :func:`consol`

    Returns:
        :class:`numpy.ndarray`: An estimate of the free concentrations
    """
    def _do_consol(_beta, _stoich, _analc, **kwargs2):
        kwconsol = kwargs.copy()
        kwconsol.update(kwargs2)
        try:
            c1 = consol(_beta, _stoich, _analc, **kwconsol)
        except excepts.TooManyIterations as err:
            c1 = err.last_value
        return c1

    def _do_intp(_conc):
        if np.any(np.isnan(_conc)):
            mc = np.ma.masked_invalid(_conc)
            interpolate_masked(mc)
            assert not np.any(np.isnan(mc))
            _conc[...] = mc

    # from libeq.damping import damping
    from libeq.pcf import pcf

    n_species = stoichiometry.shape[1]
    n_points = analyticalc.shape[0]

    # Shorten data
    # if n_points > 50:
    #     every = 20
    # elif 20 < n_points < 50:
    #     every = 10
    # else:
    #     every = 1
    # short_analc = analyticalc[::every, ...]
    # if beta.ndim == 1:
    #     short_beta = beta
    # else:
    #     short_beta = beta[::every, ...]
    # x0 = np.copy(short_analc)
    # x0[x0 < 1e-6] = 1e-6
    x0 = np.full_like(analyticalc, 1e-6)
    c1 = pcf(x0, beta, stoichiometry, analyticalc, 1e-3)

    # c1 = _do_consol(short_beta, stoichiometry, short_analc,
    #                 initial_values=x0, scaling=False, threshold=1e-5,
    #                 max_iterations=200, zero_offdiag=True, mask_flag=False,
    #                 forcer=True, step_limiter=False)

    # _do_intp(c1)
    # x0 = np.empty_like(analyticalc)
    # x0[::every, :] = c1[:, :n_species]
    # resample(x0, every)

    # damping(x0, beta, stoichiometry, analyticalc)

    # c1 = _do_consol(beta, stoichiometry, analyticalc, initial_values=x0,
    #                 scaling=False, threshold=1e-5, max_iterations=200,
    #                 zero_offdiag=True, mask_flag=False, forcer=True)
    # _do_intp(c1)
    assert not np.any(np.isnan(c1))

    c1 = _do_consol(beta, stoichiometry, analyticalc,
                    initial_values=c1[:, :n_species], max_iterations=100,
                    threshold=1e-15, scaling=False, mask_flag=False,
                    damping=False, forcer=False)
    # _do_intp(c1)
    # assert not np.any(np.isnan(c1))

    # c1 = _do_consol(beta, stoichiometry, analyticalc,
    #                 initial_values=c1[:, :n_species], max_iterations=100,
    #                 threshold=1e-15, scaling=False, mask_flag=False,
    #                 damping=False, forcer=False)

    return c1[:, :n_species]


def freeze_concentration(beta, stoichiometry, analytc, reference=-1):
    r"""Convert one component to independent variable.

    When solving the equilibrium, sometimes those concentrations are plotted
    as function of the concentration of one of them, typically the pH. That
    component, therefore, can be converted into an independent variable
    and removed from the unknowns.

    .. math::

    c_{i+S} = \beta_ic_{\omega}^{p_{i\omega}}\prod_{j\ne\omega}^Sc_j^{p_{ij}}
            = \beta_i'\prod_{j\ne\omega}^Sc_j^{p_{ij}}

    Parameters:
        beta (:class:`numpy.ndarray`): The equilibrium constants array.
        stoichiometry (:class:`numpy.ndarray`): The stoichiometric coefficient
            array
        analyticalc (:class:`numpy.ndarray`): Analytical concentrations array.
        reference (int):

    Returns:
        new_x (:class:`numpy.ndarray`): The reference component concentrations
            which is equal to the *analyltc* reference column.
        beta_prime (:class:`numpy.ndarray`): The new beta array.
        stoich_new (:class:`numpy.ndarray`): The new stoichiometry array with
            is equal to the original one with the reference component removed.
        analc_new (:class:`numpy.ndarray`): The new analytical concentrations
            array with is equal to the original one with the reference
            component removed.
    """
    new_x = analytc[:, reference]
    beta_prime = beta[np.newaxis, ...] * \
        new_x[..., np.newaxis]**stoichiometry[:, reference]
    stoich_new = np.delete(stoichiometry, reference, axis=1)
    analc_new = np.delete(analytc, reference, axis=1)
    return new_x, beta_prime, stoich_new, analc_new


def beta_uncertainty(concentration, beta_error, stoichiometry,
                     confidence=0.95):
    """Calculate errors in the free concentrations.

    Returns the error in the calculated free concentrations at a cetain
    confidence levels given an uncertanty for the equilibrium constants.

    Parameters:
        concentration (:class:`numpy.ndarray`): An (*N*, *E* + *S*)-sized
            array of floats representing the values of the free concentrations
            for all the components.
        beta_error (:class:`numpy.ndarray`): errors for equilibrium constants
        stoichiometry (:class:`numpy.ndarray`): stoichiometric coefficients
        confidence (float): level of confidence, 0.95 by default. It must
            be 0 < confidence < 1
    Returns:
        :class:`numpy.ndarray`: An (*N*, *E* + *S*)-sized array of floats
            representing the errors of concentrations.
    Raises:
        ValueError: If an invalid parameter is passed
    """
    from scipy.stats import norm
    if not 0 < confidence < 1:
        raise ValueError("Confidence must be 0 < confidence < 1")

    libaux.assert_array_dim(2, stoichiometry, concentration)
    libaux.assert_array_dim(1, beta_error)
    n_points = concentration.shape[0]
    n_equilibria, n_species = stoichiometry.shape
    libaux.assert_shape((n_points, n_equilibria+n_species), concentration)

    q = norm.interval(confidence)[1]
    return q*concentration*np.sqrt(np.sum(
        np.square(beta_error[None, None, :] * extdd_dcdb(concentration,
                                                         stoichiometry)),
        axis=-1))


def dcdb(concentration, stoichiometry):
    r"""Calculate ∂logc/∂logβ.

    This function calculates the derivative of the logarithms of
    concentrations with respect to the logarithm of betas. It solves the
    following equation for :math:`\frac{\partial\log c_k}{\partial\log\beta_b}`

    .. math:: \sum_{k=1}^S \left( \delta_{ki}c_k + \sum_{j=1}^E {
       p_{ji} p_{jk} c_{j+S} }
       \right) \frac{\partial\log c_k}{\partial \log\beta_b}
       = -p_{bi}c_{b+S}

    Parameters:
        concentration (:class:`numpy.ndarray`): The free concentrations for all
            the components (*E* + *S* components). It can be 1D or 2D. If it is
            2D it is assumed that titration points go along axis=0 and species
            go along axis=1.
        stoichiometry (:class:`numpy.ndarray`): The stoichiometric coefficient
            array.

    Returns:
        :class:`numpy.ndarray`: An (*N*, *S*, *E*) array whose values [i,j,k]
            are :math:`\frac{\partial\log c_{i,j}}{\partial\log\beta_{k}}`
            If titration points are given in **concentration** then they are
            along axis=0, and row and col are along axis=1 and 2.

    Raises:
        ValueError: If an invalid parameter is passed.

    .. note:: The result is given in natural logarithmic units.
    """
    n_species = stoichiometry.shape[1]
    n_points = concentration.shape[0]
    if concentration.ndim == 1:
        d = np.diag(concentration[:n_species])
    else:
        d = np.zeros((n_points, n_species, n_species))
        _, j, k = np.indices(d.shape)
        d[j == k] = (concentration[..., :n_species]).flat
    A = d + np.einsum('ij,ik,...i->...jk', stoichiometry, stoichiometry,
                      concentration[..., n_species:])
    B = stoichiometry[np.newaxis, ...] * \
        concentration[..., n_species:, np.newaxis]
    return np.squeeze(np.linalg.solve(A, -B.swapaxes(-2, -1)))


def extdd_dcdb(concentration, stoichiometry):
    r"""Calculate ∂logc/∂logβ.

    This function is just like :func:`dcdB` but it also returns the values
    for all the components, not just the :term:`independent components`.

    .. math::

        \frac{\partial\log c_{i+S}}{\partial\log\beta_k} =
        \delta_{ik} + \sum_{j=1}^S p_{ij}
        \frac{\partial\log c_j}{\partial\log\beta_k}

    where :math:`{\partial\log c_j}/{\partial\log\beta_k}` are calculated
    calling :func:`dcdb`.

    Parameters:
        concentration (:class:`numpy.ndarray`): The
            :term:`free concentrations array` for
            all the components (*E* + *S* components). It can be 1D or 2D. If
            it is 2D it is assumed that titration points go along axis=0 and
            species go along axis=1.
        stoichiometry (:class:`numpy.ndarray`): The
            :term:`stoichiometry array`.
    Returns:
        :class:`numpy.ndarray`: An (*N*, *S* + *E*, *E*) array whose values
            [i, j, k] are
            :math:`\frac{\partial\log c_{i,j}}{\partial\log\beta_{k}}`

    .. seealso:: :func:`dcdb`
    .. note:: The result is given in natural logarithmic units.
    """
    n_equilibria = stoichiometry.shape[0]
    x1 = dcdb(concentration, stoichiometry)
    x2 = np.eye(n_equilibria) + np.einsum("ij,...jk->...ik", stoichiometry, x1)
    return np.concatenate((x1, x2), axis=-2)   # x -> dlogc[r]/dlogB[c]


def resample(array, every):
    """Interpolate the intermediate data of an array.

    Given an array, with accurate data in intermediate positions every *every*
    data and in-place replace the middle places with linearly interpolated
    values.

    Parameters:
        array (:class:`numpy.ndarray`): The data array where the
        every (int): The span of true data
    """
    n_points = len(array)
    anchor_indices = np.arange(0, n_points, every)
    wanted_indices = [i for i in range(n_points) if i % every]
    anchor = array[::every, ...]

    import scipy.interpolate
    f = scipy.interpolate.interp1d(anchor_indices, anchor, axis=0,
                                   fill_value=anchor[-1, :],
                                   bounds_error=False)
    new = f(wanted_indices)
    array[wanted_indices, ...] = new


def interpolate_masked(data):
    """Interpolate masked values using the nearest values.
    """
    n: int = len(data)
    mask = np.any(data.mask, axis=1)
    assert len(mask) == n

    indices = np.arange(len(data))
    assert len(indices) == n

    bad_idx = indices[mask]
    valid_idx = indices[~mask]
    valid_data = data[~mask, :]

    import scipy.interpolate
    itpf = scipy.interpolate.interp1d(valid_idx, valid_data, axis=0,
                                      assume_sorted=True,
                                      fill_value="extrapolate")
    newvals = itpf(bad_idx)
    data[bad_idx, :] = newvals
    assert not np.any(np.isnan(data.data))


def compute_concentration_combined(beta, stoichiometry, analytc, initial_values):
    analytic_packed = np.concatenate(analytc)
    lims = list(itertools.accumulate([0] + [s.shape[0] for s in analytc]))
    inivals = np.concatenate(initial_values)
    concentration_packed = consol(beta, stoichiometry, analytc, inivals)
    concentration = np.vsplit(concentration_packed, lims[1:-1])
    return concentration


# def __start_profile():
#     import cProfile
#     pr = cProfile.Profile()
#     pr.enable()
#     print(20*'=')
#     return pr
# 
# 
# def __stop_profile(pr):
#     import pstats
#     import io
#     pr.disable()
#     s = io.StringIO()
#     sortby = 'cumulative'
#     my_stats = pstats.Stats(pr, stream=s).sort_stats(sortby)
#     my_stats.print_stats()
#     print(s.getvalue())
#     print(20*'=')


def __forder(array):
    "Allocate arrays in fortran ordering."
    return np.require(array, requirements='F')


def __select_fortran_call(beta):
    import libeq.mfeq
    if beta.ndim == 1:
        fcall = libeq.mfeq.newtonr1
    elif beta.ndim == 2:
        fcall = libeq.mfeq.newtonr2
    else:
        raise ValueError('Malformed beta array')
    return fcall
