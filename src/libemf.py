# -*- coding: utf-8 -*-

"""Routines for constant fitting from potentiometric data.

Module :mod:`libemf` contains the routines needed for fitting equilibrium
constants from potentiometric data. The main function to invoke here is
:func:`emffit` which handles most of the work. This is the only public
function for this module.

.. rubric:: Functions included in this module
"""

import functools
import itertools

import numpy as np
from numpy.typing import NDArray

import consts
import libaux
import libeq.consol
import libeq.jacobian
import libfit
import libmath
import excepts
# import report

__version__ = '0.5'


def hselect(array, hindices):
    """Select columns that correspond to the electroactive species.

    Given the concentrations array, selects the columns that correspond
    to the electroactive species.

    Parameters:
        array (:class:`numpy.ndarray`): The :term:`free concentrations array`
        hindices (list[int]): Indices of the electroactive specie(s).

    Returns:
        The part of C which is electroactive
    """
    return array[...,hindices]

def nernst(electroactive_conc: NDArray[float], emf0: NDArray[float],
           slope: NDArray[float] | float = 1.0, joint: NDArray[float] | float=0.0,
           temperature: float=298.15) -> NDArray[float]:
    r"""Calculate the calculated potential.

    Apply Nernst's equation to calculate potential according to 
    .. math ::

    E = E^0 + f\frac{nF}{RT}\ln[C] + J

    Parameters:
        electroactive_conc (:class:`numpy.ndarray`): a 1D array of floats
            representing the free concentrations of the electroactive species.
        emf0 (:class:`numpy.ndarray`): The :term:`standard potential`
        slope (:class:`numpy.ndarray`): The slope for :term:`Nernst's equation`
        joint (:class:`numpy.ndarray`): The liquid joint contribution for :term:`Nernst's equation`
        temperature (float): the absolute temperature
    Returns:
        :class:`numpy.ndarray`: an array of floats containing the calculated values
    """
    nernstian_slope = slope*consts.RoverF*temperature
    return emf0 + nernstian_slope*np.log(electroactive_conc) + joint


def emf_jac_beta(dlogc_dlogbeta, slope=1.0, temperature: float=298.15):
    r"""Calculate the jacobian part related to equilibrium constants.

    The calculation is done according to equation
    .. math ::

    \frac{\partial E_n}{\partial\log\beta_b}=\frac{fRT}{nF\log e}\frac{\partial\log c_{nh}}{\partial\log\beta_b} 

    Parameters:
        dlogc_dlogbeta (:class:`numpy.ndarray`): the derivative values. They can be obtained
            from :func:`libeq.jacobian.dlogcdlogbeta`.
        slope (:class:`numpy.ndarray`): The slope for :term:`Nernst's equation`
        temperature (float): the absolute temperature
    Returns:
        :class:`numpy.ndarray`: an array of floats containing the calculated values
    """
    nernstian_slope = slope*consts.RoverF*temperature
    return nernstian_slope*dlogc_dlogbeta/consts.LOGK


def emf_jac_init(dlogc_dt, slope=1.0, temperature=298.15):
    r"""Calculate the jacobian part related to the initial amount.

    The calculation is done according to equation
    .. math ::

    \frac{\partial E_n}{\partial t_i} = f\frac{RT}{nF}\frac{\partial\log c_{nh}}{\partial t_i}

    Parameters:
        dlogc_dt (:class:`numpy.ndarray`): the derivative values. They can be obtained
            from :func:`libeq.jacobian.dlogcdt`.
        slope (:class:`numpy.ndarray`): The slope for :term:`Nernst's equation`
        temperature (float): the absolute temperature
    Returns:
        :class:`numpy.ndarray`: an array of floats containing the calculated values
    """
    nernstian_slope = slope*consts.RoverF*temperature
    return nernstian_slope*dlogc_dt


def emf_jac_buret(dlogc_db, slope=1.0, temperature=298.15):
    r"""Calculate the jacobian part related to the buret concentration.

    The calculation is done according to equation
    .. math ::

    \frac{\partial E_n}{\partial b_i} = f\frac{RT}{nF}\frac{\partial\log c_{nh}}{\partial b_i}

    Parameters:
        dlogc_db (:class:`numpy.ndarray`): the derivative values. They can be obtained
            from :func:`libeq.jacobian.dlogcdb`.
        slope (:class:`numpy.ndarray`): The slope for :term:`Nernst's equation`
        temperature (float): the absolute temperature
    Returns:
        :class:`numpy.ndarray`: an array of floats containing the calculated values
    """
    nernstian_slope = slope*consts.RoverF*temperature
    return nernstian_slope*dlogc_db


def emf_jac_e0(size: int) -> NDArray[float]:
    r"""Calculate the jacobian part related to the standard potential.

    It returns ones based on the size according to equation
    .. math ::

    \frac{\partial E_n}{\partial E^0} = 1

    Parameters:
        size (int): the number of ones to return
    Returns:
        :class:`numpy.ndarray`: an array of floats containing the calculated values
    """
    return np.ones(size)


# Everything below this line can probably be deleted 


def emffit(beta, beta_flags, stoichiometry, titration_data, electrodes,
           method=consts.METHOD_LM, **kwargs):
    r"""Potentiometry data fitting.

    Parameters:
        beta (:class:`numpy.ndarray`): 1D-array of floats representing the
            initial guess for the constants in log\ :sub:`10`\ units.
        beta_flags (sequence): refinement flags for values in
            'beta'.
        stoichiometry (:class:`numpy.ndarray`): The :term:`stoichiometry array`
            in the usual format.
        titration_data (iterable of dict): Information about the titrations.
            Each member of the list is one titrations. Every element in the
            list is a dict for which the accepted keys are:

            - V (:class:`numpy.ndarray`): volume of titre in mL.
            - emf (:class:`numpy.ndarray`): containing emf measured in mV. It
                must be of the same size than :var:`V`. Each element can be a
                float or a tuple of floats if more then one electrode is being
                used.
            - weights (:class:`numpy.ndarray`): The weights for each
                experimental point. It must be of the same size than
                :var:`emf`.
            - V0 (float): the initial volume in mL
            - T0 (sequence): total amount of initial components in mmol
            - buret (sequence): concentration in the buret for
              each component in mmol/mL.
            - Tflags (sequence, optional): The refinement flag for T0
            - buret_flags (sequence, optional): The refinement flag for emf0
            - T (:class:`numpy.ndarray`): An (N, S)-array with the initial
              concentrations. If provided, V, V0, T0, buret and Tflags are
              ignored.
            - error_V (float):

        electrodes (iterable of dict):

            - E0 (:class:`numpy.ndarray`): Value for standard potential for
              each curve in mV. If more than one electrode is used, tuples are
              accepted.
            - E0flags (:class:`numpy.ndarray`, optional): refinement flags for
              **E0** in the same format.
            - hindex (list of ints or tuples of ints): Index or indices for the
              electrode sensitive species (usually, the proton H :sup:`+`).
            - fRTnF (float or list of Nd-floats): by default 25.67968mV (1 atm,
              298K, 1 electron, f=1). If float, assume it is the same value
              for all sets and internally converted into a list.
            - error_emf (float):

        method (int): A number indicating the method to be used to mimize
            the objective functions. Allowed values are 0 for the algorithm
            by Levenberg-Marquardt (nonlinear least squares, this is the
            default one) and 1 for the Nelder-Mead algorithm (the simplex
            method).
        logbeta_flag (bool, optional[True]): False by default. If true B0 is
            assumed to be the log:sup:`10` of the constants to be refined.
        free_concentration_guess (sequence, optional): The free concentrations
            used as initial guess for concentration calculations.
        debug (bool, optional): returns a text with debug information
        verbosity (int, optional): An 0-2 number indicating the level of
            verbosity to be printed. 0 for mute, 1 for normal and 2 for
            pedanting output.
        htmlout (callable, optional): A callable function that accepts an html
            snippet to be output.

    Returns:
        tuple: First element is a 1D :class:`numpy.ndarray` of floats: The
        values of the refined parameters. The second element is a list of
        2D :class:`numpy.ndarray` of floats which are the final values of
        the free concentrations.
    """
    # +-------------+
    # | input check |
    # +-------------+
    libaux.assert_array_dim(1, beta)
    libaux.assert_same_len(beta, beta_flags)
    libaux.assert_array_dim(2, stoichiometry)

    # +----------+
    # | preamble |
    # +----------+

    _check_at_least_one_parameter(beta_flags, titration_data, electrodes)

    n_equilibria, n_species = stoichiometry.shape

    # htmlout = libaux.setkw(kwargs, 'htmlout', None)
    # verbosity = libaux.setkw(kwargs, 'verbosity', 2)

    # emf = tuple(curve['emf'] if isinstance(curve['emf'], np.ndarray)
    # else np.array(curve['emf']) for curve in titration_data)
    emf = tuple(curve['emf'] for curve in titration_data)
    # v_titre = tuple(_aryfy(titration, 'V') for titration in titration_data)

    # if htmlout:
    #     htmlout(report.html_start(emf, verbosity))

    def no_refine(k, d):
        return all((k not in i) or (i[k].count(0) == len(i[k])) for i in d)

    no_refine_analconc = no_refine('Tflags', titration_data)
    no_refine_emf0 = no_refine('E0flags', electrodes)

    # TODO in the future remove the simplified version, it is not good to have two
    # routines for the same. It makes more difficult to find bugs.

    # Case for simplified version
    # If only the betas are to be refined, proceed to the simplified version
    # ----------------------------------------------------------------------

    if no_refine_analconc and no_refine_emf0:
        emf_ = [t['emf'] for t in titration_data]
        # emf_ = [np.array(t['emf']) for t in titration_data]
        emf0_ = (t['E0'] for t in electrodes)
        nernst_ = (t['fRTnF'] for t in electrodes)
        reduced_emf = [build_reduced_emf(_emf_, _emf0_, _nernst_)
                       for _emf_, _emf0_, _nernst_
                       in zip(emf_, emf0_, nernst_)]
        full_emf = np.ma.concatenate(reduced_emf, axis=0)
        if full_emf.ndim == 1:
            full_emf = np.atleast_2d(full_emf).T

        n_exp_points = full_emf.shape[0]
        full_T = np.array(libaux.build_multiT_titr3(titration_data,
                                                    concat=True))
        # import pudb
        # pudb.set_trace()
        # remember the slices before concatenating
        lims = list(itertools.accumulate([0] + [s.shape[0] for s in emf_]))
        hindex = [t['hindex'] for t in electrodes]
        fhsel = functools.partial(hselect, hindices=hindex, slices=lims[:-1])

        # TODO change this
        # weights = tuple(titr['weights'] for titr in titration_data)
        # _weights = kwargs.pop('weights', n_exp_points*[1.0/n_exp_points])
        _weights = [curve['weights'] for curve in titration_data]
        weights = (np.concatenate([np.array(a) for a in _weights])).flatten()

        # fails if masked points
        assert len(weights) == n_exp_points

        key = 'free_concentration_guess', 'c0'
        if key[0] in kwargs:
            kwargs[key[1]] = np.concatenate(kwargs.pop(key[0]))

        libaux.assert_same_len(full_emf, weights)

        args = (beta*consts.LOGK, np.array(beta_flags), stoichiometry,
                full_emf, fhsel, full_T, weights)

        if method == consts.METHOD_LM:
            func = _emffitLM2
        elif method == consts.METHOD_NM:
            func = _emffitNM1
        else:
            raise ValueError("unknown fitting method")

        try:
            ln_beta, concentration_packed, alt = func(*args, **kwargs)
        except excepts.TooManyIterations as e:
            ret = e.last_value
            ret['last_value'] = ret['last_value'] / consts.LOGK
            conc = np.vsplit(ret['concentrations'], lims[1:-1])
            residuals = np.split(ret['residuals'], lims[1:-1])
            ret['concentrations'] = conc
            ret['residuals'] = residuals
            raise e

        if 'error' in alt:
            alt['error'] /= consts.LOGK
        assert concentration_packed.shape == (n_exp_points,
                                              n_species + n_equilibria)
        concentration = np.vsplit(concentration_packed, lims[1:-1])
        residuals = np.split(alt['residuals'], lims[1:-1])
        for c, e in zip(concentration, emf):
            assert c.shape == (len(e), n_equilibria+n_species), \
                   "%d != %d" % (c.shape[0], len(e))
        alt['residuals'] = residuals
    else:
        raise NotImplementedError

    return ln_beta/consts.LOGK, concentration, alt


def _emffitLM1(ln_beta, beta_flags, stoichiometry, emf, eletrodes, titration_data,
               weights, c0=None, verbosity=1, **kwargs):
    raise NotImplementedError
    # prepare x
    for parm, flags in ((ln_beta, beta_flags), 
                        (titration_data['T0'], titration_data['Tflags']),
                        (titration_data['buret'], titration_data['buret_flags']),
                        (electrodes['E0'], electrodes['E0flags'])):
        if all(i > 0 for i in flags):
            _f_args_.append(parm)

    x0 = _ravel_parameters(*_f_args_)
    # prepare y
    y = np.vstack(emf)

    # prepare f
    # TODO
    def _func_(x, free_conc):
        _lnbeta, _analytc, _danger = _unravel_parameters(x, _f_args_)
        full_T = np.array(libaux.build_multiT_titr3(_analytc, concat=True))
        # return calculated emf
        pass

    # prepare capping
    fcapping = functools.partial(libfit.max_ratio_capping, ratio=0.05)

    # prepare free_conc
    # TODO
    # prepare jacobian
    # TODO
    ln_beta, concentration_packed, alt = \
        libfit.levenberg_marquardt(x0, y, f, free_conc, jacobian, weights,
                                   capping=fcapping, **kwargs)


def _emffitLM2(ln_beta, beta_flags, stoichiometry, remf, fhsel, analytc,
               weights, c0=None, verbosity=1, **kwargs):
    x0 = np.fromiter(libaux.unravel(ln_beta, beta_flags), dtype=np.float)
    n_species = stoichiometry.shape[1]
    y = np.ravel(remf)
    f = _fit_fobj(fhsel)
    _initial_values = c0[:, :n_species] if c0 is not None else None

    # prepare capping
    fcapping = functools.partial(libfit.max_ratio_capping, ratio=0.1)

    free_conc = _fit_free_concentration(ln_beta, beta_flags, analytc,
                                        stoichiometry,
                                        initial_values=_initial_values)
    var = np.flatnonzero(beta_flags)
    jacobian = _fit_jac(stoichiometry, fhsel, var)
    ln_beta, concentration_packed, alt = \
        libfit.levenberg_marquardt(x0, y, f, free_conc, jacobian, weights,
                                   fcapping, **kwargs)
                                   # report=_report_, **kwargs)
    error_beta, covar, correl = libfit.final_params(alt['jacobian'], alt['weights'],
                                                    alt['residuals'])
    alt.update({'error_beta': error_beta, 'covariance': covar, 'correlation': correl})
    return ln_beta, concentration_packed, alt


def _emffitNM1(ln_beta, beta_flags, stoichiometry, remf, fhsel, analytc,
               weights, report=None, c0=None, max_iterations=20, term_x=1e-5,
               term_f=1e-7, verbosity=1, **kwargs):
    r"""Potentiometric data fitting by the Simplex algorithm.

    See `Numerical Recipes in C, §10.4
    <http://www.nrbook.com/a/bookcpdf/c10-4.pdf>`_
    and `http://www.scholarpedia.org/article/Nelder-Mead_algorithm`_

    Parameters:
        ln_beta (:class:`numpy.ndarray`): initial guess for the constants in
            the form of NATURAL log. It must be an 1D array of floats.
        beta_flags (:class:`numpy.ndarray`): refinement flags for values in
            logB0.  Shape must match that of **logB0** .
        stoichiometry (:class:`numpy.ndarray`): :term:`stoichiometry array`
        remf (:class:`numpy.ndarray`): reduced emf measured (dimmensionless).
            It must be an 1D-array of floats.
        hindex (int or list of Nd-ints): Index for the electrode sensitive
          species (usually, the proton H+).
        T ((Nt, S)-arrays of floats): total amounts of species in mmol
        weights (1D-array of floats): containing the values for weighting
        c0 (:class:`numpy.ndarray`): The free concentrations used as
            initial guess for concentration calculations.
        max_iterations (int): The number of iterations at which looping will
            be interrupted.
        weights (1D-array of floats): containing the values for weighting

    Returns:
        tuple:
        - :class:`numpy.ndarray`: The refined constants in natural logarithmic
            units
        - :class:`numpy.ndarray`: The free concentrations
        - dict: Extra optional parameters

    Raises:
        ValueError: If invalid parameters are passed.
    """
    x0 = np.fromiter(libaux.unravel(ln_beta, beta_flags), dtype=np.float)
    n_species = stoichiometry.shape[1]
    y = np.ravel(remf)
    f = _fit_fobj(fhsel)
    _initial_values = c0[:, :n_species] if c0 is not None else None
    free_conc = _fit_free_concentration(ln_beta, beta_flags, analytc,
                                        stoichiometry,
                                        initial_values=_initial_values)
    var = np.flatnonzero(beta_flags)
    jacobian = _fit_jac(stoichiometry, fhsel, var)
    ln_beta_out, concentration_packed, alt = \
        libfit.simplex(x0, y, f, free_conc, weights, report=report)
    final_jac = jacobian(None, concentration_packed)
    residuals = f(ln_beta_out, concentration_packed)
    error_beta, covar, correl = libfit.final_params(final_jac,
                                                    np.diag(weights),
                                                    residuals)
    alt.update({'error_beta': error_beta, 'covariance': covar,
                'correlation': correl, 'jacobian': final_jac,
                'residuals': residuals})
    # TODO final report
    return ln_beta_out, concentration_packed, alt


def emfsim(beta, stoichiometry, titration, electrodes, **kwargs):
    r"""Simulate emf titrations. Return calculated emf.

    Parameters:
        beta (:class:`numpy.ndarray`): 1D-array of floats representing the
            stability constants.
        stoichiometry (:class:`numpy.ndarray`): The :term:`stoichiometry array`
            in the usual format.
        titration (iterable of dict): Information about the titrations. Each
            member of the list is one titrations. See structure of this dict
            at :func:`emffit`.
        electrodes (iterable of dict): See structure of this dict
            at :func:`emffit`.
    Yields:
        calculated emf.
    """
    E0 = np.array(electrodes['E0'])[np.newaxis, :]
    fRTnF = np.array(electrodes['fRTnF'])[np.newaxis, :]
    for T in libaux.build_multiT_titr3(titration):
        C = libeq.consol.consol(beta, stoichiometry, T)
        H = C[:, electrodes['hindex']]
        emf = E0 + fRTnF*np.log(H)
        yield emf


def emf_res(electroactive_conc, emf, emf0, nernst):
    """Calculate the residuals for potentiometric data.

    The residuals are calculated as :math:`\\chi = emf_{meas}-emf_{calc}`.
    and :math:`emf_{calc}` is calculated through Nernst's equation.

    Parameters:
        electroactive_conc (:class:`numpy.ndarray`): a 1D array of floats
            representing the free concentrations of the electroactive species.
        emf (:class:`numpy.ndarray`): a 1D array of floats representing the
            experimental values for the potential
        emf0 (float): The :term:`standard potential`
        nernst (float): The slope for :term:`Nernst's equation`
    Returns:
        :class:`numpy.ndarray`: a 1D array of floats containing the residuals.
    """
    return emf - (emf0 + nernst * np.log(electroactive_conc))


# def emf_jac2(stoichiometry, free_conc, fh, titration, electrode):
#     """Return the full jacobian.
# 
#     Order of parameters:
# 
#         * Equilibrium constants
#         * starting ammounts
#         * buret concentrations
#         * standard potentials
# 
#     Parameters:
#         stoichiometry (:class:`numpy.ndarray`): the :term:`stoichiometry array`
#         free_conc (:class:`numpy.ndarray`): the free concentration array
#     """
#     # Bconst, Bref, Brestr = libaux.count_flags(Bkeys)
#     # Tconst, Tref, Trestr = libaux.count_flags(Tkeys)
#     # Econst, Eref, Erestr = libaux.count_flags(E0keys)
# 
#     N = sum(c.shape[0] for c in free_conc)
#     # Nvar = Bref + len(Brestr) + Tref + len(Trestr) + Eref + len(Erestr)
#     counted_flags = libaux.count_flags(Bkeys, Tkeys, E0keys)
#     nrefs = libaux.total_refine_params(count_flags)
# 
#     J1 = emf_jac1(stoichiometry, free_conc, fh, var)  # J1 = ∂E/∂β 
#     # J2 = ?? ∂E/∂t and ∂E/∂b 
#     # J3 = ?? ∂E/∂E0
#     # TODO complete calculation of J
#     return np.vstack((J1, J2, J3))


def dcdt(concentration, stoichiometry, starting_volume, titre):
    r"""Calculate ∂c/∂T.

    T is the analytical concentration which is defined as
    :math:`T_i = \frac{a_i v_0 + b_i v}{v_0 + v}`
    where :math:`a_i` is the initial ammount for the component *i* and
    :math:`b_i` is the buret concentration for the component *i*.

    Arguments:
        concentration (:class:`numpy.ndarray`): The full free concentration
            array (N, E+S).
        stoichiometry (:class:`numpy.ndarray`): The stoichiometry array.
        starting_volume (:class:`numpy.ndarray` or float): The starting volume
        titre (:class:`numpy.ndarray`): The titre volume
    Returns:
        tuple: where the first element is an (N, S, S) array where each element
            (i,j,k) is ∂c[i,j]/∂a[k]. The second element of the tuple is also
            an (N,S,S) array each element (i,j,k) are ∂c[i,j]/∂b[k].
    """
    n_species = stoichiometry.shape[1]
    superdelta = np.eye(n_species)[np.newaxis, ...]
    A = superdelta +  \
        np.einsum('ij,ik,...i->...jk', stoichiometry, stoichiometry,
                  concentration[..., n_species:]) \
        / concentration[..., np.newaxis, :n_species]
    v1 = (starting_volume/(starting_volume + titre))[:, np.newaxis, np.newaxis]
    b1 = v1*superdelta
    v2 = (titre/(starting_volume + titre))[:, np.newaxis, np.newaxis]
    b2 = v2*superdelta
    dcdt_ = np.linalg.solve(A, b1)
    dcdb = np.linalg.solve(A, b2)
    return dcdt_, dcdb


def emf_jac1(stoichiometry, free_conc, fh, var=()):
    r"""Calculate jacobian.

    Return jacobian for potentiometric titrations when
    T and E0 are constant, based on logarithmic betas and
    dimmensionless potentials.

    Notes: It is done this way because requires less input for the function
        and the result can be conveniently transformed by the calling function
        with simple operations. :math:`E* = E / (fRT/nF)`

    .. math::

        \frac{\partial E^*}{\partial \log \beta} =
        \frac1{C_h} \frac{\partial C_h}{\partial \log \beta} =
        \frac{\partial\log C_h}{\partial \log \beta}

    where :math:`\frac{\partial C_h}{\partial \log \beta}` has to be derived
    first from the set of linear equations

    .. math::

        \sum_{k=1}^S \left(
          \delta_{ki} + \sum_{j=1}^E {
            p_{ji} p_{jk} \frac{c_{j+S}}{c_k}
          }
        \right) \frac{\partial c_k}{\partial \log\beta_b}
        = -p_{bi}c_{b+S}

    Parameters:
        stoichiometry (:class:`numpy.ndarray`): stoichiometric coefficients
        free_conc (:class:`numpy.ndarray` of floats): free concentrations
        fh (callable): Function that applied to an array returns from it the
            electroactive species only.
    Returns:
        :class:`numpy.ndarray`: 2D-array of floats containing the jacobian
            in the form ∂emf/∂logβ.
    Raises:
        :class:`RuntimeError`: if jacobian couldn't be computed
    """
    libaux.assert_array_dim(2, stoichiometry, free_conc)

    dcdb = libeq.consol.dcdb(free_conc, stoichiometry)  # N,(E+S),E
    if isinstance(fh, int):
        full_jacobian = dcdb[..., fh, var]
    else:
        full_jacobian = fh(dcdb[..., var]).T

    return np.squeeze(full_jacobian)


def hselect2(array, hindices, slices):
    """Select columns that correspond to the electroactive species.

    Given the concentrations array, selects the columns that correspond
    to the electroactive species.

    Parameters:
        array (:class:`numpy.ndarray`): The :term:`free concentrations array`
        hindices (list): List of ints or list of lists of ints with the indices
            of the electroactive species. Example: [[0,1],[1,2],[3,4],[4,5]].
            hindices are applied along axis=0
        slices (list of ints): Where to divide C. Example: [ 0, 5, 10, 15 ]
            slices are applied along axis=1

    Returns:
        The part of C which is electroactive

    >>> slices = [0, 4, 7]
    >>> hindices = [[0,1],[1,2],[3,4]]
    >>> C = np.array([[ 0.255,  0.638,  0.898,  0.503,  0.418],
    ...               [ 0.383,  0.789,  0.731,  0.713,  0.629],
    ...               [ 0.698,  0.080,  0.597,  0.503,  0.456],
    ...               [ 0.658,  0.399,  0.332,  0.700,  0.294],
    ...               [ 0.534,  0.556,  0.762,  0.493,  0.510],
    ...               [ 0.637,  0.065,  0.638,  0.770,  0.879],
    ...               [ 0.598,  0.193,  0.912,  0.263,  0.118],
    ...               [ 0.456,  0.680,  0.049,  0.381,  0.872],
    ...               [ 0.418,  0.456,  0.430,  0.842,  0.172]])
    >>> hselect(C, hindices, slices)
    array([[0.255, 0.638], [0.383, 0.789], [0.698, 0.080], [0.658, 0.399],
           [0.556, 0.762], [0.065, 0.638], [0.193, 0.912], [0.381, 0.872],
           [0.842, 0.172]])
    """
    if slices is None and isinstance(int, hindices):
        return array[:, hindices, ...]

    if len(hindices) != len(slices):
        raise TypeError("hindices and slices have wrong size")
    # libaux.assert_array_dim(2, array)

    # slices → [ 0, 5, 10, 15 ]
    #          0→4 5→9  10→14  15→end
    # hindices → [[0,1],[1,2],[3,4],[4,5]]
    # 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
    # -------------  ------------- -------------- --------
    # 0  0  0  0  0  1  1  1  1  1  3  3  3  3  3  4  4  4
    # 1  1  1  1  1  2  2  2  2  2  4  4  4  4  4  5  5  5

    num_points = array.shape[0]
    nslices = (b-a for a, b in zip(slices, slices[1:]+[array.shape[0]]))
    # y = np.array(sum([n*h for n, h in zip(nslices, hindices)], []))
    y = np.vstack([np.tile(np.array(h), (n, 1)) for h, n in zip(hindices,
                                                                nslices)])
    return np.squeeze(array[np.arange(num_points), y.T, ...].T)


def build_reduced_emf(emf, emf0, nernst):
    """Build reduced EMF.

    Arguments:
        emf (:class:`numpy.ndarray`):
        emf0 (float):
        nernst (float):
    """
    # lims = libaux.get_lims(emf)
    # emfr = np.concatenate(emf)
    # ret = np.empty_like(emfr)
    # for lim1, lim2, emf0_, nernst_ in zip(lims[:-1], lims[1:], emf0, nernst):
    #     ret[lim1:lim2] = (emfr[lim1:lim2] - emf0_) / nernst_
    # return ret
    return (emf - emf0) / nernst


def fit_final_calcs(jacobian, resids, weights):
    """Perform final calculations common to some routines.

    Parameters:
        jacobian (:class:`numpy.array`): the jacobian
        resids (:class:`numpy.array`): the residuals
        weights (:class:`numpy.array`): the weights
    Returns:
        * the error in beta
        * the correlation matrix
        * the covariance matrix
    """
    # W = np.diag(weights)
    matrix_M = np.dot(np.dot(jacobian.T, np.diag(weights)), jacobian)
    error_B = np.diag(np.linalg.inv(matrix_M))
    # F = remf - np.log(H)
    covariance = libmath.covariance(jacobian, weights, resids)
    matrix_D = np.diag(covariance)
    lenD = len(matrix_D)
    correlation = covariance/np.sqrt(np.dot(matrix_D.reshape((lenD, 1)),
                                            matrix_D.reshape((1, lenD))))
    return error_B, correlation, covariance


def _fit_free_concentration(ln_beta, beta_flags, analytc, stoichiometry,
                            initial_values=None):
    if initial_values is not None:
        _initial_values = initial_values
    else:
        _initial_values = libeq.consol.initial_guess(np.exp(ln_beta),
                                                     stoichiometry, analytc)

    def func(x):
        nonlocal _initial_values
        gen = libaux.ravel(ln_beta, x, beta_flags)
        beta = np.exp(np.fromiter(gen, dtype=np.float))
        c = libeq.consol.consol(beta, stoichiometry, analytc,
                                _initial_values)
        _initial_values = c
        return c
    return func


def _fit_free_concentration_danger(ln_beta, beta_flags, titrations, electrodes,
                                   stoichiometry, initial_values=None):
    # TODO rewrite not to include *electrodes* which is not necessary to
    # calculate the concentrations.
    def _flatten(block, key):
        return [f for t in block
                for f in t.get(key, n_species*[consts.RF_CONSTANT])]

    _initial_values = None
    if initial_values:
        _initial_values = initial_values
    else:
        analytc = np.array(libaux.build_multiT_titr3(titrations, concat=True))
        _initial_values = libeq.consol.initial_guess(np.exp(ln_beta),
                                                     stoichiometry, analytc)

    n_beta, n_species = stoichiometry.shape
    t0_flags = _flatten(titrations, 'Tflags')
    t0_orig = _flatten(titrations, 'T0')
    titre = _flatten(titrations, 'V')
    starting_volume = _flatten(titrations, 'V0')
    buret_flags = _flatten(titrations, 'buret_flags')
    buret_orig = _flatten(titrations, 'buret')
    emf0_flags = _flatten(electrodes, 'E0flags')
    emf0_orig = _flatten(electrodes, 'E0')
    n_electro = len(emf0_flags)
    flags = np.concatenate((beta_flags, t0_flags, buret_flags, emf0_flags))
    x_original = np.concatenate((ln_beta, t0_orig, buret_orig, emf0_orig))

    def func(x):
        nonlocal _initial_values
        gen = libaux.ravel(x_original, x, flags)
        _lnbeta = itertools.islice(gen, n_beta)
        _t0 = itertools.islice(gen, n_species)
        _buret = itertools.islice(gen, n_species)
        _ = itertools.islice(gen, n_electro)
        aux = libaux.build_multiT_titr2(_t0, _buret, starting_volume, titre)
        analytc = np.concatenate([np.array(t) for t in aux], axis=0)
        c = libeq.consol.consol(np.exp(_lnbeta), stoichiometry, analytc,
                                _initial_values)
        _initial_values = c
        return c

    return func


def _fit_fobj(fhsel):
    def func(_, free_conc):
        """A function that accepts the values of *x0* as well as
        the free concentrations and return the calculated values for *y*.
        """
        electroactive = fhsel(free_conc)
        calc_remf = np.log(electroactive)
        return np.ravel(calc_remf)

    return func


def _fit_fobj_dangerous(fhsel):
    def func(_, free_conc):
        """A function that accepts the values of *x0* as well as
        the free concentrations and return the calculated values for *y*.
        """
        electroactive = fhsel(free_conc)
        calc_remf = -np.log(electroactive)
        # TODO complete with other parameters
        return np.ravel(calc_remf)

    return func


def _fit_jac(stoichiometry, fhsel, var):
    def func(_, free_conc):
        j = emf_jac1(stoichiometry, free_conc, fhsel, var)
        return j
    return func


def _fit_jac_dangerous():
    raise NotImplementedError

    def func(x, free_conc):
        pass
    return func


def _capping(ratio, capping=0.1):
    """Avoid variable change bigger than 'capping'.

    Given a ratio Δx/x, calculate the number it must be multiplied by so that
    the maximum of them change by exactly 'capping' (i.e. 10%).
    """
    return min(1.0, capping/np.amax(np.abs(ratio)))


def _check_at_least_one_parameter(beta_flags, titration_data, electrodes):
    all_flags = []
    all_flags.append(beta_flags)

    if isinstance(titration_data, dict):
        all_flags.append(titration_data.get('Tflags', []))
        all_flags.append(titration_data.get('buret_flags', []))
    else:
        for titration in titration_data:
            all_flags.append(titration.get('Tflags', []))
            all_flags.append(titration.get('buret_flags', []))

    if isinstance(electrodes, dict):
        all_flags.append(electrodes.get('E0flags', []))
    else:
        for electrode in electrodes:
            all_flags.append(electrode.get('E0flags', []))

    counted_flags = libaux.count_flags(*all_flags)
    nrefs = libaux.total_refine_params(counted_flags)
    if nrefs == 0:
        raise excepts.NothingToRefine


def _ravel_parameters(*args):
    x = []
    constraints = {}
    for values_flags in args:
        for value, flag in zip(*values_flags):
            if flag == consts.RF_IGNORE:
                continue
            elif flag == consts.RF_CONSTANT:
                continue
            elif flag == consts.RF_REFINE:
                x.append(value)
            else:
                if flag in constraints:
                    continue
                else:
                    x.append(value)
                    constraints[flag] = value
    return x


def _unravel_parameters(x, *args):
    xit = iter(x.tolist())
    constraints = {}
    retval = []
    for values_flags in args:
        current_vals = []
        for value, flag in zip(*values_flags):
            if flag == consts.RF_IGNORE:
                current_vals.append(value)
            elif flag == consts.RF_CONSTANT:
                current_vals.append(value)
            elif flag == consts.RF_REFINE:
                current_vals.append(next(xit))
            else:
                if flag in constraints:
                    current_vals.append(value*constraints[flag])
                else:
                    aux = next(xit)
                    current_vals.append(aux)
                    constraints[flag] = aux/value
        # print(current_vals)
        retval.append(current_vals.copy())
    return retval

