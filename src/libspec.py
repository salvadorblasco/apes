#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Module **libspec** contains the routines needed for fitting equilibrium
constants from spectrophotometric data, either absorbance, fluorescence,
phosphorence, circular dichroism or any other that follows a linear
relationship such as 
`Beer-Lambert law <https://en.wikipedia.org/wiki/Beer%E2%80%93Lambert_law>`_.

.. glossary::

    Ns
        is the number of spectra

    Nl
        is the number of wavelengths (points per spectrum)

    Ne
        is the total number of non-transparent, refinable species

    Nc
        is the total number of non-transparent species

    Nb
        is the number of constants to be refined

    N
        is the total number of experimental points. **N** = **Ns** * **Nl**

    A
        is the absorbance matrix. It is a 2D matrix of size (**Ns**, **Nl**)  where 
        columns represent the wavelength and rows are experimental points. Very often
        this matrix must be flattened for calculations. It must be flattened column-wise.

    e
        is the epsilon matrix. It is a (Nl, E+S) array of floats representing
        the extintion coefficient.

    erefkey
        is an **E** + **S**-sized array of ints with the key to refine
        or not the epsilon.

    transp
        is an **E** + **S**-sized list of bools that indicate whether a
        specific species is transparent or not
"""

import functools

import numpy as np

import libeq

# LOGK = 2.3025851
__version__ = 'dev'


def specfit(logB0, Bflags, P, spectra, T, **kwargs):
    """Routine for spectrometry data fitting.

    Parameters:
        logB (:class:`numpy.ndarray`): 1D-array of floats representing the
            initial guess for the constants in the form of log\ :sub:`10`\ .
        Bflags (:class:`numpy.ndarray`): refinement flags for values in
            logB0. Flags are int values.
        P (:class:`numpy.ndarray`): The :term:`stoichiometry array` in the
            usual format.
        spectra (:class:`numpy.ndarray`): containing the values of the
            measured spectra in normalized format. It is an (N, Nl)-array
        T (list of :class:`numpy.ndarray`): Total amounts of the free species
            in mmol for every titration point.
        method (int): A number indicating the method to be used to mimize
            the objective functions. Allowed values are 0 for the algorithm
            by Levenberg-Marquardt (nonlinear least squares, this is the
            default one) and 1 for the Nelder-Mead algorithm (the simplex
            method).
        logB_flag (bool, optional): False by default. If true B0 is assumed
            to be the log10 of the constants to be refined
        c0 (list of :class:`numpy.ndarray`, optional): The free concentrations
            used as initial guess for concentration calculations.
        debug (bool, optional): returns a text with debug information
        mask (bool, optional): a list of boolean arrays indicating for each
            experimental point whether will be used or not.
        weights (list of 1D-array of floats, optional) containing the values
            for weighting.
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
    # ravel parameters
    #  * prepare x0
    x = libaux.unravel(logB0, Bflags)

    # transform parameters
    #  * prepare y
    y = spectra
    #  * prepare f
    #  * prepare free_conc
    #  * prepare jacobian
    #  * prepare weights

    # call libfit.levenberg_marquardt
    x, concs, ret_extra = libfit.levenberg_marquardt(x0, y, f, free_conc, jacobian, weights)

    # unravel parameters
    # return result
    return


def spec_function(free_concentrations, optical_activity, optical_path, baseline=0.0):
    """Calculate the optical interaction.

    .. math:: A_{i\lambda} = l \sum_{j=1}^{E+S} \varepsilon_{\lambda, j} c_{ij} + B
    """   
    return optical_path * free_concentrations @ optical_activity.T) + baseline


# TODO delete
def specLM1(A, e):
    """This routine performs the fitting of the constants by means of the
    Levenberg-Marquardt algorithm when only the contants are to be refined.

    The objective function is defined as

    .. math::

        f = \sum_{i,\lambda}\left( A_{i\lambda} - 
        l \sum_{j=1}^{E+S} \\varepsilon_{\lambda, j} c_{ij}
        (\\beta_1,\\beta_2,\ldots,\\beta_E) \\right)^2

    where we refine ε and β together.

    Arguments:
        A (:class:`numpy.ndarray`): The absorbance matrix. It must
            be (Ns, Nl)-sized.
        e (:class:`numpy.ndarray`): The epsilon matrix. It must
            be (Nl, Nc)-sized.

    Raises:
        ValueError: If invalid parameters are passed.
    """

    # calculate initial free concentrations
    # compute χ₂(dx)

    if 'damping' in kwargs:
        damping = kwargs.pop('damping')
    else:
        damping = 0.001

    x = libaux.unravel(logB0, Bflags)  # TODO unravel x
    J = jacobian()
    M = np.dot(np.dot(J.T, W), J)
    D = np.diag(np.diag(M))

    while True:
        try:
            dx = np.linalg.solve(M+damping*D, np.dot(np.dot(J.T, W), F.T))
        except np.linalg.linalg.LinAlgError:
            damping *= 10
            continue

        try:
            C, H = CH(libaux.ravel(logB0, x+dx, Bflags), fcalcC, hindex)
        except:
            diagn = {'iteration': iterations, 'x': x, 'dx': dx, 'J': J,
                     'chisq': chisq, 'residuals': F, 'damping': damping}
            raise excepts.FailedCalculateConcentrations(**diagn)

        # TODO calc residuals
        F = remf - np.log(H)

        if qt_handle is not None:
            qt_handle("iter %d. chisq = %f" % (iterations, chisq))

        if np.sum(F**2) >= chisq:
            damping *= 10
        else:
            if htmlout:
                htmlout(_html_iteration(iterations, x/consts.LOGK, dx/consts.LOGK, chisq))

            iterations += 1
            damping /= 10
            chisq = np.sum(F**2)
            out_chisq.append(chisq)
            x += dx
            if 'one_iter' in kwargs:
                break
            J = emf_jac1(P, C, hindex, var)
            M = np.dot(np.dot(J.T, W), J)
            D = np.diag(np.diag(M))

        if np.all(np.abs(dx)/x < threshold):
            break

        if iterations > max_iterations:
            if htmlout and verbosity:
                htmlout("Maximum number of iterations has been reached." +
                        "Fitting will stop now.")

    error_B = np.diag(np.linalg.inv(M))
    CV = covariance(J, weights, F)
    D = np.diag(CV)
    nD = len(D)
    CR = CV/np.sqrt(np.dot(D.reshape((nD, 1)), D.reshape((1, nD))))

    if htmlout and verbosity:
        htmlout(_html_finalconvg(x/consts.LOGK, error_B/consts.LOGK, CR))

    if supyquad and verbosity:
        _spyq_finalconvg(x, error_B, CR)

    ret_extra = {'error': libaux.ravel(np.zeros(E), error_B, Bflags),
                 'damping': damping,
                 'covariance': CV,
                 'correlation': CR}

    return libaux.ravel(logB0, x, Bflags), C, ret_extra


# TODO delete
def jacobian(A, e, C):
    """Returns the jacobian matrix. For constant fitting from spectroscopic
    data it is an (N×Nl,Nl×Nb)-array.

    Parameters:
        A (:class:`numpy.ndarray`): The absorbance matrix. It must
            be (Ns, Nl)-sized.
        e (:class:`numpy.ndarray`): The epsilon matrix. It must
            be (Nl, Nc)-sized.
        ntrares (:class:`numpy.ndarray`): An **E**+**S**-sized list of bools
            that indicate whether a specific species is transparent or not
        C (:class:`numpy.ndarray`): The free concentrations array. It must
            be (Ns, E+S)-sized.
    """

    Nl = e.shape[0]
    j1 = jac_sub_eps(Nl, C)
    j2 = jac_sub_beta(e, C, P)
    assert j1.shape[0] == j2.shape[0]
    return np.concatenate((j1, j2))


def jac_sub_eps(Nl, C):
    """returns the part of the jacobian that depends on the extinction
    coefficient refinement, :math:`J_{\\varepsilon}`. 

    .. math::

        \\frac{\partial A_{i\lambda}}{\partial\\varepsilon_{ab}} =
          l c_{ib} \delta_{\lambda a}

    Arguments:
        Nl (int): is the number of wavelengths used (points per spectrum)
        C (:class:`numpy.ndarray`): The free concentrations coefficient array.
            It must be an (N, Ne) array.

    Returns:
        The part of the jacobian matrix that depends on the extintion
            coefficient.
    """

    # Ne is the total number of non-transparent, refinable species
    # Nc is the total number of non-transparent species
    # N  is the total number of experimental points. N = Ns * Nl

    N, Ne = C.shape
    J = np.zeros(Nl*N, Nl*Ne)
    for wl in range(Nl):        # iterate over wavelengths
        r0 = wl*N
        r1 = (wl+1)*N
        c0 = wl*Ne
        c1 = (wl+1)*Ne
        J[r0:r1, c0:c1] = C
    return J


def jac_sub_beta(e, C, P):
    """returns the part of the jacobian that depends on the equilibrium
    constants, :math:`J_{\\beta}`. 

    .. math::

        \\frac{\partial A_{i\lambda}}{\partial\log\\beta_k} =
          l \sum_{j=1}^{E+S} \\varepsilon_{\lambda j} c_{ij}
          \\frac{\partial\log c_{ij}}{\partial\log\\beta_k}

    Arguments:
        e (:class:`numpy.ndarray`): The epsilon matrix. It must
            be (Nl, Nc)-sized.
        C (:class:`numpy.ndarray`): The free concentrations coefficient array.
            It must be an (N, E+S) array.

    Returns:
        The part of the jacobian matrix that depends on the extintion
            coefficient.
    """

    # TODO define Nb
    J = np.zeros(N, Nb)
    dcdB = numpy.stack(libeq.iterdcdB(C, P))
    assert dcdB.shape == (N, S, E)

    J = np.einsum('ij,kj,kjl->ikl', (e, C, dcdB))
    assert J.shape == (Nl, N, Nb)
    return J


def ravel_A(A):
    """This function takes all A data and flattens it.

    Parameters:
        A (:class:`numpy.ndarray`): An 2D or 3D matrix where dim 1 is the
            experimental point, dim 2 is the wavelength and dim 3 (optional)
            is the replica.

    Returns:
        :class:`numpy.ndarray`: A flattened array with all the valid data.
    """
    
    pass


def unravel_A(A):
    """This function undoes what :func:`ravel_A` does.

    Parameters:
        A (:class:`numpy.ndarray`): A flattened array with all the valid data.

    Returns:
        :class:`numpy.ndarray`: An 2D or 3D matrix where dim 1 is the
            experimental point, dim 2 is the wavelength and dim 3 (optional)
            is the replica.
    """
    pass
