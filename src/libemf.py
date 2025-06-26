# -*- coding: utf-8 -*-

"""Routines for constant fitting from potentiometric data.

Module :mod:`libemf` contains the routines needed for fitting equilibrium
constants from potentiometric data. The main function to invoke here is
:func:`emffit` which handles most of the work. This is the only public
function for this module.

.. rubric:: Functions included in this module
"""


import numpy as np
from numpy.typing import NDArray

import consts

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


def emf_jac_beta(dlogc_dlogbeta: NDArray[float], slope=1.0, temperature: float=298.15) -> NDArray[float]:
    r"""Calculate the jacobian part related to equilibrium constants.

    The calculation is done according to equation
    .. math ::

    \frac{\partial E_n}{\partial\log\beta_b}=\frac{fRT}{nF\log e}\frac{\partial\log c_{nh}}
    {\partial\log\beta_b} 

    Parameters:
        dlogc_dlogbeta (:class:`numpy.ndarray`): the derivative values. They can be obtained
            from :func:`libeq.jacobian.dlogcdlogbeta`.
        slope (:class:`numpy.ndarray`): The slope for :term:`Nernst's equation`
        temperature (float): the absolute temperature
    Returns:
        :class:`numpy.ndarray`: an array of floats containing the calculated values
    """
    nernstian_slope = slope*consts.RoverF*temperature
    return nernstian_slope*dlogc_dlogbeta


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


def emf_weights(titre: NDArray[float], titre_error: float,
                emf: NDArray[float], emf_error: float) -> NDArray[float]:
    gradient = np.gradient(emf, titre, axis=0)
    return 1/(emf_error**2 + gradient**2 * titre_error**2)


def residual_jacobian(emf: NDArray[float], calc_emf: NDArray[float],
                      weights: NDArray[float], demfdx):
    # breakpoint()
    aux = np.sqrt(weights)*(emf - calc_emf)*emf
    return -2*np.sum(aux[:,None,None]*demfdx, axis=0)
