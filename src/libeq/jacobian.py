import numpy as np


def jacobian(concentration, stoichiometry, logc=False):
    r"""Compute the jacobian array.

    This function computes the jacobian for the function :func:`fobj`,
    which is defined as

    .. math:: J = \left( \begin{array}{ccc}
       \frac{\partial f_0}{\partial c_0} & \cdots &
       \frac{\partial f_0}{\partial c_S} \\
       \vdots  & \ddots & \vdots  \\
       \frac{\partial f_N}{\partial c_0} & \cdots &
       \frac{\partial f_N}{\partial c_S} \\
      \end{array} \right)

    where :math:`f_i = c_i + \sum_{j=1}^E p_{ij}c_{j+S} - T_i`
    and therefore

    .. math:: J_{ij} = \delta_{ij} + c_j^{-1} \sum_{k=1}^E {p_{ki} p_{kj}
       c_{k+S}}

    Parameters:
        concentration (:class:`ndarray`): the :term:`free concentrations array`
            for every component.  It must be an (*N*, *E* + *S* )-sized array
            of floats.
        stoichiometry (:class:`ndarray`): The :term:`stoichiometry array`.
            It must be an (*E*, *S*)-sized array.
        log (bool): If True, the returned result will be
            :math:`J_{ij} = \frac{\partial f_i}{\partial\log c_j}`. If False
            (default) the returned result will be
            :math:`J_{ij} = \frac{\partial f_i}{\partial c_j}`
    Returns:
        :class:`ndarray`: An (*E*, *E*)-sized array which is the jacobian
            matrix.
    """
    n_species = stoichiometry.shape[1]
    aux1 = np.einsum('ij,ik,li->ljk', stoichiometry, stoichiometry,
                     concentration[:, n_species:])
    aux2 = np.eye(n_species)
    aux3 = concentration[:, np.newaxis, :n_species]
    if logc:
        return aux2 * aux3 + aux1
    else:
        return aux2 + aux1/aux3


def jacobian_solid(concentration, stoichiometry, solubility_stoich, solubility_produc):
    jf1 = jacobian(concentration, stoichiometry)
    jf2 = jacobian_f_solid(solubility_stoich)
    jg1 = jacobian_g_solid(concentration, solubility_stoich, solubility_product)
    z = np.zeros((jf1.shape[0], jf2.shape[1], jg1.shape[2]))
    return np.block([[jf1, jf2],[jg1, z]])


def jacobian_f_solid(solubility_stoich):
    return solubility_stoich.T


def jacobian_g_solid(concentration, solubility_stoich, solubility_product):
    g = fobj.gobj(concentration, solubility_stoich, solubility_product)
    return (1 + g[..., None])*solubility_stoich[None,...]/concentration[:,None,:]
