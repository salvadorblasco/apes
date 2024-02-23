#!/usr/bin/python3

"""
Variables:
- E (int) → the number of equilibria
- S (int) → the number of independent components
- N (int) → the number of experimental points
- Nn (int) → the number of nuclei per point
- d_o ((N,Nn)-array) → the experimental values for the chemical shift, δ 
- x ((N, Nn)-array) → the fractional population of nuclei
- w ((N, Nn)-array) → the weights for each point
- Pr (list) → A list of Nn ints with the correspondence nuclei/species.
    Given an index representing a nucleus, the species corresponding would
    be Pr[nucleus] and given an index representing a species, the nuclei
    corresponding to that species would be [i for i in Pr if i==species]
    NOTICE: some species might be NMR silent and no be indexed here
 
"""

import numpy as np

def f_obj(d_o):
    # d_c = 
    np.sum(w * (d_c - d_o)**2)

def jacobian(C, P):
    dcdb = libeq.dcdb(C, P)     # d log c / d log B

def nmrfit(delta_obs, T, P, B0, B_mask, d0, d_mask, X):
    """See from C. Frassineti, S. Ghelli, P. Gans, A. Sabatini,
        M. S. Moruzzi and A. Vacca, Anal. Biochem. 1995, 231, 374-382

    Function to minimise:
        U = Σ_j w_j(δ_calc - δ_obs)
        δ_calc = Σ_j(f_j.d_j) where 
            d_j chemical shift of a particular nucleus in the various species
                present.
        f_j = x_j . C_j/T_x        
            T_x is the total concentration of reagent containing the nucleus
                under consideration
            x_j is the stoichiometric coefficient X in the j_th species

    Input:
    (1) delta_obs -> array of N float numbers representing the experimental data
    """
 
def c2f(C, T, P, Pr):
    """Calculates the fractional population for the nuclei, defined as
    :math:`f_j = x_j \cdot C_j/T_x` where :math:`T_x` is the total
    concentration of reagent containing the nucleus
    under consideration, :math:`x_j` is the stoichiometric coefficient of
    reagent *x* in the *j* th species.

    Parameters:
        C (:class:`numpy.ndarray`): The :term:`free concentrations array`
        T (:class:`numpy.ndarray`): The :term:`total concentrations array`
        P (:class:`numpy.ndarray`): The :term:`stoichiometry array`

    Returns:
        :class:`numpy.ndarray`: An (N, Nn) array with the fractional
            populations of nuclei.
    """

    CC = expandC(C, Pr)
    f = np.empty(N, Nn)
    #for c in CC.T:
        


def expandC(C, Pr):
    """Convert free concentrations of species to concentration of nuclei
    """

    Nn = len(Pr)
    N = C.shape[0]

    CC = np.empty(N, Nn)
    for c, n in enumerate(Pr):
        CC[:, n] = C[:, c]

    return CC
