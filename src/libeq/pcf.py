"""
Positive Continuous Method. [Marinoni17]
"""

import numpy as np

from libeq.cpluss import cpluss

import itertools                    ##
from libeq.fobj import fobj         ##


def pcf(concentration, beta, base_stoichiometry, analyticalc, tolerance=0.25):
    theta = None
    exponent = 1.0 / _exponent(base_stoichiometry)
    stoichiometry = _morel(base_stoichiometry)

    niter = itertools.count()           ##
    # print("   n     fsq       sump        sumr     theta")

    while True:
        species = cpluss(concentration, beta, base_stoichiometry, full=True)
        f = fobj(species, base_stoichiometry, analyticalc)         ##
        fsq = np.sum(f**2)                                         ##
        if (n:=next(niter))>200: break                                             ##

        sumr, sump = _sumps(species, analyticalc, stoichiometry)
        if _converged(sumr, sump):
            break
        theta = _update_theta(theta, sump, sumr)
        ratio = (sump / sumr)**exponent[None, ...]
        concentration = theta * concentration * ratio + (1-theta) * concentration

        # print(f"{n:4}  {fsq:.4e} {np.max(sump):.4e} {np.max(sumr):.4e} {np.max(theta):.4f} {np.min(ratio):.4e} {np.max(ratio):.4e}")  ##

    return concentration


def _converged(sumr, sump, threshold=1e-3):
    factor = np.abs(sumr-sump)/(sumr+sump)
    # threshold = 1e-3  # 1e-9 if first else 0.25
    # print("\t", np.min(factor), np.max(factor))
    return np.all(factor <= threshold)


def _exponent(stoichiometry):
    return np.array([min(i for i in np.nditer(s) if i > 0)
                     for s in np.hsplit(stoichiometry, stoichiometry.shape[1])])


def _morel(array):
    """Return extended stoichiometry array."""
    n_comp = array.shape[1]
    return np.vstack((np.eye(n_comp, dtype=np.int), array))


def _sumps(species, analyticalc, stoichiometry):
    pos = stoichiometry >= 0
    pstoich = np.zeros_like(stoichiometry)
    pstoich[pos] = stoichiometry[pos]
    nstoich = np.zeros_like(stoichiometry)
    nstoich[~pos] = np.abs(stoichiometry[~pos])

    sumrp = np.dot(species, pstoich)
    sumpp = analyticalc + np.dot(species, nstoich)
    sumrn = np.abs(analyticalc) + np.dot(species, pstoich)
    sumpn = np.dot(species, nstoich)

    tpos = analyticalc >= 0.0
    sumr = np.where(tpos, sumrp, sumrn)
    sump = np.where(tpos, sumpp, sumpn)
    return sumr, sump


def _update_theta(theta, sump, sumr):
    new_theta = _weighting(sump, sumr)
    # return np.where(theta > new_theta, theta, new_theta) if theta is not None else new_theta
    return new_theta


def _weighting(sump, sumr):
    ratio = sump/sumr
    ratio[ratio==0.0] = 1e-8
    alpha, beta = 0.9, 0.8
    theta = alpha - beta*np.where(ratio < 1.0, ratio, 1/ratio)
    return theta
