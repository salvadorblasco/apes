#!/usr/bin/python
# -*- encoding: utf-8 -*-

import random
import sys
import numpy
from pyss import simulation
from scipy import stats
import math

def simulate_with_error(B, sB, P, T, pT, N=50, NB=100, confidence=0.95, x=-1):
    """
    Arguments:
        B (list): list of floats being logarithms of beta contants.
        sB (list): list of floats of same size as B being the errors of B values
        P (array): array of ints of shape (E, S) containing the 
            stoichiometric coefficients.
        T (list)
        pT (list)
        N (int): size of simulation, defaults to 50
        NB (int): defaults to 100 and is number of points for Monte Carlo 
            simulation.
        confidence (float): is the confidence level to calculate the
            uncertainty limits. Defaults to 0.95, which is 95% of conficence.
            Must lie 0.0 < confidence < 1.0
        x (int): Indicates to the species whose concentration will be the
            independent variable. Its value must be a valid column of
            array P, and defaults to -1 (the last one).

    Returns:
        (array): shape (N,) are the values of the independent variable
        (array): shape (N, E+S-1) are the values of the concentrations
            calculated for the simulation not counting the errors.
        (array): shape (N, E+S-1) are lower limits for the concentrations
            calculated.
        (array): shape (N, E+S-1) are upper limits for the concentrations
            calculated.
    """

    #TODO use value N

    assert len(B) == len(sB)
    assert isinstance(P, ndarray)
    E, S = P.shape
    assert len(B) == E
    assert len(T) == len(pT) == S
    assert isinstance(confidence, float)
    assert 0.0 < confidence < 1.0

    # I. From B and sB, generate a collection of NB Bs with a normal 
    #   distribution.
    Bs = list()
    for b, sb in zip(B, sB):
        Bs.append(
            [ random.gauss(b, sb) for i in range(NB) ])

    # II. Randomly mix the collection of number generated
    new_B = [ 
        [
            Bs[j][i] for j in range(E)
        ] for i in range(NB)
    ] 

    # III. For each Bs of the collection, calculate concentrations.
    C_sim = list()
    x0 = numpy.array([[1e-19]])

    for i, b in enumerate(new_B):
        b_ = numpy.array([10**bb for bb in b]).reshape((E,1))

        h, c = simulation(b_, P, T, pT)
        C_sim.append(c.tolist())


    h, real_c = simulation(numpy.array([10**bb for bb in B]).reshape((E,1)), P, T, pT)

    C_resh = numpy.array(C_sim)     # betas, pH, species
    assert C_resh.shape[0] == NB
    assert C_resh.shape[1] == len(h)
    assert C_resh.shape[2] == S + E

    pH = [ -math.log10(hh) for hh in h ]
    _colors = 'b', '', 'g', 'r', 'c'

    lower_c = (1.0 - confidence)/2.0
    upper_c = 1 - lower_c

    for i in range(C_resh.shape[2]):                # loop over species
        q_upper = list()
        q_lower = list()
        for h in range(C_resh.shape[1]):            # loop over pH values
            c = C_resh[:,h,i]
            par_norm = stats.norm.fit(c, moments="mvs")
            f = stats.norm(loc=par_norm[-2],scale=par_norm[-1], moments="mvs")
            q_lower.append(f.ppf(lower_c))
            q_upper.append(f.ppf(upper_c))

    return h, real_c, numpy.array(q_lower), numpy.array(q_upper)


