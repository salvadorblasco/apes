#!/usr/bin/python3

import sys
import unittest

# import hypothesis as hp
# import hypothesis.strategies as st
import numpy as np
import numpy.testing as npt

import hexaprotic

sys.path.append('../src/')

import consts
import libemf



class LibEmfTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def test_emf_jac_beta(self):
        jtest = libemf.emf_jac_beta(hexaprotic.dlogc_dlogbeta[:,1,:6])
        npt.assert_allclose(jtest, hexaprotic.emf_jac, atol=1e-3)

    def test_emf_jac_init(self):
        mtest = libemf.emf_jac_init(hexaprotic.dlogc_dt)
        npt.assert_allclose(mtest, hexaprotic.dlogc_dt*consts.NERNST, rtol=1e-6)


if __name__ == '__main__':
    unittest.main()
