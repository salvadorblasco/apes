#!/usr/bin/python3

import sys
import unittest

import hypothesis as hp
import hypothesis.strategies as st
import numpy as np

sys.path.append('../src/')

import libemf


class LibEmfTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @unittest.skip('not implemented')
    def test_emffit(self):
        ...

    # @unittest.skip('not implemented')
    def test_emffitLM1(self):
        logbeta, stoich, electrode, titrats = self.hexaprotic()
        for t in titrats:
            t['weights'] = np.ones(len(t['emf']), dtype=np.float) 
            assert len(t['weights']) == len(t['emf'])

        import consts
        electrode[0]['E0flags'] = (consts.RF_CONSTANT,)
        beta_flags = (consts.RF_REFINE, consts.RF_REFINE, consts.RF_REFINE,
                      consts.RF_REFINE, consts.RF_REFINE, consts.RF_REFINE,
                      consts.RF_CONSTANT)
        import libemf
        for noise in (0.01, 0.02, 0.05):
            with self.subTest(noise=noise):
                inbeta = logbeta + np.random.normal(loc=0.0, scale=noise, size=len(logbeta))
                rebeta, _ = libemf.emffit(inbeta, beta_flags, stoich, titrats,
                                          electrode, method=consts.METHOD_LM)
                np.testing.assert_allclose(rebeta, logbeta)


    @unittest.skip('not implemented')
    def test_emffitNM1(self):
        ...

    def test_hselect(self):
        C = np.array([[0.255, 0.638, 0.898, 0.503, 0.418],
                      [0.383, 0.789, 0.731, 0.713, 0.629],
                      [0.698, 0.080, 0.597, 0.503, 0.456],
                      [0.658, 0.399, 0.332, 0.700, 0.294],
                      [0.534, 0.556, 0.762, 0.493, 0.510],
                      [0.637, 0.065, 0.638, 0.770, 0.879],
                      [0.598, 0.193, 0.912, 0.263, 0.118],
                      [0.456, 0.680, 0.049, 0.381, 0.872],
                      [0.418, 0.456, 0.430, 0.842, 0.172]])
        slices = [0, 4, 7]
        hindices = [[0, 1], [1, 2], [3, 4]]
        out = libemf.hselect(C, hindices, slices)
        out2 = np.array([[0.255, 0.638], [0.383, 0.789], [0.698, 0.080],
                         [0.658, 0.399], [0.556, 0.762], [0.065, 0.638],
                         [0.193, 0.912], [0.381, 0.872], [0.842, 0.172]])
        np.testing.assert_array_almost_equal(out, out2)
        # ---8<---
        slices = [0]
        hindices = [0]
        out = libemf.hselect(C, hindices, slices)
        out2 = np.array([0.255, 0.383, 0.698, 0.658, 0.534, 0.637, 0.598, 0.456, 0.418])
        np.testing.assert_array_almost_equal(out, out2)


    @unittest.skip('not implemented')
    def test_emfsim(self):
        pass

    @unittest.skip('not implemented')
    def test_emf_res(self):
        pass

    @unittest.skip('not implemented')
    def test_emfjac1(self):
        pass

    @unittest.skip('not implemented')
    def test_emfjac2(self):
        pass

    def test_amatrix(self):
        C = np.array([[(5**0.5-1)/2, (5**0.5-1)/2, (3-5**0.5)/2],
                      [2**0.5-1, 2**0.5, 2-2**0.5]])
        P = np.array([[1,1]])
        B = np.array([1])
        T = np.array([[1,1],[1,2]])
        # np.einsum('ji,jk,...j->...ik', Q, Q, C[None,...])
        Areal = np.array([[[1.        , 0.38196601],
                            [0.38196601, 1.        ]],
                           [[1.        , 0.58578644],
                            [0.58578644, 2.        ]]])
        n_species = P.shape[1]
        morel = np.vstack((np.eye(n_species, dtype=int), P))
        from libemf import _jac_amatrix
        Atest = _jac_amatrix(C, morel)
        np.testing.assert_allclose(Atest, Areal)


    def test_jac_b(self):
        A = np.array([[[1.        , 0.38196601],
                        [0.38196601, 1.        ]],
                       [[1.        , 0.58578644],
                        [0.58578644, 2.        ]]])
        C = np.array([[(5**0.5-1)/2, (5**0.5-1)/2, (3-5**0.5)/2],
                      [2**0.5-1, 2**0.5, 2-2**0.5]])
        P = np.array([[1,1]])
        good = np.array([[[-0.2763932 ],
                          [-0.2763932 ]],
                         [[-0.5       ],
                          [-0.14644661]]])
        from libemf import _jac_beta
        result = _jac_beta(A, C, P)
        np.testing.assert_allclose(result, good)

    def test_jac_t(self):
        A = np.array([[[1.        , 0.38196601],
                        [0.38196601, 1.        ]],
                       [[1.        , 0.58578644],
                        [0.58578644, 2.        ]]])
        v0 = 1.0
        v = np.array([0.0, 1.0])
        good = np.array([[[ 1.17082039, -0.4472136 ],
                          [-0.4472136 ,  1.17082039]],
                         [[ 0.60355339, -0.1767767 ],
                          [-0.1767767 ,  0.3017767 ]]])
        from libemf import _jac_t
        result = _jac_t(A, v, v0)
        np.testing.assert_allclose(result, good)

    def test_jac_b(self):
        A = np.array([[[1.        , 0.38196601],
                        [0.38196601, 1.        ]],
                       [[1.        , 0.58578644],
                        [0.58578644, 2.        ]]])
        v0 = 1.0
        v = np.array([0.0, 1.0])
        good = np.array([[[ 0.        ,  0.        ],
                          [ 0.        ,  0.        ]],
                         [[ 0.60355339, -0.1767767 ],
                          [-0.1767767 ,  0.3017767 ]]])
        from libemf import _jac_b
        result = _jac_b(A, v, v0)
        np.testing.assert_allclose(result, good)

    def test_build_reduced_emf(self):
        emf = 300*(np.random.rand(50) - 0.5)
        emf0 = 386
        nernst = 25.05
        remf1 = (emf - emf0) / nernst
        remf2 = libemf.build_reduced_emf(emf, emf0, nernst)
        np.testing.assert_almost_equal(remf1, remf2)

        emf = 300*(np.random.rand(50,3) - 0.5)
        emf0 = 400*np.random.rand(3)
        nernst = 25.05
        remf1 = (emf - emf0) / nernst
        remf2 = libemf.build_reduced_emf(emf, emf0, nernst)
        np.testing.assert_almost_equal(remf1, remf2)

    @unittest.skip('not implemented')
    def test_fit_final_calcs(self):
        pass

    @unittest.skip('not implemented')
    def test_capping(self):
        pass

    @unittest.skip('FIX')
    def test_fit_free_concentration(self):
        import consts
        import equilsystems
        beta_flags = np.array((consts.RF_REFINE, consts.RF_CONSTANT))
        real_c, P, B, T = equilsystems.generic_monoprotic()
        f = libemf._fit_free_concentration(np.log(B), beta_flags, T, P)
        for x in (8.0, 9.0, 10.0):
            c = f(np.array((consts.LOGK*x,)))
            # np.testing.assert_array_almost_equal(c, real_c)
            self.assertTupleEqual(c.shape, real_c.shape)

    @unittest.skip('FIX')
    def test_fit_fobj(self):
        f = libemf._fit_fobj(lambda x: x[:,1])
        import equilsystems
        real_c, _, _, _ = equilsystems.generic_monoprotic()
        fobj = f(None, real_c)
        np.testing.assert_array_almost_equal(fobj, -np.log(real_c[:,1]))
        self.assertEqual(fobj.ndim, 1)

    @unittest.skip('FIX')
    def test_fit_jac(self):
        import functools
        import consts
        import equilsystems
        beta_flags = np.array((consts.RF_REFINE, consts.RF_CONSTANT))
        var = np.flatnonzero(beta_flags)
        real_c, P, B, T = equilsystems.generic_monoprotic()
        x = np.array((9.75,))
        fhsel = functools.partial(libemf.hselect, hindices=(1,), slices=[0,])
        jac = libemf._fit_jac(P, fhsel, var)
        j = jac(x, real_c)
        self.assertEqual(len(j), len(real_c))

    def test_ravel_parameters(self):
        import consts
        with self.subTest():
            set1 = ((1.0, 2.0, 3.0),
                    (consts.RF_REFINE, consts.RF_CONSTANT, consts.RF_REFINE))
            set2 = ((1.0, 2.0),
                    (consts.RF_CONSTANT, consts.RF_REFINE))
            xout = libemf._ravel_parameters(set1, set2)
            xres = [1.0, 3.0, 2.0]
            self.assertSequenceEqual(xout, xres)

        with self.subTest():
            set1 = ((1.0, 2.0, 3.0, 4.0),
                    (consts.RF_REFINE, consts.RF_CONSTANT, consts.RF_CONSTRAINT1, consts.RF_CONSTRAINT1))
            set2 = ((1.0, 2.0),
                    (consts.RF_CONSTRAINT2, consts.RF_REFINE))
            set3 = ((1.0, 2.0),
                    (consts.RF_CONSTRAINT2, consts.RF_REFINE))
            xout = libemf._ravel_parameters(set1, set2, set3)
            xres = [1.0, 3.0, 1.0, 2.0, 2.0]
            self.assertSequenceEqual(xout, xres)

    def test_unravel_parameters(self):
        import consts
        with self.subTest():
            set1 = ((1.0, 2.0, 3.0),
                    (consts.RF_REFINE, consts.RF_CONSTANT, consts.RF_REFINE))
            set2 = ((1.0, 2.0),
                    (consts.RF_CONSTANT, consts.RF_REFINE))
            xin = np.array([1.5, 3.5, 2.5])
            xout = libemf._unravel_parameters(xin, set1, set2)
            xres = ((1.5, 2.0, 3.5), (1.0, 2.5))
            for a, b in zip(xout, xres):
                self.assertSequenceEqual(a, b)

        with self.subTest():
            set1 = ((1.0, 2.0, 3.0, 4.0),
                    (consts.RF_REFINE, consts.RF_CONSTANT, consts.RF_CONSTRAINT1, consts.RF_CONSTRAINT1))
            set2 = ((1.0, 2.0),
                    (consts.RF_CONSTRAINT2, consts.RF_REFINE))
            set3 = ((1.0, 2.0),
                    (consts.RF_CONSTRAINT2, consts.RF_REFINE))
            xin = np.array([1.5, 3.5, 1.5, 2.5, 2.5])
            xout = libemf._unravel_parameters(xin, set1, set2, set3)
            xres = ((1.5, 2.0, 3.5, 3.5*4/3), (1.5, 2.5), (1.5, 2.5))
            for a, b in zip(xout, xres):
                self.assertSequenceEqual(a, b)

    def hexaprotic(self):
        logbeta = np.array([10., 18., 24., 28., 31., 33., -13.77])
        stoich = np.array([[1, 1], [1, 2], [1, 3], [1, 4], [1, 5], [1, 6], [0, -1]])
        electrode = {'E0': 405.0, 'fRTnF': 25.692861549, 'hindex': 1}
        titr1 = {'V': np.linspace(0.00, 12.00, 101),
                 'emf': np.array((265.6, 263.5, 261.4, 259.1, 256.7, 254.2,
                                  251.6, 248.8, 245.9, 242.9, 239.7, 236.3,
                                  232.8, 229.1, 225.2, 221.1, 216.9, 212.4,
                                  207.8, 203.0, 197.9, 192.7, 187.2, 181.5,
                                  175.5, 169.2, 162.4, 155.0, 146.7, 137.1,
                                  125.6, 112.2, 98.0, 84.9, 73.6, 63.6, 54.2,
                                  45.1, 35.7, 25.5, 14.0, 0.6 , -14.0, -27.9,
                                  -40.1, -50.6, -60.2, -69.4, -78.6, -88.4,
                                  -99.2, -111.6, -125.6, -139.7, -152.3, -163.0,
                                  -172.5, -181.1, -189.3, -197.2, -205.1, -212.9,
                                  -220.4, -227.5, -234.0, -239.7, -244.8, -249.2,
                                  -253.0, -256.5, -259.5, -262.3, -264.8, -267.1,
                                  -269.1, -271.0, -272.9, -274.5, -276.0, -277.5,
                                  -280.2, -283.7, -285.8, -287.6, -289.4, -290.2,
                                  -292.5, -294.6, -296.5)),
                  'V0': 25.0,
                  'T0': (5.0, 30.0),
                  'buret': (0.0, -0.1)}
        titr2 = {'V': np.linspace(0.00, 6.58, 95),
                 'emf': np.array((248.8, 246.4, 243.8, 241.1, 238.1, 234.9,
                                  231.5, 227.8, 223.8, 219.5, 214.8, 209.7,
                                  204.2, 198.2, 191.7, 184.6, 176.9, 168.5,
                                  159.1, 148.1, 134.5, 116.8, 96.6, 78.5, 63.5,
                                  50.0, 36.4, 21.3, 2.9, -18.2, -37.2, -52.7,
                                  -66.3, -79.7, -94.4, -111.8, -131.9, -150.3,
                                  -164.9, -177.0, -187.5, -197.2, -206.0, -214.1,
                                  -221.2, -227.5, -232.9, -237.7, -241.8, -245.4,
                                  -248.6, -251.5, -254.2, -256.5, -258.7, -260.7,
                                  -262.6, -264.3, -265.9, -267.5, -268.9, -270.3,
                                  -271.5, -272.7, -273.9, -276.0, -278.0, -280.7,
                                  -283.1, -285.3, -287.2, -290.1, -292.3, -294.6)),
                  'V0': 25.0,
                  'T0': (2.0, 12.0),
                  'buret': (0.0, -0.1)}
        titr3 = {'V': np.linspace(0.00, 9.00, 101),
                 'emf': np.array((261.6, 259.8, 257.9, 255.9, 253.8, 251.6,
                                  249.4, 247.0, 244.4, 241.8, 239.1, 236.2,
                                  233.2, 230.0, 226.7, 223.2, 219.5, 215.7,
                                  211.7, 207.6, 203.2, 198.6, 193.9, 189.0,
                                  183.9, 178.5, 172.8, 166.8, 160.4, 153.3,
                                  145.4, 136.3, 125.5, 113.0, 99.7, 87.3, 76.3,
                                  66.6, 57.7, 49.0, 40.4, 31.3, 21.3, 10.0,
                                  -3.0, -16.8, -29.5, -40.8, -50.6, -59.6, -68.2,
                                  -76.9, -85.9, -95.6, -106.7, -119.3, -132.6,
                                  -145.2, -156.2, -165.7, -174.2, -182.0, -189.5,
                                  -196.7, -203.8, -210.7, -217.3, -223.6, -229.4,
                                  -234.6, -239.4, -243.5, -247.2, -250.6, -253.6,
                                  -256.3, -258.8, -261.0, -263.1, -265.0, -266.8,
                                  -268.5, -270.1, -273.0, -275.5, -277.8, -278.9,
                                  -281.8, -283.6, -285.3, -287.5, -288.2)),
                 'V0': 25.0,
                 'T0': (4.0, 24.0),
                 'buret': (0.0, -0.1)}
        return logbeta, stoich, (electrode,), (titr1, titr2, titr3)


if __name__ == '__main__':
    unittest.main()
