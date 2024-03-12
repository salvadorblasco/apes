#!/usr/bin/python3

import sys
import unittest

# import hypothesis as hp
# import hypothesis.strategies as st
import numpy as np

sys.path.append('../src/')

import libemf

import hexaprotic


class LibEmfTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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


if __name__ == '__main__':
    unittest.main()
