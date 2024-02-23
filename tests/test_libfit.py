#!/usr/bin/python3

import unittest

import sys
import numpy as np

sys.path.append('../src')


class TestLevenberg(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.logB = np.array([10.13, 19.53, 27.80, 34.82, -13.73])
        self.B = 10**self.logB
        self.P = np.array([[1, 1], [1, 2], [1, 3], [1, 4], [0, -1]])
        with np.load('pytrenal.npz') as d:
            self.real_C = d['C']
            self.real_T = d['T']
        self.E, self.S = self.P.shape

    def test_levenberg_marquardt(self):
        pass


class TestSimplex(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.logB = np.array([10.13, 19.53, 27.80, 34.82, -13.73])
        self.B = 10**self.logB
        self.P = np.array([[1, 1], [1, 2], [1, 3], [1, 4], [0, -1]])
        with np.load('pytrenal.npz') as d:
            self.real_C = d['C']
            self.real_T = d['T']
        self.E, self.S = self.P.shape

    def test_simplex(self):
        pass


class TestAuxFunctions(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def test_trivial_capping(self):
        from libfit import trivial_capping
        a = np.random.rand(10)
        b = np.random.rand(10)
        res = trivial_capping(a, b)
        np.testing.assert_allclose(res, a+b)

    def test_ratio_capping(self):
        from libfit import max_ratio_capping
        x = np.ones(5)
        dx = np.linspace(0.1, 0.5, 5)
        ret = max_ratio_capping(x, dx, 0.25)
        tru = np.array([1.1, 1.2, 1.25, 1.25, 1.25])
        np.testing.assert_allclose(ret, tru)


if __name__ == '__main__':
    unittest.main()

