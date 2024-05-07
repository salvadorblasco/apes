#!/usr/bin/python3

import unittest
import unittest.mock
import sys

import numpy as np
from PyQt5 import QtWidgets

import _syntheticdata

sys.path.append('../src')

import consts
# import libeq.consol


class TestBridge(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.data, self.params = _syntheticdata.load_hexaprotic()

        self.pconsol = unittest.mock.patch("libeq.consol.consol", new=unittest.mock.MagicMock(return_value=self.data.free_concentration))
        self.pinitgu = unittest.mock.patch("libeq.consol.initial_guess", new=unittest.mock.MagicMock(return_value=self.data.free_concentration))

    def setUp(self):
        from bridge import Bridge
        self.bridge = Bridge(self.params)

    def test_dimmensions(self):
        self.assertTupleEqual(self.bridge.jacobian.shape, (89, 6))
        self.assertTupleEqual(self.bridge.residual.shape, (89, ))

    def test_titrations(self):
        jac, res = self.bridge.build_matrices()

        for cc in self.params.titrations.values():
            np.testing.assert_allclose(cc.free_conc, self.data.free_concentration, atol=4e-3)
            np.testing.assert_allclose(cc.amatrix, self.data.matrix_a, atol=0.2)

    def test_matrices(self):
        # import libeq.consol
        # libeq.consol.consol = unittest.mock.MagicMock(return_value=hexaprotic.free_concentration)
        # libeq.consol.initial_guess = unittest.mock.MagicMock(return_value=hexaprotic.free_concentration)
        self.pconsol.start()   
        self.pinitgu.start()
        jac, res = self.bridge.build_matrices()
        self.pconsol.stop()   
        self.pinitgu.stop()

        np.testing.assert_allclose(res, self.data.residuals, atol=0.1)
        np.testing.assert_allclose(jac, self.data.emf_jac, atol=0.001)


class TestBridge2(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        from bridge import Bridge
        self.data, self.params = _syntheticdata.load_lmh()
        self.bridge = Bridge(self.params)

    def test_dimmensions(self):
        self.assertTupleEqual(self.bridge.jacobian.shape, (404, 3))
        self.assertTupleEqual(self.bridge.residual.shape, (404, ))

    def test_titrations(self):
        # breakpoint()
        jac, res = self.bridge.build_matrices()

        for titr, creal in zip(self.params.titrations.values(), (self.data.t1_freeconc, self.data.t2_freeconc)):
            np.testing.assert_allclose(titr.free_conc, creal, atol=5e-5)


class TestParameters(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.data, self.b = _syntheticdata.load_hexaprotic()

    def test_variables(self):
        _vars = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        self.b.update_parameters(_vars)
        v = self.b.beta
        for x, y in zip(_vars, v.logbeta):
            self.assertEqual(x, y)

    def test_constraint(self):
        self.assertListEqual(self.b.constraint, 6*[None])

    def test_jacobian_part(self):
        self.assertDictEqual(self.b.jacobian_part, {'beta': slice(0,6)})

    def test_temperature(self):
        self.assertEqual(self.b.get_temp(), 298.15)

    def test_titration(self):
        for widget in self.b.titrationwidgets:
            self.assertIn(id(widget), self.b.titrations)

    def test_datawidgets(self):
        for widget in self.b.datawidgets:
            self.assertIn(id(widget), self.b.data)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)


    unittest.main()
