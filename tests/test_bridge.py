#!/usr/bin/python3

import unittest
import sys

import numpy as np
from PyQt5 import QtWidgets

sys.path.append('../src')

import hexaprotic
import consts


class TestBridge(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        from bridge import Bridge
        self.params = load_hexaprotic()
        self.bridge = Bridge(self.params)

    def test_update_titrations(self):
        beta = self.params.beta.beta()
        self.bridge.update_titrations(beta)

        tit = tuple(self.params.titrations.values())[0]
        cc = tit.free_conc
        np.testing.assert_allclose(cc, hexaprotic.free_concentration, atol=1e-2)

<<<<<<< HEAD
    def test_fobj(self):
        variables = 10**np.array([10.0, 18.0, 24.0, 28.0, 31.0, 33.0])
        f = self.bridge.generate_freeconcs()
        _ = f(variables)
        ffobj = self.bridge.generate_fobj()
        fobj = ffobj(variables)
        # np.testing.assert_allclose(jvals, jreal, atol=1e-8)

    def test_jacobian(self):
        variables = 10**np.array([10.0, 18.0, 24.0, 28.0, 31.0, 33.0])
        f = self.bridge.generate_freeconcs()
        _ = f(variables)
        fjac = self.bridge.generate_jacobian()
        jvals = fjac(variables)
        jreal = consts.NERNST*hexaprotic.dlogc_dlogbeta[:,1,:6]/variables[None,:]
        np.testing.assert_allclose(jvals, jreal, atol=1e-8)
=======
    def test_matrices(self):
        with self.subTest(t="test dimmensions"):
            self.assertTupleEqual(self.bridge.jacobian.shape, (89, 6))
            self.assertTupleEqual(self.bridge.residual.shape, (89, ))

        variables = np.array([10.0, 18.0, 24.0, 28.0, 31.0, 33.0])

        breakpoint()
        with self.subTest(t="function call"):
            jacobian, residual = self.bridge.build_matrices(variables)

        with self.subTest(t="test residual call"):
            np.testing.assert_allclose(residual, np.zeros_like(residual), atol=0.1)

        with self.subTest(t="test jacobian"):
            jreal = consts.NERNST*hexaprotic.dlogc_dlogbeta[:,1,:6]/variables[None,:]
            np.testing.assert_allclose(jacobian, jreal, atol=1e-8)
>>>>>>> 760caa2d7e14c2e1c78626ea41c9a500017975b9


class TestParameters(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.b = load_hexaprotic()

    def test_get_temp(self):
        self.assertEqual(self.b.get_temp(), 298.15)

    def test_variables(self):
        _vars = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        self.b.update_parameters(_vars)
<<<<<<< HEAD
        v = self.b.beta
        for x, y in zip(_vars, v.logbeta):
=======
        v = self.b.beta.logbeta
        for x, y in zip(_vars, v):
>>>>>>> 760caa2d7e14c2e1c78626ea41c9a500017975b9
            self.assertEqual(x, y)

    def test_constraint(self):
        self.assertListEqual(self.b.constraint, 6*[None])

    def test_jacobian_part(self):
        self.assertDictEqual(self.b.jacobian_part, {'beta': slice(0,6)})


def load_hexaprotic():
    from modelwidget import ModelWidget, ModelData
    model = ModelWidget()
    mdata = ModelData(n_equils=7, n_species=2)
    mdata.stoich = hexaprotic.stoich
    mdata.const = hexaprotic.logbeta
    mdata.const_flags = 6*[consts.RF_REFINE] + [consts.RF_CONSTANT]
    model.append(mdata)
    model.setCurrentModel(-1)

    from otherwidgets import TitrationBaseWidget
    titr = TitrationBaseWidget(model)
    titr.name = 'Titration'
    titr.initial_amount = hexaprotic.init
    titr.buret = hexaprotic.buret
    titr.starting_volume = hexaprotic.v0
    titr.titre = hexaprotic.titre.tolist()

    from emfwidget import EmfWidget
    emfw = EmfWidget(model)
    emfw.emf0 = hexaprotic.emf0
    emfw.emf = hexaprotic.emf
    emfw._titrationid = id(titr)

    from bridge import Parameters
    b = Parameters(model, [titr], [emfw])
    return b


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)

    unittest.main()
