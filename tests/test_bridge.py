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

    def test_dimmensions(self):
        self.assertTupleEqual(self.bridge.jacobian.shape, (89, 6))
        self.assertTupleEqual(self.bridge.residual.shape, (89, ))

    def test_matrices(self):
        variables = np.array([10.0, 18.0, 24.0, 28.0, 31.0, 33.0])
        jac, res = self.bridge.build_matrices(variables)

        for cc in self.params.titrations.values():
            np.testing.assert_allclose(cc.free_conc, hexaprotic.free_concentration, atol=1e-2)

        np.testing.assert_allclose(res, np.zeros_like(res), atol=0.8)

        jreal = consts.NERNST*hexaprotic.dlogc_dlogbeta[:,1,:6]/(10**variables[None,:])
        np.testing.assert_allclose(jac, jreal, atol=1e-8)


class TestParameters(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.b = load_hexaprotic()

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
