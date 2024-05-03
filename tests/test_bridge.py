#!/usr/bin/python3

import unittest
import unittest.mock
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

    def test_titrations(self):
        # breakpoint()
        jac, res = self.bridge.build_matrices()

        for cc in self.params.titrations.values():
            np.testing.assert_allclose(cc.free_conc, hexaprotic.free_concentration, atol=4e-3)
            np.testing.assert_allclose(cc.amatrix, hexaprotic.matrix_a, atol=0.2)

    def test_matrices(self):
        import libeq.consol
        libeq.consol.consol = unittest.mock.MagicMock(return_value=hexaprotic.free_concentration)
        libeq.consol.initial_guess = unittest.mock.MagicMock(return_value=hexaprotic.free_concentration)
    
        jac, res = self.bridge.build_matrices()

        np.testing.assert_allclose(res, hexaprotic.residuals, atol=0.1)
        np.testing.assert_allclose(jac, hexaprotic.emf_jac, atol=0.001)


class TestBridge2(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        from bridge import Bridge
        self.params, self.data = load_lmh()
        self.bridge = Bridge(self.params)

    def test_dimmensions(self):
        self.assertTupleEqual(self.bridge.jacobian.shape, (404, 3))
        self.assertTupleEqual(self.bridge.residual.shape, (404, ))


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


def load_lmh():
    import data_lmh

    from modelwidget import ModelWidget, ModelData
    model = ModelWidget()
    model.clear(n_equils=data_lmh.N_EQUILS, n_species=data_lmh.N_COMPON)
    mdata = ModelData(n_equils=data_lmh.N_EQUILS, n_species=data_lmh.N_COMPON)
    mdata.stoich = data_lmh.stoich
    mdata.const = data_lmh.logbeta
    mdata.const_flags = 6*[consts.RF_CONSTANT] + 3*[consts.RF_REFINE] + [consts.RF_CONSTANT]
    model.append(mdata)
    model.setCurrentModel(-1)

    from otherwidgets import TitrationBaseWidget
    titr1 = TitrationBaseWidget(model)
    titr1.name = 'Titration 1'
    titr1.initial_amount = data_lmh.t1_init
    titr1.buret = data_lmh.t1_buret
    titr1.starting_volume = data_lmh.t1_startingvol
    titr1.final_volume = data_lmh.t1_endvol

    titr2 = TitrationBaseWidget(model)
    titr2.name = 'Titration 2'
    titr2.initial_amount = data_lmh.t2_init
    titr2.buret = data_lmh.t2_buret
    titr2.starting_volume = data_lmh.t2_startingvol
    titr2.final_volume = data_lmh.t2_endvol

    from emfwidget import EmfWidget
    emfw1 = EmfWidget(model)
    emfw1.emf0 = data_lmh.t1_emf0
    emfw1.slope = (1.0, 1.0)
    emfw1.slope_flags = (0, 0)
    emfw1.nelectrons = (1, 1)
    emfw1.emf0_error = (0.01,0.01)
    emfw1.active_species = (1, 2)
    emfw1.emf0_flags = (0,0)
    emfw1.titration = titr1
    emfw1.titre = titr1.titre
    emfw1.emf = data_lmh.t1_emf

    emfw2 = EmfWidget(model)
    emfw2.emf0 = data_lmh.t2_emf0
    emfw2.slope = (1.0, 1.0)
    emfw2.slope_flags = (0, 0)
    emfw2.nelectrons = (1, 1)
    emfw2.emf0_error = (0.01,0.01)
    emfw2.active_species = (1, 2)
    emfw2.emf0_flags = (0,0)
    emfw2.titration = titr2
    emfw2.emf = data_lmh.t2_emf.T

    from bridge import Parameters
    b = Parameters(model, [titr1, titr2], [emfw1, emfw2])
    return b, data_lmh


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
    emfw.titration = titr

    from bridge import Parameters
    b = Parameters(model, [titr], [emfw])
    return b


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)

    unittest.main()
