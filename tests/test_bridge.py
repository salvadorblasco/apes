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
        self.app = QtWidgets.QApplication(sys.argv)

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

        from bridge import Bridge
        self.b = Bridge(model, [titr], [emfw])

    def test_dimmensions(self):
        self.assertTupleEqual(self.b.jacobian.shape, (89, 6))
        self.assertTupleEqual(self.b.residual.shape, (89, ))
        self.assertEqual(len(self.b.variables), 6)

    def test_variables(self):
        _vars = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        self.b.update_parameters(_vars)
        v = self.b.parameter['beta']
        for x, y in zip(_vars, v):
            self.assertEqual(x, y)

    def test_freeconcs(self):
        variables = 10**np.array([10.0, 18.0, 24.0, 28.0, 31.0, 33.0])
        f = self.b.generate_freeconcs()
        c = f(variables)
        for cc in c.values():
            np.testing.assert_allclose(cc, hexaprotic.free_concentration, atol=1e-2)

    def test_jacobian(self):
        variables = 10**np.array([10.0, 18.0, 24.0, 28.0, 31.0, 33.0])
        f = self.b.generate_freeconcs()
        _ = f(variables)
        fjac = self.b.generate_jacobian()
        ...


if __name__ == '__main__':
    unittest.main()
