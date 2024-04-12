#!/usr/bin/python3

import unittest

import sys
import numpy as np
from PyQt5 import QtWidgets

import hexaprotic

sys.path.append('../src')

import bridge
import consts
import libfit


class TestLevenberg(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.params = load_hexaprotic()
        self.bridge = bridge.Bridge(self.params)

    def test_levenberg_marquardt(self):
        initvars = np.fromiter(self.params.initial_values(), dtype=float) 
        # initvars[...] = initvars + np.random.rand(len(initvars))
        initvars[0] = initvars[0] +0.1
        weights = self.bridge.weights()
        
        x, info = libfit.levenberg_marquardt(initvars, self.bridge.build_matrices, weights)

        self.params.update_parameters(x)
        self.params.dump_to_widgets()


# class TestSimplex(unittest.TestCase):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
# 
#     def setUp(self):
#         self.logB = np.array([10.13, 19.53, 27.80, 34.82, -13.73])
#         self.B = 10**self.logB
#         self.P = np.array([[1, 1], [1, 2], [1, 3], [1, 4], [0, -1]])
#         with np.load('pytrenal.npz') as d:
#             self.real_C = d['C']
#             self.real_T = d['T']
#         self.E, self.S = self.P.shape
# 
#     def test_simplex(self):
#         pass


# class TestAuxFunctions(unittest.TestCase):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
# 
#     def test_trivial_capping(self):
#         from libfit import trivial_capping
#         a = np.random.rand(10)
#         b = np.random.rand(10)
#         res = trivial_capping(a, b)
#         np.testing.assert_allclose(res, a+b)
# 
#     def test_ratio_capping(self):
#         from libfit import max_ratio_capping
#         x = np.ones(5)
#         dx = np.linspace(0.1, 0.5, 5)
#         ret = max_ratio_capping(x, dx, 0.25)
#         tru = np.array([1.1, 1.2, 1.25, 1.25, 1.25])
#         np.testing.assert_allclose(ret, tru)


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
