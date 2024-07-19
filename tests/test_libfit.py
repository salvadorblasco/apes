#!/usr/bin/python3

import unittest

import sys
import numpy as np
import numpy.testing as npt
from PyQt5 import QtWidgets

import _syntheticdata

sys.path.append('../src')

import bridge
import consts
import libfit


class TestLevenberg(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.data, self.params = _syntheticdata.load_hexaprotic()
        self.bridge = bridge.Bridge(self.params)

    def test_levenberg(self):
        for noise in (0.0, 0.1, 0.2, 0.5, 0.75, 1.0):
            with self.subTest(noise=noise):
                self._marquardt_call(noise)
                values = np.fromiter(self.params.initial_values(), dtype=float) / consts.LOGK
                npt.assert_allclose(values, self.data.logbeta[:6], atol=1e-2)

    def _marquardt_call(self, noise_level: float):
        noise = noise_level * (np.random.rand(6) - 0.5)
        self.bridge.step_values(noise)
        self.bridge.accept_values()
        info = libfit.levenberg_marquardt(self.bridge)
        self.params.dump_to_widgets()


class TestLevenberg2(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.data, self.params = _syntheticdata.load_lmh()
        self.bridge = bridge.Bridge(self.params)

    def test_levenberg(self):
        for noise in (0.0, 0.1, 0.2, 0.5, 0.75, 1.0):
            with self.subTest(noise=noise):
                self._marquardt_call(noise)
                values = np.fromiter(self.params.initial_values(), dtype=float) / consts.LOGK
                npt.assert_allclose(values, self.data.logbeta[6:9], atol=1e-3)

    def _marquardt_call(self, noise_level: float):
        noise = noise_level * (np.random.rand(3) - 0.5)
        self.bridge.step_values(noise)
        self.bridge.accept_values()
        info = libfit.levenberg_marquardt(self.bridge)
        self.params.dump_to_widgets()


# def load_hexaprotic():
#     import hexaprotic
# 
#     from modelwidget import ModelWidget, ModelData
#     model = ModelWidget()
#     mdata = ModelData(n_equils=7, n_species=2)
#     mdata.stoich = hexaprotic.stoich
#     mdata.const = hexaprotic.logbeta
#     mdata.const_flags = 6*[consts.RF_REFINE] + [consts.RF_CONSTANT]
#     model.append(mdata)
#     model.setCurrentModel(-1)
# 
#     from otherwidgets import TitrationBaseWidget
#     titr = TitrationBaseWidget(model)
#     titr.name = 'Titration'
#     titr.initial_amount = hexaprotic.init
#     titr.buret = hexaprotic.buret
#     titr.starting_volume = hexaprotic.v0
#     titr.titre = hexaprotic.titre.tolist()
# 
#     from emfwidget import EmfWidget
#     emfw = EmfWidget(model)
#     emfw.emf0 = hexaprotic.emf0
#     emfw.emf = hexaprotic.emf
#     emfw.titration = titr
# 
#     from bridge import Parameters
#     b = Parameters(model, [titr], [emfw])
#     return hexaprotic, b


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)

    unittest.main()
