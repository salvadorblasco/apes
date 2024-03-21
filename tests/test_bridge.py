#!/usr/bin/python3

import unittest

import sys

from PyQt5 import QtWidgets

sys.path.append('../src')


class TestBridge(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.app = QtWidgets.QApplication(sys.argv)

    def test_emf_only(self):
        import hexaprotic
        import consts

        from modelwidget import ModelWidget
        model = ModelWidget()
        model.addEquilibrium(-1, (1,1), 10.00, 0.0, consts.RF_REFINE)
        model.addEquilibrium(-1, (1,2), 18.00, 0.0, consts.RF_REFINE)
        model.addEquilibrium(-1, (1,3), 24.00, 0.0, consts.RF_REFINE)
        model.addEquilibrium(-1, (1,4), 28.00, 0.0, consts.RF_REFINE)
        model.addEquilibrium(-1, (1,5), 31.00, 0.0, consts.RF_REFINE)
        model.addEquilibrium(-1, (1,6), 33.00, 0.0, consts.RF_REFINE)
        model.addEquilibrium(-1, (0,-1), -13.77, 0.0, consts.RF_CONSTANT)
        print(model.stoich)
        print(model.beta)
        print(model.beta_flags)
        sys.exit(0)

        from otherwidgets import TitrationBaseWidget
        titr = TitrationBaseWidget(model)
        titr.name = 'Titration'
        titr.init = hexaprotic.init
        titr.buret = hexaprotic.buret
        titr.starting_volume = hexaprotic.v0
        titr.final_volume = 12.00
        titr.n_points = len(hexaprotic.titre)

        from emfwidget import EmfWidget
        emfw = EmfWidget(model)
        emfw.emf0 = hexaprotic.emf0
        emfw.emf = hexaprotic.emf

        from bridge import Bridge
        b = Bridge(model, [titr], [emfw])

        print(b.jacobian.shape)

        with self.subTest("generate free concs"):
            f = b.generate_freeconcs()


if __name__ == '__main__':
    unittest.main()
