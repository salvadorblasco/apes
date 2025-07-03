"Classes related to dialogs."

from PyQt5 import QtWidgets
from PyQt5 import uic


class AboutDialog(QtWidgets.QDialog):
    '''Display 'about' dialog.'''

    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        self.ui = uic.loadUi('../forms/about.ui', self)


class Constants(QtWidgets.QDialog):
    """Dialog for calculate stepwise constants."""
    def __init__(self, beta, error_beta, labels):
        super().__init__()
        self.ui = uic.loadUi('../forms/constants.ui', self)

        combos = (self.ui.cb_reagent1, self.ui.cb_reagent2,
                  self.ui.cb_reagent3, self.ui.cb_reagent4)
        self.__ebeta = [0.0] + [b for b in beta]
        self.__eerr = [0.0] + [b for b in error_beta]

        for combo in combos:
            combo.addItem('')
            _lbl = ['log₁₀β({})'.format(i) for i in labels]
            combo.addItems(_lbl)
            combo.currentIndexChanged.connect(self._recalculate)

    def _recalculate(self):
        ucombos = (self.ui.cb_reagent1, self.ui.cb_reagent2)
        dcombos = (self.ui.cb_reagent3, self.ui.cb_reagent4)
        uind = [c.currentIndex() for c in ucombos]
        dind = [c.currentIndex() for c in dcombos]
        kval = sum(self.__ebeta[i] for i in uind) - \
               sum(self.__ebeta[i] for i in dind)
        evll = (sum(self.__eerr[i]**2 for i in uind) + \
               sum(self.__eerr[i]**2 for i in dind))**0.5

        self.ui.lbl_result.setText("logK = {:.3f}±{:.3f}".format(kval, evll))


class OptionsDialog(QtWidgets.QDialog):
    '''Dialog for setting the options.'''

    def __init__(self):
        super().__init__()
        self.ui = uic.loadUi('../forms/optionsdialog.ui', self)

    @property
    def fitparams_coarse(self):
        return float(self.ui.le_coarsetolc.text()), \
               float(self.ui.le_coarsetolv.text()), \
               self.ui.sb_coarsemaxits.value()

    @property
    def fitparams_fine(self):
        return float(self.ui.le_finetolc.text()), \
               float(self.ui.le_finetolv.text()), \
               self.ui.sb_finemaxits.value()

    @property
    def ignore_lower_than(self):
        return self.ui.dsb_ignorelwr.value()

    @property
    def unlabel_lower_than(self):
        return self.ui.dsb_unlabellwr.value()

    @property
    def verbosity(self):
        return self.ui.cb_verbosity.currentIndex()


class PropertiesDialog(QtWidgets.QDialog):
    def __init__(self, project):
        super().__init__()
        self.ui = uic.loadUi('../forms/docproperties.ui', self)


class ModelInputDialog(QtWidgets.QDialog):
    def __init__(self, project):
        super().__init__()
        self.ui = uic.loadUi('../forms/model_input.ui', self)
