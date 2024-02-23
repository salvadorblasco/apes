#!/usr/bin/python3

from datawidget import DataWidget
import ui_nmrds

class NmrWidget(DataWidget):
    def __init__(self, model):
        super().__init__(model)
        self.ui = ui_nmrds.Ui_NmrWidget()
        self.ui.setupUi(self)
