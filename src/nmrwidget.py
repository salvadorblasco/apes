#!/usr/bin/python3

from PyQt5 import uic

from datawidget import DataWidget


class NmrWidget(DataWidget):
    def __init__(self, model):
        super().__init__(model)
        self.ui = uic.loadUi('../forms/nmrds.ui', self)
