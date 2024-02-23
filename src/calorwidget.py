"""
Widget for calorimetry.

.. module:: calorwidget.py
.. moduleauthor:: Salvador Blasco <salvador.blasco@gmail.com>

"""

from PyQt5 import QtWidgets as QtGui
from PyQt5 import QtCore

import consts
import datawidget
import libqt
import ui_calords


class CalorWidget(datawidget.DataWidget):
    """Abstraction and control calorimetry titration data.

    Signals:
        - labelsChanged: when the user changes a label
    Slots:
        - updateLabel(position, new_label)
    """
    def __init__(self, model):
        """Initiate widget.

        Parameters:
            model (:class:`modelwidget.ModelWidget`): The model reference
        """
        super().__init__(model)
        self.ui = ui_calords.Ui_CalorWidget()
        self.ui.setupUi(self)
        self._connectMenu(self.ui.table_data)

    @QtCore.pyqtSlot(int, str)
    def updateLabel(self, position: int, new_label: str):
        "Slot for when a label changes"
        self.ui.table_titration.setItem(position, 0,
                                        QtGui.QTableWidgetItem(new_label))

    @property
    def amount_flags(self):
        'Flags for the starting amounts'
        return tuple(libqt.iter_column_comboindex(self.ui.table_titration,
                                                  col=3))

    @amount_flags.setter
    def amount_flags(self, flags):
        libqt.fill_column_comboindex(self.ui.table_titration, flags,
                                     consts.REFINE_LABELS, col=3)

    @property
    def buret(self):
        'Contents of the buret (in mmol/mL)'
        return list(libqt.iter_column_text(self.ui.table_params, col=2))

    @buret.setter
    def buret(self, buret):
        libqt.fill_column(self.ui.table_titration, col=2, data=buret)

    @property
    def initial_amount(self):
        'Initial amount for every component (in mmol)'
        return tuple(libqt.iter_column_text(self.ui.table_titration, col=1))

    @initial_amount.setter
    def initial_amount(self, amount):
        libqt.fill_column(self.ui.table_titration, col=1, data=amount)

    @property
    def labels(self):
        'list of labels for every component'
        return list(libqt.iter_column_text(self.ui.table_titration, col=0))

    @labels.setter
    def labels(self, labels):
        libqt.fill_column(self._tabtitr, col=self._collabels, data=labels)

    @property
    def heat(self):
        'Experimental heat values'
        return tuple(libqt.iter_column_text(self.ui.table_data, col=1))

    @heat.setter
    def heat(self, q):
        # TODO check input
        self.__recheck_table_size(len(q))
        libqt.fill_column(self.ui.table_data, col=1, data=q)

    @property
    def starting_volume(self):
        'Starting volume, in mL'
        return self.ui.dsb_V0.value()

    @starting_volume.setter
    def starting_volume(self, starting_volume: float):
        self.ui.dsb_V0.setValue(starting_volume)

    @property
    def titre(self):
        'Volume of titre used for every point.'
        return tuple(libqt.iter_column_text(self.ui.table_data, col=0))

    @titre.setter
    def titre(self, volume: float):
        # TODO check input
        self.__recheck_table_size(len(volume))
        libqt.fill_column(self.ui.table_data, col=0, data=volume)

    @property
    def volume_error(self):
        "The error in measure of volume."
        return self.ui.dsb_errV.value()

    @volume_error.setter
    def volume_error(self, verror: float):
        self.ui.dsb_errV.setValue(verror)

    def __recheck_table_size(self, ndat):
        if ndat != self.ui.table_data.rowCount():
            self.ui.table_data.setRowCount(ndat)

    def _tflgch(self):
        raise NotImplementedError
