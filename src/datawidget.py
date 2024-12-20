from PyQt5 import QtWidgets as QtGui
from PyQt5 import QtCore

import numpy as np
from numpy.typing import NDArray

import consts
import libqt
import libeq
from modelwidget import ModelWidget
from otherwidgets import TitrationBaseWidget


FloatArray = NDArray[float]


class DataWidget(QtGui.QWidget):
    '''Parent class for some data widgets.

    Manage tasks and attributes common to :class:`emfwidget.EmfWidget`,
    :class:`nmrwidget.NMRWidget`, :class:`calorwidget.CalorWidget` and
    :class:`specwidget.SpecWidget`.

    Signals:
        labelsChanged: when the user changes manually a label
    '''
    # labelsChanged = QtCore.pyqtSignal(int, str)
    titrationChanged = QtCore.pyqtSignal(object, str)

    def __init__(self, model: ModelWidget, parent=None):
        super().__init__(parent)
        self._name: str = 'No name'
        self._collabels: int = 0         # column number for table_trit
        self._temperature: float = 298.15  # Kelvin
        self._model: ModelWidget = model
        self._use: bool = True
        self._tabdata = None
        # self._tabtitr = None
        self._datafit = None    # slot for data from fitting
        # self._freeconcentration = None
        self._setweight: float = 1.0    # the weight of the whole dataset
        # self.__totalconc = None
        self._buildmenu()
        self._titrationid = None
        self._titration = None
    
    def __titration_changed(self, txt):
        # print("signal emitted", self, index)
        self.titrationChanged.emit(self, txt)

    # def calc_freec(self, altbeta=None, altstoich=None):
    #     beta = 10**np.array(self._model.beta if altbeta is None else altbeta)
    #     stoichiometry = np.array(self._model.stoich if altstoich is None else altstoich)

    #     analyticalc = np.array(self.analyticalc())
    #     if self._freeconcentration is not None:
    #         kw = {'initial_values': self._freeconcentration}
    #     else:
    #         iniguess = libeq.consol.initial_guess(beta, stoichiometry, analyticalc)
    #         kw = {'initial_values': iniguess}
    #     self._freeconcentration = libeq.consol.consol(
    #         beta, stoichiometry, analyticalc, **kw)

    def feed_data_table(self, data_by_columns):
        if self._tabdata is None:
            raise RuntimeError
        for col, data in enumerate(data_by_columns):
            libqt.fill_column(self._tabdata, col, data)

    def model(self):
        return self._model

    def reshape_data_table(self, rows, cols):
        if self._tabdata is None:
            raise RuntimeError
        self._tabdata.setColumnCount(cols)
        self._tabdata.setRowCount(rows)

    # def set_free_concentration(self, c):
    #     self._freeconcentration = c

    def populate_cb_titration(self, args):
        with libqt.signals_blocked(self.ui.cb_titration) as combo:
            combo.clear()
            combo.addItems(args)

    @property
    def initial_amount(self):
        'Initial amount for every component (in mmol)'
        stream = libqt.iter_column_text(self._tabtitr, col=1)
        return tuple(float(i) for i in stream)

    # @initial_amount.setter
    # def initial_amount(self, amount):
    #     libqt.fill_column(self._tabtitr, col=1, data=amount)

    # @property
    # def labels(self):
    #     'list of labels for every component'
    #     # return tuple(libqt.iter_column_text(self._tabtitr, col=0))
    #     return self._model.labels

    # @labels.setter
    # def labels(self, labels):
    #     libqt.fill_column(self._tabtitr, col=self._collabels, data=labels,
    #                       flags=QtCore.Qt.ItemIsSelectable)

    # @property
    # def amount_flags(self):
    #     libqt.iter_column_comboindex(self._tabtitr, col=self._colflags)

    # @amount_flags.setter
    # def amount_flags(self, flags):
    #     # TODO check input
    #     # TODO reshape table if needed
    #     libqt.fill_column_comboindex(self._tabtitr, flags,
    #                                  consts.REFINE_LABELS, col=self._colflags)

    # @property
    # def buret(self):
    #     'Contents of the buret (in mmol/mL)'
    #     stream = libqt.iter_column_text(self._tabtitr, col=2)
    #     return tuple(float(i) for i in stream)

    # @buret.setter
    # def buret(self, bdata):
    #     libqt.fill_column(self._tabtitr, col=2, data=bdata)

    @property
    def experimental_data(self):
        return libqt.tab2array(self._tabexpd, rows=self._exprow0,
                               cols=self._expcol0, masked=True)

    @experimental_data.setter
    def experimental_data(self, expdata):
        libqt.array2tab(self._tabexpd, expdata, row0=min(self._exprow0),
                        col0=min(self._expcol0))

    @property
    def free_concentration(self) -> FloatArray:
        return self._titration.free_concentration

    # @free_concentration.setter
    # def free_concentration(self, c: FloatArray):
    #     self._freeconcentration = c

    @property
    def name(self):
        'The name of the dataset.'
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def temperature(self):
        'The temperature, in Kelvin'
        return self._temperature

    @temperature.setter
    def temperature(self, temp):
        if not isinstance(temp, float):
            raise TypeError('float expected')
        if temp <= 0:
            raise ValueError('temperature must be positive')
        self._temperature = temp

    @property
    def titration(self):
        return self._titration

    @titration.setter
    def titration(self, titrationwidget):
        if type(titrationwidget) == TitrationBaseWidget:
            self._titration = titrationwidget
        else:
            self._titration = self.parent.find_titration_byname(titrationwidget)
        self._select_titration_dropdown()

    @property
    def use(self):
        "A flag indicating whether this dataset will be used or not"
        return self._use

    @use.setter
    def use(self, u: bool):
        self._use = bool(u)

    def _connectMenu(self, tab):
        assert isinstance(tab, QtGui.QTableWidget)
        self._tabdata = tab
        tab.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        tab.customContextMenuRequested.connect(self._pmtd)

    def _buildmenu(self):
        p = QtGui.QMenu(parent=self)
        # a = p.addAction("new point(s)")
        # a.triggered.connect(self._newpoint)
        a = p.addAction("use point(s)")
        a.triggered.connect(self._tdpmuse)
        a = p.addAction("ignore point(s)")
        a.triggered.connect(self._tdpmign)
        # a = p.addAction("delete point(s)")
        # a.triggered.connect(self._tdpmdel)
        self._popupm = p

    def _set_default_titration(self):
        n = self._model.number_components
        self._tabtitr.setRowCount(n)
        self.amount_flags = n*[consts.RF_CONSTANT]
        self.buret = n*[0.0]
        self.initial_amount = n*[0.001]
        self.labels = self._model.labels

    def _tdpmign(self):
        "Slot for popup menu option 'ignore points'"
        self.__ignorusecommon(True)

    def _tdpmuse(self):
        "Slot for popup menu option 'use points'"
        self.__ignorusecommon(False)

    def __ignorusecommon(self, to_cross):
        assert isinstance(self._tabdata, QtGui.QTableWidget)
        for row in {x.row() for x in self._tabdata.selectedIndexes()}:
            libqt.cross(libqt.iter_row(self._tabdata, row), to_cross)

    def _newpoint(self):
        raise NotImplementedError
        n, ok = QtGui.QInputDialog.getInt(self, 'add new titration points',
                                          'How many?', value=1, min=1,
                                          max=9999)
        if not ok:
            return
        # TODO complete

    def _tdpmdel(self):
        raise NotImplementedError

    def _pmtd(self, point):
        self._popupm.popup(self._tabdata.viewport().mapToGlobal(point))
