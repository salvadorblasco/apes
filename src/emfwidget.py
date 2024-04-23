""" ######################
    # file: emfwidget.py #
    ###################### """

import itertools

import numpy as np

from PyQt5 import QtCore, QtWidgets

import consts
import datawidget
import libaux
import libemf
import libmath
import libqt
import ui_emfds


class EmfWidget(datawidget.DataWidget):
    """Widget and encapsulation for data from potentiometric titrations.

    The data is encapsulated in three tables: (1) table_titration,
    (2) table_params and (3) table_data.
    """

    def __init__(self, model, parent):
        """Initiate widget."""
        super().__init__(model)
        self.parent = parent
        self.ui = ui_emfds.Ui_EmfDSWidget()
        self.ui.setupUi(self)

        self._tabtitr = None

        self._connectMenu(self.ui.table_data)

        self.__dangerous_allowed = False
        self.emf0_flags = (0,)      # default parameters
        self.slope_flags = (0,)      # default parameters
        self.nelectrons = (1,)
        self.active_species = (1,)
        # self.amount_flags = (0, 0)
        # self.buret_flags = (0, 0)

        for n, l in enumerate(model.labels):
            self.updateLabel(n, l)

        model.labelsChanged.connect(self.updateLabel)
        model.componentAdded.connect(self.componentAdded)
        model.componentDeleted.connect(self.componentDeleted)

        popup = self._popupm
        # action = popup.addAction('New electrode')
        # action.triggered.connect(self._newlectrode)
        # action = popup.addAction('Remove electrode')
        # action.triggered.connect(self._dellectrode)

        self.ui.cb_titration.currentTextChanged.connect(super()._DataWidget__titration_changed)

    def add_component(self, label, position=-1):
        table = self.ui.table_titration
        ncomps = self.model.number_components

        _buret = list(self.buret) + [0.1]
        _initial_amount = list(self.initial_amount) + [0.001]
        _amount_flags = list(self.amount_flags) + [consts.RF_CONSTANT]
        table.setRowCount(ncomps)
        # self.updateLabel(position, label)
        self.buret = _buret
        self.initial_amount = _initial_amount 
        self.labels = self.model.labels
        self.amount_flags = _amount_flags

    def analyticalc(self):
        'Analytical concentrations'
        return tuple(libaux.build_T_titr2(
            self.initial_amount,
            self.buret,
            self.starting_volume,
            self.titre))

    def allow_dangerous(self, allow=False):
        "Activate the refinement of dangerous parameters."
        def _touch(what):
            for dropdown in what:
                dropdown.setCurrentIndex(0)
                dropdown.setEnabled(allow)

        self.__dangerous_allowed = allow
        _touch(libqt.iter_column_widget(self.ui.table_titration, col=3))
        _touch(libqt.iter_row_widget(self.ui.table_params, row=2))

    def weight(self, scheme='auto'):
        titre = tuple(self.titre)

        def _w(x, ex):
            _tx_ = tuple(x)
            return libmath.weighting_slope(titre, _tx_,
                                           self.volume_error, ex)

        retv = [tuple(_w(_emf_, _err_))
                for _emf_, _err_ in zip(self._iter_emf(), self.emf0_error)]

        # retv = []
        # for _emf_, _err_ in zip(self._iter_emf(), self.emf0_error):
        #     _ = tuple(_w(_emf_, _err_))
        #     retv.append(_)
        #     # retv.append(tuple(_w(_emf_, _err_)))
        return tuple(zip(*retv))

    @QtCore.pyqtSlot(int, str)
    def updateLabel(self, position, new_label):
        "Slot for when a label changes"
        # item = QtWidgets.QTableWidgetItem(new_label)
        # item.setFlags(QtCore.Qt.ItemIsSelectable)
        # self.ui.table_titration.setItem(position, 0, item)
        for col in range(self.ui.table_params.columnCount()):
            combo = self.ui.table_params.cellWidget(4, col)
            combo.setItemText(position, new_label)

    @QtCore.pyqtSlot(int, str)
    def componentAdded(self, which, label):
        for col in range(self.ui.table_params.columnCount()):
            combo = self.ui.table_params.cellWidget(4, col)
            combo.insertItem(which, label)

    @QtCore.pyqtSlot(int)
    def componentDeleted(self, which):
        for col in range(self.ui.table_params.columnCount()):
            combo = self.ui.table_params.cellWidget(4, col)
            combo.removeItem(which)

    @property
    def active_species(self):
        'Indices of the active species'
        return tuple(libqt.iter_row_comboindex(self.ui.table_params, row=4))

    @active_species.setter
    def active_species(self, activ):
        # TODO check input
        # TODO reshape table if needed
        libqt.fill_row_comboindex(self.ui.table_params, self.__num2it(activ),
                                  self._model.labels, row=4)

    # @property
    # def amount_flags(self):
    #     'Flags for the starting amounts'
    #     return tuple(libqt.iter_column_comboindex(self.ui.table_titration,
    #                                               col=3))

    # @amount_flags.setter
    # def amount_flags(self, flags):
    #     libqt.fill_column_comboindex(self.ui.table_titration, flags,
    #                                  consts.REFINE_LABELS, col=3)

    @property
    def electrodes(self):
        return {'E0': np.array(self.emf0),
                'E0flags': self.emf0_flags,
                'hindex': self.active_species,
                'fRTnF': np.array(self.fRTnF)}

    # def titration(self):
    #     uf = np.fromiter(self.mask, dtype=np.bool)
    #     return {'V': np.array(np.fromiter(self.titre, dtype=np.float)),
    #             'emf': np.ma.array(self.emf, mask=uf),
    #             'V0': self.starting_volume,
    #             'T0': self.initial_amount,
    #             'buret': self.buret,
    #             'weights': self.weight(),
    #             'Tflags': self.amount_flags,
    #             'error_V': self.volume_error}

    @property
    def emf(self):
        'Measured values of emf.'
        return libqt.tab2array(
            self.ui.table_data,
            rows=range(0, self.ui.table_data.rowCount()),
            cols=range(1, self.ui.table_data.columnCount()), masked=False)

    @emf.setter
    def emf(self, emf):
        # TODO check input
        # TODO reshape table if needed
        if isinstance(emf[0], float):
            ncols = 2
            nrows = len(emf)
            self.__resize_emf_table(nrows, ncols)
            libqt.fill_column(self.ui.table_data, 1, emf,
                              formatting="{:.3f}")
        else:
            nrows = len(emf[0]) + 1
            ncols = len(emf)
            self.__resize_emf_table(nrows, ncols)
            libqt.array2tab(self.ui.table_data, emf, col0=1,
                            formatting="{:.3f}")

    @property
    def emf_fit(self):
        if self._freeconcentration is None:
            return None

        h = self._freeconcentration[:, self.active_species]
        emf0 = np.array(self.emf0)
        nernst = np.atleast_1d(np.array(self.fRTnF))
        return emf0[None, :] + nernst[None, :]*np.log(h)

    def __resize_emf_table(self, rows, cols):
        self.ui.table_data.setRowCount(rows)
        self.ui.table_data.setColumnCount(cols)

    @property
    def emf0_error(self):
        'Error(s) of emf0'
        stream = libqt.iter_row_text(self.ui.table_params, row=2)
        return tuple(float(i) for i in stream)

    @emf0_error.setter
    def emf0_error(self, erremf0):
        # TODO check input
        # TODO reshape table if needed
        libqt.fill_row(self.ui.table_params, row=2, data=self.__num2it(erremf0))

    @property
    def emf0(self):
        'Value(s) of emf0'
        stream = libqt.iter_row_text(self.ui.table_params, row=0)
        return tuple(float(i) for i in stream)

    @emf0.setter
    def emf0(self, emf0):
        # TODO check input
        # if self.nelectrodes != len(emf0):
        #     SELF._NELECTRS_CHANGED(SELF, LEN(EMF0))
        libqt.fill_row(self.ui.table_params, row=0, data=self.__num2it(emf0))

    @property
    def emf0_flags(self):
        'Refinement flags for emf0 values'
        return tuple(libqt.iter_row_comboindex(self.ui.table_params, row=1))

    @emf0_flags.setter
    def emf0_flags(self, flags):
        libqt.fill_row_comboindex(self.ui.table_params, flags, consts.REFINE_LABELS, row=1)

    @property
    def slope(self):
        stream = libqt.iter_row_text(self.ui.table_params, row=5)
        return tuple(float(i) for i in stream)

    @property
    def slope_flags(self):
        'Refinement flags for slope values'
        return tuple(libqt.iter_row_comboindex(self.ui.table_params, row=6))

    @slope_flags.setter
    def slope_flags(self, flags):
        libqt.fill_row_comboindex(self.ui.table_params, flags, consts.REFINE_LABELS, row=6)

    # deprecate - use slope
    @property
    def fRTnF(self):
        "Nernst slope factor."
        r = consts.RoverF*self.temperature
        assert isinstance(r, float)
        # if self.nelectrodes == 1:
        #     return r*self.nelectrons
        # else:
        #     return tuple(n*r for n in self.nelectrons)
        return tuple(n*r for n in self.nelectrons)

    @property
    def gradient(self):
        weights = np.array(self.weight())
        residuals = self.residuals
        jacobian = self.jacobian
        return np.sum(2*weights*residuals*jacobian, axis=0)

    # @property
    # def jacobian(self):
    #     free_conc = self._freeconcentration
    #     stoichiometry = np.array(self._model.stoich)
    #     # TODO this will not work for more than one electrode
    #     active_species = self.active_species[0]
    #     # print("<<<" , stoichiometry, free_conc, active_species)
    #     nbetas = self._model.number_equilibria
    #     jac = libemf.emf_jac1(stoichiometry, free_conc, active_species, slice(nbetas))
    #     return jac

    @property
    def mask(self):
        "Masked elements."
        return libqt.bool_crossed(libqt.iter_column(self.ui.table_data, col=0))

    @mask.setter
    def mask(self, mask):
        for col, m in enumerate(itertools.tee(mask, 1+self.nelectrodes)):
            data = libqt.iter_column(self.ui.table_data, col)
            libqt.cross(widget for widget, m_ in zip(data, m) if m_)

    @property
    def nelectrodes(self):
        'Number of electrodes used.'
        return self.ui.table_params.columnCount()

    @nelectrodes.setter
    def nelectrodes(self, n: int):
        if n < 1:
            raise ValueError(f"At least one electrode is needed. {n} given.")
        with libqt.signals_blocked(self.ui.table_params) as t:
            t.setColumnCount(n)
            labels = [f"electrode{i}" for i in range(n)]
            t.setHorizontalHeaderLabels(labels)
        with libqt.signals_blocked(self.ui.table_data) as t:
            t.setColumnCount(1+n)
            labels = ["titre"] + [f"emf{i}" for i in range(n)]
            t.setHorizontalHeaderLabels(labels)

    @property
    def nelectrons(self):
        'Number of electrons for every electrode used'
        return tuple(int(i) for i in libqt.iter_row_text(self.ui.table_params,
                                                         row=3))

    @nelectrons.setter
    def nelectrons(self, n):
        # TODO check input
        # TODO reshape table if needed
        libaux.assert_type(int, *n)
        libqt.fill_column(self.ui.table_data, col=3, data=self.__num2it(n))

    @property
    def npoints(self):
        'Number of points in the data set.'
        return self.ui.table_data.rowCount()

    @npoints.setter
    def npoints(self, number: int):
        self.ui.table_data.setRowCount(number)

    @property
    def name(self):
        'The title for this set.'
        return self.ui.le_title.text()

    @name.setter
    def name(self, title):
        self._name = title
        self.ui.le_title.setText(title)

    @property
    def residuals(self):
        emfr = self.emf
        emfc = self.emf_fit
        return emfr - emfc

    @property
    def starting_volume(self):
        'Starting volume, in mL'
        return self.ui.dsb_V0.value()

    @starting_volume.setter
    def starting_volume(self, starting_volume):
        self.ui.dsb_V0.setValue(starting_volume)

    @property
    def titration_name(self):
        return self.ui.cb_titration.currentText()

    @property
    def titre(self):
        'Volume of titre used for every point.'
        return (float(i) for i in libqt.iter_column_text(self.ui.table_data,
                                                         col=0))

    @titre.setter
    def titre(self, volume):
        # TODO check input
        # TODO reshape table if needed
        libqt.fill_column(self.ui.table_data, col=0, data=volume, formatting="{:.4f}")

    @property
    def volume_error(self):
        'Error in volume in mL.'
        return self.ui.dsb_errV.value()

    @volume_error.setter
    def volume_error(self, verror):
        self.ui.dsb_errV.setValue(verror)

    def _nelectrs_changed(self, nelectrs):
        self.ui.table_params.setColumnCount(nelectrs)
        self.ui.table_data.setColumnCount(nelectrs + 1)
        if nelectrs == 1:
            self.ui.table_params.setHorizontalHeaderLabels(('electrode',))
            self.ui.table_data.setHorizontalHeaderLabels(('titre', 'emf'))
        else:
            lbls1 = ['electrode{}'.format(i) for i in range(nelectrs)]
            lbls2 = ['titre'] + ['emf({})'.format(i) for i in range(nelectrs)]
            self.ui.table_params.setHorizontalHeaderLabels(lbls1)
            self.ui.table_data.setHorizontalHeaderLabels(lbls2)

    # def __titration_changed(self, index):
    #     self.titrationChanged.emit(self, index)

    def _newlectrode(self):
        nels = self.nelectrodes
        assert nels > 0
        assert self.ui.table_params.columnCount() == nels
        assert self.ui.table_data.columnCount() + 1 == nels
        self.ui.table_params.setColumnCount(nels + 1)
        self.ui.table_data.setColumnCount(nels + 2)

    def _dellectrode(self):
        nels = self.nelectrodes
        assert nels > 0
        assert self.ui.table_params.columnCount() == nels
        assert self.ui.table_data.columnCount() + 1 == nels
        self.ui.table_params.setColumnCount(nels - 1)
        self.ui.table_data.setColumnCount(nels)

    def _iter_emf(self):
        for col in range(1, self.ui.table_data.columnCount()):
            txts = libqt.iter_column_text(self.ui.table_data, col=col)
            yield (float(i) for i in txts)

    def __titr_params_changed(self, row, col):
        if col == 0:
            newlabel = self.ui.table_titration.item(row, col).text()
            assert isinstance(newlabel, str)
            # TODO change label in table_params
            # self.labelsChanged.emit(row, newlabel)

    @staticmethod
    def __num2it(n):
        # if isinstance(n, (int, float)):
        #     return (n,)
        # else:
        #     return n
        return (n,) if isinstance(n, (int, float)) else n
