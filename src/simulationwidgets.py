"""Module for simulation widgets.

.. module:: simulationwidgets.py
.. moduleauthor:: Salvador Blasco <salvador.blasco@gmail.com>
"""

from PyQt5 import QtCore, QtWidgets
import numpy as np

import libaux
from libeq import LIBEQ_ENGINE
import libeq.consol
import libmath
import libqt
import libplot

import ui_sdistri
import ui_titrsimu


class SimulationData(QtWidgets.QWidget):
    '''Parent class for the data for titration and speciation simulations.

    Signals:
        - labelsChanged(int, str): emitted when the user changes any label.
        - plotUpdate: emitted when the user changes any parameter that affects
            the plot.

    .. seealso:: :class:`SpeciationWidget` and :class:`TitrationWidget`.
    '''
    # labelsChanged = QtCore.pyqtSignal(int, str)
    plotUpdate = QtCore.pyqtSignal()

    def __init__(self, model):
        super().__init__()
        self._default_n = 50
        self._model = model
        self._maintable = None
        self._labelcolumn = 0
        self._combo_n = None
        self._sbconfidence = None
        self._cbploterrors = None
        self._freeconcentrations = None
        self._freeconcentrationerrors = None
        self._titration_table_label_column = 0

    # ------------------
    # ↓↓   methods    ↓↓
    # ------------------

    def clear_concentration(self):
        "Set free concentrations to None."
        self._freeconcentrations = None
        self._freeconcentrationerrors = None

    # -------------------
    # ↓↓     slots     ↓↓
    # -------------------

    @QtCore.pyqtSlot(int, str)
    def add_component(self, position, label):
        # Warning: when this routine is invoked usually self.model has changed.
        # do not use self.model.number_components
        table = self._maintable
        number_components = table.rowCount()
        print(f"position={position}, label={label}, {number_components}")
        # table = self.ui.table_titration
        if position > number_components:
            insert = number_components
        elif position < 0:
            insert = 0
        elif position == number_components:
            insert = number_components
        else:
            insert = 1 + position
        px = list(self.pX())

        with libqt.table_locked(table):
            table.insertRow(insert)
            for col in range(table.columnCount() - 1):
                self._set_main_table_defaults(insert, col)
            px.insert(position, False)
            self.set_pX(px)

    @QtCore.pyqtSlot(int, str)
    def updateLabel(self, position, new_label):
        'Slot for when a label changes'
        item = QtWidgets.QTableWidgetItem(new_label)
        item.setFlags(QtCore.Qt.ItemIsSelectable)
        self._maintable.setItem(position, self._labelcolumn, item)

    # ------------------
    # ↓↓  properties  ↓↓
    # ------------------

    @property
    def model(self):
        return self._model

    def confidence(self):
        return self._sbconfidence.value()/100.0

    def set_confidence(self, confidence):
        return self._sbconfidence.setValue(confidence)

    def free_components(self):
        return self._free_components()

    def free_concentration(self):
        'The :term:`free concentrations array`'
        return self._freeconcentrations

    def set_free_concentration(self, c):
        self._freeconcentrations = c

    def set_free_concentration_error(self, c):
        self._freeconcentrationerrors = c

    def extended_labels(self, format_):
        return libplot.make_labels(self.labels(),
                                   np.array(self._model.stoich), format_)

    def free_concentration_error(self):
        "The error in the :term:`free concentrations array`."
        return self._freeconcentrationerrors

    def groups(self):
        stoich = libaux.extend_stoich(np.array(self.model.stoich))
        xref = self.referencex()[0]
        yref = self.referencey()
        return libaux.extract_groups(stoich, xref, yref)

    def labels(self):
        return tuple(libqt.iter_column_text(self._maintable,
                                            col=self._labelcolumn))

    def n_points(self):
        return self._default_n if self._combo_n is None else self._combo_n.value()

    def plot_errors(self):
        'bool: whether the error will be plotted or not'
        return self._cbploterrors.isChecked()

    def plot_labels(self, format_='latex'):
        ref = self.referencey()
        all_labels = self.extended_labels(format_)
        if ref is not None and ref >= 0:
            n_comp = len(self._model.stoich[0])
            z = n_comp * [0]
            z[ref] = 1
            z.extend([q[ref] for q in self._model.stoich])
            retval = [a for a, w in zip(all_labels, z) if w]
        else:
            retval = all_labels
        return retval

    def referencey(self):
        'The reference species for relative concentration plot or None.'
        i = self.ui.cb_Yaxis.currentIndex()
        return i-1 if i else None

    def set_labels(self, labels):
        'Set the labels for the components.'
        libqt.fill_column(self._maintable, col=self._labelcolumn, data=labels)

    def set_n_points(self, number: int):
        'Set the number of points to plot.'
        self._combo_n.setValue(number)

    def _populate_cb_Yaxis(self):
        combo = self.ui.cb_Yaxis
        combo.clear()
        combo.addItem("conc.")
        for label in self.labels():
            combo.addItem('%' + label)

    def _recalc_free_concentration_error(self):
        # TODO use @memoize here
        ec = libeq.consol.beta_uncertainty(self._freeconcentrations,
                                           np.array(self._model.beta_error),
                                           np.array(self._model.stoich),
                                           self.confidence())
        assert self._freeconcentrations.shape == ec.shape
        self._freeconcentrationerrors = ec

    def _set_main_table_defaults(self, row, col):
        table = self._maintable
        assert -1 < col < table.columnCount(), f"col={col}"
        assert -1 < row < table.rowCount(), f"row={row}"
        if col == 0:    # labels
            labels = self.model.labels
            item = QtWidgets.QTableWidgetItem(labels[row])
            item.setFlags(QtCore.Qt.ItemIsSelectable)
        elif col in (1, 2):
            item = QtWidgets.QTableWidgetItem('0.0000')
        else:
            raise ValueError("0 <= col <=2")

        table.setItem(row, col, item)

    def _stamp_labels(self):
        libqt.fill_column(self._maintable,
                          col=0,
                          data=self._model.labels,
                          flags=QtCore.Qt.ItemIsSelectable)

    def _yindexch(self, index: int):
        assert isinstance(index, int)
        self.plotUpdate.emit()


class SpeciationWidget(SimulationData):
    '''Species distribution widget.

    Species distribution represents is a speciation diagram where the
    independent variable is one of the components.
    '''
    def __init__(self, model):
        '''Initiate widget
        '''
        super().__init__(model)
        self.ui = ui_sdistri.Ui_SpeciationWidget()
        self.ui.setupUi(self)

        self._checkboxes = QtWidgets.QButtonGroup()
        self._checkboxes.setExclusive(False)

        self._maintable = self.ui.table_species
        self._labelcolumn = 0
        self._combo_n = self.ui.sb_NPoints
        self._sbconfidence = self.ui.sbError
        self._cbploterrors = self.ui.cbPlotErrors

        self.ui.table_species.setRowCount(model.number_components)

        self.set_initial_concentration(model.number_components*(0.1,))
        self.set_final_concentration(model.number_components*(0.1,))
        self.set_pX(model.number_components*(False,))
        self._stamp_labels()

        self.ui.table_species.cellChanged.connect(self.__tableEdited)
        # self._checkboxes.buttonClicked.connect(self.__checkboxticked)
        self._checkboxes.buttonClicked.connect(libqt.checkbox_yesno)
        self.__populate_cb_Xaxis()
        self._populate_cb_Yaxis()

        self.ui.sbError.valueChanged.connect(lambda n: self.plotUpdate.emit())
        self.ui.cb_Yaxis.currentIndexChanged.connect(self._yindexch)

    def __tableEdited(self, row, col):
        'Slot for when the table is edited.'
        tx = self.ui.table_species.item(row, col).text()
        if col == 0:
            # labels have been edited
            # self.data.labels[row] = tx
            self._populate_cb_Yaxis()
            self.__populate_cb_Xaxis()
            # self.labelsChanged.emit(col, tx)
        else:
            # try:
            #     self._recalc_free_concentration()
            # except:
            #     pass
            # else:
            #     self.plotUpdate.emit()
            #     self.clear_concentration()
            self.clear_concentration()
            # self.plotUpdate.emit()

    def analyticalc(self):
        'The :term:`total concentrations array`'
        analc = tuple(zip(self.initial_concentration(), self.final_concentration()))
        return libaux.build_T(analc, self.pX(), self.n_points())

    def calc_free_concentration(self):
        xref = self.referencex()[0]
        # beta = np.array(self._model.beta)
        # stoich = np.array(self._model.stoich)
        # analc = np.array(self.analyticalc())
        kwargs = {'first': False}

        # TODO replace by frozen_params
        # x_values, frozen_beta, frozen_stoich, frozen_analc = \
        #     libeq.consol.freeze_concentration(beta, stoich, analc, xref)
        x_values, frozen_beta, frozen_stoich, frozen_analc = self.frozen_params()

        if self._freeconcentrations is None:
            initial_values = libeq.consol.initial_guess(frozen_beta, frozen_stoich,
                                                        frozen_analc, **kwargs)
            kwargs['first'] = (LIBEQ_ENGINE == 'fortran')
        else:
            actual_n = self._freeconcentrations.shape[0]
            wanted_n = self.n_points()
            _fc = self.free_components()
            assert _fc.shape[1] == self._model.number_components
            if actual_n == wanted_n:
                initial_values = _fc
            else:
                #from libmath import sample_size_change
                initial_values = libmath.sample_size_change(_fc, wanted_n)

        new_concs = libeq.consol.consol(frozen_beta, frozen_stoich, frozen_analc,
                                        initial_values, forcer=False, **kwargs)
        fullc = np.insert(new_concs, xref, x_values, axis=1)
        self.set_free_concentration(fullc)
        self._recalc_free_concentration_error()

    def frozen_params(self):
        xref = self.referencex()[0]
        beta = np.array(self._model.beta)
        stoich = np.array(self._model.stoich)
        analc = np.array(self.analyticalc())

        x_values, frozen_beta, frozen_stoich, frozen_analc = \
            libeq.consol.freeze_concentration(beta, stoich, analc, xref)
        return x_values, frozen_beta, frozen_stoich, frozen_analc

    def initial_concentration(self):
        'The total concentration at the beginning of the titration'
        txtdata = libqt.iter_column_text(self.ui.table_species, col=1)
        return tuple(float(n) for n in txtdata)

    def set_initial_concentration(self, concentration):
        libqt.fill_column(self.ui.table_species, col=1, data=concentration)

    def final_concentration(self):
        'The total concentration at the end of the titration'
        txtdata = libqt.iter_column_text(self.ui.table_species, col=2)
        return tuple(float(n) for n in txtdata)

    def set_final_concentration(self, concentration):
        libqt.fill_column(self.ui.table_species, col=2, data=concentration)

    def pX(self):
        '''Flag that indicates whether the initial/final concentrations are
        expressed in logarithmic unit or not.'''
        checkboxes = libqt.iter_column_widget(self.ui.table_species, col=3)
        return tuple(cb.isChecked() for cb in checkboxes)

    def set_pX(self, pX):
        self.__clearcheckboxes()
        tags = ('yes' if t else 'no' for t in pX)
        for row, (c, tag) in enumerate(zip(pX, tags)):
            checkbox = QtWidgets.QCheckBox(tag)
            checkbox.setChecked(c)
            self._checkboxes.addButton(checkbox)
            self.ui.table_species.setCellWidget(row, 3, checkbox)

    def referencex(self):
        """The reference species for relative concentration plot or None.

        Returns:
            tuple: first element is the index of the reference and second
                index is a bool indicating whether to use potential or not.
        """
        idx = self.ui.cb_Xaxis.currentIndex()
        return idx // 2, bool(idx % 2)

    def set_referencex(self, index: int, plog=False):
        'Set component to be *X* reference.'
        assert isinstance(index, int)
        assert index >= 0
        ext = 1 if plog else 0
        self.ui.cb_Xaxis.setCurrentIndex(index*2 + ext)

    def set_referencey(self, index: int):
        'Set component to be *Y* reference.'
        assert isinstance(index, int)
        self.ui.cb_Yaxis.setCurrentIndex(index+1)

    def xdata(self):
        'The X data to be plotted'
        reference, ispot = self.referencex()
        analc = np.array(self.analyticalc())
        return -np.log10(analc.T[reference]) if ispot else analc.T[reference]

    def xlabel(self):
        'label for the x-axis of the plot'
        return self.ui.cb_Xaxis.currentText()

    def xscalelog(self):
        'bool: True if x-scale is logarithmic or (False) linear.'
        return ((self.ui.cb_Xaxis.currentIndex() + 1) % 2) == 0

    def ydata(self):
        'The Y data to be plotted'
        if self.free_concentration() is None:
            retval = None
        if self.referencey() is None:
            retval = self.free_concentration()
        else:
            retval = libaux.percent_distribution(self.free_concentration(),
                                                 np.array(self._model.stoich),
                                                 np.array(self.analyticalc()),
                                                 self.referencey())
        return retval

    def yerrdata(self):
        'The error for Y data'
        if self.free_concentration() is None:
            return None
        if self.referencey() is None:
            return self.free_concentration_error()

        return libaux.percent_distribution(self.free_concentration_error(),
                                           np.array(self._model.stoich),
                                           np.array(self.analyticalc()),
                                           self.referencey())

    def ylabel(self):
        'label for the y-axis of the plot'
        ptype = self.ui.cb_Yaxis.currentText()
        if ptype[0] == 'c':
            return 'Concentration'
        if ptype[0] == '%':
            return '% Formation Relative to ' + ptype[1:]

        raise RuntimeError("Something is wrong with cb_Yaxis")

    def _free_components(self):
        if self._freeconcentrations is None:
            return None

        # nc = self.model.number_components
        # return np.delete(self._freeconcentrations[:, :nc],
        #                  self.referencex(), axis=1)
        return self._freeconcentrations[:, :self.model.number_components]

    def __populate_cb_Xaxis(self):
        combo = self.ui.cb_Xaxis
        combo.clear()
        for label in self.labels():
            combo.addItem(label)
            combo.addItem('p' + label)

    def __clearcheckboxes(self):
        buttons = self._checkboxes.buttons()
        for button in buttons:
            self._checkboxes.removeButton(button)


class TitrationWidget(SimulationData):
    '''Widget for a titration simulation.'''
    def __init__(self, model):
        """Initiate widget. """
        super().__init__(model)
        self.ui = ui_titrsimu.Ui_TitrationWidget()
        self.ui.setupUi(self)
        self._maintable = self.ui.table_titration

        self.ui.table_titration.setRowCount(model.number_components)
        self._stamp_labels()
        self.ui.table_titration.cellChanged.connect(self.__table_edited)

        self._maintable = self.ui.table_titration
        self._labelcolumn = 0
        self._combo_n = self.ui.sb_NPoints
        self._sbconfidence = self.ui.sbError
        self._cbploterrors = self.ui.cbPlotErrors

        self.__populate_cb_Xaxis()
        self._populate_cb_Yaxis()

    def buret(self):
        """The buret values.

        Returns:
            tuple of floats containing the values of this parameter.
        """
        ctext = libqt.iter_column_text(self.ui.table_titration, col=2)
        return tuple(float(b) for b in ctext)

    def initial_amount(self):
        """The inital amount values in millimole.

        Returns:
            tuple of floats containing the values of this parameter.
        """
        ctext = libqt.iter_column_text(self.ui.table_titration, col=1)
        return tuple(float(b) for b in ctext)

    def set_initial_amount(self, initial_amount):
        """The initial amount values in millimole.

        Parameters:
            initial_amount (sequence of floats): the values for this parameter.
                The length must be equal to the number of components.
        """
        libqt.fill_column(self.ui.table_titration, col=1, data=initial_amount)

    def set_buret(self, buret):
        """Set the buret parameter.

        Parameters:
            buret (sequence of floats): the values for buret. The length must be equal
                to the number of components.
        """
        libqt.fill_column(self.ui.table_titration, data=buret, col=2)

    def calc_free_concentration(self):
        'The :term:`free concentrations array`'
        beta = np.array(self._model.beta)
        stoich = np.array(self._model.stoich)
        analc = np.array(tuple(self.analyticalc()))
        kwargs = {'first': False}

        if self._freeconcentrations is None:
            initial_values = libeq.consol.initial_guess(beta, stoich, analc, **kwargs)
            kwargs['first'] = (LIBEQ_ENGINE == 'fortran')
        else:
            actual_n = self._freeconcentrations.shape[0]
            wanted_n = self.n_points()
            if actual_n == wanted_n:
                initial_values = self._freeconcentrations
            else:
                initial_values = libmath.sample_size_change(self._freeconcentrations, wanted_n)

        self._freeconcentrations = libeq.consol.consol(beta, stoich, analc,
                                                       initial_values, **kwargs)
        self._recalc_free_concentration_error()

    def starting_volume(self):
        'The starting volume in mL.'
        return self.ui.dsb_V0.value()

    def set_starting_volume(self, volume: float):
        'Set the starting volume in mL.'
        assert isinstance(volume, float)
        self.ui.dsb_V0.setValue(volume)

    def final_volume(self):
        'The value of the final titre in mL.'
        return self.ui.dsb_Vf.value()

    def set_final_volume(self, volume: float):
        'Set the final volume in mL.'
        assert isinstance(volume, float)
        # if volume > self.starting_volume():
        #     raise ValueError('Final volume must be bigger than start volume')
        self.ui.dsb_Vf.setValue(volume)

    def titre(self):
        """The titre values in mL.

        Yields:
            float: the value of the titre.
        """
        yield from ((self.final_volume() - self.starting_volume())
                    / self.n_points() * i for i in range(self.n_points()))

    def volume_increment(self):
        'The step volume in mL.'
        return (self.final_volume() - self.starting_volume()) / self.n_points()

    def analyticalc(self):
        "The :term:`total concentrations array`"
        return libaux.build_T_titr2(self.initial_amount(), self.buret(),
                                    self.starting_volume(), self.titre())

    def xdata(self):
        'The data for the x-axis.'
        n = self.ui.cb_Xaxis.currentIndex()
        if n == 0:
            retval = np.linspace(0, self.final_volume()-self.starting_volume(), self.n_points())
        if n == 1:
            retval = np.linspace(self.starting_volume(), self.final_volume(), self.n_points())
        else:
            xref = n // 2 - 1
            pot = bool(n % 2)
            xconc = self.free_concentration()[:, xref]
            retval = -np.log10(xconc) if pot else xconc
        return retval

    def xlabel(self):
        'The x-axis label.'
        return self.ui.cb_Xaxis.currentText()

    def ydata(self):
        'The Y data to be plotted'
        n = self.ui.cb_Yaxis.currentIndex()
        if n == 0:
            return self.free_concentration()

        reference = (n - 1) % 2
        return libaux.percent_distribution(self.free_concentration(),
                                           np.array(self.model.stoich),
                                           np.array(tuple(self.analyticalc())),
                                           reference)

    def ylabel(self):
        'The y-axis label.'
        ptype = self.ui.cb_Yaxis.currentText()
        if ptype[0] == 'c':
            return 'Concentration'
        if ptype[0] == '%':
            return '% Formation Relative to ' + ptype[1:]

        raise RuntimeError("Something is wrong with cb_Yaxis")

    @QtCore.pyqtSlot(int, str)
    def add_component(self, position, label):
        """Add a new component in the table.

        Parameters:
            position (int):
            label (str):
        """
        # Warning: when this routine is invoked usually self.model has changed.
        # do not use self.model.number_components
        number_components = self._maintable.rowCount()
        insert = libqt.clamp_range(position, number_components)

        with libqt.signals_blocked(self._maintable) as table:
            table.insertRow(insert)
            for col in range(table.columnCount()):
                self._set_main_table_defaults(insert, col)

    def __table_edited(self, row, col):
        'Slot for when the table is edited.'
        # edited_txt = self.ui.table_titration.item(row, col).text()
        if col == 0:         # labels have been edited
            # self.data.labels[row] = edited_txt
            self._populate_cb_Yaxis()
            self.__populate_cb_Xaxis()
            # self.labelsChanged.emit(col, edited_txt)
        else:
            self.plotUpdate.emit()
            self.clear_concentration()

    def __populate_cb_Xaxis(self):
        combo = self.ui.cb_Xaxis
        combo.clear()
        combo.addItem('V (titre)')
        combo.addItem('V (total)')
        for label in self.labels():
            combo.addItem(label)
            combo.addItem('p' + label)
