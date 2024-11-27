"""
Widget for the model containing stoichiometry values and beta values.

.. module:: modelwidget.py
.. moduleauthor:: Salvador Blasco <salvador.blasco@gmail.com>
"""


# from contextlib import contextmanager
import copy

from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5 import QtGui

import consts
import libaux
import libqt
import ui_model


class ModelWidget(QtWidgets.QWidget):
    '''Chemical model.

    This class contains and manages the data corresponding to a model,
    namely the equilibrium constants, the stoichiometric coefficients and
    other.

    Signals:
        - labelsChanged: when the user changes a label
        - componentAdded: when the user adds a new component
        - componentDeleted: when the user deletes a component
    Slots:
        - updateLabel(position, new_label)
    '''

    labelsChanged = QtCore.pyqtSignal(int, str)
    componentAdded = QtCore.pyqtSignal(int, str)
    componentDeleted = QtCore.pyqtSignal(int)
    equilibriaChanged = QtCore.pyqtSignal()

    def __init__(self):
        super().__init__()

        self.__badbg = QtGui.QBrush(QtGui.QColor("red"))
        self.__okbg = QtGui.QBrush(QtGui.QColor("white"))

        self.__iscalori = False                 # TODO remove
        self._temperature = 298.15              # TODO remove
        self.__currentmodel = ModelData()
        self.__models = [self.__currentmodel]
        self.__labels = ['L', 'H']

        self.ui = ui_model.Ui_ModelWidget()
        self.ui.setupUi(self)

        hheader = self.ui.table_model.horizontalHeader()
        hheader.sectionDoubleClicked.connect(self.__headerdoubleclicked)

        self.ui.table_model.cellChanged.connect(self._cell_changed)

        self.setCurrentModel(0)     # must be after connections

    def addEquilibrium(self, position=None, stoich=None, value=0.0, error=0.0, flag=consts.RF_REFINE):
        '''Add new equilibrium.

        Add a new row to the table and fills it with the default values.
        The new row is inserted in the row currently selected or in the bottom
        if no row or cell is selected. The new row will be filled with defatult
        parameters.
        '''
        with libqt.signals_blocked(self.ui.table_model):
            table = self.ui.table_model
            if position is None:
                position = table.currentRow()
            if position == -1:
                position = table.rowCount()
            table.insertRow(position)
            self.__renumber()
            for col in range(1+table.columnCount()):
                if stoich is None and col < self.number_components:
                    table.setItem(position, self.column_stoich + col, QtWidgets.QTableWidgetItem('0'))
                elif stoich is not None and col < self.number_components:
                    table.setItem(position, self.column_stoich + col, QtWidgets.QTableWidgetItem(str(stoich[col])))
                elif col == self.column_beta:
                    table.setItem(position, col, QtWidgets.QTableWidgetItem(f'{value:.4f}'))
                elif col == self.column_betaerror:
                    table.setItem(position, col, QtWidgets.QTableWidgetItem(f'{error:.4f}'))
                elif col == self.column_betaflags:
                    table.setCellWidget(position, col, libqt.create_combo(consts.REFINE_BLABELS, flag+1))
        self.equilibriaChanged.emit()

    def addComponent(self, label, position=None):
        '''Add one reagent.

        Add one column to model table with a new label.

        Parameters:
            label (str): The label for the new reagent
            position (int): position where the new reagent will be inserted.
        Returns:
            int: the position into which the component has been added.
        '''

        table = self.ui.table_model
        if position is None:
            position = table.currentColumn()
        else:
            assert isinstance(position, int)
        if position not in self.range_components:
            position = self.number_components
        self.labels.insert(position, label)

        with libqt.signals_blocked(self.ui.table_model):
            self.ui.table_model.insertColumn(position)
            libqt.fill_column(table, position,table.rowCount()*(0, ))
            self.__set_horizontal_headers()

        self.componentAdded.emit(position, label)
        return position

    def append(self, model):
        'Append a new model to the list. Do not update.'
        for mod in self.__models:
            if mod.n_species != model.n_species:
                raise ValueError('Model is ill-shaped')
        self.__models.append(model)
        if len(self.__models) == 1:
            self.setCurrentModel(0)

    def clear(self, n_equils=2, n_species=2):
        """Delete all models clear table. Create empty model.

        Parameters:
            n_equils (int, optional): the number of rows for the new cleared table.
            n_species (int, optional): the number of columns minun 3 for the new cleared table.
        """
        self.__models = []
        # self.__iscalori = False
        self.ui.table_model.clearContents()
        self.ui.table_model.setRowCount(n_equils)
        self.ui.table_model.setColumnCount(3 + n_species)
        # self.__iscalori = False
        self.__labels = _default_labels(n_species)

        self.__set_horizontal_headers()
        self.__renumber()
        for row in range(self.ui.table_model.rowCount()):
            for col in range(self.ui.table_model.columnCount()):
                self.__setdefaults(row, col)

    def current_model(self):
        """The current model.

        Saves the current model and returns it."""
        self.save_current_model()
        return self.__currentmodel

    def newModel(self):
        '''Append a new, empty model.

        The new model will be filled with default values.

        Returns:
            ModelData: the newly create model
        '''
        # save current model
        self.save_current_model()
        newmodel = copy.deepcopy(self.__currentmodel)
        newmodel.name += '(copy)'
        self.append(newmodel)
        return newmodel

    def removeEquilibrium(self, position=None, confirm=True):
        '''Remove one beta from current model.

        Parameters:
            position (optional, int): index of the equilibrium to be deleted
            confirm (bool, optional): raise dialog to ask confirmation
        Raises:
            ValueError: If there is only one equilibrium left or **position**
                is a meaningless value.
        '''
        if self.number_equilibria == 1:
            raise ValueError('Only remaining equilibrium cannot be deleted')

        if (not position) or (position is None):
            position = self.ui.table_model.currentRow()
            if position == -1:
                position = self.ui.table_model.rowCount() - 1
            self.ui.table_model.selectRow(position)
        else:
            if position < 0 or position > self.number_equilibria:
                raise ValueError("Meaningless value of 'position'")

        if confirm:
            msg = "Do you want to delete this equilibrium?"
            delete = libqt.confirm_delete(parent=self, msg=msg)
        else:
            delete = True

        if delete:
            with libqt.signals_blocked(self.ui.table_model):
                self.ui.table_model.removeRow(position)
        self.__renumber()

    def removeComponent(self, which=None):
        '''Remove one reagent from all models.

        Remove selected reagent from current table and from all the models.
        If **position** is not given, a dialog will pop to ask for one.

        Parameters:
            which (optional, int or str): the index or label of the component to be deleted
        Emits:
            componentDeleted(int)
        Raises:
            ValueError: if only two components are left and therefore neither
                can be deleted or if selected column is not a component
                column.
        .. note:: Use with care, can cause data loss.
        '''
        if self.number_components < 3:
            raise ValueError('Model must have at least two components')

        match which:
            case str():
                if which not in self.labels:
                    raise ValueError(f'component {which} does not exist')
                position = self.labels.index(which)
                confirm = False
            case int():
                if which not in range(self.number_components):
                    raise ValueError('Cannot delete this column')
                position = which
                confirm = False
            case None:
                position = self.ui.table_model.currentColumn()
                confirm = True
            case _:
                raise ValueError(f"parameter {which} not recognised")

        # table = self.ui.table_model
        if confirm:
            # c = table.currentColumn()
            # table.selectColumn(c)  # aparently this does nothing
            if not libqt.confirm_delete(parent=self,
                                    msg="Do you want to delete this reagent?"):
                return
        # if c > self.number_components:
        #     raise ValueError('Cannot delete this column')
        with libqt.signals_blocked(self.ui.table_model):
            self.ui.table_model.removeColumn(position)
        # TODO delete from all models
        self.componentDeleted.emit(position)

    def setCurrentModel(self, modelindex=0):
        '''Initialize the widget with the new values.

        Parameters:
            data (dict): The data to be loaded.
        '''
        model = self.__models[modelindex]
        self.__currentmodel = model
        columns = model.n_species + 3
        rows = model.n_equil
        table = self.ui.table_model
        with libqt.signals_blocked(self.ui.table_model):
            table.setRowCount(rows)
            table.setColumnCount(columns)
            table.clearContents()
            libqt.replace_nones(table)
            self.__set_horizontal_headers()

        self.beta_flags = model.const_flags
        self.stoich = model.stoich
        self.beta_raw = model.const
        self.beta_error = model.const_error

        # if self.__iscalori:
        #     self.enthalpy = model.enthalpy
        #     self.enthalpy_flags = model.enthalpy_flags
        #     self.enthalpy_error = model.enthalpy_error
        #     self.updateEnthropy()

    def beta_update(self, new_beta):
        raise DeprecationWarning
        with libqt.signals_blocked(self.ui.table_model):
            cells = libqt.iter_column(self.ui.table_model, self.column_beta)
            nicells = (c for i, c in zip(self.__isrowignored(), cells) if i)
            for cell, beta in zip(nicells, new_beta):
                cell.setText(str(beta))

    def updateLabel(self, position, new_label):
        "Slot for when a label changes"
        item = QtWidgets.QTableWidgetItem(new_label)
        self.ui.table_model.setHorizontalHeaderItem(position, item)

    # ------------------
    # ↓↓  properties  ↓↓
    # ------------------

    @property
    def beta(self):
        'The not ignored values of the :term:`equilibrium constants array`.'
        return tuple(10**i for i in self.__retcol(col=self.column_beta))

    @property
    def beta_log(self):
        """The not ignored values of the :term:`equilibrium constants array`
        in logarithmic units."""
        return tuple(self.__retcol(col=self.column_beta))

    @property
    def beta_raw(self):
        "The values directly read from column 0 of the table."
        itercol = libqt.iter_column_text(self.ui.table_model,
                                         col=self.column_beta)
        return tuple(float(i) for i in itercol)

    @beta_raw.setter
    def beta_raw(self, beta):
        with libqt.signals_blocked(self.ui.table_model):
            libqt.fill_column(self.ui.table_model,
                              self.column_beta,
                              beta,
                              formatting='{:.4f}')

    @property
    def beta_error(self):
        'The errors for the :term:`equilibrium constants array`.'
        # TODO use memoize here if possible
        return tuple(self.__retcol(col=self.column_betaerror))

    @beta_error.setter
    def beta_error(self, betaerr):
        with libqt.signals_blocked(self.ui.table_model):
            libqt.fill_column(self.ui.table_model, self.column_betaerror,
                              betaerr,
                              formatting='{:.4f}')

    @property
    def beta_flags(self):
        'The refinement flags for the :term:`equilibrium constants array`.'
        indices = libqt.iter_column_comboindex(self.ui.table_model,
                                               self.column_betaflags)
        return tuple(item - 1 for item in indices if item != 0)

    @beta_flags.setter
    def beta_flags(self, flags):
        assert all(tuple(map(lambda x: -2 < x < 8, flags)))
        with libqt.signals_blocked(self.ui.table_model):
            libqt.fill_column_comboindex(self.ui.table_model,
                                         (i+1 for i in flags),
                                         consts.REFINE_BLABELS,
                                         self.column_betaflags)

    @property
    def column_stoich(self) -> int:
        'The initial column for stoichiometry values.'
        return 0

    @property
    def column_beta(self) -> int:
        'The column where the beta values are suposed to be.'
        return self.number_components

    @property
    def column_betaerror(self) -> int:
        'The column where the beta error values are suposed to be.'
        return 1 + self.number_components

    @property
    def column_betaflags(self) -> int:
        'The column where the beta flags are suposed to be.'
        return 2 + self.number_components

    @property
    def extended_labels(self):
        labelstyle = 'plain'
        yield from libaux.extended_labels(self.labels, self.stoich, kind=labelstyle)

    @property
    def labels(self):
        'List of labels for every component.'
        self.__labels = [self.ui.table_model.horizontalHeaderItem(col).text()
                         for col in self.__srange()]
        return self.__labels

    @labels.setter
    def labels(self, labels):
        "Update table header when reagent labels change."
        if self.number_components != len(labels):
            columns = len(labels) + 3
            self.ui.table_model.setColumnCount(columns)
        header_labels = labels + ['Value', 'Error', 'Flag']
        self.ui.table_model.setHorizontalHeaderLabels(header_labels)
        self.__labels = labels

    @property
    def modelname(self):
        'The name of the model.'
        return self.__currentmodel.name

    @modelname.setter
    def modelname(self, name):
        self.__currentmodel.name = name

    @property
    def number_equilibria(self):
        "Number of equilibria in model"
        return self.ui.table_model.rowCount()

    @property
    def number_components(self):
        'The number of components for this model.'
        col_count = self.ui.table_model.columnCount()
        return col_count - (8 if self.__iscalori else 3)

    @property
    def range_components(self):
        'Iterable range from 0 to number_components.'
        return range(self.number_components)

    @property
    def stoich(self):
        "The :term:`stoichiometry array` of the current model."
        rows = tuple(self.__unignored_rows())
        cols = self.__srange()
        array = libqt.tab2array(self.ui.table_model, rows, cols, dtype=int)
        return array

    @stoich.setter
    def stoich(self, value):
        with libqt.signals_blocked(self.ui.table_model):
            libqt.array2tab(self.ui.table_model, value,
                            col0=self.column_stoich)

    # TODO probably not necesary
    @property
    def temperature(self):
        'Absolute temperature of the system.'
        return self._temperature

    # TODO probably not necesary
    @temperature.setter
    def temperature(self, temp):
        if temp <= 0:
            raise ValueError('Negative temperature')
        self._temperature = temp

    def _cell_changed(self, row, col):
        widget = self.ui.table_model.item(row, col)
        float_cols = (self.column_beta, self.column_betaerror)
        if col in float_cols:
            itype = float
        else:
            itype = int

        txt = widget.text()
        try:
            new_txt = str(_evaluate(txt, itype))
        except ValueError:
            new_txt = '????'

        with libqt.signals_blocked(self.ui.table_model):
            widget.setText(new_txt)

    def save_current_model(self):
        "Update info in ModelData object."
        cmod = self.__currentmodel
        # cmod.name = self.model
        cmod.const = self.beta_raw
        cmod.stoich = self.stoich
        cmod.const_flags = self.beta_flags
        cmod.const_error = self.beta_error

    def __isrowignored(self):
        "Bools listing whether row is to be used (True) or ignored (False)."
        for flag in self.__rawflags():
            yield flag != 0

    def __unignored_rows(self):
        "Bools listing whether row is to be used (True) or ignored (False)."
        return (n for n, flag in enumerate(self.__rawflags()) if flag != 0)

    def __headerdoubleclicked(self, i):
        if i >= self.number_components:
            return

        label_dialog, ok = QtWidgets.QInputDialog.getText(self, "Rename component",
                                                          "Enter new label")
        if not ok:
            return

        self.updateLabel(i, label_dialog)
        ## lbl = self.labels
        ## lbl[i] = label_dialog
        ## self.__set_horizontal_headers()
        # new_item = QtWidgets.QTableWidgetItem(label_dialog)
        # self.ui.table_model.setHorizontalHeaderItem(i, new_item)
        self.labelsChanged.emit(i, label_dialog)

    def __rawflags(self):
        '''All flags, including ignored.'''
        yield from libqt.iter_column_comboindex(self.ui.table_model,
                                                col=self.column_betaflags)

    def __renumber(self):
        table = self.ui.table_model
        rows = table.rowCount()
        table.setVerticalHeaderLabels(['%d' % (i+1) for i in range(rows)])

    def __retcol(self, col, dtype=float):
        'Return contents of a column as an iterable if not ignored'
        itercol = libqt.iter_column_text(self.ui.table_model, col)
        izip = zip(itercol, self.__isrowignored())
        return (dtype(item) for item, use in izip if use)

    def __retcolcal(self, col, dtype=float):
        if not self.isCalori():
            return None
        return self.__retcol(col, dtype)

    def __retfullcol(self, col):
        'Return contents of a column as an iterable'
        yield from libqt.iter_column_text(self.ui.table_model, col)

    def __setdefaults(self, row, col):
        table = self.ui.table_model
        if col < self.number_components:
            table.setItem(row, self.column_stoich + col,
                          QtWidgets.QTableWidgetItem('0'))
        elif col in (self.column_beta, self.column_betaerror):
            table.setItem(row, col,
                          QtWidgets.QTableWidgetItem('0.0000'))
        elif col == self.column_betaflags:
            table.setCellWidget(row, col, libqt.create_combo(consts.REFINE_LABELS, 1))

    def __set_horizontal_headers(self):
        header_labels = self.ui.table_model.columnCount() * [None]
        header_labels[self.column_beta] = 'Value'
        header_labels[self.column_betaerror] = 'Error'
        header_labels[:self.number_components] = self.__labels
        header_labels[self.column_betaflags] = 'Flag'
        self.ui.table_model.setHorizontalHeaderLabels(header_labels)

    def __srange(self):
        return range(self.column_stoich,
                     self.column_stoich+self.number_components)

    # ---------------------
    # ↓↓  magic methods  ↓↓
    # ---------------------
    # Implement list-like behaviour

    def __delitem__(self, key):
        del self.__models[key]

    def __getitem__(self, key):
        return self.__models[key]

    def __iter__(self):
        # self.__savecurrentmodel()
        self.save_current_model()
        return iter(self.__models)

    def __len__(self):
        return len(self.__models)

    def __reversed__(self):
        return reversed(self.__models)

    def __setitem__(self, key, value):
        self.__models[key] = value


class ModelData:
    "Data class for the collection of models."
    def __init__(self, n_equils=1, n_species=2):
        """Create new, empty ModelData.

        Parameters:
            n_equils (int, optional):
            n_species (int, optional):
        """
        self.name = 'No name'
        self.const = n_equils * [1.0]
        self.stoich = [[1 for n in range(n_species)] for m in range(n_equils)]
        self.const_flags = n_equils * [consts.RF_CONSTANT]
        self.const_error = n_equils * [0.0]

    @property
    def n_equil(self):
        'Number of equilibria.'
        if self.stoich is None:
            retv = None
        else:
            retv = len(self.stoich)
        return retv

    @property
    def n_species(self):
        'Number of species.'
        if self.stoich is None:
            retv = None
        else:
            retv = len(self.stoich[0])
        return retv


class SolidModelWidget(QtWidgets.QWidget):
    '''Model extension for solid equilibria.

    This class contains and manages the data corresponding to a model,
    namely the equilibrium constants, the stoichiometric coefficients and
    other.
    '''

    def __init__(self, parent_model):
        super().__init__()
        self.ui = ui_model.Ui_ModelWidget()
        self.ui.setupUi(self)
        self.__parent_model = parent_model


def _evaluate(item, itype):
    ret = None
    try:
        ret = itype(item)
    except ValueError:
        ret = itype(eval(item))
    return ret


def _default_labels(n_species):
    if n_species == 1:
        return ['L']
    if n_species == 2:
        return ['L', 'H']
    if n_species == 3:
        return ['L', 'M', 'H']
    ascii_uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if n_species < 27:
        return ascii_uppercase[:n_species]
    else:
        t = itertools.product(iter('1234'), iter(string.ascii_uppercase))
        return ["".join(reversed(next(t))) for _ in range(n_species)]
