import numpy as np
from PyQt5 import QtCore
from PyQt5 import QtGui
from PyQt5 import QtWidgets

import datawidget
import libaux
import libqt
import ui_specds


class SpecWidget(datawidget.DataWidget):
    def __init__(self, model):
        super().__init__(model)
        QtWidgets.QWidget.__init__(self)
        self.ui = ui_specds.Ui_SpecDSWidget()
        self.ui.setupUi(self)
        # self._connectMenu(self.ui.tab_data)
        self._colored_checkboxes = []

        #  self._fitdata = fitdata
        #  self._menu_location = -1
        self._islog = model.number_components * [False]

        #  self.setObjectName("Spectrum data")
        #  self.verticalLayout = QtWidgets.QVBoxLayout(self)
        #  self.verticalLayout.setObjectName("verticalLayout")
        #  self.table_data = QtWidgets.QTableWidget(self)
        #  self.table_data.setObjectName("tab_data")
        #  # TODO use model to reshape table

        self._populate_checkboxes()
        self.labels_changed()
        model.labelsChanged.connect(self.updateLabel)
        model.componentAdded.connect(self.componentAdded)
        model.equilibriaChanged.connect(self.labels_changed)
        self.ui.cb_titration.currentTextChanged.connect(super()._DataWidget__titration_changed)

        #  self.table_data.setRowCount(model.number_components+1)
        #  self.table_data.setColumnCount(2)
        #  labels = list(model.labels)
        #  self.table_data.setVerticalHeaderLabels(labels + ['200'])
        #  self.verticalLayout.addWidget(self.table_data)

        #  self._menuvh_wavel = QtWidgets.QMenu(parent=self)
        #  self._action_usewl = QtWidgets.QAction('use', parent=self._menuvh_wavel)
        #  self._action_usewl.setCheckable(True)
        #  self._action_usewl.triggered.connect(self.__use_wavelength)
        #  self._action_delwl = QtWidgets.QAction('delete', parent=self._menuvh_wavel)
        #  self._menuvh_wavel.addAction(self._action_usewl)
        #  self._menuvh_wavel.addAction(self._action_delwl)

        #  self._menuvh_specs = QtWidgets.QMenu(parent=self)
        #  self._action_makep = QtWidgets.QAction('is -log', parent=self._menuvh_specs)
        #  self._action_makep.triggered.connect(self._menuvh_makep)
        #  self._action_makep.setCheckable(True)
        #  self._action_makep.setCheckable(False)
        #  self._menuvh_specs.addAction(self._action_makep)

        #  # pop-up menu for tab_data
        #  vheader = self.table_data.verticalHeader()
        #  vheader.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #  vheader.customContextMenuRequested.connect(self._vheader_popup)

        #  # connect model signals
        #  # TODO maybe this should go on parent class
        #  # model.componentAdded.connect(self.add_component)
        #  # model.componentDeleted.connect(self.delete_component)

        #  # self.ui.tab_data.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #  # self.ui.tab_data.customContextMenuRequested.connect(self._pmtd)

    def add_component(self, position: int, label: str):
        '''Add one component to the system.

        Add one column to model table with a new label.

        Parameters:
            label (str): The label for the new reagent
            position (int): position where the new reagent will be inserted.
        Returns:
            int: the position into which the component has been added.
        '''
        self._islog.insert(position, False)
        with libqt.table_locked(self.table_data) as table:
            table.insertRow(position)
        self.__set_component_labels()

    def analyticalc(self):
        'Analytical concentrations'
        return libqt.tab2array(self.ui.tab_data,
                               rows=range(self._model.number_components),
                               cols=self.column_range_species)

    def copy_titration_data(self, titration_widget):
        """Copy concentrations from a titration widget.

        Parameters:
            titration_widget (:class:`simulationwidgets.TitrationWidget`): the widget to get
                the data from.
        """
        initial_amount = titration_widget.initial_amount()
        buret = titration_widget.buret()
        starting_volume = titration_widget.starting_volume()
        final_volume = titration_widget.final_volume()
        titration = libaux.build_T_titr(initial_amount, buret, starting_volume,
                                        final_volume-starting_volume, self.number_titration_points)
        for col, data in enumerate(titration):
            libqt.fill_column(self.ui.table_data, col, data, formatting='{:.4e}')

    def delete_component(self, position: int):
        with libqt.table_locked(self.table_data) as table:
            table.removeRow(position)
        del self._islog[position]

    def feed_data(self, wavelengths, spectra_data, analytical_concentrations=None):
        total_rows = 1 + spectra_data.shape[0] + self._model.number_components
        total_columns = spectra_data.shape[1]
        with libqt.table_locked(self.ui.table_data) as table:
            table.clear()
            table.setRowCount(total_rows)
            table.setColumnCount(total_columns)
            # libqt.fill_column(table, self.column_wavelength, wavelengths,
            #                   row0=self._model.number_components, formatting="{:.1f}")
            libqt.array2tab(table, spectra_data, row0=1+self._model.number_components,
                            col0=0, formatting="{:.4f}")
            if analytical_concentrations is not None:
                libqt.array2tab(table, analytical_concentrations,
                                col0=0, formatting="{:.4f}")
            # breakpoint()
            vertical_header_labels = self._model.labels + ["path"] + ["{:.1f}".format(w) for w in wavelengths]
            table.setVerticalHeaderLabels(vertical_header_labels)

        self._populate_checkboxes()
        # self._fitdata.set_wavelengths(wavelengths)

    def weight(self, scheme='auto'):
        pass

    # Slots
    @QtCore.pyqtSlot(int, str)
    def component_added(self, index, label):
        with libqt.table_locked(self.ui.table_data) as table:
            table.insertRow(index)

    @QtCore.pyqtSlot(int)
    def component_deleted(index):
        with libqt.table_locked(self.ui.table_data) as table:
            table.removeRow(index)

    @QtCore.pyqtSlot(int, str)
    def update_label(self, position, new_label):
        pass

    @property
    def column_range(self):
        return range(self.number_titration_points)

    @property
    def mask(self):
        "Masked elements."
        # return libqt.bool_crossed(libqt.iter_column(self.ui.table_data, col=0))
        pass

    @property
    def number_titration_points(self):
        return self.ui.table_data.columnCount()

    @property
    def row_optical_path(self) ->  int:
        return self._model.number_components

    @property
    def row_range_species(self):
        'Iterable range from 0 to number_components.'
        return range(self.row_optical_path)

    @property
    def row_range_data(self):
        return range(1+self.row_optical_path, self.ui.table_data.rowCount())

    @property
    def spectrum(self) -> np.ndarray:
        "Return whole spectrum as numpy.ndarray"
        return libqt.tab2array(self.ui.table_data, 
                               rows=self.row_range_data,
                               cols=self.column_range,
                               dtype=float)

    @property
    def wavelengths(self):
        "Return a generator containing wavelengths."
        yield from (float(self.ui.table_data.takeverticalHeaderItem(index).text())
                    for index in self.row_range_data)

    @QtCore.pyqtSlot(int, str)
    def updateLabel(self, position, new_label):
        self.labels_changed()

    @QtCore.pyqtSlot()
    def labels_changed(self):
        extended_labels = list(libaux.extended_labels(self._model.labels, self._model.stoich))
        with libqt.table_locked(self.ui.table_components) as t:
            self.__change_columns(len(extended_labels))
            t.setHorizontalHeaderLabels(extended_labels)

    @QtCore.pyqtSlot(int, str)
    def componentAdded(self, which, label):
        self.labels_changed()       # change all

    def _cchch(self, c):
        c.setText("coloured" if c.isChecked() else "transparent")

    def _pmtd(self, p):
        "pop-up menu for tab_data"
        print("context menu raised!")

    def _extpopupmenu(self):
        p = self._popupm
        a = p.addAction("use wavelength(s)")
        a.triggered.connect(self._tdwmuse)
        a = p.addAction("ignore wavelength(s)")
        a.triggered.connect(self._tdwmign)
        a = p.addAction("delete wavelength(s)")
        a.triggered.connect(self._tdwmdel)

    def _menuvh_makep(self):
        pass

    def _populate_checkboxes(self):
        table = self.ui.table_components
        for col in range(table.columnCount()):
            qchk = QtWidgets.QCheckBox('active')
            qchk.setChecked(True)
            table.setCellWidget(0, col, qchk)
    #     all_species = libaux.extended_labels(self._model.labels, self._model.stoich)        
    #     self._colored_checkboxes = [
    #         QtWidgets.QCheckBox(f"{species} is optically active")
    #         for species in all_species]
    #     layout = self.ui.layout_checkboxes
    #     for item in reversed(range(layout.count())):
    #         layout.itemAt(item).widget().setParent(None)
    #     for item in self._colored_checkboxes:
    #         layout.addWidget(item)

    def _tdwmuse(self):
        use = {x.column() for x in self._tabdata.selectedIndexes()}
        self.data.useWl(use)
        libqt.strikeOutCols(self._tabdata, use, False)

    def _tdwmign(self):
        ignore = {x.column() for x in self._tabdata.selectedIndexes()}
        self.data.ignoreWl(ignore)
        libqt.strikeOutCols(self._tabdata, ignore, True)

    def _tdwmdel(self):
        raise NotImplementedError

    def __change_columns(self, new_size):
        table = self.ui.table_components
        current_count = table.columnCount()
        table.setColumnCount(new_size)
        
        if current_count < new_size:
            for col in range(current_count, new_size):
                qchk = QtWidgets.QCheckBox('active')
                qchk.setChecked(True)
                table.setCellWidget(0, col, qchk)

    def __set_component_labels(self):
        # this function assumes the new rows have already been added/deleted
        with libqt.table_locked(self.table_data) as table:
            for position, (islog, label) in enumerate(zip(self._islog, self._model.labels)):
                lbl = "".join((['p'] if islog else []) + [label])
                item = QtWidgets.QTableWidgetItem(lbl)
                table.setVerticalHeaderItem(position, item)

    def __use_wavelength(self):
        assert self._menu_location >= 0
        # use = not self._use_wavelength[self._menu_location]
        # self._use_wavelength[self._menu_location] = use
        cells = libqt.iter_row(self.tab_data, self._menu_location)
        libqt.cross(cells)

    def _vheader_popup(self, point):
        qmodelindex = self.table_data.indexAt(point)
        row = qmodelindex.row()
        self._menu_location = row

        if row in self.column_range_species:
            menu = self._menuvh_specs
        else:
            crossed = self.table_data.item(row, 1).font().strikeOut()
            self._action_usewl.setChecked(not crossed)
            menu = self._menuvh_wavel

        menu.popup(self.table_data.viewport().mapToGlobal(point))


# class SpecFitWidget(QtWidgets.QWidget):
#     def __init__(self, model):
#         super().__init__()
# 
#         self._use_wavelength = []
#         self._coloured = []
#         self._fixed = []
#         self._menu_location = -1
#         self._model = model
# 
#         self.setObjectName("Spect. fit")
#         self.verticalLayout = QtWidgets.QVBoxLayout(self)
#         self.verticalLayout.setObjectName("verticalLayout")
#         self.tab_data = QtWidgets.QTableWidget(self)
#         self.tab_data.setObjectName("tab_data")
#         # TODO use model to reshape table
#         # self.tab_data.setRowCount(3)
#         n_compounds = model.number_components+model.number_equilibria
#         optically_active = [QtWidgets.QCheckBox("optically active", self) for _ in range(n_compounds)]
#         self.tab_data.setColumnCount(n_compounds)
#         self.tab_data.setRowCount(2)
#         labels = list(libaux.extended_labels(model.labels, model.stoich))
#         self.tab_data.setHorizontalHeaderLabels(labels)
#         for col, checkbox in enumerate(optically_active):
#             self.tab_data.setCellWidget(0, col, checkbox)
# 
#         self.verticalLayout.addWidget(self.tab_data)
# 
#         self._menuvheader = QtWidgets.QMenu(parent=self)
#         self._action_usewl = QtWidgets.QAction('use', parent=self._menuvheader)
#         self._action_usewl.setCheckable(True)
#         self._action_usewl.triggered.connect(self.__use_wavelength)
#         self._action_delwl = QtWidgets.QAction('delete', parent=self._menuvheader)
#         self._menuvheader.addAction(self._action_usewl)
#         self._menuvheader.addAction(self._action_delwl)
# 
#         self._menuhheader = QtWidgets.QMenu(parent=self)
#         self._action_coloured = self._menuhheader.addAction('coloured')
# 
#         hheader = self.tab_data.horizontalHeader()
#         hheader.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
#         hheader.customContextMenuRequested.connect(self._hheader_popup)
# 
#         vheader = self.tab_data.verticalHeader()
#         vheader.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
#         vheader.customContextMenuRequested.connect(self._vheader_popup)
# 
#         model.componentAdded.connect(self.add_component)
#         model.componentDeleted.connect(self.delete_component)
# 
#     def add_component(self, position: int, label: str):
#         '''Add one component to the system.
# 
#         Add one column to model table with a new label.
# 
#         Parameters:
#             label (str): The label for the new reagent
#             position (int): position where the new reagent will be inserted.
#         Returns:
#             int: the position into which the component has been added.
#         '''
#         with libqt.table_locked(self.tab_data) as table:
#             table.insertColumn(position)
#             table.setHorizontalHeaderLabels(self._model.labels)
#             libqt.fill_column(table, position, [0.0]*table.rowCount(), formatting="{:.4f}")
# 
#     def delete_component(self, position: int):
#         with libqt.table_locked(self.tab_data) as table:
#             table.removeColumn(position)
# 
#     def set_wavelengths(self, wavelengths):
#         vertical_header_labels = ["{:.1f}".format(w) for w in wavelengths]
#         with libqt.table_locked(self.tab_data) as table:
#             table.setRowCount(1+len(wavelengths))
#             table.setVerticalHeaderLabels([''] + vertical_header_labels)
#             libqt.array2tab(table,
#                             np.zeros((table.rowCount(), table.columnCount())),
#                             formatting="{:.4f}", row0=1)
#         self._use_wavelength = len(wavelengths) * [True]
# 
#     def _hheader_popup(self, point):
#         qmodelindex = self.tab_data.indexAt(point)
#         col = qmodelindex.column()
#         self._menu_location = col
#         self.__popup_menu(self._menuhheader, point)
# 
#     def _vheader_popup(self, point):
#         qmodelindex = self.tab_data.indexAt(point)
#         row = qmodelindex.row()
#         self._menu_location = row
#         self._action_usewl.setChecked(self._use_wavelength[row])
#         self.__popup_menu(self._menuvheader, point)
# 
#     def __popup_menu(self, menu, point):
#         menu.popup(self.tab_data.viewport().mapToGlobal(point))
# 
#     def __use_wavelength(self):
#         assert self._menu_location >= 0
#         use = not self._use_wavelength[self._menu_location]
#         self._use_wavelength[self._menu_location] = use
#         cells = libqt.iter_row(self.tab_data, self._menu_location)
#         libqt.cross(cells)
