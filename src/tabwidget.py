"""Main QTabWidget to organize tabs in mainwindow
"""

from PyQt5 import QtCore, QtGui, QtWidgets

import libio
from calorwidget import CalorWidget
from emfwidget import EmfWidget
from modelwidget import ModelWidget
from nmrwidget import NmrWidget
from otherwidgets import OutputWidget, ExternalDataWidget, IonicWidget, TitrationBaseWidget
from simulationwidgets import SpeciationWidget, TitrationWidget
from specwidget import SpecWidget
from datawidget import DataWidget


class TabWidget(QtWidgets.QTabWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self._model = None
        self.__tabdicts = {}
        self.__tabcounter = {}

    def add_external_data(self):
        widget = ExternalDataWidget()
        self.__stamp(widget, "External")
        return widget

    def add_calor(self):
        widget = CalorWidget(self._model)
        self.__stamp(widget, "Calor")
        return widget

    def add_emf(self):
        widget = EmfWidget(self._model)
        self.__stamp(widget, "EMF")
        return widget

    def add_ionic(self):
        widget = IonicWidget()
        self.__stamp(widget, "Ionic Strength")
        return widget

    def add_nmr(self):
        widget = NmrWidget(self._model)
        self.__stamp(widget, "NMR")
        return widget

    def add_model(self):
        if self._model is None:
            widget = ModelWidget()
            self.addTab(widget, "Model")
            self._model = widget
        else:
            widget = self._model
        return widget

    def add_speciation(self):
        widget = SpeciationWidget(self._model)
        self.__stamp(widget, "Speciation")
        return widget

    def add_spectrumuv(self):
        widget = SpecWidget(self._model)
        self.__stamp(widget, "UV-vis")
        return widget

    def add_titration(self):
        widget = TitrationWidget(self._model)
        self.__stamp(widget, "Titration Simulation")
        return widget

    def add_titrationbase(self):
        widget = TitrationBaseWidget(self._model)
        self.__stamp(widget, "Titration")
        self._update_titration_list()
        widget.implicitVolumeChanged.connect(self.__implicit_titre_changed)
        return widget

    def clear(self):
        "Clear all the tabs and delete everything."
        self._model = None
        self.__tabdicts = {}
        super().clear()

    def delete_current_tab(self):
        match self.currentWidget():
            case ModelWidget():
                pass
            case TitrationBaseWidget():
                self.__delete_current_tab()
                self._update_titration_list()
            case DataWidget():
                self.__delete_current_tab()

    def get_titration_names(self):
        for tab in self._itertabs():
            if isinstance(tab, TitrationBaseWidget):
                yield tab.name

    def import_txtspectrum(self, filename):
        'Import data from text file.'
        wavelength, data = libio.import_spectrum_text(filename)
        titrwidget = self.add_titrationbase()
        titrwidget.n_points = data.shape[1]
        spectwidget = self.add_spectrumuv()
        tcombo = spectwidget.ui.cb_titration
        index = tcombo.findText(titrwidget.name)
        tcombo.setCurrentIndex(index)
        spectwidget.feed_data(wavelengths=wavelength, spectra_data=data)

    def widgets_to_save(self):
        types = (ModelWidget, SpeciationWidget,
                 TitrationWidget,
                 EmfWidget, NmrWidget,
                 SpecWidget, CalorWidget,
                 ExternalDataWidget,
                 TitrationBaseWidget)
        for widget in self._itertabs():
            if isinstance(widget, types):
                yield widget

    def _freec_widgets(self):
        'Yield widgets which can calculate free concentrations.'
        for tab in self._itertabs():
            if hasattr(tab, 'analyticalc'):
                yield tab

    def _itertabs(self):
        'Iterate over widgets in the tab widget.'
        for tabn in range(self.count()):
            yield self.widget(tabn)

    def _update_titration_list(self):
        for widget in self._itertabs():
            tnames = tuple(self.get_titration_names())
            if hasattr(widget, 'populate_cb_titration'):
                widget.populate_cb_titration(tnames)

    def __delete_current_tab(self):
        idx = self.currentIndex()
        self.removeTab(idx)

    def __getitem__(self, key):
        return self.__tabdicts[key]

    def __stamp(self, widget, group):
        counter = self.__tabcounter.get(group, 0) + 1
        name = group + " " + str(counter)
        widget.name = name
        if hasattr(widget, 'populate_cb_titration'):
            tnames = tuple(self.get_titration_names())
            widget.titrationChanged.connect(self.__titration_changed)
            widget.populate_cb_titration(tnames)
        self.__tabdicts[name] = widget
        self.addTab(widget, name)
        self.__tabcounter[group] = counter

    @QtCore.pyqtSlot(object, str)
    def __titration_changed(self, widget, txt):
        # print("slot activated: ", widget, txt)
        titration = self.__tabdicts[txt]
        if titration.is_titre_implicit():
            titre = titration.titre
            widget.npoints = titration.n_points
            widget.titre = titre

    @QtCore.pyqtSlot()
    def __implicit_titre_changed(self):
        sender = self.sender()
        for widget in self._itertabs():
            if isinstance(widget, DataWidget):
                widget.npoints = sender.n_points
                widget.titre = sender.titre
