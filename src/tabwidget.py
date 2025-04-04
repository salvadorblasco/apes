"""
Main module containing the TabWidget class, which organizes and manages tab-based widgets 
for various scientific or analytical tasks within the main application window.

Dependencies:
    - PyQt5.QtCore
    - PyQt5.QtGui
    - PyQt5.QtWidgets
    - Various custom widgets and libraries (e.g., bridge, libfit, libmath)
"""

import datetime
from typing import Callable, Any

from PyQt5 import QtCore, QtGui, QtWidgets

import bridge
import consts
import libfit
import libmath
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
    """A QTabWidget subclass that organizes and manages different types of widgets (e.g., 
    analytical and simulation widgets) in the main application.

    Attributes:
        model (ModelWidget or None): The current model widget.
        output (OutputWidget or None): The output widget for displaying logs or reports.
        _specmodel (SpecModelWidget or None): The current spectral model widget.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        super().setMovable(True)
        self.model = None
        self.output = None
        self._specmodel = None

    # def add_external_data(self):
    #     widget = ExternalDataWidget()
    #     self.__stamp(widget, "External")
    #     return widget

    def add_calor(self):
        "Add a calorimetric widget to the tab."
        widget = CalorWidget(self.model)
        self.__stamp(widget, "Calor")
        return widget

    def add_emf(self):
        "Add an EMF (Electromotive Force) widget to the tab."
        widget = EmfWidget(model=self.model, parent=self)
        self.__stamp(widget, "EMF")
        return widget

    # def add_ionic(self):
    #     widget = IonicWidget()
    #     self.__stamp(widget, "Ionic Strength")
    #     return widget

    def add_nmr(self):
        "Add an NMR (Nuclear Magnetic Resonance) widget to the tab."
        widget = NmrWidget(self.model)
        self.__stamp(widget, "NMR")
        return widget

    # def addmodel(self):
    #     if self.model is None:
    #         widget = ModelWidget()
    #         self.addTab(widget, "Model")
    #         self.model = widget
    #     else:
    #         widget = self.model
    #     return widget

    # def add_speciation(self):
    #     widget = SpeciationWidget(self.model)
    #     self.__stamp(widget, "Speciation")
    #     return widget

    def add_spectrumuv(self):
        widget = SpecWidget(self.model)
        self.__stamp(widget, "UV-vis")
        return widget

    # def add_titration(self):
    #     widget = TitrationWidget(self.model)
    #     self.__stamp(widget, "Titration Simulation")
    #     return widget

    def add_titration(self):
        "Add a titration widget to the tab."
        widget = TitrationBaseWidget(self.model)
        self.__stamp(widget, "Titration")
        self._update_titration_list()
        widget.implicitVolumeChanged.connect(self.__implicit_titre_changed)
        return widget

    # def add_uvmodel(self):
    #     if self._uvmodel is None:
    #         widget = SpecModelModelWidget()
    #         self.addTab(widget, "UV Components")
    #         self._uvmodel = widget
    #     else:
    #         widget = self._uvmodel
    #     return widget

    def clear(self):
        "Clear all tabs and reset state."
        self.model = None
        super().clear()

    def delete_current_tab(self):
        "Delete the currently selected tab."
        match self.currentWidget():
            case ModelWidget():
                pass  # Prevent deletion of the model widget. 
            case TitrationBaseWidget():
                self.__delete_current_tab()
                self._update_titration_list()
            case DataWidget():
                self.__delete_current_tab()

    def get_titration_names(self):
        for tab in self._itertabs():
            if isinstance(tab, TitrationBaseWidget):
                yield tab.name

    def find_titration_byname(self, name):
        for tab in self._itertabs():
            if isinstance(tab, TitrationBaseWidget) and tab.name == name:
                return tab
        raise KeyError(f"object {name} does not exist")

    def fitting(self, option: Callable[[str], Any]) -> None:
        """Perform data fitting using the specified fitting algorithm.

        Args:
            option (Callable): A callback function to fetch fitting options.
        """
        if self.output is None:
            self.output = OutputWidget()
            self.addTab(self.output, "Output")

        start_time = datetime.datetime.now().isoformat(timespec='seconds')
        self.setCurrentWidget(self.output)
        self.output.write(f"starting fitting at {start_time} \n")

        method = option('fitting algorithm')
        if method not in (consts.METHOD_NM, consts.METHOD_LM):
            raise ValueError("Method not known")

        self.output.refresh()

        titwidgets = [w for w in self._itertabs() if isinstance(w, TitrationBaseWidget)]
        datwidgets = [w for w in self._itertabs() if isinstance(w, DataWidget)]
        params = bridge.Parameters(self.model, titwidgets, datwidgets)
        bridgeobj = bridge.Bridge(params, report_buffer=self.output.buffer)

        ffit = libfit.fitting_functions[method]
        
        info = ffit(bridgeobj, weight=bridgeobj.weights(), report=self.output.buffer, debug=False)

        covariance = libmath.covariance(info['jacobian'], bridgeobj.weights())
        errors = libmath.fitting_errors(covariance)

        self.output.refresh()
        # breakpoint()
        # params.accept_values()
        params.set_errors(errors)
        params.dump_to_widgets()

        return bridgeobj

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

    def plot(self, canvas, bridge):
        widgets_to_plot = [w for w in self._itertabs() if isinstance(w, DataWidget)]
        canvas.plot_fitting(widgets_to_plot, bridge)

    def widgets_to_save(self):
        types = (EmfWidget, NmrWidget, SpecWidget, CalorWidget, TitrationBaseWidget)
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

    def _refresh_canvas_simulate(self, canvas, option):
        widget = super().currentWidget()
        if widget.free_concentration() is None:
            return
        canvas.setStyle(option('plot style'))
        canvas.setColor(option('plot color'))
        plotoptions = {k: option(k) for k in ('unlabel_lower_than', 'ignore_lower_than',
                                              'color_by_group')}
        # externaldata = self.__filtertabs(ExternalDataWidget)
        externaldata = None
        canvas.plot_speciation(widget, externaldata, **plotoptions)

    def _update_titration_list(self):
        for widget in self._itertabs():
            tnames = tuple(self.get_titration_names())
            if hasattr(widget, 'populate_cb_titration'):
                widget.populate_cb_titration(tnames)

    def __delete_current_tab(self):
        "Delete the currently selected tab by index."
        idx = self.currentIndex()
        self.removeTab(idx)

    # def __getitem__(self, key):
    #     for tab in self._itertabs():
    #         if tab.tabText() == key:
    #             return tab.self.__tabdicts[key]

    def __stamp(self, widget, group):
        """Add a widget to the tab and assign it a unique name.

        Args:
            widget (QtWidgets.QWidget): The widget to add.
            group (str): The group name for the widget.
        """
        name = group + str(self.count())
        widget.name = name
        if hasattr(widget, 'populate_cb_titration'):
            tnames = tuple(self.get_titration_names())
            widget.titrationChanged.connect(self.__titration_changed)
            widget.populate_cb_titration(tnames)
        self.addTab(widget, name)

    @QtCore.pyqtSlot(object, str)
    def __titration_changed(self, widget, txt):
        # print("slot activated: ", widget, txt)
        # titration = self.__tabdicts[txt]
        if widget.titration.is_titre_implicit():
            titre = titration.titre
            widget.npoints = titration.n_points
            widget.titre = titre
        widget._titrationid = id(widget.titration)

    @QtCore.pyqtSlot()
    def __implicit_titre_changed(self):
        sender = self.sender()
        for widget in self._itertabs():
            if isinstance(widget, DataWidget):
                widget.npoints = sender.n_points
                widget.titre = sender.titre
