"Main window application."

import collections
import configparser
import itertools
import logging
import sys

import numpy as np

from PyQt5.Qt import QApplication, QClipboard
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5 import QtCore

from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

import csolver
import dialogs
import excepts
import libaux
import libio
import libeq
import libeq.consol
import libemf
import libqt
import canvas_pqtg as canvas
#import canvas
import consts
import ui_mainwindow as mainui
import datawidget
import otherwidgets
from emfwidget import EmfWidget
from modelwidget import ModelWidget
from calorwidget import CalorWidget
from specwidget import SpecWidget
from nmrwidget import NmrWidget
from otherwidgets import OutputWidget, ExternalDataWidget, IonicWidget, TitrationBaseWidget
from simulationwidgets import SpeciationWidget, TitrationWidget
from project import Project
from tabwidget import TabWidget


class MainWindow(QtWidgets.QMainWindow):
    "This class defines the main window."
    def __init__(self):
        super().__init__()

        self.project = Project()
        self._options_widget = dialogs.OptionsDialog()
        self._docprops_widget = dialogs.PropertiesDialog(self.project)

        self.ui = mainui.Ui_MainWindow()
        self.ui.setupUi(self)

        new_tab_main = TabWidget(self)
        self.ui.vl_1.replaceWidget(self.ui.tab_main, new_tab_main)
        self.ui.tab_main.deleteLater()
        self.ui.tab_main = new_tab_main

        self._setup_vars()          # initialize variables
        self._more_widgets()        # create extra widgets
        self._make_connections()    # pair signals/slots

        lbl = QtWidgets.QLabel()
        lbl.setPixmap(QtGui.QPixmap('../icons/idle.png').scaledToHeight(20))
        self.ui.statusbar.insertWidget(0, lbl)

        # FOR TESTING ONLY
        # _debug_fname_ =  '../data/distri_pytc8_extd.xml'
        # _debug_fname_ =  '../data/cuimpy33333_.xml'
        # self.newSpeciation()
        # self.ui.tab_main.clear()
        # _debug_fname_ = '../data/hcit1.xml'
        # _debug_fname_ = '../data/cuimpy33333_.xml'
        # _debug_fname_ = '/home/salvador/Documentos/Trabajo/datos/emf/MA_POTENCIOMETRÍA/hpytrenc8.xml'
        # _debug_fname_ = '/home/salvador/Documentos/Trabajo/datos/emf/citrate/zncit.xml'
        # _debug_fname_ = '/home/salvador/Documentos/Trabajo/datos/emf/Northover/hzn5.xml'
        # _debug_fname_ = '/home/salvador/Documentos/Trabajo/datos/emf/citrate/hcit3.xml'
        # _debug_fname_ = '/home/salvador/Documentos/Trabajo/datos/emf/EDTA/znedta.xml'
        # _debug_fname_ = '/home/salvador/Documentos/Trabajo/datos/emf/EDTA/kk.xml'
        # _debug_fname_ = '/home/salvador/Documentos/Trabajo/datos/emf/colorantes/hnbt.xml'
        # _debug_fname_ = '../data/phosphate.xml'
        # _debug_fname_ = '../data/hdtc.xml'
        # _debug_fname_ = '../data/hpytren4q.xml'
        # _debug_fname_ = '../data/distri_cudtc.xml'
        # _debug_fname_ -= '../data/external.xml'
        # _debug_fname_ = '../data/pylh.xml'
        # _debug_fname_ = '../data/universal_buffer.xml'
        # _debug_fname_ = '/home/salvador/Documentos/Trabajo/documentos/manuscritos/micelas_Mercy/distris/hpytren.xml'
        # logging.debug(f'loading {_debug_fname_}')
        # libio.loadXML(self, _debug_fname_)
        # t1 = self.ui.tab_main.add_titrationbase()
        # t1.set_volume_explicit(False)
        # m = self.ui.tab_main.add_model()
        # m.addEquilibrium(position=-1, stoich=(1, 2), value=10.0, error=0.0)
        # m.addEquilibrium(position=-1, stoich=(1, 3), value=20.0, error=0.0)
        # m.addEquilibrium(position=-1, stoich=(1, 4), value=20.0, error=0.0)
        # self.ui.tab_main.add_titrationbase()
        # self.ui.tab_main.add_spectrumuv()
        # m.removeComponent()
        # self.ui.tab_main.add_nmr()
        # self.ui.tab_main.add_ionic()
        # self.ui.tab_main.import_txtspectrum('./spec1.txt')
        # self.go()
        # self._manual_fitting()
        # self.refresh()
        # libio.importHyperquadApp(self, '/home/salvador/pztrenDoSeTy.hqd')
        # libio.importHyperquadApp(self, '/home/salvador/Documents/Trabajo/datos/emf/pdma/PDMA_0.15_25_080322.HQD')
        # libio.importSuperquadApp(self, '../data/hpytren1.sup')
        # libio.saveXML(self, '../data/hpytren1.xml')
        # libio.loadXML(self, '../data/hpytren1.xml')
        # END TESTING PART

    def go(self):
        'Start calculations.'
        # TODO for python 3.10 uses match case
        if self._is_current_tab((SpeciationWidget, TitrationWidget)):
            try:
                self.ui.tab_main.currentWidget()._recalc_free_concentration()
            except excepts.FailedCalculateConcentrations:
                libqt.popwarning('Error', 'Failed to calculate concentrations')
            except np.linalg.linalg.LinAlgError:
                libqt.popwarning('Error', 'Singular matrix')
            finally:
                self.refresh()
        elif self._is_current_tab(EmfWidget):
            self._fit_emf()
        elif self._is_current_tab(CalorWidget):
            self._fit_calor()
        elif self._is_current_tab(SpecWidget):
            self._fit_spec()
        elif self._is_current_tab(NmrWidget):
            self._fit_nmr()
        elif self._is_current_tab(IonicWidget):
            self._fit_ionic()
        else:
            self.message("Nothing to do")

    def iterate(self):
        "Perform only one iteration."
        raise NotImplementedError

    def menu_about(self):
        aboutd = dialogs.AboutDialog()
        aboutd.exec()

    def menuAddComponent(self):
        "Add a new component."
        text, ok = QtWidgets.QInputDialog.getText(self, "New component",
                                                  "Enter new label")
        if ok:
            pos = self.modelwidget.addComponent(label=text)
            for widget in self.__filtertabs(EmfWidget):
                widget.add_component(text, pos)
            # for widget in self.__filtertabs(SpeciationWidget):
            #     widget.add_component(text)

    def menuConcentrationSolver(self):
        "Open the Concentration Solver dialog."
        d = self._get_tab2dat()
        cs = csolver.CSolver(**d)
        cs.exec()

        if cs.result() == QtWidgets.QDialog.Accepted:
            widget = self.ui.tab_main.currentWidget()
            c = cs.get_full_c()
            if 'x_ref' in d:
                c = np.insert(c, d['x_ref'], d['x_values'], axis=1)
            widget.set_free_concentration(c)
            self.refresh()

    # deprecated
    def menuDeleteDs(self):
        "Delete the current dataset."
        # t = self.ui.tab_main
        # # self.data.delete(t.currentWidget().data)
        # # TODO remove signals/slots
        # t.removeTab(t.currentIndex())
        t = self.ui.tab_main.delete_current_tab()

    def menuModelRename(self):
        """Change the name of the current model.

        A dialog will pop up to introduce the new name.
        """
        t, q = QtWidgets.QInputDialog.getText(self,
                                              "Rename current model",
                                              "enter new name")
        if q:
            self.modelwidget.modelname = t
            action = self._activemodel_ag.checkedAction()
            action.setText(t)

    def menuNewProject(self):
        """Clear all the forms and start a clean project, but ask first.

        Returns True if the user agreed to wipe everything and False if
        the used decided not to do it.

        Returns:
            bool: Whether the used confirmed to clear everything.
        """
        qmsg = QtWidgets.QMessageBox(parent=self)
        qmsg.setText("Do you want to discard the current project?")
        qmsg.setStandardButtons(QtWidgets.QMessageBox.Save |
                                QtWidgets.QMessageBox.Discard |
                                QtWidgets.QMessageBox.Cancel)
        qmsg.setDefaultButton(QtWidgets.QMessageBox.Save)
        ret = qmsg.exec()

        if ret == QtWidgets.QMessageBox.Save:
            self.menuSaveFile()
        elif ret == QtWidgets.QMessageBox.Discard:
            pass
        elif ret == QtWidgets.QMessageBox.Cancel:
            return False
        else:
            raise RuntimeError

        self.ui.tab_main.clear()
        self.newModel()
        self._model_connections(renew=True)
        self._fmode = consts.FM_NONE
        self.__notifymode()
        return True

    def menuSaveFile(self):
        "Save the current file."
        if self._filename is None:
            self.menuSaveFileAs()
        else:
            libio.saveXML(self, self._filename)

    def menuSaveFileAs(self):
        "Ask for a file to save the current data."
        filters = "XML Files (*.xml);;All Files (*.*)"
        filename, ok = self.__savedialog(filters)
        if ok:
            self._filename = filename
            libio.saveXML(self, filename)

    def menuOpenFile(self):
        "Open a file and loads the project."
        if not self.menuNewProject():
            return
        filters = "XML Files (*.xml);;All Files (*.*)"
        filename, ok = self.open_dialog(filters)
        if not ok:
            return
        self._filename = filename

        try:
            libio.loadXML(self, filename)
        except AttributeError as err:
            errtxt = str(err)
            exinfo = sys.exc_info()[2]
            logging.error(f'error opening {filename}\n{errtxt}\n{exinfo}')
            self.message('Data file corrupted or incomplete.')
            import traceback
            print(traceback.format_exc())
        else:
            self._model_connections(renew=True)

    def menu_options(self):
        "Open options dialog."
        self._options_widget.exec()

    def menuPlotColor(self, qaction):
        """Select the color scheme that will be applied to the current plot.

        Several schemes are available. See also :class:`canvas.MyCanvas` and
        :func:`canvas.MyCanvas.setColor`
        """
        action_color = {self.ui.actionColorBW: canvas.MyCanvas.COLORBW,
                        self.ui.actionColorRed:canvas.MyCanvas.COLORR,
                        self.ui.actionColorBlue: canvas.MyCanvas.COLORB,
                        self.ui.actionColorRainbow: canvas.MyCanvas.COLORRB}
        assert qaction in action_color
        self.canvas.setColor(action_color[qaction])

    def menuPlotStyle(self, qaction):
        """Select the style that will be applied to the current plot.

        Several schemes are available. See also :class:`canvas.MyCanvas` and
        :func:`canvas.MyCanvas.setStyle`.
        """
        action_style = {self.ui.actionStyle1: canvas.MyCanvas.STYLE1,
                        self.ui.actionStyle2: canvas.MyCanvas.STYLE2,
                        self.ui.actionStyle3: canvas.MyCanvas.STYLE3,
                        self.ui.actionStyle4: canvas.MyCanvas.STYLE4}
        self.canvas.setStyle(action_style[qaction])

    def menuRemoveComponent(self):
        "Remove a component."
        self.modelwidget.removeComponent()

    @staticmethod
    def menuShowHelp():
        import webbrowser
        webbrowser.open("doc/build/html/index.html")

    def menu_toggle_use(self):
        q = self.ui.actionUse.isChecked()
        self.ui.tab_main.currentWidget().use = q
        self.ui.tab_main.currentWidget().setEnabled(q)

    def message(self, text):
        'Display a message in the status bar.'
        self.statusBar().showMessage(text)
        self.statusBar().repaint()

    # deprecate
    def newCalor(self):
        '''Create new calorimetry data set and return instance.

        Returns:
            :class:`CalorWidget`: The newly created widget.
        '''
        # dsw = CalorWidget(self.modelwidget)
        # # self.modelwidget.refreshWidgets()
        # self._newtab('cal. data', dsw, self.ui.actionNewCalorimetryDS)
        # return dsw
        return self.ui.tab_main.add_calor()

    # deprecate
    def newEmf(self):
        '''Create new potentiometry data set and return instance.

        Returns:
            :class:`EmfWidget`: The newly created widget.
        '''
        # dsw = EmfWidget(self.modelwidget)
        # # dsw.labelsChanged.connect(self._label_changed)
        # self.modelwidget.componentAdded.connect(dsw.add_component)
        # dsw.model = self.modelwidget
        # self._newtab('Emf data', dsw, self.ui.actionNewEmfDS)
        # self.modelwidget.labelsChanged.connect(dsw.updateLabel)
        # return dsw
        return self.ui.tab_main.add_emf()

    def newIonic(self):
        '''Create new ionis strength calculator and return instance.

        Returns:
            :class:`IonicWidget`: The newly created widget.
        '''
        # dsw = IonicWidget()
        # self._newtab('Ionic strength', dsw)
        # return dsw
        return self.ui.tab_main.add_ionic()
 
    # TODO deprecate and use tabmain.add_external_data instead
    def new_external_data(self):
        # dsw = ExternalDataWidget()
        # self._newtab('External data', dsw, self.ui.actionNewExternalData)
        # return dsw
        widget = self.ui.tab_main.add_external_data()
        return widget

    # TODO deprecate and use tabmain.add_model instead
    def newModel(self):
        '''Create new ModelWidget and return a reference to it.

        This function will create a new instance of :class:`ModelWidget` only
        if there is no ModelWidget already created. If there is one, it will
        return a reference to the existing one.

        Returns:
            :class:`ModelWidget`: The newly created widget or the currently
                existing.
        '''
        # if self.modelwidget is not None:
        #     widget = self.modelwidget
        #     widget.clear()
        # else:
        #     widget = ModelWidget()
        #     self._newtab('Model', widget)
        widget = self.ui.tab_main.add_model()
        #self.modelwidget = widget
        return widget

    # TODO deprecate - use self.ui.tab_main.add_speciation()
    def newSpeciation(self):
        # """Create new :class:`SpeciationWidget` and return instance.

        # Returns:
        #     :class:`SpeciationWidget`: the newly created widget.
        # """
        # sdw = SpeciationWidget(self.modelwidget)
        # self._newtab('Speciation', sdw, self.ui.actionNewSpeciesDist)

        # sdw.plotUpdate.connect(self.__plot_update_titr)

        # # connect w.labelsChanged con model.updateLabel
        # # sdw.labelsChanged.connect(self.modelwidget.updateLabel)
        # self.modelwidget.componentAdded.connect(sdw.add_component)
        # self.modelwidget.labelsChanged.connect(sdw.updateLabel)

        # return sdw
        return self.ui.tab_main.add_speciation()

    # deprecate
    def newSpectr(self):
        '''Create a new tab containing a spectrometry data set.

        Returns:
            :class:`SpecDSWidget`: The newly created widget.
        '''
        # if self._spect_fit is None:
        #     dswf = SpecFitWidget(self.modelwidget)
        #     self._newtab('spec fit', dswf, self.ui.actionNewSpecDS)
        # else:
        #     dswf = self._spect_fit

        # dsw = SpecWidget(self.modelwidget)
        # self._newtab('spec data', dsw, self.ui.actionNewSpecDS)
        # return dsw  #, dswf
        return self.ui.tab_main.add_spectrumuv()

    def newNmr(self):
        '''Create new :class:`NmrWidget`.

        Returns:
            :class:`NmrWidget`: the newly created widget.
        '''
        raise NotImplementedError

    def newTitration(self):
        '''Create new :class:`TitrationWidget`.

        Returns:
            :class:`TitrationWidget`: the newly created widget.
        '''
        widget = TitrationWidget(self.modelwidget)
        self._newtab('Titration', widget, self.ui.actionNewTitrSim)

        # connect w.labelsChanged con model.updateLabel
        # w.labelsChanged.connect(self.modelwidget.updateLabel)
        self.modelwidget.componentAdded.connect(widget.add_component)
        self.modelwidget.labelsChanged.connect(widget.updateLabel)
        return widget

    def new_titration_base(self):
        self.ui.tab_main.add_titrationbase()
        # widget = TitrationBaseWidget(self.modelwidget)
        # self._newtab('Titration', widget, self.ui.actionNewTitration)
        # return widget

    def open_dialog(self, filters):
        return QtWidgets.QFileDialog.getOpenFileName(
            parent=self,
            caption='Choose file to open',
            directory=self.project.default_dir,  # self.__default_dir,
            filter=filters)

    def plotEmf(self):
        '''Plot emf data.

        This routine streams emf data and fit data to canvas
        '''
        def data_stream():
            'Stream emf data.'
            # number of potentiometry datasets
            yield self._tabcounter(EmfWidget)
            # plot masked flag
            yield True
            for widget in libqt.filter_tabwidget(self.ui.tab_main, EmfWidget):
                yield (widget.use, widget.titre, widget.emf, widget.emf_fit,
                       widget.C, widget.residuals)

        self.canvas.drawEmfFit(data_stream)

    def refresh(self):
        # widget = self.ui.tab_main.currentWidget()
        if self._is_current_tab((SpeciationWidget, TitrationWidget)):
            widget = self.ui.tab_main.currentWidget()
            if widget.free_concentration() is None:
                return
            self.menuPlotStyle(self.qag_style.checkedAction())
            self.menuPlotColor(self.qag_color.checkedAction())
            plotoptions = {
                'unlabel_lower_than': self._options_widget.unlabel_lower_than,
                'ignore_lower_than': self._options_widget.ignore_lower_than,
                'color_by_group': self.ui.actionColorGroup.isChecked()}
            externaldata = self.__filtertabs(ExternalDataWidget)
            self.canvas.plot_speciation(widget, externaldata, **plotoptions)

        elif self._is_current_tab(EmfWidget):
            self.message('Resolving equilibria')
            self._compute_concentration_combined()
            self.message('Plotting')
            self.canvas.plot_emf(tuple(self.__filtertabs(EmfWidget)))
            self.message('Done')
        elif self._is_current_tab(CalorWidget):
            self.canvas.plotCalorFit()
        elif self._is_current_tab(SpecWidget):
            self.canvas.plotSpectr()
        elif self._is_current_tab(NmrWidget):
            raise NotImplementedError
        elif self._is_current_tab(otherwidgets.ManualFitWidget):
            widget = self.ui.tab_main.currentWidget()
            widget.recalc()
        else:
            # self.message("Nothing to do")
            pass

    def set_mode(self, mode: int):
        if  0 > mode > 5:
            raise ValueError()

        self._fmode = mode
        # TODO clear window and put the new widget
        self.__notifymode()

    def widgets(self):
        yield from self.__itertabs()

    @property
    def modelwidget(self):
        'ModelWidget in TabWidget.'
        widgets = tuple(libqt.filter_tabwidget(self.ui.tab_main, ModelWidget))
        if len(widgets) == 0:
            return None
        if len(widgets) == 1:
            return widgets[0]

        raise RuntimeError('There must be only one ModelWidget.')

    @property
    def temperature(self):
        "The temperature."
        return self.project.temperature  # self._temperature

    # TODO deprecate this. Use Project instead
    @temperature.setter
    def temperature(self, temp):
        if temp <= 0.0:
            raise ValueError('absolute temperature cannot be negative')
        # self._temperature = temp
        self.project.temperature = temp

    # TODO deprecate this. Use Project instead
    @property
    def title(self):
        'The title of the project.'
        # return self._title
        return self.project.title

    # TODO deprecate this. Use Project instead
    @title.setter
    def title(self, title):
        # self._title = title
        self.project.title = title

    # TODO deprecate and move to tabwidget.py
    def _compute_concentration_combined(self):
        # TODO probably move part of this to consol
        # TODO code duplicated in otherwidgets.ManualFitWidget.recalc
        widgets = [t for t in self.__itertabs() if hasattr(t, 'analyticalc')]

        try:
            aux = [np.array(t.analyticalc()) for t in widgets]  # only_analc()]
        except ValueError as e:
            #print(e)
            #print(dir(e))
            errdialog = QtWidgets.QMessageBox.critical(self, "Error", str(e))
            return

        beta = np.array(self.modelwidget.beta)
        stoichiometry = np.array(self.modelwidget.stoich)
        analytic = np.concatenate(aux)
        lims = list(itertools.accumulate([0] + [s.shape[0] for s in aux]))

        if any(w.free_concentration is None for w in widgets):
            concentration_packed = libeq.consol.initial_guess(beta,
                                                              stoichiometry,
                                                              analytic)
        else:
            inivals = np.vstack(tuple(w.free_concentration for w in widgets))
            concentration_packed = libeq.consol.consol(beta, stoichiometry,
                                                       analytic, inivals)
        concentration = np.vsplit(concentration_packed, lims[1:-1])
        for widget, conc in zip(widgets, concentration):
            assert conc.shape[0] == widget.npoints
            widget.free_concentration = conc

    def _copytdat(self):
        # 1. find out number of points
        target_widget = self.ui.tab_main.currentWidget()
        assert hasattr(target_widget, "copy_titration_data")

        # 2. find out the first titration widget
        titration_widgets = tuple(self.__filtertabs(TitrationWidget))
        if not TitrationWidget:
            # TODO warn that a TitrationWidget is needed
            raise RuntimeError
        titration_widget = titration_widgets[0]
        target_widget.copy_titration_data(titration_widget)

    def _create_new_model(self):
        # self.modelwidget.newModel()             # copy current model and append it
        # self.modelwidget.setCurrentModel(-1)    # update with new model
        # TODO update menu
        ...

    def _fit_calor(self):
        "Fit calorimetry data"
        # TODO complete
        raise NotImplementedError()

    def _fit_emf(self):
        "Start fitting for potentiometric data."
        self.__checkoutput()                             # create output widget if it does not exist
        self.ui.tab_main.setCurrentWidget(self.outputw)  # set output at forefront
        kwargs = self.__prefit()

        log_beta = np.array(self.modelwidget.beta_log)
        beta_flags = self.modelwidget.beta_flags
        stoichiometry = np.array(self.modelwidget.stoich)
        titration = []
        electrodes = []
        for emfw in self.__filtertabs(EmfWidget):
            electrodes.append(emfw.electrodes)
            titration.append(emfw.titration())
        args = log_beta, beta_flags, stoichiometry, titration, electrodes
        self.outputw.fitting_header(titration=titration, verbosity=2)  # TODO change 2 for the option set

        try:
            log_beta_out, conc, alt = libemf.emffit(*args, **kwargs)
        except excepts.FailedCalculateConcentrations:
            self.message("Failed to calculate the concentrations")
            self.outputw.appendP("Failed to calculate the concentrations")
        except excepts.TooManyIterations as e:
            self.message("Iteration limit reached.")
            msg = ("Maximum number of iterations reached",
                   "Parameters updated to value in last iteration.")
            for m in msg:
                self.outputw.appendP(m)
            self.modelwidget.beta_raw = e.last_value['last_value']
            # self.modelwidget.beta_update(e.last_value['last_value'])
            for c, tab in zip(e.last_value['concentrations'],
                              self.__filtertabs(EmfWidget)):
                tab.free_concentration = c
            # self.canvas.set_chisq(e.last_value['convergence'])
            fitdata = e.last_value
            self.canvas.plot_emf(tuple(self.__filtertabs(EmfWidget)), fitdata)
        except excepts.NothingToRefine:
            self.message("No parameters to refine. Refinement not done.")
            dialog = QtWidgets.QMessageBox.information(self,
                                                       'Nothing to refine', 
                                                       'There are no free variables to refine.')
        except ValueError as e:
            dialog = QtWidgets.QMessageBox.critical(self, 'Error in Value', str(e))
            self.outputw.appendP("Data is not correct")
            self.outputw.appendP(str(e))
        except Exception as excep:
            self.outputw.appendP("Fitting failed")
            self.message("Fitting failed")
            raise excep
        else:
            self.message("Fitting finished")

            self.outputw.save_last_result(beta=log_beta_out,
                                          beta_error=alt['error_beta'])

            # self.modelwidget.beta_raw = log_beta_out
            # self.modelwidget.beta_error = alt['error_beta']
            for c, tab in zip(conc, self.__filtertabs(EmfWidget)):
                tab.free_concentration = c

            fitdata = alt
            self.canvas.plot_emf(tuple(self.__filtertabs(EmfWidget)), fitdata)

            self.outputw.report_final(log_beta_out, errors=alt['error_beta'],
                                      **alt, verbosity=2)
            self.update()

    def _fit_ionic(self):
        widget = self.ui.tab_main.currentWidget()
        try:
            ionic = np.fromiter(widget.ionic_strength, dtype=np.float)
        except (AttributeError, ValueError):
            self.message("There is an error in the table.")

        titre = np.fromiter(widget.titre, dtype=np.float)
        self.canvas.plot_ionic(titre=titre, ionic=ionic)

    def _fit_nmr(self):
        raise NotImplementedError

    def _fit_spec(self):
        raise NotImplementedError

    def _get_fitcrit(self):
        """Fitting criteria.

        Returns:
            tuple:
        """
        # d = self.config['fitting criteria']
        if self.__fitprecactgr.checkedAction() is self.ui.actionCoarse:
            parms = self._options_widget.fitparams_coarse
        elif self.__fitprecactgr.checkedAction() is self.ui.actionFine:
            parms = self._options_widget.fitparams_fine
        else:
            raise ValueError
        return parms

    def _get_fitalg(self):
        "Return selected fitting algorithm."
        ui = self.ui
        ag = self.__fitmethactgr
        if ag.checkedAction() is ui.actionLevenberg_Marquardt:
            return consts.METHOD_LM
        if ag.checkedAction() is ui.actionNelder_Mead:
            return consts.METHOD_NM

        raise ValueError

    def _get_tab2dat(self):
        "Returns the data chunk of the current tab."
        t = self.ui.tab_main.currentWidget()

        if callable(t.free_concentration):
            c = t.free_concentration()
        else:
            c = t.free_concentration

        if isinstance(t, SpeciationWidget):
            x_values, frz_beta, frz_stoich, frz_analc = t.frozen_params()
            ret = {'beta': frz_beta,
                   'stoichiometry': frz_stoich,
                   'analytc': frz_analc,
                   'x_values': x_values,
                   'x_ref': t.referencex()[0],
                   'free_concentration': t.free_components()}
        elif isinstance(t, TitrationWidget):
            ret = {'beta': t.model.beta,
                   'stoichiometry': t.model.stoich,
                   'analytc': tuple(t.analyticalc()),
                   'free_concentration': c}
        elif isinstance(t, datawidget.DataWidget):
            ret = {'beta': t.model.beta,
                   'stoichiometry': t.model.stoich,
                   'analytc': t.analyticalc(),
                   'free_concentration': c}
        else:
            raise NotImplementedError
        return ret

    def _ignorelowerthan(self, threshold):
        self.canvas.ignorelowerthan = threshold
        # TODO update plot

    def _import_external(self):
        filters = "Text Files (*.txt);;All Files (*.*)"
        filename, ok = self.open_dialog(filters)
        if not ok:
            return
        widget = self.new_external_data()
        widget.load_from_data(filename)

    def _import_k88(self):
        raise NotImplementedError

    def _import_model(self):
        filters = "APES Files (*.xml);;All Files (*.*)"
        f, ok = self.open_dialog(filters)
        if not ok:
            return
        libio.loadModelsOnlyXML(self.modelwidget, f)

    def _import_tiamo(self):
        filters = "TIAMO Files (*.csv);;All Files (*.*)"
        f, ok = self.open_dialog(filters)
        if not ok:
            return
        #w = self.newEmf()
        w = self.ui.tab_main.add_emf()
        volume, emf = libio.import_tiamo(f)
        # always set 'emf' first because the table is reshaped accordingly
        w.emf = emf
        w.titre = volume
        w.labels = self.modelwidget.labels

    def _import_pasat(self):
        filters = "PASAT Files (*.ptr);;All Files (*.*)"
        f, ok = self.open_dialog(filters)
        if not ok:
            return
        w = self.ui.tab_main.add_emf()
        # w = self.newEmf()
        volume, emf = libio.import_pasat(f)
        # always set 'emf' first because the table is reshaped accordingly
        w.emf = emf
        w.titre = volume

    def _import_superquad(self):
        'Import data from SUPERQUAD file.'
        filters = "Superquad Files (*.sup);;All Files (*.*)"
        filename, ok = self.open_dialog(filters)
        if not ok:
            return

        # TESTING
        libio.importSuperquadApp(self, filename)
        self._fmode = consts.FM_EMF

        # PRODUCTION
        #try:
        #    libio.importSuperquadApp(self, filename)
        #except OSError:
        #    self.message('File could not be opened')
        #except Exception as e:
        #    self.message('Error importing data')
        #    print(e)
        #else:
        #    self._fmode = consts.FM_EMF

    def _import_hyperquad(self):
        'Import data from HYPERQUAD file.'
        filters = "Hyperquad Files (*.hqd);;All Files (*.*)"
        filename, ok = self.open_dialog(filters)
        if not ok:
            return

        # TESTING
        libio.importHyperquadApp(self, filename)
        self._fmode = consts.FM_EMF

    def _import_txtsp(self):
        'Import data from text file.'
        filters = "Text Files (*.txt);;All Files (*.*)"
        filename, ok = self.open_dialog(filters)
        if not ok:
            return
        self.ui.tab_main.import_txtsp(filename)

    def _is_current_tab(self, type_):
        "Check whether current tab is of a particular type."
        return isinstance(self.ui.tab_main.currentWidget(), type_)

    def _label_changed(self, index, label):
        # print("label {} is now '{}'".format(index, label))
        raise NotImplementedError

    def _manual_fitting(self):
        manfit_widgets = sum(1 for t in self.__itertabs()
                               if isinstance(t, otherwidgets.ManualFitWidget))
        assert manfit_widgets < 2
        if manfit_widgets == 0:
            model = self.modelwidget.current_model()
            labels = self.modelwidget.labels
            widgets = [t for t in self.__itertabs() if hasattr(t, 'analyticalc')]
            plotemf = lambda: self.canvas.plot_emf(tuple(self.__filtertabs(EmfWidget)))

            dsw = otherwidgets.ManualFitWidget(model, labels, widgets, plotemf)
            self._newtab("Manual fit", dsw)
        else:
            pass

    def _make_connections(self):
        """Connect signals with slots.

        Do additional preparations in the GUI, namely prepare
        QActionGroups and connect signals and slots as appropriate
        """
        def _newxag(*args):
            "Create a QActionGroup with the argument supplied and returns it"
            ag = QtWidgets.QActionGroup(self)
            for a in args:
                ag.addAction(a)
            ag.setExclusive(True)
            return ag

        ui = self.ui
        # ↓↓ action groups ↓↓
        self.__fitmethactgr = _newxag(ui.actionNelder_Mead,
                                      ui.actionLevenberg_Marquardt)
        self.__fitprecactgr = _newxag(ui.actionCoarse,
                                      ui.actionFine)
        self.__wactgr = _newxag(ui.actionWAuto, ui.actionWUnit)

        _action_fonts = (ui.actionFont6, ui.actionFont8, ui.actionFont10,
                         ui.actionFont12)
        self.__fontszactgr = _newxag(*_action_fonts)
        for act, sz in zip(_action_fonts, (8.0, 10.0, 12.0, 14.0)):
            act.setText("{:.1f}pt".format(sz))
            act.triggered.connect(lambda: self.canvas.label_size_changed(sz))

        self.qag_color = _newxag(ui.actionColorBW, ui.actionColorRed,
                                 ui.actionColorBlue,
                                 ui.actionColorRainbow)
        self.qag_color.triggered.connect(self.__plstycol)
        self.qag_style = _newxag(ui.actionStyle1, ui.actionStyle2,
                                 ui.actionStyle3, ui.actionStyle4)
        self.qag_style.triggered.connect(self.__plstych)
        # ↓↓ signals/slots ↓↓  ↑↑ action groups ↑↑
        ui.actionNew.triggered.connect(self.menuNewProject)
        ui.actionOpen.triggered.connect(self.menuOpenFile)
        ui.actionSave.triggered.connect(self.menuSaveFile)
        ui.actionSave_as.triggered.connect(self.menuSaveFileAs)
        ui.actionImportSuperquad.triggered.connect(self._import_superquad)
        ui.actionImportHyperquad.triggered.connect(self._import_hyperquad)
        ui.actionProperties.triggered.connect(self.__docprops)
        ui.actionOptions.triggered.connect(self.menu_options)
        # actionExit

        ui.actionNewModel.triggered.connect(self._create_new_model)
        ui.actionImportModel.triggered.connect(self._import_model)
        ui.actionRename.triggered.connect(self.menuModelRename)
        ui.actionAdd_reagent.triggered.connect(self.menuAddComponent)
        ui.actionNewModel.triggered.connect(self._create_new_model)
        ui.actionResetModel.triggered.connect(self._reset_model)
        self._model_connections()

        ui.actionNewEmfDS.triggered.connect(self.ui.tab_main.add_emf)
        ui.actionNewSpecDS.triggered.connect(self.ui.tab_main.add_spectrumuv)
        ui.actionNewNMRDS.triggered.connect(self.newNmr)
        ui.actionNewCalorimetryDS.triggered.connect(self.newCalor)
        ui.actionNewTitrSim.triggered.connect(self.newTitration)
        ui.actionNewSpeciesDist.triggered.connect(self.newSpeciation)
        ui.actionNewExternalData.triggered.connect(self.new_external_data)
        ui.actionDeleteDS.triggered.connect(self.menuDeleteDs)
        ui.actionNewTitration.triggered.connect(self.new_titration_base)
        ui.actionUse.triggered.connect(self.menu_toggle_use)
        ui.actionImportPASAT.triggered.connect(self._import_pasat)
        ui.actionImportTiamo.triggered.connect(self._import_tiamo)
        ui.actionImportTxtSpec.triggered.connect(self._import_txtsp)
        ui.actionCopyTitDat.triggered.connect(self._copytdat)
        ui.actionImportK88.triggered.connect(self._import_k88)
        ui.actionImportExternal.triggered.connect(self._import_external)

        ui.actionAddittional_data.triggered.connect(self.__pladdd)
        ui.actionClear_external_data.triggered.connect(self.__plclrad)
        ui.actionExportImage.triggered.connect(self.__plexp)
        ui.actionAspectRatio.triggered.connect(self.__figaspect)

        ui.actionConcentrationSolver.triggered.connect(self.menuConcentrationSolver)
        ui.actionSaveConc.triggered.connect(self.__saveconc)
        ui.actionCopyConc.triggered.connect(self.__copyconc)

        ui.actionWAuto.triggered.connect(self.project.set_weight_auto)
        ui.actionWUnit.triggered.connect(self.project.set_weight_unit)
        ui.actionAllowDP.triggered.connect(self.__dangerparams)

        ui.actionIonic_strength_calculator.triggered.connect(self.newIonic)

        ui.actionHelp.triggered.connect(self.menuShowHelp)
        ui.actionAbout_APES.triggered.connect(self.menu_about)

        ui.tab_main.currentChanged.connect(self.__tab_changed)

        ui.pb_go.clicked.connect(self.go)
        ui.pb_refresh.clicked.connect(self.refresh)
        # ui.pb_refresh.clicked.connect(self._menu_refresh)
        ui.pb_iteration.clicked.connect(self.iterate)

    def _model_connections(self, renew=False):
        unique = QtCore.Qt.AutoConnection | QtCore.Qt.UniqueConnection
        u = self.ui
        m = self.modelwidget
        u.actionAdd_beta.triggered.connect(self.__add_equilibrium)
        u.actionRemove_beta.triggered.connect(m.removeEquilibrium)
        u.actionRemove_reagent.triggered.connect(m.removeComponent)
        if not renew:
            u.actionNewModel.triggered.connect(self.__newmodel)
            u.action_ar_calorimetry.triggered.connect(self.__switch_calorimetry)
            u.action_set_temperature.triggered.connect(self.__set_temperature)
            u.action_combine_constants.triggered.connect(self.__combine_constants, type=unique)
            u.action_manual_fitting.triggered.connect(self._manual_fitting)

            self._activemodel_ag.triggered.connect(self.__change_active_model)

    def _menu_refresh(self):
        def _what(type_):
            return isinstance(self.ui.tab_main.currentWidget(), type_)

        # TODO for python 3.10 use match case
        if _what(EmfWidget):
            self._compute_concentration_combined()
            # for widget in self.__filtertabs(EmfWidget):
            #     widget.calc_freec()
            self.canvas.plot_emf(tuple(self.__filtertabs(EmfWidget)))
            self.update()
        elif _what(CalorWidget):
            raise NotImplementedError
        elif _what(SpecWidget):
            raise NotImplementedError
        elif _what(NmrWidget):
            raise NotImplementedError
        else:
            self.message("Nothing to do")

    def _more_widgets(self):
        "Create extra widgets."

        self.main_frame = QtWidgets.QWidget()
        self.canvas = canvas.MyCanvas(self.ui.widget_plot, width=8, height=6, dpi=100)
        # self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

        layout = QtWidgets.QVBoxLayout(self.ui.widget_plot)
        layout.addWidget(self.canvas)
        #layout.addWidget(self.mpl_toolbar)

        #self._newtab('Model', ModelWidget())
        self.ui.tab_main.add_model()

    # TODO deprecate
    def _newtab(self, txt, dsw, action=None):
        "Common tasks for when a new tab is opened."
        n_tabs = self._tabcounter(type(SpeciationWidget))
        if n_tabs == 0:
            append = ''
        else:
            append = str(n_tabs+1)
        i = self.ui.tab_main.addTab(dsw, txt + append)
        self.ui.tab_main.setCurrentIndex(i)

    def _reset_model(self):
        dialog = dialogs.ModelInputDialog(self)
        dialog.exec()
        if dialog.result() == QtWidgets.QDialog.Accepted:
            n_species = dialog.ui.sb_components.value()
            n_equils = dialog.ui.sb_equilibria.value()
            self.modelwidget.clear(n_equils, n_species)

    def _setup_vars(self):
        "initialize variables"
        self._filename = None
        self.outputw = None
        self._multilock = True  # if true, forbid mixing different data types
        self._fmode = consts.FM_NONE
        self._activemodel_ag = QtWidgets.QActionGroup(self)
        self._activemodel_ag.setExclusive(True)
        self._available_models = [self.ui.actionModel_1]
        self._spect_fit = None

    def _tabcounter(self, ctype=None):
        counter = collections.Counter(libqt.iter_tabwidget(self.ui.tab_main))
        if ctype is None:
            return counter

        return counter[ctype]

    def __add_equilibrium(self):
        self.modelwidget.addEquilibrium()

    def _update_available_models(self):
        self._available_models = [QtWidgets.QAction(name, self._activemodel_ag)
                                  for name in (model.name for model in self.modelwidget)]
        self.ui.menuSet_active.clear()
        for action in self._available_models:
            action.setCheckable(True)
            self.ui.menuSet_active.addAction(action)

    def __checkoutput(self):
        if self.outputw is not None:
            return
        self.outputw = OutputWidget()
        self.outputw.ui.pb_update.clicked.connect(self.__update_constants)
        ui = self.ui
        i = ui.tab_main.addTab(self.outputw, "Output")
        ui.menuOutput.setEnabled(True)
        ui.actionOutputClear.triggered.connect(self.outputw.clear)
        # self.ui.tab_main.setCurrentIndex(i)

    def __combine_constants(self):
        beta = self.modelwidget.beta_raw
        error_beta = self.modelwidget.beta_error
        labels = libaux.constant_labels(self.modelwidget.labels,
                                        self.modelwidget.stoich,
                                        kind='unicode')
        cts = dialogs.Constants(beta, error_beta, tuple(labels))
        cts.exec()

    def __configread(self):
        "Reads the config file and loads the contents"
        # TODO complete
        config = configparser.ConfigParser()
        config.read('config.ini')

    def __configwrite(self):
        "Writes the config file"
        # TODO complete
        pass

    def __dangerparams(self):
        are_allowed = self.ui.actionAllowDP.isChecked()
        for widget in self.__filtertabs(EmfWidget):
            widget.allow_dangerous(are_allowed)

    def __docprops(self):
        self._docprops_widget.exec()

    def __figaspect(self):
        asp, ok = QtWidgets.QInputDialog.getDouble(self, 
                                                   'Aspect ratio',
                                                   'select width/height ratio', 
                                                   decimals=4, min=0.1, max=10.0, value=1.33)
        if ok:
            ...

    def __newmodel(self):
        'Create new model and change to it.'
        self.modelwidget.newModel()
        self._update_available_models()
        action = self._available_models[-1]
        self.__change_active_model(action)
        # self.modelwidget.setCurrentModel(-1)

    def __notifymode(self):
        modes = {consts.FM_NONE: 'idle', consts.FM_EMF: 'potentiometry',
                 consts.FM_CALOR: 'calorimetry',
                 consts.FM_SPEC: 'spectrometry', consts.FM_NMR: 'NMR'}
        txt = "APES is now in mode %s" % modes[self._fmode]
        self.message(txt)

    def __plot_update_titr(self):
        # print("plot needs update")
        self.refresh()

    # TODO DELETE AND REPLACE WITH otherwidgets.externaldatawidget
    def __pladdd(self):
        "Add additional data to plot"

        # filters = "Text Files (*.txt);;All Files (*.*)"
        # filename, ok = self.__opendialog(filters)
        # if not ok:
        #     return

        # # TODO
        # with open(filename, 'r') as f:
        #     axttl = f.readline().strip()
        #     lbls = f.readline().split()
        #     extdat = np.loadtxt(f)

        # widget = self.ui.tab_main.currentWidget()
        # widget.set_externalDataTitle(axttl)
        # widget.setExternalData(lbls, extdat)
        # self.message("data successfully loaded")
        pass

    def __plclrad(self):
        d = self._get_tab2dat()
        # assert isinstance(d, data.SimulationData)
        d.clearExternalData()

    def __plstycol(self, style):
        "Update plot when color scheme changes."
        self.menuPlotColor(style)
        self.__plch()

    def __plstych(self, style):
        "Update plot when style changes."
        self.menuPlotStyle(style)
        self.__plch()

    # TODO deprecate this
    def __plch(self):
        "Redraw speciation plot when something changes."
        # sdw = self.ui.tab_main.currentWidget()
        # plargs = sdw.getPlot()
        # pldat = self.data.plotSpeciation(sdw.data,
        #                                  plargs.pop('ref'),
        #                                  plargs.pop('pX'))
        # self.canvas.plotSpeciation(pldat, **plargs)
        self.refresh()

    def __prefit(self):
        method = self._get_fitalg()
        kwargs = {'method': method}
        crit = self._get_fitcrit()
        if method == consts.METHOD_LM:
            kwargs['threshold'] = crit[0]
        elif method == consts.METHOD_NM:
            kwargs.update({'term_x': crit[0], 'term_f': crit[1]})
        else:
            raise RuntimeError
        kwargs['max_iterations'] = crit[2]
        self.outputw.set_report(method)
        kwargs['report'] = self.outputw.report_iteration
        return kwargs

    def __filtertabs(self, tabtype):
        'Filter tab widgets of a particular type.'
        for widget in self.__itertabs():
            if isinstance(widget, tabtype):
                yield widget

    def __itertabs(self):
        raise DeprecationWarning
        'Iterate over widgets in the tab widget.'
        for tabn in range(self.ui.tab_main.count()):
            yield self.ui.tab_main.widget(tabn)

    def __savedialog(self, filters):
        return QtWidgets.QFileDialog.getSaveFileName(
            parent=self,
            caption='Choose file to open',
            directory=self.project.default_dir,  #self.__default_dir,
            filter=filters)

    def __savefile(self):
        filters = "XML Files (*.xml);;All Files (*.*)"
        filename = self.__savedialog(filters)
        if not filename:
            return
        libio.saveXML(self, filename)

    def __concentration_matrix(self):
        current_widget = self.ui.tab_main.currentWidget()
        fail_msg = 'No concentration values available'
        if not hasattr(current_widget, 'free_concentration'):
            self.message(fail_msg)
            return None
        conc = current_widget.free_concentration()
        if conc is None:
            self.message(fail_msg)
            return None

        assert isinstance(conc, np.ndarray)
        return conc

    def __copyconc(self):
        "Copy concentration matrix to clipboard."
        conc = self.__concentration_matrix()
        labels = libaux.extended_labels(self.modelwidget.labels,
                                        self.modelwidget.stoich)
        text = "\t".join(labels) + "\n" + \
            "\n".join("\t".join(f"{num:.6e}" for num in line) for line in conc)
        clipboard = QApplication.clipboard()
        clipboard.setText(text)

    def __saveconc(self):
        "Save free concentrations to a file."
        current_widget = self.ui.tab_main.currentWidget()
        # fail_msg = 'No concentration values available'
        # if not hasattr(current_widget, 'free_concentration'):
        #     self.message(fail_msg)
        #     return
        # conc = current_widget.free_concentration()
        # if conc is None:
        #     self.message(fail_msg)
        #     return
        # assert isinstance(conc, np.ndarray)
        conc = self.__concentration_matrix()

        filters = ("Numpy files (*.npy)",
                   "Numpy compressed (*.npz)",
                   "Comma separated values (*.csv)")
        filename, type_ = self.__savedialog(';;'.join(filters))
        if not filename:
            return

        if type_ == filters[0]:
            np.save(filename, conc)
        elif type_ == filters[1]:
            er_conc = current_widget.free_concentration_error()
            labels = current_widget.extended_labels()
            beta = np.array(self.modelwidget.beta)
            analc = current_widget.analyticalc()
            stoich = np.array(self.modelwidget.stoich)
            np.savez(filename, free_concentration=conc,
                     free_concentration_error=er_conc,
                     beta=beta, analyticalc=analc,
                     stoichiometry=stoich, labels=labels)
        if type_ == filters[2]:
            np.savetxt(filename, conc)
        else:
            self.message('Unknown file type')

    def __change_active_model(self, action):
        action.setChecked(True)
        idx = self._available_models.index(action)
        self.modelwidget.save_current_model()
        self.modelwidget.setCurrentModel(idx)

    def __set_temperature(self):
        model = self.modelwidget
        get_double = QtWidgets.QInputDialog.getDouble
        temp, accept = get_double(self, 'Set Temperature', 'Temperature (K)',
                                  value=model.temperature, min=0.0, decimals=2)
        if accept:
            model.temperature = temp

    def __switch_calorimetry(self):
        if self.modelwidget.isCalori():
            self.modelwidget.removeCalorimetry()
        else:
            self.modelwidget.addCalorimetry()

    def __tab_changed(self, i):
        "Triggered when the tab changes."
        widget = self.ui.tab_main.widget(i)
        # if 'data' in dir(ctw) and 'use' in dir(ctw.data):
        if isinstance(widget, datawidget.DataWidget):
            self.ui.actionUse.setEnabled(True)
            self.ui.actionUse.setChecked(widget.use)
        else:
            self.ui.actionUse.setEnabled(False)

        buttons = (self.ui.pb_go, self.ui.pb_iteration, self.ui.pb_undo, self.ui.pb_refresh)
        actions1 = (self.ui.actionUndo, self.ui.actionRefresh,
                    self.ui.actionIteration, self.ui.actionFit)
        all_items = (*buttons, *actions1, self.ui.actionCopyTitDat)
        infodict = {
            EmfWidget: (*buttons, *actions1, self.ui.actionCopyTitDat),
            SpecWidget: (*buttons, *actions1, self.ui.actionCopyTitDat),
            TitrationWidget: (*buttons, *actions1),
            SpeciationWidget: (*buttons, *actions1),
            IonicWidget: (self.ui.pb_go,),
            otherwidgets.ManualFitWidget: (self.ui.pb_refresh,),
        }

        items = infodict.get(type(widget), [])
        for item in all_items:
            item.setEnabled(item in items)

    def __unlabellowerthan(self, threshold):
        self.canvas.unlabellowerthan = threshold
        # TODO raise plot_update

    def __update_constants(self):
        lastresult = self.outputw.get_last_result()
        self.modelwidget.beta_raw = lastresult['beta']
        self.modelwidget.beta_error = lastresult['beta_error']

    def __plexp(self):
        "export plot"
        filters = ";;".join(("Portable Network Graphics (*.png)",
                             "Portable Document Format (*.pdf)",
                             "Postscript (*.ps)",
                             "Encapsulated Postscript (*.eps)",
                             "Scalable Vector Graphics (*.svg)",
                             "All Files (*.*)"))
        filename = self.__savedialog(filters)
        if filename:
            self.canvas.saveFigure(filename)
