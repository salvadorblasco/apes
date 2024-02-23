#!/usr/bin/python3

"""
Solving equilibria.

.. module:: csolver.py
.. moduleauthor:: Salvador Blasco <salvador.blasco@gmail.com>
"""

import numpy as np

from PyQt5 import QtWidgets

from matplotlib.backends.backend_qt5agg import \
    NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import \
    FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import ui_csolver
from libeq.pcf import pcf
from libeq.fobj import fobj
from libeq.cpluss import cpluss
from libeq.damping import damping as VTdamping, _damping1 as VTdamping1
from libeq.damping import _clims, _outliers
from libeq.nr import NewtonRaphson
import libeq.consol


class CSolver(QtWidgets.QDialog):
    """This class defines the main window."""

    def __init__(self, beta, analytc, stoichiometry, free_concentration=None,
                 **kwargs):
        super().__init__()
        # print(beta)
        self.ui = ui_csolver.Ui_CSolverWidget()
        self.ui.setupUi(self)
        self.startingvals = None
        self.C = None
        self.T = np.array(analytc)
        self.B = np.array(beta)
        self.P = np.array(stoichiometry)

        n_points, n_components = self.T.shape
        self.n_points = n_points
        self.n_components = n_components
        n_equil = self.P.shape[0]
        self.n_equil = n_equil

        if free_concentration is not None:
            if free_concentration.shape[1] == n_equil + n_components:
                self.startingvals = np.copy(free_concentration[:, :n_components])
            else:
                self.startingvals = np.copy(free_concentration)
        else:
            self.ui.pbInitialVals.setEnabled(False)

        # print(beta)
        self.F = None
        self.prevC = None

        self.chisquared = []
        if self.C is not None:
            self.__dochisq()

        self.ui.sb_dampspecies.setMaximum(n_components-1)
        self.ui.sbRefFrom.setMaximum(n_points)
        self.ui.sbRefFrom.setValue(1)
        self.ui.sbRefTo.setMaximum(n_points)
        self.ui.sbRefTo.setValue(n_points)

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self)

        self.canvas_layout = QtWidgets.QVBoxLayout(self.ui.canvas)
        self.canvas_layout.addWidget(self.canvas)
        self.canvas_layout.addWidget(self.mpl_toolbar)

        self.ui.pbInitialVals.clicked.connect(self.setC0iniv)
        self.ui.pbStartAnalC.clicked.connect(self.setC0Anal)
        self.ui.pbStartFV.clicked.connect(self.set_fixed_value)
        self.ui.pbInitialGuess.clicked.connect(self.__initial_guess)
        self.ui.pbDamping.clicked.connect(self.do_damping)
        self.ui.pbDamp1.clicked.connect(self.do_damping1)
        self.ui.pb_pcf.clicked.connect(self.do_pcf)

        self.ui.pbIter1.clicked.connect(lambda: self._doiter(1))
        self.ui.pbIter10.clicked.connect(lambda: self._doiter(10))
        self.ui.pbIter100.clicked.connect(lambda: self._doiter(100))
        self.ui.pbIter1000.clicked.connect(lambda: self._doiter(1000))

        self.ui.pbSmooth.clicked.connect(self.smooth)
        self.ui.pb_zeronans.clicked.connect(self.zeronans)
        self.ui.pbSave.clicked.connect(self.save)
        self.ui.pbInterpolate.clicked.connect(self.interpolate)
        self.ui.pbFortran.clicked.connect(self.runfortran)

        self.axes_conc = self.figure.add_subplot(211)
        self.axes_concf = self.axes_conc.twinx()
        self.axes_fit = self.figure.add_subplot(212)
        self._cplotline = None
        self._chiplotline = None
        self.axes_conc.set_ylabel('concentration')

        # self.__test_setup()

    def getFitParms(self):
        pass

    def get_full_c(self):
        return cpluss(self.C, self.B, self.P, full=True)

    def interpolate(self):
        masked_c = np.ma.masked_invalid(self.C)
        libeq.consol.interpolate_masked(masked_c)
        self.C = masked_c
        self.redraw()

    def save(self):
        np.savez_compressed('save_csolver.npz', free_concentration=self.C,
                            beta=self.B, stoichiometry=self.P, analytc=self.T)

    def setC0iniv(self):
        self.C = self.startingvals.copy()
        self.reset_chisq()
        # self.__dochisq()
        # print(self.C)
        self.redraw()

    def setC0Anal(self):
        self.C = self.T.copy()
        self.C[self.C <= 0.0] = 1e-10
        self.reset_chisq()
        # self.__dochisq()
        # print(self.C)
        self.redraw()

    def set_fixed_value(self):
        number = QtWidgets.QInputDialog.getText(self, 'Enter value',
                                                'concentration')
        if number[1]:
            value = float(number[0])
            if self.C is None:
                self.C = np.empty_like(self.T)
            self.C[:] = value
            self.reset_chisq()
            # self.__dochisq()
            self.redraw()

    def do_damping(self):
        damping_factor = self.ui.dsb_dampfactor.value()
        self.prevC = self.C.copy()
        VTdamping(self.C, self.B, self.P, self.T,
                  damping_factor=damping_factor)
        self.__dochisq()
        self.redraw()

    def do_damping1(self):
        damping_factor = self.ui.dsb_dampfactor.value()
        species = self.ui.sb_dampspecies.value()

        Tcalc = VTdamping1(self.C, species, self.B, self.P, self.T,
                           damping_factor=damping_factor)

        CLIM1, CLIM2 = _clims(self.T, damping_factor)
        nout, q = _outliers(CLIM1, CLIM2, Tcalc)

        self.__dochisq()
        self.redraw()

    def do_pcf(self):
        try:
            tol = float(self.ui.le_pcf_tol.text())
        except ValueError:
            qmsg = QtWidgets.QMessageBox(parent=self)
            qmsg.setText("Tolerance value is not a number")
            qmsg.setStandardButtons(QtWidgets.QMessageBox.Ok)
            ret = qmsg.exec()
        else:
            c = pcf(self.C, self.B, self.P, self.T, tolerance=tol)
            self.C = c
            self.__dochisq()
            self.redraw()
    
    def doFixOL(self):
        "Fix outliers."
        raise NotImplementedError

    def redraw(self):
        "Update plots."
        if self.C is None:
            return

        self.axes_concf.clear()
        self.axes_conc.clear()

        if self.ui.rb_plotall.isChecked():
            conc_values = cpluss(self.C, self.B, self.P, full=True)
        else:
            conc_values = self.C

        self.axes_conc.plot(conc_values, lw=2)
        self.axes_concf.bar(np.arange(self.F.shape[0]),
                            np.max(np.abs(self.F), axis=1),
                            width=1.0, alpha=0.2, lw=0)

        if self.chisquared:
            self.axes_fit.clear()
            self.axes_fit.semilogy(self.chisquared, 'o-')

        self.axes_conc.set_ylabel('concentration')
        self.axes_concf.set_ylabel('residual')

        self.canvas.draw()
        self.canvas.flush_events()

    def reset_chisq(self):
        self.chisquared = []
        self.__dochisq()

    def runfortran(self):
        retc = libeq.consol._consol_fortran(self.B, self.P, self.T, self.C)
        self.C[:,:] = retc[:, :self.n_components]
        self.__dochisq()
        self.redraw()

    def smooth(self):
        "Smooth curves."
        import scipy.signal
        n_components = self.T.shape[1]
        for component in range(n_components):
            self.C[:, component] = scipy.signal.medfilt(self.C[:, component])
        self.__dochisq()
        self.redraw()

    def zeronans(self):
        "Replace NaN with zeroes."
        self.C[np.isnan(self.C)] = 1e-8
        assert not np.any(np.isnan(self.C))

    def _doiter(self, niterations):
        kwargs = {'damping': self.ui.cbDamping.isChecked(),
                  'forcer': self.ui.cbForcer.isChecked(),
                  'scaling': self.ui.cbScaling.isChecked(),
                  'step_limiter': self.ui.cbStepLimiter.isChecked(),
                  'zero_offdiag': self.ui.cb0OffD.isChecked(),
                  'logc': self.ui.cbLogFit.isChecked(),
                  'do_iterations': niterations,
                  'mask_flag':  self.ui.cbMask.isChecked(),
                  'panic': False}
        if self.ui.cbLogFit.isChecked():
            _c = np.log(self.C)
        else:
            _c = self.C

        if self.ui.rbRefInter.isChecked():
            _slice = slice(self.ui.sbRefFrom.value(), self.ui.sbRefTo.value())
        else:
            _slice = slice(0, self.n_points)

        if self.B.ndim == 1:
            _beta_ = self.B
        else:
            _beta_ = self.B[_slice]

        new_c = NewtonRaphson(_c[_slice], _beta_, self.P, self.T[_slice],
                              **kwargs)
        if self.ui.cbLogFit.isChecked():
            new_c = np.exp(new_c)
        self.__renew_c(new_c, _slice)
        self.__dochisq()
        self.redraw()

    def _undo(self):
        assert self.prevC is not None
        self.C = self.prevC.copy()
        self.chisquared.pop()
        self.redraw()

    def __initial_guess(self):
        import libeq.consol
        self.C = libeq.consol.initial_guess(self.B, self.P, self.T)
        self.__dochisq()
        self.redraw()

    def __dochisq(self):
        self.F = fobj(self.C, self.P, self.T, self.B)
        new_chisq = np.sum(np.ma.masked_invalid(self.F)**2)
        # print(new_chisq)
        self.chisquared.append(new_chisq)

    def __renew_c(self, new_c, _slice):
        self.prevC = self.C.copy()
        self.C[_slice] = new_c
        self.ui.pbUndo.setEnabled(True)

    def __test_setup(self):
        # TESTING ONLY
        self.setC0Anal()
        self.do_damping()
        # self.doDamping1()
        # self.ui.cbDamping.setChecked(False)
        # self.ui.cbForcer.setChecked(False)
        # self.ui.cbScaling.setChecked(False)
        # self.ui.cbStepLimiter.setChecked(False)
        # self.ui.cb0OffD.setChecked(False)
        # self.ui.cbLogFit.setChecked(True)
        # self.doIter10()


if __name__ == '__main__':
    import sys
    from PyQt5 import QtCore
    app = QtWidgets.QApplication(sys.argv)

    with np.load('consol_panic.npz') as f:
        mainwindow = CSolver(f['beta'], f['analytc'], f['stoichiometry'],
                             f['free_concentration'])
    mainwindow.show()
    QtCore.pyqtRemoveInputHook()
    sys.exit(app.exec_())
