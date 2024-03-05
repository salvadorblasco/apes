"""otherwidgets.py module

Home to miscellaneous widgets not included in other category.
"""

import enum
import itertools
import math
import re

from PyQt5 import QtCore, QtWidgets
import numpy as np

import consts
import libaux
import libqt
import report
import ui_externaldata
import ui_manualfit
import ui_output
import ui_ionicst
import ui_titration


class ManualFitWidget(QtWidgets.QWidget):
    """A class for manually fitting constants."""

    def __init__(self, model, labels, datawidgets, plotmethod):
        super().__init__()
        self.ui = ui_manualfit.Ui_ManualFitWidget()
        self.ui.setupUi(self)

        self._model = model
        self.__colbetas = model.n_species
        self.__stoich = model.stoich
        self.__datawidgets = datawidgets
        self.__plotmethod = plotmethod

        self.ui.tabdata.clear()
        self.ui.tabdata.setRowCount(model.n_equil)
        self.ui.tabdata.setColumnCount(model.n_species + 2)
        header_labels = labels  + ["Value", "Gradient"]
        self.ui.tabdata.setHorizontalHeaderLabels(header_labels)

        libqt.array2tab(self.ui.tabdata, model.stoich)
        libqt.fill_column(self.ui.tabdata, model.n_species,
                          model.const, formatting='{:.4f}')

        for widget in libqt.iter_table_widgets(self.ui.tabdata):
            if widget is None:
                continue
            widget.setTextAlignment(QtCore.Qt.AlignRight)
            widget.setFlags(QtCore.Qt.NoItemFlags)

        for col in range(model.n_species):
            self.ui.tabdata.setColumnWidth(col, 40)

        drop_down_items = [" ".join(str(_) for _ in r) + f"  →  {k}"
                           for r, k in zip(model.stoich, model.const)]
        self.ui.cb_selectconst.addItems(drop_down_items)

        self.__center_value()

        self.ui.pb_update.clicked.connect(self.accept_values)
        self.ui.slider_constant.setTracking(False)
        self.ui.slider_constant.valueChanged.connect(self.__value_changed)
        self.ui.cb_selectconst.currentIndexChanged.connect(self.__beta_changed)

    def accept_values(self):
        # TODO implement
        pass

    def recalc(self):
        beta = self._current_betas()
        stoichiometry = self.__stoich

        for widget in self.__datawidgets:
            widget.calc_freec(altbeta=beta, altstoich=stoichiometry)
        self.__plotmethod()
        self._update_betas()
        self._update_gradient()
        self.__check_out_of_bounds()

    def _current_betas(self):
        retl = list(map(float, libqt.iter_column_text(self.ui.tabdata, col=self.__colbetas)))
        idx = self.ui.cb_selectconst.currentIndex()
        value = self._current_value()
        retl[idx] = value
        # print(idx, retl, value)
        return retl

    def _current_value(self):
        return self.ui.slider_constant.value() / 100

    def _update_betas(self):
        current_betas = self._current_betas()
        libqt.fill_column(self.ui.tabdata, self._model.n_species,
                          current_betas, formatting='{:.4f}')

    def _update_gradient(self):
        grads = [widget.gradient for widget in self.__datawidgets]
        resid = sum(np.sum(np.square(widget.residuals)) for widget in self.__datawidgets)
        grad = [sum(g) for g in zip(*grads)]
        libqt.fill_column(self.ui.tabdata, 1+self._model.n_species,
                          grad, formatting='{:.3f}')
        self.ui.lbl_chisq.setText(f"{resid:.3e}")
        idx = self.ui.cb_selectconst.currentIndex()
        curr_grad = grad[idx]
        max_grad = sum(abs(_) for _ in grad)
        ratio = 100*curr_grad//max_grad
        if ratio < 0.0:
            self.ui.pb_negative.setValue(-ratio)
            self.ui.pb_positive.setValue(0)
        else:
            self.ui.pb_positive.setValue(ratio)
            self.ui.pb_negative.setValue(0)


    def __beta_changed(self, index):
        self.__center_value()
        self._update_gradient()
        # TODO insert gradient value in lbl_gradient
        # TODO modify pb_negative and pb_positive

    def __center_value(self):
        row = self.ui.cb_selectconst.currentIndex()
        qwi = self.ui.tabdata.item(row, self.__colbetas)
        value = float(qwi.text())

        minimum = 100*math.floor(value) - 25
        maximum = 100*math.ceil(value) + 25
        self._bounduaries = (math.floor(value), math.ceil(value))
        self.ui.slider_constant.setMinimum(minimum)
        self.ui.slider_constant.setMaximum(maximum)
        self.ui.slider_constant.setSliderPosition(int(100*value))
        self.ui.lbl_minconst.setText("%.1f" % ((minimum+25)/100))
        self.ui.lbl_maxconst.setText("%.1f" % ((maximum-25)/100))
        # self.ui.lbl_current_value.setText("%.1f" % value)
        self.__value_changed(100*value)

    def __check_out_of_bounds(self):
        value = self.ui.slider_constant.value()
        if self._bounduaries[0] < value < self._bounduaries[1]:
            return
        else:
            value = value / 100
            minimum = 100*math.floor(value) - 25
            maximum = 100*math.ceil(value) + 25
            self._bounduaries = (math.floor(value), math.ceil(value))
            self.ui.slider_constant.setMinimum(minimum)
            self.ui.slider_constant.setMaximum(maximum)
            self.ui.lbl_minconst.setText("%.1f" % ((minimum+25)/100))
            self.ui.lbl_maxconst.setText("%.1f" % ((maximum-25)/100))

    def __value_changed(self, value):
        self.ui.lbl_current_value.setText("%.2f" % (value/100))


class IonicWidget(QtWidgets.QWidget):
    """A class for holding text.

    This class manages the output data. This widget diplays the output of
    the internal calculations in the internal widget self.ui.textBrowser.
    """

    def __init__(self):
        super().__init__()
        self.ui = ui_ionicst.Ui_IonicWidget()
        self.ui.setupUi(self)
        self.ui.pb_plus.clicked.connect(self.add_line)
        self.ui.pb_minus.clicked.connect(self.remove_line)

    def add_line(self):
        row = self.__choose_line()
        self.ui.table_titration.insertRow(row)

    def remove_line(self):
        row = self.__choose_line()
        self.ui.table_titration.removeRow(row - 1)

    @property
    def buret(self):
        yield from (float(i) for i in libqt.iter_column_text(self.ui.table_titration, 3))

    @property
    def initial_mmoles(self):
        ivol = self.initial_volume
        yield from (float(i)*ivol for i in libqt.iter_column_text(self.ui.table_titration, 2))

    @property
    def initial_concentration(self):
        yield from (float(i) for i in libqt.iter_column_text(self.ui.table_titration, 2))

    @property
    def ionic_strength(self):
        ivol = self.initial_volume
        concs = libaux.build_T_titr2(tuple(self.initial_mmoles),
                                     tuple(self.buret),
                                     ivol,
                                     self.titre)
        charges = tuple(self.charges)
        yield from (0.5*sum(c*z**2 for c, z in zip(conc, charges))
                    for vol, conc in zip(self.titre, concs))

    @property
    def charges(self):
        yield from (int(i) for i in libqt.iter_column_text(self.ui.table_titration, 1))

    @property
    def labels(self):
        yield from libqt.iter_column_text(self.ui.table_titration, 0)

    @property
    def initial_volume(self):
        return self.ui.dsb_V0.value()

    @property
    def final_volume(self):
        return self.ui.dsb_Vf.value()

    @property
    def titre(self):
        npoints = 100
        yield from libaux.linspace(0.0, self.final_volume-self.initial_volume, npoints)

    def __choose_line(self):
        line = self.ui.table_titration.currentRow()
        lines = self.ui.table_titration.rowCount()
        which = libqt.clamp_range(line, lines)
        return which


class OutputWidget(QtWidgets.QWidget):
    """A class for holding text.

    This class manages the output data. This widget diplays the output of
    the internal calculations in the internal widget self.ui.textBrowser.
    """

    def __init__(self):
        super().__init__()
        self.ui = ui_output.Ui_OutputWidget()
        self.ui.setupUi(self)
        self.__last_result = None
        signature = "Salvador Blasco &lt;salvador.blasco@gmail.com&gt;"
        html_style = """\
        p, li {
            white-space: pre-wrap;
            margin: 10px 0px 10px 0px; }
        th, td {
            border: 1px solid black;
        }
        """
        self.clear()
        self.ui.textBrowser.append("""<hr /><p align="center">
                 <b>APES</b>, the All-Purpose Equilibrium Solver</p>
                  <p align="center">(c) 2016-2022 %s</p>
              <hr>""" % signature)

        self.ui.pb_clear.clicked.connect(self.clear)
        self.report = None

    def appendH(self, level, txt):
        self.ui.textBrowser.append("<h%s>" % level + txt + "</h%d>" % level)

    def appendP(self, txt):
        """Append a new paragraph containing a particular text.

        Parameters:
            txt (str): The HTML text that will be appended. It will be enclosed
                in <p></p>."""
        self.ui.textBrowser.append("<p>" + txt + "</p>")

    def clear(self):
        "Clears the text"
        self.ui.textBrowser.clear()

    def addRuler(self):
        "appends <hr> tag"
        self.ui.textBrowser.append("<hr>")

    def fitting_header(self, **kwargs):
        self.ui.textBrowser.append(report.html_start(**kwargs))

    def get_last_result(self):
        return self.__last_result

    def report_final(self, *args, **kwargs):
        txt = report.html_finalconvg(*args, **kwargs)
        self.ui.textBrowser.append(txt)

    def report_iteration(self, *args, **kwargs):
        txt = self.report(*args, **kwargs)
        self.ui.textBrowser.append(txt)

    def save_last_result(self, **result):
        self.__last_result = result

    def set_report(self, function=consts.METHOD_LM):
        # if function == consts.METHOD_LM:
        #     self.report = report.html_lm_iteration
        # else:
        #     self.report = report.html_nm_iteration
        self.report = report.report_function[function]


class ExternalDataWidget(QtWidgets.QWidget):
    """A class for holding data to be plotted along with concentrations.

    It is a table that holds the data with or without errors. Data is organised
    by category. Each data group contains one X column and at least one Y
    column. Optionally, error X and error Y columns can be beside their
    respective data.
    """
    def __init__(self):
        super().__init__()
        self.ui = ui_externaldata.Ui_ExternalData()
        self.ui.setupUi(self)

        self.__row_labels = 0
        self.__col_xdata = 0
        self.__with_errors = True

        self.__set_popupmenu()

    def feed_data(self, x, y, ey=None):
        """Receive data and dump it to the table.

        Parameters:
            seq_data (sequence): Each element is a tuple with the following
                data: (1) order, (2) type, (3) suborder, (4) label, (5) data.
            max_depth (int): The longest value in the data provided. It is
                used to resize the table.

        .. note: mainly called by libio.loadExternalXML
        """
        # NOTE toy code: sandbox/colsort2.py
        table = self.ui.table
        table.clear()
        table.setRowCount(len(x)+1)
        if ey is None:
            table.setColumnCount(1+len(y))
        else:
            table.setColumnCount(1+len(y)+len(ey))
        libqt.replace_nones(table)
        libqt.fill_column(table, 0, x, row0=1)

        if ey is not None:
            for n, (_y, _ey) in enumerate(zip(y, ey)):
                libqt.fill_column(table, 1+2*n, _y, row0=1)
                libqt.fill_column(table, 2+2*n, _ey, row0=1)
            headers = ['X'] + [b+a for a, b in itertools.product(map(str, range(1, 1+len(y))),
                                                                 ('Y', 'err Y'))]
        else:
            for n, _y in enumerate(y):
                libqt.fill_column(table, 1+n, _y, row0=1)
            headers = ['X'] + ['Y%d' % a for a in range(1, 1+len(y))]

        self.ui.table.setHorizontalHeaderLabels(headers)
        self.__renumber_rows()

    def iter_data(self):
        """Yield the data row-wise.

        Yields:
            tuple: data (float) with error (float) or None.
        """
        cols = range(1, self.ui.table.columnCount(),
                     2 if self.__with_errors else 1)
        for col in cols:
            data = libqt.iter_column_text(self.ui.table, col, row0=1)
            errd = libqt.iter_column_text(self.ui.table, 1+col, row0=1) \
                   if self.__with_errors else None
            yield data, errd

    def get_data(self):
        """The data in the table in an usable format.

        Return a dict in which the keys are category of the data. If there is
        only one category, then there is only one key which is None. The data
        is a tuple in which the first element is the type (X, Y, eX, eY) and
        the second is a tuple containing the data (float numbers) associated
        with it.
        """

        def __split_header(txt):
            r = re.match(r'(e?[xy])\((\d+)\)?', txt)
            return ('1', '2') if r is None else r.group(1, 2)

        cols = range(self.ui.table.columnCount())
        # TODO back to generator
        headers = [self.ui.table.horizontalHeaderItem(col) for col in cols]
        htext = (h.text() if h is not None else "" for h in headers)
        hh = (__split_header(h) for h in htext)
        out = {}
        for col, (tp, cat) in enumerate(hh):
            data = libqt.iter_column_text(self.ui.table, col, row0=1)
            item = (tp, tuple(float(d) for d in data))
            if cat in out:
                out[cat].append(item)
            else:
                out[cat] = [item]
        return out

    def labels(self):
        "The labels."
        label_cols = range(1, self.ui.table.columnCount(), 2 if self.__with_errors else 1)
        return libqt.iter_row_text(self.ui.table, row=self.__row_labels, iter_range=label_cols)

    def set_labels(self, label_sequence):
        """Set the labels in the table.

        Parameters:
            label_sequence (iterable): the labels.
        """
        lst = ['']
        if self.__with_errors:
            lst.extend(itertools.chain.from_iterable(zip(label_sequence,
                                                         len(label_sequence)*[''])))
        else:
            lst.extend(label_sequence)
        libqt.fill_row(self.ui.table, row=0, data=lst)

        itrow = libqt.iter_row(self.ui.table, row=0)
        for l, item in zip(lst, itrow):
            if l:
                item.setFlags(QtCore.Qt.ItemIsSelectable)

    def load_from_data(self, filename):
        """Load data from file.

        Parameters:
            filename (str): a valid file name to open.
        """
        with open(filename, 'r') as fhandler:
            self.set_title(fhandler.readline().strip())
            labels = fhandler.readline().split()
            extdat = np.loadtxt(fhandler)

        n_data, n_vars = extdat.shape
        table = self.ui.table
        table.clear()
        table.setRowCount(n_data+1)
        table.setColumnCount(n_vars)
        libqt.replace_nones(table)
        self.set_labels(labels)
        libqt.array2tab(table, extdat, row0=1)
        self.__renumber_rows()
        self.__renumber_cols()

    @property
    def n_variables(self):
        "The number of variables."
        # if self.__with_errors:
        #     retv = (self.ui.table.columnCount() - 1) // 2
        # else:
        #     retv = (self.ui.table.columnCount() - 1)
        # return retv
        retv = self.ui.table.columnCount() - 1
        return retv // 2 if self.__with_errors else retv

    def set_title(self, title):
        """Change the title."""
        self.ui.le_title.setText(title)

    def title(self):
        "The title."
        return self.ui.le_title.text()

    def xdata(self):
        "Data contained in the X column."
        xdtab = libqt.iter_column_text(self.ui.table, self.__col_xdata, row0=1)
        return tuple(map(float, xdtab))

    def _popupmenu(self, point):
        self.popupmenu.popup(self.ui.table.viewport().mapToGlobal(point))

    def _xcols(self):
        headers = ((i, self.ui.table.takeHorizontalHeaderItem(i))
                   for i in range(self.ui.table.columnCount()))
        rval = None
        for col, header in headers:
            text = header.text()
            if text == 'X':
                rval = col
        return rval

    def __popm_addvars(self):
        nvar, accepted = QtWidgets.QInputDialog.getInt(self, 'Add data',
                                                       'How many variables?', min=0)
        if accepted:
            if self.__with_errors:
                cols = 2*nvar + 1
            else:
                cols = nvar + 1
            self.ui.table.setColumnCount(cols)
            self.__renumber_cols()

    def __popm_addrows(self):
        rows, accepted = QtWidgets.QInputDialog.getInt(self, 'Add data',
                                                       'How many rows?', min=0)
        if accepted:
            self.ui.table.setRowCount(1+rows)
            self.__renumber_rows()

    def __popm_fromfile(self):
        'Add data from file.'

        filters = "Text Files (*.txt);;All Files (*.*)"
        filename, accepted = libqt.opendialog(self, filters)
        if not accepted:
            return

        self.load_from_data(filename)

    def __popm_show_errors(self):
        pass

    def __popm_clear_table(self):
        pass

    def __renumber_cols(self):
        table = self.ui.table
        lbls = ['X']
        for i in range(self.n_variables):
            lbls.append('Y' + str(i+1))
            if self.__with_errors:
                lbls.append('err Y' + str(i+1))
        table.setHorizontalHeaderLabels(lbls)

    def __renumber_rows(self):
        table = self.ui.table
        lbls = ['labels'] + [str(i) for i in range(1, table.rowCount())]
        table.setVerticalHeaderLabels(lbls)

    def __set_popupmenu(self):
        self.ui.table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.table.customContextMenuRequested.connect(self._popupmenu)
        self.popupmenu = QtWidgets.QMenu(self)
        tabcons = (#('Data from file', self.__popm_fromfile),
            ('Set rows size', self.__popm_addrows),
            ('Set variable number', self.__popm_addvars),
            ('Show errors', self.__popm_show_errors),
            ('Clear table', self.__popm_clear_table))

        for label, call in tabcons:
            action = QtWidgets.QAction(label, self)
            action.triggered.connect(call)
            self.popupmenu.addAction(action)

    @classmethod
    def __filter(cls, tag, *args):
        yield from filter(lambda a: a[0] == tag, args)

    # def __split_header(self, txt):
    #     import re
    #     r = re.match(r'(e?[xy])\((\d+)\)?', txt)
    #     # if r is None:
    #     #     # raise ValueError('Header is non-conforming')
    #     #     return '1', '2'
    #     # else:
    #     #     return r.group(1, 2)
    #     return ('1', '2') if r is None else r.group(1, 2)

class TitrationBaseWidget(QtWidgets.QWidget):
    def __init__(self, model):
        super().__init__()
        self.name = None
        self.ui = ui_titration.Ui_Titration()
        self.ui.setupUi(self)
        self._cols = enum.IntEnum('col', 'label init init_flags buret buret_flags', start=0)
        self._column_buret_flags = 4
        self._column_init_flags = 2
        n = model.number_components
        self.init_flags = [consts.RF_CONSTANT] * n
        self.buret_flags = [consts.RF_CONSTANT] * n

    def is_titre_implicit(self):
        return self.ui.dsb_Vf.isEnabled()

    @property
    def initial_amount(self):
        """The inital amount values in millimole.

        Returns:
            tuple of floats containing the values of this parameter.
        """
        ctext = libqt.iter_column_text(self.ui.table_titration, col=self._cols.init.value)
        return tuple(float(b) for b in ctext)

    @initial_amount.setter
    def initial_amount(self, initial_amount):
        """The initial amount values in millimole.

        Parameters:
            initial_amount (sequence of floats): the values for this parameter.
                The length must be equal to the number of components.
        """
        libqt.fill_column(self.ui.table_titration, col=self._cols.init.value, data=initial_amount)

    @property
    def n_points(self):
        return self.ui.sb_NPoints.value()

    @n_points.setter
    def n_points(self, value):
        return self.ui.sb_NPoints.setValue(value)

    @property
    def buret(self):
        """The inital amount values in millimole.

        Returns:
            tuple of floats containing the values of this parameter.
        """
        ctext = libqt.iter_column_text(self.ui.table_titration, col=self._cols.buret.value)
        return tuple(float(b) for b in ctext)

    @buret.setter
    def buret(self, buret):
        """Set the buret parameter.

        Parameters:
            buret (sequence of floats): the values for buret. The length must be equal
                to the number of components.
        """
        libqt.fill_column(self.ui.table_titration, data=buret, col=self._cols.buret.value)

    @property
    def starting_volume(self):
        'The starting volume in mL.'
        return self.ui.dsb_V0.value()

    @starting_volume.setter
    def starting_volume(self, volume: float):
        'Set the starting volume in mL.'
        assert isinstance(volume, float)
        self.ui.dsb_V0.setValue(volume)

    @property
    def final_volume(self):
        'The value of the final titre in mL.'
        return self.ui.dsb_Vf.value()

    @final_volume.setter
    def final_volume(self, volume: float):
        'Set the final volume in mL.'
        assert isinstance(volume, float)
        # if volume > self.starting_volume():
        #     raise ValueError('Final volume must be bigger than start volume')
        self.ui.dsb_Vf.setValue(volume)

    def set_labels(self, labels):
        libqt.fill_column(self.ui.table_titration, data=labels, col=self._cols.label)

    def set_volume_implicit(self, value: bool):
        for widget in (self.ui.dsb_Vf, self.ui.sb_NPoints):
            widget.setEnabled(value)

    def titre(self):
        """The titre values in mL.

        Yields:
            float: the value of the titre.
        """
        yield from ((self.final_volume - self.starting_volume)
                    / self.n_points() * i for i in range(self.n_points()))

    def volume_increment(self):
        'The step volume in mL.'
        return (self.final_volume() - self.starting_volume()) / self.n_points()

    @property
    def buret_flags(self):
        indices = libqt.iter_column_comboindex(self.ui.table_titration, self._column_buret_flags)
        return tuple(item - 1 for item in indices)

    @buret_flags.setter
    def buret_flags(self, flags):
        flagwidgets = (libqt.create_combo(consts.REFINE_LABELS, flag) for flag in flags)
        with libqt.table_locked(self.ui.table_titration):
            for row, widget in enumerate(flagwidgets):
                self.ui.table_titration.setCellWidget(row, self._column_buret_flags, widget)

    @property
    def init_flags(self):
        indices = libqt.iter_column_comboindex(self.ui.table_titration, self._column_init_flags)
        return tuple(item - 1 for item in indices)

    @init_flags.setter
    def init_flags(self, flags):
        flagwidgets = (libqt.create_combo(consts.REFINE_LABELS, flag) for flag in flags)
        with libqt.table_locked(self.ui.table_titration):
            for row, widget in enumerate(flagwidgets):
                self.ui.table_titration.setCellWidget(row, self._column_init_flags, widget)

    @property
    def volume_error(self):
        return self.ui.dsb_Verror

    @volume_error.setter
    def volume_error(self, error: float):
        self.ui.dsb_Verror.setValue(error)
