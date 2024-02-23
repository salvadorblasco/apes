"""
This module contains the canvas where the plots are drawn.

.. module:: canvas.py
.. moduleauthor:: Salvador Blasco <salvador.blasco@protonmail.com>

"""

import itertools
import dataclasses

import numpy as np
from scipy.interpolate import interp1d

from PyQt5 import QtWidgets
from PyQt5 import QtCore

import matplotlib
from matplotlib.backends.backend_qt5agg import \
    FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import libaux
import libplot


class MyCanvas(FigureCanvas):
    """Main class to do the drawing.

    This class contains the canvas where all the graphs are displayed. Also
    handles events and such.
    """

    # TODO deprecate; replace by PlotStyle
    STYLE1, STYLE2, STYLE3, STYLE4 = range(4)
    COLORBW, COLORR, COLORB, COLORRB, COLORGROUP = range(5)

    def __init__(self, parent=None, width=6, height=8, dpi=300):
        """Initialize the images for the canvas.

        Usually the default parameters are enough.
        """
        # TODO change this to uselocal .matplotlibrc file
        self.__font_size = 10.0
        matplotlib.rc('font', size=self.__font_size)
        matplotlib.rc('figure.subplot', wspace=0.0, hspace=0.30, left=0.15,
                      right=0.90)

        # Default parameters
        # self._speciation_options = {
        #     'ignore_lower_than': 4,
        #     'unlabel_lower_than': 9,
        #     'plot_labels': True
        # }
        self.__pickradius = 2
        self.title = None
        self.current_style = 0
        self.current_color = 0
        self._plotlabels = True
        self._plotlegend = False

        # the dragable labels
        self.drlbls = []

        # reference to a list of floats containing the chisq data
        self._chisq = None

        # the figure
        self._figure = Figure(figsize=(width, height), dpi=dpi)

        # the format of labels.
        # blank means use default setting
        self._labelformat = ""

        FigureCanvas.__init__(self, self._figure)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.data = None
        self.axes = {'data': [], 'conc': [], 'resid': [], 'fit': None}

        self.__setup_popupmenu()
        # cid = self.mpl_connect('button_press_event', self.on_click)

    @QtCore.pyqtSlot(float)
    def label_size_changed(self, new_size):
        """Slot for when label size changes"""
        self.__font_size = new_size
        matplotlib.rc('font', size=self.__font_size)
        print("font size is now", new_size)

    @QtCore.pyqtSlot(int)
    def style_color_changed(self, new_color_id):
        """Slot for when label size changes"""
        pass

    @QtCore.pyqtSlot(int)
    def style_scheme_changed(self, new_style_id):
        """Slot for when label size changes"""
        pass

    def save_figure(self, f):
        """Save current figure in a file.

        Arguments:
            f (file or str): The file to which the figure will be saved.
        """
        self._figure.savefig(f)

    def plot_calor(self, widget):
        """Plot calorimetry fit.

        This function plots all the information for a caloriemtry fit. It
        takes the data available from the internal *data* variable and plots
        it. The procedure is that it divides the figure in three parts. On top
        the experimental data, corresponding fit and calculated
        concentrations are plotted. In the middle, the residuals are plotted.
        On the bottom, one single plot with the variation of the residuals is
        displayed.

        Parameters:
            widget (calorwidget.CalorWidget): The calorimetry data widget
        """
        self.axes = {'data': [], 'conc': [], 'resid': [], 'fit': None}
        ax_plots = self.axes['data']
        ax_conc = self.axes['conc']
        ax_resid = self.axes['resid']

        d = self.data.plotCalor()
        NCa = next(d)
        self.figure.clear()

        for n, (vy, vn, qy, qn, vf, c, qf, ry, rn) in enumerate(zip(*d)):
            pd = self.figure.add_subplot(2, NCa, n+1)
            assert len(vy) == len(qy), "%d != %d" % (len(vy), len(qy))
            assert len(vn) == len(qn), "%d != %d" % (len(vn), len(qn))
            pd.scatter(vy, qy, c='blue')
            pd.scatter(vn, qn, c='red')
            if qf is not None:
                pd.plot(vf[1:], qf, 'r-')
            pd.set_xlabel('V / mL')
            # pd.set_ylabel('Q / kcal')
            ax_plots.append(pd)

            assert vf is not None
            if c is not None:
                pc = pd.twinx()
                pc.plot(vf, c, 'k-')
                ax_conc.append(pc)

            pr = self.figure.add_subplot(6, NCa, 3*NCa+n+1)
            pr.scatter(vy, ry, c='blue')
            pr.scatter(vn, rn, c='red')
            pr.set_ylabel('residual')
            ax_resid.append(pr)

        libplot.evenscale('y', *ax_plots)
        if len(ax_conc) > 0:
            libplot.evenscale('y', *ax_conc)
        libplot.evenscale('y', *ax_resid)

        self.axes['fit'] = self.figure.add_subplot(3, 1, 3)
        self.draw_fit_axes()
        ax_resid[0].set_ylabel("residual")
        ax_plots[0].set_ylabel("Q / kcal")

        self.draw()

    def plot_emf(self, widgets, fitdata={}):
        """Plot potentiometry fit data.

        This function plots all the information for a potentiometry fit. It
        takes the data available from the internal *data* variable and plots
        it. The procedure is that it divides the figure in three parts. On top
        the experimental data, corresponding fit and calculated
        concentrations are plotted. In the middle, the residuals are plotted.
        On the bottom, one single plot with the variation of the residuals is
        displayed.

        Parameters:
            widgets (sequence): sequence of :class:`apes.EmfWidget` type.
            fitdata (dict): additional data generated in the fiiting.
        Raises:
            ValueError: if any of the widgets if of the wrong type.
        """

        def draw_emf(ap, titre, emfr, emfc, plot_masked=True):
            # breakpoint()
            ap.scatter(titre, emfr.T, c='blue')

            if plot_masked and isinstance(emfr, np.ma.masked_array):
                ern = np.ma.array(emfr.data, mask=~emfr.mask)
                ap.scatter(titre, ern, c='red')

            if emfc is not None:
                ap.plot(titre, emfc.data, 'r-')

        def draw_fit(ax, conv):
            length = len(conv)
            ax.bar(np.arange(0.0, 0.0+length), np.array(conv), width=1)
            step = length//10
            if step == 0:
                step = 1
            if max(conv) - min(conv) > 100:
                ax.set_yscale('log')
            else:
                ax.set_yscale('linear')
                formatter = matplotlib.ticker.ScalarFormatter()
                formatter.set_powerlimits((-2, 3))
                formatter.set_scientific(True)
                ax.yaxis.set_major_formatter(formatter)
                ax.set_ylim((0.95*min(conv), 1.05*max(conv)))
            ax.set_xticks(np.arange(0, 1+length, step))

        def draw_conc(ac, v, c):
            if c is None:
                return
            ac.plot(v, c, 'k-')

        def draw_res(ar, v, r):
            if r is None:
                return

            if np.ma.isMaskedArray(r):     #  TODO stem plots do not support
                rr = r.compressed()        # masked arrays for the moment.
                vv = v[~r.mask.flatten()]  # remove this maybe in the future.
            else:
                rr = r
                vv = v

            ar.stem(vv, rr, linefmt='k-', markerfmt=' ', basefmt='k:',
                    use_line_collection=True)

        n_series = len(widgets)
        # TODO put this as function argument
        plot_masked = True

        figure = self.figure
        figure.clear()
        gs = matplotlib.gridspec.GridSpec(nrows=10, ncols=n_series,
                                          wspace=0.0, hspace=0.0)
        ax_plots = [figure.add_subplot(gs[0:5, i]) for i in range(n_series)]
        ax_conc = [ax.twinx() for ax in ax_plots]
        ax_conc[-1].set_ylabel("concentration")
        for a in ax_conc:
            a.set_navigate(False)
        ax_resid = [figure.add_subplot(gs[5:7, i], sharex=ax_plots[i]) for i in range(n_series)]
        ax_fit = figure.add_subplot(gs[8:, :])

        for ax in ax_plots[1:]:
            ax.set_yticklabels([])
        for ax in ax_resid[1:]:
            ax.set_yticklabels([])
        for ax in ax_resid:
            ax.set_xlabel("V / mL")
            ax.set_navigate(False)
        ax_plots[0].set_ylabel("emf / mV")
        ax_resid[0].set_ylabel("residual")
        ax_fit.set_xlabel('iteration')
        ax_fit.set_ylabel('χ₂')
        formatter = matplotlib.ticker.ScalarFormatter()
        formatter.set_powerlimits((-2, 3))
        formatter.set_scientific(True)
        ax_conc[-1].yaxis.set_major_formatter(formatter)

        z = zip(widgets, ax_plots, ax_conc, ax_resid)
        for n, (w, ap, ac, ar) in enumerate(z):
            # u = w.use                               # u  → the use flag
            v = np.fromiter(w.titre, dtype=np.float)  # v  → Volume of titre
            mask = np.fromiter(w.mask, np.bool)
            masked = np.any(mask)
            er = np.array(w.emf)                  # er → Experimental data
            ef = w.emf_fit                            # ef → Fitted emf
            c = w.free_concentration                  # c  → Free conc.
            r = er - ef if ef is not None else None   # r  → Residuals

            if masked:
                er = np.ma.masked_array(er, mask=mask)
                r = np.ma.masked_array(r, mask=mask)

            draw_emf(ap, v, er, ef, plot_masked)
            draw_conc(ac, v, c)
            draw_res(ar, v, r)

        libplot.evenscale('y', *ax_plots)
        libplot.evenscale('y', *ax_conc)
        libplot.evenscale('y', *ax_resid)

        for ax in zip(ax_plots, ax_resid):
            libplot.evenscale('x', *ax)
        for ax in ax_conc[:-1]:
            ax.set_yticklabels([])

        if 'convergence' in fitdata:
            draw_fit(ax_fit, fitdata['convergence'])

        self.draw()

    def plot_ionic(self, **kwargs):
        self._figure.clear()
        axes = self._figure.add_subplot(111)

        axes.set_xlabel("titre / mL")
        axes.set_ylabel("Ionic strength / mol/L")
        axes.plot(kwargs['titre'], kwargs['ionic'])

        ax2 = axes.twinx()
        ax2.set_ylabel("Percent change")
        Imin, Imax = axes.get_ylim()
        I0 = kwargs['ionic'][0]
        ax2.set_ylim(bottom=(Imin-I0)/I0*100,
                     top=(Imax-I0)/I0*100)
        self.draw()


    def plot_speciation(self, widget, externaldata=None, **plotoptions):
        """Plot Titration simulation or species ditribuction.

        Parameters:
            widget (:class:`DataWidget`): the widget to plot
            externaldata (:class:`ExternalData`, optional): 
        """
        x_value = widget.xdata()
        y_value = widget.ydata()
        species_labels = widget.plot_labels()
        x_label = widget.xlabel()
        y_label = widget.ylabel()
        if plot_errors := widget.plot_errors():
            y_error = widget.yerrdata()

        max_y_value = np.max(y_value)
        if 'ignore_lower_than' in plotoptions:
            cutvalue = plotoptions['ignore_lower_than']*max_y_value/100
            _y_value = np.ma.array(y_value, mask=y_value < cutvalue)
        else:
            _y_value = y_value

        groups = widget.groups() if plotoptions.get('color_by_group', False) else None
        what, how = self.style(num_entries=y_value.shape[1], withfill=plot_errors, groups=groups)
        # bools indicating what to plot
        plot_line, plot_error_fill, plot_error_line = what
        # generators indicating how to plot
        styp, styef, styel = how

        self._figure.clear()
        axes = self._figure.add_subplot(111)

        axes.set_xlabel(x_label)
        axes.set_ylabel(y_label)

        if 'y_error' not in plotoptions:
            y_error = itertools.cycle((None,))
        else:
            y_error = plotoptions['y_error'].T

        for value, error, label in zip(_y_value.T, y_error, species_labels):
            if error is not None:
                if plot_error_fill:
                    axes.fill_between(x_value, value-error, value+error,
                                      **next(styef), label="erf(%s)" % label)
                if plot_error_line:
                    sty = next(styel)
                    axes.plot(x_value, value-error, **sty,
                              label="erlwr(%s)" % label)
                    axes.plot(x_value, value+error, **sty,
                              label="erupr(%s)" % label)
            if plot_line:
                axes.plot(x_value, value, **next(styp),
                          pickradius=self.__pickradius,
                          label="plot(%s)" % label)

        if plotoptions.get('plot_labels', True):
            window = (axes.get_xlim(), axes.get_ylim())
            label_xy = libplot.position_labels(x_value, y_value, window)
            labels_cutoff = plotoptions.get('unlabel_lower_than', 5.0)*max_y_value/100
            self.place_labels(axes, species_labels, label_xy,
                              threshold=labels_cutoff)

        for edata in externaldata:
            ext_axes = axes.twinx()
            ext_x = edata.xdata()
            for label, (exydat, exedat) in zip(edata.labels(), edata.iter_data()):
                _y = tuple(map(float, exydat))
                ext_axes.scatter(ext_x, _y, label=label)
                if exedat is not None:
                    _ey = tuple(map(float, exedat))
                    ext_axes.errorbar(ext_x, _y, _ey, capsize=2.0)
            ext_axes.set_ylabel(edata.title())

        if 'plot_title' in plotoptions:
            axes.set_title(plotoptions['plot_title'])

        # axes.set_ylim(bottom=ymin, top=ymax)
        self.draw()
        # cid = self.mpl_connect('button_press_event', self.on_click_speciation)

    def plotSpectr(self, widget):
        """Plot all the information for a spectrometry fit.

        It takes the data available from the internal *data* variable and plots
        it. The procedure is that it divides the figure in four parts. On top
        the experimental data, corresponding fit and calculated
        concentrations are plotted. In the middle, the residuals are plotted.
        On the bottom left, the fit spectra are shown and on bottom right
        variation of the residuals is displayed.
        """
        if self.data is None:
            raise RuntimeError('No data available')

        if self.data.NSpec == 0:
            raise RuntimeError('No spectra data available')

        d = self.data.plotSpectr()

        NDs = next(d)   # number of datasets
        NSp = next(d)   # total number of spectra sets

        self.figure.clear()
        self.axes = {'data': [], 'conc': [], 'resid': [], 'fit': None}
        ax_plots = self.axes['data']
        # ax_conc = self.axes['conc']
        ax_resid = self.axes['resid']

        # I. Plot spectra
        for n, s in enumerate(next(d)):
            p = self.figure.add_subplot(2, NSp, n+1)
            p.plot(s[:, 0], s[:, 1:])
            p.set_xlabel('λ / nm')
            ax_plots.append(p)

        # II. Plot residuals
        for n, s in enumerate(next(d)):
            p = self.figure.add_subplot(6, NSp, 4*NSp+n)
            if s is not None:
                p.plot(s[:, 0], s[:, 1:])
            ax_resid.append(p)

        # III. Plot fitted spectra
        raise NotImplementedError
        # TODO implement
        # !for n, (t, s) in enumerate(zip(next(d), next(d))):
        # !    if t == data.SpecData.MODEL_LINEAR:
        # !        p = self.figure.add_subplot(3, NDs + 1, 2*(NDs+1)+n+1)
        # !        p.set_xlabel('λ / nm')
        # !        ax_conc.append(p)
        # !        if s is not None:
        # !            p.plot(s[:, 0], s[:, 1:])
        # !    else:
        # !        raise NotImplementedError

        # IV. Plot χ²
        self.axes['fit'] = self.figure.add_subplot(3, NDs + 1, 3*(NDs+1))
        self.axes['fit'].yaxis.tick_right()
        self.axes['fit'].yaxis.set_label_position("right")
        self.draw_fit_axes()
        self.draw()

    def draw_fit_axes(self):
        ax = self.axes['fit']
        ax.set_ylabel("χ²")
        ax.set_xlabel("iteration")
        ax.set_yscale('log')
        if not self._chisq:
            return
        length = len(self._chisq)

        ax.bar(np.arange(0.0, 0.0+length), np.array(self._chisq), width=1)
        step = length//10
        if step == 0:
            step = 1
        ax.set_xticks(np.arange(0, 1+length, step))
        formatter = matplotlib.ticker.ScalarFormatter()
        formatter.set_powerlimits((-2, 3))
        formatter.set_scientific(True)
        ax.yaxis.set_major_formatter(formatter)
        ax.set_ylim((0.95*min(self._chisq), 1.05*max(self._chisq)))
        return

    def __setup_popupmenu(self):
        self.menu = QtWidgets.QMenu(self)

        self.__action_set_title = QtWidgets.QAction("set title", self)
        self.menu.addAction(self.__action_set_title)

    def gca(self):
        "Return an instance of the current axes"
        return self.axes

    def on_click(self, event):
        'Handle click event.'
        # Check click on curves
        for obj in self.plot_objects:
            if isinstance(obj, Line2D):
                c, a = obj.contains(event)
                if c:
                    # print "clicked on object ", obj
                    cl = obj.get_color()
                    ls = obj.get_linestyle()
                    lw = obj.get_linewidth()

                    d = libplot.LinePropertiesDialog(self)
                    d.set_initial(cl, ls, lw)
                    d.exec_()
                    if d.result() == QtWidgets.QDialog.Accepted:
                        ans = d.get_result()
                        obj.set_linestyle(ans['linestyle'])
                        obj.set_color(ans['color'])
                        obj.set_linewidth(ans['linewidth'])
                        self.draw()
                    return

        # TODO check click on titles
        if self.title_object is not None:
            c, d = self.title_object.contains(event)
            if c:
                # print("clicked on title: ", d)
                # TODO open dialog
                pass

        xaxis_title = self.axes.get_xaxis().get_label()
        c, d = xaxis_title.contains(event)
        if c:
            # print("clicked on xaxes title")
            # TODO open dialog
            pass

        yaxis_title = self.axes.get_yaxis().get_label()
        c, d = yaxis_title.contains(event)
        if c:
            # print("clicked on yaxes title")
            # TODO open dialog
            pass

        # TODO check click on axis

        # TODO Correct coordinates
        point = QtCore.QPoint(event.x, event.y)

        # TODO If right click, open popup menu
        if event.button == 2:
            self.menu.popup(self.mapToGlobal(point))

        # TODO If hold click, drag element if dragable

        # TODO If double click, open element properties if aplicable
        if event.button == 0 and event.dblclick:
            pass

    def place_labels(self, axes, label_text, label_xy, threshold=0.0):
        """Create labels and place them.

        Parameters:
            axes (:class:`matplotlib.Axes.axes`): the axes to place the labels into.
            label_text (iterable): the text for the labels
            label_xy (iterable): the position where the labels should be.
            threshold (float): if the maximum point is lower than this, no
                label is placed.
        """
        assert len(label_text) == len(label_xy)
        label_objs = []
        self.drlbls.clear()
        for xy, l in zip(label_xy, label_text):
            assert isinstance(xy, tuple)
            assert len(xy) == 2
            assert isinstance(xy[0], float)
            assert isinstance(xy[1], float)
            if xy[1] < threshold:
                continue
            label_objs.append(
                axes.text(
                    xy[0], xy[1], l,
                    verticalalignment='bottom',
                    horizontalalignment='center',
                    color='k',
                    family="serif"
                )
            )
            drlbl = libplot.DragableLabel(label_objs[-1])
            drlbl.connect()
            drlbl.labelsChanged.connect(self.__change_labels)
            self.drlbls.append(drlbl)

    def style(self, num_entries=None, withfill=False, groups=None):
        """Generate style data.

        This function returns a generator which returns a :class:`dict` to
        be plugged into the :func:`matplotlib.pyplot.plot` function.

        The styles available are the following:

        - **style 1** Solid line
        - **style 2** solid line
        - **style 3** Solid line
        - **style 4** solid line

        Parameters:
            num_entries (int): The number of style entries to be generated. If it is
                None then an infinite generator will be returned. Ignored if
                groups is provided.
            withfill (bool):
            groups (dict): a dict that indicates how the entries are grouped.

        Returns:
            tuple: a two-element tuple. The first element is a bool tuple
                which means whether to plot (1) the line, (2) the fill
                and (3) the error line. The second element is a tuple
                of dict with information on how to plot the line, the
                fill and the error line.
        """
        if groups is None:
            # groups = {i:i for i in range(num_entries)}
            groups = list(range(num_entries))

        num_groups = len(set(groups))

        stroke_styles = {0: ((True, True, False), ('-', 1.0, ' ', 0.0), ('-', 2.0)),
                         1: ((False, True, True), (' ', 0.0, '-', 0.5), ('--', 2.0)),
                         2: ((True, True, True), ('-', 1.5, '-', 0.0), ('-', 2.0)),
                         3: ((True, True, True), ('--', 1.5, '-', 0.0), ('-', 2.0))}

        stroke = stroke_styles[self.current_style]
        if withfill:
            ret0, aux, _ = stroke
            ilsp, ilwp, ilsf, ilwf = map(itertools.repeat, aux)
        else:
            aux = stroke[2]
            ret0 = True, False, False
            ilsp, ilwp = map(itertools.repeat, aux)

        color_styles = {1: ('autumn', 'r'), 2: ('winter', 'b'), 3:('hsv', 'k')}

        if self.current_color == 0:              # black and white
            icolorp = itertools.repeat('k')
            icolorf = itertools.repeat('grey')
            icolorl = itertools.repeat('k')
        elif self.current_color in color_styles:
            cmapname, altcolor = color_styles[self.current_color]
            mycm = matplotlib.cm.get_cmap(cmapname)
            gcmap = [mycm(i/num_groups) for i in range(num_groups)]
            # gcmap = [mycm(i/num_groups, 1+num_groups) for i in range(num_groups)]
            # icolorp = itertools.repeat(altcolor) if withfill else (gcmap[v] \
            #                                      for v in groups.values())
            icolorp = itertools.repeat(altcolor) if withfill else (gcmap[v] \
                                                 for v in groups)
            # icolorf = (gcmap[i] for n, i in groups.items()) if withfill else None
            icolorf = (gcmap[i] for n, i in enumerate(groups)) if withfill else None
            icolorl = itertools.repeat(altcolor) if withfill else None
        elif self.current_color == 4:              # by group
            pass

        ret1 = ({'color': c, 'linestyle': ls, 'linewidth': lw}
                for c, ls, lw in zip(icolorp, ilsp, ilwp))

        ret2 = ({'color': c} for c in zip(icolorf)) if withfill else None
        ret3 = ({'color': c,
                 'linestyle': ls,
                 'linewidth': lw}
                for c, ls, lw in zip(icolorl, ilsp, ilwp)) if withfill else None

        return ret0, (ret1, ret2, ret3)

    def setStyle(self, s):
        """Sets the type of style to be used for plotting

        Parameters:
            s (int): a number identifying the type of style to be used
        """
        libaux.assert_type(int, s)
        if s < 0 or s > 4:
            raise ValueError("0 < s < 4")
        self.current_style = s

    def setColor(self, c):
        """Sets the type of color to be used for plotting

        Parameters:
            c (int): a number identifying the type of color scheme to be used

        Raises:
            ValueError: if the value of c is outside acceptable range
        """
        libaux.assert_type(int, c)
        if c < 0 or c > 4:
            raise ValueError("0 < c < 4")
        self.current_color = c

    def __change_labels(self, textobject):
        newlabel, ok = QtWidgets.QInputDialog.getText(self, "Label Edit", "label:", 
                                                      text=textobject.get_text())
        if not ok:
            return
        textobject.set_text(str(newlabel))
        textobject.figure.canvas.draw()

    def __line_clicked(self, event):
        # TODO Complete
        print("line clicked!!")

    def __prepare_menus(self):
        menu = QtWidgets.QMenu(parent=self)
        action_aed = QtWidgets.QAction("add external data")
        action_aed.triggered.connect(self.add_external_data)
        menu.addAction(action_aed)


def smooth_curve(xmin, xmax, ydata):
    """Smooth a curve by interpolation.

    Parameters:
        xmin (float):
        xmax (float):
        ydata (:class:`numpy.ndarray`):

    Returns: tuple where first element in the new X data and the second is the
        new Y data.
    """
    num_data = ydata.shape[0]
    new_xdata = np.linspace(xmin, xmax, num_data)
    finterp = interp1d(new_xdata, ydata, kind='cubic', axis=0)
    new_ydata = finterp(ydata)
    return new_xdata, new_ydata


@dataclasses.dataclass
class PlotStyle:
    linestyle: str
    linewidth: float
    plot_error: bool

    def kwargs(self):
        return {'linestyle': self.linestyle}
