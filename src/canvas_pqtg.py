"""
This module contains the canvas where the plots are drawn with pyqtgraph as plotter.

.. module:: canvas_pqtg.py
.. moduleauthor:: Salvador Blasco <salvador.blasco@protonmail.com>

"""

import itertools
import dataclasses

import numpy as np
from scipy.interpolate import interp1d

from PyQt5 import QtWidgets
from PyQt5 import QtCore

import pyqtgraph

import libaux
import libplot


class MyCanvas(pyqtgraph.GraphicsLayoutWidget):
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
        super().__init__(parent, size=(width, height))
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

        # the format of labels.
        # blank means use default setting
        self._labelformat = ""


        self.data = None
        self.axes = {'data': [], 'conc': [], 'resid': [], 'fit': None}

        # self.__setup_popupmenu()
        # cid = self.mpl_connect('button_press_event', self.on_click)

    @QtCore.pyqtSlot(float)
    def label_size_changed(self, new_size):
        """Slot for when label size changes"""
        ...

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
        ...

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
        ...

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
        ...

    def plot_ionic(self, **kwargs):
        ...

    def plot_speciation(self, widget, externaldata=None, **plotoptions):
        """Plot Titration simulation or species ditribuction.

        Parameters:
            widget (:class:`DataWidget`): the widget to plot
            externaldata (:class:`ExternalData`, optional): 
        """
        x_value = widget.xdata()
        y_value = widget.ydata()
        species_labels = widget.plot_labels(format_='html')
        self.clear()
        plot = self.addPlot(row=1, col=1)
        plot.setLabel("left",  widget.ylabel())
        plot.setLabel("bottom", widget.xlabel())
        for y in y_value.T:
            plot.plot(x_value, y)

        if plotoptions.get('plot_labels', True):
            window = plot.viewRange()
            max_y_value = window[1][1]
            label_xy = libplot.position_labels(x_value, y_value, window)
            labels_cutoff = plotoptions.get('unlabel_lower_than', 5.0)*max_y_value/100
            self.place_labels(plot, species_labels, label_xy,
                              threshold=labels_cutoff)

    def plotSpectr(self, widget):
        """Plot all the information for a spectrometry fit.

        It takes the data available from the internal *data* variable and plots
        it. The procedure is that it divides the figure in four parts. On top
        the experimental data, corresponding fit and calculated
        concentrations are plotted. In the middle, the residuals are plotted.
        On the bottom left, the fit spectra are shown and on bottom right
        variation of the residuals is displayed.
        """
        ...

    def draw_fit_axes(self):
        ...

    def place_labels(self, axes, label_text, label_xy, threshold=0.0):
        """Create labels and place them.

        Parameters:
            axes (:class:`matplotlib.Axes.axes`): the axes to place the labels into.
            label_text (iterable): the text for the labels
            label_xy (iterable): the position where the labels should be.
            threshold (float): if the maximum point is lower than this, no
                label is placed.
        """
        for xy, l in zip(label_xy, label_text):
            if xy[1] < threshold:
                continue
            text = pyqtgraph.TextItem()
            text.setHtml(l)
            axes.addItem(text)
            text.setPos(*xy)

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
        ...

    def setStyle(self, s):
        """Sets the type of style to be used for plotting

        Parameters:
            s (int): a number identifying the type of style to be used
        """
        ...

    def setColor(self, c):
        """Sets the type of color to be used for plotting

        Parameters:
            c (int): a number identifying the type of color scheme to be used

        Raises:
            ValueError: if the value of c is outside acceptable range
        """
        ...


    def mpl_connect(*args, **kwargs):
        pass


def smooth_curve(xmin, xmax, ydata):            # TODO move somewhere else
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

