"Auxiliary routines required for plotting"

from PyQt5 import QtWidgets as QtGui
from PyQt5 import QtCore

from ui_linepropertiesdialog import Ui_LinePropertiesDialog
import matplotlib.text

import libaux


class LinePropertiesDialog(QtGui.QDialog):
    def __init__(self, parent):
        QtGui.QDialog.__init__(self, parent=parent)
        self.setModal(True)
        self.ui = Ui_LinePropertiesDialog()
        self.ui.setupUi(self)
        self.show()

        self.__linetypes = {'solid line style': '-', 'dashed line style': '--',
                            'dash-dot line style': '-.',
                            'dotted line style': ':', 'point marker': '.',
                            'pixel marker': ',', 'circle marker': 'o',
                            'triangle_down marker': 'v',
                            'triangle_up marker': '^',
                            'triangle_left marker': '<',
                            'triangle_right marker': '>',
                            'tri_down marker': '1',
                            'tri_up marker': '2', 'tri_left marker': '3',
                            'tri_right marker': '4', 'square marker': 's',
                            'pentagon marker': 'p', 'star marker': '*',
                            'hexagon1 marker': 'h', 'hexagon2 marker': 'H',
                            'plus marker': '+', 'x marker': 'x',
                            'diamond marker': 'D', 'thin_diamond marker': 'd',
                            'vline marker': '|', 'hline marker': '_'}
        self.__reverselinetypes = dict(zip(self.__linetypes.values(),
                                           self.__linetypes.keys()))

        self.ui.cb_linetype.addItems([str(x) for x in self.__linetypes.keys()])

        self.__colors = {'blue': 'b', 'green': 'g', 'red': 'r', 'cyan': 'c',
                         'magenta': 'm', 'yellow': 'y', 'black': 'k',
                         'white': 'w'}
        self.__reversecolors = dict(zip(self.__colors.values(),
                                        self.__colors.keys()))

        self.ui.cb_color.addItems([str(x) for x in self.__colors.keys()])

    def set_initial(self, color, linestyle, linewidth):
        assert color in self.__colors.values()
        assert linestyle in self.__linetypes.values()
        assert isinstance(linewidth, float)

        self.ui.sb_thickness.setValue(linewidth)

        self.ui.cb_color.setCurrentIndex(
            self.ui.cb_color.findText(
                self.__reversecolors[color]))

        self.ui.cb_linetype.setCurrentIndex(
            self.ui.cb_linetype.findText(
                self.__reverselinetypes[linestyle]))

    def get_result(self):
        color = self.__colors[str(self.ui.cb_color.currentText())]
        ls = self.__linetypes[str(self.ui.cb_linetype.currentText())]
        lw = self.ui.sb_thickness.value()

        return {'color': color, 'linestyle': ls, 'linewidth': lw}


class DragableLabel(QtCore.QObject):
    labelsChanged = QtCore.pyqtSignal(matplotlib.text.Text)

    def __init__(self, label):
        super().__init__()
        self.label = label
        self.press = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.label.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.label.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.label.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        """Press catcher.

        On button press we will see if the mouse is over us and store some
        data.
        """
        if event.inaxes != self.label.axes:
            return

        contains, attrd = self.label.contains(event)
        if not contains:
            return

        if event.button == 1:
            x0, y0 = self.label.get_position()
            self.press = x0, y0, event.xdata, event.ydata
        elif event.button == 3:
            self.labelsChanged.emit(self.label)

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None:
            return
        if event.inaxes != self.label.axes:
            return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.label.set_x(x0+dx)
        self.label.set_y(y0+dy)

        self.label.figure.canvas.draw()

    def on_release(self, event):
        'on release we reset the press data'
        self.press = None
        self.label.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.label.figure.canvas.mpl_disconnect(self.cidpress)
        self.label.figure.canvas.mpl_disconnect(self.cidrelease)
        self.label.figure.canvas.mpl_disconnect(self.cidmotion)


def make_labels(label_seed, P):
    '''This routine takes the names of the components (label_seed) and the matrix
    of stoichiometric coefficients and makes a list of species in equilibrium.
    For instance: if primary species are 'L' and 'H' and 'L' can take up to
    three protons ('H'), the species in equilibrium are (including H-1): 'L',
    'H', 'HL', 'H_2L', 'H_3L' and 'H-1' and that is the list returned.

    Parameters:
        labels_seed (list): list of :class:`str` with the bases for making
            the labels.
        P (:class:`ndarray`): the stoichiometric coefficients.

    Returns:
        A list of strings with the constructed labels.
    '''

    # libaux.islist(label_seed)
    labels = []
    labels.extend(label_seed)

    for row in range(P.shape[0]):
        t = ""
        for col in range(P.shape[1]):
            if P[row, col] == 0:
                pass
            elif P[row, col] == 1:
                t += label_seed[col]
            else:
                t += label_seed[col]
                t += "${}_{" + str(P[row, col]) + "}$"
        labels.append(t)

    return labels


def position_labels(x, y, window):
    """Returns a list of tuples containing the optimized positions for labels
    for the given data.

    Arguments:
        x (:class:`numpy.ndarray`): 1D array containing the independent
            variable
        y (:class:`numpy.ndarray`): 2D array containing the data for the
            dependent variable
        window (:class:`tuple`): tuple containing the limits of the plot
            ((x_min, x_max), (y_min, y_max))

    Returns:
        A list of 2-sized tuples containing the (x, y) coordinates calculated
            for the labels.
    """

    libaux.assert_array_dim(1, x)
    libaux.assert_array_dim(2, y)
    assert isinstance(window, tuple)
    assert len(window) == 2 and len(window[0]) == 2 and len(window[1]) == 2

    x_min, x_max = window[0]
    y_min, y_max = window[1]
    x_window = x_max - x_min
    y_window = y_max - y_min

    # i. find maxima
    i_max = y.argmax(0)
    ret = [(x[i], y[i, n]) for n, i in enumerate(i_max)]

    # ii. prevent collisions
    # TODO
    min_dx = 0.05*x_window
    min_dy = 0.05*y_window
    for i, (x1, y1) in enumerate(ret):
        for j, (x2, y2) in enumerate(ret[i:]):
            dx = x2-x1
            dy = y2-y1
            if dx < min_dx and dy < min_dy:
                # overlapping detected
                # DO SOMETHING
                pass

    # iii. prevent outplacements
    left_margin = right_margin = 0.1 * x_window
    top_margin = bottom_margin = 0.1 * y_window
    for i, (x, y) in enumerate(ret):
        if x < (window[0][0] + left_margin):
            x = window[0][0] + left_margin
        elif x > (window[0][1] - right_margin):
            x = window[0][1] - right_margin

        if y < (window[1][0] + bottom_margin):
            y = window[1][0] + bottom_margin
        elif y > (window[1][1] - top_margin):
            y = window[1][1] - top_margin

        ret[i] = (x, y)

    return ret


def evenscale(scale, *axes):
    """Given a collections of axes, sets the same x and/or y scale for all

    Parameters:
        scale (str): The dimmension to act upon. Accepted values are 'x', 'y'
            or 'xy'.
        *axes (:class:`matplotlib.Axes`): The axes to act upon

    Raises:
        ValueError: if scale does not conform with the accepted values.
    """
    if scale not in 'xy':
        raise ValueError("%s is not valid" % scale)

    if 'x' in scale:
        xlims = [b for a in axes for b in a.get_xlim()]
        for a in axes:
            a.set_xlim((min(xlims), max(xlims)))
    if 'y' in scale:
        ylims = [b for a in axes for b in a.get_ylim()]
        for a in axes:
            a.set_ylim((min(ylims), max(ylims)))
    return
