"""Commmon tasks for Qt widgets.

This module contains some common routines to simplify interaction with
Qt libraries.
"""

from contextlib import contextmanager
from collections.abc import Iterable

from PyQt5 import QtCore
from PyQt5 import QtWidgets


# +---------------------------------------------------------------------------+
# | Functions related to QTableWidget and QTableWidgetItem                    |
# +---------------------------------------------------------------------------+


def tab2array(tab: QtWidgets.QTableWidget, rows: Iterable[int], cols: Iterable[int],
              dtype=float, masked=True) -> tuple:
    '''Write 1D or 2D array to :class:`PyQt4.QtWidgets.QTableWidget`.

    Given a :class:`PyQt4.QtWidgets.QTableWidget` and a range of rows and columns,
    the data in that area of the table is read and returned.

    Parameters:
        tab (:class:`PyQt4.QtWidgets.QTableWidget`): The table to read from
        rows (iterable): the rows to loop over and read
        cols (iterable): the columns to loop over and read
        dtype (callable): the type for the return array
        masked (bool): if set to True, values crossed over will not be
            returned. If False, all values will be returned.
    Returns:
        tuple of tuple: data to read from the table
    Raises:
        IndexError: if any item in **rows** or **cols** is out of bounds.
    '''
    def dump_col(row):
        "dummy function."
        it1 = (tab.item(row, col) for col in cols)
        if masked:
            it2 = filter_crossed(it1, invert=True)
        else:
            it2 = it1
        return tuple(dtype(i.text()) for i in it2)

    return tuple(dump_col(row) for row in rows)


def tab2maarray(tab, rows, cols, dtype=float):
    '''Write contents of table to a masked array.

    Given a :class:`PyQt4.QtWidgets.QTableWidget` and a range of rows and columns,
    the data in that area of the table is read and returned.

    Parameters:
        tab (:class:`PyQt4.QtWidgets.QTableWidget`): The table to read from
        rows (iterable): the rows to loop over and read
        cols (iterable): the columns to loop over and read
        dtype (callable): the type for the return array
    Returns:
        tuple: where the first element is the data read from the table and the
            second is the mask.
    '''
    def dump_col(row):
        it1 = (tab.item(row, col) for col in cols)
        return tuple(dtype(i.text()) for i in it1)

    def dump_mask(row):
        it1 = (tab.item(row, col) for col in cols)
        return tuple(bool_crossed(it1))

    return tuple(dump_col(row) for row in rows), \
        tuple(dump_mask(row) for row in rows)


def array2tab(tab: QtWidgets.QTableWidget, array: Iterable, row0:int=0, col0:int=0,
              formatting:str="{}") -> None:
    '''Write 2D array to QTableWidget.

    Given a :class:`PyQt4.QtWidgets.QTableWidget` and an iterable (or iterable of
    iterables), put the data from the array into the table.

    Parameters:
        tab (:class:`PyQt4.QtWidgets.QTableWidget`): The table to be written to
        array (sequence of sequences): The data to read from. The first level
            of the sequence will be interpreted as rows, and the second level
            as columns.
        row0 (int): index of first row to start writing data
        col0 (int): index of first column to start writing data

    .. note:: this function does not resize the table in any way, it only
       writes into cells.
    '''
    for row, rowcontents in enumerate(array, start=row0):
        fill_row(tab, row, rowcontents, col0, formatting=formatting)


def fill_column(tab: QtWidgets.QTableWidget, col:int, data: Iterable, row0:int=0, flags=None,
                formatting="{}") -> None:
    '''Fill column of a table with data.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): The table to write to.
        col (int): index of column to be filled with data
        data (iterable): data to fill column with
        flags (dict): extra flags to be passed to the newly created
            :class:`PyQt4.QtWidgets.QTableWidgetItem' object
        row0 (int, optional): rows to skip at beginning
    '''
    for row, item in enumerate(data, start=row0):
        widget = QtWidgets.QTableWidgetItem(formatting.format(item))
        if flags is not None:
            widget.setFlags(flags)
        tab.setItem(row, col, widget)


def fill_column_comboindex(table, indices, comboitems, col, row0=0):
    '''Fill a column with QComboBoxes.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): The table to write to.
        indices (sequence of int): The indices for the QComboBoxes to be
            initialized to corresponding to the item in **comboitems**.
        comboitems (sequence of str): The items to fill QComboBoxes.
        col (int): The index of the table column to loop over
        row0 (int): number of rows to skip at the beginning.
    '''
    for row, index in enumerate(indices, start=row0):
        table.setCellWidget(row, col, create_combo(comboitems, index))


def fill_row(table, row, data, col0=0, formatting="{}"):
    '''Fill row of a table with data.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): The table to write to.
        row (int): index of row to be filled with data
        data (iterable): data to fill row with
        col0 (int, optional): cols to skip at beginning
    '''
    for col, item in enumerate(data, start=col0):
        table.setItem(row, col, QtWidgets.QTableWidgetItem(formatting.format(item)))


def fill_row_comboindex(table, indices, comboitems, row, col0=0):
    '''Fill a column with QComboBoxes.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): The table to write to.
        indices (sequence of int): The indices for the QComboBoxes to be
            initialized to corresponding to the item in **comboitems**.
        comboitems (sequence of str): The items to fill QComboBoxes.
        col (int): The index of the table column to loop over
        row0 (int): number of rows to skip at the beginning.
    '''
    for col, index in enumerate(indices, start=col0):
        table.setCellWidget(row, col, create_combo(comboitems, index))


def cross(items, to_cross=None):
    """Cross over the items provided.

    Parameters:
        items (iterable): Provide a stream of
            :class:`PyQt4.QtWidgets.QTableWidgetItem` objects which will be
            crossed or uncrossed over.
        to_cross: If True, cross over the items; if False, uncross them over;
            if None, switch the value.
    Raises:
        TypeError: if any of the items is not
            :class:`PyQt4.QtWidgets.QTableWidgetItem`
    """
    for item in items:
        if not isinstance(item, QtWidgets.QTableWidgetItem):
            raise TypeError('not a QTableWidgetItem but a ', type(item))

        font = item.font()
        if to_cross is None:
            _cross = not font.strikeOut()
        else:
            _cross = to_cross and True
        font.setStrikeOut(_cross)
        item.setFont(font)


def count_crossed(items, invert=False):
    """Count the items crossed over.

    Parameters:
        items (iterable): must yield instances of
            :class:`PyQt4.QtWidgets.QTableWidgetItem` to be counted.
        invert (bool): if True, count the items not crossed instead.
    Returns:
        int: The number of (not) crossed items
    """
    crossed_items = filter_crossed(items, invert)
    return len(tuple(crossed_items))


def bool_crossed(items):
    '''Say wheter items are crossed over or not.

    Convert a stream of items into stream of bools where True means
    the item is crossed over and False if not.

    Parameters:
        items (iterable): must yield instances of
            :class:`PyQt4.QtWidgets.QTableWidgetItem` to be filtered.
    Yield:
        True for each item crossed over, False otherwise.
    '''
    for i in items:
        if not isinstance(i, QtWidgets.QTableWidgetItem):
            raise TypeError
        font = i.font()
        yield font.strikeOut()


def indices_crossed(items):
    '''List indices of items crossed over.

    Parameters:
        items (iterable): must yield instances of
            :class:`PyQt4.QtWidgets.QTableWidgetItem` to be filtered.
    Yield:
        indices **items** for which the text is not crossed over.
    '''
    for number, widget in enumerate(items):
        assert isinstance(widget, QtWidgets.QTableWidgetItem)
        font = widget.font()
        if font.strikeOut():
            yield number


def filter_crossed(items, invert=False):
    '''Filter items in table crossed over.

    Parameters:
        items (iterable): must yield instances of
            :class:`PyQt4.QtWidgets.QTableWidgetItem` to be filtered.
        invert (bool): if True, yield the items not crossed instead.
    Yield:
        items in **items** for which the text is crossed over.
    '''
    for widget in items:
        if not isinstance(widget, QtWidgets.QTableWidgetItem):
            raise TypeError('item must be QTableWidgetItem')

        font = widget.font()
        if invert:
            if not font.strikeOut():
                yield widget
        else:
            if font.strikeOut():
                yield widget


def filter_tabwidget(table, ctype):
    '''Filter widgets in a table of a specific type.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): The table to write to.
        ctype (type): The type to be filtered.
    Yield:
        The filtered objects.
    '''
    yield from filter(lambda w: isinstance(w, ctype),
                      iter_tabwidget(table))


def freeze_column(table, col: int):
    "Make whole column not editable."
    for row in range(table.rowCount()):
        item = table.item(row, col)
        item.setFlags(QtCore.Qt.ItemIsSelectable)


def iter_column(table, col, row0=0):
    '''Iterate over items in a column.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): A table to read from
        col (int): The index of the column to loop over.
        row0 (int): The row to start looping over.
    Yields:
        items in given column
    '''
    if row0 >= table.rowCount():
        raise ValueError("row0 must be smaller than column count")
    for row in range(row0, table.rowCount()):
        yield table.item(row, col)


def iter_column_comboindex(table, col, row0=0):
    '''Return indices of QComboBoxes in a column.

    Loop over a column from **row0** to the bottom, grab the
    :class:`PyQt4.QtWidgets.QComboBox` and return the currentIndex value

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): A table to read from
        col (int): The index of the column to loop over.
        row0 (int): The row to start looping over.
    Yields:
        int: the indices of the comboboxes along the column
    Raises:
        TypeError: if any cell on the way does not contain a QComboBox object.
    '''
    for combo in iter_column_widget(table, col, row0):
        index = combo.currentIndex()
        yield index


def iter_column_text(table, col, row0=0):
    '''Iterate over items text in a column.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): A table to read from
        col (int): The index of the column to loop over.
        row0 (int): The row to start looping over.
    Yields:
        text from items in given column
    '''
    yield from (i.text() for i in iter_column(table, col, row0))


def iter_column_widget(table, col, row0=0):
    '''Iterate over widgets in a column.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): A table to read from
        col (int): The index of the column to loop over.
        row0 (int): The row to start looping over.
    Yields:
        widgets in given column
    '''
    if row0 >= table.rowCount():
        raise ValueError("row0 must be smaller than row count")
    for row in range(row0, table.rowCount()):
        yield table.cellWidget(row, col)


def iter_row(table, row, **kwargs):
    '''Iterate over widgets in a row.

    Parameters:
        tab (:class:`PyQt4.QtWidgets.QTableWidget`): The table to iterate over
        row (int): index of the row to iterate over
        col0 (int, optional): the initial column
        columns (sequence, optional): the columns to iterate. If not given, all
            columns will be iterated upon starting at col0.
    Yield:
        QTableWidgetItem: for each column read
    '''
    for col in __range(table, **kwargs):
        yield table.item(row, col)


def iter_row_text(table, row, **kwargs):
    '''Iterate over items text in a row.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): A table to read from
        row (int): The index of the row to loop over.
        col0 (int): The column to start looping over.
    Yields:
        text from items in given row
    '''
    yield from (i.text() for i in iter_row(table, row, **kwargs))


def iter_row_widget(table, row, col0=0):
    '''Iterate over widgets in a row.

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): A table to read from
        row (int): The index of the row to loop over.
        col0 (int): The column to start looping over.
    Yields:
        widgets in given row
    '''
    if col0 >= table.columnCount():
        raise ValueError("col0 must be smaller than column count")
    for col in range(col0, table.columnCount()):
        yield table.cellWidget(row, col)


def iter_row_comboindex(table, row, col0=0):
    '''Index QComboBoxes in a column.

    Loop over a column from **col0** to the end, grab the
    :class:`PyQt4.QtWidgets.QComboBox` and return the currentIndex value

    Parameters:
        table (:class:`PyQt4.QtWidgets.QTableWidget`): A table to read from
        row (int): The index of the column to loop over.
        col0 (int): The row to start looping over.
    Yields:
        int: the indices of the comboboxes along the row
    Raises:
        TypeError: if any cell on the way does not contain a QComboBox object.
    '''
    for combo in iter_row_widget(table, row, col0):
        # libaux.assert_type(QtWidgets.QComboBox, combo)
        index = combo.currentIndex()
        yield index


def iter_table_widgets(table):
    '''Iterate over every widget in a table.

    Parameters:
        tab (:class:`PyQt4.QtWidgets.QTableWidget`): The table to iterate over
    Yield:
        QTableWidgetItem: for each cell
    '''
    for row in range(table.rowCount()):
        for col in range(table.columnCount()):
            yield table.item(row, col)


def iter_tabwidget(tabwidget):
    "Iterate over the widgets of a :class:`QtWidgets.QTabWidget`."
    if not isinstance(tabwidget, QtWidgets.QTabWidget):
        raise TypeError
    for i in range(tabwidget.count()):
        yield tabwidget.widget(i)


def replace_nones(table):
    "Fill empty table cells with QWidgetItem."
    for row in range(table.rowCount()):
        for col in range(table.columnCount()):
            item = table.item(row, col)
            if item is None:
                table.setItem(row, col, QtWidgets.QTableWidgetItem())


def remove_curr_row(tab):
    "Remove row currently selected."
    if tab.rowCount() < 2:
        return
    tab.removeRow(tab.currentRow())


def selected_rows(table):
    "Return a set with the rows currently selected in table."
    return set(i.row() for i in table.selectedIndexes())


@contextmanager
def signals_blocked(widget: QtWidgets.QWidget):
    """Suspend connections while changes are being made in table.

    The signals are blocked while editing table. Use as
    with libqt.signals_blocked(widget) as widget:
        # code that we do not want to trigger signals

    Parameters:
        widget (:class:`QtWidgets.QWidget`): the widget
    """
    widget.blockSignals(True)
    yield widget
    widget.blockSignals(False)


# +---------------------------------------------------------------------------+
# | Other functions                                                           |
# +---------------------------------------------------------------------------+

def checkbox_yesno(checkbox):
    "Quick creation of a yes/no QCheckBox."
    txt = {QtCore.Qt.Checked: "yes", QtCore.Qt.Unchecked: "no"}
    checkbox.setText(txt[checkbox.checkState()])


def clamp_range(position: int, number_components: int) -> int:
    """Constrain index within range.

    Parameters:
        position (int):
        number_components (int):

    Returns:
        int: the constrained index
    """
    return max(min(number_components, 1+position), 0)


def confirm_delete(parent, msg):
    'Pop up dialog to confirm option.'
    qmsg = QtWidgets.QMessageBox(parent=parent)
    qmsg.setText(msg)
    qmsg.setStandardButtons(QtWidgets.QMessageBox.Yes |
                            QtWidgets.QMessageBox.No)
    qmsg.setDefaultButton(QtWidgets.QMessageBox.Yes)
    ret = qmsg.exec()
    return ret == QtWidgets.QMessageBox.Yes


def create_combo(options, default_option=0):
    '''Create QComboBox from options and default option.

    Parameters:
        options (sequence): sequence of str with option for combo
        default_option (int): index of the defatult option
    Returns:
        PyQt4.QtWidgets.QComboBox
    '''
    new_combo = QtWidgets.QComboBox()
    new_combo.addItems(options)
    new_combo.setCurrentIndex(default_option)
    return new_combo


def opendialog(parent: QtWidgets.QWidget, filters: str, directory: str='.') -> QtWidgets.QDialog:
    """Open file dialog and return filename selected.

    Parameters:
        filters (str): filters in dialog
        directory (str): directory to open from.
    Returns:
        str: the file name selected

    .. seealso: :class:`PyQt4.QtWidgets.QFileDialog`
    """
    return QtWidgets.QFileDialog.getOpenFileName(
        parent=parent,
        caption='Choose file to open',
        directory=directory,
        filter=filters)


def popwarning(text: str, more_info: str) -> None:
    """Open messagebox dialog.

    Parameters:
        text (str): title of the box
        more_info (str): additional information

    .. seealso: :class:`PyQt4.QtWidgets.QMessageBox`
    """
    dialog = QtWidgets.QMessageBox()
    dialog.setIcon(QtWidgets.QMessageBox.Warning)
    dialog.setText(text)
    dialog.setInformativeText(more_info)
    dialog.addButton(QtWidgets.QMessageBox.Ok)
    dialog.exec_()


def __range(table, **kwargs):
    return kwargs.get('iter_range',
                      range(kwargs.get('col0', 0), table.columnCount()))
