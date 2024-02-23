#!/usr/bin/python3

import sys

import unittest
from PyQt5 import QtWidgets

sys.path.append('../src/')
import libqt

app = QtWidgets.QApplication(sys.argv)


class TableTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.widget = QtWidgets.QWidget()
        layout = QtWidgets.QFormLayout()
        self.qtable = QtWidgets.QTableWidget()
        layout.addWidget(self.qtable)
        self.widget.setLayout(layout)

    def tearDown(self):
        self.widget.close()

    def __fill_table_with(self, fill):
        for row in range(self.qtable.rowCount()):
            for col in range(self.qtable.columnCount()):
                if isinstance(fill, str):
                    qwa = QtWidgets.QTableWidgetItem(fill)
                else:
                    qwa = QtWidgets.QTableWidgetItem(next(fill))
                self.qtable.setItem(row, col, qwa)

    def __iter_all_cells(self):
        for row in range(self.qtable.rowCount()):
            for col in range(self.qtable.columnCount()):
                yield self.qtable.item(row, col)

    def __pattern1(self):
        return ("%d,%d" % (r, c)
                for r in range(self.qtable.rowCount())
                for c in range(self.qtable.columnCount()))

    def __random_cells(self, number):
        import random
        nrows = self.qtable.rowCount()
        ncols = self.qtable.columnCount()
        for n in range(number):
            r = random.randrange(nrows)
            c = random.randrange(ncols)
            qtwi = self.qtable.item(r, c)
            yield qtwi

    def test_tab2array(self):
        column_count = 5
        row_count = 5
        pattern = 'a'
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        self.__fill_table_with(pattern)

        for i in libqt.tab2array(self.qtable, range(row_count),
                                 range(column_count), dtype=str):
            for j in i:
                self.assertEqual(j, pattern)

    def test_array2tab(self):
        array = self.__pattern1()
        libqt.array2tab(self.qtable, array)
        array = self.__pattern1()
        for r in range(self.qtable.rowCount()):
            for c in range(self.qtable.columnCount()):
                qtwi = self.qtable.item(r, c)
                self.assertEqual(qtwi.text(), next(array))

    def test_fill_column(self):
        column_count = 5
        row_count = 5
        pattern = '1.0'
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        data = row_count * (pattern, )
        for col in range(column_count):
            # TODO add flags and row0
            libqt.fill_column(self.qtable, col, data)
            for row in range(row_count):
                t = self.qtable.item(row, col).text()
                self.assertEqual(t, pattern)

    # def test_fill_column_comboindex(self):
    #     fill_column_comboindex(table, indices, comboitems, col, row0=0)

    def test_fill_row(self):
        column_count = 5
        row_count = 5
        pattern = '1.0'
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        data = column_count * (pattern, )
        for row in range(row_count):
            # TODO add flags and col0
            libqt.fill_row(self.qtable, row, data)
            for col in range(column_count):
                t = self.qtable.item(row, col).text()
                self.assertEqual(t, pattern)

    # def test_fill_row_comboindex(self):
    #     fill_row_comboindex(table, indices, comboitems, row, col0=0)

    def test_cross(self):
        column_count = 5
        row_count = 5
        n_test_cells = 10
        pattern = 'a'
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        self.__fill_table_with(pattern)
        cells = tuple(self.__random_cells(n_test_cells))

        libqt.cross(cells, to_cross=True)

        for r in range(row_count):
            for c in range(column_count):
                cell = self.qtable.item(r, c)
                if cell in cells:
                    self.assertTrue(iscrossed(cell))
                else:
                    self.assertFalse(iscrossed(cell))

    def test_count_crossed(self):
        column_count = 5
        row_count = 5
        n_test_cells = 10
        pattern = 'a'
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        self.__fill_table_with(pattern)
        cells = tuple(self.__random_cells(n_test_cells))
        assert len(cells) == n_test_cells

        libqt.cross(cells, to_cross=True)

        counted = libqt.count_crossed(cells, invert=False)
        self.assertEqual(counted, n_test_cells)

        # TODO broken. It must iterate over all cells in the table
        # counted = libqt.count_crossed(cells, invert=True)
        # self.assertEqual(counted, column_count*row_count - n_test_cells)

    def test_bool_crossed(self):
        column_count = 5
        row_count = 5
        n_test_cells = 10
        pattern = 'a'
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        self.__fill_table_with(pattern)
        cells = tuple(self.__random_cells(n_test_cells))

        libqt.cross(cells)

        zzip = zip(self.__iter_all_cells(),
                   libqt.bool_crossed(self.__iter_all_cells()))
        for cell, crossed in zzip:
            if crossed:
                self.assertTrue(iscrossed(cell))
            else:
                self.assertFalse(iscrossed(cell))

    def test_indices_crossed(self):
        n_test_cells = 20
        cells = [QtWidgets.QTableWidgetItem("a") for _ in range(n_test_cells)]
        n_crosses = 10
        import random
        crossed_cells = random.sample(cells, n_crosses)
        libqt.cross(cells)

        indices = libqt.indices_crossed(cells)
        for n, cell in enumerate(cells):
            if iscrossed(cell):
                self.assertEqual(n, next(indices))

    def test_filter_crossed(self):
        column_count = 5
        row_count = 5
        n_test_cells = 10
        pattern = 'a'
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        self.__fill_table_with(pattern)
        cells = tuple(self.__random_cells(n_test_cells))

        libqt.cross(cells)

        items = self.__iter_all_cells()
        for cell in libqt.filter_crossed(items, invert=False):
            self.assertTrue(iscrossed(cell))

        items = self.__iter_all_cells()
        for cell in libqt.filter_crossed(items, invert=True):
            self.assertFalse(iscrossed(cell))

    # def test_filter_tabwidget(self):
    #     pass

    def test_iter_column(self):
        column_count = 5
        row_count = 5
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        self.__fill_table_with(self.__pattern1())

        for col in range(column_count):
            it = libqt.iter_column(self.qtable, col)
            for row, i in enumerate(it):
                self.assertEqual(i.text(), "%d,%d" % (row, col))

    # def test_iter_column_comboindex(self):
    #     iter_column_comboindex(table, col, row0=0)

    def test_iter_column_text(self):
        column_count = 5
        row_count = 5
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        self.__fill_table_with(self.__pattern1())

        for col in range(column_count):
            it = libqt.iter_column_text(self.qtable, col)
            for row, i in enumerate(it):
                self.assertEqual(i, "%d,%d" % (row, col))

    # def test_iter_column_widget(self):
    #     iter_column_widget(table, col, row0=0)

    def test_iter_row(self):
        column_count = 5
        row_count = 5
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        self.__fill_table_with(self.__pattern1())

        for row in range(row_count):
            it = libqt.iter_row(self.qtable, row)
            for col, i in enumerate(it):
                self.assertEqual(i.text(), "%d,%d" % (row, col))

    def test_iter_row_text(self):
        column_count = 5
        row_count = 5
        self.qtable.setColumnCount(column_count)
        self.qtable.setRowCount(row_count)
        self.__fill_table_with(self.__pattern1())

        for row in range(row_count):
            it = libqt.iter_row_text(self.qtable, row)
            for col, i in enumerate(it):
                self.assertEqual(i, "%d,%d" % (row, col))

    # def test_iter_row_widget(self):
    #     iter_row_widget(table, row, col0=0)
    # def test_iter_row_comboindex(self):
    #     iter_row_comboindex(table, row, col0=0)
    # def test_iter_tabwidget(self):
    #     iter_tabwidget(tabwidget)
    # def test_remove_curr_row(self):
    #     remove_curr_row(tab)
    # def test_selected_rows(self):
    #     selected_rows(table)
    # def test_confirm_delete(self):
    #     confirm_delete(parent, msg)
    # def test_create_combo(self):
    #     create_combo(options, default_option=0)
    # def test_opendialog(self):
    #     opendialog(parent, filters, directory='.')
    # def test_popwarning(self):
    #     popwarning(text, more_info)


def iscrossed(cell):
    return cell.font().strikeOut()


if __name__ == '__main__':
    unittest.main()
