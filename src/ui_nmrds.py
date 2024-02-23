# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/nmrds.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_NmrWidget(object):
    def setupUi(self, NmrWidget):
        NmrWidget.setObjectName("NmrWidget")
        NmrWidget.resize(611, 422)
        self.gridLayout = QtWidgets.QGridLayout(NmrWidget)
        self.gridLayout.setObjectName("gridLayout")
        self.table = QtWidgets.QTableWidget(NmrWidget)
        self.table.setObjectName("table")
        self.table.setColumnCount(2)
        self.table.setRowCount(4)
        item = QtWidgets.QTableWidgetItem()
        self.table.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.table.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.table.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.table.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table.setHorizontalHeaderItem(1, item)
        self.table.horizontalHeader().setDefaultSectionSize(80)
        self.gridLayout.addWidget(self.table, 0, 0, 1, 1)

        self.retranslateUi(NmrWidget)
        QtCore.QMetaObject.connectSlotsByName(NmrWidget)

    def retranslateUi(self, NmrWidget):
        _translate = QtCore.QCoreApplication.translate
        NmrWidget.setWindowTitle(_translate("NmrWidget", "Form"))
        item = self.table.verticalHeaderItem(0)
        item.setText(_translate("NmrWidget", "T(L)"))
        item = self.table.verticalHeaderItem(1)
        item.setText(_translate("NmrWidget", "T(H)"))
        item = self.table.verticalHeaderItem(2)
        item.setText(_translate("NmrWidget", "1"))
        item = self.table.verticalHeaderItem(3)
        item.setText(_translate("NmrWidget", "2"))
        item = self.table.horizontalHeaderItem(0)
        item.setText(_translate("NmrWidget", "1"))
        item = self.table.horizontalHeaderItem(1)
        item.setText(_translate("NmrWidget", "2"))

