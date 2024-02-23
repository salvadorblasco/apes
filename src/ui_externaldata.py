# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/externaldata.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ExternalData(object):
    def setupUi(self, ExternalData):
        ExternalData.setObjectName("ExternalData")
        ExternalData.resize(400, 300)
        self.verticalLayout = QtWidgets.QVBoxLayout(ExternalData)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(ExternalData)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.le_title = QtWidgets.QLineEdit(ExternalData)
        self.le_title.setObjectName("le_title")
        self.horizontalLayout.addWidget(self.le_title)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.table = QtWidgets.QTableWidget(ExternalData)
        self.table.setObjectName("table")
        self.table.setColumnCount(2)
        self.table.setRowCount(2)
        item = QtWidgets.QTableWidgetItem()
        self.table.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.table.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table.setHorizontalHeaderItem(1, item)
        self.verticalLayout.addWidget(self.table)

        self.retranslateUi(ExternalData)
        QtCore.QMetaObject.connectSlotsByName(ExternalData)

    def retranslateUi(self, ExternalData):
        _translate = QtCore.QCoreApplication.translate
        ExternalData.setWindowTitle(_translate("ExternalData", "Form"))
        self.label.setText(_translate("ExternalData", "Title:"))
        item = self.table.verticalHeaderItem(0)
        item.setText(_translate("ExternalData", "label"))
        item = self.table.verticalHeaderItem(1)
        item.setText(_translate("ExternalData", "1"))
        item = self.table.horizontalHeaderItem(0)
        item.setText(_translate("ExternalData", "X"))
        item = self.table.horizontalHeaderItem(1)
        item.setText(_translate("ExternalData", "Y"))

