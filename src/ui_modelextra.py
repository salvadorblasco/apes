# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/modelextra.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ModelExtra(object):
    def setupUi(self, ModelExtra):
        ModelExtra.setObjectName("ModelExtra")
        ModelExtra.resize(555, 368)
        self.verticalLayout = QtWidgets.QVBoxLayout(ModelExtra)
        self.verticalLayout.setObjectName("verticalLayout")
        self.table_emodel = QtWidgets.QTableWidget(ModelExtra)
        self.table_emodel.setObjectName("table_emodel")
        self.table_emodel.setColumnCount(0)
        self.table_emodel.setRowCount(0)
        self.verticalLayout.addWidget(self.table_emodel)

        self.retranslateUi(ModelExtra)
        QtCore.QMetaObject.connectSlotsByName(ModelExtra)

    def retranslateUi(self, ModelExtra):
        _translate = QtCore.QCoreApplication.translate
        ModelExtra.setWindowTitle(_translate("ModelExtra", "Form"))

