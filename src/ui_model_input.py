# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/model_input.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_NewModelDialog(object):
    def setupUi(self, NewModelDialog):
        NewModelDialog.setObjectName("NewModelDialog")
        NewModelDialog.resize(264, 148)
        NewModelDialog.setModal(True)
        self.verticalLayout = QtWidgets.QVBoxLayout(NewModelDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(NewModelDialog)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.sb_components = QtWidgets.QSpinBox(NewModelDialog)
        self.sb_components.setMinimum(1)
        self.sb_components.setProperty("value", 2)
        self.sb_components.setObjectName("sb_components")
        self.gridLayout.addWidget(self.sb_components, 0, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(NewModelDialog)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.sb_equilibria = QtWidgets.QSpinBox(NewModelDialog)
        self.sb_equilibria.setMinimum(1)
        self.sb_equilibria.setProperty("value", 2)
        self.sb_equilibria.setObjectName("sb_equilibria")
        self.gridLayout.addWidget(self.sb_equilibria, 1, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(NewModelDialog)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 2, 0, 1, 1)
        self.sb_solids = QtWidgets.QSpinBox(NewModelDialog)
        self.sb_solids.setEnabled(False)
        self.sb_solids.setObjectName("sb_solids")
        self.gridLayout.addWidget(self.sb_solids, 2, 1, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.buttonBox = QtWidgets.QDialogButtonBox(NewModelDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(NewModelDialog)
        self.buttonBox.accepted.connect(NewModelDialog.accept)
        self.buttonBox.rejected.connect(NewModelDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(NewModelDialog)

    def retranslateUi(self, NewModelDialog):
        _translate = QtCore.QCoreApplication.translate
        NewModelDialog.setWindowTitle(_translate("NewModelDialog", "Parameters of new model"))
        self.label.setText(_translate("NewModelDialog", "number of components"))
        self.label_2.setText(_translate("NewModelDialog", "number of solution equilibria"))
        self.label_3.setText(_translate("NewModelDialog", "number of precipitated species"))
