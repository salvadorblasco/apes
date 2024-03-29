# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/titration.ui'
#
# Created by: PyQt5 UI code generator 5.15.6
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Titration(object):
    def setupUi(self, Titration):
        Titration.setObjectName("Titration")
        Titration.resize(351, 188)
        self.verticalLayout = QtWidgets.QVBoxLayout(Titration)
        self.verticalLayout.setObjectName("verticalLayout")
        self.table_titration = QtWidgets.QTableWidget(Titration)
        self.table_titration.setObjectName("table_titration")
        self.table_titration.setColumnCount(5)
        self.table_titration.setRowCount(2)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setHorizontalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setItem(0, 0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setItem(0, 1, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setItem(0, 3, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setItem(1, 0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setItem(1, 1, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_titration.setItem(1, 3, item)
        self.table_titration.horizontalHeader().setDefaultSectionSize(65)
        self.table_titration.horizontalHeader().setMinimumSectionSize(20)
        self.table_titration.verticalHeader().setVisible(False)
        self.table_titration.verticalHeader().setCascadingSectionResizes(True)
        self.table_titration.verticalHeader().setDefaultSectionSize(30)
        self.table_titration.verticalHeader().setMinimumSectionSize(15)
        self.verticalLayout.addWidget(self.table_titration)
        self.label = QtWidgets.QLabel(Titration)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.lbl_points = QtWidgets.QLabel(Titration)
        self.lbl_points.setObjectName("lbl_points")
        self.gridLayout.addWidget(self.lbl_points, 0, 2, 1, 1)
        self.dsb_V0 = QtWidgets.QDoubleSpinBox(Titration)
        self.dsb_V0.setDecimals(3)
        self.dsb_V0.setMaximum(1000.0)
        self.dsb_V0.setSingleStep(0.1)
        self.dsb_V0.setProperty("value", 10.0)
        self.dsb_V0.setObjectName("dsb_V0")
        self.gridLayout.addWidget(self.dsb_V0, 2, 0, 1, 1)
        self.lbl_Vf = QtWidgets.QLabel(Titration)
        self.lbl_Vf.setObjectName("lbl_Vf")
        self.gridLayout.addWidget(self.lbl_Vf, 0, 1, 1, 1)
        self.lbl_Vi = QtWidgets.QLabel(Titration)
        self.lbl_Vi.setObjectName("lbl_Vi")
        self.gridLayout.addWidget(self.lbl_Vi, 0, 0, 1, 1)
        self.dsb_Vf = QtWidgets.QDoubleSpinBox(Titration)
        self.dsb_Vf.setDecimals(3)
        self.dsb_Vf.setMaximum(1000.0)
        self.dsb_Vf.setSingleStep(0.1)
        self.dsb_Vf.setProperty("value", 11.0)
        self.dsb_Vf.setObjectName("dsb_Vf")
        self.gridLayout.addWidget(self.dsb_Vf, 2, 1, 1, 1)
        self.sb_NPoints = QtWidgets.QSpinBox(Titration)
        self.sb_NPoints.setMinimum(10)
        self.sb_NPoints.setMaximum(1000)
        self.sb_NPoints.setSingleStep(5)
        self.sb_NPoints.setProperty("value", 100)
        self.sb_NPoints.setObjectName("sb_NPoints")
        self.gridLayout.addWidget(self.sb_NPoints, 2, 2, 1, 1)
        self.label_2 = QtWidgets.QLabel(Titration)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 3, 1, 1)
        self.dsb_Verror = QtWidgets.QDoubleSpinBox(Titration)
        self.dsb_Verror.setDecimals(3)
        self.dsb_Verror.setProperty("value", 0.003)
        self.dsb_Verror.setObjectName("dsb_Verror")
        self.gridLayout.addWidget(self.dsb_Verror, 2, 3, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)

        self.retranslateUi(Titration)
        QtCore.QMetaObject.connectSlotsByName(Titration)

    def retranslateUi(self, Titration):
        _translate = QtCore.QCoreApplication.translate
        Titration.setWindowTitle(_translate("Titration", "Form"))
        item = self.table_titration.verticalHeaderItem(0)
        item.setText(_translate("Titration", "1"))
        item = self.table_titration.verticalHeaderItem(1)
        item.setText(_translate("Titration", "2"))
        item = self.table_titration.horizontalHeaderItem(0)
        item.setText(_translate("Titration", "Label"))
        item = self.table_titration.horizontalHeaderItem(1)
        item.setText(_translate("Titration", "Initial"))
        item = self.table_titration.horizontalHeaderItem(2)
        item.setText(_translate("Titration", "Flag"))
        item = self.table_titration.horizontalHeaderItem(3)
        item.setText(_translate("Titration", "Buret"))
        item = self.table_titration.horizontalHeaderItem(4)
        item.setText(_translate("Titration", "Flag"))
        __sortingEnabled = self.table_titration.isSortingEnabled()
        self.table_titration.setSortingEnabled(False)
        item = self.table_titration.item(0, 0)
        item.setText(_translate("Titration", "L"))
        item = self.table_titration.item(0, 1)
        item.setText(_translate("Titration", "0.001"))
        item = self.table_titration.item(0, 3)
        item.setText(_translate("Titration", "0"))
        item = self.table_titration.item(1, 0)
        item.setText(_translate("Titration", "H"))
        item = self.table_titration.item(1, 1)
        item.setText(_translate("Titration", "0.005"))
        item = self.table_titration.item(1, 3)
        item.setText(_translate("Titration", "-0.1"))
        self.table_titration.setSortingEnabled(__sortingEnabled)
        self.label.setText(_translate("Titration", "Volume"))
        self.lbl_points.setText(_translate("Titration", "Points"))
        self.lbl_Vf.setText(_translate("Titration", "Final"))
        self.lbl_Vi.setText(_translate("Titration", "Initial"))
        self.label_2.setText(_translate("Titration", "Error"))
