# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/import.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_ImportDialog(object):
    def setupUi(self, ImportDialog):
        ImportDialog.setObjectName("ImportDialog")
        ImportDialog.resize(281, 375)
        self.verticalLayout = QtWidgets.QVBoxLayout(ImportDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(ImportDialog)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.le_filename = QtWidgets.QLineEdit(ImportDialog)
        self.le_filename.setObjectName("le_filename")
        self.horizontalLayout.addWidget(self.le_filename)
        self.tb_select = QtWidgets.QToolButton(ImportDialog)
        self.tb_select.setObjectName("tb_select")
        self.horizontalLayout.addWidget(self.tb_select)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_2 = QtWidgets.QLabel(ImportDialog)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_2.addWidget(self.label_2)
        self.cb_filetype = QtWidgets.QComboBox(ImportDialog)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cb_filetype.sizePolicy().hasHeightForWidth())
        self.cb_filetype.setSizePolicy(sizePolicy)
        self.cb_filetype.setObjectName("cb_filetype")
        self.cb_filetype.addItem("")
        self.cb_filetype.addItem("")
        self.cb_filetype.addItem("")
        self.cb_filetype.addItem("")
        self.cb_filetype.addItem("")
        self.horizontalLayout_2.addWidget(self.cb_filetype)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.pushButton = QtWidgets.QPushButton(ImportDialog)
        self.pushButton.setObjectName("pushButton")
        self.horizontalLayout_3.addWidget(self.pushButton)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.label_3 = QtWidgets.QLabel(ImportDialog)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.lw_contents = QtWidgets.QListWidget(ImportDialog)
        self.lw_contents.setObjectName("lw_contents")
        self.verticalLayout.addWidget(self.lw_contents)
        self.cb_emptyproject = QtWidgets.QCheckBox(ImportDialog)
        self.cb_emptyproject.setObjectName("cb_emptyproject")
        self.verticalLayout.addWidget(self.cb_emptyproject)
        self.buttonBox = QtWidgets.QDialogButtonBox(ImportDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(ImportDialog)
        self.buttonBox.accepted.connect(ImportDialog.accept)
        self.buttonBox.rejected.connect(ImportDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(ImportDialog)

    def retranslateUi(self, ImportDialog):
        _translate = QtCore.QCoreApplication.translate
        ImportDialog.setWindowTitle(_translate("ImportDialog", "Dialog"))
        self.label.setText(_translate("ImportDialog", "File:"))
        self.tb_select.setText(_translate("ImportDialog", "..."))
        self.label_2.setText(_translate("ImportDialog", "Type:"))
        self.cb_filetype.setItemText(0, _translate("ImportDialog", "Autodetect"))
        self.cb_filetype.setItemText(1, _translate("ImportDialog", "Superquad"))
        self.cb_filetype.setItemText(2, _translate("ImportDialog", "Hyperquad"))
        self.cb_filetype.setItemText(3, _translate("ImportDialog", "PASAT"))
        self.cb_filetype.setItemText(4, _translate("ImportDialog", "k88"))
        self.pushButton.setText(_translate("ImportDialog", "Import"))
        self.label_3.setText(_translate("ImportDialog", "Nothing to import"))
        self.cb_emptyproject.setText(_translate("ImportDialog", "Empty current project"))
