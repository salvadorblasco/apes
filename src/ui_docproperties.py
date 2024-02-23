# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/docproperties.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Properties(object):
    def setupUi(self, Properties):
        Properties.setObjectName("Properties")
        Properties.resize(400, 300)
        self.verticalLayout = QtWidgets.QVBoxLayout(Properties)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(Properties)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.le_author = QtWidgets.QLineEdit(Properties)
        self.le_author.setObjectName("le_author")
        self.gridLayout.addWidget(self.le_author, 0, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(Properties)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.le_title = QtWidgets.QLineEdit(Properties)
        self.le_title.setObjectName("le_title")
        self.gridLayout.addWidget(self.le_title, 1, 1, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.label_3 = QtWidgets.QLabel(Properties)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.pte_comment = QtWidgets.QPlainTextEdit(Properties)
        self.pte_comment.setObjectName("pte_comment")
        self.verticalLayout.addWidget(self.pte_comment)
        self.buttonBox = QtWidgets.QDialogButtonBox(Properties)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Properties)
        self.buttonBox.accepted.connect(Properties.accept)
        self.buttonBox.rejected.connect(Properties.reject)
        QtCore.QMetaObject.connectSlotsByName(Properties)

    def retranslateUi(self, Properties):
        _translate = QtCore.QCoreApplication.translate
        Properties.setWindowTitle(_translate("Properties", "Dialog"))
        self.label.setText(_translate("Properties", "Author"))
        self.label_2.setText(_translate("Properties", "Title"))
        self.label_3.setText(_translate("Properties", "Comment"))
