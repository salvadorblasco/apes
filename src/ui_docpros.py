# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/docpros.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_DocPropsWidget(object):
    def setupUi(self, DocPropsWidget):
        DocPropsWidget.setObjectName("DocPropsWidget")
        DocPropsWidget.resize(400, 300)
        self.verticalLayout = QtWidgets.QVBoxLayout(DocPropsWidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(DocPropsWidget)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.le_author = QtWidgets.QLineEdit(DocPropsWidget)
        self.le_author.setObjectName("le_author")
        self.horizontalLayout.addWidget(self.le_author)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_3 = QtWidgets.QLabel(DocPropsWidget)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_2.addWidget(self.label_3)
        self.dsb_temp = QtWidgets.QDoubleSpinBox(DocPropsWidget)
        self.dsb_temp.setMaximum(999.0)
        self.dsb_temp.setObjectName("dsb_temp")
        self.horizontalLayout_2.addWidget(self.dsb_temp)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.label_2 = QtWidgets.QLabel(DocPropsWidget)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.txt_comment = QtWidgets.QPlainTextEdit(DocPropsWidget)
        self.txt_comment.setObjectName("txt_comment")
        self.verticalLayout.addWidget(self.txt_comment)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.pbOk = QtWidgets.QPushButton(DocPropsWidget)
        self.pbOk.setObjectName("pbOk")
        self.horizontalLayout_3.addWidget(self.pbOk)
        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.retranslateUi(DocPropsWidget)
        QtCore.QMetaObject.connectSlotsByName(DocPropsWidget)

    def retranslateUi(self, DocPropsWidget):
        _translate = QtCore.QCoreApplication.translate
        DocPropsWidget.setWindowTitle(_translate("DocPropsWidget", "Form"))
        self.label.setText(_translate("DocPropsWidget", "Author:"))
        self.label_3.setText(_translate("DocPropsWidget", "Temperature:"))
        self.dsb_temp.setSuffix(_translate("DocPropsWidget", "°C"))
        self.label_2.setText(_translate("DocPropsWidget", "Comment"))
        self.pbOk.setText(_translate("DocPropsWidget", "OK"))

