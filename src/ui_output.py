# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/output.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_OutputWidget(object):
    def setupUi(self, OutputWidget):
        OutputWidget.setObjectName("OutputWidget")
        OutputWidget.resize(479, 359)
        self.verticalLayout = QtWidgets.QVBoxLayout(OutputWidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.textBrowser = QtWidgets.QTextBrowser(OutputWidget)
        font = QtGui.QFont()
        font.setFamily("Courier New")
        self.textBrowser.setFont(font)
        self.textBrowser.setObjectName("textBrowser")
        self.verticalLayout.addWidget(self.textBrowser)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.pb_update = QtWidgets.QPushButton(OutputWidget)
        self.pb_update.setObjectName("pb_update")
        self.horizontalLayout.addWidget(self.pb_update)
        self.pb_clear = QtWidgets.QPushButton(OutputWidget)
        self.pb_clear.setObjectName("pb_clear")
        self.horizontalLayout.addWidget(self.pb_clear)
        self.verticalLayout.addLayout(self.horizontalLayout)

        self.retranslateUi(OutputWidget)
        QtCore.QMetaObject.connectSlotsByName(OutputWidget)

    def retranslateUi(self, OutputWidget):
        _translate = QtCore.QCoreApplication.translate
        OutputWidget.setWindowTitle(_translate("OutputWidget", "Form"))
        self.textBrowser.setHtml(_translate("OutputWidget", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Courier New\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Noto Sans\'; font-weight:600;\">APES</span><span style=\" font-family:\'Noto Sans\';\">, the All-Purpose Equilibrium Solver</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Noto Sans\';\">(c) 2016 Salvador Blasco &lt;salvador.blasco@gmail.com&gt;</span></p>\n"
"<hr /></body></html>"))
        self.pb_update.setText(_translate("OutputWidget", "Update"))
        self.pb_clear.setText(_translate("OutputWidget", "Clear"))
