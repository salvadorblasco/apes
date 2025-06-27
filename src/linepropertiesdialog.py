# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/linepropertiesdialog.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_LinePropertiesDialog(object):
    def setupUi(self, LinePropertiesDialog):
        LinePropertiesDialog.setObjectName(_fromUtf8("LinePropertiesDialog"))
        LinePropertiesDialog.resize(210, 336)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(LinePropertiesDialog.sizePolicy().hasHeightForWidth())
        LinePropertiesDialog.setSizePolicy(sizePolicy)
        self.verticalLayout_2 = QtGui.QVBoxLayout(LinePropertiesDialog)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.label = QtGui.QLabel(LinePropertiesDialog)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout_3.addWidget(self.label)
        self.cb_linetype = QtGui.QComboBox(LinePropertiesDialog)
        self.cb_linetype.setObjectName(_fromUtf8("cb_linetype"))
        self.horizontalLayout_3.addWidget(self.cb_linetype)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.label_2 = QtGui.QLabel(LinePropertiesDialog)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout_4.addWidget(self.label_2)
        self.sb_thickness = QtGui.QDoubleSpinBox(LinePropertiesDialog)
        self.sb_thickness.setDecimals(1)
        self.sb_thickness.setObjectName(_fromUtf8("sb_thickness"))
        self.horizontalLayout_4.addWidget(self.sb_thickness)
        self.verticalLayout_2.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.label_3 = QtGui.QLabel(LinePropertiesDialog)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.horizontalLayout_5.addWidget(self.label_3)
        self.cb_color = QtGui.QComboBox(LinePropertiesDialog)
        self.cb_color.setObjectName(_fromUtf8("cb_color"))
        self.horizontalLayout_5.addWidget(self.cb_color)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.groupBox = QtGui.QGroupBox(LinePropertiesDialog)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.verticalLayout = QtGui.QVBoxLayout(self.groupBox)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.checkBox = QtGui.QCheckBox(self.groupBox)
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.verticalLayout.addWidget(self.checkBox)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.label_4 = QtGui.QLabel(self.groupBox)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.horizontalLayout.addWidget(self.label_4)
        self.comboBox = QtGui.QComboBox(self.groupBox)
        self.comboBox.setObjectName(_fromUtf8("comboBox"))
        self.horizontalLayout.addWidget(self.comboBox)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.label_5 = QtGui.QLabel(self.groupBox)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.horizontalLayout_2.addWidget(self.label_5)
        self.comboBox_2 = QtGui.QComboBox(self.groupBox)
        self.comboBox_2.setObjectName(_fromUtf8("comboBox_2"))
        self.horizontalLayout_2.addWidget(self.comboBox_2)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.checkBox_2 = QtGui.QCheckBox(self.groupBox)
        self.checkBox_2.setObjectName(_fromUtf8("checkBox_2"))
        self.verticalLayout.addWidget(self.checkBox_2)
        self.verticalLayout_2.addWidget(self.groupBox)
        self.buttonBox = QtGui.QDialogButtonBox(LinePropertiesDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout_2.addWidget(self.buttonBox)

        self.retranslateUi(LinePropertiesDialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), LinePropertiesDialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), LinePropertiesDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(LinePropertiesDialog)

    def retranslateUi(self, LinePropertiesDialog):
        LinePropertiesDialog.setWindowTitle(_translate("LinePropertiesDialog", "Dialog", None))
        self.label.setText(_translate("LinePropertiesDialog", "Linetype", None))
        self.label_2.setText(_translate("LinePropertiesDialog", "Thickness", None))
        self.label_3.setText(_translate("LinePropertiesDialog", "Color", None))
        self.groupBox.setTitle(_translate("LinePropertiesDialog", "Error", None))
        self.checkBox.setText(_translate("LinePropertiesDialog", "Show", None))
        self.label_4.setText(_translate("LinePropertiesDialog", "Color", None))
        self.label_5.setText(_translate("LinePropertiesDialog", "Linetype", None))
        self.checkBox_2.setText(_translate("LinePropertiesDialog", "Shadow", None))

