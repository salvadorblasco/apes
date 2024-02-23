#!/usr/bin/python
# -*- encoding: utf-8 -*-

import sys
from PyQt4 import QtGui, QtCore
from ui_linepropertiesdialog import Ui_LinePropertiesDialog

class LinePropertiesDialog(QtGui.QDialog):
    def __init__(self):
        QtGui.QDialog.__init__(self)
        self.ui = Ui_LinePropertiesDialog()
        self.ui.setupUi(self)

    def set_colors(self, colorlist):
        pass

def main():
    app = QtGui.QApplication(sys.argv)

    mw = LinePropertiesDialog()
    print mw.result()
    
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
