#!/usr/bin/python3
# -*- encoding: utf-8 -*-

"""
Application starter.

.. module:: apes.py
.. moduleauthor:: Salvador Blasco <salvador.blasco@gmail.com>

This module runs the main application.
"""

import logging
import sys
from typing import Final

from PyQt5 import QtWidgets, QtCore

from mainwindow import MainWindow


def main():
    """Run application."""
    app = QtWidgets.QApplication(sys.argv)
    mainwindow = MainWindow()
    mainwindow.show()
    QtCore.pyqtRemoveInputHook()
    sys.exit(app.exec_())


if __name__ == '__main__':
    logging_level: Final[int] = logging.DEBUG     # for development
    # logging_level: Final[int] = logging.INFO      # for production
    logging.basicConfig(filename='apes.log',
                        filemode='w',
                        level=logging.DEBUG,
                        format='%(asctime)s %(message)s')
    logging.info('Logging started')
    main()
    logging.info('Logging ended')
