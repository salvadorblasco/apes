#!/usr/bin/python3

import sys

import unittest
from PyQt5.QtGui import QApplication
# from PyQt4.QtTest import QTest

sys.path.append('../src/')
import mainwindow
import modelwidget

app = QApplication(sys.argv)

class ModelTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    # @classmethod
    # def setUpClass(cls):
    #     cls.app = None
    #     if QApplication.instance():
    #         cls.app = QApplication.instance()
    #     else:
    #         cls.app = QApplication([], False)

    # @classmethod
    # def tearDownClass(cls):
    #     del cls.app

    def setUp(self):
        # import pdb
        # pdb.set_trace()
        # self.app.exec()
        self.mainw = mainwindow.MainWindow()
        self.model = self.mainw.modelwidget

    def tearDown(self):
        self.mainw.close()
        # cls.app.quit()

    def test_addcalor(self):
        self.assertFalse(self.model.isCalori())
        self.model.addCalorimetry()
        self.assertTrue(self.model.isCalori())


if __name__ == '__main__':
    unittest.main()
