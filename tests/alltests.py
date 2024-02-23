#!/usr/bin/python3

# tests/runner.py
import unittest

# import your test modules
import test_gui
import test_import
import test_libaux
import test_libeq
import test_libqt
import test_libio
import test_libemf

# initialize the test suite
loader = unittest.TestLoader()
suite  = unittest.TestSuite()

# add tests to the test suite
modules = (test_gui, test_import, test_libaux, test_libeq, test_libqt)
for module in modules:
    suite.addTests(loader.loadTestsFromModule(module))

# initialize a runner, pass it your suite and run it
runner = unittest.TextTestRunner(verbosity=1)
result = runner.run(suite)
