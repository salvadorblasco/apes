#!/usr/bin/python3

import unittest

import sys
sys.path.append('../src')


class TestBridge(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        from modelwidget import ModelWidget
        model = ModelWidget()

    def test_call(self):
        ...


if __name__ == '__main__':
    unittest.main()
