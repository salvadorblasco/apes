#!/usr/bin/python3

import unittest

import sys
sys.path.append('../src')
import libio


class TestImport(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def checkis(self, what, *args):
        self.assertIsInstance(what, args[0])
        if len(args) > 1:
            for w in what:
                self.checkis(w, args[1:])

    @unittest.skip("Not implemented")
    def test_hyperquad(self):
        pass

    def test_k88(self):
        for _file in lsdir('*.k88'):
            A, B, V1, V2, idp, Q = libio.import_K88_data(_file)
            self.checkis(A, list, list, float)
            self.checkis(B, list, list, float)
            self.checkis(V1, list, float)
            self.checkis(V2, list, float)
            self.checkis(idp, list, str)
            self.checkis(Q, list, float)

    def test_pasat(self):
        for _file in lsdir('*.ptr'):
            data = libio.import_pasat_data(_file)
            self.checkis(data, tuple, list, float)

    @unittest.skip("Not implemented")
    def test_spec(self):
        pass

    def test_superquad(self):
        checkis = self.checkis
        for _file in lsdir('*.sup'):
            stream = libio.import_superquad_data(_file)
            checkis(next(stream), str)
            checkis(next(stream), tuple, int)
            checkis(next(stream), tuple, str)
            checkis(next(stream), float)
            checkis(next(stream), list, float)
            checkis(next(stream), list, list, int)
            checkis(next(stream), list, int)
            
            for s in stream:
                amm, elc, dat = s

                plot_keys, order, t0, buret, tflag = amm
                checkis(plot_keys, list, int)
                checkis(order, list, int)
                checkis(tflag, list, int)
                checkis(t0, list, float)
                checkis(buret, list, float)
                    
                V, errV, n, hindex, emf0, erremf0 = elc
                self.assertIsInstance(V, float)
                self.assertIsInstance(errV, float)
                self.assertIsInstance(n, int)
                self.assertIsInstance(hindex, int)
                self.assertIsInstance(emf0, float)
                self.assertIsInstance(erremf0, float)

                checkis(dat, tuple, tuple, float)


# def seqof(typ, seq):
#     for s in seq:
#         if not isinstance(s, typ):
#             return False
#     return True


def lsdir(pattern, folder='.'):
    import fnmatch
    import os

    for _file in os.listdir('.'):
        if fnmatch.fnmatch(_file, pattern):
            yield _file


if __name__ == '__main__':
    unittest.main()
